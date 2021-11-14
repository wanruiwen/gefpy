import gzip
import pandas as pd
import numpy as np
import logging
logging.basicConfig(level = logging.INFO,format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

""" modin can accelerate pandas 
17  Cr1l	10367	16260	1
200 Wdr77	10367	16260	1
"""



class Gem(object):
    def __init__(self, file_path, roi=None):
        self.header_count = 0
        self.header = None
        self.fd = None
        self.origin = None
        self.buffer = None
        self.buffer_shape = None
        self.gene_id_dct = None
        self.roi = roi
        self.useful_rows = None
        self._init(file_path)
    
    def _init(self, file_path):
        logger.info('Gene file path is: \'{}\'.'.format(file_path))
        # geneID  x       y       UMICount/MIDCount
        if file_path.endswith('.gz'): self.fd = gzip.open(file_path, 'rb')
        else: self.fd = open(file_path, 'rb')
        self.header = ''
        eoh = 0
        for index, line in enumerate(self.fd):
            line = line.decode("utf-8") # read in as binary, decode first
            if line.startswith('#'): # header lines always start with '#'
                self.header += line
                self.header_count += 1
                eoh = self.fd.tell() # get end-of-header position
            else:
                break
        logger.info('Header info: {}.'.format(self.header))

        # find start of expression matrix
        self.fd.seek(eoh)

        df = pd.read_csv(self.fd, sep='\t', header=0)

        x0, y0, x1, y1 = self.roi
        self.dct = dict()
        
        # find start of expression matrix
        self.fd.seek(eoh)
        for index, line in enumerate(self.fd):
            gene_id, x, y, count = str(line, encoding='utf-8').strip('\n').split('\t')
            try: 
                x, y = [int(x), int(y)]
                if ((x >= x0 and x <= x1) and (y >= y0 and y <= y1)): 
                    key = '{}_{}'.format(x, y)
                    if key not in self.dct: self.dct[key] = {}
                    self.dct[key][gene_id] = int(count)
            except: pass

        gene_ids = df['geneID'].unique()
        self.gene_id_dct = dict()
        self.origin = (df['x'].min(), df['y'].min())
        logger.info('Offset: {}.'.format(self.origin))
        df['x'] = df['x'] - df['x'].min()
        df['y'] = df['y'] - df['y'].min()
        max_x = df['x'].max() + 1
        max_y = df['y'].max() + 1
        logger.info('[GeneShape, geneIDCount]: {}, {}.'.format((max_y, max_x), len(gene_ids)))

        self.buffer_shape = (max_y, max_x)
        self.buffer = np.zeros(shape=self.buffer_shape, dtype=np.uint16)

        try: 
            new_df = df.groupby(['x', 'y']).agg(UMI_sum=('UMICount', 'sum')).reset_index()
        except: 
            new_df = df.groupby(['x', 'y']).agg(UMI_sum=('MIDCount', 'sum')).reset_index()
        self.buffer[new_df['y'], new_df['x']] = new_df['UMI_sum']
        
        logger.info('[ROI, geneIDCount]: {}, {}.'.format(self.roi, len(gene_ids)))

    def get_geneExp(self, locs):
        dct_ = dict()
        for loc in locs:
            key = '{}_{}'.format(loc[0], loc[1])
            for k, v in self.dct[key].items():
                if k in dct_: dct_[k] += v
                else: dct_[k] = v
        return dct_


def main():
    gem_file = '/jdfssz2/ST_BIGDATA/Stereomics/autoanalysis_backup/P20Z10200N0157/null/FP200000617TL_B6_web_2_backup/result/FP200000617TL_B6_web_2/02.alignment/GetExp/barcode_gene_exp.txt'
    gem = Gem(gem_file)
    gem.get_geneExp([100, 100])


if __name__ == '__main__':
    main()