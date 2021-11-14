from os import DirEntry
import h5py
from mask import Mask
from gem import Gem
from cell import Cell
import numpy as np
import cv2 as cv
import pandas as pd
import random


__cell_types__ = ['cell0', 'cell1', 'cell2', 'cell3']


def add_dataset(group, dataset_name, names, dtypes, data):
    datas = np.array(data)
    h, w = datas.shape[:2]
    ds_dt = np.dtype({'names': names,'formats': dtypes}) 
    ddd = list()
    for i in range(w):
        ddd.append(datas[:, i])
    rec_arr = np.rec.fromarrays(ddd, dtype=ds_dt)
    try:
        # normal dataset creation
        group.create_dataset(dataset_name, data=rec_arr)
    except Exception as e:
        # if dataset already exists, replace it with new data
        del group[dataset_name]
        group.create_dataset(dataset_name, data=rec_arr)


class CellExp(object):
    def __init__(self, gem_file, mask_file):
        self.cell_list = None
        self.gene_list = None
        self.cell_exp = None
        self.cell_type = None

        self.i_mask = None
        self.i_gem = None
        self._init(gem_file, mask_file)
    
    def _init(self, gem_file, mask_file):
        self.i_mask = Mask(mask_file)
        self.i_gem = Gem(gem_file, roi=self.i_mask.region)
        self.cell_list = list()
        self.cell_border_list = list()
        self.cell_type_list = list()

        for id, p in enumerate(self.i_mask.polygens):
            c_border = self.cell_border(p.border)
            c_bin = self.cell_bin(p.border)
            x, y = p.center
            count = sum(c_bin.values())
            cell_item = [x, y, p.area, len(c_bin), count]
            self.cell_list.append(cell_item)
            self.cell_border_list.append(c_border)
            self.cell_type_list.append([random.choice(__cell_types__)])
    
    def cell_border(self, border, pt_count=None):
        h, w = border.shape
        b = np.zeros((16, 2), dtype=int)
        if border is not None: 
            b[:h, :] = border
            return b
        else: 
            return None
    
    def cell_bin(self, border):
        nap = np.array(border)
        cx = nap[:, 0]
        cy = nap[:, 1]
        x0, y0, x1, y1 = [np.min(cx), np.min(cy), np.max(cx), np.max(cy)]
        w, h = [x1 - x0, y1 - y0]
        arr = np.zeros((h, w), dtype=np.uint8)
        nap[:, 0] -= x0
        nap[:, 1] -= y0
        arr = cv.fillPoly(arr, [nap], 1)
        gem_mat = arr * self.i_gem.buffer[y0: y1, x0: x1]
        y, x = np.where(gem_mat > 0)
        ind = [[x[i] + x0, y[i] + y0] for i in range(len(y))]
        return self.i_gem.get_geneExp(ind)
    
    def write(self, file_path):
        with h5py.File(file_path, 'a') as h5f:
            geneExp = h5f.require_group('geneExp')
            wholeExp = h5f.require_group('wholeExp')
            stat = h5f.require_group('stat')
            Attributes = h5f.require_group('Attributes')
            group = h5f.require_group('cellExp')

            add_dataset(group, 'cell', ['x', 'y', 'area', 'geneCount', 'expCount'], ['i4'] * 5, self.cell_list)

            ss = np.array(self.cell_border_list, dtype=int)
            
            try:
                # normal dataset creation
                group.create_dataset('cellBorder', data=np.array(self.cell_border_list, dtype=int))
            except Exception as e:
                # if dataset already exists, replace it with new data
                del group['cellBorder']
                group.create_dataset('cellBorder', data=np.array(self.cell_border_list, dtype=int))

            try:
                # normal dataset creation
                ds = group.create_dataset('cellType', data=np.array(self.cell_type_list, dtype='S32'))
            except Exception as e:
                # if dataset already exists, replace it with new data
                del group['cellType']
                ds = group.create_dataset('cellType', data=np.array(self.cell_type_list, dtype='S32'))
            group['cellType'].attrs['cellTypeCount'] = len(__cell_types__)

            # add_dataset(group, 'cellBorder', ['borderPath'], ['i4'], self.cell_border_list)





def main():
    # gem_file = '/jdfssz2/ST_BIGDATA/Stereomics/autoanalysis_backup/P20Z10200N0157/null/FP200000617TL_B6_web_2_backup/result/FP200000617TL_B6_web_2/02.alignment/GetExp/barcode_gene_exp.txt'
    # mask_file = '/zfssz3/ST_BIGDATA/stereomics/PipelineTest/data/FP200000617TL_B6/7_result/FP200000617TL_B6_mask.tif'
    gem_file = 'test_data/barCode.txt'
    mask_file = 'test_data/mask.tif'
    ce = CellExp(gem_file, mask_file)
    ce.write('./demo.gef')


if __name__ == '__main__':
    main()
        
