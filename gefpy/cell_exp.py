import h5py
from mask import Mask
import numpy as np
import pandas as pd
import cv2 as cv
from gefpy.cell_exp_cy import CellExpW

import logging
logging.basicConfig(level = logging.INFO,format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class CellExpWriterPy(object):
    cell_exp_writer: CellExpW

    def __init__(self, gef_file, mask_file, out_file):
        self.exp_matrix = None
        self.gef_file = gef_file
        self.mask_file = mask_file
        self.out_file = out_file
        self.cell_exp_writer = None

        self.x = None
        self.y = None
        self.area = None

        self.i_mask = None
        self._init()
    
    def _init(self):
        self.cell_exp_writer = CellExpW(self.out_file)
        self.cell_exp_writer.setGeneExpMap(self.gef_file)
        self.set_exp_matrix(self.gef_file)
        self.i_mask = Mask(self.mask_file)

        self.cell_num = len(self.i_mask.polygens)

        self.borders = np.zeros((self.cell_num, 16, 2), dtype=np.int8)
        self.x = np.zeros((self.cell_num,), dtype=np.uint32)
        self.y = np.zeros((self.cell_num,), dtype=np.uint32)
        self.area = np.zeros((self.cell_num, ), dtype=np.uint16)

        for id, p in enumerate(self.i_mask.polygens):
            b_len = 16 if len(p.border) >= 16 else len(p.border)
            for i in range(b_len):
                self.borders[id, i, 0] = p.border[i, 0] - p.center[0]
                self.borders[id, i, 1] = p.border[i, 1] - p.center[1]

            self.x[id] = p.center[0]
            self.y[id] = p.center[1]
            self.area[id] = p.area

            c_bin = self.cell_bin(p.border)
            self.cell_exp_writer.add_cell_bin(c_bin, len(c_bin))

    def cell_bin(self, border):
        nap = np.array(border)
        cx = nap[:, 0]
        cy = nap[:, 1]
        x0, y0, x1, y1 = [np.min(cx), np.min(cy), np.max(cx), np.max(cy)]
        w, h = [x1 - x0 + 1, y1 - y0 + 1]
        arr = np.zeros((h, w), dtype=np.uint8)
        nap[:, 0] -= x0
        nap[:, 1] -= y0

        arr = cv.fillPoly(arr, [nap], 1)

        gem_mat = arr * self.exp_matrix[y0: y1, x0: x1]
        y, x = np.where(gem_mat > 0)
        bin_index = np.zeros((len(y), 2), dtype=np.uint32)
        for i in range(len(y)):
            bin_index[i, 0] = x[i] + x0
            bin_index[i, 1] = y[i] + y0
        return bin_index

    def set_exp_matrix(self, gef_file):
        with h5py.File(gef_file, mode='r') as h5f:
            self.exp_matrix = h5f['wholeExp']['bin1']['MIDcount'].T

    def write(self):
        self.cell_exp_writer.storeGeneList()
        self.cell_exp_writer.storeCellBorder(self.borders, self.cell_num)
        self.cell_exp_writer.storeCell(self.x, self.y, self.area, self.cell_num)
        self.cell_exp_writer.storeCellExp()
        logger.info('Gef write completed : \'{}\''.format(self.out_file))


class CellExpReaderPy(object):
    def __init__(self, filepath):
        self.filepath = filepath
        self.exp_len = 0
        self.gene_num = 0
        self.cell_num = 0
        self.genes = None
        self.cells = None
        self.rows = None
        self.cols = None
        self.count = None
        self.positions = None
        self._init()

    def _init(self):
        with h5py.File(self.filepath, mode='r') as h5f:
            self.genes = pd.DataFrame(h5f['cellBin']['geneList']).values
            self.cols = pd.DataFrame(h5f['cellBin']['cellExp']['geneID']).values
            self.count = pd.DataFrame(h5f['cellBin']['cellExp']['count']).values
            self.exp_len = self.cols.shape[0]
            self.gene_num = self.genes.shape[0]
            self.positions  = pd.DataFrame(h5f['cellBin']['cell']['x', 'y']).values
            self.cell_num = self.positions.shape[0]

            self.cells = np.array(self.positions[:,0], dtype='uint64')
            self.cells = np.bitwise_or(
                np.left_shift(self.cells, 32), self.positions[:,1])

            gene_counts = pd.DataFrame(h5f['cellBin']['cell']['geneCount']).values
            self.rows = np.zeros((self.exp_len , ), dtype='uint32')
            exp_index = 0
            for i, gene_c in enumerate(gene_counts):
                for index in range(gene_c):
                    self.rows[exp_index] = i
                    exp_index += 1
