'''
Author: zhaozijian
Date: 2022-03-03 16:16:12
LastEditors: zhaozijian
LastEditTime: 2022-03-22 13:35:55
Description: file content
'''
import sys

import numpy as np
from gefpy.cgef_reader_cy import CgefR
from gefpy.cell_exp_reader import CellExpReader

from datetime import datetime

def main():
    cgef = CgefR("/ldfssz1/ST_BI/USER/zhaozijian/geftool/build/FP200000442TL_A2.cellbin.gef", True)
    # get cells
    cells = cgef.get_cells()
    
    # get genes
    #genes = cgef.get_genes()

    #a = datetime.now()
    # 限制区域
    # cgef.restrict_region(2933,5568,3933,6568)
    # cell = cgef.get_cells()
    # print(cell)
    # cgef.restrict_region(7680,9728,8703,10752)
    # cellborder = cgef.get_cellborders()
    # print(cellborder)
    #b = datetime.now()
    #print("microseconds restrict_region: ", (b - a).microseconds)
 
    cgef.cgef_close()
    return 0

from scipy.sparse import csr_matrix
def test():
    cell_bin_gef = CellExpReader("/ldfssz1/ST_BI/USER/stereopy/test/tanliwei/test/cell_correct_result/FP200000443TL_E2.adjusted.cgef")
    exp_matrix = csr_matrix((cell_bin_gef.count, (cell_bin_gef.rows, cell_bin_gef.cols)), shape=(cell_bin_gef.cell_num, cell_bin_gef.gene_num), dtype=np.uint32)
    print(cell_bin_gef.cluster)

if __name__ == "__main__":
    test()
