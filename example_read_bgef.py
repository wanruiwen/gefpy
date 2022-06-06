'''
Author: zhaozijian
Date: 2022-03-03 17:25:34
LastEditors: zhaozijian
LastEditTime: 2022-05-21 14:09:12
Description: file content
'''

import sys

import numpy as np

from gefpy.gene_exp_cy import GEF
from gefpy.bgef_reader_cy import BgefR
# from datetime import datetime
#from gefpy.bgef_writer_cy import generate_bgef
# from gefpy.bgef_writer_cy import SdDataToGef
# from gefpy.bgef_writer_cy import st_sttogef
# import h5py
# import faulthandler
# faulthandler.enable()
# import traceback
# import ctypes
# import cv2
# import matplotlib.pyplot as plt

def main():
    path="/ldfssz1/ST_BI/USER/zhaozijian/geftool/build/FP200000443TL_E2.bgef"
    #generate_bgef(path, "test.gef", 8,  region=[0,26460,0,26456])
    # a = datetime.now()
    # bgef_old = GEF(path, 200)
    # gene_num = bgef_old.get_gene_num()
    # cell_names, cell_index, count = bgef_old.get_exp_data()
    # gene_indexes, gene_names = bgef_old.get_gene_data()
    # b = datetime.now()
    # print(gene_names)
    # print(cell_names)
    # # print(bgef_old.get_cell_num())
    # # print(bgef_old.get_gene_num())

    # c = datetime.now()
    bgef = BgefR(path, 1, 4)
    gene_num1 = bgef.get_gene_num()
    uniq_cells, rows, count = bgef.get_exp_data()
    #cols, uniq_genes = bgef.get_gene_data()
    # print("end")
    # x,y = bgef.get_offset()
    print(uniq_cells)
    # cell_index, gene_index, count = bgef.get_sparse_matrix_indices2()
    # d = datetime.now()

    # print("old:", (b - a).seconds)
    # print("new:", (d - c).seconds)

    # generate gem
    #bgef.to_gem("../test_data/FP200000617TL_B6/FP200000617TL_B6.gefpy.bin1.3.gem")

    # exp = np.empty(bgef.get_expression_num()*3, dtype=np.int32)
    #exp = bgef.get_expression()
    # exp.astype(EXP_DTYPE)
    # print(bgef.get_gene_num())
    # print(bgef.get_expression_num())
    # print(bgef.get_whole_exp_matrix_shape().shape)
    # print(bgef.get_whole_exp_matrix_shape()[0])
    # print(bgef.get_whole_exp_matrix_shape()[1])
    # gene_index, gene_names = bgef.get_sparse_matrix_indices_of_gene()
    # cell_names, cell_index, count = bgef.get_sparse_matrix_indices_of_exp()
    # print(gene_names)
    # print(cell_names)

    #print(bgef.get_gene_names())
    #print(bgef.get_cell_names())
    # indices, indptr, count = bgef.get_sparse_matrix_indices()
    # print(indices)
    # print(indptr)
    # print(count)
    
    # print(cell_index)
    # print(gene_index)

    # bx1=6404-4
    # bx2=18904-4
    # by1=6011-11
    # by2=18211-11
    # generate_bgef("/jdfssz1/ST_BIGDATA/USER/zhaofuxiang/troubleShooting/tissuecut/gefpy_fail2/1.raw.gef", "out.gef", region=[bx1,bx2,by1,by2])

    # with h5py.File(path, mode='r') as h5f:
    #     h5exp = h5f['geneExp']['bin1']['expression']
    #     x_offset = h5exp.attrs['minX']
    #     y_offset = h5exp.attrs['minY']
    #     print(x_offset, y_offset)

    # img = np.zeros((9000,9000,1), dtype=np.uint8)
    # area = np.array([[7839,8784],
    #         [7835,8786],
    #         [7827,8796],
    #         [7825,8801],
    #         [7825,8807],
    #         [7830,8814],
    #         [7837,8817],
    #         [7845,8817],
    #         [7853,8810],
    #         [7855,8806],
    #         [7855,8798],
    #         [7850,8787],
    #         [7844,8784]
    #         ])
    # img = cv2.fillPoly(img, [area], 255)
    
    # plt.imshow(img)
    # plt.show()

def test1():
    path="/ldfssz1/ST_BI/USER/zhaozijian/geftool_zj/build/test"
    with h5py.File(path, mode='r') as h5f:
        ver = h5f.attrs['geftool_ver']
        print(ver)

def test2():
    path="/ldfssz1/ST_BI/USER/zhaozijian/geftool_zj/build/test"
    a = [10,30,50,14,20]
    ak = np.asarray(a,dtype=np.uint)
    SdDataToGef(path, 10, 5, ak)
    b = ["dljl","jdljfljjl","aajfjjj","bbdjkjl","cccjdlj"]
    bk = np.asarray(b, dtype=np.byte)
    st_sttogef(5,bk)

def StereoDataToGef(path, bin):
    h5f = h5py.File(path,"w")
    geneExp = h5f.create_group("geneExp")
    binsz = "bin"+str(bin)
    bing = geneExp.create_group(binsz)

    genetyp = np.dtype({'names':['gene','offset','count'], 'formats':['S32', np.uint32, np.uint32]})
    exptpy = np.dtype({'names':['x','y','count'], 'formats':[np.uint32, np.uint32,np.uint16]})

    bing["expression"] = np.array([(10,20,2), (20,40,3)], dtype=exptpy)
    bing["gene"] = np.array([("gene1",0,21), ("gene2",21,3)], dtype=genetyp)

    bing["expression"].attrs.create("minX", 100)
    bing["expression"].attrs.create("minY", 100)
    bing["expression"].attrs.create("maxX", 1000)
    bing["expression"].attrs.create("maxY", 1000)
    bing["expression"].attrs.create("minExp", 200)
    bing["expression"].attrs.create("resolution", 700)

    h5f.attrs.create("version", 2)
    h5f.close()


from gefpy.bgef_writer_cy import generate_bgef
from gefpy.cgef_writer_cy import generate_cgef
from gefpy.cgef_adjust_cy import CgefAdjust

def test():
    #bgef生成
    generate_bgef("/jdfssz2/ST_BIOINTEL/P20Z10200N0039/01.ShanghaiNeu/shanghai_mousebrain/source/20210428-T318-Z3-L-M024-02_FP200000443TL_E2/FP200000443TL_E2.gem.gz", "outbgef", 10, 1)
#     #cgef生成
#     generate_cgef(outcgef, bgef, mask, [256,256])
#     #获取修正前数据
#     cad = CgefAdjust()
#     cad.get_cell_data(bgef, cgef)
#     #写修正后的数据
#     cad.write_cgef_adjustdata(outpath, cell, dnb)

def test_segmentation_fault():
    # 对于segmentation fault并不能catch到异常，即此处try没效果
    try:
        ctypes.string_at(0)
    except Exception as e:
        print(traceback.format_exc())


if __name__ == "__main__":
    test()