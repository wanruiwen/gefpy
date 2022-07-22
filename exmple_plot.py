import sys

from gefpy.bgef_reader_cy import BgefR
from gefpy.cgef_adjust_cy import CgefAdjust
import numpy as np
from scipy.sparse import csr_matrix

def test1():
    gef = BgefR("/ldfssz1/ST_BI/USER/stereopy/test/xujunhao/data/debug/zhangxiaoyu/TestData.gef", 100, 4)
    gene_num = gef.get_gene_num()
    offset_x, offset_y = gef.get_offset()
    gef_attr = gef.get_exp_attr()

    uniq_cells, rows, count = gef.get_exp_data()
    cell_num = len(uniq_cells)

    cols, uniq_genes = gef.get_gene_data()
    position = np.array(list((zip(np.right_shift(uniq_cells, 32), np.bitwise_and(uniq_cells, 0xffff))))).astype('uint32')

    exp_matrix = csr_matrix((count, (rows, cols)), shape=(cell_num, gene_num), dtype=np.uint32)
    print("eok")

def test2():
    arry1 = [480,2680,520,2560,880,1680,920,1600,1200,1680,1320,1720,1800,1960,2120,3800,2120,3880,2000,3880,1240,3840,1000,3800,960,3760,600,3040,480,2760,480,2720]
    arry2 = [2200,2680,2680,1560,2720,1480,3400,2280,3480,3280,3480,4040,3160,3720]
    arry3 = [3960,2840,4280,1240,4360,1240,5200,1360,5360,1400,5400,1680,5520,2640,5520,2960,5480,3080,5320,3440,5120,3880,5080,3960,5040,4000,4840,4040,4280,4120,4120,4120,4080,3840,4000,3200]
    dic = [arry1,arry2,arry3]
    cg = CgefAdjust()
    cg.create_Region_Bgef("/ldfssz1/ST_BI/USER/stereopy/test/panjie/stereopy/DP8400013846TR_F5.gef", "jfkldjl", dic)

if __name__ == "__main__":
    sys.exit(test2())