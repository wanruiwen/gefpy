import sys

from gefpy.bgef_reader_cy import BgefR
from gefpy.cgef_reader_cy import CgefR
from gefpy.cgef_adjust_cy import CgefAdjust
import numpy as np
from scipy.sparse import csr_matrix

def test1():
    gef = BgefR("/ldfssz1/ST_BI/USER/stereopy/test/xujunhao/data/debug/zhangxiaoyu/TestData.gef", 1, 4)
    region = [0,2000,0,2000]
    gene=["Pharaoh_ant_13172","Pharaoh_ant_15152","ND-24","Pharaoh_ant_15157"]

    uniq_cell, gene_names, count, cell_ind, gene_ind = gef.get_filtered_data(region,gene)

    gene_num = gene_names.size
    cell_num = uniq_cell.size

    exp_matrix = csr_matrix((count, (cell_ind, gene_ind)), shape=(cell_num, gene_num), dtype=np.uint32)

    position = np.array(list((zip(np.right_shift(uniq_cell, 32), np.bitwise_and(uniq_cell, 0xffff))))).astype('uint32')

def test4():
    gef = CgefR("/ldfssz1/ST_BI/USER/stereopy/test/xujunhao/data/debug/zhangxiaoyu/TestData.gef")
    region = [0,2000,0,2000]
    #gene=["Pharaoh_ant_13172","Pharaoh_ant_15152","ND-24","Pharaoh_ant_15157"]
    gene=[]
    uniq_cell, gene_names, count, cell_ind, gene_ind = gef.get_filtered_data(region,gene)

    gene_num = gene_names.size
    cell_num = uniq_cell.size
    
    exp_matrix = csr_matrix((count, (cell_ind, gene_ind)), shape=(cell_num, gene_num), dtype=np.uint32)

    position = np.array(list((zip(np.right_shift(uniq_cell, 32), np.bitwise_and(uniq_cell, 0xffff))))).astype('uint32')



def test2():
    arry1 = [1160, 2600, 1640, 1400, 1680, 1320, 1760, 1240, 1960, 1080, 2280, 1080, 2920, 1160, 3120, 1200, 3520, 1440, 3640, 1520, 3720, 1600, 3880, 1880, 3920, 1960, 3920, 3600, 3880, 3640, 3320, 4160, 3000, 4240, 2800, 4280, 2600, 4280, 1760, 3920, 1600, 3840, 1560, 3800]
    dic = [arry1]
    cg = CgefAdjust()
    cg.create_Region_Bgef("/ldfssz1/ST_BI/USER/stereopy/test/panjie/stereopy/DP8400013846TR_F5.gef", "kktest.bgef", dic)

def test3():
    from gefpy.bgef_writer_cy import generate_bgef
    gem2 = '/ldfssz1/ST_BI/USER/stereopy/test/xujunhao/data/FP200000463BL_C4/FP200000463BL_C4.tissue.gem.gz'
    gef2_file = 'test2.gef'
    generate_bgef(gem2, gef2_file)
    gem = '/ldfssz1/ST_BI/USER/stereopy/test/xujunhao/data/mouse_embroy_e12_5/E12.5_E1S3.gem.gz'
    gef_file = 'test1.gef'
    generate_bgef(gem, gef_file)

if __name__ == "__main__":
    sys.exit(test2())