import sys

from gefpy.bgef_reader_cy import BgefR
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


if __name__ == "__main__":
    sys.exit(test1())