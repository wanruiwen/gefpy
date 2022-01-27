import sys
from gefpy.gene_exp_cy import GEF
from gefpy.bgef_reader_cy import BgefR

def main():
    # bgef_old = GEF("/Users/huangzhibo/workitems/01.github/geftools/test_data/FP200000617TL_B6/stereomics.h5", 1)
    # gene_indexes, gene_names = bgef_old.get_gene_data()
    # cell_names, cell_index, count = bgef_old.get_exp_data()
    # print(gene_names)
    # print(cell_names)

    bgef = BgefR("/Users/huangzhibo/workitems/01.github/geftools/test_data/FP200000617TL_B6/stereomics.h5", 1)
    # print(bgef.get_gene_num())
    # print(bgef.get_expression_num())
    # print(bgef.get_whole_exp_matrix_shape().shape)
    # print(bgef.get_whole_exp_matrix_shape()[0])
    # print(bgef.get_whole_exp_matrix_shape()[1])
    # gene_index, gene_names = bgef.get_sparse_matrix_indices_of_gene()
    # cell_names, cell_index, count = bgef.get_sparse_matrix_indices_of_exp()
    # print(gene_names)
    # print(cell_names)

    print(bgef.get_gene_names())
    print(bgef.get_cell_names())
    cell_index, gene_index, count = bgef.get_sparse_matrix_indices2()
    return 0


if __name__ == "__main__":
    sys.exit(main())