import sys

import numpy as np
from gefpy.cgef_reader_cy import CgefR

from datetime import datetime

def main():
    cgef = CgefR("../test_data/FP200000617TL_B6/FP200000617TL_B6.gefpy.cgef", True)
    # get cells
    cells = cgef.get_cells()

    # get genes
    genes = cgef.get_genes()

    a = datetime.now()
    # 限制区域
    cgef.restrict_region(1000, 15000, 1000, 15000)
    b = datetime.now()
    print("microseconds restrict_region: ", (b - a).microseconds)
    # print(cgef.get_expression_num())
    # print(cgef.get_gene_num())
    # print(cgef.get_cell_num())
    # print(cgef.get_cell_names())
    # print(cgef.get_cell_names().shape)
    # gene_names = cgef.get_gene_names()
    # print(gene_names)
    # print(type(gene_names[0]))

    # Delete genes with zero expression after restrict_region.
    cgef.update_gene_info()

    print("cgef.get_cell_num() 2:", cgef.get_cell_num())
    print("cgef.get_gene_num() 2:", cgef.get_gene_num())
    print("cgef.get_expression_num() 2:", cgef.get_expression_num())
    print(cgef.get_cell_names())
    print(cgef.get_cell_names().shape)
    gene_names = cgef.get_gene_names()
    print(gene_names)
    print(type(gene_names[0]))

    indices, indptr, count = cgef.get_sparse_matrix_indices(order='gene')
    print(indices)
    print(indices.shape)
    print(indptr)
    print("indptr.shape : ", indptr.shape)
    print(count)
    print(count.shape)

    from scipy import sparse
    test_csr = sparse.csr_matrix((count, indices, indptr), shape=(cgef.get_gene_num(), cgef.get_cell_num()))
     #     csr_matrix((data, indices, indptr), [shape=(M, N)])
     # |          is the standard CSR representation where the column indices for
     # |          row i are stored in ``indices[indptr[i]:indptr[i+1]]`` and their
     # |          corresponding values are stored in ``data[indptr[i]:indptr[i+1]]``.
     # |          If the shape parameter is not supplied, the matrix dimensions
     # |          are inferred from the index arrays.
    #
    # print(indices[indptr[1]:indptr[2]])
    # print(test_csr.shape)
    # print(test_csr.getrow(1))
    # (0, 16911)	1
    # 坐标， 值
    # print(test_csr.getrow(1)[0, 16911])
    #
    indices, indptr, count = cgef.get_sparse_matrix_indices(order='cell')
    test_csr = sparse.csr_matrix((count, indices, indptr), shape=(cgef.get_cell_num(), cgef.get_gene_num()))
    print(indices)
    print(max(indices))
    print(indices.shape)
    print(indptr)
    print("indptr.shape : ", indptr.shape)


    # Must call this function before reuse restrict_*.
    cgef.free_restriction()

    return 0


if __name__ == "__main__":
    sys.exit(main())
