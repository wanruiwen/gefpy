gefpy.cgef_reader_cy
===========================

Provides access to the cgef_reader interface.

.. py:module:: cgef_reader_cy

.. py:class:: CgefR(filepath)

    .. py:method:: get_expression_num(self)

        Get the number of expression.

    .. py:method:: get_cell_num(self)

        Get the number of cell.

    .. py:method:: get_gene_num(self)

        Get the number of gene.

    .. py:method:: get_gene_names(self)

        Get a list of gene names. The type of gene name is 32 chars.

    .. py:method:: get_cell_names(self)

        Get an array of cell ids. Each cell id is (cell.x <<32 | cell.y)

    .. py:method:: get_cells(self)

        Get cells, each cell include (id, x, y, offset, geneCount, expCount, dnbCount, area, cellTypeID, clusterID)

    .. py:method:: get_genes(self)

        Get genes, each gene include(geneName, offset, cellCount, expCount, maxMIDcount)

    .. py:method:: get_cellid_and_count(self)

        Get the count of each cell in each gene.

        :return:  (cell_id, count)

    .. py:method:: get_geneid_and_count(self)

        Get the count of each gene in each cell.
        
        :return:  (gene_id, count)

    .. py:method:: get_cellborders(self)

        Gets cell borders.
        
        :return: borders_list

