# -*- coding: utf-8 -*-
# distutils: language=c++
# cython: language_level=3, boundscheck=False
# cython: c_string_type=unicode, c_string_encoding=utf8
# Created by huangzhibo on 2022/01/01
"""
    Provides access to the cgef_adjust interface.
"""

from .cgef_adjust cimport *
import numpy as np
cimport numpy as np
from cython cimport view

cdef class CgefAdjust:
    cdef cellAdjust *c_instance 

    def __cinit__(self):
        self.c_instance = new cellAdjust()

    def __init__(self):
        """
        A class for reading cell bin GEF.

        :param filepath: Input cell bin GEF filepath.
        """
        pass

    def __dealloc__(self):
        del self.c_instance

    def get_cell_data(self, bgef, cgef):
        self.c_instance.readBgef(bgef)
        self.c_instance.readCgef(cgef)
        cdef vector[cellgem_label] vec_cell
        cdef vector[string] genelist
        self.c_instance.getCellLabelgem(genelist, vec_cell)
        return np.asarray(genelist), np.asarray(vec_cell)

    def write_cgef_adjustdata(self):
        cdef vector[Cell] veccell
        cdef Cell c1
        c1.cellid = 0
        c1.offset = 0
        c1.count = 2
        cdef Cell c2
        c2.cellid = 0
        c2.offset = 2
        c2.count = 3
        veccell.push_back(c1)
        veccell.push_back(c2)

        cdef vector[DnbExpression] vecdnb
        cdef DnbExpression dnb1
        dnb1.x = 112
        dnb1.y = 122
        dnb1.count = 1
        dnb1.gene_id = 0
        cdef DnbExpression dnb2
        dnb2.x = 149
        dnb2.y = 150
        dnb2.count = 2
        dnb2.gene_id = 0
        vecdnb.push_back(dnb1)
        vecdnb.push_back(dnb2)
        self.c_instance.writeCellAdjust(veccell, vecdnb)