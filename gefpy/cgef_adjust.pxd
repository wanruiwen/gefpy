# distutils: language=c++
# cython: language_level=3
import numpy as np
cimport numpy as np

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

from .gef cimport *



cdef extern from "cellAdjust.h" nogil:

    ctypedef struct cellgem_label:
        unsigned int geneid
        int x
        int y
        int midcnt
        unsigned int cellid

    ctypedef struct sapBgefData:
        unsigned short genecnt;
        unsigned int midcnt;
        int x;
        int y;

    cdef cppclass cellAdjust:
        cellAdjust()
        void readBgef(const string &strinput)
        void readCgef(const string &strinput)
        unsigned int getCellLabelgem(vector[string] &gene_list, vector[cellgem_label] &vecCellgem)
        void writeCellAdjust(const string &path, Cell *cellptr, unsigned int cellcnt, DnbExpression *dnbptr, unsigned int dnbcnt)
        void createRegionGef(const string &strout)
        void getRegionGenedata(vector[vector[int]] &vec)
        void readRawCgef(const string &strcgef)
        void getRegionCelldata(vector[vector[int]] &m_vecpos)
        void writeToCgef(const string &outpath)

        void getSapRegion(const string &strinput, int bin, int thcnt, vector[vector[int]] &vecpos, vector[sapBgefData] &vecdata)
