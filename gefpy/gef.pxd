# -*- coding: utf-8 -*-
# distutils: language=c++
# cython: language_level=3, boundscheck=False
# cython: c_string_type=unicode, c_string_encoding=utf8
# Created by huangzhibo on 2022/01/01

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool


cdef extern from "gef.h":
    ctypedef union Coordinate:
        unsigned int pos[2] # dnb coordinates x, y
        unsigned long long int pos_id

    # ctypedef struct ReadSnpInfo:
