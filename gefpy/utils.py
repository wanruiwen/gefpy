#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Provides convenient invocation methods for other modules.
"""

import h5py
import numpy as np

def gef_is_cell_bin(gef_file):
    """
    Determine if the file is a cgef file

    :param gef_file: input the gef path

    :return: the return True for cgef, False otherwise.
    :rtype: bool
    """

    h5f = h5py.File(gef_file)
    if 'cellBin' in h5f:
        h5f.close()
        return True
    else:
        h5f.close()
        return False

def StereoDataToGef(path, bin, expression, gene, attr):
    """
    Write stereodata to cgef file

    :param path: set the output path
    :param bin: set the bin size
    :param expression: input the expression data
    :param gene: input the gene data
    :param attr: input the attr data

    """
    h5f = h5py.File(path,"w")
    geneExp = h5f.create_group("geneExp")
    binsz = "bin"+str(bin)
    bing = geneExp.create_group(binsz)

    #genetype = np.dtype({'names':['gene','offset','count'], 'formats':['S32', 'i','i']})
    #exptype = np.dtype({'names':['x','y','count'], 'formats':['i', 'i','i']})

    geneExp[binsz]["expression"] = expression #np.arry([(10,20,2), (20,40,3)], dtype=exptype)
    geneExp[binsz]["gene"] = gene #np.arry([("gene1",0,21), ("gene2",21,3)], dtype=genetype)

    bing["expression"].attrs.create("minX", attr[0])
    bing["expression"].attrs.create("minY", attr[1])
    bing["expression"].attrs.create("maxX", attr[2])
    bing["expression"].attrs.create("maxY", attr[3])
    bing["expression"].attrs.create("minExp", attr[4])
    bing["expression"].attrs.create("resolution", attr[5])

    h5f.attrs.create("version", 2)
    
    h5f.close()