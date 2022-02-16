import gc

import h5py

from gefpy.bgef_writer_cy import generate_bgef
from gefpy import plot
import os
# os.environ['HDF5_USE_FILE_LOCKING']='FALSE'

# 提取sub bGEF
bgef_file = "../test_data/FP200000617TL_B6/FP200000617TL_B6.gefpy.bins.gef"
# region = [1000, 2000, 1000, 2000]
bin1_bgef_file2 = "../test_data/FP200000617TL_B6/FP200000617TL_B6.gefpy.bins.plot.sub.gef"
generate_bgef(bgef_file, bin1_bgef_file2)
# generate_bgef(bgef_file, bin1_bgef_file2, region=region)

plot.save_exp_heat_map(
        # "../test_data/FP200000617TL_B6/FP200000617TL_B6.gefpy.bin1.3.gef",
        # "../../test_data/FP200000617TL_B6/FP200000617TL_B6.bgef.h5.gef",
        bin1_bgef_file2,
        "../test_data/FP200000617TL_B6/cellbin.png")
plot.cgef_stat("../test_data/FP200000617TL_B6/FP200000617TL_B6.gefpy.cgef", "../test_data/FP200000617TL_B6/")