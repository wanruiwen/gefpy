from gefpy import plot

plot.save_exp_heat_map(
        # "../test_data/FP200000617TL_B6/FP200000617TL_B6.gefpy.bin1.3.gef",
        # "../../test_data/FP200000617TL_B6/FP200000617TL_B6.bgef.h5.gef",
        "../test_data/FP200000617TL_B6/FP200000617TL_B6.gefpy.cgef",
        "../test_data/FP200000617TL_B6/cellbin.png")
plot.cgef_stat("../test_data/FP200000617TL_B6/FP200000617TL_B6.gefpy.cgef", "../test_data/FP200000617TL_B6/")