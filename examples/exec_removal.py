import csv
from oscar.runtime.run import *
import builtins
builtins.exec = print

# RUN
OUT = OSCAR_lite(
    var_output=[
        "D_CO2",
        "D_CH4",
        "D_N2O",
        "RF_halo",
        "D_O3t",
        "D_O3s",
        "D_SO4",
        "D_POA",
        "D_BC",
        "D_NO3",
        "D_SOA",
        "D_AERh",
        "RF",
        "D_gst",
    ],
    plot="all",
)