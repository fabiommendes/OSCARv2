import csv

from oscar.config import *
from oscar.runtime.oscar import *

break_if_error = False
drive_list = []
param_list = []

# SAVE DRIVERS
drive_list.append(
    [
        data_EFF,
        data_LULCC,
        data_ECH4,
        data_EN2O,
        data_Ehalo,
        data_ENOX,
        data_ECO,
        data_EVOC,
        data_ESO2,
        data_ENH3,
        data_EOC,
        data_EBC,
        data_RFant,
        data_RFnat,
    ]
)

# SAVE PARAMETERS
param_list.append(
    [
        mod_OSNKstruct,
        mod_OSNKchem,
        mod_OSNKtrans,
        mod_LSNKnpp,
        mod_LSNKrho,
        mod_LSNKpreind,
        mod_LSNKtrans,
        mod_LSNKcover,
        mod_EFIREpreind,
        mod_EFIREtrans,
        mod_ELUCagb,
        mod_EHWPbb,
        mod_EHWPtau,
        mod_EHWPfct,
        mod_OHSNKtau,
        mod_OHSNKfct,
        mod_OHSNKtrans,
        mod_EWETpreind,
        mod_AWETtrans,
        mod_HVSNKtau,
        mod_HVSNKtrans,
        mod_HVSNKcirc,
        mod_O3Tregsat,
        mod_O3Temis,
        mod_O3Tclim,
        mod_O3Tradeff,
        mod_O3Sfracrel,
        mod_O3Strans,
        mod_O3Snitrous,
        mod_O3Sradeff,
        mod_SO4regsat,
        mod_SO4load,
        mod_SO4radeff,
        mod_POAconv,
        mod_POAregsat,
        mod_POAload,
        mod_POAradeff,
        mod_BCregsat,
        mod_BCload,
        mod_BCradeff,
        mod_BCadjust,
        mod_NO3load,
        mod_NO3radeff,
        mod_SOAload,
        mod_SOAradeff,
        mod_DUSTload,
        mod_DUSTradeff,
        mod_SALTload,
        mod_SALTradeff,
        mod_CLOUDsolub,
        mod_CLOUDerf,
        mod_CLOUDpreind,
        mod_ALBBCreg,
        mod_ALBBCrf,
        mod_ALBBCwarm,
        mod_ALBLCflux,
        mod_ALBLCalb,
        mod_ALBLCcover,
        mod_ALBLCwarm,
        mod_TEMPresp,
        mod_TEMPpattern,
        mod_PRECresp,
        mod_PRECradfact,
        mod_PRECpattern,
        mod_ACIDsurf,
        mod_SLR,
    ]
)

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
# plt.show()
writer = csv.writer(open("results/single-drive_test.csv", "w"))
writer.writerows(drive_list)
writer = csv.writer(open("results/single-param_test.csv", "w"))
writer.writerows(param_list)
writer = csv.writer(open("results/#OUT.USELESS.empty", "w"))
