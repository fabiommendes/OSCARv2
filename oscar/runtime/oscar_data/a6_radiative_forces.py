from oscar.data import load_data_and_header
from .a1_regions import ind_final, nb_regionJ, nb_kind, nb_sector, kindRF_index
from .a2_greenhouse import ind_cdiac
from ...config import dty, data_RFant, scen_RFant, data_RFnat, scen_RFnat

print("LOADING: DRIVERS")

##################################################
#   A. VECTORS
##################################################

import numpy as np

##################################################
#   5. RADIATIVE FORCING
##################################################

RFcon = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind], dtype=dty)  # {W/m2}
RFvolc = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind], dtype=dty)  # {W/m2}
RFsolar = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind], dtype=dty)  # {W/m2}

# ==================
# 5.1. Anthropogenic
# ==================

# load RF estimates from IPCC AR5
# from [IPCC, 2013] (annexe 2)
if data_RFant == "IPCC-AR5":
    path = "data/Historic_IPCC-AR5/#DATA.Historic_IPCC-AR5.1750-2011_(11for).RF.csv"
    TMP, lgd = load_data_and_header(path)

    # past estimates
    for x in range(len(lgd)):
        if lgd[x] == "Contrails":
            RFcon[51: min(ind_final, ind_cdiac) + 1, 0, 0, kindRF_index["RFcon"]] += \
                TMP[1: 1 + min(ind_final, ind_cdiac) - 51 + 1, x]

    # arbitrary extension
    if (scen_RFant == "cst") & (ind_final > ind_cdiac):
        for x in range(len(lgd)):
            if lgd[x] == "Contrails":
                RFcon[ind_cdiac + 1:, 0, 0, kindRF_index["RFcon"]] = TMP[-1, x]

# ============
# 5.2. Natural
# ============

# load RF estimates from IPCC AR5
# from [IPCC, 2013] (annexe 2)
path = "data/Historic_IPCC-AR5/#DATA.Historic_IPCC-AR5.1750-2011_(11for).RF.csv"
if data_RFnat == "IPCC-AR5":
    TMP, lgd = load_data_and_header(path)

    # past estimates
    for x in range(len(lgd)):
        if lgd[x] == "Volcano":
            RFvolc[51: min(ind_final, ind_cdiac) + 1, 0, 0, kindRF_index["RFvol"]] += \
                TMP[1: 1 + min(ind_final, ind_cdiac) - 51 + 1, x] - np.mean(TMP[1:, x])
        elif lgd[x] == "Solar":
            RFsolar[51: min(ind_final, ind_cdiac) + 1, 0, 0,
            kindRF_index["RFsol"]] += TMP[1: 1 + min(ind_final, ind_cdiac) - 51 + 1, x]

    # arbitrary extension
    if (scen_RFnat == "cst") & (ind_final > ind_cdiac):
        for x in range(len(lgd)):
            if lgd[x] == "Volcano":
                RFvolc[ind_cdiac + 1:, 0, 0, kindRF_index["RFvol"]] = 0.0
            elif lgd[x] == "Solar":
                RFsolar[ind_cdiac + 1:, 0, 0, kindRF_index["RFsol"]] = np.mean(
                    TMP[-11:, x])
