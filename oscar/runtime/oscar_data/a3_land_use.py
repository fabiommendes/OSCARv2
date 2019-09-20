import os

import numpy as np

from .a1_regions import nb_regionJ, nb_kind, nb_regionI, nb_sector, regionJ_index, regionI_index, ind_final, kLUC, biome_index, nb_biome
from .a2_greenhouse import ind_cdiac
from ...config import dty, data_LULCC, mod_LSNKcover, scen_LULCC
from ...data import load_data

##################################################
#   2. LAND-USE CHANGE
##################################################
LUC = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_biome, nb_biome], dtype=dty)  # {Mha/yr}
HARV = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_biome], dtype=dty)  # {GtC/yr}
SHIFT = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_biome, nb_biome], dtype=dty)  # {Mha/yr}

# =========
# 2.1. LUH1
# =========

# load land-use data from LUH1
# from [Hurtt et al., 2011] and updated
bio = ["des", "for", "shr", "gra", "cro", "pas", "urb"]

# LUC
if data_LULCC[:3] == "LUH":
    for b1 in range(len(bio)):
        for b2 in range(len(bio)):
            path = f"data/LandUse_{data_LULCC}/#DATA.LandUse_{data_LULCC}_{mod_LSNKcover}.1501-2015_114reg1.LUC_{bio[b1]}2{bio[b2]}.csv"
            if os.path.isfile(path):
                TMP = load_data(path)
                for i in range(1, 114 + 1):
                    LUC[1: min(ind_final, ind_cdiac) + 1, regionJ_index[i], 0, kLUC, regionI_index[i], biome_index[bio[b1]], biome_index[bio[b2]], ] += TMP[200: 200 + min(ind_final, ind_cdiac) - 1 + 1, i - 1]

# HARV
if data_LULCC[:3] == "LUH":
    for b in range(len(bio)):
        path = f"data/LandUse_{data_LULCC}/#DATA.LandUse_{data_LULCC}_{mod_LSNKcover}.1501-2015_114reg1.HARV_{bio[b]}.csv"
        if os.path.isfile(path):
            TMP = load_data(path)
            for i in range(1, 114 + 1):
                HARV[1: min(ind_final, ind_cdiac) + 1, regionJ_index[i], 0, kLUC, regionI_index[i], biome_index[bio[b]]] += TMP[200: 200 + min(ind_final, ind_cdiac) - 1 + 1, i - 1]

# SHIFT
if data_LULCC[:3] == "LUH":
    for b1 in range(len(bio)):
        for b2 in range(b1, len(bio)):
            path = f"data/LandUse_{data_LULCC}/#DATA.LandUse_{data_LULCC}_{mod_LSNKcover}.1501-2015_114reg1.SHIFT_{bio[b1]}2{bio[b2]}.csv"
            if os.path.isfile(path):
                TMP = load_data(path)
                for i in range(1, 114 + 1):
                    SHIFT[ 1: min(ind_final, ind_cdiac) + 1, regionJ_index[i], 0, kLUC, regionI_index[i], biome_index[bio[b1]], biome_index[bio[b2]], ] += TMP[200: 200 + min(ind_final, ind_cdiac) - 1 + 1, i - 1]

# ========
# 2.2. RCP
# ========

# initialization of projected drivers
LUCproj = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_biome, nb_biome], dtype=dty)
HARVproj = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_biome], dtype=dty)
SHIFTproj = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_biome, nb_biome], dtype=dty)

# projection of land-use data under RCP scenarios
# from [Hurtt et al., 2011] and [Meinshausen et al., 2011]

# LUC
if (scen_LULCC[:3] == "RCP") & (ind_final > ind_cdiac):
    for b1 in range(len(bio)):
        for b2 in range(len(bio)):
            path = f"data/LandUse_RCP/#DATA.LandUse_RCP_{mod_LSNKcover}.2006-2100_114reg1.rcp{scen_LULCC[3]}{scen_LULCC[5]}_LUC_{bio[b1]}2{bio[b2]}.csv"
            if os.path.isfile(path):
                TMP = load_data(path)
                for i in range(1, 114 + 1):
                    LUCproj[306: min(ind_final, 400) + 1, regionJ_index[i], 0,kLUC, regionI_index[i], biome_index[bio[b1]], biome_index[bio[b2]], ] += TMP[: min(ind_final, 400) - 306 + 1, i - 1]

# HARV
if (scen_LULCC[:3] == "RCP") & (ind_final > ind_cdiac):
    for b in range(len(bio)):
        path = f"data/LandUse_RCP/#DATA.LandUse_RCP_{mod_LSNKcover}.2006-2100_114reg1.rcp{scen_LULCC[3]}" f"{scen_LULCC[5]}_HARV_{bio[b]}.csv"

        if os.path.isfile(path):
            TMP = load_data(path)
            for i in range(1, 114 + 1):
                HARVproj[
                306: min(ind_final, 400) + 1, regionJ_index[i], 0, kLUC, regionI_index[i], biome_index[bio[b]]] += TMP[: min(ind_final, 400) - 306 + 1, i - 1]

# SHIFT
if (scen_LULCC[:3] == "RCP") & (ind_final > ind_cdiac):
    for b1 in range(len(bio)):
        for b2 in range(b1, len(bio)):
            path = f"data/LandUse_RCP/#DATA.LandUse_RCP_{mod_LSNKcover}.2006-2100_114reg1.rcp{scen_LULCC[3]}{scen_LULCC[5]}_SHIFT_{bio[b1]}2{bio[b2]}.csv"
            if os.path.isfile(path):
                TMP = load_data(path)
                for i in range(1, 114 + 1):
                    SHIFTproj[306: min(ind_final, 400) + 1, regionJ_index[i], 0, kLUC, regionI_index[i], biome_index[bio[b1]], biome_index[bio[b2]]] += TMP[: min(ind_final, 400) - 306 + 1, i - 1]

# ==================
# 2.A. FINAL DATASET
# ==================

# datasets mixed following various criteria
for arr, arrproj in [
    (LUC, LUCproj), (HARV, HARVproj), (SHIFT, SHIFTproj),
]:

    # stop emissions
    if (scen_LULCC == "stop") & (ind_final > ind_cdiac):
        arr[ind_cdiac + 1:, ...] = 0

    # constant emissions
    elif (scen_LULCC == "cst") & (ind_final > ind_cdiac):
        arr[ind_cdiac + 1:, ...] = arr[ind_cdiac, ...][np.newaxis, ...]

        # RCP scenarios
    # always raw discontinuity
    elif (scen_LULCC[:3] == "RCP") & (ind_final > ind_cdiac):
        arr[ind_cdiac + 1:, ...] = arrproj[ind_cdiac + 1:, ...]

# Delete individual datasets
del LUCproj
del HARVproj
del SHIFTproj
