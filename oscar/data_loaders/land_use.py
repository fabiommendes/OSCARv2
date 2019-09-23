import os

import numpy as np

from .regions import nb_regionJ, nb_kind, nb_regionI, nb_sector, regionJ_index, regionI_index, kLUC, biome_index, nb_biome
from .greenhouse import ind_cdiac
from ..data import load_data
from .. import conf

##################################################
#   2. LAND-USE CHANGE
##################################################
LUC = np.zeros([conf.ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_biome, nb_biome], dtype=conf.dty)  # {Mha/yr}
HARV = np.zeros([conf.ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_biome], dtype=conf.dty)  # {GtC/yr}
SHIFT = np.zeros([conf.ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_biome, nb_biome], dtype=conf.dty)  # {Mha/yr}

# =========
# 2.1. LUH1
# =========

# load land-use config.data_loaders from LUH1
# from [Hurtt et al., 2011] and updated
bio = ["des", "for", "shr", "gra", "cro", "pas", "urb"]

# LUC
if conf.data_LULCC[:3] == "LUH":
    for b1 in range(len(bio)):
        for b2 in range(len(bio)):
            path = f"data/LandUse_{conf.data_LULCC}/#DATA.LandUse_{conf.data_LULCC}_{conf.mod_LSNKcover}.1501-2015_114reg1.LUC_{bio[b1]}2{bio[b2]}.csv"
            if os.path.isfile(path):
                TMP = load_data(path)
                for i in range(1, 114 + 1):
                    LUC[1: min(conf.ind_final, ind_cdiac) + 1, regionJ_index[i], 0, kLUC, regionI_index[i], biome_index[bio[b1]], biome_index[bio[b2]], ] += TMP[200: 200 + min(conf.ind_final, ind_cdiac) - 1 + 1, i - 1]

# HARV
if conf.data_LULCC[:3] == "LUH":
    for b in range(len(bio)):
        path = f"data/LandUse_{conf.data_LULCC}/#DATA.LandUse_{conf.data_LULCC}_{conf.mod_LSNKcover}.1501-2015_114reg1.HARV_{bio[b]}.csv"
        if os.path.isfile(path):
            TMP = load_data(path)
            for i in range(1, 114 + 1):
                HARV[1: min(conf.ind_final, ind_cdiac) + 1, regionJ_index[i], 0, kLUC, regionI_index[i], biome_index[bio[b]]] += TMP[200: 200 + min(conf.ind_final, ind_cdiac) - 1 + 1, i - 1]

# SHIFT
if conf.data_LULCC[:3] == "LUH":
    for b1 in range(len(bio)):
        for b2 in range(b1, len(bio)):
            path = f"data/LandUse_{conf.data_LULCC}/#DATA.LandUse_{conf.data_LULCC}_{conf.mod_LSNKcover}.1501-2015_114reg1.SHIFT_{bio[b1]}2{bio[b2]}.csv"
            if os.path.isfile(path):
                TMP = load_data(path)
                for i in range(1, 114 + 1):
                    SHIFT[ 1: min(conf.ind_final, ind_cdiac) + 1, regionJ_index[i], 0, kLUC, regionI_index[i], biome_index[bio[b1]], biome_index[bio[b2]], ] += TMP[200: 200 + min(conf.ind_final, ind_cdiac) - 1 + 1, i - 1]

# ========
# 2.2. RCP
# ========

# initialization of projected drivers
LUCproj = np.zeros([conf.ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_biome, nb_biome], dtype=conf.dty)
HARVproj = np.zeros([conf.ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_biome], dtype=conf.dty)
SHIFTproj = np.zeros([conf.ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_biome, nb_biome], dtype=conf.dty)

# projection of land-use config.data_loaders under RCP scenarios
# from [Hurtt et al., 2011] and [Meinshausen et al., 2011]

# LUC
if (conf.scen_LULCC[:3] == "RCP") & (conf.ind_final > ind_cdiac):
    for b1 in range(len(bio)):
        for b2 in range(len(bio)):
            path = f"data/LandUse_RCP/#DATA.LandUse_RCP_{conf.mod_LSNKcover}.2006-2100_114reg1.rcp{conf.scen_LULCC[3]}{conf.scen_LULCC[5]}_LUC_{bio[b1]}2{bio[b2]}.csv"
            if os.path.isfile(path):
                TMP = load_data(path)
                for i in range(1, 114 + 1):
                    LUCproj[306: min(conf.ind_final, 400) + 1, regionJ_index[i], 0, kLUC, regionI_index[i], biome_index[bio[b1]], biome_index[bio[b2]], ] += TMP[: min(conf.ind_final, 400) - 306 + 1, i - 1]

# HARV
if (conf.scen_LULCC[:3] == "RCP") & (conf.ind_final > ind_cdiac):
    for b in range(len(bio)):
        path = f"data/LandUse_RCP/#DATA.LandUse_RCP_{conf.mod_LSNKcover}.2006-2100_114reg1.rcp{conf.scen_LULCC[3]}" f"{conf.scen_LULCC[5]}_HARV_{bio[b]}.csv"

        if os.path.isfile(path):
            TMP = load_data(path)
            for i in range(1, 114 + 1):
                HARVproj[
                306: min(conf.ind_final, 400) + 1, regionJ_index[i], 0, kLUC, regionI_index[i], biome_index[bio[b]]] += TMP[: min(conf.ind_final, 400) - 306 + 1, i - 1]

# SHIFT
if (conf.scen_LULCC[:3] == "RCP") & (conf.ind_final > ind_cdiac):
    for b1 in range(len(bio)):
        for b2 in range(b1, len(bio)):
            path = f"data/LandUse_RCP/#DATA.LandUse_RCP_{conf.mod_LSNKcover}.2006-2100_114reg1.rcp{conf.scen_LULCC[3]}{conf.scen_LULCC[5]}_SHIFT_{bio[b1]}2{bio[b2]}.csv"
            if os.path.isfile(path):
                TMP = load_data(path)
                for i in range(1, 114 + 1):
                    SHIFTproj[306: min(conf.ind_final, 400) + 1, regionJ_index[i], 0, kLUC, regionI_index[i], biome_index[bio[b1]], biome_index[bio[b2]]] += TMP[: min(conf.ind_final, 400) - 306 + 1, i - 1]

# ==================
# 2.A. FINAL DATASET
# ==================

# datasets mixed following various criteria
for arr, arrproj in [
    (LUC, LUCproj), (HARV, HARVproj), (SHIFT, SHIFTproj),
]:

    # stop emissions
    if (conf.scen_LULCC == "stop") & (conf.ind_final > ind_cdiac):
        arr[ind_cdiac + 1:, ...] = 0

    # constant emissions
    elif (conf.scen_LULCC == "cst") & (conf.ind_final > ind_cdiac):
        arr[ind_cdiac + 1:, ...] = arr[ind_cdiac, ...][np.newaxis, ...]

        # RCP scenarios
    # always raw discontinuity
    elif (conf.scen_LULCC[:3] == "RCP") & (conf.ind_final > ind_cdiac):
        arr[ind_cdiac + 1:, ...] = arrproj[ind_cdiac + 1:, ...]

# Delete individual datasets
del LUCproj
del HARVproj
del SHIFTproj
