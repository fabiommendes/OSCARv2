import os

import numpy as np

from .regions import nb_regionJ, nb_kind, nb_regionI, nb_sector, regionJ_index, regionI_index, ind_final, kindGHG_index
from .greenhouse import ind_cdiac, ind_edgar
from ..data import load_data_and_header, load_data
from ..constants import HFC, PFC, ODS
from ..config import dty, scen_Ehalo, mod_regionI, data_Ehalo, mod_DATAscen, mod_regionJ

nb_HFC = len(HFC)
nb_PFC = len(PFC)
nb_ODS = len(ODS)

##################################################
#   3. HALOGENATED COMPOUNDS
##################################################

EHFC = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_HFC], dtype=dty)  # {kt/yr}
EPFC = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_PFC], dtype=dty)  # {kt/yr}
EODS = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_ODS], dtype=dty)  # {kt/yr}

# ==========
# 3.1. EDGAR
# ==========

# load emissions from EDGAR v4.2
# see [JRC, 2011]

# HFCs
EHFCedgar = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_HFC], dtype=dty)
for VAR in HFC:
    path = f"data/EHaloComp_EDGAR/#DATA.EHaloComp_EDGAR.1970-{1700 + ind_edgar}_114reg0.E{VAR}.csv"
    TMP = load_data(path)
    for i in range(114 + 1):
        EHFCedgar[270: ind_edgar + 1, regionJ_index[i], 0, kindGHG_index["HFC"], regionI_index[i], HFC.index(VAR)] += TMP[: ind_edgar - 270 + 1, i]

# PFCs
EPFCedgar = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_PFC], dtype=dty)
for VAR in PFC:
    path = f"data/EHaloComp_EDGAR/#DATA.EHaloComp_EDGAR.1970-{1700 + ind_edgar}_114reg0.E{VAR}.csv"
    TMP = load_data(path)
    for i in range(114 + 1):
        EPFCedgar[270: ind_edgar + 1, regionJ_index[i], 0, kindGHG_index["PFC"], regionI_index[i], PFC.index(VAR)] += TMP[: ind_edgar - 270 + 1, i]

# =============
# 3.2. EDGAR-FT
# =============

# load emissions from EDGAR-FT v4.2-FT2010
# see [JRC, 2013]

# HFCs
EHFCeft = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_HFC], dtype=dty)
for VAR in HFC:
    path = f"data/EHaloComp_EDGAR-FT/#DATA.EHaloComp_EDGAR-FT.2008-2010_114reg0.E{VAR}.csv"
    TMP = load_data(path)
    for i in range(114 + 1):
        EHFCeft[308: 310 + 1, regionJ_index[i], 0, kindGHG_index["HFC"], regionI_index[i], HFC.index(VAR)] += TMP[: 310 - 308 + 1, i]

# PFCs
EPFCeft = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_PFC], dtype=dty)
for VAR in PFC:
    path = f"data/EHaloComp_EDGAR-FT/#DATA.EHaloComp_EDGAR-FT.2008-2010_114reg0.E{VAR}.csv"
    TMP = load_data(path)
    for i in range(114 + 1):
        EPFCeft[308: 310 + 1, regionJ_index[i], 0, kindGHG_index["PFC"], regionI_index[i], PFC.index(VAR)] += TMP[: 310 - 308 + 1, i]

# ==========
# 3.3. CMIP5
# ==========

# ODSs
EODScmip5 = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_ODS], dtype=dty)
path = "data/EHaloComp_CMIP5/#DATA.EHaloComp_CMIP5.1765-2005_(16ghg).EODS.csv"
TMP, lgd = load_data_and_header(path)
for x in range(len(lgd)):
    EODScmip5[65: 305 + 1, 0, 0, kindGHG_index["ODS"], 0, ODS.index(lgd[x])] = TMP[: min(ind_final, 305) - 65 + 1, x]

# extend dataset following RCP unique projection
path = "data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_(16ghg).rcp85_EODS.csv"
TMP, lgd = load_data_and_header(path)
for x in range(len(lgd)):
    EODScmip5[306: ind_cdiac + 1, 0, 0, kindGHG_index["ODS"], 0, ODS.index(lgd[x])] = TMP[6: min(ind_final, ind_cdiac) - 300 + 1, x]

# ========
# 3.4. RCP
# ========

# initialization of projected drivers
EHFCproj = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_HFC], dtype=dty)
EPFCproj = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_PFC], dtype=dty)
EODSproj = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_ODS], dtype=dty)

# projection of emissions under RCP scenarios
# from [Meinshausen et al., 2011]

# HFCs
if (scen_Ehalo[:3] == "RCP") & (ind_final > ind_cdiac):
    for VAR in HFC:
        path = f"data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_5reg0.rcp{scen_Ehalo[3]}{scen_Ehalo[5]}_E{VAR}.csv"

        if os.path.isfile(path):
            TMP = load_data(path)
            for i in range(4 + 1):
                if (mod_regionI == "RCP5") & (mod_regionJ == "RCP5"):
                    EHFCproj[300: min(ind_final, 400) + 1, i, 0, kindGHG_index["HFC"], i, HFC.index(VAR)] += TMP[: min(ind_final, 400) - 300 + 1, i]
                elif (mod_regionI == "RCP5") & (mod_regionJ != "RCP5"):
                    EHFCproj[300: min(ind_final, 400) + 1, 0, 0, kindGHG_index["HFC"], i, HFC.index(VAR)] += TMP[: min(ind_final, 400) - 300 + 1, i]
                elif (mod_regionI != "RCP5") & (mod_regionJ == "RCP5"):
                    EHFCproj[300: min(ind_final, 400) + 1, i, 0, kindGHG_index["HFC"], 0, HFC.index(VAR)] += TMP[: min(ind_final, 400) - 300 + 1, i]
                else:
                    EHFCproj[300: min(ind_final, 400) + 1, 0, 0, kindGHG_index["HFC"], 0, HFC.index(VAR)] += TMP[: min(ind_final, 400) - 300 + 1, i]

# PFCs
if (scen_Ehalo[:3] == "RCP") & (ind_final > ind_cdiac):
    for VAR in PFC:
        path = f"data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_5reg0.rcp{scen_Ehalo[3]}{scen_Ehalo[5]}_E{VAR}.csv"

        if os.path.isfile(path):
            TMP = load_data(path)
            for i in range(4 + 1):
                if (mod_regionI == "RCP5") & (mod_regionJ == "RCP5"):
                    EPFCproj[300: min(ind_final, 400) + 1, i, 0, kindGHG_index["PFC"], i, PFC.index(VAR)] += TMP[: min(ind_final, 400) - 300 + 1, i]
                elif (mod_regionI == "RCP5") & (mod_regionJ != "RCP5"):
                    EPFCproj[300: min(ind_final, 400) + 1, 0, 0, kindGHG_index["PFC"], i, PFC.index(VAR)] += TMP[: min(ind_final, 400) - 300 + 1, i]
                elif (mod_regionI != "RCP5") & (mod_regionJ == "RCP5"):
                    EPFCproj[300: min(ind_final, 400) + 1, i, 0, kindGHG_index["PFC"], 0, PFC.index(VAR)] += TMP[: min(ind_final, 400) - 300 + 1, i]
                else:
                    EPFCproj[300: min(ind_final, 400) + 1, 0, 0, kindGHG_index["PFC"], 0, PFC.index(VAR)] += TMP[: min(ind_final, 400) - 300 + 1, i]

# ODSs
if (scen_Ehalo[:3] == "RCP") & (ind_final > ind_cdiac):
    path = f"data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_(16ghg).rcp{scen_Ehalo[3]}{scen_Ehalo[5]}_EODS.csv"
    TMP, lgd = load_data_and_header(path)
    for x in range(len(lgd)):
        EODSproj[300: min(ind_final, 400) + 1, 0, 0, kindGHG_index["ODS"], 0, ODS.index(lgd[x])] = TMP[: min(ind_final, 400) - 300 + 1, x]

# =================
# 3.A. PAST DATASET
# =================

# datasets mixed following trends

# with EDGAR as reference
if data_Ehalo == "EDGAR":
    EHFCpast = EHFCedgar.copy()
    # follow EDGAR-FT variations after 2008
    for t in range(ind_edgar + 1, ind_cdiac + 1):
        EHFCpast[t, ...] = EHFCpast[t - 1, ...] * EHFCeft[t, ...] / EHFCeft[t - 1, ...]
        EHFCpast[np.isnan(EHFCpast) | np.isinf(EHFCpast)] = 0
        EHFCpast[t, ...] *= np.sum(EHFCpast[t - 1, ...]) / np.sum(EHFCpast[t, ...]) * np.sum(EHFCeft[t, ...]) / np.sum(EHFCeft[t - 1, ...])
    EHFCpast[np.isnan(EHFCpast) | np.isinf(EHFCpast)] = 0
    # quadratic extrapolation before 1970
    # starting year of emission based on [Meinshausen et al., 2011]
    for x in range(nb_HFC):
        if HFC[x] == "HFC23":
            for t in range(230, 270):
                EHFCpast[t, ..., x] = EHFCpast[270, ..., x] * ((t - 230) / (270 - 230.)) ** 2

VAR = "PFC"
# with EDGAR as reference
if data_Ehalo == "EDGAR":
    EPFCpast = EPFCedgar.copy()
    # follow EDGAR-FT variations after 2008
    for t in range(ind_edgar + 1, ind_cdiac + 1):
        EPFCpast[t, ...] = EPFCpast[t - 1, ...] * EPFCeft[t, ...] / EPFCeft[t - 1, ...]
        EPFCpast[np.isnan(EPFCpast) | np.isinf(EPFCpast)] = 0
        EPFCpast[t, ...] *= np.sum(EPFCpast[t - 1, ...]) / np.sum(EPFCpast[t, ...]) * np.sum(EPFCeft[t, ...]) / np.sum(EPFCeft[t - 1, ...])
    EPFCpast[np.isnan(EPFCpast) | np.isinf(EPFCpast)] = 0
    # quadratic extrapolation before 1970
    # starting year of emission based on [Meinshausen et al., 2011]
    for x in range(nb_PFC):
        if PFC[x] == "SF6":
            for t in range(250, 270):
                EPFCpast[t, ..., x] = EPFCpast[270, ..., x] * ((t - 250.) / (270 - 250.)) ** 2
        elif PFC[x] == "CF4":
            for t in range(222, 270):
                EPFCpast[t, ..., x] = EPFCpast[270, ..., x] * ((t - 222.) / (270 - 222.)) ** 2
        elif PFC[x] == "C2F6":
            for t in range(189, 270):
                EPFCpast[t, ..., x] = EPFCpast[270, ..., x] * ((t - 189) / (270 - 189.)) ** 2

# with CMIP5 as reference
if data_Ehalo == data_Ehalo:  # FIXME: is this so?
    EODSpast = EODScmip5.copy()
    # linear extrapolation before 1765
    for t in range(50, 65):
        EODSpast[t, ...] = EODSpast[65, ...] * (t - 50) / float(65 - 50)

# cut past dataset to right length
EHFC[:min(ind_cdiac, ind_final) + 1, ...] = EHFCpast[:min(ind_cdiac, ind_final) + 1, ...]
EPFC[:min(ind_cdiac, ind_final) + 1, ...] = EPFCpast[:min(ind_cdiac, ind_final) + 1, ...]
EODS[:min(ind_cdiac, ind_final) + 1, ...] = EODSpast[:min(ind_cdiac, ind_final) + 1, ...]

# ==================
# 3.B. FINAL DATASET
# ==================

# datasets mixed following various criteria
for arr, arrproj in [
    (EHFC, EHFCproj), (EPFC, EPFCproj), (EODS, EODSproj),
]:

    # stop emissions
    if (scen_Ehalo == "stop") & (ind_final > ind_cdiac):
        arr[ind_cdiac + 1:, ...] = 0

    # constant emissions
    elif (scen_Ehalo == "cst") & (ind_final > ind_cdiac):
        arr[ind_cdiac + 1:, ...] = arr[ind_cdiac, ...][np.newaxis, ...]

        # RCP scenarios
    elif (scen_Ehalo[:3] == "RCP") & (ind_final > ind_cdiac):

        # raw discontinuity
        if mod_DATAscen == "raw":
            arr[ind_cdiac + 1:, ...] = arrproj[ind_cdiac + 1:, ...]

        # offset at transition point
        elif mod_DATAscen == "offset":
            arr[ind_cdiac + 1:, ...] = arrproj[ind_cdiac + 1:, ...] - arrproj[ind_cdiac, ...] + arr[ind_cdiac, ...]
            for t in range(ind_cdiac + 1, ind_final + 1):
                def_regI = bool(np.sum(arrproj[t, :, ..., 1:]))
                def_regJ = bool(np.sum(arrproj[t, 1:, ..., :]))
                if not def_regI:
                    arr[t, :, ..., 0] += np.sum(arr[t, :, ..., 1:], -1)
                    arr[t, :, ..., 1:] = 0
                if not def_regJ:
                    arr[t, 0, ..., :] += np.sum(arr[t, 1:, ..., :], 0)
                    arr[t, 1:, ..., :] = 0

                    # linear transition over N years
        elif mod_DATAscen[:6] == "smooth":
            N = int(mod_DATAscen[6:])
            if ind_final >= ind_cdiac + N:
                for t in range(ind_cdiac + 1, ind_cdiac + N):
                    arr[t, ...] = (1 - (t - ind_cdiac) / float(N)) * arr[ind_cdiac, ...] + (t - ind_cdiac) / float(N) * arrproj[ind_cdiac + N, ...]
                    def_regI = bool(np.sum(arrproj[t, :, ..., 1:]))
                    def_regJ = bool(np.sum(arrproj[t, 1:, ..., :]))
                    if not def_regI:
                        arr[t, :, ..., 0] += np.sum(arr[t, :, ..., 1:], -1)
                        arr[t, :, ..., 1:] = 0
                    if not def_regJ:
                        arr[t, 0, ..., :] += np.sum(arr[t, 1:, ..., :], 0)
                        arr[t, 1:, ..., :] = 0
                arr[ind_cdiac + N:, ...] = arrproj[ind_cdiac + N:, ...]

        # follow trends of projection
        elif mod_DATAscen == "trends":
            for t in range(ind_cdiac + 1, ind_final + 1):
                def_regI = bool(np.sum(arrproj[t, :, ..., 1:]))
                def_regJ = bool(np.sum(arrproj[t, 1:, ..., :]))
                if def_regI and def_regJ:
                    arr[t, ...] = arr[t - 1, ...] * arrproj[t, ...] / arrproj[t - 1, ...]
                    arr[np.isnan(arr) | np.isinf(arr)] = 0
                    arr[t, ...] *= np.sum(arr[t - 1, ...]) / np.sum(arr[t, ...]) * np.sum(arrproj[t, ...]) / np.sum(arrproj[t - 1, ...])
                elif not def_regI and def_regJ:
                    arr[t, :, ..., 0] = np.sum(arr[t - 1, :, ..., :], 0) * np.sum(arrproj[t, :, ..., :], -1) / np.sum(arrproj[t - 1, :, ..., :], -1)
                    arr[np.isnan(arr) | np.isinf(arr)] = 0
                    arr[t, ...] *= np.sum(arr[t - 1, ...]) / np.sum(arr[t, ...]) * np.sum(arrproj[t, ...]) / np.sum(arrproj[t - 1, ...])
                elif def_regI and not def_regJ:
                    arr[t, 0, ..., :] = np.sum(arr[t - 1, :, ..., :], 0) * np.sum(arrproj[t, :, ..., :], 0) / np.sum(arrproj[t - 1, :, ..., :], 0)
                    arr[np.isnan(arr) | np.isinf(arr)] = 0
                    arr[t, ...] *= np.sum(arr[t - 1, ...]) / np.sum(arr[t, ...]) * np.sum(arrproj[t, ...]) / np.sum(arrproj[t - 1, ...])
                elif not def_regI and not def_regJ:
                    arr[t, 0, ..., 0] = np.sum(np.sum(arr[t - 1, :, ..., :], -1), 0) * np.sum(np.sum(arrproj[t, :, ..., :], -1), 0) / np.sum(np.sum(arrproj[t - 1, :, ..., :], -1), 0)
            arr[np.isnan(arr) | np.isinf(arr)] = 0

# delete individual datasets
del EHFCedgar, EHFCeft, EHFCpast, EHFCproj
del EPFCedgar, EPFCeft, EPFCpast, EPFCproj
del EODScmip5, EODSpast, EODSproj
