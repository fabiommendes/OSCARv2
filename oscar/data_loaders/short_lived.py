import os

import numpy as np

from .regions import nb_regionJ, nb_kind, nb_regionI, nb_sector, regionJ_index, regionI_index, ind_final, kindAER_index, kindCHI_index
from .greenhouse import sec_accmip, ind_edgar, ind_cdiac
from ..data import load_data
from ..config import dty, mod_regionI, mod_regionJ, scen_ESO2, mod_DATAscen, scen_ECO, scen_ENOX, scen_EVOC, scen_ENH3, scen_EBC, scen_EOC, data_ENOX, data_ECO, data_EVOC, data_ESO2, data_ENH3, data_EOC, data_EBC

##################################################
#   4. SHORT-LIVED SPECIES
##################################################
ENOX = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)  # {TgN/yr}
ECO = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)  # {TgC/yr}
EVOC = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)  # {Tg/yr}
ESO2 = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)  # {TgS/yr}
ENH3 = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)  # {TgN/yr}
EOC = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)  # {Tg/yr}
EBC = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)  # {Tg/yr}

# ==========
# 4.1. EDGAR
# ==========

# load emissions from EDGAR v4.2
# see [JRC, 2011]

# OzoPrec
for VAR in ["NOX", "CO", "VOC"]:
    arr = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
    if VAR == "NOX":
        ENOXedgar = arr
    elif VAR == "CO":
        ECOedgar = arr
    elif VAR == "VOC":
        EVOCedgar = arr
    else:
        raise RuntimeError

    for s in range(1, len(sec_accmip) - 2):
        path = f"data/EOzoPrec_EDGAR/#DATA.EOzoPrec_EDGAR.1970-{1700 + ind_edgar}_114reg0.E{VAR}_{sec_accmip[s]}.csv"
        if os.path.isfile(path):
            TMP = load_data(path)
            for i in range(114 + 1):
                arr[270:ind_edgar + 1, regionJ_index[i], 0, kindCHI_index[VAR], regionI_index[i]] += TMP[:ind_edgar - 270 + 1, i]

# AeroPrec
for VAR in ["SO2", "NH3"]:
    arr = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)

    if VAR == 'SO2':
        ESO2edgar = arr
    elif VAR == 'NH3':
        ENH3edgar = arr
    else:
        raise RuntimeError

    for s in range(1, len(sec_accmip) - 2):
        path = f"data/EAeroPrec_EDGAR/#DATA.EAeroPrec_EDGAR.1970-{1700 + ind_edgar}_114reg0.E{VAR}_{sec_accmip[s]}.csv"
        if os.path.isfile(path):
            TMP = load_data(path)
            for i in range(114 + 1):
                arr[270:ind_edgar + 1, regionJ_index[i], 0, kindAER_index[VAR], regionI_index[i]] += TMP[:ind_edgar - 270 + 1, i]

# PM10 (proxy for OC/BC)
EPM10edgar = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
for s in range(1, len(sec_accmip) - 2):
    path = f"data/EAeroPrec_EDGAR/#DATA.EAeroPrec_EDGAR.1970-{1700 + ind_edgar}_114reg0.EPM10_{sec_accmip[s]}.csv"
    if os.path.isfile(path):
        TMP = load_data(path)
        for i in range(114 + 1):
            EPM10edgar[270: ind_edgar + 1, regionJ_index[i], 0, kindAER_index["OC"], regionI_index[i]] += TMP[: ind_edgar - 270 + 1, i] / 2
            EPM10edgar[270: ind_edgar + 1, regionJ_index[i], 0, kindAER_index["BC"], regionI_index[i]] += TMP[: ind_edgar - 270 + 1, i] / 2

# ===============
# 4.2. EDGAR-HTAP
# ===============

# load emissions from EDGAR-HTAP v2
# see [JRC, 2013]

# OzoPrec
for VAR in ["NOX", "CO", "VOC"]:
    arr = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
    if VAR == 'NOX':
        ENOXehtap = arr
    elif VAR == 'CO':
        ECOehtap = arr
    elif VAR == 'VOC':
        EVOCehtap = arr

    for s in range(1, len(sec_accmip) - 2):
        path = f"data/EOzoPrec_EDGAR-HTAP/#DATA.EOzoPrec_EDGAR-HTAP.2008-2010_114reg0.E{VAR}_{sec_accmip[s]}.csv"
        if os.path.isfile(path):
            TMP = load_data(path)
            for i in range(114 + 1):
                arr[308:310 + 1, regionJ_index[i], 0, kindCHI_index[VAR], regionI_index[i]] += TMP[:ind_edgar - 270 + 1, i]

# AeroPrec
for VAR in ["SO2", "NH3", "OC", "BC"]:
    arr = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)

    if VAR == "SO2":
        ESO2ehtap = arr
    elif VAR == "NH3":
        ENH3ehtap = arr
    elif VAR == "OC":
        EOCehtap = arr
    elif VAR == "BC":
        EBCehtap = arr

    for s in range(1, len(sec_accmip) - 2):
        path = f"data/EAeroPrec_EDGAR-HTAP/#DATA.EAeroPrec_EDGAR-HTAP.2008-2010_114reg0.E{VAR}_{sec_accmip[s]}.csv"
        if os.path.isfile(path):
            TMP = load_data(path)
            for i in range(114 + 1):
                arr[308:310 + 1, regionJ_index[i], 0, kindAER_index[VAR], regionI_index[i]] += TMP[:ind_edgar - 270 + 1, i]

# ===========
# 4.3. ACCMIP
# ===========

# load emissions from ACCMIP
# from [Lamarque et al., 2010]

# OzoPrec
for VAR in ["NOX", "CO", "VOC"]:
    accmip = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
    bio = np.zeros([nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)

    if VAR == "NOX":
        ENOXaccmip = accmip
        p_ENOX_bio = bio
    if VAR == "CO":
        ECOaccmip = accmip
        p_ECO_bio = bio
    if VAR == "VOC":
        EVOCaccmip = accmip
        p_EVOC_bio = bio

    for s in range(1, len(sec_accmip) - 2):
        path = f"data/EOzoPrec_ACCMIP/#DATA.EOzoPrec_ACCMIP.1850-2000_114reg0.E{VAR}_{sec_accmip[s]}.csv"
        if os.path.isfile(path):
            TMP = load_data(path)
            for i in range(114 + 1):
                accmip[150:300 + 1, regionJ_index[i], 0, kindCHI_index[VAR], regionI_index[i]] += TMP[:300 - 150 + 1, i]
                if sec_accmip[s] in ["agr", "awb", "wst"]:
                    bio[regionJ_index[i], 0, kindCHI_index[VAR], regionI_index[i]] += TMP[0, i]
    bio /= accmip[150]
    bio[np.isnan(bio) | np.isinf(bio)] = 0

# AeroPrec
for VAR in ["SO2", "NH3", "OC", "BC"]:
    accmip = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
    bio = np.zeros([nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)

    if VAR == "SO2":
        ESO2accmip = accmip
        p_ESO2_bio = bio
    elif VAR == "NH3":
        ENH3accmip = accmip
        p_ENH3_bio = bio
    elif VAR == "OC":
        EOCaccmip = accmip
        p_EOC_bio = bio
    elif VAR == "BC":
        EBCaccmip = accmip
        p_EBC_bio = bio

    for s in range(1, len(sec_accmip) - 2):
        path = f"data/EAeroPrec_ACCMIP/#DATA.EAeroPrec_ACCMIP.1850-2000_114reg0.E{VAR}_{sec_accmip[s]}.csv"
        if os.path.isfile(path):
            TMP = load_data(path)
            for i in range(114 + 1):
                accmip[150:300 + 1, regionJ_index[i], 0, kindAER_index[VAR], regionI_index[i]] += TMP[:300 - 150 + 1, i]
                if sec_accmip[s] in ["agr", "awb", "wst"]:
                    bio[regionJ_index[i], 0, kindAER_index[VAR], regionI_index[i]] += TMP[0, i]
    bio /= accmip[150]
    bio[np.isnan(bio) | np.isinf(bio)] = 0

# =========
# 4.4. SRES
# =========

# initialization of projected drivers
ENOXproj = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
ECOproj = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
EVOCproj = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
ESO2proj = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
ENH3proj = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
EOCproj = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
EBCproj = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)

# projection of emissions under SRES scenarios
# from [IPCC, 2000]

# OzoPrec
for VAR, scen, arr in [
    ("NOX", scen_ENOX, ENOXproj), ("CO", scen_ECO, ECOproj), ("VOC", scen_EVOC, EVOCproj),
]:
    if (scen[:4] == "SRES") & (ind_final > ind_cdiac):
        path = f"data/EOzoPrec_SRES/#DATA.EOzoPrec_SRES.2000-2100_4reg0.{scen[5:]}_E{VAR}.csv"
        TMP = load_data(path)

        for i in range(4 + 1):
            if (mod_regionI == "SRES4") & (mod_regionJ == "SRES4"):
                arr[300:min(ind_final, 400) + 1, i, 0, kindCHI_index[VAR], i] += TMP[:min(ind_final, 400) - 300 + 1, i]
            elif (mod_regionI == "SRES4") & (mod_regionJ != "SRES4"):
                arr[300:min(ind_final, 400) + 1, 0, 0, kindCHI_index[VAR], i] += TMP[:min(ind_final, 400) - 300 + 1, i]
            elif (mod_regionI != "SRES4") & (mod_regionJ == "SRES4"):
                arr[300:min(ind_final, 400) + 1, i, 0, kindCHI_index[VAR], 0] += TMP[:min(ind_final, 400) - 300 + 1, i]
            elif (mod_regionI != "SRES4") & (mod_regionJ != "SRES4"):
                arr[300:min(ind_final, 400) + 1, 0, 0, kindCHI_index[VAR], 0] += TMP[:min(ind_final, 400) - 300 + 1, i]

# AeroPrec (SO2 only)
if (scen_ESO2[:4] == "SRES") & (ind_final > ind_cdiac):
    path = f"data/EAeroPrec_SRES/#DATA.EAeroPrec_SRES.2000-2100_4reg0.{scen_ESO2[5:]}_ESO2.csv"
    TMP = load_data(path)
    for i in range(4 + 1):
        if (mod_regionI == "SRES4") & (mod_regionJ == "SRES4"):
            ESO2proj[300: min(ind_final, 400) + 1, i, 0, kindAER_index["SO2"], i] += TMP[: min(ind_final, 400) - 300 + 1, i]
        elif (mod_regionI == "SRES4") & (mod_regionJ != "SRES4"):
            ESO2proj[300: min(ind_final, 400) + 1, 0, 0, kindAER_index["SO2"], i] += TMP[: min(ind_final, 400) - 300 + 1, i]
        elif (mod_regionI != "SRES4") & (mod_regionJ == "SRES4"):
            ESO2proj[300: min(ind_final, 400) + 1, i, 0, kindAER_index["SO2"], 0] += TMP[: min(ind_final, 400) - 300 + 1, i]
        elif (mod_regionI != "SRES4") & (mod_regionJ != "SRES4"):
            ESO2proj[300: min(ind_final, 400) + 1, 0, 0, kindAER_index["SO2"], 0] += TMP[: min(ind_final, 400) - 300 + 1, i]

# ========
# 4.5. RCP
# ========

# projection of emissions under RCP scenarios
# from [Meinshausen et al., 2011]

# OzoPrec
for VAR, scen, arr in [
    ("NOX", scen_ENOX, ENOXproj), ("CO", scen_ECO, ECOproj), ("VOC", scen_EVOC, EVOCproj),
]:
    if (scen[:3] == "RCP") & (ind_final > ind_cdiac):
        for s in range(1, len(sec_accmip) - 2):
            path = f"data/EOzoPrec_RCP/#DATA.EOzoPrec_RCP.2000-2100_5reg0.rcp{scen[3]}{scen[5]}_E{VAR}_{sec_accmip[s]}.csv"
            if os.path.isfile(path):
                TMP = load_data(path)
                for i in range(5 + 1):
                    if (mod_regionI == "RCP5") & (mod_regionJ == "RCP5"):
                        arr[300:min(ind_final, 400) + 1, i, 0, kindCHI_index[VAR], i] += TMP[:min(ind_final, 400) - 300 + 1, i]
                    elif (mod_regionI == "RCP5") & (mod_regionJ != "RCP5"):
                        arr[300:min(ind_final, 400) + 1, 0, 0, kindCHI_index[VAR], i] += TMP[:min(ind_final, 400) - 300 + 1, i]
                    elif (mod_regionI != "RCP5") & (mod_regionJ == "RCP5"):
                        arr[300:min(ind_final, 400) + 1, i, 0, kindCHI_index[VAR], 0] += TMP[:min(ind_final, 400) - 300 + 1, i]
                    elif (mod_regionI != "RCP5") & (mod_regionJ != "RCP5"):
                        arr[300:min(ind_final, 400) + 1, 0, 0, kindCHI_index[VAR], 0] += TMP[:min(ind_final, 400) - 300 + 1, i]

# AeroPrec
for VAR, scen, arr in [
    ("SO2", scen_ESO2, ESO2proj), ("NH3", scen_ENH3, ENH3proj), ("OC", scen_EOC, EOCproj), ("BC", scen_EBC, EBCproj),
]:
    if (scen[:3] == "RCP") & (ind_final > ind_cdiac):
        for s in range(1, len(sec_accmip) - 2):
            path = f"data/EAeroPrec_RCP/#DATA.EAeroPrec_RCP.2000-2100_5reg0.rcp{scen[3]}{scen[5]}_E{VAR}_" f"{sec_accmip[s]}.csv"

            if os.path.isfile(path):
                TMP = load_data(path)
                for i in range(5 + 1):
                    if (mod_regionI == "RCP5") & (mod_regionJ == "RCP5"):
                        arr[300:min(ind_final, 400) + 1, i, 0, kindAER_index[VAR], i] += TMP[:min(ind_final, 400) - 300 + 1, i]
                    elif (mod_regionI == "RCP5") & (mod_regionJ != "RCP5"):
                        arr[300:min(ind_final, 400) + 1, 0, 0, kindAER_index[VAR], i] += TMP[:min(ind_final, 400) - 300 + 1, i]
                    elif (mod_regionI != "RCP5") & (mod_regionJ == "RCP5"):
                        arr[300:min(ind_final, 400) + 1, i, 0, kindAER_index[VAR], 0] += TMP[:min(ind_final, 400) - 300 + 1, i]
                    elif (mod_regionI != "RCP5") & (mod_regionJ != "RCP5"):
                        arr[300:min(ind_final, 400) + 1, 0, 0, kindAER_index[VAR], 0] += TMP[:min(ind_final, 400) - 300 + 1, i]

# =================
# 4.A. PAST DATASET
# =================

# datasets mixed following trends
for VAR, data, edgar, raw, ehtap, accmip, bio in [
    ("NOX", data_ENOX, ENOXedgar, ENOXedgar if data_ENOX == "EDGAR" else ENOXaccmip, ENOXehtap, ENOXaccmip, p_ENOX_bio), ("CO", data_ECO, ECOedgar, ECOedgar if data_ECO == "EDGAR" else ECOaccmip, ECOehtap, ECOaccmip, p_ECO_bio), ("VOC", data_EVOC, EVOCedgar, EVOCedgar if data_EVOC == "EDGAR" else EVOCaccmip, EVOCehtap, EVOCaccmip, p_EVOC_bio), ("SO2", data_ESO2, ESO2edgar, ESO2edgar if data_ESO2 == "EDGAR" else ESO2accmip, ESO2ehtap, ESO2accmip, p_ESO2_bio), ("NH3", data_ENH3, ENH3edgar, ENH3edgar if data_ENH3 == "EDGAR" else ENH3accmip, ENH3ehtap, ENH3accmip, p_ENH3_bio),
]:
    # with EDGAR as reference
    past = raw.copy()
    if data == "EDGAR":
        # follow EDGAR-HTAP variations after 2008
        for t in range(ind_edgar + 1, ind_cdiac + 1):
            past[t, ...] = past[t - 1, ...] * ehtap[t, ...] / ehtap[t - 1, ...]
            past[np.isnan(past) | np.isinf(past)] = 0
            past[t, ...] *= np.sum(past[t - 1, ...]) / np.sum(past[t, ...]) * np.sum(ehtap[t, ...]) / np.sum(ehtap[t - 1, ...])
        past[np.isnan(past) | np.isinf(past)] = 0

        # follow ACCMIP variations before 1970
        for t in range(150, 270)[::-1]:
            past[t, ...] = past[t + 1, ...] * accmip[t, ...] / accmip[t + 1, ...]
            past[np.isnan(past) | np.isinf(past)] = 0
            past[t, ...] *= np.sum(past[t + 1, ...]) / np.sum(past[t, ...]) * np.sum(accmip[t, ...]) / np.sum(accmip[t + 1, ...])
        past[np.isnan(past) | np.isinf(past)] = 0

        # linear extrapolation before 1850
        for t in range(50, 150):
            past[t, ...] = (1 - bio) * past[150, ...] * (t - 50) / float(150 - 50)
        for t in range(0, 150):
            past[t, ...] += bio * past[150, ...] * (t - -200) / float(150 - -200)
        arr_0 = past[50, ...]

    # with ACCMIP as reference
    elif data == "ACCMIP":
        # follow EDGAR variations after 2000
        for t in range(300 + 1, ind_edgar + 1):
            past[np.isnan(past) | np.isinf(past)] = 0
            past[t, ...] = past[t - 1, ...] * edgar[t, ...] / edgar[t - 1, ...]
            past[t, ...] *= np.sum(past[t - 1, ...]) / np.sum(past[t, ...]) * np.sum(edgar[t, ...]) / np.sum(edgar[t - 1, ...])
        past[np.isnan(past) | np.isinf(past)] = 0

        # follow EDGAR-HTAP variations after 2008
        for t in range(ind_edgar + 1, ind_cdiac + 1):
            past[t, ...] = past[t - 1, ...] * ehtap[t, ...] / ehtap[t - 1, ...]
            past[np.isnan(past) | np.isinf(past)] = 0
            past[t, ...] *= np.sum(past[t - 1, ...]) / np.sum(past[t, ...]) * np.sum(ehtap[t, ...]) / np.sum(ehtap[t - 1, ...])
        past[np.isnan(past) | np.isinf(past)] = 0

        # linear extrapolation before 1850
        for t in range(50, 150):
            past[t, ...] = (1 - bio) * past[150, ...] * (t - 50) / float(150 - 50)
        for t in range(0, 150):
            past[t, ...] += bio * past[150, ...] * (t - -200) / float(150 + 200)
        arr_0 = arr[50, ...]
    else:
        raise ValueError(data)

    if VAR == "NOX":
        ENOXpast = past
        ENOX_0 = arr_0
    elif VAR == "CO":
        ECOpast = past
        ECO_0 = arr_0
    elif VAR == "VOC":
        EVOCpast = past
        EVOC_0 = arr_0
    elif VAR == "SO2":
        ESO2past = past
        ESO2_0 = arr_0
    elif VAR == "NH3":
        ENH3past = past
        ENH3_0 = arr_0

for VAR, data, bio, accmip, ehtap in [
    ("OC", data_EOC, p_EOC_bio, EOCaccmip, EOCehtap), ("BC", data_EBC, p_EBC_bio, EBCaccmip, EBCehtap),
]:
    # with ACCMIP as reference
    if data == "ACCMIP":
        past = accmip.copy()

        # follow EDGAR (PM10) variations after 2000
        for t in range(300 + 1, ind_edgar + 1):
            past[t, ...] = past[t - 1, ...] * EPM10edgar[t, ...] / EPM10edgar[t - 1, ...]
            past[np.isnan(past) | np.isinf(past)] = 0
            past[t, ...] *= np.sum(past[t - 1, ...]) / np.sum(past[t, ...]) * np.sum(EPM10edgar[t, ...]) / np.sum(EPM10edgar[t - 1, ...])
        past[np.isnan(past) | np.isinf(past)] = 0

        # follow EDGAR-HTAP variations after 2008
        for t in range(ind_edgar + 1, ind_cdiac + 1):
            past[t, ...] = past[t - 1, ...] * ehtap[t, ...] / ehtap[t - 1, ...]
            past[np.isnan(past) | np.isinf(past)] = 0
            past[t, ...] *= np.sum(past[t - 1, ...]) / np.sum(past[t, ...]) * np.sum(ehtap[t, ...]) / np.sum(ehtap[t - 1, ...])
        past[np.isnan(past) | np.isinf(past)] = 0

        # linear extrapolation before 1850
        for t in range(50, 150):
            past[t, ...] = (1 - bio) * past[150, ...] * (t - 50) / float(150 - 50)
        for t in range(0, 150):
            past[t, ...] += bio * past[150, ...] * (t - -200) / float(150 - -200)
        arr_0 = past[50, ...]
    else:
        past = arr_0 = None

    if VAR == "OC":
        EOCpast = past
        EOC_0 = arr_0
    elif VAR == "BC":
        EBCpast = past
        EBC_0 = arr_0

# cut past dataset to right length
ENOX[:min(ind_cdiac, ind_final) + 1, ...] = ENOXpast[:min(ind_cdiac, ind_final) + 1, ...]
ECO[:min(ind_cdiac, ind_final) + 1, ...] = ECOpast[:min(ind_cdiac, ind_final) + 1, ...]
EVOC[:min(ind_cdiac, ind_final) + 1, ...] = EVOCpast[:min(ind_cdiac, ind_final) + 1, ...]
ESO2[:min(ind_cdiac, ind_final) + 1, ...] = ESO2past[:min(ind_cdiac, ind_final) + 1, ...]
ENH3[:min(ind_cdiac, ind_final) + 1, ...] = ENH3past[:min(ind_cdiac, ind_final) + 1, ...]
EOC[:min(ind_cdiac, ind_final) + 1, ...] = EOCpast[:min(ind_cdiac, ind_final) + 1, ...]
EBC[:min(ind_cdiac, ind_final) + 1, ...] = EBCpast[:min(ind_cdiac, ind_final) + 1, ...]

# ==================
# 4.B. FINAL DATASET
# ==================

# datasets mixed following various criteria
for VAR, scen, arr, arrproj, arr_0 in [
    ("NOX", scen_ENOX, ENOX, ENOXproj, ENOX_0), ("CO", scen_ECO, ECO, ECOproj, ECO_0), ("VOC", scen_EVOC, EVOC, EVOCproj, EVOC_0), ("SO2", scen_ESO2, ESO2, ESO2proj, ESO2_0), ("NH3", scen_ENH3, ENH3, ENH3proj, ENH3_0), ("OC", scen_EOC, EOC, EOCproj, EOC_0), ("BC", scen_EBC, EBC, EBCproj, EBC_0),
]:
    # stop emissions
    if (scen == "stop") & (ind_final > ind_cdiac):
        arr[ind_cdiac + 1:, ...] = arr_0[np.newaxis, ...]

    # constant emissions
    elif (scen == "cst") & (ind_final > ind_cdiac):
        arr[ind_cdiac + 1:, ...] = arr[ind_cdiac, ...][np.newaxis, ...]

    # RCP or SRES scenarios
    elif ((scen[:4] == "SRES") | (scen[:3] == "RCP")) & (ind_final > ind_cdiac):

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
                    arr[t, :, ..., 0] = np.sum(arr[t - 1, :, ..., :], -1) * np.sum(arrproj[t, :, ..., :], -1) / np.sum(arrproj[t - 1, :, ..., :], -1)
                    arr[np.isnan(arr) | np.isinf(arr)] = 0
                    arr[t, ...] *= np.sum(arr[t - 1, ...]) / np.sum(arr[t, ...]) * np.sum(arrproj[t, ...]) / np.sum(arrproj[t - 1, ...])
                elif def_regI and not def_regJ:
                    arr[t, 0, ..., :] = np.sum(arr[t - 1, :, ..., :], 0) * np.sum(arrproj[t, :, ..., :], 0) / np.sum(arrproj[t - 1, :, ..., :], 0)
                    arr[np.isnan(arr) | np.isinf(arr)] = 0
                    arr[t, ...] *= np.sum(arr[t - 1, ...]) / np.sum(arr[t, ...]) * np.sum(arrproj[t, ...]) / np.sum(arrproj[t - 1, ...])
                elif not def_regI and not def_regJ:
                    arr[t, 0, ..., 0] = np.sum(np.sum(arr[t - 1, :, ..., :], -1), 0) * np.sum(np.sum(arrproj[t, :, ..., :], -1), 0) / np.sum(np.sum(arrproj[t - 1, :, ..., :], -1), 0)
            arr[np.isnan(arr) | np.isinf(arr)] = 0

# Delete individual datasets
del ENOXedgar, ENOXehtap, ENOXaccmip, ENOXpast, ENOXproj, p_ENOX_bio
del ECOedgar, ECOehtap, ECOaccmip, ECOpast, ECOproj, p_ECO_bio
del EVOCedgar, EVOCehtap, EVOCaccmip, EVOCpast, EVOCproj, p_EVOC_bio
del ESO2edgar, ESO2ehtap, ESO2accmip, ESO2past, ESO2proj, p_ESO2_bio
del ENH3edgar, ENH3ehtap, ENH3accmip, ENH3past, ENH3proj, p_ENH3_bio
del EOCehtap, EOCaccmip, EOCpast, EOCproj, p_EOC_bio
del EBCehtap, EBCaccmip, EBCpast, EBCproj, p_EBC_bio
del EPM10edgar
