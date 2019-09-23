import os

import numpy as np

from .regions import nb_regionJ, nb_kind, nb_regionI, nb_sector, regionJ_index, kFF, regionI_index, kindGHG_index
from ..data import load_data, load_data_and_header
from .. import conf

##################################################
#   1. GREENHOUSE GASES
##################################################

init = lambda: np.zeros([conf.ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=conf.dty)
EFF = init()  # {GtC/yr}
ECH4 = init()  # {TgC/yr}
EN2O = init()  # {TgN/yr}

# ==========
# 1.1. CDIAC
# ==========

# load CDIAC emissions
# from [Boden et al., 2012]
ind_cdiac = 310
kin = ["sol", "liq", "gas", "cem", "fla"]

# total emissions
EFFcdiac = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=conf.dty)

for k in range(len(kin)):
    path = f"data/EFossil_CDIAC/#DATA.EFossil_CDIAC.1751-2010_114reg0.EFF_{kin[k]}.csv"
    TMP = load_data(path)
    for i in range(114 + 1):
        EFFcdiac[51: ind_cdiac + 1, regionJ_index[i], 0, kFF, regionI_index[i]] += TMP[: ind_cdiac - 51 + 1, i]

# distribution among fuels
# TODO

# ========
# 1.2. EPA
# ========

# load EPA emissions
# see [EPA, 2012]
ind_epa = 310
sec_epa = ["agr_ferm", "agr_manu", "agr_othr", "agr_rice", "agr_soil", "ene_burn", "ene_coal", "ene_comb", "ene_ngos", "ene_othr", "ind_acid"]
sec_epa1 = ["agr_ferm", "agr_manu", "agr_othr", "agr_rice", "agr_soil", "ene_coal", "ene_comb", "ene_ngos", "ene_othr", "ind_acid"]

# CH4
ECH4epa = np.zeros([conf.ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=conf.dty)
for s in range(len(sec_epa1)):
    path = f"data/EMethane_EPA/#DATA.EMethane_EPA.1990-{1700 + ind_epa}_114reg0.ECH4_{sec_epa1[s]}.csv"
    if os.path.isfile(path):
        TMP = load_data(path)
        for i in range(114 + 1):
            ECH4epa[290: ind_epa + 1, regionJ_index[i], 0, kindGHG_index["CH4"], regionI_index[i]] += TMP[: ind_epa - 290 + 1, i]

# N2O
EN2Oepa = np.zeros([conf.ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=conf.dty)
for s in range(len(sec_epa1)):
    path = f"data/ENitrousOx_EPA/#DATA.ENitrousOx_EPA.1990-{1700 + ind_epa}_114reg0.EN2O_{sec_epa1[s]}.csv"
    if os.path.isfile(path):
        TMP = load_data(path)
        for i in range(114 + 1):
            EN2Oepa[290: ind_epa + 1, regionJ_index[i], 0, kindGHG_index["N2O"], regionI_index[i]] += TMP[: ind_epa - 290 + 1, i]

# ==========
# 1.3. EDGAR
# ==========

# load emissions from EDGAR v4.2
# see [JRC, 2011]
ind_edgar = 308
sec_ehyde = ["oo", "fc", "fp", "bc", "in", "al", "an", "aw", "lf", "sb", "df"]
sec_accmip = ["ooo", "ene", "ind", "tra", "shp", "air", "dom", "slv", "agr", "awb", "wst", "for", "gra"]

# FF
EFFedgar = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=conf.dty)
for s in range(1, len(sec_ehyde) - 2):
    path = f"data/ECarbon_EDGAR/#DATA.ECarbon_EDGAR.1970-{1700 + ind_edgar}_114reg0.ECO2_{sec_ehyde[s]}.csv"
    if os.path.isfile(path):
        TMP = load_data(path)
        for i in range(114 + 1):
            EFFedgar[270: ind_edgar + 1, regionJ_index[i], 0, kFF, regionI_index[i]] += TMP[: ind_edgar - 270 + 1, i]

# CH4
ECH4edgar = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=conf.dty)
for s in range(1, len(sec_accmip) - 2):
    path = f"data/EMethane_EDGAR/#DATA.EMethane_EDGAR.1970-{1700 + ind_edgar}_114reg0.ECH4_{sec_accmip[s]}.csv"
    if os.path.isfile(path):
        TMP = load_data(path)
        for i in range(114 + 1):
            ECH4edgar[270: ind_edgar + 1, regionJ_index[i], 0, kindGHG_index["CH4"], regionI_index[i]] += TMP[: ind_edgar - 270 + 1, i]

# N2O
EN2Oedgar = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=conf.dty)
for s in range(1, len(sec_ehyde) - 2):
    path = f"data/ENitrousOx_EDGAR/#DATA.ENitrousOx_EDGAR.1970-{1700 + ind_edgar}_114reg0.EN2O_{sec_ehyde[s]}.csv"
    if os.path.isfile(path):
        TMP = load_data(path)
        for i in range(114 + 1):
            EN2Oedgar[270: ind_edgar + 1, regionJ_index[i], 0, kindGHG_index["N2O"], regionI_index[i]] += TMP[: ind_edgar - 270 + 1, i]

# =============
# 1.4. EDGAR-FT
# =============

# load emissions from EDGAR-FT v4.2-FT2010
# see [JRC, 2013]

# FF
EFFeft = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=conf.dty)
for s in range(1, len(sec_ehyde) - 2):
    path = f"data/ECarbon_EDGAR-FT/#DATA.ECarbon_EDGAR-FT.2008-2010_114reg0.ECO2_{sec_ehyde[s]}.csv"
    if os.path.isfile(path):
        TMP = load_data(path)
        for i in range(114 + 1):
            EFFeft[308: 310 + 1, regionJ_index[i], 0, kFF, regionI_index[i]] += TMP[: 310 - 308 + 1, i]

# CH4
ECH4eft = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=conf.dty)
for s in range(1, len(sec_accmip) - 2):
    path = f"data/EMethane_EDGAR-FT/#DATA.EMethane_EDGAR-FT.2008-2010_114reg0.ECH4_{sec_accmip[s]}.csv"
    if os.path.isfile(path):
        TMP = load_data(path)
        for i in range(114 + 1):
            ECH4eft[308: 310 + 1, regionJ_index[i], 0, kindGHG_index["CH4"], regionI_index[i]] += TMP[: 310 - 308 + 1, i]

# N2O
EN2Oeft = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=conf.dty)
for s in range(1, len(sec_ehyde) - 2):
    path = f"data/ENitrousOx_EDGAR-FT/#DATA.ENitrousOx_EDGAR-FT.2008-2010_114reg0.EN2O_{sec_ehyde[s]}.csv"

    if os.path.isfile(path):
        TMP = load_data(path)
        for i in range(114 + 1):
            EN2Oeft[308: 310 + 1, regionJ_index[i], 0, kindGHG_index["N2O"], regionI_index[i]] += TMP[: 310 - 308 + 1, i]

# ===========
# 1.5. ACCMIP
# ===========

# load emissions from ACCMIP
# from [Lamarque et al., 2010]

# CH4
ECH4accmip = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=conf.dty)
p_ECH4_bio = np.zeros([nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=conf.dty)
for s in range(1, len(sec_accmip) - 2):
    path = f"data/EMethane_ACCMIP/#DATA.EMethane_ACCMIP.1850-2000_114reg0.ECH4_{sec_accmip[s]}.csv"

    if os.path.isfile(path):
        TMP = load_data(path)
        for i in range(114 + 1):
            ECH4accmip[150: 300 + 1, regionJ_index[i], 0, kindGHG_index["CH4"], regionI_index[i]] += TMP[: 300 - 150 + 1, i]
            if sec_accmip[s] in ["agr", "awb", "wst"]:
                p_ECH4_bio[regionJ_index[i], 0, kindGHG_index["CH4"], regionI_index[i]] += TMP[0, i]
p_ECH4_bio /= ECH4accmip[150]
p_ECH4_bio[np.isnan(p_ECH4_bio) | np.isinf(p_ECH4_bio)] = 0

# ===============
# 1.6. EDGAR-HYDE
# ===============

# load emissions from EDGAR-HYDE v1.3
# from [van Aardenne et al., 2001]

# N2O
EN2Oehyde = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=conf.dty)
p_EN2O_bio = np.zeros([nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=conf.dty)
for s in range(1, len(sec_ehyde) - 2):
    path = f"data/ENitrousOx_EDGAR-HYDE/#DATA.ENitrousOx_EDGAR-HYDE.1890-1990_114reg0.EN2O_{sec_ehyde[s]}.csv"

    if os.path.isfile(path):
        TMP = load_data(path)
        for i in range(114 + 1):
            EN2Oehyde[190: 290 + 1, regionJ_index[i], 0, kindGHG_index["N2O"], regionI_index[i]] += TMP[: 290 - 190 + 1, i]
            if sec_ehyde[s] in ["al", "an", "aw", "lf"]:
                p_EN2O_bio[regionJ_index[i], 0, kindGHG_index["N2O"], regionI_index[i]] += TMP[0, i]
p_EN2O_bio /= EN2Oehyde[190]
p_EN2O_bio[np.isnan(p_EN2O_bio) | np.isinf(p_EN2O_bio)] = 0

# ==============
# 1.7. Stern1998
# ==============

# load emissions from [Stern et al., 1998]
ECH4stern = np.zeros([ind_cdiac + 1], dtype=conf.dty)
path = "data/EMethane_Stern1998/#DATA.EMethane_Stern1998.1860-1994_(7sec).ECH4.csv"
TMP, lgd = load_data_and_header(path)
for s in range(len(lgd)):
    if not lgd[s] in ["Biomass Burning"]:
        ECH4stern[160: 294 + 1] += TMP[: 294 - 160 + 1, s]

# =================
# 1.8. Davidson2009
# =================

# load emissions from [Davidson et al., 2009]
EN2Odavidson = np.zeros([ind_cdiac + 1], dtype=conf.dty)
path = "data/ENitrousOx_Davidson2009/#DATA.ENitrousOx_Davidson2009.1860-2005_(5sec).EN2O.csv"
TMP, lgd = load_data_and_header(path)
for s in range(len(lgd)):
    if lgd[s] in ["nyl_prod", "ff_burn", "biogen"]:
        EN2Odavidson[160: 305 + 1] += TMP[: 305 - 160 + 1, s]

# global rescaling of EDGAR-HYDE emissions
EN2OehydeR = EN2Oehyde.copy()
if True:
    EN2OehydeR[190: 290 + 1, ...] *= (EN2Odavidson[190: 290 + 1] / np.sum(np.sum(np.sum(np.sum(EN2Oehyde[190: 290 + 1, ...], 4), 3), 2), 1))[:, np.newaxis, np.newaxis, np.newaxis, np.newaxis]

# =========
# 1.9. SRES
# =========

# initialization of projected drivers
EFFproj = np.zeros([conf.ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=conf.dty)
ECH4proj = np.zeros([conf.ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=conf.dty)
EN2Oproj = np.zeros([conf.ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=conf.dty)

# projection of emissions under SRES scenarios
# from [IPCC, 2000]

# FF
if (conf.scen_EFF[:4] == "SRES") & (conf.ind_final > ind_cdiac):
    path = f"data/EFossil_SRES/#DATA.EFossil_SRES.2000-2100_4reg0.{conf.scen_EFF[5:]}_EFF.csv"
    TMP = load_data(path)
    for i in range(4 + 1):
        if (conf.mod_regionI == "SRES4") & (conf.mod_regionJ == "SRES4"):
            EFFproj[300: min(conf.ind_final, 400) + 1, i, 0, kFF, i] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]
        elif (conf.mod_regionI == "SRES4") & (conf.mod_regionJ != "SRES4"):
            EFFproj[300: min(conf.ind_final, 400) + 1, 0, 0, kFF, i] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]
        elif (conf.mod_regionI != "SRES4") & (conf.mod_regionJ == "SRES4"):
            EFFproj[300: min(conf.ind_final, 400) + 1, i, 0, kFF, 0] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]
        elif (conf.mod_regionI != "SRES4") & (conf.mod_regionJ != "SRES4"):
            EFFproj[300: min(conf.ind_final, 400) + 1, 0, 0, kFF, 0] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]

# CH4
if (conf.scen_ECH4[:4] == "SRES") & (conf.ind_final > ind_cdiac):
    path = f"data/EMethane_SRES/#DATA.EMethane_SRES.2000-2100_4reg0.{conf.scen_ECH4[5:]}_ECH4.csv"
    TMP = load_data(path)
    for i in range(4 + 1):
        if (conf.mod_regionI == "SRES4") & (conf.mod_regionJ == "SRES4"):
            ECH4proj[300: min(conf.ind_final, 400) + 1, i, 0, kindGHG_index["CH4"], i] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]
        elif (conf.mod_regionI == "SRES4") & (conf.mod_regionJ != "SRES4"):
            ECH4proj[300: min(conf.ind_final, 400) + 1, 0, 0, kindGHG_index["CH4"], i] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]
        elif (conf.mod_regionI != "SRES4") & (conf.mod_regionJ == "SRES4"):
            ECH4proj[300: min(conf.ind_final, 400) + 1, i, 0, kindGHG_index["CH4"], 0] += TMP[
                                                                                     : min(conf.ind_final, 400) - 300 + 1, i]
        elif (conf.mod_regionI != "SRES4") & (conf.mod_regionJ != "SRES4"):
            ECH4proj[300: min(conf.ind_final, 400) + 1, 0, 0, kindGHG_index["CH4"], 0] += TMP[
                                                                                     : min(conf.ind_final, 400) - 300 + 1, i]

# N2O
if (conf.scen_EN2O[:4] == "SRES") & (conf.ind_final > ind_cdiac):
    path = f"data/ENitrousOx_SRES/#DATA.ENitrousOx_SRES.2000-2100_4reg0.{conf.scen_EN2O[5:]}_EN2O.csv"
    TMP = load_data(path)
    for i in range(4 + 1):
        if (conf.mod_regionI == "SRES4") & (conf.mod_regionJ == "SRES4"):
            EN2Oproj[300: min(conf.ind_final, 400) + 1, i, 0, kindGHG_index["N2O"], i] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]
        elif (conf.mod_regionI == "SRES4") & (conf.mod_regionJ != "SRES4"):
            EN2Oproj[300: min(conf.ind_final, 400) + 1, 0, 0, kindGHG_index["N2O"], i] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]
        elif (conf.mod_regionI != "SRES4") & (conf.mod_regionJ == "SRES4"):
            EN2Oproj[300: min(conf.ind_final, 400) + 1, i, 0, kindGHG_index["N2O"], 0] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]
        elif (conf.mod_regionI != "SRES4") & (conf.mod_regionJ != "SRES4"):
            EN2Oproj[300: min(conf.ind_final, 400) + 1, 0, 0, kindGHG_index["N2O"], 0] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]

# =========
# 1.10. RCP
# =========

# projection of emissions under RCP scenarios
# from [Meinshausen et al., 2011]

# FF
if (conf.scen_EFF[:3] == "RCP") & (conf.ind_final > ind_cdiac):
    path = f"data/EFossil_RCP/#DATA.EFossil_RCP.2000-2100_5reg0.rcp{conf.scen_EFF[3]}{conf.scen_EFF[5]}_EFF.csv"
    TMP = load_data(path)
    for i in range(5 + 1):
        if (conf.mod_regionI == "RCP5") & (conf.mod_regionJ == "RCP5"):
            EFFproj[300: min(conf.ind_final, 400) + 1, i, 0, kFF, i] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]
        elif (conf.mod_regionI == "RCP5") & (conf.mod_regionJ != "RCP5"):
            EFFproj[300: min(conf.ind_final, 400) + 1, 0, 0, kFF, i] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]
        elif (conf.mod_regionI != "RCP5") & (conf.mod_regionJ == "RCP5"):
            EFFproj[300: min(conf.ind_final, 400) + 1, i, 0, kFF, 0] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]
        elif (conf.mod_regionI != "RCP5") & (conf.mod_regionJ != "RCP5"):
            EFFproj[300: min(conf.ind_final, 400) + 1, 0, 0, kFF, 0] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]

# CH4
if (conf.scen_ECH4[:3] == "RCP") & (conf.ind_final > ind_cdiac):
    for s in range(1, len(sec_accmip) - 2):
        path = f"data/EMethane_RCP/#DATA.EMethane_RCP.2000-2100_5reg0.rcp{conf.scen_ECH4[3]}{conf.scen_ECH4[5]}_ECH4_" f"{sec_accmip[s]}.csv"

        if os.path.isfile(path):
            TMP = load_data(path)
            for i in range(5 + 1):
                if (conf.mod_regionI == "RCP5") & (conf.mod_regionJ == "RCP5"):
                    ECH4proj[300: min(conf.ind_final, 400) + 1, i, 0, kindGHG_index["CH4"], i] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]
                elif (conf.mod_regionI == "RCP5") & (conf.mod_regionJ != "RCP5"):
                    ECH4proj[300: min(conf.ind_final, 400) + 1, 0, 0, kindGHG_index["CH4"], i] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]
                elif (conf.mod_regionI != "RCP5") & (conf.mod_regionJ == "RCP5"):
                    ECH4proj[300: min(conf.ind_final, 400) + 1, i, 0, kindGHG_index["CH4"], 0] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]
                elif (conf.mod_regionI != "RCP5") & (conf.mod_regionJ != "RCP5"):
                    ECH4proj[300: min(conf.ind_final, 400) + 1, 0, 0, kindGHG_index["CH4"], 0] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]

# N2O
if (conf.scen_EN2O[:3] == "RCP") & (conf.ind_final > ind_cdiac):
    path = f"data/ENitrousOx_RCP/#DATA.ENitrousOx_RCP.2000-2100_5reg0.rcp{conf.scen_EN2O[3]}{conf.scen_EN2O[5]}_EN2O.csv"
    TMP = load_data(path)
    for i in range(5 + 1):
        if (conf.mod_regionI == "RCP5") & (conf.mod_regionJ == "RCP5"):
            EN2Oproj[300: min(conf.ind_final, 400) + 1, i, 0, kindGHG_index["N2O"], i] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]
        elif (conf.mod_regionI == "RCP5") & (conf.mod_regionJ != "RCP5"):
            EN2Oproj[300: min(conf.ind_final, 400) + 1, 0, 0, kindGHG_index["N2O"], i] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]
        elif (conf.mod_regionI != "RCP5") & (conf.mod_regionJ == "RCP5"):
            EN2Oproj[300: min(conf.ind_final, 400) + 1, i, 0, kindGHG_index["N2O"], 0] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]
        elif (conf.mod_regionI != "RCP5") & (conf.mod_regionJ != "RCP5"):
            EN2Oproj[300: min(conf.ind_final, 400) + 1, 0, 0, kindGHG_index["N2O"], 0] += TMP[: min(conf.ind_final, 400) - 300 + 1, i]

# =================
# 1.A. PAST DATASET
# =================

# datasets mixed following trends


# with CDIAC as reference
if conf.data_EFF == "CDIAC":
    EFFpast = EFFcdiac.copy()

# with EDGAR as reference
elif conf.data_EFF == "EDGAR":
    EFFpast = EFFedgar.copy()

    # follow EDGAR-FT variations after 2008
    for t in range(ind_edgar + 1, ind_cdiac + 1):
        EFFpast[t, ...] = EFFpast[t - 1, ...] * EFFeft[t, ...] / EFFeft[t - 1, ...]
        EFFpast[np.isnan(EFFpast) | np.isinf(EFFpast)] = 0
        EFFpast[t, ...] *= np.sum(EFFpast[t - 1, ...]) / np.sum(EFFpast[t, ...]) * np.sum(EFFeft[t, ...]) / np.sum(EFFeft[t - 1, ...])

    EFFpast[np.isnan(EFFpast) | np.isinf(EFFpast)] = 0

    # follow CDIAC variations before 1970
    for t in range(50, 270)[::-1]:
        EFFpast[t, ...] = EFFpast[t + 1, ...] * EFFcdiac[t, ...] / EFFcdiac[t + 1, ...]
        EFFpast[np.isnan(EFFpast) | np.isinf(EFFpast)] = 0
        EFFpast[t, ...] *= np.sum(EFFpast[t + 1, ...]) / np.sum(EFFpast[t, ...]) * np.sum(EFFcdiac[t, ...]) / np.sum(EFFcdiac[t + 1, ...])
    EFFpast[np.isnan(EFFpast) | np.isinf(EFFpast)] = 0

# with EDGAR as reference
if conf.data_ECH4 == "EDGAR":
    ECH4past = ECH4edgar.copy()
    # follow EDGAR-FT variations after 2008
    for t in range(ind_edgar + 1, ind_cdiac + 1):
        ECH4past[t, ...] = ECH4past[t - 1, ...] * ECH4eft[t, ...] / ECH4eft[t - 1, ...]
        ECH4past[np.isnan(ECH4past) | np.isinf(ECH4past)] = 0
        ECH4past[t, ...] *= np.sum(ECH4past[t - 1, ...]) / np.sum(ECH4past[t, ...]) * np.sum(ECH4eft[t, ...]) / np.sum(ECH4eft[t - 1, ...])
    ECH4past[np.isnan(ECH4past) | np.isinf(ECH4past)] = 0

    # follow ACCMIP variations before 1970
    for t in range(150, 270)[::-1]:
        ECH4past[t, ...] = ECH4past[t + 1, ...] * ECH4accmip[t, ...] / ECH4accmip[
            t + 1, ...]
        ECH4past[np.isnan(ECH4past) | np.isinf(ECH4past)] = 0
        ECH4past[t, ...] *= np.sum(ECH4past[t + 1, ...]) / np.sum(ECH4past[t, ...]) * np.sum(ECH4accmip[t, ...]) / np.sum(ECH4accmip[t + 1, ...])
    ECH4past[np.isnan(ECH4past) | np.isinf(ECH4past)] = 0
    # linear extrapolation before 1850
    for t in range(50, 150):
        ECH4past[t, ...] = (1 - p_ECH4_bio) * ECH4past[150, ...] * (t - 50) / float(150 - 50)
    for t in range(0, 150):
        ECH4past[t, ...] += p_ECH4_bio * ECH4past[150, ...] * (t - -200) / float(150 - -200)
    ECH4_0 = ECH4past[50, ...]

# with ACCMIP as reference
elif conf.data_ECH4 == "ACCMIP":
    ECH4past = ECH4accmip.copy()
    # follow EDGAR variations after 2000
    for t in range(300 + 1, ind_edgar + 1):
        ECH4past[t, ...] = ECH4past[t - 1, ...] * ECH4edgar[t, ...] / ECH4edgar[
            t - 1, ...]
        ECH4past[np.isnan(ECH4past) | np.isinf(ECH4past)] = 0
        ECH4past[t, ...] *= np.sum(ECH4past[t - 1, ...]) / np.sum(ECH4past[t, ...]) * np.sum(ECH4edgar[t, ...]) / np.sum(ECH4edgar[t - 1, ...])
    ECH4past[np.isnan(ECH4past) | np.isinf(ECH4past)] = 0

    # follow EDGAR-FT variations after 2008
    for t in range(ind_edgar + 1, ind_cdiac + 1):
        ECH4past[t, ...] = ECH4past[t - 1, ...] * ECH4eft[t, ...] / ECH4eft[t - 1, ...]
        ECH4past[np.isnan(ECH4past) | np.isinf(ECH4past)] = 0
        ECH4past[t, ...] *= np.sum(ECH4past[t - 1, ...]) / np.sum(ECH4past[t, ...]) * np.sum(ECH4eft[t, ...]) / np.sum(ECH4eft[t - 1, ...])
    ECH4past[np.isnan(ECH4past) | np.isinf(ECH4past)] = 0

    # linear extrapolation before 1850
    for t in range(50, 150):
        ECH4past[t, ...] = (1 - p_ECH4_bio) * ECH4past[150, ...] * (t - 50) / float(150 - 50)
    for t in range(0, 150):
        ECH4past[t, ...] += p_ECH4_bio * ECH4past[150, ...] * (t - -200) / float(150 + 200)
    ECH4_0 = ECH4past[50, ...]

# with EPA as reference
elif conf.data_ECH4 == "EPA":
    ECH4past = ECH4epa.copy()
    # follow ACCMIP variations before 1990
    for t in range(150, 290)[::-1]:
        ECH4past[t, ...] = ECH4past[t + 1, ...] * ECH4accmip[t, ...] / ECH4accmip[
            t + 1, ...]
        ECH4past[np.isnan(ECH4past) | np.isinf(ECH4past)] = 0
        ECH4past[t, ...] *= np.sum(ECH4past[t + 1, ...]) / np.sum(ECH4past[t, ...]) * np.sum(ECH4accmip[t, ...]) / np.sum(ECH4accmip[t + 1, ...])
    ECH4past[np.isnan(ECH4past) | np.isinf(ECH4past)] = 0

    # linear extrapolation before 1850
    for t in range(50, 150):
        ECH4past[t, ...] = (1 - p_ECH4_bio) * ECH4past[150, ...] * (t - 50) / float(150 - 50)
    for t in range(0, 150):
        ECH4past[t, ...] += p_ECH4_bio * ECH4past[150, ...] * (t - -200) / float(150 - -200)
    ECH4_0 = ECH4past[50, ...]

# with EDGAR as reference
if conf.data_EN2O == "EDGAR":
    EN2Opast = EN2Oedgar.copy()
    # follow EDGAR-FT variations after 2008
    for t in range(ind_edgar + 1, ind_cdiac + 1):
        EN2Opast[t, ...] = EN2Opast[t - 1, ...] * EN2Oeft[t, ...] / EN2Oeft[t - 1, ...]
        EN2Opast[np.isnan(EN2Opast) | np.isinf(EN2Opast)] = 0
        EN2Opast[t, ...] *= np.sum(EN2Opast[t - 1, ...]) / np.sum(EN2Opast[t, ...]) * np.sum(EN2Oeft[t, ...]) / np.sum(EN2Oeft[t - 1, ...])
    EN2Opast[np.isnan(EN2Opast) | np.isinf(EN2Opast)] = 0

    # follow EDGAR-HYDE variations before 1970
    for t in range(190, 270)[::-1]:
        EN2Opast[t, ...] = EN2Opast[t + 1, ...] * EN2OehydeR[t, ...] / EN2OehydeR[t + 1, ...]
        EN2Opast[np.isnan(EN2Opast) | np.isinf(EN2Opast)] = 0
        EN2Opast[t, ...] *= np.sum(EN2Opast[t + 1, ...]) / np.sum(EN2Opast[t, ...]) * np.sum(EN2OehydeR[t, ...]) / np.sum(EN2OehydeR[t + 1, ...])
    EN2Opast[np.isnan(EN2Opast) | np.isinf(EN2Opast)] = 0

    # linear extrapolation before 1890
    for t in range(50, 190):
        EN2Opast[t, ...] = (1 - p_EN2O_bio) * EN2Opast[190, ...] * (t - 50) / float(190 - 50)
    for t in range(0, 190):
        EN2Opast[t, ...] += p_EN2O_bio * EN2Opast[190, ...] * (t - -200) / float(190 - -200)
    EN2O_0 = EN2Opast[50, ...]

# with EPA as reference
elif conf.data_EN2O == "EPA":
    EN2Opast = EN2Oepa.copy()
    # follow EDGAR-HYDE variations before 1990
    for t in range(190, 290)[::-1]:
        EN2Opast[t, ...] = EN2Opast[t + 1, ...] * EN2OehydeR[t, ...] / EN2OehydeR[
            t + 1, ...]
        EN2Opast[np.isnan(EN2Opast) | np.isinf(EN2Opast)] = 0
        EN2Opast[t, ...] *= np.sum(EN2Opast[t + 1, ...]) / np.sum(EN2Opast[t, ...]) * np.sum(EN2OehydeR[t, ...]) / np.sum(EN2OehydeR[t + 1, ...])
    EN2Opast[np.isnan(EN2Opast) | np.isinf(EN2Opast)] = 0
    # linear extrapolation before 1850
    for t in range(50, 190):
        EN2Opast[t, ...] = (1 - p_EN2O_bio) * EN2Opast[190, ...] * (t - 50) / float(190 - 50)
    for t in range(0, 190):
        EN2Opast[t, ...] += p_EN2O_bio * EN2Opast[190, ...] * (t - -200) / float(190 - -200)
    EN2O_0 = EN2Opast[50, ...]

# Cut past dataset to right length
for arr, past in [(EFF, EFFpast), (ECH4, ECH4past), (EN2O, EN2Opast)]:
    arr[:min(ind_cdiac, conf.ind_final) + 1, ...] = past[:min(ind_cdiac, conf.ind_final) + 1, ...]

# ==================
# 1.B. FINAL DATASET
# ==================

# datasets mixed following various criteria
for VAR, scen, arr, arr_0, arrproj in [
    ("FF", conf.scen_EFF, EFF, None, EFFproj), ("CH4", conf.scen_ECH4, ECH4, ECH4_0, ECH4proj), ("N2O", conf.scen_EN2O, EN2O, EN2O_0, EN2Oproj),
]:
    # stop emissions
    if (scen == "stop") & (conf.ind_final > ind_cdiac):
        if VAR in ["FF"]:
            EFF[ind_cdiac + 1:, ...] = 0
        else:
            arr[ind_cdiac + 1:, ...] = arr_0[np.newaxis, ...]

    # constant emissions
    elif (scen == "cst") & (conf.ind_final > ind_cdiac):
        arr[ind_cdiac + 1:, ...] = arr[ind_cdiac, ...][np.newaxis, ...]

    # RCP or SRES scenarios
    elif ((scen[:4] == "SRES") | (scen[:3] == "RCP")) & (conf.ind_final > ind_cdiac):

        # raw discontinuity
        if conf.mod_DATAscen == "raw":
            arr[ind_cdiac + 1:, ...] = arrproj[ind_cdiac + 1:, ...]

        # offset at transition point
        elif conf.mod_DATAscen == "offset":
            arr[ind_cdiac + 1:, ...] = arrproj[ind_cdiac + 1:, ...] - arrproj[
                ind_cdiac, ...] + arr[ind_cdiac, ...]
            for t in range(ind_cdiac + 1, conf.ind_final + 1):
                def_regI = bool(np.sum(arrproj[t, :, ..., 1:]))
                def_regJ = bool(np.sum(arrproj[t, 1:, ..., :]))
                if not def_regI:
                    arr[t, :, ..., 0] += np.sum(arr[t, :, ..., 1:], -1)
                    arr[t, :, ..., 1:] = 0
                if not def_regJ:
                    arr[t, 0, ..., :] += np.sum(arr[t, 1:, ..., :], 0)
                    arr[t, 1:, ..., :] = 0

                    # linear transition over N years
        elif conf.mod_DATAscen[:6] == "smooth":
            N = int(conf.mod_DATAscen[6:])
            if conf.ind_final >= ind_cdiac + N:
                for t in range(ind_cdiac + 1, ind_cdiac + N):
                    arr[t, ...] = (1 - (t - ind_cdiac) / float(N)) * arr[
                        ind_cdiac, ...] + (t - ind_cdiac) / float(N) * arrproj[
                                      ind_cdiac + N, ...]
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
        elif conf.mod_DATAscen == "trends":
            for t in range(ind_cdiac + 1, conf.ind_final + 1):
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

# delete individual datasets
del EFFcdiac, EFFedgar, EFFeft, EFFpast, EFFproj
del ECH4epa, ECH4edgar, ECH4eft, ECH4accmip, ECH4past, ECH4proj, ECH4stern, p_ECH4_bio
del EN2Oepa, EN2Oedgar, EN2Oeft, EN2Oehyde, EN2Opast, EN2Oproj, EN2Odavidson, p_EN2O_bio

# ===========
# 1.11. PETERS
# ===========

# TODO
