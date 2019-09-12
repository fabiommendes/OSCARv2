import csv
import os

import numpy as np

from .a1_regions import nb_regionJ, nb_kind, nb_regionI, \
    nb_sector, regionJ_index, kFF, regionI_index, ind_final, kindGHG_index
from ...config import dty, scen_EFF, mod_regionI, scen_EN2O, scen_ECH4, data_ECH4, \
    data_EN2O, mod_DATAscen, mod_regionJ, data_EFF

# _exec = exec
# _glob = globals()
#
#
# def exec(x: int):
#     print(x)
#     _exec(x, _glob, _glob)


##################################################
#   1. GREENHOUSE GASES
##################################################

EFF = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
               dtype=dty)  # {GtC/yr}
ECH4 = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                dtype=dty)  # {TgC/yr}
EN2O = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                dtype=dty)  # {TgN/yr}

# ==========
# 1.1. CDIAC
# ==========

# load CDIAC emissions
# from [Boden et al., 2012]
ind_cdiac = 310
kin = ["sol", "liq", "gas", "cem", "fla"]

# total emissions
EFFcdiac = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                    dtype=dty)
for k in range(len(kin)):
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open("data/EFossil_CDIAC/#DATA.EFossil_CDIAC.1751-2010_114reg0.EFF_" +
                 kin[k] + ".csv", "r")
        )],
        dtype=dty,
    )
    for i in range(114 + 1):
        EFFcdiac[51: ind_cdiac + 1, regionJ_index[i], 0, kFF, regionI_index[i]] += \
            TMP[: ind_cdiac - 51 + 1, i]

# distribution among fuels
# TODO

# ========
# 1.2. EPA
# ========

# load EPA emissions
# see [EPA, 2012]
ind_epa = 310
sec_epa = [
    "agr_ferm",
    "agr_manu",
    "agr_othr",
    "agr_rice",
    "agr_soil",
    "ene_burn",
    "ene_coal",
    "ene_comb",
    "ene_ngos",
    "ene_othr",
    "ind_acid",
]
sec_epa1 = [
    "agr_ferm",
    "agr_manu",
    "agr_othr",
    "agr_rice",
    "agr_soil",
    "ene_coal",
    "ene_comb",
    "ene_ngos",
    "ene_othr",
    "ind_acid",
]

# CH4
ECH4epa = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
for s in range(len(sec_epa1)):
    if os.path.isfile(
            "data/EMethane_EPA/#DATA.EMethane_EPA.1990-" + str(
                1700 + ind_epa) + "_114reg0.ECH4_" + sec_epa1[s] + ".csv"
    ):
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/EMethane_EPA/#DATA.EMethane_EPA.1990-"
                    + str(1700 + ind_epa)
                    + "_114reg0.ECH4_"
                    + sec_epa1[s]
                    + ".csv",
                    "r",
                )
            )],
            dtype=dty,
        )
        for i in range(114 + 1):
            ECH4epa[290: ind_epa + 1, regionJ_index[i], 0, kindGHG_index["CH4"],
            regionI_index[i]] += TMP[: ind_epa - 290 + 1, i]

# N2O
EN2Oepa = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
for s in range(len(sec_epa1)):
    if os.path.isfile(
            "data/ENitrousOx_EPA/#DATA.ENitrousOx_EPA.1990-" + str(
                1700 + ind_epa) + "_114reg0.EN2O_" + sec_epa1[s] + ".csv"
    ):
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/ENitrousOx_EPA/#DATA.ENitrousOx_EPA.1990-"
                    + str(1700 + ind_epa)
                    + "_114reg0.EN2O_"
                    + sec_epa1[s]
                    + ".csv",
                    "r",
                )
            )],
            dtype=dty,
        )
        for i in range(114 + 1):
            EN2Oepa[290: ind_epa + 1, regionJ_index[i], 0, kindGHG_index["N2O"],
            regionI_index[i]] += TMP[: ind_epa - 290 + 1, i]

# ==========
# 1.3. EDGAR
# ==========

# load emissions from EDGAR v4.2
# see [JRC, 2011]
ind_edgar = 308
sec_ehyde = ["oo", "fc", "fp", "bc", "in", "al", "an", "aw", "lf", "sb", "df"]
sec_accmip = ["ooo", "ene", "ind", "tra", "shp", "air", "dom", "slv", "agr", "awb", "wst",
              "for", "gra"]

# FF
EFFedgar = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                    dtype=dty)
for s in range(1, len(sec_ehyde) - 2):
    if os.path.isfile(
            "data/ECarbon_EDGAR/#DATA.ECarbon_EDGAR.1970-"
            + str(1700 + ind_edgar)
            + "_114reg0.ECO2_"
            + sec_ehyde[s]
            + ".csv"
    ):
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/ECarbon_EDGAR/#DATA.ECarbon_EDGAR.1970-"
                    + str(1700 + ind_edgar)
                    + "_114reg0.ECO2_"
                    + sec_ehyde[s]
                    + ".csv",
                    "r",
                )
            )],
            dtype=dty,
        )
        for i in range(114 + 1):
            EFFedgar[270: ind_edgar + 1, regionJ_index[i], 0, kFF,
            regionI_index[i]] += TMP[: ind_edgar - 270 + 1, i]

# CH4
ECH4edgar = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                     dtype=dty)
for s in range(1, len(sec_accmip) - 2):
    if os.path.isfile(
            "data/EMethane_EDGAR/#DATA.EMethane_EDGAR.1970-"
            + str(1700 + ind_edgar)
            + "_114reg0.ECH4_"
            + sec_accmip[s]
            + ".csv"
    ):
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/EMethane_EDGAR/#DATA.EMethane_EDGAR.1970-"
                    + str(1700 + ind_edgar)
                    + "_114reg0.ECH4_"
                    + sec_accmip[s]
                    + ".csv",
                    "r",
                )
            )],
            dtype=dty,
        )
        for i in range(114 + 1):
            ECH4edgar[270: ind_edgar + 1, regionJ_index[i], 0, kindGHG_index["CH4"],
            regionI_index[i]] += TMP[: ind_edgar - 270 + 1, i]

# N2O
EN2Oedgar = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                     dtype=dty)
for s in range(1, len(sec_ehyde) - 2):
    if os.path.isfile(
            "data/ENitrousOx_EDGAR/#DATA.ENitrousOx_EDGAR.1970-"
            + str(1700 + ind_edgar)
            + "_114reg0.EN2O_"
            + sec_ehyde[s]
            + ".csv"
    ):
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/ENitrousOx_EDGAR/#DATA.ENitrousOx_EDGAR.1970-"
                    + str(1700 + ind_edgar)
                    + "_114reg0.EN2O_"
                    + sec_ehyde[s]
                    + ".csv",
                    "r",
                )
            )],
            dtype=dty,
        )
        for i in range(114 + 1):
            EN2Oedgar[270: ind_edgar + 1, regionJ_index[i], 0, kindGHG_index["N2O"],
            regionI_index[i]] += TMP[: ind_edgar - 270 + 1, i]

# =============
# 1.4. EDGAR-FT
# =============

# load emissions from EDGAR-FT v4.2-FT2010
# see [JRC, 2013]

# FF
EFFeft = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
for s in range(1, len(sec_ehyde) - 2):
    if os.path.isfile(
            "data/ECarbon_EDGAR-FT/#DATA.ECarbon_EDGAR-FT.2008-2010_114reg0.ECO2_" +
            sec_ehyde[s] + ".csv"):
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/ECarbon_EDGAR-FT/#DATA.ECarbon_EDGAR-FT.2008-2010_114reg0.ECO2_" +
                    sec_ehyde[s] + ".csv",
                    "r",
                )
            )],
            dtype=dty,
        )
        for i in range(114 + 1):
            EFFeft[308: 310 + 1, regionJ_index[i], 0, kFF, regionI_index[i]] += \
                TMP[: 310 - 308 + 1, i]

# CH4
ECH4eft = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
for s in range(1, len(sec_accmip) - 2):
    if os.path.isfile(
            "data/EMethane_EDGAR-FT/#DATA.EMethane_EDGAR-FT.2008-2010_114reg0.ECH4_" +
            sec_accmip[s] + ".csv"
    ):
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/EMethane_EDGAR-FT/#DATA.EMethane_EDGAR-FT.2008-2010_114reg0.ECH4_"
                    + sec_accmip[s]
                    + ".csv",
                    "r",
                )
            )],
            dtype=dty,
        )
        for i in range(114 + 1):
            ECH4eft[308: 310 + 1, regionJ_index[i], 0, kindGHG_index["CH4"],
            regionI_index[i]] += TMP[: 310 - 308 + 1, i]

# N2O
EN2Oeft = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
for s in range(1, len(sec_ehyde) - 2):
    if os.path.isfile(
            "data/ENitrousOx_EDGAR-FT/#DATA.ENitrousOx_EDGAR-FT.2008-2010_114reg0.EN2O_" +
            sec_ehyde[s] + ".csv"
    ):
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/ENitrousOx_EDGAR-FT/#DATA.ENitrousOx_EDGAR-FT.2008-2010_114reg0.EN2O_"
                    + sec_ehyde[s]
                    + ".csv",
                    "r",
                )
            )],
            dtype=dty,
        )
        for i in range(114 + 1):
            EN2Oeft[308: 310 + 1, regionJ_index[i], 0, kindGHG_index["N2O"],
            regionI_index[i]] += TMP[: 310 - 308 + 1, i]

# ===========
# 1.5. ACCMIP
# ===========

# load emissions from ACCMIP
# from [Lamarque et al., 2010]

# CH4
ECH4accmip = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                      dtype=dty)
p_ECH4_bio = np.zeros([nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
for s in range(1, len(sec_accmip) - 2):
    if os.path.isfile(
            "data/EMethane_ACCMIP/#DATA.EMethane_ACCMIP.1850-2000_114reg0.ECH4_" +
            sec_accmip[s] + ".csv"):
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/EMethane_ACCMIP/#DATA.EMethane_ACCMIP.1850-2000_114reg0.ECH4_" +
                    sec_accmip[s] + ".csv",
                    "r",
                )
            )],
            dtype=dty,
        )
        for i in range(114 + 1):
            ECH4accmip[150: 300 + 1, regionJ_index[i], 0, kindGHG_index["CH4"],
            regionI_index[i]] += TMP[: 300 - 150 + 1, i]
            if sec_accmip[s] in ["agr", "awb", "wst"]:
                p_ECH4_bio[regionJ_index[i], 0, kindGHG_index["CH4"], regionI_index[i]] += \
                    TMP[0, i]
p_ECH4_bio /= ECH4accmip[150]
p_ECH4_bio[np.isnan(p_ECH4_bio) | np.isinf(p_ECH4_bio)] = 0

# ===============
# 1.6. EDGAR-HYDE
# ===============

# load emissions from EDGAR-HYDE v1.3
# from [van Aardenne et al., 2001]

# N2O
EN2Oehyde = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                     dtype=dty)
p_EN2O_bio = np.zeros([nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
for s in range(1, len(sec_ehyde) - 2):
    if os.path.isfile(
            "data/ENitrousOx_EDGAR-HYDE/#DATA.ENitrousOx_EDGAR-HYDE.1890-1990_114reg0.EN2O_" +
            sec_ehyde[s] + ".csv"
    ):
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/ENitrousOx_EDGAR-HYDE/#DATA.ENitrousOx_EDGAR-HYDE.1890-1990_114reg0.EN2O_"
                    + sec_ehyde[s]
                    + ".csv",
                    "r",
                )
            )],
            dtype=dty,
        )
        for i in range(114 + 1):
            EN2Oehyde[190: 290 + 1, regionJ_index[i], 0, kindGHG_index["N2O"],
            regionI_index[i]] += TMP[: 290 - 190 + 1, i]
            if sec_ehyde[s] in ["al", "an", "aw", "lf"]:
                p_EN2O_bio[regionJ_index[i], 0, kindGHG_index["N2O"], regionI_index[i]] += \
                    TMP[0, i]
p_EN2O_bio /= EN2Oehyde[190]
p_EN2O_bio[np.isnan(p_EN2O_bio) | np.isinf(p_EN2O_bio)] = 0

# ==============
# 1.7. Stern1998
# ==============

# load emissions from [Stern et al., 1998]
ECH4stern = np.zeros([ind_cdiac + 1], dtype=dty)
TMP = np.array(
    [
        line
        for line in csv.reader(
        open("data/EMethane_Stern1998/#DATA.EMethane_Stern1998.1860-1994_(7sec).ECH4.csv",
             "r"))][1:],
    dtype=dty,
)
lgd = [
    line for line in csv.reader(
        open("data/EMethane_Stern1998/#DATA.EMethane_Stern1998.1860-1994_(7sec).ECH4.csv",
             "r"))
][0]
for s in range(len(lgd)):
    if not lgd[s] in ["Biomass Burning"]:
        ECH4stern[160: 294 + 1] += TMP[: 294 - 160 + 1, s]

# =================
# 1.8. Davidson2009
# =================

# load emissions from [Davidson et al., 2009]
EN2Odavidson = np.zeros([ind_cdiac + 1], dtype=dty)
TMP = np.array(
    [
        line
        for line in csv.reader(
        open(
            "data/ENitrousOx_Davidson2009/#DATA.ENitrousOx_Davidson2009.1860-2005_(5sec).EN2O.csv",
            "r")
    )][1:],
    dtype=dty,
)
lgd = [
    line
    for line in csv.reader(
        open(
            "data/ENitrousOx_Davidson2009/#DATA.ENitrousOx_Davidson2009.1860-2005_(5sec).EN2O.csv",
            "r")
    )
][0]
for s in range(len(lgd)):
    if lgd[s] in ["nyl_prod", "ff_burn", "biogen"]:
        EN2Odavidson[160: 305 + 1] += TMP[: 305 - 160 + 1, s]

# global rescaling of EDGAR-HYDE emissions
EN2OehydeR = EN2Oehyde.copy()
if True:
    EN2OehydeR[190: 290 + 1, ...] *= (
                                             EN2Odavidson[190: 290 + 1] / np.sum(np.sum(
                                         np.sum(np.sum(EN2Oehyde[190: 290 + 1, ...], 4),
                                                3), 2), 1)
                                     )[:, np.newaxis, np.newaxis, np.newaxis, np.newaxis]

# =========
# 1.9. SRES
# =========

# initialization of projected drivers
EFFproj = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
ECH4proj = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                    dtype=dty)
EN2Oproj = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                    dtype=dty)

# projection of emissions under SRES scenarios
# from [IPCC, 2000]

# FF
if (scen_EFF[:4] == "SRES") & (ind_final > ind_cdiac):
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open("data/EFossil_SRES/#DATA.EFossil_SRES.2000-2100_4reg0." + scen_EFF[
                                                                           5:] + "_EFF.csv",
                 "r")
        )],
        dtype=dty,
    )
    for i in range(4 + 1):
        if (mod_regionI == "SRES4") & (mod_regionJ == "SRES4"):
            EFFproj[300: min(ind_final, 400) + 1, i, 0, kFF, i] += TMP[: min(ind_final,
                                                                             400) - 300 + 1,
                                                                   i]
        elif (mod_regionI == "SRES4") & (mod_regionJ != "SRES4"):
            EFFproj[300: min(ind_final, 400) + 1, 0, 0, kFF, i] += TMP[: min(ind_final,
                                                                             400) - 300 + 1,
                                                                   i]
        elif (mod_regionI != "SRES4") & (mod_regionJ == "SRES4"):
            EFFproj[300: min(ind_final, 400) + 1, i, 0, kFF, 0] += TMP[: min(ind_final,
                                                                             400) - 300 + 1,
                                                                   i]
        elif (mod_regionI != "SRES4") & (mod_regionJ != "SRES4"):
            EFFproj[300: min(ind_final, 400) + 1, 0, 0, kFF, 0] += TMP[: min(ind_final,
                                                                             400) - 300 + 1,
                                                                   i]

# CH4
if (scen_ECH4[:4] == "SRES") & (ind_final > ind_cdiac):
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open("data/EMethane_SRES/#DATA.EMethane_SRES.2000-2100_4reg0." + scen_ECH4[
                                                                             5:] + "_ECH4.csv",
                 "r")
        )],
        dtype=dty,
    )
    for i in range(4 + 1):
        if (mod_regionI == "SRES4") & (mod_regionJ == "SRES4"):
            ECH4proj[300: min(ind_final, 400) + 1, i, 0, kindGHG_index["CH4"], i] += TMP[
                                                                                     : min(
                                                                                         ind_final,
                                                                                         400) - 300 + 1,
                                                                                     i]
        elif (mod_regionI == "SRES4") & (mod_regionJ != "SRES4"):
            ECH4proj[300: min(ind_final, 400) + 1, 0, 0, kindGHG_index["CH4"], i] += TMP[
                                                                                     : min(
                                                                                         ind_final,
                                                                                         400) - 300 + 1,
                                                                                     i]
        elif (mod_regionI != "SRES4") & (mod_regionJ == "SRES4"):
            ECH4proj[300: min(ind_final, 400) + 1, i, 0, kindGHG_index["CH4"], 0] += TMP[
                                                                                     : min(
                                                                                         ind_final,
                                                                                         400) - 300 + 1,
                                                                                     i]
        elif (mod_regionI != "SRES4") & (mod_regionJ != "SRES4"):
            ECH4proj[300: min(ind_final, 400) + 1, 0, 0, kindGHG_index["CH4"], 0] += TMP[
                                                                                     : min(
                                                                                         ind_final,
                                                                                         400) - 300 + 1,
                                                                                     i]

# N2O
if (scen_EN2O[:4] == "SRES") & (ind_final > ind_cdiac):
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/ENitrousOx_SRES/#DATA.ENitrousOx_SRES.2000-2100_4reg0." + scen_EN2O[
                                                                                5:] + "_EN2O.csv",
                "r")
        )],
        dtype=dty,
    )
    for i in range(4 + 1):
        if (mod_regionI == "SRES4") & (mod_regionJ == "SRES4"):
            EN2Oproj[300: min(ind_final, 400) + 1, i, 0, kindGHG_index["N2O"], i] += TMP[
                                                                                     : min(
                                                                                         ind_final,
                                                                                         400) - 300 + 1,
                                                                                     i]
        elif (mod_regionI == "SRES4") & (mod_regionJ != "SRES4"):
            EN2Oproj[300: min(ind_final, 400) + 1, 0, 0, kindGHG_index["N2O"], i] += TMP[
                                                                                     : min(
                                                                                         ind_final,
                                                                                         400) - 300 + 1,
                                                                                     i]
        elif (mod_regionI != "SRES4") & (mod_regionJ == "SRES4"):
            EN2Oproj[300: min(ind_final, 400) + 1, i, 0, kindGHG_index["N2O"], 0] += TMP[
                                                                                     : min(
                                                                                         ind_final,
                                                                                         400) - 300 + 1,
                                                                                     i]
        elif (mod_regionI != "SRES4") & (mod_regionJ != "SRES4"):
            EN2Oproj[300: min(ind_final, 400) + 1, 0, 0, kindGHG_index["N2O"], 0] += TMP[
                                                                                     : min(
                                                                                         ind_final,
                                                                                         400) - 300 + 1,
                                                                                     i]

# =========
# 1.10. RCP
# =========

# projection of emissions under RCP scenarios
# from [Meinshausen et al., 2011]

# FF
if (scen_EFF[:3] == "RCP") & (ind_final > ind_cdiac):
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/EFossil_RCP/#DATA.EFossil_RCP.2000-2100_5reg0.rcp" + scen_EFF[3] +
                scen_EFF[5] + "_EFF.csv",
                "r",
            )
        )],
        dtype=dty,
    )
    for i in range(5 + 1):
        if (mod_regionI == "RCP5") & (mod_regionJ == "RCP5"):
            EFFproj[300: min(ind_final, 400) + 1, i, 0, kFF, i] += \
                TMP[: min(ind_final, 400) - 300 + 1, i]
        elif (mod_regionI == "RCP5") & (mod_regionJ != "RCP5"):
            EFFproj[300: min(ind_final, 400) + 1, 0, 0, kFF, i] += \
                TMP[: min(ind_final, 400) - 300 + 1, i]
        elif (mod_regionI != "RCP5") & (mod_regionJ == "RCP5"):
            EFFproj[300: min(ind_final, 400) + 1, i, 0, kFF, 0] += \
                TMP[: min(ind_final, 400) - 300 + 1, i]
        elif (mod_regionI != "RCP5") & (mod_regionJ != "RCP5"):
            EFFproj[300: min(ind_final, 400) + 1, 0, 0, kFF, 0] += \
                TMP[: min(ind_final, 400) - 300 + 1, i]

# CH4
if (scen_ECH4[:3] == "RCP") & (ind_final > ind_cdiac):
    for s in range(1, len(sec_accmip) - 2):
        if os.path.isfile(
                "data/EMethane_RCP/#DATA.EMethane_RCP.2000-2100_5reg0.rcp"
                + scen_ECH4[3]
                + scen_ECH4[5]
                + "_ECH4_"
                + sec_accmip[s]
                + ".csv"
        ):
            TMP = np.array(
                [
                    line
                    for line in csv.reader(
                    open(
                        "data/EMethane_RCP/#DATA.EMethane_RCP.2000-2100_5reg0.rcp"
                        + scen_ECH4[3]
                        + scen_ECH4[5]
                        + "_ECH4_"
                        + sec_accmip[s]
                        + ".csv",
                        "r",
                    )
                )],
                dtype=dty,
            )
            for i in range(5 + 1):
                if (mod_regionI == "RCP5") & (mod_regionJ == "RCP5"):
                    ECH4proj[300: min(ind_final, 400) + 1, i, 0, kindGHG_index["CH4"],
                    i] += TMP[: min(ind_final, 400) - 300 + 1, i]
                elif (mod_regionI == "RCP5") & (mod_regionJ != "RCP5"):
                    ECH4proj[300: min(ind_final, 400) + 1, 0, 0, kindGHG_index["CH4"],
                    i] += TMP[: min(ind_final, 400) - 300 + 1, i]
                elif (mod_regionI != "RCP5") & (mod_regionJ == "RCP5"):
                    ECH4proj[300: min(ind_final, 400) + 1, i, 0, kindGHG_index["CH4"],
                    0] += TMP[: min(ind_final, 400) - 300 + 1, i]
                elif (mod_regionI != "RCP5") & (mod_regionJ != "RCP5"):
                    ECH4proj[300: min(ind_final, 400) + 1, 0, 0, kindGHG_index["CH4"],
                    0] += TMP[: min(ind_final, 400) - 300 + 1, i]

# N2O
if (scen_EN2O[:3] == "RCP") & (ind_final > ind_cdiac):
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/ENitrousOx_RCP/#DATA.ENitrousOx_RCP.2000-2100_5reg0.rcp"
                + scen_EN2O[3]
                + scen_EN2O[5]
                + "_EN2O.csv",
                "r",
            )
        )],
        dtype=dty,
    )
    for i in range(5 + 1):
        if (mod_regionI == "RCP5") & (mod_regionJ == "RCP5"):
            EN2Oproj[300: min(ind_final, 400) + 1, i, 0, kindGHG_index["N2O"], i] += TMP[
                                                                                     : min(
                                                                                         ind_final,
                                                                                         400) - 300 + 1,
                                                                                     i]
        elif (mod_regionI == "RCP5") & (mod_regionJ != "RCP5"):
            EN2Oproj[300: min(ind_final, 400) + 1, 0, 0, kindGHG_index["N2O"], i] += TMP[
                                                                                     : min(
                                                                                         ind_final,
                                                                                         400) - 300 + 1,
                                                                                     i]
        elif (mod_regionI != "RCP5") & (mod_regionJ == "RCP5"):
            EN2Oproj[300: min(ind_final, 400) + 1, i, 0, kindGHG_index["N2O"], 0] += TMP[
                                                                                     : min(
                                                                                         ind_final,
                                                                                         400) - 300 + 1,
                                                                                     i]
        elif (mod_regionI != "RCP5") & (mod_regionJ != "RCP5"):
            EN2Oproj[300: min(ind_final, 400) + 1, 0, 0, kindGHG_index["N2O"], 0] += TMP[
                                                                                     : min(
                                                                                         ind_final,
                                                                                         400) - 300 + 1,
                                                                                     i]

# =================
# 1.A. PAST DATASET
# =================

# datasets mixed following trends
for VAR in ["FF", "CH4", "N2O"]:

    if VAR in ["FF"]:

        # with CDIAC as reference
        if data_EFF == "CDIAC":
            exec("E" + VAR + "past = E" + VAR + "cdiac.copy()")

        # with EDGAR as reference
        elif data_EFF == "EDGAR":
            exec("E" + VAR + "past = E" + VAR + "edgar.copy()")
            # follow EDGAR-FT variations after 2008
            for t in range(ind_edgar + 1, ind_cdiac + 1):
                exec(
                    "E"
                    + VAR
                    + "past[t,...] = E"
                    + VAR
                    + "past[t-1,...] * E"
                    + VAR
                    + "eft[t,...]/E"
                    + VAR
                    + "eft[t-1,...]"
                )
                exec(
                    "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
                exec(
                    "E"
                    + VAR
                    + "past[t,...] *= np.sum(E"
                    + VAR
                    + "past[t-1,...])/np.sum(E"
                    + VAR
                    + "past[t,...]) * np.sum(E"
                    + VAR
                    + "eft[t,...])/np.sum(E"
                    + VAR
                    + "eft[t-1,...])"
                )
            exec(
                "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
            # follow CDIAC variations before 1970
            for t in range(50, 270)[::-1]:
                exec(
                    "E"
                    + VAR
                    + "past[t,...] = E"
                    + VAR
                    + "past[t+1,...] * E"
                    + VAR
                    + "cdiac[t,...]/E"
                    + VAR
                    + "cdiac[t+1,...]"
                )
                exec(
                    "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
                exec(
                    "E"
                    + VAR
                    + "past[t,...] *= np.sum(E"
                    + VAR
                    + "past[t+1,...])/np.sum(E"
                    + VAR
                    + "past[t,...]) * np.sum(E"
                    + VAR
                    + "cdiac[t,...])/np.sum(E"
                    + VAR
                    + "cdiac[t+1,...])"
                )
            exec(
                "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")

    elif VAR in ["CH4"]:

        # with EDGAR as reference
        if data_ECH4 == "EDGAR":
            exec("E" + VAR + "past = E" + VAR + "edgar.copy()")
            # follow EDGAR-FT variations after 2008
            for t in range(ind_edgar + 1, ind_cdiac + 1):
                exec(
                    "E"
                    + VAR
                    + "past[t,...] = E"
                    + VAR
                    + "past[t-1,...] * E"
                    + VAR
                    + "eft[t,...]/E"
                    + VAR
                    + "eft[t-1,...]"
                )
                exec(
                    "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
                exec(
                    "E"
                    + VAR
                    + "past[t,...] *= np.sum(E"
                    + VAR
                    + "past[t-1,...])/np.sum(E"
                    + VAR
                    + "past[t,...]) * np.sum(E"
                    + VAR
                    + "eft[t,...])/np.sum(E"
                    + VAR
                    + "eft[t-1,...])"
                )
            exec(
                "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
            # follow ACCMIP variations before 1970
            for t in range(150, 270)[::-1]:
                exec(
                    "E"
                    + VAR
                    + "past[t,...] = E"
                    + VAR
                    + "past[t+1,...] * E"
                    + VAR
                    + "accmip[t,...]/E"
                    + VAR
                    + "accmip[t+1,...]"
                )
                exec(
                    "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
                exec(
                    "E"
                    + VAR
                    + "past[t,...] *= np.sum(E"
                    + VAR
                    + "past[t+1,...])/np.sum(E"
                    + VAR
                    + "past[t,...]) * np.sum(E"
                    + VAR
                    + "accmip[t,...])/np.sum(E"
                    + VAR
                    + "accmip[t+1,...])"
                )
            exec(
                "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
            # linear extrapolation before 1850
            for t in range(50, 150):
                exec(
                    "E"
                    + VAR
                    + "past[t,...] = (1-p_E"
                    + VAR
                    + "_bio) * E"
                    + VAR
                    + "past[150,...] * (t-50)/float(150-50)"
                )
            for t in range(0, 150):
                exec(
                    "E"
                    + VAR
                    + "past[t,...] += p_E"
                    + VAR
                    + "_bio * E"
                    + VAR
                    + "past[150,...] * (t--200)/float(150--200)"
                )
            exec("E" + VAR + "_0 = E" + VAR + "past[50,...]")

        # with ACCMIP as reference
        elif data_ECH4 == "ACCMIP":
            exec("E" + VAR + "past = E" + VAR + "accmip.copy()")
            # follow EDGAR variations after 2000
            for t in range(300 + 1, ind_edgar + 1):
                exec(
                    "E"
                    + VAR
                    + "past[t,...] = E"
                    + VAR
                    + "past[t-1,...] * E"
                    + VAR
                    + "edgar[t,...]/E"
                    + VAR
                    + "edgar[t-1,...]"
                )
                exec(
                    "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
                exec(
                    "E"
                    + VAR
                    + "past[t,...] *= np.sum(E"
                    + VAR
                    + "past[t-1,...])/np.sum(E"
                    + VAR
                    + "past[t,...]) * np.sum(E"
                    + VAR
                    + "edgar[t,...])/np.sum(E"
                    + VAR
                    + "edgar[t-1,...])"
                )
            exec(
                "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
            # follow EDGAR-FT variations after 2008
            for t in range(ind_edgar + 1, ind_cdiac + 1):
                exec(
                    "E"
                    + VAR
                    + "past[t,...] = E"
                    + VAR
                    + "past[t-1,...] * E"
                    + VAR
                    + "eft[t,...]/E"
                    + VAR
                    + "eft[t-1,...]"
                )
                exec(
                    "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
                exec(
                    "E"
                    + VAR
                    + "past[t,...] *= np.sum(E"
                    + VAR
                    + "past[t-1,...])/np.sum(E"
                    + VAR
                    + "past[t,...]) * np.sum(E"
                    + VAR
                    + "eft[t,...])/np.sum(E"
                    + VAR
                    + "eft[t-1,...])"
                )
            exec(
                "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
            # linear extrapolation before 1850
            for t in range(50, 150):
                exec(
                    "E"
                    + VAR
                    + "past[t,...] = (1-p_E"
                    + VAR
                    + "_bio) * E"
                    + VAR
                    + "past[150,...] * (t-50)/float(150-50)"
                )
            for t in range(0, 150):
                exec(
                    "E"
                    + VAR
                    + "past[t,...] += p_E"
                    + VAR
                    + "_bio * E"
                    + VAR
                    + "past[150,...] * (t--200)/float(150--200)"
                )
            exec("E" + VAR + "_0 = E" + VAR + "past[50,...]")

        # with EPA as reference
        elif data_ECH4 == "EPA":
            exec("E" + VAR + "past = E" + VAR + "epa.copy()")
            # follow ACCMIP variations before 1990
            for t in range(150, 290)[::-1]:
                exec(
                    "E"
                    + VAR
                    + "past[t,...] = E"
                    + VAR
                    + "past[t+1,...] * E"
                    + VAR
                    + "accmip[t,...]/E"
                    + VAR
                    + "accmip[t+1,...]"
                )
                exec(
                    "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
                exec(
                    "E"
                    + VAR
                    + "past[t,...] *= np.sum(E"
                    + VAR
                    + "past[t+1,...])/np.sum(E"
                    + VAR
                    + "past[t,...]) * np.sum(E"
                    + VAR
                    + "accmip[t,...])/np.sum(E"
                    + VAR
                    + "accmip[t+1,...])"
                )
            exec(
                "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
            # linear extrapolation before 1850
            for t in range(50, 150):
                exec(
                    "E"
                    + VAR
                    + "past[t,...] = (1-p_E"
                    + VAR
                    + "_bio) * E"
                    + VAR
                    + "past[150,...] * (t-50)/float(150-50)"
                )
            for t in range(0, 150):
                exec(
                    "E"
                    + VAR
                    + "past[t,...] += p_E"
                    + VAR
                    + "_bio * E"
                    + VAR
                    + "past[150,...] * (t--200)/float(150--200)"
                )
            exec("E" + VAR + "_0 = E" + VAR + "past[50,...]")

    elif VAR in ["N2O"]:

        # with EDGAR as reference
        if data_EN2O == "EDGAR":
            exec("E" + VAR + "past = E" + VAR + "edgar.copy()")
            # follow EDGAR-FT variations after 2008
            for t in range(ind_edgar + 1, ind_cdiac + 1):
                exec(
                    "E"
                    + VAR
                    + "past[t,...] = E"
                    + VAR
                    + "past[t-1,...] * E"
                    + VAR
                    + "eft[t,...]/E"
                    + VAR
                    + "eft[t-1,...]"
                )
                exec(
                    "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
                exec(
                    "E"
                    + VAR
                    + "past[t,...] *= np.sum(E"
                    + VAR
                    + "past[t-1,...])/np.sum(E"
                    + VAR
                    + "past[t,...]) * np.sum(E"
                    + VAR
                    + "eft[t,...])/np.sum(E"
                    + VAR
                    + "eft[t-1,...])"
                )
            exec(
                "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
            # follow EDGAR-HYDE variations before 1970
            for t in range(190, 270)[::-1]:
                exec(
                    "E"
                    + VAR
                    + "past[t,...] = E"
                    + VAR
                    + "past[t+1,...] * E"
                    + VAR
                    + "ehydeR[t,...]/E"
                    + VAR
                    + "ehydeR[t+1,...]"
                )
                exec(
                    "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
                exec(
                    "E"
                    + VAR
                    + "past[t,...] *= np.sum(E"
                    + VAR
                    + "past[t+1,...])/np.sum(E"
                    + VAR
                    + "past[t,...]) * np.sum(E"
                    + VAR
                    + "ehydeR[t,...])/np.sum(E"
                    + VAR
                    + "ehydeR[t+1,...])"
                )
            exec(
                "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
            # linear extrapolation before 1890
            for t in range(50, 190):
                exec(
                    "E"
                    + VAR
                    + "past[t,...] = (1-p_E"
                    + VAR
                    + "_bio) * E"
                    + VAR
                    + "past[190,...] * (t-50)/float(190-50)"
                )
            for t in range(0, 190):
                exec(
                    "E"
                    + VAR
                    + "past[t,...] += p_E"
                    + VAR
                    + "_bio * E"
                    + VAR
                    + "past[190,...] * (t--200)/float(190--200)"
                )
            exec("E" + VAR + "_0 = E" + VAR + "past[50,...]")

        # with EPA as reference
        elif data_EN2O == "EPA":
            exec("E" + VAR + "past = E" + VAR + "epa.copy()")
            # follow EDGAR-HYDE variations before 1990
            for t in range(190, 290)[::-1]:
                exec(
                    "E"
                    + VAR
                    + "past[t,...] = E"
                    + VAR
                    + "past[t+1,...] * E"
                    + VAR
                    + "ehydeR[t,...]/E"
                    + VAR
                    + "ehydeR[t+1,...]"
                )
                exec(
                    "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
                exec(
                    "E"
                    + VAR
                    + "past[t,...] *= np.sum(E"
                    + VAR
                    + "past[t+1,...])/np.sum(E"
                    + VAR
                    + "past[t,...]) * np.sum(E"
                    + VAR
                    + "ehydeR[t,...])/np.sum(E"
                    + VAR
                    + "ehydeR[t+1,...])"
                )
            exec(
                "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
            # linear extrapolation before 1850
            for t in range(50, 190):
                exec(
                    "E"
                    + VAR
                    + "past[t,...] = (1-p_E"
                    + VAR
                    + "_bio) * E"
                    + VAR
                    + "past[190,...] * (t-50)/float(190-50)"
                )
            for t in range(0, 190):
                exec(
                    "E"
                    + VAR
                    + "past[t,...] += p_E"
                    + VAR
                    + "_bio * E"
                    + VAR
                    + "past[190,...] * (t--200)/float(190--200)"
                )
            exec("E" + VAR + "_0 = E" + VAR + "past[50,...]")

    # cut past dataset to right length
    exec(
        "E" + VAR + "[:min(ind_cdiac,ind_final)+1,...] = E" + VAR + "past[:min(ind_cdiac,ind_final)+1,...]")

# ==================
# 1.B. FINAL DATASET
# ==================

# datasets mixed following various criteria
for VAR in ["FF", "CH4", "N2O"]:
    exec("scen = scen_E" + VAR)

    # stop emissions
    if (scen == "stop") & (ind_final > ind_cdiac):
        if VAR in ["FF"]:
            exec("E" + VAR + "[ind_cdiac+1:,...] = 0")
        else:
            exec("E" + VAR + "[ind_cdiac+1:,...] = E" + VAR + "_0[np.newaxis,...]")

    # constant emissions
    elif (scen == "cst") & (ind_final > ind_cdiac):
        exec(
            "E" + VAR + "[ind_cdiac+1:,...] = E" + VAR + "[ind_cdiac,...][np.newaxis,...]")

        # RCP or SRES scenarios
    elif ((scen[:4] == "SRES") | (scen[:3] == "RCP")) & (ind_final > ind_cdiac):

        # raw discontinuity
        if mod_DATAscen == "raw":
            exec("E" + VAR + "[ind_cdiac+1:,...] = E" + VAR + "proj[ind_cdiac+1:,...]")

        # offset at transition point
        elif mod_DATAscen == "offset":
            exec(
                "E"
                + VAR
                + "[ind_cdiac+1:,...] = E"
                + VAR
                + "proj[ind_cdiac+1:,...] - E"
                + VAR
                + "proj[ind_cdiac,...] + E"
                + VAR
                + "[ind_cdiac,...]"
            )
            for t in range(ind_cdiac + 1, ind_final + 1):
                def_regI = eval("bool(np.sum(E" + VAR + "proj[t,:,...,1:]))")
                def_regJ = eval("bool(np.sum(E" + VAR + "proj[t,1:,...,:]))")
                if not def_regI:
                    exec("E" + VAR + "[t,:,...,0] += np.sum(E" + VAR + "[t,:,...,1:],-1)")
                    exec("E" + VAR + "[t,:,...,1:] = 0")
                if not def_regJ:
                    exec("E" + VAR + "[t,0,...,:] += np.sum(E" + VAR + "[t,1:,...,:],0)")
                    exec("E" + VAR + "[t,1:,...,:] = 0")

                    # linear transition over N years
        elif mod_DATAscen[:6] == "smooth":
            N = int(mod_DATAscen[6:])
            if ind_final >= ind_cdiac + N:
                for t in range(ind_cdiac + 1, ind_cdiac + N):
                    exec(
                        "E"
                        + VAR
                        + "[t,...] = (1-(t-ind_cdiac)/float(N)) * E"
                        + VAR
                        + "[ind_cdiac,...] + (t-ind_cdiac)/float(N) * E"
                        + VAR
                        + "proj[ind_cdiac+N,...]"
                    )
                    def_regI = exec("bool(np.sum(E" + VAR + "proj[t,:,...,1:]))")
                    def_regJ = exec("bool(np.sum(E" + VAR + "proj[t,1:,...,:]))")
                    if not def_regI:
                        exec(
                            "E" + VAR + "[t,:,...,0] += np.sum(E" + VAR + "[t,:,...,1:],-1)")
                        exec("E" + VAR + "[t,:,...,1:] = 0")
                    if not def_regJ:
                        exec(
                            "E" + VAR + "[t,0,...,:] += np.sum(E" + VAR + "[t,1:,...,:],0)")
                        exec("E" + VAR + "[t,1:,...,:] = 0")
                exec(
                    "E" + VAR + "[ind_cdiac+N:,...] = E" + VAR + "proj[ind_cdiac+N:,...]")

        # follow trends of projection
        elif mod_DATAscen == "trends":
            for t in range(ind_cdiac + 1, ind_final + 1):
                exec("def_regI = bool(np.sum(E" + VAR + "proj[t,:,...,1:]))")
                exec("def_regJ = bool(np.sum(E" + VAR + "proj[t,1:,...,:]))")
                if def_regI and def_regJ:
                    exec(
                        "E"
                        + VAR
                        + "[t,...] = E"
                        + VAR
                        + "[t-1,...] * E"
                        + VAR
                        + "proj[t,...]/E"
                        + VAR
                        + "proj[t-1,...]"
                    )
                    exec(
                        "E" + VAR + "[np.isnan(E" + VAR + ")|np.isinf(E" + VAR + ")] = 0")
                    exec(
                        "E"
                        + VAR
                        + "[t,...] *= np.sum(E"
                        + VAR
                        + "[t-1,...])/np.sum(E"
                        + VAR
                        + "[t,...]) * np.sum(E"
                        + VAR
                        + "proj[t,...])/np.sum(E"
                        + VAR
                        + "proj[t-1,...])"
                    )
                elif not def_regI and def_regJ:
                    exec(
                        "E"
                        + VAR
                        + "[t,:,...,0] = np.sum(E"
                        + VAR
                        + "[t-1,:,...,:],-1) * np.sum(E"
                        + VAR
                        + "proj[t,:,...,:],-1)/np.sum(E"
                        + VAR
                        + "proj[t-1,:,...,:],-1)"
                    )
                    exec(
                        "E" + VAR + "[np.isnan(E" + VAR + ")|np.isinf(E" + VAR + ")] = 0")
                    exec(
                        "E"
                        + VAR
                        + "[t,...] *= np.sum(E"
                        + VAR
                        + "[t-1,...])/np.sum(E"
                        + VAR
                        + "[t,...]) * np.sum(E"
                        + VAR
                        + "proj[t,...])/np.sum(E"
                        + VAR
                        + "proj[t-1,...])"
                    )
                elif def_regI and not def_regJ:
                    exec(
                        "E"
                        + VAR
                        + "[t,0,...,:] = np.sum(E"
                        + VAR
                        + "[t-1,:,...,:],0) * np.sum(E"
                        + VAR
                        + "proj[t,:,...,:],0)/np.sum(E"
                        + VAR
                        + "proj[t-1,:,...,:],0)"
                    )
                    exec(
                        "E" + VAR + "[np.isnan(E" + VAR + ")|np.isinf(E" + VAR + ")] = 0")
                    exec(
                        "E"
                        + VAR
                        + "[t,...] *= np.sum(E"
                        + VAR
                        + "[t-1,...])/np.sum(E"
                        + VAR
                        + "[t,...]) * np.sum(E"
                        + VAR
                        + "proj[t,...])/np.sum(E"
                        + VAR
                        + "proj[t-1,...])"
                    )
                elif not def_regI and not def_regJ:
                    exec(
                        "E"
                        + VAR
                        + "[t,0,...,0] = np.sum(np.sum(E"
                        + VAR
                        + "[t-1,:,...,:],-1),0) * np.sum(np.sum(E"
                        + VAR
                        + "proj[t,:,...,:],-1),0)/np.sum(np.sum(E"
                        + VAR
                        + "proj[t-1,:,...,:],-1),0)"
                    )
            exec("E" + VAR + "[np.isnan(E" + VAR + ")|np.isinf(E" + VAR + ")] = 0")

# delete individual datasets
for VAR in ["FF"]:
    exec(
        "del E" + VAR + "cdiac,E" + VAR + "edgar,E" + VAR + "eft,E" + VAR + "past,E" + VAR + "proj")
for VAR in ["CH4"]:
    exec(
        "del E"
        + VAR
        + "epa,E"
        + VAR
        + "edgar,E"
        + VAR
        + "eft,E"
        + VAR
        + "accmip,E"
        + VAR
        + "past,E"
        + VAR
        + "proj,E"
        + VAR
        + "stern,p_E"
        + VAR
        + "_bio"
    )
for VAR in ["N2O"]:
    exec(
        "del E"
        + VAR
        + "epa,E"
        + VAR
        + "edgar,E"
        + VAR
        + "eft,E"
        + VAR
        + "ehyde,E"
        + VAR
        + "past,E"
        + VAR
        + "proj,E"
        + VAR
        + "davidson,p_E"
        + VAR
        + "_bio"
    )

# ===========
# 1.11. PETERS
# ===========

# TODO
