import csv
import os

import numpy as np

from .a1_regions import nb_regionJ, nb_kind, nb_regionI, nb_sector, regionJ_index, \
    regionI_index, ind_final, kindAER_index, kindCHI_index
from .a2_greenhouse import sec_accmip, ind_edgar, ind_cdiac
from ...config import dty, mod_regionI, mod_regionJ, scen_ESO2, mod_DATAscen, scen_ECO, \
    scen_ENOX, scen_EVOC, scen_ENH3, scen_EBC, scen_EOC, data_ENOX, data_ECO, data_EVOC, \
    data_ESO2, data_ENH3, data_EBC, data_EOC

##################################################
#   4. SHORT-LIVED SPECIES
##################################################

ENOX = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                dtype=dty)  # {TgN/yr}
ECO = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
               dtype=dty)  # {TgC/yr}
EVOC = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                dtype=dty)  # {Tg/yr}
ESO2 = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                dtype=dty)  # {TgS/yr}
ENH3 = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                dtype=dty)  # {TgN/yr}
EOC = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
               dtype=dty)  # {Tg/yr}
EBC = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
               dtype=dty)  # {Tg/yr}

# ==========
# 4.1. EDGAR
# ==========

# load emissions from EDGAR v4.2
# see [JRC, 2011]

# OzoPrec
for VAR in ["NOX", "CO", "VOC"]:
    exec(
        "E" + VAR + "edgar = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)")
    for s in range(1, len(sec_accmip) - 2):
        if os.path.isfile(
                "data/EOzoPrec_EDGAR/#DATA.EOzoPrec_EDGAR.1970-"
                + str(1700 + ind_edgar)
                + "_114reg0.E"
                + VAR
                + "_"
                + sec_accmip[s]
                + ".csv"
        ):
            TMP = np.array(
                [
                    line
                    for line in csv.reader(
                    open(
                        "data/EOzoPrec_EDGAR/#DATA.EOzoPrec_EDGAR.1970-"
                        + str(1700 + ind_edgar)
                        + "_114reg0.E"
                        + VAR
                        + "_"
                        + sec_accmip[s]
                        + ".csv",
                        "r",
                    )
                )],
                dtype=dty,
            )
            for i in range(114 + 1):
                exec(
                    "E"
                    + VAR
                    + "edgar[270:ind_edgar+1,regionJ_index[i],0,kindCHI_index[VAR],regionI_index[i]] += TMP[:ind_edgar-270+1,i]"
                )

# AeroPrec
for VAR in ["SO2", "NH3"]:
    exec(
        "E" + VAR + "edgar = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)")
    for s in range(1, len(sec_accmip) - 2):
        if os.path.isfile(
                "data/EAeroPrec_EDGAR/#DATA.EAeroPrec_EDGAR.1970-"
                + str(1700 + ind_edgar)
                + "_114reg0.E"
                + VAR
                + "_"
                + sec_accmip[s]
                + ".csv"
        ):
            TMP = np.array(
                [
                    line
                    for line in csv.reader(
                    open(
                        "data/EAeroPrec_EDGAR/#DATA.EAeroPrec_EDGAR.1970-"
                        + str(1700 + ind_edgar)
                        + "_114reg0.E"
                        + VAR
                        + "_"
                        + sec_accmip[s]
                        + ".csv",
                        "r",
                    )
                )],
                dtype=dty,
            )
            for i in range(114 + 1):
                exec(
                    "E"
                    + VAR
                    + "edgar[270:ind_edgar+1,regionJ_index[i],0,kindAER_index[VAR],regionI_index[i]] += TMP[:ind_edgar-270+1,i]"
                )

# PM10 (proxy for OC/BC)
EPM10edgar = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                      dtype=dty)
for s in range(1, len(sec_accmip) - 2):
    if os.path.isfile(
            "data/EAeroPrec_EDGAR/#DATA.EAeroPrec_EDGAR.1970-"
            + str(1700 + ind_edgar)
            + "_114reg0.EPM10_"
            + sec_accmip[s]
            + ".csv"
    ):
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/EAeroPrec_EDGAR/#DATA.EAeroPrec_EDGAR.1970-"
                    + str(1700 + ind_edgar)
                    + "_114reg0.EPM10_"
                    + sec_accmip[s]
                    + ".csv",
                    "r",
                )
            )],
            dtype=dty,
        )
        for i in range(114 + 1):
            EPM10edgar[270: ind_edgar + 1, regionJ_index[i], 0, kindAER_index["OC"],
            regionI_index[i]] += (
                    TMP[: ind_edgar - 270 + 1, i] / 2
            )
            EPM10edgar[270: ind_edgar + 1, regionJ_index[i], 0, kindAER_index["BC"],
            regionI_index[i]] += (
                    TMP[: ind_edgar - 270 + 1, i] / 2
            )

# ===============
# 4.2. EDGAR-HTAP
# ===============

# load emissions from EDGAR-HTAP v2
# see [JRC, 2013]

# OzoPrec
for VAR in ["NOX", "CO", "VOC"]:
    exec(
        "E" + VAR + "ehtap = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)")
    for s in range(1, len(sec_accmip) - 2):
        if os.path.isfile(
                "data/EOzoPrec_EDGAR-HTAP/#DATA.EOzoPrec_EDGAR-HTAP.2008-2010_114reg0.E"
                + VAR
                + "_"
                + sec_accmip[s]
                + ".csv"
        ):
            TMP = np.array(
                [
                    line
                    for line in csv.reader(
                    open(
                        "data/EOzoPrec_EDGAR-HTAP/#DATA.EOzoPrec_EDGAR-HTAP.2008-2010_114reg0.E"
                        + VAR
                        + "_"
                        + sec_accmip[s]
                        + ".csv",
                        "r",
                    )
                )],
                dtype=dty,
            )
            for i in range(114 + 1):
                kindCHI_index
                exec(
                    "E"
                    + VAR
                    + "ehtap[308:310+1,regionJ_index[i],0,kindCHI_index[VAR],regionI_index[i]] += TMP[:ind_edgar-270+1,i]"
                )

# AeroPrec
for VAR in ["SO2", "NH3", "OC", "BC"]:
    exec(
        "E" + VAR + "ehtap = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)")
    for s in range(1, len(sec_accmip) - 2):
        if os.path.isfile(
                "data/EAeroPrec_EDGAR-HTAP/#DATA.EAeroPrec_EDGAR-HTAP.2008-2010_114reg0.E"
                + VAR
                + "_"
                + sec_accmip[s]
                + ".csv"
        ):
            TMP = np.array(
                [
                    line
                    for line in csv.reader(
                    open(
                        "data/EAeroPrec_EDGAR-HTAP/#DATA.EAeroPrec_EDGAR-HTAP.2008-2010_114reg0.E"
                        + VAR
                        + "_"
                        + sec_accmip[s]
                        + ".csv",
                        "r",
                    )
                )],
                dtype=dty,
            )
            for i in range(114 + 1):
                exec(
                    "E"
                    + VAR
                    + "ehtap[308:310+1,regionJ_index[i],0,kindAER_index[VAR],regionI_index[i]] += TMP[:ind_edgar-270+1,i]"
                )

# ===========
# 4.3. ACCMIP
# ===========

# load emissions from ACCMIP
# from [Lamarque et al., 2010]

# OzoPrec
for VAR in ["NOX", "CO", "VOC"]:
    exec(
        "E" + VAR + "accmip = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)")
    exec(
        "p_E" + VAR + "_bio = np.zeros([nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)")
    for s in range(1, len(sec_accmip) - 2):
        if os.path.isfile(
                "data/EOzoPrec_ACCMIP/#DATA.EOzoPrec_ACCMIP.1850-2000_114reg0.E" + VAR + "_" +
                sec_accmip[s] + ".csv"
        ):
            TMP = np.array(
                [
                    line
                    for line in csv.reader(
                    open(
                        "data/EOzoPrec_ACCMIP/#DATA.EOzoPrec_ACCMIP.1850-2000_114reg0.E"
                        + VAR
                        + "_"
                        + sec_accmip[s]
                        + ".csv",
                        "r",
                    )
                )],
                dtype=dty,
            )
            for i in range(114 + 1):
                exec(
                    "E"
                    + VAR
                    + "accmip[150:300+1,regionJ_index[i],0,kindCHI_index[VAR],regionI_index[i]] += TMP[:300-150+1,i]"
                )
                if sec_accmip[s] in ["agr", "awb", "wst"]:
                    exec(
                        "p_E" + VAR + "_bio[regionJ_index[i],0,kindCHI_index[VAR],regionI_index[i]] += TMP[0,i]")
    exec("p_E" + VAR + "_bio /= E" + VAR + "accmip[150]")
    exec(
        "p_E" + VAR + "_bio[np.isnan(p_E" + VAR + "_bio)|np.isinf(p_E" + VAR + "_bio)] = 0")

# AeroPrec
for VAR in ["SO2", "NH3", "OC", "BC"]:
    exec(
        "E" + VAR + "accmip = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)")
    exec(
        "p_E" + VAR + "_bio = np.zeros([nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)")
    for s in range(1, len(sec_accmip) - 2):
        if os.path.isfile(
                "data/EAeroPrec_ACCMIP/#DATA.EAeroPrec_ACCMIP.1850-2000_114reg0.E" + VAR + "_" +
                sec_accmip[s] + ".csv"
        ):
            TMP = np.array(
                [
                    line
                    for line in csv.reader(
                    open(
                        "data/EAeroPrec_ACCMIP/#DATA.EAeroPrec_ACCMIP.1850-2000_114reg0.E"
                        + VAR
                        + "_"
                        + sec_accmip[s]
                        + ".csv",
                        "r",
                    )
                )],
                dtype=dty,
            )
            for i in range(114 + 1):
                exec(
                    "E"
                    + VAR
                    + "accmip[150:300+1,regionJ_index[i],0,kindAER_index[VAR],regionI_index[i]] += TMP[:300-150+1,i]"
                )
                if sec_accmip[s] in ["agr", "awb", "wst"]:
                    exec(
                        "p_E" + VAR + "_bio[regionJ_index[i],0,kindAER_index[VAR],regionI_index[i]] += TMP[0,i]")
    exec("p_E" + VAR + "_bio /= E" + VAR + "accmip[150]")
    exec(
        "p_E" + VAR + "_bio[np.isnan(p_E" + VAR + "_bio)|np.isinf(p_E" + VAR + "_bio)] = 0")

# =========
# 4.4. SRES
# =========

# initialization of projected drivers
for VAR in ["NOX", "CO", "VOC"] + ["SO2", "NH3", "OC", "BC"]:
    exec(
        "E" + VAR + "proj = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)")

# projection of emissions under SRES scenarios
# from [IPCC, 2000]

# OzoPrec
scen_ENOX, scen_ECO, scen_EVOC
for VAR in ["NOX", "CO", "VOC"]:
    exec("scen = scen_E" + VAR)
    if (scen[:4] == "SRES") & (ind_final > ind_cdiac):
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/EOzoPrec_SRES/#DATA.EOzoPrec_SRES.2000-2100_4reg0." + scen[
                                                                                5:] + "_E" + VAR + ".csv",
                    "r"
                )
            )],
            dtype=dty,
        )
        for i in range(4 + 1):
            if (mod_regionI == "SRES4") & (mod_regionJ == "SRES4"):
                exec(
                    "E"
                    + VAR
                    + "proj[300:min(ind_final,400)+1,i,0,kindCHI_index[VAR],i] += TMP[:min(ind_final,400)-300+1,i]"
                )
            elif (mod_regionI == "SRES4") & (mod_regionJ != "SRES4"):
                exec(
                    "E"
                    + VAR
                    + "proj[300:min(ind_final,400)+1,0,0,kindCHI_index[VAR],i] += TMP[:min(ind_final,400)-300+1,i]"
                )
            elif (mod_regionI != "SRES4") & (mod_regionJ == "SRES4"):
                exec(
                    "E"
                    + VAR
                    + "proj[300:min(ind_final,400)+1,i,0,kindCHI_index[VAR],0] += TMP[:min(ind_final,400)-300+1,i]"
                )
            elif (mod_regionI != "SRES4") & (mod_regionJ != "SRES4"):
                exec(
                    "E"
                    + VAR
                    + "proj[300:min(ind_final,400)+1,0,0,kindCHI_index[VAR],0] += TMP[:min(ind_final,400)-300+1,i]"
                )

# AeroPrec (SO2 only)
if (scen_ESO2[:4] == "SRES") & (ind_final > ind_cdiac):
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open("data/EAeroPrec_SRES/#DATA.EAeroPrec_SRES.2000-2100_4reg0." + scen_ESO2[
                                                                               5:] + "_ESO2.csv",
                 "r")
        )],
        dtype=dty,
    )
    for i in range(4 + 1):
        if (mod_regionI == "SRES4") & (mod_regionJ == "SRES4"):
            ESO2proj[300: min(ind_final, 400) + 1, i, 0, kindAER_index["SO2"], i] += TMP[
                                                                                     : min(
                                                                                         ind_final,
                                                                                         400) - 300 + 1,
                                                                                     i]
        elif (mod_regionI == "SRES4") & (mod_regionJ != "SRES4"):
            ESO2proj[300: min(ind_final, 400) + 1, 0, 0, kindAER_index["SO2"], i] += TMP[
                                                                                     : min(
                                                                                         ind_final,
                                                                                         400) - 300 + 1,
                                                                                     i]
        elif (mod_regionI != "SRES4") & (mod_regionJ == "SRES4"):
            ESO2proj[300: min(ind_final, 400) + 1, i, 0, kindAER_index["SO2"], 0] += TMP[
                                                                                     : min(
                                                                                         ind_final,
                                                                                         400) - 300 + 1,
                                                                                     i]
        elif (mod_regionI != "SRES4") & (mod_regionJ != "SRES4"):
            ESO2proj[300: min(ind_final, 400) + 1, 0, 0, kindAER_index["SO2"], 0] += TMP[
                                                                                     : min(
                                                                                         ind_final,
                                                                                         400) - 300 + 1,
                                                                                     i]

# ========
# 4.5. RCP
# ========

# projection of emissions under RCP scenarios
# from [Meinshausen et al., 2011]

# OzoPrec
for VAR in ["NOX", "CO", "VOC"]:
    exec("scen = scen_E" + VAR)
    if (scen[:3] == "RCP") & (ind_final > ind_cdiac):
        for s in range(1, len(sec_accmip) - 2):
            if os.path.isfile(
                    "data/EOzoPrec_RCP/#DATA.EOzoPrec_RCP.2000-2100_5reg0.rcp"
                    + scen[3]
                    + scen[5]
                    + "_E"
                    + VAR
                    + "_"
                    + sec_accmip[s]
                    + ".csv"
            ):
                TMP = np.array(
                    [
                        line
                        for line in csv.reader(
                        open(
                            "data/EOzoPrec_RCP/#DATA.EOzoPrec_RCP.2000-2100_5reg0.rcp"
                            + scen[3]
                            + scen[5]
                            + "_E"
                            + VAR
                            + "_"
                            + sec_accmip[s]
                            + ".csv",
                            "r",
                        )
                    )],
                    dtype=dty,
                )
                for i in range(5 + 1):
                    if (mod_regionI == "RCP5") & (mod_regionJ == "RCP5"):
                        exec(
                            "E"
                            + VAR
                            + "proj[300:min(ind_final,400)+1,i,0,kindCHI_index[VAR],i] += TMP[:min(ind_final,400)-300+1,i]"
                        )
                    elif (mod_regionI == "RCP5") & (mod_regionJ != "RCP5"):
                        exec(
                            "E"
                            + VAR
                            + "proj[300:min(ind_final,400)+1,0,0,kindCHI_index[VAR],i] += TMP[:min(ind_final,400)-300+1,i]"
                        )
                    elif (mod_regionI != "RCP5") & (mod_regionJ == "RCP5"):
                        exec(
                            "E"
                            + VAR
                            + "proj[300:min(ind_final,400)+1,i,0,kindCHI_index[VAR],0] += TMP[:min(ind_final,400)-300+1,i]"
                        )
                    elif (mod_regionI != "RCP5") & (mod_regionJ != "RCP5"):
                        exec(
                            "E"
                            + VAR
                            + "proj[300:min(ind_final,400)+1,0,0,kindCHI_index[VAR],0] += TMP[:min(ind_final,400)-300+1,i]"
                        )

# AeroPrec
scen_ENH3
scen_EOC
scen_EBC

for VAR in ["SO2", "NH3", "OC", "BC"]:
    exec("scen = scen_E" + VAR)
    if (scen[:3] == "RCP") & (ind_final > ind_cdiac):
        for s in range(1, len(sec_accmip) - 2):
            if os.path.isfile(
                    "data/EAeroPrec_RCP/#DATA.EAeroPrec_RCP.2000-2100_5reg0.rcp"
                    + scen[3]
                    + scen[5]
                    + "_E"
                    + VAR
                    + "_"
                    + sec_accmip[s]
                    + ".csv"
            ):
                TMP = np.array(
                    [
                        line
                        for line in csv.reader(
                        open(
                            "data/EAeroPrec_RCP/#DATA.EAeroPrec_RCP.2000-2100_5reg0.rcp"
                            + scen[3]
                            + scen[5]
                            + "_E"
                            + VAR
                            + "_"
                            + sec_accmip[s]
                            + ".csv",
                            "r",
                        )
                    )],
                    dtype=dty,
                )
                for i in range(5 + 1):
                    if (mod_regionI == "RCP5") & (mod_regionJ == "RCP5"):
                        exec(
                            "E"
                            + VAR
                            + "proj[300:min(ind_final,400)+1,i,0,kindAER_index[VAR],i] += TMP[:min(ind_final,400)-300+1,i]"
                        )
                    elif (mod_regionI == "RCP5") & (mod_regionJ != "RCP5"):
                        exec(
                            "E"
                            + VAR
                            + "proj[300:min(ind_final,400)+1,0,0,kindAER_index[VAR],i] += TMP[:min(ind_final,400)-300+1,i]"
                        )
                    elif (mod_regionI != "RCP5") & (mod_regionJ == "RCP5"):
                        exec(
                            "E"
                            + VAR
                            + "proj[300:min(ind_final,400)+1,i,0,kindAER_index[VAR],0] += TMP[:min(ind_final,400)-300+1,i]"
                        )
                    elif (mod_regionI != "RCP5") & (mod_regionJ != "RCP5"):
                        exec(
                            "E"
                            + VAR
                            + "proj[300:min(ind_final,400)+1,0,0,kindAER_index[VAR],0] += TMP[:min(ind_final,400)-300+1,i]"
                        )

# =================
# 4.A. PAST DATASET
# =================

# datasets mixed following trends
series = {
    "NOX": data_ENOX,
    "CO": data_ECO,
    "VOC": data_EVOC,
    "SO2": data_ESO2,
    "NH3": data_ENH3,
    "BC": data_EBC,
    "OC": data_EOC,
}
for VAR in ["NOX", "CO", "VOC"] + ["SO2", "NH3", "OC", "BC"]:

    if VAR in ["NOX", "CO", "VOC"] + ["SO2", "NH3"]:
        exec("data = data_E" + VAR)

        # with EDGAR as reference
        if data == "EDGAR":
            exec("E" + VAR + "past = E" + VAR + "edgar.copy()")
            # follow EDGAR-HTAP variations after 2008
            for t in range(ind_edgar + 1, ind_cdiac + 1):
                exec(
                    "E"
                    + VAR
                    + "past[t,...] = E"
                    + VAR
                    + "past[t-1,...] * E"
                    + VAR
                    + "ehtap[t,...]/E"
                    + VAR
                    + "ehtap[t-1,...]"
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
                    + "ehtap[t,...])/np.sum(E"
                    + VAR
                    + "ehtap[t-1,...])"
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
        elif data == "ACCMIP":
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
            # follow EDGAR-HTAP variations after 2008
            for t in range(ind_edgar + 1, ind_cdiac + 1):
                exec(
                    "E"
                    + VAR
                    + "past[t,...] = E"
                    + VAR
                    + "past[t-1,...] * E"
                    + VAR
                    + "ehtap[t,...]/E"
                    + VAR
                    + "ehtap[t-1,...]"
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
                    + "ehtap[t,...])/np.sum(E"
                    + VAR
                    + "ehtap[t-1,...])"
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

    elif VAR in ["OC", "BC"]:
        exec("data = data_E" + VAR)

        # with ACCMIP as reference
        if data == "ACCMIP":
            exec("E" + VAR + "past = E" + VAR + "accmip.copy()")
            # follow EDGAR (PM10) variations after 2000
            for t in range(300 + 1, ind_edgar + 1):
                exec(
                    "E" + VAR + "past[t,...] = E" + VAR + "past[t-1,...] * EPM10edgar[t,...]/EPM10edgar[t-1,...]")
                exec(
                    "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
                exec(
                    "E"
                    + VAR
                    + "past[t,...] *= np.sum(E"
                    + VAR
                    + "past[t-1,...])/np.sum(E"
                    + VAR
                    + "past[t,...]) * np.sum(EPM10edgar[t,...])/np.sum(EPM10edgar[t-1,...])"
                )
            exec(
                "E" + VAR + "past[np.isnan(E" + VAR + "past)|np.isinf(E" + VAR + "past)] = 0")
            # follow EDGAR-HTAP variations after 2008
            for t in range(ind_edgar + 1, ind_cdiac + 1):
                exec(
                    "E"
                    + VAR
                    + "past[t,...] = E"
                    + VAR
                    + "past[t-1,...] * E"
                    + VAR
                    + "ehtap[t,...]/E"
                    + VAR
                    + "ehtap[t-1,...]"
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
                    + "ehtap[t,...])/np.sum(E"
                    + VAR
                    + "ehtap[t-1,...])"
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

    # cut past dataset to right length
    exec(
        "E" + VAR + "[:min(ind_cdiac,ind_final)+1,...] = E" + VAR + "past[:min(ind_cdiac,ind_final)+1,...]")

# ==================
# 4.B. FINAL DATASET
# ==================

# datasets mixed following various criteria
for VAR in ["NOX", "CO", "VOC"] + ["SO2", "NH3", "OC", "BC"]:
    exec("scen = scen_E" + VAR)

    # stop emissions
    if (scen == "stop") & (ind_final > ind_cdiac):
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
                exec("def_regI = bool(np.sum(E" + VAR + "proj[t,:,...,1:]))")
                exec("def_regJ = bool(np.sum(E" + VAR + "proj[t,1:,...,:]))")
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
                    exec("def_regI = bool(np.sum(E" + VAR + "proj[t,:,...,1:]))")
                    exec("def_regJ = bool(np.sum(E" + VAR + "proj[t,1:,...,:]))")
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
for VAR in ["NOX", "CO", "VOC"] + ["SO2", "NH3"]:
    exec(
        "del E"
        + VAR
        + "edgar,E"
        + VAR
        + "ehtap,E"
        + VAR
        + "accmip,E"
        + VAR
        + "past,E"
        + VAR
        + "proj,p_E"
        + VAR
        + "_bio"
    )
for VAR in ["OC", "BC"]:
    exec(
        "del E" + VAR + "ehtap,E" + VAR + "accmip,E" + VAR + "past,E" + VAR + "proj,p_E" + VAR + "_bio")
for VAR in ["PM10"]:
    exec("del E" + VAR + "edgar")
