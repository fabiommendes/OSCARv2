import csv
import os

import numpy as np

from .a1_regions import nb_regionJ, nb_kind, nb_regionI, \
    nb_sector, regionJ_index, regionI_index, ind_final, kindGHG_index, nb_HFC, HFC, PFC, \
    ODS, nb_PFC, nb_ODS
from .a2_greenhouse import ind_cdiac, ind_edgar
from ...config import dty, scen_Ehalo, mod_regionI, data_Ehalo, mod_DATAscen, mod_regionJ

##################################################
#   3. HALOGENATED COMPOUNDS
##################################################

EHFC = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_HFC],
                dtype=dty)  # {kt/yr}
EPFC = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_PFC],
                dtype=dty)  # {kt/yr}
EODS = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_ODS],
                dtype=dty)  # {kt/yr}

# ==========
# 3.1. EDGAR
# ==========

# load emissions from EDGAR v4.2
# see [JRC, 2011]

# HFCs
EHFCedgar = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_HFC],
                     dtype=dty)
for VAR in HFC:
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/EHaloComp_EDGAR/#DATA.EHaloComp_EDGAR.1970-"
                + str(1700 + ind_edgar)
                + "_114reg0.E"
                + VAR
                + ".csv",
                "r",
            )
        )],
        dtype=dty,
    )
    for i in range(114 + 1):
        EHFCedgar[
        270: ind_edgar + 1, regionJ_index[i], 0, kindGHG_index["HFC"], regionI_index[i],
        HFC.index(VAR)] += TMP[: ind_edgar - 270 + 1, i]

# PFCs
EPFCedgar = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_PFC],
                     dtype=dty)
for VAR in PFC:
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/EHaloComp_EDGAR/#DATA.EHaloComp_EDGAR.1970-"
                + str(1700 + ind_edgar)
                + "_114reg0.E"
                + VAR
                + ".csv",
                "r",
            )
        )],
        dtype=dty,
    )
    for i in range(114 + 1):
        EPFCedgar[
        270: ind_edgar + 1, regionJ_index[i], 0, kindGHG_index["PFC"], regionI_index[i],
        PFC.index(VAR)] += TMP[: ind_edgar - 270 + 1, i]

# =============
# 3.2. EDGAR-FT
# =============

# load emissions from EDGAR-FT v4.2-FT2010
# see [JRC, 2013]

# HFCs
EHFCeft = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_HFC],
                   dtype=dty)
for VAR in HFC:
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/EHaloComp_EDGAR-FT/#DATA.EHaloComp_EDGAR-FT.2008-2010_114reg0.E" + VAR + ".csv",
                "r")
        )],
        dtype=dty,
    )
    for i in range(114 + 1):
        EHFCeft[308: 310 + 1, regionJ_index[i], 0, kindGHG_index["HFC"], regionI_index[i],
        HFC.index(VAR)] += TMP[: 310 - 308 + 1, i]

# PFCs
EPFCeft = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_PFC],
                   dtype=dty)
for VAR in PFC:
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/EHaloComp_EDGAR-FT/#DATA.EHaloComp_EDGAR-FT.2008-2010_114reg0.E" + VAR + ".csv",
                "r")
        )],
        dtype=dty,
    )
    for i in range(114 + 1):
        EPFCeft[308: 310 + 1, regionJ_index[i], 0, kindGHG_index["PFC"], regionI_index[i],
        PFC.index(VAR)] += TMP[: 310 - 308 + 1, i]

# ==========
# 3.3. CMIP5
# ==========

# ODSs
EODScmip5 = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_ODS],
                     dtype=dty)
TMP = np.array(
    [line for line in csv.reader(
        open("data/EHaloComp_CMIP5/#DATA.EHaloComp_CMIP5.1765-2005_(16ghg).EODS.csv",
             "r"))][
    1:],
    dtype=dty,
)
lgd = [line for line in csv.reader(
    open("data/EHaloComp_CMIP5/#DATA.EHaloComp_CMIP5.1765-2005_(16ghg).EODS.csv", "r"))][
    0
]
for x in range(len(lgd)):
    EODScmip5[65: 305 + 1, 0, 0, kindGHG_index["ODS"], 0, ODS.index(lgd[x])] = TMP[: min(
        ind_final, 305) - 65 + 1, x]

# extend dataset following RCP unique projection
TMP = np.array(
    [line for line in csv.reader(
        open("data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_(16ghg).rcp85_EODS.csv",
             "r"))][
    1:],
    dtype=dty,
)
lgd = [
    line for line in csv.reader(
        open("data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_(16ghg).rcp85_EODS.csv",
             "r"))
][0]
for x in range(len(lgd)):
    EODScmip5[306: ind_cdiac + 1, 0, 0, kindGHG_index["ODS"], 0, ODS.index(lgd[x])] = TMP[
                                                                                      6: min(
                                                                                          ind_final,
                                                                                          ind_cdiac) - 300 + 1,
                                                                                      x]

# ========
# 3.4. RCP
# ========

# initialization of projected drivers
for VAR in ["HFC", "PFC", "ODS"]:
    exec(
        "E" + VAR + "proj = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI,nb_" + VAR + "], dtype=dty)")

# projection of emissions under RCP scenarios
# from [Meinshausen et al., 2011]

# HFCs
if (scen_Ehalo[:3] == "RCP") & (ind_final > ind_cdiac):
    for VAR in HFC:
        if os.path.isfile(
                "data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_5reg0.rcp"
                + scen_Ehalo[3]
                + scen_Ehalo[5]
                + "_E"
                + VAR
                + ".csv"
        ):
            TMP = np.array(
                [
                    line
                    for line in csv.reader(
                    open(
                        "data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_5reg0.rcp"
                        + scen_Ehalo[3]
                        + scen_Ehalo[5]
                        + "_E"
                        + VAR
                        + ".csv",
                        "r",
                    )
                )],
                dtype=dty,
            )
            for i in range(4 + 1):
                if (mod_regionI == "RCP5") & (mod_regionJ == "RCP5"):
                    EHFCproj[300: min(ind_final, 400) + 1, i, 0, kindGHG_index["HFC"], i,
                    HFC.index(VAR)] += TMP[: min(ind_final, 400) - 300 + 1, i]
                elif (mod_regionI == "RCP5") & (mod_regionJ != "RCP5"):
                    EHFCproj[300: min(ind_final, 400) + 1, 0, 0, kindGHG_index["HFC"], i,
                    HFC.index(VAR)] += TMP[: min(ind_final, 400) - 300 + 1, i]
                elif (mod_regionI != "RCP5") & (mod_regionJ == "RCP5"):
                    EHFCproj[300: min(ind_final, 400) + 1, i, 0, kindGHG_index["HFC"], 0,
                    HFC.index(VAR)] += TMP[: min(ind_final, 400) - 300 + 1, i]
                else:
                    EHFCproj[300: min(ind_final, 400) + 1, 0, 0, kindGHG_index["HFC"], 0,
                    HFC.index(VAR)] += TMP[: min(ind_final, 400) - 300 + 1, i]

# PFCs
if (scen_Ehalo[:3] == "RCP") & (ind_final > ind_cdiac):
    for VAR in PFC:
        if os.path.isfile(
                "data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_5reg0.rcp"
                + scen_Ehalo[3]
                + scen_Ehalo[5]
                + "_E"
                + VAR
                + ".csv"
        ):
            TMP = np.array(
                [
                    line
                    for line in csv.reader(
                    open(
                        "data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_5reg0.rcp"
                        + scen_Ehalo[3]
                        + scen_Ehalo[5]
                        + "_E"
                        + VAR
                        + ".csv",
                        "r",
                    )
                )],
                dtype=dty,
            )
            for i in range(4 + 1):
                if (mod_regionI == "RCP5") & (mod_regionJ == "RCP5"):
                    EPFCproj[300: min(ind_final, 400) + 1, i, 0, kindGHG_index["PFC"], i,
                    PFC.index(VAR)] += TMP[: min(ind_final, 400) - 300 + 1, i]
                elif (mod_regionI == "RCP5") & (mod_regionJ != "RCP5"):
                    EPFCproj[300: min(ind_final, 400) + 1, 0, 0, kindGHG_index["PFC"], i,
                    PFC.index(VAR)] += TMP[: min(ind_final, 400) - 300 + 1, i]
                elif (mod_regionI != "RCP5") & (mod_regionJ == "RCP5"):
                    EPFCproj[300: min(ind_final, 400) + 1, i, 0, kindGHG_index["PFC"], 0,
                    PFC.index(VAR)] += TMP[: min(ind_final, 400) - 300 + 1, i]
                else:
                    EPFCproj[300: min(ind_final, 400) + 1, 0, 0, kindGHG_index["PFC"], 0,
                    PFC.index(VAR)] += TMP[: min(ind_final, 400) - 300 + 1, i]

# ODSs
if (scen_Ehalo[:3] == "RCP") & (ind_final > ind_cdiac):
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_(16ghg).rcp"
                + scen_Ehalo[3]
                + scen_Ehalo[5]
                + "_EODS.csv",
                "r",
            )
        )][1:],
        dtype=dty,
    )
    lgd = [
        line
        for line in csv.reader(
            open(
                "data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_(16ghg).rcp"
                + scen_Ehalo[3]
                + scen_Ehalo[5]
                + "_EODS.csv",
                "r",
            )
        )][0]
    for x in range(len(lgd)):
        EODSproj[300: min(ind_final, 400) + 1, 0, 0, kindGHG_index["ODS"], 0,
        ODS.index(lgd[x])] = TMP[: min(ind_final, 400) - 300 + 1, x]

# =================
# 3.A. PAST DATASET
# =================

# datasets mixed following trends
for VAR in ["HFC", "PFC", "ODS"]:

    if VAR in ["HFC"]:

        # with EDGAR as reference
        if data_Ehalo == "EDGAR":
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
            # quadratic extrapolation before 1970
            # starting year of emission based on [Meinshausen et al., 2011]
            for x in range(nb_HFC):
                if HFC[x] == "HFC23":
                    for t in range(230, 270):
                        exec(
                            "E" + VAR + "past[t,...,x] = E" + VAR + "past[270,...,x] * ((t-230)/(270-230.))**2")

    elif VAR in ["PFC"]:

        # with EDGAR as reference
        if data_Ehalo == "EDGAR":
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
            # quadratic extrapolation before 1970
            # starting year of emission based on [Meinshausen et al., 2011]
            for x in range(nb_PFC):
                if PFC[x] == "SF6":
                    for t in range(250, 270):
                        exec(
                            "E" + VAR + "past[t,...,x] = E" + VAR + "past[270,...,x] * ((t-250.)/(270-250.))**2")
                elif PFC[x] == "CF4":
                    for t in range(222, 270):
                        exec(
                            "E" + VAR + "past[t,...,x] = E" + VAR + "past[270,...,x] * ((t-222.)/(270-222.))**2")
                elif PFC[x] == "C2F6":
                    for t in range(189, 270):
                        exec(
                            "E" + VAR + "past[t,...,x] = E" + VAR + "past[270,...,x] * ((t-189)/(270-189.))**2")

    elif VAR in ["ODS"]:

        # with CMIP5 as reference
        if data_Ehalo == data_Ehalo:  # FIXME: is this so?
            exec("E" + VAR + "past = E" + VAR + "cmip5.copy()")
            # linear extrapolation before 1765
            for t in range(50, 65):
                exec(
                    "E" + VAR + "past[t,...] = E" + VAR + "past[65,...] * (t-50)/float(65-50)")

    # cut past dataset to right length
    exec("E" + VAR + "[:min(ind_cdiac,ind_final)+1,...] = E" + VAR + "past[:min(ind_cdiac,ind_final)+1,...]")

# ==================
# 3.B. FINAL DATASET
# ==================

# datasets mixed following various criteria
for VAR in ["HFC", "PFC", "ODS"]:

    # stop emissions
    if (scen_Ehalo == "stop") & (ind_final > ind_cdiac):
        exec("E" + VAR + "[ind_cdiac+1:,...] = 0")

    # constant emissions
    elif (scen_Ehalo == "cst") & (ind_final > ind_cdiac):
        exec(
            "E" + VAR + "[ind_cdiac+1:,...] = E" + VAR + "[ind_cdiac,...][np.newaxis,...]")

        # RCP scenarios
    elif (scen_Ehalo[:3] == "RCP") & (ind_final > ind_cdiac):

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
                        + "[t-1,:,...,:],0) * np.sum(E"
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
for VAR in ["HFC", "PFC"]:
    exec("del E" + VAR + "edgar,E" + VAR + "eft,E" + VAR + "past,E" + VAR + "proj")
for VAR in ["ODS"]:
    exec("del E" + VAR + "cmip5,E" + VAR + "past,E" + VAR + "proj")
