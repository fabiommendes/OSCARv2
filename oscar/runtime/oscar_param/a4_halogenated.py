import csv
import os

import numpy as np

from ..oscar_data import *
from oscar.data import load_data
from .a2_methane import scale_OH
from ..oscar_data import HFC, PFC, ODS, nb_HFC, nb_PFC, nb_ODS
from ...config import dty

##################################################
#   4. HALOGENATED COMPOUNDS
##################################################

# ===============
# 4.1. ATMOSPHERE
# ===============

# conversion of HaloC from {ppt} to {kt}
alpha_HFC = 0.1765 * np.array(
    [70.0, 52.0, 120.0, 102.0, 84.0, 66.0, 170.0, 152.0, 134.0, 148.1, 252.1], dtype=dty)
alpha_PFC = 0.1765 * np.array(
    [146.1, 71.0, 88.0, 138.0, 188.0, 200.0, 238.0, 288.0, 338.0, 388.1], dtype=dty)
alpha_ODS = 0.1765 * np.array(
    [137.4, 120.9, 187.4, 170.9, 154.5, 153.8, 133.4, 86.5, 117.0, 100.5, 165.4, 209.8,
     148.9, 259.8, 94.9, 50.5],
    dtype=dty,
)

# historic HaloC from IPCC-AR5 {ppt}
# from [IPCC WG1, 2013] (annexe 2)
for VAR in ["HFC", "PFC", "ODS"]:
    HFC_ipcc = np.ones([311+1,nb_HFC], dtype=dty) * np.nan
    PFC_ipcc = np.ones([311+1,nb_PFC], dtype=dty) * np.nan
    ODS_ipcc = np.ones([311+1,nb_ODS], dtype=dty) * np.nan

for VAR in HFC:
    path = f"data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1940-2011.{VAR}.csv"
    if os.path.isfile(path):
        TMP = load_data(path)
        HFC_ipcc[240:, HFC.index(VAR)] = TMP[:, 0]

for VAR in PFC:
    if os.path.isfile(f"data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1900-2011.{VAR}.csv"):
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1900-2011." + VAR + ".csv",
                    "r")
            )],
            dtype=dty,
        )
        PFC_ipcc[200:, PFC.index(VAR)] = TMP[:, 0]

for VAR in ODS:
    if os.path.isfile(f"data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1960-2011.{VAR}.csv"):
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1960-2011." + VAR + ".csv",
                    "r")
            )],
            dtype=dty,
        )
        ODS_ipcc[260:, ODS.index(VAR)] = TMP[:, 0]

# historic HaloC from CMIP5 {ppt}
# from [Meinshausen et al., 2011]
HFC_cmip5 = np.ones([305+1,nb_HFC], dtype=dty) * np.nan
PFC_cmip5 = np.ones([305+1,nb_PFC], dtype=dty) * np.nan
ODS_cmip5 = np.ones([305+1,nb_ODS], dtype=dty) * np.nan

for VAR in HFC:
    if os.path.isfile(f"data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005.{VAR}.csv"):
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open("data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005." + VAR + ".csv",
                     "r"))],
            dtype=dty,
        )
        HFC_cmip5[65:, HFC.index(VAR)] = TMP[:, 0]

for VAR in PFC:
    if os.path.isfile(f"data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005.{VAR}.csv"):
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open("data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005." + VAR + ".csv",
                     "r"))],
            dtype=dty,
        )
        PFC_cmip5[65:, PFC.index(VAR)] = TMP[:, 0]

for VAR in ODS:
    if os.path.isfile(f"data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005.{VAR}.csv"):
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open("data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005." + VAR + ".csv",
                     "r"))],
            dtype=dty,
        )
        ODS_cmip5[65:, ODS.index(VAR)] = TMP[:, 0]

# preindustrial HaloC concentrations {ppt}
# from [IPCC WG1, 2013] (annexe 2) and [Meinshausen et al., 2011]
HFC_0 = np.zeros([nb_HFC], dtype=dty)
PFC_0 = np.zeros([nb_PFC], dtype=dty)
ODS_0 = np.zeros([nb_ODS], dtype=dty)
PFC_0[PFC.index("CF4")] = 35.0
ODS_0[ODS.index("CH3Br")] = 5.8
ODS_0[ODS.index("CH3Cl")] = 480.0

# ==============
# 4.2. CHEMISTRY
# ==============

# atmospheric OH lifetime {yr}
# from [WMO, 2011] (table 1-3)
# rescaled to follow the arbitrary rescaling of tau_CH4_OH
tau_HFC_OH = scale_OH * np.array([245, 5.5, 32, 14.3, 55, 1.6, 44.5, 253, 8.2, 9.3, 17.9],
                                 dtype=dty)
tau_PFC_OH = scale_OH * np.array(
    [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf],
    dtype=dty
)
tau_ODS_OH = scale_OH * np.array(
    [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, 6.1, 12.8, 10.7, 19.3, np.inf,
     np.inf, np.inf, np.inf, 1.5, 1.9],
    dtype=dty,
)

# stratospheric lifetime {yr}
# idem
# rescaled to follow [Prather et al., 2015]
tau_HFC_hv = 1.06 * np.array([2347, 89, 246, 232, 327, 45.4, 310, 5676, 116, 125, 157],
                             dtype=dty)
tau_PFC_hv = 1.06 * np.array(
    [3200, 500, 50000, 10000, 2600, 3200, 2600, 4100, 3100, 3000], dtype=dty)
tau_ODS_hv = 1.06 * np.array(
    [45, 100, 85, 190, 1020, 35, 39, 186, 64.9, 160, np.inf, np.inf, np.inf, np.inf,
     np.inf, np.inf], dtype=dty
)

# other sinks lifetime {yr}
# idem
tau_HFC_othr = np.array(
    [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
     np.inf], dtype=dty
)
tau_PFC_othr = np.array(
    [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf],
    dtype=dty)
tau_ODS_othr = np.array(
    [np.inf, np.inf, np.inf, np.inf, np.inf, 101, 89, np.inf, np.inf, np.inf, 16, 2.9, 65,
     20, 3, 1.4], dtype=dty
)
