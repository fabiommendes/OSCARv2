import csv

import numpy as np

from .a1_carbon import CO2_0, load_data
from .a2_methane import CH4_0
from .a3_nitrous_oxide import N2O_0
from ..oscar_data import nb_regionI, nb_biome, regionI_index, biome_index
from ...config import dty, mod_O3Tradeff, mod_SO4radeff, mod_POAradeff, mod_BCradeff, \
    mod_NO3radeff, mod_SOAradeff, mod_BCadjust, mod_CLOUDsolub, mod_CLOUDerf, \
    mod_CLOUDpreind, mod_ALBBCrf, mod_ALBBCreg, mod_ALBBCwarm, mod_ALBLCalb, \
    mod_ALBLCflux, mod_ALBLCcover, mod_ALBLCwarm, mod_O3Sradeff

##################################################
#   7. RADIATIVE FORCING
##################################################

# ====================
# 7.A. Reconstructions
# ====================

# historic RF from IPCC-AR5 {W/m2}
# from [IPCC WG1, 2013] annexe 2
RF_ipcc = np.zeros([311 + 1], dtype=dty)
RF_ipcc[:50] = np.nan

RF_WMGHG_ipcc = RF_ipcc.copy()
RF_O3_ipcc = RF_ipcc.copy()
RF_AER_ipcc = RF_ipcc.copy()
RF_Alb_ipcc = RF_ipcc.copy()

TMP = load_data("data/Historic_IPCC-AR5/#DATA.Historic_IPCC-AR5.1750-2011_(11for).RF.csv", slice=1)
path = "data/Historic_IPCC-AR5/#DATA.Historic_IPCC-AR5.1750-2011_(11for).RF.csv"
lgd = [line for line in csv.reader(open(path, "r"))][0]

for x in range(len(lgd)):
    if lgd[x] in ["CO2", "GHG Other", "H2O (Strat)"]:
        RF_WMGHG_ipcc[51:] += TMP[1:, x]
        RF_ipcc[51:] += TMP[1:, x]
    elif lgd[x] in ["O3 (Trop)", "O3 (Strat)"]:
        RF_O3_ipcc[51:] += TMP[1:, x]
        RF_ipcc[51:] += TMP[1:, x]
    elif lgd[x] in ["Aerosol (Total)"]:
        RF_AER_ipcc[51:] += TMP[1:, x]
        RF_ipcc[51:] += TMP[1:, x]
    elif lgd[x] in ["LUC", "BC Snow"]:
        RF_Alb_ipcc[51:] += TMP[1:, x]
        RF_ipcc[51:] += TMP[1:, x]
    elif lgd[x] in ["Volcano"]:
        RF_ipcc[51:] += TMP[1:, x] - np.mean(TMP[1:, x])
    else:
        RF_ipcc[51:] += TMP[1:, x]

# historic RF from CMIP5 {W/m2}
# from [Meinshausen et al., 2011]
RF_cmip5 = np.zeros([305 + 1], dtype=dty)
RF_cmip5[:65] = np.nan

path = "data/Historic_CMIP5/#DATA.Historic_CMIP5.1765-2005_(19for).RF.csv"
TMP = load_data(path, slice=1)

path = "data/Historic_CMIP5/#DATA.Historic_CMIP5.1765-2005_(19for).RF.csv"
lgd = [line for line in csv.reader(open(path, "r"))][0]

for x in range(len(lgd)):
    if lgd[x] == "VOLC":
        RF_cmip5[66:] += TMP[1:, x] - np.mean(TMP[1:, x])
    else:
        RF_cmip5[66:] += TMP[1:, x]

# load RCP radiative forcing {W/m2}
# from [Meinshausen et al., 2011]
RF_rcp = np.zeros([800 + 1, 6], dtype=dty)
RF_rcp[:300] = np.nan
n = -1
for rcp in ["rcp26", "rcp45", "rcp60", "rcp85", "rcp45to26", "rcp60to45"]:
    n += 1
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/Scenario_ECP/#DATA.Scenario_ECP.2000-2500_(19for)." + rcp + "_RF.csv",
                "r")
        )][1:],
        dtype=dty,
    )
    TMP2 = np.array(
        [line for line in csv.reader(
            open("data/Historic_CMIP5/#DATA.Historic_CMIP5.1765-2005_(19for).RF.csv",
                 "r"))][
        1:],
        dtype=dty,
    )
    lgd = [
        line
        for line in csv.reader(open(
            "data/Scenario_ECP/#DATA.Scenario_ECP.2000-2500_(19for)." + rcp + "_RF.csv",
            "r"))][0]
    for x in range(len(lgd)):
        RF_rcp[300:, n] += TMP[:, x]
        if lgd[x] == "VOLC":
            RF_rcp[:7] -= np.mean(TMP2[1:, x])


# =====================
# 7.1. GREENHOUSE GASES
# =====================

# radiative forcing functions {W/m2}
# from IPCC-TAR [Ramaswamy et al., 2001]
def f_RF_CO2(D_CO2):
    RF = 5.35 * np.log(1 + D_CO2 / CO2_0)
    return np.array(RF, dtype=dty)


def f_RF_CH4(D_CH4):
    RF = 0.036 * (np.sqrt(CH4_0 + D_CH4) - np.sqrt(CH4_0))
    return np.array(RF, dtype=dty)


def f_RF_H2Os(D_CH4_lag):
    RF = 0.15 * 0.036 * (np.sqrt(CH4_0 + D_CH4_lag) - np.sqrt(CH4_0))
    return np.array(RF, dtype=dty)


def f_RF_N2O(D_N2O):
    RF = 0.12 * (np.sqrt(N2O_0 + D_N2O) - np.sqrt(N2O_0))
    return np.array(RF, dtype=dty)


def f_RF_overlap(D_CH4, D_N2O):
    RF = 0.47 * np.log(
        1
        + 2.01e-5 * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 0.75
        + 5.31e-15 * (CH4_0 + D_CH4) * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 1.52
    )
    RF -= 0.47 * np.log(1 + 2.01e-5 * (CH4_0 * N2O_0) ** 0.75 + 5.31e-15 * CH4_0 * (
            CH4_0 * N2O_0) ** 1.52)
    return np.array(RF, dtype=dty)


def df_RF_overlap_dCH4(D_CH4, D_N2O):
    RF = 0.47 * (
            2.01e-5 * 0.75 * (CH4_0 + D_CH4) ** (0.75 - 1) * (N2O_0 + D_N2O) ** (0.75)
            + 5.31e-15 * (1.52 + 1) * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 1.52
    )
    RF /= (
            1
            + 2.01e-5 * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 0.75
            + 5.31e-15 * (CH4_0 + D_CH4) * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 1.52
    )
    return np.array(RF, dtype=dty)


def df_RF_overlap_dN2O(D_CH4, D_N2O):
    RF = 0.47 * (
            2.01e-5 * 0.75 * (CH4_0 + D_CH4) ** (0.75) * (N2O_0 + D_N2O) ** (0.75 - 1)
            + 5.31e-15 * 1.52 * (CH4_0 + D_CH4) ** (1.52 + 1) * (N2O_0 + D_N2O) ** (
                    1.52 - 1)
    )
    RF /= (
            1
            + 2.01e-5 * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 0.75
            + 5.31e-15 * (CH4_0 + D_CH4) * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 1.52
    )
    return np.array(RF, dtype=dty)


def df_RF_overlap(D_CH4, D_N2O):
    df_1 = df_RF_overlap_dCH4(D_CH4, D_N2O) * D_CH4
    df_2 = df_RF_overlap_dN2O(D_CH4, D_N2O) * D_N2O
    df_tot = df_1 + df_2
    return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot)]


# radiative efficiency of halogenated compounds {{W/m2}/ppt}
# from IPCC-AR5 [Myrhe et al., 2013] (table 8.A.1)
radeff_HFC = 1e-3 * np.array(
    [0.18, 0.11, 0.23, 0.16, 0.16, 0.10, 0.26, 0.24, 0.24, 0.22, 0.42], dtype=dty)
radeff_PFC = 1e-3 * np.array([0.57, 0.20, 0.09, 0.25, 0.28, 0.32, 0.36, 0.41, 0.44, 0.50],
                             dtype=dty)
radeff_ODS = 1e-3 * np.array(
    [0.26, 0.32, 0.30, 0.31, 0.20, 0.17, 0.07, 0.21, 0.16, 0.19, 0.29, 0.27, 0.30, 0.31,
     0.004, 0.01], dtype=dty
)

# ==========
# 7.2. OZONE
# ==========

# radiative efficiency of tropospheric O3 {{W/m2}/DU}
# from IPCC-AR5 [Myhre et al., 2013]
if mod_O3Tradeff == "IPCC-AR5":
    radeff_O3t = np.array([0.042], dtype=dty)
# from IPCC-AR4 [Forster et al., 2007]
elif mod_O3Tradeff == "IPCC-AR4":
    radeff_O3t = np.array([0.032], dtype=dty)
# from ACCMIP [Stevenson et al., 2013] (table 3)
elif mod_O3Tradeff == "mean-ACCMIP":
    radeff_O3t = np.array([0.377 / 8.9], dtype=dty)
elif mod_O3Tradeff == "CESM-CAM-superfast":
    radeff_O3t = np.array([0.446 / 10.0], dtype=dty)
elif mod_O3Tradeff == "CICERO-OsloCTM2":
    radeff_O3t = np.array([0.401 / 9.3], dtype=dty)
elif mod_O3Tradeff == "CMAM":
    radeff_O3t = np.array([0.322 / 7.6], dtype=dty)
elif mod_O3Tradeff == "EMAC":
    radeff_O3t = np.array([0.460 / 10.8], dtype=dty)
elif mod_O3Tradeff == "GEOSCCM":
    radeff_O3t = np.array([0.387 / 8.7], dtype=dty)
elif mod_O3Tradeff == "GFDL-AM3":
    radeff_O3t = np.array([0.423 / 10.3], dtype=dty)
elif mod_O3Tradeff == "GISS-E2-R":
    radeff_O3t = np.array([0.314 / 8.3], dtype=dty)
elif mod_O3Tradeff == "GISS-E2-R-TOMAS":
    radeff_O3t = np.array([0.333 / 8.7], dtype=dty)
elif mod_O3Tradeff == "HadGEM2":
    radeff_O3t = np.array([0.303 / 7.3], dtype=dty)
elif mod_O3Tradeff == "LMDzORINCA":
    radeff_O3t = np.array([0.351 / 8.2], dtype=dty)
elif mod_O3Tradeff == "MIROC-CHEM":
    radeff_O3t = np.array([0.402 / 9.2], dtype=dty)
elif mod_O3Tradeff == "MOCAGE":
    radeff_O3t = np.array([0.219 / 4.8], dtype=dty)
elif mod_O3Tradeff == "NCAR-CAM-35":
    radeff_O3t = np.array([0.433 / 10.2], dtype=dty)
elif mod_O3Tradeff == "STOC-HadAM3":
    radeff_O3t = np.array([0.437 / 10.5], dtype=dty)
elif mod_O3Tradeff == "UM-CAM":
    radeff_O3t = np.array([0.376 / 8.7], dtype=dty)
elif mod_O3Tradeff == "TM5":
    radeff_O3t = np.array([0.422 / 10.0], dtype=dty)

# radiative efficiency of stratospheric O3 {{W/m2}/DU}
# from IPCC-AR4 [Forster et al., 2007]
if mod_O3Sradeff == "IPCC-AR4":
    radeff_O3s = np.array([0.004], dtype=dty)
# from ACCENT [Gauss et al., 2006] (tables 4 & 6)
elif mod_O3Sradeff == "mean-ACCENT":
    radeff_O3s = np.array([-0.058 / -13.9], dtype=dty)
elif mod_O3Sradeff == "ULAQ":
    radeff_O3s = np.array([-0.059 / -12.6], dtype=dty)
elif mod_O3Sradeff == "DLR-E39C":
    radeff_O3s = np.array([-0.027 / -16.1], dtype=dty)
elif mod_O3Sradeff == "NCAR-MACCM":
    radeff_O3s = np.array([-0.019 / -12.7], dtype=dty)
elif mod_O3Sradeff == "CHASER":
    radeff_O3s = np.array([-0.126 / -14.1], dtype=dty)

# =============
# 7.3. AEROSOLS
# =============

# -------------
# 7.3.1. Direct
# -------------

# radiative efficiency of sulfate aerosols {{W/m2}/Tg}
# from AeroCom2 [Myhre et al., 2013] (table 4)
if mod_SO4radeff == "mean-AeroCom2":
    radeff_SO4 = 1e12 / 510_072e9 * np.array([-185], dtype=dty)
elif mod_SO4radeff == "BCC":
    radeff_SO4 = 1e12 / 510_072e9 * np.array([-108], dtype=dty)
elif mod_SO4radeff == "CAM4-Oslo":
    radeff_SO4 = 1e12 / 510_072e9 * np.array([-173], dtype=dty)
elif mod_SO4radeff == "CAM-51":
    radeff_SO4 = 1e12 / 510_072e9 * np.array([-104], dtype=dty)
elif mod_SO4radeff == "GEOS-CHEM":
    radeff_SO4 = 1e12 / 510_072e9 * np.array([-123], dtype=dty)
elif mod_SO4radeff == "GISS-MATRIX":
    radeff_SO4 = 1e12 / 510_072e9 * np.array([-196], dtype=dty)
elif mod_SO4radeff == "GISS-modelE":
    radeff_SO4 = 1e12 / 510_072e9 * np.array([-307], dtype=dty)
elif mod_SO4radeff == "GMI":
    radeff_SO4 = 1e12 / 510_072e9 * np.array([-195], dtype=dty)
elif mod_SO4radeff == "GOCART":
    radeff_SO4 = 1e12 / 510_072e9 * np.array([-238], dtype=dty)
elif mod_SO4radeff == "HadGEM2":
    radeff_SO4 = 1e12 / 510_072e9 * np.array([-193], dtype=dty)
elif mod_SO4radeff == "IMPACT-Umich":
    radeff_SO4 = 1e12 / 510_072e9 * np.array([-113], dtype=dty)
elif mod_SO4radeff == "INCA":
    radeff_SO4 = 1e12 / 510_072e9 * np.array([-180], dtype=dty)
elif mod_SO4radeff == "MPIHAM":
    radeff_SO4 = 1e12 / 510_072e9 * np.array([-125], dtype=dty)
elif mod_SO4radeff == "NCAR-CAM-35":
    radeff_SO4 = 1e12 / 510_072e9 * np.array([-354], dtype=dty)
elif mod_SO4radeff == "OsloCTM2":
    radeff_SO4 = 1e12 / 510_072e9 * np.array([-192], dtype=dty)
elif mod_SO4radeff == "SPRINTARS":
    radeff_SO4 = 1e12 / 510_072e9 * np.array([-172], dtype=dty)

# radiative efficiency of primary organic aerosols {{W/m2}/Tg}
# from AeroCom2 [Myhre et al., 2013] (table 6)
if mod_POAradeff == "mean-AeroCom2":
    radeff_POA = 1e12 / 510_072e9 * np.array([-113], dtype=dty)
elif mod_POAradeff == "BCC":
    radeff_POA = 1e12 / 510_072e9 * np.array([-97], dtype=dty)
elif mod_POAradeff == "CAM4-Oslo":
    radeff_POA = 1e12 / 510_072e9 * np.array([-118], dtype=dty)
elif mod_POAradeff == "CAM-51":
    radeff_POA = 1e12 / 510_072e9 * np.array([-69], dtype=dty)
elif mod_POAradeff == "GEOS-CHEM":
    radeff_POA = 1e12 / 510_072e9 * np.array([-95], dtype=dty)
elif mod_POAradeff == "GISS-MATRIX":
    radeff_POA = 1e12 / 510_072e9 * np.array([-129], dtype=dty)
elif mod_POAradeff == "GISS-modelE":
    radeff_POA = 1e12 / 510_072e9 * np.array([-76], dtype=dty)
elif mod_POAradeff == "GMI":
    radeff_POA = 1e12 / 510_072e9 * np.array([-189], dtype=dty)
elif mod_POAradeff == "GOCART":
    radeff_POA = 1e12 / 510_072e9 * np.array([-144], dtype=dty)
elif mod_POAradeff == "HadGEM2":
    radeff_POA = 1e12 / 510_072e9 * np.array([-145], dtype=dty)
elif mod_POAradeff == "IMPACT-Umich":
    radeff_POA = 1e12 / 510_072e9 * np.array([-141], dtype=dty)
elif mod_POAradeff == "INCA":
    radeff_POA = 1e12 / 510_072e9 * np.array([-76], dtype=dty)
elif mod_POAradeff == "MPIHAM":
    radeff_POA = 1e12 / 510_072e9 * np.array([-41], dtype=dty)
elif mod_POAradeff == "NCAR-CAM-35":
    radeff_POA = 1e12 / 510_072e9 * np.array([-48], dtype=dty)
elif mod_POAradeff == "OsloCTM2":
    radeff_POA = 1e12 / 510_072e9 * np.array([-165], dtype=dty)
elif mod_POAradeff == "SPRINTARS":
    radeff_POA = 1e12 / 510_072e9 * np.array([-102], dtype=dty)

# radiative efficiency of black carbon aerosols {{W/m2}/Tg}
# from AeroCom2 [Myhre et al., 2013] (table 5)
if mod_BCradeff == "mean-AeroCom2":
    radeff_BC = 1e12 / 510_072e9 * np.array([1438], dtype=dty)
elif mod_BCradeff == "BCC":
    radeff_BC = 1e12 / 510_072e9 * np.array([650], dtype=dty)
elif mod_BCradeff == "CAM4-Oslo":
    radeff_BC = 1e12 / 510_072e9 * np.array([1763], dtype=dty)
elif mod_BCradeff == "CAM-51":
    radeff_BC = 1e12 / 510_072e9 * np.array([2661], dtype=dty)
elif mod_BCradeff == "GEOS-CHEM":
    radeff_BC = 1e12 / 510_072e9 * np.array([1067], dtype=dty)
elif mod_BCradeff == "GISS-MATRIX":
    radeff_BC = 1e12 / 510_072e9 * np.array([2484], dtype=dty)
elif mod_BCradeff == "GISS-modelE":
    radeff_BC = 1e12 / 510_072e9 * np.array([1253], dtype=dty)
elif mod_BCradeff == "GMI":
    radeff_BC = 1e12 / 510_072e9 * np.array([1208], dtype=dty)
elif mod_BCradeff == "GOCART":
    radeff_BC = 1e12 / 510_072e9 * np.array([874], dtype=dty)
elif mod_BCradeff == "HadGEM2":
    radeff_BC = 1e12 / 510_072e9 * np.array([612], dtype=dty)
elif mod_BCradeff == "IMPACT-Umich":
    radeff_BC = 1e12 / 510_072e9 * np.array([1467], dtype=dty)
elif mod_BCradeff == "INCA":
    radeff_BC = 1e12 / 510_072e9 * np.array([1160], dtype=dty)
elif mod_BCradeff == "MPIHAM":
    radeff_BC = 1e12 / 510_072e9 * np.array([1453], dtype=dty)
elif mod_BCradeff == "NCAR-CAM-35":
    radeff_BC = 1e12 / 510_072e9 * np.array([1364], dtype=dty)
elif mod_BCradeff == "OsloCTM2":
    radeff_BC = 1e12 / 510_072e9 * np.array([2161], dtype=dty)
elif mod_BCradeff == "SPRINTARS":
    radeff_BC = 1e12 / 510_072e9 * np.array([1322], dtype=dty)

# radiative efficiency of nitrate aerosols {{W/m2}/Tg}
# from AeroCom2 [Myhre et al., 2013] (table 8)
if mod_NO3radeff == "mean-AeroCom2":
    radeff_NO3 = 1e12 / 510_072e9 * np.array([-166], dtype=dty)
elif mod_NO3radeff == "GEOS-CHEM":
    radeff_NO3 = 1e12 / 510_072e9 * np.array([-136], dtype=dty)
elif mod_NO3radeff == "GISS-MATRIX":
    radeff_NO3 = 1e12 / 510_072e9 * np.array([-240], dtype=dty)
elif mod_NO3radeff == "GMI":
    radeff_NO3 = 1e12 / 510_072e9 * np.array([-103], dtype=dty)
elif mod_NO3radeff == "HadGEM2":
    radeff_NO3 = 1e12 / 510_072e9 * np.array([-249], dtype=dty)
elif mod_NO3radeff == "IMPACT-Umich":
    radeff_NO3 = 1e12 / 510_072e9 * np.array([-155], dtype=dty)
elif mod_NO3radeff == "INCA":
    radeff_NO3 = 1e12 / 510_072e9 * np.array([-110], dtype=dty)
elif mod_NO3radeff == "NCAR-CAM-35":
    radeff_NO3 = 1e12 / 510_072e9 * np.array([-91], dtype=dty)
elif mod_NO3radeff == "OsloCTM2":
    radeff_NO3 = 1e12 / 510_072e9 * np.array([-173], dtype=dty)

# radiative efficiency of secondary organic aerosols {{W/m2}/Tg}
# from AeroCom2 [Myhre et al., 2013] (table 7)
if mod_SOAradeff == "mean-AeroCom2":
    radeff_SOA = 1e12 / 510_072e9 * np.array([-122], dtype=dty)
elif mod_SOAradeff == "CAM-51":
    radeff_SOA = 1e12 / 510_072e9 * np.array([-45], dtype=dty)
elif mod_SOAradeff == "GEOS-CHEM":
    radeff_SOA = 1e12 / 510_072e9 * np.array([-45], dtype=dty)
elif mod_SOAradeff == "IMPACT-Umich":
    radeff_SOA = 1e12 / 510_072e9 * np.array([-218], dtype=dty)
elif mod_SOAradeff == "MPIHAM":
    radeff_SOA = 1e12 / 510_072e9 * np.array([-139], dtype=dty)
elif mod_SOAradeff == "OsloCTM2":
    radeff_SOA = 1e12 / 510_072e9 * np.array([-161], dtype=dty)

# radiative efficiency of natural aerosols {{W/m2}/Tg}
# set to zero in this version
radeff_DUST = 1e12 / 510_072e9 * np.array([0], dtype=dty)
radeff_SALT = 1e12 / 510_072e9 * np.array([0], dtype=dty)

# ---------------
# 7.3.2. Indirect
# ---------------

# semi-direct effect (adjustements induced by the direct RF of BC)
# best-guess from IPCC AR5 [Boucher et al., 2013]
if mod_BCadjust == "Boucher2013":
    k_BC_adjust = np.array([-0.1 / 0.6], dtype=dty)
# variations from [Lohmann et al., 2010] (figure 2; data provided by author)
elif mod_BCadjust == "CSIRO":
    k_BC_adjust = np.array([-0.37], dtype=dty) * (-0.1 / 0.6) / -0.111
elif mod_BCadjust == "GISS":
    k_BC_adjust = np.array([-0.225], dtype=dty) * (-0.1 / 0.6) / -0.111
elif mod_BCadjust == "HadGEM2":
    k_BC_adjust = np.array([-0.13], dtype=dty) * (-0.1 / 0.6) / -0.111
elif mod_BCadjust == "ECHAM5":
    k_BC_adjust = np.array([0.05], dtype=dty) * (-0.1 / 0.6) / -0.111
elif mod_BCadjust == "ECMWF":
    k_BC_adjust = np.array([0.12], dtype=dty) * (-0.1 / 0.6) / -0.111

# solubility of aerosols for the aerosol-cloud interaction
# from [Hansen et al., 2005]
if mod_CLOUDsolub == "Hansen2005":
    solub_SO4 = np.array([1.0], dtype=dty)
    solub_POA = np.array([0.8], dtype=dty)
    solub_BC = np.array([0.7], dtype=dty)  # average of FF (0.6) and BB (0.8)
    solub_NO3 = np.array([1.0], dtype=dty)
    solub_SOA = np.array([0.8], dtype=dty)  # assumed
    solub_DUST = np.array([0.0], dtype=dty)
    solub_SALT = np.array([1.0], dtype=dty)
# from [Lamarque et al., 2011] (RCP database)
elif mod_CLOUDsolub == "Lamarque2011":
    solub_SO4 = np.array([1.0], dtype=dty)
    solub_POA = np.array([0.86], dtype=dty)
    solub_BC = np.array([0.80], dtype=dty)
    solub_NO3 = np.array([1.0], dtype=dty)
    solub_SOA = np.array([1.0], dtype=dty)
    solub_DUST = np.array([0.12], dtype=dty)
    solub_SALT = np.array([0.05], dtype=dty)

# ERF over 1850-2000 for the aerosol-cloud interaction
# from ACCMIP [Shindell et al., 2013] (table 7)
# rescaled to the best guess of IPCC AR5 [Boucher et al., 2013]
if mod_CLOUDerf == "mean-ACCMIP":
    RF_ref1 = np.array([-0.84], dtype=dty) * -0.45 / -0.84
elif mod_CLOUDerf == "CSIRO-Mk360":
    RF_ref1 = np.array([-0.99], dtype=dty) * -0.45 / -0.84
elif mod_CLOUDerf == "GFDL-AM3":
    RF_ref1 = np.array([-0.82], dtype=dty) * -0.45 / -0.84
elif mod_CLOUDerf == "GISS-E2-R":
    RF_ref1 = np.array([-0.61], dtype=dty) * -0.45 / -0.84
elif mod_CLOUDerf == "HadGEM2":
    RF_ref1 = np.array([-0.89], dtype=dty) * -0.45 / -0.84
elif mod_CLOUDerf == "LMDzORINCA":
    RF_ref1 = np.array([-0.21], dtype=dty) * -0.45 / -0.84
elif mod_CLOUDerf == "MIROC-CHEM":
    RF_ref1 = np.array([-1.12], dtype=dty) * -0.45 / -0.84
elif mod_CLOUDerf == "NCAR-CAM-51":
    RF_ref1 = np.array([-1.22], dtype=dty) * -0.45 / -0.84

# soluble aerosol load for the aerosol-cloud interaction
# load pre-processed ACCMIP data
TMP = np.array(
    [
        line
        for line in csv.reader(
        open(
            "data/AeroCloud_ACCMIP/#DATA.AeroCloud_" + mod_CLOUDerf + ".(2yr)_(7aer).LOAD.csv",
            "r")
    )][1:],
    dtype=dty,
)[:, 1:]
path = f"data/AeroCloud_ACCMIP/#DATA.AeroCloud_{mod_CLOUDerf}.(2yr)_(7aer).LOAD.csv"
lgd = [line for line in csv.reader(open(path, "r"))][0][1:]
AER_ref0 = 0
AER_ref1 = 0
for n in range(len(lgd)):
    if not np.isnan(TMP[0, n]):
        AER_ref0 += TMP[0,n] * globals()["solub_" + lgd[n]]
        AER_ref1 += TMP[1,n] * globals()["solub_" + lgd[n]]

# calculate parameters for aerosol-cloud interaction
# intensity of indirect effect {W/m2}
# based on a logarithmic formulation [e.g. Gultepe and Isaac, 1999]
Phi_0 = RF_ref1 / np.log(AER_ref1 / AER_ref0)

# preindustrial load of soluble aerosols {Tg}
# reduced by a factor from [Carslaw et al., 2013] and two arbitrary variations
if mod_CLOUDpreind == "median":
    AERh_0 = AER_ref0 * np.exp(1 * (1.42 - 1.30) / Phi_0)
elif mod_CLOUDpreind == "high":
    AERh_0 = AER_ref0 * np.exp(0 * (1.42 - 1.30) / Phi_0)
elif mod_CLOUDpreind == "low":
    AERh_0 = AER_ref0 * np.exp(2 * (1.42 - 1.30) / Phi_0)

# ----------------
# 7.3.3. Volcanoes
# ----------------

# warming efficacy of volcano forcing {.}
# based on [Gregory et al., 2016]
warmeff_volc = np.array([0.6], dtype=dty)

# ===========
# 7.4. ALBEDO
# ===========

# -------------------
# 7.4.1. Black Carbon
# -------------------

# read region distribution
TMP = np.array(
    [line for line in csv.reader(
        open("data/RegDiv_Reddy2007/#DATA.RegDiv_Reddy2007.114reg1_(9reg0).AREA.csv",
             "r"))][
    1:],
    dtype=dty,
)
p_reg9 = np.zeros([nb_regionI, 9 + 1], dtype=dty)
for i in range(1, 114 + 1):
    p_reg9[regionI_index[i], :] += TMP[i - 1, :]
p_reg9 /= np.sum(p_reg9, 1)[:, np.newaxis]
p_reg9[np.isnan(p_reg9) | np.isinf(p_reg9)] = 0

# regional effect coefficient {.}
# from [Reddy and Boucher, 2007] (table 1)
if mod_ALBBCreg == "Reddy2007":
    w_reg_BCsnow = np.array(
        [1, 1 / 6.0, 11 / 11.0, 1 / 10.0, 63 / 12.0, 2 / 3.0, 2 / 13.0, 17 / 43.0,
         1 / 1.0, 2 / 1.0], dtype=dty
    )

# global normalized radiative effect {{W/m2}/{Tg/yr}}
# from ACCMIP [Lee et al., 2013] (tab. 3 & fig. 15)
# rescaled to the best guess of IPCC AR5 [Boucher et al., 2013] (using table 7.1a)
if mod_ALBBCrf == "mean-ACCMIP":
    radeff_BCsnow = np.array([0.0146 / (7.9 - 3.2)], dtype=dty) * (0.04 / 4.8) / (
            0.0146 / (7.9 - 3.2))
elif mod_ALBBCrf == "CICERO-OsloCTM2":
    radeff_BCsnow = np.array([0.0131 / (7.8 - 3.1)], dtype=dty) * (0.04 / 4.8) / (
            0.0146 / (7.9 - 3.2))
elif mod_ALBBCrf == "GFDL-AM3":
    radeff_BCsnow = np.array([0.0130 / (7.8 - 3.1)], dtype=dty) * (0.04 / 4.8) / (
            0.0146 / (7.9 - 3.2))
elif mod_ALBBCrf == "GISS-E2-R":
    radeff_BCsnow = np.array([0.0142 / (8.8 - 4.0)], dtype=dty) * (0.04 / 4.8) / (
            0.0146 / (7.9 - 3.2))
elif mod_ALBBCrf == "GISS-E2-R-TOMAS":
    radeff_BCsnow = np.array([0.0175 / (7.8 - 3.1)], dtype=dty) * (0.04 / 4.8) / (
            0.0146 / (7.9 - 3.2))
elif mod_ALBBCrf == "HadGEM2":
    radeff_BCsnow = np.array([0.0133 / (7.8 - 3.1)], dtype=dty) * (0.04 / 4.8) / (
            0.0146 / (7.9 - 3.2))
elif mod_ALBBCrf == "MIROC-CHEM":
    radeff_BCsnow = np.array([0.0173 / (7.7 - 3.0)], dtype=dty) * (0.04 / 4.8) / (
            0.0146 / (7.9 - 3.2))
elif mod_ALBBCrf == "NCAR-CAM-35":
    radeff_BCsnow = np.array([0.0143 / (7.8 - 3.1)], dtype=dty) * (0.04 / 4.8) / (
            0.0146 / (7.9 - 3.2))
elif mod_ALBBCrf == "NCAR-CAM-51":
    radeff_BCsnow = np.array([0.0141 / (7.8 - 3.1)], dtype=dty) * (0.04 / 4.8) / (
            0.0146 / (7.9 - 3.2))

# warming efficacy of black carbon on snow {.}
# from IPCC AR5 [Boucher et al., 2013] (sect. 7.5.2.3)
if mod_ALBBCwarm == "median":
    warmeff_BCsnow = np.array([3.0], dtype=dty)
elif mod_ALBBCwarm == "low":
    warmeff_BCsnow = np.array([2.0], dtype=dty)
elif mod_ALBBCwarm == "high":
    warmeff_BCsnow = np.array([4.0], dtype=dty)

# -----------------
# 7.4.2. Land-Cover
# -----------------

# basic biomes of aggregation
bio = ["des", "for", "shr", "gra", "cro", "pas", "urb"]

# load pre-processed albedo climatology
TMP = np.array(
    [
        line
        for line in csv.reader(
        open(
            "data/Albedo_"
            + mod_ALBLCalb
            + "/#DATA.Albedo_"
            + mod_ALBLCalb
            + "_"
            + mod_ALBLCflux
            + "_"
            + mod_ALBLCcover
            + ".114reg1_7bio.alb.csv",
            "r",
        )
    )],
    dtype=dty,
)
TMP2 = np.array(
    [
        line
        for line in csv.reader(
        open(
            "data/Albedo_"
            + mod_ALBLCalb
            + "/#DATA.Albedo_"
            + mod_ALBLCalb
            + "_"
            + mod_ALBLCflux
            + "_"
            + mod_ALBLCcover
            + ".114reg1_7bio.RSDS.csv",
            "r",
        )
    )],
    dtype=dty,
)
alpha_alb = np.zeros([nb_regionI, nb_biome])
RSDS_alb = np.zeros([nb_regionI, nb_biome])
for i in range(1, 114 + 1):
    for b in range(len(bio)):
        alpha_alb[regionI_index[i], biome_index[bio[b]]] += TMP[i - 1, b] * TMP2[i - 1, b]
        RSDS_alb[regionI_index[i], biome_index[bio[b]]] += TMP2[i - 1, b]
alpha_alb /= RSDS_alb
alpha_alb[np.isnan(alpha_alb) | np.isinf(alpha_alb)] = 0

# set pasture albedo
if np.sum(alpha_alb[:, biome_index["pas"]]) == 0:
    for i in range(1, 114 + 1):
        alpha_alb[regionI_index[i], biome_index["pas"]] += (
                0.6 * TMP[i - 1, bio.index("gra")] * TMP2[i - 1, bio.index("gra")]
                + 0.4 * TMP[i - 1, bio.index("des")] * TMP2[i - 1, bio.index("des")]
        )
        RSDS_alb[regionI_index[i], biome_index["pas"]] += (
                0.6 * TMP2[i - 1, bio.index("gra")] + 0.4 * TMP2[i - 1, bio.index("des")]
        )
    alpha_alb[:, biome_index["pas"]] /= RSDS_alb[:, biome_index["pas"]]
    alpha_alb[np.isnan(alpha_alb) | np.isinf(alpha_alb)] = 0

# load pre-processed radiation climatology
TMP = np.array(
    [
        line
        for line in csv.reader(
        open(
            "data/RadFlux_" + mod_ALBLCflux + "/#DATA.RadFlux_" + mod_ALBLCflux + ".114reg1.rsds.csv",
            "r")
    )],
    dtype=dty,
)
TMP2 = np.array(
    [
        line
        for line in csv.reader(
        open(
            "data/RadFlux_" + mod_ALBLCflux + "/#DATA.RadFlux_" + mod_ALBLCflux + ".114reg1.AREA.csv",
            "r")
    )],
    dtype=dty,
)
rsds_alb = np.zeros([nb_regionI])
AREA_alb = np.zeros([nb_regionI])
for i in range(1, 114 + 1):
    rsds_alb[regionI_index[i]] += TMP[i - 1, 0] * TMP2[i - 1, 0]
    AREA_alb[regionI_index[i]] += TMP2[i - 1, 0]
rsds_alb /= AREA_alb
rsds_alb[np.isnan(rsds_alb) | np.isinf(rsds_alb)] = 0

# upward transmittance from [Lenton and Vaughan, 2009]
p_trans = np.array([-0.854], dtype=dty)

# final albedo parameters {{W/m2}/Mha}
alpha_LCC = p_trans * alpha_alb * rsds_alb[:, np.newaxis] / (510_072e9 / 1e10)

# warming efficacy of land-cover change albedo effect {.}
# from [Bright et al., 2015] (tab. 7)
if mod_ALBLCwarm == "Hansen2005":
    warmeff_LCC = np.array([1.02], dtype=dty)
elif mod_ALBLCwarm == "Davin2007":
    warmeff_LCC = np.array([0.5], dtype=dty)
elif mod_ALBLCwarm == "Davin2010":
    warmeff_LCC = np.array([0.78], dtype=dty)
elif mod_ALBLCwarm == "Jones2013":
    warmeff_LCC = np.array([0.79], dtype=dty)
