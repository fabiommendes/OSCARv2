import csv
import os

import numpy as np
from scipy.optimize import fmin, fsolve
from scipy.special import gammainc

from ..oscar_data import *
from ..oscar_data import nb_regionI, nb_biome, regionI_index, biome_index
from ...config import dty, mod_OSNKstruct, mod_OSNKchem, mod_OSNKtrans, PI_1750, \
    data_LULCC, mod_LSNKcover, ind_final, mod_EHWPspeed, mod_EHWPtau, mod_EHWPbb, \
    mod_biomeURB, mod_biomeV3, mod_biomeSHR, mod_ELUCagb, mod_EPFmain, mod_EFIREpreind, \
    mod_LSNKpreind, mod_LSNKnpp, mod_LSNKrho, mod_LSNKtrans, mod_EFIREtrans

##################################################
#   1. CARBON DIOXIDE
##################################################

# ===============
# 1.1. ATMOSPHERE
# ===============

# conversion of CO2 from {ppm} to {GtC}
alpha_CO2 = 0.1765 * np.array([12.0], dtype=dty)

# historic CO2 from IPCC-AR5 {ppm}
# from [IPCC WG1, 2013] annexe 2
CO2_ipcc = np.ones([311 + 1], dtype=dty) * np.nan
TMP = np.array(
    [line for line in csv.reader(
        open("data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1750-2011.CO2.csv", "r"))],
    dtype=dty,
)
CO2_ipcc[50:] = TMP[:, 0]
CO2_0 = np.array([CO2_ipcc[50]], dtype=dty)

# historic CO2 from CMIP5 {ppm}
# from [Meinshausen et al., 2011]
CO2_cmip5 = np.ones([305 + 1], dtype=dty) * np.nan
TMP = np.array(
    [line for line in
     csv.reader(open("data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005.CO2.csv", "r"))],
    dtype=dty
)
CO2_cmip5[65:] = TMP[:, 0]

# historic CO2 from NOAA {ppm}
# from the website
CO2_noaa_ml = np.ones([314 + 1], dtype=dty) * np.nan
TMP = np.array(
    [
        line
        for line in csv.reader(open(
        "data/HistAtmo_NOAA-ESRL/#DATA.HistAtmo_NOAA-ESRL.1959-2014.CO2_maunaloa.csv",
        "r"))],
    dtype=dty,
)
CO2_noaa_ml[259:] = TMP[:, 0]
CO2_noaa_gl = np.ones([314 + 1], dtype=dty) * np.nan
TMP = np.array(
    [
        line
        for line in csv.reader(
        open("data/HistAtmo_NOAA-ESRL/#DATA.HistAtmo_NOAA-ESRL.1980-2014.CO2_global.csv",
             "r"))],
    dtype=dty,
)
CO2_noaa_gl[280:] = TMP[:, 0]

# historic CO2 from Law Dome ice cores {ppm}
# from [Etheridge et al., 1996] and [MacFarling Meure et al., 2006]
CO2_lawdome = np.array(
    [
        line
        for line in csv.reader(open(
        "data/HistAtmo_NOAA-NCDC/#DATA.HistAtmo_NOAA-NCDC.(IceCores).CO2_lawdome.csv",
        "r"))][1:],
    dtype=dty,
)

# load RCP concentrations {ppm}
# from [Meinshausen et al., 2011]
CO2_rcp = np.ones([800 + 1, 6], dtype=dty) * np.nan
n = -1
for rcp in ["rcp26", "rcp45", "rcp60", "rcp85", "rcp45to26", "rcp60to45"]:
    n += 1
    TMP = np.array(
        [line for line in csv.reader(
            open("data/Scenario_ECP/#DATA.Scenario_ECP.2000-2500." + rcp + "_CO2.csv",
                 "r"))],
        dtype=dty,
    )
    CO2_rcp[300:, n] = TMP[:, 0]

# global CO2 historic flux {GtC/yr}
# from [Le Quere et al., 2015]
TMP = np.array(
    [line for line in csv.reader(
        open("data/Historic_GCP/#DATA.Historic_GCP.1959-2014_(5flx).budget.csv", "r"))][
    1:],
    dtype=dty,
)
n = -1
for VAR in ["EFF", "ELUC", "d_CO2", "OSNK", "LSNK"]:
    n += 1
    exec(VAR + "_gcp = np.ones([314+1], dtype=dty) * np.nan")
    exec(VAR + "_gcp[259:] = TMP[:,n]")
OSNK_gcp *= -1
LSNK_gcp *= -1

# ==========
# 1.2. OCEAN
# ==========

# ----------------
# 1.2.1. Structure
# ----------------

# parameters of the surface layer
# from HILDA and other models comparison [Joos et al., 1996]
# gas-exchange coefficient {/yr}, area {m2}, depth {m} and temperature {degC}
if mod_OSNKstruct == "HILDA":
    v_fg = np.array([1 / 9.06], dtype=dty)
    A_ocean = np.array([3.62e14], dtype=dty)
    mld_0 = np.array([75], dtype=dty)
    sst_0 = np.array([18.2], dtype=dty)
elif mod_OSNKstruct == "BD-model":
    v_fg = np.array([1 / 7.80], dtype=dty)
    A_ocean = np.array([3.62e14], dtype=dty)
    mld_0 = np.array([75], dtype=dty)
    sst_0 = np.array([17.7], dtype=dty)
elif mod_OSNKstruct == "2D-model":
    v_fg = np.array([1 / 7.46], dtype=dty)
    A_ocean = np.array([3.54e14], dtype=dty)
    mld_0 = np.array([50], dtype=dty)
    sst_0 = np.array([18.3], dtype=dty)
elif mod_OSNKstruct == "3D-model":
    v_fg = np.array([1 / 7.66], dtype=dty)
    A_ocean = np.array([3.55e14], dtype=dty)
    mld_0 = np.array([50.9], dtype=dty)
    sst_0 = np.array([17.7], dtype=dty)

# sizes and timescales of surface subdivisions for transportation to deep ocean
# adapted from the models comparison by [Joos et al., 1996]
nb_obox = 8
p_circ = np.zeros([nb_obox], dtype=dty)
tau_circ = np.ones([nb_obox], dtype=dty) * np.inf
if mod_OSNKstruct == "HILDA":
    p_circ[1: 1 + 6] = np.array(
        [0.24278, 0.13963, 0.089_318, 0.037_820, 0.035_549, 0.022_936], dtype=dty)
    p_circ[0] = 1 - np.sum(p_circ[1:])
    tau_circ[1: 1 + 6] = np.array([1.2679, 5.2528, 18.601, 68.736, 232.30, np.inf],
                                  dtype=dty)
    tau_circ[0] = 1 / 3.0
elif mod_OSNKstruct == "BD-model":
    p_circ[1: 1 + 7] = np.array(
        [0.16851, 0.11803, 0.076_817, 0.050_469, 0.010_469, 0.031_528, 0.019_737],
        dtype=dty)
    p_circ[0] = 1 - np.sum(p_circ[:])
    tau_circ[1: 1 + 7] = np.array(
        [1.6388, 4.8702, 14.172, 43.506, 148.77, 215.71, np.inf], dtype=dty)
    tau_circ[0] = 1 / 3.0
elif mod_OSNKstruct == "2D-model":
    p_circ[1: 1 + 6] = np.array(
        [0.067_380, 0.036_608, 0.026_994, 0.026_933, 0.012_456, 0.013_691], dtype=dty)
    p_circ[0] = 1 - np.sum(p_circ[:])
    tau_circ[1: 1 + 6] = np.array([10.515, 11.677, 38.946, 107.57, 331.54, np.inf],
                                  dtype=dty)
    tau_circ[0] = 1 / 2.0
elif mod_OSNKstruct == "3D-model":
    p_circ[1: 1 + 6] = np.array(
        [0.70367, 0.24966, 0.066_485, 0.038_344, 0.019_439, 0.014_819], dtype=dty)
    p_circ[1] = 1 - np.sum(p_circ[2:])
    tau_circ[1: 1 + 6] = np.array([0.70177, 2.3488, 15.281, 65.359, 347.55, np.inf],
                                  dtype=dty)

# ----------------
# 1.2.2. Chemistry
# ----------------

# conversion factor for surface ocean{{mumol/kgsol}{m3/ppm}}
# after [Joos et al., 1996]
alpha_sol = np.array([1.722e17], dtype=dty)
alpha_dic = alpha_sol / alpha_CO2 / A_ocean / mld_0


# fonction for converting concentration to mass
def f_dic(D_CSURF, D_mld):
    D_dic = alpha_dic * D_CSURF / (1 + D_mld / mld_0)
    return np.array(D_dic, dtype=dty)


def df_dic(D_CSURF, D_mld):
    df_1 = alpha_dic / (1 + D_mld / mld_0) * D_CSURF
    df_2 = -alpha_dic * D_CSURF / mld_0 / (1 + D_mld / mld_0) ** 2 * D_mld
    df_tot = df_1 + df_2
    return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot)]


# emulation of the carbonate chemistry {ppm}
# after [Lewis and Wallas, 1998] and the CSIRO CTR047
if mod_OSNKchem == "CO2SysPade":

    dic_0 = (
            30015.6
            * (1 - 0.022_653_6 * (sst_0 - 15) + 0.000_167_105 * (sst_0 - 15) ** 2)
            * (CO2_0 / 380)
            / (
                    1
                    + 13.4574 * (1 - 0.019_829 * (sst_0 - 15) + 0.000_113_872 * (
                    sst_0 - 15) ** 2) * (CO2_0 / 380)
                    - 0.243_121 * (1 + 0.000_443_511 * (sst_0 - 15) - 0.000_473_227 * (
                    sst_0 - 15) ** 2) * (CO2_0 / 380) ** 2
            )
    )
    CSURF_0 = p_circ * dic_0 / alpha_dic


    def f_pCO2(D_dic, D_sst):
        a0_0 = 30015.6 * (
                1 - 0.022_653_6 * (sst_0 - 15) + 0.000_167_105 * (sst_0 - 15) ** 2)
        a1_0 = 13.4574 * (
                1 - 0.019_829 * (sst_0 - 15) + 0.000_113_872 * (sst_0 - 15) ** 2)
        a2_0 = -0.243_121 * (
                1 + 0.000_443_511 * (sst_0 - 15) - 0.000_473_227 * (sst_0 - 15) ** 2)
        a0 = 30015.6 * (1 - 0.022_653_6 * (sst_0 + D_sst - 15) + 0.000_167_105 * (
                sst_0 + D_sst - 15) ** 2)
        a1 = 13.4574 * (1 - 0.019_829 * (sst_0 + D_sst - 15) + 0.000_113_872 * (
                sst_0 + D_sst - 15) ** 2)
        a2 = -0.243_121 * (1 + 0.000_443_511 * (sst_0 + D_sst - 15) - 0.000_473_227 * (
                sst_0 + D_sst - 15) ** 2)
        D_pCO2 = 380 * (
                a0 - a1 * (dic_0 + D_dic) - np.sqrt(
            (a0 - a1 * (dic_0 + D_dic)) ** 2 - 4 * a2 * (dic_0 + D_dic) ** 2)
        ) / (2 * a2 * (dic_0 + D_dic)) - 380 * (
                         a0_0 - a1_0 * dic_0 - np.sqrt(
                     (a0_0 - a1_0 * dic_0) ** 2 - 4 * a2_0 * dic_0 ** 2)
                 ) / (
                         2 * a2_0 * dic_0
                 )
        return np.array(D_pCO2, dtype=dty)


    def df_pCO2_ddic(D_dic, D_sst):
        a0 = 30015.6 * (1 - 0.022_653_6 * (sst_0 + D_sst - 15) + 0.000_167_105 * (
                sst_0 + D_sst - 15) ** 2)
        a1 = 13.4574 * (1 - 0.019_829 * (sst_0 + D_sst - 15) + 0.000_113_872 * (
                sst_0 + D_sst - 15) ** 2)
        a2 = -0.243_121 * (1 + 0.000_443_511 * (sst_0 + D_sst - 15) - 0.000_473_227 * (
                sst_0 + D_sst - 15) ** 2)
        rac = -np.sqrt((a0 - a1 * (dic_0 + D_dic)) ** 2 - 4 * a2 * (dic_0 + D_dic) ** 2)
        D_pCO2 = (
                380
                * (
                        (-2 * a0 * a1 * a2 * (dic_0 + D_dic) + (
                                2 * a2 * a1 ** 2 - 8 * a2 ** 2) * (
                                 dic_0 + D_dic) ** 2) / rac
                        - 2 * a0 * a2
                        - 2 * a2 * rac
                )
                / (2 * a2 * (dic_0 + D_dic)) ** 2
        )
        return np.array(D_pCO2, dtype=dty)


    def df_pCO2_dsst(D_dic, D_sst):
        a0 = 30015.6 * (1 - 0.022_653_6 * (sst_0 + D_sst - 15) + 0.000_167_105 * (
                sst_0 + D_sst - 15) ** 2)
        a1 = 13.4574 * (1 - 0.019_829 * (sst_0 + D_sst - 15) + 0.000_113_872 * (
                sst_0 + D_sst - 15) ** 2)
        a2 = -0.243_121 * (1 + 0.000_443_511 * (sst_0 + D_sst - 15) - 0.000_473_227 * (
                sst_0 + D_sst - 15) ** 2)
        da0 = 30015.6 * (-0.022_653_6 + 2 * 0.000_167_105 * (sst_0 + D_sst - 15))
        da1 = 13.4574 * (-0.019_829 + 2 * 0.000_113_872 * (sst_0 + D_sst - 15))
        da2 = -0.243_121 * (0.000_443_511 - 2 * 0.000_473_227 * (sst_0 + D_sst - 15))
        rac = -np.sqrt((a0 - a1 * (dic_0 + D_dic)) ** 2 - 4 * a2 * (dic_0 + D_dic) ** 2)
        D_pCO2 = (
                380
                * (
                        2 * (a2 * da0 - a0 * da2) * (dic_0 + D_dic)
                        + 2 * (a1 * da2 - da1 * a2) * (dic_0 + D_dic) ** 2
                        + (
                                2 * a2 * da0 * a0 * (dic_0 + D_dic)
                                - 2 * a2 * (a0 * da1 + a1 * da0) * (dic_0 + D_dic) ** 2
                                + 2 * a2 * (da1 * a1 - 2 * da2) * (dic_0 + D_dic) ** 3
                        )
                        / rac
                        - 2 * da2 * (dic_0 + D_dic) * rac
                )
                / (2 * a2 * (dic_0 + D_dic)) ** 2
        )
        return np.array(D_pCO2, dtype=dty)


    def df_pCO2(D_dic, D_sst):
        df_1 = df_pCO2_ddic(D_dic, D_sst) * D_dic
        df_2 = df_pCO2_dsst(D_dic, D_sst) * D_sst
        df_tot = df_1 + df_2
        return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot)]


# after [Lewis and Wallas, 1998] and the CSIRO CTR047
elif mod_OSNKchem == "CO2SysPower":

    dic_0 = (
            2160.156
            * (1 - 0.003_450_63 * (sst_0 - 15) - 0.000_025_001_6 * (sst_0 - 15) ** 2)
            * ((CO2_0 / 380) - 0.318_665 * (
            1 - 0.001_512_92 * (sst_0 - 15) - 0.000_198_978 * (sst_0 - 15) ** 2))
            ** (0.059_596_1 * (
            1 + 0.020_032_8 * (sst_0 - 15) + 0.000_192_084 * (sst_0 - 15) ** 2))
    )
    CSURF_0 = p_circ * dic_0 / alpha_dic


    def f_pCO2(D_dic, D_sst):
        p0_0 = 2160.156 * (
                1 - 0.003_450_63 * (sst_0 - 15) - 0.000_025_001_6 * (sst_0 - 15) ** 2)
        p1_0 = 0.059_596_1 * (
                1 + 0.020_032_8 * (sst_0 - 15) + 0.000_192_084 * (sst_0 - 15) ** 2)
        p2_0 = 0.318_665 * (
                1 - 0.001_512_92 * (sst_0 - 15) - 0.000_198_978 * (sst_0 - 15) ** 2)
        p0 = 2160.156 * (1 - 0.003_450_63 * (sst_0 + D_sst - 15) - 0.000_025_001_6 * (
                sst_0 + D_sst - 15) ** 2)
        p1 = 0.059_596_1 * (1 + 0.020_032_8 * (sst_0 + D_sst - 15) + 0.000_192_084 * (
                sst_0 + D_sst - 15) ** 2)
        p2 = 0.318_665 * (1 - 0.001_512_92 * (sst_0 + D_sst - 15) - 0.000_198_978 * (
                sst_0 + D_sst - 15) ** 2)
        D_pCO2 = 380 * (p2 + ((dic_0 + D_dic) / p0) ** (1 / p1)) - 380 * (
                p2_0 + (dic_0 / p0_0) ** (1 / p1_0))
        return np.array(D_pCO2, dtype=dty)


    def df_pCO2_ddic(D_dic, D_sst):
        p0 = 2160.156 * (1 - 0.003_450_63 * (sst_0 + D_sst - 15) - 0.000_025_001_6 * (
                sst_0 + D_sst - 15) ** 2)
        p1 = 0.059_596_1 * (1 + 0.020_032_8 * (sst_0 + D_sst - 15) + 0.000_192_084 * (
                sst_0 + D_sst - 15) ** 2)
        p2 = 0.318_665 * (1 - 0.001_512_92 * (sst_0 + D_sst - 15) - 0.000_198_978 * (
                sst_0 + D_sst - 15) ** 2)
        D_pCO2 = 380 * (1 / p0 / p1) * ((dic_0 + D_dic) / p0) ** (1 / p1 - 1)
        return np.array(D_pCO2, dtype=dty)


    def df_pCO2_dsst(D_dic, D_sst):
        p0 = 2160.156 * (1 - 0.003_450_63 * (sst_0 + D_sst - 15) - 0.000_025_001_6 * (
                sst_0 + D_sst - 15) ** 2)
        p1 = 0.059_596_1 * (1 + 0.020_032_8 * (sst_0 + D_sst - 15) + 0.000_192_084 * (
                sst_0 + D_sst - 15) ** 2)
        p2 = 0.318_665 * (1 - 0.001_512_92 * (sst_0 + D_sst - 15) - 0.000_198_978 * (
                sst_0 + D_sst - 15) ** 2)
        dp0 = 2160.156 * (-0.003_450_63 - 2 * 0.000_025_001_6 * (sst_0 + D_sst - 15))
        dp1 = 0.059_596_1 * (0.020_032_8 + 2 * 0.000_192_084 * (sst_0 + D_sst - 15))
        dp2 = 0.318_665 * (-0.001_512_92 - 2 * 0.000_198_978 * (sst_0 + D_sst - 15))
        D_pCO2 = 380 * (
                dp2
                + ((-dp1 / p1 / p1) * np.log((dic_0 + D_dic) / p0) - (dp0 / p0 / p1)) * (
                        (dic_0 + D_dic) / p0) ** (1 / p1)
        )
        return np.array(D_pCO2, dtype=dty)


    def df_pCO2(D_dic, D_sst):
        df_1 = df_pCO2_ddic(D_dic, D_sst) * D_dic
        df_2 = df_pCO2_dsst(D_dic, D_sst) * D_sst
        df_tot = df_1 + df_2
        return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot)]

# ------------
# 1.2.3. CMIP5
# ------------

# period of CMIP5 simulations and depth of MLD per model
prd = {"ctrl": "451yr", "hist": "1850-2005", "rcp85": "2006-2300"}
lng = {"ctrl": 451, "hist": 156, "rcp85": 295}
mld = {"HILDA": "75m", "BD-model": "75m", "2D-model": "50m", "3D-model": "50m9"}

# load pre-processed CMIP5 results for specified model
if mod_OSNKtrans != "":
    for sim in ["ctrl", "hist", "rcp85"]:
        for VAR in ["FCO2", "FEXP"]:
            if os.path.isfile(
                    "data/Ocean_CMIP5/#DATA.Ocean_"
                    + mod_OSNKtrans
                    + "."
                    + prd[sim]
                    + "_18x10lat."
                    + sim
                    + "_"
                    + VAR
                    + ".csv"
            ):
                TMP = np.array(
                    [
                        line
                        for line in csv.reader(
                        open(
                            "data/Ocean_CMIP5/#DATA.Ocean_"
                            + mod_OSNKtrans
                            + "."
                            + prd[sim]
                            + "_18x10lat."
                            + sim
                            + "_"
                            + VAR
                            + ".csv",
                            "r",
                        )
                    )],
                    dtype=dty,
                )
                exec(VAR + "_" + sim + " = np.sum(TMP,1)")
            else:
                exec(VAR + "_" + sim + " = np.zeros([lng[sim]], dtype=dty)")
        for var in ["tos", "dpCO2", "sic", "mld"]:
            if os.path.isfile(
                    "data/Ocean_CMIP5/#DATA.Ocean_"
                    + mod_OSNKtrans
                    + "."
                    + prd[sim]
                    + "_18x10lat."
                    + sim
                    + "_"
                    + var
                    + ".csv"
            ):
                TMP = np.array(
                    [
                        line
                        for line in csv.reader(
                        open(
                            "data/Ocean_CMIP5/#DATA.Ocean_"
                            + mod_OSNKtrans
                            + "."
                            + prd[sim]
                            + "_18x10lat."
                            + sim
                            + "_"
                            + var
                            + ".csv",
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
                            "data/Ocean_CMIP5/#DATA.Ocean_" + mod_OSNKtrans + ".451yr_18x10lat.SURF.csv",
                            "r")
                    )],
                    dtype=dty,
                )
                exec(
                    var + "_" + sim + " = np.sum(TMP*TMP2[:len(TMP)],1)/np.sum(TMP2[:len(TMP)],1)")
            else:
                exec(var + "_" + sim + " = np.zeros([lng[sim]], dtype=dty)")
        for var in ["amoc"]:
            if os.path.isfile(
                    "data/Ocean_CMIP5/#DATA.Ocean_" + mod_OSNKtrans + "." + prd[
                        sim] + "." + sim + "_" + var + ".csv"
            ):
                TMP = np.array(
                    [
                        line
                        for line in csv.reader(
                        open(
                            "data/Ocean_CMIP5/#DATA.Ocean_"
                            + mod_OSNKtrans
                            + "."
                            + prd[sim]
                            + "."
                            + sim
                            + "_"
                            + var
                            + ".csv",
                            "r",
                        )
                    )],
                    dtype=dty,
                )
                exec(var + "_" + sim + " = TMP[:,0]")
            else:
                exec(var + "_" + sim + " = np.zeros([lng[sim]], dtype=dty)")
    # aggregate all experiments
    if mod_OSNKtrans != "":
        for VAR in ["FCO2", "FEXP"] + ["tos", "dpCO2", "sic", "mld"] + ["amoc"]:
            exec(VAR + "_all = np.array(list(" + VAR + "_hist)+list(" + VAR + "_rcp85))")

# definition of parameters
# for sea ice concentration reduction {pp}&{./K}
Alpha_sic = np.array([0], dtype=dty)
gamma_sic = np.array([0], dtype=dty)
# for sea surface exchange rate {.}&{./K}
gamma_fg = np.array([0], dtype=dty)
# for mixing layer stratification {.}&{./K}
alpha_mld = np.array([0], dtype=dty)
gamma_mld = np.array([0], dtype=dty)
# for overturning circulation slowdown {.}&{./K}
alpha_amoc = np.array([0], dtype=dty)
gamma_amoc = np.array([0], dtype=dty)
# for biological pump reduction decrease {GtC/yr}&{./K}
Alpha_exp = np.array([0], dtype=dty)
gamma_exp = np.array([0], dtype=dty)

# fit of parameters
# sic
if mod_OSNKtrans != "":
    obj = sic_all


    def err(var):
        carb = 1
        clim = np.mean(sic_ctrl) * np.exp(
            1 - np.exp(var[0] * (tos_all - np.mean(tos_ctrl))))
        return np.sum((obj - carb * clim) ** 2)


    [gamma_sic[0]] = fmin(err, [0], disp=False)
    Alpha_sic[0] = np.mean(sic_ctrl)
# fg
if mod_OSNKtrans != "":
    diff = FCO2_all - np.mean(FCO2_ctrl)


    def err(var):
        carb = 1
        clim = var[0] * (1 - sic_all) * (dpCO2_all - np.mean(dpCO2_ctrl)) * (
                1 + var[1] * (tos_all - np.mean(tos_ctrl)))
        return np.sum((diff - carb * clim) ** 2)


    [TMP, gamma_fg[0]] = fmin(err, [-1, 0], disp=False)
# mld
if mod_OSNKtrans != "":
    ratio = mld_all / np.mean(mld_ctrl)


    def err(var):
        carb = 1
        clim = (1 - var[0]) + var[0] * np.exp(var[1] * (tos_all - np.mean(tos_ctrl)))
        return np.sum((ratio - carb * clim) ** 2)


    [alpha_mld[0], gamma_mld[0]] = fmin(err, [0.5, 0], disp=False)
# amoc
if mod_OSNKtrans != "":
    ratio = amoc_all / np.mean(amoc_ctrl)


    def err(var):
        carb = 1
        clim = (1 - var[0]) + var[0] * np.exp(
            1 - np.exp(var[1] * (tos_all - np.mean(tos_ctrl))))
        return np.sum((ratio - carb * clim) ** 2)


    [alpha_amoc[0], gamma_amoc[0]] = fmin(err, [0.5, 0], disp=False)
# exp
if mod_OSNKtrans != "":
    diff = FEXP_all - np.mean(FEXP_ctrl)


    def err(var):
        carb = 1
        clim = var[0] * (1 - np.exp(var[1] * (tos_all - np.mean(tos_ctrl))))
        return np.sum((diff - carb * clim) ** 2)


    [Alpha_exp[0], gamma_exp[0]] = fmin(err, [-1, 0], disp=False)

# =========
# 1.3. LAND
# =========

# ----------------
# 1.3.1. Functions
# ----------------

# reference beta-factor of NPP {.}
# from FACE experiments [Norby et al., 2005]
beta_ref = np.array([0.60], dtype=dty)

# reference Q10 of HR {.}
# from [Mahecha et al., 2010]
Q10_ref = np.array([1.4], dtype=dty)

# expression of NPP function
# logarithmic formulation [Friedlingstein et al., 1995]
if mod_LSNKnpp == "log":

    def f_npp(D_CO2, D_lst, D_lyp):
        D_npp = (1 + beta_npp * np.log(1 + D_CO2 / CO2_0)) * (
                1 + gamma_nppT * D_lst + gamma_nppP * D_lyp) - 1
        return np.array(D_npp, dtype=dty)


    def df_npp_dCO2(D_CO2, D_lst, D_lyp):
        D_npp = beta_npp * (1 + gamma_nppT * D_lst + gamma_nppP * D_lyp) / (CO2_0 + D_CO2)
        return np.array(D_npp, dtype=dty)


    def df_npp_dlst(D_CO2, D_lst, D_lyp):
        D_npp = gamma_nppT * (1 + beta_npp * np.log(1 + D_CO2 / CO2_0))
        return np.array(D_npp, dtype=dty)


    def df_npp_dlyp(D_CO2, D_lst, D_lyp):
        D_npp = gamma_nppP * (1 + beta_npp * np.log(1 + D_CO2 / CO2_0))
        return np.array(D_npp, dtype=dty)


    def df_npp(D_CO2, D_lst, D_lyp):
        df_1 = df_npp_dCO2(D_CO2, D_lst, D_lyp) * D_CO2
        df_2 = df_npp_dlst(D_CO2, D_lst, D_lyp) * D_lst
        df_3 = df_npp_dlyp(D_CO2, D_lst, D_lyp) * D_lyp
        df_tot = df_1 + df_2 + df_3
        return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot),
                np.nan_to_num(df_3 / df_tot)]


# hyperbolic formulation [Friedlingstein et al., 1995]
elif mod_LSNKnpp == "hyp":

    def f_npp(D_CO2, D_lst, D_lyp):
        K2 = ((2 * CO2_0 - CO2_comp) - beta_npp0 * (CO2_0 - CO2_comp)) / (
                (beta_npp0 - 1) * (CO2_0 - CO2_comp) * (2 * CO2_0 - CO2_comp)
        )
        K1 = (1 + K2 * (CO2_0 - CO2_comp)) / (CO2_0 - CO2_comp)
        D_npp = (
                K1
                * (CO2_0 + D_CO2 - CO2_comp)
                / (1 + K2 * (CO2_0 + D_CO2 - CO2_comp))
                * (1 + gamma_nppT * D_lst + gamma_nppP * D_lyp)
                - 1
        )
        return np.array(D_npp, dtype=dty)


    def df_npp_dCO2(D_CO2, D_lst, D_lyp):
        K2 = ((2 * CO2_0 - CO2_comp) - beta_npp0 * (CO2_0 - CO2_comp)) / (
                (beta_npp0 - 1) * (CO2_0 - CO2_comp) * (2 * CO2_0 - CO2_comp)
        )
        K1 = (1 + K2 * (CO2_0 - CO2_comp)) / (CO2_0 - CO2_comp)
        D_npp = (1 + gamma_nppT * D_lst + gamma_nppP * D_lyp) * K1 / (
                1 + K2 * (CO2_0 + D_CO2 - CO2_comp)) ** 2
        return np.array(D_npp, dtype=dty)


    def df_npp_dlst(D_CO2, D_lst, D_lyp):
        K2 = ((2 * CO2_0 - CO2_comp) - beta_npp0 * (CO2_0 - CO2_comp)) / (
                (beta_npp0 - 1) * (CO2_0 - CO2_comp) * (2 * CO2_0 - CO2_comp)
        )
        K1 = (1 + K2 * (CO2_0 - CO2_comp)) / (CO2_0 - CO2_comp)
        D_npp = gamma_nppT * K1 * (CO2_0 + D_CO2 - CO2_comp) / (
                1 + K2 * (CO2_0 + D_CO2 - CO2_comp))
        return np.array(D_npp, dtype=dty)


    def df_npp_dlyp(D_CO2, D_lst, D_lyp):
        K2 = ((2 * CO2_0 - CO2_comp) - beta_npp0 * (CO2_0 - CO2_comp)) / (
                (beta_npp0 - 1) * (CO2_0 - CO2_comp) * (2 * CO2_0 - CO2_comp)
        )
        K1 = (1 + K2 * (CO2_0 - CO2_comp)) / (CO2_0 - CO2_comp)
        D_npp = gamma_nppP * K1 * (CO2_0 + D_CO2 - CO2_comp) / (
                1 + K2 * (CO2_0 + D_CO2 - CO2_comp))
        return np.array(D_npp, dtype=dty)


    def df_npp(D_CO2, D_lst, D_lyp):
        df_1 = df_npp_dCO2(D_CO2, D_lst, D_lyp) * D_CO2
        df_2 = df_npp_dlst(D_CO2, D_lst, D_lyp) * D_lst
        df_3 = df_npp_dlyp(D_CO2, D_lst, D_lyp) * D_lyp
        df_tot = df_1 + df_2 + df_3
        return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot),
                np.nan_to_num(df_3 / df_tot)]

# expression of RH function
# exponential formulation [Tuomi et al., 2008]
if mod_LSNKrho == "exp":

    def f_rho(D_lst, D_lyp):
        D_rho = np.exp(gamma_rhoT * D_lst + gamma_rhoP * D_lyp) - 1
        return np.array(D_rho, dtype=dty)


    def df_rho_dlst(D_lst, D_lyp):
        D_rho = gamma_rhoT * np.exp(gamma_rhoT * D_lst + gamma_rhoP * D_lyp)
        return np.array(D_rho, dtype=dty)


    def df_rho_dlyp(D_lst, D_lyp):
        D_rho = gamma_rhoP * np.exp(gamma_rhoT * D_lst + gamma_rhoP * D_lyp)
        return np.array(D_rho, dtype=dty)


    def df_rho(D_lst, D_lyp):
        df_1 = df_rho_dlst(D_lst, D_lyp) * D_lst
        df_2 = df_rho_dlyp(D_lst, D_lyp) * D_lyp
        df_tot = df_1 + df_2
        return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot)]


# gaussian formulation [Tuomi et al., 2008]
elif mod_LSNKrho == "gauss":

    def f_rho(D_lst, D_lyp):
        D_rho = np.exp(
            gamma_rhoT1 * D_lst + gamma_rhoT2 * D_lst ** 2 + gamma_rhoP * D_lyp) - 1
        return np.array(D_rho, dtype=dty)


    def df_rho_dlst(D_lst, D_lyp):
        D_rho = (gamma_rhoT1 + 2 * gamma_rhoT2 * D_lst) * np.exp(
            gamma_rhoT1 * D_lst + gamma_rhoT2 * D_lst ** 2 + gamma_rhoP * D_lyp
        )
        return np.array(D_rho, dtype=dty)


    def df_rho_dlyp(D_lst, D_lyp):
        D_rho = gamma_rhoP * np.exp(
            gamma_rhoT1 * D_lst + gamma_rhoT2 * D_lst ** 2 + gamma_rhoP * D_lyp)
        return np.array(D_rho, dtype=dty)


    def df_rho(D_lst, D_lyp):
        df_1 = df_rho_dlst(D_lst, D_lyp) * D_lst
        df_2 = df_rho_dlyp(D_lst, D_lyp) * D_lyp
        df_tot = df_1 + df_2
        return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot)]

# -------------
# 1.3.2. TRENDY
# -------------

# basic biomes of aggregation
bio = ["des", "for", "shr", "gra", "cro", "pas", "urb"]

# load pre-processed TRENDY results for specified model
# data related to preindustrial C-cycle
for VAR in ["AREA", "CSOIL", "CLITTER", "CVEG", "NPP", "RH"]:
    exec(VAR + "_pi = np.zeros([nb_regionI,len(bio)], dtype=dty)")
    if os.path.isfile(
            "data/Land_TRENDYv2/#DATA.Land_" + mod_LSNKpreind + ".1910s_114reg1_7bio." + VAR + ".csv"):
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/Land_TRENDYv2/#DATA.Land_" + mod_LSNKpreind + ".1910s_114reg1_7bio." + VAR + ".csv",
                    "r")
            )],
            dtype=dty,
        )
        for i in range(1, 114 + 1):
            exec(VAR + "_pi[regionI_index[i],:] += TMP[i-1,:]")

# complete with arbitrary values
# if no litter given by model
if np.sum(CLITTER_pi) == 0:
    CLITTER_pi = 0.05 * CSOIL_pi
    CSOIL_pi = 0.95 * CSOIL_pi
# if no shrubs in the model
if (nb_biome > 1) & (np.sum(AREA_pi[:, bio.index("shr")]) == 0) & (mod_biomeSHR == "SHR"):
    AREA_pi[:, bio.index("shr")] = AREA_pi[:, bio.index("gra")]
    for VAR in ["CSOIL", "CLITTER", "CVEG", "NPP", "RH"]:
        exec(
            VAR
            + '_pi[:,bio.index("shr")] = 0.85*'
            + VAR
            + '_pi[:,bio.index("gra")] + 0.15*'
            + VAR
            + '_pi[:,bio.index("for")]*AREA_pi[:,bio.index("gra")]/AREA_pi[:,bio.index("for")]'
        )
    TMP = AREA_pi[:, :, bio.index("for")] == 0
    exec(
        VAR + '_pi[:,:,bio.index("shr")][TMP] = ' + VAR + '_pi[:,:,bio.index("gra")][TMP]')
# if no crops in the model
if (nb_biome > 1) & (np.sum(AREA_pi[:, bio.index("cro")]) == 0):
    for VAR in ["AREA", "CSOIL", "CLITTER", "CVEG", "NPP", "RH"]:
        exec(VAR + '_pi[:,bio.index("cro")] = ' + VAR + '_pi[:,bio.index("gra")]')
# if no pastures in the model
if (nb_biome > 1) & (np.sum(AREA_pi[:, bio.index("pas")]) == 0):
    AREA_pi[:, bio.index("pas")] = AREA_pi[:, bio.index("gra")]
    for VAR in ["CSOIL", "CLITTER", "CVEG", "NPP", "RH"]:
        exec(
            VAR
            + '_pi[:,bio.index("pas")] = 0.60*'
            + VAR
            + '_pi[:,bio.index("gra")] + 0.40*'
            + VAR
            + '_pi[:,bio.index("des")]*AREA_pi[:,bio.index("gra")]/AREA_pi[:,bio.index("des")]'
        )
        TMP = AREA_pi[:, bio.index("des")] == 0
        exec(
            VAR + '_pi[:,bio.index("pas")][TMP] = ' + VAR + '_pi[:,bio.index("gra")][TMP]')
# if no urban in the model
if (nb_biome > 1) & (np.sum(AREA_pi[:, bio.index("urb")]) == 0) & (
        (mod_biomeURB == "URB") | mod_biomeV3):
    for VAR in ["AREA", "CSOIL", "CLITTER", "CVEG", "NPP", "RH"]:
        exec(VAR + '_pi[:,bio.index("urb")] = ' + VAR + '_pi[:,bio.index("des")]')

# data related to preindustrial fires
for VAR in ["EFIRE", "CBURNT"]:
    exec(VAR + "_pi = np.zeros([nb_regionI,len(bio)], dtype=dty)")
    if mod_EFIREpreind != "":
        if os.path.isfile(
                "data/BioBurn_TRENDYv2/#DATA.BioBurn_" + mod_EFIREpreind + ".1910s_114reg1_7bio." + VAR + ".csv"
        ):
            TMP = np.array(
                [
                    line
                    for line in csv.reader(
                    open(
                        "data/BioBurn_TRENDYv2/#DATA.BioBurn_"
                        + mod_EFIREpreind
                        + ".1910s_114reg1_7bio."
                        + VAR
                        + ".csv",
                        "r",
                    )
                )],
                dtype=dty,
            )
        for i in range(1, 114 + 1):
            exec(VAR + "_pi[regionI_index[i],:] += TMP[i-1,:]")

# complete with arbitrary values
# if no shrubs in the model
if (nb_biome > 1) & (np.sum(CBURNT_pi[:, bio.index("shr")]) == 0) & (
        mod_biomeSHR == "SHR"):
    CBURNT_pi[:, bio.index("shr")] = CBURNT_pi[:, bio.index("gra")]
    for VAR in ["EFIRE"]:
        exec(
            VAR
            + '_pi[:,bio.index("shr")] = 0.85*'
            + VAR
            + '_pi[:,bio.index("gra")] + 0.15*'
            + VAR
            + '_pi[:,bio.index("for")]*CBURNT_pi[:,bio.index("gra")]/CBURNT_pi[:,bio.index("for")]'
        )
        TMP = CBURNT_pi[:, :, bio.index("for")] == 0
        exec(
            VAR + '_pi[:,:,bio.index("shr")][TMP] = ' + VAR + '_pi[:,:,bio.index("gra")][TMP]')
# if no crops in the model
if (nb_biome > 1) & (np.sum(CBURNT_pi[:, bio.index("cro")]) == 0):
    for VAR in ["EFIRE", "CBURNT"]:
        exec(VAR + '_pi[:,bio.index("cro")] = ' + VAR + '_pi[:,bio.index("gra")]')
# if no pastures in the model
if (nb_biome > 1) & (np.sum(CBURNT_pi[:, bio.index("pas")]) == 0):
    CBURNT_pi[:, bio.index("pas")] = CBURNT_pi[:, bio.index("gra")]
    for VAR in ["EFIRE"]:
        exec(
            VAR
            + '_pi[:,bio.index("pas")] = 0.60*'
            + VAR
            + '_pi[:,bio.index("gra")] + 0.40*'
            + VAR
            + '_pi[:,bio.index("des")]*CBURNT_pi[:,bio.index("gra")]/CBURNT_pi[:,bio.index("des")]'
        )
        TMP = CBURNT_pi[:, bio.index("des")] == 0
        exec(
            VAR + '_pi[:,bio.index("pas")][TMP] = ' + VAR + '_pi[:,bio.index("gra")][TMP]')
# if no urban in the model
if (nb_biome > 1) & (np.sum(CBURNT_pi[:, bio.index("urb")]) == 0) & (
        (mod_biomeURB == "URB") | mod_biomeV3):
    for VAR in ["EFIRE", "CBURNT"]:
        exec(VAR + '_pi[:,bio.index("urb")] = ' + VAR + '_pi[:,bio.index("des")]')

# aggregate biomes
for VAR in ["AREA", "CSOIL", "CLITTER", "CVEG", "NPP", "RH"] + ["EFIRE", "CBURNT"]:
    exec("TMP = " + VAR + "_pi.copy()")
    exec(VAR + "_pi = np.zeros([nb_regionI,nb_biome], dtype=dty)")
    for b in range(len(bio)):
        exec(VAR + "_pi[:,biome_index[bio[b]]] += TMP[:,b]")

# ratio for metabolization of litter to soil {.}
# from [Foley, 1995]
k_met = np.array([0.3 / 0.7], dtype=dty)

# calculate preindustrial flux parameters
# net primary productivity {GtC/Mha/yr}
npp_0 = NPP_pi / AREA_pi
# rate of biomass mortality {./yr}
mu_0 = RH_pi / CVEG_pi
# rates of heterotrophic respiration {./yr}
rho_0 = RH_pi / (CLITTER_pi + CSOIL_pi)
rho1_0 = 1 / (1 + k_met) * RH_pi / CLITTER_pi
rho2_0 = k_met / (1 + k_met) * RH_pi / CSOIL_pi
# rate of fire ignition {./yr}
igni_0 = EFIRE_pi / CBURNT_pi
# [NaN]
for var in ["npp", "mu", "rho", "rho1", "rho2", "igni"]:
    exec(var + "_0[np.isnan(" + var + "_0)|np.isinf(" + var + "_0)] = 0")

# arbitrary adjustment of parameters
# lower preindustrial npp
npp_0 *= CO2_0 / np.mean(CO2_cmip5[201:231])
# crops: yearly harvest of 80% & protection against wildfires of 100%
if nb_biome > 1:
    npp_0[:, biome_index["cro"]] *= 0.2
    mu_0[:, biome_index["cro"]] = 1.0
    igni_0[:, biome_index["cro"]] *= 0.0
# pasts: harvest of 0% & protection against wildfires of 0%
if nb_biome > 1:
    npp_0[:, biome_index["pas"]] *= 1.0
    igni_0[:, biome_index["pas"]] *= 1.0
# urban: protection against wildfires of 100%
if (nb_biome > 1) & ((mod_biomeURB == "URB") | mod_biomeV3):
    igni_0[:, biome_index["urb"]] = 0.0

# calculate preindustrial stocks {GtC/Mha}
cveg_0 = npp_0 / (mu_0 + igni_0)
csoil_0 = mu_0 / rho_0 * cveg_0
csoil1_0 = 1 / (1 + k_met) * mu_0 / rho1_0 * cveg_0
csoil2_0 = k_met * (rho1_0 / rho2_0) * csoil1_0
# [NaN]
for var in ["cveg", "csoil", "csoil1", "csoil2"]:
    exec(var + "_0[np.isnan(" + var + "_0)|np.isinf(" + var + "_0)] = 0")

# -----------------
# 1.3.3. Land-Cover
# -----------------

# load preindustrial land-cover {Mha}
AREA_0 = np.zeros([nb_regionI, nb_biome], dtype=dty)
TMP = np.array(
    [
        line
        for line in csv.reader(
        open(
            "data/LandCover_"
            + data_LULCC
            + "/#DATA.LandCover_"
            + data_LULCC
            + "_"
            + mod_LSNKcover
            + ".1700_114reg1_7bio.AREA.csv",
            "r",
        )
    )],
    dtype=dty,
)
for i in range(1, 114 + 1):
    for b in range(len(bio)):
        AREA_0[regionI_index[i], biome_index[bio[b]]] += TMP[i - 1, b]

# force true 1750 preindustrial (land-cover)
if PI_1750:
    LUC_tmp = np.zeros([50 + 1, nb_regionI, nb_biome, nb_biome], dtype=dty)
    bio = ["des", "for", "shr", "gra", "cro", "pas", "urb"]
    # load LUH data
    if data_LULCC[:3] == "LUH":
        for b1 in range(len(bio)):
            for b2 in range(len(bio)):
                if os.path.isfile(
                        "data/LandUse_"
                        + data_LULCC
                        + "/#DATA.LandUse_"
                        + data_LULCC
                        + "_"
                        + mod_LSNKcover
                        + ".1501-2015_114reg1.LUC_"
                        + bio[b1]
                        + "2"
                        + bio[b2]
                        + ".csv"
                ):
                    TMP = np.array(
                        [
                            line
                            for line in csv.reader(
                            open(
                                "data/LandUse_LUH1/#DATA.LandUse_"
                                + data_LULCC
                                + "_"
                                + mod_LSNKcover
                                + ".1501-2015_114reg1.LUC_"
                                + bio[b1]
                                + "2"
                                + bio[b2]
                                + ".csv",
                                "r",
                            )
                        )],
                        dtype=dty,
                    )
                    for i in range(1, 114 + 1):
                        LUC_tmp[1: 50 + 1, regionI_index[i], biome_index[bio[b1]],
                        biome_index[bio[b2]]] += TMP[
                                                 200: 200 + 50, i - 1]
    # correct land-cover
    AREA_0 += np.sum(np.sum(LUC_tmp, 2), 0) - np.sum(np.sum(LUC_tmp, 3), 0)
    del LUC_tmp

# ------------
# 1.3.4. CMIP5
# ------------

# basic biomes of aggregation
bio = ["des", "for", "shr", "gra", "cro", "pas", "urb"]

# load pre-processed CMIP5 results for specified model
# data related to C-cycle
for sim in ["ctrl", "upct", "fxcl", "fdbk"]:
    for VAR in ["AREA", "CSOIL0", "NPP", "RH", "FINPUT"]:
        exec(VAR + "_" + sim + " = np.zeros([140,nb_regionI,len(bio)], dtype=dty)")
        for b in range(len(bio)):
            if os.path.isfile(
                    "data/Land_CMIP5/#DATA.Land_"
                    + mod_LSNKtrans
                    + ".140yr_114reg1."
                    + sim
                    + "_"
                    + VAR
                    + "_"
                    + bio[b]
                    + ".csv"
            ):
                TMP = np.array(
                    [
                        line
                        for line in csv.reader(
                        open(
                            "data/Land_CMIP5/#DATA.Land_"
                            + mod_LSNKtrans
                            + ".140yr_114reg1."
                            + sim
                            + "_"
                            + VAR
                            + "_"
                            + bio[b]
                            + ".csv",
                            "r",
                        )
                    )],
                    dtype=dty,
                )
                for i in range(1, 114 + 1):
                    exec(VAR + "_" + sim + "[:,regionI_index[i],b] += TMP[:,i-1]")
    for var in ["tas", "pr"]:
        exec(var + "_" + sim + " = np.zeros([140,nb_regionI], dtype=dty)")
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/Land_CMIP5/#DATA.Land_" + mod_LSNKtrans + ".140yr_114reg1." + sim + "_" + var + ".csv",
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
                    "data/Land_CMIP5/#DATA.Land_" + mod_LSNKtrans + ".140yr_114reg1.AREA.csv",
                    "r")
            )],
            dtype=dty,
        )
        for i in range(1, 114 + 1):
            exec(var + "_" + sim + "[:,regionI_index[i]] += TMP[:,i-1]*TMP2[:,i-1]")
        TMP = np.zeros([140, nb_regionI], dtype=dty)
        for i in range(1, 114 + 1):
            exec("TMP[:,regionI_index[i]] += TMP2[:,i-1]")
        exec(var + "_" + sim + " /= TMP")
        exec(
            var + "_" + sim + "[np.isnan(" + var + "_" + sim + ")|np.isinf(" + var + "_" + sim + ")] = 0")

# complete with arbitrary values
# if no shrubs in the model
if (nb_biome > 1) & (np.sum(AREA_ctrl[:, :, bio.index("shr")]) == 0) & (
        mod_biomeSHR == "SHR"):
    for sim in ["ctrl", "upct", "fxcl", "fdbk"]:
        exec(
            "AREA_" + sim + '[:,:,bio.index("shr")] = AREA_' + sim + '[:,:,bio.index("gra")]')
        for VAR in ["CSOIL0", "NPP", "RH", "FINPUT"]:
            exec(
                VAR
                + "_"
                + sim
                + '[:,:,bio.index("shr")] = 0.85*'
                + VAR
                + "_"
                + sim
                + '[:,:,bio.index("gra")] + 0.15*'
                + VAR
                + "_"
                + sim
                + '[:,:,bio.index("for")]*AREA_'
                + sim
                + '[:,:,bio.index("gra")]/AREA_'
                + sim
                + '[:,:,bio.index("for")]'
            )
            exec("TMP = (AREA_" + sim + '[:,:,bio.index("for")] == 0)')
            exec(
                VAR + "_" + sim + '[:,:,bio.index("shr")][TMP] = ' + VAR + "_" + sim + '[:,:,bio.index("gra")][TMP]')
# if no crops in the model
if (nb_biome > 1) & (np.sum(AREA_ctrl[:, :, bio.index("cro")]) == 0):
    for sim in ["ctrl", "upct", "fxcl", "fdbk"]:
        for VAR in ["AREA", "CSOIL0", "NPP", "RH", "FINPUT"]:
            exec(
                VAR + "_" + sim + '[:,:,bio.index("cro")] = ' + VAR + "_" + sim + '[:,:,bio.index("gra")]')
# if no pastures in the model
if (nb_biome > 1) & (np.sum(AREA_ctrl[:, :, bio.index("pas")]) == 0):
    for sim in ["ctrl", "upct", "fxcl", "fdbk"]:
        exec(
            "AREA_" + sim + '[:,:,bio.index("pas")] = AREA_' + sim + '[:,:,bio.index("gra")]')
        for VAR in ["CSOIL0", "NPP", "RH", "FINPUT"]:
            exec(
                VAR
                + "_"
                + sim
                + '[:,:,bio.index("pas")] = 0.60*'
                + VAR
                + "_"
                + sim
                + '[:,:,bio.index("gra")] + 0.40*'
                + VAR
                + "_"
                + sim
                + '[:,:,bio.index("des")]*AREA_'
                + sim
                + '[:,:,bio.index("gra")]/AREA_'
                + sim
                + '[:,:,bio.index("des")]'
            )
            exec("TMP = (AREA_" + sim + '[:,:,bio.index("des")] == 0)')
            exec(
                VAR + "_" + sim + '[:,:,bio.index("pas")][TMP] = ' + VAR + "_" + sim + '[:,:,bio.index("gra")][TMP]')
# if no urban in the model
if (nb_biome > 1) & (np.sum(AREA_ctrl[:, :, bio.index("urb")]) == 0) & (
        (mod_biomeURB == "URB") | mod_biomeV3):
    for sim in ["ctrl", "upct", "fxcl", "fdbk"]:
        for VAR in ["AREA", "CSOIL0", "NPP", "RH", "FINPUT"]:
            exec(
                VAR + "_" + sim + '[:,:,bio.index("urb")] = ' + VAR + "_" + sim + '[:,:,bio.index("des")]')

# aggregate biomes
for sim in ["ctrl", "upct", "fxcl", "fdbk"]:
    for VAR in ["AREA", "CSOIL0", "NPP", "RH", "FINPUT"]:
        exec("TMP = " + VAR + "_" + sim + ".copy()")
        exec(VAR + "_" + sim + " = np.zeros([140,nb_regionI,nb_biome], dtype=dty)")
        for b in range(len(bio)):
            exec(VAR + "_" + sim + "[:,:,biome_index[bio[b]]] += TMP[:,:,b]")
# aggregate all experiments
CO2_all = np.array(
    list(CO2_cmip5[150] * 1.01 ** np.arange(140))
    + list(CO2_cmip5[150] * 1.01 ** np.arange(140))
    + list(CO2_cmip5[150] * np.ones(140)),
    dtype=dty,
)
for VAR in ["AREA", "CSOIL0", "NPP", "RH", "FINPUT"] + ["tas", "pr"]:
    exec(
        VAR + "_all = np.array(list(" + VAR + "_upct)+list(" + VAR + "_fxcl)+list(" + VAR + "_fdbk), dtype=dty)")
# decadal means
for VAR in ["AREA", "CSOIL0", "NPP", "RH", "FINPUT"] + ["tas", "pr"] + ["CO2"]:
    exec(
        VAR + "_dec= np.zeros([3*(140-10)]+list(np.shape(" + VAR + "_all)[1:]), dtype=dty)")
    for t in range(130):
        exec(VAR + "_dec[t,...] = np.mean(" + VAR + "_all[t:t+10,...],0)")
        exec(VAR + "_dec[t+140-10,...] = np.mean(" + VAR + "_all[140+t:140+t+10,...],0)")
        exec(
            VAR + "_dec[t+2*140-2*10,...] = np.mean(" + VAR + "_all[2*140+t:2*140+t+10,...],0)")

# definition of parameters
# sensitivity of NPP to CO2 {.} or {.}&{ppm}
if mod_LSNKnpp == "log":
    beta_npp = np.zeros([nb_regionI, nb_biome], dtype=dty)
elif mod_LSNKnpp == "hyp":
    beta_npp0 = np.ones([nb_regionI, nb_biome], dtype=dty) * 1.000_001
    CO2_comp = np.ones([nb_regionI, nb_biome], dtype=dty) * 1e18
# sensitivity of NPP to climate {/K}&{/mm}
gamma_nppT = np.zeros([nb_regionI, nb_biome], dtype=dty)
gamma_nppP = np.zeros([nb_regionI, nb_biome], dtype=dty)
# sensitivity of RH rate to soil input {.} (not used)
prim_rho = np.zeros([nb_regionI, nb_biome], dtype=dty)
# sensitivity of RH rate to climate {/K}&{/mm} or {/K}&{/K**2}&{/mm}
if mod_LSNKrho == "exp":
    gamma_rhoT = np.zeros([nb_regionI, nb_biome], dtype=dty)
    gamma_rhoP = np.zeros([nb_regionI, nb_biome], dtype=dty)
elif mod_LSNKrho == "gauss":
    gamma_rhoT1 = np.zeros([nb_regionI, nb_biome], dtype=dty)
    gamma_rhoT2 = np.zeros([nb_regionI, nb_biome], dtype=dty)
    gamma_rhoP = np.zeros([nb_regionI, nb_biome], dtype=dty)

# fit of parameters
for i in range(1, nb_regionI):
    for b in range(nb_biome):

        # for NPP
        ratio = (NPP_all / AREA_all)[:, i, b] / np.mean((NPP_ctrl / AREA_ctrl)[:, i, b])
        ratio_dec = (NPP_dec / AREA_dec)[:, i, b] / np.mean(
            (NPP_ctrl / AREA_ctrl)[:, i, b])
        if mod_LSNKnpp == "log":

            def err(var):
                carb = 1 + np.abs(var[0]) * np.log(CO2_dec[:] / CO2_cmip5[150])
                clim = 1 + var[1] * (tas_dec[:, i] - np.mean(tas_ctrl[:, i]))
                return np.sum((ratio_dec - carb * clim) ** 2)


            [beta_npp[i, b], gamma_nppT[i, b]] = fmin(err, [beta_ref[0], 0], disp=False)
            beta_npp[i, b] = np.abs(beta_npp[i, b])


            def err(var):
                carb = 1 + beta_npp[i, b] * np.log(CO2_all[:] / CO2_cmip5[150])
                clim = (
                        1
                        + gamma_nppT[i, b] * (tas_all[:, i] - np.mean(tas_ctrl[:, i]))
                        + var[0] * (pr_all[:, i] - np.mean(pr_ctrl[:, i]))
                )
                return np.sum((ratio - carb * clim) ** 2)


            [gamma_nppP[i, b]] = fmin(err, [0], disp=False)
        elif mod_LSNKnpp == "hyp":

            def err(var):
                K2 = (
                             (2 * CO2_cmip5[150] - np.abs(var[1])) - (
                             1 + np.abs(var[0])) * (CO2_cmip5[150] - np.abs(var[1]))
                     ) / (
                             (1 + np.abs(var[0]) - 1) * (
                             CO2_cmip5[150] - np.abs(var[1])) * (
                                     2 * CO2_cmip5[150] - np.abs(var[1]))
                     )
                K1 = (1 + K2 * (CO2_cmip5[150] - np.abs(var[1]))) / (
                        CO2_cmip5[150] - np.abs(var[1]))
                carb = (K1 * (CO2_dec[:] - np.abs(var[1]))) / (
                        1 + K2 * (CO2_dec[:] - np.abs(var[1])))
                clim = 1 + var[2] * (tas_dec[:, i] - np.mean(tas_ctrl[:, i]))
                return np.sum((ratio_dec - carb * clim) ** 2)


            [beta_npp0[i, b], CO2_comp[i, b], gamma_nppT[i, b]] = fmin(err, [
                np.log(2) * beta_ref[0], 0, 0], disp=False)
            beta_npp0[i, b] = 1 + np.abs(beta_npp0[i, b])
            CO2_comp[i, b] = np.abs(CO2_comp[i, b])


            def err(var):
                K2 = ((2 * CO2_cmip5[150] - CO2_comp[i, b]) - beta_npp0[i, b] * (
                        CO2_cmip5[150] - CO2_comp[i, b])) / (
                             (beta_npp0[i, b] - 1) * (CO2_cmip5[150] - CO2_comp[i, b]) * (
                             2 * CO2_cmip5[150] - CO2_comp[i, b])
                     )
                K1 = (1 + K2 * (CO2_cmip5[150] - CO2_comp[i, b])) / (
                        CO2_cmip5[150] - CO2_comp[i, b])
                carb = (K1 * (CO2_all[:] - CO2_comp[i, b])) / (
                        1 + K2 * (CO2_all[:] - CO2_comp[i, b]))
                clim = (
                        1
                        + gamma_nppT[i, b] * (tas_all[:, i] - np.mean(tas_ctrl[:, i]))
                        + var[0] * (pr_all[:, i] - np.mean(pr_ctrl[:, i]))
                )
                return np.sum((ratio - carb * clim) ** 2)


            [gamma_nppP[i, b]] = fmin(err, [0], disp=False)

        # for RH rate
        ratio = (RH_all / CSOIL0_all)[:, i, b] / np.mean((RH_ctrl / CSOIL0_ctrl)[:, i, b])
        ratio_dec = (RH_dec / CSOIL0_dec)[:, i, b] / np.mean(
            (RH_ctrl / CSOIL0_ctrl)[:, i, b])
        if mod_LSNKrho == "exp":

            def err(var):
                carb = 1 + np.abs(var[0]) * (
                        FINPUT_dec[:, i, b] - np.mean(FINPUT_ctrl[:, i, b]))
                clim = np.exp(var[1] * (tas_dec[:, i] - np.mean(tas_ctrl[:, i])))
                return np.sum((ratio_dec - carb * clim) ** 2)


            first_guess = (ratio_dec[-131] - 1) / (
                    FINPUT_dec[-131, i, b] - np.mean(FINPUT_ctrl[:, i, b]))
            [prim_rho[i, b], gamma_rhoT[i, b]] = fmin(err, [first_guess,
                                                            np.log(Q10_ref[0]) / 10],
                                                      disp=False)
            prim_rho[i, b] = np.abs(prim_rho[i, b])


            def err(var):
                carb = 1 + prim_rho[i, b] * (
                        FINPUT_all[:, i, b] - np.mean(FINPUT_ctrl[:, i, b]))
                clim = np.exp(
                    gamma_rhoT[i, b] * (tas_all[:, i] - np.mean(tas_ctrl[:, i]))) * (
                               1 + var[0] * (pr_all[:, i] - np.mean(pr_ctrl[:, i]))
                       )
                return np.sum((ratio - carb * clim) ** 2)


            [gamma_rhoP[i, b]] = fmin(err, [0], disp=False)
        elif mod_LSNKrho == "gauss":

            def err(var):
                carb = 1 + np.abs(var[0]) * (
                        FINPUT_dec[:, i, b] - np.mean(FINPUT_ctrl[:, i, b]))
                clim = np.exp(
                    np.abs(var[1]) * (tas_dec[:, i] - np.mean(tas_ctrl[:, i]))
                    - np.abs(var[2]) * (tas_dec[:, i] - np.mean(tas_ctrl[:, i])) ** 2
                )
                return np.sum((ratio_dec - carb * clim) ** 2)


            first_guess = (ratio_dec[-131] - 1) / (
                    FINPUT_dec[-131, i, b] - np.mean(FINPUT_ctrl[:, i, b]))
            [prim_rho[i, b], gamma_rhoT1[i, b], gamma_rhoT2[i, b]] = fmin(
                err, [first_guess, np.log(Q10_ref[0]) / 10, 0], disp=False
            )
            prim_rho[i, b] = np.abs(prim_rho[i, b])
            gamma_rhoT1[i, b] = np.abs(gamma_rhoT1[i, b])
            gamma_rhoT2[i, b] = -np.abs(gamma_rhoT2[i, b])


            def err(var):
                carb = 1 + prim_rho[i, b] * (
                        FINPUT_all[:, i, b] - np.mean(FINPUT_ctrl[:, i, b]))
                clim = np.exp(
                    gamma_rhoT1[i, b] * (tas_all[:, i] - np.mean(tas_ctrl[:, i]))
                    + gamma_rhoT2[i, b] * (tas_all[:, i] - np.mean(tas_ctrl[:, i])) ** 2
                ) * (1 + var[0] * (pr_all[:, i] - np.mean(pr_ctrl[:, i])))
                return np.sum((ratio - carb * clim) ** 2)


            [gamma_rhoP[i, b]] = fmin(err, [0], disp=False)

# arbitrary values of parameters
# avoid diverging value of sigma
if mod_LSNKnpp == "hyp":
    beta_npp0[beta_npp0 == 1] = 1.000_001
# urban: no change to preindustrial
if (nb_biome > 1) & ((mod_biomeURB == "URB") | mod_biomeV3):
    if mod_LSNKnpp == "log":
        beta_npp[:, biome_index["urb"]] = 0
    elif mod_LSNKnpp == "hyp":
        beta_npp0[:, biome_index["urb"]] = 1.000_001
        CO2_comp[:, biome_index["urb"]] = 1e18
    gamma_nppT[:, biome_index["urb"]] = 0
    gamma_nppP[:, biome_index["urb"]] = 0
    prim_rho[:, biome_index["urb"]] = 0
    if mod_LSNKrho == "exp":
        gamma_rhoT[:, biome_index["urb"]] = 0
        gamma_rhoP[:, biome_index["urb"]] = 0
    elif mod_LSNKrho == "gauss":
        gamma_rhoT1[:, biome_index["urb"]] = 0
        gamma_rhoT2[:, biome_index["urb"]] = 0
        gamma_rhoP[:, biome_index["urb"]] = 0

# load pre-processed CMIP5 results for specified model
# data related to wildfires
for sim in ["ctrl", "upct", "fxcl", "fdbk"]:
    for VAR in ["EFIRE", "CBURNT"]:
        exec(VAR + "_" + sim + " = np.zeros([140,nb_regionI,len(bio)], dtype=dty)")
        if mod_EFIREtrans != "":
            for b in range(len(bio)):
                if os.path.isfile(
                        "data/BioBurn_CMIP5/#DATA.BioBurn_"
                        + mod_EFIREtrans
                        + ".140yr_114reg1."
                        + sim
                        + "_"
                        + VAR
                        + "_"
                        + bio[b]
                        + ".csv"
                ):
                    TMP = np.array(
                        [
                            line
                            for line in csv.reader(
                            open(
                                "data/BioBurn_CMIP5/#DATA.BioBurn_"
                                + mod_EFIREtrans
                                + ".140yr_114reg1."
                                + sim
                                + "_"
                                + VAR
                                + "_"
                                + bio[b]
                                + ".csv",
                                "r",
                            )
                        )],
                        dtype=dty,
                    )
                    for i in range(1, 114 + 1):
                        exec(VAR + "_" + sim + "[:,regionI_index[i],b] += TMP[:,i-1]")
    for var in ["tas2", "pr2"]:
        exec(var + "_" + sim + " = np.zeros([140,nb_regionI], dtype=dty)")
        if mod_EFIREtrans != "":
            TMP = np.array(
                [
                    line
                    for line in csv.reader(
                    open(
                        "data/BioBurn_CMIP5/#DATA.BioBurn_"
                        + mod_EFIREtrans
                        + ".140yr_114reg1."
                        + sim
                        + "_"
                        + var
                        + ".csv",
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
                        "data/BioBurn_CMIP5/#DATA.BioBurn_" + mod_EFIREtrans + ".140yr_114reg1.AREA2.csv",
                        "r")
                )],
                dtype=dty,
            )
            for i in range(1, 114 + 1):
                exec(var + "_" + sim + "[:,regionI_index[i]] += TMP[:,i-1]*TMP2[:,i-1]")
            TMP = np.zeros([140, nb_regionI], dtype=dty)
            for i in range(1, 114 + 1):
                exec("TMP[:,regionI_index[i]] += TMP2[:,i-1]")
            exec(var + "_" + sim + " /= TMP")
            exec(
                var + "_" + sim + "[np.isnan(" + var + "_" + sim + ")|np.isinf(" + var + "_" + sim + ")] = 0")

# complete with arbitrary values
# if no shrubs in the model
if (nb_biome > 1) & (np.sum(CBURNT_ctrl[:, :, bio.index("shr")]) == 0) & (
        mod_biomeSHR == "SHR"):
    for sim in ["ctrl", "upct", "fxcl", "fdbk"]:
        exec(
            "CBURNT_" + sim + '[:,:,bio.index("shr")] = CBURNT_' + sim + '[:,:,bio.index("gra")]')
        for VAR in ["EFIRE"]:
            exec(
                VAR
                + "_"
                + sim
                + '[:,:,bio.index("shr")] = 0.85*'
                + VAR
                + "_"
                + sim
                + '[:,:,bio.index("gra")] + 0.15*'
                + VAR
                + "_"
                + sim
                + '[:,:,bio.index("for")]*CBURNT_'
                + sim
                + '[:,:,bio.index("gra")]/CBURNT_'
                + sim
                + '[:,:,bio.index("for")]'
            )
            exec("TMP = (CBURNT_" + sim + '[:,:,bio.index("for")] == 0)')
            exec(
                VAR + "_" + sim + '[:,:,bio.index("shr")][TMP] = ' + VAR + "_" + sim + '[:,:,bio.index("gra")][TMP]')
# if no crops in the model
if (nb_biome > 1) & (np.sum(CBURNT_ctrl[:, :, bio.index("cro")]) == 0):
    for sim in ["ctrl", "upct", "fxcl", "fdbk"]:
        for VAR in ["CBURNT", "EFIRE"]:
            exec(
                VAR + "_" + sim + '[:,:,bio.index("cro")] = ' + VAR + "_" + sim + '[:,:,bio.index("gra")]')
# if no pastures in the model
if (nb_biome > 1) & (np.sum(CBURNT_ctrl[:, :, bio.index("pas")]) == 0):
    for sim in ["ctrl", "upct", "fxcl", "fdbk"]:
        exec(
            "CBURNT_" + sim + '[:,:,bio.index("pas")] = CBURNT_' + sim + '[:,:,bio.index("gra")]')
        for VAR in ["EFIRE"]:
            exec(
                VAR
                + "_"
                + sim
                + '[:,:,bio.index("pas")] = 0.60*'
                + VAR
                + "_"
                + sim
                + '[:,:,bio.index("gra")] + 0.40*'
                + VAR
                + "_"
                + sim
                + '[:,:,bio.index("des")]*CBURNT_'
                + sim
                + '[:,:,bio.index("gra")]/CBURNT_'
                + sim
                + '[:,:,bio.index("des")]'
            )
            exec("TMP = (CBURNT_" + sim + '[:,:,bio.index("des")] == 0)')
            exec(
                VAR + "_" + sim + '[:,:,bio.index("pas")][TMP] = ' + VAR + "_" + sim + '[:,:,bio.index("gra")][TMP]')
# if no urban in the model
if (nb_biome > 1) & (np.sum(CBURNT_ctrl[:, :, bio.index("urb")]) == 0) & (
        (mod_biomeURB == "URB") | mod_biomeV3):
    for sim in ["ctrl", "upct", "fxcl", "fdbk"]:
        for VAR in ["CBURNT", "EFIRE"]:
            exec(
                VAR + "_" + sim + '[:,:,bio.index("urb")] = ' + VAR + "_" + sim + '[:,:,bio.index("des")]')

# aggregate biomes
for sim in ["ctrl", "upct", "fxcl", "fdbk"]:
    for VAR in ["EFIRE", "CBURNT"]:
        exec("TMP = " + VAR + "_" + sim + ".copy()")
        exec(VAR + "_" + sim + " = np.zeros([140,nb_regionI,nb_biome], dtype=dty)")
        for b in range(len(bio)):
            exec(VAR + "_" + sim + "[:,:,biome_index[bio[b]]] += TMP[:,:,b]")
# aggregate all experiments
CO2_all = np.array(
    list(CO2_cmip5[150] * 1.01 ** np.arange(140))
    + list(CO2_cmip5[150] * 1.01 ** np.arange(140))
    + list(CO2_cmip5[150] * np.ones(140)),
    dtype=dty,
)
for VAR in ["EFIRE", "CBURNT"] + ["tas2", "pr2"]:
    exec(
        VAR + "_all = np.array(list(" + VAR + "_upct)+list(" + VAR + "_fxcl)+list(" + VAR + "_fdbk), dtype=dty)")
# decadal means
for VAR in ["EFIRE", "CBURNT"] + ["tas2", "pr2"]:
    exec(
        VAR + "_dec= np.zeros([3*(140-10)]+list(np.shape(" + VAR + "_all)[1:]), dtype=dty)")
    for t in range(130):
        exec(VAR + "_dec[t,...] = np.mean(" + VAR + "_all[t:t+10,...],0)")
        exec(VAR + "_dec[t+140-10,...] = np.mean(" + VAR + "_all[140+t:140+t+10,...],0)")
        exec(
            VAR + "_dec[t+2*140-2*10,...] = np.mean(" + VAR + "_all[2*140+t:2*140+t+10,...],0)")

# definition of parameters
# sensitivity of ignition rate to climate {/ppm}&{/K}&{/mm}
gamma_igniC = np.zeros([nb_regionI, nb_biome], dtype=dty)
gamma_igniT = np.zeros([nb_regionI, nb_biome], dtype=dty)
gamma_igniP = np.zeros([nb_regionI, nb_biome], dtype=dty)

# fit of parameters
for i in range(1, nb_regionI):
    for b in range(nb_biome):

        # for ignition rate
        if mod_EFIREtrans != "":
            ratio = (EFIRE_all / CBURNT_all)[:, i, b] / np.mean(
                (EFIRE_ctrl / CBURNT_ctrl)[:, i, b])
            ratio_dec = (EFIRE_dec / CBURNT_dec)[:, i, b] / np.mean(
                (EFIRE_ctrl / CBURNT_ctrl)[:, i, b])


            def err(var):
                carb = 1
                clim = (
                        1
                        + var[0] * (CO2_dec[:] - CO2_cmip5[150])
                        + np.abs(var[1]) * (tas2_dec[:, i] - np.mean(tas2_ctrl[:, i]))
                        - np.abs(var[2]) * (pr2_dec[:, i] - np.mean(pr2_ctrl[:, i]))
                )
                return np.sum((ratio_dec - carb * clim) ** 2)


            first_guess = (ratio_dec[-131] - 1) / (CO2_dec[-131] - CO2_cmip5[150])
            [gamma_igniC[i, b], gamma_igniT[i, b], gamma_igniP[i, b]] = fmin(err,
                                                                             [first_guess,
                                                                              0, 0],
                                                                             disp=False)
            gamma_igniT[i, b] = np.abs(gamma_igniT[i, b])
            gamma_igniP[i, b] = -np.abs(gamma_igniP[i, b])

# arbitrary values of parameters
# crops: no wildfire
if nb_biome > 1:
    gamma_igniC[:, biome_index["cro"]] = 0.0
    gamma_igniT[:, biome_index["cro"]] = 0.0
    gamma_igniP[:, biome_index["cro"]] = 0.0
# urban: no wildfire
if (nb_biome > 1) & ((mod_biomeURB == "URB") | mod_biomeV3):
    gamma_igniC[:, biome_index["urb"]] = 0.0
    gamma_igniT[:, biome_index["urb"]] = 0.0
    gamma_igniP[:, biome_index["urb"]] = 0.0

# -----------------
# 1.3.5. Permafrost
# -----------------

# permafrost regions
regionPF = ["EUAS", "NAM"]
regionPF_name = ["Eurasia", "North America"]
nb_regionPF = len(regionPF)


# function for speed of thawing/freezing
def f_v_PF(pthaw_bar, pthaw):
    v = v_thaw * (pthaw_bar - pthaw > 0) + v_froz * (pthaw_bar - pthaw < 0)
    return np.array(v, dtype=dty)


# load other calibrated parameters
if mod_EPFmain == "JSBACH":

    # regional weighting for local temperature {.}
    w_reg_lstPF = np.array([1.856_835_9, 1.947_070_3], dtype=dty)

    # sensitivities of local heterotrophic respiration {/K}&{/K2}&{.}
    gamma_rhoPF1 = np.array([0.096_905_43, 0.100_760_28], dtype=dty)
    gamma_rhoPF2 = -np.array([0.002_063_56, 0.002_141_74], dtype=dty)
    w_rhoPF = np.array([0.938_127_85, 0.734_540_44], dtype=dty)

    # parameters for thawing function {.}&{/K}&{.}
    pthaw_min = np.array([0.253_198_50, 0.195_031_77], dtype=dty)
    k_pthaw = np.array([2.301_288_7, 1.417_431_2], dtype=dty)
    gamma_pthaw = np.array([0.175_823_92, 0.286_380_24], dtype=dty)

    # speed of thawing/freezing {/yr}
    v_thaw = np.array([4.894_087_3, 0.499_786_20], dtype=dty)
    v_froz = np.array([0.0, 0.0], dtype=dty)

    # weights and turnover times for the three emitting pools {.}&{yr}
    p_PF1 = np.array([0.069_769_70, 0.099_873_57], dtype=dty)
    p_PF2 = np.array([0.263_969_75, 0.234_255_85], dtype=dty)
    p_PF3 = np.array([0.666_260_55, 0.665_870_58], dtype=dty)
    tau_PF1 = np.array([7.961_624_9, 11.632_080], dtype=dty)
    tau_PF2 = np.array([88.489_543, 103.06668], dtype=dty)
    tau_PF3 = np.array([1355.6291, 1240.6807], dtype=dty)

    # initial pool of frozen carbon {GtC}
    CFROZ_0 = np.array([413.553, 277.483], dtype=dty)

elif mod_EPFmain == "ORCHIDEE-MICT":

    # regional weighting for local temperature {.}
    w_reg_lstPF = np.array([1.867_871_0, 1.782_812_5], dtype=dty)

    # sensitivities of local heterotrophic respiration {/K}&{/K2}&{.}
    gamma_rhoPF1 = np.array([0.114_007_56, 0.109_569_53], dtype=dty)
    gamma_rhoPF2 = -np.array([0.003_663_66, 0.003_454_03], dtype=dty)
    w_rhoPF = np.array([3.380_518_3, 3.707_680_2], dtype=dty)

    # parameters for thawing function {.}&{/K}&{.}
    pthaw_min = np.array([0.080_653_49, 0.875_467_73], dtype=dty)
    k_pthaw = np.array([0.356_892_19, 8.442_261_9], dtype=dty)
    gamma_pthaw = np.array([1.618_303_2, 0.097_257_53], dtype=dty)

    # speed of thawing/freezing {/yr}
    v_thaw = np.array([0.286_963_67, 0.499_533_14], dtype=dty)
    v_froz = np.array([0.052_892_77, 0.049_027_89], dtype=dty)

    # weights and turnover times for the three emitting pools {.}&{yr}
    p_PF1 = np.array([0.243_892_66, 0.288_349_40], dtype=dty)
    p_PF2 = np.array([0.756_107_34, 0.711_650_60], dtype=dty)
    p_PF3 = np.array([0.0, 0.0], dtype=dty)
    tau_PF1 = np.array([1332.6992, 2950.1099], dtype=dty)
    tau_PF2 = np.array([14953.331, 23437.897], dtype=dty)
    tau_PF3 = np.array([np.inf, np.inf], dtype=dty)

    # initial pool of frozen carbon {GtC}
    CFROZ_0 = np.array([270.56425, 118.19786], dtype=dty)

elif mod_EPFmain == "JULES-DeepResp":

    # regional weighting for local temperature {.}
    w_reg_lstPF = np.array([1.960_351_5, 2.025_585_9], dtype=dty)

    # sensitivities of local heterotrophic respiration {/K}&{/K2}&{.}
    gamma_rhoPF1 = np.array([0.168_177_39, 0.143_768_88], dtype=dty)
    gamma_rhoPF2 = -np.array([0.005_312_66, 0.003_796_30], dtype=dty)
    w_rhoPF = np.array([1.374_876_2, 1.582_355_5], dtype=dty)

    # parameters for thawing function {.}&{/K}&{.}
    pthaw_min = np.array([0.772_467_63, 0.958_939_20], dtype=dty)
    k_pthaw = np.array([2.052_901_6, 1.914_339_4], dtype=dty)
    gamma_pthaw = np.array([0.153_014_8, 0.144_579_36], dtype=dty)

    # speed of thawing/freezing {/yr}
    v_thaw = np.array([0.266_194_01, 0.598_990_92], dtype=dty)
    v_froz = np.array([0.041_027_76, 0.032_966_77], dtype=dty)

    # weights and turnover times for the three emitting pools {.}&{yr}
    p_PF1 = np.array([0.033_667_69, 0.651_148_42], dtype=dty)
    p_PF2 = np.array([0.966_332_31, 0.348_851_58], dtype=dty)
    p_PF3 = np.array([0.0, 0.0], dtype=dty)
    tau_PF1 = np.array([156.38595, 1969.9072], dtype=dty)
    tau_PF2 = np.array([3157.9293, 5370.1557], dtype=dty)
    tau_PF3 = np.array([np.inf, np.inf], dtype=dty)

    # initial pool of frozen carbon {GtC}
    CFROZ_0 = np.array([480.92224, 175.53549], dtype=dty)

elif mod_EPFmain == "JULES-SuppressResp":

    # regional weighting for local temperature {.}
    w_reg_lstPF = np.array([1.959_082_0, 2.024_218_7], dtype=dty)

    # sensitivities of local heterotrophic respiration {/K}&{/K2}&{.}
    gamma_rhoPF1 = np.array([0.134_720_92, 0.111_936_08], dtype=dty)
    gamma_rhoPF2 = -np.array([0.004_988_00, 0.003_394_37], dtype=dty)
    w_rhoPF = np.array([0.661_802_79, 2.185_352_3], dtype=dty)

    # parameters for thawing function {.}&{/K}&{.}
    pthaw_min = np.array([0.870_200_23, 1.164_973_9], dtype=dty)
    k_pthaw = np.array([1.630_524_6, 1.683_605_4], dtype=dty)
    gamma_pthaw = np.array([0.179_786_47, 0.190_059_63], dtype=dty)

    # speed of thawing/freezing {/yr}
    v_thaw = np.array([0.373_391_28, 0.492_236_81], dtype=dty)
    v_froz = np.array([0.061_284_38, 0.042_925_02], dtype=dty)

    # weights and turnover times for the three emitting pools {.}&{yr}
    p_PF1 = np.array([1.0, 1.0], dtype=dty)
    p_PF2 = np.array([0.0, 0.0], dtype=dty)
    p_PF3 = np.array([0.0, 0.0], dtype=dty)
    tau_PF1 = np.array([9433.4705, 4322.9007], dtype=dty)
    tau_PF2 = np.array([np.inf, np.inf], dtype=dty)
    tau_PF3 = np.array([np.inf, np.inf], dtype=dty)

    # initial pool of frozen carbon {GtC}
    CFROZ_0 = np.array([372.75814, 121.20251], dtype=dty)

elif mod_EPFmain == "":

    # regional weighting for local temperature {.}
    w_reg_lstPF = np.array([1.0, 1.0], dtype=dty)

    # sensitivities of local heterotrophic respiration {/K}&{/K2}&{.}
    gamma_rhoPF1 = np.array([0.0, 0.0], dtype=dty)
    gamma_rhoPF2 = -np.array([0.0, 0.0], dtype=dty)
    w_rhoPF = np.array([1.0, 1.0], dtype=dty)

    # parameters for thawing function {.}&{/K}&{.}
    pthaw_min = np.array([0.0, 0.0], dtype=dty)
    k_pthaw = np.array([1.0, 1.0], dtype=dty)
    gamma_pthaw = np.array([0.0, 0.0], dtype=dty)

    # speed of thawing/freezing {/yr}
    v_thaw = np.array([1.0, 1.0], dtype=dty)
    v_froz = np.array([0.0, 0.0], dtype=dty)

    # weights and turnover times for the three emitting pools {.}&{yr}
    p_PF1 = np.array([0.0, 0.0], dtype=dty)
    p_PF2 = np.array([0.0, 0.0], dtype=dty)
    p_PF3 = np.array([0.0, 0.0], dtype=dty)
    tau_PF1 = np.array([np.inf, np.inf], dtype=dty)
    tau_PF2 = np.array([np.inf, np.inf], dtype=dty)
    tau_PF3 = np.array([np.inf, np.inf], dtype=dty)

    # initial pool of frozen carbon {GtC}
    CFROZ_0 = np.array([0.0, 0.0], dtype=dty)

# =============
# 1.4. LAND-USE
# =============

# -------------
# 1.4.1. TRENDY
# -------------

# basic biomes of aggregation
bio = ["des", "for", "shr", "gra", "cro", "pas", "urb"]

# load pre-processed TRENDY results for specified model
# data related to preindustrial C-cycle
for VAR in ["AGB", "BGB"]:
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/WoodUse_TRENDYv2/#DATA.WoodUse_" + mod_ELUCagb + ".1910s_114reg1_7bio." + VAR + ".csv",
                "r")
        )],
        dtype=dty,
    )
    exec(VAR + "_pi = np.zeros([nb_regionI,len(bio)], dtype=dty)")
    for i in range(1, 114 + 1):
        exec(VAR + "_pi[regionI_index[i],:] += TMP[i-1,:]")

# complete with arbitrary values
# if no shrubs in the model
if (nb_biome > 1) & (np.sum((AGB_pi + BGB_pi)[:, bio.index("shr")]) == 0) & (
        mod_biomeSHR == "SHR"):
    for VAR in ["AGB", "BGB"]:
        exec(
            VAR
            + '_pi[:,bio.index("shr")] = 0.85*'
            + VAR
            + '_pi[:,bio.index("gra")] + 0.15*'
            + VAR
            + '_pi[:,bio.index("for")]'
        )
# if no crops in the model
if (nb_biome > 1) & (np.sum((AGB_pi + BGB_pi)[:, bio.index("cro")]) == 0):
    for VAR in ["AGB", "BGB"]:
        exec(VAR + '_pi[:,bio.index("cro")] = ' + VAR + '_pi[:,bio.index("gra")]')
# if no pastures in the model
if (nb_biome > 1) & (np.sum((AGB_pi + BGB_pi)[:, bio.index("pas")]) == 0):
    for VAR in ["AGB", "BGB"]:
        exec(
            VAR
            + '_pi[:,bio.index("pas")] = 0.60*'
            + VAR
            + '_pi[:,bio.index("gra")] + 0.40*'
            + VAR
            + '_pi[:,bio.index("des")]'
        )
# if no urban in the model
if (nb_biome > 1) & (np.sum((AGB_pi + BGB_pi)[:, bio.index("urb")]) == 0) & (
        (mod_biomeURB == "URB") | mod_biomeV3):
    for VAR in ["AGB", "BGB"]:
        exec(VAR + '_pi[:,bio.index("urb")] = ' + VAR + '_pi[:,bio.index("des")]')

# aggregate biomes
for VAR in ["AGB", "BGB"]:
    exec("TMP = " + VAR + "_pi.copy()")
    exec(VAR + "_pi = np.zeros([nb_regionI,nb_biome], dtype=dty)")
    for b in range(len(bio)):
        exec(VAR + "_pi[:,biome_index[bio[b]]] += TMP[:,b]")

# calculate preindustrial flux parameters
p_AGB = AGB_pi / (AGB_pi + BGB_pi)
p_AGB[np.isnan(p_AGB)] = 0

# ---------------
# 1.4.2. Wood-Use
# ---------------

# characteristic time for shifting cultivation {yr}
# adjusted; 15yr is from [Hurtt et al., 2006]
tau_shift = np.array([15], dtype=dty)

# aggregation of fractions for wood-use {.}
# (n=0) for slash; (n=1) for burnt; (n=2,3) for wood decay
# base data from [Earles et al., 2012]
for n in range(4):
    exec("p_HWP" + str(n) + " = np.zeros([nb_regionI,nb_biome], dtype=dty)")
TMP = np.array(
    [line for line in csv.reader(
        open("data/WoodUse_Earles/#DATA.WoodUse_Earles.2000s_114reg1_(4use).HWP.csv",
             "r"))][
    1:],
    dtype=dty,
)
for i in range(1, 114 + 1):
    for n in range(4):
        exec("p_HWP" + str(n) + "[regionI_index[i],:] += TMP[i-1,n]")
# biomass-burning of non-merchantable AGB
if mod_EHWPbb == "high":
    p_HWP1 += 0.5 * p_HWP0
    p_HWP0 -= 0.5 * p_HWP0
elif mod_EHWPbb == "low":
    p_HWP1 += 0.0 * p_HWP0
    p_HWP0 -= 0.0 * p_HWP0
# remove some products for some biomes
if nb_biome > 1:
    for bio in ["des", "shr", "gra"]:
        p_HWP2[:, biome_index[bio]] = p_HWP3[:, biome_index[bio]] = 0
    for bio in ["cro", "pas"]:
        p_HWP1[:, biome_index[bio]] = p_HWP2[:, biome_index[bio]] = p_HWP3[:,
                                                                    biome_index[bio]] = 0
    if (mod_biomeURB == "URB") | mod_biomeV3:
        p_HWP1[:, biome_index["urb"]] = p_HWP2[:, biome_index["urb"]] = p_HWP3[:,
                                                                        biome_index[
                                                                            "urb"]] = 0
# normalize
TMP = p_HWP0 + p_HWP1 + p_HWP2 + p_HWP3
for n in range(4):
    exec("p_HWP" + str(n) + " /= TMP")
    exec("p_HWP" + str(n) + "[np.isinf(p_HWP" + str(n) + ")|np.isnan(p_HWP" + str(
        n) + ")] = 0")
p_HWP0 = 1 - p_HWP1 - p_HWP2 - p_HWP3

# fraction of actually burnt HWP1 {}
# used to adjust non-CO2 anthropogenic BB emissions
p_HWP1_bb = 0.5

# wood-use decay constants {yr}
# based roughly on [Houghton et Hackler, 2001]
if mod_EHWPtau == "Houghton2001":
    tau_HWP1 = np.array([1], dtype=dty)
    tau_HWP2 = np.array([10], dtype=dty)
    tau_HWP3 = np.array([100], dtype=dty)
# based on [Earles et al., 2012]
elif mod_EHWPtau == "Earles2012":
    tau_HWP1 = np.array([0.5], dtype=dty)
    tau_HWP2 = np.array([2], dtype=dty)
    tau_HWP3 = np.array([30], dtype=dty)

# decay times variations
# fast: 20% remaining after 0.8*tau
if mod_EHWPspeed == "fast":
    tau_HWP1 *= 0.8 / -np.log(0.2)
    tau_HWP2 *= 0.8 / -np.log(0.2)
    tau_HWP3 *= 0.8 / -np.log(0.2)
# fast: 30% remaining after 1.5*tau
elif mod_EHWPspeed == "slow":
    tau_HWP1 *= 1.5 / -np.log(0.3)
    tau_HWP2 *= 1.5 / -np.log(0.3)
    tau_HWP3 *= 1.5 / -np.log(0.3)

# IRF for wood product pools {.}
# forced to be exponential in this version
mod_EHWPfct = "exp"
if mod_EHWPfct == "gamma":

    def f_HWP(tau):
        err = lambda a0: gammainc(a0, (a0 + 1) / 2.5) - 0.05
        a0 = fsolve(err, 0.1)
        HWP = gammainc(a0, tau * (a0 + 1) / np.arange(ind_final + 1))
        HWP[0] = 1.0
        return np.array(HWP, dtype=dty)


elif mod_EHWPfct == "lin":

    def f_HWP(tau):
        HWP = (1 - np.arange(ind_final + 1) / tau) * (np.arange(ind_final + 1) < tau)
        return np.array(HWP, dtype=dty)


elif mod_EHWPfct == "exp":

    def f_HWP(tau):
        HWP = np.exp(-np.arange(ind_final + 1) / tau)
        return np.array(HWP, dtype=dty)

# IRF for wood product fluxes ratio {/yr}
for N in ["1", "2", "3"]:
    exec("r_HWP" + N + " = np.zeros([ind_final+1], dtype=dty)")
    exec(
        "r_HWP"
        + N
        + "[1:] = -(f_HWP(tau_HWP"
        + N
        + ")[1:]-f_HWP(tau_HWP"
        + N
        + ")[:-1])/(f_HWP(tau_HWP"
        + N
        + ")[1:]+f_HWP(tau_HWP"
        + N
        + ")[:-1])/0.5"
    )
    exec("r_HWP" + N + "[np.isnan(r_HWP" + N + ")|np.isinf(r_HWP" + N + ")] = 1")
    # exec('r_HWP'+N+'[r_HWP'+N+'>1] = 1')
    exec("r_HWP" + N + "[0] = 0")

# -----------
# 1.4.3. GFED
# -----------

# biomes of GFED and correspondence
bio = ["def", "for", "woo", "sav", "agr", "pea"]
BIO = {"def": "for", "for": "for", "woo": "shr", "sav": "gra", "agr": "cro", "pea": "pea"}

# load GFED v3.1 BB emissions {TgX/GtC}&{Tg/GtC}
# from [Randerson et al., 2013]
for VAR in ["CO2"] + ["CH4", "NOX", "CO", "VOC", "N2O", "SO2", "NH3", "BC", "OC"]:
    exec("alpha_BB_" + VAR + " = np.zeros([nb_regionI,nb_biome], dtype=dty)")
    for b in range(len(bio) - 1):
        if os.path.isfile(
                "data/BioBurn_GFED/#DATA.BioBurn_GFED.1997-2011_114reg1.E" + VAR + "_" +
                bio[b] + ".csv"):
            TMP = np.array(
                [
                    line
                    for line in csv.reader(
                    open(
                        "data/BioBurn_GFED/#DATA.BioBurn_GFED.1997-2011_114reg1.E" + VAR + "_" +
                        bio[b] + ".csv",
                        "r",
                    )
                )],
                dtype=dty,
            )
            for i in range(1, 114 + 1):
                exec(
                    "alpha_BB_" + VAR + "[regionI_index[i],biome_index[BIO[bio[b]]]] += np.sum(TMP[:,i-1])")
                # set pastures and baresoil as grasslands
                if bio[b] == "sav":
                    exec(
                        "alpha_BB_" + VAR + '[regionI_index[i],biome_index["des"]] += np.sum(TMP[:,i-1])')
                    exec(
                        "alpha_BB_" + VAR + '[regionI_index[i],biome_index["pas"]] += np.sum(TMP[:,i-1])')
for VAR in ["CH4", "NOX", "CO", "VOC", "N2O", "SO2", "NH3", "BC", "OC"] + ["CO2"]:
    exec("alpha_BB_" + VAR + " /= alpha_BB_CO2")
    exec(
        "alpha_BB_" + VAR + "[np.isnan(alpha_BB_" + VAR + ")|np.isinf(alpha_BB_" + VAR + ")] = 0")
    for b in range(nb_biome):
        exec(
            "alpha_BB_"
            + VAR
            + "[:,b][alpha_BB_"
            + VAR
            + "[:,b]==0] = np.sum(alpha_BB_"
            + VAR
            + "[:,b])/np.sum(alpha_BB_"
            + VAR
            + "[:,b]!=0)"
        )
    exec(
        "alpha_BB_" + VAR + "[np.isnan(alpha_BB_" + VAR + ")|np.isinf(alpha_BB_" + VAR + ")] = 0")
