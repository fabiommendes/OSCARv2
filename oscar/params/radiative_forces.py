import numpy as np

from .carbon import CO2_0
from .methane import CH4_0
from .nitrous_oxide import N2O_0
from ..data_loaders import nb_regionI, nb_biome, regionI_index, biome_index
from ..config import dty, mod_O3Tradeff, mod_SO4radeff, mod_POAradeff, mod_BCradeff, mod_NO3radeff, mod_SOAradeff, mod_BCadjust, mod_CLOUDsolub, mod_CLOUDerf, mod_CLOUDpreind, mod_ALBBCrf, mod_ALBBCreg, mod_ALBBCwarm, mod_ALBLCalb, mod_ALBLCflux, mod_ALBLCcover, mod_ALBLCwarm, mod_O3Sradeff
from ..data import load_data, load_data_and_header

##################################################
#   7. RADIATIVE FORCING
##################################################

# load RCP radiative forcing {W/m2}
# from [Meinshausen et al., 2011]
RF_rcp = np.zeros([800 + 1, 6], dtype=dty)
RF_rcp[:300] = np.nan
n = -1
for rcp in ["rcp26", "rcp45", "rcp60", "rcp85", "rcp45to26", "rcp60to45"]:
    n += 1
    path = f"data/Scenario_ECP/#DATA.Scenario_ECP.2000-2500_(19for).{rcp}_RF.csv"
    path2 = "data/Historic_CMIP5/#DATA.Historic_CMIP5.1765-2005_(19for).RF.csv"
    TMP, lgd = load_data_and_header(path)
    TMP2 = load_data(path, start=1)
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
    RF = 0.47 * np.log(1
        + 2.01e-5 * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 0.75
        + 5.31e-15 * (CH4_0 + D_CH4) * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 1.52)
    RF -= 0.47 * np.log(1 + 2.01e-5 * (CH4_0 * N2O_0) ** 0.75 + 5.31e-15 * CH4_0 * (CH4_0 * N2O_0) ** 1.52)
    return np.array(RF, dtype=dty)


def df_RF_overlap_dCH4(D_CH4, D_N2O):
    RF = 0.47 * (2.01e-5 * 0.75 * (CH4_0 + D_CH4) ** (0.75 - 1) * (N2O_0 + D_N2O) ** (0.75) + 5.31e-15 * (1.52 + 1) * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 1.52)
    RF /= (1 + 2.01e-5 * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 0.75 + 5.31e-15 * (CH4_0 + D_CH4) * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 1.52)
    return np.array(RF, dtype=dty)


def df_RF_overlap_dN2O(D_CH4, D_N2O):
    RF = 0.47 * (2.01e-5 * 0.75 * (CH4_0 + D_CH4) ** (0.75) * (N2O_0 + D_N2O) ** (0.75 - 1) + 5.31e-15 * 1.52 * (CH4_0 + D_CH4) ** (1.52 + 1) * (N2O_0 + D_N2O) ** (1.52 - 1))
    RF /= (1 + 2.01e-5 * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 0.75 + 5.31e-15 * (CH4_0 + D_CH4) * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 1.52)
    return np.array(RF, dtype=dty)


def df_RF_overlap(D_CH4, D_N2O):
    df_1 = df_RF_overlap_dCH4(D_CH4, D_N2O) * D_CH4
    df_2 = df_RF_overlap_dN2O(D_CH4, D_N2O) * D_N2O
    df_tot = df_1 + df_2
    return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot)]


# radiative efficiency of halogenated compounds {{W/m2}/ppt}
# from IPCC-AR5 [Myrhe et al., 2013] (table 8.A.1)
radeff_HFC = 1e-3 * np.array([0.18, 0.11, 0.23, 0.16, 0.16, 0.10, 0.26, 0.24, 0.24, 0.22, 0.42], dtype=dty)
radeff_PFC = 1e-3 * np.array([0.57, 0.20, 0.09, 0.25, 0.28, 0.32, 0.36, 0.41, 0.44, 0.50], dtype=dty)
radeff_ODS = 1e-3 * np.array([0.26, 0.32, 0.30, 0.31, 0.20, 0.17, 0.07, 0.21, 0.16, 0.19, 0.29, 0.27, 0.30, 0.31, 0.004, 0.01], dtype=dty)

# ==========
# 7.2. OZONE
# ==========

# radiative efficiency of tropospheric O3 {{W/m2}/DU}
# from IPCC-AR5 [Myhre et al., 2013]
radeff_O3t_map = {
    "IPCC-AR5": 0.042,
    # from IPCC-AR4 [Forster et al., 2007]
    "IPCC-AR4": 0.032,
    # from ACCMIP [Stevenson et al., 2013] (table 3)
    "mean-ACCMIP": 0.377 / 8.9,
    "CESM-CAM-superfast": 0.446 / 10.0,
    "CICERO-OsloCTM2": 0.401 / 9.3,
    "CMAM": 0.322 / 7.6,
    "EMAC": 0.460 / 10.8,
    "GEOSCCM": 0.387 / 8.7,
    "GFDL-AM3": 0.423 / 10.3,
    "GISS-E2-R": 0.314 / 8.3,
    "GISS-E2-R-TOMAS": 0.333 / 8.7,
    "HadGEM2": 0.303 / 7.3,
    "LMDzORINCA": 0.351 / 8.2,
    "MIROC-CHEM": 0.402 / 9.2,
    "MOCAGE": 0.219 / 4.8,
    "NCAR-CAM-35": 0.433 / 10.2,
    "STOC-HadAM3": 0.437 / 10.5,
    "UM-CAM": 0.376 / 8.7,
    "TM5": 0.422 / 10.0,
}
radeff_O3t = radeff_O3t_map[mod_O3Tradeff]

# radiative efficiency of stratospheric O3 {{W/m2}/DU}
# from IPCC-AR4 [Forster et al., 2007]
radeff_O3s_map = {
    "IPCC-AR4":0.004,
    # ACCENT [Gauss et al., 2006] (tables 4 & 6)
    "mean-ACCENT":-0.058 / -13.9,
    "ULAQ":-0.059 / -12.6,
    "DLR-E39C":-0.027 / -16.1,
    "NCAR-MACCM":-0.019 / -12.7,
    "CHASER":-0.126 / -14.1,
}
radeff_O3s = radeff_O3s_map[mod_O3Sradeff]


# =============
# 7.3. AEROSOLS
# =============

# -------------
# 7.3.1. Direct
# -------------

#: Radiative efficiency of sulfate aerosols {{W/m2}/Tg}, AeroCom2 [Myhre et al., 2013] (table 4)
radeff_SO4_map = {
    "mean-AeroCom2": (-185) * 1e12 / 510_072e9,
    "BCC": (-108) * 1e12 / 510_072e9,
    "CAM4-Oslo": (-173) * 1e12 / 510_072e9,
    "CAM-51": (-104) * 1e12 / 510_072e9,
    "GEOS-CHEM": (-123) * 1e12 / 510_072e9,
    "GISS-MATRIX": (-196) * 1e12 / 510_072e9,
    "GISS-modelE": (-307) * 1e12 / 510_072e9,
    "GMI": (-195) * 1e12 / 510_072e9,
    "GOCART": (-238) * 1e12 / 510_072e9,
    "HadGEM2": (-193) * 1e12 / 510_072e9,
    "IMPACT-Umich": (-113) * 1e12 / 510_072e9,
    "INCA": (-180) * 1e12 / 510_072e9,
    "MPIHAM": (-125) * 1e12 / 510_072e9,
    "NCAR-CAM-35": (-354) * 1e12 / 510_072e9,
    "OsloCTM2": (-192) * 1e12 / 510_072e9,
    "SPRINTARS": (-172) * 1e12 / 510_072e9,
}
radeff_SO4 = radeff_SO4_map[mod_SO4radeff]

# radiative efficiency of primary organic aerosols {{W/m2}/Tg}
# from AeroCom2 [Myhre et al., 2013] (table 6)
radeff_POA_map = {
    "mean-AeroCom2": (-113) * 1e12 / 510_072e9,
    "BCC": (-97) * 1e12 / 510_072e9,
    "CAM4-Oslo": (-118) * 1e12 / 510_072e9,
    "CAM-51": (-69) * 1e12 / 510_072e9,
    "GEOS-CHEM": (-95) * 1e12 / 510_072e9,
    "GISS-MATRIX": (-129) * 1e12 / 510_072e9,
    "GISS-modelE": (-76) * 1e12 / 510_072e9,
    "GMI": (-189) * 1e12 / 510_072e9,
    "GOCART": (-144) * 1e12 / 510_072e9,
    "HadGEM2": (-145) * 1e12 / 510_072e9,
    "IMPACT-Umich": (-141) * 1e12 / 510_072e9,
    "INCA": (-76) * 1e12 / 510_072e9,
    "MPIHAM": (-41) * 1e12 / 510_072e9,
    "NCAR-CAM-35": (-48) * 1e12 / 510_072e9,
    "OsloCTM2": (-165) * 1e12 / 510_072e9,
    "SPRINTARS": (-102) * 1e12 / 510_072e9,
}
radeff_POA = radeff_POA_map[mod_POAradeff]

# radiative efficiency of black carbon aerosols {{W/m2}/Tg}
# from AeroCom2 [Myhre et al., 2013] (table 5)
radeff_BC_map = {
    "mean-AeroCom2": 1e12 / 510_072e9 * 1438,
    "BCC": 1e12 / 510_072e9 * 650,
    "CAM4-Oslo": 1e12 / 510_072e9 * 1763,
    "CAM-51": 1e12 / 510_072e9 * 2661,
    "GEOS-CHEM": 1e12 / 510_072e9 * 1067,
    "GISS-MATRIX": 1e12 / 510_072e9 * 2484,
    "GISS-modelE": 1e12 / 510_072e9 * 1253,
    "GMI": 1e12 / 510_072e9 * 1208,
    "GOCART": 1e12 / 510_072e9 * 874,
    "HadGEM2": 1e12 / 510_072e9 * 612,
    "IMPACT-Umich": 1e12 / 510_072e9 * 1467,
    "INCA": 1e12 / 510_072e9 * 1160,
    "MPIHAM": 1e12 / 510_072e9 * 1453,
    "NCAR-CAM-35": 1e12 / 510_072e9 * 1364,
    "OsloCTM2": 1e12 / 510_072e9 * 2161,
    "SPRINTARS": 1e12 / 510_072e9 * 1322,
}
radeff_BC = radeff_BC_map[mod_BCradeff]

# radiative efficiency of nitrate aerosols {{W/m2}/Tg}
# from AeroCom2 [Myhre et al., 2013] (table 8)
radeff_NO3_map = {
    "mean-AeroCom2": (-166) * 1e12 / 510_072e9,
    "GEOS-CHEM": (-136) * 1e12 / 510_072e9,
    "GISS-MATRIX": (-240) * 1e12 / 510_072e9,
    "GMI": (-103) * 1e12 / 510_072e9,
    "HadGEM2": (-249) * 1e12 / 510_072e9,
    "IMPACT-Umich": (-155) * 1e12 / 510_072e9,
    "INCA": (-110) * 1e12 / 510_072e9,
    "NCAR-CAM-35": (-91) * 1e12 / 510_072e9,
    "OsloCTM2": (-173) * 1e12 / 510_072e9,
}
radeff_NO3 = radeff_NO3_map[mod_NO3radeff]

# radiative efficiency of secondary organic aerosols {{W/m2}/Tg}
# from AeroCom2 [Myhre et al., 2013] (table 7)
radeff_SOA_map = {
    "mean-AeroCom2": (-122) * 1e12 / 510_072e9,
    "CAM-51": (-45) * 1e12 / 510_072e9,
    "GEOS-CHEM": (-45) * 1e12 / 510_072e9,
    "IMPACT-Umich": (-218) * 1e12 / 510_072e9,
    "MPIHAM": (-139) * 1e12 / 510_072e9,
    "OsloCTM2": (-161) * 1e12 / 510_072e9,
}
radeff_SOA = radeff_SOA_map[mod_SOAradeff]

# radiative efficiency of natural aerosols {{W/m2}/Tg}
# set to zero in this version
radeff_DUST = 1e12 / 510_072e9 * 0
radeff_SALT = 1e12 / 510_072e9 * 0

# ---------------
# 7.3.2. Indirect
# ---------------

#: Semi-direct effect (adjustements induced by the direct RF of BC), best-guess from IPCC AR5 [Boucher et al., 2013]
k_BC_adjust_map = {
    "Boucher2013": -0.1 / 0.6,
    # variations from [Lohmann et al., 2010] (figure 2; data_loaders provided by author)
    "CSIRO": -0.37 * (-0.1 / 0.6) / -0.111,
    "GISS": -0.225 * (-0.1 / 0.6) / -0.111,
    "HadGEM2": -0.13 * (-0.1 / 0.6) / -0.111,
    "ECHAM5": 0.05 * (-0.1 / 0.6) / -0.111,
    "ECMWF": 0.12 * (-0.1 / 0.6) / -0.111,
}
k_BC_adjust = k_BC_adjust_map[mod_BCadjust]

# solubility of aerosols for the aerosol-cloud interaction
# from [Hansen et al., 2005]
if mod_CLOUDsolub == "Hansen2005":
    solub_SO4 = 1.0
    solub_POA = 0.8
    solub_BC = 0.7  # average of FF (0.6) and BB (0.8)
    solub_NO3 = 1.0
    solub_SOA = 0.8  # assumed
    solub_DUST = 0.0
    solub_SALT = 1.0
# from [Lamarque et al., 2011] (RCP database)
elif mod_CLOUDsolub == "Lamarque2011":
    solub_SO4 = 1.0
    solub_POA = 0.86
    solub_BC = 0.80
    solub_NO3 = 1.0
    solub_SOA = 1.0
    solub_DUST = 0.12
    solub_SALT = 0.05
else:
    raise RuntimeError

# ERF over 1850-2000 for the aerosol-cloud interaction
# from ACCMIP [Shindell et al., 2013] (table 7)
# rescaled to the best guess of IPCC AR5 [Boucher et al., 2013]
RF_ref1_map = {
    "mean-ACCMIP": -0.84 * -0.45 / -0.84,
    "CSIRO-Mk360": -0.99 * -0.45 / -0.84,
    "GFDL-AM3": -0.82 * -0.45 / -0.84,
    "GISS-E2-R": -0.61 * -0.45 / -0.84,
    "HadGEM2": -0.89 * -0.45 / -0.84,
    "LMDzORINCA": -0.21 * -0.45 / -0.84,
    "MIROC-CHEM": -1.12 * -0.45 / -0.84,
    "NCAR-CAM-51": -1.22 * -0.45 / -0.84,
}
RF_ref1 = RF_ref1_map[mod_CLOUDerf]

# soluble aerosol load for the aerosol-cloud interaction
# load pre-processed ACCMIP data_loaders
path = f"data/AeroCloud_ACCMIP/#DATA.AeroCloud_{mod_CLOUDerf}.(2yr)_(7aer).LOAD.csv"
TMP, lgd = load_data_and_header(path)
del lgd[0]
AER_ref0 = 0
AER_ref1 = 0
for n in range(len(lgd)):
    if not np.isnan(TMP[0, n]):
        AER_ref0 += TMP[0, n] * globals()["solub_" + lgd[n]]
        AER_ref1 += TMP[1, n] * globals()["solub_" + lgd[n]]

#: Aerosol-cloud interaction intensity of indirect effect {W/m2}
#: based on a logarithmic formulation [e.g. Gultepe and Isaac, 1999]
Phi_0 = RF_ref1 / np.log(AER_ref1 / AER_ref0)

#: Preindustrial load of soluble aerosols {Tg}
#: reduced by a factor from [Carslaw et al., 2013] and two arbitrary variations
AERh_0_map = {
    "median": AER_ref0 * np.exp(1 * (1.42 - 1.30) / Phi_0),
    "high": AER_ref0 * np.exp(0 * (1.42 - 1.30) / Phi_0),
    "low": AER_ref0 * np.exp(2 * (1.42 - 1.30) / Phi_0),
}
AERh_0 = AERh_0_map[mod_CLOUDpreind]

# ----------------
# 7.3.3. Volcanoes
# ----------------

# warming efficacy of volcano forcing {.}
# based on [Gregory et al., 2016]
warmeff_volc: float = 0.6

# ===========
# 7.4. ALBEDO
# ===========

# -------------------
# 7.4.1. Black Carbon
# -------------------

# read region distribution
path = "data/RegDiv_Reddy2007/#DATA.RegDiv_Reddy2007.114reg1_(9reg0).AREA.csv"
TMP = load_data(path, start=1)
p_reg9 = np.zeros([nb_regionI, 9 + 1], dtype=dty)
for i in range(1, 114 + 1):
    p_reg9[regionI_index[i], :] += TMP[i - 1, :]
p_reg9 /= np.sum(p_reg9, 1)[:, np.newaxis]
p_reg9[np.isnan(p_reg9) | np.isinf(p_reg9)] = 0

# regional effect coefficient {.}
# from [Reddy and Boucher, 2007] (table 1)
if mod_ALBBCreg == "Reddy2007":
    w_reg_BCsnow = np.array([1, 1 / 6.0, 11 / 11.0, 1 / 10.0, 63 / 12.0, 2 / 3.0, 2 / 13.0, 17 / 43.0, 1 / 1.0, 2 / 1.0], dtype=dty)

# global normalized radiative effect {{W/m2}/{Tg/yr}}
# from ACCMIP [Lee et al., 2013] (tab. 3 & fig. 15)
# rescaled to the best guess of IPCC AR5 [Boucher et al., 2013] (using table 7.1a)
radeff_BCsnow_map = {
    "mean-ACCMIP":  0.0146 / (7.9 - 3.2) * (0.04 / 4.8) / (0.0146 / (7.9 - 3.2)),
    "CICERO-OsloCTM2":  0.0131 / (7.8 - 3.1) * (0.04 / 4.8) / (0.0146 / (7.9 - 3.2)),
    "GFDL-AM3":  0.0130 / (7.8 - 3.1) * (0.04 / 4.8) / (0.0146 / (7.9 - 3.2)),
    "GISS-E2-R":  0.0142 / (8.8 - 4.0) * (0.04 / 4.8) / (0.0146 / (7.9 - 3.2)),
    "GISS-E2-R-TOMAS":  0.0175 / (7.8 - 3.1) * (0.04 / 4.8) / (0.0146 / (7.9 - 3.2)),
    "HadGEM2":  0.0133 / (7.8 - 3.1) * (0.04 / 4.8) / (0.0146 / (7.9 - 3.2)),
    "MIROC-CHEM":  0.0173 / (7.7 - 3.0) * (0.04 / 4.8) / (0.0146 / (7.9 - 3.2)),
    "NCAR-CAM-35":  0.0143 / (7.8 - 3.1) * (0.04 / 4.8) / (0.0146 / (7.9 - 3.2)),
    "NCAR-CAM-51":  0.0141 / (7.8 - 3.1) * (0.04 / 4.8) / (0.0146 / (7.9 - 3.2)),
}
radeff_BCsnow = radeff_BCsnow_map[mod_ALBBCrf]

# warming efficacy of black carbon on snow {.}
# from IPCC AR5 [Boucher et al., 2013] (sect. 7.5.2.3)
warmeff_BCsnow_map = {"median": 3.0, "low": 2.0, "high": 4.0}
warmeff_BCsnow = warmeff_BCsnow_map[mod_ALBBCwarm]

# -----------------
# 7.4.2. Land-Cover
# -----------------

# basic biomes of aggregation
bio = ["des", "for", "shr", "gra", "cro", "pas", "urb"]

# load pre-processed albedo climatology
path = f"data/Albedo_{mod_ALBLCalb}/#DATA.Albedo_{mod_ALBLCalb}_{mod_ALBLCflux}_{mod_ALBLCcover}.114reg1_7bio.alb.csv"
path2 = f"data/Albedo_{mod_ALBLCalb}/#DATA.Albedo_{mod_ALBLCalb}_{mod_ALBLCflux}_{mod_ALBLCcover}.114reg1_7bio.RSDS.csv"
TMP = load_data(path)
TMP2 = load_data(path2)
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
        alpha_alb[regionI_index[i], biome_index["pas"]] += (0.6 * TMP[i - 1, bio.index("gra")] * TMP2[i - 1, bio.index("gra")]
                + 0.4 * TMP[i - 1, bio.index("des")] * TMP2[i - 1, bio.index("des")])
        RSDS_alb[regionI_index[i], biome_index["pas"]] += (0.6 * TMP2[i - 1, bio.index("gra")] + 0.4 * TMP2[i - 1, bio.index("des")])
    alpha_alb[:, biome_index["pas"]] /= RSDS_alb[:, biome_index["pas"]]
    alpha_alb[np.isnan(alpha_alb) | np.isinf(alpha_alb)] = 0

# load pre-processed radiation climatology
path = f"data/RadFlux_{mod_ALBLCflux}/#DATA.RadFlux_{mod_ALBLCflux}.114reg1.rsds.csv"
path2 = f"data/RadFlux_{mod_ALBLCflux}/#DATA.RadFlux_{mod_ALBLCflux}.114reg1.AREA.csv"
TMP = load_data(path)
TMP2 = load_data(path2)
rsds_alb = np.zeros([nb_regionI])
AREA_alb = np.zeros([nb_regionI])
for i in range(1, 114 + 1):
    rsds_alb[regionI_index[i]] += TMP[i - 1, 0] * TMP2[i - 1, 0]
    AREA_alb[regionI_index[i]] += TMP2[i - 1, 0]
rsds_alb /= AREA_alb
rsds_alb[np.isnan(rsds_alb) | np.isinf(rsds_alb)] = 0

# upward transmittance from [Lenton and Vaughan, 2009]
p_trans = -0.854

# final albedo parameters {{W/m2}/Mha}
alpha_LCC = p_trans * alpha_alb * rsds_alb[:, np.newaxis] / (510_072e9 / 1e10)

# warming efficacy of land-cover change albedo effect {.}
# from [Bright et al., 2015] (tab. 7)
warmeff_LCC_map = {"Hansen2005": 1.02, "Davin2007": 0.5, "Davin2010": 0.78, "Jones2013": 0.79}
warmeff_LCC = warmeff_LCC_map[mod_ALBLCwarm]