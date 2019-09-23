import numpy as np

from .carbon import CO2_0, CO2_cmip5
from ..data import load_data, load_header
from ..data_loaders import nb_regionI, regionI_index, nb_biome, biome_index
from .. import historical
from .. import conf

CH4_0 = historical.CH4_0

##################################################
#   2. METHANE
##################################################

# ===============
# 2.1. ATMOSPHERE
# ===============

#: Conversion of CH4 from {ppb} to {TgC}
alpha_CH4: float = 0.1765 * 12.0

# ==============
# 2.2. CHEMISTRY
# ==============

# ----------------
# 2.2.1. Lifetimes
# ----------------

#: Preindustrial lifetime of CH4 {yr} [Prather et al., 2012]
#: Rescale of strato lifetime from [Prather et al., 2015]
tau_CH4_hv: float = 1.06 * 120.0
tau_CH4_soil = 150.0
tau_CH4_ocean = 200.0
tau_CH4_OH_map = {# Average from [Prather et al., 2012]
    "Prather2012": 11.2,

    # Variations from [Naik et al., 2013] (table 1)
    "CESM-CAM-superfast": 8.4 * 11.2 / 9.7, "CICERO-OsloCTM2": 10.0 * 11.2 / 9.7, "CMAM": 9.4 * 11.2 / 9.7, "EMAC": 9.1 * 11.2 / 9.7, "GEOSCCM": 9.6 * 11.2 / 9.7, "GDFL-AM3": 9.4 * 11.2 / 9.7, "GISS-E2-R": 10.6 * 11.2 / 9.7,
    "GISS-E2-R-TOMAS": 9.2 * 11.2 / 9.7, "HadGEM2": 11.6 * 11.2 / 9.7, "LMDzORINCA": 10.5 * 11.2 / 9.7, "MIROC-CHEM": 8.7 * 11.2 / 9.7, "MOCAGE": 7.1 * 11.2 / 9.7, "NCAR-CAM-35": 9.2 * 11.2 / 9.7, "STOC-HadAM3": 9.1 * 11.2 / 9.7,
    "TM5": 9.9 * 11.2 / 9.7, "UM-CAM": 14.0 * 11.2 / 9.7, }
tau_CH4_OH = tau_CH4_OH_map[conf.mod_OHSNKtau]

# arbitrary rescaling (down) of OH lifetime
scale_OH: float = 0.80
tau_CH4_OH *= scale_OH

# -----------------
# 2.2.2. Holmes2013
# -----------------

#: Sensitivity of OH sink to climate, ozone and methane {.} [Holmes et al., 2013] and [Ehhalt et al., 2001] (TAR, mean-OxComp)
chi_OH_Tatm_map = {"Holmes2013": 3.0, "UCI-CTM": 3.9, "Oslo-CTM3": 2.8, "GEOS-Chem": 2.2, "mean-OxComp": 0.0}
chi_OH_Qatm_map = {"Holmes2013": 0.32, "UCI-CTM": 0.32, "Oslo-CTM3": 0.29, "GEOS-Chem": 0.34, "mean-OxComp": 0.0}
chi_OH_O3_map = {"Holmes2013": -0.55, "UCI-CTM": -0.66, "Oslo-CTM3": -0.43, "GEOS-Chem": -0.61, "mean-OxComp": 0.0}
chi_OH_CH4_map = {"Holmes2013": -0.31, "UCI-CTM": -0.363, "Oslo-CTM3": -0.307, "GEOS-Chem": -0.274, "mean-OxComp": -0.32}
chi_OH_Tatm = chi_OH_Tatm_map[conf.mod_OHSNKtrans]
chi_OH_Qatm = chi_OH_Qatm_map[conf.mod_OHSNKtrans]
chi_OH_O3 = chi_OH_O3_map[conf.mod_OHSNKtrans]
chi_OH_CH4 = chi_OH_CH4_map[conf.mod_OHSNKtrans]

#: Sensitivity of OH sink to ozone precursors {.}; from [Holmes et al., 2013]
chi_OH_map = {# Logarithmic sensitivity, kept logarithmic
    ('log', 'NOX'): {"Holmes2013": -0.14, "UCI-CTM": -0.15, "Oslo-CTM3": -0.10, "GEOS-Chem": -0.16, "mean-OxComp": -0.137}, ('log', 'CO'): {"Holmes2013": 0.06, "UCI-CTM": 0.066, "Oslo-CTM3": 0.05, "GEOS-Chem": 0.065, "mean-OxComp": 0.11},
    ('log', 'VOC'): {"Holmes2013": 0.04, "UCI-CTM": 0.04, "Oslo-CTM3": 0.04, "GEOS-Chem": 0.04, "mean-OxComp": 0.047},

    # Linear sensitivity {./{TgN/yr}}&{./{TgC/yr}}&{./{Tg/yr}}, rescaled using the TAR [Ehhalt et al., 2001]
    ('lin', 'NOX'): {"Holmes2013": 0.0042 * -0.14 / -0.137, "UCI-CTM": 0.0042 * -0.15 / -0.137, "Oslo-CTM3": 0.0042 * -0.10 / -0.137, "GEOS-Chem": 0.0042 * -0.16 / -0.137, "mean-OxComp": 0.0042},
    ('lin', 'CO'): {"Holmes2013": -1.05e-4 * 0.06 / 0.11 * 28 / 12.0, "UCI-CTM": -1.05e-4 * 0.066 / 0.11 * 28 / 12.0, "Oslo-CTM3": -1.05e-4 * 0.05 / 0.11 * 28 / 12.0, "GEOS-Chem": -1.05e-4 * 0.065 / 0.11 * 28 / 12.0,
                    "mean-OxComp": -1.05e-4 * 28 / 12.0},
    ('lin', 'VOC'): {"Holmes2013": -3.14e-4 * 0.04 / 0.047, "UCI-CTM": -3.14e-4 * 0.04 / 0.047, "Oslo-CTM3": -3.14e-4 * 0.04 / 0.047, "GEOS-Chem": -3.14e-4 * 0.04 / 0.047, "mean-OxComp": -3.14e-4}, }

chi_OH_NOX = chi_OH_map[conf.mod_OHSNKfct, 'NOX'][conf.mod_OHSNKtrans]
chi_OH_CO = chi_OH_map[conf.mod_OHSNKfct, 'CO'][conf.mod_OHSNKtrans]
chi_OH_VOC = chi_OH_map[conf.mod_OHSNKfct, 'VOC'][conf.mod_OHSNKtrans]

#: Constant parameters for OH sink function {.}&{.}&{K}&{DU} [Holmes et al., 2013] (supp.)
k_Tatm: float = 0.94  # +/- 0.1
k_Qatm = 1.5  # +/- 0.1
Tatm_0 = 251.0  # +/- 1
# from [Jacobson, 2005] (eq. 2.62)
k_svp: float = 17.67
T_svp = 243.5 - 273.15
# from [Cionni et al., 2010]
O3s_0: float = 280.0
# from [Skeie et al., 2011] (tab. 1)
ENOX_oh: float = 13.0
ECO_oh = 180.0 * 12 / 28.0
EVOC_oh = 39.0 + 220.0 + 175.0

# expression of OH sink function {.}
# from [Holmes et al., 2013] (supp.)
if conf.mod_OHSNKfct == "log":

    def f_kOH(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(1 + D_O3s / O3s_0) + chi_OH_NOX * np.log(1 + ENOX / ENOX_oh) + chi_OH_CO * np.log(1 + ECO / ECO_oh) + chi_OH_VOC * np.log(1 + EVOC / EVOC_oh) + chi_OH_Tatm * np.log(
                1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))) - 1)
        return np.array(D_kOH, dtype=conf.dty)


    def df_kOH_dCH4(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (chi_OH_CH4 / (CH4_0 + D_CH4) * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(1 + D_O3s / O3s_0) + chi_OH_NOX * np.log(1 + ENOX / ENOX_oh) + chi_OH_CO * np.log(1 + ECO / ECO_oh) + chi_OH_VOC * np.log(1 + EVOC / EVOC_oh) + chi_OH_Tatm * np.log(
                1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))))
        return np.array(D_kOH, dtype=conf.dty)


    def df_kOH_dO3s(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (chi_OH_O3 / (O3s_0 + D_O3s) * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(1 + D_O3s / O3s_0) + chi_OH_NOX * np.log(1 + ENOX / ENOX_oh) + chi_OH_CO * np.log(1 + ECO / ECO_oh) + chi_OH_VOC * np.log(1 + EVOC / EVOC_oh) + chi_OH_Tatm * np.log(
                1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))))
        return np.array(D_kOH, dtype=conf.dty)


    def df_kOH_dgst(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (chi_OH_Tatm * k_Qatm * k_svp * k_Tatm / (Tatm_0 + T_svp) * np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) / (1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1)) * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(1 + D_O3s / O3s_0) + chi_OH_NOX * np.log(1 + ENOX / ENOX_oh) + chi_OH_CO * np.log(1 + ECO / ECO_oh) + chi_OH_VOC * np.log(1 + EVOC / EVOC_oh) + chi_OH_Tatm * np.log(
                1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))))
        return np.array(D_kOH, dtype=conf.dty)


    def df_kOH_dENOX(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (chi_OH_NOX / (ENOX_oh + ENOX) * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(1 + D_O3s / O3s_0) + chi_OH_NOX * np.log(1 + ENOX / ENOX_oh) + chi_OH_CO * np.log(1 + ECO / ECO_oh) + chi_OH_VOC * np.log(1 + EVOC / EVOC_oh) + chi_OH_Tatm * np.log(
                1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))))
        return np.array(D_kOH, dtype=conf.dty)


    def df_kOH_dECO(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (chi_OH_CO / (ECO_oh + ECO) * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(1 + D_O3s / O3s_0) + chi_OH_NOX * np.log(1 + ENOX / ENOX_oh) + chi_OH_CO * np.log(1 + ECO / ECO_oh) + chi_OH_VOC * np.log(1 + EVOC / EVOC_oh) + chi_OH_Tatm * np.log(
                1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))))
        return np.array(D_kOH, dtype=conf.dty)


    def df_kOH_dEVOC(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (chi_OH_VOC / (EVOC_oh + EVOC) * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(1 + D_O3s / O3s_0) + chi_OH_NOX * np.log(1 + ENOX / ENOX_oh) + chi_OH_CO * np.log(1 + ECO / ECO_oh) + chi_OH_VOC * np.log(1 + EVOC / EVOC_oh) + chi_OH_Tatm * np.log(
                1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))))
        return np.array(D_kOH, dtype=conf.dty)


    def df_kOH(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        df_1 = df_kOH_dCH4(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * D_CH4
        df_2 = df_kOH_dO3s(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * D_O3s
        df_3 = df_kOH_dgst(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * D_gst
        df_4 = df_kOH_dENOX(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * ENOX
        df_5 = df_kOH_dECO(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * ECO
        df_6 = df_kOH_dEVOC(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * EVOC
        df_tot = df_1 + df_2 + df_3 + df_4 + df_5 + df_6
        return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot), np.nan_to_num(df_3 / df_tot), np.nan_to_num(df_4 / df_tot), np.nan_to_num(df_5 / df_tot), np.nan_to_num(df_6 / df_tot)]


elif conf.mod_OHSNKfct == "lin":

    def f_kOH(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (np.exp(chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(1 + D_O3s / O3s_0) + chi_OH_NOX * ENOX + chi_OH_CO * ECO + chi_OH_VOC * EVOC + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
            1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))) - 1)
        return np.array(D_kOH, dtype=conf.dty)


    def df_kOH_dCH4(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (chi_OH_CH4 / (CH4_0 + D_CH4) * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(1 + D_O3s / O3s_0) + chi_OH_NOX * ENOX + chi_OH_CO * ECO + chi_OH_VOC * EVOC + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))))
        return np.array(D_kOH, dtype=conf.dty)


    def df_kOH_dO3s(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (chi_OH_O3 / (O3s_0 + D_O3s) * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(1 + D_O3s / O3s_0) + chi_OH_NOX * ENOX + chi_OH_CO * ECO + chi_OH_VOC * EVOC + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))))
        return np.array(D_kOH, dtype=conf.dty)


    def df_kOH_dgst(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (chi_OH_Tatm * k_Qatm * k_svp * k_Tatm / (Tatm_0 + T_svp) * np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) / (1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1)) * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(1 + D_O3s / O3s_0) + chi_OH_NOX * ENOX + chi_OH_CO * ECO + chi_OH_VOC * EVOC + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))))
        return np.array(D_kOH, dtype=conf.dty)


    def df_kOH_dENOX(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = chi_OH_NOX * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(1 + D_O3s / O3s_0) + chi_OH_NOX * ENOX + chi_OH_CO * ECO + chi_OH_VOC * EVOC + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1)))
        return np.array(D_kOH, dtype=conf.dty)


    def df_kOH_dECO(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = chi_OH_CO * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(1 + D_O3s / O3s_0) + chi_OH_NOX * ENOX + chi_OH_CO * ECO + chi_OH_VOC * EVOC + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1)))
        return np.array(D_kOH, dtype=conf.dty)


    def df_kOH_dEVOC(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = chi_OH_VOC * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(1 + D_O3s / O3s_0) + chi_OH_NOX * ENOX + chi_OH_CO * ECO + chi_OH_VOC * EVOC + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1)))
        return np.array(D_kOH, dtype=conf.dty)


    def df_kOH(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        df_1 = df_kOH_dCH4(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * D_CH4
        df_2 = df_kOH_dO3s(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * D_O3s
        df_3 = df_kOH_dgst(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * D_gst
        df_4 = df_kOH_dENOX(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * ENOX
        df_5 = df_kOH_dECO(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * ECO
        df_6 = df_kOH_dEVOC(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * EVOC
        df_tot = df_1 + df_2 + df_3 + df_4 + df_5 + df_6
        return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot), np.nan_to_num(df_3 / df_tot), np.nan_to_num(df_4 / df_tot), np.nan_to_num(df_5 / df_tot), np.nan_to_num(df_6 / df_tot)]

# =========
# 2.3. LAND
# =========

# ---------------
# 2.3.1. WETCHIMP
# ---------------

# basic biomes of aggregation
bio = ["des", "for", "shr", "gra", "cro", "pas"]

# load pre-processed WETCHIMP results for specified model
# preindustrial partition of wetlands
conf.mod_LSNKcover_save = conf.mod_LSNKcover
cover_options = ["ESA-CCI", "MODIS", "Ramankutty1999", "Levavasseur2012", "mean-TRENDYv2", "CLM-45", "JSBACH", "JULES", "LPJ", "LPJ-GUESS", "LPX-Bern", "OCN", "ORCHIDEE", "VISIT"]
AREA_wet = np.zeros([nb_regionI, nb_biome], dtype=conf.dty)
if conf.mod_EWETpreind != "":
    for conf.mod_LSNKcover in cover_options:
        path = f"data/Wetlands_WETCHIMP/#DATA.Wetlands_{conf.mod_EWETpreind}_{conf.mod_LSNKcover}.1910s_114reg1_4bio.AREA.csv"
        TMP = load_data(path)
        for i in range(1, 114 + 1):
            for b in range(len(bio) - 2):
                AREA_wet[regionI_index[i], biome_index[bio[b]]] += TMP[i - 1, b]

conf.mod_LSNKcover = conf.mod_LSNKcover_save
del conf.mod_LSNKcover_save

# preindustrial area and emissions
AREA_wet0 = np.zeros([nb_regionI], dtype=conf.dty)
ECH4_wet0 = np.zeros([nb_regionI], dtype=conf.dty)
for VAR, arr in [("AREA", AREA_wet0), ("ECH4", ECH4_wet0)]:
    if conf.mod_EWETpreind != "":
        TMP = load_data(f"data/Wetlands_WETCHIMP/#DATA.Wetlands_{conf.mod_EWETpreind}.1910s_114reg1_(1exp).{VAR}.csv", start=1)
        for i in range(1, 114 + 1):
            arr[regionI_index[i]] += TMP[i - 1, 0]

# calculate preindustrial parameters {TgC/Mha}&{Mha}&{.}
# with arbitrary adjustment of areal emissions
ewet_0 = ECH4_wet0 / AREA_wet0 * CO2_0 / np.mean(CO2_cmip5[201:231])
ewet_0[np.isnan(ewet_0) | np.isinf(ewet_0)] = 0
AWET_0 = AREA_wet0.copy()
p_wet = AREA_wet / np.sum(AREA_wet, 1)[:, np.newaxis]

# ensure NaN and zeros removed
ewet_0[np.isnan(ewet_0) | np.isinf(ewet_0)] = 0
p_wet[np.isnan(p_wet) | np.isinf(p_wet)] = 0
p_wet[p_wet == 0] = 1e-18
p_wet /= np.sum(p_wet, 1)[:, np.newaxis]


# load pre-processed WETCHIMP results for specified model
# response to perturbation experiments
def wetchimp(sim):
    arr = np.zeros([nb_regionI], dtype=conf.dty)
    if conf.mod_AWETtrans != "":
        TMP = load_data(f"data/Wetlands_WETCHIMP/#DATA.Wetlands_{conf.mod_AWETtrans}.1910s_114reg1_(4exp).AREA.csv", start=1)
        path = f"data/Wetlands_WETCHIMP/#DATA.Wetlands_{conf.mod_AWETtrans}.1910s_114reg1_(4exp).AREA.csv"
        lgd = load_header(path)
        for i in range(1, 114 + 1):
            arr[regionI_index[i]] += TMP[i - 1, lgd.index(sim)]
    return arr


AREA_exp0 = wetchimp("exp0")
AREA_expC = wetchimp("expC")
AREA_expT = wetchimp("expT")
AREA_expP = wetchimp("expP")

# value of perturbations
# see [Melton et al., 2013]
CO2_expC: float = 857.0 - 303.0
T_expT = 3.4
P_expP = np.zeros([nb_regionI], dtype=conf.dty)
S_tmp = np.zeros([nb_regionI], dtype=conf.dty)
TMP = load_data("data/HistLand_CRU/#DATA.HistLand_CRU.1901-2014_114reg1.lyp.csv")
TMP2 = load_data("data/HistLand_CRU/#DATA.HistLand_CRU.1901-2014_114reg1.AREA.csv")
for i in range(1, 114 + 1):
    P_expP[regionI_index[i]] += np.mean(TMP[:30, i - 1] * TMP2[:30, i - 1], 0)
    S_tmp[regionI_index[i]] += np.mean(TMP2[:30, i - 1], 0)
P_expP /= S_tmp
P_expP[np.isnan(P_expP) | np.isinf(P_expP)] = 0
P_expP *= 0.039

# calculation of parameters
# sensitivity of wetland extent to CO2 {./ppm}
gamma_wetC = (AREA_expC / AREA_exp0 - 1) / CO2_expC
# sensitivity of wetland extent to temperature {./K}
gamma_wetT = (AREA_expT / AREA_exp0 - 1) / T_expT
# sensitivity of wetland extent to precipitation {./mm}
gamma_wetP = (AREA_expP / AREA_exp0 - 1) / P_expP
# [NaN]
gamma_wetT[np.isnan(gamma_wetT) | np.isinf(gamma_wetT)] = 0
gamma_wetP[np.isnan(gamma_wetP) | np.isinf(gamma_wetP)] = 0
gamma_wetC[np.isnan(gamma_wetC) | np.isinf(gamma_wetC)] = 0

# -----------------
# 2.3.2. Permafrost
# -----------------

# fraction of methane in for PF emissions {.}
# best guess from [Schuur et al., 2015]
p_PF_CH4_map = {"zero": 0.000, "best": 0.023, "twice": 0.046}
p_PF_CH4 = p_PF_CH4_map[conf.mod_EPFmethane]

# fraction of instantaneous PF emission {.}
p_PF_inst: float = 0.0
