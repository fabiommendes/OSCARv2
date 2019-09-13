from .a1_carbon import CO2_0, CO2_cmip5, load_data
from ..oscar_data import *
from ..oscar_data import nb_regionI, regionI_index, nb_biome, biome_index
from ...config import dty, mod_LSNKcover, mod_OHSNKtau, mod_OHSNKtrans, mod_OHSNKfct, \
    mod_EWETpreind, mod_AWETtrans, mod_EPFmethane

##################################################
#   2. METHANE
##################################################

# ===============
# 2.1. ATMOSPHERE
# ===============

# conversion of CH4 from {ppb} to {TgC}
alpha_CH4 = 0.1765 * np.array([12.0], dtype=dty)

# historic CH4 from IPCC-AR5 {ppb}
# from [IPCC WG1, 2013] annexe 2
CH4_ipcc = np.ones([311 + 1], dtype=dty) * np.nan
TMP = load_data("data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1750-2011.CH4.csv")
CH4_ipcc[50:] = TMP[:, 0]
CH4_0 = np.array([CH4_ipcc[50]], dtype=dty)

# historic CH4 from CMIP5 {ppb}
# from [Meinshausen et al., 2011]
CH4_cmip5 = np.ones([305 + 1], dtype=dty) * np.nan
TMP = load_data("data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005.CH4.csv")
CH4_cmip5[65:] = TMP[:, 0]

# historic CH4 from AGAGE {ppb}
# from [Prinn et al., 2013] updated from the website
CH4_agage = np.ones([313 + 1], dtype=dty) * np.nan
TMP = load_data("data/HistAtmo_AGAGE/#DATA.HistAtmo_AGAGE.1987-2013.CH4_global.csv")
CH4_agage[287:] = TMP[:, 0]

# historic CH4 from Law Dome ice cores {ppm}
# from [Etheridge et al., 1998] and [MacFarling Meure et al., 2006]
CH4_lawdome = load_data("data/HistAtmo_NOAA-NCDC/#DATA.HistAtmo_NOAA-NCDC.(IceCores).CH4_lawdome.csv", slice=1)

# load RCP concentrations {ppb}
# from [Meinshausen et al., 2011]
CH4_rcp = np.ones([800 + 1, 6], dtype=dty) * np.nan
n = -1
for rcp in ["rcp26", "rcp45", "rcp60", "rcp85", "rcp45to26", "rcp60to45"]:
    n += 1
    TMP = load_data(f"data/Scenario_ECP/#DATA.Scenario_ECP.2000-2500.{rcp}_CH4.csv")
    CH4_rcp[300:, n] = TMP[:, 0]

# ==============
# 2.2. CHEMISTRY
# ==============

# ----------------
# 2.2.1. Lifetimes
# ----------------

# preindustrial lifetime of CH4 {yr}
# from [Prather et al., 2012]
# rescale of strato lifetime from [Prather et al., 2015]
tau_CH4_hv = 1.06 * np.array([120.0], dtype=dty)
tau_CH4_soil = np.array([150.0], dtype=dty)
tau_CH4_ocean = np.array([200.0], dtype=dty)
# average from [Prather et al., 2012]
if mod_OHSNKtau == "Prather2012":
    tau_CH4_OH = np.array([11.2], dtype=dty)
# variations from [Naik et al., 2013] (table 1)
elif mod_OHSNKtau == "CESM-CAM-superfast":
    tau_CH4_OH = np.array([8.4], dtype=dty) * 11.2 / 9.7
elif mod_OHSNKtau == "CICERO-OsloCTM2":
    tau_CH4_OH = np.array([10.0], dtype=dty) * 11.2 / 9.7
elif mod_OHSNKtau == "CMAM":
    tau_CH4_OH = np.array([9.4], dtype=dty) * 11.2 / 9.7
elif mod_OHSNKtau == "EMAC":
    tau_CH4_OH = np.array([9.1], dtype=dty) * 11.2 / 9.7
elif mod_OHSNKtau == "GEOSCCM":
    tau_CH4_OH = np.array([9.6], dtype=dty) * 11.2 / 9.7
elif mod_OHSNKtau == "GDFL-AM3":
    tau_CH4_OH = np.array([9.4], dtype=dty) * 11.2 / 9.7
elif mod_OHSNKtau == "GISS-E2-R":
    tau_CH4_OH = np.array([10.6], dtype=dty) * 11.2 / 9.7
elif mod_OHSNKtau == "GISS-E2-R-TOMAS":
    tau_CH4_OH = np.array([9.2], dtype=dty) * 11.2 / 9.7
elif mod_OHSNKtau == "HadGEM2":
    tau_CH4_OH = np.array([11.6], dtype=dty) * 11.2 / 9.7
elif mod_OHSNKtau == "LMDzORINCA":
    tau_CH4_OH = np.array([10.5], dtype=dty) * 11.2 / 9.7
elif mod_OHSNKtau == "MIROC-CHEM":
    tau_CH4_OH = np.array([8.7], dtype=dty) * 11.2 / 9.7
elif mod_OHSNKtau == "MOCAGE":
    tau_CH4_OH = np.array([7.1], dtype=dty) * 11.2 / 9.7
elif mod_OHSNKtau == "NCAR-CAM-35":
    tau_CH4_OH = np.array([9.2], dtype=dty) * 11.2 / 9.7
elif mod_OHSNKtau == "STOC-HadAM3":
    tau_CH4_OH = np.array([9.1], dtype=dty) * 11.2 / 9.7
elif mod_OHSNKtau == "TM5":
    tau_CH4_OH = np.array([9.9], dtype=dty) * 11.2 / 9.7
elif mod_OHSNKtau == "UM-CAM":
    tau_CH4_OH = np.array([14.0], dtype=dty) * 11.2 / 9.7

# arbitrary rescaling (down) of OH lifetime
scale_OH = 0.80
tau_CH4_OH *= scale_OH

# -----------------
# 2.2.2. Holmes2013
# -----------------

# sensitivity of OH sink to climate, ozone and methane {.}
# from [Holmes et al., 2013]
if mod_OHSNKtrans == "Holmes2013":
    chi_OH_Tatm = np.array([3.0], dtype=dty)
    chi_OH_Qatm = np.array([0.32], dtype=dty)
    chi_OH_O3 = np.array([-0.55], dtype=dty)
    chi_OH_CH4 = np.array([-0.31], dtype=dty)
elif mod_OHSNKtrans == "UCI-CTM":
    chi_OH_Tatm = np.array([3.9], dtype=dty)
    chi_OH_Qatm = np.array([0.32], dtype=dty)
    chi_OH_O3 = np.array([-0.66], dtype=dty)
    chi_OH_CH4 = np.array([-0.363], dtype=dty)
elif mod_OHSNKtrans == "Oslo-CTM3":
    chi_OH_Tatm = np.array([2.8], dtype=dty)
    chi_OH_Qatm = np.array([0.29], dtype=dty)
    chi_OH_O3 = np.array([-0.43], dtype=dty)
    chi_OH_CH4 = np.array([-0.307], dtype=dty)
elif mod_OHSNKtrans == "GEOS-Chem":
    chi_OH_Tatm = np.array([2.2], dtype=dty)
    chi_OH_Qatm = np.array([0.34], dtype=dty)
    chi_OH_O3 = np.array([-0.61], dtype=dty)
    chi_OH_CH4 = np.array([-0.274], dtype=dty)
# but also from the TAR [Ehhalt et al., 2001]
elif mod_OHSNKtrans == "mean-OxComp":
    chi_OH_Tatm = np.array([0.0], dtype=dty)
    chi_OH_Qatm = np.array([0.0], dtype=dty)
    chi_OH_O3 = np.array([0.0], dtype=dty)
    chi_OH_CH4 = np.array([-0.32], dtype=dty)

# logarithmic sensitivity of OH sink to ozone precursors {.}
# kept logarithmic; from [Holmes et al., 2013]
if mod_OHSNKfct == "log":
    if mod_OHSNKtrans == "Holmes2013":
        chi_OH_NOX = np.array([-0.14], dtype=dty)
        chi_OH_CO = np.array([0.06], dtype=dty)
        chi_OH_VOC = np.array([0.04], dtype=dty)
    elif mod_OHSNKtrans == "UCI-CTM":
        chi_OH_NOX = np.array([-0.15], dtype=dty)
        chi_OH_CO = np.array([0.066], dtype=dty)
        chi_OH_VOC = np.array([0.04], dtype=dty)
    elif mod_OHSNKtrans == "Oslo-CTM3":
        chi_OH_NOX = np.array([-0.10], dtype=dty)
        chi_OH_CO = np.array([0.05], dtype=dty)
        chi_OH_VOC = np.array([0.04], dtype=dty)
    elif mod_OHSNKtrans == "GEOS-Chem":
        chi_OH_NOX = np.array([-0.16], dtype=dty)
        chi_OH_CO = np.array([0.065], dtype=dty)
        chi_OH_VOC = np.array([0.04], dtype=dty)
    elif mod_OHSNKtrans == "mean-OxComp":
        chi_OH_NOX = np.array([-0.137], dtype=dty)
        chi_OH_CO = np.array([0.11], dtype=dty)
        chi_OH_VOC = np.array([0.047], dtype=dty)

    # linear sensitivity {./{TgN/yr}}&{./{TgC/yr}}&{./{Tg/yr}}
# rescaled using the TAR [Ehhalt et al., 2001]
elif mod_OHSNKfct == "lin":
    if mod_OHSNKtrans == "Holmes2013":
        chi_OH_NOX = np.array([0.0042 * -0.14 / -0.137], dtype=dty)
        chi_OH_CO = np.array([-1.05e-4 * 0.06 / 0.11], dtype=dty) * 28 / 12.0
        chi_OH_VOC = np.array([-3.14e-4 * 0.04 / 0.047], dtype=dty)
    elif mod_OHSNKtrans == "UCI-CTM":
        chi_OH_NOX = np.array([0.0042 * -0.15 / -0.137], dtype=dty)
        chi_OH_CO = np.array([-1.05e-4 * 0.066 / 0.11], dtype=dty) * 28 / 12.0
        chi_OH_VOC = np.array([-3.14e-4 * 0.04 / 0.047], dtype=dty)
    elif mod_OHSNKtrans == "Oslo-CTM3":
        chi_OH_NOX = np.array([0.0042 * -0.10 / -0.137], dtype=dty)
        chi_OH_CO = np.array([-1.05e-4 * 0.05 / 0.11], dtype=dty) * 28 / 12.0
        chi_OH_VOC = np.array([-3.14e-4 * 0.04 / 0.047], dtype=dty)
    elif mod_OHSNKtrans == "GEOS-Chem":
        chi_OH_NOX = np.array([0.0042 * -0.16 / -0.137], dtype=dty)
        chi_OH_CO = np.array([-1.05e-4 * 0.065 / 0.11], dtype=dty) * 28 / 12.0
        chi_OH_VOC = np.array([-3.14e-4 * 0.04 / 0.047], dtype=dty)
    elif mod_OHSNKtrans == "mean-OxComp":
        chi_OH_NOX = np.array([0.0042], dtype=dty)
        chi_OH_CO = np.array([-1.05e-4], dtype=dty) * 28 / 12.0
        chi_OH_VOC = np.array([-3.14e-4], dtype=dty)

    # constant parameters for OH sink function {.}&{.}&{K}&{DU}
# from [Holmes et al., 2013] (supp.)
k_Tatm = 0.94  # +/- 0.1
k_Qatm = 1.5  # +/- 0.1
Tatm_0 = 251.0  # +/- 1
# from [Jacobson, 2005] (eq. 2.62)
k_svp = 17.67
T_svp = 243.5 - 273.15
# from [Cionni et al., 2010]
O3s_0 = 280.0
# from [Skeie et al., 2011] (tab. 1)
ENOX_oh = 13.0
ECO_oh = 180.0 * 12 / 28.0
EVOC_oh = 39.0 + 220.0 + 175.0

# expression of OH sink function {.}
# from [Holmes et al., 2013] (supp.)
if mod_OHSNKfct == "log":

    def f_kOH(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (
                np.exp(
                    chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0)
                    + chi_OH_O3 * np.log(1 + D_O3s / O3s_0)
                    + chi_OH_NOX * np.log(1 + ENOX / ENOX_oh)
                    + chi_OH_CO * np.log(1 + ECO / ECO_oh)
                    + chi_OH_VOC * np.log(1 + EVOC / EVOC_oh)
                    + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0)
                    + chi_OH_Qatm * np.log(1 + k_Qatm * (
                            np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))
                )
                - 1
        )
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dCH4(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (
                chi_OH_CH4
                / (CH4_0 + D_CH4)
                * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0)
            + chi_OH_O3 * np.log(1 + D_O3s / O3s_0)
            + chi_OH_NOX * np.log(1 + ENOX / ENOX_oh)
            + chi_OH_CO * np.log(1 + ECO / ECO_oh)
            + chi_OH_VOC * np.log(1 + EVOC / EVOC_oh)
            + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0)
            + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))
        )
        )
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dO3s(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (
                chi_OH_O3
                / (O3s_0 + D_O3s)
                * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0)
            + chi_OH_O3 * np.log(1 + D_O3s / O3s_0)
            + chi_OH_NOX * np.log(1 + ENOX / ENOX_oh)
            + chi_OH_CO * np.log(1 + ECO / ECO_oh)
            + chi_OH_VOC * np.log(1 + EVOC / EVOC_oh)
            + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0)
            + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))
        )
        )
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dgst(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (
                chi_OH_Tatm
                * k_Qatm
                * k_svp
                * k_Tatm
                / (Tatm_0 + T_svp)
                * np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp))
                / (1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))
                * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0)
            + chi_OH_O3 * np.log(1 + D_O3s / O3s_0)
            + chi_OH_NOX * np.log(1 + ENOX / ENOX_oh)
            + chi_OH_CO * np.log(1 + ECO / ECO_oh)
            + chi_OH_VOC * np.log(1 + EVOC / EVOC_oh)
            + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0)
            + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))
        )
        )
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dENOX(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (
                chi_OH_NOX
                / (ENOX_oh + ENOX)
                * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0)
            + chi_OH_O3 * np.log(1 + D_O3s / O3s_0)
            + chi_OH_NOX * np.log(1 + ENOX / ENOX_oh)
            + chi_OH_CO * np.log(1 + ECO / ECO_oh)
            + chi_OH_VOC * np.log(1 + EVOC / EVOC_oh)
            + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0)
            + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))
        )
        )
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dECO(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (
                chi_OH_CO
                / (ECO_oh + ECO)
                * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0)
            + chi_OH_O3 * np.log(1 + D_O3s / O3s_0)
            + chi_OH_NOX * np.log(1 + ENOX / ENOX_oh)
            + chi_OH_CO * np.log(1 + ECO / ECO_oh)
            + chi_OH_VOC * np.log(1 + EVOC / EVOC_oh)
            + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0)
            + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))
        )
        )
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dEVOC(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (
                chi_OH_VOC
                / (EVOC_oh + EVOC)
                * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0)
            + chi_OH_O3 * np.log(1 + D_O3s / O3s_0)
            + chi_OH_NOX * np.log(1 + ENOX / ENOX_oh)
            + chi_OH_CO * np.log(1 + ECO / ECO_oh)
            + chi_OH_VOC * np.log(1 + EVOC / EVOC_oh)
            + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0)
            + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))
        )
        )
        return np.array(D_kOH, dtype=dty)


    def df_kOH(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        df_1 = df_kOH_dCH4(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * D_CH4
        df_2 = df_kOH_dO3s(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * D_O3s
        df_3 = df_kOH_dgst(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * D_gst
        df_4 = df_kOH_dENOX(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * ENOX
        df_5 = df_kOH_dECO(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * ECO
        df_6 = df_kOH_dEVOC(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * EVOC
        df_tot = df_1 + df_2 + df_3 + df_4 + df_5 + df_6
        return [
            np.nan_to_num(df_1 / df_tot),
            np.nan_to_num(df_2 / df_tot),
            np.nan_to_num(df_3 / df_tot),
            np.nan_to_num(df_4 / df_tot),
            np.nan_to_num(df_5 / df_tot),
            np.nan_to_num(df_6 / df_tot), ]


elif mod_OHSNKfct == "lin":

    def f_kOH(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (
                np.exp(
                    chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0)
                    + chi_OH_O3 * np.log(1 + D_O3s / O3s_0)
                    + chi_OH_NOX * ENOX
                    + chi_OH_CO * ECO
                    + chi_OH_VOC * EVOC
                    + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0)
                    + chi_OH_Qatm * np.log(1 + k_Qatm * (
                            np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))
                )
                - 1
        )
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dCH4(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (
                chi_OH_CH4
                / (CH4_0 + D_CH4)
                * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0)
            + chi_OH_O3 * np.log(1 + D_O3s / O3s_0)
            + chi_OH_NOX * ENOX
            + chi_OH_CO * ECO
            + chi_OH_VOC * EVOC
            + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0)
            + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))
        )
        )
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dO3s(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (
                chi_OH_O3
                / (O3s_0 + D_O3s)
                * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0)
            + chi_OH_O3 * np.log(1 + D_O3s / O3s_0)
            + chi_OH_NOX * ENOX
            + chi_OH_CO * ECO
            + chi_OH_VOC * EVOC
            + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0)
            + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))
        )
        )
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dgst(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = (
                chi_OH_Tatm
                * k_Qatm
                * k_svp
                * k_Tatm
                / (Tatm_0 + T_svp)
                * np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp))
                / (1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))
                * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0)
            + chi_OH_O3 * np.log(1 + D_O3s / O3s_0)
            + chi_OH_NOX * ENOX
            + chi_OH_CO * ECO
            + chi_OH_VOC * EVOC
            + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0)
            + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))
        )
        )
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dENOX(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = chi_OH_NOX * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0)
            + chi_OH_O3 * np.log(1 + D_O3s / O3s_0)
            + chi_OH_NOX * ENOX
            + chi_OH_CO * ECO
            + chi_OH_VOC * EVOC
            + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0)
            + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))
        )
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dECO(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = chi_OH_CO * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0)
            + chi_OH_O3 * np.log(1 + D_O3s / O3s_0)
            + chi_OH_NOX * ENOX
            + chi_OH_CO * ECO
            + chi_OH_VOC * EVOC
            + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0)
            + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))
        )
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dEVOC(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = chi_OH_VOC * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0)
            + chi_OH_O3 * np.log(1 + D_O3s / O3s_0)
            + chi_OH_NOX * ENOX
            + chi_OH_CO * ECO
            + chi_OH_VOC * EVOC
            + chi_OH_Tatm * np.log(1 + k_Tatm * D_gst / Tatm_0)
            + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))
        )
        return np.array(D_kOH, dtype=dty)


    def df_kOH(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        df_1 = df_kOH_dCH4(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * D_CH4
        df_2 = df_kOH_dO3s(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * D_O3s
        df_3 = df_kOH_dgst(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * D_gst
        df_4 = df_kOH_dENOX(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * ENOX
        df_5 = df_kOH_dECO(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * ECO
        df_6 = df_kOH_dEVOC(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * EVOC
        df_tot = df_1 + df_2 + df_3 + df_4 + df_5 + df_6
        return [
            np.nan_to_num(df_1 / df_tot),
            np.nan_to_num(df_2 / df_tot),
            np.nan_to_num(df_3 / df_tot),
            np.nan_to_num(df_4 / df_tot),
            np.nan_to_num(df_5 / df_tot),
            np.nan_to_num(df_6 / df_tot), ]

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
mod_LSNKcover_save = mod_LSNKcover
cover_options = [
    "ESA-CCI",
    "MODIS",
    "Ramankutty1999",
    "Levavasseur2012",
    "mean-TRENDYv2",
    "CLM-45",
    "JSBACH",
    "JULES",
    "LPJ",
    "LPJ-GUESS",
    "LPX-Bern",
    "OCN",
    "ORCHIDEE",
    "VISIT",
]
AREA_wet = np.zeros([nb_regionI, nb_biome], dtype=dty)
if mod_EWETpreind != "":
    for mod_LSNKcover in cover_options:
        path = f"data/Wetlands_WETCHIMP/#DATA.Wetlands_{mod_EWETpreind}_{mod_LSNKcover}.1910s_114reg1_4bio.AREA.csv"
        TMP = load_data(path)
        for i in range(1, 114 + 1):
            for b in range(len(bio) - 2):
                AREA_wet[regionI_index[i], biome_index[bio[b]]] += TMP[i - 1, b]

mod_LSNKcover = mod_LSNKcover_save
del mod_LSNKcover_save

# preindustrial area and emissions
AREA_wet0 = np.zeros([nb_regionI], dtype=dty)
ECH4_wet0 = np.zeros([nb_regionI], dtype=dty)
for VAR, arr in [("AREA", AREA_wet0), ("ECH4", ECH4_wet0)]:
    if mod_EWETpreind != "":
        TMP = load_data(f"data/Wetlands_WETCHIMP/#DATA.Wetlands_{mod_EWETpreind}.1910s_114reg1_(1exp).{VAR}.csv",
                        slice=1)
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
    arr = np.zeros([nb_regionI], dtype=dty)
    if mod_AWETtrans != "":
        TMP = load_data(f"data/Wetlands_WETCHIMP/#DATA.Wetlands_{mod_AWETtrans}.1910s_114reg1_(4exp).AREA.csv", slice=1)
        path = f"data/Wetlands_WETCHIMP/#DATA.Wetlands_{mod_AWETtrans}.1910s_114reg1_(4exp).AREA.csv"
        lgd = [line for line in csv.reader(open(path, "r"))][0]
        for i in range(1, 114 + 1):
            arr[regionI_index[i]] += TMP[i - 1, lgd.index(sim)]
    return arr


AREA_exp0 = wetchimp("exp0")
AREA_expC = wetchimp("expC")
AREA_expT = wetchimp("expT")
AREA_expP = wetchimp("expP")

# value of perturbations
# see [Melton et al., 2013]
CO2_expC = 857.0 - 303.0
T_expT = 3.4
P_expP = np.zeros([nb_regionI], dtype=dty)
S_tmp = np.zeros([nb_regionI], dtype=dty)
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
if mod_EPFmethane == "zero":
    p_PF_CH4 = np.array([0.000], dtype=dty)
elif mod_EPFmethane == "best":
    p_PF_CH4 = np.array([0.023], dtype=dty)
elif mod_EPFmethane == "twice":
    p_PF_CH4 = np.array([0.046], dtype=dty)

# fraction of instantaneous PF emission {.}
p_PF_inst = np.array([0.0], dtype=dty)
