from scipy.optimize import fmin

from ..oscar_data import *
from ..oscar_data import nb_ODS, ODS
from ...config import dty, mod_HVSNKtau, mod_HVSNKtrans, mod_HVSNKcirc
from .a1_carbon import load_data

##################################################
#   3. NITROUS OXIDE
##################################################

# ===============
# 3.1. ATMOSPHERE
# ===============

# conversion of N2O from {ppb} to {TgN}
alpha_N2O = 0.1765 * np.array([28.0], dtype=dty)

# historic N2O from IPCC-AR5 {ppb}
# from [IPCC WG1, 2013] annexe 2
N2O_ipcc = np.ones([311 + 1], dtype=dty) * np.nan
TMP = np.array(
    [line for line in csv.reader(
        open("data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1750-2011.N2O.csv", "r"))],
    dtype=dty,
)
N2O_ipcc[50:] = TMP[:, 0]
N2O_0 = np.array([N2O_ipcc[50]], dtype=dty)

# historic N2O from CMIP5 {ppb}
# from [Meinshausen et al., 2011]
N2O_cmip5 = np.ones([305 + 1], dtype=dty) * np.nan
TMP = np.array(
    [line for line in
     csv.reader(open("data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005.N2O.csv", "r"))],
    dtype=dty
)
N2O_cmip5[65:] = TMP[:, 0]

# historic N2O from AGAGE {ppb}
# from [Prinn et al., 2013] updated from the website
N2O_agage = np.ones([313 + 1], dtype=dty) * np.nan
TMP = np.array(
    [line for line in csv.reader(
        open("data/HistAtmo_AGAGE/#DATA.HistAtmo_AGAGE.1979-2013.N2O_global.csv", "r"))],
    dtype=dty,
)
N2O_agage[279:] = TMP[:, 0]

# historic N2O from Law Dome ice cores {ppm}
# from [MacFarling Meure et al., 2006]
N2O_lawdome = np.array(
    [
        line
        for line in csv.reader(open(
        "data/HistAtmo_NOAA-NCDC/#DATA.HistAtmo_NOAA-NCDC.(IceCores).N2O_lawdome.csv",
        "r"))][1:],
    dtype=dty,
)

# load RCP concentrations {ppb}
# from [Meinshausen et al., 2011]
N2O_rcp = np.ones([800 + 1, 6], dtype=dty) * np.nan
n = -1
for rcp in ["rcp26", "rcp45", "rcp60", "rcp85", "rcp45to26", "rcp60to45"]:
    n += 1
    TMP = np.array(
        [line for line in csv.reader(
            open("data/Scenario_ECP/#DATA.Scenario_ECP.2000-2500." + rcp + "_N2O.csv",
                 "r"))],
        dtype=dty,
    )
    N2O_rcp[300:, n] = TMP[:, 0]

# ==============
# 3.2. CHEMISTRY
# ==============

# ----------------
# 3.2.1. Lifetimes
# ----------------

# preindustrial lifetime of N2O {yr}
# average from [Prather et al., 2012]
if mod_HVSNKtau == "Prather2015":
    tau_N2O_hv = np.array([123.0], dtype=dty)
# variations also from [Prather et al., 2015] (table 2)
elif mod_HVSNKtau == "GMI":
    tau_N2O_hv = np.array([137.4], dtype=dty) * 123.0 / 132.5
elif mod_HVSNKtau == "GEOSCCM":
    tau_N2O_hv = np.array([120.2], dtype=dty) * 123.0 / 132.5
elif mod_HVSNKtau == "G2d-M":
    tau_N2O_hv = np.array([127.0], dtype=dty) * 123.0 / 132.5
elif mod_HVSNKtau == "G2d":
    tau_N2O_hv = np.array([129.5], dtype=dty) * 123.0 / 132.5
elif mod_HVSNKtau == "Oslo-c29":
    tau_N2O_hv = np.array([126.1], dtype=dty) * 123.0 / 132.5
elif mod_HVSNKtau == "Oslo-c36":
    tau_N2O_hv = np.array([146.7], dtype=dty) * 123.0 / 132.5
elif mod_HVSNKtau == "UCI-c29":
    tau_N2O_hv = np.array([126.2], dtype=dty) * 123.0 / 132.5
elif mod_HVSNKtau == "UCI-c36":
    tau_N2O_hv = np.array([146.2], dtype=dty) * 123.0 / 132.5

# lag used for lagged concentrations {yr}
# adjusted; 3yr is typical WMO value
tau_lag = np.array([3.0], dtype=dty)

# ------------------
# 3.2.2. Prather2015
# ------------------

# estimate EESC for use to deduce strato sink sensitivity
# ODS concentration from [IPCC WG1, 2013] (annex 2) in year 2005
EESC_hv = np.zeros([nb_ODS], dtype=dty)
EESC_hv0 = np.zeros([nb_ODS], dtype=dty)
for VAR in ODS:
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
    EESC_hv[ODS.index(VAR)] = TMP[45, :]
EESC_hv0[ODS.index("CH3Br")] = 5.8
EESC_hv0[ODS.index("CH3Cl")] = 480.0

#  fractional releases from [Newman et al., 2007]
for arr in [EESC_hv, EESC_hv0]:
    arr *= np.array([0.47, 0.23, 0.29, 0.12, 0.05, 0.56, 0.67, 0.13, 0.08, 0.01, 0.62, 0.62, 0.28, 0.65, 0.60, 0.44],
                    dtype=dty)
    arr *= np.array([3, 2, 3, 2, 1, 4, 3, 1, 2, 1, 1 + 60 * 1, 0 + 60 * 2, 0 + 60 * 1, 0 + 60 * 2, 0 + 60 * 1, 1],
                    dtype=dty)

EESC_hv = np.sum(EESC_hv)
EESC_hv0 = np.sum(EESC_hv0)

# sensitivity of strato sink to N2O, ODSs and age-of-air
# from [Prather et al., 2015] (table 3)
# change in age-of-air based on [Fleming et al., 2011] (figure 12)
if mod_HVSNKtrans == "Prather2015":
    chi_hv_N2O = np.array([0.065], dtype=dty)
    chi_hv_EESC = np.array([0.04 / np.log(EESC_hv / EESC_hv0)], dtype=dty)
    chi_hv_age = np.array([0], dtype=dty)
elif mod_HVSNKtrans == "G2d":
    chi_hv_N2O = np.array([0.018 / np.log(321 / 270.0)], dtype=dty)
    chi_hv_EESC = np.array([0.048 / np.log(EESC_hv / EESC_hv0)], dtype=dty)
    chi_hv_age = np.array([0.011 / np.log(4.0 / 4.5)], dtype=dty)
elif mod_HVSNKtrans == "Oslo-c29":
    chi_hv_N2O = np.array([0.010 / np.log(321 / 270.0)], dtype=dty)
    chi_hv_EESC = np.array([0.033 / np.log(EESC_hv / EESC_hv0)], dtype=dty)
    chi_hv_age = np.array([0], dtype=dty)
elif mod_HVSNKtrans == "UCI-c29":
    chi_hv_N2O = np.array([0.012 / np.log(321 / 270.0)], dtype=dty)
    chi_hv_EESC = np.array([0.029 / np.log(EESC_hv / EESC_hv0)], dtype=dty)
    chi_hv_age = np.array([0], dtype=dty)
# but also from [Prather et al., 2012]
elif mod_HVSNKtrans == "Prather2012":
    chi_hv_N2O = np.array([0.08], dtype=dty)
    chi_hv_EESC = np.array([0], dtype=dty)
    chi_hv_age = np.array([0], dtype=dty)


# expression of hv sink function {.}
# adapted from [Prather et al., 2015]
def f_hv(D_N2O, D_EESC, D_gst):
    from .a5_ozone import EESC_0
    D_hv = (
            np.exp(
                chi_hv_N2O * np.log(1 + D_N2O / N2O_0)
                + chi_hv_EESC * np.log(1 + D_EESC / EESC_0)
                - chi_hv_age * np.log(1 + gamma_age * D_gst)
            )
            - 1
    )
    return np.array(D_hv, dtype=dty)


def df_hv_dN2O(D_N2O, D_EESC, D_gst):
    from .a5_ozone import EESC_0
    D_hv = (
            chi_hv_N2O
            / (N2O_0 + D_N2O)
            * np.exp(
        chi_hv_N2O * np.log(1 + D_N2O / N2O_0)
        + chi_hv_EESC * np.log(1 + D_EESC / EESC_0)
        - chi_hv_age * np.log(1 + gamma_age * D_gst)
    )
    )
    return np.array(D_hv, dtype=dty)


def df_hv_dEESC(D_N2O, D_EESC, D_gst):
    from .a5_ozone import EESC_0
    D_hv = (
            chi_hv_EESC
            / (EESC_0 + D_EESC)
            * np.exp(
        chi_hv_N2O * np.log(1 + D_N2O / N2O_0)
        + chi_hv_EESC * np.log(1 + D_EESC / EESC_0)
        - chi_hv_age * np.log(1 + gamma_age * D_gst)
    )
    )
    return np.array(D_hv, dtype=dty)


def df_hv_dgst(D_N2O, D_EESC, D_gst):
    from .a5_ozone import EESC_0
    D_hv = (
            chi_hv_age
            * gamma_age
            / (1 + gamma_age * D_gst)
            * np.exp(
        chi_hv_N2O * np.log(1 + D_N2O / N2O_0)
        + chi_hv_EESC * np.log(1 + D_EESC / EESC_0)
        - chi_hv_age * np.log(1 + gamma_age * D_gst)
    )
    )
    return np.array(D_hv, dtype=dty)


def df_hv(D_N2O, D_EESC, D_gst):
    df_1 = df_hv_dN2O(D_N2O, D_EESC, D_gst) * D_N2O
    df_2 = df_hv_dEESC(D_N2O, D_EESC, D_gst) * D_EESC
    df_3 = df_hv_dgst(D_N2O, D_EESC, D_gst) * D_gst
    df_tot = df_1 + df_2 + df_3
    return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot),
            np.nan_to_num(df_3 / df_tot)]


# --------------
# 3.2.3. CCMVAL2
# --------------

# load pre-processed CCMVal2 results for specified model
age_atm = np.zeros([2099 - 1961 + 1], dtype=dty)
ta_atm = np.zeros([2099 - 1961 + 1], dtype=dty)
for VAR, arr in [('age', age_atm), ('ta', ta_atm)]:
    path = f"data/Atmosphere_CCMVal2/#DATA.Atmosphere_{mod_HVSNKcirc}.1961-2099_(1lvl).{VAR}.csv"
    TMP = load_data(path, slice=1)
    arr[:] = TMP[:, 0]

# definition of parameter
# sensitivity of stratospheric lag to temperature {./K}
gamma_age = np.array([0], dtype=dty)

# fit of parameter
ratio = np.mean(age_atm[:10]) / age_atm[:]


def err(var):
    clim = 1 + var[0] * (ta_atm[:] - np.mean(ta_atm[:10]))
    return np.sum((ratio - clim) ** 2)


[gamma_age[0]] = fmin(err, [0], disp=False)
