import csv
import os

import numpy as np
from scipy.optimize import fmin

from .a1_carbon import CO2_0, load_data
from .a7_radiative_forces import f_RF_CO2
from ..oscar_data import nb_regionI, regionI_index
from ...config import dty, mod_TEMPresp, mod_PRECradfact, mod_TEMPpattern, mod_ACIDsurf, \
    mod_PRECpattern, mod_PRECresp


##################################################
#   8. CLIMATE
##################################################

# ==========
# 8.1. GLOBE
# ==========

# ----------------------
# 8.1.A. Reconstructions
# ----------------------


# historical global temperatures from NOAA/NCDC {degC}
# from [Smith et al., 2008] updated from website
def historical_temperatures(model):
    arr = np.ones([314 + 1], dtype=dty) * np.nan
    path = f"data/HistClim_NOAA-NCDC/#DATA.HistClim_NOAA-NCDC.1880-2014.{model}.csv"
    arr[180:] = load_data(path)[:, 0]
    return arr


gst_ncdc = historical_temperatures("gst")
lst_ncdc = historical_temperatures("lst")
sst_ncdc = historical_temperatures("sst")

# historical global temperature from GISTEMP {degC}
# from [Hansen et al., 2010] updated from website
gst_giss = np.ones([314 + 1], dtype=dty) * np.nan
path = f"data/HistClim_GISTEMP/#DATA.HistClim_GISTEMP.1880-2014.gst.csv"
TMP = load_data(path)
gst_giss[180:] = TMP[:, 0]

# historical global temperature from HadCRUT4 {degC}
# from [Morice et al., 2012] updated from website
gst_had = np.ones([314 + 1], dtype=dty) * np.nan
path = "data/HistClim_HadCRUT4/#DATA.HistClim_HadCRUT4.1850-2014.gst.csv"
TMP = load_data(path)
gst_had[150:] = TMP[:, 0]

# ------------
# 8.1.1. CMIP5
# ------------

# length of CMIP5 simulations per model
lng = {
    "mean-CMIP5": 140,
    "ACCESS-10": 150,
    "ACCESS-13": 151,
    "BCC-CSM-11": 150,
    "BCC-CSM-11m": 150,
    "CanESM2": 150,
    "CCSM4": 151,
    "CNRM-CM5": 150,
    "CNRM-CM5-2": 140,
    "CSIRO-Mk360": 150,
    "GFDL-CM3": 150,
    "GFDL-ESM2G": 300,
    "GFDL-ESM2M": 300,
    "GISS-E2-H": 151,
    "GISS-E2-R": 151,
    "HadGEM2-ES": 151,
    "IPSL-CM5A-LR": 260,
    "IPSL-CM5A-MR": 140,
    "IPSL-CM5B-LR": 160,
    "MIROC5": 151,
    "MIROC-ESM": 150,
    "MPI-ESM-LR": 150,
    "MPI-ESM-MR": 150,
    "MPI-ESM-P": 150,
    "MRI-CGCM3": 150,
    "NorESM1-M": 150,
}


# load pre-processed CMIP5 results for specified model
# data related to temperature change
def load_cmip5_temperature_change(sim):
    path = f"data/Climate_CMIP5/#DATA.Climate_{mod_TEMPresp}.{lng[mod_TEMPresp]}yr_(7var).{sim}_global.csv"
    TMP = load_data(path, slice=1)

    path = f"data/Climate_CMIP5/#DATA.Climate_{mod_TEMPresp}.{lng[mod_TEMPresp]}yr_(7var).{sim}_global.csv"
    lgd = [line for line in csv.reader(open(path, "r"))][0]

    gst = TMP[:, lgd.index("tas")]
    erb = TMP[:, lgd.index("rsdt")] - TMP[:, lgd.index("rsut")] - TMP[:, lgd.index("rlut")]
    return gst, erb


gst_ctrlT, erb_ctrlT = load_cmip5_temperature_change('ctrl')
gst_quadT, erb_quadT = load_cmip5_temperature_change('quad')


# data related to precipitations change
def load_precipitation_change(sim):
    path = f"data/Climate_CMIP5/#DATA.Climate_{mod_PRECresp}.{lng[mod_PRECresp]}yr_(7var).{sim}_global.csv"
    TMP = load_data(path, slice=1)
    lgd = [line for line in csv.reader(open(path, 'r'))][0]
    gst = TMP[:, lgd.index("tas")]
    gyp = TMP[:, lgd.index("pr")]
    return gst, gyp


gst_ctrlP, gyp_ctrlP = load_precipitation_change("ctrl")
gst_quadP, gyp_quadP = load_precipitation_change("quad")

# definition of parameters
# equilibrium climate sensitivity {K}&{K/{W/m2}}
lambda_0 = np.array([0], dtype=dty)
# dynamical parameters for temperature change {.}&{yr}&{yr}
theta_0 = np.array([0], dtype=dty)
tau_gst = np.array([0], dtype=dty)
tau_gst0 = np.array([0], dtype=dty)
# precipitation sensitivities to temperature and radiative forcing {mm/K}&{mm/{W/m2}}
alpha_gyp = np.array([0], dtype=dty)
beta_gyp = np.array([0], dtype=dty)

# calculation of parameters
# equilibrium temperature change
# following the approach by [Gregory et al., 2004]
diff = gst_quadT - np.mean(gst_ctrlT, 0)


def err(var):
    clim = var[0] * (erb_quadT - np.mean(erb_ctrlT, 0)) + var[1]
    return np.sum((diff - clim) ** 2)


[TMP, T4x] = fmin(err, [-1.0, 5.0], disp=False)
lambda_0[0] = T4x / f_RF_CO2(3 * CO2_0)

# impulse response function for temperature change
diff = gst_quadT - np.mean(gst_ctrlT, 0)


def err(var):
    clim = T4x * (
            1
            - var[2] * np.exp(-np.arange(0.5, lng[mod_TEMPresp] + 0.5, 1) / var[0])
            - (1 - var[2]) * np.exp(-np.arange(0.5, lng[mod_TEMPresp] + 0.5, 1) / var[1])
    )
    return np.sum((diff - clim) ** 2)


[tau_Tfast, tau_Tslow, p_Tfast] = fmin(err, [4.0, 400.0, 0.5], disp=False)
p_Tslow = 1 - p_Tfast
if tau_Tfast > tau_Tslow:
    TMP = tau_Tfast
    tau_Tfast = tau_Tslow
    tau_Tslow = TMP
    TMP = p_Tfast
    p_Tfast = p_Tslow
    p_Tslow = TMP

# parameters of corresponding two-box model
# formulation based on [Geoffroy et al., 2013]
tau_gst[0] = tau_Tfast * tau_Tslow / (tau_Tfast * p_Tslow + tau_Tslow * p_Tfast)
tau_gst0[0] = tau_Tfast * p_Tfast + tau_Tslow * p_Tslow - tau_gst[0]
theta_0[0] = tau_gst0[0] / (tau_Tfast * p_Tslow + tau_Tslow * p_Tfast)

# sensitivities of precipitation change
diff = gyp_quadP - np.mean(gyp_ctrlP, 0)


def err(var):
    clim = np.abs(var[0]) * (gst_quadP - np.mean(gst_ctrlP, 0)) - np.abs(var[1])
    return np.sum((diff - clim) ** 2)


[alpha_gyp[0], beta_gyp[0]] = fmin(err, [10.0, -50.0], disp=False)
alpha_gyp[0] = np.abs(alpha_gyp[0])
beta_gyp[0] = -np.abs(beta_gyp[0]) / f_RF_CO2(3 * CO2_0)

# ---------------------
# 8.1.2. Precipitations
# ---------------------

# global precipitation response differentiated by forcing agent
# formulation based on [Allan et al., 2014]

# factors from [Andrews et al., 2010] (tab. 3)
# assumptions following [Allan et al., 2014]
if mod_PRECradfact == "Andrews2010":
    p_atm_CO2 = np.array([0.8], dtype=dty)
    p_atm_noCO2 = np.array([0.5], dtype=dty)
    p_atm_O3t = np.array([-0.3], dtype=dty)
    p_atm_strat = np.array([0.0], dtype=dty)  # assumed
    p_atm_scatter = np.array([0.0], dtype=dty)
    p_atm_absorb = np.array([2.5], dtype=dty)
    p_atm_cloud = np.array([0.0], dtype=dty)  # assumed
    p_atm_alb = np.array([0.0], dtype=dty)
    p_atm_solar = np.array([0.2], dtype=dty)

# factors from [Kvalevag et al., 2013] (tab. 2; highest perturbation)
# assumptions following [Allan et al., 2014] and nil if not given
elif mod_PRECradfact == "Kvalevag2013":
    p_atm_CO2 = np.array([2.0 / 3.6], dtype=dty)
    p_atm_noCO2 = np.array([0.5 / 1.8], dtype=dty)
    p_atm_O3t = np.array([0.0], dtype=dty)  # not given
    p_atm_strat = np.array([0.0], dtype=dty)  # assumed
    p_atm_scatter = np.array([0.6 / -1.4], dtype=dty)
    p_atm_absorb = np.array([3.1 / 0.5], dtype=dty)
    p_atm_cloud = np.array([0.0], dtype=dty)  # assumed
    p_atm_alb = np.array([0.0], dtype=dty)  # not given
    p_atm_solar = np.array([0.4 / 3.5], dtype=dty)

# normalization of radiative factor of precipitations
beta_gyp /= p_atm_CO2


# ==========
# 8.2. OCEAN
# ==========

# ----------------------
# 8.2.A. Reconstructions
# ----------------------

# historical sea climate from HadISST1 {degC}&{Mha}
# from [Rayner et al., 2003] updated from website
def load_sea_climate(var):
    arr = np.ones([314 + 1], dtype=dty) * np.nan
    path = f"data/HistOcean_HadISST1/#DATA.HistOcean_HadISST1.1870-2014_18x10lat.{var}.csv"
    path2 = "data/HistOcean_HadISST1/#DATA.HistOcean_HadISST1.1870-2014_18x10lat.SURF.csv"
    TMP = load_data(path)
    TMP2 = load_data(path2)
    arr[170:] = np.sum(TMP[:, :] * TMP2[:, :], 1)
    arr[170:] /= np.sum(TMP2[:, :], 1)
    return arr


sst_had = load_sea_climate("sst")
sic_had = load_sea_climate("sic")

# historical sea climate from ERSST4 {degC}
# from [Huang et al., 2015] updated from website
sst_ncdc2 = np.ones([314 + 1], dtype=dty) * np.nan
path = "data/HistOcean_ERSST4/#DATA.HistOcean_ERSST4.1854-2014_18x10lat.sst.csv"
path2 = "data/HistOcean_ERSST4/#DATA.HistOcean_ERSST4.1854-2014_18x10lat.SURF.csv"
TMP = load_data(path)
TMP2 = load_data(path2)
sst_ncdc2[154:] = np.sum(TMP[:, :] * TMP2[:, :], 1)
sst_ncdc2[154:] /= np.sum(TMP2[:, :], 1)

# historical ocean heat content
# from [Levitus et al., 2012] updated from website
# reference period is whole period
D_OHC_nodc = np.ones([311 + 1], dtype=dty) * np.nan
path = "data/HistOcean_NOAA-NODC/#DATA.HistOcean_NOAA-NODC.1955-2011.D_OHC.csv"
TMP = load_data(path)
D_OHC_nodc[255:] = TMP[:, 0]

# ------------
# 8.2.1. CMIP5
# ------------

# length of CMIP5 simulations per model or simulation
lng = {
    "mean-CMIP5": 140,
    "ACCESS-10": 150,
    "ACCESS-13": 151,
    "BCC-CSM-11": 150,
    "BCC-CSM-11m": 150,
    "CanESM2": 150,
    "CCSM4": 151,
    "CNRM-CM5": 150,
    "CNRM-CM5-2": 140,
    "CSIRO-Mk360": 150,
    "GFDL-CM3": 150,
    "GFDL-ESM2G": 300,
    "GFDL-ESM2M": 300,
    "GISS-E2-H": 151,
    "GISS-E2-R": 151,
    "HadGEM2-ES": 151,
    "IPSL-CM5A-LR": 260,
    "IPSL-CM5A-MR": 140,
    "IPSL-CM5B-LR": 160,
    "MIROC5": 151,
    "MIROC-ESM": 150,
    "MPI-ESM-LR": 150,
    "MPI-ESM-MR": 150,
    "MPI-ESM-P": 150,
    "MRI-CGCM3": 150,
    "NorESM1-M": 150,
}
prd = {
    "ctrl": "251yr",
    "hist": "1850-2005",
    "rcp26": "2006-2100",
    "rcp45": "2006-2100",
    "rcp60": "2006-2100",
    "rcp85": "2006-2100",
}

# load pre-processed CMIP5 results for specified model
# temperature patterns based on 'abrupt4xCO2'
if mod_TEMPpattern == "4xCO2":
    def load_abrupt4xCO2(sim):
        # global
        path = f"data/Climate_CMIP5/#DATA.Climate_{mod_TEMPresp}.{lng[mod_TEMPresp]}yr_(7var).{sim}_global.csv"
        TMP = load_data(path, slice=1)
        lgd = [line for line in csv.reader(open(path, "r"))][0]
        gst = TMP[:, lgd.index("tas")]

        # local
        path = f"data/Climate_CMIP5/#DATA.Climate_{mod_TEMPresp}.{lng[mod_TEMPresp]}yr_18x10lat.{sim}_tos.csv"
        path2 = f"data/Climate_CMIP5/#DATA.Climate_{mod_TEMPresp}.{lng[mod_TEMPresp]}yr_18x10lat.SURF.csv"
        TMP = load_data(path)
        TMP2 = load_data(path2)
        tos = np.sum(TMP * TMP2, 1) / np.sum(TMP2, 1)
        return gst, tos


    gst_ctrlTR, tos_ctrlTR = load_abrupt4xCO2("ctrl")
    gst_quadTR, tos_quadTR = load_abrupt4xCO2("quad")

# temperature patterns based on 'historical' and 'rcp'
elif mod_TEMPpattern == "hist&RCPs":
    def load_temp_historical(sim):
        # global
        path = f"data/ClimReg_CMIP5/#DATA.ClimReg_{mod_TEMPresp}.{prd[sim]}_(3var).{sim}_global.csv"
        if os.path.isfile(path):
            path = f"data/ClimReg_CMIP5/#DATA.ClimReg_{mod_TEMPresp}.{prd[sim]}_(3var).{sim}_global.csv"
            TMP = load_data(path, slice=1)
            lgd = [line for line in csv.reader(open(path, "r"))][0]
            gst = TMP[:, lgd.index("tas")]

        # local
        path = f"data/ClimReg_CMIP5/#DATA.ClimReg_{mod_TEMPresp}.{prd[sim]}_18x10lat.{sim}_tos.csv"
        if os.path.isfile(path):
            TMP = load_data(path)
            TMP2 = load_data(f"data/ClimReg_CMIP5/#DATA.ClimReg_{mod_TEMPresp}.251yr_18x10lat.SURF.csv")
            tos = np.sum(TMP * TMP2[:len(TMP)], 1) / np.sum(TMP2[:len(TMP)], 1)
        return gst, tos


    gst_ctrlTR, tos_ctrlTR = load_temp_historical("ctrl")
    gst_histTR, tos_histTR = load_temp_historical("hist")
    gst_rcp26TR, tos_rcp26TR = load_temp_historical("rcp26")
    gst_rcp45TR, tos_rcp45TR = load_temp_historical("rcp45")
    gst_rcp60TR, tos_rcp60TR = load_temp_historical("rcp60")
    gst_rcp85TR, tos_rcp85TR = load_temp_historical("rcp85")


def histTR(var):
    mapping = {"tos": tos_histTR, "gst": gst_histTR}
    try:
        return mapping[var]
    except KeyError:
        mapping = {"tas": tas_histTR}
        return mapping[var]


def quadTR(var):
    mapping = {"tos": tos_quadTR, "gst": gst_quadTR}
    try:
        return mapping[var]
    except KeyError:
        mapping = {"tas": tas_quadTR}
        return mapping[var]


def get_TR(var, sim):
    mapping = {
        ("tos", "rcp26"): tos_rcp26TR,
        ("gst", "rcp26"): gst_rcp26TR,
        ("tos", "rcp45"): tos_rcp45TR,
        ("gst", "rcp45"): gst_rcp45TR,
        ("tos", "rcp60"): tos_rcp60TR,
        ("gst", "rcp60"): gst_rcp60TR,
        ("tos", "rcp85"): tos_rcp85TR,
        ("gst", "rcp85"): gst_rcp85TR,
    }
    try:
        return mapping[var, sim]
    except KeyError:
        mapping = {
            ("tas", "rcp26"): tas_rcp26TR,
            ("tas", "rcp45"): tas_rcp45TR,
            ("tas", "rcp60"): tas_rcp60TR,
            ("tas", "rcp85"): tas_rcp85TR,
        }
        return mapping[var, sim]


# aggregate all experiments
def aggregate(VAR):
    if mod_TEMPpattern == "4xCO2":
        return np.array(list(quadTR(VAR)), dtype=dty)
    elif mod_TEMPpattern == "hist&RCPs":
        res = []
        for sim in ["rcp26", "rcp45", "rcp60", "rcp85"]:
            path = f"data/ClimReg_CMIP5/#DATA.ClimReg_{mod_TEMPresp}.{prd[sim]}_(3var).{sim}_global.csv"
            if os.path.isfile(path):
                res += list(histTR(VAR)) + list(get_TR(VAR, sim))
        if res == []:
            res.extend(list(histTR(VAR)))
        return np.array(res, dtype=dty)
    else:
        raise ValueError


if mod_TEMPpattern in ("4xCO2", "hist&RCPs"):
    tos_allTR = aggregate("tos")
    gst_allTR = aggregate("gst")

# definition of parameter
# scaling of sea surface temperature over gst {.}
w_reg_sst = np.array([0], dtype=dty)

# fit of parameter
diff = tos_allTR - np.mean(tos_ctrlTR, 0)


def err(var):
    clim = var[0] * (gst_allTR - np.mean(gst_ctrlTR, 0))
    return np.sum((diff - clim) ** 2)


[w_reg_sst[0]] = fmin(err, [1], disp=False)

# -------------------
# 8.2.2. Heat content
# -------------------

# conversion factor from radiative forcing {W/m2} to heat {ZJ}
alpha_OHC = np.array([510_072e9 * 3600 * 24 * 365.25 / 1e21], dtype=dty)

# fraction of energy used to heat the land, the atmosphere, and to melt ice {.}
# from [Otto et al., 2013] (SI)
p_OHC = np.array([1 - 0.03 - 0.01 - 0.02], dtype=dty)

# --------------------
# 8.2.3. Acidification
# --------------------

# fonction relating surface pH and atmospheric CO2
# from [Tans, 2009] based on [Tans, 2008]
if mod_ACIDsurf == "Tans2009":

    def f_pH(D_CO2):
        D_pH = -0.85 * np.log(1 + D_CO2 / CO2_0)
        return np.array(D_pH, dtype=dty)


# from [Bernie et al., 2010]
elif mod_ACIDsurf == "Bernie2010":

    def f_pH(D_CO2):
        D_pH = (
                -0.00173 * D_CO2
                + 1.3264e-6 * (2 * CO2_0 * D_CO2 + D_CO2 ** 2)
                - 4.4943e-10 * (
                        3 * D_CO2 * CO2_0 ** 2 + 3 * CO2_0 * D_CO2 ** 2 + D_CO2 ** 3)
        )
        return np.array(D_pH, dtype=dty)


# ----------------
# 8.2.4. Sea-level
# ----------------

# TODO

# =========
# 8.3. LAND
# =========

# ----------------------
# 8.3.A. Reconstructions
# ----------------------


# historical local climate from CRU {degC}&{mm}
# from [Harris et al., 2014]
def load_historical_climate(var):
    arr = np.zeros([314 + 1, nb_regionI], dtype=dty)
    path = f"data/HistLand_CRU/#DATA.HistLand_CRU.1901-2014_114reg1.{var}.csv"
    TMP = load_data(path)

    path = "data/HistLand_CRU/#DATA.HistLand_CRU.1901-2014_114reg1.AREA.csv"
    TMP2 = load_data(path)

    for i in range(1, 114 + 1):
        arr[201:, regionI_index[i]] += TMP[:, i - 1] * TMP2[:, i - 1]
    TMP = np.zeros([314 + 1, nb_regionI], dtype=dty)
    for i in range(1, 114 + 1):
        TMP[201:, regionI_index[i]] += TMP2[:, i - 1]
    arr /= TMP[:, :]
    return arr


lst_cru = load_historical_climate("lst")
lyp_cru = load_historical_climate("lyp")

# historical local climate from GHCN+CAMS {degC}
# from [Fan and van den Dool, 2008] updated from website
lst_ghcn = np.zeros([314 + 1, nb_regionI], dtype=dty)
path = f"data/HistLand_GHCN-CAMS/#DATA.HistLand_GHCN-CAMS.1948-2014_114reg1.lst.csv"
path2 = "data/HistLand_GHCN-CAMS/#DATA.HistLand_GHCN-CAMS.1948-2014_114reg1.AREA.csv"
TMP = load_data(path)
TMP2 = load_data(path2)
for i in range(1, 114 + 1):
    lst_ghcn[248:, regionI_index[i]] += TMP[:, i - 1] * TMP2[:, i - 1]
TMP = np.zeros([314 + 1, nb_regionI], dtype=dty)
for i in range(1, 114 + 1):
    TMP[248:, regionI_index[i]] += TMP2[:, i - 1]
lst_ghcn /= TMP[:, :]

# ------------
# 8.3.1. CMIP5
# ------------

# length of CMIP5 simulations per model or simulation
lng = {
    "mean-CMIP5": 140,
    "ACCESS-10": 150,
    "ACCESS-13": 151,
    "BCC-CSM-11": 150,
    "BCC-CSM-11m": 150,
    "CanESM2": 150,
    "CCSM4": 151,
    "CNRM-CM5": 150,
    "CNRM-CM5-2": 140,
    "CSIRO-Mk360": 150,
    "GFDL-CM3": 150,
    "GFDL-ESM2G": 300,
    "GFDL-ESM2M": 300,
    "GISS-E2-H": 151,
    "GISS-E2-R": 151,
    "HadGEM2-ES": 151,
    "IPSL-CM5A-LR": 260,
    "IPSL-CM5A-MR": 140,
    "IPSL-CM5B-LR": 160,
    "MIROC5": 151,
    "MIROC-ESM": 150,
    "MPI-ESM-LR": 150,
    "MPI-ESM-MR": 150,
    "MPI-ESM-P": 150,
    "MRI-CGCM3": 150,
    "NorESM1-M": 150,
}
prd = {
    "ctrl": "251yr",
    "hist": "1850-2005",
    "rcp26": "2006-2100",
    "rcp45": "2006-2100",
    "rcp60": "2006-2100",
    "rcp85": "2006-2100",
}

# load pre-processed CMIP5 results for specified model
# temperature patterns based on 'abrupt4xCO2'
if mod_TEMPpattern == "4xCO2":
    def load_temp_patterns(sim):
        # global
        path = f"data/Climate_CMIP5/#DATA.Climate_{mod_TEMPresp}.{lng[mod_TEMPresp]}yr_(7var).{sim}global.csv"
        TMP = load_data(path, slice=1)
        lgd = [line for line in csv.reader(open(path, "r"))][0]
        gst = TMP[:, lgd.index("tas")]

        # local
        # tas = np.zeros([lng[mod_TEMPresp],nb_regionI], dtype=dty)
        path = f"data/Climate_CMIP5/#DATA.Climate_{mod_TEMPresp}.{lng[mod_TEMPresp]}yr_114reg1.{sim}_tas.csv"
        path2 = f"data/Climate_CMIP5/#DATA.Climate_{mod_TEMPresp}.{lng[mod_TEMPresp]}yr_114reg1.AREA.csv"
        TMP = load_data(path)
        TMP2 = load_data(path2)
        tas = np.zeros([len(TMP), nb_regionI], dtype=dty)
        area = np.zeros([len(TMP), nb_regionI], dtype=dty)
        for i in range(1, 114 + 1):
            tas[:, regionI_index[i]] += TMP[:, i - 1] * TMP2[:len(TMP), i - 1]
            area[:, regionI_index[i]] += TMP2[:len(TMP), i - 1]
        tas /= area
        tas[np.isnan(tas) | np.isinf(tas)] = 0

        return gst, tas, area


    gst_ctrlTR, tas_ctrlTR, AREA_ctrlTR = load_temp_patterns("ctrl")
    gst_quadTR, tas_quadTR, AREA_quadTR = load_temp_patterns("quad")


# temperature patterns based on 'historical' and 'rcp'
elif mod_TEMPpattern == "hist&RCPs":
    def load_temperature_patterns(sim):
        # global
        path = f"data/ClimReg_CMIP5/#DATA.ClimReg_{mod_TEMPresp}.{prd[sim]}_(3var).{sim}_global.csv"
        if os.path.isfile(path):
            TMP = load_data(path, slice=1)
            lgd = [line for line in csv.reader(open(path, "r"))][0]
            gst = TMP[:, lgd.index("tas")]

        # local
        path = f"data/ClimReg_CMIP5/#DATA.ClimReg_{mod_TEMPresp}.{prd[sim]}_114reg1.{sim}_tas.csv"
        if os.path.isfile(path):
            path2 = f"data/ClimReg_CMIP5/#DATA.ClimReg_{mod_TEMPresp}.251yr_114reg1.AREA.csv"
            TMP = load_data(path)
            TMP2 = load_data(path2)
            tas = np.zeros([len(TMP), nb_regionI], dtype=dty)
            area = np.zeros([len(TMP), nb_regionI], dtype=dty)
            for i in range(1, 114 + 1):
                tas[:, regionI_index[i]] += TMP[:, i - 1] * TMP2[:len(TMP), i - 1]
                area[:, regionI_index[i]] += TMP2[:len(TMP), i - 1]
            tas /= area
            tas[np.isnan(tas) | np.isinf(tas)] = 0

            return gst, tas, area


    gst_ctrlTR, tas_ctrlTR, AREA_ctrlTR = load_temperature_patterns("ctrl")
    gst_histTR, tas_histTR, AREA_histTR = load_temperature_patterns("hist")
    gst_rcp26TR, tas_rcp26TR, AREA_rcp26TR = load_temperature_patterns("rcp26")
    gst_rcp45TR, tas_rcp45TR, AREA_rcp45TR = load_temperature_patterns("rcp45")
    gst_rcp60TR, tas_rcp60TR, AREA_rcp60TR = load_temperature_patterns("rcp60")
    gst_rcp85TR, tas_rcp85TR, AREA_rcp85TR = load_temperature_patterns("rcp85")

# aggregate all experiments
# FIXME: this seems to be repeated. Does anything change from first definition
# of <VAR>_allTR to now?
nrcp = 0


def aggregate(VAR):
    global nrcp
    if mod_TEMPpattern == "4xCO2":
        return np.array(list(quadTR(VAR)), dtype=dty)
    elif mod_TEMPpattern == "hist&RCPs":
        res = []
        nrcp = 0
        for sim in ["rcp26", "rcp45", "rcp60", "rcp85"]:
            path = f"data/ClimReg_CMIP5/#DATA.ClimReg_{mod_TEMPresp}.{prd[sim]}_(3var).{sim}_global.csv"
            if os.path.isfile(path):
                nrcp += 1
                res += list(histTR(VAR)) + list(get_TR(VAR, sim))
        if res == []:
            res += list(histTR(VAR))
        return np.array(res, dtype=dty)


gst_allTR = aggregate("gst")
tas_allTR = aggregate("tas")


# decadal means
def decadal_means(base):
    if mod_TEMPpattern == "4xCO2":
        arr = np.zeros([lng[mod_TEMPresp] - 10] + list(np.shape(base)[1:]), dtype=dty)
        for t in range(lng[mod_TEMPresp] - 10):
            arr[t, ...] = np.mean(base[t:t + 10, ...], 0)
    elif mod_TEMPpattern == "hist&RCPs":
        arr = np.zeros([(156 - 10) + nrcp * (95 - 10)] + list(np.shape(base)[1:]), dtype=dty)
        for t in range(156 - 10):
            arr[t, ...] = np.mean(base[t:t + 10, ...], 0)
        for n in range(nrcp):
            for t in range(95 - 10):
                arr[(156 - 10) + n * (95 - 10) + t, ...] = np.mean(
                    base[(156 - 10) + n * (95 - 10) + t:(156 - 10) + n * (95 - 10) + t + 10, ...], 0)
    else:
        raise ValueError('invalid module!')
    return arr


gst_decTR = decadal_means(gst_allTR)
tas_decTR = decadal_means(tas_allTR)

# definition of parameter
# scaling of local surface temperature over gst {.}
w_reg_lst = np.zeros([nb_regionI], dtype=dty)

# fit of parameter
for i in range(nb_regionI):
    diff = tas_decTR[:, i] - np.mean(tas_ctrlTR[:, i], 0)


    def err(var):
        clim = var[0] * (gst_decTR - np.mean(gst_ctrlTR, 0))
        return np.sum((diff - clim) ** 2)


    [w_reg_lst[i]] = fmin(err, [1], disp=False)

# load pre-processed CMIP5 results for specified model
# precipitation patterns based on 'abrupt4xCO2'
if mod_PRECpattern == "4xCO2":
    def load_mod_prec(sim):
        # global
        path = f"data/Climate_CMIP5/#DATA.Climate_{mod_PRECresp}.{lng[mod_PRECresp]}yr_(7var).{sim}_global.csv"
        TMP = load_data(path, slice=1)
        lgd = [line for line in csv.reader(open(path, "r"))][0]
        gyp = TMP[:, lgd.index("pr")]

        # local
        path = f"data/Climate_CMIP5/#DATA.Climate_{mod_PRECresp}.{lng[mod_PRECresp]}yr_114reg1.{sim}_pr.csv"
        path2 = f"data/Climate_CMIP5/#DATA.Climate_{mod_PRECresp}.{lng[mod_PRECresp]}yr_114reg1.AREA.csv"
        TMP = load_data(path)
        TMP2 = load_data(path2)
        pr = p.zeros([len(TMP), nb_regionI], dtype=dty)
        area = np.zeros([len(TMP), nb_regionI], dtype=dty)
        for i in range(1, 114 + 1):
            pr[:, regionI_index[i]] += TMP[:, i - 1] * TMP2[:len(TMP), i - 1]
            area[:, regionI_index[i]] += TMP2[:len(TMP), i - 1]
        pr /= area
        pr[np.isnan(pr) | np.isinf(pr)] = 0
        return gyp, pr, area


    gyp_ctrlPR, pr_ctrlPR, AREA_ctrlPR = load_mod_prec("ctrl")
    gyp_quadPR, pr_quadPR, AREA_quadPR = load_mod_prec("quad")


# precipitation patterns based on 'historical' and 'rcp'
elif mod_PRECpattern == "hist&RCPs":
    def load_precipitation(sim):
        # global
        path = f"data/ClimReg_CMIP5/#DATA.ClimReg_{mod_PRECresp}.{prd[sim]}_(3var).{sim}_global.csv"
        if os.path.isfile(path):
            TMP = load_data(path, slice=1)
            lgd = [line for line in csv.reader(open(path, "r"))][0]
            gyp = TMP[:, lgd.index("pr")]

        # local
        path = f"data/ClimReg_CMIP5/#DATA.ClimReg_{mod_PRECresp}.{prd[sim]}_114reg1.{sim}_pr.csv"
        path2 = f"data/ClimReg_CMIP5/#DATA.ClimReg_{mod_PRECresp}.251yr_114reg1.AREA.csv"
        if os.path.isfile(path):
            TMP = load_data(path)
            TMP2 = load_data(path2)
            pr = np.zeros([len(TMP), nb_regionI], dtype=dty)
            area = np.zeros([len(TMP), nb_regionI], dtype=dty)
            for i in range(1, 114 + 1):
                pr[:, regionI_index[i]] += TMP[:, i - 1] * TMP2[:len(TMP), i - 1]
                area[:, regionI_index[i]] += TMP2[:len(TMP), i - 1]
            pr /= area
            pr[np.isnan(pr) | np.isinf(pr)] = 0

        return gyp, pr, area


    gyp_ctrlPR, pr_ctrlPR, AREA_ctrlPR = load_precipitation("ctrl")
    gyp_histPR, pr_histPR, AREA_histPR = load_precipitation("hist")
    gyp_rcp26PR, pr_rcp26PR, AREA_rcp26PR = load_precipitation("rcp26")
    gyp_rcp45PR, pr_rcp45PR, AREA_rcp45PR = load_precipitation("rcp45")
    gyp_rcp60PR, pr_rcp60PR, AREA_rcp60PR = load_precipitation("rcp60")
    gyp_rcp85PR, pr_rcp85PR, AREA_rcp85PR = load_precipitation("rcp85")


# aggregate all experiments
def quadPR(var):
    return {"gyp": gyp_quadPR, "pr": pr_quadPR}[var]


def histPR(var):
    return {"gyp": gyp_histPR, "pr": pr_histPR}[var]


def getPR(var, sim):
    return {
        ("gyp", "rcp26"): gyp_rcp26PR,
        ("gyp", "rcp45"): gyp_rcp45PR,
        ("gyp", "rcp60"): gyp_rcp60PR,
        ("gyp", "rcp85"): gyp_rcp85PR,
        ("pr", "rcp26"): pr_rcp26PR,
        ("pr", "rcp45"): pr_rcp45PR,
        ("pr", "rcp60"): pr_rcp60PR,
        ("pr", "rcp85"): pr_rcp85PR,
    }[var, sim]


def aggregate(VAR):
    global nrcp

    if mod_PRECpattern == "4xCO2":
        return np.array(list(quadPR(VAR)), dtype=dty)
    elif mod_PRECpattern == "hist&RCPs":
        res = []
        nrcp = 0
        for sim in ["rcp26", "rcp45", "rcp60", "rcp85"]:
            path = f"data/ClimReg_CMIP5/#DATA.ClimReg_{mod_PRECresp}.{prd[sim]}_(3var).{sim}_global.csv"
            if os.path.isfile(path):
                nrcp += 1
                res += list(histPR(VAR)) + list(getPR(VAR, sim))
        if res == []:
            res += list(histPR(VAR))
        return np.array(res, dtype=dty)


gyp_allPR = aggregate("gyp")
pr_allPR = aggregate("pr")


# decadal means
def decadal_means(base):
    if mod_PRECpattern == "4xCO2":
        dec = np.zeros([lng[mod_PRECresp] - 10] + list(np.shape(base)[1:]), dtype=dty)
        for t in range(lng[mod_PRECresp] - 10):
            dec[t, ...] = np.mean(base[t:t + 10, ...], 0)
    elif mod_PRECpattern == "hist&RCPs":
        dec = np.zeros([(156 - 10) + nrcp * (95 - 10)] + list(np.shape(base)[1:]), dtype=dty)
        for t in range(156 - 10):
            dec[t, ...] = np.mean(base[t:t + 10, ...], 0)
        for n in range(nrcp):
            for t in range(95 - 10):
                dec[(156 - 10) + n * (95 - 10) + t, ...] = np.mean(
                    base[(156 - 10) + n * (95 - 10) + t:(156 - 10) + n * (95 - 10) + t + 10, ...], 0)
    else:
        raise ValueError(mod_PRECpattern)
    return dec


gyp_decPR = decadal_means(gyp_allPR)
pr_decPR = decadal_means(pr_allPR)

# definition of parameter
# scaling of local yearly precipitations over gst {.}
w_reg_lyp = np.zeros([nb_regionI], dtype=dty)

# fit of parameter
for i in range(nb_regionI):
    diff = pr_decPR[:, i] - np.mean(pr_ctrlPR[:, i], 0)


    def err(var):
        clim = var[0] * (gyp_decPR - np.mean(gyp_ctrlPR, 0))
        return np.sum((diff - clim) ** 2)


    [w_reg_lyp[i]] = fmin(err, [1], disp=False)
