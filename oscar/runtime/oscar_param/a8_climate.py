import csv
import os

import numpy as np
from scipy.optimize import fmin

from .a1_carbon import CO2_0
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
for var in ["gst", "lst", "sst"]:
    exec(var + "_ncdc = np.ones([314+1], dtype=dty) * np.nan")
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/HistClim_NOAA-NCDC/#DATA.HistClim_NOAA-NCDC.1880-2014." + var + ".csv",
                "r")
        )],
        dtype=dty,
    )
    exec(var + "_ncdc[180:] = TMP[:,0]")

# historical global temperature from GISTEMP {degC}
# from [Hansen et al., 2010] updated from website
for var in ["gst"]:
    exec(var + "_giss = np.ones([314+1], dtype=dty) * np.nan")
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open("data/HistClim_GISTEMP/#DATA.HistClim_GISTEMP.1880-2014." + var + ".csv",
                 "r"))],
        dtype=dty,
    )
    exec(var + "_giss[180:] = TMP[:,0]")

# historical global temperature from HadCRUT4 {degC}
# from [Morice et al., 2012] updated from website
for var in ["gst"]:
    exec(var + "_had = np.ones([314+1], dtype=dty) * np.nan")
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/HistClim_HadCRUT4/#DATA.HistClim_HadCRUT4.1850-2014." + var + ".csv",
                "r")
        )],
        dtype=dty,
    )
    exec(var + "_had[150:] = TMP[:,0]")

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
for sim in ["ctrl", "quad"]:
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/Climate_CMIP5/#DATA.Climate_"
                + mod_TEMPresp
                + "."
                + str(lng[mod_TEMPresp])
                + "yr_(7var)."
                + sim
                + "_global.csv",
                "r",
            )
        )][1:],
        dtype=dty,
    )
    lgd = [
        line
        for line in csv.reader(
            open(
                "data/Climate_CMIP5/#DATA.Climate_"
                + mod_TEMPresp
                + "."
                + str(lng[mod_TEMPresp])
                + "yr_(7var)."
                + sim
                + "_global.csv",
                "r",
            )
        )][0]
    exec("gst_" + sim + 'T = TMP[:,lgd.index("tas")]')
    exec(
        "erb_" + sim + 'T = TMP[:,lgd.index("rsdt")] - TMP[:,lgd.index("rsut")] - TMP[:,lgd.index("rlut")]')
# data related to precipitations change
for sim in ["ctrl", "quad"]:
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/Climate_CMIP5/#DATA.Climate_"
                + mod_PRECresp
                + "."
                + str(lng[mod_PRECresp])
                + "yr_(7var)."
                + sim
                + "_global.csv",
                "r",
            )
        )][1:],
        dtype=dty,
    )
    lgd = [
        line
        for line in csv.reader(
            open(
                "data/Climate_CMIP5/#DATA.Climate_"
                + mod_PRECresp
                + "."
                + str(lng[mod_PRECresp])
                + "yr_(7var)."
                + sim
                + "_global.csv",
                "r",
            )
        )][0]
    exec("gst_" + sim + 'P = TMP[:,lgd.index("tas")]')
    exec("gyp_" + sim + 'P = TMP[:,lgd.index("pr")]')

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
for var in ["sst", "sic"]:
    exec(var + "_had = np.ones([314+1], dtype=dty) * np.nan")
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/HistOcean_HadISST1/#DATA.HistOcean_HadISST1.1870-2014_18x10lat." + var + ".csv",
                "r")
        )],
        dtype=dty,
    )
    TMP2 = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/HistOcean_HadISST1/#DATA.HistOcean_HadISST1.1870-2014_18x10lat.SURF.csv",
                "r")
        )],
        dtype=dty,
    )
    exec(var + "_had[170:] = np.sum(TMP[:,:]*TMP2[:,:],1)")
    exec(var + "_had[170:] /= np.sum(TMP2[:,:],1)")

# historical sea climate from ERSST4 {degC}
# from [Huang et al., 2015] updated from website
for var in ["sst"]:
    exec(var + "_ncdc2 = np.ones([314+1], dtype=dty) * np.nan")
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/HistOcean_ERSST4/#DATA.HistOcean_ERSST4.1854-2014_18x10lat." + var + ".csv",
                "r")
        )],
        dtype=dty,
    )
    TMP2 = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/HistOcean_ERSST4/#DATA.HistOcean_ERSST4.1854-2014_18x10lat.SURF.csv",
                "r")
        )],
        dtype=dty,
    )
    exec(var + "_ncdc2[154:] = np.sum(TMP[:,:]*TMP2[:,:],1)")
    exec(var + "_ncdc2[154:] /= np.sum(TMP2[:,:],1)")

# historical ocean heat content
# from [Levitus et al., 2012] updated from website
# reference period is whole period
for var in ["D_OHC"]:
    exec(var + "_nodc = np.ones([311+1], dtype=dty) * np.nan")
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open("data/HistOcean_NOAA-NODC/#DATA.HistOcean_NOAA-NODC.1955-2011.D_OHC.csv",
                 "r"))],
        dtype=dty,
    )
    exec(var + "_nodc[255:] = TMP[:,0]")

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
    for sim in ["ctrl", "quad"]:
        # global
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/Climate_CMIP5/#DATA.Climate_"
                    + mod_TEMPresp
                    + "."
                    + str(lng[mod_TEMPresp])
                    + "yr_(7var)."
                    + sim
                    + "_global.csv",
                    "r",
                )
            )][1:],
            dtype=dty,
        )
        lgd = [
            line
            for line in csv.reader(
                open(
                    "data/Climate_CMIP5/#DATA.Climate_"
                    + mod_TEMPresp
                    + "."
                    + str(lng[mod_TEMPresp])
                    + "yr_(7var)."
                    + sim
                    + "_global.csv",
                    "r",
                )
            )][0]
        exec("gst_" + sim + 'TR = TMP[:,lgd.index("tas")]')
        # local
        for var in ["tos"]:
            TMP = np.array(
                [
                    line
                    for line in csv.reader(
                    open(
                        "data/Climate_CMIP5/#DATA.Climate_"
                        + mod_TEMPresp
                        + "."
                        + str(lng[mod_TEMPresp])
                        + "yr_18x10lat."
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
                        "data/Climate_CMIP5/#DATA.Climate_"
                        + mod_TEMPresp
                        + "."
                        + str(lng[mod_TEMPresp])
                        + "yr_18x10lat.SURF.csv",
                        "r",
                    )
                )],
                dtype=dty,
            )
            exec(var + "_" + sim + "TR = np.sum(TMP*TMP2,1)/np.sum(TMP2,1)")

# temperature patterns based on 'historical' and 'rcp'
elif mod_TEMPpattern == "hist&RCPs":
    for sim in ["ctrl", "hist", "rcp26", "rcp45", "rcp60", "rcp85"]:
        # global
        if os.path.isfile(
                "data/ClimReg_CMIP5/#DATA.ClimReg_" + mod_TEMPresp + "." + prd[
                    sim] + "_(3var)." + sim + "_global.csv"
        ):
            TMP = np.array(
                [
                    line
                    for line in csv.reader(
                    open(
                        "data/ClimReg_CMIP5/#DATA.ClimReg_"
                        + mod_TEMPresp
                        + "."
                        + prd[sim]
                        + "_(3var)."
                        + sim
                        + "_global.csv",
                        "r",
                    )
                )][1:],
                dtype=dty,
            )
            lgd = [
                line
                for line in csv.reader(
                    open(
                        "data/ClimReg_CMIP5/#DATA.ClimReg_"
                        + mod_TEMPresp
                        + "."
                        + prd[sim]
                        + "_(3var)."
                        + sim
                        + "_global.csv",
                        "r",
                    )
                )][0]
            exec("gst_" + sim + 'TR = TMP[:,lgd.index("tas")]')
        # local
        for var in ["tos"]:
            if os.path.isfile(
                    "data/ClimReg_CMIP5/#DATA.ClimReg_"
                    + mod_TEMPresp
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
                            "data/ClimReg_CMIP5/#DATA.ClimReg_"
                            + mod_TEMPresp
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
                            "data/ClimReg_CMIP5/#DATA.ClimReg_" + mod_TEMPresp + ".251yr_18x10lat.SURF.csv",
                            "r")
                    )],
                    dtype=dty,
                )
                exec(
                    var + "_" + sim + "TR = np.sum(TMP*TMP2[:len(TMP)],1)/np.sum(TMP2[:len(TMP)],1)")

# aggregate all experiments
for VAR in ["gst", "tos"]:
    if mod_TEMPpattern == "4xCO2":
        exec(VAR + "_allTR = np.array(list(" + VAR + "_quadTR), dtype=dty)")
    elif mod_TEMPpattern == "hist&RCPs":
        exec(VAR + "_allTR = []")
        for sim in ["rcp26", "rcp45", "rcp60", "rcp85"]:
            if os.path.isfile(
                    "data/ClimReg_CMIP5/#DATA.ClimReg_" + mod_TEMPresp + "." + prd[
                        sim] + "_(3var)." + sim + "_global.csv"
            ):
                exec(
                    VAR + "_allTR += list(" + VAR + "_histTR)+list(" + VAR + "_" + sim + "TR)")
        exec("test = " + VAR + "_allTR")
        if test == []:
            exec(VAR + "_allTR += list(" + VAR + "_histTR)")
        exec(VAR + "_allTR = np.array(" + VAR + "_allTR, dtype=dty)")

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
        D_pH = D_pH = (
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
for var in ["lst", "lyp"]:
    exec(var + "_cru = np.zeros([314+1,nb_regionI], dtype=dty)")
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open("data/HistLand_CRU/#DATA.HistLand_CRU.1901-2014_114reg1." + var + ".csv",
                 "r"))],
        dtype=dty,
    )
    TMP2 = np.array(
        [line for line in csv.reader(
            open("data/HistLand_CRU/#DATA.HistLand_CRU.1901-2014_114reg1.AREA.csv",
                 "r"))],
        dtype=dty,
    )
    for i in range(1, 114 + 1):
        exec(var + "_cru[201:,regionI_index[i]] += TMP[:,i-1]*TMP2[:,i-1]")
    TMP = np.zeros([314 + 1, nb_regionI], dtype=dty)
    for i in range(1, 114 + 1):
        TMP[201:, regionI_index[i]] += TMP2[:, i - 1]
    exec(var + "_cru /= TMP[:,:]")

# historical local climate from GHCN+CAMS {degC}
# from [Fan and van den Dool, 2008] updated from website
for var in ["lst"]:
    exec(var + "_ghcn = np.zeros([314+1,nb_regionI], dtype=dty)")
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/HistLand_GHCN-CAMS/#DATA.HistLand_GHCN-CAMS.1948-2014_114reg1." + var + ".csv",
                "r")
        )],
        dtype=dty,
    )
    TMP2 = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/HistLand_GHCN-CAMS/#DATA.HistLand_GHCN-CAMS.1948-2014_114reg1.AREA.csv",
                "r")
        )],
        dtype=dty,
    )
    for i in range(1, 114 + 1):
        exec(var + "_ghcn[248:,regionI_index[i]] += TMP[:,i-1]*TMP2[:,i-1]")
    TMP = np.zeros([314 + 1, nb_regionI], dtype=dty)
    for i in range(1, 114 + 1):
        TMP[248:, regionI_index[i]] += TMP2[:, i - 1]
    exec(var + "_ghcn /= TMP[:,:]")

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
    for sim in ["ctrl", "quad"]:
        # global
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/Climate_CMIP5/#DATA.Climate_"
                    + mod_TEMPresp
                    + "."
                    + str(lng[mod_TEMPresp])
                    + "yr_(7var)."
                    + sim
                    + "_global.csv",
                    "r",
                )
            )][1:],
            dtype=dty,
        )
        lgd = [
            line
            for line in csv.reader(
                open(
                    "data/Climate_CMIP5/#DATA.Climate_"
                    + mod_TEMPresp
                    + "."
                    + str(lng[mod_TEMPresp])
                    + "yr_(7var)."
                    + sim
                    + "_global.csv",
                    "r",
                )
            )][0]
        exec("gst_" + sim + 'TR = TMP[:,lgd.index("tas")]')
        # local
        for var in ["tas"]:
            exec(
                var + "_" + sim + "TR = np.zeros([lng[mod_TEMPresp],nb_regionI], dtype=dty)")
            TMP = np.array(
                [
                    line
                    for line in csv.reader(
                    open(
                        "data/Climate_CMIP5/#DATA.Climate_"
                        + mod_TEMPresp
                        + "."
                        + str(lng[mod_TEMPresp])
                        + "yr_114reg1."
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
                        "data/Climate_CMIP5/#DATA.Climate_"
                        + mod_TEMPresp
                        + "."
                        + str(lng[mod_TEMPresp])
                        + "yr_114reg1.AREA.csv",
                        "r",
                    )
                )],
                dtype=dty,
            )
            exec(var + "_" + sim + "TR = np.zeros([len(TMP),nb_regionI], dtype=dty)")
            exec("AREA_" + sim + "TR = np.zeros([len(TMP),nb_regionI], dtype=dty)")
            for i in range(1, 114 + 1):
                exec(
                    var + "_" + sim + "TR[:,regionI_index[i]] += TMP[:,i-1]*TMP2[:len(TMP),i-1]")
                exec("AREA_" + sim + "TR[:,regionI_index[i]] += TMP2[:len(TMP),i-1]")
            exec(var + "_" + sim + "TR /= AREA_" + sim + "TR")
            exec(
                var + "_" + sim + "TR[np.isnan(" + var + "_" + sim + "TR)|np.isinf(" + var + "_" + sim + "TR)] = 0")

# temperature patterns based on 'historical' and 'rcp'
elif mod_TEMPpattern == "hist&RCPs":
    for sim in ["ctrl", "hist", "rcp26", "rcp45", "rcp60", "rcp85"]:
        # global
        if os.path.isfile(
                "data/ClimReg_CMIP5/#DATA.ClimReg_" + mod_TEMPresp + "." + prd[
                    sim] + "_(3var)." + sim + "_global.csv"
        ):
            TMP = np.array(
                [
                    line
                    for line in csv.reader(
                    open(
                        "data/ClimReg_CMIP5/#DATA.ClimReg_"
                        + mod_TEMPresp
                        + "."
                        + prd[sim]
                        + "_(3var)."
                        + sim
                        + "_global.csv",
                        "r",
                    )
                )][1:],
                dtype=dty,
            )
            lgd = [
                line
                for line in csv.reader(
                    open(
                        "data/ClimReg_CMIP5/#DATA.ClimReg_"
                        + mod_TEMPresp
                        + "."
                        + prd[sim]
                        + "_(3var)."
                        + sim
                        + "_global.csv",
                        "r",
                    )
                )][0]
            exec("gst_" + sim + 'TR = TMP[:,lgd.index("tas")]')
        # local
        for var in ["tas"]:
            if os.path.isfile(
                    "data/ClimReg_CMIP5/#DATA.ClimReg_"
                    + mod_TEMPresp
                    + "."
                    + prd[sim]
                    + "_114reg1."
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
                            "data/ClimReg_CMIP5/#DATA.ClimReg_"
                            + mod_TEMPresp
                            + "."
                            + prd[sim]
                            + "_114reg1."
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
                            "data/ClimReg_CMIP5/#DATA.ClimReg_" + mod_TEMPresp + ".251yr_114reg1.AREA.csv",
                            "r")
                    )],
                    dtype=dty,
                )
                exec(var + "_" + sim + "TR = np.zeros([len(TMP),nb_regionI], dtype=dty)")
                exec("AREA_" + sim + "TR = np.zeros([len(TMP),nb_regionI], dtype=dty)")
                for i in range(1, 114 + 1):
                    exec(
                        var + "_" + sim + "TR[:,regionI_index[i]] += TMP[:,i-1]*TMP2[:len(TMP),i-1]")
                    exec("AREA_" + sim + "TR[:,regionI_index[i]] += TMP2[:len(TMP),i-1]")
                exec(var + "_" + sim + "TR /= AREA_" + sim + "TR")
                exec(
                    var + "_" + sim + "TR[np.isnan(" + var + "_" + sim + "TR)|np.isinf(" + var + "_" + sim + "TR)] = 0"
                )

# aggregate all experiments
for VAR in ["gst", "tas"]:
    if mod_TEMPpattern == "4xCO2":
        exec(VAR + "_allTR = np.array(list(" + VAR + "_quadTR), dtype=dty)")
    elif mod_TEMPpattern == "hist&RCPs":
        exec(VAR + "_allTR = []")
        nrcp = 0
        for sim in ["rcp26", "rcp45", "rcp60", "rcp85"]:
            if os.path.isfile(
                    "data/ClimReg_CMIP5/#DATA.ClimReg_" + mod_TEMPresp + "." + prd[
                        sim] + "_(3var)." + sim + "_global.csv"
            ):
                nrcp += 1
                exec(
                    VAR + "_allTR += list(" + VAR + "_histTR)+list(" + VAR + "_" + sim + "TR)")
        exec("test = " + VAR + "_allTR")
        if test == []:
            exec(VAR + "_allTR += list(" + VAR + "_histTR)")
        exec(VAR + "_allTR = np.array(" + VAR + "_allTR, dtype=dty)")

# decadal means
for VAR in ["gst", "tas"]:
    if mod_TEMPpattern == "4xCO2":
        exec(
            VAR + "_decTR= np.zeros([lng[mod_TEMPresp]-10]+list(np.shape(" + VAR + "_allTR)[1:]), dtype=dty)")
        for t in range(lng[mod_TEMPresp] - 10):
            exec(VAR + "_decTR[t,...] = np.mean(" + VAR + "_allTR[t:t+10,...],0)")
    elif mod_TEMPpattern == "hist&RCPs":
        exec(
            VAR + "_decTR= np.zeros([(156-10)+nrcp*(95-10)]+list(np.shape(" + VAR + "_allTR)[1:]), dtype=dty)")
        for t in range(156 - 10):
            exec(VAR + "_decTR[t,...] = np.mean(" + VAR + "_allTR[t:t+10,...],0)")
        for n in range(nrcp):
            for t in range(95 - 10):
                exec(
                    VAR
                    + "_decTR[(156-10)+n*(95-10)+t,...] = np.mean("
                    + VAR
                    + "_allTR[(156-10)+n*(95-10)+t:(156-10)+n*(95-10)+t+10,...],0)"
                )

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
    for sim in ["ctrl", "quad"]:
        # global
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/Climate_CMIP5/#DATA.Climate_"
                    + mod_PRECresp
                    + "."
                    + str(lng[mod_PRECresp])
                    + "yr_(7var)."
                    + sim
                    + "_global.csv",
                    "r",
                )
            )][1:],
            dtype=dty,
        )
        lgd = [
            line
            for line in csv.reader(
                open(
                    "data/Climate_CMIP5/#DATA.Climate_"
                    + mod_PRECresp
                    + "."
                    + str(lng[mod_PRECresp])
                    + "yr_(7var)."
                    + sim
                    + "_global.csv",
                    "r",
                )
            )][0]
        exec("gyp_" + sim + 'PR = TMP[:,lgd.index("pr")]')
        # local
        for var in ["pr"]:
            TMP = np.array(
                [
                    line
                    for line in csv.reader(
                    open(
                        "data/Climate_CMIP5/#DATA.Climate_"
                        + mod_PRECresp
                        + "."
                        + str(lng[mod_PRECresp])
                        + "yr_114reg1."
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
                        "data/Climate_CMIP5/#DATA.Climate_"
                        + mod_PRECresp
                        + "."
                        + str(lng[mod_PRECresp])
                        + "yr_114reg1.AREA.csv",
                        "r",
                    )
                )],
                dtype=dty,
            )
            exec(var + "_" + sim + "PR = np.zeros([len(TMP),nb_regionI], dtype=dty)")
            exec("AREA_" + sim + "PR = np.zeros([len(TMP),nb_regionI], dtype=dty)")
            for i in range(1, 114 + 1):
                exec(
                    var + "_" + sim + "PR[:,regionI_index[i]] += TMP[:,i-1]*TMP2[:len(TMP),i-1]")
                exec("AREA_" + sim + "PR[:,regionI_index[i]] += TMP2[:len(TMP),i-1]")
            exec(var + "_" + sim + "PR /= AREA_" + sim + "PR")
            exec(
                var + "_" + sim + "PR[np.isnan(" + var + "_" + sim + "PR)|np.isinf(" + var + "_" + sim + "PR)] = 0")

# precipitation patterns based on 'historical' and 'rcp'
elif mod_PRECpattern == "hist&RCPs":
    for sim in ["ctrl", "hist", "rcp26", "rcp45", "rcp60", "rcp85"]:
        # global
        if os.path.isfile(
                "data/ClimReg_CMIP5/#DATA.ClimReg_" + mod_PRECresp + "." + prd[
                    sim] + "_(3var)." + sim + "_global.csv"
        ):
            TMP = np.array(
                [
                    line
                    for line in csv.reader(
                    open(
                        "data/ClimReg_CMIP5/#DATA.ClimReg_"
                        + mod_PRECresp
                        + "."
                        + prd[sim]
                        + "_(3var)."
                        + sim
                        + "_global.csv",
                        "r",
                    )
                )][1:],
                dtype=dty,
            )
            lgd = [
                line
                for line in csv.reader(
                    open(
                        "data/ClimReg_CMIP5/#DATA.ClimReg_"
                        + mod_PRECresp
                        + "."
                        + prd[sim]
                        + "_(3var)."
                        + sim
                        + "_global.csv",
                        "r",
                    )
                )][0]
            exec("gyp_" + sim + 'PR = TMP[:,lgd.index("pr")]')
        # local
        for var in ["pr"]:
            if os.path.isfile(
                    "data/ClimReg_CMIP5/#DATA.ClimReg_"
                    + mod_PRECresp
                    + "."
                    + prd[sim]
                    + "_114reg1."
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
                            "data/ClimReg_CMIP5/#DATA.ClimReg_"
                            + mod_PRECresp
                            + "."
                            + prd[sim]
                            + "_114reg1."
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
                            "data/ClimReg_CMIP5/#DATA.ClimReg_" + mod_PRECresp + ".251yr_114reg1.AREA.csv",
                            "r")
                    )],
                    dtype=dty,
                )
                exec(var + "_" + sim + "PR = np.zeros([len(TMP),nb_regionI], dtype=dty)")
                exec("AREA_" + sim + "PR = np.zeros([len(TMP),nb_regionI], dtype=dty)")
                for i in range(1, 114 + 1):
                    exec(
                        var + "_" + sim + "PR[:,regionI_index[i]] += TMP[:,i-1]*TMP2[:len(TMP),i-1]")
                    exec("AREA_" + sim + "PR[:,regionI_index[i]] += TMP2[:len(TMP),i-1]")
                exec(var + "_" + sim + "PR /= AREA_" + sim + "PR")
                exec(
                    var + "_" + sim + "PR[np.isnan(" + var + "_" + sim + "PR)|np.isinf(" + var + "_" + sim + "PR)] = 0"
                )

# aggregate all experiments
for VAR in ["gyp", "pr"]:
    if mod_PRECpattern == "4xCO2":
        exec(VAR + "_allPR = np.array(list(" + VAR + "_quadPR), dtype=dty)")
    elif mod_PRECpattern == "hist&RCPs":
        exec(VAR + "_allPR = []")
        nrcp = 0
        for sim in ["rcp26", "rcp45", "rcp60", "rcp85"]:
            if os.path.isfile(
                    "data/ClimReg_CMIP5/#DATA.ClimReg_" + mod_PRECresp + "." + prd[
                        sim] + "_(3var)." + sim + "_global.csv"
            ):
                nrcp += 1
                exec(
                    VAR + "_allPR += list(" + VAR + "_histPR)+list(" + VAR + "_" + sim + "PR)")
        exec("test = " + VAR + "_allPR")
        if test == []:
            exec(VAR + "_allPR += list(" + VAR + "_histPR)")
        exec(VAR + "_allPR = np.array(" + VAR + "_allPR, dtype=dty)")
# decadal means
for VAR in ["gyp", "pr"]:
    if mod_PRECpattern == "4xCO2":
        exec(
            VAR + "_decPR= np.zeros([lng[mod_PRECresp]-10]+list(np.shape(" + VAR + "_allPR)[1:]), dtype=dty)")
        for t in range(lng[mod_PRECresp] - 10):
            exec(VAR + "_decPR[t,...] = np.mean(" + VAR + "_allPR[t:t+10,...],0)")
    elif mod_PRECpattern == "hist&RCPs":
        exec(
            VAR + "_decPR= np.zeros([(156-10)+nrcp*(95-10)]+list(np.shape(" + VAR + "_allPR)[1:]), dtype=dty)")
        for t in range(156 - 10):
            exec(VAR + "_decPR[t,...] = np.mean(" + VAR + "_allPR[t:t+10,...],0)")
        for n in range(nrcp):
            for t in range(95 - 10):
                exec(
                    VAR
                    + "_decPR[(156-10)+n*(95-10)+t,...] = np.mean("
                    + VAR
                    + "_allPR[(156-10)+n*(95-10)+t:(156-10)+n*(95-10)+t+10,...],0)"
                )

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
