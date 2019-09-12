import csv

import numpy as np
from scipy.optimize import fmin

from ..oscar_data import *
from ...config import dty, mod_POAconv, mod_SO4regsat, mod_POAregsat, mod_BCregsat, \
    mod_POAload, mod_BCload, mod_NO3load, mod_SOAload, mod_DUSTload, mod_SALTload, \
    mod_SO4load

##################################################
#   6. AEROSOLS
##################################################

# ===============
# 6.1. ATMOSPHERE
# ===============

# conversion of SO4 from {TgS} to {Tg(SO4)}
alpha_SO4 = np.array([96 / 32.0], dtype=dty)
# conversion of POM from {Tg(OC)} to {Tg(OM)}
if mod_POAconv == "default":
    alpha_POM = np.array([1.4], dtype=dty)
elif mod_POAconv == "GFDL":
    alpha_POM = np.array([1.6], dtype=dty)
elif mod_POAconv == "CSIRO":
    alpha_POM = np.array([1.3], dtype=dty)
# conversion of NO3 from {TgN} to {Tg(NO3)}
alpha_NO3 = np.array([62 / 14.0], dtype=dty)

# ==============
# 6.2. CHEMISTRY
# ==============

# -----------
# 6.2.1. HTAP
# -----------

# regional saturation coefficient {.}
# from HTAP experiments [Yu et al., 2013] (table 6 extended; provided by author)
# sulfate
if mod_SO4regsat == "mean-HTAP":
    w_reg_SO2 = np.array([-3.51, -3.87, -3.92, -2.85, -3.92], dtype=dty) / -3.51
elif mod_SO4regsat == "CAMCHEM":
    w_reg_SO2 = np.array([-4.26, -5.01, -4.86, -3.29, -4.55], dtype=dty) / -4.26
elif mod_SO4regsat == "GISS-PUCCINI":
    w_reg_SO2 = np.array([-2.73, -3.88, -2.92, -1.78, -3.44], dtype=dty) / -2.73
elif mod_SO4regsat == "GMI":
    w_reg_SO2 = np.array([-3.61, -3.55, -3.96, -3.30, -3.59], dtype=dty) / -3.61
elif mod_SO4regsat == "GOCART":
    w_reg_SO2 = np.array([-4.05, -3.86, -4.68, -3.79, -3.06], dtype=dty) / -4.05
elif mod_SO4regsat == "INCA2":
    w_reg_SO2 = np.array([-4.03, -4.52, -4.27, -3.33, -5.45], dtype=dty) / -4.03
elif mod_SO4regsat == "LLNL-IMPACT":
    w_reg_SO2 = np.array([-3.70, -3.74, -4.27, -2.84, -4.82], dtype=dty) / -3.70
elif mod_SO4regsat == "SPRINTARS":
    w_reg_SO2 = np.array([-2.16, -2.51, -2.53, -1.64, -2.56], dtype=dty) / -2.16
else:
    w_reg_SO2 = np.array([1, 1, 1, 1, 1], dtype=dty)

# primary organic aerosols
if mod_POAregsat == "mean-HTAP":
    w_reg_OC = np.array([-4.00, -4.39, -4.30, -3.68, -4.09], dtype=dty) / -4.00
elif mod_POAregsat == "CAMCHEM":
    w_reg_OC = np.array([-3.30, -3.90, -3.23, -2.80, -3.71], dtype=dty) / -3.30
elif mod_POAregsat == "GISS-PUCCINI":
    w_reg_OC = np.array([-6.26, -5.86, -6.07, -6.44, -6.38], dtype=dty) / -6.26
elif mod_POAregsat == "GMI":
    w_reg_OC = np.array([-4.62, -5.30, -6.13, -4.22, -4.16], dtype=dty) / -4.62
elif mod_POAregsat == "GOCART":
    w_reg_OC = np.array([-3.57, -3.65, -3.39, -3.28, -4.00], dtype=dty) / -3.57
elif mod_POAregsat == "INCA2":
    w_reg_OC = np.array([-3.07, -4.27, -3.33, -2.88, -2.66], dtype=dty) / -3.07
elif mod_POAregsat == "LLNL-IMPACT":
    w_reg_OC = np.array([-1.31, -1.41, -1.97, -0.99, -1.29], dtype=dty) / -1.31
elif mod_POAregsat == "SPRINTARS":
    w_reg_OC = np.array([-5.86, -6.32, -5.97, -5.12, -6.45], dtype=dty) / -5.86
else:
    w_reg_OC = np.array([1, 1, 1, 1, 1], dtype=dty)

# black carbon
if mod_BCregsat == "mean-HTAP":
    w_reg_BC = np.array([29.51, 27.31, 37.36, 28.36, 25.31], dtype=dty) / 29.51
elif mod_BCregsat == "CAMCHEM":
    w_reg_BC = np.array([27.56, 28.00, 35.71, 25.24, 24.08], dtype=dty) / 27.56
elif mod_BCregsat == "GISS-PUCCINI":
    w_reg_BC = np.array([60.41, 51.67, 69.53, 65.06, 45.46], dtype=dty) / 60.41
elif mod_BCregsat == "GMI":
    w_reg_BC = np.array([26.68, 25.80, 42.13, 24.81, 15.00], dtype=dty) / 26.68
elif mod_BCregsat == "GOCART":
    w_reg_BC = np.array([46.20, 42.30, 52.21, 45.87, 43.71], dtype=dty) / 46.20
elif mod_BCregsat == "INCA2":
    w_reg_BC = np.array([17.32, 16.88, 20.37, 14.14, 23.69], dtype=dty) / 17.32
elif mod_BCregsat == "LLNL-IMPACT":
    w_reg_BC = np.array([7.25, 6.63, 12.99, 5.66, 5.56], dtype=dty) / 7.25
elif mod_BCregsat == "SPRINTARS":
    w_reg_BC = np.array([21.16, 19.93, 28.58, 17.75, 19.69], dtype=dty) / 21.16
else:
    w_reg_BC = np.array([1, 1, 1, 1, 1], dtype=dty)

# -------------
# 6.2.2. ACCMIP
# -------------

# period of ACCMIP/CMIP5 simulations per simulation
prd = {
    "ctrl": "251yr",
    "hist": "1850-2005",
    "rcp26": "2006-2100",
    "rcp45": "2006-2100",
    "rcp60": "2006-2100",
    "rcp85": "2006-2100",
}

# load pre-processed ACCMIP/CMIP5 results for specified model
# sulfate
for sim in ["hist", "rcp26", "rcp45", "rcp60", "rcp85"]:
    for VAR in ["ESO2", "EDMS", "SO4"] + ["tas"]:
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/AeroChem_ACCMIP/#DATA.AeroChem_"
                    + mod_SO4load
                    + "."
                    + prd[sim]
                    + "."
                    + sim
                    + "_"
                    + VAR
                    + ".csv",
                    "r",
                )
            )],
            dtype=dty,
        )
        exec(VAR + "_" + sim + "S = TMP[:,0].copy()")
for VAR in ["ESO2", "EDMS", "SO4"] + ["tas"]:
    exec(
        VAR
        + "_allS = np.array(list("
        + VAR
        + "_histS)+list("
        + VAR
        + "_histS)+list("
        + VAR
        + "_histS)+list("
        + VAR
        + "_histS)+list("
        + VAR
        + "_rcp26S)+list("
        + VAR
        + "_rcp45S)+list("
        + VAR
        + "_rcp60S)+list("
        + VAR
        + "_rcp85S), dtype=dty)"
    )

# primary organic aerosols
for sim in ["hist", "rcp26", "rcp45", "rcp60", "rcp85"]:
    for VAR in ["EOM", "EOMBB", "POA"] + ["tas"]:
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/AeroChem_ACCMIP/#DATA.AeroChem_"
                    + mod_POAload
                    + "."
                    + prd[sim]
                    + "."
                    + sim
                    + "_"
                    + VAR
                    + ".csv",
                    "r",
                )
            )],
            dtype=dty,
        )
        exec(VAR + "_" + sim + "P = TMP[:,0].copy()")
for VAR in ["EOM", "EOMBB", "POA"] + ["tas"]:
    exec(
        VAR
        + "_allP = np.array(list("
        + VAR
        + "_histP)+list("
        + VAR
        + "_histP)+list("
        + VAR
        + "_histP)+list("
        + VAR
        + "_histP)+list("
        + VAR
        + "_rcp26P)+list("
        + VAR
        + "_rcp45P)+list("
        + VAR
        + "_rcp60P)+list("
        + VAR
        + "_rcp85P), dtype=dty)"
    )

# black carbon
for sim in ["hist", "rcp26", "rcp45", "rcp60", "rcp85"]:
    for VAR in ["EBC", "EBCBB", "BC"] + ["tas"]:
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/AeroChem_ACCMIP/#DATA.AeroChem_"
                    + mod_BCload
                    + "."
                    + prd[sim]
                    + "."
                    + sim
                    + "_"
                    + VAR
                    + ".csv",
                    "r",
                )
            )],
            dtype=dty,
        )
        exec(VAR + "_" + sim + "B = TMP[:,0].copy()")
for VAR in ["EBC", "EBCBB", "BC"] + ["tas"]:
    exec(
        VAR
        + "_allB = np.array(list("
        + VAR
        + "_histB)+list("
        + VAR
        + "_histB)+list("
        + VAR
        + "_histB)+list("
        + VAR
        + "_histB)+list("
        + VAR
        + "_rcp26B)+list("
        + VAR
        + "_rcp45B)+list("
        + VAR
        + "_rcp60B)+list("
        + VAR
        + "_rcp85B), dtype=dty)"
    )

# nitrate
if not mod_NO3load in ["Bellouin2011", "Hauglustaine2014"]:
    for sim in ["hist", "rcp26", "rcp45", "rcp60", "rcp85"]:
        for VAR in ["ENOX", "ENH3", "NO3"] + ["tas2"]:
            TMP = np.array(
                [
                    line
                    for line in csv.reader(
                    open(
                        "data/AeroChem_ACCMIP/#DATA.AeroChem_"
                        + mod_NO3load
                        + "."
                        + prd[sim]
                        + "."
                        + sim
                        + "_"
                        + VAR
                        + ".csv",
                        "r",
                    )
                )],
                dtype=dty,
            )
            exec(VAR + "_" + sim + "N = TMP[:,0].copy()")
    for VAR in ["ENOX", "ENH3", "NO3"] + ["tas2"]:
        exec(
            VAR
            + "_allN = np.array(list("
            + VAR
            + "_histN)+list("
            + VAR
            + "_histN)+list("
            + VAR
            + "_histN)+list("
            + VAR
            + "_histN)+list("
            + VAR
            + "_rcp26N)+list("
            + VAR
            + "_rcp45N)+list("
            + VAR
            + "_rcp60N)+list("
            + VAR
            + "_rcp85N), dtype=dty)"
        )

# secondary organic aerosols
if mod_SOAload != "":
    for sim in ["hist", "rcp26", "rcp45", "rcp60", "rcp85"]:
        for VAR in ["EVOC", "EBVOC", "SOA"] + ["tas2"]:
            TMP = np.array(
                [
                    line
                    for line in csv.reader(
                    open(
                        "data/AeroChem_ACCMIP/#DATA.AeroChem_"
                        + mod_SOAload
                        + "."
                        + prd[sim]
                        + "."
                        + sim
                        + "_"
                        + VAR
                        + ".csv",
                        "r",
                    )
                )],
                dtype=dty,
            )
            exec(VAR + "_" + sim + "Q = TMP[:,0].copy()")
    for VAR in ["EVOC", "EBVOC", "SOA"] + ["tas2"]:
        exec(
            VAR
            + "_allQ = np.array(list("
            + VAR
            + "_histQ)+list("
            + VAR
            + "_histQ)+list("
            + VAR
            + "_histQ)+list("
            + VAR
            + "_histQ)+list("
            + VAR
            + "_rcp26Q)+list("
            + VAR
            + "_rcp45Q)+list("
            + VAR
            + "_rcp60Q)+list("
            + VAR
            + "_rcp85Q), dtype=dty)"
        )

# mineral dusts
for sim in ["hist", "rcp26", "rcp45", "rcp60", "rcp85"]:
    for VAR in ["EDUST", "DUST"] + ["tas"]:
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/AeroChem_ACCMIP/#DATA.AeroChem_"
                    + mod_DUSTload
                    + "."
                    + prd[sim]
                    + "."
                    + sim
                    + "_"
                    + VAR
                    + ".csv",
                    "r",
                )
            )],
            dtype=dty,
        )
        exec(VAR + "_" + sim + "D = TMP[:,0].copy()")
for VAR in ["EDUST", "DUST"] + ["tas"]:
    exec(
        VAR
        + "_allD = np.array(list("
        + VAR
        + "_histD)+list("
        + VAR
        + "_histD)+list("
        + VAR
        + "_histD)+list("
        + VAR
        + "_histD)+list("
        + VAR
        + "_rcp26D)+list("
        + VAR
        + "_rcp45D)+list("
        + VAR
        + "_rcp60D)+list("
        + VAR
        + "_rcp85D), dtype=dty)"
    )

# sea salts
for sim in ["hist", "rcp26", "rcp45", "rcp60", "rcp85"]:
    for VAR in ["ESALT", "SALT"] + ["tas3"]:
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/AeroChem_ACCMIP/#DATA.AeroChem_"
                    + mod_SALTload
                    + "."
                    + prd[sim]
                    + "."
                    + sim
                    + "_"
                    + VAR
                    + ".csv",
                    "r",
                )
            )],
            dtype=dty,
        )
        exec(VAR + "_" + sim + "T = TMP[:,0].copy()")
for VAR in ["ESALT", "SALT"] + ["tas3"]:
    exec(
        VAR
        + "_allT = np.array(list("
        + VAR
        + "_histT)+list("
        + VAR
        + "_histT)+list("
        + VAR
        + "_histT)+list("
        + VAR
        + "_histT)+list("
        + VAR
        + "_rcp26T)+list("
        + VAR
        + "_rcp45T)+list("
        + VAR
        + "_rcp60T)+list("
        + VAR
        + "_rcp85T), dtype=dty)"
    )

# definition of parameters
# lifetimes of sulfate precursors {yr}
tau_SO2 = np.array([0], dtype=dty)
tau_DMS = np.array([0], dtype=dty)
# sensitivity of sulfate to climate change {Tg/K}
Gamma_SO4 = np.array([0], dtype=dty)
# lifetimes of primary organic aerosols {yr}
tau_OMff = np.array([0], dtype=dty)
tau_OMbb = np.array([0], dtype=dty)
# sensitivity of primary organic aerosols to climate change {Tg/K}
Gamma_POA = np.array([0], dtype=dty)
# lifetimes of black carbon {yr}
tau_BCff = np.array([0], dtype=dty)
tau_BCbb = np.array([0], dtype=dty)
# sensitivity of black carbon to climate change {Tg/K}
Gamma_BC = np.array([0], dtype=dty)
# lifetimes of nitrate precursors {yr}
tau_NOX = np.array([0], dtype=dty)
tau_NH3 = np.array([0], dtype=dty)
# sensitivity of nitrate to climate change {Tg/K}
Gamma_NO3 = np.array([0], dtype=dty)
# lifetimes of secondary organic aerosols {yr}
tau_VOC = np.array([0], dtype=dty)
tau_BVOC = np.array([0], dtype=dty)
# sensitivity of secondary organic aerosols to climate change {Tg/K}
Gamma_SOA = np.array([0], dtype=dty)
# lifetime of mineral dusts {yr}
tau_DUST = np.array([0], dtype=dty)
# sensitivity of mineral dusts to climate change {Tg/K}
Gamma_DUST = np.array([0], dtype=dty)
# lifetime of sea salts {yr}
tau_SALT = np.array([0], dtype=dty)
# sensitivity of sea salts to climate change {Tg/K}
Gamma_SALT = np.array([0], dtype=dty)

# fit of parameters
# sulfate
diff = SO4_allS - np.mean(SO4_allS[:10])


def err(var):
    conc = np.abs(var[0]) * alpha_SO4 * (ESO2_allS - np.mean(ESO2_allS[:10])) + np.abs(
        var[1]) * alpha_SO4 * (
                   EDMS_allS - np.mean(EDMS_allS[:10])
           )
    clim = var[2] * (tas_allS - np.mean(tas_allS[:10]))
    return np.sum((diff - (conc + clim)) ** 2)


[tau_SO2[0], tau_DMS[0], Gamma_SO4[0]] = fmin(err, [0.01, 0.01, 0], disp=False)
tau_SO2 = np.abs(tau_SO2)
tau_DMS = np.abs(tau_DMS)

# primary organic aerosols
diff = POA_allP - np.mean(POA_allP[:10])


def err(var):
    conc = np.abs(var[0]) * (EOM_allP - np.mean(EOM_allP[:10])) + np.abs(var[1]) * (
            EOMBB_allP - np.mean(EOMBB_allP[:10])
    )
    clim = var[2] * (tas_allP - np.mean(tas_allP[:10]))
    return np.sum((diff - (conc + clim)) ** 2)


[tau_OMff[0], tau_OMbb[0], Gamma_POA[0]] = fmin(err, [0.01, 0.01, 0], disp=False)
tau_OMff = np.abs(tau_OMff)
tau_OMbb = np.abs(tau_OMbb)

# black carbon
diff = BC_allB - np.mean(BC_allB[:10])


def err(var):
    conc = np.abs(var[0]) * (EBC_allB - np.mean(EBC_allB[:10])) + np.abs(var[1]) * (
            EBCBB_allB - np.mean(EBCBB_allB[:10])
    )
    clim = var[2] * (tas_allB - np.mean(tas_allB[:10]))
    return np.sum((diff - (conc + clim)) ** 2)


[tau_BCff[0], tau_BCbb[0], Gamma_BC[0]] = fmin(err, [0.01, 0.01, 0], disp=False)
tau_BCff = np.abs(tau_BCff)
tau_BCbb = np.abs(tau_BCbb)

# nitrate
if not mod_NO3load in ["Bellouin2011", "Hauglustaine2014"]:
    diff = NO3_allN - np.mean(NO3_allN[:10])


    def err(var):
        conc = np.abs(var[0]) * alpha_NO3 * (
                ENOX_allN - np.mean(ENOX_allN[:10])) + np.abs(var[1]) * alpha_NO3 * (
                       ENH3_allN - np.mean(ENH3_allN[:10])
               )
        clim = var[2] * (tas2_allN - np.mean(tas2_allN[:10]))
        return np.sum((diff - (conc + clim)) ** 2)


    [tau_NOX[0], tau_NH3[0], Gamma_NO3[0]] = fmin(err, [0.01, 0.01, 0], disp=False)
    tau_NOX = np.abs(tau_NOX)
    tau_NH3 = np.abs(tau_NH3)

# secondary organic aerosols
if mod_SOAload != "":
    diff = SOA_allQ - np.mean(SOA_allQ[:10])


    def err(var):
        conc = np.abs(var[0]) * (EVOC_allQ - np.mean(EVOC_allQ[:10])) + np.abs(var[1]) * (
                EBVOC_allQ - np.mean(EBVOC_allQ[:10])
        )
        clim = var[2] * (tas2_allQ - np.mean(tas2_allQ[:10]))
        return np.sum((diff - (conc + clim)) ** 2)


    [tau_VOC[0], tau_BVOC[0], Gamma_SOA[0]] = fmin(err, [0.01, 0.01, 0], disp=False)
    tau_VOC = np.abs(tau_VOC)
    tau_BVOC = np.abs(tau_BVOC)

# mineral dusts
diff = DUST_allD - np.mean(DUST_allD[:10])


def err(var):
    conc = np.abs(var[0]) * (EDUST_allD - np.mean(EDUST_allD[:10]))
    clim = var[1] * (tas_allD - np.mean(tas_allD[:10]))
    return np.sum((diff - (conc + clim)) ** 2)


[tau_DUST[0], Gamma_DUST[0]] = fmin(err, [0.01, 0], disp=False)
tau_DUST = np.abs(tau_DUST)

# sea salts
diff = SALT_allT - np.mean(SALT_allT[:10])


def err(var):
    conc = np.abs(var[0]) * (ESALT_allT - np.mean(ESALT_allT[:10]))
    clim = var[1] * (tas3_allT - np.mean(tas3_allT[:10]))
    return np.sum((diff - (conc + clim)) ** 2)


[tau_SALT[0], Gamma_SALT[0]] = fmin(err, [0.01, 0], disp=False)
tau_SALT = np.abs(tau_SALT)

# ---------------
# 6.2.3. Nitrates
# ---------------

# load data for nitrate aerosols
# from HadGEM2 [Bellouin et al., 2011] (also RCP and CMIP5 data)
if mod_NO3load == "Bellouin2011":
    ENOX_nitrate = np.array([5.7, 37.4, 18.4, 16.3, 16.6, 23.8, 18.4, 16.3, 16.6, 23.8],
                            dtype=dty)
    ENH3_nitrate = np.array([16.6, 41.4, 67.2, 49.2, 63.0, 70.0, 67.2, 49.2, 63.0, 70.0],
                            dtype=dty)
    tas_nitrate = np.array(
        [13.55, 14.07, 15.39, 16.46, 17.07, 18.66, 13.55, 13.55, 13.55, 13.55], dtype=dty)
    NO3_nitrate = np.array([0.05, 0.34, 0.56, 0.29, 0.41, 0.52, 0.63, 0.36, 0.50, 0.68],
                           dtype=dty)
# from LMDz4-INCA3 [Hauglustaine et al., 2014] (tables 1,5)
elif mod_NO3load == "Hauglustaine2014":
    ENOX_nitrate = np.array(
        [10, 36, 29, 26, 14, 32, 26, 14, 30, 27, 13, 38, 30, 21, 21, 14], dtype=dty)
    ENH3_nitrate = np.array(
        [21, 29, 41, 46, 58, 35, 36, 33, 36, 43, 51, 42, 48, 57, 33, 57], dtype=dty)
    tas_nitrate = np.array(
        [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
         0.00, 0.00, 0.00], dtype=dty
    )
    # NO3_nitrate = np.array([0.23,0.48,0.46,0.47,0.37,0.48,0.46,0.42,0.47,0.48,0.40,0.54,0.52,0.52,0.47,0.43], dtype=dty) # HNO3 and NO3-
    NO3_nitrate = np.array(
        [0.09, 0.18, 0.21, 0.23, 0.21, 0.2, 0.2, 0.18, 0.19, 0.21, 0.21, 0.23, 0.24, 0.25,
         0.2, 0.22], dtype=dty
    )  # NO3- only

# fit of parameters (defined in previous section)
diff = NO3_nitrate - NO3_nitrate[0]


def err(var):
    conc = np.abs(var[0]) * alpha_NO3 * (ENOX_nitrate - ENOX_nitrate[0]) + np.abs(
        var[1]) * alpha_NO3 * (
                   ENH3_nitrate - ENH3_nitrate[0]
           )
    clim = var[2] * (tas_nitrate - tas_nitrate[0])
    if mod_NO3load == "Hauglustaine2014":
        clim = 0 * clim
    return np.sum((diff - (conc + clim)) ** 2)


[tau_NOX[0], tau_NH3[0], Gamma_NO3[0]] = fmin(err, [0.01, 0.01, 0], disp=False)
tau_NOX = np.abs(tau_NOX)
tau_NH3 = np.abs(tau_NH3)
if mod_NO3load == "Hauglustaine2014":
    Gamma_NO3 = 0 * Gamma_NO3
