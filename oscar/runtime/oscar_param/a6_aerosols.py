import numpy as np
from ...data import load_data
from scipy.optimize import fmin
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
def load_accmip(var, sim):
    path = f"data/AeroChem_ACCMIP/#DATA.AeroChem_{mod_SO4load}.{prd[sim]}.{sim}_{var}.csv"
    return load_data(path).copy()

ESO2_histS = load_accmip('ESO2', 'hist')
EDMS_histS = load_accmip('EDMS', 'hist')
SO4_histS = load_accmip('SO4', 'hist')
tas_histS = load_accmip('tas', 'hist')
ESO2_rcp26S = load_accmip('ESO2', 'rcp26')
EDMS_rcp26S = load_accmip('EDMS', 'rcp26')
SO4_rcp26S = load_accmip('SO4', 'rcp26')
tas_rcp26S = load_accmip('tas', 'rcp26')
ESO2_rcp45S = load_accmip('ESO2', 'rcp45')
EDMS_rcp45S = load_accmip('EDMS', 'rcp45')
SO4_rcp45S = load_accmip('SO4', 'rcp45')
tas_rcp45S = load_accmip('tas', 'rcp45')
ESO2_rcp60S = load_accmip('ESO2', 'rcp60')
EDMS_rcp60S = load_accmip('EDMS', 'rcp60')
SO4_rcp60S = load_accmip('SO4', 'rcp60')
tas_rcp60S = load_accmip('tas', 'rcp60')
ESO2_rcp85S = load_accmip('ESO2', 'rcp85')
EDMS_rcp85S = load_accmip('EDMS', 'rcp85')
SO4_rcp85S = load_accmip('SO4', 'rcp85')
tas_rcp85S = load_accmip('tas', 'rcp85')

ESO2_allS = np.array(list(ESO2_histS)+list(ESO2_histS)+list(ESO2_histS)+list(ESO2_histS)+list(ESO2_rcp26S)+list(ESO2_rcp45S)+list(ESO2_rcp60S)+list(ESO2_rcp85S), dtype=dty)
EDMS_allS = np.array(list(EDMS_histS)+list(EDMS_histS)+list(EDMS_histS)+list(EDMS_histS)+list(EDMS_rcp26S)+list(EDMS_rcp45S)+list(EDMS_rcp60S)+list(EDMS_rcp85S), dtype=dty)
SO4_allS = np.array(list(SO4_histS)+list(SO4_histS)+list(SO4_histS)+list(SO4_histS)+list(SO4_rcp26S)+list(SO4_rcp45S)+list(SO4_rcp60S)+list(SO4_rcp85S), dtype=dty)
tas_allS = np.array(list(tas_histS)+list(tas_histS)+list(tas_histS)+list(tas_histS)+list(tas_rcp26S)+list(tas_rcp45S)+list(tas_rcp60S)+list(tas_rcp85S), dtype=dty)

# primary organic aerosols
def load_organic_aero_chem(sim, VAR):
    path = f"data/AeroChem_ACCMIP/#DATA.AeroChem_{mod_POAload}.{prd[sim]}.{sim}_{VAR}.csv"
    TMP = load_data(path)
    return TMP[:,0].copy()


EOM_histP = load_organic_aero_chem('hist', 'EOM')
EOMBB_histP = load_organic_aero_chem('hist', 'EOMBB')
POA_histP = load_organic_aero_chem('hist', 'POA')
tas_histP = load_organic_aero_chem('hist', 'tas')
EOM_rcp26P = load_organic_aero_chem('rcp26', 'EOM')
EOMBB_rcp26P = load_organic_aero_chem('rcp26', 'EOMBB')
POA_rcp26P = load_organic_aero_chem('rcp26', 'POA')
tas_rcp26P = load_organic_aero_chem('rcp26', 'tas')
EOM_rcp45P = load_organic_aero_chem('rcp45', 'EOM')
EOMBB_rcp45P = load_organic_aero_chem('rcp45', 'EOMBB')
POA_rcp45P = load_organic_aero_chem('rcp45', 'POA')
tas_rcp45P = load_organic_aero_chem('rcp45', 'tas')
EOM_rcp60P = load_organic_aero_chem('rcp60', 'EOM')
EOMBB_rcp60P = load_organic_aero_chem('rcp60', 'EOMBB')
POA_rcp60P = load_organic_aero_chem('rcp60', 'POA')
tas_rcp60P = load_organic_aero_chem('rcp60', 'tas')
EOM_rcp85P = load_organic_aero_chem('rcp85', 'EOM')
EOMBB_rcp85P = load_organic_aero_chem('rcp85', 'EOMBB')
POA_rcp85P = load_organic_aero_chem('rcp85', 'POA')
tas_rcp85P = load_organic_aero_chem('rcp85', 'tas')

EOM_allP = np.array(list(EOM_histP)+list(EOM_histP)+list(EOM_histP)+list(EOM_histP)+list(EOM_rcp26P)+list(EOM_rcp45P)+list(EOM_rcp60P)+list(EOM_rcp85P), dtype=dty)
EOMBB_allP = np.array(list(EOMBB_histP)+list(EOMBB_histP)+list(EOMBB_histP)+list(EOMBB_histP)+list(EOMBB_rcp26P)+list(EOMBB_rcp45P)+list(EOMBB_rcp60P)+list(EOMBB_rcp85P), dtype=dty)
POA_allP = np.array(list(POA_histP)+list(POA_histP)+list(POA_histP)+list(POA_histP)+list(POA_rcp26P)+list(POA_rcp45P)+list(POA_rcp60P)+list(POA_rcp85P), dtype=dty)
tas_allP = np.array(list(tas_histP)+list(tas_histP)+list(tas_histP)+list(tas_histP)+list(tas_rcp26P)+list(tas_rcp45P)+list(tas_rcp60P)+list(tas_rcp85P), dtype=dty)

# black carbon
# for sim in ["hist", "rcp26", "rcp45", "rcp60", "rcp85"]:
#     for VAR in ["EBC", "EBCBB", "BC"] + ["tas"]:
#         print(f'{VAR}_{sim}B = load_black_carbon({sim!r}, {VAR!r})')


def load_black_carbon(sim, VAR):
    path = f"data/AeroChem_ACCMIP/#DATA.AeroChem_{mod_BCload}.{prd[sim]}.{sim}_{VAR}.csv"
    TMP = load_data(path)
    return TMP[:,0].copy()

EBC_histB = load_black_carbon('hist', 'EBC')
EBCBB_histB = load_black_carbon('hist', 'EBCBB')
BC_histB = load_black_carbon('hist', 'BC')
tas_histB = load_black_carbon('hist', 'tas')
EBC_rcp26B = load_black_carbon('rcp26', 'EBC')
EBCBB_rcp26B = load_black_carbon('rcp26', 'EBCBB')
BC_rcp26B = load_black_carbon('rcp26', 'BC')
tas_rcp26B = load_black_carbon('rcp26', 'tas')
EBC_rcp45B = load_black_carbon('rcp45', 'EBC')
EBCBB_rcp45B = load_black_carbon('rcp45', 'EBCBB')
BC_rcp45B = load_black_carbon('rcp45', 'BC')
tas_rcp45B = load_black_carbon('rcp45', 'tas')
EBC_rcp60B = load_black_carbon('rcp60', 'EBC')
EBCBB_rcp60B = load_black_carbon('rcp60', 'EBCBB')
BC_rcp60B = load_black_carbon('rcp60', 'BC')
tas_rcp60B = load_black_carbon('rcp60', 'tas')
EBC_rcp85B = load_black_carbon('rcp85', 'EBC')
EBCBB_rcp85B = load_black_carbon('rcp85', 'EBCBB')
BC_rcp85B = load_black_carbon('rcp85', 'BC')
tas_rcp85B = load_black_carbon('rcp85', 'tas')

EBC_allB = np.array(list(EBC_histB)+list(EBC_histB)+list(EBC_histB)+list(EBC_histB)+list(EBC_rcp26B)+list(EBC_rcp45B)+list(EBC_rcp60B)+list(EBC_rcp85B), dtype=dty)
EBCBB_allB = np.array(list(EBCBB_histB)+list(EBCBB_histB)+list(EBCBB_histB)+list(EBCBB_histB)+list(EBCBB_rcp26B)+list(EBCBB_rcp45B)+list(EBCBB_rcp60B)+list(EBCBB_rcp85B), dtype=dty)
BC_allB = np.array(list(BC_histB)+list(BC_histB)+list(BC_histB)+list(BC_histB)+list(BC_rcp26B)+list(BC_rcp45B)+list(BC_rcp60B)+list(BC_rcp85B), dtype=dty)
tas_allB = np.array(list(tas_histB)+list(tas_histB)+list(tas_histB)+list(tas_histB)+list(tas_rcp26B)+list(tas_rcp45B)+list(tas_rcp60B)+list(tas_rcp85B), dtype=dty)

# nitrate
if not mod_NO3load in ["Bellouin2011", "Hauglustaine2014"]:
    # for sim in ["hist", "rcp26", "rcp45", "rcp60", "rcp85"]:
    #     for VAR in ["ENOX", "ENH3", "NO3"] + ["tas2"]:
    #         print(f'{VAR}_{sim}N = load_nitrate({sim!r}, {VAR!r})')


    def load_nitrate(sim, VAR):
        path = f"data/AeroChem_ACCMIP/#DATA.AeroChem_{mod_NO3load}.{prd[sim]}.{sim}_{VAR}.csv"
        TMP = load_data(path)
        return TMP[:,0].copy()

    ENOX_histN = load_nitrate('hist', 'ENOX')
    ENH3_histN = load_nitrate('hist', 'ENH3')
    NO3_histN = load_nitrate('hist', 'NO3')
    tas2_histN = load_nitrate('hist', 'tas2')
    ENOX_rcp26N = load_nitrate('rcp26', 'ENOX')
    ENH3_rcp26N = load_nitrate('rcp26', 'ENH3')
    NO3_rcp26N = load_nitrate('rcp26', 'NO3')
    tas2_rcp26N = load_nitrate('rcp26', 'tas2')
    ENOX_rcp45N = load_nitrate('rcp45', 'ENOX')
    ENH3_rcp45N = load_nitrate('rcp45', 'ENH3')
    NO3_rcp45N = load_nitrate('rcp45', 'NO3')
    tas2_rcp45N = load_nitrate('rcp45', 'tas2')
    ENOX_rcp60N = load_nitrate('rcp60', 'ENOX')
    ENH3_rcp60N = load_nitrate('rcp60', 'ENH3')
    NO3_rcp60N = load_nitrate('rcp60', 'NO3')
    tas2_rcp60N = load_nitrate('rcp60', 'tas2')
    ENOX_rcp85N = load_nitrate('rcp85', 'ENOX')
    ENH3_rcp85N = load_nitrate('rcp85', 'ENH3')
    NO3_rcp85N = load_nitrate('rcp85', 'NO3')
    tas2_rcp85N = load_nitrate('rcp85', 'tas2')

    ENOX_allN = np.array(list(ENOX_histN)+list(ENOX_histN)+list(ENOX_histN)+list(ENOX_histN)+list(ENOX_rcp26N)+list(ENOX_rcp45N)+list(ENOX_rcp60N)+list(ENOX_rcp85N), dtype=dty)
    ENH3_allN = np.array(list(ENH3_histN)+list(ENH3_histN)+list(ENH3_histN)+list(ENH3_histN)+list(ENH3_rcp26N)+list(ENH3_rcp45N)+list(ENH3_rcp60N)+list(ENH3_rcp85N), dtype=dty)
    NO3_allN = np.array(list(NO3_histN)+list(NO3_histN)+list(NO3_histN)+list(NO3_histN)+list(NO3_rcp26N)+list(NO3_rcp45N)+list(NO3_rcp60N)+list(NO3_rcp85N), dtype=dty)
    tas2_allN = np.array(list(tas2_histN)+list(tas2_histN)+list(tas2_histN)+list(tas2_histN)+list(tas2_rcp26N)+list(tas2_rcp45N)+list(tas2_rcp60N)+list(tas2_rcp85N), dtype=dty)

# secondary organic aerosols
if mod_SOAload != "":
    # for sim in ["hist", "rcp26", "rcp45", "rcp60", "rcp85"]:
    #     for VAR in ["EVOC", "EBVOC", "SOA"] + ["tas2"]:
    #         print(f'{VAR}_{sim}Q = load_secondary_aerosols({sim!r}, {VAR!r})')

    def load_secondary_aerosols(sim, VAR):
        path = f"data/AeroChem_ACCMIP/#DATA.AeroChem_{mod_SOAload}.{prd[sim]}.{sim}_{VAR}.csv"
        TMP = load_data(path)
        return TMP[:,0].copy()

    EVOC_histQ = load_secondary_aerosols('hist', 'EVOC')
    EBVOC_histQ = load_secondary_aerosols('hist', 'EBVOC')
    SOA_histQ = load_secondary_aerosols('hist', 'SOA')
    tas2_histQ = load_secondary_aerosols('hist', 'tas2')
    EVOC_rcp26Q = load_secondary_aerosols('rcp26', 'EVOC')
    EBVOC_rcp26Q = load_secondary_aerosols('rcp26', 'EBVOC')
    SOA_rcp26Q = load_secondary_aerosols('rcp26', 'SOA')
    tas2_rcp26Q = load_secondary_aerosols('rcp26', 'tas2')
    EVOC_rcp45Q = load_secondary_aerosols('rcp45', 'EVOC')
    EBVOC_rcp45Q = load_secondary_aerosols('rcp45', 'EBVOC')
    SOA_rcp45Q = load_secondary_aerosols('rcp45', 'SOA')
    tas2_rcp45Q = load_secondary_aerosols('rcp45', 'tas2')
    EVOC_rcp60Q = load_secondary_aerosols('rcp60', 'EVOC')
    EBVOC_rcp60Q = load_secondary_aerosols('rcp60', 'EBVOC')
    SOA_rcp60Q = load_secondary_aerosols('rcp60', 'SOA')
    tas2_rcp60Q = load_secondary_aerosols('rcp60', 'tas2')
    EVOC_rcp85Q = load_secondary_aerosols('rcp85', 'EVOC')
    EBVOC_rcp85Q = load_secondary_aerosols('rcp85', 'EBVOC')
    SOA_rcp85Q = load_secondary_aerosols('rcp85', 'SOA')
    tas2_rcp85Q = load_secondary_aerosols('rcp85', 'tas2')

    EVOC_allQ = np.array(list(EVOC_histQ)+list(EVOC_histQ)+list(EVOC_histQ)+list(EVOC_histQ)+list(EVOC_rcp26Q)+list(EVOC_rcp45Q)+list(EVOC_rcp60Q)+list(EVOC_rcp85Q), dtype=dty)
    EBVOC_allQ = np.array(list(EBVOC_histQ)+list(EBVOC_histQ)+list(EBVOC_histQ)+list(EBVOC_histQ)+list(EBVOC_rcp26Q)+list(EBVOC_rcp45Q)+list(EBVOC_rcp60Q)+list(EBVOC_rcp85Q), dtype=dty)
    SOA_allQ = np.array(list(SOA_histQ)+list(SOA_histQ)+list(SOA_histQ)+list(SOA_histQ)+list(SOA_rcp26Q)+list(SOA_rcp45Q)+list(SOA_rcp60Q)+list(SOA_rcp85Q), dtype=dty)
    tas2_allQ = np.array(list(tas2_histQ)+list(tas2_histQ)+list(tas2_histQ)+list(tas2_histQ)+list(tas2_rcp26Q)+list(tas2_rcp45Q)+list(tas2_rcp60Q)+list(tas2_rcp85Q), dtype=dty)

# mineral dusts
# for sim in ["hist", "rcp26", "rcp45", "rcp60", "rcp85"]:
#     for VAR in ["EDUST", "DUST"] + ["tas"]:
#         print(f'{VAR}_{sim}D = load_mineral_dusts({sim!r}, {VAR!r})')


def load_mineral_dusts(sim, VAR):
    path = f"data/AeroChem_ACCMIP/#DATA.AeroChem_{mod_DUSTload}.{prd[sim]}.{sim}_{VAR}.csv"
    return load_data(path)[:, 0].copy()

EDUST_histD = load_mineral_dusts('hist', 'EDUST')
DUST_histD = load_mineral_dusts('hist', 'DUST')
tas_histD = load_mineral_dusts('hist', 'tas')
EDUST_rcp26D = load_mineral_dusts('rcp26', 'EDUST')
DUST_rcp26D = load_mineral_dusts('rcp26', 'DUST')
tas_rcp26D = load_mineral_dusts('rcp26', 'tas')
EDUST_rcp45D = load_mineral_dusts('rcp45', 'EDUST')
DUST_rcp45D = load_mineral_dusts('rcp45', 'DUST')
tas_rcp45D = load_mineral_dusts('rcp45', 'tas')
EDUST_rcp60D = load_mineral_dusts('rcp60', 'EDUST')
DUST_rcp60D = load_mineral_dusts('rcp60', 'DUST')
tas_rcp60D = load_mineral_dusts('rcp60', 'tas')
EDUST_rcp85D = load_mineral_dusts('rcp85', 'EDUST')
DUST_rcp85D = load_mineral_dusts('rcp85', 'DUST')
tas_rcp85D = load_mineral_dusts('rcp85', 'tas')

EDUST_allD = np.array(list(EDUST_histD)+list(EDUST_histD)+list(EDUST_histD)+list(EDUST_histD)+list(EDUST_rcp26D)+list(EDUST_rcp45D)+list(EDUST_rcp60D)+list(EDUST_rcp85D), dtype=dty)
DUST_allD = np.array(list(DUST_histD)+list(DUST_histD)+list(DUST_histD)+list(DUST_histD)+list(DUST_rcp26D)+list(DUST_rcp45D)+list(DUST_rcp60D)+list(DUST_rcp85D), dtype=dty)
tas_allD = np.array(list(tas_histD)+list(tas_histD)+list(tas_histD)+list(tas_histD)+list(tas_rcp26D)+list(tas_rcp45D)+list(tas_rcp60D)+list(tas_rcp85D), dtype=dty)

# sea salts
# for sim in ["hist", "rcp26", "rcp45", "rcp60", "rcp85"]:
#     for VAR in ["ESALT", "SALT"] + ["tas3"]:
#         print(f'{VAR}_{sim}T = load_sea_salts({sim!r}, {VAR!r})')


def load_sea_salts(sim, VAR):
    path = f"data/AeroChem_ACCMIP/#DATA.AeroChem_{mod_SALTload}.{prd[sim]}.{sim}_{VAR}.csv"
    return load_data(path)[:, 0].copy()

ESALT_histT = load_sea_salts('hist', 'ESALT')
SALT_histT = load_sea_salts('hist', 'SALT')
tas3_histT = load_sea_salts('hist', 'tas3')
ESALT_rcp26T = load_sea_salts('rcp26', 'ESALT')
SALT_rcp26T = load_sea_salts('rcp26', 'SALT')
tas3_rcp26T = load_sea_salts('rcp26', 'tas3')
ESALT_rcp45T = load_sea_salts('rcp45', 'ESALT')
SALT_rcp45T = load_sea_salts('rcp45', 'SALT')
tas3_rcp45T = load_sea_salts('rcp45', 'tas3')
ESALT_rcp60T = load_sea_salts('rcp60', 'ESALT')
SALT_rcp60T = load_sea_salts('rcp60', 'SALT')
tas3_rcp60T = load_sea_salts('rcp60', 'tas3')
ESALT_rcp85T = load_sea_salts('rcp85', 'ESALT')
SALT_rcp85T = load_sea_salts('rcp85', 'SALT')
tas3_rcp85T = load_sea_salts('rcp85', 'tas3')

ESALT_allT = np.array(list(ESALT_histT)+list(ESALT_histT)+list(ESALT_histT)+list(ESALT_histT)+list(ESALT_rcp26T)+list(ESALT_rcp45T)+list(ESALT_rcp60T)+list(ESALT_rcp85T), dtype=dty)
SALT_allT = np.array(list(SALT_histT)+list(SALT_histT)+list(SALT_histT)+list(SALT_histT)+list(SALT_rcp26T)+list(SALT_rcp45T)+list(SALT_rcp60T)+list(SALT_rcp85T), dtype=dty)
tas3_allT = np.array(list(tas3_histT)+list(tas3_histT)+list(tas3_histT)+list(tas3_histT)+list(tas3_rcp26T)+list(tas3_rcp45T)+list(tas3_rcp60T)+list(tas3_rcp85T), dtype=dty)

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
