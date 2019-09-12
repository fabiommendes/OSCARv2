from __future__ import absolute_import
from .default import *
from random import choice

# Drivers
data_EFF = choice(["CDIAC", "EDGAR"])
data_LULCC = choice(["LUH1"])
data_ECH4 = choice(["EDGAR", "ACCMIP", "EPA"])
data_EN2O = choice(["EDGAR", "EPA"])
data_Ehalo = choice(["EDGAR"])
data_ENOX = choice(["EDGAR", "ACCMIP"])
data_ECO = choice(["EDGAR", "ACCMIP"])
data_EVOC = choice(["EDGAR", "ACCMIP"])
data_ESO2 = choice(["EDGAR", "ACCMIP"])
data_ENH3 = choice(["EDGAR", "ACCMIP"])
data_EOC = choice(["ACCMIP"])
data_EBC = choice(["ACCMIP"])
data_RFant = choice(["IPCC-AR5"])
data_RFnat = choice(["IPCC-AR5"])

# Ocean
mod_OSNKstruct = choice(["HILDA", "BD-model", "2D-model", "3D-model"])
mod_OSNKchem = choice(["CO2SysPade", "CO2SysPower"])
mod_OSNKtrans = choice(["CESM1-BGC", "IPSL-CM5A-LR", "MPI-ESM-LR"])

# Biosphere
mod_LSNKnpp = choice(["log", "hyp"])
mod_LSNKrho = choice(["exp", "gauss"])
mod_LSNKpreind = choice(["CLM-45", "JSBACH", "JULES", "LPJ", "LPJ-GUESS", "LPX-Bern", "OCN", "ORCHIDEE", "VISIT"])
mod_LSNKtrans = choice(["BCC-CSM-11", "CESM1-BGC", "CanESM2", "HadGEM2-ES", "IPSL-CM5A-LR", "MPI-ESM-LR", "NorESM1-ME"])
mod_LSNKcover = choice(
    [
        "ESA-CCI",
        "MODIS",
        "Ramankutty1999",
        "Levavasseur2012",
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
)

# Wildfires
mod_EFIREpreind = choice(["", "CLM-45", "JSBACH", "LPJ", "LPJ-GUESS", "ORCHIDEE", "VISIT"])
mod_EFIREtrans = choice(["", "CESM1-BGC", "IPSL-CM5A-LR", "MPI-ESM-LR", "NorESM1-ME"])

# Permafrost
mod_EPFmain = choice(["", "JSBACH", "ORCHIDEE-MICT", "JULES-DeepResp", "JULES-SuppressResp"])
mod_EPFmethane = choice(["zero", "best", "twice"])

# Land-use
mod_ELUCagb = choice(["CLM-45", "LPJ-GUESS", "ORCHIDEE"])
mod_EHWPbb = choice(["high", "low"])
mod_EHWPtau = choice(["Houghton2001", "Earles2012"])
mod_EHWPspeed = choice(["normal", "fast", "slow"])

# Hydroxyl
mod_OHSNKtau = choice(
    [
        "Prather2012",
        "CESM-CAM-superfast",
        "CICERO-OsloCTM2",
        "CMAM",
        "EMAC",
        "GEOSCCM",
        "GDFL-AM3",
        "GISS-E2-R",
        "GISS-E2-R-TOMAS",
        "HadGEM2",
        "LMDzORINCA",
        "MIROC-CHEM",
        "MOCAGE",
        "NCAR-CAM-35",
        "STOC-HadAM3",
        "TM5",
        "UM-CAM",
    ]
)
mod_OHSNKfct = choice(["lin", "log"])
mod_OHSNKtrans = choice(["mean-OxComp", "Holmes2013", "GEOS-Chem", "Oslo-CTM3", "UCI-CTM"])

# Wetlands
mod_EWETpreind = choice(["", "CLM-4Me", "DLEM", "IAP-RAS", "LPJ-Bern", "LPJ-WSL", "ORCHIDEE", "SDGVM"])
mod_AWETtrans = choice(["", "CLM-4Me", "DLEM", "LPJ-Bern", "ORCHIDEE", "SDGVM", "UVic-ESCM"])

# Photolysis
mod_HVSNKtau = choice(["Prather2015", "GMI", "GEOSCCM", "G2d-M", "G2d", "Oslo-c29", "Oslo-c36", "UCI-c29", "UCI-c36"])
mod_HVSNKtrans = choice(["Prather2012", "Prather2015", "G2d", "Oslo-c29", "UCI-c29"])
mod_HVSNKcirc = choice(["CAM-35", "CMAM", "Niwa-SOCOL", "SOCOL", "ULAQ", "UMUKCA-UCAM"])

# Ozone tropo
mod_O3Tregsat = choice(
    [
        "",
        "CAMCHEM",
        "FRSGCUCI",
        "GISS-modelE",
        "GMI",
        "INCA",
        "LLNL-IMPACT",
        "MOZART-GFDL",
        "MOZECH",
        "STOC-HadAM3",
        "TM5-JRC",
        "UM-CAM",
    ]
)
mod_O3Temis = choice(["mean-OxComp", "CICERO-OsloCTM2", "NCAR-CAM-35", "STOC-HadAM3", "UM-CAM"])
mod_O3Tclim = choice(
    ["", "CESM-CAM-superfast", "GFDL-AM3", "GISS-E2-R", "MIROC-CHEM", "MOCAGE", "NCAR-CAM-35", "STOC-HadAM3", "UM-CAM"]
)
mod_O3Tradeff = choice(
    [
        "IPCC-AR5",
        "IPCC-AR4",
        "CESM-CAM-superfast",
        "CICERO-OsloCTM2",
        "CMAM",
        "EMAC",
        "GEOSCCM",
        "GFDL-AM3",
        "GISS-E2-R",
        "GISS-E2-R-TOMAS",
        "HadGEM2",
        "LMDzORINCA",
        "MIROC-CHEM",
        "MOCAGE",
        "NCAR-CAM-35",
        "STOC-HadAM3",
        "UM-CAM",
        "TM5",
    ]
)

# Ozone strato
mod_O3Sfracrel = choice(["Newman2006", "Laube2013"])
mod_O3Strans = choice(
    [
        "AMTRAC",
        "CCSR-NIES",
        "CMAM",
        "CNRM-ACM",
        "LMDZrepro",
        "MRI",
        "Niwa-SOCOL",
        "SOCOL",
        "ULAQ",
        "UMSLIMCAT",
        "UMUKCA-UCAM",
    ]
)
mod_O3Snitrous = choice(["", "Daniel2010"])
mod_O3Sradeff = choice(["IPCC-AR4", "ULAQ", "DLR-E39C", "NCAR-MACCM", "CHASER"])

# Sulfate
mod_SO4regsat = choice(["", "CAMCHEM", "GISS-PUCCINI", "GMI", "GOCART", "INCA2", "LLNL-IMPACT", "SPRINTARS"])
mod_SO4load = choice(["CSIRO-Mk360", "GFDL-AM3", "GISS-E2-R", "MIROC-CHEM"])
mod_SO4radeff = choice(
    [
        "BCC",
        "CAM4-Oslo",
        "CAM-51",
        "GEOS-CHEM",
        "GISS-MATRIX",
        "GISS-modelE",
        "GMI",
        "GOCART",
        "HadGEM2",
        "IMPACT-Umich",
        "INCA",
        "MPIHAM",
        "NCAR-CAM-35",
        "OsloCTM2",
        "SPRINTARS",
    ]
)

# POA
mod_POAconv = choice(["default", "GFDL", "CSIRO"])
mod_POAregsat = choice(["", "CAMCHEM", "GISS-PUCCINI", "GMI", "GOCART", "INCA2", "LLNL-IMPACT", "SPRINTARS"])
mod_POAload = choice(["CSIRO-Mk360", "GFDL-AM3", "GISS-E2-R", "MIROC-CHEM"])
mod_POAradeff = choice(
    [
        "BCC",
        "CAM4-Oslo",
        "CAM-51",
        "GEOS-CHEM",
        "GISS-MATRIX",
        "GISS-modelE",
        "GMI",
        "GOCART",
        "HadGEM2",
        "IMPACT-Umich",
        "INCA",
        "MPIHAM",
        "NCAR-CAM-35",
        "OsloCTM2",
        "SPRINTARS",
    ]
)

# BC
mod_BCregsat = choice(["", "CAMCHEM", "GISS-PUCCINI", "GMI", "GOCART", "INCA2", "LLNL-IMPACT", "SPRINTARS"])
mod_BCload = choice(["CSIRO-Mk360", "GFDL-AM3", "GISS-E2-R", "MIROC-CHEM"])
mod_BCradeff = choice(
    [
        "BCC",
        "CAM4-Oslo",
        "CAM-51",
        "GEOS-CHEM",
        "GISS-MATRIX",
        "GISS-modelE",
        "GMI",
        "GOCART",
        "HadGEM2",
        "IMPACT-Umich",
        "INCA",
        "MPIHAM",
        "NCAR-CAM-35",
        "OsloCTM2",
        "SPRINTARS",
    ]
)
mod_BCadjust = choice(["Boucher2013", "CSIRO", "GISS", "HadGEM2", "ECHAM5", "ECMWF"])

# Nitrate
mod_NO3load = choice(["Bellouin2011", "Hauglustaine2014"])
mod_NO3radeff = choice(
    ["GEOS-CHEM", "GISS-MATRIX", "GMI", "HadGEM2", "IMPACT-Umich", "INCA", "NCAR-CAM-35", "OsloCTM2"]
)

# SOA
mod_SOAload = choice(["", "GFDL-AM3", "GISS-E2-R"])
mod_SOAradeff = choice(["CAM-51", "GEOS-CHEM", "IMPACT-Umich", "MPIHAM", "OsloCTM2"])

# Dust
mod_DUSTload = choice(["CSIRO-Mk360", "GFDL-AM3", "GISS-E2-R", "MIROC-CHEM"])
mod_DUSTradeff = choice([""])  # [no value for now]

# Salt
mod_SALTload = choice(["GFDL-AM3", "GISS-E2-R", "MIROC-CHEM"])
mod_SALTradeff = choice([""])  # [no value for now]

# Cloud
mod_CLOUDsolub = choice(["Hansen2005", "Lamarque2011"])
mod_CLOUDerf = choice(
    ["mean-ACCMIP", "CSIRO-Mk360", "GFDL-AM3", "GISS-E2-R", "HadGEM2", "LMDzORINCA", "MIROC-CHEM", "NCAR-CAM-51"]
)
mod_CLOUDpreind = choice(["low", "median", "high"])

# Albedo bc
mod_ALBBCreg = choice(["Reddy2007"])
mod_ALBBCrf = choice(
    [
        "CICERO-OsloCTM2",
        "GFDL-AM3",
        "GISS-E2-R",
        "GISS-E2-R-TOMAS",
        "HadGEM2",
        "MIROC-CHEM",
        "NCAR-CAM-35",
        "NCAR-CAM-51",
    ]
)
mod_ALBBCwarm = choice(["low", "median", "high"])

# Albedo lc
mod_ALBLCflux = choice(["CERES", "GEWEX", "MERRA"])
mod_ALBLCalb = choice(["GlobAlbedo", "MODIS"])
mod_ALBLCcover = choice(["ESA-CCI", "MODIS"])
mod_ALBLCwarm = choice(["Hansen2005", "Davin2007", "Davin2010", "Jones2013"])

# Temperature
mod_TEMPresp = choice(
    [
        "ACCESS-10",
        "ACCESS-13",
        "BCC-CSM-11",
        "BCC-CSM-11m",
        "CanESM2",
        "CCSM4",
        "CNRM-CM5",
        "CNRM-CM5-2",
        "CSIRO-Mk360",
        "GFDL-CM3",
        "GFDL-ESM2G",
        "GFDL-ESM2M",
        "GISS-E2-H",
        "GISS-E2-R",
        "HadGEM2-ES",
        "IPSL-CM5A-LR",
        "IPSL-CM5A-MR",
        "IPSL-CM5B-LR",
        "MIROC5",
        "MIROC-ESM",
        "MPI-ESM-LR",
        "MPI-ESM-MR",
        "MPI-ESM-P",
        "MRI-CGCM3",
        "NorESM1-M",
    ]
)
mod_TEMPpattern = choice(["4xCO2", "hist&RCPs"])

# Precipitations
mod_PRECresp = choice(
    [
        "ACCESS-10",
        "ACCESS-13",
        "BCC-CSM-11",
        "BCC-CSM-11m",
        "CanESM2",
        "CCSM4",
        "CNRM-CM5",
        "CNRM-CM5-2",
        "CSIRO-Mk360",
        "GFDL-CM3",
        "GFDL-ESM2G",
        "GFDL-ESM2M",
        "GISS-E2-H",
        "GISS-E2-R",
        "HadGEM2-ES",
        "IPSL-CM5A-LR",
        "IPSL-CM5A-MR",
        "IPSL-CM5B-LR",
        "MIROC5",
        "MIROC-ESM",
        "MPI-ESM-LR",
        "MPI-ESM-MR",
        "MPI-ESM-P",
        "MRI-CGCM3",
        "NorESM1-M",
    ]
)
mod_PRECradfact = choice(["Andrews2010", "Kvalevag2013"])
mod_PRECpattern = choice(["4xCO2", "hist&RCPs"])

# Acidification
mod_ACIDsurf = choice(["Tans2009", "Bernie2010"])

# SLR
mod_SLR = choice([""])  # [no value for now]
