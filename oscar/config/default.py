"""
Copyright: IIASA (International Institute for Applied Systems Analysis), 2016-2018; CEA (Commissariat a L'Energie Atomique) & UVSQ (Universite de Versailles et Saint-Quentin), 2016
Contributor(s): Thomas Gasser (gasser@iiasa.ac.at)

This software is a computer program whose purpose is to simulate the behavior of the Earth system, with a specific but not exclusive focus on anthropogenic climate change.

This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software.  You can use, modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 

As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability. 

In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the software by the user in light of its specific status of free software, that may mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data to be ensured and,  more generally, to use and operate it in the same conditions as regards security. 

The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
"""

import numpy as np

# time step (p-th year; must be >4)
p = 6
# carbon feedback (0 or 1)
fC = 1
# climate feedback (0 or 1)
fT = 1
# precision (float32 or float64)
dty = np.float32
# if False simulates the 1700-1750 period
PI_1750 = True

# ending year of run (+1700)
ind_final = 310
# starting year of attribution
ind_attrib = 0

# [deprecated]
attrib_DRIVERS = "producers"
# [deprecated]
attrib_FEEDBACKS = "emitters"
# [deprecated]
attrib_ELUCdelta = "causal"
# [deprecated]
attrib_ELUCampli = "causal"

# SRES4 | SRES11 | RECCAP* | Raupach* | Houghton | IMACLIM | Kyoto | RCP5 | RCP10*
mod_regionI = "Houghton"
# SRES4 | SRES11 | RECCAP* | Raupach* | Houghton | IMACLIM | Kyoto | RCP5 | RCP10*
mod_regionJ = "RCP5"

# '' | Time | TimeRCP
mod_sector = ""
# one | CDIAC
mod_kindFF = "one"
# one | all
mod_kindLUC = "one"
# one | RCP
mod_kindGHG = "one"
# one | all
mod_kindCHI = "one"
# one | all
mod_kindAER = "one"
# one | two | all
mod_kindRF = "one"
# '' | PUP
mod_kindGE = ""
# SHR | w/GRA | w/FOR
mod_biomeSHR = "w/GRA"
# URB | w/DES
mod_biomeURB = "w/DES"
# True | False
mod_biomeV3 = True

# CDIAC | EDGAR
data_EFF = "CDIAC"
# LUH1
data_LULCC = "LUH1"
# EDGAR | ACCMIP | EPA
data_ECH4 = "EDGAR"
# EDGAR | EPA
data_EN2O = "EDGAR"
# EDGAR
data_Ehalo = "EDGAR"
# EDGAR | ACCMIP
data_ENOX = "EDGAR"
# EDGAR | ACCMIP
data_ECO = "EDGAR"
# EDGAR | ACCMIP
data_EVOC = "EDGAR"
# EDGAR | ACCMIP
data_ESO2 = "EDGAR"
# EDGAR | ACCMIP
data_ENH3 = "EDGAR"
# ACCMIP
data_EOC = "ACCMIP"
# ACCMIP
data_EBC = "ACCMIP"

# '' | IPCC-AR5
data_RFant = "IPCC-AR5"
# '' | IPCC-AR5
data_RFnat = "IPCC-AR5"
# raw | offset | smoothX (X in yr) | trends
mod_DATAscen = "trends"

# '' | stop | cst | RCP8.5 | RCP6.0 | RCP4.5 | RCP2.6
scen_ALL = ""
# stop | cst | SRES-A1B | SRES-A1FI | SRES-A1T | SRES-A2 | SRES-B1 | SRES-B2
# | RCP8.5 | RCP6.0 | RCP4.5 | RCP2.6
scen_EFF = "stop"
# stop | cst | RCP8.5 | RCP6.0 | RCP4.5 | RCP2.6
scen_LULCC = "stop"
# stop | cst | SRES-A1B | SRES-A1FI | SRES-A1T | SRES-A2 | SRES-B1 | SRES-B2
# | RCP8.5 | RCP6.0 | RCP4.5 | RCP2.6
scen_ECH4 = "stop"
# stop | cst | SRES-A1B | SRES-A1FI | SRES-A1T | SRES-A2 | SRES-B1 | SRES-B2
# | RCP8.5 | RCP6.0 | RCP4.5 | RCP2.6
scen_EN2O = "stop"
# stop | cst | RCP8.5 | RCP6.0 | RCP4.5 | RCP2.6
scen_Ehalo = "stop"
# stop | cst | SRES-A1B | SRES-A1FI | SRES-A1T | SRES-A2 | SRES-B1 | SRES-B2
# | RCP8.5 | RCP6.0 | RCP4.5 | RCP2.6
scen_ENOX = "stop"
# stop | cst | SRES-A1B | SRES-A1FI | SRES-A1T | SRES-A2 | SRES-B1 | SRES-B2
# | RCP8.5 | RCP6.0 | RCP4.5 | RCP2.6
scen_ECO = "stop"
# stop | cst | SRES-A1B | SRES-A1FI | SRES-A1T | SRES-A2 | SRES-B1 | SRES-B2
# | RCP8.5 | RCP6.0 | RCP4.5 | RCP2.6
scen_EVOC = "stop"
# stop | cst | SRES-A1B | SRES-A1FI | SRES-A1T | SRES-A2 | SRES-B1 | SRES-B2
# | RCP8.5 | RCP6.0 | RCP4.5 | RCP2.6
scen_ESO2 = "stop"
# stop | cst | RCP8.5 | RCP6.0 | RCP4.5 | RCP2.6
scen_ENH3 = "stop"
# stop | cst | RCP8.5 | RCP6.0 | RCP4.5 | RCP2.6
scen_EOC = "stop"
# stop | cst | RCP8.5 | RCP6.0 | RCP4.5 | RCP2.6
scen_EBC = "stop"
# stop | cst
scen_RFant = "stop"
# stop | cst
scen_RFnat = "stop"

# HILDA | BD-model | 2D-model | 3D-model
mod_OSNKstruct = "HILDA"
# CO2SysPade | CO2SysPower
mod_OSNKchem = "CO2SysPower"
# mean-CMIP5 | CESM1-BGC | IPSL-CM5A-LR | MPI-ESM-LR
mod_OSNKtrans = "mean-CMIP5"
# log | hyp
mod_LSNKnpp = "hyp"
# exp | gauss
mod_LSNKrho = "exp"
# mean-TRENDYv2 | CLM-45 | JSBACH | JULES | LPJ | LPJ-GUESS | LPX-Bern | OCN
# | ORCHIDEE | VISIT
mod_LSNKpreind = "mean-TRENDYv2"
# mean-CMIP5 | BCC-CSM-11 | CESM1-BGC | CanESM2 | HadGEM2-ES | IPSL-CM5A-LR
# | MPI-ESM-LR | NorESM1-ME
mod_LSNKtrans = "mean-CMIP5"
# ESA-CCI | MODIS | Ramankutty1999 | Levavasseur2012 | mean-TRENDYv2 | CLM-45
# | JSBACH | JULES | LPJ | LPJ-GUESS | LPX-Bern | OCN | ORCHIDEE | VISIT
mod_LSNKcover = "mean-TRENDYv2"

# '' | mean-TRENDYv2 | CLM-45 | JSBACH | LPJ | LPJ-GUESS | ORCHIDEE | VISIT
mod_EFIREpreind = "mean-TRENDYv2"
# '' | mean-CMIP5 | CESM1-BGC | IPSL-CM5A-LR | MPI-ESM-LR | NorESM1-ME
mod_EFIREtrans = "mean-CMIP5"
# '' | JSBACH | ORCHIDEE-MICT | JULES-DeepResp | JULES-SuppressResp
mod_EPFmain = ""
# zero | best | twice
mod_EPFmethane = "best"
# mean-TRENDYv2 | CLM-45 | LPJ-GUESS | ORCHIDEE
mod_ELUCagb = "mean-TRENDYv2"

# high | low
mod_EHWPbb = "high"
# Houghton2001 | Earles2012
mod_EHWPtau = "Earles2012"
# normal | fast | slow
mod_EHWPspeed = "normal"

# Prather2012 | CESM-CAM-superfast | CICERO-OsloCTM2 | CMAM | EMAC | GEOSCCM
# | GDFL-AM3 | GISS-E2-R | GISS-E2-R-TOMAS | HadGEM2 | LMDzORINCA | MIROC-CHEM
# | MOCAGE | NCAR-CAM-35 | STOC-HadAM3 | TM5 | UM-CAM
mod_OHSNKtau = "Prather2012"
# lin | log
mod_OHSNKfct = "lin"
# mean-OxComp | Holmes2013 | GEOS-Chem | Oslo-CTM3 | UCI-CTM
mod_OHSNKtrans = "Holmes2013"
# '' | mean-WETCHIMP | CLM-4Me | DLEM | IAP-RAS | LPJ-Bern | LPJ-WSL | ORCHIDEE
# | SDGVM
mod_EWETpreind = "mean-WETCHIMP"
# '' | mean-WETCHIMP | CLM-4Me | DLEM | LPJ-Bern | ORCHIDEE | SDGVM | UVic-ESCM
mod_AWETtrans = "mean-WETCHIMP"

# Prather2015 | GMI | GEOSCCM | G2d-M | G2d | Oslo-c29 | Oslo-c36 | UCI-c29 | UCI-c36
mod_HVSNKtau = "Prather2015"
# Prather2012 | Prather2015 | G2d | Oslo-c29 | UCI-c29
mod_HVSNKtrans = "Prather2015"
# mean-CCMVal2 | AMTRAC | CAM-35 | CMAM | Niwa-SOCOL | SOCOL | ULAQ | UMUKCA-UCAM
mod_HVSNKcirc = "mean-CCMVal2"

# '' | mean-HTAP | CAMCHEM | FRSGCUCI | GISS-modelE | GMI | INCA | LLNL-IMPACT
# | MOZART-GFDL | MOZECH | STOC-HadAM3 | TM5-JRC | UM-CAM
mod_O3Tregsat = "mean-HTAP"
# mean-OxComp | mean-ACCMIP | CICERO-OsloCTM2 | NCAR-CAM-35 | STOC-HadAM3 | UM-CAM
mod_O3Temis = "mean-ACCMIP"
# '' | mean-ACCMIP | CESM-CAM-superfast | GFDL-AM3 | GISS-E2-R | MIROC-CHEM
# | MOCAGE | NCAR-CAM-35 | STOC-HadAM3 | UM-CAM
mod_O3Tclim = "mean-ACCMIP"
# IPCC-AR5 | IPCC-AR4 | mean-ACCMIP | CESM-CAM-superfast | CICERO-OsloCTM2 | CMAM
# | EMAC | GEOSCCM | GFDL-AM3 | GISS-E2-R | GISS-E2-R-TOMAS | HadGEM2 | LMDzORINCA
# | MIROC-CHEM | MOCAGE | NCAR-CAM-35 | STOC-HadAM3 | UM-CAM | TM5
mod_O3Tradeff = "IPCC-AR5"
# Newman2006 | Laube2013
mod_O3Sfracrel = "Newman2006"
# mean-CCMVal2 | AMTRAC | CCSR-NIES | CMAM | CNRM-ACM | LMDZrepro | MRI
# | Niwa-SOCOL | SOCOL | ULAQ | UMSLIMCAT | UMUKCA-UCAM
mod_O3Strans = "mean-CCMVal2"
# '' | Daniel2010
mod_O3Snitrous = "Daniel2010"
# IPCC-AR4 | mean-ACCENT | ULAQ | DLR-E39C | NCAR-MACCM | CHASER
mod_O3Sradeff = "IPCC-AR4"

# '' | mean-HTAP | CAMCHEM | GISS-PUCCINI | GMI | GOCART | INCA2 | LLNL-IMPACT
# | SPRINTARS
mod_SO4regsat = "mean-HTAP"
# mean-ACCMIP | CSIRO-Mk360 | GFDL-AM3 | GISS-E2-R | MIROC-CHEM
mod_SO4load = "mean-ACCMIP"
# mean-AeroCom2 | BCC | CAM4-Oslo | CAM-51 | GEOS-CHEM | GISS-MATRIX | GISS-modelE
# | GMI | GOCART | HadGEM2 | IMPACT-Umich | INCA | MPIHAM | NCAR-CAM-35
# | OsloCTM2 | SPRINTARS
mod_SO4radeff = "mean-AeroCom2"

# default | GFDL | CSIRO
mod_POAconv = "default"
# '' | mean-HTAP | CAMCHEM | GISS-PUCCINI | GMI | GOCART | INCA2 | LLNL-IMPACT
# | SPRINTARS
mod_POAregsat = "mean-HTAP"
# mean-ACCMIP | CSIRO-Mk360 | GFDL-AM3 | GISS-E2-R | MIROC-CHEM
mod_POAload = "mean-ACCMIP"
# mean-AeroCom2 | BCC | CAM4-Oslo | CAM-51 | GEOS-CHEM | GISS-MATRIX
# | GISS-modelE | GMI | GOCART | HadGEM2 | IMPACT-Umich | INCA | MPIHAM
# | NCAR-CAM-35 | OsloCTM2 | SPRINTARS
mod_POAradeff = "mean-AeroCom2"

# '' | mean-HTAP | CAMCHEM | GISS-PUCCINI | GMI | GOCART | INCA2 | LLNL-IMPACT
# | SPRINTARS
mod_BCregsat = "mean-HTAP"
# mean-ACCMIP | CSIRO-Mk360 | GFDL-AM3 | GISS-E2-R | MIROC-CHEM
mod_BCload = "mean-ACCMIP"
# mean-AeroCom2 | BCC | CAM4-Oslo | CAM-51 | GEOS-CHEM | GISS-MATRIX
# | GISS-modelE | GMI | GOCART | HadGEM2 | IMPACT-Umich | INCA | MPIHAM
# | NCAR-CAM-35 | OsloCTM2 | SPRINTARS
mod_BCradeff = "mean-AeroCom2"
# Boucher2013 | CSIRO | GISS | HadGEM2 | ECHAM5 | ECMWF
mod_BCadjust = "Boucher2013"

# Bellouin2011 | Hauglustaine2014
mod_NO3load = "Bellouin2011"
# mean-AeroCom2 | GEOS-CHEM | GISS-MATRIX | GMI | HadGEM2 | IMPACT-Umich
# | INCA | NCAR-CAM-35 | OsloCTM2
mod_NO3radeff = "mean-AeroCom2"

# '' | mean-ACCMIP | GFDL-AM3 | GISS-E2-R
mod_SOAload = "mean-ACCMIP"
# mean-AeroCom2 | CAM-51 | GEOS-CHEM | IMPACT-Umich | MPIHAM | OsloCTM2
mod_SOAradeff = "mean-AeroCom2"

# mean-ACCMIP | CSIRO-Mk360 | GFDL-AM3 | GISS-E2-R | MIROC-CHEM
mod_DUSTload = "mean-ACCMIP"
# [no value for now]
mod_DUSTradeff = ""

# mean-ACCMIP | GFDL-AM3 | GISS-E2-R | MIROC-CHEM
mod_SALTload = "mean-ACCMIP"
# [no value for now]
mod_SALTradeff = ""

# Hansen2005 | Lamarque2011
mod_CLOUDsolub = "Lamarque2011"
# mean-ACCMIP | CSIRO-Mk360 | GFDL-AM3 | GISS-E2-R | HadGEM2 | LMDzORINCA
# | MIROC-CHEM | NCAR-CAM-51
mod_CLOUDerf = "mean-ACCMIP"
# low | median | high
mod_CLOUDpreind = "median"

# Reddy2007
mod_ALBBCreg = "Reddy2007"
# mean-ACCMIP | CICERO-OsloCTM2 | GFDL-AM3 | GISS-E2-R | GISS-E2-R-TOMAS
# | HadGEM2 | MIROC-CHEM | NCAR-CAM-35 | NCAR-CAM-51
mod_ALBBCrf = "mean-ACCMIP"
# low | median | high
mod_ALBBCwarm = "median"
# CERES | GEWEX | MERRA
mod_ALBLCflux = "CERES"
# GlobAlbedo | MODIS
mod_ALBLCalb = "GlobAlbedo"
# ESA-CCI | MODIS
mod_ALBLCcover = "ESA-CCI"
# Hansen2005 | Davin2007 | Davin2010 | Jones2013
mod_ALBLCwarm = "Jones2013"

# mean-CMIP5 | ACCESS-10 | ACCESS-13 | BCC-CSM-11 | BCC-CSM-11m | CanESM2
# | CCSM4 | CNRM-CM5 | CNRM-CM5-2 | CSIRO-Mk360 | GFDL-CM3 | GFDL-ESM2G
# | GFDL-ESM2M | GISS-E2-H | GISS-E2-R | HadGEM2-ES | IPSL-CM5A-LR | IPSL-CM5A-MR
# | IPSL-CM5B-LR | MIROC5 | MIROC-ESM | MPI-ESM-LR | MPI-ESM-MR | MPI-ESM-P
# | MRI-CGCM3 | NorESM1-M
mod_TEMPresp = "mean-CMIP5"
# 4xCO2 | hist&RCPs
mod_TEMPpattern = "hist&RCPs"

# mean-CMIP5 | ACCESS-10 | ACCESS-13 | BCC-CSM-11 | BCC-CSM-11m | CanESM2
# | CCSM4 | CNRM-CM5 | CNRM-CM5-2 | CSIRO-Mk360 | GFDL-CM3 | GFDL-ESM2G
# | GFDL-ESM2M | GISS-E2-H | GISS-E2-R | HadGEM2-ES | IPSL-CM5A-LR | IPSL-CM5A-MR
# | IPSL-CM5B-LR | MIROC5 | MIROC-ESM | MPI-ESM-LR | MPI-ESM-MR | MPI-ESM-P
# | MRI-CGCM3 | NorESM1-M
mod_PRECresp = "mean-CMIP5"
# Andrews2010 | Kvalevag2013
mod_PRECradfact = "Andrews2010"
# 4xCO2 | hist&RCPs
mod_PRECpattern = "hist&RCPs"
# Tans2009 | Bernie2010
mod_ACIDsurf = "Bernie2010"
# [no value for now]
mod_SLR = ""

##################################################
#   2. OSCAR
##################################################

if scen_ALL != "":
    scen_LULCC = (
        scen_EFF
    ) = (
        scen_ECH4
    ) = (
        scen_EN2O
    ) = (
        scen_Ehalo
    ) = (
        scen_ENOX
    ) = scen_ECO = scen_EVOC = scen_ESO2 = scen_ENH3 = scen_EOC = scen_EBC = scen_RFant = scen_RFnat = scen_ALL
    if scen_RFant[:3] == "RCP":
        scen_RFant = "cst"
    if scen_RFnat[:3] == "RCP":
        scen_RFnat = "cst"
