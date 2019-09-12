##################################################
#   B. FINAL DRIVERS
##################################################
import numpy as np

from .a2_greenhouse import ECH4, EN2O, EFF, ECH4_0, EN2O_0
from .a3_land_use import LUC, HARV, SHIFT
from .a4_halogenated import EHFC, EPFC, EODS
from .a5_short_lived import ENOX, ECO, EVOC, ESO2, EOC, EBC, ENH3, \
    ENOX_0, ECO_0, EVOC_0, ESO2_0, EOC_0, EBC_0, ENH3_0
from .a6_radiative_forces import RFcon, RFvolc, RFsolar
from ...config import PI_1750, mod_sector, ind_final, ind_attrib

var_map = {
    "EFF": EFF,
    "ECH4": ECH4,
    "EN2O": EN2O,
    "LUC": LUC,
    "HARV": HARV,
    "SHIFT": SHIFT,
    "EHFC": EHFC,
    "EPFC": EPFC,
    "EODS": EODS,
    "ENOX": ENOX,
    "ECO": ECO,
    "EVOC": EVOC,
    "ESO2": ESO2,
    "ENH3": ENH3,
    "EOC": EOC,
    "EBC": EBC,
    "RFcon": RFcon,
    "RFvolc": RFvolc,
    "RFsolar": RFsolar,
    "ECH4_0": ECH4_0,
    "EN2O_0": EN2O_0,
    "ENOX_0": ENOX_0,
    "ECO_0": ECO_0,
    "EVOC_0": EVOC_0,
    "ESO2_0": ESO2_0,
    "ENH3_0": ENH3_0,
    "EOC_0": EOC_0,
    "EBC_0": EBC_0,
}

# ==================
# B.1. PREINDUSTRIAL
# ==================

# reference emissions for preindustrial
for VAR in ["ECH4", "EN2O"] + ["ENOX", "ECO", "EVOC", "ESO2", "ENH3", "EOC", "EBC"]:
    exec(VAR + "[:,...] -= " + VAR + "_0[np.newaxis,...]")

# force true 1750 preindustrial (drivers)
if PI_1750:
    for VAR in (
            ["EFF", "ECH4", "EN2O"]
            + ["LUC", "HARV", "SHIFT"]
            + ["EHFC", "EPFC", "EODS"]
            + ["ENOX", "ECO", "EVOC", "ESO2", "ENH3", "EOC", "EBC"]
            + ["RFcon", "RFvolc", "RFsolar"]
    ):
        exec(VAR + "[:50+1] *= 0")

# ================
# B.2. ATTRIBUTION
# ================


# starting year of regional attribution
for VAR in (
        ["EFF", "ECH4", "EN2O"]
        + ["LUC", "HARV", "SHIFT"]
        + ["EHFC", "EPFC", "EODS"]
        + ["ENOX", "ECO", "EVOC", "ESO2", "ENH3", "EOC", "EBC"]
        + ["RFcon", "RFvolc", "RFsolar"]
):
    ind_attrib
    exec(VAR + "[:ind_attrib,0,:,:,...] = np.sum(" + VAR + "[:ind_attrib,:,:,:,...],1)")
    exec(VAR + "[:ind_attrib,1:,:,:,...] = 0")

# timeframed sectoral attribution
if (mod_sector == "Time") & (300 < ind_final <= 310):
    for VAR in (
            ["EFF", "ECH4", "EN2O"]
            + ["LUC", "HARV", "SHIFT"]
            + ["EHFC", "EPFC", "EODS"]
            + ["ENOX", "ECO", "EVOC", "ESO2", "ENH3", "EOC", "EBC"]
            + ["RFcon", "RFvolc", "RFsolar"]
    ):
        exec(VAR + "[150:200,:,:,1,:,...] = " + VAR + "[150:200,:,:,0,:,...]")
        exec(VAR + "[200:210,:,:,2,:,...] = " + VAR + "[200:210,:,:,0,:,...]")
        exec(VAR + "[210:220,:,:,3,:,...] = " + VAR + "[210:220,:,:,0,:,...]")
        exec(VAR + "[220:230,:,:,4,:,...] = " + VAR + "[220:230,:,:,0,:,...]")
        exec(VAR + "[230:240,:,:,5,:,...] = " + VAR + "[230:240,:,:,0,:,...]")
        exec(VAR + "[240:250,:,:,6,:,...] = " + VAR + "[240:250,:,:,0,:,...]")
        exec(VAR + "[250:260,:,:,7,:,...] = " + VAR + "[250:260,:,:,0,:,...]")
        exec(VAR + "[260:265,:,:,8,:,...] = " + VAR + "[260:265,:,:,0,:,...]")
        exec(VAR + "[265:270,:,:,9,:,...] = " + VAR + "[265:270,:,:,0,:,...]")
        exec(VAR + "[270:275,:,:,10,:,...] = " + VAR + "[270:275,:,:,0,:,...]")
        exec(VAR + "[275:280,:,:,11,:,...] = " + VAR + "[275:280,:,:,0,:,...]")
        exec(VAR + "[280:285,:,:,12,:,...] = " + VAR + "[280:285,:,:,0,:,...]")
        exec(VAR + "[285:290,:,:,13,:,...] = " + VAR + "[285:290,:,:,0,:,...]")
        exec(VAR + "[290:292,:,:,14,:,...] = " + VAR + "[290:292,:,:,0,:,...]")
        exec(VAR + "[292:294,:,:,15,:,...] = " + VAR + "[292:294,:,:,0,:,...]")
        exec(VAR + "[294:296,:,:,16,:,...] = " + VAR + "[294:296,:,:,0,:,...]")
        exec(VAR + "[296:298,:,:,17,:,...] = " + VAR + "[296:298,:,:,0,:,...]")
        exec(VAR + "[298:300,:,:,18,:,...] = " + VAR + "[298:300,:,:,0,:,...]")
        for t in range(300, ind_final + 1):
            exec(VAR + "[t,:,:,19+t-300,:,...] = " + VAR + "[t,:,:,0,:,...]")
        exec(VAR + "[150:,:,:,0,:,...] = 0")
elif (mod_sector == "TimeRCP") & (ind_final == 400):
    for VAR in (
            ["EFF", "ECH4", "EN2O"]
            + ["LUC", "HARV", "SHIFT"]
            + ["EHFC", "EPFC", "EODS"]
            + ["ENOX", "ECO", "EVOC", "ESO2", "ENH3", "EOC", "EBC"]
            + ["RFcon", "RFvolc", "RFsolar"]
    ):
        exec(VAR + "[311:316,:,1,:,...] = " + VAR + "[311:316,:,0,:,...]")
        exec(VAR + "[316:321,:,2,:,...] = " + VAR + "[316:321,:,0,:,...]")
        exec(VAR + "[321:326,:,3,:,...] = " + VAR + "[321:326,:,0,:,...]")
        exec(VAR + "[326:331,:,4,:,...] = " + VAR + "[326:331,:,0,:,...]")
        exec(VAR + "[331:336,:,5,:,...] = " + VAR + "[331:336,:,0,:,...]")
        exec(VAR + "[336:341,:,6,:,...] = " + VAR + "[336:341,:,0,:,...]")
        exec(VAR + "[341:346,:,7,:,...] = " + VAR + "[341:346,:,0,:,...]")
        exec(VAR + "[346:351,:,8,:,...] = " + VAR + "[346:351,:,0,:,...]")
        exec(VAR + "[351:356,:,9,:,...] = " + VAR + "[351:356,:,0,:,...]")
        exec(VAR + "[356:361,:,10,:,...] = " + VAR + "[356:361,:,0,:,...]")
        exec(VAR + "[361:366,:,11,:,...] = " + VAR + "[361:366,:,0,:,...]")
        exec(VAR + "[366:371,:,12,:,...] = " + VAR + "[366:371,:,0,:,...]")
        exec(VAR + "[371:376,:,13,:,...] = " + VAR + "[371:376,:,0,:,...]")
        exec(VAR + "[376:381,:,14,:,...] = " + VAR + "[376:381,:,0,:,...]")
        exec(VAR + "[381:386,:,15,:,...] = " + VAR + "[381:386,:,0,:,...]")
        exec(VAR + "[386:391,:,16,:,...] = " + VAR + "[386:391,:,0,:,...]")
        exec(VAR + "[391:396,:,17,:,...] = " + VAR + "[391:396,:,0,:,...]")
        exec(VAR + "[396:401,:,18,:,...] = " + VAR + "[396:401,:,0,:,...]")
        exec(VAR + "[311:,:,0,:,...] = 0")
