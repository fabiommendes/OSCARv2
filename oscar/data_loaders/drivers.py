import numpy as np

from .greenhouse import ECH4, EN2O, EFF, ECH4_0, EN2O_0
from .land_use import LUC, HARV, SHIFT
from .halogenated import EHFC, EPFC, EODS
from .short_lived import ENOX, ECO, EVOC, ESO2, EOC, EBC, ENH3, ENOX_0, ECO_0, EVOC_0, ESO2_0, EOC_0, EBC_0, ENH3_0
from .radiative_forces import RFcon, RFvolc, RFsolar
from .. import conf


var_map = {
    "EFF": EFF, "ECH4": ECH4, "EN2O": EN2O, "LUC": LUC, "HARV": HARV, "SHIFT": SHIFT, "EHFC": EHFC, "EPFC": EPFC, "EODS": EODS, "ENOX": ENOX, "ECO": ECO, "EVOC": EVOC, "ESO2": ESO2, "ENH3": ENH3, "EOC": EOC, "EBC": EBC, "RFcon": RFcon, "RFvolc": RFvolc, "RFsolar": RFsolar, "ECH4_0": ECH4_0, "EN2O_0": EN2O_0, "ENOX_0": ENOX_0, "ECO_0": ECO_0, "EVOC_0": EVOC_0, "ESO2_0": ESO2_0, "ENH3_0": ENH3_0, "EOC_0": EOC_0, "EBC_0": EBC_0,
}

# ==================
# B.1. PREINDUSTRIAL
# ==================

# reference emissions for preindustrial
for VAR, VAR_0 in [
    (ECH4, ECH4_0), (EN2O, EN2O_0), (ENOX, ENOX_0), (ECO, ECO_0), (EVOC, EVOC_0), (ESO2, ESO2_0), (ENH3, ENH3_0), (EOC, EOC_0), (EBC, EBC_0),
]:
    VAR[:, ...] -= VAR_0[np.newaxis, ...]

# force true 1750 preindustrial (drivers)
if conf.PI_1750:
    EFF[:50 + 1] *= 0
    EFF[:conf.ind_attrib, 0, :, :, ...] = np.sum(EFF[:conf.ind_attrib, :, :, :, ...], 1)
    EFF[:conf.ind_attrib, 1, :, :, ...] = 0

    ECH4[:50 + 1] *= 0
    ECH4[:conf.ind_attrib, 0, :, :, ...] = np.sum(ECH4[:conf.ind_attrib, :, :, :, ...], 1)
    ECH4[:conf.ind_attrib, 1, :, :, ...] = 0

    EN2O[:50 + 1] *= 0
    EN2O[:conf.ind_attrib, 0, :, :, ...] = np.sum(EN2O[:conf.ind_attrib, :, :, :, ...], 1)
    EN2O[:conf.ind_attrib, 1, :, :, ...] = 0

    LUC[:50 + 1] *= 0
    LUC[:conf.ind_attrib, 0, :, :, ...] = np.sum(LUC[:conf.ind_attrib, :, :, :, ...], 1)
    LUC[:conf.ind_attrib, 1, :, :, ...] = 0

    HARV[:50 + 1] *= 0
    HARV[:conf.ind_attrib, 0, :, :, ...] = np.sum(HARV[:conf.ind_attrib, :, :, :, ...], 1)
    HARV[:conf.ind_attrib, 1, :, :, ...] = 0

    SHIFT[:50 + 1] *= 0
    SHIFT[:conf.ind_attrib, 0, :, :, ...] = np.sum(SHIFT[:conf.ind_attrib, :, :, :, ...], 1)
    SHIFT[:conf.ind_attrib, 1, :, :, ...] = 0

    EHFC[:50 + 1] *= 0
    EHFC[:conf.ind_attrib, 0, :, :, ...] = np.sum(EHFC[:conf.ind_attrib, :, :, :, ...], 1)
    EHFC[:conf.ind_attrib, 1, :, :, ...] = 0

    EPFC[:50 + 1] *= 0
    EPFC[:conf.ind_attrib, 0, :, :, ...] = np.sum(EPFC[:conf.ind_attrib, :, :, :, ...], 1)
    EPFC[:conf.ind_attrib, 1, :, :, ...] = 0

    EODS[:50 + 1] *= 0
    EODS[:conf.ind_attrib, 0, :, :, ...] = np.sum(EODS[:conf.ind_attrib, :, :, :, ...], 1)
    EODS[:conf.ind_attrib, 1, :, :, ...] = 0

    ENOX[:50 + 1] *= 0
    ENOX[:conf.ind_attrib, 0, :, :, ...] = np.sum(ENOX[:conf.ind_attrib, :, :, :, ...], 1)
    ENOX[:conf.ind_attrib, 1, :, :, ...] = 0

    ECO[:50 + 1] *= 0
    ECO[:conf.ind_attrib, 0, :, :, ...] = np.sum(ECO[:conf.ind_attrib, :, :, :, ...], 1)
    ECO[:conf.ind_attrib, 1, :, :, ...] = 0

    EVOC[:50 + 1] *= 0
    EVOC[:conf.ind_attrib, 0, :, :, ...] = np.sum(EVOC[:conf.ind_attrib, :, :, :, ...], 1)
    EVOC[:conf.ind_attrib, 1, :, :, ...] = 0

    ESO2[:50 + 1] *= 0
    ESO2[:conf.ind_attrib, 0, :, :, ...] = np.sum(ESO2[:conf.ind_attrib, :, :, :, ...], 1)
    ESO2[:conf.ind_attrib, 1, :, :, ...] = 0

    ENH3[:50 + 1] *= 0
    ENH3[:conf.ind_attrib, 0, :, :, ...] = np.sum(ENH3[:conf.ind_attrib, :, :, :, ...], 1)
    ENH3[:conf.ind_attrib, 1, :, :, ...] = 0

    EOC[:50 + 1] *= 0
    EOC[:conf.ind_attrib, 0, :, :, ...] = np.sum(EOC[:conf.ind_attrib, :, :, :, ...], 1)
    EOC[:conf.ind_attrib, 1, :, :, ...] = 0

    EBC[:50 + 1] *= 0
    EBC[:conf.ind_attrib, 0, :, :, ...] = np.sum(EBC[:conf.ind_attrib, :, :, :, ...], 1)
    EBC[:conf.ind_attrib, 1, :, :, ...] = 0

    RFcon[:50 + 1] *= 0
    RFcon[:conf.ind_attrib, 0, :, :, ...] = np.sum(RFcon[:conf.ind_attrib, :, :, :, ...], 1)
    RFcon[:conf.ind_attrib, 1, :, :, ...] = 0

    RFvolc[:50 + 1] *= 0
    RFvolc[:conf.ind_attrib, 0, :, :, ...] = np.sum(RFvolc[:conf.ind_attrib, :, :, :, ...], 1)
    RFvolc[:conf.ind_attrib, 1, :, :, ...] = 0

    RFsolar[:50 + 1] *= 0
    RFsolar[:conf.ind_attrib, 0, :, :, ...] = np.sum(RFsolar[:conf.ind_attrib, :, :, :, ...], 1)
    RFsolar[:conf.ind_attrib, 1, :, :, ...] = 0

# timeframed sectoral attribution
if (conf.mod_sector == "Time") & (300 < conf.ind_final <= 310):
    for VAR in (EFF, ECH4, EN2O, LUC, HARV, SHIFT, EHFC, EPFC, EODS, ENOX, ECO, EVOC, ESO2, ENH3, EOC, EBC, RFcon, RFvolc, RFsolar, ):
        VAR[150:200, :, :, 1, :, ...] = VAR[150:200, :, :, 0, :, ...]
        VAR[200:210, :, :, 2, :, ...] = VAR[200:210, :, :, 0, :, ...]
        VAR[210:220, :, :, 3, :, ...] = VAR[210:220, :, :, 0, :, ...]
        VAR[220:230, :, :, 4, :, ...] = VAR[220:230, :, :, 0, :, ...]
        VAR[230:240, :, :, 5, :, ...] = VAR[230:240, :, :, 0, :, ...]
        VAR[240:250, :, :, 6, :, ...] = VAR[240:250, :, :, 0, :, ...]
        VAR[250:260, :, :, 7, :, ...] = VAR[250:260, :, :, 0, :, ...]
        VAR[260:265, :, :, 8, :, ...] = VAR[260:265, :, :, 0, :, ...]
        VAR[265:270, :, :, 9, :, ...] = VAR[265:270, :, :, 0, :, ...]
        VAR[270:275, :, :, 10, :, ...] = VAR[270:275, :, :, 0, :, ...]
        VAR[275:280, :, :, 11, :, ...] = VAR[275:280, :, :, 0, :, ...]
        VAR[280:285, :, :, 12, :, ...] = VAR[280:285, :, :, 0, :, ...]
        VAR[285:290, :, :, 13, :, ...] = VAR[285:290, :, :, 0, :, ...]
        VAR[290:292, :, :, 14, :, ...] = VAR[290:292, :, :, 0, :, ...]
        VAR[292:294, :, :, 15, :, ...] = VAR[292:294, :, :, 0, :, ...]
        VAR[294:296, :, :, 16, :, ...] = VAR[294:296, :, :, 0, :, ...]
        VAR[296:298, :, :, 17, :, ...] = VAR[296:298, :, :, 0, :, ...]
        VAR[298:300, :, :, 18, :, ...] = VAR[298:300, :, :, 0, :, ...]
        for t in range(300, conf.ind_final + 1):
            VAR[t, :, :, 19 + t - 300, :, ...] = VAR[t, :, :, 0, :, ...]
        VAR[150:, :, :, 0, :, ...] = 0

elif (conf.mod_sector == "TimeRCP") & (conf.ind_final == 400):
    for VAR in (EFF, ECH4, EN2O, LUC, HARV, SHIFT, EHFC, EPFC, EODS, ENOX, ECO, EVOC, ESO2, ENH3, EOC, EBC, RFcon, RFvolc, RFsolar, ):
        VAR[311:316, :, 1, :, ...] = VAR[311:316, :, 0, :, ...]
        VAR[316:321, :, 2, :, ...] = VAR[316:321, :, 0, :, ...]
        VAR[321:326, :, 3, :, ...] = VAR[321:326, :, 0, :, ...]
        VAR[326:331, :, 4, :, ...] = VAR[326:331, :, 0, :, ...]
        VAR[331:336, :, 5, :, ...] = VAR[331:336, :, 0, :, ...]
        VAR[336:341, :, 6, :, ...] = VAR[336:341, :, 0, :, ...]
        VAR[341:346, :, 7, :, ...] = VAR[341:346, :, 0, :, ...]
        VAR[346:351, :, 8, :, ...] = VAR[346:351, :, 0, :, ...]
        VAR[351:356, :, 9, :, ...] = VAR[351:356, :, 0, :, ...]
        VAR[356:361, :, 10, :, ...] = VAR[356:361, :, 0, :, ...]
        VAR[361:366, :, 11, :, ...] = VAR[361:366, :, 0, :, ...]
        VAR[366:371, :, 12, :, ...] = VAR[366:371, :, 0, :, ...]
        VAR[371:376, :, 13, :, ...] = VAR[371:376, :, 0, :, ...]
        VAR[376:381, :, 14, :, ...] = VAR[376:381, :, 0, :, ...]
        VAR[381:386, :, 15, :, ...] = VAR[381:386, :, 0, :, ...]
        VAR[386:391, :, 16, :, ...] = VAR[386:391, :, 0, :, ...]
        VAR[391:396, :, 17, :, ...] = VAR[391:396, :, 0, :, ...]
        VAR[396:401, :, 18, :, ...] = VAR[396:401, :, 0, :, ...]
        VAR[311:, :, 0, :, ...] = 0
