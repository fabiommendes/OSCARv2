import numpy as np

from .methane import scale_OH
from .. import config

##################################################
#   4. HALOGENATED COMPOUNDS
##################################################

# ===============
# 4.1. ATMOSPHERE
# ===============

# conversion of HaloC from {ppt} to {kt}
alpha_HFC = 0.1765 * np.array([70.0, 52.0, 120.0, 102.0, 84.0, 66.0, 170.0, 152.0, 134.0, 148.1, 252.1], dtype=config.dty)
alpha_PFC = 0.1765 * np.array([146.1, 71.0, 88.0, 138.0, 188.0, 200.0, 238.0, 288.0, 338.0, 388.1], dtype=config.dty)
alpha_ODS = 0.1765 * np.array([137.4, 120.9, 187.4, 170.9, 154.5, 153.8, 133.4, 86.5, 117.0, 100.5, 165.4, 209.8, 148.9, 259.8, 94.9, 50.5], dtype=config.dty)


# ==============
# 4.2. CHEMISTRY
# ==============

# atmospheric OH lifetime {yr}
# from [WMO, 2011] (table 1-3)
# rescaled to follow the arbitrary rescaling of tau_CH4_OH
infs = lambda n: [np.inf] * n
tau_HFC_OH = scale_OH * np.array([245, 5.5, 32, 14.3, 55, 1.6, 44.5, 253, 8.2, 9.3, 17.9], dtype=config.dty)
tau_PFC_OH = scale_OH * np.array(infs(10), dtype=config.dty)
tau_ODS_OH = scale_OH * np.array([*infs(6), 6.1, 12.8, 10.7, 19.3, *infs(4), 1.5, 1.9], dtype=config.dty)

# stratospheric lifetime {yr}
# idem
# rescaled to follow [Prather et al., 2015]
tau_HFC_hv = 1.06 * np.array([2347, 89, 246, 232, 327, 45.4, 310, 5676, 116, 125, 157], dtype=config.dty)
tau_PFC_hv = 1.06 * np.array([3200, 500, 50000, 10000, 2600, 3200, 2600, 4100, 3100, 3000], dtype=config.dty)
tau_ODS_hv = 1.06 * np.array([45, 100, 85, 190, 1020, 35, 39, 186, 64.9, 160, *infs(6)], dtype=config.dty)

# other sinks lifetime {yr}
# idem
tau_HFC_othr = np.array(infs(11), dtype=config.dty)
tau_PFC_othr = np.array(infs(10), dtype=config.dty)
tau_ODS_othr = np.array([*infs(5), 101, 89, *infs(3), 16, 2.9, 65, 20, 3, 1.4], dtype=config.dty)
