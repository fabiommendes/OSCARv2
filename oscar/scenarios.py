import numpy as np

from .config import dty
from .data import load_data

__all__ = ['CO2_rcp', 'CH4_rcp', 'N2O_rcp']


def _load_rcp_scenario(var):
    arr = np.ones([800 + 1, 6], dtype=dty) * np.nan
    for n, rcp in enumerate(["rcp26", "rcp45", "rcp60", "rcp85", "rcp45to26", "rcp60to45"]):
        aux = load_data(f"data/Scenario_ECP/#DATA.Scenario_ECP.2000-2500.{rcp}_{var}.csv")
        arr[300:, n] = aux[:, 0]
    return arr


#: RCP concentrations {ppm} [Meinshausen et al., 2011]
CO2_rcp = _load_rcp_scenario('CO2')

#: RCP concentrations {ppb} [Meinshausen et al., 2011]
CH4_rcp = _load_rcp_scenario('CH4')

#: RPC concentrations {ppb} [Meinshausen et al., 2011]
N2O_rcp = _load_rcp_scenario('N2O')
