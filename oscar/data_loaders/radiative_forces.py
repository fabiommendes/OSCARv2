"""
Radiative forces
"""
import numpy as np

from .regions import nb_regionJ, nb_kind, nb_sector, kindRF_index
from .greenhouse import ind_cdiac
from .. import config
from ..data import load_data_and_header

_kind_idx = {"Contrails": kindRF_index["RFcon"], "Solar": kindRF_index["RFsol"], "Volcano": kindRF_index["RFvol"]}


def _load_rf(which):
    data = np.zeros([config.ind_final + 1, nb_regionJ, nb_sector, nb_kind], dtype=config.dty)
    if config.data_RFant != "IPCC-AR5":
        return data

    path = "data/Historic_IPCC-AR5/#DATA.Historic_IPCC-AR5.1750-2011_(11for).RF.csv"
    aux, lgd = load_data_and_header(path)

    # Past estimates
    i = lgd.index(which)
    data[51: min(config.ind_final, ind_cdiac) + 1, 0, 0, _kind_idx[which]] += aux[1: 1 + min(config.ind_final, ind_cdiac) - 51 + 1, i]

    # Arbitrary extension
    if (config.scen_RFant == "cst") & (config.ind_final > ind_cdiac):
        data[ind_cdiac + 1:, 0, 0, _kind_idx[which]] = aux[-1, i]
    return data


#: RF contrails estimates from IPCC AR5 {W/m2} [IPCC, 2013] (annexe 2)
RFcon = _load_rf("Contrails")

#: RF solar radiation estimates from IPCC AR5 {W/m2} [IPCC, 2013] (annexe 2)
RFsolar = _load_rf("Solar")

#: RF volcanic forcings estimates from IPCC AR5 {W/m2} [IPCC, 2013] (annexe 2)
RFvolc = _load_rf("Volcano")
