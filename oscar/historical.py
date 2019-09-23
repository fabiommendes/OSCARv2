import os
from functools import partial

import numpy as np

from .constants import HFC, PFC, ODS
from .data import load_data, load_data_and_header
from . import conf

__all__ = ['CO2_ipcc', 'CO2_lawdome', 'CO2_noaa_gl', 'CO2_noaa_ml', 'EFF_gcp', 'ELUC_gcp', 'LSNK_gcp', 'OSNK_gcp', 'd_CO2_gcp', 'CH4_0', 'CH4_agage', 'CH4_cmip5', 'CH4_ipcc', 'CH4_lawdome', 'N2O_0', 'N2O_agage', 'N2O_cmip5', 'N2O_ipcc',
    'N2O_lawdome', 'ODS_0', 'ODS_cmip5', 'ODS_ipcc', 'PFC', 'PFC_0', 'PFC_cmip5', 'PFC_ipcc', 'HFC', 'HFC_0', 'HFC_cmip5', 'HFC_ipcc', 'RF_AER_ipcc', 'RF_Alb_ipcc', 'RF_O3_ipcc', 'RF_WMGHG_ipcc', 'RF_cmip5', 'RF_ipcc', ]


def _load_historical(path, size, fill):
    data = np.ones([size], dtype=conf.dty) * np.nan
    aux = load_data(path)
    data[fill:] = aux[:, 0]
    return data


def _load_co2_flux(n, tmp):
    arr = np.ones([314 + 1], dtype=conf.dty) * np.nan
    arr[259:] = tmp[:, n]
    return arr


def _load_halogenated(kind, path, size, fill):
    arr = np.ones([size, len(kind)], dtype=conf.dty) * np.nan
    for i, var in enumerate(kind):
        path = path.format(VAR=var)
        if os.path.isfile(path):
            aux = load_data(path)
            arr[fill:, i] = aux[:, 0]
    return arr


def _empty(size, nan_start=None, nan_end=None, dtype=conf.dty):
    arr = np.zeros(size, dtype=dtype)
    if nan_start is not None:
        arr[:nan_start] = np.nan
    if nan_end is not None:
        arr[-nan_end:] = np.nan
    return arr


def _update_rf(path, total, start, mapping, volcanic='Volcano'):
    """
    Update radiative forces for all arrays in given mapping.

    Args:
         path: Path to config.data_loaders file.
         total: Array holding total radiative forcing.
         mapping: A map between column names to arrays holding the given RFs.
         volcanic: Volcanic RF are handled differently. Sets the column name for volcanic RFs.
    """
    aux, header = load_data_and_header(path)
    print(aux.shape)
    for i, col in enumerate(header):
        total[start:] += aux[1:, i]
        if col == volcanic:
            total[start:] -= np.mean(aux[1:, i])
        elif col in mapping:
            mapping[col][start:] += aux[1:, i]


_load_ipcc = partial(_load_historical, size=311 + 1, fill=50)
_load_cmip5 = partial(_load_historical, size=305 + 1, fill=65)

# ==============================================================================
# CARBON
# ==============================================================================

#: Historic CO2 from IPCC-AR5 {ppm} [IPCC WG1, 2013] annexe 2
CO2_ipcc = _load_ipcc("data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1750-2011.CO2.csv")
CO2_0 = CO2_ipcc[50]

#: Historic CO2 from CMIP5 {ppm} [Meinshausen et al., 2011]
CO2_cmip5 = _load_cmip5("data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005.CO2.csv")

#: Historic CO2 from NOAA {ppm}, from the website
CO2_noaa_ml = _load_historical("data/HistAtmo_NOAA-ESRL/#DATA.HistAtmo_NOAA-ESRL.1959-2014.CO2_maunaloa.csv", 314 + 1, fill=259)
CO2_noaa_gl = _load_historical("data/HistAtmo_NOAA-ESRL/#DATA.HistAtmo_NOAA-ESRL.1980-2014.CO2_global.csv", 314 + 1, fill=280)

#: Historic CO2 from Law Dome ice cores {ppm} [Etheridge et al., 1996] and [MacFarling Meure et al., 2006]
CO2_lawdome = load_data("data/HistAtmo_NOAA-NCDC/#DATA.HistAtmo_NOAA-NCDC.(IceCores).CO2_lawdome.csv", start=1)

#: Global CO2 historic flux {GtC/yr} [Le Quere et al., 2015]
_aux = load_data("data/Historic_GCP/#DATA.Historic_GCP.1959-2014_(5flx).budget.csv", start=1)
EFF_gcp = _load_co2_flux(0, _aux)
ELUC_gcp = _load_co2_flux(1, _aux)
d_CO2_gcp = _load_co2_flux(2, _aux)
OSNK_gcp = _load_co2_flux(3, _aux)
LSNK_gcp = _load_co2_flux(4, _aux)
del _aux, _load_co2_flux

# ==============================================================================
# METHANE
# ==============================================================================

#: Historic CH4 from IPCC-AR5 {ppb} [IPCC WG1, 2013] annexe 2
CH4_ipcc = _load_ipcc("data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1750-2011.CH4.csv")
CH4_0 = CH4_ipcc[50]

#: Historic CH4 from CMIP5 {ppb} [Meinshausen et al., 2011]
CH4_cmip5 = _load_cmip5("data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005.CH4.csv")

#: Historic CH4 from AGAGE {ppb} [Prinn et al., 2013] updated from the website
CH4_agage = _load_historical("data/HistAtmo_AGAGE/#DATA.HistAtmo_AGAGE.1987-2013.CH4_global.csv", size=313 + 1, fill=287)

#: Historic CH4 from Law Dome ice cores {ppm} [Etheridge et al., 1998] and [MacFarling Meure et al., 2006]
CH4_lawdome = load_data("data/HistAtmo_NOAA-NCDC/#DATA.HistAtmo_NOAA-NCDC.(IceCores).CH4_lawdome.csv", start=1)

# ==============================================================================
# NITROUS OXIDE
# ==============================================================================

#: Historic N2O from IPCC-AR5 {ppb} [IPCC WG1, 2013] annexe 2
N2O_ipcc = _load_ipcc("data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1750-2011.N2O.csv")
N2O_0 = N2O_ipcc[50]

#: Historic N2O from CMIP5 {ppb} [Meinshausen et al., 2011]
N2O_cmip5 = _load_cmip5("data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005.N2O.csv")

#: Historic N2O from AGAGE {ppb} [Prinn et al., 2013] updated from the website
N2O_agage = _load_historical("data/HistAtmo_AGAGE/#DATA.HistAtmo_AGAGE.1979-2013.N2O_global.csv", size=313 + 1, fill=279)

#: Historic N2O from Law Dome ice cores {ppm} [MacFarling Meure et al., 2006]
N2O_lawdome = load_data("data/HistAtmo_NOAA-NCDC/#DATA.HistAtmo_NOAA-NCDC.(IceCores).N2O_lawdome.csv", start=1)

# ==============================================================================
# HALOGENATED
# ==============================================================================

_path_hfc = "data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1940-2011.{VAR}.csv"
_path_pfc = "data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1900-2011.{VAR}.csv"
_path_ods = "data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1960-2011.{VAR}.csv"
_path_cmip5 = "data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005.{VAR}.csv"

#: Historic HaloC from IPCC-AR5 {ppt} [IPCC WG1, 2013] (annexe 2)
HFC_ipcc = _load_halogenated(HFC, _path_hfc, size=311 + 1, fill=240)
PFC_ipcc = _load_halogenated(PFC, _path_pfc, size=311 + 1, fill=200)
ODS_ipcc = _load_halogenated(ODS, _path_ods, size=311 + 1, fill=260)

#: Historic HaloC from CMIP5 {ppt} [Meinshausen et al., 2011]
HFC_cmip5 = _load_halogenated(HFC, _path_cmip5, size=305 + 1, fill=65)
PFC_cmip5 = _load_halogenated(PFC, _path_cmip5, size=305 + 1, fill=65)
ODS_cmip5 = _load_halogenated(ODS, _path_cmip5, size=305 + 1, fill=65)

#: Preindustrial HaloC concentrations {ppt} [IPCC WG1, 2013] (annexe 2) and [Meinshausen et al., 2011]
HFC_0 = np.zeros([len(HFC)], dtype=conf.dty)
PFC_0 = np.zeros([len(PFC)], dtype=conf.dty)
ODS_0 = np.zeros([len(ODS)], dtype=conf.dty)
PFC_0[PFC.index("CF4")] = 35.0
ODS_0[ODS.index("CH3Br")] = 5.8
ODS_0[ODS.index("CH3Cl")] = 480.0

# ==============================================================================
# RADIATIVE FORCINGS
# ==============================================================================

#: Historic RF from IPCC-AR5 {W/m2} [IPCC WG1, 2013] (annexe 2)
RF_ipcc = _empty(311 + 1, nan_start=50)
RF_WMGHG_ipcc = _empty(311 + 1, nan_start=50)
RF_O3_ipcc = _empty(311 + 1, nan_start=50)
RF_AER_ipcc = _empty(311 + 1, nan_start=50)
RF_Alb_ipcc = _empty(311 + 1, nan_start=50)

_path = "data/Historic_IPCC-AR5/#DATA.Historic_IPCC-AR5.1750-2011_(11for).RF.csv"
_update_rf(_path, RF_ipcc, 51, {"CO2": RF_WMGHG_ipcc, "GHG Other": RF_WMGHG_ipcc, "H2O (Strat)": RF_WMGHG_ipcc, "O3 (Trop)": RF_O3_ipcc, "O3 (Strat)": RF_O3_ipcc, "Aerosol (Total)": RF_AER_ipcc, "LUC": RF_Alb_ipcc, "BC Snow": RF_Alb_ipcc})

#: Historic RF from CMIP5 {W/m2} [Meinshausen et al., 2011]
RF_cmip5 = _empty(305 + 1, nan_start=65)

_path = "data/Historic_CMIP5/#DATA.Historic_CMIP5.1765-2005_(19for).RF.csv"
_update_rf(_path, RF_cmip5, 66, {}, volcanic='VOLC')
