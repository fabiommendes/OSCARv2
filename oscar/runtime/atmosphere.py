import numpy as np

from ..constants import HFC, PFC, ODS
from .drivers import EFF, ECH4, EN2O, EHFC, EPFC, EODS
from .oscar_param import alpha_CO2, alpha_CH4, alpha_N2O, alpha_HFC, alpha_PFC, alpha_ODS, tau_lag
from .util import mk_scalar_var, mk_linear_var, mk_scalar

nb_HFC = len(HFC)
nb_PFC = len(PFC)
nb_ODS = len(ODS)


class AtmosphereSimulatorMixin:
    dt: float
    D_OHSNK_HFC: np.ndarray
    D_OHSNK_PFC: np.ndarray
    D_OHSNK_ODS: np.ndarray
    D_HVSNK_HFC: np.ndarray
    D_HVSNK_PFC: np.ndarray
    D_HVSNK_ODS: np.ndarray
    D_XSNK_HFC: np.ndarray
    D_XSNK_PFC: np.ndarray
    D_XSNK_ODS: np.ndarray
    D_OHSNK_CH4: np.ndarray
    D_HVSNK_CH4: np.ndarray
    D_XSNK_CH4: np.ndarray
    D_HVSNK_N2O: np.ndarray
    ELUC: np.ndarray
    LSNK: np.ndarray
    OSNK: np.ndarray
    D_FOXI_CH4: np.ndarray
    EPF_CO2: np.ndarray
    D_EBB_CH4: np.ndarray
    D_EWET: np.ndarray
    EPF_CH4: np.ndarray
    D_EBB_N2O: np.ndarray

    D_CO2 = mk_scalar_var()
    D_CH4 = mk_scalar_var()
    D_CH4_lag = mk_scalar_var()
    D_N2O = mk_scalar_var()
    D_N2O_lag = mk_scalar_var()
    D_HFC = mk_linear_var(nb_HFC)
    D_HFC_lag = mk_linear_var(nb_HFC)
    D_PFC = mk_linear_var(nb_PFC)
    D_PFC_lag = mk_linear_var(nb_PFC)
    D_ODS = mk_linear_var(nb_ODS)
    D_ODS_lag = mk_linear_var(nb_ODS)

    D_CO2_t = mk_scalar()
    D_CH4_t = mk_scalar()
    D_CH4_lag_t = mk_scalar()
    D_N2O_t = mk_scalar()
    D_N2O_lag_t = mk_scalar()

    def step_atmosphere(self, t):
        dt = self.dt

        # stocks
        self.D_CO2 += dt * (1 / alpha_CO2) * (np.sum(EFF[t]) + np.sum(self.ELUC) + self.LSNK + self.OSNK + self.D_FOXI_CH4 + np.sum(self.EPF_CO2))
        self.D_CH4 += dt * (1 / alpha_CH4) * (np.sum(ECH4[t]) + np.sum(self.D_EBB_CH4) + np.sum(self.D_EWET) + np.sum(self.EPF_CH4) + self.D_OHSNK_CH4 + self.D_HVSNK_CH4 + self.D_XSNK_CH4)
        self.D_N2O += dt * (1 / alpha_N2O) * (np.sum(EN2O[t]) + np.sum(self.D_EBB_N2O) + self.D_HVSNK_N2O)
        self.D_HFC += dt * (1 / alpha_HFC) * (np.sum(EHFC[t], 0) + self.D_OHSNK_HFC + self.D_HVSNK_HFC + self.D_XSNK_HFC)
        self.D_PFC += dt * (1 / alpha_PFC) * (np.sum(EPFC[t], 0) + self.D_OHSNK_PFC + self.D_HVSNK_PFC + self.D_XSNK_PFC)
        self.D_ODS += dt * (1 / alpha_ODS) * (np.sum(EODS[t], 0) + self.D_OHSNK_ODS + self.D_HVSNK_ODS + self.D_XSNK_ODS)
        self.D_CH4_lag += dt * ((1 / tau_lag) * self.D_CH4 - (1 / tau_lag) * self.D_CH4_lag)
        self.D_N2O_lag += dt * ((1 / tau_lag) * self.D_N2O - (1 / tau_lag) * self.D_N2O_lag)
        self.D_HFC_lag += dt * ((1 / tau_lag) * self.D_HFC - (1 / tau_lag) * self.D_HFC_lag)
        self.D_PFC_lag += dt * ((1 / tau_lag) * self.D_PFC - (1 / tau_lag) * self.D_PFC_lag)
        self.D_ODS_lag += dt * ((1 / tau_lag) * self.D_ODS - (1 / tau_lag) * self.D_ODS_lag)