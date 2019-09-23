import numpy as np

from .util import scalar_timeseries, obox_timeseries
from ..params import nb_obox, mld_0, alpha_mld, gamma_mld, p_circ, v_fg, alpha_CO2, f_pCO2, tau_circ, alpha_dic
from .. import conf


class OceanSimulationMixin:
    # Missing variables
    dt: float
    fT: int
    D_sst: np.ndarray
    D_CO2: np.ndarray

    # (obox) variables
    D_mld_t = scalar_timeseries()
    D_dic_t = scalar_timeseries()
    D_FIN_t = obox_timeseries()
    D_FOUT_t = obox_timeseries()
    D_FCIRC_t = obox_timeseries()
    D_CSURF_t = obox_timeseries()

    def __init__(self):
        super().__init__()
        self.D_dic = 0
        self.D_CSURF = np.zeros([nb_obox], dtype=conf.dty)

    def step_ocean(self, t):
        dt = self.dt

        # structure
        self.D_mld = mld_0 * alpha_mld * (np.exp(gamma_mld * self.fT * self.D_sst) - 1)

        # fluxes
        self.D_FIN = p_circ * v_fg * alpha_CO2 * self.D_CO2
        self.D_FOUT = p_circ * v_fg * alpha_CO2 * f_pCO2(self.D_dic, self.fT * self.D_sst)
        self.D_FCIRC = self.D_CSURF * (1 / tau_circ)
        self.OSNK = np.sum(self.D_FOUT - self.D_FIN)

        # stocks
        self.D_CSURF += dt * (self.D_FIN - self.D_FOUT - self.D_FCIRC)
        self.D_dic = alpha_dic * np.sum(self.D_CSURF) / (1 + self.D_mld / mld_0)