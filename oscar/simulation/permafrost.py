import numpy as np

from ..params import nb_regionPF, w_rhoPF, gamma_rhoPF1, w_reg_lstPF, gamma_rhoPF2, pthaw_min, k_pthaw, gamma_pthaw, f_v_PF, CFROZ_0, tau_PF1, tau_PF2, tau_PF3, p_PF_CH4, p_PF_inst, p_PF1, p_PF2, p_PF3
from .util import linear_var, regionPF_timeseries


class PermafrostSimulatorMixin:
    dt: float
    fT: np.ndarray
    D_gst: np.ndarray

    pthaw = linear_var(nb_regionPF)
    CTHAW1 = linear_var(nb_regionPF)
    CTHAW2 = linear_var(nb_regionPF)
    CTHAW3 = linear_var(nb_regionPF)
    D_CFROZ = linear_var(nb_regionPF)

    pthaw_t = regionPF_timeseries()
    CTHAW1_t = regionPF_timeseries()
    CTHAW2_t = regionPF_timeseries()
    CTHAW3_t = regionPF_timeseries()
    D_CFROZ_t = regionPF_timeseries()
    pthaw_bar_t = regionPF_timeseries()
    FTHAW_t = regionPF_timeseries()
    ETHAW1_t = regionPF_timeseries()
    ETHAW2_t = regionPF_timeseries()
    ETHAW3_t = regionPF_timeseries()
    EPF_CO2_t = regionPF_timeseries()
    EPF_CH4_t = regionPF_timeseries()
    EPF_t = regionPF_timeseries()

    def step_permafrost(self, t):
        dt = self.dt

        # factors
        self.rD_rhoPF = np.exp(w_rhoPF * gamma_rhoPF1 * w_reg_lstPF * self.fT * self.D_gst + w_rhoPF * gamma_rhoPF2 * (w_reg_lstPF * self.fT * self.D_gst) ** 2) - 1

        # fraction thawed
        self.pthaw_bar = -pthaw_min + (1 + pthaw_min) / (1 + ((1 / pthaw_min + 1) ** k_pthaw - 1) * np.exp(-gamma_pthaw * k_pthaw * w_reg_lstPF * self.fT * self.D_gst)) ** (1 / k_pthaw)
        self.d_pthaw = f_v_PF(self.pthaw_bar, self.pthaw) * (self.pthaw_bar - self.pthaw)
        self.pthaw += dt * self.d_pthaw

        # fluxes
        self.FTHAW = CFROZ_0 * self.d_pthaw
        self.ETHAW1 = 1 / tau_PF1 * (1 + self.rD_rhoPF) * self.CTHAW1
        self.ETHAW2 = 1 / tau_PF2 * (1 + self.rD_rhoPF) * self.CTHAW2
        self.ETHAW3 = 1 / tau_PF3 * (1 + self.rD_rhoPF) * self.CTHAW3
        self.EPF_CO2 = (1 - p_PF_CH4) * (self.ETHAW1 + self.ETHAW2 + self.ETHAW3 + p_PF_inst * self.FTHAW)
        self.EPF_CH4 = 1000.0 * p_PF_CH4 * (self.ETHAW1 + self.ETHAW2 + self.ETHAW3 + p_PF_inst * self.FTHAW)
        self.EPF = self.EPF_CO2 + 0.001 * self.EPF_CH4

        # stocks
        self.D_CFROZ -= dt * self.FTHAW
        self.CTHAW1 += dt * (p_PF1 * (1 - p_PF_inst) * self.FTHAW - self.ETHAW1)
        self.CTHAW2 += dt * (p_PF2 * (1 - p_PF_inst) * self.FTHAW - self.ETHAW2)
        self.CTHAW3 += dt * (p_PF3 * (1 - p_PF_inst) * self.FTHAW - self.ETHAW3)