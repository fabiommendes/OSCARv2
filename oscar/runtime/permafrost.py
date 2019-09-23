import numpy as np

from oscar.runtime.oscar_param import nb_regionPF, w_rhoPF, gamma_rhoPF1, w_reg_lstPF, gamma_rhoPF2, pthaw_min, k_pthaw, gamma_pthaw, f_v_PF, CFROZ_0, tau_PF1, tau_PF2, tau_PF3, p_PF_CH4, p_PF_inst, p_PF1, p_PF2, p_PF3
from oscar.runtime.util import mk_linear_var, mk_regionPF_var


class PermafrostSimulatorMixin:
    dt: float
    fT: np.ndarray
    D_gst: np.ndarray

    pthaw = mk_linear_var(nb_regionPF)
    CTHAW1 = mk_linear_var(nb_regionPF)
    CTHAW2 = mk_linear_var(nb_regionPF)
    CTHAW3 = mk_linear_var(nb_regionPF)
    D_CFROZ = mk_linear_var(nb_regionPF)

    pthaw_t = mk_regionPF_var()
    CTHAW1_t = mk_regionPF_var()
    CTHAW2_t = mk_regionPF_var()
    CTHAW3_t = mk_regionPF_var()
    D_CFROZ_t = mk_regionPF_var()
    pthaw_bar_t = mk_regionPF_var()
    FTHAW_t = mk_regionPF_var()
    ETHAW1_t = mk_regionPF_var()
    ETHAW2_t = mk_regionPF_var()
    ETHAW3_t = mk_regionPF_var()
    EPF_CO2_t = mk_regionPF_var()
    EPF_CH4_t = mk_regionPF_var()
    EPF_t = mk_regionPF_var()

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