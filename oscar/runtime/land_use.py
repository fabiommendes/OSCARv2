import numpy as np

from oscar.config import ind_final, dty
from oscar.runtime.drivers import LUC, HARV, SHIFT
from oscar.runtime.oscar_data import nb_biome, nb_regionI
from oscar.runtime.oscar_param import cveg_0, csoil1_0, csoil2_0, p_AGB, p_HWP0, p_HWP1, p_HWP2, p_HWP3, mu_0, tau_shift, igni_0, rho1_0, k_met, rho2_0, r_HWP1, r_HWP2, r_HWP3, AREA_0, p_HWP1_bb, alpha_BB_CO2, alpha_BB_CH4, alpha_BB_N2O, \
    alpha_BB_NOX, alpha_BB_CO, alpha_BB_VOC, alpha_BB_SO2, alpha_BB_NH3, alpha_BB_OC, alpha_BB_BC
from oscar.runtime.util import mk_land_use_var, mk_region_biome2_age_var, mk_region_var


class LandUseSimulatorMixin:
    dt: float

    # Land variables
    D_AREA: np.ndarray
    D_cveg: np.ndarray
    D_csoil1: np.ndarray
    D_csoil2: np.ndarray
    D_efire: np.ndarray
    D_k_igni: np.ndarray
    D_k_rho: np.ndarray

    CVEG_luc = mk_land_use_var()
    CSOIL1_luc = mk_land_use_var()
    CSOIL2_luc = mk_land_use_var()
    CHWP1_luc = mk_land_use_var()
    CHWP2_luc = mk_land_use_var()
    CHWP3_luc = mk_land_use_var()

    EFIRE_luc_t = mk_region_biome2_age_var()
    FMORT_luc_t = mk_region_biome2_age_var()
    RH1_luc_t = mk_region_biome2_age_var()
    FMET_luc_t = mk_region_biome2_age_var()
    RH2_luc_t = mk_region_biome2_age_var()
    EHWP1_luc_t = mk_region_biome2_age_var()
    EHWP2_luc_t = mk_region_biome2_age_var()
    EHWP3_luc_t = mk_region_biome2_age_var()
    CVEG_luc_t = mk_region_biome2_age_var()
    CSOIL1_luc_t = mk_region_biome2_age_var()
    CSOIL2_luc_t = mk_region_biome2_age_var()
    CHWP1_luc_t = mk_region_biome2_age_var()
    CHWP2_luc_t = mk_region_biome2_age_var()
    CHWP3_luc_t = mk_region_biome2_age_var()

    ELUC_t = mk_region_var()
    D_EBB_CO2_t = mk_region_var()
    D_EBB_CH4_t = mk_region_var()
    D_EBB_N2O_t = mk_region_var()
    D_EBB_NOX_t = mk_region_var()
    D_EBB_CO_t = mk_region_var()
    D_EBB_VOC_t = mk_region_var()
    D_EBB_SO2_t = mk_region_var()
    D_EBB_NH3_t = mk_region_var()
    D_EBB_OC_t = mk_region_var()
    D_EBB_BC_t = mk_region_var()

    def step_land_use(self, t):
        dt = self.dt
        self.__initialization(t)
        self.__harvest(t)
        self.__shifting_cultivation(t)
        self.__book_keeping_model(t)
        self.__biomass_burning(t)

        # stocks
        self.CVEG_luc += dt * (self.NPP_luc - self.FMORT_luc - self.EFIRE_luc)
        self.CSOIL1_luc += dt * (self.FMORT_luc - self.FMET_luc - self.RH1_luc)
        self.CSOIL2_luc += dt * (self.FMET_luc - self.RH2_luc)
        self.CHWP1_luc += dt * -self.EHWP1_luc
        self.CHWP2_luc += dt * -self.EHWP2_luc
        self.CHWP3_luc += dt * -self.EHWP3_luc

    def __initialization(self, t):
        dt = self.dt
        for b1 in range(nb_biome):
            for b2 in range(nb_biome):
                self.CVEG_luc[:, b1, b2, t] += (-dt) * (cveg_0 + self.D_cveg)[:, b2] * LUC[t, :, b1, b2]
                self.CSOIL1_luc[:, b1, b2, t] += (dt * ((csoil1_0 + self.D_csoil1)[:, b1] - (csoil1_0 + self.D_csoil1)[:, b2]) * LUC[t, :, b1, b2])
                self.CSOIL2_luc[:, b1, b2, t] += (dt * ((csoil2_0 + self.D_csoil2)[:, b1] - (csoil2_0 + self.D_csoil2)[:, b2]) * LUC[t, :, b1, b2])
                self.CSOIL1_luc[:, b1, b2, t] += (dt * (cveg_0 + self.D_cveg)[:, b1] * (p_AGB[:, b1] * p_HWP0[:, b1] + (1 - p_AGB[:, b1])) * LUC[t, :, b1, b2])
                self.CHWP1_luc[:, b1, b2, t] += (dt * (cveg_0 + self.D_cveg)[:, b1] * p_AGB[:, b1] * p_HWP1[:, b1] * LUC[t, :, b1, b2])
                self.CHWP2_luc[:, b1, b2, t] += (dt * (cveg_0 + self.D_cveg)[:, b1] * p_AGB[:, b1] * p_HWP2[:, b1] * LUC[t, :, b1, b2])
                self.CHWP3_luc[:, b1, b2, t] += (dt * (cveg_0 + self.D_cveg)[:, b1] * p_AGB[:, b1] * p_HWP3[:, b1] * LUC[t, :, b1, b2])

    def __harvest(self, t):
        dt = self.dt
        for b in range(nb_biome):
            self.CVEG_luc[:, b, b, t] += dt * -HARV[t, :, b]
            self.CSOIL1_luc[:, b, b, t] += dt * p_HWP0[:, b] * HARV[t, :, b]
            self.CHWP1_luc[:, b, b, t] += dt * p_HWP1[:, b] * HARV[t, :, b]
            self.CHWP2_luc[:, b, b, t] += dt * p_HWP2[:, b] * HARV[t, :, b]
            self.CHWP3_luc[:, b, b, t] += dt * p_HWP3[:, b] * HARV[t, :, b]

    def __shifting_cultivation(self, t):
        dt = self.dt
        for b1 in range(nb_biome):
            for b2 in range(b1, nb_biome):
                self.CVEG_luc[:, b1, b2, t] -= (dt * (cveg_0 + self.D_cveg)[:, b2] * (1 - np.exp(-mu_0[:, b2] * tau_shift)) * SHIFT[t, :, b1, b2])
                self.CVEG_luc[:, b2, b1, t] -= (dt * (cveg_0 + self.D_cveg)[:, b1] * (1 - np.exp(-mu_0[:, b1] * tau_shift)) * SHIFT[t, :, b1, b2])

                self.CSOIL1_luc[:, b1, b2, t] += (dt * (cveg_0 + self.D_cveg)[:, b1] * (1 - np.exp(-mu_0[:, b1] * tau_shift)) * (p_AGB[:, b1] * p_HWP0[:, b1] + (1 - p_AGB[:, b1])) * SHIFT[t, :, b1, b2])
                self.CSOIL1_luc[:, b2, b1, t] += (dt * (cveg_0 + self.D_cveg)[:, b2] * (1 - np.exp(-mu_0[:, b2] * tau_shift)) * (p_AGB[:, b2] * p_HWP0[:, b2] + (1 - p_AGB[:, b2])) * SHIFT[t, :, b1, b2])

                self.CHWP1_luc[:, b1, b2, t] += (dt * (cveg_0 + self.D_cveg)[:, b1] * (1 - np.exp(-mu_0[:, b1] * tau_shift)) * p_AGB[:, b1] * p_HWP1[:, b1] * SHIFT[t, :, b1, b2])
                self.CHWP1_luc[:, b2, b1, t] += (dt * (cveg_0 + self.D_cveg)[:, b2] * (1 - np.exp(-mu_0[:, b2] * tau_shift)) * p_AGB[:, b2] * p_HWP1[:, b2] * SHIFT[t, :, b1, b2])

                self.CHWP2_luc[:, b1, b2, t] += (dt * (cveg_0 + self.D_cveg)[:, b1] * (1 - np.exp(-mu_0[:, b1] * tau_shift)) * p_AGB[:, b1] * p_HWP2[:, b1] * SHIFT[t, :, b1, b2])
                self.CHWP2_luc[:, b2, b1, t] += (dt * (cveg_0 + self.D_cveg)[:, b2] * (1 - np.exp(-mu_0[:, b2] * tau_shift)) * p_AGB[:, b2] * p_HWP2[:, b2] * SHIFT[t, :, b1, b2])

                self.CHWP3_luc[:, b1, b2, t] += (dt * (cveg_0 + self.D_cveg)[:, b1] * (1 - np.exp(-mu_0[:, b1] * tau_shift)) * p_AGB[:, b1] * p_HWP3[:, b1] * SHIFT[t, :, b1, b2])
                self.CHWP3_luc[:, b2, b1, t] += (dt * (cveg_0 + self.D_cveg)[:, b2] * (1 - np.exp(-mu_0[:, b2] * tau_shift)) * p_AGB[:, b2] * p_HWP3[:, b2] * SHIFT[t, :, b1, b2])

    def __book_keeping_model(self, t):
        self.NPP_luc = 0 * self.CVEG_luc
        self.EFIRE_luc = (igni_0 * (1 + self.D_k_igni))[:, np.newaxis, :, np.newaxis] * self.CVEG_luc
        self.FMORT_luc = mu_0[:, np.newaxis, :, np.newaxis] * self.CVEG_luc
        self.RH1_luc = (rho1_0 * (1 + self.D_k_rho))[:, np.newaxis, :, np.newaxis] * self.CSOIL1_luc
        self.FMET_luc = k_met * self.RH1_luc
        self.RH2_luc = (rho2_0 * (1 + self.D_k_rho))[:, np.newaxis, :, np.newaxis] * self.CSOIL2_luc
        self.EHWP1_luc = np.zeros([nb_regionI, nb_biome, nb_biome, ind_final + 1], dtype=dty)
        self.EHWP1_luc[:, :, :, : t + 1] = r_HWP1[np.newaxis, np.newaxis, np.newaxis, t::-1] * self.CHWP1_luc[:, :, :, : t + 1]
        self.EHWP2_luc = np.zeros([nb_regionI, nb_biome, nb_biome, ind_final + 1], dtype=dty)
        self.EHWP2_luc[:, :, :, : t + 1] = r_HWP2[np.newaxis, np.newaxis, np.newaxis, t::-1] * self.CHWP2_luc[:, :, :, : t + 1]
        self.EHWP3_luc = np.zeros([nb_regionI, nb_biome, nb_biome, ind_final + 1], dtype=dty)
        self.EHWP3_luc[:, :, :, : t + 1] = r_HWP3[np.newaxis, np.newaxis, np.newaxis, t::-1] * self.CHWP3_luc[:, :, :, : t + 1]
        self.ELUC = np.sum(np.sum(np.sum(self.RH1_luc + self.RH2_luc + self.EFIRE_luc + self.EHWP1_luc + self.EHWP2_luc + self.EHWP3_luc - self.NPP_luc, 3), 2), 1)

    def __biomass_burning(self, t):
        self.D_EBB_CO2 = self._biomass_burning_diff(alpha_BB_CO2)
        self.D_EBB_CH4 = self._biomass_burning_diff(alpha_BB_CH4)
        self.D_EBB_N2O = self._biomass_burning_diff(alpha_BB_N2O)
        self.D_EBB_NOX = self._biomass_burning_diff(alpha_BB_NOX)
        self.D_EBB_CO = self._biomass_burning_diff(alpha_BB_CO)
        self.D_EBB_VOC = self._biomass_burning_diff(alpha_BB_VOC)
        self.D_EBB_SO2 = self._biomass_burning_diff(alpha_BB_SO2)
        self.D_EBB_NH3 = self._biomass_burning_diff(alpha_BB_NH3)
        self.D_EBB_OC = self._biomass_burning_diff(alpha_BB_OC)
        self.D_EBB_BC = self._biomass_burning_diff(alpha_BB_BC)

    def _biomass_burning_diff(self, alpha_BB):
        new = np.newaxis
        diff = np.sum(alpha_BB * (igni_0 * cveg_0 * self.D_AREA + self.D_efire * AREA_0 + self.D_efire * self.D_AREA), 1)
        diff += p_HWP1_bb * np.sum(np.sum(np.sum(alpha_BB[:, :, new, new] * self.EHWP1_luc, 3), 2), 1)
        diff += np.sum(np.sum(np.sum(alpha_BB[:, new, :, new] * self.EFIRE_luc, 3), 2), 1)
        return diff
