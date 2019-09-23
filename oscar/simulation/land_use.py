import numpy as np

from .drivers import LUC, HARV, SHIFT
from .util import land_use_var, region_biome2_age_timeseries, region_timeseries
from ..data_loaders import nb_biome, nb_regionI
from ..params import cveg_0, csoil1_0, csoil2_0, p_AGB, p_HWP0, p_HWP1, p_HWP2, p_HWP3, mu_0, tau_shift, igni_0, rho1_0, k_met, rho2_0, r_HWP1, r_HWP2, r_HWP3, AREA_0, p_HWP1_bb, alpha_BB_CO2, alpha_BB_CH4, alpha_BB_N2O, \
    alpha_BB_NOX, alpha_BB_CO, alpha_BB_VOC, alpha_BB_SO2, alpha_BB_NH3, alpha_BB_OC, alpha_BB_BC
from .. import config


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

    CVEG_luc = land_use_var()
    CSOIL1_luc = land_use_var()
    CSOIL2_luc = land_use_var()
    CHWP1_luc = land_use_var()
    CHWP2_luc = land_use_var()
    CHWP3_luc = land_use_var()

    EFIRE_luc_t = region_biome2_age_timeseries()
    FMORT_luc_t = region_biome2_age_timeseries()
    RH1_luc_t = region_biome2_age_timeseries()
    FMET_luc_t = region_biome2_age_timeseries()
    RH2_luc_t = region_biome2_age_timeseries()
    EHWP1_luc_t = region_biome2_age_timeseries()
    EHWP2_luc_t = region_biome2_age_timeseries()
    EHWP3_luc_t = region_biome2_age_timeseries()
    CVEG_luc_t = region_biome2_age_timeseries()
    CSOIL1_luc_t = region_biome2_age_timeseries()
    CSOIL2_luc_t = region_biome2_age_timeseries()
    CHWP1_luc_t = region_biome2_age_timeseries()
    CHWP2_luc_t = region_biome2_age_timeseries()
    CHWP3_luc_t = region_biome2_age_timeseries()

    ELUC_t = region_timeseries()
    D_EBB_CO2_t = region_timeseries()
    D_EBB_CH4_t = region_timeseries()
    D_EBB_N2O_t = region_timeseries()
    D_EBB_NOX_t = region_timeseries()
    D_EBB_CO_t = region_timeseries()
    D_EBB_VOC_t = region_timeseries()
    D_EBB_SO2_t = region_timeseries()
    D_EBB_NH3_t = region_timeseries()
    D_EBB_OC_t = region_timeseries()
    D_EBB_BC_t = region_timeseries()

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
        self.EHWP1_luc = np.zeros([nb_regionI, nb_biome, nb_biome, config.ind_final + 1], dtype=config.dty)
        self.EHWP1_luc[:, :, :, : t + 1] = r_HWP1[np.newaxis, np.newaxis, np.newaxis, t::-1] * self.CHWP1_luc[:, :, :, : t + 1]
        self.EHWP2_luc = np.zeros([nb_regionI, nb_biome, nb_biome, config.ind_final + 1], dtype=config.dty)
        self.EHWP2_luc[:, :, :, : t + 1] = r_HWP2[np.newaxis, np.newaxis, np.newaxis, t::-1] * self.CHWP2_luc[:, :, :, : t + 1]
        self.EHWP3_luc = np.zeros([nb_regionI, nb_biome, nb_biome, config.ind_final + 1], dtype=config.dty)
        self.EHWP3_luc[:, :, :, : t + 1] = r_HWP3[np.newaxis, np.newaxis, np.newaxis, t::-1] * self.CHWP3_luc[:, :, :, : t + 1]
        self.ELUC = np.sum(np.sum(np.sum(self.RH1_luc + self.RH2_luc + self.EFIRE_luc + self.EHWP1_luc + self.EHWP2_luc + self.EHWP3_luc - self.NPP_luc, 3), 2), 1)

    def __biomass_burning(self, t):
        self.D_EBB_CO2 = _biomass_burning_diff(alpha_BB_CO2, self.D_AREA, self.D_efire, self.EFIRE_luc, self.EHWP1_luc)
        self.D_EBB_CH4 = _biomass_burning_diff(alpha_BB_CH4, self.D_AREA, self.D_efire, self.EFIRE_luc, self.EHWP1_luc)
        self.D_EBB_N2O = _biomass_burning_diff(alpha_BB_N2O, self.D_AREA, self.D_efire, self.EFIRE_luc, self.EHWP1_luc)
        self.D_EBB_NOX = _biomass_burning_diff(alpha_BB_NOX, self.D_AREA, self.D_efire, self.EFIRE_luc, self.EHWP1_luc)
        self.D_EBB_CO = _biomass_burning_diff(alpha_BB_CO, self.D_AREA, self.D_efire, self.EFIRE_luc, self.EHWP1_luc)
        self.D_EBB_VOC = _biomass_burning_diff(alpha_BB_VOC, self.D_AREA, self.D_efire, self.EFIRE_luc, self.EHWP1_luc)
        self.D_EBB_SO2 = _biomass_burning_diff(alpha_BB_SO2, self.D_AREA, self.D_efire, self.EFIRE_luc, self.EHWP1_luc)
        self.D_EBB_NH3 = _biomass_burning_diff(alpha_BB_NH3, self.D_AREA, self.D_efire, self.EFIRE_luc, self.EHWP1_luc)
        self.D_EBB_OC = _biomass_burning_diff(alpha_BB_OC, self.D_AREA, self.D_efire, self.EFIRE_luc, self.EHWP1_luc)
        self.D_EBB_BC = _biomass_burning_diff(alpha_BB_BC, self.D_AREA, self.D_efire, self.EFIRE_luc, self.EHWP1_luc)


def _biomass_burning_diff(alpha_BB, area, efire, efire_luc, ehwp1_luc):
    new = np.newaxis
    diff = np.sum(alpha_BB * (igni_0 * cveg_0 * area + efire * AREA_0 + efire * area), 1)
    diff += p_HWP1_bb * np.sum(np.sum(np.sum(alpha_BB[:, :, new, new] * ehwp1_luc, 3), 2), 1)
    diff += np.sum(np.sum(np.sum(alpha_BB[:, new, :, new] * efire_luc, 3), 2), 1)
    return diff
