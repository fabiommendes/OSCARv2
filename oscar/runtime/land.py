import numpy as np

from oscar.runtime.drivers import LUC
from oscar.runtime.oscar_param import AWET_0, gamma_wetT, gamma_wetP, gamma_wetC, gamma_igniT, gamma_igniP, gamma_igniC, f_rho, npp_0, f_npp, igni_0, cveg_0, mu_0, rho1_0, csoil1_0, k_met, rho2_0, csoil2_0, ewet_0, p_wet, AREA_0
from oscar.runtime.util import mk_land_var, mk_region_biome_var, mk_region_var, mk_scalar


class LandSimulatorMixin:
    dt: float
    fT: np.ndarray
    D_CO2: np.ndarray
    D_lst: np.ndarray
    D_lyp: np.ndarray

    D_AREA = mk_land_var()
    D_cveg = mk_land_var()
    D_csoil1 = mk_land_var()
    D_csoil2 = mk_land_var()

    D_AREA_t = mk_region_biome_var()
    D_cveg_t = mk_region_biome_var()
    D_csoil1_t = mk_region_biome_var()
    D_csoil2_t = mk_region_biome_var()

    D_AWET_t = mk_region_var()
    D_npp_t = mk_region_biome_var()
    D_efire_t = mk_region_biome_var()
    D_fmort_t = mk_region_biome_var()
    D_rh1_t = mk_region_biome_var()
    D_fmet_t = mk_region_biome_var()
    D_rh2_t = mk_region_biome_var()
    D_EWET_t = mk_region_var()
    D_ewet_t = mk_region_var()
    LSNK_t = mk_scalar()

    def step_land(self, t):
        dt = self.dt

        # land-cover
        self.D_AREA += dt * (np.sum(LUC[t], 1) - np.sum(LUC[t], 2))
        self.D_AWET = AWET_0 * (gamma_wetT * self.fT * self.D_lst + gamma_wetP * self.fT * self.D_lyp + gamma_wetC * self.fT * self.D_CO2)

        # factors
        self.D_k_igni = (gamma_igniT * self.fT * self.D_lst[:, np.newaxis] + gamma_igniP * self.fT * self.D_lyp[:, np.newaxis] + gamma_igniC * self.fT * self.D_CO2)
        self.D_k_rho = f_rho(self.fT * self.D_lst[:, np.newaxis], self.fT * self.D_lyp[:, np.newaxis])

        # fluxes
        self.D_npp = npp_0 * f_npp(self.D_CO2, self.fT * self.D_lst[:, np.newaxis], self.fT * self.D_lyp[:, np.newaxis])
        self.D_efire = igni_0 * ((1 + self.D_k_igni) * (cveg_0 + self.D_cveg) - cveg_0)
        self.D_fmort = mu_0 * self.D_cveg
        self.D_rh1 = rho1_0 * ((1 + self.D_k_rho) * (csoil1_0 + self.D_csoil1) - csoil1_0)
        self.D_fmet = k_met * self.D_rh1
        self.D_rh2 = rho2_0 * ((1 + self.D_k_rho) * (csoil2_0 + self.D_csoil2) - csoil2_0)
        self.D_ewet = ewet_0 * np.nan_to_num(np.sum(p_wet * self.D_csoil1, 1) / np.sum(p_wet * csoil1_0, 1))
        self.LSNK = np.sum((self.D_rh1 + self.D_rh2 + self.D_efire - self.D_npp) * (AREA_0 + self.D_AREA))
        self.D_EWET = ewet_0 * self.D_AWET + self.D_ewet * AWET_0 + self.D_ewet * self.D_AWET

        # stocks
        self.D_cveg += dt * (self.D_npp - self.D_fmort - self.D_efire)
        self.D_csoil1 += dt * (self.D_fmort - self.D_fmet - self.D_rh1)
        self.D_csoil2 += dt * (self.D_fmet - self.D_rh2)