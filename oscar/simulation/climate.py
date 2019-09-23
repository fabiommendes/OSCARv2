import numpy as np

from .drivers import EBC, RFcon, RFvolc, RFsolar
from ..data_loaders import nb_regionI
from ..params import f_RF_CO2, f_RF_CH4, f_RF_overlap, f_RF_H2Os, f_RF_N2O, radeff_HFC, radeff_PFC, radeff_ODS, radeff_O3t, radeff_O3s, radeff_SO4, radeff_POA, radeff_BC, radeff_NO3, radeff_SOA, radeff_DUST, radeff_SALT, \
    k_BC_adjust, Phi_0, AERh_0, radeff_BCsnow, w_reg_BCsnow, p_reg9, alpha_LCC, warmeff_BCsnow, warmeff_LCC, warmeff_volc, p_atm_CO2, p_atm_noCO2, p_atm_O3t, p_atm_strat, p_atm_scatter, p_atm_absorb, p_atm_cloud, p_atm_alb, p_atm_solar, \
    tau_gst, lambda_0, theta_0, tau_gst0, w_reg_sst, w_reg_lst, alpha_gyp, beta_gyp, w_reg_lyp, p_OHC, alpha_OHC, f_pH
from .util import scalar_var, linear_var, scalar_timeseries


class ClimateSimulatorMixin:
    dt: float
    D_O3t: np.ndarray
    D_O3s: np.ndarray
    D_SO4: np.ndarray
    D_POA: np.ndarray
    D_BC: np.ndarray
    D_NO3: np.ndarray
    D_SOA: np.ndarray
    D_DUST: np.ndarray
    D_SALT: np.ndarray
    D_CO2: np.ndarray
    D_CH4: np.ndarray
    D_CH4_lag: np.ndarray
    D_N2O: np.ndarray
    D_HFC: np.ndarray
    D_PFC: np.ndarray
    D_ODS: np.ndarray
    D_AREA: np.ndarray
    D_AERh: np.ndarray
    D_EBB_BC: np.ndarray

    D_gst = scalar_var()
    D_gst0 = scalar_var()
    D_sst = scalar_var()
    D_gyp = scalar_var()
    D_OHC = scalar_var()
    D_lst = linear_var(nb_regionI)
    D_lyp = linear_var(nb_regionI)

    RF_t = scalar_timeseries()
    RF_warm_t = scalar_timeseries()
    RF_atm_t = scalar_timeseries()
    RF_CO2_t = scalar_timeseries()
    RF_CH4_t = scalar_timeseries()
    RF_H2Os_t = scalar_timeseries()
    RF_N2O_t = scalar_timeseries()
    RF_halo_t = scalar_timeseries()
    RF_O3t_t = scalar_timeseries()
    RF_O3s_t = scalar_timeseries()
    RF_SO4_t = scalar_timeseries()
    RF_POA_t = scalar_timeseries()
    RF_BC_t = scalar_timeseries()
    RF_NO3_t = scalar_timeseries()
    RF_SOA_t = scalar_timeseries()
    RF_cloud_t = scalar_timeseries()
    RF_BCsnow_t = scalar_timeseries()
    RF_LCC_t = scalar_timeseries()
    D_gst_t = scalar_timeseries()
    D_sst_t = scalar_timeseries()
    D_gyp_t = scalar_timeseries()
    D_OHC_t = scalar_timeseries()

    def step_climate(self, t):
        dt = self.dt

        # fluxes
        # per component
        self.RF_CO2 = f_RF_CO2(self.D_CO2)
        self.RF_CH4 = f_RF_CH4(self.D_CH4) - (f_RF_overlap(self.D_CH4, self.D_N2O) - f_RF_overlap(0, self.D_N2O))
        self.RF_H2Os = f_RF_H2Os(self.D_CH4_lag)
        self.RF_N2O = f_RF_N2O(self.D_N2O) - (f_RF_overlap(self.D_CH4, self.D_N2O) - f_RF_overlap(self.D_CH4, 0))
        self.RF_halo = np.sum(radeff_HFC * self.D_HFC) + np.sum(radeff_PFC * self.D_PFC) + np.sum(radeff_ODS * self.D_ODS)
        self.RF_O3t = radeff_O3t * self.D_O3t
        self.RF_O3s = radeff_O3s * self.D_O3s
        self.RF_SO4 = radeff_SO4 * self.D_SO4
        self.RF_POA = radeff_POA * self.D_POA
        self.RF_BC = radeff_BC * self.D_BC
        self.RF_NO3 = radeff_NO3 * self.D_NO3
        self.RF_SOA = radeff_SOA * self.D_SOA
        self.RF_DUST = radeff_DUST * self.D_DUST
        self.RF_SALT = radeff_SALT * self.D_SALT
        self.RF_cloud = k_BC_adjust * self.RF_BC + Phi_0 * np.log(1 + self.D_AERh / AERh_0)
        self.RF_BCsnow = radeff_BCsnow * np.sum(w_reg_BCsnow * np.sum(p_reg9 * (EBC[t] + self.D_EBB_BC)[:, np.newaxis], 0))
        self.RF_LCC = np.sum(alpha_LCC * self.D_AREA)

        # totals
        self.RF = (
                self.RF_CO2
                + self.RF_CH4
                + self.RF_H2Os
                + self.RF_N2O
                + self.RF_halo
                + self.RF_O3t
                + self.RF_O3s
                + self.RF_SO4
                + self.RF_POA
                + self.RF_BC
                + self.RF_NO3
                + self.RF_SOA
                + self.RF_DUST
                + self.RF_SALT
                + self.RF_cloud
                + self.RF_BCsnow
                + self.RF_LCC
                + RFcon[t]
                + RFvolc[t]
                + RFsolar[t])
        self.RF_warm = (
                self.RF_CO2
                + self.RF_CH4
                + self.RF_H2Os
                + self.RF_N2O
                + self.RF_halo
                + self.RF_O3t
                + self.RF_O3s
                + self.RF_SO4
                + self.RF_POA
                + self.RF_BC
                + self.RF_NO3
                + self.RF_SOA
                + self.RF_DUST
                + self.RF_SALT
                + self.RF_cloud
                + warmeff_BCsnow * self.RF_BCsnow
                + warmeff_LCC * self.RF_LCC
                + warmeff_volc * RFvolc[t]
                + RFcon[t]
                + RFsolar[t])
        self.RF_atm = (
                p_atm_CO2 * self.RF_CO2
                + p_atm_noCO2 * (self.RF_CH4 + self.RF_N2O + self.RF_halo)
                + p_atm_O3t * self.RF_O3t
                + p_atm_strat * (self.RF_O3s + self.RF_H2Os)
                + p_atm_scatter * (self.RF_SO4 + self.RF_POA + self.RF_NO3 + self.RF_SOA + self.RF_DUST + self.RF_SALT + RFvolc[t])
                + p_atm_absorb * self.RF_BC
                + p_atm_cloud * (self.RF_cloud + RFcon[t])
                + p_atm_alb * (self.RF_BCsnow + self.RF_LCC)
                + p_atm_solar * RFsolar[t])

        # stocks
        # temperatures
        self.D_gst += dt * (1 / tau_gst) * (lambda_0 * self.RF_warm - self.D_gst - theta_0 * (self.D_gst - self.D_gst0))
        self.D_gst0 += dt * (1 / tau_gst0) * theta_0 * (self.D_gst - self.D_gst0)
        self.D_sst = w_reg_sst * self.D_gst
        self.D_lst = w_reg_lst * self.D_gst

        # precipitations
        self.D_gyp = alpha_gyp * self.D_gst + beta_gyp * self.RF_atm
        self.D_lyp = w_reg_lyp * self.D_gyp

        # ocean
        self.D_OHC += dt * p_OHC * alpha_OHC * (self.RF - self.D_gst / lambda_0)
        self.D_pH = f_pH(self.D_CO2)