import numpy as np

from ..constants import HFC, PFC, ODS
from ..historical import HFC_0, PFC_0, ODS_0
from .drivers import ENOX, ECO, EVOC, ECH4, ESO2, EOC, EBC, ENH3
from ..params import f_kOH, f_hv, alpha_CH4, tau_CH4_OH, CH4_0, tau_CH4_hv, tau_CH4_soil, tau_CH4_ocean, alpha_N2O, tau_N2O_hv, N2O_0, alpha_HFC, tau_HFC_OH, alpha_PFC, tau_PFC_OH, alpha_ODS, tau_ODS_OH, tau_HFC_hv, \
    tau_PFC_hv, tau_ODS_hv, tau_HFC_othr, tau_PFC_othr, tau_ODS_othr, chi_O3t_CH4, Gamma_O3t, chi_O3t_NOX, w_reg_NOX, p_reg4, chi_O3t_CO, w_reg_CO, chi_O3t_VOC, w_reg_VOC, f_fracrel, tau_lag, n_Cl, alpha_Br, n_Br, chi_O3s_EESC, \
    chi_O3s_N2O, \
    EESC_x, Gamma_O3s, alpha_SO4, tau_SO2, w_reg_SO2, tau_DMS, Gamma_SO4, tau_OMff, alpha_POM, w_reg_OC, tau_OMbb, Gamma_POA, tau_BCff, w_reg_BC, tau_BCbb, Gamma_BC, alpha_NO3, tau_NOX, tau_NH3, Gamma_NO3, tau_VOC, tau_BVOC, Gamma_SOA, \
    tau_DUST, Gamma_DUST, tau_SALT, Gamma_SALT, solub_SO4, solub_POA, solub_BC, solub_NO3, solub_SOA, solub_DUST, solub_SALT
from .util import scalar_var, species_timeseries, scalar_timeseries

nb_HFC = len(HFC)
nb_PFC = len(PFC)
nb_ODS = len(ODS)


# noinspection PyAttributeOutsideInit
class ChemestrySimulatorMixin:
    fT: np.ndarray
    D_EWET: np.ndarray
    D_CH4: np.ndarray
    D_CH4_lag: np.ndarray
    D_gst: np.ndarray
    D_N2O_lag: np.ndarray
    D_HFC: np.ndarray
    D_HFC_lag: np.ndarray
    D_PFC: np.ndarray
    D_PFC_lag: np.ndarray
    D_ODS: np.ndarray
    D_ODS_lag: np.ndarray
    D_EBB_NOX: np.ndarray
    D_EBB_CO: np.ndarray
    D_EBB_VOC: np.ndarray
    D_EBB_CH4: np.ndarray
    D_EBB_SO2: np.ndarray
    D_EBB_BC: np.ndarray
    D_EBB_NH3: np.ndarray

    # FIXME: are those necessary?
    D_EESC = scalar_var()
    D_O3s = scalar_var()

    D_OHSNK_HFC_t = species_timeseries(nb_HFC)
    D_HVSNK_HFC_t = species_timeseries(nb_HFC)
    D_XSNK_HF_t = species_timeseries(nb_HFC)
    D_OHSNK_PFC_t = species_timeseries(nb_PFC)
    D_HVSNK_PFC_t = species_timeseries(nb_PFC)
    D_XSNK_PFC_t = species_timeseries(nb_PFC)
    D_OHSNK_ODS_t = species_timeseries(nb_ODS)
    D_HVSNK_ODS_t = species_timeseries(nb_ODS)
    D_XSNK_ODS_t = species_timeseries(nb_ODS)
    D_OHSNK_CH4_t = scalar_timeseries()
    D_HVSNK_CH4_t = scalar_timeseries()
    D_XSNK_CH4_t = scalar_timeseries()
    D_HVSNK_N2O_t = scalar_timeseries()
    D_kOH_t = scalar_timeseries()
    D_hv_t = scalar_timeseries()
    D_O3t_t = scalar_timeseries()
    D_EESC_t = scalar_timeseries()
    D_O3s_t = scalar_timeseries()
    D_SO4_t = scalar_timeseries()
    D_POA_t = scalar_timeseries()
    D_BC_t = scalar_timeseries()
    D_NO3_t = scalar_timeseries()
    D_SOA_t = scalar_timeseries()
    D_DUST_t = scalar_timeseries()
    D_SALT_t = scalar_timeseries()
    D_AERh_t = scalar_timeseries()

    def step_chemestry(self, t):
        # factors
        self.D_kOH = f_kOH(self.D_CH4, self.D_O3s, self.fT * self.D_gst, np.sum(ENOX[t] + self.D_EBB_NOX), np.sum(ECO[t] + self.D_EBB_CO), np.sum(EVOC[t] + self.D_EBB_VOC))
        self.D_hv = f_hv(self.D_N2O_lag, self.D_EESC, self.fT * self.D_gst)

        # fluxes
        self.D_OHSNK_CH4 = -alpha_CH4 / tau_CH4_OH * (CH4_0 * self.D_kOH + self.D_CH4 + self.D_kOH * self.D_CH4)
        self.D_HVSNK_CH4 = -alpha_CH4 / tau_CH4_hv * (CH4_0 * self.D_hv + self.D_CH4_lag + self.D_hv * self.D_CH4_lag)
        self.D_XSNK_CH4 = -alpha_CH4 * (1 / tau_CH4_soil + 1 / tau_CH4_ocean) * self.D_CH4
        self.D_FOXI_CH4 = -0.001 * (1.0 * np.sum(ECH4[t]) + np.sum(self.D_EBB_CH4) + np.sum(self.D_EWET) + self.D_OHSNK_CH4 + self.D_HVSNK_CH4 + self.D_XSNK_CH4)
        self.D_HVSNK_N2O = -alpha_N2O / tau_N2O_hv * (N2O_0 * self.D_hv + self.D_N2O_lag + self.D_hv * self.D_N2O_lag)
        self.D_OHSNK_HFC = -alpha_HFC / tau_HFC_OH * (HFC_0 * self.D_kOH + self.D_HFC + self.D_kOH * self.D_HFC)
        self.D_OHSNK_PFC = -alpha_PFC / tau_PFC_OH * (PFC_0 * self.D_kOH + self.D_PFC + self.D_kOH * self.D_PFC)
        self.D_OHSNK_ODS = -alpha_ODS / tau_ODS_OH * (ODS_0 * self.D_kOH + self.D_ODS + self.D_kOH * self.D_ODS)
        self.D_HVSNK_HFC = -alpha_HFC / tau_HFC_hv * (HFC_0 * self.D_hv + self.D_HFC_lag + self.D_hv * self.D_HFC_lag)
        self.D_HVSNK_PFC = -alpha_PFC / tau_PFC_hv * (PFC_0 * self.D_hv + self.D_PFC_lag + self.D_hv * self.D_PFC_lag)
        self.D_HVSNK_ODS = -alpha_ODS / tau_ODS_hv * (ODS_0 * self.D_hv + self.D_ODS_lag + self.D_hv * self.D_ODS_lag)
        self.D_XSNK_HFC = -alpha_HFC / tau_HFC_othr * self.D_HFC
        self.D_XSNK_PFC = -alpha_PFC / tau_PFC_othr * self.D_PFC
        self.D_XSNK_ODS = -alpha_ODS / tau_ODS_othr * self.D_ODS

        # stocks
        self.D_O3t = chi_O3t_CH4 * np.log(1 + self.D_CH4 / CH4_0) + Gamma_O3t * self.fT * self.D_gst
        self.D_O3t += chi_O3t_NOX * np.sum(w_reg_NOX * np.sum(p_reg4 * (ENOX[t] + self.D_EBB_NOX)[:, np.newaxis], 0))
        self.D_O3t += chi_O3t_CO * np.sum(w_reg_CO * np.sum(p_reg4 * (ECO[t] + self.D_EBB_CO)[:, np.newaxis], 0))
        self.D_O3t += chi_O3t_VOC * np.sum(w_reg_VOC * np.sum(p_reg4 * (EVOC[t] + self.D_EBB_VOC)[:, np.newaxis], 0))
        self.D_EESC = np.sum(f_fracrel(tau_lag) * (n_Cl + alpha_Br * n_Br) * self.D_ODS_lag)
        self.D_O3s = chi_O3s_EESC * self.D_EESC + chi_O3s_N2O * self.D_N2O_lag * (1 - self.D_EESC / EESC_x) + Gamma_O3s * self.fT * self.D_gst
        self.D_SO4 = alpha_SO4 * tau_SO2 * np.sum(w_reg_SO2 * np.sum(p_reg4 * (ESO2[t] + self.D_EBB_SO2)[:, np.newaxis], 0)) + alpha_SO4 * tau_DMS * 0 + Gamma_SO4 * self.fT * self.D_gst
        self.D_POA = tau_OMff * alpha_POM * np.sum(w_reg_OC * np.sum(p_reg4 * (EOC[t])[:, np.newaxis], 0)) + tau_OMbb * alpha_POM * np.sum(self.D_EBB_OC) + Gamma_POA * self.fT * self.D_gst
        self.D_BC = tau_BCff * np.sum(w_reg_BC * np.sum(p_reg4 * (EBC[t])[:, np.newaxis], 0)) + tau_BCbb * np.sum(self.D_EBB_BC) + Gamma_BC * self.fT * self.D_gst
        self.D_NO3 = alpha_NO3 * tau_NOX * np.sum(ENOX[t] + self.D_EBB_NOX) + alpha_NO3 * tau_NH3 * np.sum(ENH3[t] + self.D_EBB_NH3) + Gamma_NO3 * self.fT * self.D_gst
        self.D_SOA = tau_VOC * np.sum(EVOC[t] + self.D_EBB_VOC) + tau_BVOC * 0 + Gamma_SOA * self.fT * self.D_gst
        self.D_DUST = 0 * (tau_DUST * 0 + Gamma_DUST * self.fT * self.D_gst)
        self.D_SALT = 0 * (tau_SALT * 0 + Gamma_SALT * self.fT * self.D_gst)
        self.D_AERh = solub_SO4 * self.D_SO4 + solub_POA * self.D_POA + solub_BC * self.D_BC + solub_NO3 * self.D_NO3 + solub_SOA * self.D_SOA + solub_DUST * self.D_DUST + solub_SALT * self.D_SALT