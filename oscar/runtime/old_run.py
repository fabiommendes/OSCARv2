import matplotlib.pyplot as plt
import numpy as np

from ..historical import ODS_0, HFC_0, PFC_0
from ..constants import HFC, ODS, PFC
from .oscar_param import AERh_0, AREA_0, AWET_0, CFROZ_0, CH4_0, EESC_x, Gamma_BC, Gamma_DUST, Gamma_NO3, Gamma_O3s, Gamma_O3t, Gamma_POA, Gamma_SALT, Gamma_SO4, Gamma_SOA, N2O_0, Phi_0, alpha_BB_BC, alpha_BB_CH4, alpha_BB_CO, alpha_BB_CO2, alpha_BB_N2O, alpha_BB_NH3, alpha_BB_NOX, alpha_BB_OC, alpha_BB_SO2, alpha_BB_VOC, alpha_Br, alpha_CH4, alpha_CO2, alpha_HFC, alpha_LCC, alpha_N2O, alpha_NO3, alpha_ODS, alpha_OHC, alpha_PFC, alpha_POM, alpha_SO4, alpha_dic, alpha_gyp, alpha_mld, beta_gyp, chi_O3s_EESC, chi_O3s_N2O, chi_O3t_CH4, chi_O3t_CO, chi_O3t_NOX, chi_O3t_VOC, csoil1_0, csoil2_0, cveg_0, dty, ewet_0, f_RF_CH4, f_RF_CO2, f_RF_H2Os, f_RF_N2O, f_RF_overlap, f_fracrel, f_hv, f_kOH, f_npp, f_pCO2, f_pH, f_rho, f_v_PF, gamma_igniC, gamma_igniP, gamma_igniT, gamma_mld, gamma_pthaw, gamma_rhoPF1, gamma_rhoPF2, gamma_wetC, gamma_wetP, gamma_wetT, igni_0, ind_final, k_BC_adjust, k_met, k_pthaw, lambda_0, mld_0, mu_0, n_Br, n_Cl, nb_biome, nb_obox, nb_regionI, nb_regionPF, npp_0, p_AGB, p_HWP0, p_HWP1, p_HWP1_bb, p_HWP2, p_HWP3, p_OHC, p_PF1, p_PF2, p_PF3, p_PF_CH4, p_PF_inst, p_atm_CO2, p_atm_O3t, p_atm_absorb, p_atm_alb, p_atm_cloud, p_atm_noCO2, p_atm_scatter, p_atm_solar, p_atm_strat, p_circ, p_reg4, p_reg9, p_wet, pthaw_min, r_HWP1, r_HWP2, r_HWP3, radeff_BC, radeff_BCsnow, radeff_DUST, radeff_HFC, radeff_NO3, radeff_O3s, radeff_O3t, radeff_ODS, radeff_PFC, radeff_POA, radeff_SALT, radeff_SO4, radeff_SOA, rho1_0, rho2_0, solub_BC, solub_DUST, solub_NO3, solub_POA, solub_SALT, solub_SO4, solub_SOA, tau_BCbb, tau_BCff, tau_BVOC, tau_CH4_OH, tau_CH4_hv, tau_CH4_ocean, tau_CH4_soil, tau_DMS, tau_DUST, tau_HFC_OH, tau_HFC_hv, tau_HFC_othr, tau_N2O_hv, tau_NH3, tau_NOX, tau_ODS_OH, tau_ODS_hv, tau_ODS_othr, tau_OMbb, tau_OMff, tau_PF1, tau_PF2, tau_PF3, tau_PFC_OH, tau_PFC_hv, tau_PFC_othr, tau_SALT, tau_SO2, tau_VOC, tau_circ, tau_gst, tau_gst0, tau_lag, tau_shift, theta_0, v_fg, w_reg_BC, w_reg_BCsnow, w_reg_CO, w_reg_NOX, w_reg_OC, w_reg_SO2, w_reg_VOC, w_reg_lst, w_reg_lstPF, w_reg_lyp, w_reg_sst, w_rhoPF, warmeff_BCsnow, warmeff_LCC, warmeff_volc
from ..config import fT, p
from .drivers import EFF, ECH4, EN2O, LUC, HARV, SHIFT, EHFC, EPFC, EODS, ENOX, ECO, EVOC, ESO2, ENH3, EOC, EBC, RFcon, RFvolc, RFsolar, ECH4_0, EN2O_0, ENOX_0, ECO_0, EVOC_0, ESO2_0, ENH3_0, EOC_0, EBC_0
from .util import mk_scalar, mk_region_var, mk_region_biome_var, mk_region_biome2_age_var, mk_obox_var, mk_species_var, mk_regionPF_var, mk_land_var, mk_land_use_var, mk_scalar_var, mk_linear_var

nb_HFC = len(HFC)
nb_PFC = len(PFC)
nb_ODS = len(ODS)


# noinspection PyAttributeOutsideInit
class Simulation:
    # Global scalar variables
    D_mld_t = mk_scalar()
    D_dic_t = mk_scalar()
    D_pH_t = mk_scalar()
    OSNK_t = mk_scalar()
    LSNK_t = mk_scalar()
    D_FOXI_CH4_t = mk_scalar()
    D_OHSNK_CH4_t = mk_scalar()
    D_HVSNK_CH4_t = mk_scalar()
    D_XSNK_CH4_t = mk_scalar()
    D_HVSNK_N2O_t = mk_scalar()
    D_kOH_t = mk_scalar()
    D_hv_t = mk_scalar()
    D_O3t_t = mk_scalar()
    D_EESC_t = mk_scalar()
    D_O3s_t = mk_scalar()
    D_SO4_t = mk_scalar()
    D_POA_t = mk_scalar()
    D_BC_t = mk_scalar()
    D_NO3_t = mk_scalar()
    D_SOA_t = mk_scalar()
    D_AERh_t = mk_scalar()
    D_CO2_t = mk_scalar()
    D_CH4_t = mk_scalar()
    D_CH4_lag_t = mk_scalar()
    D_N2O_t = mk_scalar()
    D_N2O_lag_t = mk_scalar()
    RF_t = mk_scalar()
    RF_warm_t = mk_scalar()
    RF_atm_t = mk_scalar()
    RF_CO2_t = mk_scalar()
    RF_CH4_t = mk_scalar()
    RF_H2Os_t = mk_scalar()
    RF_N2O_t = mk_scalar()
    RF_halo_t = mk_scalar()
    RF_O3t_t = mk_scalar()
    RF_O3s_t = mk_scalar()
    RF_SO4_t = mk_scalar()
    RF_POA_t = mk_scalar()
    RF_BC_t = mk_scalar()
    RF_NO3_t = mk_scalar()
    RF_SOA_t = mk_scalar()
    RF_cloud_t = mk_scalar()
    RF_BCsnow_t = mk_scalar()
    RF_LCC_t = mk_scalar()
    D_gst_t = mk_scalar()
    D_sst_t = mk_scalar()
    D_gyp_t = mk_scalar()
    D_OHC_t = mk_scalar()

    # (region) variables
    ELUC_t = mk_region_var()
    D_AWET_t = mk_region_var()
    D_EWET_t = mk_region_var()
    D_ewet_t = mk_region_var()
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
    D_lst_t = mk_region_var()
    D_lyp_t = mk_region_var()

    # (region)*(biome) variables
    D_AREA_t = mk_region_biome_var()
    D_npp_t = mk_region_biome_var()
    D_efire_t = mk_region_biome_var()
    D_fmort_t = mk_region_biome_var()
    D_rh1_t = mk_region_biome_var()
    D_fmet_t = mk_region_biome_var()
    D_rh2_t = mk_region_biome_var()
    D_cveg_t = mk_region_biome_var()
    D_csoil1_t = mk_region_biome_var()
    D_csoil2_t = mk_region_biome_var()

    # (region)*(biome)*(biome)*(age) variables
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

    # (obox) variables
    D_FIN_t = mk_obox_var()
    D_FOUT_t = mk_obox_var()
    D_FCIRC_t = mk_obox_var()
    D_CSURF_t = mk_obox_var()

    # (species) variables
    D_HFC_t = mk_species_var(nb_HFC)
    D_HFC_lag_t = mk_species_var(nb_HFC)
    D_OHSNK_HFC_t = mk_species_var(nb_HFC)
    D_HVSNK_HFC_t = mk_species_var(nb_HFC)
    D_XSNK_HF_t = mk_species_var(nb_HFC)

    D_PFC_t = mk_species_var(nb_PFC)
    D_PFC_lag_t = mk_species_var(nb_PFC)
    D_OHSNK_PFC_t = mk_species_var(nb_PFC)
    D_HVSNK_PFC_t = mk_species_var(nb_PFC)
    D_XSNK_PFC_t = mk_species_var(nb_PFC)

    D_ODS_t = mk_species_var(nb_ODS)
    D_ODS_lag_t = mk_species_var(nb_ODS)
    D_OHSNK_ODS_t = mk_species_var(nb_ODS)
    D_HVSNK_ODS_t = mk_species_var(nb_ODS)
    D_XSNK_ODS_t = mk_species_var(nb_ODS)

    # (regionPF) variables
    pthaw_t = mk_regionPF_var()
    pthaw_bar_t = mk_regionPF_var()
    FTHAW_t = mk_regionPF_var()
    ETHAW1_t = mk_regionPF_var()
    ETHAW2_t = mk_regionPF_var()
    ETHAW3_t = mk_regionPF_var()
    EPF_CO2_t = mk_regionPF_var()
    EPF_CH4_t = mk_regionPF_var()
    EPF_t = mk_regionPF_var()
    CTHAW1_t = mk_regionPF_var()
    CTHAW2_t = mk_regionPF_var()
    CTHAW3_t = mk_regionPF_var()
    D_CFROZ_t = mk_regionPF_var()

    # Land
    D_AREA = mk_land_var()
    D_cveg = mk_land_var()
    D_csoil1 = mk_land_var()
    D_csoil2 = mk_land_var()

    # Land use
    CVEG_luc = mk_land_use_var()
    CSOIL1_luc = mk_land_use_var()
    CSOIL2_luc = mk_land_use_var()
    CHWP1_luc = mk_land_use_var()
    CHWP2_luc = mk_land_use_var()
    CHWP3_luc = mk_land_use_var()

    # atmosphere
    D_CO2 = mk_scalar_var()
    D_CH4 = mk_scalar_var()
    D_CH4_lag = mk_scalar_var()
    D_N2O = mk_scalar_var()
    D_N2O_lag = mk_scalar_var()
    D_EESC = mk_scalar_var()
    D_O3s = mk_scalar_var()

    D_HFC = mk_linear_var(nb_HFC)
    D_HFC_lag = mk_linear_var(nb_HFC)
    D_PFC = mk_linear_var(nb_PFC)
    D_PFC_lag = mk_linear_var(nb_PFC)
    D_ODS = mk_linear_var(nb_ODS)
    D_ODS_lag = mk_linear_var(nb_ODS)

    # Climate
    D_gst = mk_scalar_var()
    D_gst0 = mk_scalar_var()
    D_sst = mk_scalar_var()
    D_gyp = mk_scalar_var()
    D_OHC = mk_scalar_var()
    D_lst = mk_linear_var(nb_regionI)
    D_lyp = mk_linear_var(nb_regionI)

    # Permafrost
    pthaw = mk_linear_var(nb_regionPF)
    CTHAW1 = mk_linear_var(nb_regionPF)
    CTHAW2 = mk_linear_var(nb_regionPF)
    CTHAW3 = mk_linear_var(nb_regionPF)
    D_CFROZ = mk_linear_var(nb_regionPF)

    def __init__(self, p=p, fT=fT, outputs=["ELUC", "OSNK", "LSNK", "D_CO2", "RF", "D_gst"], plot='all'):
        var_plot = self.plot_vars(plot)
        self.dt = dt = 1 / p

        # ocean
        self.D_dic = np.array([0], dtype=dty)
        self.D_CSURF = np.zeros([nb_obox], dtype=dty)

        # save variables
        self.var_timeseries = list(set(outputs) | set(var_plot))

        # =======
        # B. RUN
        # =======
        print('STARTING SIMULATION')
        print('\n'.join(self.var_timeseries))

        for t in range(1, ind_final + 1):
            for tt in range(p):
                self.step_ocean(t)
                self.step_land(t)
                self.step_land_use(t)
                self.step_permafrost(t)
                self.step_chemestry(t)
                self.step_atmosphere(t)
                self.step_climate(t)

                # save time series variables
                loc = locals()
                for name in self.var_timeseries:
                    arr = getattr(self, name + '_t')
                    try:
                        run = loc[name]
                    except KeyError:
                        run = getattr(self, name)
                    arr[t] += dt * run
                del loc

        # FIGURES
        self.plot(plot)

        # OUTPUTS
        output = []
        for name in outputs:
            output.append(getattr(self, name + '_t'))
        self.output = output

    def step_ocean(self, t):
        dt = self.dt

        # structure
        D_mld = mld_0 * alpha_mld * (np.exp(gamma_mld * fT * self.D_sst) - 1)

        # fluxes
        self.D_FIN = p_circ * v_fg * alpha_CO2 * self.D_CO2
        self.D_FOUT = p_circ * v_fg * alpha_CO2 * f_pCO2(self.D_dic, fT * self.D_sst)
        self.D_FCIRC = self.D_CSURF * (1 / tau_circ)
        self.OSNK = np.sum(self.D_FOUT - self.D_FIN)

        # stocks
        self.D_CSURF += dt * (self.D_FIN - self.D_FOUT - self.D_FCIRC)
        self.D_dic = alpha_dic * np.sum(self.D_CSURF) / (1 + D_mld / mld_0)

    def step_land(self, t):
        dt = self.dt

        # land-cover
        self.D_AREA += dt * (np.sum(LUC[t], 1) - np.sum(LUC[t], 2))
        self.D_AWET = AWET_0 * (gamma_wetT * fT * self.D_lst + gamma_wetP * fT * self.D_lyp + gamma_wetC * fT * self.D_CO2)

        # factors
        self.D_k_igni = (gamma_igniT * fT * self.D_lst[:, np.newaxis] + gamma_igniP * fT * self.D_lyp[:, np.newaxis] + gamma_igniC * fT * self.D_CO2)
        self.D_k_rho = f_rho(fT * self.D_lst[:, np.newaxis], fT * self.D_lyp[:, np.newaxis])

        # fluxes
        self.D_npp = npp_0 * f_npp(self.D_CO2, fT * self.D_lst[:, np.newaxis], fT * self.D_lyp[:, np.newaxis])
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

    def step_land_use(self, t):
        dt = self.dt
        # initialization
        # land-use change
        for b1 in range(nb_biome):
            for b2 in range(nb_biome):
                self.CVEG_luc[:, b1, b2, t] += (-dt) * (cveg_0 + self.D_cveg)[:, b2] * LUC[t, :, b1, b2]
                self.CSOIL1_luc[:, b1, b2, t] += (dt * ((csoil1_0 + self.D_csoil1)[:, b1] - (csoil1_0 + self.D_csoil1)[:, b2]) * LUC[t, :, b1, b2])
                self.CSOIL2_luc[:, b1, b2, t] += (dt * ((csoil2_0 + self.D_csoil2)[:, b1] - (csoil2_0 + self.D_csoil2)[:, b2]) * LUC[t, :, b1, b2])
                self.CSOIL1_luc[:, b1, b2, t] += (dt * (cveg_0 + self.D_cveg)[:, b1] * (p_AGB[:, b1] * p_HWP0[:, b1] + (1 - p_AGB[:, b1])) * LUC[t, :, b1, b2])
                self.CHWP1_luc[:, b1, b2, t] += (dt * (cveg_0 + self.D_cveg)[:, b1] * p_AGB[:, b1] * p_HWP1[:, b1] * LUC[t, :, b1, b2])
                self.CHWP2_luc[:, b1, b2, t] += (dt * (cveg_0 + self.D_cveg)[:, b1] * p_AGB[:, b1] * p_HWP2[:, b1] * LUC[t, :, b1, b2])
                self.CHWP3_luc[:, b1, b2, t] += (dt * (cveg_0 + self.D_cveg)[:, b1] * p_AGB[:, b1] * p_HWP3[:, b1] * LUC[t, :, b1, b2])

        # harvest
        for b in range(nb_biome):
            self.CVEG_luc[:, b, b, t] += dt * -HARV[t, :, b]
            self.CSOIL1_luc[:, b, b, t] += dt * p_HWP0[:, b] * HARV[t, :, b]
            self.CHWP1_luc[:, b, b, t] += dt * p_HWP1[:, b] * HARV[t, :, b]
            self.CHWP2_luc[:, b, b, t] += dt * p_HWP2[:, b] * HARV[t, :, b]
            self.CHWP3_luc[:, b, b, t] += dt * p_HWP3[:, b] * HARV[t, :, b]

        # shifting cultivation
        for b1 in range(nb_biome):
            for b2 in range(b1, nb_biome):
                self.CVEG_luc[:, b1, b2, t] += (dt * -(cveg_0 + self.D_cveg)[:, b2] * (1 - np.exp(-mu_0[:, b2] * tau_shift)) * SHIFT[t, :, b1, b2])
                self.CSOIL1_luc[:, b1, b2, t] += (dt * (cveg_0 + self.D_cveg)[:, b1] * (1 - np.exp(-mu_0[:, b1] * tau_shift)) * (p_AGB[:, b1] * p_HWP0[:, b1] + (1 - p_AGB[:, b1])) * SHIFT[t, :, b1, b2])
                self.CHWP1_luc[:, b1, b2, t] += (dt * (cveg_0 + self.D_cveg)[:, b1] * (1 - np.exp(-mu_0[:, b1] * tau_shift)) * p_AGB[:, b1] * p_HWP1[:, b1] * SHIFT[t, :, b1, b2])
                self.CHWP2_luc[:, b1, b2, t] += (dt * (cveg_0 + self.D_cveg)[:, b1] * (1 - np.exp(-mu_0[:, b1] * tau_shift)) * p_AGB[:, b1] * p_HWP2[:, b1] * SHIFT[t, :, b1, b2])
                self.CHWP3_luc[:, b1, b2, t] += (dt * (cveg_0 + self.D_cveg)[:, b1] * (1 - np.exp(-mu_0[:, b1] * tau_shift)) * p_AGB[:, b1] * p_HWP3[:, b1] * SHIFT[t, :, b1, b2])
                self.CVEG_luc[:, b2, b1, t] += (dt * -(cveg_0 + self.D_cveg)[:, b1] * (1 - np.exp(-mu_0[:, b1] * tau_shift)) * SHIFT[t, :, b1, b2])
                self.CSOIL1_luc[:, b2, b1, t] += (dt * (cveg_0 + self.D_cveg)[:, b2] * (1 - np.exp(-mu_0[:, b2] * tau_shift)) * (p_AGB[:, b2] * p_HWP0[:, b2] + (1 - p_AGB[:, b2])) * SHIFT[t, :, b1, b2])
                self.CHWP1_luc[:, b2, b1, t] += (dt * (cveg_0 + self.D_cveg)[:, b2] * (1 - np.exp(-mu_0[:, b2] * tau_shift)) * p_AGB[:, b2] * p_HWP1[:, b2] * SHIFT[t, :, b1, b2])
                self.CHWP2_luc[:, b2, b1, t] += (dt * (cveg_0 + self.D_cveg)[:, b2] * (1 - np.exp(-mu_0[:, b2] * tau_shift)) * p_AGB[:, b2] * p_HWP2[:, b2] * SHIFT[t, :, b1, b2])
                self.CHWP3_luc[:, b2, b1, t] += (dt * (cveg_0 + self.D_cveg)[:, b2] * (1 - np.exp(-mu_0[:, b2] * tau_shift)) * p_AGB[:, b2] * p_HWP3[:, b2] * SHIFT[t, :, b1, b2])

        # fluxes
        # book-keeping model
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

        # biomass burning
        def biomass_burning_diff(alpha_BB):
            new = np.newaxis
            diff = np.sum(alpha_BB * (igni_0 * cveg_0 * self.D_AREA + self.D_efire * AREA_0 + self.D_efire * self.D_AREA), 1)
            diff += p_HWP1_bb * np.sum(np.sum(np.sum(alpha_BB[:, :, new, new] * self.EHWP1_luc, 3), 2), 1)
            diff += np.sum(np.sum(np.sum(alpha_BB[:, new, :, new] * self.EFIRE_luc, 3), 2), 1)
            return diff

        self.D_EBB_CO2 = biomass_burning_diff(alpha_BB_CO2)
        self.D_EBB_CH4 = biomass_burning_diff(alpha_BB_CH4)
        self.D_EBB_N2O = biomass_burning_diff(alpha_BB_N2O)
        self.D_EBB_NOX = biomass_burning_diff(alpha_BB_NOX)
        self.D_EBB_CO = biomass_burning_diff(alpha_BB_CO)
        self.D_EBB_VOC = biomass_burning_diff(alpha_BB_VOC)
        self.D_EBB_SO2 = biomass_burning_diff(alpha_BB_SO2)
        self.D_EBB_NH3 = biomass_burning_diff(alpha_BB_NH3)
        self.D_EBB_OC = biomass_burning_diff(alpha_BB_OC)
        self.D_EBB_BC = biomass_burning_diff(alpha_BB_BC)

        # stocks
        self.CVEG_luc += dt * (self.NPP_luc - self.FMORT_luc - self.EFIRE_luc)
        self.CSOIL1_luc += dt * (self.FMORT_luc - self.FMET_luc - self.RH1_luc)
        self.CSOIL2_luc += dt * (self.FMET_luc - self.RH2_luc)
        self.CHWP1_luc += dt * -self.EHWP1_luc
        self.CHWP2_luc += dt * -self.EHWP2_luc
        self.CHWP3_luc += dt * -self.EHWP3_luc

    def step_permafrost(self, t):
        dt = self.dt

        # factors
        self.rD_rhoPF = np.exp(w_rhoPF * gamma_rhoPF1 * w_reg_lstPF * fT * self.D_gst + w_rhoPF * gamma_rhoPF2 * (w_reg_lstPF * fT * self.D_gst) ** 2) - 1

        # fraction thawed
        self.pthaw_bar = -pthaw_min + (1 + pthaw_min) / (1 + ((1 / pthaw_min + 1) ** k_pthaw - 1) * np.exp(-gamma_pthaw * k_pthaw * w_reg_lstPF * fT * self.D_gst)) ** (1 / k_pthaw)
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

    def step_chemestry(self, t):
        dt = self.dt

        # factors
        self.D_kOH = f_kOH(self.D_CH4, self.D_O3s, fT * self.D_gst, np.sum(ENOX[t] + self.D_EBB_NOX), np.sum(ECO[t] + self.D_EBB_CO), np.sum(EVOC[t] + self.D_EBB_VOC))
        self.D_hv = f_hv(self.D_N2O_lag, self.D_EESC, fT * self.D_gst)

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
        self.D_O3t = chi_O3t_CH4 * np.log(1 + self.D_CH4 / CH4_0) + Gamma_O3t * fT * self.D_gst
        self.D_O3t += chi_O3t_NOX * np.sum(w_reg_NOX * np.sum(p_reg4 * (ENOX[t] + self.D_EBB_NOX)[:, np.newaxis], 0))
        self.D_O3t += chi_O3t_CO * np.sum(w_reg_CO * np.sum(p_reg4 * (ECO[t] + self.D_EBB_CO)[:, np.newaxis], 0))
        self.D_O3t += chi_O3t_VOC * np.sum(w_reg_VOC * np.sum(p_reg4 * (EVOC[t] + self.D_EBB_VOC)[:, np.newaxis], 0))
        self.D_EESC = np.sum(f_fracrel(tau_lag) * (n_Cl + alpha_Br * n_Br) * self.D_ODS_lag)
        self.D_O3s = chi_O3s_EESC * self.D_EESC + chi_O3s_N2O * self.D_N2O_lag * (1 - self.D_EESC / EESC_x) + Gamma_O3s * fT * self.D_gst
        self.D_SO4 = alpha_SO4 * tau_SO2 * np.sum(w_reg_SO2 * np.sum(p_reg4 * (ESO2[t] + self.D_EBB_SO2)[:, np.newaxis], 0)) + alpha_SO4 * tau_DMS * 0 + Gamma_SO4 * fT * self.D_gst
        self.D_POA = tau_OMff * alpha_POM * np.sum(w_reg_OC * np.sum(p_reg4 * (EOC[t])[:, np.newaxis], 0)) + tau_OMbb * alpha_POM * np.sum(self.D_EBB_OC) + Gamma_POA * fT * self.D_gst
        self.D_BC = tau_BCff * np.sum(w_reg_BC * np.sum(p_reg4 * (EBC[t])[:, np.newaxis], 0)) + tau_BCbb * np.sum(self.D_EBB_BC) + Gamma_BC * fT * self.D_gst
        self.D_NO3 = alpha_NO3 * tau_NOX * np.sum(ENOX[t] + self.D_EBB_NOX) + alpha_NO3 * tau_NH3 * np.sum(ENH3[t] + self.D_EBB_NH3) + Gamma_NO3 * fT * self.D_gst
        self.D_SOA = tau_VOC * np.sum(EVOC[t] + self.D_EBB_VOC) + tau_BVOC * 0 + Gamma_SOA * fT * self.D_gst
        self.D_DUST = 0 * (tau_DUST * 0 + Gamma_DUST * fT * self.D_gst)
        self.D_SALT = 0 * (tau_SALT * 0 + Gamma_SALT * fT * self.D_gst)
        self.D_AERh = solub_SO4 * self.D_SO4 + solub_POA * self.D_POA + solub_BC * self.D_BC + solub_NO3 * self.D_NO3 + solub_SOA * self.D_SOA + solub_DUST * self.D_DUST + solub_SALT * self.D_SALT

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
        self.RF = self.RF_CO2 + self.RF_CH4 + self.RF_H2Os + self.RF_N2O + self.RF_halo + self.RF_O3t + self.RF_O3s + self.RF_SO4 + self.RF_POA + self.RF_BC + self.RF_NO3 + self.RF_SOA + self.RF_DUST + self.RF_SALT + self.RF_cloud + self.RF_BCsnow + self.RF_LCC + RFcon[t] + RFvolc[t] + RFsolar[t]
        self.RF_warm = self.RF_CO2 + self.RF_CH4 + self.RF_H2Os + self.RF_N2O + self.RF_halo + self.RF_O3t + self.RF_O3s + self.RF_SO4 + self.RF_POA + self.RF_BC + self.RF_NO3 + self.RF_SOA + self.RF_DUST + self.RF_SALT + self.RF_cloud + warmeff_BCsnow * self.RF_BCsnow + warmeff_LCC * self.RF_LCC + RFcon[t] + warmeff_volc * RFvolc[t] + RFsolar[t]
        self.RF_atm = p_atm_CO2 * self.RF_CO2 + p_atm_noCO2 * (self.RF_CH4 + self.RF_N2O + self.RF_halo) + p_atm_O3t * self.RF_O3t + p_atm_strat * (self.RF_O3s + self.RF_H2Os) + p_atm_scatter * (self.RF_SO4 + self.RF_POA + self.RF_NO3 + self.RF_SOA + self.RF_DUST + self.RF_SALT + RFvolc[t]) + p_atm_absorb * self.RF_BC + p_atm_cloud * (self.RF_cloud + RFcon[t]) + p_atm_alb * (self.RF_BCsnow + self.RF_LCC) + p_atm_solar * RFsolar[t]

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

    def plot_vars(self, plot):
        # plot variables
        var_plot = []
        if plot == "all" or plot == "CO2" or "CO2" in plot:
            var_plot += ["D_CO2", "OSNK", "LSNK", "ELUC", "D_AREA", "D_npp", "D_efire", "D_fmort", "D_rh1", "D_fmet", "D_rh2", "D_FIN", "D_FOUT", "D_FCIRC", "EFIRE_luc", "FMORT_luc", "RH1_luc", "FMET_luc", "RH2_luc", "EHWP1_luc", "EHWP2_luc", "EHWP3_luc", ]
        if plot == "all" or plot == "CH4" or "CH4" in plot:
            var_plot += ["D_CH4", "D_OHSNK_CH4", "D_HVSNK_CH4", "D_XSNK_CH4", "D_EWET", "D_EBB_CH4"]
        if plot == "all" or plot == "N2O" or "N2O" in plot:
            var_plot += ["D_N2O", "D_HVSNK_N2O", "D_EBB_N2O"]
        if plot == "all" or plot == "O3" or "O3" in plot:
            var_plot += ["D_O3t", "D_O3s", "D_EESC", "D_N2O_lag", "D_gst"]
        if plot == "all" or plot == "AER" or "AER" in plot:
            var_plot += ["D_SO4", "D_POA", "D_BC", "D_NO3", "D_SOA", "D_AERh", "RF_SO4", "RF_POA", "RF_BC", "RF_NO3", "RF_SOA", "RF_cloud", ]
        if plot == "all" or plot == "clim" or "clim" in plot:
            var_plot += ["RF", "D_gst", "D_gyp", "RF_CO2", "RF_CH4", "RF_H2Os", "RF_N2O", "RF_halo", "RF_O3t", "RF_O3s", "RF_SO4", "RF_POA", "RF_BC", "RF_NO3", "RF_SOA", "RF_cloud", "RF_BCsnow", "RF_LCC"]
        return var_plot

    def plot(self, plot):
        from ..plot import plot_AER, plot_CH4, plot_clim, plot_CO2, plot_N2O, plot_O3

        if plot == "all" or plot == "CO2" or "CO2" in plot:
            plot_CO2(self.D_CO2_t, self.OSNK_t, self.LSNK_t, self.ELUC_t, EFF, self.D_AREA_t, self.D_npp_t, self.D_efire_t, self.D_fmort_t, self.D_rh1_t, self.D_fmet_t, self.D_rh2_t, self.D_FIN_t, self.D_FOUT_t, self.D_FCIRC_t, self.EFIRE_luc_t, self.FMORT_luc_t, self.RH1_luc_t, self.RH2_luc_t, self.EHWP1_luc_t, self.EHWP2_luc_t, self.EHWP3_luc_t)
            plt.savefig("results/plot-CO2.svg")
        if plot == "all" or plot == "CH4" or "CH4" in plot:
            plot_CH4(self.D_CH4_t, self.D_OHSNK_CH4_t, self.D_HVSNK_CH4_t, self.D_XSNK_CH4_t, self.D_EWET_t, self.D_EBB_CH4_t, ECH4)
            plt.savefig("results/plot-CH4.svg")
        if plot == "all" or plot == "N2O" or "N2O" in plot:
            plot_N2O(self.D_N2O_t, self.D_HVSNK_N2O_t, self.D_EBB_N2O_t, EN2O)
            plt.savefig("results/plot-N20.svg")
        if plot == "all" or plot == "O3" or "O3" in plot:
            plot_O3(self.D_O3t_t, self.D_O3s_t, self.D_EESC_t, self.D_N2O_lag_t, self.D_gst_t)
            plt.savefig("results/plot-03.svg")
        if plot == "all" or plot == "AER" or "AER" in plot:
            plot_AER(self.D_SO4_t, self.D_POA_t, self.D_BC_t, self.D_NO3_t, self.D_SOA_t, self.D_AERh_t, self.RF_SO4_t, self.RF_POA_t, self.RF_BC_t, self.RF_NO3_t, self.RF_SOA_t, self.RF_cloud_t)
            plt.savefig("results/plot-AER.svg")
        if plot == "all" or plot == "clim" or "clim" in plot:
            plot_clim(self.RF_t, self.D_gst_t, self.D_gyp_t, self.RF_CO2_t, self.RF_CH4_t, self.RF_H2Os_t, self.RF_N2O_t, self.RF_halo_t, self.RF_O3t_t, self.RF_O3s_t, self.RF_SO4_t, self.RF_POA_t, self.RF_BC_t, self.RF_NO3_t, self.RF_SOA_t, self.RF_cloud_t, self.RF_BCsnow_t, self.RF_LCC_t, RFcon, RFvolc, RFsolar)
            plt.savefig("results/plot-clim.svg")


def OSCAR_lite(var_output=["D_CO2", "D_CH4", "D_N2O", "RF_halo", "D_O3t", "D_O3s", "D_SO4", "D_POA", "D_BC", "D_NO3", "D_SOA", "D_AERh", "RF", "D_gst"], plot="all"):
    sim = Simulation(outputs=var_output, plot=plot)
    # all_vars = set(var_output) | set(sim.plot_variables(plot))
    # print('STARTING SIMULATION')
    # sim.track(all_vars)
    # sim.run()
    # print('PLOTING VARIABLES', plot)
    # sim.plot(plot)
    return sim.output