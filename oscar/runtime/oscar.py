import matplotlib.pyplot as plt
import numpy as np

from .oscar_param import *
from .oscar_param import nb_regionI, nb_biome, nb_HFC, radeff_NO3, tau_OMff, f_RF_CO2, nb_obox
from ..config import fT, p, dty, ind_final


def reduce_driver(driver):
    return np.sum(np.sum(np.sum(driver, 3), 2), 1)


def reduce_param(param):
    return np.sum(np.sum(np.sum(param, 2), 1), 0)


def OSCAR_lite(
        p=p,
        fT=fT,
        force_CO2=False,
        force_GHG=False,
        force_halo=False,
        force_RF=False,
        force_RFs=False,
        force_clim=False,
        var_output=["ELUC", "OSNK", "LSNK", "D_CO2", "RF", "D_gst"],
        plot=[],
        D_CO2_force=None,
        D_CH4_force=None,
        D_N2O_force=None,
        D_HFC_force=None,
        D_PFC_force=None,
        D_ODS_force=None,
        RF_CO2_force=None,
        RF_CH4_force=None,
        RF_H2Os_force=None,
        RF_N2O_force=None,
        RF_halo_force=None,
        RF_O3t_force=None,
        RF_O3s_force=None,
        RF_SO4_force=None,
        RF_POA_force=None,
        RF_BC_force=None,
        RF_NO3_force=None,
        RF_SOA_force=None,
        RF_DUST_force=None,
        RF_SALT_force=None,
        RF_cloud_force=None,
        RF_BCsnow_force=None,
        RF_LCC_force=None,
        D_gst_force=None,
        D_sst_force=None,
        D_lst_force=None,
        D_lyp_force=None,
        RF_force=None,
):
    # ===============
    # A. DEFINITIONS
    # ===============

    # plot variables
    var_plot = []
    if plot == "all" or plot == "CO2" or "CO2" in plot:
        var_plot += ["D_CO2", "OSNK", "LSNK", "ELUC", "D_AREA", "D_npp", "D_efire", "D_fmort", "D_rh1", "D_fmet",
                     "D_rh2", "D_FIN", "D_FOUT", "D_FCIRC", "EFIRE_luc", "FMORT_luc", "RH1_luc", "FMET_luc",
                     "RH2_luc", "EHWP1_luc", "EHWP2_luc", "EHWP3_luc", ]
    if plot == "all" or plot == "CH4" or "CH4" in plot:
        var_plot += ["D_CH4", "D_OHSNK_CH4", "D_HVSNK_CH4", "D_XSNK_CH4", "D_EWET", "D_EBB_CH4"]
    if plot == "all" or plot == "N2O" or "N2O" in plot:
        var_plot += ["D_N2O", "D_HVSNK_N2O", "D_EBB_N2O"]
    if plot == "all" or plot == "O3" or "O3" in plot:
        var_plot += ["D_O3t", "D_O3s", "D_EESC", "D_N2O_lag", "D_gst"]
    if plot == "all" or plot == "AER" or "AER" in plot:
        var_plot += ["D_SO4", "D_POA", "D_BC", "D_NO3", "D_SOA", "D_AERh", "RF_SO4", "RF_POA", "RF_BC", "RF_NO3",
                     "RF_SOA", "RF_cloud", ]
    if plot == "all" or plot == "clim" or "clim" in plot:
        var_plot += ["RF", "D_gst", "D_gyp", "RF_CO2", "RF_CH4", "RF_H2Os", "RF_N2O", "RF_halo", "RF_O3t", "RF_O3s",
                     "RF_SO4", "RF_POA", "RF_BC", "RF_NO3", "RF_SOA", "RF_cloud", "RF_BCsnow", "RF_LCC"]

    #
    # Global scalar variables
    #
    mk_scalar = lambda: np.zeros([ind_final + 1], dtype=dty)
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

    #
    # (region) variables
    #
    mk_region_var = lambda: np.zeros([ind_final + 1, nb_regionI], dtype=dty)
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
    mk_region_biome_var = lambda: np.zeros([ind_final + 1, nb_regionI, nb_biome], dtype=dty)
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

    #
    # (region)*(biome)*(biome)*(age) variables
    #
    mk_region_biome2_age_var = lambda: np.zeros([ind_final + 1, nb_regionI, nb_biome, nb_biome, ind_final + 1],
                                                dtype=dty)
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
    mk_obox_var = lambda: np.zeros([ind_final + 1, nb_obox], dtype=dty)
    D_FIN_t = mk_obox_var()
    D_FOUT_t = mk_obox_var()
    D_FCIRC_t = mk_obox_var()
    D_CSURF_t = mk_obox_var()

    # (species) variables
    mk_species_var = lambda n: np.zeros([ind_final + 1, n], dtype=dty)
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
    mk_regionPF_var = lambda: np.zeros([ind_final + 1, nb_regionPF], dtype=dty)
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

    # Run variables
    # ocean
    D_dic = np.array([0], dtype=dty)
    D_CSURF = np.zeros([nb_obox], dtype=dty)

    # land
    mk_land_var = lambda: np.zeros([nb_regionI, nb_biome], dtype=dty)
    D_AREA = mk_land_var()
    D_cveg = mk_land_var()
    D_csoil1 = mk_land_var()
    D_csoil2 = mk_land_var()

    # land use
    mk_land_use_var = lambda: np.zeros([nb_regionI, nb_biome, nb_biome, ind_final + 1], dtype=dty)
    CVEG_luc = mk_land_use_var()
    CSOIL1_luc = mk_land_use_var()
    CSOIL2_luc = mk_land_use_var()
    CHWP1_luc = mk_land_use_var()
    CHWP2_luc = mk_land_use_var()
    CHWP3_luc = mk_land_use_var()

    # atmosphere
    mk_scalar_var = lambda value=0: np.array([value], dtype=dty)
    D_CO2 = mk_scalar_var()
    D_CH4 = mk_scalar_var()
    D_CH4_lag = mk_scalar_var()
    D_N2O = mk_scalar_var()
    D_N2O_lag = mk_scalar_var()
    D_EESC = mk_scalar_var()
    D_O3s = mk_scalar_var()

    mk_linear_var = lambda n: np.zeros([n], dtype=dty)
    D_HFC = mk_linear_var(nb_HFC)
    D_HFC_lag = mk_linear_var(nb_HFC)
    D_PFC = mk_linear_var(nb_PFC)
    D_PFC_lag = mk_linear_var(nb_PFC)
    D_ODS = mk_linear_var(nb_ODS)
    D_ODS_lag = mk_linear_var(nb_ODS)

    # climate
    D_gst = mk_scalar_var()
    D_gst0 = mk_scalar_var()
    D_sst = mk_scalar_var()
    D_gyp = mk_scalar_var()
    D_OHC = mk_scalar_var()

    D_lst = mk_linear_var(nb_regionI)
    D_lyp = mk_linear_var(nb_regionI)

    # permafrost
    pthaw = mk_linear_var(nb_regionPF)
    CTHAW1 = mk_linear_var(nb_regionPF)
    CTHAW2 = mk_linear_var(nb_regionPF)
    CTHAW3 = mk_linear_var(nb_regionPF)
    D_CFROZ = mk_linear_var(nb_regionPF)

    # Forces
    D_CO2_force = D_CO2 * 0 if D_CO2_force is None else D_CO2_force
    D_CO2_force = D_CO2 * 0 if D_CO2_force is None else D_CO2_force
    D_CH4_force = D_CH4 * 0 if D_CH4_force is None else D_CH4_force
    D_N2O_force = D_N2O * 0 if D_N2O_force is None else D_N2O_force
    D_HFC_force = D_HFC * 0 if D_HFC_force is None else D_HFC_force
    D_PFC_force = D_PFC * 0 if D_PFC_force is None else D_PFC_force
    D_ODS_force = D_ODS * 0 if D_ODS_force is None else D_ODS_force

    RF_force = mk_scalar_var() if RF_force is None else RF_force
    RF_CO2_force = D_CO2 * 0 if RF_CO2_force is None else RF_CO2_force
    RF_CH4_force = D_CH4 * 0 if RF_CH4_force is None else RF_CH4_force
    RF_H2Os_force = mk_scalar_var() if RF_H2Os_force is None else RF_H2Os_force
    RF_N2O_force = D_N2O * 0 if RF_N2O_force is None else RF_N2O_force
    RF_halo_force = mk_scalar_var() if RF_halo_force is None else RF_halo_force
    RF_O3t_force = mk_scalar_var() if RF_O3t_force is None else RF_O3t_force
    RF_O3s_force = D_O3s * 0 if RF_O3s_force is None else RF_O3s_force
    RF_SO4_force = mk_scalar_var() if RF_SO4_force is None else RF_SO4_force
    RF_POA_force = mk_scalar_var() if RF_POA_force is None else RF_POA_force
    RF_BC_force = mk_scalar_var() if RF_BC_force is None else RF_BC_force
    RF_NO3_force = mk_scalar_var() if RF_NO3_force is None else RF_NO3_force
    RF_SOA_force = mk_scalar_var() if RF_SOA_force is None else RF_SOA_force
    RF_DUST_force = mk_scalar_var() if RF_DUST_force is None else RF_DUST_force
    RF_SALT_force = mk_scalar_var() if RF_SALT_force is None else RF_SALT_force
    RF_cloud_force = mk_scalar_var() if RF_cloud_force is None else RF_cloud_force
    RF_BCsnow_force = mk_scalar_var() if RF_BCsnow_force is None else RF_BCsnow_force
    RF_LCC_force = mk_scalar_var() if RF_LCC_force is None else RF_LCC_force

    D_gst_force = D_gst if D_gst_force is None else D_gst_force
    D_sst_force = D_sst if D_sst_force is None else D_sst_force
    D_lst_force = D_lst if D_lst_force is None else D_lst_force
    D_lyp_force = D_lyp if D_lyp_force is None else D_lyp_force

    # save variables
    loc = locals()
    var_timeseries = list(set(var_output) | set(var_plot))
    values_timeseries = [loc[name + '_t'] for name in var_timeseries]
    del loc

    # =======
    # B. RUN
    # =======
    dt = 1 / p
    
    for t in range(1, ind_final + 1):
        for tt in range(p):

            # ---------
            # 1. OCEAN
            # ---------

            # structure
            D_mld = mld_0 * alpha_mld * (np.exp(gamma_mld * fT * D_sst) - 1)
            # fluxes
            D_FIN = p_circ * v_fg * alpha_CO2 * D_CO2
            D_FOUT = p_circ * v_fg * alpha_CO2 * f_pCO2(D_dic, fT * D_sst)
            D_FCIRC = D_CSURF * (1 / tau_circ)
            OSNK = np.sum(D_FOUT - D_FIN)
            # stocks
            D_CSURF += dt * (D_FIN - D_FOUT - D_FCIRC)
            D_dic = alpha_dic * np.sum(D_CSURF) / (1 + D_mld / mld_0)

            # --------
            # 2. LAND
            # --------

            # land-cover
            D_AREA += dt * (np.sum(LUC[t], 1) - np.sum(LUC[t], 2))
            D_AWET = AWET_0 * (
                    gamma_wetT * fT * D_lst + gamma_wetP * fT * D_lyp + gamma_wetC * fT * D_CO2)
            # factors
            D_k_igni = (
                    gamma_igniT * fT * D_lst[:, np.newaxis]
                    + gamma_igniP * fT * D_lyp[:, np.newaxis]
                    + gamma_igniC * fT * D_CO2
            )
            D_k_rho = f_rho(fT * D_lst[:, np.newaxis], fT * D_lyp[:, np.newaxis])
            # fluxes
            D_npp = npp_0 * f_npp(D_CO2, fT * D_lst[:, np.newaxis],
                                  fT * D_lyp[:, np.newaxis])
            D_efire = igni_0 * ((1 + D_k_igni) * (cveg_0 + D_cveg) - cveg_0)
            D_fmort = mu_0 * D_cveg
            D_rh1 = rho1_0 * ((1 + D_k_rho) * (csoil1_0 + D_csoil1) - csoil1_0)
            D_fmet = k_met * D_rh1
            D_rh2 = rho2_0 * ((1 + D_k_rho) * (csoil2_0 + D_csoil2) - csoil2_0)
            D_ewet = ewet_0 * np.nan_to_num(
                np.sum(p_wet * D_csoil1, 1) / np.sum(p_wet * csoil1_0, 1))
            LSNK = np.sum((D_rh1 + D_rh2 + D_efire - D_npp) * (AREA_0 + D_AREA))
            D_EWET = ewet_0 * D_AWET + D_ewet * AWET_0 + D_ewet * D_AWET

            # stocks
            D_cveg += dt * (D_npp - D_fmort - D_efire)
            D_csoil1 += dt * (D_fmort - D_fmet - D_rh1)
            D_csoil2 += dt * (D_fmet - D_rh2)

            # ------------
            # 3. LAND-USE
            # ------------

            # initialization
            # land-use change
            for b1 in range(nb_biome):
                for b2 in range(nb_biome):
                    CVEG_luc[:, b1, b2, t] += dt * -(cveg_0 + D_cveg)[:, b2] * LUC[t, :, b1, b2]
                    CSOIL1_luc[:, b1, b2, t] += (
                            dt
                            * ((csoil1_0 + D_csoil1)[:, b1] - (csoil1_0 + D_csoil1)[:, b2])
                            * LUC[t, :, b1, b2]
                    )
                    CSOIL2_luc[:, b1, b2, t] += (
                            dt
                            * ((csoil2_0 + D_csoil2)[:, b1] - (csoil2_0 + D_csoil2)[:, b2])
                            * LUC[t, :, b1, b2]
                    )
                    CSOIL1_luc[:, b1, b2, t] += (
                            dt
                            * (cveg_0 + D_cveg)[:, b1]
                            * (p_AGB[:, b1] * p_HWP0[:, b1] + (1 - p_AGB[:, b1]))
                            * LUC[t, :, b1, b2]
                    )
                    CHWP1_luc[:, b1, b2, t] += (
                            dt * (cveg_0 + D_cveg)[:, b1] * p_AGB[:, b1] * p_HWP1[:, b1] * LUC[t, :, b1, b2]
                    )
                    CHWP2_luc[:, b1, b2, t] += (
                            dt * (cveg_0 + D_cveg)[:, b1] * p_AGB[:, b1] * p_HWP2[:, b1] * LUC[t, :, b1, b2]
                    )
                    CHWP3_luc[:, b1, b2, t] += (
                            dt * (cveg_0 + D_cveg)[:, b1] * p_AGB[:, b1] * p_HWP3[:, b1] * LUC[t, :, b1, b2]
                    )

            # harvest
            for b in range(nb_biome):
                CVEG_luc[:, b, b, t] += dt * -HARV[t, :, b]
                CSOIL1_luc[:, b, b, t] += dt * p_HWP0[:, b] * HARV[t, :, b]
                CHWP1_luc[:, b, b, t] += dt * p_HWP1[:, b] * HARV[t, :, b]
                CHWP2_luc[:, b, b, t] += dt * p_HWP2[:, b] * HARV[t, :, b]
                CHWP3_luc[:, b, b, t] += dt * p_HWP3[:, b] * HARV[t, :, b]

            # shifting cultivation
            for b1 in range(nb_biome):
                for b2 in range(b1, nb_biome):
                    CVEG_luc[:, b1, b2, t] += (
                            dt
                            * -(cveg_0 + D_cveg)[:, b2]
                            * (1 - np.exp(-mu_0[:, b2] * tau_shift))
                            * SHIFT[t, :, b1, b2]
                    )
                    CSOIL1_luc[:, b1, b2, t] += (
                            dt
                            * (cveg_0 + D_cveg)[:, b1]
                            * (1 - np.exp(-mu_0[:, b1] * tau_shift))
                            * (p_AGB[:, b1] * p_HWP0[:, b1] + (1 - p_AGB[:, b1]))
                            * SHIFT[t, :, b1, b2]
                    )
                    CHWP1_luc[:, b1, b2, t] += (
                            dt
                            * (cveg_0 + D_cveg)[:, b1]
                            * (1 - np.exp(-mu_0[:, b1] * tau_shift))
                            * p_AGB[:, b1]
                            * p_HWP1[:, b1]
                            * SHIFT[t, :, b1, b2]
                    )
                    CHWP2_luc[:, b1, b2, t] += (
                            dt
                            * (cveg_0 + D_cveg)[:, b1]
                            * (1 - np.exp(-mu_0[:, b1] * tau_shift))
                            * p_AGB[:, b1]
                            * p_HWP2[:, b1]
                            * SHIFT[t, :, b1, b2]
                    )
                    CHWP3_luc[:, b1, b2, t] += (
                            dt
                            * (cveg_0 + D_cveg)[:, b1]
                            * (1 - np.exp(-mu_0[:, b1] * tau_shift))
                            * p_AGB[:, b1]
                            * p_HWP3[:, b1]
                            * SHIFT[t, :, b1, b2]
                    )
                    CVEG_luc[:, b2, b1, t] += (
                            dt
                            * -(cveg_0 + D_cveg)[:, b1]
                            * (1 - np.exp(-mu_0[:, b1] * tau_shift))
                            * SHIFT[t, :, b1, b2]
                    )
                    CSOIL1_luc[:, b2, b1, t] += (
                            dt
                            * (cveg_0 + D_cveg)[:, b2]
                            * (1 - np.exp(-mu_0[:, b2] * tau_shift))
                            * (p_AGB[:, b2] * p_HWP0[:, b2] + (1 - p_AGB[:, b2]))
                            * SHIFT[t, :, b1, b2]
                    )
                    CHWP1_luc[:, b2, b1, t] += (
                            dt
                            * (cveg_0 + D_cveg)[:, b2]
                            * (1 - np.exp(-mu_0[:, b2] * tau_shift))
                            * p_AGB[:, b2]
                            * p_HWP1[:, b2]
                            * SHIFT[t, :, b1, b2]
                    )
                    CHWP2_luc[:, b2, b1, t] += (
                            dt
                            * (cveg_0 + D_cveg)[:, b2]
                            * (1 - np.exp(-mu_0[:, b2] * tau_shift))
                            * p_AGB[:, b2]
                            * p_HWP2[:, b2]
                            * SHIFT[t, :, b1, b2]
                    )
                    CHWP3_luc[:, b2, b1, t] += (
                            dt
                            * (cveg_0 + D_cveg)[:, b2]
                            * (1 - np.exp(-mu_0[:, b2] * tau_shift))
                            * p_AGB[:, b2]
                            * p_HWP3[:, b2]
                            * SHIFT[t, :, b1, b2]
                    )

            # fluxes
            # book-keeping model
            NPP_luc = 0 * CVEG_luc
            EFIRE_luc = (igni_0 * (1 + D_k_igni))[:, np.newaxis, :, np.newaxis] * CVEG_luc
            FMORT_luc = mu_0[:, np.newaxis, :, np.newaxis] * CVEG_luc
            RH1_luc = (rho1_0 * (1 + D_k_rho))[:, np.newaxis, :, np.newaxis] * CSOIL1_luc
            FMET_luc = k_met * RH1_luc
            RH2_luc = (rho2_0 * (1 + D_k_rho))[:, np.newaxis, :, np.newaxis] * CSOIL2_luc
            EHWP1_luc = np.zeros([nb_regionI, nb_biome, nb_biome, ind_final + 1], dtype=dty)
            EHWP1_luc[:, :, :, : t + 1] = \
                r_HWP1[np.newaxis, np.newaxis, np.newaxis, t::-1] * CHWP1_luc[:, :, :, : t + 1]
            EHWP2_luc = np.zeros([nb_regionI, nb_biome, nb_biome, ind_final + 1], dtype=dty)
            EHWP2_luc[:, :, :, : t + 1] = \
                r_HWP2[np.newaxis, np.newaxis, np.newaxis, t::-1] * CHWP2_luc[:, :, :, : t + 1]
            EHWP3_luc = np.zeros([nb_regionI, nb_biome, nb_biome, ind_final + 1], dtype=dty)
            EHWP3_luc[:, :, :, : t + 1] = \
                r_HWP3[np.newaxis, np.newaxis, np.newaxis, t::-1] * CHWP3_luc[:, :, :, : t + 1]
            ELUC = \
                np.sum(
                    np.sum(np.sum(RH1_luc + RH2_luc + EFIRE_luc + EHWP1_luc + EHWP2_luc + EHWP3_luc - NPP_luc, 3), 2),
                    1)

            # biomass burning
            def biomass_burning_diff(alpha_BB):
                new = np.newaxis
                diff = np.sum(alpha_BB * (igni_0 * cveg_0 * D_AREA + D_efire * AREA_0 + D_efire * D_AREA), 1)
                diff += p_HWP1_bb * np.sum(np.sum(np.sum(alpha_BB[:, :, new, new] * EHWP1_luc, 3), 2), 1)
                diff += np.sum(np.sum(np.sum(alpha_BB[:, new, :, new] * EFIRE_luc, 3), 2), 1)
                return diff

            D_EBB_CO2 = biomass_burning_diff(alpha_BB_CO2)
            D_EBB_CH4 = biomass_burning_diff(alpha_BB_CH4)
            D_EBB_N2O = biomass_burning_diff(alpha_BB_N2O)
            D_EBB_NOX = biomass_burning_diff(alpha_BB_NOX)
            D_EBB_CO = biomass_burning_diff(alpha_BB_CO)
            D_EBB_VOC = biomass_burning_diff(alpha_BB_VOC)
            D_EBB_SO2 = biomass_burning_diff(alpha_BB_SO2)
            D_EBB_NH3 = biomass_burning_diff(alpha_BB_NH3)
            D_EBB_OC = biomass_burning_diff(alpha_BB_OC)
            D_EBB_BC = biomass_burning_diff(alpha_BB_BC)

            # stocks
            CVEG_luc += dt * (NPP_luc - FMORT_luc - EFIRE_luc)
            CSOIL1_luc += dt * (FMORT_luc - FMET_luc - RH1_luc)
            CSOIL2_luc += dt * (FMET_luc - RH2_luc)
            CHWP1_luc += dt * -EHWP1_luc
            CHWP2_luc += dt * -EHWP2_luc
            CHWP3_luc += dt * -EHWP3_luc

            # ---------------
            # 3b. PERMAFROST
            # ---------------

            # factors
            rD_rhoPF = \
                np.exp(w_rhoPF * gamma_rhoPF1 * w_reg_lstPF * fT * D_gst
                       + w_rhoPF * gamma_rhoPF2 * (w_reg_lstPF * fT * D_gst) ** 2) - 1

            # fraction thawed
            pthaw_bar = \
                -pthaw_min + (1 + pthaw_min) / (1 + ((1 / pthaw_min + 1) ** k_pthaw - 1)
                                                * np.exp(-gamma_pthaw * k_pthaw * w_reg_lstPF * fT * D_gst)) ** (
                        1 / k_pthaw)
            d_pthaw = f_v_PF(pthaw_bar, pthaw) * (pthaw_bar - pthaw)
            pthaw += dt * d_pthaw

            # fluxes
            FTHAW = CFROZ_0 * d_pthaw
            ETHAW1 = 1 / tau_PF1 * (1 + rD_rhoPF) * CTHAW1
            ETHAW2 = 1 / tau_PF2 * (1 + rD_rhoPF) * CTHAW2
            ETHAW3 = 1 / tau_PF3 * (1 + rD_rhoPF) * CTHAW3
            EPF_CO2 = (1 - p_PF_CH4) * (ETHAW1 + ETHAW2 + ETHAW3 + p_PF_inst * FTHAW)
            EPF_CH4 = 1000.0 * p_PF_CH4 * (
                    ETHAW1 + ETHAW2 + ETHAW3 + p_PF_inst * FTHAW)
            EPF = EPF_CO2 + 0.001 * EPF_CH4

            # stocks
            D_CFROZ -= dt * FTHAW
            CTHAW1 += dt * (p_PF1 * (1 - p_PF_inst) * FTHAW - ETHAW1)
            CTHAW2 += dt * (p_PF2 * (1 - p_PF_inst) * FTHAW - ETHAW2)
            CTHAW3 += dt * (p_PF3 * (1 - p_PF_inst) * FTHAW - ETHAW3)

            # -------------
            # 4. CHEMISTRY
            # -------------

            # factors
            D_kOH = f_kOH(
                D_CH4,
                D_O3s,
                fT * D_gst,
                np.sum(ENOX[t] + D_EBB_NOX),
                np.sum(ECO[t] + D_EBB_CO),
                np.sum(EVOC[t] + D_EBB_VOC),
            )
            D_hv = f_hv(D_N2O_lag, D_EESC, fT * D_gst)

            # fluxes
            D_OHSNK_CH4 = -alpha_CH4 / tau_CH4_OH * (CH4_0 * D_kOH + D_CH4 + D_kOH * D_CH4)
            D_HVSNK_CH4 = -alpha_CH4 / tau_CH4_hv * (CH4_0 * D_hv + D_CH4_lag + D_hv * D_CH4_lag)
            D_XSNK_CH4 = -alpha_CH4 * (1 / tau_CH4_soil + 1 / tau_CH4_ocean) * D_CH4
            D_FOXI_CH4 = -0.001 * (1.0 * np.sum(ECH4[t]) + np.sum(D_EBB_CH4) + np.sum(D_EWET) + D_OHSNK_CH4 + D_HVSNK_CH4 + D_XSNK_CH4)
            D_HVSNK_N2O = -alpha_N2O / tau_N2O_hv * (N2O_0 * D_hv + D_N2O_lag + D_hv * D_N2O_lag)
            D_OHSNK_HFC = -alpha_HFC / tau_HFC_OH * (HFC_0 * D_kOH + D_HFC + D_kOH * D_HFC)
            D_OHSNK_PFC = -alpha_PFC / tau_PFC_OH * (PFC_0 * D_kOH + D_PFC + D_kOH * D_PFC)
            D_OHSNK_ODS = -alpha_ODS / tau_ODS_OH * (ODS_0 * D_kOH + D_ODS + D_kOH * D_ODS)

            D_HVSNK_HFC = -alpha_HFC / tau_HFC_hv * (HFC_0 * D_hv + D_HFC_lag + D_hv * D_HFC_lag)
            D_HVSNK_PFC = -alpha_PFC / tau_PFC_hv * (PFC_0 * D_hv + D_PFC_lag + D_hv * D_PFC_lag)
            D_HVSNK_ODS = -alpha_ODS / tau_ODS_hv * (ODS_0 * D_hv + D_ODS_lag + D_hv * D_ODS_lag)

            D_XSNK_HFC = -alpha_HFC / tau_HFC_othr * D_HFC
            D_XSNK_PFC = -alpha_PFC / tau_PFC_othr * D_PFC
            D_XSNK_ODS = -alpha_ODS / tau_ODS_othr * D_ODS

            # stocks
            D_O3t = chi_O3t_CH4 * np.log(1 + D_CH4 / CH4_0) + Gamma_O3t * fT * D_gst
            D_O3t += chi_O3t_NOX * np.sum(w_reg_NOX * np.sum(p_reg4 * (ENOX[t] + D_EBB_NOX)[:, np.newaxis], 0))
            D_O3t += chi_O3t_CO * np.sum(w_reg_CO * np.sum(p_reg4 * (ECO[t] + D_EBB_CO)[:, np.newaxis], 0))
            D_O3t += chi_O3t_VOC * np.sum(w_reg_VOC * np.sum(p_reg4 * (EVOC[t] + D_EBB_VOC)[:, np.newaxis], 0))
            D_EESC = np.sum(f_fracrel(tau_lag) * (n_Cl + alpha_Br * n_Br) * D_ODS_lag)
            D_O3s = chi_O3s_EESC * D_EESC + chi_O3s_N2O * D_N2O_lag \
                    * (1 - D_EESC / EESC_x) + Gamma_O3s * fT * D_gst
            D_SO4 = (alpha_SO4 * tau_SO2 * np.sum(
                w_reg_SO2 * np.sum(p_reg4 * (ESO2[t] + D_EBB_SO2)[:, np.newaxis], 0))
                     + alpha_SO4 * tau_DMS * 0
                     + Gamma_SO4 * fT * D_gst
                     )
            D_POA = (
                    tau_OMff * alpha_POM * np.sum(
                w_reg_OC * np.sum(p_reg4 * (EOC[t])[:, np.newaxis], 0))
                    + tau_OMbb * alpha_POM * np.sum(D_EBB_OC)
                    + Gamma_POA * fT * D_gst
            )
            D_BC = (
                    tau_BCff * np.sum(
                w_reg_BC * np.sum(p_reg4 * (EBC[t])[:, np.newaxis], 0))
                    + tau_BCbb * np.sum(D_EBB_BC)
                    + Gamma_BC * fT * D_gst
            )
            D_NO3 = (
                    alpha_NO3 * tau_NOX * np.sum(ENOX[t] + D_EBB_NOX)
                    + alpha_NO3 * tau_NH3 * np.sum(ENH3[t] + D_EBB_NH3)
                    + Gamma_NO3 * fT * D_gst
            )
            D_SOA = tau_VOC * np.sum(EVOC[t] + D_EBB_VOC) + tau_BVOC * 0 + Gamma_SOA * fT * D_gst
            D_DUST = 0 * (tau_DUST * 0 + Gamma_DUST * fT * D_gst)
            D_SALT = 0 * (tau_SALT * 0 + Gamma_SALT * fT * D_gst)
            D_AERh = (
                    solub_SO4 * D_SO4
                    + solub_POA * D_POA
                    + solub_BC * D_BC
                    + solub_NO3 * D_NO3
                    + solub_SOA * D_SOA
                    + solub_DUST * D_DUST
                    + solub_SALT * D_SALT
            )

            # --------------
            # 5. ATMOSPHERE
            # --------------

            # stocks
            D_CO2 += (
                    dt
                    * (1 / alpha_CO2)
                    * (np.sum(EFF[t]) + np.sum(
                ELUC) + LSNK + OSNK + D_FOXI_CH4 + np.sum(EPF_CO2))
            )
            D_CH4 += (
                    dt
                    * (1 / alpha_CH4)
                    * (
                            np.sum(ECH4[t])
                            + np.sum(D_EBB_CH4)
                            + np.sum(D_EWET)
                            + np.sum(EPF_CH4)
                            + D_OHSNK_CH4
                            + D_HVSNK_CH4
                            + D_XSNK_CH4
                    )
            )
            D_N2O += dt * (1 / alpha_N2O) * (np.sum(EN2O[t]) + np.sum(D_EBB_N2O) + D_HVSNK_N2O)
            D_HFC += dt * (1 / alpha_HFC) * (np.sum(EHFC[t], 0) + D_OHSNK_HFC + D_HVSNK_HFC + D_XSNK_HFC)
            D_PFC += dt * (1 / alpha_PFC) * (np.sum(EPFC[t], 0) + D_OHSNK_PFC + D_HVSNK_PFC + D_XSNK_PFC)
            D_ODS += dt * (1 / alpha_ODS) * (np.sum(EODS[t], 0) + D_OHSNK_ODS + D_HVSNK_ODS + D_XSNK_ODS)

            D_CH4_lag += dt * ((1 / tau_lag) * D_CH4 - (1 / tau_lag) * D_CH4_lag)
            D_N2O_lag += dt * ((1 / tau_lag) * D_N2O - (1 / tau_lag) * D_N2O_lag)
            D_HFC_lag += dt * ((1 / tau_lag) * D_HFC - (1 / tau_lag) * D_HFC_lag)
            D_PFC_lag += dt * ((1 / tau_lag) * D_PFC - (1 / tau_lag) * D_PFC_lag)
            D_ODS_lag += dt * ((1 / tau_lag) * D_ODS - (1 / tau_lag) * D_ODS_lag)

            # FORCE
            if force_CO2:
                D_CO2 = D_CO2_force[t]

            if force_GHG:
                D_CO2 = D_CO2_force[t]
                D_CH4 = D_CH4_force[t]
                D_N2O = D_N2O_force[t]

            if force_halo:
                D_HFC[:] = D_HFC_force[t]
                D_PFC[:] = D_PFC_force[t]
                D_ODS[:] = D_ODS_force[t]

            # -----------
            # 6. CLIMATE
            # -----------

            # fluxes
            # per component
            RF_CO2 = f_RF_CO2(D_CO2)
            RF_CH4 = f_RF_CH4(D_CH4) - (
                    f_RF_overlap(D_CH4, D_N2O) - f_RF_overlap(0, D_N2O))
            RF_H2Os = f_RF_H2Os(D_CH4_lag)
            RF_N2O = f_RF_N2O(D_N2O) - (
                    f_RF_overlap(D_CH4, D_N2O) - f_RF_overlap(D_CH4, 0))
            RF_halo = np.sum(radeff_HFC * D_HFC) + np.sum(
                radeff_PFC * D_PFC) + np.sum(radeff_ODS * D_ODS)

            RF_O3t = radeff_O3t * D_O3t
            RF_O3s = radeff_O3s * D_O3s
            RF_SO4 = radeff_SO4 * D_SO4
            RF_POA = radeff_POA * D_POA
            RF_BC = radeff_BC * D_BC
            RF_NO3 = radeff_NO3 * D_NO3
            RF_SOA = radeff_SOA * D_SOA
            RF_DUST = radeff_DUST * D_DUST
            RF_SALT = radeff_SALT * D_SALT
            RF_cloud = k_BC_adjust * RF_BC + Phi_0 * np.log(1 + D_AERh / AERh_0)
            RF_BCsnow = radeff_BCsnow * np.sum(
                w_reg_BCsnow * np.sum(p_reg9 * (EBC[t] + D_EBB_BC)[:, np.newaxis], 0)
            )
            RF_LCC = np.sum(alpha_LCC * D_AREA)

            # FORCE
            if force_RFs:
                RF_CO2 = RF_CO2_force[t]
                RF_CH4 = RF_CH4_force[t]
                RF_H2Os = RF_H2Os_force[t]
                RF_N2O = RF_N2O_force[t]
                RF_halo = RF_halo_force[t]
                RF_O3t = RF_O3t_force[t]
                RF_O3s = RF_O3s_force[t]
                RF_SO4 = RF_SO4_force[t]
                RF_POA = RF_POA_force[t]
                RF_BC = RF_BC_force[t]
                RF_NO3 = RF_NO3_force[t]
                RF_SOA = RF_SOA_force[t]
                RF_DUST = RF_DUST_force[t]
                RF_SALT = RF_SALT_force[t]
                RF_cloud = RF_cloud_force[t]
                RF_BCsnow = RF_BCsnow_force[t]
                RF_LCC = RF_LCC_force[t]

            # totals
            RF = (
                    RF_CO2
                    + RF_CH4
                    + RF_H2Os
                    + RF_N2O
                    + RF_halo
                    + RF_O3t
                    + RF_O3s
                    + RF_SO4
                    + RF_POA
                    + RF_BC
                    + RF_NO3
                    + RF_SOA
                    + RF_DUST
                    + RF_SALT
                    + RF_cloud
                    + RF_BCsnow
                    + RF_LCC
                    + RFcon[t]
                    + RFvolc[t]
                    + RFsolar[t]
            )
            RF_warm = (
                    RF_CO2
                    + RF_CH4
                    + RF_H2Os
                    + RF_N2O
                    + RF_halo
                    + RF_O3t
                    + RF_O3s
                    + RF_SO4
                    + RF_POA
                    + RF_BC
                    + RF_NO3
                    + RF_SOA
                    + RF_DUST
                    + RF_SALT
                    + RF_cloud
                    + warmeff_BCsnow * RF_BCsnow
                    + warmeff_LCC * RF_LCC
                    + RFcon[t]
                    + warmeff_volc * RFvolc[t]
                    + RFsolar[t]
            )
            RF_atm = (
                    p_atm_CO2 * RF_CO2
                    + p_atm_noCO2 * (RF_CH4 + RF_N2O + RF_halo)
                    + p_atm_O3t * RF_O3t
                    + p_atm_strat * (RF_O3s + RF_H2Os)
                    + p_atm_scatter * (RF_SO4 + RF_POA + RF_NO3 + RF_SOA + RF_DUST + RF_SALT + RFvolc[t])
                    + p_atm_absorb * RF_BC
                    + p_atm_cloud * (RF_cloud + RFcon[t])
                    + p_atm_alb * (RF_BCsnow + RF_LCC)
                    + p_atm_solar * RFsolar[t]
            )

            # FORCE
            if force_RF:
                RF_warm = RF_force[t] * (RF_warm / RF)
                RF_atm = RF_force[t] * (RF_atm / RF)
                RF = RF_force[t]

            # stocks
            # temperatures
            D_gst += dt * (1 / tau_gst) * (
                    lambda_0 * RF_warm - D_gst - theta_0 * (D_gst - D_gst0))
            D_gst0 += dt * (1 / tau_gst0) * theta_0 * (D_gst - D_gst0)
            D_sst = w_reg_sst * D_gst
            D_lst = w_reg_lst * D_gst
            # precipitations
            D_gyp = alpha_gyp * D_gst + beta_gyp * RF_atm
            D_lyp = w_reg_lyp * D_gyp
            # ocean
            D_OHC += dt * p_OHC * alpha_OHC * (RF - D_gst / lambda_0)
            D_pH = f_pH(D_CO2)

            # FORCE
            if force_clim:
                D_gst = D_gst_force[t]
                D_sst = D_sst_force[t]
                D_lst = D_lst_force[t]
                D_lyp = D_lyp_force[t]

            # -----------
            # Y. SAVE
            # -----------

            loc = locals()
            values_run = [loc[name] for name in var_timeseries]
            for arr, run in zip(values_timeseries, values_run):
                arr[t] += dt * run
            del loc

            # ---------
            # Z. TESTS
            # ---------

            if np.isnan(np.sum(D_CO2)):
                print("D_CO2 = NaN at t = " + str(t) + " and tt = " + str(tt))
                print("OSNK = " + str(np.sum(OSNK)))
                print("LSNK = " + str(np.sum(LSNK)))
                print("ELUC = " + str(np.sum(ELUC)))
                break
            if np.isnan(np.sum(D_CH4)):
                print("D_CH4 = NaN at t = " + str(t) + " and tt = " + str(tt))
                print("D_EWET = " + str(np.sum(D_EWET)))
                print("D_OHSNK = " + str(np.sum(D_OHSNK_CH4)))
                print("D_HVSNK = " + str(np.sum(D_HVSNK_CH4)))
                break
            if np.isnan(np.sum(D_gst)):
                print("D_gst = NaN at t = " + str(t) + " and tt = " + str(tt))
                print("RF_CO2 = " + str(np.sum(RF_CO2)))
                print("RF_CH4 = " + str(np.sum(RF_CH4)))
                print("RF_H2Os = " + str(np.sum(RF_H2Os)))
                print("RF_N2O = " + str(np.sum(RF_N2O)))
                print("RF_halo = " + str(np.sum(RF_halo)))
                print("RF_O3t = " + str(np.sum(RF_O3t)))
                print("RF_O3s = " + str(np.sum(RF_O3s)))
                print("RF_SO4 = " + str(np.sum(RF_SO4)))
                print("RF_POA = " + str(np.sum(RF_POA)))
                print("RF_BC = " + str(np.sum(RF_BC)))
                print("RF_NO3 = " + str(np.sum(RF_NO3)))
                print("RF_SOA = " + str(np.sum(RF_SOA)))
                print("RF_DUST = " + str(np.sum(RF_DUST)))
                print("RF_SALT = " + str(np.sum(RF_SALT)))
                print("RF_cloud = " + str(np.sum(RF_cloud)))
                print("RF_BCsnow = " + str(np.sum(RF_BCsnow)))
                print("RF_LCC = " + str(np.sum(RF_LCC)))
                break

        if np.isnan(np.sum(D_CO2)) | np.isnan(np.sum(D_CH4)) | np.isnan(np.sum(D_gst)):
            for arr in values_timeseries:
                if t < ind_final:
                    arr[t + 1:] = np.nan
            break

    # ===========
    # C. FIGURES
    # ===========

    print("PLOTTING", plot)
    from ..plot import plot_AER, plot_CH4, plot_clim, plot_CO2, plot_N2O, plot_O3

    if plot == "all" or plot == "CO2" or "CO2" in plot:
        plot_CO2(
            D_CO2_t,
            OSNK_t,
            LSNK_t,
            ELUC_t,
            EFF,
            D_AREA_t,
            D_npp_t,
            D_efire_t,
            D_fmort_t,
            D_rh1_t,
            D_fmet_t,
            D_rh2_t,
            D_FIN_t,
            D_FOUT_t,
            D_FCIRC_t,
            EFIRE_luc_t,
            FMORT_luc_t,
            RH1_luc_t,
            RH2_luc_t,
            EHWP1_luc_t,
            EHWP2_luc_t,
            EHWP3_luc_t,
        )
        plt.savefig("results/plot-CO2.svg")
    if plot == "all" or plot == "CH4" or "CH4" in plot:
        plot_CH4(D_CH4_t, D_OHSNK_CH4_t, D_HVSNK_CH4_t, D_XSNK_CH4_t, D_EWET_t,
                 D_EBB_CH4_t, ECH4)
        plt.savefig("results/plot-CH4.svg")
    if plot == "all" or plot == "N2O" or "N2O" in plot:
        plot_N2O(D_N2O_t, D_HVSNK_N2O_t, D_EBB_N2O_t, EN2O)
        plt.savefig("results/plot-N20.svg")
    if plot == "all" or plot == "O3" or "O3" in plot:
        plot_O3(D_O3t_t, D_O3s_t, D_EESC_t, D_N2O_lag_t, D_gst_t)
        plt.savefig("results/plot-03.svg")
    if plot == "all" or plot == "AER" or "AER" in plot:
        plot_AER(
            D_SO4_t,
            D_POA_t,
            D_BC_t,
            D_NO3_t,
            D_SOA_t,
            D_AERh_t,
            RF_SO4_t,
            RF_POA_t,
            RF_BC_t,
            RF_NO3_t,
            RF_SOA_t,
            RF_cloud_t,
        )
        plt.savefig("results/plot-AER.svg")
    if plot == "all" or plot == "clim" or "clim" in plot:
        plot_clim(
            RF_t,
            D_gst_t,
            D_gyp_t,
            RF_CO2_t,
            RF_CH4_t,
            RF_H2Os_t,
            RF_N2O_t,
            RF_halo_t,
            RF_O3t_t,
            RF_O3s_t,
            RF_SO4_t,
            RF_POA_t,
            RF_BC_t,
            RF_NO3_t,
            RF_SOA_t,
            RF_cloud_t,
            RF_BCsnow_t,
            RF_LCC_t,
            RFcon,
            RFvolc,
            RFsolar,
        )
        plt.savefig("results/plot-clim.svg")

    # ===========
    # D. OUTPUTS
    # ===========

    output = []
    loc = locals()
    for var in var_output:
        output.append(loc[var + '_t'])
    return output


#
# REMOVE ATTRIBUTION AXIS
#

# Drivers
EFF = reduce_driver(EFF)
ECH4 = reduce_driver(ECH4)
EN2O = reduce_driver(EN2O)
LUC = reduce_driver(LUC)
HARV = reduce_driver(HARV)
SHIFT = reduce_driver(SHIFT)
EHFC = reduce_driver(EHFC)
EPFC = reduce_driver(EPFC)
EODS = reduce_driver(EODS)
ENOX = reduce_driver(ENOX)
ECO = reduce_driver(ECO)
EVOC = reduce_driver(EVOC)
ESO2 = reduce_driver(ESO2)
ENH3 = reduce_driver(ENH3)
EOC = reduce_driver(EOC)
EBC = reduce_driver(EBC)
RFcon = reduce_driver(RFcon)
RFvolc = reduce_driver(RFvolc)
RFsolar = reduce_driver(RFsolar)

# parameters
ECH4_0 = reduce_param(ECH4_0)
EN2O_0 = reduce_param(EN2O_0)
ENOX_0 = reduce_param(ENOX_0)
ECO_0 = reduce_param(ECO_0)
EVOC_0 = reduce_param(EVOC_0)
ESO2_0 = reduce_param(ESO2_0)
ENH3_0 = reduce_param(ENH3_0)
EOC_0 = reduce_param(EOC_0)
EBC_0 = reduce_param(EBC_0)
