import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties

from ..data_loaders import ind_cdiac
from ..params import alpha_CO2, AREA_0, alpha_CH4, tau_CH4_OH, alpha_N2O, tau_N2O_hv, gamma_age, tau_lag, tau_CH4_hv, tau_CH4_soil, tau_CH4_ocean, chi_O3s_EESC, EESC_x, gst_giss, gst_had, gst_ncdc, chi_O3s_N2O
from .. import historical
from .. import scenarios
from .. import config


def plot_CO2(D_CO2, OSNK, LSNK, ELUC, EFF, D_AREA, D_npp, D_efire, D_fmort, D_rh1, D_fmet, D_rh2, D_FIN, D_FOUT, D_FCIRC, D_MORT_luc, D_EFIRE_luc, D_RH1_luc, D_RH2_luc, EHWP1_luc, EHWP2_luc, EHWP3_luc,):
    plt.figure()

    # atmospheric CO2
    ax = plt.subplot(2, 3, 1)
    plt.plot(1700 + np.arange(config.ind_final + 1), D_CO2, color="k", lw=2, label="OSCAR")
    plt.plot(1700 + np.arange(len(historical.CO2_ipcc)), historical.CO2_ipcc - historical.CO2_0, color="r", lw=2, ls="--", label="IPCC")
    if config.ind_final > ind_cdiac:
        plt.plot(1700 + np.arange(min(len(scenarios.CO2_rcp), config.ind_final + 1)), scenarios.CO2_rcp[: min(len(scenarios.CO2_rcp), config.ind_final + 1), 0] - historical.CO2_0, color="0.8", lw=2, ls=":", label="RCP2.6", )
        plt.plot(1700 + np.arange(min(len(scenarios.CO2_rcp), config.ind_final + 1)), scenarios.CO2_rcp[: min(len(scenarios.CO2_rcp), config.ind_final + 1), 1] - historical.CO2_0, color="0.6", lw=2, ls=":", label="RCP4.5", )
        plt.plot(1700 + np.arange(min(len(scenarios.CO2_rcp), config.ind_final + 1)), scenarios.CO2_rcp[: min(len(scenarios.CO2_rcp), config.ind_final + 1), 2] - historical.CO2_0, color="0.4", lw=2, ls=":", label="RCP6.0", )
        plt.plot(1700 + np.arange(min(len(scenarios.CO2_rcp), config.ind_final + 1)), scenarios.CO2_rcp[: min(len(scenarios.CO2_rcp), config.ind_final + 1), 3] - historical.CO2_0, color="0.2", lw=2, ls=":", label="RCP8.5", )
    plt.title(r"$\Delta$CO2 (ppm)$", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + config.ind_final])

    # budget fluxes
    ax = plt.subplot(2, 3, 2)
    plt.plot([1700, 1700 + config.ind_final + 1], [0, 0], "k-")
    plt.plot(1700 + np.arange(config.ind_final + 1), np.sum(EFF, 1), color="#666666", lw=2, label="EFF")
    plt.plot(1700 + np.arange(config.ind_final + 1), np.sum(ELUC, 1), color="#993300", lw=2, label="ELUC")
    plt.plot(1700 + np.arange(config.ind_final + 1), OSNK, color="#000099", lw=2, label="OSNK")
    plt.plot(1700 + np.arange(config.ind_final + 1), LSNK, color="#009900", lw=2, label="LSNK")
    plt.plot(1700 + np.arange(config.ind_final) + 1, alpha_CO2 * (D_CO2[1:] - D_CO2[:-1]), color="#FFCC00", lw=2, label="d_CO2")
    plt.plot(1700 + np.arange(len(historical.EFF_gcp)), historical.EFF_gcp, color="#666666", ls="--")
    plt.plot(1700 + np.arange(len(historical.ELUC_gcp)), historical.ELUC_gcp, color="#CC3300", ls="--")
    plt.plot(1700 + np.arange(len(historical.OSNK_gcp)), historical.OSNK_gcp, color="#000099", ls="--")
    plt.plot(1700 + np.arange(len(historical.LSNK_gcp)), historical.LSNK_gcp, color="#009900", ls="--")
    plt.plot(1700 + np.arange(len(historical.d_CO2_gcp)), historical.d_CO2_gcp, color="#FFCC00", ls="--")
    plt.plot([1700, 1700], [0, 0], "k--", label="GCP")
    plt.title("CO2 fluxes (GtC/yr)", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + config.ind_final])

    # airborne fraction
    ax = plt.subplot(2, 3, 3)
    plt.plot([1700, 1700 + config.ind_final + 1], [0, 0], "k-")
    plt.plot(1700 + np.arange(config.ind_final) + 1, alpha_CO2 * (D_CO2[1:] - D_CO2[:-1]) / np.sum(EFF + ELUC, 1)[1:], color="#FFCC00", lw=1, label="AF", )
    plt.plot(1700 + np.arange(config.ind_final + 1), -OSNK / np.sum(EFF + ELUC, 1), color="#000099", lw=1, label="OF")
    plt.plot(1700 + np.arange(config.ind_final + 1), -LSNK / np.sum(EFF + ELUC, 1), color="#009900", lw=1, label="LF")
    plt.plot(np.arange(1959, 1700 + ind_cdiac + 1), np.ones([ind_cdiac - 259 + 1])
        * np.mean((alpha_CO2 * (D_CO2[1:] - D_CO2[:-1]) / np.sum(EFF + ELUC, 1)[1:])[
                  259 - 1: ind_cdiac]), color="k", lw=2, label="OSCAR", )
    plt.plot(np.arange(1959, 1700 + ind_cdiac + 1), np.ones([ind_cdiac - 259 + 1]) * np.mean((historical.d_CO2_gcp / (historical.EFF_gcp + historical.ELUC_gcp))[259: ind_cdiac + 1]), color="r", lw=2, ls="--", label="GCP", )
    plt.title("airborne fraction", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + config.ind_final])
    ax.set_ylim([-0.2, 1.2])

    # ELUC details
    ax = plt.subplot(2, 3, 4)
    plt.plot([1700, 1700 + config.ind_final + 1], [0, 0], "k-")
    plt.plot(1700 + np.arange(config.ind_final + 1), np.sum(ELUC, 1), color="k", ls="-.", lw=2, label="ELUC")
    plt.plot(1700 + np.arange(config.ind_final + 1), np.sum(np.sum(np.sum(np.sum(D_EFIRE_luc + D_RH1_luc + D_RH2_luc, 4), 3), 2), 1), color="#009900", lw=2, label="ELUC_bio", )
    plt.plot(1700 + np.arange(config.ind_final + 1), np.sum(np.sum(np.sum(np.sum(EHWP1_luc + EHWP2_luc + EHWP3_luc, 4), 3), 2), 1), color="#993300", lw=2, label="ELUC_hwp", )
    plt.plot(1700 + np.arange(config.ind_final + 1), np.sum(np.sum(np.sum(np.sum(EHWP1_luc, 4), 3), 2), 1), color="#FF3300", lw=1, label="EHWP1", )
    plt.plot(1700 + np.arange(config.ind_final + 1), np.sum(np.sum(np.sum(np.sum(EHWP2_luc, 4), 3), 2), 1), color="#CC9900", lw=1, label="EHWP2", )
    plt.plot(1700 + np.arange(config.ind_final + 1), np.sum(np.sum(np.sum(np.sum(EHWP3_luc, 4), 3), 2), 1), color="#663300", lw=1, label="EHWP3", )
    plt.title("ELUC fluxes (GtC/yr)", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))

    # LSNK details
    ax = plt.subplot(2, 3, 5)
    plt.plot([1700, 1700 + config.ind_final + 1], [0, 0], "k-")
    plt.plot(1700 + np.arange(config.ind_final + 1), -LSNK, color="k", lw=2, ls="-.", label="$-$LSNK")
    plt.plot(1700 + np.arange(config.ind_final + 1), np.sum(np.sum(D_npp * (AREA_0 + D_AREA), 2), 1), color="#009900", lw=2, label="D_NPP", )
    plt.plot(1700 + np.arange(config.ind_final + 1), np.sum(np.sum(D_efire * (AREA_0 + D_AREA), 2), 1), color="#FF3300", lw=2, label="D_EFIRE", )
    plt.plot(1700 + np.arange(config.ind_final + 1), np.sum(np.sum(D_fmort * (AREA_0 + D_AREA), 2), 1), color="#336633", lw=2, label="D_FMORT", )
    plt.plot(1700 + np.arange(config.ind_final + 1), np.sum(np.sum((D_rh1 + D_rh2) * (AREA_0 + D_AREA), 2), 1), color="#663300", lw=2, label="D_RH", )
    plt.title("LSNK fluxes (GtC/yr)", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))

    # OSNK details
    ax = plt.subplot(2, 3, 6)
    plt.plot([1700, 1700 + config.ind_final + 1], [0, 0], "k-")
    plt.plot(1700 + np.arange(config.ind_final + 1), -OSNK, color="k", lw=2, ls="-.", label="$-$OSNK")
    plt.plot(1700 + np.arange(config.ind_final + 1), np.sum(D_FIN, 1), color="#000099", lw=2, label="D_FIN")
    plt.plot(1700 + np.arange(config.ind_final + 1), np.sum(D_FOUT, 1), color="#0099FF", lw=2, label="D_FOUT")
    plt.plot(1700 + np.arange(config.ind_final + 1), np.sum(D_FCIRC, 1), color="#663399", lw=2, label="D_FCIRC")
    plt.title("OSNK fluxes (GtC/yr)", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))


# =========
# 2.2. CH4
# =========


def plot_CH4(D_CH4, D_OHSNK_CH4, D_HVSNK_CH4, D_XSNK_CH4, D_EWET, D_EBB_CH4, ECH4):
    plt.figure()

    # atmospheric CH4
    ax = plt.subplot(2, 3, 1)
    plt.plot(1700 + np.arange(config.ind_final + 1), D_CH4, color="k", lw=2, label="OSCAR")
    plt.plot(1700 + np.arange(len(historical.CH4_ipcc)), historical.CH4_ipcc - historical.CH4_0, color="r", lw=2, ls="--", label="IPCC")
    if config.ind_final > ind_cdiac:
        plt.plot(1700 + np.arange(min(len(scenarios.CH4_rcp), config.ind_final + 1)), scenarios.CH4_rcp[: min(len(scenarios.CH4_rcp), config.ind_final + 1), 0] - historical.CH4_0, color="0.8", lw=2, ls=":", label="RCP2.6", )
        plt.plot(1700 + np.arange(min(len(scenarios.CH4_rcp), config.ind_final + 1)), scenarios.CH4_rcp[: min(len(scenarios.CH4_rcp), config.ind_final + 1), 1] - historical.CH4_0, color="0.6", lw=2, ls=":", label="RCP4.5", )
        plt.plot(1700 + np.arange(min(len(scenarios.CH4_rcp), config.ind_final + 1)), scenarios.CH4_rcp[: min(len(scenarios.CH4_rcp), config.ind_final + 1), 2] - historical.CH4_0, color="0.4", lw=2, ls=":", label="RCP6.0", )
        plt.plot(1700 + np.arange(min(len(scenarios.CH4_rcp), config.ind_final + 1)), scenarios.CH4_rcp[: min(len(scenarios.CH4_rcp), config.ind_final + 1), 3] - historical.CH4_0, color="0.2", lw=2, ls=":", label="RCP8.5", )
    plt.title("$\Delta$CH4 (ppb)", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + config.ind_final])

    # budget fluxes
    ax = plt.subplot(2, 3, 2)
    plt.plot([1700, 1700 + config.ind_final + 1], [0, 0], "k-")
    plt.plot(1700 + np.arange(config.ind_final + 1), np.sum(ECH4, 1), color="#666666", lw=2, label="ECH4")
    plt.plot(1700 + np.arange(config.ind_final + 1), np.sum(D_EBB_CH4, 1), color="#993300", lw=2, label="D_EBB")
    plt.plot(1700 + np.arange(config.ind_final + 1), np.sum(D_EWET, 1), color="#006666", lw=2, label="D_EWET")
    plt.plot(1700 + np.arange(config.ind_final + 1), (D_OHSNK_CH4 + D_HVSNK_CH4 + D_XSNK_CH4), color="#990066", lw=2, label="D_SNK")
    plt.plot(1700 + np.arange(config.ind_final) + 1, alpha_CH4 * (D_CH4[1:] - D_CH4[:-1]), color="#FFCC00", lw=2, label="d_CH4")
    plt.title("CH4 fluxes (MtC/yr)", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + config.ind_final])

    # lifetime
    ax = plt.subplot(2, 3, 3)
    plt.plot(1700 + np.arange(config.ind_final + 1), alpha_CH4
        * (historical.CH4_0 + D_CH4)
        / (alpha_CH4 * historical.CH4_0 * (1 / tau_CH4_OH + 1 / tau_CH4_hv + 1 / tau_CH4_soil + 1 / tau_CH4_ocean)
                - D_OHSNK_CH4
                - D_HVSNK_CH4
                - D_XSNK_CH4), color="k", lw=2, label="OSCAR", )
    plt.plot(1700 + np.arange(config.ind_final + 1), alpha_CH4 * (historical.CH4_0 + D_CH4) / (alpha_CH4 * historical.CH4_0 / tau_CH4_OH - D_OHSNK_CH4), color="k", lw=1, label="OH only", )
    plt.title("CH4 lifetime (yr)", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + config.ind_final])

    # wetlands
    ax = plt.subplot(2, 3, 4)
    plt.title("wetlands", fontsize="medium")

    # biomass burning
    ax = plt.subplot(2, 3, 5)
    plt.title("biomass burning", fontsize="medium")


# =========
# 2.3. N2O
# =========


def plot_N2O(D_N2O, D_HVSNK_N2O, D_EBB_N2O, EN2O):
    plt.figure()

    # atmospheric N2O
    ax = plt.subplot(2, 3, 1)
    plt.plot(1700 + np.arange(config.ind_final + 1), D_N2O, color="k", lw=2, label="OSCAR")
    plt.plot(1700 + np.arange(len(historical.N2O_ipcc)), historical.N2O_ipcc - historical.N2O_0, color="r", lw=2, ls="--", label="IPCC")
    if config.ind_final > ind_cdiac:
        plt.plot(1700 + np.arange(min(len(scenarios.N2O_rcp), config.ind_final + 1)), scenarios.N2O_rcp[: min(len(scenarios.N2O_rcp), config.ind_final + 1), 0] - historical.N2O_0, color="0.8", lw=2, ls=":", label="RCP2.6", )
        plt.plot(1700 + np.arange(min(len(scenarios.N2O_rcp), config.ind_final + 1)), scenarios.N2O_rcp[: min(len(scenarios.N2O_rcp), config.ind_final + 1), 1] - historical.N2O_0, color="0.6", lw=2, ls=":", label="RCP4.5", )
        plt.plot(1700 + np.arange(min(len(scenarios.N2O_rcp), config.ind_final + 1)), scenarios.N2O_rcp[: min(len(scenarios.N2O_rcp), config.ind_final + 1), 2] - historical.N2O_0, color="0.4", lw=2, ls=":", label="RCP6.0", )
        plt.plot(1700 + np.arange(min(len(scenarios.N2O_rcp), config.ind_final + 1)), scenarios.N2O_rcp[: min(len(scenarios.N2O_rcp), config.ind_final + 1), 3] - historical.N2O_0, color="0.2", lw=2, ls=":", label="RCP8.5", )
    plt.title("$\Delta$N2O (ppb)", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + config.ind_final])

    # budget fluxes
    ax = plt.subplot(2, 3, 2)
    plt.plot([1700, 1700 + config.ind_final + 1], [0, 0], "k-")
    plt.plot(1700 + np.arange(config.ind_final + 1), np.sum(EN2O, 1), color="#666666", lw=2, label="EN2O")
    plt.plot(1700 + np.arange(config.ind_final + 1), np.sum(D_EBB_N2O, 1), color="#993300", lw=2, label="D_EBB")
    plt.plot(1700 + np.arange(config.ind_final + 1), D_HVSNK_N2O, color="#990066", lw=2, label="D_SNK")
    plt.plot(1700 + np.arange(config.ind_final) + 1, alpha_N2O * (D_N2O[1:] - D_N2O[:-1]), color="#FFCC00", lw=2, label="d_N2O")
    plt.title("N2O fluxes (MtN/yr)", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + config.ind_final])

    # lifetime
    ax = plt.subplot(2, 3, 3)
    plt.plot(1700 + np.arange(config.ind_final + 1), alpha_N2O * (historical.N2O_0 + D_N2O) / (alpha_N2O * historical.N2O_0 / tau_N2O_hv - D_HVSNK_N2O), color="k", lw=2, label="OSCAR", )
    plt.title("N2O lifetime (yr)", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + config.ind_final])


# ========
# 2.4. O3
# ========


def plot_O3(D_O3t, D_O3s, D_EESC, D_N2O_lag, D_gst):
    plt.figure()

    # tropospheric O3
    ax = plt.subplot(2, 3, 1)
    plt.plot(1700 + np.arange(config.ind_final + 1), D_O3t, color="k", lw=2, label="OSCAR")
    plt.title("$\Delta$O3 trop. (DU)", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + config.ind_final])

    # stratospheric O3
    ax = plt.subplot(2, 3, 2)
    plt.plot(1700 + np.arange(config.ind_final + 1), D_O3s, color="k", lw=2, label="OSCAR")
    plt.title("$\Delta$O3 strat. (DU)", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + config.ind_final])

    # EESC
    ax = plt.subplot(2, 3, 3)
    plt.plot(1700 + np.arange(config.ind_final + 1), D_EESC, color="k", lw=2, label="OSCAR")
    plt.plot(1700 + np.arange(config.ind_final + 1), (chi_O3s_N2O * D_N2O_lag * (1 - D_EESC / EESC_x) / chi_O3s_EESC), color="k", lw=1, label="N2O effect", )
    plt.title("$\Delta$EESC (ppt)", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + config.ind_final])

    # age-of-air
    ax = plt.subplot(2, 3, 4)
    plt.plot(1700 + np.arange(config.ind_final + 1), tau_lag / (1 + gamma_age * D_gst), color="k", lw=2, label="OSCAR")
    plt.title("mean age-of-air (yr)", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + config.ind_final])


# ==============
# 2.5. Aerosols
# ==============


def plot_AER(D_SO4, D_POA, D_BC, D_NO3, D_SOA, D_AERh, RF_SO4, RF_POA, RF_BC, RF_NO3, RF_SOA, RF_cloud):
    plt.figure()

    # atmospheric burden
    ax = plt.subplot(2, 3, 1)
    plt.plot(1700 + np.arange(config.ind_final + 1), D_SO4, color="b", lw=2, label="D_SO4")
    plt.plot(1700 + np.arange(config.ind_final + 1), D_POA, color="m", lw=2, label="D_POA")
    plt.plot(1700 + np.arange(config.ind_final + 1), D_BC, color="r", lw=2, label="D_BC")
    plt.plot(1700 + np.arange(config.ind_final + 1), D_NO3, color="g", lw=2, label="D_NO3")
    plt.plot(1700 + np.arange(config.ind_final + 1), D_SOA, color="y", lw=2, label="D_SOA")
    plt.plot(1700 + np.arange(config.ind_final + 1), D_AERh, color="c", lw=2, label="D_AERh")
    plt.title("$\Delta$ burdens (Tg)", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + config.ind_final])

    # radiative forcing
    ax = plt.subplot(2, 3, 4)
    plt.plot([1700, 1700 + config.ind_final], [0, 0], "k-")
    plt.plot(1700 + np.arange(config.ind_final + 1), RF_SO4, color="b", lw=2, label="RF_SO4")
    plt.errorbar([2010], [-0.40], yerr=[[0.20], [0.20]], marker="o", mfc="b", color="k")
    plt.plot(1700 + np.arange(config.ind_final + 1), RF_POA, color="m", lw=2, label="RF_POA")
    plt.errorbar([2010], [-0.29], yerr=[[-0.29 * 0.63], [-0.29 * 0.72]], marker="o", mfc="m", color="k")
    plt.plot(1700 + np.arange(config.ind_final + 1), RF_BC, color="r", lw=2, label="RF_BC")
    plt.errorbar([2010], [+0.60], yerr=[[+0.60 * 0.61], [+0.60 * 0.70]], marker="o", mfc="r", color="k")
    plt.plot(1700 + np.arange(config.ind_final + 1), RF_NO3, color="g", lw=2, label="RF_NO3")
    plt.errorbar([2010], [-0.11], yerr=[[0.19], [0.08]], marker="o", mfc="g", color="k")
    plt.plot(1700 + np.arange(config.ind_final + 1), RF_SOA, color="y", lw=2, label="RF_SOA")
    plt.errorbar([2010], [-0.03], yerr=[[0.24], [0.23]], marker="o", mfc="y", color="k")
    plt.plot(1700 + np.arange(config.ind_final + 1), RF_cloud, color="c", lw=2, label="RF_cloud")
    plt.errorbar([2010], [-0.45], yerr=[[0.75], [0.45]], marker="o", mfc="c", color="k")
    # plt.errorbar([2010],[-0.10],yerr=[[0.20],[0.20]],marker='o',mfc='0.5',color='k')
    plt.title("RF (W/m2)", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, max(1700 + config.ind_final, 2010 + 10)])


# =============
# 2.6. Climate
# =============


def plot_clim(RF, D_gst, D_gyp, RF_CO2, RF_CH4, RF_H2Os, RF_N2O, RF_halo, RF_O3t, RF_O3s, RF_SO4, RF_POA, RF_BC, RF_NO3, RF_SOA, RF_cloud, RF_BCsnow, RF_LCC, RFcon, RFvolc, RFsolar,):
    plt.figure()

    # radiative forcing
    ax = plt.subplot(2, 3, 1)
    plt.plot([1700, 1700 + config.ind_final], [0, 0], "k-")
    plt.plot(1700 + np.arange(config.ind_final + 1), RF, color="k", lw=2, label="OSCAR")
    plt.plot(1700 + np.arange(len(historical.RF_ipcc)), historical.RF_ipcc, color="r", lw=2, ls="--", label="IPCC")
    plt.title("RF (W/m2)", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + config.ind_final])

    # global temperature
    ax = plt.subplot(2, 3, 2)
    plt.plot([1700, 1700 + config.ind_final], [0, 0], "k-")
    plt.plot(1700 + np.arange(config.ind_final + 1), D_gst - np.mean(D_gst[200:230]), color="k", lw=2, label="OSCAR")
    plt.plot(1700 + np.arange(len(gst_giss)), gst_giss - np.mean(gst_giss[200:230]), color="b", ls="--", label="GISS")
    plt.plot(1700 + np.arange(len(gst_had)), gst_had - np.mean(gst_had[200:230]), color="g", ls="--", label="Hadley")
    plt.plot(1700 + np.arange(len(gst_ncdc)), gst_ncdc - np.mean(gst_ncdc[200:230]), color="m", ls="--", label="NCDC")
    plt.title("$\Delta$ temp. (K)", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + config.ind_final])

    # global precipitations
    ax = plt.subplot(2, 3, 3)
    plt.plot(1700 + np.arange(config.ind_final + 1), D_gyp, color="k", lw=2, label="OSCAR")
    plt.title("$\Delta$ precip. (mm)", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + config.ind_final])

    # RF details
    ax = plt.subplot(2, 3, 4)
    plt.plot([1700, 1700 + config.ind_final], [0, 0], "k-")
    plt.plot(1700 + np.arange(config.ind_final + 1), RF_CO2 + RF_CH4 + RF_N2O + RF_halo + RF_H2Os, color="r", lw=2, label="WMGHG")
    plt.plot(1700 + np.arange(len(historical.RF_WMGHG_ipcc)), historical.RF_WMGHG_ipcc, color="r", ls="--")
    plt.plot(1700 + np.arange(config.ind_final + 1), RF_O3t + RF_O3s, color="y", lw=2, label="O3")
    plt.plot(1700 + np.arange(len(historical.RF_O3_ipcc)), historical.RF_O3_ipcc, color="y", ls="--")
    plt.plot(1700 + np.arange(config.ind_final + 1), RF_SO4 + RF_POA + RF_BC + RF_NO3 + RF_SOA + RF_cloud, color="b", lw=2, label="AER", )
    plt.plot(1700 + np.arange(len(historical.RF_AER_ipcc)), historical.RF_AER_ipcc, color="b", ls="--")
    plt.plot(1700 + np.arange(config.ind_final + 1), RF_BCsnow + RF_LCC, color="g", lw=2, label="Alb.")
    plt.plot(1700 + np.arange(len(historical.RF_Alb_ipcc)), historical.RF_Alb_ipcc, color="g", ls="--")
    plt.plot(1700 + np.arange(config.ind_final + 1), RFcon, color="k", ls="--", label="Ant.")
    plt.plot(1700 + np.arange(config.ind_final + 1), RFvolc + RFsolar, color="0.5", ls="--", label="Nat.")
    plt.title("RF (W/m2)", fontsize="medium")
    plt.legend(loc=0, ncol=2, prop=FontProperties(size="small"))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + config.ind_final])
