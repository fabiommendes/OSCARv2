from scipy.optimize import fmin

from .a4_halogenated import *
from .a3_nitrous_oxide import tau_lag
from ..oscar_data import nb_regionI, regionI_index, nb_ODS, ODS
from ...config import dty, mod_O3Tregsat, mod_O3Temis, mod_O3Tclim, mod_O3Sfracrel, \
    mod_O3Strans, mod_O3Snitrous

##################################################
#   5. OZONE
##################################################

# ================
# 5.1. TROPOSPHERE
# ================

# -----------
# 5.1.1. HTAP
# -----------

# read region distribution
TMP = np.array(
    [line for line in csv.reader(
        open("data/RegDiv_HTAP/#DATA.RegDiv_HTAP.114reg1_(4reg0).AREA.csv", "r"))][1:],
    dtype=dty,
)
p_reg4 = np.zeros([nb_regionI, 4 + 1], dtype=dty)

for i in range(1, 114 + 1):
    p_reg4[regionI_index[i], :] += TMP[i - 1, :]
p_reg4 /= np.sum(p_reg4, 1)[:, np.newaxis]
p_reg4[np.isnan(p_reg4) | np.isinf(p_reg4)] = 0

# regional weight of ozone precursors {.}
# from HTAP experiments [Fry et al., 2013] (table S5) & [Fiore et al., 2013] (table S1)
if mod_O3Tregsat == "mean-HTAP":
    w_reg_NOX = (
            np.array([-3.31, -1.32, -0.67, -0.91, -0.41], dtype=dty)
            / np.array([27.5, 8.7, 8.4, 7.0, 3.3], dtype=dty)
            / -0.2
    )
    w_reg_CO = (
            np.array([-1.52, -0.42, -0.27, -0.52, -0.31], dtype=dty) / np.array(
        [461, 123, 88, 151, 98], dtype=dty) / -0.2
    )
    w_reg_VOC = (
            np.array([-2.53, -0.36, -0.42, -0.36, -2.53], dtype=dty)
            / np.array([191.8, 67.8, 39.1, 49.8, 35.2], dtype=dty)
            / -0.2
    )
elif mod_O3Tregsat == "CAMCHEM":
    w_reg_NOX = (
            np.array([-1.91, -0.62, -0.44, -0.55, -0.30], dtype=dty)
            / np.array([27.6, 8.2, 9.2, 6.8, 3.4], dtype=dty)
            / -0.2
    )
    w_reg_CO = (
            np.array([-1.49, -0.36, -0.22, -0.58, -0.33], dtype=dty) / np.array(
        [626, 154, 106, 221, 145], dtype=dty) / -0.2
    )
    w_reg_VOC = (
            np.array([-1.25, -0.29, -0.59, -0.21, -0.16], dtype=dty)
            / np.array([245.2, 92.4, 60.6, 55.5, 36.7], dtype=dty)
            / -0.2
    )
elif mod_O3Tregsat == "FRSGCUCI":
    w_reg_NOX = (
            np.array([-1.58, -0.50, -0.22, -0.62, -0.24], dtype=dty)
            / np.array([27.1, 8.8, 8.4, 6.9, 3.0], dtype=dty)
            / -0.2
    )
    w_reg_CO = (
            np.array([-1.54, -0.47, -0.28, -0.49, -0.30], dtype=dty) / np.array(
        [408, 131, 74, 129, 74], dtype=dty) / -0.2
    )
    w_reg_VOC = (
            np.array([-2.42, -0.60, -0.78, -0.68, -0.36], dtype=dty)
            / np.array([175.1, 57.3, 33.5, 50.5, 33.8], dtype=dty)
            / -0.2
    )
elif mod_O3Tregsat == "GISS-modelE":
    w_reg_NOX = (
            np.array([-2.84, -1.14, -0.55, -0.71, -0.44], dtype=dty)
            / np.array([31.4, 8.9, 8.6, 10.8, 3.1], dtype=dty)
            / -0.2
    )
    w_reg_CO = (
            np.array([-2.23, -0.61, -0.36, -0.90, -0.36], dtype=dty) / np.array(
        [417, 111, 69, 170, 67], dtype=dty) / -0.2
    )
    w_reg_VOC = (
            np.array([-0.41, -0.08, -0.08, -0.14, -0.11], dtype=dty)
            / np.array([111.7, 29.0, 26.9, 31.7, 24.1], dtype=dty)
            / -0.2
    )
elif mod_O3Tregsat == "GMI":
    w_reg_NOX = (
            np.array([-2.61, -0.96, -0.54, -0.75, -0.36], dtype=dty)
            / np.array([24.4, 8.0, 7.5, 5.8, 3.1], dtype=dty)
            / -0.2
    )
    w_reg_CO = (
            np.array([-2.19, -0.52, -0.38, -0.87, -0.42], dtype=dty) / np.array(
        [528, 134, 101, 194, 99], dtype=dty) / -0.2
    )
    w_reg_VOC = (
            np.array([-1.08, -0.31, -0.15, -0.38, -0.24], dtype=dty)
            / np.array([165.7, 60.1, 32.0, 42.0, 31.6], dtype=dty)
            / -0.2
    )
elif mod_O3Tregsat == "INCA":
    w_reg_NOX = (
            np.array([-3.08, -1.13, -0.56, -0.67, -0.72], dtype=dty)
            / np.array([26.9, 8.6, 7.3, 7.2, 3.8], dtype=dty)
            / -0.2
    )
    w_reg_CO = (
            np.array([-0.91, -0.18, -0.17, -0.33, -0.23], dtype=dty) / np.array(
        [376, 74, 67, 127, 108], dtype=dty) / -0.2
    )
    w_reg_VOC = (
            np.array([-1.05, -0.40, -0.35, -0.30, 0.00], dtype=dty)
            / np.array([279.1, 100.1, 56.7, 77.0, 45.3], dtype=dty)
            / -0.2
    )
elif mod_O3Tregsat == "LLNL-IMPACT":
    w_reg_NOX = (
            np.array([-1.99, -0.82, -0.36, -0.50, -0.31], dtype=dty)
            / np.array([30.0, 8.8, 9.7, 7.3, 4.2], dtype=dty)
            / -0.2
    )
    w_reg_CO = (
            np.array([-1.41, -0.31, -0.24, -0.40, -0.46], dtype=dty) / np.array(
        [518, 130, 111, 153, 124], dtype=dty) / -0.2
    )
    w_reg_VOC = (
            np.array([-0.22, -0.05, -0.10, -0.04, -0.03], dtype=dty)
            / np.array([143.4, 53.2, 23.3, 36.2, 30.7], dtype=dty)
            / -0.2
    )
elif mod_O3Tregsat == "MOZART-GFDL":
    w_reg_NOX = (
            np.array([-3.00, -1.18, -0.57, -0.81, -0.44], dtype=dty)
            / np.array([25.8, 9.4, 8.6, 5.2, 2.6], dtype=dty)
            / -0.2
    )
    w_reg_CO = (
            np.array([-1.29, -0.34, -0.35, -0.36, -0.24], dtype=dty) / np.array(
        [493, 124, 130, 134, 105], dtype=dty) / -0.2
    )
    w_reg_VOC = (
            np.array([-0.84, -0.19, -0.28, -0.20, -0.17], dtype=dty)
            / np.array([186.8, 68.7, 32.4, 48.9, 36.8], dtype=dty)
            / -0.2
    )
elif mod_O3Tregsat == "MOZECH":
    w_reg_NOX = (
            np.array([-2.24, -0.81, -0.58, -0.50, -0.35], dtype=dty)
            / np.array([26.9, 8.9, 7.4, 7.0, 3.6], dtype=dty)
            / -0.2
    )
    w_reg_CO = (
            np.array([-1.61, -0.46, -0.31, -0.48, -0.36], dtype=dty) / np.array(
        [470, 107, 85, 154, 124], dtype=dty) / -0.2
    )
    w_reg_VOC = (
            np.array([-1.94, -0.69, -0.55, -0.45, -0.25], dtype=dty)
            / np.array([250.6, 107.0, 44.7, 62.0, 36.9], dtype=dty)
            / -0.2
    )
elif mod_O3Tregsat == "STOC-HadAM3":
    w_reg_NOX = (
            np.array([-1.43, -0.59, -0.27, -0.46, -0.11], dtype=dty)
            / np.array([27.9, 9.0, 8.6, 7.1, 3.2], dtype=dty)
            / -0.2
    )
    w_reg_CO = (
            np.array([-1.29, -0.39, -0.25, -0.40, -0.25], dtype=dty) / np.array(
        [406, 129, 74, 127, 76], dtype=dty) / -0.2
    )
    w_reg_VOC = (
            np.array([-1.76, -0.40, -0.64, -0.45, -0.27], dtype=dty)
            / np.array([216.9, 70.9, 52.2, 50.0, 43.8], dtype=dty)
            / -0.2
    )
elif mod_O3Tregsat == "TM5-JRC":
    w_reg_NOX = (
            np.array([-3.80, -1.32, -0.72, -0.98, -0.78], dtype=dty)
            / np.array([26.6, 8.7, 8.5, 6.4, 3.0], dtype=dty)
            / -0.2
    )
    w_reg_CO = (
            np.array([-1.34, -0.54, -0.16, -0.42, -0.22], dtype=dty) / np.array(
        [399, 127, 75, 123, 74], dtype=dty) / -0.2
    )
    w_reg_VOC = (
            np.array([-1.91, -0.56, -0.49, -0.58, -0.28], dtype=dty)
            / np.array([164.3, 50.3, 34.6, 47.0, 32.4], dtype=dty)
            / -0.2
    )
elif mod_O3Tregsat == "UM-CAM":
    w_reg_NOX = (
            np.array([-2.53, -0.94, -0.50, -0.68, -0.41], dtype=dty)
            / np.array([27.9, 8.9, 8.7, 7.0, 3.3], dtype=dty)
            / -0.2
    )
    w_reg_CO = (
            np.array([-1.42, -0.41, -0.28, -0.44, -0.29], dtype=dty) / np.array(
        [428, 137, 81, 131, 79], dtype=dty) / -0.2
    )
    w_reg_VOC = (
            np.array([-1.91, -0.38, -0.64, -0.58, -0.31], dtype=dty)
            / np.array([171.1, 56.4, 32.7, 47.4, 34.6], dtype=dty)
            / -0.2
    )
elif mod_O3Tregsat == "":
    w_reg_NOX = np.array([1, 1, 1, 1, 1], dtype=dty)
    w_reg_CO = np.array([1, 1, 1, 1, 1], dtype=dty)
    w_reg_VOC = np.array([1, 1, 1, 1, 1], dtype=dty)
w_reg_NOX /= w_reg_NOX[0]
w_reg_CO /= w_reg_CO[0]
w_reg_VOC /= w_reg_VOC[0]

# -------------
# 5.1.2. ACCMIP
# -------------

# load pre-processed ACCMIP results for specified model
# for sensitivity to emissions
if mod_O3Temis != "mean-OxComp":
    def load_ozoe(var):
        TMP = np.array(
            [
                line
                for line in csv.reader(
                open(
                    "data/OzoChem_ACCMIP/#DATA.OzoChem_" + mod_O3Temis + ".2000s_(6exp)." + var + ".csv",
                    "r")
            )][1:],
            dtype=dty,
        )
        lgd = [
            line
            for line in csv.reader(
                open(
                    "data/OzoChem_ACCMIP/#DATA.OzoChem_" + mod_O3Temis + ".2000s_(6exp)." + var + ".csv",
                    "r")
            )][0]
        return TMP[0,:]

    O3t_ozoe = load_ozoe("O3t")
    CH4_ozoe = load_ozoe("CH4")
    ECO_ozoe = load_ozoe("ECO")
    EVOC_ozoe = load_ozoe("EVOC")
    ENOX_ozoe = load_ozoe("ENOX")

# setting of parameters
# based on ACCMIP
# to be used with formula by [Ehhalt et al., 2001]
if mod_O3Temis != "mean-OxComp":
    # sensitivity of tropospheric O3 to atmospheric CH4 {DU/.}
    chi_O3t_CH4 = np.array(
        [(O3t_ozoe[2] - O3t_ozoe[1]) / np.log(CH4_ozoe[2] / CH4_ozoe[1])], dtype=dty)
    # sensitivity of tropospheric O3 to ozoesions of ozone precursors {DU/{TgC/yr}}&{DU/{Tg/yr}}&{DU/{TgN/yr}}
    chi_O3t_CO = np.array([(O3t_ozoe[3] - O3t_ozoe[1]) / (ECO_ozoe[3] - ECO_ozoe[1])],
                          dtype=dty)
    chi_O3t_VOC = np.array([(O3t_ozoe[4] - O3t_ozoe[1]) / (EVOC_ozoe[4] - EVOC_ozoe[1])],
                           dtype=dty)
    chi_O3t_NOX = np.array([(O3t_ozoe[5] - O3t_ozoe[1]) / (ENOX_ozoe[5] - ENOX_ozoe[1])],
                           dtype=dty)
    # normalization of the non-linearity
    chi_O3t_CH4 *= (O3t_ozoe[1] - O3t_ozoe[0]) / (
            4 * O3t_ozoe[1] - O3t_ozoe[2] - O3t_ozoe[3] - O3t_ozoe[4] - O3t_ozoe[5]
    )
    chi_O3t_CO *= (O3t_ozoe[1] - O3t_ozoe[0]) / (
            4 * O3t_ozoe[1] - O3t_ozoe[2] - O3t_ozoe[3] - O3t_ozoe[4] - O3t_ozoe[5]
    )
    chi_O3t_VOC *= (O3t_ozoe[1] - O3t_ozoe[0]) / (
            4 * O3t_ozoe[1] - O3t_ozoe[2] - O3t_ozoe[3] - O3t_ozoe[4] - O3t_ozoe[5]
    )
    chi_O3t_NOX *= (O3t_ozoe[1] - O3t_ozoe[0]) / (
            4 * O3t_ozoe[1] - O3t_ozoe[2] - O3t_ozoe[3] - O3t_ozoe[4] - O3t_ozoe[5]
    )

# taken from the TAR [Ehhalt et al., 2001]
elif mod_O3Temis == "mean-OxComp":
    chi_O3t_CH4 = np.array([5.0], dtype=dty)
    chi_O3t_CO = np.array([0.0011], dtype=dty) * 28 / 12.0
    chi_O3t_VOC = np.array([0.0033], dtype=dty)
    chi_O3t_NOX = np.array([0.125], dtype=dty)

# load pre-processed ACCMIP results for specified model
# for sensitivity to climate
if mod_O3Tclim != "":
    TMP = load_data(f"data/OzoChem_ACCMIP/#DATA.OzoChem_{mod_O3Tclim}.2000s_(4exp).O3t.csv", slice=1)
    O3t_ozoc = TMP[0,:]

    TMP = load_data(f"data/OzoChem_ACCMIP/#DATA.OzoChem_{mod_O3Tclim}.2000s_(4exp).tas.csv", slice=1)
    tas_ozoc = TMP[0, :]

# definition of parameter
# sensitivity of tropospheric O3 to climate change {DU/K}
Gamma_O3t = np.array([0], dtype=dty)

# fit of parameter
if mod_O3Tclim != "":
    diff = O3t_ozoc - O3t_ozoc[0]


    def err(var):
        clim = var[0] * (tas_ozoc[:] - tas_ozoc[0])
        return np.sum((diff - clim) ** 2)


    [Gamma_O3t[0]] = fmin(err, [0], disp=False)

# =================
# 5.2. STRATOSPHERE
# =================

# ---------------
# 5.2.1. Chlorine
# ---------------

# fractional release factors {.}
# fits based on [Newman et al., 2006]
if mod_O3Sfracrel in ["Newman2006"]:

    def f_fracrel(tau):
        fracrel = np.zeros([nb_ODS], dtype=dty)
        fracrel[ODS.index(
            "CFC11")] = 6.53976e-02 * tau + 4.18938e-02 * tau ** 2 + -3.68985e-03 * tau ** 3
        fracrel[ODS.index(
            "CFC12")] = 4.06080e-02 * tau + 6.74585e-04 * tau ** 2 + 3.70165e-03 * tau ** 3
        fracrel[ODS.index("CFC113")] = (
                                               5.36055e-02 * tau + 7.16185e-08 * tau ** 2) * np.exp(
            tau / 4.95237e00)
        fracrel[ODS.index("CFC114")] = (
                                               2.01017e-02 * tau + -4.67409e-08 * tau ** 2) * np.exp(
            tau / 4.28272e00)
        fracrel[ODS.index("CFC115")] = (
                                               3.55000e-03 * tau + 7.67310e-04 * tau ** 2) * np.exp(
            tau / 4.33197e00)
        fracrel[ODS.index("CCl4")] = (
                                             8.50117e-02 * tau + 8.33504e-02 * tau ** 2) * np.exp(
            -tau / 5.12976e00)
        fracrel[ODS.index("CH3CCl3")] = (
                                                1.63270e-01 * tau + 1.12119e-01 * tau ** 2) * np.exp(
            -tau / 3.74978e00)
        fracrel[ODS.index(
            "HCFC22")] = 3.02660e-02 * tau + 2.49148e-04 * tau ** 2 + 1.38022e-03 * tau ** 3
        fracrel[ODS.index("HCFC141b")] = (
                                                 2.89263e-03 * tau + 1.14686e-08 * tau ** 2) * np.exp(
            tau / 1.36374e00)
        fracrel[ODS.index("HCFC142b")] = (
                                                 1.22909e-05 * tau + 1.06509e-04 * tau ** 2) * np.exp(
            tau / 1.23520e00)
        fracrel[ODS.index("Halon1211")] = (
                                                  -1.51872e-02 * tau + 1.68718e-01 * tau ** 2) * np.exp(
            -tau / 3.43897e00)
        fracrel[ODS.index("Halon1202")] = (
                                                  -1.51872e-02 * tau + 1.68718e-01 * tau ** 2) * np.exp(
            -tau / 3.43897e00)
        fracrel[ODS.index("Halon1301")] = (
                                                  5.46396e-02 * tau + 3.03982e-08 * tau ** 2) * np.exp(
            tau / 5.63826e00)
        fracrel[ODS.index(
            "Halon2402")] = 2.25980e-01 * tau + 2.23266e-04 * tau ** 2 + -1.13882e-03 * tau ** 3
        fracrel[ODS.index("CH3Br")] = (
                                              7.60148e-02 * tau + 1.12771e-01 * tau ** 2) * np.exp(
            -tau / 4.09494e00)
        fracrel[ODS.index(
            "CH3Cl")] = 1.39444e-01 * tau + 1.46219e-04 * tau ** 2 + 8.40557e-04 * tau ** 3
        return fracrel


    def df_fracrel_dtau(tau):
        fracrel = np.zeros([nb_ODS], dtype=dty)
        fracrel[ODS.index(
            "CFC11")] = 6.53976e-02 + 2 * 4.18938e-02 * tau + 3 * -3.68985e-03 * tau ** 2
        fracrel[ODS.index(
            "CFC12")] = 4.06080e-02 + 2 * 6.74585e-04 * tau + 3 * 3.70165e-03 * tau ** 2
        fracrel[ODS.index("CFC113")] = (
                                               5.36055e-02
                                               + 2 * 7.16185e-08 * tau
                                               + (5.36055e-02 / 4.95237e00) * tau
                                               + (7.16185e-08 / 4.95237e00) * tau ** 2
                                       ) * np.exp(tau / 4.95237e00)
        fracrel[ODS.index("CFC114")] = (
                                               2.01017e-02
                                               + 2 * -4.67409e-08 * tau
                                               + (2.01017e-02 / 4.28272e00) * tau
                                               + (-4.67409e-08 / 4.28272e00) * tau ** 2
                                       ) * np.exp(tau / 4.28272e00)
        fracrel[ODS.index("CFC115")] = (
                                               3.55000e-03
                                               + 2 * 7.67310e-04 * tau
                                               + (3.55000e-03 / 4.33197e00) * tau
                                               + (7.67310e-04 / 4.33197e00) * tau ** 2
                                       ) * np.exp(tau / 4.33197e00)
        fracrel[ODS.index("CCl4")] = (
                                             8.50117e-02
                                             + 2 * 8.33504e-02 * tau
                                             - (8.50117e-02 / 5.12976e00) * tau
                                             - (8.33504e-02 / 5.12976e00) * tau ** 2
                                     ) * np.exp(-tau / 5.12976e00)
        fracrel[ODS.index("CH3CCl3")] = (
                                                1.63270e-01
                                                + 2 * 1.12119e-01 * tau
                                                - (1.63270e-01 / 3.74978e00) * tau
                                                - (1.12119e-01 / 3.74978e00) * tau ** 2
                                        ) * np.exp(-tau / 3.74978e00)
        fracrel[ODS.index(
            "HCFC22")] = 3.02660e-02 + 2 * 2.49148e-04 * tau + 3 * 1.38022e-03 * tau ** 2
        fracrel[ODS.index("HCFC141b")] = (
                                                 2.89263e-03
                                                 + 2 * 1.14686e-08 * tau
                                                 + (2.89263e-03 / 1.36374e00) * tau
                                                 + (1.14686e-08 / 1.36374e00) * tau ** 2
                                         ) * np.exp(tau / 1.36374e00)
        fracrel[ODS.index("HCFC142b")] = (
                                                 1.22909e-05
                                                 + 2 * 1.06509e-04 * tau
                                                 + (1.22909e-05 / 1.23520e00) * tau
                                                 + (1.06509e-04 / 1.23520e00) * tau ** 2
                                         ) * np.exp(tau / 1.23520e00)
        fracrel[ODS.index("Halon1211")] = (
                                                  -1.51872e-02
                                                  + 2 * 1.68718e-01 * tau
                                                  - (-1.51872e-02 / 3.43897e00) * tau
                                                  - (1.68718e-01 / 3.43897e00) * tau ** 2
                                          ) * np.exp(-tau / 3.43897e00)
        fracrel[ODS.index("Halon1202")] = (
                                                  -1.51872e-02
                                                  + 2 * 1.68718e-01 * tau
                                                  - (-1.51872e-02 / 3.43897e00) * tau
                                                  - (1.68718e-01 / 3.43897e00) * tau ** 2
                                          ) * np.exp(-tau / 3.43897e00)
        fracrel[ODS.index("Halon1301")] = (
                                                  5.46396e-02
                                                  + 2 * 3.03982e-08 * tau
                                                  + (5.46396e-02 / 5.63826e00) * tau
                                                  + (3.03982e-08 / 5.63826e00) * tau ** 2
                                          ) * np.exp(tau / 5.63826e00)
        fracrel[ODS.index(
            "Halon2402")] = 2.25980e-01 + 2 * 2.23266e-04 * tau + 3 * -1.13882e-03 * tau ** 2
        fracrel[ODS.index("CH3Br")] = (
                                              7.60148e-02
                                              + 2 * 1.12771e-01 * tau
                                              + (7.60148e-02 / 4.09494e00) * tau
                                              + (1.12771e-01 / 4.09494e00) * tau ** 2
                                      ) * np.exp(-tau / 4.09494e00)
        fracrel[ODS.index(
            "CH3Cl")] = 1.39444e-01 + 2 * 1.46219e-04 * tau + 3 * 8.40557e-04 * tau ** 2
        return fracrel

    # alternative values (mid-latitudes) from [Laube et al., 2013]
# missing values still from [Newman et al., 2006]
elif mod_O3Sfracrel in ["Laube2013"]:

    def f_fracrel(tau):
        fracrel = np.zeros([nb_ODS], dtype=dty)
        fracrel[ODS.index("CFC11")] = -0.0173 + 0.098_666 * tau + 0.008_169_55 * tau ** 2
        fracrel[ODS.index("CFC12")] = -0.0154 + 0.046_244 * tau + 0.007_073_56 * tau ** 2
        fracrel[ODS.index("CFC113")] = -0.0059 + 0.049_669 * tau + 0.008_624_13 * tau ** 2
        fracrel[ODS.index("CFC114")] = (
                                               2.01017e-02 * tau + -4.67409e-08 * tau ** 2) * np.exp(
            tau / 4.28272e00
        )  # not given
        fracrel[ODS.index("CFC115")] = (
                                               3.55000e-03 * tau + 7.67310e-04 * tau ** 2) * np.exp(
            tau / 4.33197e00
        )  # not given
        fracrel[ODS.index("CCl4")] = -0.0139 + 0.131_338 * tau + 0.004_648_06 * tau ** 2
        fracrel[
            ODS.index("CH3CCl3")] = -0.0227 + 0.254_820 * tau + -0.015_059_46 * tau ** 2
        fracrel[ODS.index("HCFC22")] = -0.0190 + 0.021_203 * tau + 0.002_294_34 * tau ** 2
        fracrel[
            ODS.index("HCFC141b")] = -0.0635 + 0.050_362 * tau + 0.009_399_38 * tau ** 2
        fracrel[
            ODS.index("HCFC142b")] = -0.0032 + 0.010_130 * tau + 0.002_331_85 * tau ** 2
        fracrel[
            ODS.index("Halon1211")] = -0.0535 + 0.204_371 * tau + -0.004_646_44 * tau ** 2
        fracrel[ODS.index("Halon1202")] = (
                                                  -1.51872e-02 * tau + 1.68718e-01 * tau ** 2) * np.exp(
            -tau / 3.43897e00
        )  # not given
        fracrel[
            ODS.index("Halon1301")] = -0.0185 + 0.061_608 * tau + 0.010_518_28 * tau ** 2
        fracrel[ODS.index("Halon2402")] = (
                2.25980e-01 * tau + 2.23266e-04 * tau ** 2 + -1.13882e-03 * tau ** 3
        )  # not given
        fracrel[ODS.index("CH3Br")] = (
                                              7.60148e-02 * tau + 1.12771e-01 * tau ** 2) * np.exp(
            -tau / 4.09494e00
        )  # not given
        fracrel[ODS.index(
            "CH3Cl")] = 1.39444e-01 * tau + 1.46219e-04 * tau ** 2 + 8.40557e-04 * tau ** 3  # not given
        return fracrel


    def df_fracrel_dtau(tau):
        fracrel = np.zeros([nb_ODS], dtype=dty)
        fracrel[ODS.index("CFC11")] = 0.098_666 + 2 * 0.008_169_55 * tau
        fracrel[ODS.index("CFC12")] = 0.046_244 + 2 * 0.007_073_56 * tau
        fracrel[ODS.index("CFC113")] = 0.049_669 + 2 * 0.008_624_13 * tau
        fracrel[ODS.index("CFC114")] = (
                                               2.01017e-02
                                               + 2 * -4.67409e-08 * tau
                                               + (2.01017e-02 / 4.28272e00) * tau
                                               + (-4.67409e-08 / 4.28272e00) * tau ** 2
                                       ) * np.exp(
            tau / 4.28272e00
        )  # not given
        fracrel[ODS.index("CFC115")] = (
                                               3.55000e-03
                                               + 2 * 7.67310e-04 * tau
                                               + (3.55000e-03 / 4.33197e00) * tau
                                               + (7.67310e-04 / 4.33197e00) * tau ** 2
                                       ) * np.exp(
            tau / 4.33197e00
        )  # not given
        fracrel[ODS.index("CCl4")] = 0.131_338 + 2 * 0.004_648_06 * tau
        fracrel[ODS.index("CH3CCl3")] = 0.254_820 + 2 * -0.015_059_46 * tau
        fracrel[ODS.index("HCFC22")] = 0.021_203 + 2 * 0.002_294_34 * tau
        fracrel[ODS.index("HCFC141b")] = 0.050_362 + 2 * 0.009_399_38 * tau
        fracrel[ODS.index("HCFC142b")] = 0.010_130 + 2 * 0.002_331_85 * tau
        fracrel[ODS.index("Halon1211")] = 0.204_371 + 2 * -0.004_646_44 * tau
        fracrel[ODS.index("Halon1202")] = (
                                                  -1.51872e-02
                                                  + 2 * 1.68718e-01 * tau
                                                  - (-1.51872e-02 / 3.43897e00) * tau
                                                  - (1.68718e-01 / 3.43897e00) * tau ** 2
                                          ) * np.exp(
            -tau / 3.43897e00
        )  # not given
        fracrel[ODS.index("Halon1301")] = 0.061_608 + 2 * 0.010_518_28 * tau
        fracrel[ODS.index(
            "Halon2402")] = 2.25980e-01 + 2 * 2.23266e-04 * tau + 3 * -1.13882e-03 * tau ** 2  # not given
        fracrel[ODS.index("CH3Br")] = (
                                              7.60148e-02
                                              + 2 * 1.12771e-01 * tau
                                              + (7.60148e-02 / 4.09494e00) * tau
                                              + (1.12771e-01 / 4.09494e00) * tau ** 2
                                      ) * np.exp(
            -tau / 4.09494e00
        )  # not given
        fracrel[ODS.index(
            "CH3Cl")] = 1.39444e-01 + 2 * 1.46219e-04 * tau + 3 * 8.40557e-04 * tau ** 2  # not given
        return fracrel

    # others values for EESC {.}
# relative strength of bromine from [Daniel et al., 2007]
alpha_Br = np.array([60.0], dtype=dty)
n_Cl = np.array([3, 2, 3, 2, 1, 4, 3, 1, 2, 1, 1, 0, 0, 0, 0, 1], dtype=dty)
n_Br = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 2, 1, 0], dtype=dty)
EESC_0 = np.sum(f_fracrel(tau_lag) * (n_Cl + alpha_Br * n_Br) * ODS_0)

# --------------
# 5.2.2. CCMVal2
# --------------

# surface temperature trend from CCMVal2
# load pre-processed CCMVal2 results for specified model
for var in ["ta2"]:
    ta2_atm = np.zeros([2099-1961+1], dtype=dty)
    TMP = np.array(
        [
            line
            for line in csv.reader(
            open(
                "data/Atmosphere_CCMVal2/#DATA.Atmosphere_" + mod_O3Strans + ".1961-2099_(1lvl)." + var + ".csv",
                "r",
            )
        )][1:],
        dtype=dty,
    )
    ta2_atm[:] = TMP[:,0]

# fit of linear yearly trend
ta_trend = np.array([0], dtype=dty)
diff = ta2_atm[:] - np.mean(ta2_atm[:10])


def err(var):
    return np.sum((diff - var[0] * np.arange(len(ta2_atm))) ** 2)


[ta_trend[0]] = fmin(err, [0], disp=False)

# sensitivity of strato ozone to chlorine and climate {DU/ppt}&{DU/K}
# from [Douglass et al., 2014] (figure 2; data provided by author)
if mod_O3Strans == "mean-CCMVal2":
    chi_O3s_EESC = np.array([-12.5e-3], dtype=dty)
    Gamma_O3s = np.array([0.012], dtype=dty) / ta_trend
elif mod_O3Strans == "AMTRAC":
    chi_O3s_EESC = np.array([-13.7e-3], dtype=dty)
    Gamma_O3s = np.array([-0.015], dtype=dty) / ta_trend
elif mod_O3Strans == "CCSR-NIES":
    chi_O3s_EESC = np.array([-7.9e-3], dtype=dty)
    Gamma_O3s = np.array([0.013], dtype=dty) / ta_trend
elif mod_O3Strans == "CMAM":
    chi_O3s_EESC = np.array([-8.1e-3], dtype=dty)
    Gamma_O3s = np.array([-0.15], dtype=dty) / ta_trend
elif mod_O3Strans == "CNRM-ACM":
    chi_O3s_EESC = np.array([-23.2e-3], dtype=dty)
    Gamma_O3s = np.array([0.0066], dtype=dty) / ta_trend
elif mod_O3Strans == "LMDZrepro":
    chi_O3s_EESC = np.array([-10.6e-3], dtype=dty)
    Gamma_O3s = np.array([0.042], dtype=dty) / ta_trend
elif mod_O3Strans == "MRI":
    chi_O3s_EESC = np.array([-20.3e-3], dtype=dty)
    Gamma_O3s = np.array([-0.030], dtype=dty) / ta_trend
elif mod_O3Strans == "Niwa-SOCOL":
    chi_O3s_EESC = np.array([-5.2e-3], dtype=dty)
    Gamma_O3s = np.array([0.087], dtype=dty) / ta_trend
elif mod_O3Strans == "SOCOL":
    chi_O3s_EESC = np.array([-10.5e-3], dtype=dty)
    Gamma_O3s = np.array([0.067], dtype=dty) / ta_trend
elif mod_O3Strans == "ULAQ":
    chi_O3s_EESC = np.array([-15.4e-3], dtype=dty)
    Gamma_O3s = np.array([0.062], dtype=dty) / ta_trend
elif mod_O3Strans == "UMSLIMCAT":
    chi_O3s_EESC = np.array([-10.6e-3], dtype=dty)
    Gamma_O3s = np.array([0.0071], dtype=dty) / ta_trend
elif mod_O3Strans == "UMUKCA-UCAM":
    chi_O3s_EESC = np.array([-11.5e-3], dtype=dty)
    Gamma_O3s = np.array([0.04], dtype=dty) / ta_trend

# sensitivity of stratospheric O3 to N2O {DU/ppb}&{DU/ppb/ppt}
# formulation and values from [Daniel et al., 2010]
if mod_O3Snitrous == "Daniel2010":
    chi_O3s_N2O = chi_O3s_EESC * f_fracrel(3)[0] * 6.4 * (1.53 + 0.53 * 240 / 1400.0)
    EESC_x = np.array([1400 / 0.53], dtype=dty)
else:
    chi_O3s_N2O = 0 * chi_O3s_EESC
    EESC_x = np.array([np.inf], dtype=dty)
