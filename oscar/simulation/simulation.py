import matplotlib.pyplot as plt

from ..constants import HFC, ODS, PFC
from .. import conf
from .drivers import EFF, ECH4, EN2O, RFcon, RFvolc, RFsolar
from .atmosphere import AtmosphereSimulatorMixin
from .chemestry import ChemestrySimulatorMixin
from .climate import ClimateSimulatorMixin
from .land import LandSimulatorMixin
from .land_use import LandUseSimulatorMixin
from .ocean import OceanSimulationMixin
from .permafrost import PermafrostSimulatorMixin
from .util import scalar_timeseries, region_timeseries, species_timeseries

nb_HFC = len(HFC)
nb_PFC = len(PFC)
nb_ODS = len(ODS)


# noinspection PyAttributeOutsideInit
class Simulation(OceanSimulationMixin,
                 LandSimulatorMixin,
                 LandUseSimulatorMixin,
                 PermafrostSimulatorMixin,
                 ChemestrySimulatorMixin,
                 AtmosphereSimulatorMixin,
                 ClimateSimulatorMixin):

    # Global scalar variables
    D_pH_t = scalar_timeseries()
    OSNK_t = scalar_timeseries()
    D_FOXI_CH4_t = scalar_timeseries()

    # (region) variables
    D_lst_t = region_timeseries()
    D_lyp_t = region_timeseries()

    # (species) variables
    D_HFC_t = species_timeseries(nb_HFC)
    D_HFC_lag_t = species_timeseries(nb_HFC)
    D_PFC_t = species_timeseries(nb_PFC)
    D_PFC_lag_t = species_timeseries(nb_PFC)
    D_ODS_t = species_timeseries(nb_ODS)
    D_ODS_lag_t = species_timeseries(nb_ODS)

    def __init__(self, p=conf.p, fT=conf.fT, track=["ELUC", "OSNK", "LSNK", "D_CO2", "RF", "D_gst"], plot='all'):
        super().__init__()
        var_plot = self.plot_vars(plot)
        self.fT = fT
        self.dt = 1 / p
        self.p = p

        # save variables
        self.tracking = list(set(track) | set(var_plot))
        self.time = 1

    def run(self):
        """
        Run simulation to final time.
        """
        for t in range(1, conf.ind_final + 1):
            self.time = t
            self.step(t)

    def step(self, t):
        """
        Advance 1 step in simulation.
        """
        for tt in range(self.p):
            self.step_ocean(t)
            self.step_land(t)
            self.step_land_use(t)
            self.step_permafrost(t)
            self.step_chemestry(t)
            self.step_atmosphere(t)
            self.step_climate(t)

        # Save tracked time series variables
        for name in self.tracking:
            arr = getattr(self, name + '_t')
            run = getattr(self, name)
            arr[t] = run

    def output(self, variables=None):
        variables = self.tracking if variables is None else variables
        return [getattr(self, name + '_t') for name in variables]

    def plot_vars(self, plot):
        # plot variables
        var_plot = []
        if plot == "all" or plot == "CO2" or "CO2" in plot:
            var_plot += ["D_CO2", "OSNK", "LSNK", "ELUC", "D_AREA", "D_npp", "D_efire", "D_fmort", "D_rh1", "D_fmet", "D_rh2", "D_FIN", "D_FOUT", "D_FCIRC", "EFIRE_luc", "FMORT_luc", "RH1_luc", "FMET_luc", "RH2_luc", "EHWP1_luc",
                         "EHWP2_luc", "EHWP3_luc", ]
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
            plot_CO2(self.D_CO2_t, self.OSNK_t, self.LSNK_t, self.ELUC_t, EFF, self.D_AREA_t, self.D_npp_t, self.D_efire_t, self.D_fmort_t, self.D_rh1_t, self.D_fmet_t, self.D_rh2_t, self.D_FIN_t, self.D_FOUT_t, self.D_FCIRC_t,
                     self.EFIRE_luc_t, self.FMORT_luc_t, self.RH1_luc_t, self.RH2_luc_t, self.EHWP1_luc_t, self.EHWP2_luc_t, self.EHWP3_luc_t)
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
            plot_clim(self.RF_t, self.D_gst_t, self.D_gyp_t, self.RF_CO2_t, self.RF_CH4_t, self.RF_H2Os_t, self.RF_N2O_t, self.RF_halo_t, self.RF_O3t_t, self.RF_O3s_t, self.RF_SO4_t, self.RF_POA_t, self.RF_BC_t, self.RF_NO3_t,
                      self.RF_SOA_t, self.RF_cloud_t, self.RF_BCsnow_t, self.RF_LCC_t, RFcon, RFvolc, RFsolar)
            plt.savefig("results/plot-clim.svg")


def run(track=(), plot="all"):
    sim = Simulation(track=track, plot=plot)
    print('STARTING SIMULATION')
    sim.run()
    print('PLOTING VARIABLES')
    sim.plot(plot)
    return sim.output