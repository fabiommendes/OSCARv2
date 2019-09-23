import numpy as np

from ..data_loaders import nb_regionI, nb_biome
from ..params import nb_obox, nb_regionPF
from .. import config


def reduce_driver(driver):
    return np.sum(np.sum(np.sum(driver, 3), 2), 1)


def reduce_param(param):
    return np.sum(np.sum(np.sum(param, 2), 1), 0)


class Lazy:
    def __init__(self, fn):
        self.fn = fn
        self.name = None

    def __get__(self, instance, cls=None):
        if instance is None:
            return self
        instance.__dict__[self.name] = res = self.fn(instance)
        return res

    def __set_name__(self, owner, name):
        if self.name is None:
            self.name = name


scalar_timeseries = lambda: Lazy(lambda _: np.zeros([config.ind_final + 1], dtype=config.dty))
region_timeseries = lambda: Lazy(lambda _: np.zeros([config.ind_final + 1, nb_regionI], dtype=config.dty))
region_biome_timeseries = lambda: Lazy(lambda _: np.zeros([config.ind_final + 1, nb_regionI, nb_biome], dtype=config.dty))
region_biome2_age_timeseries = lambda: Lazy(lambda _: np.zeros([config.ind_final + 1, nb_regionI, nb_biome, nb_biome, config.ind_final + 1], dtype=config.dty))
obox_timeseries = lambda: Lazy(lambda _: np.zeros([config.ind_final + 1, nb_obox], dtype=config.dty))
species_timeseries = lambda n: Lazy(lambda _: np.zeros([config.ind_final + 1, n], dtype=config.dty))
regionPF_timeseries = lambda: Lazy(lambda _: np.zeros([config.ind_final + 1, nb_regionPF], dtype=config.dty))
land_var = lambda: Lazy(lambda _: np.zeros([nb_regionI, nb_biome], dtype=config.dty))
land_use_var = lambda: Lazy(lambda _: np.zeros([nb_regionI, nb_biome, nb_biome, config.ind_final + 1], dtype=config.dty))
scalar_var = lambda value=0.0: Lazy(lambda _: value)
linear_var = lambda n: Lazy(lambda _: np.zeros([n], dtype=config.dty))
