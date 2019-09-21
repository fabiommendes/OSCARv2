import numpy as np

from ..config import ind_final, dty
from ..runtime.oscar_data import nb_regionI, nb_biome
from ..runtime.oscar_param import nb_obox, nb_regionPF


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


mk_scalar = lambda: Lazy(lambda _: np.zeros([ind_final + 1], dtype=dty))
mk_region_var = lambda: Lazy(lambda _: np.zeros([ind_final + 1, nb_regionI], dtype=dty))
mk_region_biome_var = lambda: Lazy(lambda _: np.zeros([ind_final + 1, nb_regionI, nb_biome], dtype=dty))
mk_region_biome2_age_var = lambda: Lazy(lambda _: np.zeros([ind_final + 1, nb_regionI, nb_biome, nb_biome, ind_final + 1], dtype=dty))
mk_obox_var = lambda: Lazy(lambda _: np.zeros([ind_final + 1, nb_obox], dtype=dty))
mk_species_var = lambda n: Lazy(lambda _: np.zeros([ind_final + 1, n], dtype=dty))
mk_regionPF_var = lambda: Lazy(lambda _: np.zeros([ind_final + 1, nb_regionPF], dtype=dty))
mk_land_var = lambda: Lazy(lambda _: np.zeros([nb_regionI, nb_biome], dtype=dty))
mk_land_use_var = lambda: Lazy(lambda _: np.zeros([nb_regionI, nb_biome, nb_biome, ind_final + 1], dtype=dty))
mk_scalar_var = lambda value=0: Lazy(lambda _: np.array([value], dtype=dty))
mk_linear_var = lambda n: Lazy(lambda _: np.zeros([n], dtype=dty))