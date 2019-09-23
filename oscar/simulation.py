import numpy as np


class external:
    __slots__ = 'name', 'varname', 'default'

    def __init__(self, varname, default=None):
        self.varname = varname
        self.default = default
        self.name = None

    def __get__(self, instance, cls=None):
        pass

    def __set__(self, instance, value):
        pass

    def __set_name__(self, owner, name):
        if self.name is None:
            self.name = name


class Simulation:
    """
    State of simulation.
    """
    __slots__ = '_variables',

    def __init__(self, variables):
        self._variables = {}
        for v in variables:
            self.add_variable(v)

    def add_variable(self, v):
        """"""

    def step(self):
        """
        Advance simulation by a single step.
        """
        forces = [var.force() for var in self._variables]
        for var, force in zip(self._variables.values(), forces):
            var.step(self.dt, force=force)


class Parameter:
    """
    Any parameter that depends on the a simulation.
    """
    __slots__ = 'simulation', 'name'
    _index = 0

    def __init__(self, simulation, name=None):
        self.simulation = simulation
        if name is None:
            Parameter._index += 1
            name = f'v{Parameter._index}'
        self.name = name


class Emissions:
    """
    Represent emissions including historical data_loaders.
    """


class Variable:
    __slots__ = 'config', 'state', 'value', 'timeseries', '__dict__'

    def __init__(self, config, state, value, dtype=float):
        self.value = value
        self.diff = value * 0
        self.config = config
        self.state = state
        self.timeseries = np.array([value], dtype=dtype)
        self.__dict__ = {}

    def update(self, value):
        pass

    def step(self, dt):
        pass

    def force(self):
        pass

