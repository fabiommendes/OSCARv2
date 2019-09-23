class State:
    """
    Simulation state.
    """
    __slots__ = 'time',

    def __init__(self, time=0):
        self.time = time
