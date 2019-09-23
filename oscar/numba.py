try:
    from numba import jit, autojit, jitclass, prange
except ImportError:
    def jit(func):
        return func
    autojit = jitclass = autojit
    prange = range