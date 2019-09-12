"""
Same as executing the OSCAR script with random values.

This modules exposes the OSCAR_lite that runs simulations.
"""
import sys

from ..config import load_oscar_script

globals().update(load_oscar_script(random=True))
_glob = globals()


def reload_params():
    frame = sys._getframe(1)
    new = load_oscar_script(random=True)
    _glob.update(new)
    frame.f_globals.update(new)
    return new
