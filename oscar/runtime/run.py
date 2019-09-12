"""
Same as executing the OSCAR script with default values.

This modules exposes the OSCAR_lite that runs simulations.
"""
from ..config import load_oscar_script

globals().update(load_oscar_script(random=False))
