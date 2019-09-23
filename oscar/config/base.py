import os
import sys

path = os.path.dirname(os.path.dirname(__file__))
script_path = os.path.join(path, "simulation", "simulation.py")


def random_config_module():
    """
    Return a config module with random values for initialization variables.
    """

    sys.modules.pop("oscar.config.random", "")
    from oscar.config import random as mod

    return mod


def random_config_vars():
    """
    Return a dictionary mapping configuration variables to their corresponding
    values.
    """

    mod = random_config_module()
    return {k: v for k, v, in vars(mod).items() if not k.startswith("_")}


def load_oscar_script(random=False):
    """
    Load Oscar script and return a dictionary with all computed variables.
    """
    # Load config
    if random:
        ns = random_config_vars()
    else:
        from oscar.config import default

        ns = vars(default)

    # Load script
    with open(script_path) as fd:
        src = fd.read()
    exec(src, ns)
    return ns
