import numpy as np
import csv
from functools import lru_cache
from .config import dty


@lru_cache(512)
def load_data(path, slice=None, dtype=dty, use='csv'):
    """
    Return an array with data loaded from given path.
    """
    if use == 'numpy' and slice in (0, None):
        return np.genfromtxt(path, dtype=dtype)
    elif use == 'numpy' and slice == 1:
        return np.genfromtxt(path, dtype=dtype, skip_header=True)
    else:
        with open(path, "r") as fd:
            data = [line for line in csv.reader(fd)]
        if slice is not None:
            data = data[slice:]
        return np.array(data, dtype=dtype)


def load_data_and_header(path):
    """
    Return (header, array) with a list of column names and numeric data loaded
    from given path.
    """
    header, *data = csv.reader(open(path, "r"))
    return np.array(data, dtype=dty), header


def load_header(path):
    """
    Return a list with the header from a table stored as CSV data.
    """
    return next(iter(csv.reader(open(path, "r"))))
