import numpy as np
import csv
import os
from functools import lru_cache
from .config import dty


@lru_cache(128)
def load_data(path, start=None, dtype=dty, use='csv'):
    """
    Return an array with data loaded from given path.
    """

    cached = os.path.join('cache', path) + '.npy'
    try:
        return np.load(cached)
    except FileNotFoundError:
        pass

    if use == 'numpy' and start in (0, None):
        data = np.genfromtxt(path, dtype=dtype)
    elif use == 'numpy' and start == 1:
        data = np.genfromtxt(path, dtype=dtype, skip_header=True)
    else:
        with open(path, "r") as fd:
            data = [line for line in csv.reader(fd)]
        if start is not None:
            data = data[start:]
        data = np.array(data, dtype=dtype)

    if dtype is not object:
        os.makedirs(os.path.dirname(cached), exist_ok=True)
        np.save(cached, data)
    return data


@lru_cache(128)
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
