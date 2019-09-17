import csv

import numpy as np

from .config import dty


def load_data(path, slice=None):
    """
    Return an array with data loaded from given path.
    """
    data = [line for line in csv.reader(open(path, "r"))]
    if slice is not None:
        data = data[slice:]
    return np.array(data, dtype=dty)


def load_header_and_table(path):
    """
    Return (header, array) with a list of column names and numeric data loaded
    from given path.
    """
    header, *data = csv.reader(open(path, "r"))
    return header, np.array(data, dtype=dty)