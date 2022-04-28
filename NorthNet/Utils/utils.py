import hashlib
import numpy as np


def sha1_hash(text, num_chars=7):
    """
    Generate a SHA1 hash of a string with a defined length (<= 40 chars).

    Parameters
    ----------
    text: string

    length: int

    Returns
    -------
    sha1_hash: str
    """

    full_hash = hashlib.sha1(str.encode(text)).hexdigest()

    first_n_chars = full_hash[:num_chars]

    return first_n_chars


def is_float(thing):
    """
    Test if an thing (e.g. str) can be converted to a float.

    Parameters
    ----------
    x: any type

    Returns
    -------
    bool
    """

    try:
        float(thing)
        return True
    except ValueError:
        return False


def is_int(x):
    """
    Test if variable can be converted to an integer.

    Parameters
    ----------
    x: any type

    Returns
    -------
    bool
    """

    try:
        int(x)
        return True
    except ValueError:
        return False


def indices_from_boundary(data, start, end):
    """
    Get the indices of elements of the data array between start and end.

    Parameters
    ----------
    data: ndarray
        Array from which indice will be found
    start: float
        Lower boundary.
    end: float
        Upper boundary.

    Returns
    -------
    indices: ndarray
        Indices where start < data < end.
    """

    pre_indices = np.where((data >= start) & (data <= end))

    indices = pre_indices[0]

    return indices
