import glob
import os

from parse import parse


def match_paths(base, pattern, sort=True):
    full_pattern = os.path.join(base, pattern)
    paths = glob.glob(os.path.join(base, "*"))
    out = {}
    for path in paths:
        parsed = parse(full_pattern, path)
        if parsed is not None:
            out[path] = parsed.named

    if sort:
        out = dict(sorted(out.items(), key=lambda item: tuple(item[1].values())))

    return out


import numpy as np


def detect_threshold_crossings(signal: np.ndarray, threshold: float) -> np.ndarray:
    """
    Detects threshold crossings in a given signal where the signal transitions
    from below the threshold to above the threshold.

    Parameters
    ----------
    signal : np.ndarray
        1D array containing the data signal.
    threshold : float
        The threshold value to detect crossings.

    Returns
    -------
    np.ndarray
        Array of indices where the signal crosses the threshold from below.

    Examples
    --------
    >>> signal = np.array([0.1, 0.3, 0.2, 0.5, 0.7, 0.4, 1.0])
    >>> threshold = 0.4
    >>> detect_threshold_crossings(signal, threshold)
    array([3])
    """
    # Identify where signal is above the threshold
    above_threshold = signal >= threshold
    # Shift 'above_threshold' by one to align for comparison
    shift_above_threshold = np.roll(above_threshold, 1)
    # The first element of shifted version should not be True because it wraps around
    shift_above_threshold[0] = False
    # Find crossings: False to True transition
    crossings = (above_threshold & ~shift_above_threshold)
    # Get indices where crossings occur
    return np.where(crossings)[0]
