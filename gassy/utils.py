import numpy as np
from typing import Optional


def smooth(a:np.array,WSZ:Optional[int]=5)->np.array:
    """Python implementation of Matlab's `smooth` function

    Ref: https://stackoverflow.com/questions/40443020/

    Args:
        a (1d np.array): data to be smoothed
        WSZ (int): smoothing window size (must be odd number)

    Returns:
        1d np.array: smoothed data

    """
    avg = np.convolve(a,np.ones(WSZ,dtype=int),'valid')/WSZ
    r = np.arange(1,WSZ-1,2)
    start = np.cumsum(a[:WSZ-1])[::2]/r
    stop = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
    return np.concatenate((  start , avg, stop  ))
