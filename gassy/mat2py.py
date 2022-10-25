"""Matlab functions in python"""
from typing import Optional

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d




def interp(vector, factor):
    """
    Interpolate a given 1D vector by a specific interpolation factor.
    :param vector: 1D data vector
    :param factor: factor for interpolation (must be integer)
    :return: interpolated 1D vector by a given factor

    Ref: https://gist.github.com/piotrdzwiniel/64b2dff08ae629bc73909e6ceb6d7f35
    """
    x = np.arange(np.size(vector))
    y = vector
    f = interp1d(x, y)

    x_extended_by_factor = np.linspace(x[0], x[-1], np.size(x) * factor)
    y_interpolated = np.zeros(np.size(x_extended_by_factor))

    i = 0
    for x in x_extended_by_factor:
        y_interpolated[i] = f(x)
        i += 1

    return y_interpolated


def interp1(x, v, xq, method="linear", extrapolation="extrap"):
    """
    vq = interp1(x,v,xq,method,extrapolation) specifies a strategy for
    evaluating points that lie outside the domain of x.
    Set extrapolation to 'extrap' when you want to use the method algorithm
    for extrapolation.

    Alternatively, you can specify a scalar value,
    in which case, interp1 returns that value for
    all points outside the domain of x.
    """
    f = interp1d(x, v, kind="linear", fill_value="extrapolate")
    return f(xq)


def integral(fun, xmin, xmax):
    """
    https://au.mathworks.com/help/matlab/ref/integral.html
    https://docs.scipy.org/doc/scipy/tutorial/integrate.html
    """
    val, err = quad(fun, xmin, xmax)  # , args=fun_args)
    return val
