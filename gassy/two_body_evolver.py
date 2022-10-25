"""
[adapted from spiralling113MESA15]
"""
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np

from .constants import G, Msol, Rsol, pi, c
from .conversions import get_mu, get_period, mag

from gassy.stellar_profiles import get_smooth_profile_functions
from gassy.ode_driver import ode_driver


def get_evolution_initial_conditions(M, m, a, e, mesa_profile_name, ):
    rho_func, c_func, M_e_func = get_smooth_profile_functions(mesa_profile_name)
    M_e = M_e_func(a)
    T = get_period(M, m, a)
    L = np.sqrt(G * (M + m + M_e * a * (1 - np.power(e, 2))))
    y0 = [(1 + e), 0, 0, L * T / (np.power(a, 2) * (1 + e))]
    return T, y0


def evolve_bodies(M, m, a, e, Tend, mesa_profile_name, num_points=1000) -> Tuple:
    """
    Solves the equations of motion for a two-body system in a common envelope
    with polytropic index p, i.e. P ~ rho^(1+1/p).

    Args:
        M (float): description

    Returns:
        X (type): description
    """
    T, y0 = get_evolution_initial_conditions(M, m, a, e, mesa_profile_name)
    rho_func, c_func, M_e_func = get_smooth_profile_functions(mesa_profile_name)

    tspan = np.array([0, Tend]) / T
    args = (M, m, a, c_func, rho_func, M_e_func)
    # Integrator
    t = np.linspace(0, Tend, 1000000)
    y = ode_driver(twobody_in_gas_dydt_ode, t=t, y0=y0, args=args, rtol=1e-8, atol=1e-8)

    # unpack and return data
    y = y.T
    X = a * y[0, :]
    v_x = a * y[1, :] / T
    Y = a * y[2, :]
    v_y = a * y[3, :] / T
    t = t * T
    return X, Y, v_x, v_y, t


def get_N_steps(M: float, v: float, a: float, R: float, T: float) -> int:
    # TODO: what is this?

    b_90 = max([G * M / np.power(v * a / T, 2), 1e-1 * Rsol])
    N = R / b_90
    return N


def get_term_A(v, rdot, M, m, r, nu):
    c2, c4, c5 = np.power(c, 2), np.power(c, 4), np.power(c, 5)
    v2, v4 = np.power(v, 2), np.power(v, 4)
    rdot2, rdot3, rdot4 = np.power(rdot, 2), np.power(rdot, 3), np.power(rdot, 4)
    nu2 = np.power(nu, 2)
    GMmr = G * (M + m) / r
    GMmr2 = np.power(GMmr, 2)
    nuv2 = nu * v2
    nu2v2 = nu2 * v2
    nuv4 = nu * v4
    nu2v4 = nu2 * v4

    GMmTerm = (
            GMmr * (-2 * rdot2 - 25 * rdot2 * nu - 2 * rdot2 * nu2 - 13 * nuv2 / 2 + 2 * nu2v2) +
            GMmr2 * (9 + 87 * nu / 4)
    )

    line1 = (1 / c2) * (-3 * rdot2 * (nu / 2) + v2 + 3 * nuv2 - GMmr * (4 + 2 * nu))
    line2 = (1 / c4) * (
            15 * rdot4 * (nu / 8) - 45 * rdot4 * (nu2 / 8) - 9 * rdot2 * (nuv2 / 2) +
            6 * rdot2 * nu2v2 + 3 * nuv4 - 4 * nu2v4 + GMmTerm
    )
    line3 = (1 / c5) * (-24 * nuv2 * rdot * GMmr / 5 - (136 * rdot * nu / 15) * GMmr2)

    # TODO: what is this? where are these numbers coming from?
    A = line1 + line2 + line3
    return A


def get_term_B(v, rdot, M, m, r, nu):
    c2, c4, c5 = np.power(c, 2), np.power(c, 4), np.power(c, 5)
    v2, v4 = np.power(v, 2), np.power(v, 4)
    rdot2, rdot3, rdot4 = np.power(rdot, 2), np.power(rdot, 3), np.power(rdot, 4)
    nu2 = np.power(nu, 2)
    GMmr = G * (M + m) / r
    GMmr2 = np.power(GMmr, 2)
    B = (
            1 / c2 * (-4 * rdot + 2 * rdot * nu)
            + 1
            / c4
            * (
                    9 * rdot3 * nu / 2
                    + 3 * rdot3 * nu2
                    - 15 * rdot * nu * v2 / 2
                    - 2 * rdot * nu2 * v2
                    + GMmr * (2 * rdot + 41 * rdot * nu / 2 + 4 * rdot * nu2)
            )
            + 1 / c5 * ((8 * nu * v2 / 5) * GMmr + 24 * nu / 5 * GMmr2)
    )
    return B


def compute_drag_coef(v, a, T, c_s, N, M, rho, constant_drag_coef=0.0):
    if constant_drag_coef != None:
        return constant_drag_coef

    # TODO: what is this?
    Mach = (v * a / T) / c_s

    # TODO: what is this?
    # Note f is cut at 1 so as not to diverge
    h1 = max([1 / N, 1e-2])
    N2 = np.power(N, 2)
    h = np.power((1 - (2 - h1) * np.exp(-2 + 2 * h1) / (N2 * h1)), (-1 / 2)) - 1

    # Ostriker # TODO: what is this?
    if Mach < 1 - h1:
        I = 0.5 * np.log((1 + Mach) / (1. - Mach)) - Mach
    else:
        if (Mach >= 1 - h1) and (Mach < 1 + h):
            I = 0.5 * np.log((2 - h1) / h1) - 1 + h1
            print("Speed of sound reached")
        else:
            I = 0.5 * np.log(1 - 1 / np.power(Mach, 2)) + np.log(N)
            if Mach < 1:
                print("Mach < 1 !")  # should this be an error?

    assert not np.isnan(I), f"I is a nan!"
    assert I >= 0, "I < 0!"

    # TODO: what is this?
    f = (I * 4 * pi * np.power(G * M, 2) * rho) / ((a * v / T) ** 3)

    return f


def twobody_in_gas_dydt_ode(data, t, M, m, a, c_s_func, rho_func, M_e_func):
    dydt = np.zeros(4)
    R = 111.3642 * Rsol
    a = 0.7 * R  # why is this hardcoded?
    M = 15 * Msol
    m = 1 * Msol
    mu = get_mu(M, m)
    T = get_period(M, m, a)

    if np.isnan(data).any():
        raise ValueError("NaN in data")

    x, y = data[0], data[2]
    v_x, v_y = data[1], data[3]

    r = mag([x, y])
    v = mag([v_x, v_y])
    c_s = c_s_func(a * r)
    rho = rho_func(a * r)

    # TODO: i think this can be simplified? We should have the mathematica code that auto-generates alot of this...
    # Store the mathematica code as C++

    N = get_N_steps(M, v, a, R, T)
    f = compute_drag_coef(v, a, T, c_s, N, M, rho)

    rdot = (x * v_x + y * v_y) / r
    rdot = rdot * a / T
    nu = mu / (M + m)

    # TODO: wut why?
    v = a * v / T
    r = a * r

    # TODO: what is this? where are these numbers coming from?
    A = get_term_A(v, rdot, M, m, r, nu)
    B = get_term_B(v, rdot, M, m, r, nu)

    r = r / a
    r2, r3 = np.power(r, 2), np.power(r, 3)
    pi2 = np.power(pi, 2)

    M_e = M_e_func(a * r)

    # wut mate?
    dydt[0] = y
    dydt[1] = (
            -4 * pi2 * ((1 + A) * x / r + B * v_x * a / T) / r2
            - T * f * v_x / M
            - 4 * pi2 * ((M_e - m) / (M + m)) * x / r3
    )

    dydt[2] = v_y
    dydt[3] = (
            -4 * pi2 * ((1 + A) * y / r + B * v_y * a / T) / r2
            - T * f * v_y / M
            - 4 * pi2 * ((M_e - m) / (M + m)) * y / r3
    )

    Rt = 1e-1 * Rsol * max([(M / m) ** (1 / 3), (m / M) ** (1 / 3)])

    if r * a <= Rt:
        dydt[0] = 0
        dydt[1] = 0
        dydt[2] = 0
        dydt[3] = 0
        print(t * T)

    if np.isnan(dydt).any():
        raise ValueError("NaN in dydt")

    return dydt
