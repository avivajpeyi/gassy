"""
[adapted from spiralling113MESA15]
"""
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np

from .constants import G, Msol, Rsol, pi
from .conversions import get_mu, get_period, mag
from .mat2py import integral, interp1, smooth
from .stellar_profiles import read_profile


def smooth_rho(r, R, rho, q) -> np.float64:
    """what is this doing?
    rho = @(r) (r<=R).*real(interp1(q,Rho,r,'linear','extrap')) + (r > R).*0;
    """
    if r <= R:
        return np.real(interp1(q, rho, r, "linear", "extrap"))
    else:
        return 0


def smooth_c(r, q, c_s) -> np.float64:
    """what is this doing?
    C = @(r) abs(interp1(q,c_s,r,'linear','extrap'));
    """
    return np.abs(interp1(q, c_s, r, "linear", "extrap"))


def smooth_M_e(x, R, rho, q) -> np.float64:
    """what is this doing?
    M_e = @(x)
    (x<=R).*integral(@(r) 4*pi*r.^2.*rho(r),0,x) +
    (x > R).*integral(@(r) 4*pi*r.^2.*rho(r),0,R);
    """
    integrand = lambda r: 4 * pi * r**2 * smooth_rho(r, R, rho, q)
    if x <= R:
        return x * integral(integrand, 0, x)
    else:
        return x * integral(integrand, 0, R)


def get_evolution_initial_conditions(M, m, a, e, mesa_profile_name):
    profile = read_profile(mesa_profile_name)
    rho = profile.rho
    q = profile.q * Rsol  # conversion to cgs
    c_s = profile.c_s * 100  # conversion to cgs (i think?)
    R = max(q)
    q = smooth(q)
    M_e = smooth_M_e(a, R, rho, q)

    T = get_period(M, m, a)
    L = np.sqrt(G * (M + m + M_e * a * (1 - e**2)))
    y0 = [(1 + e), 0, 0, L * T / (a**2 * (1 + e))]
    return T, y0


def evolve_bodies(M, m, a, e, dt, Tend, mesa_profile_name) -> Tuple:
    """
    Solves the equations of motion for a two-body system in a common envelope
    with polytropic index p, i.e. P ~ rho^(1+1/p).

    Args:
        M (float): description
        m (float): description
        a (float): description
        e (float): description
        dt (float): description
        Tend (float): description

    Returns:
        X (type): description
        Y (type): description
        v_x (type): description
        v_y (type): description
        t (type): description
    """

    X, Y, v_x, v_y, t = None, None, None, None, None
    T, y0 = get_evolution_initial_conditions(M, m, a, e, mesa_profile_name)

    tspan = [0, Tend] / T
    opts = odeset("RelTol", 1e-10, "Stats", "on")

    # Integrator
    # [t,y] = ode113(@odefun,tspan,y0,opts);

    # unpack and return data
    X = a * data[1, :]
    v_x = a * data[2, :] / T
    Y = a * data[3, :]
    v_y = a * data[4, :] / T
    t = t * T
    return X, Y, v_x, v_y, t


def odefun(t, data, M, m, a, orig_c_s, q, orig_rho, R):
    dydt = np.zeros(4)

    M = M * Msol
    m = m * Msol
    mu = get_mu(M, m)
    T = get_period(M, m, a)

    x, y = data[0], data[2]
    v_x, v_y = data[1], data[3]

    r = mag([x, y])
    v = mag([v_x, v_y])
    c_s = smooth_c(a * r, q, orig_c_s)
    rho = smooth_rho(a * r, R, orig_rho, q)
    M_e = smooth_M_e(a * r, R, orig_rho, q)

    # TODO: what is this?
    b_90 = max([G * M / (v * a / T) ** 2, 1e-1 * Rsol])
    N = R / b_90

    # TODO: what is this?
    Mach = (v * a / T) / c_s

    # TODO: what is this?
    # Note f is cut at 1 so as not to diverge
    h1 = max([1 / N, 1e-2])
    h = (1 - (2 - h1) * np.exp(-2 + 2 * h1) / (N**2 * h1)) ** (-1 / 2) - 1

    # Ostriker # TODO: what is this?
    if Mach < 1 - h1:
        I = 0.5 * log((1 + Mach) / (1 - Mach)) - Mach
    else:
        if (Mach >= 1 - h1) and (Mach < 1 + h):
            I = 0.5 * np.log((2 - h1) / h1) - 1 + h1
            print("Speed of sound reached")
        else:
            I = 0.5 * log(1 - 1 / Mach ^ 2) + log(N)
            if Mach < 1:
                print("Wrong condition")

    # TODO: what is this?
    f = I * 4 * pi * G**2 * M**2 * rho / ((a * v / T) ** 3)
    if I < 0:
        print("boo")  # TODO: am i raising an error?

    rdot = (x * v_x + y * v_y) / r
    rdot = rdot * a / T
    nu = mu / (M + m)

    v = a * v / T
    r = a * r

    c2 = c**2
    c4 = c**4
    c5 = c**5
    v2 = v**2
    rdot2 = rdot**2
    rdot3 = rdot**3
    rdot4 = rdot**4
    nu2 = nu**2
    GMmr = G * (M + m) / r
    GMmr2 = (GMm / r) ** 2
    pi2 = pi**2

    # TODO: what is this? where are these numbers coming from?
    A = (
        1
        / c2
        * (-3 * rdot2 * nu / 2 + v2 + 3 * nu * v2 - G * (M + m) * (4 + 2 * nu) / r)
        + 1
        / c4
        * (
            15 * rdot4 * nu / 8
            - 45 * rdot4 * nu2 / 8
            - 9 * rdot2 * nu * v2 / 2
            + 6 * rdot2 * nu2 * v2
            + 3 * nu * v4
            - 4 * nu2 * v4
            + GMmr
            * (
                -2 * rdot2
                - 25 * rdot2 * nu
                - 2 * rdot2 * nu2
                - 13 * nu * v2 / 2
                + 2 * nu2 * v2
            )
            + GMmr2 * (9 + 87 * nu / 4)
        )
        + 1 / c5 * (-24 * nu * rdot * v2 * GMmr / 5 - (136 * rdot * nu / 15) * GMmr2)
    )

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

    r = r / a

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

    return dydt
