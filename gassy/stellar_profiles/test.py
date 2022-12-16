import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.constants import G, M_sun, R_sun
from scipy.integrate import odeint, quad
from scipy.interpolate import interp1d
from scipy.optimize import root

MSol = M_sun.cgs.value
R_sun = R_sun.cgs.value
G = G.cgs.value
pi = np.pi


def get_rho(Pi, gamma, S):
    return np.power(Pi / S, 1.0 / gamma)


def get_P(gamma, S, rhoi):
    return S * np.power(rhoi, gamma)


def get_Pj(Pi, mi, ri, rj, rhoi):
    dr = rj - ri
    dP_dr = -(G * mi * rhoi) / np.power(ri, 2)
    Pj = dP_dr * dr + Pi
    return Pj


def get_mj(mi, ri, rj, rhoi):
    dm_dr = (4.0 / 3.0) * pi * (rj**3 - ri**3) * rhoi
    dr = rj - ri
    return mi + dm_dr * dr


def get_jth_values(Pi, mi, ri, rhoi, dr, S0, gamma):
    rj = ri + dr
    Pj = get_Pj(Pi, mi, ri, rj, rhoi)
    rhoj = get_rho(Pi=Pj, S=S0, gamma=gamma)
    mj = get_mj(mi, ri, rj, rhoi)
    return rj, Pj, rhoj, mj


class Star:
    def __init__(self, r, P, rho, m):
        self.r = r
        self.P = P
        self.rho = rho
        self.m = m

    @classmethod
    def build(cls, S0, rho0, gamma, m0):

        P0 = get_P(rhoi=rho0, S=S0, gamma=gamma)

        r_guess = cls.get_system_radii(P0, m0, 0.01, rho0, S0, gamma)

        r, dr = np.linspace(0.01, r_guess, 1000, retstep=True)
        P = np.zeros(len(r))
        rho = np.zeros(len(r))
        m = np.zeros(len(r))

        rho[0], m[0], P[0] = rho0, m0, P0

        for j in range(1, len(r) - 1):
            Pi, mi, ri, rhoi = P[j - 1], m[j - 1], r[j - 1], rho[j - 1]
            _, P[j], rho[j], m[j] = get_jth_values(Pi, mi, ri, rhoi, dr, S0, gamma)

        # mask data to where rho > 0
        rmask = rho > 0
        r = r[rmask]
        P = P[rmask]
        rho = rho[rmask]
        m = m[rmask]

        # if len(r) < 10:
        #     raise ValueError("Too few points where ddensity is above 0")
        # save only 1000 points

        return cls(r, P, rho, m)

    @staticmethod
    def get_system_radii(P, m, r, rho, S0, gamma):
        num_it, dr = 0, R_sun / 1e5

        r = dr
        m = 4 / 3 * pi * r**3
        while P > 0:
            num_it += 1
            r, P, rho, m = get_jth_values(P, m, r, rho, dr, S0, gamma)
        return r - dr

    def plot_profile(self):
        """plot star rho, mass, pressure"""
        plt.close("all")
        fig, ax = plt.subplots(1, 3, figsize=(12, 4))
        ax[0].loglog(self.r / R_sun, self.rho)
        ax[0].set_xlabel("r [R_sun]")
        ax[0].set_ylabel("density")
        ax[1].loglog(self.r / R_sun, self.m / MSol)
        ax[1].set_xlabel("r [R_sun]")
        ax[1].set_ylabel("mass [M_sun]")
        ax[2].loglog(self.r / R_sun, self.P)
        ax[2].set_xlabel("r [m]")
        ax[2].set_ylabel("Pressure")
        # use sf in x axis

        plt.tight_layout()
        plt.show()

    def to_dataframe(self):
        return pd.DataFrame({"r": self.r, "P": self.P, "rho": self.rho, "m": self.m})


star = Star.build(S0=1e15, rho0=1, gamma=5.0 / 3.0, m0=1)
star.plot_profile()
print(star.to_dataframe())
