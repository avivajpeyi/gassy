from typing import Tuple, Callable, Optional
from scipy.interpolate import interp1d
from gassy.stellar_profiles.parser import read_profile
from gassy.stellar_profiles.plotter import plot_1d_profile_data, plot_desnity_grid
import numpy as np


class StellarProfile:
    def __init__(self, q, rho, c_s):
        self.q = smooth(q)
        self.minr, self.R = min(self.q), max(self.q)

        self._rho_interp, self.rho_range = self._build_rho_interp(self.q, rho)
        self.rho = np.vectorize(self._rho)

        self._c_s_interp = interp(self.q, c_s)
        self.c_s = np.vectorize(self._c_s)

        self._M_e_interp, self.M_e_range = self._build_M_e_interp(self.q)
        self.M_e = np.vectorize(self._M_e)

    @classmethod
    def load_matlab_profile(cls, profile_name):
        profile = read_profile(profile_name)
        return cls(profile.q, profile.rho, profile.c_s)

    def _build_rho_interp(self, q, rho):
        rho_interp = interp(q, rho)
        rho_range = (rho_interp(self.minr), rho_interp(self.R))
        return rho_interp, rho_range

    def _build_M_e_interp(self, q):
        integrand = 4 * np.pi * np.power(q, 2) * self.rho(q)
        m = np.array([np.trapz(integrand[0:i], x=self.q[0:i]) for i in range(2, len(integrand) + 1)])
        m_interp = interp(q[1:], m)
        m_range = (m_interp(self.minr), m_interp(self.R))
        return m_interp, m_range

    def _rho(self, r):
        rho = np.real(self._rho_interp(r))
        if rho > self.R:
            rho = 0
        return rho

    def _M_e(self, r):
        m = np.real(self._M_e_interp(r))
        return np.clip(m, *self.M_e_range)

    def _c_s(self, r):
        return np.abs(self._c_s_interp(r))

    def plot_grid(self, ax=None, savefig=False, **kwargs):
        return plot_desnity_grid(self, ax=ax, savefig=savefig,**kwargs)

    def plot_1d_profile_data(self):
        return plot_1d_profile_data(self)


def get_smooth_profile_functions(profile_name: str) -> Tuple[Callable, Callable, Callable]:
    profile = StellarProfile.load_matlab_profile(profile_name)
    return profile.rho, profile.c_s, profile.M_e


def smooth(a: np.array, WSZ: Optional[int] = 5) -> np.array:
    """Python implementation of Matlab's `smooth` function

    Ref: https://stackoverflow.com/questions/40443020/

    Args:
        a (1d np.array): data to be smoothed
        WSZ (int): smoothing window size (must be odd number)

    Returns:
        1d np.array: smoothed data

    """
    avg = np.convolve(a, np.ones(WSZ, dtype=int), "valid") / WSZ
    r = np.arange(1, WSZ - 1, 2)
    start = np.cumsum(a[: WSZ - 1])[::2] / r
    stop = (np.cumsum(a[:-WSZ:-1])[::2] / r)[::-1]
    return np.concatenate((start, avg, stop))


def interp(x, y):
    return interp1d(x, y, kind="cubic", fill_value="extrapolate", bounds_error=False)
