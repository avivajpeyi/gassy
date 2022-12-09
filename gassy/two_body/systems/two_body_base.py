import warnings
from typing import Optional

import numpy as np

from gassy.constants import G, Msol, Rsol, kgm_s2, m_per_s, m_per_s2, pi

from .orbit_type import OrbitType


class TwoBodyBase:
    def __init__(
        self,
        m: float,
        M: float,
        r: float,
        init_vy: Optional[float] = None,
        continue_on_error: Optional[bool] = False,
    ):
        self.m = m * Msol
        self.M = M * Msol
        self.r = np.array([r, 0]) * Rsol

        if init_vy is None:
            init_vy = self._orbital_vel(self.rmag)

        self.v = np.array([0, init_vy])
        self.init_y = self.pack_data()
        self.bound_orbit_check(continue_on_error)

    def bound_orbit_check(self, continue_on_error=False):
        msg = (
            f"{OrbitType} orbit: "
            f"|Init Vel| {mag(self.v):.3f} > |Escape Vel| ({self.escape_vel:.3f}). "
            f"Are the init conditions correct?"
        )
        if mag(self.v) > self.escape_vel or self.orbit_type == OrbitType.UNBOUND:
            if continue_on_error:
                warnings.warn(msg)
            else:
                raise ValueError(msg)

    @property
    def rmag(self):
        return mag(self.r)

    @property
    def rhat(self):
        return self.r / self.rmag

    @property
    def vmag(self):
        return mag(self.v)

    @property
    def vhat(self):
        return self.v / self.vmag

    @property
    def Me(self):
        """Mass of inner body inside from 0 to a"""
        return self.M

    @property
    def escape_vel(self):
        return np.sqrt(2 * G * self.Me / self.rmag)

    def _orbital_vel(self, radius: float) -> float:
        return np.sqrt(np.abs(G * self.Me / radius))

    @property
    def period(self) -> float:
        m, M, a = (self.m, self.Me, self.rmag)
        return np.sqrt((4 * pi**2) * (a**3) / (G * (M + m)))

    @property
    def gravitational_force(self) -> np.ndarray:
        f_mag = -G * self.m * self.M / self.rmag**2
        return f_mag * self.rhat

    @property
    def drag_force(self) -> np.ndarray:
        return np.array([0, 0])

    @property
    def net_force(self) -> np.ndarray:
        """Returns net [Fx, Fy] on m"""
        return self.gravitational_force + self.drag_force

    @property
    def accel(self) -> np.ndarray:
        return self.net_force / self.m

    @property
    def mu(self) -> float:
        return self.Me * self.m / (self.Me + self.m)

    @property
    def mass_moment(self):
        """ddot M"""
        x, y, vx, vy, ax, ay = np.concatenate([self.r, self.v, self.accel])
        M00 = 2 * (x * ax + vx**2)
        M01 = x * ay + y * ax + 2 * vx * vy
        M10 = M01
        M11 = 2 * (y * ay + vy**2)
        return np.array([M00, M01, M10, M11]) * self.mu

    @property
    def orbit_type(self) -> OrbitType:
        if self.Egpe + self.Ek < 0:
            return OrbitType.BOUND
        return OrbitType.UNBOUND

    def update(self, pos_vel_array: np.ndarray):
        """Updates the system's position and velocity (pos_vel_array: [x, y, vx, vy])"""
        self.r = pos_vel_array[0:2]
        self.v = pos_vel_array[2:4]

    def pack_data(self, all_data: Optional[bool] = False) -> np.ndarray:
        """pack two-body data into an array

        if all_data:
            Returns [x, y, vx, vy, Ek, Egpe, L, M00 M01 M10, M11]
        else:
            Returns [x, y, vx, vy]
        """
        if all_data:
            data = np.array(
                [*self.r, *self.v, self.Ek, self.Egpe, self.L, *self.mass_moment]
            )
        else:
            data = np.array([*self.r, *self.v])
        return data

    @property
    def data_dim(self):
        return self.pack_data(all_data=True).shape[0]

    @property
    def label(self):
        return f"NO DRAG Point particles: {self._m_label}"

    @property
    def _m_label(self):
        msun = r"$M_{\odot}$"
        return f"m={self.m / Msol} {msun}, M={self.M / Msol} {msun}"

    def __repr__(self):
        return f"{self.label} [{self.orbit_type}]"

    @property
    def Egpe(self) -> float:
        """Gravitational Potential Energy"""
        return -G * self.m * self.Me / self.rmag

    @property
    def Ek(self) -> float:
        """Kinetic Energy"""
        return 0.5 * np.power(self.vmag, 2)

    @property
    def Etot(self) -> float:
        return self.Egpe + self.Ek

    @property
    def L(self) -> float:
        return self.m * np.cross(self.r, self.v)

    @property
    def eccentricity(self):
        """Eccentricity"""
        return np.sqrt(
            1 + 2 * self.Etot * self.L**2 / (G**2 * self.m**2 * self.Me**2)
        )

    @property
    def a(self):
        """Semi-major axis"""
        return self.rmag / (1 - self.eccentricity)


def mag(x, axis=None):
    return np.linalg.norm(x, axis=axis)
