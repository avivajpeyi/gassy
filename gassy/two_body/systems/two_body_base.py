import warnings
from .orbit_type import OrbitType
from typing import Optional

import numpy as np

from gassy.constants import G, pi





class TwoBodyBase:
    def __init__(
            self,
            m: float,
            M: float,
            init_x: float,
            init_vy: float,
    ):
        """
        :param m: mass of the smaller body (kg)
        Args:
            m:
            M:
            init_x:
            init_vy:
        """
        self.m = m
        self.M = M
        self.r = np.array([init_x, 0])
        self.v = np.array([0, init_vy])
        self.init_y = self.pack_data()
        self.bound_orbit_check()

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

    @property
    def period(self) -> float:
        m, M, a = (self.m, self.Me, self.rmag)
        return np.sqrt((4 * pi ** 2) * (a ** 3) / (G * (M + m)))

    @property
    def gravitational_force(self) -> np.ndarray:
        f_mag = -G * self.m * self.M / self.rmag ** 2
        return f_mag * self.rhat

    @property
    def drag_force(self) -> np.ndarray:
        return np.array([0,0])

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
        q = np.zeros(4)
        x, y, vx, vy, ax, ay = np.concatenate([self.r, self.v, self.accel])
        q[0] = 2 * (x * ax + vx ** 2)
        q[1] = x * ay + y * ax + 2 * vx * vy
        q[2] = q[1]
        q[3] = 2 * (y * ay + vy ** 2)
        return q * self.mu

    @property
    def orbit_type(self) -> OrbitType:
        if self.Egpe + self.Ek < 0:
            return OrbitType.BOUND
        return OrbitType.UNBOUND

    def update(self, pos_vel_array: np.ndarray) -> np.ndarray:
        """Updates the system's position and velocity (pos_vel_array: [x, y, vx, vy])"""
        self.r = pos_vel_array[0:2]
        self.v = pos_vel_array[2:4]
        return self.pack_data()

    def pack_data(self) -> np.ndarray:
        """Returns [x, y, vx, vy, Ek, Egpe, L, M00 M01 M10, M11]"""
        return np.array(
            [*self.r, *self.v, self.Ek, self.Egpe, self.L, *self.mass_moment]
        )

    @property
    def label(self):
        return f"Point particles: {self._m_label}"

    @property
    def _m_label(self):
        msun = r"$M_{\odot}$"
        return f"m={self.m} {msun}, M={self.M} {msun}"

    def __repr__(self):
        return f"{self.label} [{self.orbit_type}]"

    @property
    def Egpe(self) -> float:
        """Gravitational Potential Energy"""
        return -G * self.m * self.Me / self.rmag

    @property
    def Ek(self) -> float:
        """Kinetic Energy"""
        return 0.5 * self.m * self.vmag ** 2

    @property
    def L(self) -> float:
        return self.m * np.cross(self.r, self.v)


def mag(x, axis=None):
    return np.linalg.norm(x, axis=axis)
