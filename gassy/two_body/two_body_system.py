import warnings

from typing import Optional, Tuple, Union
import numpy as np
from gassy.conversions import mag, get_period
from gassy.constants import G, Msol, Rsol, pi, c
from enum import Enum


class OrbitType(Enum):
    BOUND = 0
    UNBOUND = 1


class TwoBodySystem:
    def __init__(self, m: float, M: float, init_x: float, init_vy: float, drag_coeff: Optional[float] = None):
        self.m = m
        self.M = M
        self.pos = np.array([init_x, 0])
        self.v = np.array([0, init_vy])
        self.init_y = self.get_pos_vel()
        self.drag_coeff = drag_coeff if drag_coeff is not None else 0
        if self.orbit_type == OrbitType.UNBOUND:
            warnings.warn("Unbound orbit -- are the init conditions correct?")

    @property
    def a(self) -> float:
        return mag(self.pos)

    @property
    def init_a(self) -> float:
        return mag(self.init_y[:2])

    @property
    def period(self, use_init=False) -> float:
        a = self.a if not use_init else self.init_a
        return get_period(self.m, self.M, a)

    @property
    def gravitational_energy(self) -> float:
        return self.__get_grav_potential_energy(self.m, self.M, self.pos)

    @property
    def kinetic_energy(self) -> float:
        return self.__get_kinetic_energy(self.m, self.v)

    @property
    def gravitational_force(self) -> np.ndarray:
        f_mag = -G * self.m * self.M / mag(self.pos)
        return f_mag * self.pos / self.a

    @property
    def drag_force(self) -> np.ndarray:
        return -self.drag_coeff * self.v / mag(self.v)

    @property
    def net_force(self) -> np.ndarray:
        """Returns net [Fx, Fy] on m"""
        return self.gravitational_force + self.drag_force

    @property
    def net_acceleration(self) -> np.ndarray:
        return self.net_force / self.m

    @property
    def orbit_type(self) -> OrbitType:
        if self.kinetic_energy > self.gravitational_energy:
            return OrbitType.BOUND
        return OrbitType.UNBOUND

    def update(self, pos_vel_array: np.ndarray) -> None:
        """Updates the system's position and velocity (pos_vel_array: [x, y, vx, vy])"""
        self.pos = pos_vel_array[:2]
        self.v = pos_vel_array[2:]

    def get_pos_vel(self) -> np.ndarray:
        """Returns [x, y, vx, vy]"""
        return np.concatenate((self.pos, self.v))

    def get_energies(self, cached_pos, cached_v) -> Tuple[np.ndarray, np.ndarray]:
        kinetic_energies = self.__get_kinetic_energy(self.m,  cached_v, axis=1)
        gravitational_energies = self.__get_grav_potential_energy(self.m, self.M, cached_pos, axis=1)
        return kinetic_energies, gravitational_energies

    @property
    def label(self):
        return f"m={self.m}, M={self.M}, drag={self.drag_coeff}"

    @staticmethod
    def __get_kinetic_energy(m: float,  v: np.ndarray, axis=None) -> Union[float, np.ndarray]:
        return 0.5 * m * mag(v, axis=axis) ** 2

    @staticmethod
    def __get_grav_potential_energy(m: float, M: float, r: np.ndarray, axis=None) -> Union[float, np.ndarray]:
        return -G * m * M / mag(r, axis=axis)
