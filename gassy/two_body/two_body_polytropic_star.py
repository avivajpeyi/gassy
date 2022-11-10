import numpy as np

from ..constants import G, pi
from ..stellar_profiles.polytropic_star import PolytropicStar
from .two_body_system import TwoBodySystem


class TwoBodyPolytriopicStar(TwoBodySystem):
    def __init__(self, m: float, M: float, init_x: float, init_vy: float, n: float):
        super().__init__(m, M, init_x, init_vy)
        self.stellar_profile = PolytropicStar(n=n, mass=M, radius=1)

    @property
    def drag_force(self) -> np.ndarray:
        """Bondi-Hoyle-Lyttleton drag force"""
        star = self.stellar_profile
        c_s, rho = star.c_s(self.r), star.rho(self.r)
        Fd = 4 * pi * G**2
        Fd *= self.m**2 * rho * self.vmag
        Fd /= np.pow(c_s**2 + self.vmag**2, 1.5)
        return -self.drag_coeff * self.v / self.vmag
