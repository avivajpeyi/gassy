from typing import Optional

import numpy as np

from gassy.constants import G, pi
from gassy.stellar_profiles.polytropic_star import PolytropicStar

from .two_body_base import TwoBodyBase


class TwoBodyInStellarProfile(TwoBodyBase):
    def __init__(
        self,
        m: float,
        M: float,
        r: float,
        init_vy: float,
        n: Optional[float] = None,
        mesa_profile_fname: Optional[str] = None,
        continue_on_error: bool = False,
    ):
        super().__init__(m, M, r, init_vy, continue_on_error)
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
