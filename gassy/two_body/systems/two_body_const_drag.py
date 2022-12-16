import numpy as np

from .two_body_base import TwoBodyBase


class TwoBodyConstDrag(TwoBodyBase):
    def __init__(
        self,
        m: float,
        M: float,
        r: float,
        init_vy: float,
        drag_coeff: float,
        continue_on_error: bool = False,
    ):
        self.drag_coeff = drag_coeff
        super().__init__(m, M, r, init_vy, continue_on_error)

    @property
    def drag_force(self) -> np.ndarray:
        return -self.drag_coeff * self.vhat * self.vmag**2

    @property
    def label(self):
        return (
            f"Point particles: {self._m_label}, drag 'force'= -{self.drag_coeff}*(v^2)"
        )
