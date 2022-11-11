from .two_body_base import TwoBodyBase
import numpy as np


class TwoBodyConstDrag(TwoBodyBase):

    def __init__(self, m: float, M: float, init_x: float, init_vy: float, drag_coeff: float):
        self.drag_coeff = drag_coeff
        super().__init__(m, M, init_x, init_vy)


    @property
    def drag_force(self) -> np.ndarray:
        return -self.drag_coeff * self.vhat

    @property
    def label(self):
        return f"Point particles: {self._m_label}, drag 'force'= -{self.drag_coeff}*vhat"
