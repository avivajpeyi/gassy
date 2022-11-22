import numpy as np

from .two_body_base import TwoBodyBase


class TwoBodyConstDrag(TwoBodyBase):
    def __init__(
        self,
        m: float,
        M: float,
        init_x: float,
        init_vy: float,
        drag_coeff: float,
        continue_on_error: bool = False,
    ):
        self.drag_coeff = drag_coeff
        super().__init__(m, M, init_x, init_vy, continue_on_error)

    @property
    def drag_force(self) -> np.ndarray:
        # bondi-hoyle-lyttleton drag force
        return -self.drag_coeff * self.vhat * self.vmag**2

    @property
    def label(self):
        return (
            f"Point particles: {self._m_label}, drag 'force'= -{self.drag_coeff}*vhat"
        )
