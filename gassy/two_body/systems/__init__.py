from typing import Optional

from .orbit_type import OrbitType
from .two_body_base import TwoBodyBase
from .two_body_const_drag import TwoBodyConstDrag
from .two_body_in_stellar_profile import TwoBodyInStellarProfile


def create_two_body_system(
    m: float,
    M: float,
    init_x: float,
    init_vy: Optional[float] = None,
    drag_coeff: Optional[float] = None,
    stellar_polytropic_index: Optional[float] = None,
    mesa_profile_fname: Optional[str] = None,
    continue_on_error: Optional[bool] = False,
) -> TwoBodyBase:
    kwargs = dict(
        m=m, M=M, init_x=init_x, init_vy=init_vy, continue_on_error=continue_on_error
    )
    if drag_coeff is not None:
        return TwoBodyConstDrag(drag_coeff=drag_coeff, **kwargs)
    elif mesa_profile_fname is not None:
        return TwoBodyInStellarProfile(mesa_profile_fname=mesa_profile_fname, **kwargs)
    elif stellar_polytropic_index is not None:
        return TwoBodyInStellarProfile(n=stellar_polytropic_index, **kwargs)
    else:
        return TwoBodyBase(**kwargs)
