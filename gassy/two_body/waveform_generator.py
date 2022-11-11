"""
Calculates the gravitational waveform for a binary system with a given mass ratio, initial separation, and initial velocity.
"""

from . import TwoBodyBase, Evolver, Strain, create_two_body_system
from typing import Optional
from .plotter import plot_diagnostic


class WaveformGenerator:
    def __init__(self,
                 m: float, M: float, init_x: float, init_vy: float,
                 num_periods: float,
                 drag_coeff: Optional[float] = None, stellar_polytropic_index: Optional[float] = None,
                 mesa_profile_fname: Optional[str] = None):
        kwgs = dict(m=m, M=M, init_x=init_x, init_vy=init_vy, drag_coeff=drag_coeff,
                    stellar_polytropic_index=stellar_polytropic_index, mesa_profile_fname=mesa_profile_fname)
        self.two_body_sys = create_two_body_system(**kwgs)
        self.history = Evolver(self.two_body_sys, num_periods=num_periods).history
        self.strain = Strain(self.history.mass_moment)
        self.label = self.two_body_sys.label

    def __call__(self, distance, theta=0, phi=0, save_plot_fname=""):
        h = self.strain.h(distance=distance, theta=theta, phi=phi)
        if save_plot_fname:
            plot_diagnostic(
                pos=self.history.pos,
                ke=self.history.Ek,
                gpe=self.history.Egpe,
                t=self.history.time,
                h=h,
                label=self.two_body_sys.label,
                save_fname=save_plot_fname
            )
        return self.history.time, h
