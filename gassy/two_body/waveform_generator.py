"""
Calculates the gravitational waveform for a binary system with a given mass ratio, initial separation, and initial velocity.
"""

from typing import Optional

from . import Evolver, Strain, TwoBodyBase, create_two_body_system
from .plotter import plot_diagnostic


class WaveformGenerator:
    def __init__(
        self,
        m: float,
        M: float,
        r: float,
        num_periods: float,
        drag_coeff: Optional[float] = None,
        stellar_polytropic_index: Optional[float] = None,
        mesa_profile_fname: Optional[str] = None,
    ):
        kwgs = dict(
            m=m,
            M=M,
            init_x=r,
            drag_coeff=drag_coeff,
            stellar_polytropic_index=stellar_polytropic_index,
            mesa_profile_fname=mesa_profile_fname,
        )
        self.two_body_sys = create_two_body_system(**kwgs)
        self.history = Evolver(self.two_body_sys, num_periods=num_periods).history
        self.strain = Strain(self.history.mass_moment)
        self.label = self.two_body_sys.label

    def __call__(self, distance, theta=0, phi=0):
        h = self.strain.h(distance=distance, theta=theta, phi=phi)
        return self.history.time, h

    def plot(self, distance, theta=0, phi=0, save_dir=""):
        t, h = self(distance, theta, phi)
        save_fname = f"{save_dir}/{self.label}.png" if len(save_dir) > 0 else ""
        plot_diagnostic(
            pos=self.history.pos,
            ke=self.history.Ek,
            gpe=self.history.Egpe,
            t=self.history.time,
            vel=self.history.vel,
            h=h,
            label=self.two_body_sys.label,
            save_fname=save_fname,
        )


def orbital_velocity(m, M, r):
    return np.sqrt(G * (m + M) / r)
