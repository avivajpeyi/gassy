"""
Calculates the gravitational waveform for a binary system with a given mass ratio, initial separation, and initial velocity.
"""

import os
from typing import Optional

from . import Evolver, History, Strain, create_two_body_system
from .plotter import plot_diagnostic


class WaveformGenerator:
    def __init__(
        self,
        evolution_history: History,
        label: str,
    ):
        self.history = evolution_history
        self.strain = Strain(self.history.mass_moment)
        self.label = label

    @classmethod
    def from_evol_inital_conditions(
        cls,
        m: float,
        M: float,
        r: float,
        num_periods: float,
        drag_coeff: Optional[float] = None,
        stellar_polytropic_index: Optional[float] = None,
        mesa_profile_fname: Optional[str] = None,
        cache_dir: Optional[str] = ".",
    ):
        two_body_kwgs = dict(
            m=m,
            M=M,
            r=r,
            drag_coeff=drag_coeff,
            stellar_polytropic_index=stellar_polytropic_index,
            mesa_profile_fname=mesa_profile_fname,
        )
        two_body_sys = create_two_body_system(**two_body_kwgs)
        history = Evolver(two_body_sys, num_periods=num_periods).history
        if not os.path.isdir(cache_dir):
            os.makedirs(cache_dir)
        history.save(f"{cache_dir}/{two_body_sys.label}.npz")
        return cls(history, two_body_sys.label)

    @classmethod
    def from_cache(cls, cache_fname):
        history = History.load(cache_fname)
        label = os.path.basename(cache_fname).split(".")[-1]
        return cls(history, label)

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
            label=self.label,
            save_fname=save_fname,
        )
