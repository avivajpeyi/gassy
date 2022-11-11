from typing import Optional

import numpy as np
from scipy.integrate import odeint

from .history import History
from .systems import TwoBodyBase


class Evolver:
    def __init__(
            self,
            two_body_system: TwoBodyBase,
            dt: Optional[float] = None,
            num_periods: Optional[float] = None,
            max_steps: Optional[int] = 10000,
    ):
        self.two_body = two_body_system
        T = self.two_body.period
        self.dt = T / 1000 if dt is None else dt
        self.N_T = 3.0 if num_periods is None else num_periods
        self.n_timesteps = int(self.N_T * T / self.dt)
        self.t = np.linspace(0, T * self.N_T, self.n_timesteps)
        if self.n_timesteps > max_steps:
            print(
                f"Too many steps requested {self.n_timesteps}. Reducing to first {max_steps} steps"
            )
        self.t = self.t[0:max_steps]
        print(f"{self.n_timesteps} evolutions!")
        self.history = self.evolve()

    def evolve(self):
        y0, t = self.two_body.pack_data(), self.t
        ode_out = odeint(
            self.dydt, y0=y0, t=t, args=(self.two_body,), atol=1e-12, rtol=1e-12
        )
        return History.from_ode_out(ode_out, t)

    def dydt(self, y, t, two_body_system: TwoBodyBase):
        pos, vel = y[:2], y[2:4]
        dpos_dt, dvel_dt = vel, two_body_system.accel
        extra_data = two_body_system.update(
            y
        )  # update the system's position and velocity for next step
        dydt = np.concatenate((dpos_dt, dvel_dt, extra_data[4:11]))
        return dydt
