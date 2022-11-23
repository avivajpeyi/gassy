from typing import Optional

import numpy as np

from .history import History
from .ode_driver import ode_driver
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
        self.history = self.evolve()

    def evolve(self) -> History:
        evolved_pos_vel = ode_driver(
            fun=self.dydt,
            y0=self.two_body.pack_data(),
            t=self.t,
            args=(self.two_body,),
            # type='odeint'
        )
        return self._generate_two_body_history(evolved_pos_vel)

    def _generate_two_body_history(self, pos_vel: np.ndarray) -> History:
        num_points = len(pos_vel)
        two_body_data = np.zeros((num_points, self.two_body.data_dim))
        for i in range(num_points):
            self.two_body.update(pos_vel[i])
            two_body_data[i] = self.two_body.pack_data(all_data=True)
        return History.from_ode_out(two_body_data, self.t)

    def dydt(self, t, y, two_body_system: TwoBodyBase) -> np.ndarray:
        pos, vel = y[0:2], y[2:4]
        dpos_dt, dvel_dt = vel, two_body_system.accel
        two_body_system.update(
            y
        )  # update the system's position and velocity for next step
        dydt = np.concatenate((dpos_dt, dvel_dt))
        return dydt
