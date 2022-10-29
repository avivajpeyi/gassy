import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

from gassy.evolver.ode_driver import ode_driver
from gassy.two_body.two_body_system import TwoBodySystem


class TwoBodyEvolver:
    def __init__(self, two_body_system, dt):
        self.two_body_system = two_body_system
        self.pos, self.vel = self.evolve(dt)

    def evolve(self, dt):
        y0 = self.two_body_system.get_pos_vel()
        t = np.linspace(0, two_body_system.period * 3, 1001)
        ode_out = ode_driver(self.dydt, y0=y0, t=t, args=(self.two_body_system,))
        pos_cache, vel_cache = ode_out[:, 0:2], ode_out[:, 2:4]
        return pos_cache, vel_cache

    def dydt(self, y, t, two_body_system: TwoBodySystem):
        pos, vel = y[:2], y[2:]
        new_y = np.concatenate((vel, two_body_system.net_acceleration))
        two_body_system.update(y)  # update the system's position and velocity for next step
        return new_y

    def plot_energy(self, ax=None, save_fname=None):
        if ax is None:
            fig, ax = plt.subplots(1, 1)
        ax1, ax2 = ax, ax.twinx()

        ke, gpe = self.two_body_system.get_energies(
            cached_pos=self.pos, cached_v=self.vel
        )
        ax1.plot(ke - np.mean(ke), label="Kinetic", color="C1")
        ax1.plot(gpe - np.mean(gpe), label="Gravitational", color="C2")
        tot = ke + gpe
        ax2.plot(tot - np.mean(tot), label="Total", color="C3", linestyle="--")
        ax1.legend(frameon=True, loc='upper right')
        ax1.set_xlabel("t [??]")
        ax1.set_ylabel("E - avg(E) [J?]")
        ax2.set_ylabel("Tot (E - avg(E)) [J?]", color="C3")
        ax2.tick_params(axis='y', labelcolor="C3")
        if save_fname is not None:
            plt.savefig(save_fname)
        return ax

    def plot_orbit(self, ax=None, save_fname=None):
        if ax is None:
            fig, ax = plt.subplots()
        ax.plot(self.pos[:, 0], self.pos[:, 1], '-o', markersize=0.1)
        ax.scatter(0, 0, marker="*", color="red")
        max_xy = self.two_body_system.init_a * 1.1
        # TODO: add orbit Type annotation on top left

        ax.set_xlim(-max_xy, max_xy)
        ax.set_ylim(-max_xy, max_xy)
        ax.set_xlabel("x [??]")
        ax.set_ylabel("y [??]")
        ax.set_aspect('equal', adjustable='box')
        if save_fname is not None:
            plt.savefig(save_fname)
        return ax

    def plot_diagnostic(self, save_fname="orbit_diagnostic.png"):
        fig = plt.figure( figsize=(6, 6))
        gs = GridSpec(3, 2)
        ax1 = fig.add_subplot(gs[0:2, 0:2])
        ax2 = fig.add_subplot(gs[2, 0:2])
        self.plot_orbit(ax=ax1)
        self.plot_energy(ax=ax2)

        plt.suptitle(f"{self.two_body_system.label}")
        plt.tight_layout()
        plt.savefig(save_fname)


if __name__ == '__main__':
    two_body_system = TwoBodySystem(m=1, M=1.5e7, init_x=1, init_vy=1)
    two_body_evolver = TwoBodyEvolver(two_body_system, 0.1)
    two_body_evolver.plot_diagnostic("no_drag.png")

    two_body_system = TwoBodySystem(m=1, M=1.5e7, init_x=1, init_vy=1, drag_coeff=0.01)
    two_body_evolver = TwoBodyEvolver(two_body_system, 0.1)
    two_body_evolver.plot_diagnostic("drag_0_1.png")

    two_body_system = TwoBodySystem(m=1, M=1.5e7, init_x=1, init_vy=1, drag_coeff=0.07)
    two_body_evolver = TwoBodyEvolver(two_body_system, 0.1)
    two_body_evolver.plot_diagnostic("drag_0_7.png")