import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

from gassy.constants import Rsol

VEL_CMAP = "Blues_r"


def plot_energy(ke, gpe, t, ax=None, save_fname=None):
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    ax1, ax2 = ax, ax.twinx()

    tot = ke + gpe
    ax1.plot(t, ke / np.abs(ke[0]), label="Kinetic", color="C1")
    ax1.plot(t, gpe / np.abs(gpe[0]), label="Gravitational", color="C2")
    ax2.plot(t, (tot - tot[0]) / tot[0], label="Total", color="C3", linestyle="--")

    ax1.legend(frameon=True, loc="upper left")
    ax1.set_xlabel("t [s]")
    ax1.set_ylabel("$E/|E_0|$")

    ax2.set_ylabel("Total $\Delta E / E_0$", color="C3")
    ax2.tick_params(axis="y", labelcolor="C3")

    ax1.axhline(0, linestyle="-", color="black", linewidth=0.5)
    ax2.axhline(0, linestyle="-", color="red", linewidth=0.5)

    if save_fname is not None:
        plt.savefig(save_fname)
    return ax


def plot_angular_momentum(angular_momentum, save_fname=None):
    pass


def plot_waveform(h, t, ax=None, save_fname=None):
    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(t, h[0, :], label="+")
    ax.plot(t, h[1, :], label="x")
    ax.set_xlabel("t [s]")
    ax.set_ylabel(r"$h_{\rm GW}$")
    ax.legend(frameon=True, loc="upper left")
    if save_fname is not None:
        plt.savefig(save_fname)
    return ax


def plot_velocity(vel, t, ax=None, save_fname=None, escape_vel=None):
    if ax is None:
        fig, ax = plt.subplots()
    vel = get_mag(vel)
    rel_vel = vel / vel[0]
    ax.scatter(t, rel_vel, c=rel_vel, cmap=VEL_CMAP, s=0.75)
    ax.set_xlabel("t [s]")
    ax.set_ylabel(r"$|\vec{v}/\vec{v}_0|$")

    if escape_vel is not None:
        ax.axhline(escape_vel, linestyle="--", color="red", linewidth=0.5)
        ## anotate axhline with text at y=escape_vel
        ax.text(
            0.1,
            escape_vel,
            "escape velocity",
            color="red",
            fontsize=8,
            ha="center",
            va="bottom",
        )

    if save_fname is not None:
        plt.savefig(save_fname)
    return ax


def plot_orbit(pos, vel=[], ax=None, save_fname=None):
    if ax is None:
        fig, ax = plt.subplots()
    ax.set_facecolor("black")

    valid_idx = np.isfinite(pos[:, 0])
    pos = pos[valid_idx, :] / Rsol
    if len(vel) > 0:
        vel = get_mag(vel[valid_idx, :])
        rel_vel = vel / vel[0]
        ax.scatter(
            pos[valid_idx, 0], pos[valid_idx, 1], s=0.75, c=rel_vel, cmap=VEL_CMAP
        )
    else:
        ax.plot(pos[valid_idx, 0], pos[valid_idx, 1], linewidth=0.5, c="white")

    ax.scatter(0, 0, marker="*", color="red")
    ax.scatter(pos[0, 0], pos[0, 1], marker="o", color="white")

    # max pos
    rs = get_mag(pos)
    r = np.max(rs[np.isfinite(rs)]) * 1.1
    ax.set_xlim(-r, r)
    ax.set_ylim(-r, r)

    ax.set_xlabel(r"x [$R_{\odot}$]")
    ax.set_ylabel(r"y [$R_{\odot}$]")
    ax.set_aspect("equal", adjustable="box")

    if save_fname is not None:
        plt.savefig(save_fname)
    return ax


def plot_diagnostic(
    pos, ke, gpe, t, h=[], vel=[], save_fname="orbit_diagnostic.png", label=None
):
    fig = plt.figure(figsize=(12, 6))
    gs = GridSpec(3, 4)
    ax_orbit = fig.add_subplot(gs[:, 0:2])
    ax_energy = fig.add_subplot(gs[0, 2:4])

    plot_orbit(pos, vel, ax=ax_orbit)
    plot_energy(ke, gpe, t, ax=ax_energy)

    if len(vel) > 0:
        ax_vel = fig.add_subplot(gs[1, 2:4], sharex=ax_energy)
        plot_velocity(vel, t, ax=ax_vel)

    if len(h) > 0:
        ax_strain = fig.add_subplot(gs[2, 2:4], sharex=ax_energy)
        plot_waveform(h, t, ax=ax_strain)

    if label is not None:
        plt.suptitle(label)

    plt.tight_layout()
    plt.savefig(save_fname)


def get_mag(x: np.ndarray):
    """Get magnitude of a 2d vector"""
    return np.sqrt(x[:, 0] ** 2 + x[:, 1] ** 2)
