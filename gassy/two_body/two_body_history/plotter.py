import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec


def plot_energy(ke, gpe, ax=None, save_fname=None):
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    ax1, ax2 = ax, ax.twinx()

    tot = ke + gpe
    ax1.plot(ke / ke[0], label="Kinetic", color="C1")
    ax1.plot(gpe / gpe[0], label="Gravitational", color="C2")
    ax2.plot(tot / tot[0], label="Total", color="C3", linestyle="--")

    ax1.legend(frameon=True, loc="upper right")
    ax1.set_xlabel("t [??]")
    ax1.set_ylabel("E/E[0]")
    ax2.set_ylabel("TotE/TotE[0]", color="C3")
    ax2.tick_params(axis="y", labelcolor="C3")

    if save_fname is not None:
        plt.savefig(save_fname)
    return ax


def plot_angular_momentum(angular_momentum, save_fname=None):
    pass


def plot_orbit(pos, ax=None, save_fname=None):
    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(pos[:, 0], pos[:, 1], "-o", markersize=0.1)
    ax.scatter(0, 0, marker="*", color="red")

    # max pos
    r = np.max(np.sqrt(np.einsum("ij,ij->i", pos, pos)))
    ax.set_xlim(-r, r)
    ax.set_ylim(-r, r)

    ax.set_xlabel("x [??]")
    ax.set_ylabel("y [??]")
    ax.set_aspect("equal", adjustable="box")

    if save_fname is not None:
        plt.savefig(save_fname)
    return ax


def plot_diagnostic(pos, ke, gpe, save_fname="orbit_diagnostic.png", label=None):
    fig = plt.figure(figsize=(6, 6))
    gs = GridSpec(3, 2)
    ax1 = fig.add_subplot(gs[0:2, 0:2])
    ax2 = fig.add_subplot(gs[2, 0:2])
    plot_orbit(pos, ax=ax1)
    plot_energy(ke, gpe, ax=ax2)

    if label is not None:
        plt.suptitle(label)

    plt.tight_layout()
    plt.savefig(save_fname)
