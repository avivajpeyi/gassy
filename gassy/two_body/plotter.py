import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec


def plot_energy(ke, gpe, t, ax=None, save_fname=None):
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    ax1, ax2 = ax, ax.twinx()

    tot = ke + gpe
    ax1.plot(t, ke / ke[0], label="Kinetic", color="C1")
    ax1.plot(t, gpe / gpe[0], label="Gravitational", color="C2")
    ax2.plot(t, tot / tot[0], label="Total", color="C3", linestyle="--")

    ax1.legend(frameon=True, loc="upper left")
    ax1.set_xlabel("t [s]")
    ax1.set_ylabel("E/E[0]")
    ax2.set_ylabel("TotE/TotE[0]", color="C3")
    ax2.tick_params(axis="y", labelcolor="C3")

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


def plot_orbit(pos, ax=None, save_fname=None):
    if ax is None:
        fig, ax = plt.subplots()

    valid_idx = np.isfinite(pos[:,0])
    ax.plot(pos[valid_idx, 0], pos[valid_idx, 1], "-o", markersize=0.1)
    ax.scatter(0, 0, marker="*", color="red")

    # max pos
    rs = np.sqrt(np.einsum("ij,ij->i", pos, pos))
    r = np.max(rs[np.isfinite(rs)])
    ax.set_xlim(-r, r)
    ax.set_ylim(-r, r)

    ax.set_xlabel(r"x [$R_{\odot}$]")
    ax.set_ylabel(r"y [$R_{\odot}$]")
    ax.set_aspect("equal", adjustable="box")

    if save_fname is not None:
        plt.savefig(save_fname)
    return ax


def plot_diagnostic(pos, ke, gpe, t, h=[], save_fname="orbit_diagnostic.png", label=None):
    fig = plt.figure(figsize=(12, 6))
    gs = GridSpec(2, 4)
    ax1 = fig.add_subplot(gs[0:2, 0:2])
    ax2 = fig.add_subplot(gs[0, 2:4])

    plot_orbit(pos, ax=ax1)
    plot_energy(ke, gpe, t, ax=ax2)
    if len(h) > 0:
        ax3 = fig.add_subplot(gs[1, 2:4])
        plot_waveform(h, t, ax=ax3)

    if label is not None:
        plt.suptitle(label)

    plt.tight_layout()
    plt.savefig(save_fname)
