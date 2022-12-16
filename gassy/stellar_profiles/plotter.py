import warnings

import matplotlib.pyplot as plt
import numpy as np

from gassy.constants import M_sun, R_sun

warnings.simplefilter("ignore", UserWarning)


def plot_desnity_grid(profile, ax=None, fname=None, **kwargs):
    if ax == None:
        fig, ax = plt.subplots(figsize=(5, 5))
    y = profile.q
    p = np.linspace(0, 2 * np.pi, 50)
    R, P = np.meshgrid(y, p)
    X, Y = R * np.cos(P), R * np.sin(P)
    Z = profile.rho(R)
    ax.pcolormesh(X / R_sun, Y / R_sun, Z, cmap="hot", **kwargs)
    # plt.colorbar(c)
    ax.set_facecolor("black")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"x $[R_{\odot}]$")
    ax.set_ylabel(r"y $[R_{\odot}]$")
    if fname is not None:
        plt.tight_layout()
        fig.savefig(fname)
    return ax


def plot_1d_profile_data(profile, fname="1dprofile.png", ax=[], save=True):
    if len(ax) == 0:
        fig, axes = plt.subplots(1, 1, figsize=(5, 3))
        ax = [axes, axes.twinx(), axes.twinx()]
        ax[2].tick_params(axis="y", direction="in", labelleft="on", pad=-22)
        ax[2].yaxis.set_label_position("left")
        ax[2].yaxis.set_ticks_position("left")
    for i, a in enumerate(ax):
        a.tick_params(axis="y", labelcolor=f"C{i}", colors=f"C{i}")

    ax[0].semilogx(profile.q / R_sun, profile.c_s(profile.q), label="c_s", color="C0")
    ax[1].semilogx(
        profile.q / R_sun,
        profile.rho(profile.q),
        label="rho",
        color="C1",
        ls="--",
        alpha=0.3,
    )
    ax[2].semilogx(
        profile.q / R_sun,
        profile.M_e(profile.q) / M_sun,
        label="M_e",
        color="C2",
        alpha=0.3,
    )
    ax[0].set_xlabel(r"$R\ [R_{\odot}]$")

    ax[0].set_ylabel("$c_s$ [km/s]", color="C0")
    ax[1].set_ylabel(r"$\rho$ [kg/m$^3$]", color="C1")
    ax[2].set_ylabel(r"$M_e [M_{\odot}]$", color="C2", labelpad=-35)
    ax[0].text(
        0.825,
        0.925,
        f"${profile.MaxM / M_sun:.2f} M_{{\\odot}}$",
        transform=ax[0].transAxes,
        color="C2",
        zorder=100,
    )

    if save:
        fig = ax[0].get_figure()
        plt.tight_layout()
        fig.savefig(fname)
