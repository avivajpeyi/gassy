import matplotlib.pyplot as plt
import numpy as np

from gassy.constants import Msol, Rsol


def plot_desnity_grid(profile, ax=None, fname=None, **kwargs):
    if ax == None:
        fig, ax = plt.subplots(figsize=(5, 5))
    y = profile.q
    p = np.linspace(0, 2 * np.pi, 50)
    R, P = np.meshgrid(y, p)
    X, Y = R * np.cos(P), R * np.sin(P)
    Z = profile.rho(R)
    ax.pcolormesh(X / Rsol, Y / Rsol, Z, cmap="hot", **kwargs)
    # plt.colorbar(c)
    ax.set_facecolor("black")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"x $[R_{\odot}]$")
    ax.set_ylabel(r"y $[R_{\odot}]$")
    if fname is not None:
        plt.tight_layout()
        fig.savefig(fname)
    return ax


def plot_1d_profile_data(profile, fname="1dprofile.png"):
    fig, ax = plt.subplots(3, 1, figsize=(5, 8), sharex=True)
    ax[0].plot(profile.q / Rsol, profile.c_s(profile.q), label="c_s")
    ax[0].set_ylabel("c_s [km/s]")
    ax[0].set_yscale("log")
    ax[1].plot(profile.q / Rsol, profile.rho(profile.q), label="rho")
    ax[1].set_ylabel(r"$\rho$ [g/cm$^3$]")
    ax[1].set_yscale("log")
    ax[2].plot(profile.q / Rsol, profile.M_e(profile.q) / Msol, label="M_e")
    ax[2].set_ylabel(r"$M_e [M_{\odot}]$")
    ax[2].set_xlabel(r"$R\ [R_{\odot}]$")
    ax[2].set_yscale("log")
    ax[2].set_xscale("log")
    # remove whitespace between subplots and save without cutting axes labels
    fig.subplots_adjust(hspace=0, wspace=0)
    plt.tight_layout()
    fig.savefig(fname)
