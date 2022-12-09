import warnings
from functools import cached_property
from typing import List

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from scipy.optimize import root

from gassy.constants import G, Msol, Rsol

MAX_XI = 20


class PolytropicStar:
    """
    The main sequence profile for a Lane-emden polytropic star (given the mass, radius and n)

    Profile consists of the following:
     - central pressure, density, the constant of proportionality K and the scale length ($r_{n}$).

    The formulas listed are derived and discussed in this
    [book](https://link.springer.com/book/10.1007/978-1-4419-9110-2):
    $$P_c = \frac{8.952e+14}{(n+1)(\theta'_{n})^2_{\\xi_1}}\left(\frac{M}{M_\odot}\right)^2\left(\frac{R}{R_\odot}\right)^{-4}dyne.cm^{-2}$$
    $$K = \frac{G}{n+1}M^{1-1/n}R^{-1+3/n}\left(\frac{4\pi}{\\xi_1^{n+1}(-\theta'_n\Bigr|_{\\xi_1})^{n-1}}\right)^{1/n}$$
    $$\rho_c = \left(\frac{P_c}{K}\right)^{\frac{n}{n+1}} g/cc$$
    $$r_n = \sqrt{\frac{(n+1)P_c}{4\pi G\rho_c^2}}$$

    We can use the dimensionless definitions used to derive the Lane-Emden equation to determine the density and pressure profiles
    $$P = P_{c} \theta_{n}^{n+1}$$
    $$ \rho = \rho_{c}\theta_{n}^{n}$$

    To get the mass profile in terms of the polytrope, we separate variables in the mass continuity equation and
    substitute the polytropic properties where applicable
    $$\frac{dM}{dr} = 4\pi r^{2}\rho$$
    $$dM = 4\pi r^{2} \rho dr$$
    $$M(r) = 4\pi \int_{0}^{R} r^{2}\rho dr$$
    $$M(\\xi) = 4\pi r_n^3 \rho_c (-\\xi^2\theta'_n\Bigr|_\\xi)$$

    """

    def __init__(self, n, mass=1, radius=1):
        """
        n: polytropic index
        mass: in Solar Masses
        radius: in Solar radii
        """
        self.n = n
        self.xi_root, self.xi, self.theta = self.lane_emden_solver(n)
        if self.xi_root > MAX_XI:
            warnings.warn(
                f"xi_root {self.xi_root} is greater than {MAX_XI}. This may result in an invalid Profile."
            )
        self.d_theta = np.gradient(self.theta, self.xi)
        self.d_theta_xi_root = np.interp(self.xi_root, self.xi, self.d_theta)
        self.MaxM = mass
        self.MaxR = radius

    @cached_property
    def K(self) -> float:
        n, mass, radius = self.n, self.MaxM, self.MaxR
        term1 = G / (n + 1)
        term2 = np.power(mass * Msol, 1.0 - 1.0 / n, dtype=complex) * np.power(
            radius * Rsol, -1.0 + 3.0 / n, dtype=complex
        )
        term3 = 4 * np.pi
        term4 = np.power(self.xi_root, n + 1.0, dtype=complex) * np.power(
            self.d_theta_xi_root, n - 1.0, dtype=complex
        )
        K = term1 * term2 * np.power(term3 / term4, 1.0 / n, dtype=complex)
        return np.abs(K)

    @cached_property
    def P_c(self) -> float:
        n, mass, radius = self.n, self.MaxM, self.MaxR
        term1 = 8.952e14 * np.power(mass, 2.0) * np.power(radius, -4.0)
        term2 = (n + 1) * np.power(self.d_theta_xi_root, 2.0, dtype=complex)
        P_c = term1 / term2  # dyne/cm^2
        return np.abs(P_c)

    @cached_property
    def P(self) -> List[float]:
        P_c, theta, n = self.P_c, self.theta, self.n
        P = P_c * np.power(theta, n + 1, dtype=complex)
        return np.abs(P)

    @cached_property
    def rho_c(self) -> float:
        P_c, K, n = self.P_c, self.K, self.n
        rho_c = np.power(P_c / K, n / (n + 1.0), dtype=complex)
        return np.abs(rho_c)

    @cached_property
    def rho(self) -> np.ndarray:
        rho_c, theta, n = self.rho_c, self.theta, self.n
        rho = rho_c * np.power(theta, n, dtype=complex)
        return np.abs(rho)

    @cached_property
    def r_n(self) -> float:
        P_c, rho_c, n = self.P_c, self.rho_c, self.n
        term1 = (n + 1) * P_c
        term2 = 4 * np.pi * G * np.power(rho_c, 2.0, dtype=complex)
        r_n = np.sqrt(term1 / term2)
        return np.abs(r_n)

    @cached_property
    def r(self) -> np.ndarray:
        r_n, xi = self.r_n, self.xi
        r = r_n * xi
        return np.abs(r)

    @cached_property
    def m_e(self) -> np.ndarray:
        r_n, xi, rho_c = self.r_n, self.xi, self.rho_c
        dxi = np.interp(xi, xi, self.d_theta)
        m = (
            -4
            * np.pi
            * np.power(r_n, 3.0, dtype=complex)
            * rho_c
            * np.power(xi, 2.0, dtype=complex)
            * dxi
        )
        return np.abs(m)

    @cached_property
    def c_s(self) -> np.ndarray:
        gamma = 1.0 + 1.0 / self.n
        c_s = np.sqrt(gamma * self.P / self.rho)
        assert len(c_s[np.isnan(c_s)]) == 0, c_s
        return c_s

    def plot_profile(self, ax=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots(1, 5, figsize=(18, 4))

        ax[0].plot(self.r, self.rho)
        ax[0].set_xlabel(r"$\xi$")
        ax[0].set_ylabel(r"$\Theta(\xi)$")
        ax[0].plot(self.xi, self.theta, **kwargs)
        col = ax[0].lines[-1].get_color()
        ax[0].scatter(self.xi_root, 0, color=col)
        ax[0].set_xlim(0, MAX_XI)
        if self.xi_root > MAX_XI:
            ax[0].set_xlim(0.5, self.xi_root * 1.5)
            ax[0].set_xscale("log")
        ax[0].set_ylim(0, 1)

        ax[1].plot(self.r / Rsol, self.m_e / Msol, **kwargs)
        ax[1].set_ylabel(r"$M\ [M_{\odot}]$")
        ax[1].set_xlabel("$r\ [R_{\odot}]$")
        ax[1].set_xlim(0, self.MaxR)
        ax[1].set_ylim(0, self.MaxM)

        ax[2].plot(self.r / Rsol, self.P, **kwargs)
        ax[2].set_ylabel(r"$P\ [dyne/m^{2}]$ (?)")
        ax[2].set_xlabel("$r\ [R_{\odot}]$")
        ax[2].set_yscale("log")
        ax[2].set_xlim(0, self.MaxR)

        ax[3].plot(self.r / Rsol, self.rho, **kwargs)
        ax[3].set_ylabel(r"$\rho\ [kg/m^3]$ (?)")
        ax[3].set_xlabel("$r\ [R_{\odot}]$")
        ax[3].set_yscale("log")
        ax[3].set_xlim(0, self.MaxR)

        ax[4].plot(self.r / Rsol, self.c_s, **kwargs)
        ax[4].set_ylabel(r"$c_s\ [m/s]$ (?)")
        ax[4].set_xlabel("$r\ [R_{\odot}]$")
        ax[4].set_yscale("log")
        ax[4].set_xlim(0, self.MaxR)

        return ax

    @staticmethod
    def lane_emden_solver(n, numpts=5000, max_xi=MAX_XI):
        """
        Solve the Lane-Emden equation for a polytropic star

        In simple terms, the polytropic equation of state is a power-law relation between the
        pressure (P) and the density $(\rho)$. Algebraically, the relation can be written down as
        $$P = K\rho^{1+1/n}$$ where n is the polytropic index.
        The above equation of state is substituted in the equation of hydrostatic equilibrium
        $$\frac{dP}{dr}=-\frac{GM}{r^{2}}\rho$$ to arrive at the Lane-Emden equation
        $$\frac{1}{\\xi^2}\frac{d}{d\\xi}\left(\\xi^2\frac{d\theta_n}{d\\xi}\right) = -\theta_n^n$$
        Full derivation of the Lane-Emden equation is provided
        [here](https://github.com/jaadt7/Lane_Emden/blob/master/lane_emden_derivation.pdf).

        The $\theta_n(\\xi)$ is known as the Lane-Emden solution of index n, or polytrope of index n.
        The differential equation has a mixed boundary condition at the center where
        $$\theta_{n}(0)=1 ; \frac{d\theta_{n}}{d\\xi}\Bigr|_{\\xi=0}=0\ .$$

        Their physical significance reflects the existance of central densities and pressures,
        and no exchange to a "negative" radius.

        Another property of polytropes is that they admit roots for indices between 0 and 5; with 5 exclusive.
        Since $\\xi$ is defined as the dimensionless radius, the axis of the polytrope should only extend to the root
        since it would be the surface of the star.

        To numerically solve the Lane-Emden equation, we split it into 2 first order differential equations:
        $$\frac{d\theta}{d\\xi} = -\frac{\phi}{\\xi^2}$$
        $$\frac{d\phi}{d\\xi}=\theta^n\\xi^2$$ and feed it into the solver with the aforementioned boundary conditions.
        """
        dxi, theta0 = max_xi / numpts, np.array([1, 0])
        xi = np.arange(dxi, dxi * numpts, dxi)
        if n > 4:
            xi = np.geomspace(dxi, max_xi * 5, numpts)
        dtheta_dxi = lambda t, xi: np.real(
            np.array([-t[1] / (xi**2), np.power(t[0], n, dtype=complex) * (xi**2)])
        )
        theta = odeint(dtheta_dxi, theta0, xi)[:, 0]

        # finding the roots
        xi_root_guess = min(xi[theta < 0]) if np.any(theta < 0) else np.infty
        interp_theta = interp1d(
            xi, theta, kind="cubic", fill_value=np.infty, bounds_error=False
        )
        xi_root = root(lambda xi_: interp_theta(xi_), xi_root_guess).x[0]

        # truncate the solution to the root
        theta, xi = theta[xi <= xi_root], xi[xi <= xi_root]
        xi, theta = xi[~np.isnan(theta)], theta[~np.isnan(theta)]

        return (
            xi_root,
            xi,
            theta,
        )
