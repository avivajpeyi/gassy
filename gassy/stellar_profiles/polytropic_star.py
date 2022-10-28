from scipy.integrate import odeint
from scipy.interpolate import interp1d
from scipy.optimize import root
import numpy as np
import matplotlib.pyplot as plt
from gassy.constants import Rsol, G, Msol
import warnings

MAX_XI = 20


class PolytropicStar:
    def __init__(self, n, mass=1, radius=1):
        self.n = n
        self.xi, self.theta = self.lane_emden_solver(n)
        self.xi_root, self.theta_root = self.xi[-1], self.theta[-1]
        if self.xi_root > MAX_XI:
            warnings.warn(f"xi_root {self.xi_root} is greater than {MAX_XI}. This may result in an invalid Profile.")
        self.d_theta = np.gradient(self.theta, self.xi)
        self.d_theta_xi_root = np.interp(self.xi_root, self.xi, self.d_theta)
        self.MaxM = mass
        self.MaxR = radius
        self.r, self.rho, self.P, self.m = self.get_main_sequence_profile(mass, radius)

    def _calc_K(self, mass, radius):
        n = self.n
        term1 = G / (n + 1)
        term2 = np.power(mass * Msol, 1. - 1. / n, dtype=complex) * np.power(radius * Rsol, -1. + 3. / n, dtype=complex)
        term3 = 4 * np.pi
        term4 = np.power(self.xi_root, n + 1., dtype=complex) * np.power(self.d_theta_xi_root, n - 1., dtype=complex)
        K = term1 * term2 * np.power(term3 / term4, 1. / n, dtype=complex)
        return K

    def get_main_sequence_profile(self, mass, radius):
        """
        Given a mass and radius, calculate the main sequence profile for a polytropic star.

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

        K = self._calc_K(mass, radius)
        n, xi, theta = self.n, self.xi, self.theta

        # pressure
        term1 = 8.952e+14 * np.power(mass, 2.) * np.power(radius, -4.)
        term2 = (n + 1) * np.power(self.d_theta_xi_root, 2., dtype=complex)
        P_c = term1 / term2  # dyne/cm^2
        P = P_c * np.power(theta, n + 1, dtype=complex)

        # density
        rho_c = np.power(P_c / K, n / (n + 1.), dtype=complex)
        rho = rho_c * np.power(theta, n, dtype=complex)

        # length scale and radius
        term1 = (n + 1) * P_c
        term2 = 4 * np.pi * G * np.power(rho_c, 2., dtype=complex)
        r_n = np.sqrt(term1 / term2)
        r = r_n * xi

        # calculating the mass profile
        dxi = np.interp(xi, xi, self.d_theta)
        m = (-4 * np.pi *
             np.power(r_n, 3., dtype=complex) *
             rho_c * np.power(xi, 2., dtype=complex) *
             dxi
             )

        return np.abs(r), np.abs(rho), np.abs(P), np.abs(m)


    def c_s(self):
        gamma = 1. + 1. / self.n
        return np.sqrt(gamma * self.P / self.rho)




    def plot_profile(self, ax=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots(1, 5, figsize=(18, 4))

        ax[0].plot(self.r, self.rho)
        ax[0].set_xlabel(r'$\xi$')
        ax[0].set_ylabel(r'$\Theta(\xi)$')
        ax[0].plot(self.xi, self.theta, **kwargs)
        col = ax[0].lines[-1].get_color()
        ax[0].scatter(self.xi_root, self.theta_root, color=col)
        ax[0].set_xlim(0, MAX_XI)
        if self.xi_root > MAX_XI:
            ax[0].set_xlim(0.5, self.xi_root * 1.5)
            ax[0].set_xscale('log')
        ax[0].set_ylim(0, 1)

        ax[1].plot(self.r / Rsol, self.m / Msol, **kwargs)
        ax[1].set_ylabel(r'$M\ [M_{\odot}]$')
        ax[1].set_xlabel("$r\ [R_{\odot}]$")
        ax[1].set_xlim(0, self.MaxR)
        ax[1].set_ylim(0, self.MaxM)

        ax[2].plot(self.r / Rsol, self.P, **kwargs)
        ax[2].set_ylabel(r'$P\ [dyne/cm^{2}]$')
        ax[2].set_xlabel("$r\ [R_{\odot}]$")
        ax[2].set_yscale('log')
        ax[2].set_xlim(0, self.MaxR)

        ax[3].plot(self.r / Rsol, self.rho, **kwargs)
        ax[3].set_ylabel(r'$\rho\ [g/cm^3]$')
        ax[3].set_xlabel("$r\ [R_{\odot}]$")
        ax[3].set_yscale('log')
        ax[3].set_xlim(0, self.MaxR)

        ax[4].plot(self.r / Rsol, self.c_s(), **kwargs)
        ax[4].set_ylabel(r'$c_s\ [cm/s]$')
        ax[4].set_xlabel("$r\ [R_{\odot}]$")
        ax[4].set_yscale('log')
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
        dtheta_dxi = lambda t, xi: np.real(np.array([
            -t[1] / (xi ** 2),
            np.power(t[0], n, dtype=complex) * (xi ** 2)
        ]))
        theta = odeint(dtheta_dxi, theta0, xi)[:, 0]

        # finding the roots
        xi_root_guess = min(xi[theta < 0]) if np.any(theta < 0) else np.infty
        interp_theta = interp1d(xi, theta, kind='cubic', fill_value=np.infty, bounds_error=False)
        xi_root = root(lambda xi_: interp_theta(xi_), xi_root_guess).x[0]

        # truncate the solution to the root
        theta, xi = theta[xi <= xi_root], xi[xi <= xi_root]
        theta, xi = np.append(theta, 0), np.append(xi, xi_root)
        xi, theta = xi[~np.isnan(theta)], theta[~np.isnan(theta)]

        return xi, theta


if __name__ == '__main__':
    fig, ax = plt.subplots(1, 5, figsize=(16, 4))
    for i, n in enumerate(np.arange(1, 5, 0.25)):
        profile = PolytropicStar(n, 1, 1)
        profile.plot_profile(ax=ax, label=f'n={n:.2f}', color=f'C{i}')
    ax[0].legend(loc='upper right', frameon=False, borderaxespad=0., labelspacing=0.25, handlelength=0.5,
                 fontsize='small')

    plt.suptitle("Lane-Emden Polytropic Star")
    plt.tight_layout()
    plt.savefig("profiles.png")
