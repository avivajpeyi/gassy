from functools import partial

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate._ivp.base import \
    OdeSolver  # this is the class we will monkey patch
from tqdm import tqdm

### monkey patching the ode solvers with a progress bar

# save the old methods - we still need them
old_init = OdeSolver.__init__
old_step = OdeSolver.step


# define our own methods
def new_init(self, fun, t0, y0, t_bound, vectorized, support_complex=False):
    # define the progress bar
    self.pbar = tqdm(total=t_bound - t0, unit="ut", initial=t0, ascii=True, desc="IVP")
    self.last_t = t0

    # call the old method - we still want to do the old things too!
    old_init(self, fun, t0, y0, t_bound, vectorized, support_complex)


def new_step(self):
    # call the old method
    old_step(self)

    # update the bar
    tst = self.t - self.last_t
    self.pbar.update(tst)
    self.last_t = self.t

    # close the bar if the end is reached
    if self.t >= self.t_bound:
        self.pbar.close()


# overwrite the old methods with our customized ones
OdeSolver.__init__ = new_init
OdeSolver.step = new_step


# dydt function for the Lorenz Equations
def lorenz_dydt(_t, y, s=10, r=28, b=2.667):

    xp = s * (y[1] - y[0])
    yp = y[0] * (r - y[2]) - y[1]
    zp = y[0] * y[1] - b * y[2]

    return np.asarray([xp, yp, zp])


if __name__ == "__main__":

    # parameters
    sigma = 10
    rho = 28
    beta = 2.667

    # fix the parameters
    lorenz_dydt_ode = partial(lorenz_dydt, s=sigma, r=rho, b=beta)

    # solve the system
    tspan, y0 = [0, 4000], np.asarray([0.2, 0.3, 0.4])
    sol = solve_ivp(lorenz_dydt_ode, tspan, y0, method="RK45", rtol=1e-12)

    # optional plotting
    ax = plt.axes(projection="3d")
    ax.plot3D(*[column for column in sol.y], "blue")
    ax.set_title(rf"Lorenz Equations ($\sigma={sigma}, \rho={rho}, \beta={beta}$)")
    plt.show()
