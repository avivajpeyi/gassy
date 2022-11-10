"""
Monkeypatching the ODE solver to add a progress bar

#TODO: doenst work for ODEInt
#TODO: add option for progress bar
"""
from scipy.integrate import odeint, solve_ivp
from scipy.integrate._ivp.base import \
    OdeSolver  # this is the class we will monkey patch
from tqdm import tqdm

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


def ode_driver(fun, y0, t, args=(), **kwargs):
    """
    Wrapper for scipy.integrate.solve_ivp that adds a progress bar.
    """
    if "rtol" not in kwargs:
        kwargs["rtol"] = 1e-12
    if "atol" not in kwargs:
        kwargs["atol"] = 1e-12
    return odeint(fun, y0, t, args=args, **kwargs)
