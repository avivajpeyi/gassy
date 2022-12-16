import warnings

from scipy.integrate import odeint, solve_ivp
from scipy.integrate._ivp.base import OdeSolver
from scipy.integrate.odepack import ODEintWarning
from tqdm import tqdm

warnings.filterwarnings("ignore", category=ODEintWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


class OdeSolverWithProgressBar(OdeSolver):
    """
    Monkeypatch the solve_ivp ODE solver to add a progress bar
    """

    def __init__(self, fun, t0, y0, t_bound, vectorized, support_complex=False):
        self.pbar = tqdm(
            total=t_bound - t0, unit="ut", initial=t0, ascii=True, desc="ODE"
        )
        self.last_t = t0

        # call the super init method to set up the rest of the class
        super().__init__(fun, t0, y0, t_bound, vectorized, support_complex)

    def step(self):
        # call the old method
        super().step()

        # update the bar
        tst = self.t - self.last_t
        self.pbar.update(tst)
        self.last_t = self.t

        # close the bar if the end is reached
        if self.t >= self.t_bound:
            self.pbar.close()


# monkey patch the ivp_solve ODE solver
OdeSolver = OdeSolverWithProgressBar


def ode_driver(
    fun,
    t,
    y0,
    args=(),
    ode_type="odeint",
    **kwargs,
):
    """Wrapper for ODE solvers"""
    kwargs["atol"] = kwargs.get("atol", 1e-12)
    kwargs["rtol"] = kwargs.get("rtol", 1e-12)

    if ode_type == "odeint":
        kwargs["tfirst"] = kwargs.get("tfirst", True)
        sol = odeint(fun, y0, t, args=args, **kwargs)
    elif ode_type == "ivp":  # this is much slower than odeint
        t_span = (t[0], t[-1])
        kwargs["method"] = kwargs.get("method", "LSODA")
        sol = solve_ivp(fun, t_span, y0, args=args, t_eval=t, **kwargs)
        sol = sol.y.T
    else:
        raise ValueError("Unknown ODE solver type")

    return sol
