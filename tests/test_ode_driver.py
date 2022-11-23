from unittest import TestCase

import numpy as np

from gassy.two_body.ode_driver import ode_driver


def fun(t, y):
    return -y


class TestOdeDriver(TestCase):
    def setUp(self) -> None:
        self.t = np.linspace(0, 10, 100)
        self.y0 = np.array([1])

    def test_odeint_driver(self):
        y = ode_driver(fun, self.t, self.y0, ode_type="odeint")
        self.assertEqual(len(y), len(self.t))

    def test_solve_ivp_driver(self):
        y = ode_driver(fun, self.t, self.y0, ode_type="ivp")
        self.assertEqual(len(y), len(self.t))
