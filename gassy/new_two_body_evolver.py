"""
[adapted from spiralling113MESA15]
"""
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np

from .constants import G, Msol, Rsol, pi, c
from .conversions import get_mu, get_period, mag

from gassy.stellar_profiles import get_smooth_profile_functions
from gassy.ode_driver import ode_driver


class TwoBodySystem:
    def __init__(self, m1, m2, a):
        self.m1 = m1
        self.m2 = m2
        self.a = a
        self.T = get_period(m1, m2, a)

    def evolve(self):
        pass

    def

    def plot_orbit(self):
        pass



if __name__ == '__main__':
    two_body_evolver =  TwoBodySystem(1, 1, 1)
    two_body_evolver.evolve()
    two_body_evolver.plot_orbit()