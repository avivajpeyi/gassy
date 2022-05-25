import unittest

import matplotlib.pyplot as plt
import numpy as np
import pytest
from deepdiff import DeepDiff
from scipy.io import loadmat

from gassy.constants import G, Msol, Rsol, pi
from gassy.mat2py import smooth
from gassy.stellar_profiles import read_profile
from gassy.two_body_evolver import (evolve_bodies,
                                    get_evolution_initial_conditions, odefun,
                                    smooth_M_e, smooth_rho)

VISUAL_COMPARISONS = False


class TestStringMethods(unittest.TestCase):
    def setUp(self):
        self.M = 100
        self.m = 10
        self.a = 1
        self.e = 0
        self.dt = 0.5
        self.Tend = 10
        self.profile_name = "15Msol"

    @pytest.mark.skip("cant run atm")
    def test_evolver(self):
        """
        [X,Y,v_x,v_y,t] = spiralling113MESA15(100,10,1,0,0.5,10)
        out = struct('X', X, 'Y', Y, 'v_x', v_x, 'v_y', v_y, 't', t)
        save('spiral_out' out)
        """
        r = loadmat("tests/data/spiral_out.mat")
        matlab_results = dict(X=r[0], Y=r[1], v_x=r[2], v_y=r[3], t=r[4])
        X, Y, v_x, v_y, t = evolve_bodies(
            M=self.M,
            m=self.m,
            a=self.a,
            e=self.e,
            dt=self.dt,
            Tend=self.Tend,
            mesa_profile_name=self.profile_name,
        )
        python_results = dict(X=X, Y=Y, v_x=v_x, v_y=v_y, t=t)
        diffs = DeepDiff(matlab_results, python_results)
        diff = DeepDiff([T, y0], [true_T, true_y0])
        self.assertTrue(len(diff) == 1, f"Difference bw true and calc: {diff}")

    def test_smoothning(self):
        # test inputs
        r = 2
        profile = read_profile(self.profile_name)
        rho = profile.rho
        q = profile.q * Rsol  # conversion to cgs
        R = max(q)
        q = smooth(q)

        # test 1
        true_rho = 1.2861e03
        calc_rho = smooth_rho(r, R, rho, q)
        assert np.isclose(true_rho, calc_rho, rtol=1)
        if VISUAL_COMPARISONS:
            x = np.logspace(1, 13, 100)
            plt.plot(x, [smooth_rho(xi, R, rho, q) for xi in x])
            plt.xscale("log")
            plt.savefig("tests/visual_comparisons/smooth_rho/python.png")

        # test 2
        true_M_e = 5.3871e03  # from matlab
        calc_M_e = smooth_M_e(self.a, R, rho, q)
        print(f"{true_M_e}, {calc_M_e}")
        assert np.isclose(true_M_e, calc_M_e)

    def test_inital_conditions(self):
        T, y0 = get_evolution_initial_conditions(
            self.M, self.m, self.a, self.e, self.profile_name
        )
        true_T = 2.3196e03
        true_y0 = [1, 0, 0, 44.4169]
        diff = DeepDiff([T, y0], [true_T, true_y0])
        self.assertTrue(len(diff) == 1, f"Difference bw true and calc: {diff}")

    def test_ode_step(self):
        profile = read_profile(self.profile_name)
        rho = profile.rho
        q = profile.q * Rsol  # conversion to cgs
        R = max(q)
        q = smooth(q)
        dydx = odefun(
            t=0,
            data=[1, 0, 0, 44.4169],
            M=15,
            m=1,
            a=0.7 * R,
            orig_c_s=profile.c_s,
            q=q,
            orig_rho=rho,
            R=R,
        )
        print(dydx)
