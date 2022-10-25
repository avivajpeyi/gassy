import unittest

import matplotlib.pyplot as plt
import numpy as np
import pytest
from deepdiff import DeepDiff
from scipy.io import loadmat
import os

from gassy.constants import G, Msol, Rsol, pi

from gassy.stellar_profiles import read_profile, get_smooth_profile_functions, StellarProfile
from gassy.two_body_evolver import (evolve_bodies,
                                    get_evolution_initial_conditions, twobody_in_gas_dydt_ode,
                                    get_N_steps,
                                    get_term_A, get_term_B, compute_drag_coef
                                    )



VISUAL_COMPARISONS = False

DIR = os.path.abspath(os.path.dirname(__file__))

fmt = lambda x: f"{x:.4e}"


class TestTwoBodyEvolver(unittest.TestCase):
    def setUp(self):
        self.M = 100
        self.m = 10
        self.a = 1
        self.e = 0
        self.Tend = 10
        self.profile_name = "15Msol"
        self.tol = 1e-6

    def test_inital_conditions(self):
        T, y0 = get_evolution_initial_conditions(
            self.M, self.m, self.a, self.e, self.profile_name
        )
        true_T = 2.3196e03
        true_y0 = [1, 0, 0, 44.4169]  # from matlab (lines 38-41, spiral113MESA15.m)
        diff = DeepDiff([T, y0], [true_T, true_y0])
        self.assertTrue(len(diff) == 1, f"Difference bw true and calc: {diff}")

    def test_odefun_get_N(self):
        N = get_N_steps(M=2.982825e+34, v=44, a=5428329882948, R=7754756975640, T=1.724996e+06)
        expected_N = 74.7268
        self.assertTrue(
            np.allclose(N, expected_N, rtol=self.tol),
            f"Expected N: {expected_N}, got: {N}"
        )

        N = get_N_steps(M=2.982825e+34, v=44, a=5428329882948, R=7754756975640, T=1.724996e+06)

    def test_get_terms(self):
        """
        fprintf("c_s=%d,v=%d,rdot=%d,M=%d,m=%d,r=%d,nu=%d", c_s, v, rdot, M, m, r, nu)
        """
        kwgs = dict(
            v=1.384621e+08, rdot=0,
            M=2.982825e+34, m=1.988550e+33, r=5428329882948,
            nu=5.859375e-02
        )
        A = get_term_A(**kwgs)
        expected_A = 2.3290e-05
        self.assertTrue(np.allclose(A, expected_A, rtol=self.tol), f"Expected A: {expected_A}, got: {A}")

        B = get_term_B(**kwgs)
        expected_B = 3.0792e-23
        self.assertTrue(np.allclose(B, expected_B, rtol=self.tol), f"Expected B: {expected_B}, got: {B}")

    def test_ode_step_zero_case(self):
        q, rho_func, c_s_func, M_e_func = get_smooth_profile_functions(self.profile_name)
        dydt = twobody_in_gas_dydt_ode(
            t=0,
            data=[0, 0, 0, 0],
            M=15,
            m=1,
            a=0.7,
            rho_func=rho_func,
            c_s_func=c_s_func,
            M_e_func=M_e_func,
        )
        self.assertTrue(np.allclose(dydt, [0, 0, 0, 0]))

    def test_ode_step(self):
        q, rho_func, c_s_func, M_e_func = get_smooth_profile_functions(self.profile_name)
        dydt = twobody_in_gas_dydt_ode(
            t=0,
            data=[1, 0, 0, 44],
            M=15,
            m=1,
            a=0.7,
            rho_func=rho_func,
            c_s_func=c_s_func,
            M_e_func=M_e_func,
        )
        expected_dydt = np.array([0, -74.8371, 44, -0.0187])  # obtained from matlab odefun(0, [1 0 0 44])
        diff = DeepDiff(expected_dydt, dydt)
        self.assertTrue(np.allclose(expected_dydt, dydt, rtol=1e-2), f"Difference bw true and calc: {diff}")

        dydt = twobody_in_gas_dydt_ode(
            t=0,
            data=[1, 0, 0, 6.15135515],
            M=15,
            m=1,
            a=0.7,
            rho_func=rho_func,
            c_s_func=c_s_func,
            M_e_func=M_e_func,
        )
        expected_dydt = [0,-74.8362,6.1514, -0.0742]
        diff = DeepDiff(expected_dydt, dydt)
        self.assertTrue(np.allclose(expected_dydt, dydt, rtol=1e-3), f"Difference bw true and calc: {diff}")

    @staticmethod
    def load_matlab_evolution_data():
        r = loadmat(f"{DIR}/data/spiral_out.mat")['out']
        fmt = lambda x: r[x][0][0].flatten()
        matlab_results = dict(X=fmt("X"), Y=fmt("Y"), v_x=fmt("v_x"), v_y=fmt("v_y"), t=fmt("t"))
        return matlab_results

    @pytest.mark.slow
    def test_evolver(self):
        """
        Runs the evolver and compares the results to the matlab version

        matlab code to generate test data:
        [X,Y,v_x,v_y,t] = spiralling113MESA15(100,10,1,0,0.5,10)
        out = struct('X', X, 'Y', Y, 'v_x', v_x, 'v_y', v_y, 't', t)
        save('spiral_out' out)
        """
        R = 111.3642 * Rsol
        X, Y, v_x, v_y, t = evolve_bodies(
            M=15 * Msol,
            m=1 * Msol,
            a=0.7 * R,
            e=0,
            Tend=1e10,
            mesa_profile_name=self.profile_name,
        )
        profile = StellarProfile.load_matlab_profile(self.profile_name)
        python_results = dict(X=X, Y=Y, v_x=v_x, v_y=v_y, t=t)

        fig, ax = plt.subplots(1,1, figsize=(5, 5))
        # profile.plot_grid(ax)
        ax.plot(python_results['X'], python_results['Y'], 'o-', label="Model")

        plt.savefig(f"{DIR}/test_plots/spiral_comparison.png")


    def test_drag_coeff(self):
        vs = np.linspace(0.01, 2)
        vec_coef = np.vectorize(compute_drag_coef)
        plt.plot(vs, vec_coef(v=vs,a=0.1,T=1,c_s=1,N=1000,M=1,rho=1))
        plt.plot(vs, vec_coef(v=vs, a=0.1, T=1, c_s=1, N=0.00001, M=1, rho=1))
        plt.plot(vs, vec_coef(v=vs, a=0.1, T=1000, c_s=1, N=0.00001, M=1, rho=1))
        plt.plot(vs, vec_coef(v=vs, a=0.1, T=1000, c_s=1, N=0.00001, M=15, rho=1))
        plt.plot(vs, vec_coef(v=vs, a=0.1, T=1000, c_s=1, N=0.00001, M=15, rho=0.00001))
        plt.savefig("coef.png")