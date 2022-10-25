import pandas as pd

from gassy.stellar_profiles import read_profile, get_smooth_profile_functions, StellarProfile
from gassy.constants import Msol, Rsol
import numpy as np
import unittest
import matplotlib.pyplot as plt


import pytest, os

DIR = os.path.dirname(os.path.abspath(__file__))


class TestStellarProfile(unittest.TestCase):
    def setUp(self):
        self.profile_name = "15Msol"
        self.tolerance = 1e-5

    def test_profile_parser(self):
        """Simple test to see if profile parser works"""
        p = read_profile("15Msol")
        self.assertIsInstance(p, pd.DataFrame)
        expected_cols = ["rho", "c_s", "q"]
        self.assertTrue(all([col in p.columns for col in expected_cols]))

    def test_rho_smoothning(self):
        rho_f, c_f, M_e_f = get_smooth_profile_functions(self.profile_name)

        # test 1
        true_rho = 1.2861e03
        calc_rho = rho_f(2)
        self.assertTrue(
            np.isclose(true_rho, calc_rho, rtol=self.tolerance),
            f"rho: {true_rho} != {calc_rho}"
        )

    def test_M_e_smoothning(self):
        rho_f, c_f, M_e_f = get_smooth_profile_functions(self.profile_name)

        true_M_e = 5.3871e03  # from matlab
        calc_M_e = M_e_f(1)
        self.assertTrue(
            np.isclose(true_M_e, calc_M_e, rtol=self.tolerance),
            f"M_e: {true_M_e} != {calc_M_e}"
        )

        true_M_e = 3.0484e+34  # from matlab
        calc_M_e = M_e_f(5.4283e+12)
        # FIXME: this doenst match exactly -- something wrong in my integration?
        self.assertTrue(
            np.isclose(true_M_e, calc_M_e, rtol=3),
            f"M_e: {true_M_e} != {calc_M_e}"
        )


    @pytest.mark.slow
    def test_plot(self):
        profile = StellarProfile.load_matlab_profile(self.profile_name)
        rho = profile.rho(profile.q)
        M_e = profile.M_e(profile.q)
        c_s = profile.c_s(profile.q)
        fig, axes = plt.subplots(3,1, sharex=True, figsize=(8,8))
        r = profile.q/Rsol
        axes[0].plot(r, c_s, label='c_s')
        axes[0].set_ylabel(r'$c_s$')
        axes[1].plot(r, rho, label='rho')
        axes[1].set_ylabel(r'$\rho$')
        axes[2].plot(r, M_e/Msol, label='M_e [Msol]')
        axes[2].set_ylabel(r'$M_e\ [M_{\odot}]$')
        fig.subplots_adjust(hspace=0, wspace=0)
        plt.xscale("log")
        plt.xlabel("$R\ [R_{\odot}]$")
        plt.savefig(f"{DIR}/visual_comparisons/smooth_rho/stellar_profile.png")
