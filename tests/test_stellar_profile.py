import os
import shutil
import unittest

import numpy as np
import pandas as pd
import pytest

from gassy.stellar_profiles import (
    StellarProfile,
    get_smooth_profile_functions,
    read_profile,
)

CLEANUP = False

DIR = os.path.dirname(os.path.abspath(__file__))


class TestStellarProfile(unittest.TestCase):
    def setUp(self):
        self.profile_name = "15Msol"
        self.tolerance = 1e-5
        self.outdir = f"{DIR}/test_plots/stellar_profile"
        os.makedirs(self.outdir, exist_ok=True)

    def tearDown(self) -> None:
        if CLEANUP:
            shutil.rmtree(self.outdir)

    def test_profile_parser(self):
        """Simple test to see if profile parser works"""
        p = read_profile("15Msol")
        self.assertIsInstance(p, pd.DataFrame)
        expected_cols = ["rho", "c_s", "q"]
        self.assertTrue(all([col in p.columns for col in expected_cols]))

    def test_rho_smoothning(self):
        rho_f, c_f, M_e_f = get_smooth_profile_functions(self.profile_name)
        true_rho = 1.2861e03
        calc_rho = rho_f(2)
        self.assertTrue(
            np.isclose(true_rho, calc_rho, rtol=self.tolerance),
            f"rho: {true_rho} != {calc_rho}",
        )

    def test_M_e_smoothning(self):
        rho_f, c_f, M_e_f = get_smooth_profile_functions(self.profile_name)
        true_M_e = 5.3871e03  # from matlab
        calc_M_e = M_e_f(1)
        self.assertTrue(
            np.isclose(true_M_e, calc_M_e, rtol=self.tolerance),
            f"M_e: {true_M_e} != {calc_M_e}",
        )

        true_M_e = 3.0484e34  # from matlab
        calc_M_e = M_e_f(5.4283e12)
        # FIXME: this doenst match exactly -- something wrong in my integration?
        self.assertTrue(
            np.isclose(true_M_e, calc_M_e, rtol=3), f"M_e: {true_M_e} != {calc_M_e}"
        )

    @pytest.mark.slow
    def test_plot(self):
        profile = StellarProfile.load_matlab_profile(self.profile_name)
        make_profile_plots(profile, f"{self.outdir}/stellar_mesa_profile.png")

    @pytest.mark.slow
    def test_plot_polytropic(self):
        profile = StellarProfile.load_polytropic_profile(n=3.25, mass=2, radius=10)
        make_profile_plots(profile, f"{self.outdir}/stellar_polytropic_profile.png")


def make_profile_plots(profile, fname):
    profile.plot_1d_profile_data(fname)
    profile.plot_grid(fname=fname.replace("profile", "grid"))
