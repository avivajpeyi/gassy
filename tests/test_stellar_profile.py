import os
import shutil
import unittest

import pandas as pd
import pytest

from gassy.stellar_profiles import StellarProfile, read_profile

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
    profile.plot_grid(fname=fname.replace("profile.png", "grid.png"))
