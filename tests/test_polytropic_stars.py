import os
import shutil
import unittest

import matplotlib.pyplot as plt
import numpy as np
import pytest

from gassy.stellar_profiles.polytropic_star import PolytropicStar

CLEANUP = False
DIR = os.path.dirname(os.path.abspath(__file__))


class TestPolytropicProfile(unittest.TestCase):
    def setUp(self):
        self.outdir = f"{DIR}/test_plots/polytropic"
        os.makedirs(self.outdir, exist_ok=True)

    def tearDown(self):
        if CLEANUP:
            if os.path.isdir(self.outdir):
                shutil.rmtree(self.outdir)

    def test_polytropic_stars_roots(self):
        star = PolytropicStar(n=3.5)
        self.assertFalse(np.isnan(star.xi_root))

    @pytest.mark.slow
    def test_variety_of_stars_plot(self):
        params = [dict(mass=1, radius=1), dict(mass=2, radius=1)]

        for p in params:
            fig, ax = plt.subplots(1, 5, figsize=(16, 4))
            # for i, n in enumerate(np.arange(1, 5, 0.25)):
            for i, n in enumerate([1, 1.5, 2, 2.5, 3, 3.5, 4]):
                profile = PolytropicStar(n, **p)
                profile.plot_profile(ax=ax, label=f"n={n:.2f}", color=f"C{i}")
            ax[0].legend(
                loc="upper right",
                frameon=False,
                borderaxespad=0.0,
                labelspacing=0.25,
                handlelength=0.5,
                fontsize="small",
            )
            plt.suptitle(f"Lane-Emden Polytropic Star ({p})")
            plt.tight_layout()
            plt.savefig(f"{self.outdir}/lane_emden_polytropic_stars_{p}.png")
