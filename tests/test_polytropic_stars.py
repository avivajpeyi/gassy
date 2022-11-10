import os
import shutil
import unittest

import matplotlib.pyplot as plt
import numpy as np

from gassy.stellar_profiles.polytropic_star import PolytropicStar

CLEANUP = False
DIR = os.path.dirname(os.path.abspath(__file__))


class TestPolytropicProfile(unittest.TestCase):
    def setUp(self):
        self.outdir = f"{DIR}/test_plots/polytropic"
        os.makedirs(self.outdir, exist_ok=True)

    def test_polytropic_stars_roots(self):
        star = PolytropicStar(n=3.5)
        self.assertFalse(np.isnan(star.xi_root))

    def test_variety_of_stars_plot(self):
        fig, ax = plt.subplots(1, 5, figsize=(16, 4))
        for i, n in enumerate(np.arange(1, 5, 0.25)):
            profile = PolytropicStar(n, 1, 1)
            profile.plot_profile(ax=ax, label=f"n={n:.2f}", color=f"C{i}")
        ax[0].legend(
            loc="upper right",
            frameon=False,
            borderaxespad=0.0,
            labelspacing=0.25,
            handlelength=0.5,
            fontsize="small",
        )
        plt.suptitle("Lane-Emden Polytropic Star")
        plt.tight_layout()
        plt.savefig(f"{self.outdir}/lane_emden_polytropic_stars.png")
