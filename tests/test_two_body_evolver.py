import os
import unittest

from gassy.evolver.two_body_evolver import TwoBodyEvolver, TwoBodySystem

DIR = os.path.abspath(os.path.dirname(__file__))

import shutil

CLEANUP = False


class TestTwoBodyEvolver(unittest.TestCase):
    def setUp(self) -> None:
        self.body_kwgs = dict(m=1, M=1e4, init_x=-1, init_vy=0.025)
        self.evol_kwgs = dict(num_periods=3)
        self.outdir = f"{DIR}/test_plots/evolver"
        os.makedirs(self.outdir, exist_ok=True)

    def tearDown(self) -> None:
        if CLEANUP:
            shutil.rmtree(self.outdir)

    def test_no_drag_evol(self):
        two_body_system = TwoBodySystem(**self.body_kwgs)
        two_body_evolver = TwoBodyEvolver(two_body_system, **self.evol_kwgs)
        two_body_evolver.cache.plot(f"{self.outdir}/no_drag.png")

    def test_small_drag(self):
        two_body_system = TwoBodySystem(drag_coeff=0.0001, **self.body_kwgs)
        two_body_evolver = TwoBodyEvolver(two_body_system, **self.evol_kwgs)
        two_body_evolver.cache.plot(f"{self.outdir}/small_drag.png")

    def test_large_drag(self):
        two_body_system = TwoBodySystem(drag_coeff=0.00001, **self.body_kwgs)
        two_body_evolver = TwoBodyEvolver(two_body_system, **self.evol_kwgs)
        two_body_evolver.cache.plot(f"{self.outdir}/large_drag.png")
