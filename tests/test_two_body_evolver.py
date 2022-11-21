import os
import shutil
import unittest

from gassy.two_body import Evolver, create_two_body_system

DIR = os.path.abspath(os.path.dirname(__file__))


CLEANUP = False


class TestTwoBodyEvolver(unittest.TestCase):
    def setUp(self) -> None:
        self.body_kwgs = dict(m=1, M=1, init_x=-1)
        self.evol_kwgs = dict(num_periods=2, max_steps=10000)
        self.outdir = f"{DIR}/test_plots/evolver"
        os.makedirs(self.outdir, exist_ok=True)

    def tearDown(self) -> None:
        if CLEANUP:
            shutil.rmtree(self.outdir)

    def test_no_drag_evol(self):
        two_body_system = create_two_body_system(**self.body_kwgs)
        two_body_evolver = Evolver(two_body_system, **self.evol_kwgs)
        two_body_evolver.history.plot(f"{self.outdir}/no_drag.png")

    def test_small_drag(self):
        two_body_system = create_two_body_system(drag_coeff=1e-5, **self.body_kwgs)
        two_body_evolver = Evolver(two_body_system, **self.evol_kwgs)
        two_body_evolver.history.plot(f"{self.outdir}/small_drag.png")

    def test_large_drag(self):
        two_body_system = create_two_body_system(drag_coeff=1e-4, **self.body_kwgs)
        two_body_evolver = Evolver(two_body_system, **self.evol_kwgs)
        two_body_evolver.history.plot(f"{self.outdir}/large_drag.png")
