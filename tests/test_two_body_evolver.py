import os
import shutil
import unittest

import pytest

from gassy.two_body import Evolver, History, create_two_body_system

DIR = os.path.abspath(os.path.dirname(__file__))

CLEANUP = True


class TestTwoBodyEvolver(unittest.TestCase):
    def setUp(self) -> None:
        self.body_kwgs = dict(m=1, M=10, r=-1)
        self.evol_kwgs = dict(num_periods=0.1, max_steps=200)
        self.evol_kwgs_slow = dict(num_periods=10, max_steps=10000)
        self.outdir = f"{DIR}/test_plots/evolver"
        os.makedirs(self.outdir, exist_ok=True)

    def tearDown(self) -> None:
        if CLEANUP:
            shutil.rmtree(self.outdir)

    def test_evols(self):
        two_body_system = create_two_body_system(**self.body_kwgs)
        Evolver(two_body_system, **self.evol_kwgs)
        two_body_system = create_two_body_system(drag_coeff=1e-5, **self.body_kwgs)
        Evolver(two_body_system, **self.evol_kwgs)

    def test_history_save_and_load(self):
        two_body_system = create_two_body_system(**self.body_kwgs)
        hist = Evolver(two_body_system, **self.evol_kwgs).history
        fname = f"{self.outdir}/history.npz"
        hist.save(fname)
        loaded_hist = History.load(fname)
        self.assertTrue((hist.pos == loaded_hist.pos).all())

    @pytest.mark.slow
    def test_history_plots(self):
        global CLEANUP
        CLEANUP = False

        # zero_drag_system = create_two_body_system(**self.body_kwgs)
        # two_body_evolver = Evolver(zero_drag_system, **self.evol_kwgs_slow)
        # two_body_evolver.history.plot(f"{self.outdir}/zero_drag.png")

        drag_system = create_two_body_system(drag_coeff=1e20, **self.body_kwgs)
        two_body_evolver = Evolver(drag_system, **self.evol_kwgs_slow)
        two_body_evolver.history.plot(f"{self.outdir}/some_drag.png")
