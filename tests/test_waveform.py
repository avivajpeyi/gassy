import os
import shutil
import unittest

import pytest

from gassy.two_body import WaveformGenerator, create_two_body_system

DIR = os.path.abspath(os.path.dirname(__file__))

CLEANUP = True


class TestWaveformGenerator(unittest.TestCase):
    def setUp(self) -> None:
        self.outdir = f"{DIR}/test_plots/waveform"
        os.makedirs(self.outdir, exist_ok=True)
        self.plot_kwrgs = dict(distance=1000, save_dir=self.outdir)
        self.kwgs = dict(m=1, M=10, r=-1, num_periods=0.1, cache_dir=self.outdir)
        self.slow_kwgs = dict(m=1, M=10, r=-1, num_periods=10)

    def tearDown(self) -> None:
        if CLEANUP and os.path.exists(self.outdir):
            shutil.rmtree(self.outdir)

    def test_waveform_generator(self):
        two_body_system = create_two_body_system(
            m=self.kwgs["m"], M=self.kwgs["M"], r=self.kwgs["r"]
        )
        t, h = WaveformGenerator.from_evol_inital_conditions(**self.kwgs)(
            distance=10, theta=0, phi=0
        )
        self.assertTrue(h.shape == (2, len(t)))
        cache_fname = os.path.join(self.outdir, f"{two_body_system.label}.npz")
        t1, h1 = WaveformGenerator.from_cache(cache_fname)(distance=10, theta=0, phi=0)
        self.assertTrue((t == t1).all())
        self.assertTrue((h == h1).all())

    @pytest.mark.slow
    def test_waveform_generator_plots(self):
        global CLEANUP
        CLEANUP = False

        waveform_generator = WaveformGenerator.from_evol_inital_conditions(
            **self.slow_kwgs
        )
        waveform_generator.plot(**self.plot_kwrgs)

        waveform_generator = WaveformGenerator.from_evol_inital_conditions(
            **self.slow_kwgs, drag_coeff=1e20
        )
        waveform_generator.plot(**self.plot_kwrgs)
