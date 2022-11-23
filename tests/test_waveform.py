import os
import shutil
import unittest

import pytest

from gassy.two_body import WaveformGenerator

DIR = os.path.abspath(os.path.dirname(__file__))

CLEANUP = True


class TestWaveformGenerator(unittest.TestCase):
    def setUp(self) -> None:
        self.kwgs = dict(m=1, M=20, r=-1, num_periods=0.1)
        self.slow_kwgs = dict(m=1, M=20, r=-1, num_periods=10)
        self.outdir = f"{DIR}/test_plots/waveform"
        os.makedirs(self.outdir, exist_ok=True)
        self.plot_kwrgs = dict(distance=1000, save_dir=self.outdir)

    def tearDown(self) -> None:
        if CLEANUP and os.path.exists(self.outdir):
            shutil.rmtree(self.outdir)

    def test_waveform_generator(self):
        t, h = WaveformGenerator(**self.kwgs)(distance=10, theta=0, phi=0)
        self.assertTrue(h.shape == (2, len(t)))

    @pytest.mark.slow
    def test_waveform_generator_plots(self):
        waveform_generator = WaveformGenerator(**self.slow_kwgs)
        waveform_generator.plot(**self.plot_kwrgs)

        waveform_generator = WaveformGenerator(**self.slow_kwgs, drag_coeff=1e-8)
        waveform_generator.plot(**self.plot_kwrgs)
