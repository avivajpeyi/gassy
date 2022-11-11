import os
import shutil
import unittest

from gassy.two_body import WaveformGenerator

DIR = os.path.abspath(os.path.dirname(__file__))

CLEANUP = False


class TestWaveformGenerator(unittest.TestCase):
    def setUp(self) -> None:
        self.kwgs = dict(m=1, M=1e4, init_x=-1, init_vy=0.025, num_periods=0.5)
        self.outdir = f"{DIR}/test_plots/waveform"
        os.makedirs(self.outdir, exist_ok=True)

    def tearDown(self) -> None:
        if CLEANUP:
            shutil.rmtree(self.outdir)

    def test_waveform_generator_no_drag(self):
        waveform = WaveformGenerator(**self.kwgs)
        waveform(distance=1000, save_plot_fname=f"{self.outdir}/{waveform.label}.png")

    def test_waveform_generator_small_drag(self):
        kwgs = self.kwgs.copy()
        kwgs["num_periods"] = kwgs["num_periods"] * 2
        waveform = WaveformGenerator(**kwgs, drag_coeff=1e-5)
        waveform(distance=1000, save_plot_fname=f"{self.outdir}/{waveform.label}.png")

    def test_waveform_generator_large_drag(self):
        waveform = WaveformGenerator(**self.kwgs, drag_coeff=1e-4)
        waveform(distance=1000, save_plot_fname=f"{self.outdir}/{waveform.label}.png")
