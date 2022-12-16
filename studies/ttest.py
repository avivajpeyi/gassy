from gassy.two_body import WaveformGenerator

waveform_generator = WaveformGenerator.from_evol_inital_conditions(
    m=1, M=10, r=-100, num_periods=0.1, stellar_polytropic_index=1
)
waveform_generator.plot(distance=100)
