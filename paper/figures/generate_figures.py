from gassy.constants import kpc
from gassy.two_body import WaveformGenerator

OUTDIR = "."


def plot_system_with_constant_drag(outdir):
    drag_amount = [0, 1e20]
    for drag in drag_amount:
        waveform_generator = WaveformGenerator.from_evol_inital_conditions(
            m=1, M=2, r=-1, num_periods=10, drag_coeff=drag
        )
        fig = waveform_generator.plot(distance=100 * kpc, save_dir=None)
        fig.tight_layout()
        fig.savefig(f"{outdir}/drag_{drag}.png")


def plot_common_envelope_example(outdir):
    waveform_generator = WaveformGenerator.from_evol_inital_conditions(
        m=1, M=10, r=-1, num_periods=10, stellar_polytropic_index=0.9
    )
    fig = waveform_generator.plot(distance=100 * kpc, save_dir=None)
    fig.tight_layout()
    fig.savefig(f"{outdir}/common_envelope_example.png")


def main(outdir=OUTDIR):
    # plot_system_with_constant_drag(outdir)
    plot_common_envelope_example(outdir)


if __name__ == "__main__":
    main()
