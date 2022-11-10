import numpy as np

from .plotter import plot_diagnostic


class TwoBodyHistory:
    """
    Stores the data for the two body system at each timestep.
    """

    def __init__(self, pos, vel, time, mass_moment, Ek, Egpe, L):
        self.pos = pos
        self.vel = vel
        self.time = time
        self.mass_moment = mass_moment
        self.Ek = Ek
        self.Egpe = Egpe
        self.L = L

    @staticmethod
    def __check_fname(fname):
        if not fname.endswith(".npz"):
            raise ValueError("File must be .npz")

    def save(self, fname):
        """Save cache to file"""
        self.__check_fname(fname)
        np.savez(
            fname,
            pos=self.pos,
            vel=self.vel,
            time=self.time,
            mass_moment=self.mass_moment,
            kinetic_energy=self.Ek,
            gravitational_energy=self.Egpe,
            angular_momentum=self.L,
        )

    @classmethod
    def load(cls, fname):
        """Load cache from file"""
        cls.__check_fname(fname)
        data = np.load(fname)
        return cls(
            pos=data["pos"],
            vel=data["vel"],
            time=data["time"],
            mass_moment=data["mass_moment"],
            Ek=data["kinetic_energy"],
            Egpe=data["gravitational_energy"],
            L=data["angular_momentum"],
        )

    @classmethod
    def from_ode_out(cls, y: np.ndarray, t: np.ndarray):
        assert y.shape == (len(t), 11)
        return cls(
            pos=y[:, 0:2],
            vel=y[:, 2:4],
            Ek=y[:, 4],
            Egpe=y[:, 5],
            L=y[:, 6],
            mass_moment=y[:, 7:11],
            time=t,
        )

    def plot(self, save_fname=""):
        plot_diagnostic(self.pos, self.Ek, self.Egpe, save_fname=save_fname)
