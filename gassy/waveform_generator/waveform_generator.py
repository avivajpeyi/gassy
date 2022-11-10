"""
Calculates the gravitational waveform for a binary system with a given mass ratio, initial separation, and initial velocity.
"""
import numpy as np
from numpy import cos, sin

from gassy.constants import G, c


class WaveformGenerator:
    def __init__(self, mass_moment):
        self.mass_moment = mass_moment
        self.M00 = self.mass_moment[:, 0]
        self.M01 = self.mass_moment[:, 1]
        self.M10 = self.mass_moment[:, 2]
        self.M11 = self.mass_moment[:, 3]

    def __h_cross(self, distance, phi, theta):
        """
        h_cross =  G/(D*c^4)*( ( M00 - M11).*sin(2*phi)*cos(i) + 2*M01.*cos(2*phi)*cos(i));
        """
        p, t, D = phi, theta, distance
        Gc4_D = G * c**4 / D
        M00, M11, M01 = self.M00, self.M11, self.M01
        return Gc4_D * (
            (M00 - M11) * sin(2 * phi) * cos(theta)
            + 2 * M01 * cos(2 * phi) * cos(theta)
        )

    def __h_plus(self, distance, phi, theta):
        """
        h_plus = G/(D*c^4)*(
           M00.*(cos(phi)^2 - sin(phi)^2*cos(i)^2) +
           M11.*(sin(phi)^2 - cos(phi)^2*cos(i)^2) - M01.*(sin(2*phi)*(1+cos(i)^2) );
        """
        p, t, D = phi, theta, distance
        Gc4_D = G * c**4 / D
        M00, M11, M01 = self.M00, self.M11, self.M01
        R00 = (cos(phi) ** 2) - (sin(phi) ** 2) * (cos(theta) ** 2)
        R11 = (sin(phi) ** 2) - (cos(phi) ** 2) * (cos(theta) ** 2)
        R01 = sin(2 * phi) * (1 + cos(theta) ** 2)
        return Gc4_D * (M00 * R00 + M11 * R11 - M01 * R01)

    def h(self, phi, theta, distance):
        return np.concatenate(
            self.__h_cross(phi, theta, distance), self.__h_plus(phi, theta, distance)
        )
