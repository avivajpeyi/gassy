"""
Calculates the gravitational waveform for a binary system with a given mass ratio, initial separation, and initial velocity.
"""
import numpy as np


class WaveformGenerator:
    def __init__(self, , ):
        self.two_body_system = two_body_system

    def __h_cross(self, phi, theta, distance):
        pass

    def __h_plus(self, phi, theta, distance):
        pass

    def h(self, phi, theta, distance):
        return np.concatenate(
            self.__h_cross(phi, theta, distance),
            self.__h_plus(phi, theta, distance)
        )




