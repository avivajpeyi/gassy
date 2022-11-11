import unittest

import numpy as np

from gassy.two_body import OrbitType, TwoBodyBase


class TestTwoBodyBasics(unittest.TestCase):
    def test_two_body_config(self):
        bodies = TwoBodyBase(m=1, M=100, init_x=-1, init_vy=0)
        Fg = bodies.gravitational_force
        self.assertEqual(Fg[1], 0)  # no y component for force
        self.assertEqual(bodies.orbit_type, OrbitType.BOUND)
        self.assertTrue(Fg[0] > 0)  # should be towards +ive x dir

        # give object a 'kick'
        bodies.update(np.array([-1, 0, 100, 0]))
        self.assertEqual(bodies.orbit_type, OrbitType.UNBOUND)
        self.assertTrue(bodies.escape_vel < np.linalg.norm(bodies.v))
        with self.assertWarns(Warning):
            bodies.bound_orbit_check(continue_on_error=True)

    def test_two_body_raise_error(self):
        with self.assertRaises(ValueError):
            TwoBodyBase(m=1, M=100, init_x=-1, init_vy=1000)
