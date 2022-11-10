from gassy.two_body.two_body_system import TwoBodySystem, OrbitType
import numpy as np


def test_force():
    bodies = TwoBodySystem(m=1, M=100, init_x=-1, init_vy=0)
    Fg = bodies.gravitational_force  # should be towards +ive x dir
    assert Fg[1] == 0
    assert Fg[0] > 0
    assert bodies.orbit_type == OrbitType.BOUND
    # give object a 'kick'
    bodies.update(np.array([-1, 0, 100, 0]))
    assert bodies.orbit_type == OrbitType.UNBOUND
    assert bodies.escape_vel < np.linalg.norm(bodies.v)
