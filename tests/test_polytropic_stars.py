from gassy.stellar_profiles.polytropic_star import PolytropicStar
import numpy as np


def test_polytropic_stars_roots():
    star = PolytropicStar(n=3.5)
    assert not np.isnan(star.xi_root)