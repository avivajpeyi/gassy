from gassy.stellar_profiles import read_profile

import pandas as pd

def test_profile_parser():
    p = read_profile("5Msol")
    assert isinstance(p, pd.DataFrame)
