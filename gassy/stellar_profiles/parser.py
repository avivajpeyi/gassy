import glob
import os
import re

import pandas as pd
from scipy.io import loadmat

from gassy import constants

HERE = os.path.abspath(os.path.dirname(__file__))
RE = "/profile(.*?).mat"


def get_profile_paths():
    paths = glob.glob(f"{HERE}/data/*.mat")
    profiles = {}
    for p in paths:
        name = re.search(RE, p).group(1)
        profiles[name] = p
    return profiles


def read_profile(profile_name: str) -> pd.DataFrame:
    profile_paths = get_profile_paths()
    available_profiles = list(profile_paths.keys())
    if profile_name not in available_profiles:
        raise ValueError(
            f"The profile {profile_name} is not available."
            f"The available profiles are {available_profiles}."
        )

    profile = loadmat(profile_paths[profile_name])
    params = ["c_s", "q", "rho"]
    profile = {p: profile[p].flatten() for p in params}
    data = pd.DataFrame(profile)
    data["c_s"] = data["c_s"] * 100  # convert to cgs (i think?) #TODO: check with evgeni
    data["q"] = data['q'] * constants.Rsol  # convert to cgs
    return data


