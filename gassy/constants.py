"""constants in cgs"""
import numpy as np
from astropy.constants import G, M_sun, R_sun, au, c, kpc

G = G.cgs.value
c = c.cgs.value
Rsol = R_sun.cgs.value
Msol = M_sun.cgs.value
AU = au.cgs.value
pi = np.pi
kpc = kpc.cgs.value
