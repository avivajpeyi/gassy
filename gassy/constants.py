"""constants in SI"""
import numpy as np
from astropy import units
from astropy.constants import G, M_sun, R_sun, au, c, kpc

m_per_s = units.m / units.s
m_per_s2 = units.m / units.s**2
kgm_s2 = units.kg * m_per_s2
G = G.si.value
c = c.si.value
R_sun = R_sun.si.value
M_sun = M_sun.si.value
AU = au.si.value
kpc = kpc.si.value
pi = np.pi
