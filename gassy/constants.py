"""constants in cgs"""
import numpy as np
from astropy.constants import G, M_sun, R_sun, au, c

G = G.cgs.value
G = 6.67e-8 # ginat value
c = c.cgs.value
c = 29979245800.0 # ginat value
Rsol = R_sun.cgs.value
Rsol = 696342*1e5 # ginat value
Msol = M_sun.cgs.value
Msol = 1.98855*1e33 # ginat value
AU = au.cgs.value
pi = np.pi

