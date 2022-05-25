import numpy as np
from .constants import pi, G

def get_period(M,m,a):
    return np.sqrt((4*pi**2)*(a**3)/(G*(M+m)))

def get_mu(M,m):
    return M*m/(M+m)

def mag(x):
    return np.sqrt([xi**2 for xi in x])
