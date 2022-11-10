import matplotlib.pyplot as plt

from gassy.constants import Msol, Rsol
from gassy.two_body_evolver import evolve_bodies

R = 111.3642 * Rsol
X, Y, v_x, v_y, t = evolve_bodies(
    M=15 * Msol,
    m=1 * Msol,
    a=0.7 * R,
    e=0,
    Tend=3e6,
    mesa_profile_name="15Msol",
)
python_results = dict(X=X, Y=Y, v_x=v_x, v_y=v_y, t=t)

plt.figure(figsize=(5, 5))

plt.plot(python_results["X"], python_results["Y"], "o-", label="Python")
plt.legend()
plt.savefig(f"spiral_orbit.png")
