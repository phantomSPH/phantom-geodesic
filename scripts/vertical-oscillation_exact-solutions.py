"""
Plot the vertical-oscillation frequencies as a function of radius,
for difference black hole spins. Starting from the ISCO (innermost stable circular orbits.)
Also calculates the number of radii to use for spins when simulating.
(This is done relative to spin a=1, which extends all the way in, so as to keep sampling density the same.)
"""

import numpy as np
import matplotlib.pyplot as plt
from oscillation_frequencies import omega_z

spin=np.linspace(-1,1,11)
risco=np.array([ 9.     ,  8.43176,  7.85069,  7.25427,  6.63904,  6. ,5.32944,  4.61434,  3.82907,  2.90664,  1.     ])

for i in range(len(spin)):
   a = spin[i]
   r = np.linspace(risco[i],14,500)
   plt.plot(r,omega_z(r,a),label=a)

plt.legend()
plt.show()

rmax = 14

n_radii_max = 50
n_radii=(rmax-risco)/(rmax-risco[-1])*n_radii_max

for i in range(len(n_radii)):
    n_radii[i] = int(n_radii[i])

print(n_radii)
