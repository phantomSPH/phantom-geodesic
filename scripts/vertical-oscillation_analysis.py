"""
This script calculates the vertical oscillation frequency, and can compare it to the exact solution for a given radius.
Call it within the same directory as your positions.dat file, as well as the spin.value file.
It will output a file containing the raddii and the corresponding measured frequencies, as well as the exact values.
(omega-z.dat)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import periodogram
from oneDpowerspec import powerspec,peak_freq,plot_powerspec
from oscillation_frequencies import omega_z

data=np.loadtxt('positions.dat',skiprows=3)
a   =np.loadtxt('spin.value')
time=data[:,0]
s = np.shape(data)
n = int((s[1]-1)/3)
radius=data[:,1:s[1]:3]
theta =data[:,2:s[1]:3]

r0      = radius[0,:]
omegaz  = np.zeros(n)
for i in range(n):
    freqs,power=powerspec(time,theta[:,i])
    omegaz[i] = peak_freq(freqs,power)
    omegazapprox = omegaz[i]
    omegazexact  = omega_z(r0[i],a)
    print(omegazapprox,omegazexact,omegazapprox-omegazexact)


fs=30

plt.figure(figsize=(10,7.5))
plt.plot(r0,omegaz,'bx',label='Simulation')
plt.plot(r0,omega_z(r0,a),'g-',label='Exact')
plt.legend(loc='best')
# plt.plot(r_smooth,k_smooth,'k-')
plt.ylabel(r'$\Omega_z$',fontsize=fs)
plt.xlabel(r'$r$',fontsize=fs)
plt.title('a = '+str(a))
plt.savefig('python_plot.pdf')
plt.show()

df = freqs[1]-freqs[0]
print(df)
np.savetxt('omega-z.dat',np.column_stack((r0,omegaz,omega_z(r0,a))),header=str(a)+'\n'+str(df))

# plt.figure()
# plt.plot(r0,abs(omegaz-omega_z(r0,a)))
# plt.yscale('log')
# plt.show()
