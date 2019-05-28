"""
This script calculates the epicyclic frequency, and can compare it to the exact solution for a given radius.
Call it within the same directory as your positions.dat file, as well as the spin.value file.
It will output a file containing the raddii and the corresponding measured frequencies, as well as the exact values.

"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import periodogram
from oneDpowerspec import powerspec,peak_freq,plot_powerspec
from oscillation_frequencies import kappa
from scipy.interpolate import interp1d

data=np.loadtxt('positions.dat',skiprows=3)
a   =np.loadtxt('spin.value')
time=data[:,0]
s = np.shape(data)
n = int((s[1]-1)/3)
radius=data[:,1:s[1]:3]

r0 = np.zeros(n)
k  = np.zeros(n)

t_fine = np.linspace(time[0],time[-1],len(time)*10)

for i in range(n):
    rad_i = radius[:,i]

    # half = int(len(rad_i)/2)
    # rad_i = rad_i[:half]

    freqs,power=powerspec(time,rad_i)
    # rad_i_interp = interp1d(time,rad_i)
    # freqs,power  = powerspec(t_fine,rad_i_interp(t_fine))
    r0[i] = rad_i[0]
    k[i]  = peak_freq(freqs,power)
    k_approx = k[i]
    k_exact  = kappa(rad_i[0],a)

    amp = (max(rad_i)-r0[i])/2
    rexact = -amp*np.cos(k_exact*time) + r0[i] + amp

    err = np.sum((rexact-rad_i)**2)
    err = err/(n*np.max(rad_i))
    err = np.sqrt(err)
    print('Simulation  kappa: ',k_approx)
    print('Exact       kappa: ',k_exact)
    print('kappa       error: ',k_approx-k_exact)
    print('kappa frac. error: ',(k_approx-k_exact)/k_exact)
    print('L2 theta(t) error: ',err)
    print('------------------------------------------------------------------------------------------------------')

# plt.show()

fs=30

plt.figure(figsize=(10,7.5))
plt.plot(r0,k,'bx',label='Simulation')
plt.plot(r0,kappa(r0,a),'g-',label='Exact')
plt.legend(loc='best')
# plt.plot(r_smooth,k_smooth,'k-')
plt.ylabel(r'$\kappa$',fontsize=fs)
plt.xlabel(r'$r$',fontsize=fs)
plt.title('a = '+str(a))
plt.savefig('python_plot.pdf')
# plt.show()

df = freqs[1]-freqs[0]
print(df)
np.savetxt('kappa.dat',np.column_stack((r0,k,kappa(r0,a))),header=str(a)+'\n'+str(df))
