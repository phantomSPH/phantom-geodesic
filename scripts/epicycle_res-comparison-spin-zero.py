"""
Plots a comparison of two simulations (different length of time) of the epicyclic frequency.
Requires the kappa.dat files that are made by epicyclic_analysis.py.
"""

import numpy as np
import matplotlib.pyplot as plt

r,kal,ke = np.loadtxt('spin-zero0.0/kappa.dat',unpack=True)
r,kah,k  = np.loadtxt('high-res_spin-zero0.0/kappa.dat',unpack=True)

fs = 25
plt.figure(figsize=(15,10))
plt.subplot(211)
plt.xlabel(r'$r$',fontsize=fs)
plt.ylabel(r'$\kappa$',fontsize=fs)
plt.plot(r,kal,'-x',label='Low res.')
plt.plot(r,kah,'-x',label='High res.')
plt.legend()

plt.subplot(212)
plt.xlabel(r'$r$',fontsize=fs)
plt.ylabel(r'Error ($| \kappa - \kappa_{\rm{exact}} |$)',fontsize=fs)
plt.plot(r,abs(kal-ke),'-x',label='Low res.')
plt.plot(r,abs(kah-ke),'-x',label='High res.')
plt.yscale('log')
plt.legend()

plt.savefig('python_plot.pdf')
plt.show()
