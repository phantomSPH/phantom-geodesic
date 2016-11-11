"""
This script combines all the output files made by vertical-oscillation_analysis.py (i.e. omega-z.dat)
for different values of spin and plots them on the same graph.
Call it from the directory containing the directories of each simulation with different spin.
"""

from bash import BASH
import numpy as np
import matplotlib.pyplot as plt
from oscillation_frequencies import omega_z

directories = BASH('ls -d spin-*/').split('\n')

colormap = plt.cm.nipy_spectral #I suggest to use nipy_spectral, Set1,Paired
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111)
ax.set_color_cycle([colormap(i) for i in np.linspace(0, 1,len(directories))])

for d in directories:
    try:
        f = d+'omega-z.dat'
        print('Loading: '+f)
        r,omegaz,omegaz_exact=np.loadtxt(f,unpack=True)
        rfine=np.linspace(r[0],r[-1],250)
        i = 0
        for line in open(f):
            if i==0: a=float(line.strip("#").strip())
            if i==1:
                df=float(line.strip("#").strip())
                break
            i+=1
        print('')
        # p = plt.plot(r,omegaz,'x')
        p = plt.errorbar(r,omegaz,yerr=df/2,fmt='.')
        c = p[0].get_color()
        # plt.plot(r,omegaz_exact,label="a = "+str(a),color=c)
        plt.plot(rfine,omega_z(rfine,a),label="a = "+str(a),color=c)
    except:
        print('  |--> Failed')

fs=25
plt.xlabel(r'$r$',fontsize=fs)
plt.ylabel(r'$\Omega_z$',fontsize=fs)
plt.legend()

plt.savefig('python_plot.pdf')

plt.show()
