"""
This script combines all the output files made by vertical-oscillation_analysis.py (i.e. omega-z.dat)
for different values of spin and plots them on the same graph.
Call it from the directory containing the directories of each simulation with different spin.
"""

from bash import BASH
import numpy as np
import matplotlib.pyplot as plt
from oscillation_frequencies import omega_z,isco

d_neg = BASH('ls -dr spin-neg*/').split('\n')
d_zer = BASH('ls -d spin-zero*/').split('\n')
d_pos = BASH('ls -d spin-pos*/').split('\n')
directories = d_neg+d_zer+d_pos
plt.style.use('paper')

# colormap = plt.cm.nipy_spectral #I suggest to use nipy_spectral, Set1,Paired
colormap = plt.cm.terrain
fig1 = plt.figure(1)
fig2 = plt.figure(2)
ax = fig1.add_subplot(111)
ax2= fig2.add_subplot(111)
ax.set_color_cycle([colormap(i) for i in np.linspace(0, 0.25,len(directories))])
ax2.set_color_cycle([colormap(i) for i in np.linspace(0, 0.25,len(directories))])

for d in directories:
    try:
        f = d+'omega-z.dat'
        print('Loading: '+f)
        r,omegaz,omegaz_exact=np.loadtxt(f,unpack=True)
        print('  |--> Success')
    except:
        print('  |--> Failed')
        continue
    i = 0
    for line in open(f):
        if i==0: a=float(line.strip("#").strip())
        if i==1:
            df=float(line.strip("#").strip())
            break
        i+=1
    plt.figure(1)
    # p = plt.plot(r,omegaz,'x')
    p = plt.errorbar(r,omegaz,yerr=df/2,fmt='.')
    print(df/2/omegaz_exact)
    c = p[0].get_color()
    # rfine=np.linspace(r[0],r[-1],250)
    r0 = isco(a)
    rmax = 13.
    rfine = np.linspace(r0,rmax,500)
    # plt.plot(r,omegaz_exact,label="a = "+str(a),color=c)
    plt.plot(rfine,omega_z(rfine,a),label="a = "+str(a),color=c,lw=1)
    plt.figure(2)
    plt.plot(r,np.abs(omega_z(r,a)-omegaz)/omega_z(r,a),'o-')

plt.yscale('log')
plt.figure(1)
# fs=25
plt.xlabel(r'$r$')
plt.ylabel(r'$\Omega_z$')
plt.xlim(xmax=rmax)
plt.ylim(ymin=0,ymax=0.18)
# plt.legend(frameon=False)

plt.savefig('python_plot.pdf',bbox_inches='tight')

plt.show()
