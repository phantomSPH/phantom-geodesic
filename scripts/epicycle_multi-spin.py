"""
This script combines all the output files made by epicycle_analysis.py (i.e. kappa.dat)
for different values of spin and plots them on the same graph.
Call it from the directory containing the directories of each simulation with different spin.
"""

from bash import BASH
import numpy as np
import matplotlib.pyplot as plt
from oscillation_frequencies import kappa, kappa_r0

d_neg = BASH('ls -dr spin-neg*/').split('\n')
d_zer = BASH('ls -d spin-zero*/').split('\n')
d_pos = BASH('ls -d spin-pos*/').split('\n')
directories = d_neg+d_zer+d_pos

fig_width_pt  = 240*3                       # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27                   # Convert pt to inches
golden_mean   = (np.sqrt(5)-1.0)/2.0        # Aesthetic ratio
fig_width     = fig_width_pt*inches_per_pt  # width in inches
fig_height    = fig_width*golden_mean       # height in inches
fig_size      = [fig_width,fig_width]
# fig_size = [3.3208800332088004, 2.0524167330839185]
plt.style.use('paper')

# colormap = plt.cm.nipy_spectral #I suggest to use nipy_spectral, Set1,Paired
colormap = plt.cm.terrain
# fig1 = plt.figure(1,figsize=(15,10))
# fig2 = plt.figure(2,figsize=(15,10))
fig1 = plt.figure(1)
fig2 = plt.figure(2)
ax = fig1.add_subplot(111)
ax2= fig2.add_subplot(111)
ax.set_color_cycle([colormap(i) for i in np.linspace(0, 0.25,len(directories))])
ax2.set_color_cycle([colormap(i) for i in np.linspace(0, 0.25,len(directories))])

for d in directories:
    try:
        f = d+'kappa.dat'
        print('Loading: '+f)
        data = np.loadtxt(f)
        extraD = np.loadtxt(d+'extra-data-point/kappa.dat')
        data   = np.concatenate((np.array([extraD]),data),axis=0)
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
    r = data[:,0]
    k = data[:,1]
    k_exact = data[:,2]
    plt.figure(1)
    # p = plt.plot(r,k,'x')
    p = plt.errorbar(r,k,yerr=df/2,fmt='.')
    print(df/2./k_exact)
    c = p[0].get_color()
    # plt.plot(r,k_exact,label="a = "+str(a),color=c)
    r0 = kappa_r0(a)
    rfine=np.linspace(r0,20.,500)
    plt.plot(rfine,kappa(rfine,a),label=r"$a = $"+str(a).rjust(4),color=c,lw=1)
    # plt.plot(rfine,kappa(rfine,a),label="a = "+str(a),color='k')
    plt.figure(2)
    plt.plot(r,np.abs(kappa(r,a)-k)/kappa(r,a),'o-')

plt.yscale('log')
plt.figure(1)

# plt.tick_params(axis='both', which='major', labelsize=lbfs)
plt.xlabel(r'$r$')
plt.ylabel(r'$\kappa$')
plt.xlim(xmax=20)
plt.ylim(ymin=0)
plt.legend(frameon=False,loc='upper right')
# plt.tight_layout()

plt.savefig('python_plot.pdf',bbox_inches='tight')

plt.show()
