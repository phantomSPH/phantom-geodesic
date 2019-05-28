"""
Plot the angular momentum conservation
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.style.use('paper')

colormap = plt.cm.terrain
cblue = colormap(0.)
cgreen= colormap(0.25)
mpl.rc('font', size=14)

L_init = 4.7500155230447643
T_orbit = 2390.

t1,E1,L1 = np.loadtxt('ev.dat-precDP-spo10k-tol1e-7',unpack=True)
t2,E2,L2 = np.loadtxt('ev.dat-precDP-spo10k-tol1e-10',unpack=True)
t3,E3,L3 = np.loadtxt('ev.dat-precDP-spo10k-tol1e-15',unpack=True)

t4,E4,L4 = np.loadtxt('ev.dat-precQP-spo10k-tol1e-15',unpack=True)
t5,E5,L5 = np.loadtxt('ev.dat-precQP-spo10k-tol1e-22',unpack=True)
t6,E6,L6 = np.loadtxt('ev.dat-precQP-spo10k-tol1e-30',unpack=True)

L=np.array([L1,L2,L3,L4,L5,L6])
T=np.array([t1,t2,t3,t4,t5,t6])
T=T/T_orbit
L=np.abs(L/L_init)
L[L==0.] = np.nan # Convert to nan when error is zero, so that there's no infinite lines on a log scale

# plt.plot(T[0,:],L[0,:],T[1,:],L[1,:],T[2,:],L[2,:],color=cblue)
# plt.plot(T[3,:],L[3,:],T[4,:],L[4,:],T[5,:],L[5,:],color=cgreen)

plt.plot(T[0,:],L[0,:],color=cblue,ls='solid',label=r'tol=1e-7')
plt.plot(T[1,:],L[1,:],color=cblue,ls='dashed',label=r'tol=1e-10')
plt.plot(T[2,:],L[2,:],color=cblue,ls='dotted',label=r'tol=1e-15')

plt.plot(T[3,:],L[3,:],color=cgreen,ls='solid',label=r'tol=1e-15')
plt.plot(T[4,:],L[4,:],color=cgreen,ls='dashed',label=r'tol=1e-22')
plt.plot(T[5,:],L[5,:],color=cgreen,ls='dotted',label=r'tol=1e-30')

# plt.plot(0.,1.e-30,color='black',ls='solid',label='1k steps/orbit')
# plt.plot(0.,1.e-30,color='black',ls='dashed',label='10k steps/orbit')
# plt.plot(0.,1.e-30,color='black',ls='dotted',label='100k steps/orbit')

# plt.grid()
plt.xlim(0.,np.max(T))
plt.ylim(1e-34,1e1-1.)
plt.yscale('log')
plt.xlabel('Number of orbits')
plt.ylabel(r'$\Delta L/L_0$')
plt.legend(loc='best',framealpha=0.,ncol=2)
# plt.tight_layout()
plt.savefig('python_plot.pdf',bbox_inches='tight')
plt.show()
