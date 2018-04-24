"""
Plots the precessing orbit when positions.dat is in spherical coordinates
"""

import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl

mpl.rc('font', size=19)
mpl.rc('legend', fontsize=15)

t,r,theta,phi=np.loadtxt('positions.dat',skiprows=3,unpack=True)
circle1=plt.Circle((0,0),2,color='Black')

colormap = plt.cm.terrain
cblue = colormap(0.)
cgreen= colormap(0.25)

x2,y2,z2=np.loadtxt('clement_positions.dat',unpack=True)

x = r*np.sin(theta)*np.cos(phi)
y = r*np.sin(theta)*np.sin(phi)
z = r*np.cos(theta)

fig = plt.figure(1)
ax = fig.add_subplot(111,aspect='equal', adjustable='box-forced')
# ax.axis('equal')
ax.add_artist(circle1)
# plt.plot(x2,y2,'r',lw=2,label='Simulation')
# plt.plot(x,y,'k',label='Geodesic')
plt.plot(x2,y2,lw=2,color=cgreen,label='Geodesic')
plt.plot(x,y,lw=1,color=cblue,label='Sim.')
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.ylim([-110,110])
plt.xlim([-110,110])
plt.yticks([-100,0,100])
# plt.legend(frameon=False,bbox_to_anchor=(0,1.035),loc='upper left')
plt.legend(frameon=False,loc='upper left')
plt.tight_layout()
plt.savefig('python_plot.pdf',bbox_inches='tight')
plt.show()
