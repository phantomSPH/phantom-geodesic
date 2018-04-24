"""
Plots a circular orbit at radius of 2 when positions.dat is in spherical coordinates
for the kerr metric with spin +1.0, and therefore event horizon = 1 + sqrt(1 - a^2)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rc('font', size=18)

colormap = plt.cm.terrain
cblue = colormap(0.)
cgreen= colormap(0.25)

spin = 1.0
rh = 1. + np.sqrt(1.-spin**2)

t,r,theta,phi=np.loadtxt('positions.dat',skiprows=3,unpack=True)
circle1=plt.Circle((0,0),rh,color='Black')

x = r*np.sin(theta)*np.cos(phi)
y = r*np.sin(theta)*np.sin(phi)
z = r*np.cos(theta)

phi_exact = np.linspace(0,2.*np.pi,100)
x_exact = 2.*np.cos(phi_exact)
y_exact = 2.*np.sin(phi_exact)

fig = plt.figure(1)
ax = fig.add_subplot(111,aspect='equal', adjustable='box-forced')
ax.add_artist(circle1)
# ax.axis('equal')
# plt.plot(x_exact,y_exact,'r',lw=2,label='Exact')
# plt.plot(x,y,'k',label='Simulation')
plt.plot(x_exact,y_exact,lw=2,color=cgreen,label='Exact')
plt.plot(x,y,lw=0.5,color=cblue,label='Simulation')
plt.xlim([-11,11])
plt.ylim([-11,11])
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.tight_layout()
plt.legend(frameon=False,loc='upper right',fontsize='x-small')
plt.savefig('python_plot.pdf',bbox_inches='tight')
plt.show()
