"""
Plots a circular orbit when positions.dat is in spherical coordinates
for the Schwarzschild metric where event horizon = 2M
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rc('font', size=18)

colormap = plt.cm.terrain
cblue = colormap(0.)
cgreen= colormap(0.25)


t,r,theta,phi = np.loadtxt('positions.dat',skiprows=3,unpack=True)
circle1 = plt.Circle((0,0),2,color='Black')

x = r*np.sin(theta)*np.cos(phi)
y = r*np.sin(theta)*np.sin(phi)
z = r*np.cos(theta)

phi_exact = np.linspace(0,2.*np.pi,100)
x_exact = 10.*np.cos(phi_exact)
y_exact = 10.*np.sin(phi_exact)

fig = plt.figure(1)
ax = fig.add_subplot(111,aspect='equal', adjustable='box-forced')
# ax.axis('equal')
ax.add_artist(circle1)
# plt.plot(x_exact,y_exact,'r',lw=2,label='Exact')
# plt.plot(x,y,'k',label='Simulation')
plt.plot(x_exact,y_exact,lw=2,color=cgreen,label='Exact')
plt.plot(x,y,lw=0.5,color=cblue,label='Simulation')
plt.xlim([-11,11])
plt.ylim([-11,11])
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.tight_layout()
# plt.legend(frameon=False,loc='upper right',fontsize='smaller')
plt.savefig('python_plot.pdf',bbox_inches='tight')
plt.show()
