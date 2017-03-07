"""
Plots the precessing orbit when positions.dat is in spherical coordinates
"""

import numpy as np
import matplotlib.pyplot as plt
plt.style.use('paper')

colormap = plt.cm.terrain
cblue = colormap(0.)
cgreen= colormap(0.25)

tzero,xzero,yzero,zzero = np.loadtxt('positions-zero0.0.dat',skiprows=3,unpack=True)
tpos ,xpos ,ypos ,zpos  = np.loadtxt('positions-pos0.1.dat' ,skiprows=3,unpack=True)
tneg ,xneg ,yneg ,zneg  = np.loadtxt('positions-neg0.1.dat' ,skiprows=3,unpack=True)

circle1=plt.Circle((0,0),2,color='Black')
circle2=plt.Circle((0,0),2,color='Black')

fig1 = plt.figure(1)
fig2 = plt.figure(2)
ax1 = fig1.add_subplot(111,aspect='equal', adjustable='box-forced')
ax2 = fig2.add_subplot(111,aspect='equal', adjustable='box-forced')
# ax1.axis('equal')
# ax2.axis('equal')
ax1.add_artist(circle1)
ax2.add_artist(circle2)

plt.figure(1)
# plt.plot(xzero,yzero,'r',label=r'$a=0.0$')
# plt.plot(xpos,ypos,'k',label=r'$a=0.1$')
plt.plot(xzero,yzero,color=cblue,lw=2,label=r'$a = $'+str(0.0))
plt.plot(xpos,ypos,color=cgreen,lw=1,label=r'$a = $'+str(0.1))
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.ylim([-110,110])
plt.xlim([-110,110])
plt.legend(frameon=False,loc='upper right')
plt.savefig('python_plot_a.pdf',bbox_inches='tight')
# plt.tight_layout()

plt.figure(2)
# plt.plot(xzero,yzero,'r',label=r'$a=\quad \, 0.0$')
# plt.plot(xneg,yneg,'k',label=r'$a=-0.1$')
plt.plot(xzero,yzero,color=cblue,lw=2,label=r'$a = $'+str(0.0))  #r'$a=\quad \, 0.0$'
plt.plot(xneg,yneg,color=cgreen,lw=1,label=r'$a = $'+str(-0.1))
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.ylim([-110,110])
plt.xlim([-110,110])
plt.legend(frameon=False,loc='upper right')
plt.savefig('python_plot_b.pdf',bbox_inches='tight')

plt.show()
