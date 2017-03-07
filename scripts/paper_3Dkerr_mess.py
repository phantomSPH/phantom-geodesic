import numpy as np
import matplotlib.pyplot as plt
plt.style.use('paper')

colormap = plt.cm.terrain
cblue = colormap(0.)
cgreen= colormap(0.25)

tzero,x,y,z = np.loadtxt('positions.dat',skiprows=3,unpack=True)

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
plt.plot(x,y,color=cblue)
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.ylim([-110,110])
plt.xlim([-110,110])
# plt.savefig('python_plot_a.pdf',bbox_inches='tight')
# plt.tight_layout()

plt.figure(2)
# plt.plot(xzero,yzero,'r',label=r'$a=0.0$')
# plt.plot(xpos,ypos,'k',label=r'$a=0.1$')
plt.plot(x,z,color=cblue)
plt.xlabel(r'$x$')
plt.ylabel(r'$z$')
plt.ylim([-110,110])
plt.xlim([-110,110])
# plt.savefig('python_plot_b.pdf',bbox_inches='tight')
# plt.tight_layout()

plt.show()
