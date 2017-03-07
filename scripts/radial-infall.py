"""
Plots the veloctity as a function of time for a radially infalling particle.
Also computes and compares this to the exact value.
"""

import numpy as np
import matplotlib.pyplot as plt
plt.style.use('paper')

colormap = plt.cm.terrain
cblue = colormap(0.)
cgreen= colormap(0.25)

def vr(r,r0):
   vel=(1-2/r)/np.sqrt(1-2/r0)*np.sqrt(2*(1/r-1/r0))
   return vel

# fs = 40
# plt.figure(figsize=(15,10))
#
# for r0 in range(5,40,5):
#     r=np.linspace(2,r0,1000)
#     plt.plot(r,vr(r,r0))
#
# plt.ylim(ymin=0)
# plt.ylabel(r'$v_r$',fontsize=fs)
# plt.xlabel(r'$r$',fontsize=fs)
# plt.show()

positions=np.loadtxt('positions.dat',skiprows=3)
velocities=np.loadtxt('velocities.dat',skiprows=3)

s = np.shape(positions)
n = int((s[1]-1)/3)
radius = positions[:,1:s[1]:3]
velr   = velocities[:,1:s[1]:3]

# fs = 40
fig1 = plt.figure()
colormap = plt.cm.terrain
ax = fig1.add_subplot(111)
plt.style.use('paper')
# plt.subplot(211)
plt.ylabel(r'$v_r$')
plt.xlabel(r'$r$')
# ax.set_color_cycle([colormap(i) for i in np.linspace(0, 0.25,n)])
# plt.ylim(ymin=0,ymax=np.nanmax(np.abs(velr)))
# plt.subplot(212)
# plt.ylabel(r'$\Delta v_r / v_r$',fontsize=fs)
# plt.ylabel(r'$\Delta v_r$',fontsize=fs)
# plt.xlabel(r'$r$',fontsize=fs)
# plt.yscal e('log')

for i in range(n):
    r = radius[:,i]
    v = abs(velr[:,i])
    vexact = vr(r,r[0])
    # print(r[0])

    if i==0:
        # plt.subplot(211)
        # plt.plot(r,vexact,'r-',lw=2,label='Exact')
        # plt.plot(r,v,'k-',label='Simulation')
        plt.plot(r,vexact,color=cgreen,lw=2,label='Exact')
        plt.plot(r,v,color=cblue,lw=1,label='Simulation')

        # plt.subplot(212)
        # plt.plot(r,abs(v-vexact),'r-',label='Residual')
    else:
        # plt.subplot(211)
        # plt.plot(r,vexact,'r-',lw=2)
        # plt.plot(r,v,'k-')
        plt.plot(r,vexact,color=cgreen,lw=2)
        plt.plot(r,v,color=cblue,lw=1)


        # plt.subplot(212)
        # plt.plot(r,abs(v-vexact),'r-')

plt.xlim(xmin=0,xmax=42)
plt.ylim(ymin=0,ymax=0.43)
plt.legend(frameon=False,loc='upper right')
# plt.subplot(211)
# plt.legend()
# plt.tight_layout()
plt.savefig('python_plot.pdf',bbox_inches='tight')
plt.show()
