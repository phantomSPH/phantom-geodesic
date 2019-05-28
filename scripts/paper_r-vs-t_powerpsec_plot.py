import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import periodogram
from oneDpowerspec import powerspec,peak_freq#,plot_powerspec
from oscillation_frequencies import kappa
from scipy.interpolate import interp1d
import matplotlib as mpl
plt.style.use('paper')

mpl.rc('font', size=15)

colormap = plt.cm.terrain
cblue = colormap(0.)
cgreen= colormap(0.25)

data = np.loadtxt('positions.dat',skiprows=3)
a    = np.loadtxt('spin.value')
time = data[:,0]
s    = np.shape(data)
n    = int((s[1]-1)/3)
radius = data[:,1:s[1]:3]

r0 = np.zeros(n)
k  = np.zeros(n)

t_fine = np.linspace(time[0],time[-1],len(time)*10)

i = 0
rad_i = radius[:,i]

fig1 = plt.figure(1)
ax1  = fig1.add_subplot(111)
# ax.ticklabel_format(useOffset=False)
ax1.plot(time,rad_i,lw=1,color=cblue)
ax1.set_xlabel(r'$t$')
ax1.set_ylabel(r'$r$')
plt.xlim([0,max(time)])
plt.ylim([rad_i[0]-0.001,rad_i[0]+0.002])
plt.savefig('python_plot_a.pdf',bbox_inches='tight')

freqs,power = powerspec(time,rad_i)
r0[i]       = rad_i[0]
k[i]        = peak_freq(freqs,power)
k_approx = k[i]
k_exact  = kappa(rad_i[0],a)
# print(k_approx,k_exact,k_approx-k_exact)

fig2 = plt.figure(2)
ax2  = fig2.add_subplot(111)
# ax.ticklabel_format(useOffset=False)
ax2.plot(freqs,power,lw=1,color=cgreen)
ax2.set_xlabel(r'Frequency ( s$^{-1}$)')
ax2.set_ylabel(r'Power')

plt.xscale('log')
plt.xlim([min(freqs),max(freqs)])
plt.savefig('python_plot_b.pdf',bbox_inches='tight')
offset = ax2.get_yaxis().get_offset_text()
ax2.set_ylabel('{0} ({1})'.format(ax2.get_ylabel(), offset.get_text()) )
# ax.set_ylabel( str(ax.get_ylabel()) + '('+str(offset.get_text())+')' )
offset.set_visible(False)
plt.savefig('python_plot_b.pdf',bbox_inches='tight')
plt.show()

#
# fs=30
#
# plt.figure(figsize=(10,7.5))
# plt.plot(r0,k,'bx',label='Simulation')
# plt.plot(r0,kappa(r0,a),'g-',label='Exact')
# plt.legend(loc='best')
# # plt.plot(r_smooth,k_smooth,'k-')
# plt.ylabel(r'$\kappa$',fontsize=fs)
# plt.xlabel(r'$r$',fontsize=fs)
# plt.title('a = '+str(a))
# plt.savefig('python_plot.pdf')
# # plt.show()
#
# df = freqs[1]-freqs[0]
# print(df)
# np.savetxt('kappa.dat',np.column_stack((r0,k,kappa(r0,a))),header=str(a)+'\n'+str(df))
