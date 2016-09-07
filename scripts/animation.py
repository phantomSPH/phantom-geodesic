"""
ORIGINALLY:
    Matplotlib Animation Example

    author: Jake Vanderplas
    email: vanderplas@astro.washington.edu
    website: http://jakevdp.github.com
    license: BSD
    Please feel free to use and modify this, but keep the above information. Thanks!

MODIFIED by: David Liptai, Monash University, 2016.

USE:


"""
import matplotlib
# matplotlib.use('TkAgg')
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import sys
from mpl_toolkits.mplot3d import Axes3D

try:
    datafilename1 = sys.argv[1]
except:
    quit('Quitting... No input file give.')

print('Loading files:',sys.argv[1:])
xyz = np.loadtxt(datafilename1,skiprows=2)
n   = int(np.genfromtxt(datafilename1,max_rows=1))
tnmax  = len(xyz)
x=xyz[:,1:n+1]
y=xyz[:,n+1:2*n+1]
z=xyz[:,2*n+1:]

# Choose how many lines of the data file to skip between each frame.
# i.e. Set the speed of the animation
# Use this if you have a lot of small steps.
speed = 2
n_frames = int(n/speed)

# First set up the figure, the axis, and the plot element we want to animate
markersize = 5
# fig = plt.figure(figsize=plt.figaspect(1)*1.5)
fig = plt.figure()
# ax = plt.axes()
ax = fig.gca(projection='3d')

# ax.axis('equal')
# ax.set_xlim([-100,100])
# ax.set_ylim([-100,100])
# ax.set_zlim([-10,10])

phi   = np.linspace(0, 2 * np.pi, 100)
theta = np.linspace(0, np.pi, 100)

x_bh = 2 * np.outer(np.cos(phi), np.sin(theta))
y_bh = 2 * np.outer(np.sin(phi), np.sin(theta))
z_bh = 2 * np.outer(np.ones(np.size(phi)), np.cos(theta))

body=[]
for i in range(n):
    body += ax.plot([],[],[],color='red',marker='.',ms=markersize)

# initialization function: plot the background of each frame
def init():
    body[0].set_data([],[])
    body[0].set_3d_properties([])
    ax.plot_surface(x_bh, y_bh, z_bh,color='black')
    return body[0],

# animation function.  This is called sequentially
def animate(i):
    for j in range(n):
        body[j].set_data(x[i*speed,j],y[i*speed,j])
        body[j].set_3d_properties(z[i*speed,j])

    return body[0],

# call the animator.  blit=True means only re-draw the parts that have changed.
# Note: when using the Mac OS X Backend, blit=True will not work!!
#       Need to manually set matplotlib.use('TkAgg') first....
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=n_frames, interval=1, blit=False)

plt.show()
