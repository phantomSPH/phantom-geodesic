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
    3D animation of particle orbits.
    Call this script with 'positions.dat' file as an input
    [Should contain: (time,xyz,xyz...) as columns]

"""
import matplotlib
# matplotlib.use('TkAgg')
import numpy as np
from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation
import sys

try:
    datafilename1 = sys.argv[1]
except:
    quit('Quitting... No input file give.')

print('Loading files:',sys.argv[1:])
xyz = np.loadtxt(datafilename1,skiprows=2)
n   = int(np.genfromtxt(datafilename1,max_rows=1))
tnmax  = len(xyz)
# x=xyz[:,1:n+1]
# y=xyz[:,n+1:2*n+1]
# z=xyz[:,2*n+1:]
x=xyz[:,1:tnmax:3]
y=xyz[:,2:tnmax:3]
z=xyz[:,3:tnmax:3]

# Choose how many lines of the data file to skip between each frame.
# i.e. Set the speed of the animation
# Use this if you have a lot of small steps.
speed = 2
n_frames = int(tnmax/speed)

# First set up the figure, the axis, and the plot element we want to animate
markersize = 5
# fig = plt.figure(figsize=plt.figaspect(1)*1.5)
# ax = plt.axes()
fig = plt.figure(figsize=(15*2,11*2))
ax  = p3.Axes3D(fig)
# ax = fig.gca(projection='3d')

# ax.axis('equal')
# ax.set_xlim([-100,100])
# ax.set_ylim([-100,100])
# ax.set_zlim([-100,100])
ax.set_xlim3d([-100.0, 100.0])
ax.set_xlabel('X')
ax.set_ylim3d([-100.0, 100.0])
ax.set_ylabel('Y')
ax.set_zlim3d([-100.0, 100.0])
ax.set_zlabel('Z')

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
    ax.plot_surface(x_bh, y_bh, z_bh,color='black')
    return body[0],

# animation function.  This is called sequentially
def animate(i):
    for j in range(n):
        body[j].set_data(x[i*speed,j],y[i*speed,j])
        body[j].set_3d_properties(z[i*speed,j])
    return body,

#### call the animator.  blit=True means only re-draw the parts that have changed.
#### Note: when using the Mac OS X Backend, blit=True will not work!!
####    Need to manually set matplotlib.use('TkAgg') first....
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=n_frames, interval=1, blit=False)
plt.show()

# #### This is for creating image files to make a movie
# init()
# for i in range(tnmax):
#     animate(i)
#     png_name = 'pythonplot_'+str(i).zfill(5)+'.png'
#     print('Writing '+png_name)
#     # plt.show()
#     plt.savefig(png_name)
