"""
Code to create a plot of an orbit in 3D.

INPUT: file of x,y,z data in 3 seperate columns

USE: python plot3Dorbit.py filename

Written by:
David Liptai, Monash University, 2016.

(Python3)
"""

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import sys

try:
    datafilename = sys.argv[1]
except:
    quit('Quitting... No input file give.')
print('Loading file:',datafilename)


phi   = np.linspace(0, 2 * np.pi, 100)
theta = np.linspace(0, np.pi, 100)

x_bh = 2 * np.outer(np.cos(phi), np.sin(theta))
y_bh = 2 * np.outer(np.sin(phi), np.sin(theta))
z_bh = 2 * np.outer(np.ones(np.size(phi)), np.cos(theta))


xyz = np.loadtxt(datafilename,skiprows=4)
x   = xyz[:,1]
y   = xyz[:,2]
z   = xyz[:,3]
mpl.rcParams['legend.fontsize'] = 10
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x,y,z, label='orbit',color='b')
ax.plot_surface(x_bh, y_bh, z_bh,color='black')
ax.legend()
plt.show()

# make_figure = False
# while make_figure is False:
#     usrinpt = input('Save figure? (y/n): ')
#     yes=['y','Y','yes','Yes','YES','true','True','TRUE']
#     no =['n','N','no','No','NO','false','False','FALSE']
#     if usrinpt in yes:
#         make_figure=True
#     if usrinpt in no:
#         print('Bye.')
#         quit()
#
# if make_figure:
#     fig.savefig('figure.pdf')
