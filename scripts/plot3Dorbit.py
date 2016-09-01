"""
Custom code to create a 3D plot of an orbit.

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

xyz = np.loadtxt(datafilename)
x   = xyz[:,0]
y   = xyz[:,1]
z   = xyz[:,2]
mpl.rcParams['legend.fontsize'] = 10
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x,y,z, label='orbit',color='b')
ax.legend()
plt.show()

make_figure = False
while make_figure is False:
    usrinpt = input('Save figure? (y/n): ')
    yes=['y','Y','yes','Yes','YES','true','True','TRUE']
    no =['n','N','no','No','NO','false','False','FALSE']
    if usrinpt in yes:
        make_figure=True
    if usrinpt in no:
        print('Bye.')
        quit()

if make_figure:
    fig.savefig('figure.pdf')
    