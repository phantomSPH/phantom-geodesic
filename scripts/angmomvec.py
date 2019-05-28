import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import sys

print()
parser = argparse.ArgumentParser(description='files to analyse')
parser.add_argument('position', nargs=1, help='xyz positions as function of time')
parser.add_argument('velocities', nargs=1, help='vxyz velocities as function of time')
args = parser.parse_args()
file_xyz  = args.position[0]
file_vxyz = args.velocities[0]

# # Pandas is faster to load file
skiprows = 4
txyz  = np.array(pd.read_csv(file_xyz, delim_whitespace=True,skip_blank_lines=True,skiprows=skiprows,header=None,dtype=np.float64))
tvxyz = np.array(pd.read_csv(file_vxyz,delim_whitespace=True,skip_blank_lines=True,skiprows=skiprows,header=None,dtype=np.float64))

# Total angular momentum vector as function of time
Ltot  = np.cross(txyz[:,1:],tvxyz[:,1:])

# Magnitude of angmomvec at each time
# Lmag  = np.sqrt(np.sum(np.square(Ltot), axis=1))
Lmag = np.linalg.norm(Ltot,axis=1)

# Unit angular momentum vector at each time
Lhat  = Ltot/Lmag[:,None]       # makes some imaginary axes when dividing

# Projection of unit angmomvec in the x-y plane
Lxy = np.copy(Lhat)
Lxy[:,2] = 0.

# Convert to a unit vector
Lxymag = np.linalg.norm(Lxy, axis=1)
Lhatxy = Lxy/Lxymag[:,None]

# Find the angle between the angular momentum vector and the x-y plane
inc = np.arccos(np.sum(Lhat*Lhatxy,axis=1))*180./np.pi

# Find the angle between the angular momentum
rot = np.arccos(np.sum(Lhat*np.array([1,0,0]),axis=1))*180./np.pi

time = txyz[:,0]

r     = np.linalg.norm(txyz[:,1:],axis=1)
theta = np.arccos(txyz[:,3]/r)*180/np.pi
phi   = np.arctan2(txyz[:,2],txyz[:,1])*180/np.pi

# plt.plot(time,inc-inc[0],label='change in inclination')
# plt.plot(time,rot-rot[0],label='change in twist')
plt.figure()
plt.subplot(311)
plt.plot(time,inc,label='inclination')
plt.legend(loc='best')
plt.subplot(312)
plt.plot(time,rot,label='twist')
plt.legend(loc='best')
plt.subplot(313)
plt.plot(time,r,label='radius')
plt.legend(loc='best')
plt.show()

# plt.plot(time,theta)
# plt.plot(phi,txyz[:,3])
plt.plot(theta,phi)
plt.plot(theta[0],phi[0],color='blue',marker='.',ms=10)
plt.show()
