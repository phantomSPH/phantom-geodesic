"""
Module containing functions that compute the exact frequencies of epicycles and vertical-oscillations,
given a radius and black hole spin value.
"""
import numpy as np

def kappa(r,a=0):
    omega  = 1./(r**1.5 + a)
    kappa2 = omega**2*(1. - (6./r - 8.*a*r**(-1.5) + 3.*a**2/r**2))
    k = np.sqrt(kappa2)
    return k

def omega_z(r,a):
    omega   = 1./(r**1.5 + a)
    omegaz2 = omega**2*(1. - (4.*a*r**(-1.5) - 3.*a**2/r**2))
    omegaz = np.sqrt(omegaz2)
    return omegaz
