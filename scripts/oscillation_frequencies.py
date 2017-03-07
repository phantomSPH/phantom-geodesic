"""
Module containing functions that compute the exact frequencies of epicycles and vertical-oscillations,
given a radius and black hole spin value.
"""
import numpy as np
from scipy.optimize import brentq as root_find
import sys

# Epicyclic frequency
def kappa(r,a=0):
    k = None
    omega  = 1./(r**1.5 + a)
    kappa2 = omega**2*(1. - (6./r - 8.*a*r**(-1.5) + 3.*a**2/r**2))
    k = np.sqrt(kappa2)
    return k

# Given a value of spin, find the radius at which kappa=0. (i.e. the ISCO)
def kappa_r0(a):
    def term(r):
        return (1. - (6./r - 8.*a*r**(-1.5) + 3.*a**2/r**2))
    r0 = root_find(term,1.e-10,100.)
    # If the computed root is slightly to the left of the true root,
    # (i.e. gives kappa nan or neg or anything other than positive)
    # then increase the value of the root by small increments, until
    # kappa(r0) returns a real value that is >=0.
    count = 0
    while not kappa(r0,a)>=0.:
        count+=1
        r0 += 1.e-10
        if count>=100:
            print('WARNING: A real value of kappa was not returned, given the root found.')
            break
    return r0

# Wrapper to call kappa_r0, which actually also returns the ISCO
def isco(a):
    return kappa_r0(a)

# Vertical-oscillation frequency
def omega_z(r,a):
    omega   = 1./(r**1.5 + a)
    omegaz2 = omega**2*(1. - (4.*a*r**(-1.5) - 3.*a**2/r**2))
    omegaz = np.sqrt(omegaz2)
    return omegaz
