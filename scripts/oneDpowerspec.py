"""
Module containing functions to compute/plot power spectra, as well as find the peak value.

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import periodogram
from oscillation_frequencies import kappa

# Wrapper to call the desired way of computing the power spectrum.
def powerspec(time,radius):
    freqs,power=powerspec_periodogram(time,radius)
    # freqs,power=powerspec_fft(time,radius)
    return freqs,power

def powerspec_periodogram(time,radius):
    dt = time[1]-time[0]
    sampling_frequency=1/dt
    freqs,power=periodogram(radius,sampling_frequency,scaling='spectrum')
    #Subtract the mean values
    # dr = radius-np.mean(radius)
    # freqs,power=periodogram(dr,sampling_frequency,scaling='spectrum')
    # Turn frequencies into angular frequencies
    freqs = 2*np.pi*freqs
    return freqs,power

def powerspec_fft(time,radius):
    dt = time[1]-time[0]

    dr = (radius-np.mean(radius))
    dr = 1./np.max(dr)*dr
    power = np.abs(np.fft.fft(dr))**2
    freqs = abs(np.fft.fftfreq(dr.size,dt))
    # power = np.abs(np.fft.fft(radius))**2
    # freqs = np.fft.fftfreq(radius.size,dt)

    idx   = np.argsort(freqs)
    freqs = freqs[idx]
    power = power[idx]
    freqs = 2*np.pi*freqs
    return freqs,power


def plot_powerspec(freqs,power):
    # plt.figure()
    plt.plot(freqs,power)
    plt.grid()
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlabel('Angular frequency')
    plt.ylabel('Power')
    # plt.show()

def peak_freq(freqs,power):
   maxpower = max(power)
   maxp_freq= freqs[power==maxpower]
   return maxp_freq[0]
