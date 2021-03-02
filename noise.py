import numpy as np
import scipy.signal as ss

def make_noise(num_samples=10000, mean_fr = 1):
    fv = np.linspace(0, 1, 20);                                # Normalised Frequencies
    a = 1/(1 + fv*2);                                      # Amplitudes Of '1/f'
    b = ss.firls(43, fv, a);                                   # Filter Numerator Coefficients
    wn = np.random.rand(1,num_samples) + mean_fr - 0.5
    invfn = ss.filtfilt(b, 1, wn);                             # Create '1/f' Noise
    return invfn[0]