# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 09:13:34 2023

Functions to design digital filters

@author: Ryan Tsai
"""

import math
import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
from scipy import signal
from scipy import fft
from typing import Literal

def firls_rate_change(updn: Literal["up", "down"], ntaps, obw, fs_in, R, passband_ripple_spec_db=0.1, stopband_rej_spec_db=50, en_plot=False):
    """
    Description
    -----------
    Generates AAF/interpolation filter coefficients

    Parameters
    ----------
    updn : 'up' or 'down' to specify upsampling or downsampling
    ntaps : filter length
    obw : occupied BW of desired lowpass signal
    fs_in : input sampling rate
    R : integer rate change
    passband_ripple_spec_db : passband ripple spec in dB
    stopband_rej_spec_db : stopband rejection spec in dB

    Returns
    -------
    b : AAF/interpolation filter coefficients

    """
    
    # Determine filter sampling rate
    if updn == 'up':
        fs = fs_in*R 
    elif updn == 'down':
        fs = fs_in
    
    # Determine passband and stopband
    passband = (obw/2)/(fs/2)
    stopband = (fs/R-obw/2)/(fs/2)
    bands = [0, passband, stopband, 1]
    amps = [1, 1, 0, 0]
    
    # Determine firls weighting
    passband_lin_dev = 1-10**(-passband_ripple_spec_db/20)
    stopband_lin_dev = 10**(-stopband_rej_spec_db/20)
    weights = np.array([passband_lin_dev, stopband_lin_dev])
    # weights = max(weights)/weights
    # weights = weights**2
    weights = (weights.max()/weights)**2
    
    # Generate filter coefficients
    b = signal.firls(ntaps, bands, amps, weight=weights, fs=2)
    w, h = signal.freqz(b, fs=2)
    # h_pb = 20*np.log10(abs(h[w <= passband]))
    # h_sb = -20*np.log10(abs(h[w >= stopband]))
    # wc_pb_ripple = max(abs(h_pb))
    # wc_sb_rej = min(h_sb)
    h_pb = 20*np.log10(np.abs(h[w <= passband]))
    h_sb = -20*np.log10(np.abs(h[w >= stopband]))
    wc_pb_ripple = np.abs(h_pb).max() # dB
    wc_sb_rej = h_sb.min() # dB
    
    print('digital_filter_design.firls_rate_change()')
    print('Largest passband ripple (dB) = ' + str(round(wc_pb_ripple,3)))
    print('Smallest stopband rejection (dB) = ' + str(round(wc_sb_rej,1)))
    print('\n\n')
    
    if en_plot:
        plt.figure()
        plt.plot(w, 20*np.log10(abs(h)))
        plt.title("digital_filter_design.firls_rate_change()", {'fontsize':40})
        plt.xlabel("Normalized Digital Frequency", {'fontsize':30})
        plt.ylabel("Magnitude Response (dB)", {'fontsize':30})
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.grid()

    return b