# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 08:55:31 2023

Functions to generate tones

@author: Ryan Tsai
"""

import math
import numpy as np
from scipy import signal
from scipy import fft
import matplotlib.pyplot as plt
import calc

def tonegen(nsamp,fs,fc,cossin='cos',theta0=0):
    """
    Generates real or complex tone at fc with starting phase of theta0
    x = cos(wn + theta0) or sin(wn + theta0) or exp(j(wn + theta0))

    Parameters
    ----------
    nsamp : number of samples
    fs : sampling rate
    fc : tone frequency
    cossin : 'cos', 'sin', or 'exp'
    theta0 : starting phase in degrees

    Returns
    -------
    x : real or complex tone

    """
    
    # Convert theta0 from degrees to radians
    theta0 = theta0*np.pi/180
    
    wc = np.pi*fc/(fs/2) # digital LO frequency in rad/s
    n = np.array(range(nsamp))
    if cossin == 'cos':
        x = np.cos(wc*n + theta0)
    elif cossin == 'sin':
        x = np.sin(wc*n + theta0)
    elif cossin == 'exp':
        x = np.exp(1j*(wc*n + theta0))

    return x