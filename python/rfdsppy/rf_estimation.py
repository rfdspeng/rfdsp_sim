# -*- coding: utf-8 -*-
"""
Created on Tue Nov 11 2025

Functions and classes for modeling RF estimation algorithms (FW or HW)

@author: Ryan Tsai
"""

import numpy as np

def est_tone_amp_phase(x: np.ndarray, fs, f0):
    """
    Estimate amplitude and phase of a sinusoid

    x: sinusoid (real or complex)
    fs = sampling rate (MHz)
    f0 = tone frequency (MHz)
    
    """

    w0 = f0*2*np.pi/fs

    c = (x * np.exp(-1j*w0*np.arange(x.size, dtype="float"))).sum()/x.size

    A = np.abs(c)
    phi = np.angle(c)

    return (A, phi)