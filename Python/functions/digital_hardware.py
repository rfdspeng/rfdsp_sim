# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 09:55:56 2023

Functions to mimic digital hardware

@author: Ryan Tsai
"""

import math
import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
from scipy import signal
from scipy import fft

def polyphase_downsampler(x,b,frac_bits,R):
    """
    Description
    -----------
    Implements cycle-based polyphase downsampler

    Parameters
    ----------
    x : input signal
    b : prototype filter coefficients
    frac_bits : fractional bitwidth for normalization (0 means floating point)
    R : integer downsampling ratio

    Returns
    -------
    y : output signal

    """
    
    # Generate branches
    branch_len = math.ceil(len(b)/R)
    b_pp = np.zeros((R,branch_len))
    gd = int((len(b)-1)/2)
    b = np.concatenate((b,np.zeros(b_pp.size-len(b)+1)))
    for branch in range(R):
        b_pp[branch,:] = b[branch:-1:R]
    
    reg_in = np.zeros(branch_len*R) + 1j*np.zeros(branch_len*R) # input tapped delay line
    
    y = np.zeros(math.floor(len(x)/R)) + 1j*np.zeros(math.floor(len(x)/R))
    
    # Zero pad the input signal
    x = np.concatenate((np.zeros(len(reg_in)-1),x,np.zeros(gd)))
    
    branch = 0 # branch index
    y_br = np.zeros(R) + 1j*np.zeros(R) # branch outputs
    ydx = 0 # output index
    for x_end in range(len(reg_in)-1,len(x)):
        # Update input tapped delay line
        reg_in[:] = np.flip(x[x_end-len(reg_in)+1:x_end+1])
        x_tapped = reg_in[0:-1:R]
        
        b_br = b_pp[branch,:] # branch coefficients
        
        # Compute branch output
        conv_out = sum(b_br*x_tapped)
        #if frac_bits > 0:
        #    conv_out = round(conv_out/2**frac_bits)
    
        y_br[branch] = conv_out
        
        # Compute output sample
        if branch == 0:
            y[ydx] = sum(y_br)
            if frac_bits > 0:
                y[ydx] = round(y[ydx]/2**frac_bits)
            ydx += 1
        
        # Update branch
        if branch == 0:
            branch = R-1
        else:
            branch -= 1
        """
        if branch == 2:
            branch = 0
        else:
            branch += 1
        """
        
        if ydx >= len(y):
            break

    return y