# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 07:02:46 2023

@author: Ryan Tsai
"""

import math
import numpy as np
from scipy import signal
from scipy import fft
import matplotlib.pyplot as plt
from IPython import get_ipython
import sys
sys.path.append("functions")
import ofdm
import calc
import digital_filter_design
import digital_hardware

if __name__ == '__main__':
    plt.close('all')
    get_ipython().magic('reset -sf')
    
    N = 10 # period
    
    x_base = np.zeros(N)
    x_base[0] = 1
    
    # Need to have multiple periods to emulate DTFT
    # Otherwise, taking DFT directly just yields a vector of ones
    x = x_base
    for idx in range(9):
        x = np.concatenate((x,x_base))
    
    x = x/10
    x_fft = fft.fft(x)
    
    w = np.arange(x_fft.size)*2*np.pi/x_fft.size
    
    plt.stem(w,abs(x_fft))    
    plt.title("100-point FFT of DT Impulse Train of Period 10",{'fontsize':40})
    plt.xlabel("Frequency (0 to 2pi)",{'fontsize':30})
    plt.ylabel("abs(FFT)",{'fontsize':30})
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    #plt.autoscale(enable=True,axis='both',tight=True)
    plt.grid()