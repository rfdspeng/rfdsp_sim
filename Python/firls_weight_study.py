# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 15:18:23 2023

@author: tsair
"""

import math
import numpy as np
from scipy import signal
from scipy import fft
import matplotlib.pyplot as plt
from IPython import get_ipython
import sys
sys.path.append("tools")
sys.path.append("models")

if __name__ == '__main__':
    plt.close('all')
    get_ipython().magic('reset -sf')
    
    fs = 2 # Normalized digital frequency
    bands = [0,0.4,0.6,1]
    amps = [1,1,0,0]
    weights = np.array([1,1])
    
    passband_ripple_spec_db = 0.1 # +/- dB
    stopband_rej_spec_db = 50
    stopband_rej_spec_db = 100
    ntaps = 15 # scipy firls forces ntaps to be odd - is this generally required for a linear-phase filter?
    
    passband_lin_dev = 1-10**(-passband_ripple_spec_db/20)
    stopband_lin_dev = 10**(-stopband_rej_spec_db/20)
    
    b = signal.firls(ntaps,bands,amps,weight=weights,fs=fs)
    w,h = signal.freqz(b,fs=fs)

    plt.figure()
    plt.plot(w,20*np.log10(abs(h)))
    
    #plt.title("Kaiser Window FFT in Log Domain",{'fontsize':40})
    #plt.xlabel("Normalized Digital Frequency",{'fontsize':30})
    #plt.ylabel("FFT (dB)",{'fontsize':30})
    #plt.legend(loc="upper right",fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    #plt.xlim([0, 0.015])
    #plt.ylim([-300, 0])
    #plt.autoscale(enable=True,axis='both',tight=True)
    plt.grid()
    
    weights = np.array([passband_lin_dev,stopband_lin_dev])
    weights = max(weights)/weights
    
    b = signal.firls(ntaps,bands,amps,weight=weights,fs=fs)
    w,h = signal.freqz(b,fs=fs)
    
    plt.plot(w,20*np.log10(abs(h)))
    
    weights = weights**2
    
    b = signal.firls(ntaps,bands,amps,weight=weights,fs=fs)
    w,h = signal.freqz(b,fs=fs)
    
    plt.plot(w,20*np.log10(abs(h)))
    
    """ Print LSE to desired response """