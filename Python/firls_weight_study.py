# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 15:18:23 2023

Study the effect of weights on firls performance

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
    passband = 0.4
    stopband = 0.6
    bands = [0,passband,stopband,1]
    amps = [1,1,0,0]
    weights = np.array([1,1])
    
    passband_ripple_spec_db = 0.1 # +/- dB
    stopband_rej_spec_db = 50
    stopband_rej_spec_db = 80
    ntaps = 31 # scipy firls forces ntaps to be odd - is this generally required for a linear-phase filter?
    
    passband_lin_dev = 1-10**(-passband_ripple_spec_db/20)
    stopband_lin_dev = 10**(-stopband_rej_spec_db/20)
    
    plt.figure()
    
    # Flat weighting ---------------------------------------------------------------------------------------------------------------------------------
    
    print('_____Flat weighting_____')
    
    b = signal.firls(ntaps,bands,amps,weight=weights,fs=fs)
    w,h = signal.freqz(b,fs=fs)
    h_pb = 20*np.log10(abs(h[w <= passband]));
    h_sb = -20*np.log10(abs(h[w >= stopband]));
    wc_pb_ripple = max(abs(h_pb))
    wc_sb_rej = min(h_sb)
    
    print('Largest passband ripple (dB) = ' + str(round(wc_pb_ripple,3)))
    print('Smallest stopband rejection (dB) = ' + str(round(wc_sb_rej,1)))
    print('\n\n')
    
    plt.plot(w,20*np.log10(abs(h)),label='Flat weighting')
    
    # Weighting using linear deviations --------------------------------------------------------------------------------------------------------------
    
    print('_____Linear deviation weighting_____')
    
    weights = np.array([passband_lin_dev,stopband_lin_dev])
    weights = max(weights)/weights
    
    b = signal.firls(ntaps,bands,amps,weight=weights,fs=fs)
    w,h = signal.freqz(b,fs=fs)
    h_pb = 20*np.log10(abs(h[w <= passband]));
    h_sb = -20*np.log10(abs(h[w >= stopband]));
    wc_pb_ripple = max(abs(h_pb))
    wc_sb_rej = min(h_sb)
    
    print('Weights: ' + str(weights))
    print('Largest passband ripple (dB) = ' + str(round(wc_pb_ripple,3)))
    print('Smallest stopband rejection (dB) = ' + str(round(wc_sb_rej,1)))
    print('\n\n')
    
    plt.plot(w,20*np.log10(abs(h)),label='Linear deviation weighting')
    
    # Weighting using squared deviations ------------------------------------------------------------------------------------------------------------
    
    print('_____Squared deviation weighting_____')
    
    weights = weights**2
    
    b = signal.firls(ntaps,bands,amps,weight=weights,fs=fs)
    w,h = signal.freqz(b,fs=fs)
    h_pb = 20*np.log10(abs(h[w <= passband]));
    h_sb = -20*np.log10(abs(h[w >= stopband]));
    wc_pb_ripple = max(abs(h_pb))
    wc_sb_rej = min(h_sb)
    
    plt.plot(w,20*np.log10(abs(h)))
    
    print('Weights: ' + str(weights))
    print('Largest passband ripple (dB) = ' + str(round(wc_pb_ripple,3)))
    print('Smallest stopband rejection (dB) = ' + str(round(wc_sb_rej,1)))
    print('\n\n')
    
    plt.plot(w,20*np.log10(abs(h)),label='Squared deviation weighting')
    
    plt.title("FIRLS Weighting Study",{'fontsize':40})
    plt.xlabel("Normalized Digital Frequency",{'fontsize':30})
    plt.ylabel("Magnitude Response (dB)",{'fontsize':30})
    plt.legend(loc="upper right",fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    #plt.autoscale(enable=True,axis='both',tight=True)
    plt.grid()