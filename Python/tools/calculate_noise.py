# -*- coding: utf-8 -*-
"""
Created on Sat May 21 08:56:25 2022

@author: Ryan Tsai
"""

import math
import numpy as np
from scipy import signal
from scipy import fft
import matplotlib.pyplot as plt
from calculate_psd import calculate_psd

def calculate_noise(x,fs,rbw,sigf,noisef,cfg={}):
    wintype = cfg['wintype'] if 'wintype' in cfg else 'kaiser'
    [p,f] = calculate_psd(x,fs,rbw,wintype)
    
    sigfl = sigf[0]; sigfh = sigf[1]
    psig = p[(f >= sigfl) & (f <= sigfh)]
    fsig = f[(f >= sigfl) & (f <= sigfh)]
    
    noisefl = noisef[0]; noisefh = noisef[1]
    pnoise = p[(f >= noisefl) & (f <= noisefh)]
    fnoise = f[(f >= noisefl) & (f <= noisefh)]
    
    psig_sum = sum(psig)
    pnoise_sum = sum(pnoise)
    noise_dbc = 10*np.log10(pnoise_sum/psig_sum)
    
    en_plot = cfg['en_plot'] if 'en_plot' in cfg else 0
    if en_plot:
        fig = plt.figure()
        titlestr = cfg['title'] if 'title' in cfg else 'PSD'
        plt.plot(f,10*np.log10(p),linewidth=2.5,label='Full PSD')
        plt.plot(fsig,10*np.log10(psig),label='Signal')
        plt.plot(fnoise,10*np.log10(pnoise),label='Noise')
        plt.title(titlestr,{'fontsize':40})
        plt.xlabel("Frequency (MHz)",{'fontsize':30})
        plt.ylabel("PSD (dBm)",{'fontsize':30})
        plt.legend(loc="lower center",fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.grid()
    
    return noise_dbc