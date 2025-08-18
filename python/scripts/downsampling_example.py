# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 11:04:59 2023

@author: Ryan Tsai
"""

import math
import numpy as np
from scipy import signal
from scipy import fft
from scipy import linalg
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
    
    R = 3
    pdelta = 3 # Power delta in dB for downsampling example (blocker minus desired)
    
    # Generate waveforms for downsampling example ************************************************************************************************************************
    
    nsym = 14
    bw = 20
    scs = 15
    modorder = 4
    en_tprecode = 0
    
    # Desired signal
    num_sc_sig = 100*12
    start_sc_sig = 0
    x_signal,x_standard,cfg_ofdm = ofdm.ofdm_wavgen(nsym,bw,scs,num_sc_sig,start_sc_sig,modorder,en_tprecode,osr=R)
    fs_in = cfg_ofdm['fs']
    
    # Blocker
    num_sc_block = 20
    start_sc_block = round((100*12-num_sc_block)/2)
    x_blocker,_,_ = ofdm.ofdm_wavgen(nsym,bw,scs,num_sc_block,start_sc_block,modorder,en_tprecode,osr=R)
    w_blocker = 2*np.pi/R
    x_blocker = x_blocker*np.exp(1j*w_blocker*np.array(range(len(x_blocker)))) # Frequency shift blocker to correct location
    
    # Add signal and blocker
    x_signal = x_signal/calc.rms(x_signal)
    x_blocker = x_blocker/calc.rms(x_blocker)*calc.rms(x_signal)*10**(pdelta/20)
    x = x_signal+x_blocker
    
    [p,f] = calc.psd(x,fs_in,fs_in/2048)
    plt.figure()
    plt.plot(f,10*np.log10(p))
    
    # Generate AAF *********************************************************************************************************************************************************
    
    """
    passband = num_sc_sig*scs/1000/fs_in
    stopband = (fs_in/R-num_sc_block*scs/1000/2)/(fs_in/2)
    bands = [0,passband,stopband,1]
    amps = [1,1,0,0]
    passband_ripple_spec_db = 0.1 # +/- dB
    stopband_rej_spec_db = 50
    ntaps = 15
    
    passband_lin_dev = 1-10**(-passband_ripple_spec_db/20)
    stopband_lin_dev = 10**(-stopband_rej_spec_db/20)
    
    weights = np.array([passband_lin_dev,stopband_lin_dev])
    weights = max(weights)/weights
    weights = weights**2
    
    b = signal.firls(ntaps,bands,amps,weight=weights,fs=2)
    w,h = signal.freqz(b,fs=2)
    h_pb = 20*np.log10(abs(h[w <= passband]));
    h_sb = -20*np.log10(abs(h[w >= stopband]));
    wc_pb_ripple = max(abs(h_pb))
    wc_sb_rej = min(h_sb)
    
    plt.figure()
    plt.plot(w,20*np.log10(abs(h)))
    
    print('AAF weights: ' + str(weights))
    print('Largest passband ripple (dB) = ' + str(round(wc_pb_ripple,3)))
    print('Smallest stopband rejection (dB) = ' + str(round(wc_sb_rej,1)))
    print('\n\n')
    """
    
    ntaps = 31
    obw = num_sc_sig*scs/1000
    b = digital_filter_design.firls_rate_change('down',ntaps,obw,fs_in,R,en_plot=True)
    
    # Downsample *********************************************************************************************************************************************************
    
    y_aaf = signal.lfilter(b,1,x)
    #y_aaf = x
    y_down = y_aaf[list(range(0,len(y_aaf),R))]
    [p,f] = calc.psd(y_down,fs_in/R,fs_in/R/2048)
    plt.figure()
    plt.plot(f,10*np.log10(p))