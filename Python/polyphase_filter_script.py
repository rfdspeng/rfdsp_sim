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
    
    R = 3
    pdelta = 3 # Power delta in dB for downsampling example (blocker minus desired)
    
    # Generate waveforms for downsampling example ************************************************************************************************************************
    
    sig_frac_bits = 16
    
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
    
    if sig_frac_bits > 0:
        x = (x*2**sig_frac_bits).round()
    
    [p,f] = calc.psd(x,fs_in,fs_in/2048)
    plt.figure()
    plt.plot(f,10*np.log10(p))
    
    # Generate AAF *********************************************************************************************************************************************************
    
    ntaps = 31
    obw = num_sc_sig*scs/1000
    b = digital_filter_design.firls_rate_change('down',ntaps,obw,fs_in,R,en_plot=True)
    frac_bits = 0
    frac_bits = 16
    b = (b*2**frac_bits).round()
    b = b.astype(int)
    
    # Downsample *********************************************************************************************************************************************************
    
    y_prototype = signal.lfilter(b,1,x)
    if frac_bits > 0:
        y_prototype = (y_prototype/2**frac_bits).round()
    y_prototype_down = y_prototype[0:-1:R]
    
    y_polyphase = digital_hardware.polyphase_downsampler(x,b,frac_bits,R)
    
    plt.figure()
    
    [p,f] = calc.psd(y_prototype_down,fs_in/R,fs_in/R/2048)
    plt.plot(f,10*np.log10(p))
    
    [p,f] = calc.psd(y_polyphase,fs_in/R,fs_in/R/2048)
    plt.plot(f,10*np.log10(p))
    
    sum_delta = sum(y_polyphase-y_prototype_down)
    print('Sum delta = ' + str(sum_delta))