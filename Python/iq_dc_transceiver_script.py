# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 08:46:21 2023

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
import tonegen
import iq_dc_transceiver as trx

if __name__ == '__main__':
    plt.close('all')
    get_ipython().magic('reset -sf')
    
    fs = 307.2
    fc = 61.44
    fbb = 3.84
    
    # Generate Tx complex baseband tone at fbb and upconvert to fc
    tx_bb = tonegen.tonegen(2**16,fs,fbb,cossin='exp')
    tx_rf = trx.upconvert(tx_bb,fs,fc)
    
    # Downconvert
    rx_bb_raw = trx.downconvert(tx_rf,fs,fc)
    
    # Filter 2fc component
    passband = fbb*2/(fs/2)
    stopband = (fc*2-fbb*2)/(fs/2)
    bands = [0,passband,stopband,1]
    amps = [1,1,0,0]
    weights = [1,1]
    ntaps = 31
    b = signal.firls(ntaps,bands,amps,weight=weights,fs=2)
    rx_bb = signal.lfilter(b,1,rx_bb_raw)
    
    # Plot filter response
    w,h = signal.freqz(b,fs=fs)
    plt.figure()
    plt.plot(w,20*np.log10(abs(h)))
    
    # Plot
    """
    plt.figure()
    [p,f] = calc.psd(tx_bb,fs,fs/2048)
    plt.plot(f,10*np.log10(p),label='Tx BB')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.autoscale(enable=True,axis='both',tight=True)
    plt.grid()
    
    plt.figure()
    [p,f] = calc.psd(tx_rf,fs,fs/2048)
    plt.plot(f,10*np.log10(p),label='Tx RF')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.autoscale(enable=True,axis='both',tight=True)
    plt.grid()
    
    plt.figure()
    [p,f] = calc.psd(rx_bb,fs,fs/2048)
    plt.plot(f,10*np.log10(p),label='Rx BB')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.autoscale(enable=True,axis='both',tight=True)
    plt.grid()
    """
    
    plt.figure()
    [p,f] = calc.psd(tx_bb,fs,fs/2048)
    plt.plot(f,10*np.log10(p),color='blue',linewidth=5,label='Tx BB')
    [p,f] = calc.psd(tx_rf,fs,fs/2048)
    plt.plot(f,10*np.log10(p),color='purple',label='Tx RF')
    [p,f] = calc.psd(rx_bb,fs,fs/2048)
    plt.plot(f,10*np.log10(p),color='orange',label='Rx BB')
    
    plt.title("Simple IQ Up/Downconversion",{'fontsize':40})
    plt.xlabel("Frequency (MHz)",{'fontsize':30})
    plt.ylabel("PSD (dB)",{'fontsize':30})
    plt.legend(loc="upper right",fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.autoscale(enable=True,axis='both',tight=True)
    plt.grid()
    
    # Downconvert with LO IQ mismatch
    common_phase = np.random.rand()*360
    rx_bb_iqmm_raw = trx.downconvert(tx_rf,fs,fc,common_phase=common_phase,lo_phase_mismatch=3,lo_gain_mismatch=0.05)
    
    # Filter 2fc component
    rx_bb_iqmm = signal.lfilter(b,1,rx_bb_iqmm_raw)
    
    plt.figure()
    [p,f] = calc.psd(tx_bb,fs,fs/2048)
    plt.plot(f,10*np.log10(p),color='blue',linewidth=5,label='Tx BB')
    [p,f] = calc.psd(tx_rf,fs,fs/2048)
    plt.plot(f,10*np.log10(p),color='purple',label='Tx RF')
    [p,f] = calc.psd(rx_bb_iqmm,fs,fs/2048)
    plt.plot(f,10*np.log10(p),color='orange',label='Rx BB w/LO IQ Mismatch')
    
    plt.title("Rx LO IQ Mismatch",{'fontsize':40})
    plt.xlabel("Frequency (MHz)",{'fontsize':30})
    plt.ylabel("PSD (dB)",{'fontsize':30})
    plt.legend(loc="upper right",fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.autoscale(enable=True,axis='both',tight=True)
    plt.grid()
    
    # Estimate IQ mismatch
    tone_pfbb = tonegen.tonegen(2**16,fs,fbb,cossin='exp') # Generate tone at +fbb
    tone_nfbb = tonegen.tonegen(2**16,fs,-fbb,cossin='exp') # Generate tone at -fbb
    rsb = np.mean(tone_pfbb*rx_bb_iqmm) # Correlate Rx signal with +fbb tone to get RSB information
    sig = np.mean(tone_nfbb*rx_bb_iqmm) # Correlate Rx signal with -fbb tone to get signal information
    rsb_angle = np.angle(rsb)
    rsb_abs = abs(rsb)
    
    sig1 = sig*np.exp(1j*rsb_angle)
    sig2 = sig1/np.real(sig1)
    theta = np.imag(sig2)*-2
    sig1_abs = abs(sig1)
    r = sig1_abs/np.sqrt(1+theta**2/4)
    ep = rsb_abs*2/r
    