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
    #common_phase = 0
    lo_phase_mismatch = 2
    lo_gain_mismatch = 0.1
    rx_bb_iqmm_raw = trx.downconvert(tx_rf,fs,fc,common_phase=common_phase,lo_phase_mismatch=lo_phase_mismatch,lo_gain_mismatch=lo_gain_mismatch)
    
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
    spdft_freq = round(fbb/(fs/len(rx_bb_iqmm)))*(fs/len(rx_bb_iqmm))
    tone_pfbb = tonegen.tonegen(2**16,fs,spdft_freq,cossin='exp') # Generate tone at +fbb
    tone_nfbb = tonegen.tonegen(2**16,fs,-spdft_freq,cossin='exp') # Generate tone at -fbb
    plt.figure()
    [p,f] = calc.psd(tone_pfbb,fs,fs/2048)
    plt.plot(f,10*np.log10(p))
    [p,f] = calc.psd(tone_nfbb,fs,fs/2048)
    plt.plot(f,10*np.log10(p))
    rsb = np.mean(tone_pfbb*rx_bb_iqmm) # Correlate Rx signal with +fbb tone to get RSB information
    sig = np.mean(tone_nfbb*rx_bb_iqmm) # Correlate Rx signal with -fbb tone to get signal information
    
    theoretical_rsb_dbc = -10*np.log10((lo_gain_mismatch**2+(lo_phase_mismatch*np.pi/180)**2)/4)
    print('Theoretical RSB (dBc) = ' + str(theoretical_rsb_dbc))
    
    est_rsb_dbc = 20*np.log10(abs(sig)/abs(rsb))
    print('Estimated RSB (dBc) = ' + str(est_rsb_dbc))
    
    r = np.abs(sig)
    theta0 = np.angle(sig)
    
    rsb_processed = rsb*np.exp(1j*theta0)/r
    ep = np.real(rsb_processed)*2
    theta = np.imag(rsb_processed)*2
    
    print('r = ' + str(r))
    print('theta (deg) = ' + str(theta*180/np.pi))
    print('epsilon = ' + str(ep))
    
    # Compensate IQ mismatch
    i_mm = np.real(rx_bb_iqmm)
    q_mm = np.imag(rx_bb_iqmm)
    i_mm_comp = i_mm/(1+ep/2)/np.cos(theta/2) - q_mm*np.tan(theta/2)
    q_mm_comp = q_mm/(1-ep/2)/np.cos(theta/2) - i_mm*np.tan(theta/2)
    rx_bb_iqmm_comp = i_mm_comp + 1j*q_mm_comp
    
    rsb_comp = np.mean(tone_pfbb*rx_bb_iqmm_comp) # Correlate Rx signal with +fbb tone to get RSB information
    sig_comp = np.mean(tone_nfbb*rx_bb_iqmm_comp) # Correlate Rx signal with -fbb tone to get signal information
    est_rsb_comp_dbc = 20*np.log10(abs(sig_comp)/abs(rsb_comp))
    print('Estimated RSB after compensation (dBc) = ' + str(est_rsb_comp_dbc))
    
    plt.figure()
    [p,f] = calc.psd(tx_bb,fs,fs/2048)
    plt.plot(f,10*np.log10(p),color='blue',linewidth=10,label='Tx BB')
    [p,f] = calc.psd(tx_rf,fs,fs/2048)
    plt.plot(f,10*np.log10(p),color='purple',label='Tx RF')
    [p,f] = calc.psd(rx_bb_iqmm,fs,fs/2048)
    plt.plot(f,10*np.log10(p),color='orange',linewidth=5,label='Rx BB w/LO IQ Mismatch')
    [p,f] = calc.psd(rx_bb_iqmm_comp,fs,fs/2048)
    plt.plot(f,10*np.log10(p),label='Rx BB w/LO IQ Mismatch + Compensation')
    
    plt.title("Rx LO IQ Mismatch",{'fontsize':40})
    plt.xlabel("Frequency (MHz)",{'fontsize':30})
    plt.ylabel("PSD (dB)",{'fontsize':30})
    plt.legend(loc="upper right",fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.autoscale(enable=True,axis='both',tight=True)
    plt.grid()