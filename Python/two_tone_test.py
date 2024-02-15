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
    
    
    wavtype = 'tones' # 'tones' or 'ofdm'
    #wavtype = 'ofdm'
    
    pin = 0 # dBm @ RF. For tones, this is power of one tone.
    piip2 = 20 # dBm @ RF
    piip3 = 20 # dBm @ RF
    
    
    fc = 61.44 # LO frequency
    
    
    # Generate waveforms
    if wavtype == 'tones': # Generate real baseband tone at fbb and upconvert to fc
        fs = 307.2 # sampling rate
        fbb = 3.84 # tones at +/- fbb
        
        # Generate RF input signal for RF simulation
        # r*cos((wc - wbb)t) + r*cos((wc + wbb)t)
        tx_bb = tonegen.tonegen(2**16,fs,fbb)
        r = math.sqrt(2*50*1e-3)*10**(pin/20) # RF tone amplitude
        tx_bb = tx_bb*r*2 # RF tone amplitude gets halved relative to BB amplitude due to upconversion
        tx_rf = trx.upconvert(tx_bb,fs,fc)
        
        # Generate baseband input signal for baseband simulation
        # r*exp(1j*wbb*t) + r*exp(-1j*wbb*t)
        tx_bb_m = tonegen.tonegen(2**16,fs,-fbb,cossin='exp')
        tx_bb_p = tonegen.tonegen(2**16,fs,+fbb,cossin='exp')
        tx_bb_m = tx_bb_m*r
        tx_bb_p = tx_bb_p*r
    elif wavtype == 'ofdm': # Generate OFDM baseband waveform and upconvert to fc
        nsym = 14; bw = 5; scs = 15; num_sc = 300; start_sc = 0; ncp = 7;
        modorder = 4; en_tprecode = 0; wola = 10;
        x,x_standard,cfg_evm = ofdm.ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,osr=40,ncp=ncp,wola=wola)
        fs = cfg_evm['fs']
        fbb = bw/2
        tx_bb = x/calc.rms(x)*calc.dbm2v(pin,'dBm')
        tx_rf = trx.upconvert(tx_bb,fs,fc)
    
    # Nonlinear intercepts
    aiip2 = math.sqrt(2*50*1e-3)*10**(piip2/20)
    aiip3 = math.sqrt(2*50*1e-3)*10**(piip3/20)
    
    a1 = 1
    
    # Apply RF nonlinear model
    rx_rf = a1*tx_rf + (a1/aiip2)*tx_rf**2 + (-4/3*a1/aiip3**2)*tx_rf**3
    
    # Apply baseband nonlinear model
    rx_bb = a1*(tx_bb_p + tx_bb_m) + (a1/aiip2)*(tx_bb_p**2 + tx_bb_m**2) + (-a1/aiip3**2)*(tx_bb_p**3 + tx_bb_m**3)
    
    """
    # Apply RF nonlinear model
    a1 = 1
    a2 = a1/aiip2
    a3 = -4/3*a1/aiip3**2
    rx_rf = a1*tx_rf + a2*tx_rf**2 + a3*tx_rf**3
    
    # Apply baseband nonlinear model
    if wavtype == 'tones':
        a1 = 1
        a2 = 2*a1/aiip2
        a3 = -4/3*a1/aiip3**2
        rx_bb = a1*(tx_bb_m + tx_bb_p) + a2/2*(tx_bb_m**2 + tx_bb_p**2) + 3/4*a3*(tx_bb_m**3 + tx_bb_p**3)
    elif wavtype == 'ofdm':
        a1 = 1
        a2 = a1/aiip2
        a3 = -4/3*a1/aiip3**2
        rx_bb = a1*tx_bb + a2/2*abs(tx_bb)**2 + 3/4*a3*abs(tx_bb)**2*tx_bb
    """
    
    # Isolate IM2 (RF model only)
    passband = fbb*5/(fs/2)
    stopband = (fc - fbb*5)/(fs/2)
    bands = [0,passband,stopband,1]
    amps = [1,1,0,0]
    weights = [1,1]
    ntaps = 31
    b = signal.firls(ntaps,bands,amps,weight=weights,fs=2)
    b = b/sum(b)
    rx_rf_im2 = signal.lfilter(b,1,rx_rf)
    
    # Plot filter response
    w,h = signal.freqz(b,fs=fs)
    plt.figure()
    plt.plot(w,20*np.log10(abs(h)))
    plt.title("Filter Response",{'fontsize':40})
    plt.xlabel("Frequency (MHz)",{'fontsize':30})
    plt.ylabel("Magnitude Response (dB)",{'fontsize':30})
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.autoscale(enable=True,axis='both',tight=True)
    plt.grid()
    
    # Downconvert (RF model only)
    rx_rf_down_raw = trx.downconvert(rx_rf,fs,fc)
    
    # Isolate desired signal
    rx_rf_down = signal.lfilter(b,1,rx_rf_down_raw)
    rx_rf_down = rx_rf_down + rx_rf_im2
    
    """
    Important note on using SPDFT to estimate RF signal power
    
    In the RF simulation, the RF input signal is r*cos((wc - wbb)t) + r*cos((wc + wbb)t)
    r is set s.t. the RF power in one tone is equal to Pin
    
    Let's say we use SPDFT to calculate power at +(wc - wbb). The real signal is r*cos((wc - wbb)t), so the signal component at +(wc - wbb) is r/2*exp(1j(wc - wbb)t).
    
    Therefore, SPDFT estimate yields r/2. In order to calculate r, we need to multiply SPDFT output by 2.
    
    Of course, the estimated relative IMD (dBc) won't change, but this is important if you want to confirm that the expected IMD matches the simulation (since this affects the signal power estimation).
    
    In the baseband simulation, the input signal is r*exp(1j*wbb*t) + r*exp(-1j*wbb*t), so the SPDFT directly estimates r. There is no need to multiply by 2. 
    
    """
    
    # Calculate IM2 using SPDFTs (RF model)
    if wavtype == 'tones':
        f_sig = (fc-fbb)
        f_imd = fbb*2
        tone_sig = tonegen.tonegen(2**16,fs,f_sig,cossin='exp') # Generate tone at f_sig
        tone_imd = tonegen.tonegen(2**16,fs,f_imd,cossin='exp') # Generate tone at f_imd
        v_sig = 2*np.mean(rx_rf*tone_sig) # Multiply by 2 to get r (see note above)
        v_imd = 2*np.mean(rx_rf*tone_imd)
        
        sig_dbm = 10*np.log10(abs(v_sig)**2/2/50/1e-3)
        imd_dbm = 10*np.log10(abs(v_imd)**2/2/50/1e-3)
        imd_dbc = imd_dbm - sig_dbm
        imd_dbc_calc = pin - piip2
 
    # Print sim outputs for IM2
    if wavtype == 'tones':
        print('RF two-tone simulation for IM2')
        print('---------------------------')
        print('IIP2 (dBm) = ' + str(np.round(piip2,2)))
        print('Sig (dBm) = ' + str(np.round(sig_dbm,2)))
        print('IM2 (dBm) = ' + str(np.round(imd_dbm,2)))
        print('IM2 (dBc) = ' + str(np.round(imd_dbc,2)))
        print('Expected IM2 (dBc) = ' + str(np.round(imd_dbc_calc,2)))
        print('---------------------------')
    
    # Calculate IM2 using SPDFTs (baseband model)
    if wavtype == 'tones':
        f_sig = fbb
        f_imd = fbb*2
        tone_sig = tonegen.tonegen(2**16,fs,f_sig,cossin='exp') # Generate tone at f_sig
        tone_imd = tonegen.tonegen(2**16,fs,f_imd,cossin='exp') # Generate tone at f_imd
        v_sig = np.mean(rx_bb*tone_sig) # No need to multiply by 2 (see note above)
        v_imd = np.mean(rx_bb*tone_imd)
        
        sig_dbm = 10*np.log10(abs(v_sig)**2/2/50/1e-3)
        imd_dbm = 10*np.log10(abs(v_imd)**2/2/50/1e-3)
        imd_dbc = imd_dbm - sig_dbm
        imd_dbc_calc = pin - piip2
    
    # Print sim outputs for IM2
    if wavtype == 'tones':
        print('Baseband two-tone simulation for IM2')
        print('---------------------------')
        print('IIP2 (dBm) = ' + str(np.round(piip2,2)))
        print('Sig (dBm) = ' + str(np.round(sig_dbm,2)))
        print('IM2 (dBm) = ' + str(np.round(imd_dbm,2)))
        print('IM2 (dBc) = ' + str(np.round(imd_dbc,2)))
        print('Expected IM2 (dBc) = ' + str(np.round(imd_dbc_calc,2)))
        print('---------------------------')
        
    # Calculate IM3 using SPDFTs (RF model)
    if wavtype == 'tones':
        f_sig = (fc-fbb)
        f_imd = fc-3*fbb
        tone_sig = tonegen.tonegen(2**16,fs,f_sig,cossin='exp') # Generate tone at f_sig
        tone_imd = tonegen.tonegen(2**16,fs,f_imd,cossin='exp') # Generate tone at f_imd
        v_sig = 2*np.mean(rx_rf*tone_sig) # Multiply by 2 to get r (see note above)
        v_imd = 2*np.mean(rx_rf*tone_imd)
        
        sig_dbm = 10*np.log10(abs(v_sig)**2/2/50/1e-3)
        imd_dbm = 10*np.log10(abs(v_imd)**2/2/50/1e-3)
        imd_dbc = imd_dbm - sig_dbm
        imd_dbc_calc = 2*pin - 2*piip3
    
    # Print sim outputs for IM3
    if wavtype == 'tones':
        print('RF two-tone simulation for IM3')
        print('---------------------------')
        print('IIP3 (dBm) = ' + str(np.round(piip3,2)))
        print('Sig (dBm) = ' + str(np.round(sig_dbm,2)))
        print('IM3 (dBm) = ' + str(np.round(imd_dbm,2)))
        print('IM3 (dBc) = ' + str(np.round(imd_dbc,2)))
        print('Expected IM3 (dBc) = ' + str(np.round(imd_dbc_calc,2)))
        print('---------------------------')
    
    # Calculate IM3 using SPDFTs (baseband model)
    if wavtype == 'tones':
        f_sig = fbb
        f_imd = fbb*3
        tone_sig = tonegen.tonegen(2**16,fs,f_sig,cossin='exp') # Generate tone at f_sig
        tone_imd = tonegen.tonegen(2**16,fs,f_imd,cossin='exp') # Generate tone at f_imd
        v_sig = np.mean(rx_bb*tone_sig) # No need to multiply by 2 (see note above)
        v_imd = np.mean(rx_bb*tone_imd)
        
        sig_dbm = 10*np.log10(abs(v_sig)**2/2/50/1e-3)
        imd_dbm = 10*np.log10(abs(v_imd)**2/2/50/1e-3)
        imd_dbc = imd_dbm - sig_dbm
        imd_dbc_calc = 2*pin - 2*piip3
    
    # Print sim outputs for IM3
    if wavtype == 'tones':
        print('Baseband two-tone simulation for IM3')
        print('---------------------------')
        print('IIP3 (dBm) = ' + str(np.round(piip3,2)))
        print('Sig (dBm) = ' + str(np.round(sig_dbm,2)))
        print('IM3 (dBm) = ' + str(np.round(imd_dbm,2)))
        print('IM3 (dBc) = ' + str(np.round(imd_dbc,2)))
        print('Expected IM3 (dBc) = ' + str(np.round(imd_dbc_calc,2)))
        print('---------------------------')

    # Plots    
    plt.figure()
    if wavtype == 'tones':
        [p,f] = calc.psd(tx_bb_m+tx_bb_p,fs,fs/2048)
    elif wavtype == 'ofdm':
        [p,f] = calc.psd(tx_bb,fs,fs/2048)
    plt.plot(f,10*np.log10(p),color='blue',linewidth=5,label='Tx BB')
    [p,f] = calc.psd(rx_bb,fs,fs/2048)
    plt.plot(f,10*np.log10(p),color='orange',label='Rx BB (Post-NL)')
    
    plt.title("Baseband Nonlinear Model",{'fontsize':40})
    plt.xlabel("Frequency (MHz)",{'fontsize':30})
    plt.ylabel("PSD (dB)",{'fontsize':30})
    plt.legend(loc="lower right",fontsize=10)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.autoscale(enable=True,axis='both',tight=True)
    plt.xlim((-5*fbb,+5*fbb))
    plt.ylim(bottom=-100)
    plt.grid()
    
    plt.figure()
    [p,f] = calc.psd(tx_bb,fs,fs/2048)
    plt.plot(f,10*np.log10(p),color='blue',linewidth=10,label='Tx BB')
    [p,f] = calc.psd(tx_rf,fs,fs/2048)
    plt.plot(f,10*np.log10(p),color='red',linewidth=10,label='Tx RF (Upconverted)')
    [p,f] = calc.psd(rx_rf,fs,fs/2048)
    plt.plot(f,10*np.log10(p),color='green',linewidth=5,label='Rx RF (Post-NL)')
    [p,f] = calc.psd(rx_rf_down,fs,fs/2048)
    plt.plot(f,10*np.log10(p),color='orange',label='Rx BB (Post-NL and Downconverted)')
    
    plt.title("RF Nonlinear Model",{'fontsize':40})
    plt.xlabel("Frequency (MHz)",{'fontsize':30})
    plt.ylabel("PSD (dB)",{'fontsize':30})
    plt.legend(loc="lower right",fontsize=10)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.autoscale(enable=True,axis='both',tight=True)
    plt.xlim((-fc-5*fbb,+fc+5*fbb))
    plt.xlim((-5*fbb,+5*fbb))
    plt.ylim(bottom=-100)
    plt.grid()