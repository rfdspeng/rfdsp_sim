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
    
    
    
    
    fc = 61.44
        
    if wavtype == 'tones':
        # Generate real baseband tone at fbb and upconvert to fc
        fs = 307.2
        fbb = 3.84
        tx_bb = tonegen.tonegen(2**16,fs,fbb)
        tx_bb = tx_bb*2 # Play with signal level
        tx_rf = trx.upconvert(tx_bb,fs,fc)
        p_avg,p_pk = calc.power_dbm(tx_bb)
    elif wavtype == 'ofdm':  
        # Generate OFDM baseband waveform and upconvert to fc
        nsym = 14; bw = 5; scs = 15; num_sc = 300; start_sc = 0; ncp = 7;
        modorder = 4; en_tprecode = 0; wola = 10;
        x,x_standard,cfg_evm = ofdm.ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,osr=40,ncp=ncp,wola=wola)
        fs = cfg_evm['fs']
        fbb = bw/2
        p_avg = 13
        tx_bb = x/calc.rms(x)*calc.dbm2v(p_avg,'dBm')
        p_avg,p_pk = calc.power_dbm(tx_bb)
        tx_rf = trx.upconvert(tx_bb,fs,fc)
    
    # Nonlinear model specs
    piip2 = 40 # dBm @ RF
    piip3 = 30 # dBm @ RF
    
    #aiip2 = math.sqrt(2)*10**(piip2/20)
    #aiip3 = math.sqrt(2)*10**(piip3/20)
    aiip2 = math.sqrt(2*50*1e-3)*10**(piip2/20)
    aiip3 = math.sqrt(2*50*1e-3)*10**(piip3/20)

    a1 = 1
    a2 = a1/aiip2
    a3 = -4/3*a1/aiip3**2
    #a3 = 0
    
    # Apply nonlinear models (baseband and RF)
    rx_bb = a1*tx_bb + 1/2*a2*abs(tx_bb)**2 + 3/4*a3*abs(tx_bb)**2*tx_bb
    rx_rf = a1*tx_rf + a2*tx_rf**2 + a3*tx_rf**3
    
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
    
    # Calculate IM2 using SPDFTs
    if wavtype == 'tones':
        #"""
        f_sig = (fc-fbb)
        f_imd = fbb*2
        tone_sig = tonegen.tonegen(2**16,fs,f_sig,cossin='exp') # Generate tone at f_sig
        tone_imd = tonegen.tonegen(2**16,fs,f_imd,cossin='exp') # Generate tone at f_imd
        v_sig = np.mean(rx_rf*tone_sig)
        v_imd = np.mean(rx_rf*tone_imd)
        
        sig_dbm = calc.dbm2v(abs(v_sig),'V')
        imd_dbm = calc.dbm2v(abs(v_imd),'V')
        imd_dbc = imd_dbm - sig_dbm
        imd_dbc_calc = sig_dbm + 6 - piip2
        #"""
        """
        f_sig = fbb
        f_imd = fbb*2
        tone_sig = tonegen.tonegen(2**16,fs,f_sig,cossin='exp') # Generate tone at f_sig
        tone_imd = tonegen.tonegen(2**16,fs,f_imd,cossin='exp') # Generate tone at f_imd
        v_sig = np.mean(rx_bb*tone_sig)
        v_imd = np.mean(rx_bb*tone_imd)
        
        sig_dbm = calc.dbm2v(abs(v_sig),'V')
        imd_dbm = calc.dbm2v(abs(v_imd),'V')
        imd_dbc = imd_dbm - sig_dbm
        imd_dbc_calc = sig_dbm - piip2
        """
        
        
        
        
        
    # Print sim outputs for IM2
    if wavtype == 'tones':
        print('Two-tone simulation for IM2')
        print('---------------------------')
        print('IIP2 (dBm) = ' + str(np.round(piip2,2)))
        print('Sig (dBm) = ' + str(np.round(sig_dbm,2)))
        print('IM2 (dBm) = ' + str(np.round(imd_dbm,2)))
        print('IM2 (dBc) = ' + str(np.round(imd_dbc,2)))
        print('Expected IM2 (dBc) = ' + str(np.round(imd_dbc_calc,2)))
        print('---------------------------')
        
    # Calculate IM3 using SPDFTs
    if wavtype == 'tones':
        f_sig = (fc-fbb)
        f_imd = fc-3*fbb
        tone_sig = tonegen.tonegen(2**16,fs,f_sig,cossin='exp') # Generate tone at f_sig
        tone_imd = tonegen.tonegen(2**16,fs,f_imd,cossin='exp') # Generate tone at f_imd
        v_sig = np.mean(rx_rf*tone_sig)
        v_imd = np.mean(rx_rf*tone_imd)
        
        sig_dbm = calc.dbm2v(abs(v_sig),'V')
        imd_dbm = calc.dbm2v(abs(v_imd),'V')
        imd_dbc = imd_dbm - sig_dbm
        imd_dbc_calc = 2*sig_dbm - 2*piip3
    
    # Print sim outputs for IM2
    if wavtype == 'tones':
        print('Two-tone simulation for IM3')
        print('---------------------------')
        print('IIP3 (dBm) = ' + str(np.round(piip3,2)))
        print('Sig (dBm) = ' + str(np.round(sig_dbm,2)))
        print('IM3 (dBm) = ' + str(np.round(imd_dbm,2)))
        print('IM3 (dBc) = ' + str(np.round(imd_dbc,2)))
        print('Expected IM3 (dBc) = ' + str(np.round(imd_dbc_calc,2)))
        print('---------------------------')
        


    # Plots    
    plt.figure()
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
    plt.ylim(bottom=-50)
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
    plt.ylim(bottom=-50)
    plt.grid()