# -*- coding: utf-8 -*-
"""
Created on Fri May 20 10:06:24 2022

Test WOLA implementation in OFDM waveform generator

@author: Ryan Tsai
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from IPython import get_ipython
import sys
sys.path.append("tools")
sys.path.append("models")

from ofdm_wavgen import ofdm_wavgen
from ofdm_evm_calculator import ofdm_evm_calculator
from calculate_psd import calculate_psd
from scale_psd import scale_psd

if __name__ == '__main__':
    plt.close('all'); plt.ion()
    get_ipython().magic('reset -sf')
    
    # User options
    en_plot = 1
    wolas = np.array(range(0,24,4)) # WOLA percentage sweep
    en_tprecode = 1
    
    # WOLA study
    wola_lens = np.zeros(len(wolas)) # store WOLA lengths
    wavs = {} # store waveforms for plotting
    evms = np.zeros(len(wolas)) # store EVM
    for wdx in range(0,len(wolas)):
        # Generate waveform
        nsym = 14; bw = 20; scs = 15; num_sc = 1200; start_sc = 600-round(num_sc/2)
        modorder = 4; ncp = 144/2048*100
        wola = wolas[wdx]
        x,x_standard,cfg_evm = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola)
        wavs[wdx] = x
        
        wola_len = round(cfg_evm['wola_len'])
        fs = cfg_evm['fs']
        wola_lens[wdx] = wola_len
        
        # EVM
        cfg_evm['en_plot'] = en_plot; cfg_evm['title'] = "WOLA % = " + str(wola)
        evm = ofdm_evm_calculator(cfg_evm,x_standard,x[round(wola_len/2):])
        evms[wdx] = evm
    
    max_wola_len = max(wola_lens)
    # Plot samples and PSDs
    if en_plot:
        fig_samp = plt.figure()
        fig_psd = plt.figure()
        for wdx in range(0,len(wola_lens)):
            wola = wolas[wdx]
            wola_len = round(wola_lens[wdx])
            delta = max_wola_len - wola_len
            x = wavs[wdx]
            
            # Samples
            x = np.concatenate((np.zeros(round(delta/2)), x, np.zeros(round(delta/2))))
            plt.figure(fig_samp)
            plt.plot(abs(x),label=str(wola))
            
            # PSD
            rbw = scs/1000/2**2
            p,f = calculate_psd(x,fs,rbw)
            p = scale_psd(p,f,bw,scs,start_sc,num_sc)
            plt.figure(fig_psd)
            plt.plot(f,10*np.log10(p),label=str(wola))
    
    plt.figure(fig_samp)
    plt.title("Envelope",{'fontsize':40})
    plt.xlabel("Sample Index",{'fontsize':30})
    plt.ylabel("Envelope",{'fontsize':30})
    plt.legend(loc="upper right",fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.autoscale(enable=True,axis='both',tight=True)
    plt.grid()
    
    plt.figure(fig_psd)
    plt.title("PSD",{'fontsize':40})
    plt.xlabel("Frequency (MHz)",{'fontsize':30})
    plt.ylabel("PSD (dBm/Bin)",{'fontsize':30})
    plt.legend(loc="lower center",fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.autoscale(enable=True,axis='both',tight=True)
    plt.grid()
    
    outarray = np.concatenate((wolas.reshape(len(wola_lens),1), evms.reshape(len(evms),1).round(3)), axis=1)
    print("WOLA (%) / EVM (%)")
    print(outarray)