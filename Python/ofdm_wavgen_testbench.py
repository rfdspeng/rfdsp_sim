# -*- coding: utf-8 -*-
"""
Created on Fri May 20 10:06:24 2022

Test ofdm_wavgen.py and ofdm_evm_calculator.py

@author: Ryan Tsai
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython import get_ipython
import sys
sys.path.append("tools")
sys.path.append("models")

from ofdm_wavgen import ofdm_wavgen
from ofdm_evm_calculator import ofdm_evm_calculator
from calculate_psd import calculate_psd

if __name__ == '__main__':
    plt.close('all')
    get_ipython().magic('reset -sf')
    
    en_td_plot = 0
    en_psd_plot = 0
    en_const_plot = 1
    
    # Waveform params
    en_tprecodes = [1,0] # DFT-s-OFDM, CP-OFDM
    modorders = [4,16,64,256]
    wola = 10
    
    # Sweep waveforms
    out = {}
    out['Waveform'] = []
    out['EVM (%)'] = []
    out['SNR (dB)'] = []
    for precode_idx in range(len(en_tprecodes)):
        en_tprecode = en_tprecodes[precode_idx]
        
        for mod_idx in range(len(modorders)):
            modorder = modorders[mod_idx]
    
            # Generate waveform
            nsym = 14; bw = 20; scs = 15; num_sc = 1200; start_sc = 600-round(num_sc/2); ncp = 7;
            x,x_standard,cfg_evm = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola)
            
            wola_len = cfg_evm['wola_len']
            fs = cfg_evm['fs']
            
            # Plot time-domain samples
            if en_td_plot:
                x_standard_plt = np.concatenate((np.zeros(round(wola_len/2)), x_standard, np.zeros(round(wola_len/2))))
                fig = plt.figure()
                plt.plot(abs(x_standard_plt))
                plt.plot(abs(x))
            
            # Plot PSDs
            if en_psd_plot:
                rbw = scs/1000/2**2
                p,f = calculate_psd(x,fs,rbw)
                ps,f = calculate_psd(x_standard,fs,rbw)
                fig = plt.figure()
                plt.plot(f,10*np.log10(ps))
                plt.plot(f,10*np.log10(p))
                plt.title("PSD",{'fontsize':40})
                plt.xlabel("Frequency (MHz)",{'fontsize':30})
                plt.ylabel("PSD (dB/Bin)",{'fontsize':30})
                plt.xticks(fontsize=20)
                plt.yticks(fontsize=20)
                plt.autoscale(enable=True,axis='both',tight=True)
                plt.grid()
            
            # EVM
            if en_const_plot:
                cfg_evm['en_plot'] = 1
            evm = ofdm_evm_calculator(cfg_evm,x_standard,x[round(wola_len/2):])
            snr = round(-20*np.log10(evm/100),2)
            
            wavstr = 'DFT-s-OFDM ' if en_tprecode == 1 else 'CP-OFDM '
            wavstr = wavstr + 'QPSK' if modorder == 4 else wavstr + str(modorder) + 'QAM'
            out['Waveform'].append(wavstr)
            out['EVM (%)'].append(evm)
            out['SNR (dB)'].append(snr)
    
    out = pd.DataFrame(data=out)
    print(out)