# -*- coding: utf-8 -*-
"""
Created on Fri May 20 10:06:24 2022

Test calculate_papr.py

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
from calculate_papr import calculate_papr

if __name__ == '__main__':
    plt.close('all')
    get_ipython().magic('reset -sf')
    
    # Waveform params
    en_tprecodes = [1,0] # DFT-s-OFDM, CP-OFDM
    modorders = [4,16,64,256]
    wola = 0
    
    # Sweep waveforms
    out = {}
    out['Waveform'] = []
    out['99.99% PAPR (dB)'] = []
    for precode_idx in range(len(en_tprecodes)):
        en_tprecode = en_tprecodes[precode_idx]
        
        for mod_idx in range(len(modorders)):
            modorder = modorders[mod_idx]
    
            # Generate waveform
            nsym = 14*20; bw = 20; scs = 15; num_sc = 1200; start_sc = 600-round(num_sc/2); ncp = 7;
            x,x_standard,cfg_evm = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola)
            
            papr = calculate_papr(x,99.99)
            papr = round(papr,2)
            
            wavstr = 'DFT-s-OFDM ' if en_tprecode == 1 else 'CP-OFDM '
            wavstr = wavstr + 'QPSK' if modorder == 4 else wavstr + str(modorder) + 'QAM'
            out['Waveform'].append(wavstr)
            out['99.99% PAPR (dB)'].append(papr)
    
    out = pd.DataFrame(data=out)
    print(out)