# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 17:34:40 2022

Test bench for generating a "realistic" PA model

Assume APT mode needs to support up to 18dBm without saturating for a max PA bias of 3.5V

@author: tsair
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
from calculate_compression import calculate_compression
import pa_model

if __name__ == '__main__':
    plt.close('all')
    get_ipython().magic('reset -sf')
    
    # Generate waveform
    nsym = 14; bw = 20; scs = 15; num_sc = 1200; start_sc = 600-round(num_sc/2)
    modorder = 4; en_tprecode = 1; ncp = 7; wola = 2
    x,x_standard,cfg_evm = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola)
    
    x = x/max(abs(x))
    y = pa_model.pa_model(x)
    
    cfg_comp = {'en_plot':1}
    comp,nlse = calculate_compression(x,y,cfg=cfg_comp)