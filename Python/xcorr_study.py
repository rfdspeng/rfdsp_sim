# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 20:59:58 2022

Mess around with linear algebra concepts

@author: Ryan Tsai
"""

import numpy as np
import matplotlib.pyplot as plt
from IPython import get_ipython
import sys
sys.path.append("tools")
sys.path.append("models")

from ofdm_wavgen import ofdm_wavgen
from numpy import linalg as nplinalg
from scipy import linalg as sclinalg

if __name__ == '__main__':
    plt.close('all'); plt.ion()
    get_ipython().magic('reset -sf')
    
    # User options
    en_plot = 1
    en_tprecode = 1
    
    # Generate waveform
    nsym = 14; bw = 20; scs = 15; num_sc = 1200; start_sc = 600-round(num_sc/2)
    modorder = 4; ncp = 144/2048*100
    wola = 2
    x,x_standard,cfg_evm = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola)
    
    # Cross-correlation
    xcorr = np.correlate(x,x,mode='full')
    xcorr_abs = np.absolute(xcorr)
    xcorr_max = np.amax(xcorr_abs)
    lags = np.array(range(-len(x)+1,len(x)))
    plt.plot(lags,xcorr_abs/xcorr_max)
    
    # Generate random waveform