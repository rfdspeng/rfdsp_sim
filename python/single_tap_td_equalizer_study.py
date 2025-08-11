# -*- coding: utf-8 -*-
"""
Created on Tue May 24 18:30:38 2022

Project vector x onto vector y.

For real vectors, the dot product is commutative, so the equalizer coefficient is the same regardless of order: y'x/y'y or x'y/y'y
However, for complex vectors, the order matters because y'x is not equal to x'y. The coefficient is calculated as y'x/y'y.

@author: tsair
"""

import numpy as np
from IPython import get_ipython
import sys
sys.path.append("tools")
sys.path.append("models")

from ofdm_wavgen import ofdm_wavgen

if __name__ == '__main__':
    get_ipython().magic('reset -sf')
    
    # Generate waveform
    nsym = 14; bw = 20; scs = 15; num_sc = 1200; start_sc = 600-round(num_sc/2)
    modorder = 4; en_tprecode = 1; ncp = 7; wola = 2
    x,x_standard,cfg_evm = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola)
    
    # Random gain
    rng = np.random.default_rng()
    g = rng.random() + 1j*rng.random()
    y = g*x
    
    # TD equalizer
    c1 = np.vdot(y,x)/np.vdot(y,y)
    c2 = np.vdot(x,y)/np.vdot(y,y)
    y1 = c1*y
    y2 = c2*y
    
    # SNR
    e = y-x
    snr = 10*np.log10(np.vdot(x,x)/np.vdot(e,e)).real
    print('Baseline SNR (dB): ' + str(snr))
    e = y1-x
    snr = 10*np.log10(np.vdot(x,x)/np.vdot(e,e)).real
    print('SNR with correct equalizer (dB): ' + str(snr))
    e = y2-x
    snr = 10*np.log10(np.vdot(x,x)/np.vdot(e,e)).real
    print('SNR with incorrect equalizer (dB): ' + str(snr))