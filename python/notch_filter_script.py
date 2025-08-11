# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 20:55:42 2023

@author: Ryan Tsai
"""

import math
import numpy as np
from scipy import signal
from scipy import fft
from scipy import linalg
import matplotlib.pyplot as plt
from IPython import get_ipython
import sys
sys.path.append("tools")
sys.path.append("models")
import ofdm
import calc

def notch_filter(x,wo):
    return 1

if __name__ == '__main__':
    plt.close('all')
    get_ipython().magic('reset -sf')
    
    nsym = 14
    bw = 20
    scs = 15
    num_sc = 100*12
    start_sc = 0
    modorder = 4
    en_tprecode = 0
    x,x_standard,cfg_ofdm = ofdm.ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode)
    
    wo = np.pi/2
    r = 0.9
    
    b = [1,-np.exp(-1j*wo)]
    a = [1,-r*np.exp(-1j*wo)]
    
    y = signal.lfilter(b,a,x)
    
    [px,fx] = calc.psd(x,cfg_ofdm['fs'],cfg_ofdm['fs']/2048)
    [py,fy] = calc.psd(y,cfg_ofdm['fs'],cfg_ofdm['fs']/2048)
    
    fig = plt.figure()
    plt.plot(fx,10*np.log10(px))
    plt.plot(fy,10*np.log10(py))