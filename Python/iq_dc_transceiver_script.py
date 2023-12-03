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
    tx_bb = tonegen.tonegen(2**12,fs,fbb,cossin='exp')
    tx_rf = trx.upconvert(tx_bb,fs,fc)
    
    plt.figure()
    plt.plot(np.real(tx_bb))
    plt.plot(np.imag(tx_bb))
    
    tx_bb_fft = fft.fft(tx_bb)
    plt.figure()
    plt.plot(20*np.log10(abs(tx_bb_fft)))
    
    # Downconvert
    rx_bb = trx.downconvert(tx_rf,fs,fc)
    
    
    
    
    # Plot
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
    [p,f] = calc.psd(tx_rf,fs,fs/2048)
    plt.plot(f,10*np.log10(p),label='Tx RF')
    [p,f] = calc.psd(rx_bb,fs,fs/2048)
    plt.plot(f,10*np.log10(p),label='Rx BB')
    
    
    plt.legend(loc="upper right",fontsize=20)
    """