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
from power_voltage_conversion import power_voltage_conversion
from rms import rms
import pa_model

def pa_params():
    # Return PA params to target various compression points

if __name__ == '__main__':
    plt.close('all')
    get_ipython().magic('reset -sf')
    
    en_tprecode = 1; modorder = 4
    clipped_papr = 4.5; mpr = 1
    unclipped_papr = 6.5
    p_avg = 18 # dBm @ PA output (RF power) for MPR0
    p_avg = p_avg-mpr
    p_peak = p_avg+clipped_papr
    p_avg = p_peak-unclipped_papr
    v_rms = power_voltage_conversion(p_avg,'dBm')
    pa_gain = 30 # dB
    pa_gain_lin = 10**(pa_gain/20)
    v_rms_in = v_rms/pa_gain_lin
    
    # Generate waveform
    nsym = 14; bw = 20; scs = 15; num_sc = 1200; start_sc = 600-round(num_sc/2)
    ncp = 7; wola = 1
    x,x_standard,cfg_evm = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola)
    wola_len = cfg_evm['wola_len']
    fs = cfg_evm['fs']
    
    x = x/rms(x)*v_rms_in
    
    # Rapp params
    cfg = {'g':pa_gain_lin}
    cfg['smoothness'] = 2
    cfg['osat'] = 4
    
    # Saleh params
    cfg['a'] = 2
    cfg['b'] = 10
    
    cfg['en_plot'] = 1
    
    # PA model
    y = pa_model.rapp_saleh_model(cfg,x)
    #x = x/max(abs(x))
    #y = pa_model.pa_model(x)
    
    # Calculate peak compression
    cfg_comp = {'en_plot':1}
    comp,nlse = calculate_compression(x,y,cfg=cfg_comp)
    print(comp)
    print(nlse)
    
    # Calculate EVM
    cfg_evm['en_plot'] = 1
    evm = ofdm_evm_calculator(cfg_evm,x_standard,y[round(wola_len/2):])
    snr = round(-20*np.log10(evm/100),2)
    print(evm)
    print(snr)
    
    # Calculate ACLR