# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 17:34:40 2022

Test bench for generating a "realistic" PA model

Assume APT mode needs to support up to 18dBm without saturating for a max PA bias of 3.5V

@author: Ryan Tsai
"""

from pathlib import Path
Path().resolve().parent
import numpy as np
import matplotlib.pyplot as plt
from IPython import get_ipython
import sys
sys.path.append('tools')
sys.path.append('models')
sys.path.append('algos')

import ofdm
import calc
import pa_model

def get_pa_params(target_comp):
    # Return PA params to target various compression points
    # PA gain = 30dB, PA output power = 27dBm
    
    pa_gain = 30 # dB
    pa_gain_lin = 10**(pa_gain/20)
    
    pa_params = {'g':pa_gain_lin}
    pa_params['smoothness'] = 2
    pa_params['a'] = 1
    pa_params['b'] = 30
    
    if target_comp == 1:
        pa_params['osat'] = 13
    elif target_comp == 2:
        pa_params['osat'] = 10
    elif target_comp == 3:
        pa_params['osat'] = 8.5
    else:
        pa_params['osat'] = target_comp
        
    return pa_params

if __name__ == '__main__':
    plt.close('all')
    get_ipython().magic('reset -sf')
    
    en_tprecode = 1; modorder = 4
    clipped_papr = 4.5; mpr = 1
    unclipped_papr = 6.5
    p_avg = 27 # dBm @ PA output (RF power) for MPR0
    p_avg = p_avg-mpr
    p_peak = p_avg+(clipped_papr+3)
    p_avg = p_peak-(unclipped_papr+3)
    v_rms = calc.power_voltage_conversion(p_avg,'dBm')
    pa_gain = 30 # dB
    pa_gain_lin = 10**(pa_gain/20)
    v_rms_in = v_rms/pa_gain_lin
    target_comp = 3
    
    # Generate waveform
    nsym = 14; bw = 20; scs = 15; num_sc = 1200; start_sc = 600-round(num_sc/2)
    ncp = 7; wola = 1; osr = 4; seed = 1;
    [x,x_standard,cfg_evm] = ofdm.ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp=ncp,wola=wola,osr=osr,seed=seed)
    wola_len = cfg_evm['wola_len']
    fs = cfg_evm['fs']
    
    x = x/calc.rms(x)*v_rms_in

    # PA model
    cfg = get_pa_params(target_comp)
    cfg['en_plot'] = 1
    y = pa_model.rapp_saleh_model(cfg,x)
    #y = pa_model.pa_model(x)
    #x = x/max(abs(x))
    #y = pa_model.pa_model(x)
    
    # Plot input and output PSDs
    rbw = scs/1000
    [px,f] = calc.calculate_psd(x,fs,rbw)
    [py,_] = calc.calculate_psd(y,fs,rbw)
    px = calc.scale_psd(px,f,bw,scs,start_sc,num_sc)
    py = calc.scale_psd(py,f,bw,scs,start_sc,num_sc)
    fig = plt.figure()
    plt.plot(f,10*np.log10(px),linewidth=2.5,label='PA Input')
    plt.plot(f,10*np.log10(py),label='PA Output')
    plt.title('PSD',{'fontsize':40})
    plt.xlabel("Frequency (MHz)",{'fontsize':30})
    plt.ylabel("PSD (dBm)",{'fontsize':30})
    plt.legend(loc="lower center",fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid()
    
    # Calculate power
    p_pa = calc.calculate_power(y)
    print('RF power @ PA output (dBm): ' + str(p_pa))
    
    # Calculate peak compression
    [comp,nlse] = calc.calculate_compression(x,y,cfg={'en_plot':1})
    print('Compression (dB): ' + str(comp))
    print('Forward model NLSE (dB): ' + str(nlse))
    
    # Calculate EVM
    cfg_evm['en_plot'] = 1
    evm = ofdm.ofdm_evm_calculator(cfg_evm,x_standard,y[round(wola_len/2):])
    snr = round(-20*np.log10(evm/100),2)
    print('EVM (%): ' + str(evm))
    print('SNR (dB): ' + str(snr))
    
    # Calculate ACLR
    [aclrm,aclrp] = calc.calculate_aclr(y,fs,bw,scs,en_plot=1)
    print('ACLR- (dB): ' + str(aclrm))
    print('ACLR+ (dB): ' + str(aclrp))