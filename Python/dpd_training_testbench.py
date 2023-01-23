# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 17:34:40 2022

Test bench for generating a "realistic" PA model

Assume APT mode needs to support up to 18dBm without saturating for a max PA bias of 3.5V

@author: tsair
"""

import math
import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
from IPython import get_ipython
import sys
sys.path.append('tools')
sys.path.append('models')
sys.path.append('algos')

import ofdm
import calc
import pa_model
import dpd
import linalg_custom

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

def get_ktups(idx):
    if idx == 0:
        ktups = [('GMP',1,0,0), \
                 ('GMP',3,0,0), \
                 ('GMP',5,0,0), \
                 ('GMP',7,0,0), \
                 ('GMP',9,0,0)  ]
    elif idx == 1:
        ktups = [('GMP',1,0,0), \
                 ('GMP',2,0,0), \
                 ('GMP',3,0,0), \
                 ('GMP',4,0,0), \
                 ('GMP',5,0,0), \
                 ('GMP',6,0,0), \
                 ('GMP',7,0,0), \
                 ('GMP',8,0,0), \
                 ('GMP',9,0,0)  ]
        """  
        ktups = [('GMP',1,0,0), \
                 ('GMP',2,0,0), \
                 ('GMP',3,0,0), \
                 ('GMP',4,0,0), \
                 ('GMP',5,0,0), ]
            """
    
    
    
    
    
    
    
    
    
    
    return ktups

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
    x_dig = x/max(abs(x))
    dig_gain = max(abs(x))

    # PA model
    cfg = get_pa_params(target_comp)
    cfg['en_plot'] = 1
    y = pa_model.rapp_saleh_model(cfg,x)
    y_dig = y/max(abs(y))

    # Kernel matrix
    ktups = get_ktups(0)
    [ykmat,ykstr] = dpd.generate_kernel_matrix(y_dig,ktups)
    
    # LU decomposition
    [L,U,elim_idx] = linalg_custom.lu(ykmat.T.conj() @ ykmat,40)
    
    # Ersatz linear algebra solution
    c = linalg.pinv(ykmat) @ x_dig
    [xkmat,xkstr] = dpd.generate_kernel_matrix(x_dig,ktups)
    x_dpd_dig = np.matmul(xkmat,c)
    x_dpd = x_dpd_dig*dig_gain
    cfg = get_pa_params(target_comp)
    cfg['en_plot'] = 1
    y_dpd = pa_model.rapp_saleh_model(cfg,x_dpd)
    
    xe = np.matmul(ykmat,c)
    error = xe-x_dig
    error2 = sum(abs(error)**2)
    power = sum(abs(x_dig)**2)
    nlse = 10*np.log10(error2/power)
    print('Reverse model NLSE (dB): ' + str(nlse))
    fig = plt.figure()
    plt.plot(x_dig.real,xe.real,'.',label='Real Part')
    plt.plot(x_dig.imag,xe.imag,'.',label='Imag Part')
    plt.title('Reverse Model',{'fontsize':40})
    plt.xlabel('x',{'fontsize':30})
    plt.ylabel('x estimated from y',{'fontsize':30})
    plt.legend(loc='lower right',fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid()
    
    #x = x_dpd; y = y_dpd
    # Plot input and output PSDs
    rbw = scs/1000
    [px,f] = calc.calculate_psd(x_dpd,fs,rbw)
    [py,_] = calc.calculate_psd(y_dpd,fs,rbw)
    px = calc.scale_psd(px,f,bw,scs,start_sc,num_sc)
    py = calc.scale_psd(py,f,bw,scs,start_sc,num_sc)
    fig = plt.figure()
    plt.plot(f,10*np.log10(px),linewidth=2.5,label='PA Input')
    plt.plot(f,10*np.log10(py),label='PA Output')
    plt.title('PSD',{'fontsize':40})
    plt.xlabel('Frequency (MHz)',{'fontsize':30})
    plt.ylabel('PSD (dBm)',{'fontsize':30})
    plt.legend(loc='lower center',fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid()
    
    # Calculate power
    [pout_avg,pout_peak] = calc.calculate_power(y)
    [pout_avg_dpd,pout_peak_dpd] = calc.calculate_power(y_dpd)
    print('Average RF power (dBm): ' + str(round(pout_avg,2)))
    print('Peak RF power (dBm): ' + str(round(pout_peak,2)))
    print('Average RF power with DPD (dBm): ' + str(round(pout_avg_dpd,2)))
    print('Peak RF power with DPD (dBm): ' + str(round(pout_peak_dpd,2)))
    
    # Calculate peak compression
    [comp,nlse] = calc.calculate_compression(x_dpd,y_dpd,cfg={'en_plot':1})
    print('Compression (dB): ' + str(comp))
    print('Forward model NLSE (dB): ' + str(nlse))
    
    # Calculate EVM
    cfg_evm['en_plot'] = 1
    evm = ofdm.ofdm_evm_calculator(cfg_evm,x_standard,y_dpd[round(wola_len/2):])
    snr = round(-20*np.log10(evm/100),2)
    print('EVM (%): ' + str(evm))
    print('SNR (dB): ' + str(snr))
    
    # Calculate ACLR
    [aclrm,aclrp] = calc.calculate_aclr(y_dpd,fs,bw,scs,en_plot=1)
    print('ACLR- (dB): ' + str(aclrm))
    print('ACLR+ (dB): ' + str(aclrp))