# -*- coding: utf-8 -*-
"""
Created on Mon May 23 09:20:17 2022

If you assume all symbols are the same length, then you can also vectorize the symbol slicing
And if you assume no dynamic power change, you can also vectorize the single-tap TD equalizer

The standard FD equalizer and my equalizer are the same - how?

Calculates EVM (%) for OFDM waveform

Input arguments
cfg: dictionary containing demod params
x: reference
y: signal under test
Signals must be at baseband sampling rate and time-aligned

@author: tsair
"""

import math
import numpy as np
from scipy import fft
import matplotlib.pyplot as plt
from rms import rms

def ofdm_evm_calculator(cfg,x,y):
    en_plot = 0 if not('en_plot' in cfg) else cfg['en_plot']
    en_fd_eq = 0 if not('en_fd_eq' in cfg) else cfg['en_fd_eq']
    
    # Convert from TD symbols to FD subcarriers
    n = min(len(x),len(y))
    nfft = cfg['nfft']; ncp = cfg['ncp']; sym_len = nfft+ncp
    fft_start = round(ncp/2)
    tone_idx = cfg['tone_idx']
    tone_idx = tone_idx == 2
    for sdx in range(0,n,sym_len):
        sym_x = x[sdx:sdx+sym_len]
        sym_y = y[sdx:sdx+sym_len]
        
        # Single-tap TD equalizer: Project x onto y
        # y'x / y'y * y
        sym_y = np.vdot(sym_y,sym_x)/np.vdot(sym_y,sym_y)*sym_y
        
        # Remove CP
        sym_x = sym_x[fft_start:fft_start+nfft]
        sym_y = sym_y[fft_start:fft_start+nfft]
        nshift = ncp-fft_start
        sym_x = np.roll(sym_x,-nshift)
        sym_y = np.roll(sym_y,-nshift)
        
        # FFT
        tones_x = fft.fft(sym_x); tones_y = fft.fft(sym_y)
        tones_x = tones_x[tone_idx]; tones_y = tones_y[tone_idx]
        tones_x = tones_x.reshape(1,len(tones_x))
        tones_y = tones_y.reshape(1,len(tones_y))
        if sdx == 0:
            tones_xx = tones_x # nsym x ntone (occupied tones only)
            tones_yy = tones_y
        else:
            tones_xx = np.concatenate((tones_xx,tones_x))
            tones_yy = np.concatenate((tones_yy,tones_y))
    
    # FD equalizer
    fd_eqs = np.ones(tones_x.size) + 1j*np.zeros(tones_x.size)
    if en_fd_eq:
        for edx in range(0,len(fd_eqs)):
            xx = tones_xx[:,edx]; yy = tones_yy[:,edx]
            
            # Standard definition
            fd_eq = np.vdot(xx,xx)/np.vdot(xx,yy)
            tones_yy[:,edx] = fd_eq*yy
            fd_eqs[edx] = fd_eq
            
            # My definition - this also works
            # Project x onto y
            # y'x / y'y * y
            'fd_eq = np.vdot(yy,xx)/np.vdot(yy,yy)'
            
    # Inverse transform precoding
    if cfg['en_tprecode']:
        for sdx in range(0,tones_xx.shape[0]):
            tones_xx[sdx,:] = fft.ifft(tones_xx[sdx,:])
            tones_yy[sdx,:] = fft.ifft(tones_yy[sdx,:])
    
    # Calculate EVM
    err = tones_xx-tones_yy
    errv = err.reshape(1,err.size)
    tonesv = tones_xx.reshape(1,tones_xx.size)
    errp = np.vdot(errv,errv).real
    tonesp = np.vdot(tonesv,tonesv).real
    evm = math.sqrt(errp/tonesp)*100
    
    # Plot constellation
    if en_plot:
        titlestr = cfg['title'] if 'title' in cfg else 'Constellation'
        tonesvx = tones_xx.reshape(1,tones_xx.size)
        tonesvy = tones_yy.reshape(1,tones_yy.size)
        c = rms(tonesvx)
        tonesvx = tonesvx/c; tonesvy = tonesvy/c
        tonesvx = tonesvx[0]; tonesvy = tonesvy[0]
        
        fig = plt.figure()
        plt.plot(tonesvx.real,tonesvx.imag,'x',markersize=50)
        plt.plot(tonesvy.real,tonesvy.imag,'x',markersize=25)
        plt.title(titlestr,{'fontsize':40})
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.grid()
    
    # Return
    return evm
        