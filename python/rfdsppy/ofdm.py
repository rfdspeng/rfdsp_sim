# -*- coding: utf-8 -*-
"""
Created on Tue May 17 08:50:42 2022

@author: Ryan Tsai
"""

import math
import numpy as np
from scipy import signal
from scipy import fft
import matplotlib.pyplot as plt
from rfdsppy import calc
from rfdsppy.digital_modulation import modulation_mapper

def ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,osr=1,ncp=0,wola=0,seed=[]):
    """
    OFDM waveform generator
    """
    
    if seed == []:
        rng = np.random.default_rng()
    else:
        rng = np.random.default_rng(seed)
    
    nrb = round(bw*5*15/scs) # Max number of resource blocks
    nfft = 2**math.ceil(math.log2(nrb*12)) # Natural FFT size (osr=1)
    if osr > 1:
        nfft = nfft*osr
    ncp = round(nfft*ncp/100)
    sym_len = nfft+ncp
    fs = nfft*scs/1000 # Waveform sampling rate
    bps = round(math.log2(modorder))
    wola_len = round(nfft*wola/100)
    wola_len = wola_len + (wola_len % 2) # Needs to be even (half outside of symbol, half inside of symbol)
    b = generate_wola_window(wola_len,sym_len)
    
    x_standard = np.zeros(nsym*sym_len) + 1j*np.zeros(nsym*sym_len) # No WOLA
    x = np.zeros(nsym*sym_len+wola_len) + 1j*np.zeros(nsym*sym_len+wola_len) # With WOLA
    for sdx in range(nsym):
        bitstream = rng.random(bps*num_sc).round() # Bits for current symbol
        tones = modulation_mapper(bitstream,modorder) # Convert from bits to modulated tones
        if en_tprecode: tones = fft.fft(tones) # Transform precoding
        tones_nrb = np.zeros(nrb*12) + 1j*np.zeros(nrb*12) # Available tones
        tones_nrb[start_sc:start_sc+num_sc] = tones # Allocated tones
        tones_all = np.concatenate((np.zeros(round((nfft-nrb*12)/2)), tones_nrb, np.zeros(round((nfft-nrb*12)/2)))) # Freq-domain zero-padding
        tones_all = np.concatenate((tones_all[round(nfft/2):], tones_all[0:round(nfft/2)])) # FFT shift prior to IFFT
        sym = fft.ifft(tones_all) # OFDM symbol
        sym_cp = np.concatenate((sym[nfft-ncp:], sym)) # Add CP
        x_standard[sdx*sym_len:(sdx+1)*sym_len] = sym_cp # Standard waveform
        
        # Apply WOLA
        sym_wola = np.concatenate((sym[round(nfft-ncp-wola_len/2):], sym, sym[0:round(wola_len/2)]))
        sym_wola = np.multiply(sym_wola,b)
        x[sdx*sym_len:(sdx+1)*sym_len+wola_len] = sym_wola + x[sdx*sym_len:(sdx+1)*sym_len+wola_len] # Waveform with WOLA
        
        # Save allocated tone indices for demodulation
        if sdx == 0:
            tone_idx = np.ones(nrb*12) # 1 = Unallocated tones (IBE)
            tone_idx[start_sc:start_sc+num_sc] = 2 # 2 = Allocated tones (EVM, equalizer)
            tone_idx = np.concatenate((np.zeros(round((nfft-nrb*12)/2)), tone_idx, np.zeros(round((nfft-nrb*12)/2)))) # 0 = Zero-padded tones
            tone_idx = np.concatenate((tone_idx[round(nfft/2):], tone_idx[0:round(nfft/2)]))
    
    # Output dictionary for demodulation
    cfg = {}
    cfg['fs'] = fs; cfg['nfft'] = nfft; cfg['ncp'] = ncp; cfg['nrb'] = nrb; cfg['nsym'] = nsym
    cfg['bw'] = bw; cfg['scs'] = scs; cfg['num_sc'] = num_sc; cfg['start_sc'] = start_sc
    cfg['modorder'] = modorder; cfg['en_tprecode'] = en_tprecode; cfg['tone_idx'] = tone_idx
    cfg['wola_len'] = wola_len

    return (x,x_standard,cfg)

def generate_wola_window(wola_len,sym_len):
    # sym_len = nfft+ncp
    
    if wola_len > 0:
        b = signal.hann(2*wola_len)
        b = np.concatenate((b[0:wola_len], np.ones(sym_len-wola_len), b[wola_len:]))
    else:
        b = np.ones(sym_len)
        
    return b

def ofdm_evm_calculator(cfg,x,y):
    """
    If you assume all symbols are the same length, then you can also vectorize the symbol slicing
    And if you assume no dynamic power change, you can also vectorize the single-tap TD equalizer

    The standard FD equalizer and my equalizer are the same - how?

    Calculates EVM (%) for OFDM waveform

    Input arguments
    cfg: dictionary containing demod params
    x: reference
    y: signal under test
    Signals must be at baseband sampling rate and time-aligned
    """
    
    en_plot = 0 if not('en_plot' in cfg) else cfg['en_plot']
    en_fd_eq = 0 if not('en_fd_eq' in cfg) else cfg['en_fd_eq']
    
    # Convert from TD symbols to FD subcarriers
    nfft = cfg['nfft']; ncp = cfg['ncp']; sym_len = nfft+ncp
    fft_start = round(ncp/2)
    tone_idx = cfg['tone_idx']
    tone_idx = tone_idx == 2
    n = min(len(x),len(y))-1 - sym_len # last possible start index (sdx)
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
        c = calc.rms(tonesvx)
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