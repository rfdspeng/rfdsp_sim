# -*- coding: utf-8 -*-
"""
Created on Tue May 17 08:50:42 2022

OFDM waveform generator

@author: Ryan Tsai
"""

import math
import numpy as np
from scipy import signal
from scipy import fft

def ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp=0,wola=0,seed=[]):
    if seed == []:
        rng = np.random.default_rng()
    else:
        rng = np.random.default_rng(seed)
    
    nrb = bw*5 # max number of resource blocks
    nfft = 2**math.ceil(math.log2(nrb*12))
    ncp = round(nfft*ncp/100)
    sym_len = nfft+ncp
    fs = nfft*scs/1000 # waveform sampling rate
    bps = round(math.log2(modorder))
    wola_len = round(nfft*wola/100)
    wola_len = wola_len + (wola_len % 2) # needs to be even (half outside of symbol, half inside of symbol)
    b = generate_wola_window(wola_len,sym_len)
    
    x_standard = np.zeros(nsym*sym_len) + 1j*np.zeros(nsym*sym_len) # no WOLA
    x = np.zeros(nsym*sym_len+wola_len) + 1j*np.zeros(nsym*sym_len+wola_len) # with WOLA
    for sdx in range(nsym):
        bitstream = rng.random(bps*num_sc).round() # bits for current symbol
        tones = modulation_mapper(bitstream,modorder) # convert from bits to modulated tones
        if en_tprecode: tones = fft.fft(tones) # transform precoding
        tones_nrb = np.zeros(nrb*12) + 1j*np.zeros(nrb*12) # available tones
        tones_nrb[start_sc:start_sc+num_sc] = tones # allocated tones
        tones_all = np.concatenate((np.zeros(round((nfft-nrb*12)/2)), tones_nrb, np.zeros(round((nfft-nrb*12)/2)))) # freq-domain zero-padding
        tones_all = np.concatenate((tones_all[round(nfft/2):], tones_all[0:round(nfft/2)])) # FFT shift prior to IFFT
        sym = fft.ifft(tones_all) # OFDM symbol
        sym_cp = np.concatenate((sym[nfft-ncp:], sym)) # add CP
        x_standard[sdx*sym_len:(sdx+1)*sym_len] = sym_cp # standard waveform
        
        # Apply WOLA
        sym_wola = np.concatenate((sym[round(nfft-ncp-wola_len/2):], sym, sym[0:round(wola_len/2)]))
        sym_wola = np.multiply(sym_wola,b)
        x[sdx*sym_len:(sdx+1)*sym_len+wola_len] = sym_wola + x[sdx*sym_len:(sdx+1)*sym_len+wola_len] # waveform with WOLA
        
        # Save allocated tone indices for demodulation
        if sdx == 0:
            tone_idx = np.ones(nrb*12) # 1 = unallocated tones (IBE)
            tone_idx[start_sc:start_sc+num_sc] = 2 # 2 = allocated tones (EVM, equalizer)
            tone_idx = np.concatenate((np.zeros(round((nfft-nrb*12)/2)), tone_idx, np.zeros(round((nfft-nrb*12)/2)))) # 0 = zero-padded tones
            tone_idx = np.concatenate((tone_idx[round(nfft/2):], tone_idx[0:round(nfft/2)]))
    
    # Output dictionary for demodulation
    cfg = {}
    cfg['fs'] = fs; cfg['nfft'] = nfft; cfg['ncp'] = ncp; cfg['nrb'] = nrb; cfg['nsym'] = nsym
    cfg['bw'] = bw; cfg['scs'] = scs; cfg['num_sc'] = num_sc; cfg['start_sc'] = start_sc
    cfg['modorder'] = modorder; cfg['en_tprecode'] = en_tprecode; cfg['tone_idx'] = tone_idx
    cfg['wola_len'] = wola_len

    return (x,x_standard,cfg)
            
def modulation_mapper(bitstream,modorder):
    if modorder == 4:
        tones = 1-2*bitstream[0::2] + 1j*(1-2*bitstream[1::2])
        tones = tones/math.sqrt(2)
    elif modorder == 16:
        tones = (1-2*bitstream[0::4]) * (2-(1-2*bitstream[2::4])) + 1j * (1-2*bitstream[1::4]) * (2-(1-2*bitstream[3::4]))
        tones = tones/math.sqrt(10)
    elif modorder == 64:
        tones = (1-2*bitstream[0::6]) * (4-(1-2*bitstream[2::6])*(2-(1-2*bitstream[4::6]))) + 1j * (1-2*bitstream[1::6]) * (4-(1-2*bitstream[3::6])*(2-(1-2*bitstream[5::6])))
        tones = tones/math.sqrt(42)
    elif modorder == 256:
        tones = (1-2*bitstream[0::8]) * (8-(1-2*bitstream[2::8])*(4-(1-2*bitstream[4::8])*(2-(1-2*bitstream[6::8])))) + 1j * (1-2*bitstream[1::8]) * (8-(1-2*bitstream[3::8])*(4-(1-2*bitstream[5::8])*(2-(1-2*bitstream[7::8]))))
        tones = tones/math.sqrt(170)
        
    return tones

def generate_wola_window(wola_len,sym_len):
    # sym_len = nfft+ncp
    
    if wola_len > 0:
        b = signal.hann(2*wola_len)
        b = np.concatenate((b[0:wola_len], np.ones(sym_len-wola_len), b[wola_len:]))
    else:
        b = np.ones(sym_len)
        
    return b
        