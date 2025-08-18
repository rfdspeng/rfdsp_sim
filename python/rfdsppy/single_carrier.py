# -*- coding: utf-8 -*-
"""
Created on Tue May 17 08:50:42 2022

Functions for generating a single-carrier waveform

@author: Ryan Tsai
"""

from typing import Literal
import numpy as np
import math

def dsss_wavgen(nsym, bw, modorder, osr=1, seed=None):
    """
    Direct spread-sequence waveform generator
    """
    
    if not seed:
        rng = np.random.default_rng()
    else:
        rng = np.random.default_rng(seed)
    
    fs = bw*osr # sampling rate
    bps = round(math.log2(modorder)) # bits per symbol
    
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