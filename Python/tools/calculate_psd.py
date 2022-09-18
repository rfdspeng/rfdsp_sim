# -*- coding: utf-8 -*-
"""
Created on Sat May 21 08:56:25 2022

@author: Ryan Tsai
"""

import math
from scipy import signal
from scipy import fft

def calculate_psd(x,fs,rbw,wintype='kaiser'):
    nfft = math.ceil(fs/rbw)
    
    if wintype == 'kaiser':
        taps = signal.windows.kaiser(nfft,25)
    elif wintype == 'blackmanharris':
        taps = signal.windows.blackmanharris(nfft)
    elif wintype == 'hann':
        taps = signal.windows.hann(nfft)
        
    f,p = signal.welch(x,fs,taps,nperseg=None,noverlap=None,nfft=None,
                       detrend='constant',return_onesided=False,scaling='spectrum')
    f = fft.fftshift(f)
    p = fft.fftshift(p)
    return (p,f)