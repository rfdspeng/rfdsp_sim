# -*- coding: utf-8 -*-
"""
Created on Sat May 21 08:56:25 2022

@author: Ryan Tsai
"""

import math
from scipy import signal
from scipy import fft

def calculate_psd(x,fs,rbw):
    nfft = math.ceil(fs/rbw)
    hannwin = signal.windows.hann(nfft)
    f,p = signal.welch(x,fs,hannwin,nperseg=None,noverlap=None,nfft=None,
                       detrend='constant',return_onesided=False,scaling='spectrum')
    f = fft.fftshift(f)
    p = fft.fftshift(p)
    return (p,f)