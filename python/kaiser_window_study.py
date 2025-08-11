# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 11:46:33 2022

https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.kaiser.html

@author: Ryan Tsai
"""

import math
import numpy as np
from scipy import signal
from scipy import fft
import matplotlib.pyplot as plt
from IPython import get_ipython
import sys
sys.path.append("tools")
sys.path.append("models")

if __name__ == '__main__':
    plt.close('all')
    get_ipython().magic('reset -sf')
    
    winlen = 2**12
    fftmult = 2**6
    f = np.linspace(-1,1,winlen*fftmult+1); f = f[0:-1]
    
    # 0 = rectangular
    # 5 = Hamming-like
    # 6 = Hann-like
    # 8.6 = Blackman-like
    # 14 = recommended as a starting point by SciPy documentation
    # 25 = for PSD estimation
    betas = [0, 5, 6, 8.6, 14, 25]
    betas_str = list(map(str,betas))
    betas_str = ["beta="+beta for beta in betas_str]
    
    fig_td = plt.figure()
    fig_fd = plt.figure()
    for bdx in range(0,len(betas)):
        beta = betas[bdx]
        win = signal.windows.kaiser(winlen,beta)
        winfft = fft.fft(win,n=winlen*fftmult)
        winfft = abs(fft.fftshift(winfft))
        
        plt.figure(fig_td)
        plt.plot(win,label=betas_str[bdx])
        plt.figure(fig_fd)
        plt.plot(f,20*np.log10(winfft/max(winfft)),label=betas_str[bdx])
        
    plt.figure(fig_td)
    plt.title("Kaiser Window",{'fontsize':40})
    plt.xlabel("Sample Index",{'fontsize':30})
    plt.legend(loc="upper right",fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.autoscale(enable=True,axis='both',tight=True)
    plt.grid()
    
    plt.figure(fig_fd)
    plt.title("Kaiser Window FFT in Log Domain",{'fontsize':40})
    plt.xlabel("Normalized Digital Frequency",{'fontsize':30})
    plt.ylabel("FFT (dB)",{'fontsize':30})
    plt.legend(loc="upper right",fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim([0, 0.015])
    plt.ylim([-300, 0])
    #plt.autoscale(enable=True,axis='both',tight=True)
    plt.grid()