# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 17:41:44 2023

Mimic C FIR filter implementation

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

def sample_based_fir_filter(h,x):
    L = len(h)
    x = np.concatenate((x,np.zeros(round((L-1)/2))))
    
    y = np.zeros(len(x))
    for n in range(len(x)):
        accum = 0
        for k in range(L):
            if n-k >= 0:
                accum += h[k]*x[n-k] # y(n) = sum{h(k)x(n-k)}
        
        y[n] = accum
    
    return y

if __name__ == '__main__':
    L = 5
    N = 10
    
    #b = (np.random.rand(L)*10).round()
    
    #b = np.zeros(L)
    #b[round((L-1)/2)] = 1
    
    b = np.zeros(L)
    b[0:math.ceil(L/2)] = (np.random.rand(math.ceil(L/2))*10).round()
    b[math.ceil(L/2):len(b)] = np.flip(b[0:math.ceil(L/2)-1])
    
    x = (np.random.rand(N)*10).round()
    
    y = signal.lfilter(b,1,np.concatenate((x,np.zeros(round((L-1)/2)))))
    y_scipy = y
    print('Using scipy.signal.lfilter():')
    print('x = ' + str(x))
    print('y = ' + str(y))
    
    y = sample_based_fir_filter(b,x)
    y_samp = y
    print('Using sample_based_fir_filter():')
    print('x = ' + str(x))
    print('y = ' + str(y))
    
    print(sum(y_scipy - y_samp))