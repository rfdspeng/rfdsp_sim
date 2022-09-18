# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 17:35:53 2022

Calculate compression by using polynomial curve fitting

cfg['polyorders'] to specify the polynomial terms, e.g. [1,3,5,7,9] for odd orders up to 9. By default, 1 through 9.

@author: tsair
"""

import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt

def calculate_compression(x,y,cfg={}):
    # x and y are input and output signals (numpy array)
    
    x = x.reshape(x.size,1)
    y = y.reshape(y.size,1)
    
    x = x-np.mean(x)
    y = y-np.mean(y)
    
    x = abs(x)/max(abs(x))
    y = abs(y)/max(abs(y))
    
    # Generate kernel matrix, poly order 1 through 9
    X = np.empty(0)
    ks = cfg['polyorders'] if 'polyorders' in cfg else range(1,10)
    for k in ks:
        X = x**k if k == 1 else np.hstack((X,x**k))
        
    c = np.matmul(linalg.pinv(X),y)
    
    xe = np.linspace(0,1,2^16+1)
    ye = np.zeros_like(xe)
    ye2 = np.zeros_like(x)
    for k in range(1,10):
        ye = ye+c[k-1]*xe**k
        ye2 = ye2+c[k-1]*x**k # To estimate LSE
        
    g = ye[1:]/xe[1:]
    comp = 20*np.log10(g[0]/g[-1])
    
    error = ye2-y
    error2 = sum(error*error)
    power = sum(y*y)
    nlse = 10*np.log10(error2/power)
    nlse = nlse[0]
    
    en_plot = cfg['en_plot'] if 'en_plot' in cfg else 0
    if en_plot:
        fig = plt.figure()
        titlestr = cfg['amam_fitting_title'] if 'amam_fitting_title' in cfg else 'AMAM'
        plt.plot(x,y,'.',label='Data')
        plt.plot(xe,ye,label='Fitted')
        plt.title(titlestr,{'fontsize':40})
        plt.xlabel("X (Normalized)",{'fontsize':30})
        plt.ylabel("Y (Normalized)",{'fontsize':30})
        plt.legend(loc="lower right",fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.grid()
        
        fig = plt.figure()
        titlestr = cfg['gain_title'] if 'gain_title' in cfg else 'Gain'
        plt.plot(x,20*np.log10(y/x),'.',label='Data')
        plt.plot(xe[1:],20*np.log10(g),label='Fitted')
        plt.title(titlestr,{'fontsize':40})
        plt.xlabel("X (Normalized)",{'fontsize':30})
        plt.ylabel("Gain (dB)",{'fontsize':30})
        plt.legend(loc="lower left",fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.grid()
        
    return (comp,nlse)