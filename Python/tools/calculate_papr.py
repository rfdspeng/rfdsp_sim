# -*- coding: utf-8 -*-
"""
Created on Sat Sep 17 10:40:14 2022

Calculate PAPR of signal x (complex ndarray)
p = percentile in %. For example, if you wish to calculate 99.99% PAPR, set p = 99.99

@author: Ryan Tsai
"""

import numpy as np

def calculate_papr(x,p):
    env2 = x*np.conjugate(x)
    env2 = env2.real
    env2 = np.sort(env2)
    avg_power = np.mean(env2)
    peak_power = env2[round(np.ceil(env2.size*p/100))]
    papr = 10*np.log10(peak_power/avg_power)
    
    return papr