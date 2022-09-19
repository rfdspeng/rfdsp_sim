# -*- coding: utf-8 -*-
"""
Created on Sun Sep 18 18:08:21 2022

Converts dBm to V and vice versa

Assumes dBm is @ RF (-3dB relative to complex baseband model)

x is the value, zo is characteristic impedance

@author: Ryan Tsai
"""

import math
import numpy as np

def power_voltage_conversion(x,unit,zo=50):
    if unit == 'dBm':
        x = x+3 # Complex baseband is +3dB relative to RF
        V = math.sqrt(10**(x/10)*zo*1e-3)
        return V
    elif unit == 'V':
        P = 10*np.log10(x**2/zo/1e-3)
        P = P-3 # RF is -3dB relative to complex baseband
        return P