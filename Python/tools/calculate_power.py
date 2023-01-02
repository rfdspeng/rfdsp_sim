# -*- coding: utf-8 -*-
"""
Created on Sun Sep 18 18:08:21 2022

Calculates the power of complex signal x (power @ RF in dBm)

@author: Ryan Tsai
"""

import math
import numpy as np
from rms import rms

def calculate_power(x,zo=50):
        P = 10*np.log10(rms(x)**2/zo/1e-3)
        P = P-3 # RF is -3dB relative to complex baseband
        return P