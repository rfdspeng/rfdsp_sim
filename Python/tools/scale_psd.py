# -*- coding: utf-8 -*-
"""
Created on Sat May 28 09:01:35 2022

Scale PSD so average signal bin power is 0dBm

@author: tsair
"""

import numpy as np

def scale_psd(p,f,bw,scs,start_sc,num_sc):
    nrb = bw*5
    sigl = -nrb*12*scs/1000/2 + start_sc*scs/1000
    sigh = sigl + (num_sc-1)*scs/1000
    psig = p[np.logical_and(f >= sigl, f <= sigh)]
    psig = sum(psig)/len(psig)
    p = p/psig
    return p