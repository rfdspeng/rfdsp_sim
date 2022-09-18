# -*- coding: utf-8 -*-
"""
Created on Wed May 25 14:46:47 2022

Calculate rms of a complex vector

@author: tsair
"""

import math
import numpy as np

def rms(x):
    return math.sqrt(np.vdot(x,x).real/x.size)