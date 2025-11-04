# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 09:55:56 2023

Functions and classes for modeling RF Rx/FBRx digital front end

@author: Ryan Tsai
"""

import math
import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
from scipy import signal
from scipy import fft
from typing import Literal

class NotchFilter:
    def __init__(self, w0, r, sim_type: Literal["vectorized", "hardware"]="vectorized", Nest: int | float | bool=False, bitwidth: int | float | bool = False):
        self.w0 = w0
        self.r = r
        self.sim_type = sim_type
        self.Nest = Nest
        self.bitwidth = bitwidth # Not currently supported

        self.b_ = [1, -np.exp(1j*w0)]
        self.a_ = [1, -r*np.exp(1j*w0)]
    
    def hw_process(self, x: np.ndarray) -> np.ndarray:
        y = np.zeros_like(x)
        reg = 0 + 0j
        # if self.Nest:
        #     acc = 0 + 0j
        #     for idx in range(self.Nest):
        #         acc = x[idx] + acc*np.exp(1j*self.w0)
        #     acc = acc/self.Nest
        #     reg = acc/(1-self.r)
        if self.Nest:
            for idx in range(self.Nest):
                reg = x[idx] + reg*np.exp(1j*self.w0)
            reg = reg/self.Nest/(1-self.r)

        start_idx = 0 if self.Nest == False else self.Nest
        for idx in range(start_idx, len(y)):
            y[idx] = x[idx] + reg*self.r*np.exp(1j*self.w0) + reg*-1*np.exp(1j*self.w0)
            reg = x[idx] + reg*self.r*np.exp(1j*self.w0)
        
        return y
    
    def transform(self, x: np.ndarray) -> np.ndarray:
        if self.sim_type == "vectorized":
            return signal.lfilter(self.b_, self.a_, x)
        elif self.sim_type == "hardware":
            return self.hw_process(x)