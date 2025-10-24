# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 2025

Functions and classes for modeling RF Tx analog blocks

@author: Ryan Tsai
"""

import numpy as np

class IQUpconverter:
    def __init__(self, theta=0, ep=0):
        """
        theta = IQ phase mismatch in radians
        ep = IQ gain mismatch (linear)

        """
        self.theta = theta
        self.ep = ep
    
    # def fit(self):
        self.LO_I_ = (1 + self.ep/2)*np.exp(1j*self.theta/2)
        self.LO_Q_ = -1j*(1 - self.ep/2)*np.exp(-1j*self.theta/2)
    
    def transform(self, x: np.ndarray) -> np.ndarray:
        return x.real*self.LO_I_ - x.imag*self.LO_Q_