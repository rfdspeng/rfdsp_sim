# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 2025

Functions and classes for modeling RF Tx digital front end

@author: Ryan Tsai
"""

import numpy as np

class TxIQMC:
    """
    class TxIQMC

    Tx frequency-independent IQ mismatch compensation

    """
    
    def __init__(self, ep, theta):
        """
        ep = gain mismatch (linear)
        theta = phase mismatch (radians)

        IQ mismatch model:
        LO_I = (1 + ep/2)*cos(wc*t + theta/2)
        LO_Q = (1 - ep/2)*sin(wc*t - theta/2)
        """

        self.ep = ep
        self.theta = theta
    
    # def fit(self):
        self.gii_ = 1/((1+self.ep/2)*np.cos(self.theta/2)*(1-np.tan(self.theta/2)**2))
        self.gqi_ = self.gii_*-np.tan(self.theta/2)
        self.gqq_ = 1/((1-self.ep/2)*np.cos(self.theta/2)*(1-np.tan(self.theta/2)**2))
        self.giq_ = self.gqq_*-np.tan(self.theta/2)
        
    def transform(self, x: np.ndarray) -> np.ndarray:
        I = self.gii_*x.real + self.gqi_*x.imag
        Q = self.gqq_*x.imag + self.giq_*x.real

        return I + 1j*Q