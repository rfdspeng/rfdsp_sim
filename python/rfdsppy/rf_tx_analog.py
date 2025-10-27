# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 2025

Functions and classes for modeling RF Tx analog blocks

@author: Ryan Tsai
"""

import numpy as np
import matplotlib.pyplot as plt

class IQUpconverter:
    """
    class IQUpconverter

    Lowpass equivalent model of quadrature upconverter
    
    """
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
    
def rapp_saleh_model(x: np.ndarray, cfg: dict | None = None):
    """
    Rapp-Saleh model (baseband equivalent of RF model)

    https://www.mathworks.com/help/comm/ref/memorylessnonlinearity.html
    - Rapp model for AMAM
    - Saleh model for AMPM

    Parameters
    ----------
    x: complex baseband time-domain signal (V)
    cfg: dictionary
        - 'g': gain in the linear region (Rapp)
        - 'smoothness': smoothness factor (Rapp)
        - 'osat': output saturation level (Rapp)
        - 'a': AMPM alpha (Saleh)
        - 'b': AMPM beta (Saleh)
        - 'en_plot'
    
    Returns
    -------
    y: the output of the Rapp-Saleh model

    """
    
    cfg = cfg if cfg else {}

    g = cfg.get("g", 10**(30/20))
    s = 2*cfg.get("smoothness", 2)
    osat = cfg.get("osat", 25)
    a = cfg.get("a", 0.5)
    b = cfg.get("b", 10)
    
    # Apply model to x
    env = np.abs(x)
    ph = np.angle(x)
    env_y = g*env/(1+(g*env/osat)**s)**(1/s)
    ampm = (a*env**2)/(1+b*env**2)
    y = env_y*np.exp(1j*ph)*np.exp(1j*ampm)
    
    if cfg.get('en_plot', False):
        plt.figure()
        plt.plot(env, env_y, '.')
        plt.title('Rapp-Saleh AMAM', {'fontsize':40})
        plt.xlabel("X (V)", {'fontsize':30})
        plt.ylabel("Y (V)", {'fontsize':30})
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.grid()
        
        plt.figure()
        plt.plot(env, ampm*180/np.pi, '.')
        plt.title('Rapp-Saleh AMPM', {'fontsize':40})
        plt.xlabel("X (V)", {'fontsize':30})
        plt.ylabel("AMPM (Degrees)", {'fontsize':30})
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.grid()
    
    return y

class RappSaleh:
    """
    Rapp-Saleh model (lowpass equivalent amplifier model)

    https://www.mathworks.com/help/comm/ref/memorylessnonlinearity.html
    - Rapp model for AMAM
    - Saleh model for AMPM
    
    Constructor parameters
    ----------------------
    cfg: dictionary
        - 'g': gain in the linear region (Rapp)
        - 'smoothness': smoothness factor (Rapp)
        - 'osat': output saturation level (Rapp)
        - 'a': AMPM alpha (Saleh)
        - 'b': AMPM beta (Saleh)

    Methods
    -------
    transform(x)

    """
    
    def __init__(self, cfg: dict | None = None):
        cfg = cfg if cfg else {}

        self.g = cfg.get("g", 10**(30/20))
        self.s = 2*cfg.get("smoothness", 2)
        self.osat = cfg.get("osat", 25)
        self.a = cfg.get("a", 0.5)
        self.b = cfg.get("b", 10)
    
    def transform(self, x: np.ndarray) -> np.ndarray:
        """
        Apply model to input waveform

        Parameters
        ----------
        x: complex baseband time-domain signal (V)

        Returns
        -------
        y: transformed signal (V)

        """
        env = np.abs(x)
        ph = np.angle(x)

        env_y = self.g*env/(1+(self.g*env/self.osat)**self.s)**(1/self.s)
        ampm = (self.a*env**2)/(1+self.b*env**2)
        y = env_y*np.exp(1j*ph)*np.exp(1j*ampm)

        return y