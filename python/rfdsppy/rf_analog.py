# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 2025

Functions and classes for modeling RF analog blocks and impairments

@author: Ryan Tsai
"""

import numpy as np
import matplotlib.pyplot as plt
from rfdsppy import calc
import scipy.fft
import math

class AWGN:
    """
    class AWGN

    Can generate either real or complex noise
    
    """

    def __init__(self, power, bw, fs):
        """
        power: power in a 50Ohm system (dBm). This is power within the signal bandwidth.
        bw: signal bandwidth (MHz)
        fs: sampling rate (MHz)
        
        If the input signal to transform() is real, then the noise generated is real.
        If the input signal to transform() is complex, then the noise generated is complex, and the power is divided between I and Q.
        
        """

        self.rng = np.random.default_rng()
        self.vrms = calc.dbm2v(power, unit="dBm")*np.sqrt(fs/bw)

    def transform(self, x: np.ndarray) -> np.ndarray:
        if np.iscomplexobj(x):
            n = self.rng.normal(loc=0, scale=self.vrms/np.sqrt(2), size=x.size) + \
                1j*self.rng.normal(loc=0, scale=self.vrms/np.sqrt(2), size=x.size)
        else:
            n = self.rng.normal(loc=0, scale=self.vrms, size=x.size)

        return x + n

class FlickerNoise:
    ""

class PhaseNoise:
    """
    Lowpass equivalent model
    
    """

    def __init__(self, l0, lfloor, bpll, fcorner, fs):
        """
        l0 = inband noise floor (dBc/Hz)
        lfloor = out-of-band noise floor (dBc/Hz)
        bpll = PLL 3dB BW (MHz)
        fcorner = flicker noise corner frequency (MHz)
        fs = sampling rate (MHz)
        
        """
        self.l0 = l0
        # self.l0_lin = 10**(l0/10)
        self.lfloor = lfloor
        # self.lfloor_lin = 10**(lfloor/10)
        self.bpll = bpll
        self.fcorner = fcorner
        self.fs = fs
        self.rng = np.random.default_rng()
    
    def transform(self, x: np.ndarray) -> np.ndarray:
        # One-sided PSD
        # f = np.linspace(0, self.fs/2, x.size+2 + (x.size+1)%2)[1:]
        f = np.linspace(0, self.fs/2, math.ceil(x.size/2)+2 + (math.ceil(x.size/2)+1)%2)[1:]
        l0_lin = 10**(self.l0/10)
        lfloor_lin = 10**(self.lfloor/10)
        fbin = (f[1]-f[0])*1e6 # Hz
        # l0_lin = 10**(self.l0/10)*fbin
        # lfloor_lin = 10**(self.lfloor/10)*fbin
        P = self.bpll**2*l0_lin/(self.bpll**2 + f**2)*(1+self.fcorner/f) + lfloor_lin
        phi = self.rng.uniform(low=-np.pi, high=np.pi, size=P.size)

        # FT
        F = np.concatenate(([0], np.sqrt(fbin*P[:-1]/2), np.sqrt(fbin*P[::-1]/2)))
        phi2 = np.concatenate(([0], np.exp(1j*phi[:-1]), np.exp(-1j*phi[::-1])))
        F = F*phi2*F.size

        # IFFT
        theta = scipy.fft.ifft(F)
        theta = theta.real
        # theta = theta*4343*np.sqrt(1000)

        self.theta_ = theta
        self.P_ = P
        self.Pf_ = f
        self.F_ = F
        self.Ff_ = np.linspace(0, self.fs, F.size)

        return x*np.exp(1j*theta[:x.size])

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

class IP2_IP3_Nonlinearity:
    """
    Receiver model
    
    """

    def __init__(self, IIP2, IIP3):
        """
        IIP2: in dBm (power of one tone)
        IIP3: in dBm (power of one tone)
        
        
        """
        

        self.IIP2 = IIP2
        self.IIP3 = IIP3