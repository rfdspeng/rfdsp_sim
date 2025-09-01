# -*- coding: utf-8 -*-
"""
Created on Sat Aug 30 10:06 2025

Classes and functions for mixed-signal processing

@author: Ryan Tsai
"""

import numpy as np
from scipy.signal import firls, firwin2, lfilter

class DAC:
    def __init__(self, fs_in, fs_eq: float | None=None, en_eq: bool=True, osr: int=5):
        self.fs_in = fs_in # DAC sampling rate
        self.fs_eq = fs_eq
        self.en_eq = en_eq
        self.osr = osr

    def fit(self, n_taps_eq: int | None=None):
        """
        Fitted parameters:
            fs_out_: continuous-time "sampling rate" - sampling rate of the DAC output, which is, in reality, a continuous-time waveform
            zoh_: coefficients of the ZOH "filter" - model for a realizable DAC
            g_zoh_: ZOH frequency response
            eq_taps_: equalizer coefficients - for equalizing ZOH response in the signal bandwidth
            fs_eq_: equalization bandwidth - attempts to equalize up to fs_eq_ frequency
            f_, g_: desired equalization response - passed to filter design function

        """
        self.fs_out_ = self.fs_in*self.osr # CT "sampling rate"
        self.zoh_ = np.ones(self.osr)
        if self.en_eq:
            n_taps_eq = 4096 if n_taps_eq is None else n_taps_eq
            fs_eq = self.fs_eq if self.fs_eq is not None else self.fs_in/2
            f = np.linspace(0, 1, 4096+1) # 0 to fs_in/2
            T = 1/self.fs_in
            g_zoh = T*np.sinc(f/2)
            g = 1/g_zoh # self.fs_in/2*f*T = f/2
            g[f > fs_eq/(self.fs_in/2)] = 0
            self.eq_taps_ = firwin2(n_taps_eq, f, g)
            self.fs_eq_ = fs_eq
            self.f_ = f*self.fs_in/2
            self.g_ = g
            self.g_zoh_ = g_zoh

        return self

    def transform(self, x: np.ndarray):
        assert x.ndim == 1

        y = np.zeros(x.size*self.osr, dtype=x.dtype)
        y[::self.osr] = x
        y = lfilter(self.zoh_, 1, y)
        y = lfilter(self.eq_taps_, 1, y)

        return y