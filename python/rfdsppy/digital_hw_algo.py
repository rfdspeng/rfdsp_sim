# -*- coding: utf-8 -*-
"""
Created on Mon Nov 3 2025

Functions and classes for modeling generic digital hardware algorithms
and operations

@author: Ryan Tsai
"""

import numpy as np

def round(x: float | int | np.ndarray):
    pass

def sat(x: float | int | np.ndarray):
    pass

class RoundSat:
    def __init__(self):
        pass

def upsample(x: np.ndarray, M: int | float) -> np.ndarray:
    """
    Integer upsampling

    Arguments
        x: numpy array
        M: upsampling ratio

    """
    assert x.ndim == 1, "Input vector should be 1-dimensional"

    M = int(M)
    y = np.zeros(x.size*M, dtype=x.dtype)
    y[::M] = x
    
    return y

def downsample(x: np.ndarray, M: int | float) -> np.ndarray:
    """
    Integer downsampling
    
    Arguments
        x: numpy array
        M: downsampling ratio

    """

    assert x.ndim == 1, "Input vector should be 1-dimensional"

    return x[::M].copy()

class CORDIC:
    def __init__(self, niter):
        pass