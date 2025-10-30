# -*- coding: utf-8 -*-
"""
Created on Thurs Sep 4 09:42 2025

Classes and functions for working with digital signals and systems

@author: Ryan Tsai
"""

import numpy as np

def round(x: float | int | np.ndarray):
    pass

def sat(x: float | int | np.ndarray):
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