# -*- coding: utf-8 -*-
"""
Created on Mon Nov 3 2025

Functions and classes for modeling generic digital hardware algorithms
and operations

@author: Ryan Tsai
"""

import numpy as np
from typing import Literal

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
    """
    CORDIC
        
    """

    def __init__(self, N, mode: Literal["vectoring", "rotation"], sim_type: Literal["vectorized", "hardware"]="vectorized"):
        """
        N = number of iterations
        mode = vectoring or rotation

        """

        self.N = N
        self.mode = mode
        self.sim_type = sim_type

        self.K_ = np.prod(1/np.sqrt(1+2**(-2*np.arange(N))))
        self.theta_i_ = CORDIC.get_rotation_angles(N)

        self.rot_mat_p = [np.array([[1, -1/2**idx],[1/2**idx, 1]]) for idx in range(N)]
        self.rot_mat_m = [np.array([[1, 1/2**idx],[-1/2**idx, 1]]) for idx in range(N)]
    
    @staticmethod
    def get_rotation_angles(N):
        """
        Get the angles of rotation for each iteration (radians)

        """
        return np.arctan(1/2**(np.arange(N)))
    
    def transform(self, I: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """
        In vectoring mode, the input are the rectangular representation I and Q
        In rotation mode, the inputs are the polar representation env and phase
        
        """

        if self.mode == "vectoring":
            I = I.copy()
            Q = Q.copy()
            self.theta_acc_ = np.zeros(I.size)

            self.theta_acc_[I < 0] = np.pi
            Q[I < 0] = -1*Q[I < 0]
            Q[I < 0] = -1*I[I < 0]

            env, ph = self.vectoring(I, Q)

        elif self.mode == "rotation":
            env = I.copy()
            ph = Q.copy()
            
            ph = np.mod(ph, 2*np.pi) # wrap to 0 to 2*pi
            ph[ph >= np.pi] = ph[ph >= np.pi] - 2*np.pi # wrap to -pi, pi
            self.final_multiply_ = np.ones(ph.size)
            self.final_multiply_[(ph < -np.pi/2) and (ph > np.pi/2)] = -1
            ph[ph < -np.pi/2] = ph[ph < -np.pi/2] + np.pi
            ph[ph > np.pi/2] = ph[ph > np.pi/2] - np.pi
            self.err_ = ph

            I, Q = self.rotation(env, ph)

    def vectoring(self, I, Q):
        if self.sim_type == "vectorized":
            for i in range(self.N):
                y_sign = np.sign(Q)
                
    
    def rotation(self, env, ph):
        if self.sim_type == "vectorized":
            pass