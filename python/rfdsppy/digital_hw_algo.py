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

def polyphase_downsampler(x: np.ndarray, b, R, frac_bits: int | float | bool = False):
    """
    Description
    -----------
    Implements cycle-based polyphase downsampler

    Parameters
    ----------
    x : input signal
    b : prototype filter coefficients
    frac_bits : fractional bitwidth for normalization (0/False means floating point)
    R : integer downsampling ratio

    Returns
    -------
    y : output signal

    """
    
    # Generate branches
    branch_len = math.ceil(len(b)/R)
    b_pp = np.zeros((R, branch_len))
    gd = int((len(b)-1)/2)
    b = np.concatenate((b, np.zeros(b_pp.size-len(b)+1)))
    for branch in range(R):
        b_pp[branch,:] = b[branch:-1:R]
    
    reg_in = np.zeros(branch_len*R) + 1j*np.zeros(branch_len*R) # input tapped delay line
    
    y = np.zeros(math.floor(len(x)/R)) + 1j*np.zeros(math.floor(len(x)/R))
    
    # Zero pad the input signal
    x = np.concatenate((np.zeros(len(reg_in)-1), x, np.zeros(gd)))
    
    branch = 0 # branch index
    y_br = np.zeros(R) + 1j*np.zeros(R) # branch outputs
    ydx = 0 # output index
    for x_end in range(len(reg_in)-1, len(x)):
        # Update input tapped delay line
        reg_in[:] = np.flip(x[x_end-len(reg_in)+1:x_end+1])
        x_tapped = reg_in[0:-1:R]
        
        b_br = b_pp[branch,:] # branch coefficients
        
        # Compute branch output
        conv_out = sum(b_br*x_tapped)
        #if frac_bits > 0:
        #    conv_out = round(conv_out/2**frac_bits)
    
        y_br[branch] = conv_out
        
        # Compute output sample
        if branch == 0:
            y[ydx] = sum(y_br)
            if frac_bits > 0:
                # y[ydx] = round(y[ydx]/2**frac_bits)
                y[ydx] = (y[ydx]/2**frac_bits).round()
            ydx += 1
        
        # Update branch
        if branch == 0:
            branch = R-1
        else:
            branch -= 1
        """
        if branch == 2:
            branch = 0
        else:
            branch += 1
        """
        
        if ydx >= len(y):
            break

    return y

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

        self.rot_mat_p_ = [np.array([[1, -1/2**idx],[1/2**idx, 1]]) for idx in range(N)]
        self.rot_mat_m_ = [np.array([[1, 1/2**idx],[-1/2**idx, 1]]) for idx in range(N)]
    
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

            I = I.reshape((1, I.size))
            Q = Q.reshape((1, Q.size))
            IQ = np.vstack((I, Q))
            env, ph = self.vectoring(IQ)

            return (env, ph)

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

            I = env.reshape((1, env.size))
            Q = np.zeros_like(I)
            IQ = np.vstack((I, Q))
            I, Q = self.rotation(IQ)

            return (I, Q)

    def vectoring(self, IQ):
        """
        IQ is a 2 x nsamp matrix
        First row is I
        Second row is Q

        """

        IQ = IQ.copy()
        if self.sim_type == "vectorized":
            for i in range(self.N):
                y_sign = np.sign(IQ[1,:])
                IQp = self.rot_mat_p_ @ IQ
                IQm = self.rot_mat_m_ @ IQ
                IQ[:, y_sign == -1] = IQp[:, y_sign == -1]
                IQ[:, y_sign == +1] = IQm[:, y_sign == +1]
                # IQ[:, y_sign == 0] Do nothing, you've converged

                self.theta_acc_ = self.theta_acc_ + -1*y_sign*self.theta_i_[i]
            
            env = IQ[0, :]
            ph = self.theta_acc_.copy()
        elif self.sim_type == "hardware":
            pass
        
        return (env, ph)
    
    def rotation(self, IQ):
        """
        IQ is a 2 x nsamp matrix
        First row is I (starting value is env)
        Second row is Q (starting value is 0)

        """

        IQ = IQ.copy()
        if self.sim_type == "vectorized":
            for i in range(self.N):
                err_sign = np.sign(self.err_)
                IQp = self.rot_mat_p_ @ IQ
                IQm = self.rot_mat_m_ @ IQ
                IQ[:, err_sign == 1] = IQp[:, err_sign == 1]
                IQ[:, err_sign == -1] = IQm[:, err_sign == -1]
                # IQ[:, err_sign == 0] Do nothing, you've converged

                self.err_ = self.err_ - err_sign*self.theta_i_[i]
                        
            I = IQ[0,:]*self.final_multiply_
            Q = IQ[1,:]*self.final_multiply_
        elif self.sim_type == "hardware":
            pass

        return (I, Q)