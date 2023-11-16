# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 18:39:34 2023

Mimic C matrix multiplication

@author: Ryan Tsai
"""

import math
import numpy as np
from scipy import signal
from scipy import fft
from scipy import linalg
import matplotlib.pyplot as plt
from IPython import get_ipython
import sys
sys.path.append("tools")
sys.path.append("models")

def matmul_custom(A,B):
    # Sanity check dimensions
    if len(A.shape) != 2 or len(B.shape) != 2 or A.shape[1] != B.shape[0]:
        raise Exception('matmul_custom(): A and B dimensions are wrong')
    
    C = np.zeros((A.shape[0],B.shape[1]))
    for b_col_idx in range(B.shape[1]):
        b_col = B[:,b_col_idx]
        
        for a_row_idx in range(A.shape[0]):
            a_row = A[a_row_idx,:]
            
            accum = 0
            for ele_idx in range(len(a_row)):
                accum += b_col[ele_idx]*a_row[ele_idx]
            
            C[a_row_idx,b_col_idx] = accum
    
    return C

if __name__ == '__main__':
    N = 3
    
    A = (np.random.rand(N,N)*10).round()
    B = (np.random.rand(N,N)*10).round()
    
    C = np.matmul(A,B)
    C_builtin = C
    print('Using numpy.matmul(): ')
    print(C)
    
    C = matmul_custom(A,B)
    C_custom = C
    print('Using matmul_custom(): ')
    print(C)
    
    print(sum(C_builtin-C_custom))