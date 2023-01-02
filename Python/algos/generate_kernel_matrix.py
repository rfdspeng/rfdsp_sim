# -*- coding: utf-8 -*-
"""
Created on Sun Jan  1 08:27:45 2023

Generate kernel matrix for DPD training

ktups is a list of tuples. Each tuple consists of (identifier string, params)
The identifier string is 'GMP', 'DDR', etc. The params depend on the kernel type.

@author: Ryan Tsai
"""

import numpy as np

def generate_kernel_matrix(x,ktups):
    x = x.reshape(x.size,1)
    x_abs = abs(x)
    
    kmat = np.empty(0)
    kstr = []
    for kdx in range(len(ktups)):
        ktup = ktups[kdx]
        ktype = ktup[0]
    
        if ktype == 'GMP':
            p = ktup[1] # Total order
            m = ktup[2] # Complex delay
            l = ktup[3] # Envelope delay (relative to complex delay)
            
            # Complex term
            if m == 0:
                x_iq = x
            elif m > 0:
                x_iq = np.vstack( (np.zeros((m,1)), x[0:-m]) )
            elif m < 0:
                x_iq = np.vstack( (x[-m:], np.zeros((-m,1))) )
            
            # Envelope term
            if m+l == 0:
                x_env = x_abs
            elif m+l > 0:
                x_env = np.vstack( (np.zeros((m+l,1)), x_abs[0:-(m+l)]) )
            elif m+l < 0:
                x_env = np.vstack( (x_abs[-(m+l):], np.zeros((-(m+l),1))) )
                
            kernel = x_iq*x_env**(p-1)
            kid = 'x' + str(m) if p == 1 else 'x' + str(m) + '|x' + str(l) + '|^' + str(p-1)
        elif ktype == 'DDR':
            'dummy'
        
        kmat = kernel if kdx == 0 else np.hstack((kmat,kernel))
        kstr.append(kid) 
        
    return (kmat,kstr)