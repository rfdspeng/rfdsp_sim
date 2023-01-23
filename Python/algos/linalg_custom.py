# -*- coding: utf-8 -*-
"""
Created on Sat Jan 21 09:00:29 2023

Custom linear algebra functions

@author: Ryan Tsai
"""

import numpy as np

def identity_matrix(n):
    """
    Generate nxn identity matrix
    
    """
    I = np.zeros((n,n))
    for idx in range(n):
        I[idx,idx] = 1
    
    return I

def lu(A,prune_thr_db):
    """
    A is a square matrix
    prune_thr_db is the pruning threshold in dB (a positive value)
    When a pivot is less than (first pivot)*10^(-prune_thr_db/20), zero out the row and column and move to the next one
    
    Returns L, U, and a vector that contains the indices (starting from 0) of the non-pruned rows and columns
    
    """
    
    if not(len(A.shape) == 2):
        raise Exception("lu: A must be 2-dimensional")
    elif not(A.shape[0] == A.shape[1]):
        raise Exception("lu: A must be a square matrix")
    elif prune_thr_db < 0:
        raise Exception("lu: prune_thr_db must be positive")
        
    ndim = A.shape[0] # Matrix size
    E = identity_matrix(ndim)
    U = A.copy()
    elim_idx = []
    for idx in range(ndim):
        pivot = U[idx,idx]
        if idx == 0:
            first_pivot = abs(pivot)
        
        if abs(pivot) <= first_pivot*10**(-prune_thr_db/20):
            # Pruning condition
            U[idx,:] = 0
            U[:,idx] = 0
            elim_idx.append(idx)
        else:
            # Elimination
            Ei = identity_matrix(ndim)
            for jdx in range(idx+1,ndim):
                c = -U[jdx,idx]/pivot
                U[jdx,:] = c*U[idx,:] + U[jdx,:]
                Ei[jdx,idx] = c
                
            E = Ei @ E
    
    L = np.linalg.inv(E) # Write function to invert triangular matrix
    return (L,U,np.array(elim_idx))

def chol():
    
    
    
    
    
    
    
    
    
    return 1

def solve_ls(A,b,decomp_option="lu",prune_thr_db=999):
    """
    Solve Ax = b
    
    decomp_option: Choose the algorithm to use for decomposition of A
    
    """
    
    if not(len(A.shape) == 2):
        raise Exception("solve_ls: A must be 2-dimensional")
    elif not(len(b.shape) == 2):
        raise Exception("solve_ls: b must be 2-dimensional")
    elif not(b.shape[1] == 1):
        raise Exception("solve_ls: b must have 1 column")
    elif not(A.shape[0] == b.shape[0]):
        raise Exception("solve_ls: A and b must have the same number of rows")
        
    if decomp_option == "lu":
        # Ax = b
        # AHAx = AHb
        AHA = A.T.conj() @ A
        AHb = A.T.conj() @ b
        
        # LUx = AHb
        [L,U,elim_idx] = lu(AHA,prune_thr_db=prune_thr_db)
        
        # Solve Ly = AHb for y, where y = Ux
        y = np.zeros((L.shape[0],1))
        for idx in range(L.shape[0]):
            if not(idx in elim_idx):
                y[idx] = (AHb[idx] - L[idx,:]*y.T)/L[idx,idx]
        
        # Solve Ux = y for x
        
    
    return 1