# -*- coding: utf-8 -*-
"""
Created on Sun Jun 12 07:22:28 2022

University of Utah
Math 2270
Patrick Dylan Zwick
https://www.math.utah.edu/~zwick/Classes/Fall2012_2270/

Lecture 2

@author: Ryan Tsai
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from IPython import get_ipython
import sys

if __name__ == '__main__':
    plt.close('all')
    get_ipython().magic('reset -sf')
    
    "Dot product example in 2D space"
    "Generate a random vector x"
    "Generate multiple random vectors y by rotating x"
    "Take the dot products of x and y"
    ""
    
    # Generate random 2D vector x
    x1 = np.random.rand()*2-1
    x2 = np.random.rand()*2-1
    x = np.array((x1,x2)) # 2D vector
    
    # Convert x to complex number in order to perform rotation operation
    # Also normalize x to be unit length
    x_iq = x[0] + 1j*x[1]
    x_iq = x_iq/abs(x_iq)
    
    # Rotate x to generate y
    angles = np.linspace(0,2*math.pi,2**4+1)
    angles = angles[0:-1]
    
    rot = np.exp(1j*angles)
    y_iq = x_iq*rot
    y_i = y_iq.real.reshape((16,1))
    y_q = y_iq.imag.reshape((16,1))
    
    # Final x and y vectors
    x = np.array((x_iq.real, x_iq.imag))
    x = x.reshape((2,1))
    y = np.concatenate((y_i,y_q),axis=1)
    
    # Take dot product of x and y
    z = y @ x
    
    # Plot y
    fig = plt.figure()
    for vdx in range(0,y.shape[0]):
        plt.plot(np.array((0,y[vdx,0])),np.array((0,y[vdx,1])))
    plt.title("Y Vectors",{'fontsize':40})
    plt.xlabel("dim1",{'fontsize':30})
    plt.ylabel("dim2",{'fontsize':30})
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim([-1,1])
    plt.ylim([-1,1])
    #plt.autoscale(enable=True,axis='both',tight=True)
    plt.grid()
    
    # Plot dot product vs. angle of rotation
    fig = plt.figure()
    plt.plot(angles/math.pi,z)
    plt.title("Dot Product vs. Angle",{'fontsize':40})
    plt.xlabel("Normalized Angle (Rad/PI)",{'fontsize':30})
    plt.ylabel("Dot Product",{'fontsize':30})
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.autoscale(enable=True,axis='both',tight=True)
    plt.grid()
    
    "For real vectors, the dot product returns a real number that depends on the angle b/w the vectors"
    "For complex vectors, the dot product returns a complex number"