# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 17:00:56 2022

rapp_saleh_model(): Generic model
pa_model(): Rapp-Saleh model with parameters wrapped in

@author: tsair
"""

import numpy as np

def rapp_saleh_model(cfg,x):
    """
    https://www.mathworks.com/help/comm/ref/memorylessnonlinearity.html
    Rapp model for AMAM
    Saleh model for AMPM

    Rapp parameters
    cfg['g']: Gain in the linear region
    cfg['smoothness']: Smoothness factor
    cfg['osat']: Output saturation level

    Saleh parameters
    cfg['a']: AMPM alpha
    cfg['b']: AMPM beta
    
    x = complex input data (numpy array)
    """
    
    g = cfg['g']; s = 2*cfg['smoothness']; osat = cfg['osat'];
    a = cfg['a']; b = cfg['b'];
    
    env = abs(x)
    ph = np.angle(x)
    env_y = g*env/(1+(g*env/osat)**s)**(1/s)
    ampm = (a*env**2)/(1+b*env**2)
    y = env_y*np.exp(1j*ph)*np.exp(1j*ampm)
    
    return y

def pa_model(x):
    cfg = {}
    
    # Rapp params
    cfg['g'] = 10**(30/20)
    cfg['smoothness'] = 2
    cfg['osat'] = 25
    
    # Saleh params
    cfg['a'] = 0.5
    cfg['b'] = 10
    
    y = rapp_saleh_model(cfg,x)
    
    return y