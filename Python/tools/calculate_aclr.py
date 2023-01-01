# -*- coding: utf-8 -*-
"""
Created on Sat Dec 31 15:20:39 2022

@author: Ryan Tsai
"""

from calculate_noise import calculate_noise

def calculate_aclr(x,fs,bw,scs,en_plot=0):
    rbw = scs/1000
    nrb = round(bw*5*15/scs)
    obw = nrb*12*scs/1000
    sigl = -obw/2; sigh = obw/2-scs/1000
    sigf = [sigl,sigh]
    noisef = [bw-obw/2,bw+obw/2]
    aclrp = calculate_noise(x,fs,rbw,sigf,noisef,cfg={'en_plot':en_plot,'title':"ACLR+"})
    noisef = [-bw-obw/2,-bw+obw/2]
    aclrm = calculate_noise(x,fs,rbw,sigf,noisef,cfg={'en_plot':en_plot,'title':"ACLR-"})
        
    return (aclrm,aclrp)