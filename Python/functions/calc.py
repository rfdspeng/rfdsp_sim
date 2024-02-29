# -*- coding: utf-8 -*-
"""
Created on Wed May 25 14:46:47 2022

@author: tsair
"""

import math
import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
from scipy import signal
from scipy import fft

def aclr(x,fs,bw,scs,en_plot=0):
    rbw = scs/1000
    nrb = round(bw*5*15/scs)
    obw = nrb*12*scs/1000
    sigl = -obw/2; sigh = obw/2-scs/1000
    sigf = [sigl,sigh]
    noisef = [bw-obw/2,bw+obw/2]
    aclrp = noise_dbc(x,fs,rbw,sigf,noisef,cfg={'en_plot':en_plot,'title':"ACLR+"})
    noisef = [-bw-obw/2,-bw+obw/2]
    aclrm = noise_dbc(x,fs,rbw,sigf,noisef,cfg={'en_plot':en_plot,'title':"ACLR-"})
        
    return (aclrm,aclrp)

def comp_db(x,y,cfg={}):
    """
    Calculate compression by using polynomial curve fitting

    x and y are input and output signals (numpy array)
    cfg['polyorders'] to specify the polynomial terms, e.g. [1,3,5,7,9] for odd orders up to 9. By default, 1 through 9.
    """
    
    x = x.reshape(x.size,1)
    y = y.reshape(y.size,1)
    
    x = x-np.mean(x)
    y = y-np.mean(y)
    
    x = abs(x)/max(abs(x))
    y = abs(y)/max(abs(y))
    
    # Generate kernel matrix, poly order 1 through 9
    X = np.empty(0)
    ks = cfg['polyorders'] if 'polyorders' in cfg else range(1,10)
    for k in ks:
        X = x**k if k == 1 else np.hstack((X,x**k))
        
    c = np.matmul(linalg.pinv(X),y)
    
    xe = np.linspace(0,1,2^16+1)
    ye = np.zeros_like(xe)
    ye2 = np.zeros_like(x)
    for k in range(1,10):
        ye = ye+c[k-1]*xe**k
        ye2 = ye2+c[k-1]*x**k # To estimate LSE
        
    g = ye[1:]/xe[1:]
    comp = 20*np.log10(g[0]/g[-1])
    
    error = ye2-y
    error2 = sum(error*error)
    power = sum(y*y)
    nlse = 10*np.log10(error2/power)
    nlse = nlse[0]
    
    en_plot = cfg['en_plot'] if 'en_plot' in cfg else 0
    if en_plot:
        fig = plt.figure()
        titlestr = cfg['amam_fitting_title'] if 'amam_fitting_title' in cfg else 'AMAM'
        plt.plot(x,y,'.',label='Data')
        plt.plot(xe,ye,label='Fitted')
        plt.title(titlestr,{'fontsize':40})
        plt.xlabel("X (Normalized)",{'fontsize':30})
        plt.ylabel("Y (Normalized)",{'fontsize':30})
        plt.legend(loc="lower right",fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.grid()
        
        fig = plt.figure()
        titlestr = cfg['gain_title'] if 'gain_title' in cfg else 'Gain'
        plt.plot(x,20*np.log10(y/x),'.',label='Data')
        plt.plot(xe[1:],20*np.log10(g),label='Fitted')
        plt.title(titlestr,{'fontsize':40})
        plt.xlabel("X (Normalized)",{'fontsize':30})
        plt.ylabel("Gain (dB)",{'fontsize':30})
        plt.legend(loc="lower left",fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.grid()
        
    return (comp,nlse)

def noise_dbc(x,fs,rbw,sigf,noisef,cfg={}):
    wintype = cfg['wintype'] if 'wintype' in cfg else 'kaiser'
    [p,f] = psd(x,fs,rbw,wintype)
    
    sigfl = sigf[0]; sigfh = sigf[1]
    psig = p[(f >= sigfl) & (f <= sigfh)]
    fsig = f[(f >= sigfl) & (f <= sigfh)]
    
    noisefl = noisef[0]; noisefh = noisef[1]
    pnoise = p[(f >= noisefl) & (f <= noisefh)]
    fnoise = f[(f >= noisefl) & (f <= noisefh)]
    
    psig_sum = sum(psig)
    pnoise_sum = sum(pnoise)
    noise_dbc = 10*np.log10(pnoise_sum/psig_sum)
    
    en_plot = cfg['en_plot'] if 'en_plot' in cfg else 0
    if en_plot:
        fig = plt.figure()
        titlestr = cfg['title'] if 'title' in cfg else 'PSD'
        plt.plot(f,10*np.log10(p),linewidth=2.5,label='Full PSD')
        plt.plot(fsig,10*np.log10(psig),label='Signal')
        plt.plot(fnoise,10*np.log10(pnoise),label='Noise')
        plt.title(titlestr,{'fontsize':40})
        plt.xlabel("Frequency (MHz)",{'fontsize':30})
        plt.ylabel("PSD (dBm)",{'fontsize':30})
        plt.legend(loc="lower center",fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.grid()
    
    return noise_dbc

def papr(x,p=99.99):
    """
    Calculate PAPR of signal x (complex ndarray)
    p = percentile in %. For example, if you wish to calculate 99.99% PAPR, set p = 99.99
    """

    env2 = x*np.conjugate(x)
    env2 = env2.real
    env2 = np.sort(env2)
    avg_power = np.mean(env2)
    peak_power = env2[round(np.ceil(env2.size*p/100))]
    papr = 10*np.log10(peak_power/avg_power)
    
    return papr

def power_dbm(x,zo=50):
    """
    Calculates the power of complex signal x (power @ RF in dBm)
    """
    
    # Average power
    p_avg = 10*np.log10(rms(x)**2/zo/1e-3)
    p_avg = p_avg-10*np.log10(2) # RF is -3dB relative to complex baseband
    
    # Peak power
    p_peak = 10*np.log10(max(abs(x))**2/zo/1e-3)
    return (p_avg,p_peak)
    
def psd(x,fs,rbw,wintype='kaiser'):
    nfft = math.ceil(fs/rbw)
    
    if wintype == 'kaiser':
        taps = signal.windows.kaiser(nfft,25)
    elif wintype == 'blackmanharris':
        taps = signal.windows.blackmanharris(nfft)
    elif wintype == 'hann':
        taps = signal.windows.hann(nfft)
        
    #f,p = signal.welch(x,fs,taps,nperseg=None,noverlap=None,nfft=None,
    #                   detrend='constant',return_onesided=False,scaling='spectrum')
    #f,p = signal.welch(x,fs,taps,nperseg=None,noverlap=None,nfft=None,
    #                   detrend=False,return_onesided=False,scaling='spectrum')
    f,p = signal.welch(x,fs,taps,nperseg=None,noverlap=None,nfft=None,
                       detrend=False,return_onesided=False)
    f = fft.fftshift(f)
    p = fft.fftshift(p)
    return (p,f)

def psd_dbm(x,fs,rbw,low,high,zo=50):
    """
    Calculates PSD of x and then computes power from PSD
    fs and rbw are in MHz
    low and high are the frequency integration limits in MHz
    zo is impedance of the system
    """
    
    # For accurate PSD (scaling='density', see below), fs must be in Hz. Convert everything from MHz to Hz.
    fs = fs*1e6; rbw = rbw*1e6; low = low*1e6; high = high*1e6
    nfft = math.ceil(fs/rbw)
    taps = signal.windows.kaiser(nfft,25)
    
    # Using scaling='spectrum' seems to lead to inaccurate power estimation
    # scaling='density' returns PSD in V^2/Hz if fs is in Hz
    f,p = signal.welch(x,fs,taps,nperseg=None,noverlap=None,nfft=None,
                       detrend=False,return_onesided=False,scaling='density')
    f = fft.fftshift(f)
    p = fft.fftshift(p)
    
    psig = p[(f >= low) & (f <= high)] # integration limits
    fsig = f[(f >= low) & (f <= high)] # for debug only
    p_lin = sum(psig*rbw) # convert from V^2/Hz to V^2
    p_dbm = 10*np.log10(p_lin/zo/1e-3)
    
    return p_dbm
    

def dbm2v(x,unit,zo=50):
    """
    Converts dBm to V and vice versa

    Assumes dBm is @ RF (-3dB relative to complex baseband model)

    x is the value, zo is characteristic impedance
    
    
    xrf = r*cos(wt + phi)
    xbb = r*exp(1j*phi)
    
    
    
    """

    if unit == 'dBm':
        x = x+10*np.log10(2) # Complex baseband is +3dB relative to RF
        V = math.sqrt(10**(x/10)*zo*1e-3)
        return V
    elif unit == 'V':
        P = 10*np.log10(x**2/zo/1e-3)
        P = P-10*np.log10(2) # RF is -3dB relative to complex baseband
        return P
    
def rms(x):
    """
    Calculate rms of a complex vector
    """

    return math.sqrt(np.vdot(x,x).real/x.size)

def scale_psd(p,f,bw,scs,start_sc,num_sc):
    """
    Scale PSD so average signal bin power is 0dBm
    """

    nrb = bw*5
    sigl = -nrb*12*scs/1000/2 + start_sc*scs/1000
    sigh = sigl + (num_sc-1)*scs/1000
    psig = p[(f >= sigl) & (f <= sigh)]
    psig = sum(psig)/len(psig)
    p = p/psig
    return p