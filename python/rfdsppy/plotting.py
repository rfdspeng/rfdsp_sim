# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 08:55:00 2025

Functions for plotting

@author: Ryan Tsai
"""

import numpy as np
import matplotlib.pyplot as plt

def plot_ct_sig(t: np.ndarray, x: np.ndarray, labels: list[str] | None=None, title: str | None=None):
    """
    To plot multiple signals, let t and x have shape nsig * nsamp
    
    """
    t = t.copy()
    x = x.copy()

    if t.ndim == 1:
        t = t[np.newaxis, :]
    if x.ndim == 1:
        x = x[np.newaxis, :]
    
    assert t.shape == x.shape

    labels = labels.copy() if labels is not None else [f"x{i}" for i in range(x.shape[0])]
    assert x.shape[0] == len(labels)

    fig, ax = plt.subplots(dpi=150)
    for idx in range(x.shape[0]):
        ax.plot(t[idx, :], x[idx, :], label=labels[idx]) #, linewidth = 10/2**(idx+1))

    ax.set_xlim(left=0, right=t.max())
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Signal Amplitude")
    ax.grid()
    ax.legend(fontsize=8)
    if title is None:
        title = "CT Signals"
    fig.suptitle(title)

    return (fig, ax)

def plot_ctfs(f: np.ndarray, Xk: np.ndarray, title: str | None=None):
    f = f.copy()
    Xk = Xk.copy()

    fig, axs = plt.subplots(nrows=2, sharex=True, dpi=150)
    axs[0].stem(f, np.abs(Xk))
    axs[0].set_ylabel("|Xk|")
    axs[0].tick_params(labelbottom=True)
    axs[0].grid()
    axs[1].stem(f, np.angle(Xk))
    axs[1].set_xlim(left=f.min(), right=f.max())
    axs[1].set_ylabel("Angle Xk (rad)")
    axs[1].set_xlabel("Frequency (Hz)")
    axs[1].grid()
    if title is None:
        title = "CT Fourier Series"
    fig.suptitle(title)

    return (fig, axs)