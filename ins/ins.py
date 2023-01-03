#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 15:37:55 2022

@author: raoni
"""

import math
import numpy as np
from scipy.signal import hilbert
from . import utils

def statio_test_funct(x, timeFrequencyRepresentation, numberOfSurrogates, Nh0):
    
    x = x/np.std(x)
    Nx = len(x)
    
    if timeFrequencyRepresentation == "mtfr":
        Nh = int((2*np.round(Nx*Nh0/2)) - 1)
        nFFT = int(np.power(2, np.ceil(np.log2(Nh))))
        dt = np.floor((Nh+1)/8)
        sides = int((Nh+1)/2)
        tt = np.arange(sides, Nx - sides, dt, dtype=int)
        ttred = tt
        Mh = 5
        tm = 5
        opt1 = 1
    else:
        raise Exception(timeFrequencyRepresentation + " Time-frequency representation is not implemented")
    
    opt_dist = 8
    doBS = 0
    fa_rate = 0.05
    Nhist = 20
    
    
    if timeFrequencyRepresentation == "mtfr":
        MSp = tfrsp_hm(hilbert(x), tt, nFFT, Nh, Mh, tm)
        print("aqui")
    
    INS = 1
    return INS

def tfrsp_hm(x, t, nFFT, Nh, M, tm):
    
    h, Dh, tt = hermf(Nh, M, tm)
    
    S = np.zeros((Nh, M, tm))

    for k in range(M):
        spt = tfrsp_h(x, t, nFFT, h[k, :], Dh[k, :])
        S[:, :, k] = spt
    
    return S
    
    
def hermf(N, M, tm):
    dt = 2*tm/(N-1)
    tt = np.linspace(-tm, tm, N)
    g = np.exp(-np.power(tt,2)/2)
    
    P = np.ones((M+1, N))
    P[1,:] = 2*tt
    
    Htemp = np.zeros((M+1, N))
    
    Dh = np.zeros((M,N))
    
    for k in np.arange(2, M+1):
        P[k,:] = np.multiply(2*tt, P[k-1,:]) - 2*(k-1)*P[k-2, :]
        
    for k in np.arange(M+1):
        Htemp[k,:] = np.multiply(P[k, :], g/np.sqrt(np.sqrt(np.pi) * np.power(2, k) * math.gamma(k+1)) * np.sqrt(dt))
        
    h = Htemp[0:M, :]
    
    for k in range(M):
        Dh[k, :] = (np.multiply(tt, Htemp[k, :]) - np.sqrt(2*(k+1))*Htemp[k+1, :])*dt
        
    return h, Dh, tt
        
def tfrsp_h(x, t, nFFT, h, Dh):
    
    xrow = x.shape[0]
    hlength = int(np.floor(nFFT/4))
    hlength = int(hlength + 1 - math.remainder(hlength, 2))
    tcol = len(t)
    hrow = len(h)
    Lh = (hrow-1)/2
    
    if tcol == 1:
        Dt = 1
    else:
        Deltat = t[1:tcol] - t[0:tcol-1]
        Mini = np.min(Deltat)
        Maxi = np.max(Deltat)
        
        if Mini != Maxi:
            raise Exception("The time instants must be regularly sampled.")
        else:
            Dt = Mini
            
    S = np.zeros((nFFT, tcol))
    tf2 = np.zeros((nFFT, tcol))
    tf3 = np.zeros((nFFT, tcol))
    
    for icol in range(tcol):
        ti = t[icol]
        tau = np.arange(-min([np.round(nFFT/2)-1, Lh, ti-1]), min([np.round(nFFT/2)-1, Lh, xrow-ti]) + 1, dtype=int)
        indices = np.remainder(nFFT+tau, nFFT) + 1
        norm_h = np.norm(h[int(Lh) + tau])
        pass
    
    
        
        
    
    
    return 0


def ins(inputData, **kwargs):
    
    vNh0 = [0.015, 0.03, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]
    vINS = np.zeros((len(vNh0), 1))
    numberOfSurrogates = 50
    data, sampleRate = utils.getData(inputData)
    timeFrequencyRepresentation = "mtfr"
    
    for i in range(len(vNh0)):
        vINS[i] = statio_test_funct(data, timeFrequencyRepresentation,
                                    numberOfSurrogates, vNh0[i])
    
    return vINS