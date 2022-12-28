#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 15:45:32 2022

@author: raoni
"""

import numpy as np
from scipy.io import wavfile

def getData(input):
    if isinstance(input, str):
        samplerate, data = wavfile.read(input)
        if np.max(data) > 1:
            data = data/np.power(2,16)
    elif isinstance(input, np.ndarray):
        samplerate = 0
        data = input
    return data, samplerate