#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 28 11:08:06 2022

@author: raoni
"""

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import ins
import time
from scipy.io import wavfile

t = time.time()
A = ins.ins("../bin/fasw0_0.wav")
elapsed = time.time() - t
print("Tempo = " + str(elapsed) + " segundos")

'''
samplerate, dataNoisy = wavfile.read("../bin/ruidoso.wav")
samplerate, dataClean = wavfile.read("../bin/voz.wav")
samplerate, dataOutput = wavfile.read("../bin/output.wav")

stoiNP = stoi(dataClean, dataNoisy, samplerate, extended=False)
stoiIBM = stoi(dataClean, dataOutput, samplerate, extended=False)

print("STOI  NP = " + str(stoiNP))
print("STOI IBM = " + str(stoiIBM))
'''