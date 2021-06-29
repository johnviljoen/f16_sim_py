#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 27 22:59:50 2021

@author: johnviljoen
"""

from ctypes import CDLL
import ctypes
import os
import numpy as np
import matplotlib.pyplot as plt

so_file = os.getcwd() + "/C/nlplant_xcg35.so"
nlplant = CDLL(so_file)

xu = np.zeros(18)
xdot = np.zeros(18)
fi_flag = 1

alt = 10000
vt = 700
coeff = np.zeros(3)

nlplant.atmos(ctypes.c_double(alt),ctypes.c_double(vt),ctypes.c_void_p(coeff.ctypes.data))
atmos_out = coeff[1]/coeff[2] * 9.05

alpha_deg = 1.0721
LF_state = -alpha_deg
time_step = 0.001
lef = 0.25012
time_start = 0
time_end = 3

rng = np.linspace(time_start, time_end, int((time_end-time_start)/time_step))
lef_storage = np.zeros([len(rng),1])
lef_cmd_storage = np.zeros([len(rng),1])

for idx, val in enumerate(rng):
    
    LF_err = alpha_deg - (LF_state + (2 * alpha_deg))
    
    LF_state += LF_err*7.25*time_step
    
    LF_out = (LF_state + (2 * alpha_deg)) * 1.38
    
    lef_cmd = LF_out + 1.45 - atmos_out
    
    ##############THEPROBLEM##############
    # command saturation
    lef_cmd = np.clip(lef_cmd,0,25)
    
    # rate saturation
    lef_err = np.clip((1/0.136) * (lef_cmd - lef),-25,25)
    
    # integrate
    lef += lef_err*time_step
    ######################################
    lef_storage[idx,:] = lef
    lef_cmd_storage[idx,:] = lef_cmd
    

fig, axs = plt.subplots(2, 1)
axs[0].plot(rng, lef_storage[:,0])
axs[0].set_ylabel('npos (ft)')

axs[1].plot(rng, lef_cmd_storage[:,0])
axs[1].set_ylabel('npos (ft)')