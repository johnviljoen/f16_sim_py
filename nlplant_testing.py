#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 17:18:09 2021

@author: johnviljoen
"""

from ctypes import CDLL
import ctypes
import os
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from scipy.integrate import ode

so_file = os.getcwd() + "/C/nlplant_xcg25.so"
nlplant = CDLL(so_file)

xu = np.zeros(18)
xdot = np.zeros(18)
fi_flag = 1

time_step = 0.001
time_start = 0
time_end = 10
xu[2] = 10000
xu[6] = 700

rng = np.linspace(time_start, time_end, int((time_end-time_start)/time_step))
xu_storage = np.zeros([len(rng),len(xu)])
xdot_storage = np.zeros([len(rng),len(xdot)])

def f(t,xu):
    
    
    
    return 

for idx, val in enumerate(rng):

    nlplant.Nlplant(ctypes.c_void_p(xu.ctypes.data), ctypes.c_void_p(xdot.ctypes.data), ctypes.c_int(fi_flag))
    xu[0:11] += xdot[0:11]*time_step
    xu_storage[idx,:] = xu
    xdot_storage[idx,:] = xdot
    print(idx)
    
# In[]
    
fig, axs = plt.subplots(12, 1)
#fig.suptitle('Vertically stacked subplots')
axs[0].plot(rng, xu_storage[:,0])
axs[0].set_ylabel('npos (ft)')

axs[1].plot(rng, xu_storage[:,1])
axs[1].set_ylabel('epos (ft)')

axs[2].plot(rng, xu_storage[:,2])
axs[2].set_ylabel('h (ft)')

axs[3].plot(rng, xu_storage[:,3])
axs[3].set_ylabel('$\phi$ (rad)')

axs[4].plot(rng, xu_storage[:,4])
axs[4].set_ylabel('$\theta$ (rad)')

axs[5].plot(rng, xu_storage[:,5])
axs[5].set_ylabel('$\psi$ (rad)')

axs[6].plot(rng, xu_storage[:,6])
axs[6].set_ylabel("V_t (ft/s)")

axs[7].plot(rng, xu_storage[:,7]*180/pi)
axs[7].set_ylabel('alpha (deg)')

axs[8].plot(rng, xu_storage[:,8]*180/pi)
axs[8].set_ylabel('beta (deg)')

axs[9].plot(rng, xu_storage[:,9]*180/pi)
axs[9].set_ylabel('p (deg/s)')

axs[10].plot(rng, xu_storage[:,10]*180/pi)
axs[10].set_ylabel('q (deg/s)')

axs[11].plot(rng, xu_storage[:,11]*180/pi)
axs[11].set_ylabel('r (deg/s)')
axs[11].set_xlabel('time (s)')

fig2, axs2 = plt.subplots(5,1)

axs2[0].plot(rng, xu_storage[:,12])
axs2[0].set_ylabel('P3')

axs2[1].plot(rng, xu_storage[:,13])
axs2[1].set_ylabel('dh')

axs2[2].plot(rng, xu_storage[:,14])
axs2[2].set_ylabel('da')

axs2[3].plot(rng, xu_storage[:,15])
axs2[3].set_ylabel('dr')

axs2[4].plot(rng, xu_storage[:,16])
axs2[4].set_ylabel('lef')