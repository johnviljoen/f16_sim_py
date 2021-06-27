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

so_file = os.getcwd() + "/C/nlplant_xcg35.so"
nlplant = CDLL(so_file)

xu = np.zeros(18)
xdot = np.zeros(18)
fi_flag = 1

nlplant.Nlplant(ctypes.c_void_p(xu.ctypes.data), ctypes.c_void_p(xdot.ctypes.data), ctypes.c_int(fi_flag))

# In[]

# two doubles of alt and vt followed by an array of pointers of [mach, qbar, ps]

coeff = np.zeros(3)

alt = 10000
vt = 700

nlplant.atmos(ctypes.c_double(alt),ctypes.c_double(vt),ctypes.c_void_p(coeff.ctypes.data))

atmos_out = coeff[1]/coeff[2] * 9.05

# In[]

alpha_deg = 3


