#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 21:37:08 2021

@author: johnviljoen
"""

from scipy.integrate import ode

y0, t0 = [1.0j, 2.0], 0

def f(t, y):
    return [1j*y[0] + y[1], -y[1]**2]

def jac(t, y):
    return [[1j, 1], [0, -2*y[1]]]

r = ode(f, jac).set_integrator('zvode', method='bdf')
r.set_initial_value(y0, t0).set_f_params(2.0).set_jac_params(2.0)
t1 = 10
dt = 1

while r.successful() and r.t < t1:
    print(r.t+dt, r.integrate(r.t+dt))