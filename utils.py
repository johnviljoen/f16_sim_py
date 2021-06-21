#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 16:18:57 2020

@author: johnviljoen
"""

import numpy as np
from numpy import cos, sin, pi

def create_2d_grid_points(x_min, x_max, x_spacing, y_min, y_max, y_spacing):
    x_vals = np.linspace(x_min, x_max, np.int(np.abs(np.abs(x_max-x_min)/x_spacing + 1)))
    y_vals = np.linspace(y_min, y_max, np.int(np.abs(np.abs(y_max-y_min)/y_spacing + 1)))
    points = (x_vals, y_vals)
    return points

def angle2quat(psi, theta, phi, input_units, orientation_order):
    
    if input_units == 'rad':
        pass
        
    elif input_units == 'deg':
        phi, theta, psi = phi*pi/180, theta*pi/180, psi*pi/180
    
    cphi = cos(phi/2)
    cthe = cos(theta/2)
    cpsi = cos(psi/2)
    sphi = sin(phi/2)
    sthe = sin(theta/2)
    spsi = sin(psi/2)
    
    if orientation_order == 'ZYX': # YAW PITCH ROLL quarternion!!!!
        q_a2b = np.array([[cphi*cthe*cpsi + sphi*sthe*spsi],
                          [sphi*cthe*cpsi - cphi*sthe*spsi],
                          [cphi*sthe*cpsi + sphi*cthe*spsi],
                          [cphi*cthe*spsi - sphi*sthe*cpsi]])
        
    return q_a2b

def quat2angle(q, output_units, orientation_order): # room for optimisation!
    
    DCM = np.array([[q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2, 2*(q[1]*q[2] + q[0]*q[3]), 2*(q[1]*q[3] - q[0]*q[2])],
                    [2*(q[1]*q[2] - q[0]*q[3]), q[0]**2 - q[1]**2 +q[2]**2 - q[3]**2, 2*(q[2]*q[3] + q[0]*q[1])],
                    [2*(q[1]*q[3] + q[0]*q[2]), q*(q[2]*q[3] - q[0]*q[1]), q[0]**2 - q[1]**2 -q[2]**2 + q[3]**2]])
    
    phi = np.arctan2(DCM[0,1], DCM[0,0]) 
    theta = np.arcsin(-DCM[0,2])
    psi = np.arctan2(DCM[1,2], DCM[2,2])
    
    if output_units == 'deg':
        phi = phi*180/pi
        theta = theta * 180/pi
        psi = psi*180/pi
        
    return psi, theta, phi
    
# testing SUCCESS!

# psi = 25 * pi/180
# theta = 4 * pi/180
# phi = -20 * pi/180

# #print(psi, theta, phi)
# q = angle2quat(phi, theta, psi, input_units = 'rad', orientation_order='ZYX')

# psi, theta, phi = quat2angle(q, output_units='deg', orientation_order='ZYX')
#print(psi, theta, phi)

def qdot_calc(qt, p, q, r):
    
    omega = np.matrix([[0, -p, -q, -r],
                      [p, 0, r, -q],
                      [q, -r, 0, p],
                      [r, q, -p, 0]])
    
    # lbda = 1 - np.linalg.norm(q)
    
    # drift_correction = np.array([[lbda * qt[0]],
    #                              [lbda * qt[1]],
    #                              [lbda * qt[2]],
    #                              [lbda * qt[3]]])
    
    # qt = qt.reshape([4,])
    # drift_correction = 1-np.linalg.norm(drift_correction)
    # print(drift_correction.transpose().dot(qt))
    
    
    qdot = 0.5 * np.matmul(omega, qt) # + drift_correction.transpose().dot(qt)
    qdot = qdot.reshape([4,1])
    
    # print('omega:')
    # print(omega)
    # print('quaternion:')
    # print(qt)
    # print('quaternion dot:')
    # print(qdot)
    
    return qdot

# qdot = qdot_calc(q, 0.02, 0.01, 0.005)
# #print(qdot)

# print(q)
# q = q + qdot*0.001
# print(q)

# print(psi, theta, phi)
# psi, theta, phi = quat2angle(q, output_units='deg', orientation_order='ZYX')
# print(psi, theta, phi)
    
    