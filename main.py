# In[] imports

# from ctypes import *
from ctypes import CDLL
import ctypes
import os

# import numpy and sin, cos for convenience
import numpy as np
from numpy import pi

# import matplotlib for visualisation
import matplotlib.pyplot as plt

# import cython for C library calls
# import cython

# import progressbar for convenience
import progressbar

# import parameters
from parameters import aircraft_properties, initial_state_vector, simulation_parameters

# import utility functions
# from utils import angle2quat, quat2angle, qdot_calc

# import scipy integrator RK45
# from scipy.integrate import RK45

# import exit() function for debugging
from sys import exit

# control systems library for actuator modelling
from control.matlab import *

# In[]

#----------------------------------------------------------------------------#
#-------------------------prepare data for nlplant.c-------------------------#
#----------------------------------------------------------------------------#

# unwrap initial inputs for debugging
npos, epos, h, phi, theta, psi, vt, alpha, beta, P, Q, R, P3, dh, da, dr, lef, fi_flag = initial_state_vector
# Ixx, Iyy, Izz, Ixz, weight, b, S, cbar, He, x_cg, x_cg_ref = aircraft_properties
time_step, time_start, time_end, g, stability_flag = simulation_parameters
fi_flag = 1

# convert inputs to correct units for nlplant.c
m2f = 3.28084 # metres to feet conversion
f2m = 1/m2f # feet to metres conversion
initial_state_vector_ft_rad = np.array([npos*m2f, epos*m2f, h*m2f, phi, theta, psi, vt*m2f, alpha, beta, P, Q, R, P3, dh, da, dr, lef, fi_flag])

# create interface with c shared library .so file in folder "C"
if stability_flag == 1:
    so_file = os.getcwd() + "/C/nlplant_xcg35.so"
elif stability_flag == 0:
    so_file = os.getcwd() + "/C/nlplant_xcg25.so"
    
nlplant = CDLL(so_file)

# initialise xu and xdot
xu = initial_state_vector_ft_rad
xdot = np.zeros(18)

# In[]

#----------------------------------------------------------------------------#
#-------------------------functions for simulation---------------------------#
#----------------------------------------------------------------------------#

def update_xdot(xu, xdot, fi_flag):

    nlplant.Nlplant(ctypes.c_void_p(xu.ctypes.data), ctypes.c_void_p(xdot.ctypes.data), ctypes.c_int(fi_flag))


def update_xu(xu, xdot, time_step):

    xu[0:11] = xu[0:11] + xdot[0:11]*time_step

    # lets upgrade to quaternions!

    return xu

#----------------------------------------------------------------------------#
#---------------------------------Simulate-----------------------------------#
#----------------------------------------------------------------------------#

rng = np.linspace(time_start, time_end, int((time_end-time_start)/time_step))
bar = progressbar.ProgressBar(maxval=len(rng)).start()
# create storage
xu_storage = np.zeros([len(rng),len(xu)])
xdot_storage = np.zeros([len(rng),len(xdot)])

for idx, val in enumerate(rng):
    
    #--------------Take Action---------------#
    T_cmd = 2000.
    dstab_cmd = 0.
    ail_cmd = 0.
    rud_cmd = 0.
    lef_cmd = 0.
    
    #------------Actuator Models-------------#
    
    ### thrust
    # command saturation
    T_cmd = np.clip(T_cmd,1000,19000)
    
    # rate saturation
    T_err = np.clip(T_cmd - xu[12], -10000, 10000)
    
    # integrate
    xu[12] += T_err*time_step
    
    ### Dstab
    # command saturation
    dstab_cmd = np.clip(dstab_cmd,-25,25)
    
    # rate saturation
    dstab_err = np.clip(20.2*(dstab_cmd - xu[13]), -60, 60)
    
    # integrate
    xu[13] += dstab_err*time_step
    
    ### aileron
    # command saturation
    ail_cmd = np.clip(ail_cmd,-21.5,21.5)
    
    # rate saturation
    ail_err = np.clip(20.2*(ail_cmd - xu[14]), -80, 80)
    
    # integrate
    xu[14] += ail_err*time_step
    
    ### rudder
    # command saturation
    rud_cmd = np.clip(rud_cmd,-30,30)
    
    # rate saturation
    rud_err = np.clip(20.2*(rud_cmd - xu[15]), -120, 120)
    
    # integrate
    xu[15] += rud_err*time_step
    
    ### leading edge flap
    # find command from current states
    xu[2] # altitude
    xu[6] # velocity
    
    # command saturation
    lef_cmd = np.clip(lef_cmd,0,25)
    
    # rate saturation
    lef_err = np.clip((1/0.136) * (lef_cmd - xu[16]),-25,25)
    
    # integrate
    xu[16] += lef_err*time_step
        
    #--------------Integrator----------------#
    update_xdot(xu, xdot, fi_flag)
    xu = update_xu(xu, xdot, time_step)

    #------------Store History---------------#
    xu_storage[idx,:] = xu
    xdot_storage[idx,:] = xdot

    bar.update(idx)

#print(xu_storage.shape)

# In[]

#----------------------------------------------------------------------------#
#---------------------------------Visualise----------------------------------#
#----------------------------------------------------------------------------#

#%matplotlib qt

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

axs[7].plot(rng, xu_storage[:,7])
axs[7].set_ylabel('alpha (rad)')

axs[8].plot(rng, xu_storage[:,8])
axs[8].set_ylabel('beta (rad)')

axs[9].plot(rng, xu_storage[:,9])
axs[9].set_ylabel('p (rad/s)')

axs[10].plot(rng, xu_storage[:,10])
axs[10].set_ylabel('q (rad/s)')

axs[11].plot(rng, xu_storage[:,11])
axs[11].set_ylabel('r (rad/s)')
axs[11].set_xlabel('time (s)')



# %%
