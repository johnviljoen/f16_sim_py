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

# import progressbar for convenience
import progressbar

# import parameters
from parameters import initial_state_vector, simulation_parameters

# import exit() function for debugging
from sys import exit

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

# initialise Mach, qbar, ps storage
coeff = np.zeros(3)

# initialise LF_state
LF_state = -xu[7] * 180/pi

# In[]

#----------------------------------------------------------------------------#
#---------------------------------Simulate-----------------------------------#
#----------------------------------------------------------------------------#

rng = np.linspace(time_start, time_end, int((time_end-time_start)/time_step))
bar = progressbar.ProgressBar(maxval=len(rng)).start()

# create storage
xu_storage = np.zeros([len(rng),len(xu)])
xdot_storage = np.zeros([len(rng),len(xdot)])

for idx, val in enumerate(rng):
    
    #----------------------------------------#
    #--------------Take Action---------------#
    #----------------------------------------#
    
    T_cmd = 2886.6468
    dstab_cmd = -2.0385
    ail_cmd = -0.087577
    rud_cmd = -0.03877
    
    #----------------------------------------#
    #------------Actuator Models-------------#
    #----------------------------------------#
    
    #--------------Thrust Model--------------#
    # command saturation
    T_cmd = np.clip(T_cmd,1000,19000)
    # rate saturation
    T_err = np.clip(T_cmd - xu[12], -10000, 10000)
    # integrate
    xu[12] += T_err*time_step
    #--------------Dstab Model---------------#
    # command saturation
    dstab_cmd = np.clip(dstab_cmd,-25,25)
    # rate saturation
    dstab_err = np.clip(20.2*(dstab_cmd - xu[13]), -60, 60)
    # integrate
    xu[13] += dstab_err*time_step
    #-------------aileron model--------------#
    # command saturation
    ail_cmd = np.clip(ail_cmd,-21.5,21.5)
    # rate saturation
    ail_err = np.clip(20.2*(ail_cmd - xu[14]), -80, 80)
    # integrate
    xu[14] += ail_err*time_step
    #--------------rudder model--------------#
    # command saturation
    rud_cmd = np.clip(rud_cmd,-30,30)
    # rate saturation
    rud_err = np.clip(20.2*(rud_cmd - xu[15]), -120, 120)
    # integrate
    xu[15] += rud_err*time_step
    #--------leading edge flap model---------#
    # find command from current states
    nlplant.atmos(ctypes.c_double(xu[2]),ctypes.c_double(xu[6]),ctypes.c_void_p(coeff.ctypes.data))
    atmos_out = coeff[1]/coeff[2] * 9.05
    alpha_deg = xu[7]*180/pi
    LF_err = alpha_deg - (LF_state + (2 * alpha_deg))
    LF_state += LF_err*7.25*time_step
    LF_out = (LF_state + (2 * alpha_deg)) * 1.38
    lef_cmd = LF_out + 1.45 - atmos_out
    # command saturation
    lef_cmd = np.clip(lef_cmd,0,25)
    # rate saturation
    lef_err = np.clip((1/0.136) * (lef_cmd - xu[16]),-25,25)
    # integrate
    xu[16] += lef_err*time_step
    
    #----------run nlplant for xdot----------#
    nlplant.Nlplant(ctypes.c_void_p(xu.ctypes.data), ctypes.c_void_p(xdot.ctypes.data), ctypes.c_int(fi_flag))
    
    #----------------------------------------#
    #--------------Integrator----------------#
    #----------------------------------------#
    
    # update xu
    xu[0:12] += xdot[0:12]*time_step
    
    #----------------------------------------#
    #------------Store History---------------#
    #----------------------------------------#
    
    xu_storage[idx,:] = xu
    xdot_storage[idx,:] = xdot

    bar.update(idx)

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


# %%
