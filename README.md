# Description

This is a simulation of an F16 based upon two works:

1. Nguyen "Simulator study of stall/post-stall characteristics of a fighter airplane with relaxed longitudinal static stability" - https://core.ac.uk/download/pdf/42866809.pdf

2. Stevens and Lewis "Aircraft Control and Simulation: Dynamics, Controls Design, and Autonomous Systems, 3rd Edition" - https://www.wiley.com/en-gb/Aircraft+Control+and+Simulation%3A+Dynamics%2C+Controls+Design%2C+and+Autonomous+Systems%2C+3rd+Edition-p-9781118870983

The former is the higher fidelity model, with accurate angles of attack between -20 -> +90 degrees, thanks to some crazy pilots and a lot of empirical testing. The latter is a reduced version of this model correct between angles of attacked between -10 -> +45 degrees. Both are accurate in a sideslip range between +-30 degrees.

This work is based upon the C implementation of both of these codes for MATLAB found here "https://dept.aem.umn.edu/~balas/darpa_sec/SEC.Software.html#F16Manual". This was originally designed to be a MEX file for MATLAB compilation and use, but I have cut those bits and am now using it as a shared object (.so) file to be called from Python.

This C implementation produces the instantaneous time derivatives of the 12 primary aircraft states (p5 of this document describes the IO: https://dept.aem.umn.edu/~balas/darpa_sec/software/F16Manual.pdf) these are:

x = {npos,epos,h,phi,theta,psi,Vt,alpha,beta,p,q,r}

Therefore it does not complete the simulation, an integrator and an attitude determination system is required. Due to my desire to use this simulation for reinforcement learning I have kept the integrator as a simple time step to ensure every agent interaction occurs for the same amount of time. As for the attitude determination this has been kept as a simple euler angle formulation for now as the nonlinear manouvres which shall be investigated may not encounter singularities which would require a quaternion formulation. Therefore this has been left to future work.

# Technicalities

The C .so file has been compiled twice for two primary configurations of the F16 -> stable and unstable. These occur when xcg = 25, and 35 respectively. Should other parameters be changed the flags for gcc compilation of "nlplant.c" were:

gcc -fPIC -shared -lm -o nlplant.so nlplant.c

These two files are selected in the Python main.py using the "stability_flag" which is 1 for unstable and 0 for stable.

main.py calls the .so file using ctypes.CDLL. This allows access for individual functions in nlplant.so, two of which are called by main.py -> Nlplant and atmos. Nlplant is the main aforementioned one which generates the state derivates, and atmos is used for calculating the leading edge flap (lef) deflection, which is a function of a number of variables.


