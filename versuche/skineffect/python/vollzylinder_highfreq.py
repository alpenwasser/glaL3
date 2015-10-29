#!/usr/bin/env python3

from sympy import *
from sympy.external import import_module
#from numpy import *
from mpmath import *
from matplotlib.pyplot import *
from math import copysign

init_printing()

# ---------------------------------------------------------#
# Measurement in Center, depending on Frequency            #
# ---------------------------------------------------------#


# ---------------------------------------------------------#
# Init, Define Variables and Constants                     #
# ---------------------------------------------------------#
var('mu0 B_abs B_arg B B0 j0 k r r0 w f sigma denom enum')
mu0   = 4*pi*1e-7
#sigma = 37.7e6                 # conductivity of aluminium (de.wikipedia.org)
#sigma = 24e6
sigma = 21.5e6
#r     = 0
r0    = 45e-3
#B0    = 6.9e-2                               # adjust this as needed for scaling
B0    = 6.2e-2                               # adjust this as needed for scaling
w     = 2*np.pi*450                        # frequency was fixed at 30 Hz


# ---------------------------------------------------------#
# Function for magnetic Field B                            #
# ---------------------------------------------------------#
# See formula 21 on p.11 of script for experiment.

#k = lambda w: sqrt((w*mu0*sigma)/2)*(mpc(1,-1))                     # rad/sec
k = lambda r: sqrt((w*mu0*sigma)/2)*(mpc(1,-1))          # degrees/sec

# Enumerator:
enum  = lambda r: besselj(0,k(w)*r)
denom = lambda r: besselj(0,k(w)*r0)

B = lambda r: enum(r) / denom(r) * B0

B_abs = lambda r: abs(B(r))
B_arg = lambda r: arg(B(r))

# Generate points for omega axis of B, store in Bw
var('npts expufunc n B_r Br Babs Barg fmax fmin')
npts = 1e3
rmin=25e-3
#rmax=50e-3
rmax=45e-3
n = np.linspace(rmin,rmax,npts)

# Generate B-Field values for that frequency vector
var('Bargufunc Babsufunc Babs')
Babsufunc = np.frompyfunc(B_abs,1,1)
Babs      = Babsufunc(n)
Bargufunc = np.frompyfunc(B_arg,1,1)
Barg      = Bargufunc(n)


# ---------------------------------------------------------#
# Unfortunately, the  arg() function only  delivers values #
# between -pi and  +pi for the angle of  a complex number, #
# which,  while  correct,  is   not  suitable  for  pretty #
# plotting, so we  will shift the values  larger then zero #
# accordingly for a continuous curve.                      #
# ---------------------------------------------------------#

Barg = np.unwrap(Barg)
Barg = 180/np.pi*Barg-360                          # degrees

# ---------------------------------------------------------#
# Measurement Values from the actual experiment            #
# ---------------------------------------------------------#
radii          = 1e-3 * np.array([    25,  27.5,    30,  32.5,    35,   37.5,    40,  42.5,    45,   47.5,    50])
voltages       =        np.array([1.5e-3,2.2e-3,3.6e-3,5.9e-3,9.5e-3,1.55e-2,2.5e-2,3.9e-2,5.5e-2,5.75e-2,3.8e-2])
phases_degrees =        np.array([   215,   183,   152,   125,   100,     73,    47,    24,   5.2,    0.2,     0])
#phases_rad     = np.pi/180*phases_degrees


# ---------------------------------------------------------#
# Plot the Things                                          #
# ---------------------------------------------------------#
subplot(2,1,1)
plot(n,Babs)
scatter(radii,voltages)
#xscale('log')
subplot(2,1,2)
plot(n,Barg)
scatter(radii,-phases_degrees)
#xscale('log')
show()
