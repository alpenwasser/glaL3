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
w     = 2*np.pi*30                         # frequency was fixed at 30 Hz


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
rmin=0
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
Barg = 180/np.pi*Barg                              # degrees

# ---------------------------------------------------------#
# Measurement Values from the actual experiment            #
# ---------------------------------------------------------#
radii          = 1e-3 * np.array([      0,      5,     10,    15,  20,    25,    30,    35,    40,    45,    50])
voltages       =        np.array([2.86e-2,2.85e-2,2.87e-2,2.9e-2,3e-2,3.3e-2,3.8e-2,4.5e-2,5.4e-2,6.2e-2,3.7e-2])
phases_degrees =        np.array([    111,    109,    104,    94,  81,    65,  48.5,    32,    16,   2.7,     0])
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
