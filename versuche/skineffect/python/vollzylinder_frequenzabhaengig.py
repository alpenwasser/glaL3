#!/usr/bin/env python3

from sympy import *
from sympy.external import import_module
#import numpy as numpy
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
sigma = 24e6
r     = 0
r0    = 45e-3
B0    = 6.9e-2                             # adjust this as needed for scaling


# ---------------------------------------------------------#
# Function for magnetic Field B                            #
# ---------------------------------------------------------#
# See formula 21 on p.11 of script for experiment.

# TODO: Switch from w to f!!!!!

#k = lambda w: sqrt((w*mu0*sigma)/2)*(mpc(1,-1))                     # rad/sec
k = lambda w: sqrt((2*np.pi*w*mu0*sigma)/2)*(mpc(1,-1))          # degrees/sec

# Enumerator:
enum  = lambda w: besselj(0,k(w)*r)
denom = lambda w: besselj(0,k(w)*r0)

B = lambda w: enum(w) / denom(w) * B0

B_abs = lambda w: abs(B(w))
B_arg = lambda w: arg(B(w))

# Generate points for omega axis of B, store in Bw
var('npts expufunc n B_w Bw Babs Barg fmax fmin')
npts = 1e3
fmin=1
fmax=250
n = np.linspace(1,npts,npts)
expufunc = np.frompyfunc(exp,1,1)
Bw = 1*expufunc(n*log(fmax-1)/npts)


# Generate B-Field values for that frequency vector
var('Bargufunc Babsufunc Babs')
Babsufunc = np.frompyfunc(B_abs,1,1)
Babs      = Babsufunc(Bw)
Bargufunc = np.frompyfunc(B_arg,1,1)
Barg      = Bargufunc(Bw) # radians
#Barg      = 180/np.pi*Bargufunc(Bw) # degrees


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
var('frequencies voltages phases_degrees phases_rad')
frequencies    = np.array([     1,     5,    10,    15,  20,     30,    40,     60,  80,   100,   120,   160, 200, 250])
phases_degrees = np.array([   5.4,    26,    50,    69,  85,    111,   132,    166, 196,   220,   243,   283, 320, 350])
#phases_degrees = np.array([   5.4,    26,    50,    69,  85,    111,   132,    166, 196-360,   220-360,   243-360,   283-360, 320-360, 350-360])
voltages       = np.array([6.9e-2,6.5e-2,5.7e-2,4.8e-2,4e-2,2.85e-2,2.1e-2,1.25e-2,8e-3,5.4e-3,3.6e-3,1.9e-3,1e-3,6e-4])
phases_rad     = np.pi/180*phases_degrees



# ---------------------------------------------------------#
# Plot the Things                                          #
# ---------------------------------------------------------#
subplot(2,1,1)
plot(Bw,Babs)
scatter(frequencies,voltages)
xscale('log')
subplot(2,1,2)
plot(Bw,-Barg)
scatter(frequencies,phases_degrees)
xscale('log')
show()
