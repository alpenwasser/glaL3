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
#sigma = 24e6
sigma = 23.75e6
r     = 0
dsp   = 98e-3                                              # inner dia of coil
rsp   = dsp / 2
r0    = 45e-3
B0    = 6.9e-2                             # adjust this as needed for scaling
N0    = 574                                         # number of turns for coil
l     = 500e-3                                                # length of coil


# ---------------------------------------------------------#
# Function for magnetic Field B                            #
# ---------------------------------------------------------#
# See formula 21 on p.11 of script for experiment.

# TODO: Switch from w to f!!!!!

#k = lambda w: sqrt((w*mu0*sigma)/2)*(mpc(1,-1))                     # rad/sec
k = lambda w: sqrt((2*np.pi*w*mu0*sigma)/2)*(mpc(1,-1))          # degrees/sec

#LRand = #TODO
LRand = (mu0 * 2 * pi * r0 * (rsp - r0) * N0**2) / l

L = lambda w: (mu0*2*pi*r0*N0**2) / l * re(besselj(1,k(w)*r0) / (k(w) * besselj(0,k(w)*r0))) + LRand


# Generate points for omega axis of B, store in Bw
var('npts expufunc n B_w Bw Babs Barg fmax fmin')
npts = 1e3
fmin=1
fmax=250
n = np.linspace(1,npts,npts)
expufunc         = np.frompyfunc(exp,1,1)
frequency_vector = 1*expufunc(n*log(fmax-1)/npts)


# Generate B-Field values for that frequency vector
var('Bargufunc Babsufunc Babs')
Lufunc = np.frompyfunc(L,1,1)
Lnum   = Lufunc(frequency_vector)                    # numerical results for L

# ---------------------------------------------------------#
# Calculated Values Based on Measurements                  #
# ---------------------------------------------------------#
frequencies    = np.array([     1,     5,    10,    15,  20,     30,    40,     60,  80,   100,   120,   160, 200, 250])
L_measured = Lufunc(frequencies)


# ---------------------------------------------------------#
# Plot the Things                                          #
# ---------------------------------------------------------#
font = {
        #'family' : 'monospace',
        'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 16,
        }

#subplot(2,1,1)
plot(frequency_vector,Lnum,color='blue')
scatter(frequencies,L_measured,color='black',s=64,label='Messwerte')
xscale('log')
#xlabel('Frequenz (Hertz)',fontdict=font)
#ylabel('Spannung (Volt)',fontdict=font)
#title('Betrag des Magnetfelds in Zylinderspule mit Vollzylinder aus Aluminium, (Messpunkt: Zylinderachse, horizontal zentriert)',fontdict=font)
#legend(fontsize=16)
#subplot(2,1,2)
#plot(Bw,-Barg,color='blue',label='Fitfunktion')
#xscale('log')
#xlabel('Frequenz (Hertz)',fontdict=font)
#ylabel('Phase (Grad))',fontdict=font)
#title('Phase des Magnetfelds in Zylinderspule mit Vollzylinder aus Aluminium (Messpunkt: Zylinderachse, horizontal zentriert)',fontdict=font)
#legend(fontsize=16)
show()
