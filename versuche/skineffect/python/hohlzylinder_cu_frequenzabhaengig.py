#!/usr/bin/env python3

from sympy import *
from sympy.external import import_module
#from numpy import *
from mpmath import *
from matplotlib.pyplot import *
from math import copysign

init_printing()

#mp.dps=27  # precision in decimal places, prec ~= 3.33*dps
mp.prec=80 # precision in bits

# ---------------------------------------------------------#
# Measurement in Center, depending on Frequency            #
# ---------------------------------------------------------#


# ---------------------------------------------------------#
# Init, Define Variables and Constants                     #
# ---------------------------------------------------------#
var('mu0 B_abs B_arg B B0 j0 k r r0 w f sigma denom enum')
mu0   = 4*pi*1e-7
#sigma = 58.1e6                  # de.wikipedia.org/wiki/Kupfer
sigma = 52e6
r     = 30e-3
r0    = 45e-3
r1    = 30e-3
r2    = 35e-3
B0    = 6.9e-2                             # adjust this as needed for scaling


# ---------------------------------------------------------#
# Function for magnetic Field B                            #
# ---------------------------------------------------------#
# See formula 21 on p.11 of script for experiment.

# TODO: Switch from w to f!!!!!

#k = lambda w: sqrt((w*mu0*sigma)/2)*(mpc(1,-1))                     # rad/sec
k = lambda w: sqrt((2*np.pi*w*mu0*sigma)/2)*(mpc(1,-1))          # degrees/sec

# Enumerator:
enum  = lambda w: besselj(0,k(w)*r)  * bessely(2,k(w)*r1) - besselj(2,k(w)*r1) * bessely(0,k(w)*r)
denom = lambda w: besselj(0,k(w)*r2) * bessely(2,k(w)*r1) - besselj(2,k(w)*r1) * bessely(0,k(w)*r2)

B = lambda w: enum(w) / denom(w) * B0

B_abs = lambda w: abs(B(w))
B_arg = lambda w: arg(B(w))

# Generate points for omega axis of B, store in Bw
var('npts expufunc n B_w Bw Babs Barg fmax fmin')
npts = 1e3
fmin=1
#fmax=1943 # maximum calculable frequency
fmax = 2500
n = np.linspace(1,npts,npts)
expufunc = np.frompyfunc(exp,1,1)
Bw = 1*expufunc(n*log(fmax-1)/npts)


# Generate B-Field values for that frequency vector
var('Bargufunc Babsufunc Babs')
Babsufunc = np.frompyfunc(B_abs,1,1)
Babs      = Babsufunc(Bw)
Bargufunc = np.frompyfunc(B_arg,1,1)
#Barg      = Bargufunc(Bw) # radians
Barg      = Bargufunc(Bw) # radians


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
frequencies    = np.array([    1,     10,      20,      40,      80,     120,     160,   200,    400,    600,    800,   1000, 1200, 1500])
phases_degrees = np.array([    2,   19.2,    35.2,    56.7,    76.7,      87,      94,   100,    121,    140,    155,    170,  180,  200])
voltages       = np.array([ 7e-2, 6.6e-2, 5.78e-2, 4.18e-2, 2.44e-2, 1.69e-2, 1.27e-2,  1e-2, 4.8e-3, 2.9e-3, 1.9e-3, 1.4e-3, 1e-3, 7e-4])
phases_rad     = np.pi/180*phases_degrees


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

subplot(2,1,1)
subplot(2,1,1)
plot(Bw,Babs,color='blue',label='Fitfunktion')
scatter(frequencies,voltages,color='black',s=64,label='Messwerte')
xlabel('Frequenz (Hz)',fontdict=font)
ylabel('Spannung (Volt)',fontdict=font)
title('Betrag des Magnetfelds in Zylinderspule mit Hohlzylinder aus Kupfer, (Messpunkt: auf Zylinderachse, horizontal zentriert)',fontdict=font)
legend(fontsize=16)
xscale('log')
subplot(2,1,2)
plot(Bw,Barg,color='blue',label='Fitfunktion')
scatter(frequencies,-phases_degrees,color='black',s=64,label='Messwerte')
xlabel('Frequenz (Hz)',fontdict=font)
ylabel('Phase (Grad))',fontdict=font)
title('Phase des Magnetfelds in Zylinderspule mit Hohlzylinder aus Kupfer (Messpunkt: auf Zylinderachse, horizontal zentriert)',fontdict=font)
legend(fontsize=16)
xscale('log')
show()
