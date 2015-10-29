#!/usr/bin/env python3

from sympy import *
from sympy.external import import_module
#from numpy import *
from mpmath import *
from matplotlib.pyplot import *
from math import copysign

init_printing()

#mp.dps=27  # precision in decimal places, prec ~= 3.33*dps
#mp.prec=80 # precision in bits
mp.prec=64 # precision in bits

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
dsp   = 98e-3                                               # diameter of coil
rsp   = dsp / 2
r0    = 45e-3
r1    = 30e-3
r2    = 35e-3
B0    = 6.9e-2                             # adjust this as needed for scaling
N0    = 574
l     = 500e-3
R_0   = 6.2                               # resistance of coil


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

enum1     = lambda w: besselj(0,k(w)*r1) * bessely(2,k(w)*r1) - besselj(2,k(w)*r1) * bessely(0,k(w) * r1)
denom1    = lambda w: besselj(0,k(w)*r2) * bessely(2,k(w)*r1) - besselj(2,k(w)*r1) * bessely(0,k(w) * r2)
enum2     = lambda w: r2 * (besselj(1,k(w)*r2) * bessely(2,k(w)*r1) - besselj(2,k(w)*r1) * bessely(1,k(w) * r2)) - r1 * (besselj(1,k(w)*r1) * bessely(2,k(w)*r1) - besselj(2,k(w)*r1) * bessely(1,k(w) * r1))
denom2    = lambda w: besselj(0,k(w)*r2) * bessely(2,k(w)*r1) - besselj(2,k(w)*r1) * bessely(0,k(w) * r2)
term3     = rsp ** 2 - r2**2
prefactor = mu0 * pi * N0**2 / l

phi_norm = lambda w: prefactor * (r1**2 * enum1(w)/denom1(w) + 2/k(w) * enum2(w)/denom2(w) + term3)

R = lambda w: -w * im(phi_norm(w)) + R_0

# ---------------------------------------------------------#
# Generate points for omega axis                           #
# ---------------------------------------------------------#
var('npts expufunc n B_w Bw Babs Barg fmax fmin')
npts = 1e2
fmin=1
#fmax=1943 # maximum calculable frequency
#fmax = 2500
fmax = 1500
n = np.linspace(1,npts,npts)
expufunc = np.frompyfunc(exp,1,1)
frequency_vector = 1*expufunc(n*log(fmax-1)/npts)


# Generate B-Field values for that frequency vector
R_ufunc = np.frompyfunc(R,1,1)
R_num   = R_ufunc(frequency_vector)

# ---------------------------------------------------------#
# Measurement Values from the actual experiment            #
# ---------------------------------------------------------#
frequencies = np.array([    1,     10,      20,      40,      80,     120,     160,   200,    400,    600,    800,   1000, 1200, 1500])
R_measured  = R_ufunc(frequencies)

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
plot(frequency_vector,R_num,color='blue',label='Fitfunktion')
scatter(frequencies,R_measured,color='black',s=64,label='Messwerte')
#xlabel('Frequenz (Hz)',fontdict=font)
#ylabel('Spannung (Volt)',fontdict=font)
#title('Betrag des Magnetfelds in Zylinderspule mit Hohlzylinder aus Kupfer, (Messpunkt: auf Zylinderachse, horizontal zentriert)',fontdict=font)
#legend(fontsize=16)
xscale('log')
#subplot(2,1,2)
#plot(frequency_vector,phi_norm_arg_num,color='blue',label='Fitfunktion')
#scatter(frequencies,phi_norm_arg_measured,color='black',s=64,label='Messwerte')
#xlabel('Frequenz (Hz)',fontdict=font)
#ylabel('Phase (Grad))',fontdict=font)
#title('Phase des Magnetfelds in Zylinderspule mit Hohlzylinder aus Kupfer (Messpunkt: auf Zylinderachse, horizontal zentriert)',fontdict=font)
#legend(fontsize=16)
#xscale('log')
show()