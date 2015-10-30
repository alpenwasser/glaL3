#!/usr/bin/env python3

from sympy import *
from sympy.external import import_module
#from numpy import *
from mpmath import *
from matplotlib.pyplot import *
from math import copysign
init_printing()


# ************************************************************************** #
# Self-Inductance L of copper coil with hollow copper cylinder inserted.     #
# ************************************************************************** #

# All values are in standard SI units unless otherwise noted.

# -------------------------------------------------------- #
# Default precision is insufficient, therefore we increase #
# precision.   One  can  increase the  number  of  decimal #
# places or bits, where the number of bits places is ~3.33 #
# times the number of decimal places.                      #
# -------------------------------------------------------- #
#mp.dps=25  # decimal places
mp.prec=80 # precision in bits


# ---------------------------------------------------------#
# Init, Define Variables and Constants                     #
# ---------------------------------------------------------#
var('mu0 B_abs B_arg B B0 j0 k r w f sigma denom enum')
mu0   = 4*pi*1e-7                                        # vacuum permeability
sigma = 52e6                            # de.wikipedia.org/wiki/Kupfer: 58.1e6
r     = 0             # radial position of measurement probe. Centered on axis
dsp   = 98e-3                                               # diameter of coil
rsp   = dsp / 2                                               # radius of coil
r1    = 30e-3                                # inner radius of copper cylinder
r2    = 35e-3                                # outer radius of copper cylinder
B0    = 6.9e-2                             # adjust this as needed for scaling
N0    = 574                                   # number of turns of copper coil
l     = 500e-3                                         # length of copper coil


# ---------------------------------------------------------#
# Function for magnetic Field B                            #
#                                                          #
# See formula 28 on p.15 of script for experiment.         #
#                                                          #
# NOTE: We use  frequency f  instead of  angular frequency #
# omega since that is what we actually set on the function #
# generator.                                               #
# ---------------------------------------------------------#

k = lambda f: sqrt((2*np.pi*f*mu0*sigma)/2)*(mpc(1,-1))          # degrees/sec

enum1 = lambda f: besselj(0,k(f)*r1) * bessely(2,k(f)*r1) - besselj(2,k(f)*r1) * bessely(0,k(f) * r1)
denom1    = lambda f: besselj(0,k(f)*r2) * bessely(2,k(f)*r1) - besselj(2,k(f)*r1) * bessely(0,k(f) * r2)
enum2     = lambda f: r2 * (besselj(1,k(f)*r2) * bessely(2,k(f)*r1) - besselj(2,k(f)*r1) * bessely(1,k(f) * r2)) - r1 * (besselj(1,k(f)*r1) * bessely(2,k(f)*r1) - besselj(2,k(f)*r1) * bessely(1,k(f) * r1))
denom2    = lambda f: besselj(0,k(f)*r2) * bessely(2,k(f)*r1) - besselj(2,k(f)*r1) * bessely(0,k(f) * r2)
term3     = rsp ** 2 - r2**2
prefactor = mu0 * pi * N0**2 / l

phi_norm = lambda f: prefactor * (r1**2 * enum1(f)/denom1(f) + 2/k(f) * enum2(f)/denom2(f) + term3)

L = lambda f: re(phi_norm(f))


# ---------------------------------------------------------#
# Generate points for omega axis                           #
# ---------------------------------------------------------#
npts = 1e3
fmin=1
#fmax=1943 # maximum calculable frequency for default precision
#fmax = 2500
fmax = 2500
n = np.linspace(1,npts,npts)
expufunc = np.frompyfunc(exp,1,1)
frequency_vector = 1*expufunc(n*log(fmax-1)/npts)


# Generate B-Field values for that frequency vector
L_ufunc = np.frompyfunc(L,1,1)
L_num   = L_ufunc(frequency_vector)


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

plot(frequency_vector,L_num,color='blue',label='Fitfunktion')
#scatter(frequencies,L_measured,color='black',s=64,label='Messwerte')
#xlabel('Frequenz (Hz)',fontdict=font)
#ylabel('Spannung (Volt)',fontdict=font)
#title('Betrag des Magnetfelds in Zylinderspule mit Hohlzylinder aus Kupfer, (Messpunkt: auf Zylinderachse, horizontal zentriert)',fontdict=font)
#legend(fontsize=16)
xscale('log')
show()
