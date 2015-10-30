#!/usr/bin/env python3

from sympy import *
from sympy.external import import_module
from mpmath import *
from matplotlib.pyplot import *
#init_printing()     # make things prettier when we print stuff for debugging.


# ************************************************************************** #
# Magnetic Flow normed  by current, copper coil with  hollow copper cylinder #
# inserted.                                                                  #
# ************************************************************************** #

# All values are in standard SI units unless otherwise noted.

# -------------------------------------------------------- #
# Default precision is insufficient, therefore we increase #
# precision.   One  can  increase the  number  of  decimal #
# places or bits, where the number of bits places is ~3.33 #
# times the number of decimal places.                      #
# The highest calculable  frequency with default precision #
# was determined to be 1943 Hz                             #
# -------------------------------------------------------- #
#mp.dps=25  # decimal places
mp.prec=80 # precision in bits

# ---------------------------------------------------------#
# Init, Define Variables and Constants                     #
# ---------------------------------------------------------#
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

    # -----------------------------------------------------#
    # NOTE: According to  formula 26 on p.14,  the B-Field #
    # inside the  copper cylinder  (r<r1) is equal  to the #
    # B-Field at the inner boundary of the copper cylinder #
    # (B(r1)),  therefore  we  set  r to  r1  for  further #
    # calculations.                                        #
    # -----------------------------------------------------#
r     = 30e-3


# ---------------------------------------------------------#
# Functions                                                #
#                                                          #
# See formula 26 on p.14 of script for experiment.         #
#                                                          #
# NOTE: We use  frequency f  instead of  angular frequency #
# omega since that is what we actually set on the function #
# generator.                                               #
# ---------------------------------------------------------#
var('f')

k = lambda f: sqrt((2*np.pi*f*mu0*sigma)/2)*(mpc(1,-1))

enum  = lambda f:(
          besselj(0,k(f)*r)
        * bessely(2,k(f)*r1)
        - besselj(2,k(f)*r1)
        * bessely(0,k(f)*r)
    )
denom = lambda f:(
          besselj(0,k(f)*r2)
        * bessely(2,k(f)*r1)
        - besselj(2,k(f)*r1)
        * bessely(0,k(f)*r2)
    )

B = lambda f: enum(f) / denom(f) * B0

B_abs = lambda f: abs(B(f))
B_arg = lambda f: arg(B(f))


# ---------------------------------------------------------#
# Generate points for omega axis                           #
# ---------------------------------------------------------#
npts = 1e3
fmin=1
fmax = 2500
n = np.linspace(1,npts,npts)
expufunc = np.frompyfunc(exp,1,1)
frequency_vector = 1*expufunc(n*log(fmax-1)/npts)


# ---------------------------------------------------------#
# Numerically evaluate functions                           #
# ---------------------------------------------------------#
Babsufunc = np.frompyfunc(B_abs,1,1)
B_abs_num = Babsufunc(frequency_vector)
Bargufunc = np.frompyfunc(B_arg,1,1)
B_arg_num = Bargufunc(frequency_vector)


# ---------------------------------------------------------#
# Unfortunately, the  arg() function only  delivers values #
# between -pi and  +pi for the angle of  a complex number, #
# which,  while  correct,  is   not  suitable  for  pretty #
# plotting, so we  will shift the values  larger then zero #
# accordingly for a continuous curve.                      #
# ---------------------------------------------------------#
B_arg_num = np.unwrap(B_arg_num)
B_arg_num = 180/np.pi*B_arg_num


# ---------------------------------------------------------#
# Measurement Values from the experiment                   #
# ---------------------------------------------------------#
frequencies_measured = np.array([    1,     10,      20,      40,      80,     120,     160,   200,    400,    600,    800,   1000, 1200, 1500])
phases_degrees       = np.array([    2,   19.2,    35.2,    56.7,    76.7,      87,      94,   100,    121,    140,    155,    170,  180,  200])
voltages             = np.array([ 7e-2, 6.6e-2, 5.78e-2, 4.18e-2, 2.44e-2, 1.69e-2, 1.27e-2,  1e-2, 4.8e-3, 2.9e-3, 1.9e-3, 1.4e-3, 1e-3, 7e-4])


# ---------------------------------------------------------#
# Plot the Things                                          #
# ---------------------------------------------------------#
font = {
        'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 16,
        }

subplot(2,1,1)
plot(frequency_vector,B_abs_num,color='blue',label='Fitfunktion')
scatter(frequencies_measured,voltages,color='black',s=64,label='Messwerte')
xlabel('Frequenz (Hz)',fontdict=font)
ylabel('Spannung (Volt)',fontdict=font)
title('Betrag des Magnetfelds in Zylinderspule mit Hohlzylinder aus Kupfer, (Messpunkt: auf Zylinderachse, horizontal zentriert)',fontdict=font)
legend(fontsize=16)
xscale('log')
subplot(2,1,2)
plot(frequency_vector,B_arg_num,color='blue',label='Fitfunktion')
scatter(frequencies_measured,-phases_degrees,color='black',s=64,label='Messwerte')
xlabel('Frequenz (Hz)',fontdict=font)
ylabel('Phase (Grad)',fontdict=font)
title('Phase des Magnetfelds in Zylinderspule mit Hohlzylinder aus Kupfer (Messpunkt: auf Zylinderachse, horizontal zentriert)',fontdict=font)
legend(fontsize=16)
xscale('log')
show()
