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
# -------------------------------------------------------- #
#mp.dps=25  # decimal places
mp.prec=512 # precision in bits


# ---------------------------------------------------------#
# Init, Define Variables and Constants                     #
# ---------------------------------------------------------#
mu0   = 4*pi*1e-7                                        # vacuum permeability
sigma = 1.25e6                               # fit parameter: adjust as needed
r     = 0             # radial position of measurement probe. Centered on axis
dsp   = 98e-3                                               # diameter of coil
rsp   = dsp / 2                                               # radius of coil
r1    = 30e-3                                # inner radius of copper cylinder
r2    = 35e-3                                # outer radius of copper cylinder
B0    = 6.9e-2                          # fit parameter: adjust this as needed
N0    = 574                                   # number of turns of copper coil
l     = 500e-3                                         # length of copper coil
npts  = 1e3
fmin  = 8e1
fmax  = 5e4
font = {
        'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 16,
        }
plot_color_fit          = 'blue'
plot_linewidth          = 1
plot_scale_x            = 'log'
plot_label_x            = 'Frequenz (Hz)'
plot_label_y            = 'Selbstinduktion L (mH)'
plot_title              = """
Selbstinduktionskoeffizient, Kupferspule mit Hohlzylinder \
aus rostfreiem Stahl
"""


# ---------------------------------------------------------#
# Functions                                                #
#                                                          #
# See formula 28 on p.15 of script for experiment.         #
#                                                          #
# NOTE: We use  frequency f  instead of  angular frequency #
# omega since that is what we actually set on the function #
# generator.                                               #
# ---------------------------------------------------------#

var('f')

k = lambda f: sqrt((2*np.pi*f*mu0*sigma)/2)*(mpc(1,-1))

enum1 = lambda f:(
          besselj(0,k(f)*r1)
        * bessely(2,k(f)*r1)
        - besselj(2,k(f)*r1)
        * bessely(0,k(f) * r1)
    )
denom1 = lambda f: (
          besselj(0,k(f)*r2)
        * bessely(2,k(f)*r1)
        - besselj(2,k(f)*r1)
        * bessely(0,k(f) * r2)
    )
enum2 = lambda f:(
          r2 * (
            besselj(1,k(f)*r2)
            * bessely(2,k(f)*r1)
            - besselj(2,k(f)*r1)
            * bessely(1,k(f) * r2)
        )
        - r1 * (
            besselj(1,k(f)*r1)
            * bessely(2,k(f)*r1)
            - besselj(2,k(f)*r1)
            * bessely(1,k(f) * r1)
        )
    )
denom2 = lambda f: (
          besselj(0,k(f)*r2)
        * bessely(2,k(f)*r1)
        - besselj(2,k(f)*r1)
        * bessely(0,k(f) * r2)
    )
term3 = rsp ** 2 - r2**2
prefactor = mu0 * pi * N0**2 / l

phi_norm = lambda f:(
        prefactor * (
            r1**2    * enum1(f)/denom1(f)
            + 2/k(f) * enum2(f)/denom2(f)
            + term3
        )
    )

phi_norm_abs = lambda w: abs(phi_norm(w))
phi_norm_arg = lambda w: arg(phi_norm(w))

L = lambda f: re(phi_norm(f))


# ---------------------------------------------------------#
# Generate points for omega axis                           #
# ---------------------------------------------------------#
n                = np.linspace(0,npts,npts)
expufunc         = np.frompyfunc(exp,1,1)
logufunc         = np.frompyfunc(log,1,1)

# Old version:
#frequency_vector = fmin*expufunc(n*log(fmax-fmin)/npts)

# New version:
frequency_vector = expufunc((1-n/npts)*log(fmin)) * expufunc(n*log(fmax)/npts)

# Debugging section:
#scatter(frequency_vector,n)
#xscale('log')
#frequency_vector1 = expufunc((1-n/npts)*log(fmin))
#frequency_vector2 = expufunc(n*log(fmax)/npts)
#print(frequency_vector)
#scatter(frequency_vector1,n,color='red')
#scatter(frequency_vector2,n,color='magenta')
#show()
#exit()


# ---------------------------------------------------------#
# Numerically evaluate functions                           #
# ---------------------------------------------------------#
L_ufunc = np.frompyfunc(L,1,1)
L_num   = L_ufunc(frequency_vector)
L_num   = 1e3 * L_num           # improve legibility of plot


# ---------------------------------------------------------#
# Plot the Things                                          #
# ---------------------------------------------------------#
fig  = figure(1)
axes = fig.add_subplot(111)
axes.plot(frequency_vector,L_num,linewidth=plot_linewidth,color=plot_color_fit)
axes.set_xscale(plot_scale_x)
axes.set_xlim([fmin*0.9,fmax*1.1])
axes.set_xlabel(plot_label_x,fontdict=font)
axes.set_ylabel(plot_label_y,fontdict=font)
axes.set_title(plot_title,fontdict=font)

show()
