#!/usr/bin/env python3

from sympy import *
from mpmath import *
from matplotlib.pyplot import *
#init_printing()     # make things prettier when we print stuff for debugging.


# ************************************************************************** #
# Self-Inductance L of copper coil with hollow copper cylinder inserted.     #
# ************************************************************************** #

# All values are in standard SI units unless otherwise noted.

# -------------------------------------------------------- #
# Default precision is insufficient, therefore we increase #
# precision.   One  can  increase the  number  of  decimal #
# places or bits, where the number of bits places is ~3.33 #
# times the number of decimal places.                      #
# For a maximum  frequency of 50 KHz,  a minimum precision #
# of  273 bits  is required,  otherwise there  will be  an #
# error regarding division by  zero. We therefore set this #
# to 320 bits for  a comfortable margin, without incurring #
# a gigantic  penalty in  execution speed (although  it is #
# already significant at this precision).                  #
# -------------------------------------------------------- #
#mp.dps=25  # decimal places
#mp.prec=320 # precision in bits


# ---------------------------------------------------------#
# Init, Define Variables and Constants                     #
# ---------------------------------------------------------#
mu0   = 4*pi*1e-7                                        # vacuum permeability
sigma = 1.25e6                          # de.wikipedia.org/wiki/Kupfer: 58.1e6
dsp   = 98e-3                                               # diameter of coil
rsp   = dsp / 2                                               # radius of coil
r1    = 30e-3                                # inner radius of copper cylinder
r2    = 35e-3                                # outer radius of copper cylinder
B0    = 6.9e-2                             # adjust this as needed for scaling
N0    = 574                                   # number of turns of copper coil
l     = 500e-3                                         # length of copper coil
R_0   = 5                                  # resistance of coil, assumed value
npts  = 1e3                                        # number of points for plot
fmin  = 8e1                                                # minimum frequency
fmax  = 5e4                                                # maximum frequency
font = {
        'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 11,
        }
plot_color_fit = 'blue'
plot_linewidth = 1
plot_scale_x   = 'log'
plot_label_x   = 'Frequenz (Hz)'
plot_label_y   = 'Widerstand (Ohm)'
plot_title     = r"Ohm'scher Widerstand, Spule mit Stahlrohr"


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

R = lambda f: -2 * pi * f * im(phi_norm(f)) + R_0


# ---------------------------------------------------------#
# Generate points for frequency axis                       #
# ---------------------------------------------------------#
n                = np.linspace(0,npts,npts)
expufunc         = np.frompyfunc(exp,1,1)
#frequency_vector = 1*expufunc(n*log(fmax-1)/npts)
frequency_vector = expufunc((1-n/npts)*log(fmin)) * expufunc(n*log(fmax)/npts)


# ---------------------------------------------------------#
# Numerically evaluate function                            #
# ---------------------------------------------------------#
R_ufunc = np.frompyfunc(R,1,1)
R_num   = R_ufunc(frequency_vector)


# ---------------------------------------------------------#
# Plot the Things                                          #
# ---------------------------------------------------------#
matplotlib.pyplot.rc('text', usetex=True)
matplotlib.pyplot.rc('font', family='serif')

fig  = figure(1)
axes = fig.add_subplot(111)
axes.plot(frequency_vector,R_num,linewidth=plot_linewidth,color=plot_color_fit)
axes.set_xscale(plot_scale_x)
axes.set_xlim([fmin*0.9,fmax*1.1])
axes.set_xlabel(plot_label_x,fontdict=font)
axes.set_ylabel(plot_label_y,fontdict=font)
axes.set_title(plot_title,fontdict=font)

fig.subplots_adjust(bottom=0.1,left=0.1,right=0.9,top=0.95,hspace=0.5)

fig.savefig('plots-pgf/hollow--st--R.pgf')
fig.savefig('plots-pdf/hollow--st--R.pdf')
