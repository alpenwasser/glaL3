#!/usr/bin/env python3

from sympy import *
from mpmath import *
from matplotlib.pyplot import *
#init_printing()     # make things prettier when we print stuff for debugging.


# ************************************************************************** #
# Magnetic Flow normed  by current, copper coil with  hollow copper cylinder #
# inserted.                                                                  #
# ************************************************************************** #

# All values are in standard SI units unless otherwise noted.


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
        'size'   : 11,
        }
plot_legend_fontsize    = 11
plot_color_fit          = 'blue'
plot_color_measurements = 'black'
plot_label_measurements = 'Messwerte'
plot_size_measurements  = 32
plot_scale_x            = 'log'
plot_label_fit          = 'Fitfunktion'
plot_label_x            = 'Frequenz (Hz)'
plot_1_label_y          = r"$\displaystyle \biggl| \frac{\Phi}{I} \biggr|$ $\biggl( \displaystyle \frac{Vs}{A} \biggr)$"
plot_2_label_y          = r"$\displaystyle arg\biggl( \frac{\Phi}{I} \biggr)$ (Grad)"
plot_1_title            = r"Betrag Magn. Fluss normiert auf Spulenstrom, Spule mit Stahlrohr"
plot_2_title            = r"Phase Magn. Fluss normiert auf Spulenstrom, Spule mit Stahlrohr"


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

phi_norm_abs = lambda f: abs(phi_norm(f))
phi_norm_arg = lambda f: arg(phi_norm(f))


# ---------------------------------------------------------#
# Generate points for omega axis                           #
# ---------------------------------------------------------#
# See also separate file stuetzpunkte.py
n                = np.linspace(0,npts,npts)
expufunc         = np.frompyfunc(exp,1,1)
frequency_vector = expufunc((1-n/npts)*log(fmin)) * expufunc(n*log(fmax)/npts)


# ---------------------------------------------------------#
# Numerically evaluate functions                           #
# ---------------------------------------------------------#
phi_norm_abs_ufunc = np.frompyfunc(phi_norm_abs,1,1)
phi_norm_abs_num   = phi_norm_abs_ufunc(frequency_vector)
phi_norm_arg_ufunc = np.frompyfunc(phi_norm_arg,1,1)
phi_norm_arg_num   = phi_norm_arg_ufunc(frequency_vector)


# ---------------------------------------------------------#
# Unfortunately, the  arg() function only  delivers values #
# between -pi and  +pi for the angle of  a complex number, #
# which,  while  correct,  is   not  suitable  for  pretty #
# plotting, so we  will shift the values  larger then zero #
# accordingly for a continuous curve.                      #
# ---------------------------------------------------------#
phi_norm_arg_num   = np.unwrap(phi_norm_arg_num)


# ---------------------------------------------------------#
# Scale values for improved legibility in plot             #
# ---------------------------------------------------------#
phi_norm_abs_num   = phi_norm_abs_ufunc(frequency_vector)
phi_norm_abs_num   = 1e3 * phi_norm_abs_num
phi_norm_arg_num   = 180 / pi * phi_norm_arg_num


# ---------------------------------------------------------#
# Plot the Things                                          #
# ---------------------------------------------------------#
matplotlib.pyplot.rc('text', usetex=True)
matplotlib.pyplot.rc('font', family='serif')

fig   = figure(1)
axes1 = fig.add_subplot(211)
axes1.plot(frequency_vector,phi_norm_abs_num,color=plot_color_fit,label=plot_label_fit)
axes1.set_xlim([fmin*0.9,fmax*1.1])
axes1.set_xscale(plot_scale_x)
axes1.set_xlabel(plot_label_x,fontdict=font)
axes1.set_ylabel(plot_1_label_y,fontdict=font)
axes1.set_title(plot_1_title,fontdict=font)

axes2 = fig.add_subplot(212)
axes2.plot(frequency_vector,phi_norm_arg_num,color=plot_color_fit,label=plot_label_fit)
axes2.set_xlim([fmin*0.9,fmax*1.1])
axes2.set_xscale(plot_scale_x)
axes2.set_xlabel(plot_label_x,fontdict=font)
axes2.set_ylabel(plot_2_label_y,fontdict=font)
axes2.set_title(plot_2_title,fontdict=font)

fig.subplots_adjust(bottom=0.1,left=0.125,right=0.925,top=0.95,hspace=0.5)

fig.savefig('plots-pgf/hollow--st--freq--phi-norm.pgf')
fig.savefig('plots-pdf/hollow--st--freq--phi-norm.pdf')
