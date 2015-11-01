#!/usr/bin/env python3

from sympy import *
from mpmath import *
from matplotlib.pyplot import *
#init_printing()     # make things prettier when we print stuff for debugging.


# ************************************************************************** #
# Self-Inductance L of copper coil with massive aluminium cylinder inserted  #
# ************************************************************************** #

# All values are in standard SI units unless otherwise noted.

# ---------------------------------------------------------#
# Define Variables and Constants                           #
# ---------------------------------------------------------#
mu0   = 4*pi*1e-7                                        # vacuum permeability
#sigma = 37.7e6                 # conductivity of aluminium (de.wikipedia.org)
sigma = 23.75e6                         # de.wikipedia.org/wiki/Kupfer: 58.1e6
r     = 0             # radial position of measurement probe. Centered on axis
dsp   = 98e-3                                               # diameter of coil
rsp   = dsp / 2                                               # radius of coil
r0    = 45e-3                                         # radius of alu cylinder
B0    = 6.9e-2                             # adjust this as needed for scaling
N0    = 574                                   # number of turns of copper coil
l     = 500e-3                                         # length of copper coil
npts  = 1e3
fmin  = 1
fmax  = 250
font = {
        'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 11,
        }
plot_color_fit = 'blue'
plot_linewidth = 1
plot_label_fit = 'Fitfunktion'
plot_scale_x   = 'log'
plot_label_x   = 'Frequenz (Hz)'
plot_1_label_y = r"$\displaystyle \biggl| \frac{\Phi}{I} \biggr|$ $\biggl( \displaystyle \frac{Vs}{A} \biggr)$"
plot_2_label_y = r"$\displaystyle arg\biggl( \frac{\Phi}{I} \biggr)$ (Grad)"
plot_1_title   = r"Betrag Magn. Fluss normiert auf Spulenstrom, Spule mit Kupferrohr"
plot_2_title   = r"Phase Magn. Fluss normiert auf Spulenstrom, Spule mit Kupferrohr"


# ---------------------------------------------------------#
# Functions                                                #
#                                                          #
# See formula 23 on p.12 of script for experiment.         #
#                                                          #
# NOTE: We use  frequency f  instead of  angular frequency #
# omega since that is what we actually set on the function #
# generator.                                               #
# ---------------------------------------------------------#
var('f')

k = lambda f: sqrt((2*np.pi*f*mu0*sigma)/2)*(mpc(1,-1))

phi_norm = lambda f:(
          (mu0 * 2 * pi * r0 * N0**2) / l
        * (
              besselj(1,k(f)*r0) / (k(f) * besselj(0,k(f)*r0))
            + rsp - r0
        )
    )
phi_norm_abs = lambda f: abs(phi_norm(f))
phi_norm_arg = lambda f: arg(phi_norm(f))


# ---------------------------------------------------------#
# Generate points for frequency axis                       #
# ---------------------------------------------------------#
n                = np.linspace(1,npts,npts)
expufunc         = np.frompyfunc(exp,1,1)
frequency_vector = fmin*expufunc(n*log(fmax-fmin)/npts)


# ---------------------------------------------------------#
# Numerically evaluate function                            #
# ---------------------------------------------------------#
phi_norm_abs_ufunc = np.frompyfunc(phi_norm_abs,1,1)
phi_norm_abs_num   = phi_norm_abs_ufunc(frequency_vector)
phi_norm_arg_ufunc = np.frompyfunc(phi_norm_arg,1,1)
phi_norm_arg_num   = phi_norm_arg_ufunc(frequency_vector)


# ---------------------------------------------------------#
# Correct Phase                                            #
# ---------------------------------------------------------#
phi_norm_arg_num = np.unwrap(phi_norm_arg_num)


# ---------------------------------------------------------#
# Scale values for improved legibility in plot             #
# ---------------------------------------------------------#
phi_norm_arg_num   = 180 / pi * phi_norm_arg_num
phi_norm_abs_num   = phi_norm_abs_ufunc(frequency_vector)
phi_norm_abs_num   = 1e3 * phi_norm_abs_num


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

fig.savefig('plots-pgf/massive--alu--freq--phi-norm.pgf')
fig.savefig('plots-pdf/massive--alu--freq--phi-norm.pdf')
