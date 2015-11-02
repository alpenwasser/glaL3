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
# Init, Define Variables and Constants                     #
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
fmin=1
fmax=250
font = {
        'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 11,
        }
plot_color_fit          = 'blue'
plot_color_measurements = 'black'
plot_linewidth          = 1
plot_scale_x            = 'log'
plot_label_x            = 'Frequenz (Hz)'
plot_label_y            = 'Selbstinduktion L (mH)'
plot_title              = "Selbstinduktionskoeffizient, Spule mit Vollzylinder"


# ---------------------------------------------------------#
# Functions                                                #
#                                                          #
# See formula 22 on p.12 of script for experiment.         #
#                                                          #
# NOTE: We use  frequency f  instead of  angular frequency #
# omega since that is what we actually set on the function #
# generator.                                               #
# ---------------------------------------------------------#
var('f')

k = lambda f: sqrt((2*np.pi*f*mu0*sigma)/2)*(mpc(1,-1))

LRand = (mu0 * 2 * pi * r0 * (rsp - r0) * N0**2) / l

L = lambda f:(
          (mu0*2*pi*r0*N0**2) / l
        * re(besselj(1,k(f)*r0) / (k(f)
        * besselj(0,k(f)*r0)))
        + LRand
    )


# ---------------------------------------------------------#
# Generate points for frequency axis                       #
# ---------------------------------------------------------#
n                = np.linspace(1,npts,npts)
expufunc         = np.frompyfunc(exp,1,1)
frequency_vector = fmin*expufunc(n*log(fmax-fmin)/npts)


# ---------------------------------------------------------#
# Numerically evaluate function                            #
# ---------------------------------------------------------#
L_ufunc = np.frompyfunc(L,1,1)
L_num   = L_ufunc(frequency_vector)
L_num = 1e3 * L_num                     # improve legibility


# ---------------------------------------------------------#
# Plot the Things                                          #
# ---------------------------------------------------------#
matplotlib.pyplot.rc('text', usetex=True)
matplotlib.pyplot.rc('font', family='serif')

figwidth = 8.27 # in inches
fig  = figure(1,figsize=(figwidth,figwidth*0.4))
axes = fig.add_subplot(111)
axes.plot(frequency_vector,L_num,linewidth=plot_linewidth,color=plot_color_fit)
axes.set_xscale(plot_scale_x)
axes.set_xlim([fmin*0.9,fmax*1.1])
axes.set_xlabel(plot_label_x,fontdict=font)
axes.set_ylabel(plot_label_y,fontdict=font)
axes.set_title(plot_title,fontdict=font)

fig.subplots_adjust(bottom=0.15,left=0.1,right=0.9,top=0.9,hspace=0.5)

fig.savefig('plots-pgf/massive--alu--L.pgf')
fig.savefig('plots-pdf/massive--alu--L.pdf')
