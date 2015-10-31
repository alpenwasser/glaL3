#!/usr/bin/env python3

from sympy import *
from sympy.external import import_module
#from numpy import *
from mpmath import *
from matplotlib.pyplot import *
from math import copysign
import matplotlib.ticker as plticker
#init_printing()     # make things prettier when we print stuff for debugging.


# ************************************************************************** #
# B-Field, Cylinder Coild with Massive Alu Cylinder                          #
# ************************************************************************** #

# All values are in standard SI units unless otherwise noted.

# ---------------------------------------------------------#
# Measurement in Center, depending on Frequency            #
# ---------------------------------------------------------#


# ---------------------------------------------------------#
# Init, Define Variables and Constants                     #
# ---------------------------------------------------------#
var('mu0 B_abs B_arg B B0 j0 k r r0 w f sigma denom enum')
mu0   = 4*pi*1e-7
#sigma = 37.7e6                 # conductivity of aluminium (de.wikipedia.org)
sigma = 18e6                                                   # affects phase
B0    = 5.5e-2                 # does not affect phase, use for scaling abs(B)
r0    = 45e-3
freq = 450                                      # frequency was fixed at 30 Hz
npts = 1e2
rmin=25e-3
#rmax=50e-3
rmax=45e-3
font = {
        'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 16,
        }
plot_legend_fontsize    = 16
plot_color_fit          = 'blue'
plot_color_measurements = 'black'
plot_label_measurements = 'Messwerte'
plot_size_measurements  = 64
plot_scale_x            = 'linear'
plot_label_fit          = 'Fitfunktion'
plot_label_x            = 'radiale Position bezogen auf Zylinderachse (mm)'
plot_1_label_y          = 'gemessene Spannung (mV)'
plot_2_label_y          = 'Phase (Grad)'
plot_1_title            = """
Betrag des Magnetfelds in Zylinderspule mit Vollzylinder aus \
Aluminium, (Frequenz: 450 Hz, horizontal zentriert)
"""
plot_2_title            = """
Phase des Magnetfelds in Zylinderspule mit Vollzylinder aus \
Aluminium (Frequenz: 450 Hz, horizontal zentriert)
"""
# Set ticker intervals for plots (in millimeters)
loc = plticker.MultipleLocator(base=2.5)


# ---------------------------------------------------------#
# Function for magnetic Field B                            #
# ---------------------------------------------------------#
# See formula 21 on p.11 of script for experiment.

#k = lambda w: sqrt((w*mu0*sigma)/2)*(mpc(1,-1))                     # rad/sec
var('f')
k = lambda f: sqrt((2*pi*f*mu0*sigma)/2)*(mpc(1,-1))          # degrees/sec

# Enumerator:
enum  = lambda r: besselj(0,k(freq)*r)
denom =           besselj(0,k(freq)*r0)

B = lambda r: enum(r) / denom * B0

B_abs = lambda r: abs(B(r))
B_arg = lambda r: arg(B(r))

# Generate points for omega axis of B, store in Bw
radii = np.linspace(rmin,rmax,npts)

# Generate B-Field values for that frequency vector
Babsufunc = np.frompyfunc(B_abs,1,1)
B_abs_num = Babsufunc(radii)
Bargufunc = np.frompyfunc(B_arg,1,1)
B_arg_num = Bargufunc(radii)


# ---------------------------------------------------------#
# Unfortunately, the  arg() function only  delivers values #
# between -pi and  +pi for the angle of  a complex number, #
# which,  while  correct,  is   not  suitable  for  pretty #
# plotting, so we  will shift the values  larger then zero #
# accordingly for a continuous curve.                      #
# ---------------------------------------------------------#
B_arg_num = np.unwrap(B_arg_num)


# ---------------------------------------------------------#
# Measurement Values from the actual experiment            #
# ---------------------------------------------------------#
radii_measured  = np.array([    25,  27.5,    30,  32.5,    35,   37.5,    40,  42.5,    45,   47.5,    50])
voltages        = np.array([1.5e-3,2.2e-3,3.6e-3,5.9e-3,9.5e-3,1.55e-2,2.5e-2,3.9e-2,5.5e-2,5.75e-2,3.8e-2])
phases_degrees  = np.array([   215,   183,   152,   125,   100,     73,    47,    24,   5.2,    0.2,     0])


# ---------------------------------------------------------#
# Scale values for improved legibility in plot             #
# ---------------------------------------------------------#
# We scale from meters to millimeters,  from rad to degress.
B_abs_num = 1e3 * B_abs_num
radii     = 1e3 * radii
voltages  = 1e3 * voltages
B_arg_num = 180/pi*B_arg_num - 360
rmin      = 1e3 * rmin
rmax      = 1e3 * rmax


# ---------------------------------------------------------#
# Plot the Things                                          #
# ---------------------------------------------------------#

fig   = figure(1)
axes1 = fig.add_subplot(211)
axes1.plot(radii,B_abs_num,color=plot_color_fit,label=plot_label_fit)
axes1.scatter(radii_measured,
        voltages,
        color=plot_color_measurements,
        s=plot_size_measurements,
        label=plot_label_measurements
        )
axes1.set_xlim([rmin*0.9,rmax*1.1])
axes1.set_xscale(plot_scale_x)
axes1.set_xlabel(plot_label_x,fontdict=font)
axes1.set_ylabel(plot_1_label_y,fontdict=font)
axes1.set_title(plot_1_title,fontdict=font)
axes1.legend(fontsize=plot_legend_fontsize,loc='upper left')
axes1.xaxis.set_major_locator(loc)

axes2 = fig.add_subplot(212)
axes2.plot(radii,B_arg_num,color=plot_color_fit,label=plot_label_fit)
axes2.scatter(radii_measured,
        -phases_degrees,
        color=plot_color_measurements,
        s=plot_size_measurements,
        label=plot_label_measurements
        )
axes2.set_xlim([rmin*0.9,rmax*1.1])
axes2.set_xscale(plot_scale_x)
axes2.set_xlabel(plot_label_x,fontdict=font)
axes2.set_ylabel(plot_2_label_y,fontdict=font)
axes2.set_title(plot_2_title,fontdict=font)
axes2.legend(fontsize=plot_legend_fontsize,loc='upper left')
axes2.xaxis.set_major_locator(loc)

show()
