#!/usr/bin/env python3

from sympy import *
from mpmath import *
from matplotlib.pyplot import *
import matplotlib.ticker as plticker
#init_printing()     # make things prettier when we print stuff for debugging.


# ************************************************************************** #
# B-Field,   Cylinder  Coil   with  Massive   Alu  Cylinder                  #
#                                                                            #
# High-frequency approximation, applied to both  high frequency (450 Hz) and #
# low  frequency (30  Hz),  to show  that  it is  indeed  suitable for  high #
# frequencies and unsuitable for low frequencies.                            #
# ************************************************************************** #

# All values are in standard SI units unless otherwise noted.


# ---------------------------------------------------------#
# Define Variables and Constants                           #
# ---------------------------------------------------------#
var('mu0 B_abs B_arg B B0 j0 k r r0 w f sigma denom enum')
mu0   = 4*pi*1e-7
#sigma = 37.7e6                 # conductivity of aluminium (de.wikipedia.org)
sigma = 18e6                                                   # affects phase
B0    = 5.5e-2                 # does not affect phase, use for scaling abs(B)
r0    = 45e-3
freq = 450                                     # frequency was fixed at 450 Hz
npts = 1e2
rmin_high=25e-3
#rmax_high=50e-3
rmax_high=45e-3
rmin_low=0
rmax_low=45e-3
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
plot_11_label_y          = 'gemessene Spannung (mV)'
plot_21_label_y          = 'gemessene Spannung (mV)'
plot_12_label_y          = 'Phase (Grad)'
plot_11_title            = r"Hochfrequenzn\"aherung: Betrag des Magnetfelds in Zylinderspule mit Vollzylinder aus Aluminium, (Frequenz: 450 Hz, horizontal zentriert)"
plot_12_title            = r"Hochfrequenzn\"aherung: Phase des Magnetfelds in Zylinderspule mit Vollzylinder aus Aluminium (Frequenz: 450 Hz, horizontal zentriert)"
plot_21_title            = r"Hochfrequenzn\"aherung: Betrag des Magnetfelds in Zylinderspule mit Vollzylinder aus Aluminium, (Frequenz: 30 Hz, horizontal zentriert)"
plot_22_title            = r"Hochfrequenzn\"aherung: Phase des Magnetfelds in Zylinderspule mit Vollzylinder aus Aluminium (Frequenz: 30 Hz, horizontal zentriert)"

# Set ticker intervals for plots (in millimeters)
loc1 = plticker.MultipleLocator(base=2.5)
loc2 = plticker.MultipleLocator(base=5)


# ---------------------------------------------------------#
# Function for magnetic Field B                            #
# ---------------------------------------------------------#
# See formula 11 on p.8 of script for experiment.

s_skin = sqrt(2/(2*pi*freq*mu0*sigma))
x = lambda r: r0-r
B = lambda r: B0 * exp(-x(r)/s_skin) * exp(mpc(0,-x(r)/s_skin))

B_abs = lambda r: abs(B(r))
B_arg = lambda r: arg(B(r))


# ---------------------------------------------------------#
# Generate points for radius axis                          #
# ---------------------------------------------------------#
radii_high = np.linspace(rmin_high,rmax_high,npts)
radii_low = np.linspace(rmin_low,rmax_low,npts)


# ---------------------------------------------------------#
# Numerically evaluate function                            #
# ---------------------------------------------------------#
Babsufunc = np.frompyfunc(B_abs,1,1)
Bargufunc = np.frompyfunc(B_arg,1,1)
B_abs_num_high = Babsufunc(radii_high)
B_arg_num_high = Bargufunc(radii_high)
B_abs_num_low = Babsufunc(radii_low)
B_arg_num_low = Bargufunc(radii_low)


# ---------------------------------------------------------#
# Unfortunately, the  arg() function only  delivers values #
# between -pi and  +pi for the angle of  a complex number, #
# which,  while  correct,  is   not  suitable  for  pretty #
# plotting, so we  will shift the values  larger then zero #
# accordingly for a continuous curve.                      #
# ---------------------------------------------------------#
B_arg_num_high = np.unwrap(B_arg_num_high)
B_arg_num_low = np.unwrap(B_arg_num_low)


# ---------------------------------------------------------#
# Measurement Values from the actual experiment            #
# ---------------------------------------------------------#
radii_measured_high  = np.array([    25,  27.5,    30,  32.5,    35,   37.5,    40,  42.5,    45,   47.5,    50])
voltages_high        = np.array([1.5e-3,2.2e-3,3.6e-3,5.9e-3,9.5e-3,1.55e-2,2.5e-2,3.9e-2,5.5e-2,5.75e-2,3.8e-2])
phases_degrees_high  = np.array([   215,   183,   152,   125,   100,     73,    47,    24,   5.2,    0.2,     0])

radii_measured_low = np.array([      0,      5,     10,    15,  20,    25,    30,    35,    40,    45,    50])
voltages_low       = np.array([2.86e-2,2.85e-2,2.87e-2,2.9e-2,3e-2,3.3e-2,3.8e-2,4.5e-2,5.4e-2,6.2e-2,3.7e-2])
phases_degrees_low = np.array([    111,    109,    104,    94,  81,    65,  48.5,    32,    16,   2.7,     0])

# ---------------------------------------------------------#
# Scale values for improved legibility in plot             #
# ---------------------------------------------------------#
# We scale from meters to millimeters,  from rad to degress.
B_abs_num_high = 1e3 * B_abs_num_high
B_abs_num_low = 1e3 * B_abs_num_low
radii_high     = 1e3 * radii_high
radii_low     = 1e3 * radii_low
voltages_high  = 1e3 * voltages_high
voltages_low  = 1e3 * voltages_low
B_arg_num_high = 180/pi*B_arg_num_high - 360
B_arg_num_low = 180/pi*B_arg_num_low - 360
rmin_high      = 1e3 * rmin_high
rmax_high      = 1e3 * rmax_high
rmin_low      = 1e3 * rmin_low
rmax_low      = 1e3 * rmax_low


# ---------------------------------------------------------#
# Plot the Things                                          #
# ---------------------------------------------------------#

matplotlib.pyplot.rc('text', usetex=True)
matplotlib.pyplot.rc('font', family='serif')

fig1   = figure(1)
axes11 = fig1.add_subplot(211)
axes11.plot(radii_high,B_abs_num_high,color=plot_color_fit,label=plot_label_fit)
axes11.scatter(radii_measured_high,
        voltages_high,
        color=plot_color_measurements,
        s=plot_size_measurements,
        label=plot_label_measurements
        )
axes11.set_xlim([rmin_high*0.9,rmax_high*1.1])
axes11.set_xscale(plot_scale_x)
axes11.set_xlabel(plot_label_x,fontdict=font)
axes11.set_ylabel(plot_11_label_y,fontdict=font)
axes11.set_title(plot_11_title,fontdict=font)
axes11.legend(fontsize=plot_legend_fontsize,loc='upper left')
axes11.xaxis.set_major_locator(loc1)

axes12 = fig1.add_subplot(212)
axes12.plot(radii_high,B_arg_num_high,color=plot_color_fit,label=plot_label_fit)
axes12.scatter(radii_measured_high,
        -phases_degrees_high,
        color=plot_color_measurements,
        s=plot_size_measurements,
        label=plot_label_measurements
        )
axes12.set_xlim([rmin_high*0.9,rmax_high*1.1])
axes12.set_xscale(plot_scale_x)
axes12.set_xlabel(plot_label_x,fontdict=font)
axes12.set_ylabel(plot_12_label_y,fontdict=font)
axes12.set_title(plot_12_title,fontdict=font)
axes12.legend(fontsize=plot_legend_fontsize,loc='upper left')
axes12.xaxis.set_major_locator(loc1)

fig2   = figure(2)
axes21 = fig2.add_subplot(211)
axes21.plot(radii_low,B_abs_num_low,color=plot_color_fit,label=plot_label_fit)
axes21.scatter(radii_measured_low,
        voltages_low,
        color=plot_color_measurements,
        s=plot_size_measurements,
        label=plot_label_measurements
        )
axes21.set_xlim([rmin_low*0.9,rmax_low*1.1])
axes21.set_xscale(plot_scale_x)
axes21.set_xlabel(plot_label_x,fontdict=font)
axes21.set_ylabel(plot_21_label_y,fontdict=font)
axes21.set_title(plot_21_title,fontdict=font)
axes21.legend(fontsize=plot_legend_fontsize,loc='upper left')
axes21.xaxis.set_major_locator(loc2)

axes22 = fig2.add_subplot(212)
axes22.plot(radii_low,B_arg_num_low,color=plot_color_fit,label=plot_label_fit)
axes22.scatter(radii_measured_low,
        -phases_degrees_low,
        color=plot_color_measurements,
        s=plot_size_measurements,
        label=plot_label_measurements
        )
axes22.set_xlim([rmin_low*0.9,rmax_low*1.1])
axes22.set_xscale(plot_scale_x)
axes22.set_xlabel(plot_label_x,fontdict=font)
axes22.set_ylabel(plot_12_label_y,fontdict=font)
axes22.set_title(plot_12_title,fontdict=font)
axes22.legend(fontsize=plot_legend_fontsize,loc='upper left')
axes22.xaxis.set_major_locator(loc2)

show()
