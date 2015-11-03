#!/usr/bin/env python3

from sympy import *
from mpmath import *
from matplotlib.pyplot import *
#init_printing()     # make things prettier when we print stuff for debugging.


# ************************************************************************** #
# Magnetic field inside copper coil with massive alu cylinder                #
# ************************************************************************** #

# All values are in standard SI units unless otherwise noted.


# ---------------------------------------------------------#
# Define Variables and Constants                           #
# ---------------------------------------------------------#
mu0   = 4*pi*1e-7
#sigma = 37.7e6                 # conductivity of aluminium (de.wikipedia.org)
sigma = 23.75e6
r     = 0
r0    = 45e-3
B0    = 6.9e-2                             # adjust this as needed for scaling
npts  = 1e3
fmin  = 1
fmax  = 250
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
plot_size_measurements  = 16
plot_scale_x            = 'log'
plot_label_fit          = 'Fitfunktion'
plot_label_x            = 'Frequenz (Hz)'
plot_1_label_y          = 'gemessene Spannung (mV)'
plot_2_label_y          = 'Phase (Grad)'
plot_1_title            = r"Betrag Magnetfeld, Spule mit Vollzylinder"
plot_2_title            = r"Phase Magnetfeld, Spule mit Vollzylinder"


# ---------------------------------------------------------#
# Function for magnetic Field B                            #
# ---------------------------------------------------------#
# See formula 21 on p.11 of script for experiment.

var('f')

k = lambda f: sqrt((2*np.pi*f*mu0*sigma)/2)*(mpc(1,-1))

# Enumerator:
enum  = lambda f: besselj(0,k(f)*r)
denom = lambda f: besselj(0,k(f)*r0)

B = lambda f: enum(f) / denom(f) * B0

B_abs = lambda f: abs(B(f))
B_arg = lambda f: arg(B(f))

# ---------------------------------------------------------#
# Generate points for frequency axis                       #
# ---------------------------------------------------------#
n                = np.linspace(0,npts,npts)
expufunc         = np.frompyfunc(exp,1,1)
frequency_vector = fmin*expufunc(n*log(fmax-fmin)/npts)


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


# ---------------------------------------------------------#
# Measurement Values from the actual experiment            #
# ---------------------------------------------------------#
frequencies_measured = np.array([     1,     5,    10,    15,  20,     30,    40,     60,  80,   100,   120,   160, 200, 250])
phases_degrees       = np.array([   5.4,    26,    50,    69,  85,    111,   132,    166, 196,   220,   243,   283, 320, 350])
voltages             = np.array([6.9e-2,6.5e-2,5.7e-2,4.8e-2,4e-2,2.85e-2,2.1e-2,1.25e-2,8e-3,5.4e-3,3.6e-3,1.9e-3,1e-3,6e-4])


# ---------------------------------------------------------#
# Scale values for improved legibility in plot             #
# ---------------------------------------------------------#
B_abs_num = 1e3 * B_abs_num
voltages  = 1e3 * voltages
B_arg_num = 180/pi*B_arg_num

# ---------------------------------------------------------#
# Plot the Things                                          #
# ---------------------------------------------------------#
matplotlib.pyplot.rc('text', usetex=True)
matplotlib.pyplot.rc('font', family='serif')

fig   = figure(1)
#axes1 = fig.add_subplot(211,axisbg='black')
axes1 = fig.add_subplot(211)
axes1.plot(frequency_vector,B_abs_num,color=plot_color_fit,label=plot_label_fit)
axes1.scatter(frequencies_measured,
        voltages,
        color=plot_color_measurements,
        s=plot_size_measurements,
        label=plot_label_measurements
        )
axes1.set_xlim([fmin*0.9,fmax*1.1])
axes1.set_xscale(plot_scale_x)
axes1.set_xlabel(plot_label_x,fontdict=font)
axes1.set_ylabel(plot_1_label_y,fontdict=font)
axes1.set_title(plot_1_title,fontdict=font)
axes1.legend(fontsize=plot_legend_fontsize)

axes2 = fig.add_subplot(212)
axes2.plot(frequency_vector,B_arg_num,color=plot_color_fit,label=plot_label_fit)
axes2.scatter(frequencies_measured,
        -phases_degrees,
        color=plot_color_measurements,
        s=plot_size_measurements,
        label=plot_label_measurements
        )
axes2.set_xlim([fmin*0.9,fmax*1.1])
axes2.set_xscale(plot_scale_x)
axes2.set_xlabel(plot_label_x,fontdict=font)
axes2.set_ylabel(plot_2_label_y,fontdict=font)
axes2.set_title(plot_2_title,fontdict=font)
axes2.legend(fontsize=plot_legend_fontsize)

fig.subplots_adjust(bottom=0.1,left=0.1,right=0.9,top=0.95,hspace=0.5)

fig.savefig('plots-pgf/massive--alu--freq.pgf')
fig.savefig('plots-pdf/massive--alu--freq.pdf')
