#!/usr/bin/env python3

from sympy import *
from mpmath import *
from matplotlib.pyplot import *
#init_printing()     # make things prettier when we print stuff for debugging.


# ************************************************************************** #
# Magnetic field inside copper coil with  hollow stainless steel cylinder    #
# Low-frequency approximation                                                #
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
#mp.prec=80 # precision in bits

# ---------------------------------------------------------#
# Init, Define Variables and Constants                     #
# ---------------------------------------------------------#
mu0   = 4*pi*1e-7                                        # vacuum permeability
#sigma = 52e6                            # de.wikipedia.org/wiki/Kupfer: 58.1e6
sigma = 1.25e6                            # de.wikipedia.org/wiki/Kupfer: 58.1e6
dsp   = 98e-3                                               # diameter of coil
rsp   = dsp / 2                                               # radius of coil
r1    = 30e-3                                # inner radius of copper cylinder
r2    = 35e-3                                # outer radius of copper cylinder
r_avg = (r1+r2)/2                                 # average radius of cylinder
d_rohr = r2 - r1                           # wall thickness of copper cylinder
N0    = 574                                   # number of turns of copper coil
l     = 500e-3                                         # length of copper coil
npts  = 1e2
fmin  = 1
fmax  = 7500
font = {
        'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 16,
        }
plot_legend_fontsize    = 16
plot_color_fit          = 'blue'
plot_color_ratio        = 'magenta'
plot_color_measurements = 'black'
plot_label_measurements = 'Messwerte'
plot_size_measurements  = 64
plot_scale_x            = 'log'
plot_label_fit          = r"Fitfunktion (N\"aherung)"
plot_label_ratio        = r"$\displaystyle \frac{d_{Rohr}}{s_{skin}}$: Sollte $<1$ sein f\"ur G\"ultigkeit der Approximation."
plot_label_x            = 'Frequenz (Hz)'
plot_1_label_y          = 'gemessene Spannung (mV)'
plot_2_label_y          = 'Phase (Grad)'
plot_1_title            = r"N\"aherungsl\"osung: Betrag des Magnetfelds in Zylinderspule mit Hohlzylinder aus rostfreiem Stahl, (Messpunkt: auf Zylinderachse, horizontal zentriert)"
plot_2_title            = r"N\"aherungsl\"osung: Phase des Magnetfelds in Zylinderspule mit Hohlzylinder aus rostfreiem Stahl (Messpunkt: auf Zylinderachse, horizontal zentriert)"

    # ---------------------------------------------------- #
    # current in copper coil. This is a scaling parameter, #
    # not the  measured value. measured value  was: 200 mA #
    # This is due to the  fact that the measurement values #
    # are  voltages  representing  the  B-Field,  not  the #
    # actual B-Field itself.                               #
    # ---------------------------------------------------- #
I0    = 48.5

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

enum1  = mu0*N0*I0
denom1 = l
enum2  = 2
denom2 = lambda f: mpc(2,2*pi*f*mu0*r_avg*d_rohr*sigma)

B = lambda f: enum1 / denom1 * enum2 / denom2(f)

B_abs = lambda f: abs(B(f))
B_arg = lambda f: arg(B(f))

s_skin = lambda f: sqrt(2/(2*pi*f*mu0*sigma))

# ---------------------------------------------------------#
# Generate points for frequency axis                       #
# ---------------------------------------------------------#
n = np.linspace(1,npts,npts)
expufunc = np.frompyfunc(exp,1,1)
frequency_vector = fmin*expufunc(n*log(fmax-fmin)/npts)


# ---------------------------------------------------------#
# Numerically evaluate functions                           #
# ---------------------------------------------------------#
Babsufunc        = np.frompyfunc(B_abs,1,1)
B_abs_num        = Babsufunc(frequency_vector)
Bargufunc        = np.frompyfunc(B_arg,1,1)
B_arg_num        = Bargufunc(frequency_vector)
s_skin_ufunc     = np.frompyfunc(s_skin,1,1)
s_skin_num       = s_skin_ufunc(frequency_vector)
s_skin_ratio_num = d_rohr / s_skin_num # should be < 1 for validity
#print(s_skin_ratio_num)
#print(B_abs_num)
#exit()


# ---------------------------------------------------------#
# Unfortunately, the  arg() function only  delivers values #
# between -pi and  +pi for the angle of  a complex number, #
# which,  while  correct,  is   not  suitable  for  pretty #
# plotting, so we  will shift the values  larger then zero #
# accordingly for a continuous curve.                      #
# ---------------------------------------------------------#
B_arg_num = np.unwrap(B_arg_num)


# ---------------------------------------------------------#
# Measurement Values from experiment                       #
# ---------------------------------------------------------#
frequencies_measured = np.array([       1,      10,      20,      40,      80,     120,     160,     200,     400,     600,    800,    1000,    1200,   1500,   1750,    2000,   2500,   3500,   5000,   7500])
phases_degrees       = np.array([       0,    0.45,    0.95,     1.8,     3.6,     5.4,     7.2,       9,    17.5,    25.4,   32.4,    38.4,    43.5,     50,     54,      58,     64,     71,     78,     88])
voltages             = np.array([ 6.96e-2, 6.97e-2, 6.97e-2, 6.97e-2, 6.92e-2, 6.91e-2, 6.87e-2, 6.62e-2, 6.27e-2, 6.27e-2, 5.9e-2, 5.45e-2, 5.05e-2, 4.5e-2, 4.1e-2, 3.72e-2, 3.2e-2, 2.4e-2, 1.8e-2, 1.2e-2])


# ---------------------------------------------------------#
# Scale values for improved legibility in plot             #
# ---------------------------------------------------------#
B_abs_num = 1e3 * B_abs_num
voltages  = 1e3 * voltages
B_arg_num = 180/np.pi*B_arg_num


# ---------------------------------------------------------#
# Plot the Things                                          #
# ---------------------------------------------------------#
matplotlib.pyplot.rc('text', usetex=True)
matplotlib.pyplot.rc('font', family='serif')

fig   = figure(1)
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
axes2.legend(fontsize=plot_legend_fontsize,loc='center left')

axes3 = axes2.twinx()
axes3.plot(frequency_vector,s_skin_ratio_num,color=plot_color_ratio,label=plot_label_ratio)
axes3.legend(fontsize=plot_legend_fontsize,loc='upper right')
axes3.set_xlim([fmin*0.9,fmax*1.1])

show()
