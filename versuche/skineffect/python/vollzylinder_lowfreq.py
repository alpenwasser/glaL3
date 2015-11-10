#!/usr/bin/env python3

from sympy import *
from mpmath import *
from matplotlib.pyplot import *
import matplotlib.ticker as plticker
import numpy as np
#init_printing()     # make things prettier when we print stuff for debugging.


# ************************************************************************** #
# B-Field, Cylinder Coild with Massive Alu Cylinder                          #
# ************************************************************************** #

# All values are in standard SI units unless otherwise noted.


# ---------------------------------------------------------#
# Init, Define Variables and Constants                     #
# ---------------------------------------------------------#
mu0   = 4*pi*1e-7
#sigma = 37.7e6                 # conductivity of aluminium (de.wikipedia.org)
rho_kuchling   = 0.027e-6  # resistivity Kuchling 17th edition, p.649, tab. 45
sigma_kuchling = 1/rho_kuchling
#sigma = 21.5e6
sigma_abs = 22.5e6 # great fit for phase
sigma_arg = 23.5e6 # great fit for phase
r0    = 45e-3
#B0    = 6.9e-2                               # adjust this as needed for scaling
B0    = 6.6e-2                               # adjust this as needed for scaling
freq = 30                                     # frequency was fixed at 450 Hz
npts = 1e3
rmin=0
#rmax=50e-3
rmax=45e-3
    # -----------------------------------------------------#
    # Create  a list  for convenient  printing of  vars to #
    # file, add LaTeX where necessary.                     #
    # -----------------------------------------------------#
params = [
        '        ' + r'\textcolor{red}{$\sigma_{Fit,|\hat{B}|}'      + r'$} & \textcolor{red}{$' +  '\SI{'   + str(sigma_abs)       + r'}{\ampere\per\volt\per\meter}' + r'$}\\' + "\n",
        '        ' + r'\textcolor{red}{$\sigma_{Fit,\angle\hat{B}}'  + r'$} & \textcolor{red}{$' +  '\SI{'   + str(sigma_arg)       + r'}{\ampere\per\volt\per\meter}' + r'$}\\' + "\n",
        '        ' + r'\textcolor{red}{$\sigma_{Kuch}' + r'$} & \textcolor{red}{$' +  '\SI{'   + str(sigma_kuchling)  + r'}{\ampere\per\volt\per\meter}' + r'$}\\' + "\n",
        '        ' + '$\mu_0'   + '$ & $' +  '\SI{'   + str(mu0)    + r'}{\newton\per\ampere\squared}' + r'$\\' + "\n",
        '        ' + '$r_0'     + '$ & $' +  '\SI{'   + str(r0)     + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$r_{max}' + '$ & $' +  '\SI{'   + str(rmax)   + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$r_{min}' + '$ & $' +  '\SI{'   + str(rmin)   + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$B_0'     + '$ & $' +  '\SI{'   + str(B0)     + r'}{\tesla}'                     + r'$\\' + "\n",
        '        ' + '$NPTS'    + '$ & $' +  r'\num{' + str(npts)   + '}'                              + r'$\\' + "\n",
        '        ' + '$f'       + '$ & $' +  '\SI{'   + str(freq)   + r'}{\hertz}'                     + r'$\\' + "\n",
        ]
font = {
        'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 9,
        }
titlefont = {
        'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 10,
        }
plot_legend_fontsize    = 9
plot_color_fit          = 'blue'
plot_color_measurements = 'black'
plot_label_measurements = 'Messwerte'
plot_size_measurements  = 16
plot_scale_x            = 'linear'
plot_label_fit          = 'Fit-Funktion'
plot_label_x            = 'radiale Position bezogen auf Zylinderachse (mm)'
plot_1_label_y          = 'gemessene Spannung (mV)'
plot_2_label_y          = 'Phase (Grad)'
plot_1_title            = r"Exakte L\"osung: Betrag Magnetfeld Spule mit Vollzylinder (30 Hz)"
plot_2_title            = r"Exakte L\"osung: Phase Magnetfeld Spule mit Vollzylinder (30 Hz)"

loc = plticker.MultipleLocator(base=5)


# ---------------------------------------------------------#
# Functions                                                #
#                                                          #
# See formula 21 on p.11 of script for experiment.         #
#                                                          #
# NOTE: We use  frequency f  instead of  angular frequency #
# omega since that is what we actually set on the function #
# generator.                                               #
# NOTE: We evaluate B_abs and B_arg based on two different #
# values for sigma, which allows to fit each of the curves #
# more accurately.                                         #
# ---------------------------------------------------------#

k_abs = lambda f: sqrt((2*pi*f*mu0*sigma_abs)/2)*(mpc(1,-1))
k_arg = lambda f: sqrt((2*pi*f*mu0*sigma_arg)/2)*(mpc(1,-1))

# Enumerator:
enum_abs  = lambda r: besselj(0,k_abs(freq)*r)
denom_abs =           besselj(0,k_abs(freq)*r0)
enum_arg  = lambda r: besselj(0,k_arg(freq)*r)
denom_arg =           besselj(0,k_arg(freq)*r0)

B_abs = lambda r: abs(enum_abs(r) / denom_abs * B0)
B_arg = lambda r: arg(enum_arg(r) / denom_arg * B0)


# ---------------------------------------------------------#
# Generate points for radius axis                          #
# ---------------------------------------------------------#
radii = np.linspace(rmin,rmax,npts)


# ---------------------------------------------------------#
# Numerically evaluate function                            #
# ---------------------------------------------------------#
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
radii_measured = np.array([      0,      5,     10,    15,  20,    25,    30,    35,    40,    45,    50])
voltages       = np.array([2.86e-2,2.85e-2,2.87e-2,2.9e-2,3e-2,3.3e-2,3.8e-2,4.5e-2,5.4e-2,6.2e-2,3.7e-2])
phases_degrees = np.array([    111,    109,    104,    94,  81,    65,  48.5,    32,    16,   2.7,     0])


# ---------------------------------------------------------#
# Scale values for improved legibility in plot             #
# ---------------------------------------------------------#
# We scale from meters to millimeters,  from rad to degress.
B_abs_num = 1e3 * B_abs_num
radii     = 1e3 * radii
voltages  = 1e3 * voltages
B_arg_num = 180/pi*B_arg_num
rmin      = 1e3 * rmin
rmax      = 1e3 * rmax


# ---------------------------------------------------------#
# Plot the Things                                          #
# ---------------------------------------------------------#
matplotlib.pyplot.rc('text', usetex=True)
matplotlib.pyplot.rc('font', family='serif')

fig   = figure(1)
axes1 = fig.add_subplot(211)
axes1.plot(radii,B_abs_num,color=plot_color_fit,label=plot_label_fit)
axes1.scatter(radii_measured,
        voltages,
        color=plot_color_measurements,
        s=plot_size_measurements,
        label=plot_label_measurements
        )
axes1.set_xlim([rmin-1,rmax*1.1+5])
axes1.set_xscale(plot_scale_x)
axes1.set_xlabel(plot_label_x,fontdict=font)
axes1.set_ylabel(plot_1_label_y,fontdict=font)
axes1.set_title(plot_1_title,fontdict=titlefont)
axes1.legend(fontsize=plot_legend_fontsize,loc='upper left')
axes1.xaxis.set_major_locator(loc)
axes1.tick_params(labelsize=9)

axes2 = fig.add_subplot(212)
axes2.plot(radii,B_arg_num,color=plot_color_fit,label=plot_label_fit)
axes2.scatter(radii_measured,
        -phases_degrees,
        color=plot_color_measurements,
        s=plot_size_measurements,
        label=plot_label_measurements
        )
axes2.set_xlim([rmin-1,rmax*1.1+5])
axes2.set_xscale(plot_scale_x)
axes2.set_xlabel(plot_label_x,fontdict=font)
axes2.set_ylabel(plot_2_label_y,fontdict=font)
axes2.set_title(plot_2_title,fontdict=titlefont)
axes2.legend(fontsize=plot_legend_fontsize,loc='upper left')
axes2.xaxis.set_major_locator(loc)
axes2.tick_params(labelsize=9)

fig.subplots_adjust(bottom=0.1,left=0.1,right=0.9,top=0.95,hspace=0.5)

fig.savefig('plots-pgf/massive--alu--low-freq--exact.pgf')
fig.savefig('plots-pdf/massive--alu--low-freq--exact.pdf')


# ---------------------------------------------------------#
# Save listing to file                                     #
# ---------------------------------------------------------#
dumpfile = open('listings/massive--alu--low-freq--exact.tex', 'w')

table_opening = r"""
{%
    \begin{center}
    \captionof{table}{%
        Parameter f\"ur Fit-Funktion aus Abbildung~\ref{fig:alu:rad:low:sensor}
    }
    \label{tab:fitparams:alu:freq:low:exact}
    \sisetup{%
        %math-rm=\mathtt,
        scientific-notation=engineering,
        table-format = +3.2e+2,
        round-precision = 3,
        round-mode = figures,
    }
    \begin{tabular}{lr}
    \toprule
"""
table_closing = r"""
    \bottomrule
    \end{tabular}
    \end{center}
}

"""

dumpfile.writelines(table_opening)

for line in params:
    dumpfile.writelines(line)

dumpfile.writelines(table_closing)
dumpfile.close()


# ---------------------------------------------------------#
# Save Value of sigma to file for error analysis           #
# ---------------------------------------------------------#
np.savetxt('numpy-txt/massive--alu--low-freq--exact.txt',([sigma_abs,sigma_arg]))
