#!/usr/bin/env python3

from sympy import *
from mpmath import *
from matplotlib.pyplot import *
import matplotlib.ticker as plticker
#init_printing()     # make things prettier when we print stuff for debugging.


# ************************************************************************** #
# B-Field, Cylinder Coild with Massive Alu Cylinder                          #
# ************************************************************************** #

# All values are in standard SI units unless otherwise noted.


# ---------------------------------------------------------#
# Define Variables and Constants                           #
# ---------------------------------------------------------#
mu0   = 4*pi*1e-7
#sigma = 37.7e6                 # conductivity of aluminium (de.wikipedia.org)
sigma = 18e6                                                   # affects phase
B0    = 5.5e-2                 # does not affect phase, use for scaling abs(B)
r0    = 45e-3
freq = 450                                     # frequency was fixed at 450 Hz
npts = 1e3
rmin=25e-3
#rmax=50e-3
rmax=45e-3
    # -----------------------------------------------------#
    # Create  a list  for convenient  printing of  vars to #
    # file, add LaTeX where necessary.                     #
    # -----------------------------------------------------#
params = [
        '        ' + '$\mu_0'   + '$ & $' +  '\SI{'   + str(mu0)    + r'}{\newton\per\ampere\squared}' + r'$\\' + "\n",
        '        ' + '$\sigma'  + '$ & $' +  '\SI{'   + str(sigma)  + r'}{\ampere\per\volt\per\meter}' + r'$\\' + "\n",
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
        'size'   : 11,
        }
plot_legend_fontsize    = 11
plot_color_fit          = 'blue'
plot_color_measurements = 'black'
plot_label_measurements = 'Messwerte'
plot_size_measurements  = 32
plot_scale_x            = 'linear'
plot_label_fit          = 'Fitfunktion'
plot_label_x            = 'radiale Position bezogen auf Zylinderachse (mm)'
plot_1_label_y          = 'gemessene Spannung (mV)'
plot_2_label_y          = 'Phase (Grad)'
plot_1_title            = r"Exakte L\"osung: Betrag Magnetfeld Spule mit Vollzylinder (450 Hz)"
plot_2_title            = r"Exakte L\"osung: Phase Magnetfeld Spule mit Vollzylinder (450 Hz)"
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

fig.subplots_adjust(bottom=0.1,left=0.1,right=0.9,top=0.95,hspace=0.5)

fig.savefig('plots-pgf/massive--alu--high-freq--exact.pgf')
fig.savefig('plots-pdf/massive--alu--high-freq--exact.pdf')


# ---------------------------------------------------------#
# Save listing to file                                     #
# ---------------------------------------------------------#
dumpfile = open('listings/massive--alu--high-freq--exact.tex', 'w')

table_opening = r"""
{%
    \begin{center}
    \captionof{table}{%
        Parameter        f\"ur         Fitfunktion        aus        Abbildung
        \ref{fig:alu:rad:exact:high}.
    }
    \label{tab:fitparams:alu:freq:high:exact}
    \sisetup{%
        %math-rm=\mathtt,
        scientific-notation=engineering,
        table-format = +3.2e+2,
        round-precision = 2,
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
