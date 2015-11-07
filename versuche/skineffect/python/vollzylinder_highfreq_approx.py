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
mu0   = 4*pi*1e-7
#sigma = 37.7e6                 # conductivity of aluminium (de.wikipedia.org)
sigma = 18e6                                                   # affects phase
B0    = 5.5e-2                 # does not affect phase, use for scaling abs(B)
r0    = 45e-3
freq = 450                                     # frequency was fixed at 450 Hz
npts = 1e3
rmin=25e-3
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
plot_11_label_y          = 'gemessene Spannung (mV)'
plot_21_label_y          = 'gemessene Spannung (mV)'
plot_12_label_y          = 'Phase (Grad)'
plot_11_title            = r"Hochfrequenzn\"aherung: Betrag Magnetfeld Spule mit Vollzylinder (450 Hz)"
plot_12_title            = r"Hochfrequenzn\"aherung: Phase Magnetfeld Spule mit Vollzylinder (450 Hz)"

# Set ticker intervals for plots (in millimeters)
loc1 = plticker.MultipleLocator(base=2.5)


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
radii = np.linspace(rmin,rmax,npts)


# ---------------------------------------------------------#
# Numerically evaluate function                            #
# ---------------------------------------------------------#
Babsufunc = np.frompyfunc(B_abs,1,1)
Bargufunc = np.frompyfunc(B_arg,1,1)
B_abs_num= Babsufunc(radii)
B_arg_num= Bargufunc(radii)


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

fig1   = figure(1)
axes11 = fig1.add_subplot(211)
axes11.plot(radii,B_abs_num,color=plot_color_fit,label=plot_label_fit)
axes11.scatter(radii_measured,
        voltages,
        color=plot_color_measurements,
        s=plot_size_measurements,
        label=plot_label_measurements
        )
axes11.set_xlim([rmin*0.9,rmax*1.1])
axes11.set_xscale(plot_scale_x)
axes11.set_xlabel(plot_label_x,fontdict=font)
axes11.set_ylabel(plot_11_label_y,fontdict=font)
axes11.set_title(plot_11_title,fontdict=font)
axes11.legend(fontsize=plot_legend_fontsize,loc='upper left')
axes11.xaxis.set_major_locator(loc1)

axes12 = fig1.add_subplot(212)
axes12.plot(radii,B_arg_num,color=plot_color_fit,label=plot_label_fit)
axes12.scatter(radii_measured,
        -phases_degrees,
        color=plot_color_measurements,
        s=plot_size_measurements,
        label=plot_label_measurements
        )
axes12.set_xlim([rmin*0.9,rmax*1.1])
axes12.set_xscale(plot_scale_x)
axes12.set_xlabel(plot_label_x,fontdict=font)
axes12.set_ylabel(plot_12_label_y,fontdict=font)
axes12.set_title(plot_12_title,fontdict=font)
axes12.legend(fontsize=plot_legend_fontsize,loc='upper left')
axes12.xaxis.set_major_locator(loc1)

fig1.subplots_adjust(bottom=0.1,left=0.1,right=0.9,top=0.95,hspace=0.5)

fig1.savefig('plots-pgf/massive--alu--high-freq-approx--high.pgf')
fig1.savefig('plots-pdf/massive--alu--high-freq-approx--high.pdf')


# ---------------------------------------------------------#
# Save listing to file                                     #
# ---------------------------------------------------------#
dumpfile = open('listings/massive--alu--high-freq-approx--high.tex', 'w')

table_opening = r"""
{%
    \begin{center}
    \captionof{table}{%
        Paramaterwerte f\"ur Fitfunktion basierend auf der N\"aherungsl\"osung
        f\"ur hohe Frequenzen.
    }
    \label{tab:fitparams:alu:freq:approx:high}
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
