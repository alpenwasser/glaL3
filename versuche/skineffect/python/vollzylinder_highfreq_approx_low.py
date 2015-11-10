#!/usr/bin/env python3

from sympy import *
from mpmath import *
from matplotlib.pyplot import *
import matplotlib.ticker as plticker
#init_printing()     # make things prettier when we print stuff for debugging.


# ************************************************************************** #
# B-Field,   Cylinder  Coil   with  Massive   Alu  Cylinder                  #
#                                                                            #
# High frequenzy  approximation applied to low  frequency measurements. This #
# should not give an accurate match between curve and measurement points.    #
# ************************************************************************** #

# All values are in standard SI units unless otherwise noted.


# ---------------------------------------------------------#
# Define Variables and Constants                           #
# ---------------------------------------------------------#
mu0   = 4*pi*1e-7
rho_kuchling   = 0.027e-6  # resistivity Kuchling 17th edition, p.649, tab. 45
sigma_kuchling = 1/rho_kuchling
#sigma = 37.7e6                 # conductivity of aluminium (de.wikipedia.org)
#sigma = 18e6                                                   # affects phase
sigma = 17e6                                                   # affects phase
B0    = 6.2e-2                 # does not affect phase, use for scaling abs(B)
r0    = 45e-3
freq = 30                                       # frequency was fixed at 30 Hz
npts = 1e3
rmin=0
rmax=45e-3
    # -----------------------------------------------------#
    # Create  a list  for convenient  printing of  vars to #
    # file, add LaTeX where necessary.                     #
    # -----------------------------------------------------#
params = [
        '        ' + r'\textcolor{red}{$\sigma_{Fit}'  + r'$} & \textcolor{red}{$' +  '\SI{'   + str(sigma)           + r'}{\ampere\per\volt\per\meter}' + r'$}\\' + "\n",
        '        ' + r'\textcolor{red}{$\sigma_{Kuch}' + r'$} & \textcolor{red}{$' +  '\SI{'   + str(sigma_kuchling)  + r'}{\ampere\per\volt\per\meter}' + r'$}\\' + "\n",
        '        ' + '$\mu_0'   + '$ & $' +  '\SI{'   + str(mu0)    + r'}{\newton\per\ampere\squared}' + r'$\\' + "\n",
        '        ' + '$\sigma'  + '$ & $' +  '\SI{'   + str(sigma)  + r'}{\ampere\per\volt\per\meter}' + r'$\\' + "\n",
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
plot_11_label_y          = 'gemessene Spannung (mV)'
plot_21_label_y          = 'gemessene Spannung (mV)'
plot_12_label_y          = 'Phase (Grad)'
plot_21_title            = r"Hochfrequenzn\"aherung: Betrag Magnetfeld Spule mit Vollzylinder (30 Hz)"
plot_22_title            = r"Hochfrequenzn\"aherung: Phase Magnetfeld Spule mit Vollzylinder (30 Hz)"

# Set ticker intervals for plots (in millimeters)
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
radii = np.linspace(rmin,rmax,npts)


# ---------------------------------------------------------#
# Numerically evaluate function                            #
# ---------------------------------------------------------#
Babsufunc = np.frompyfunc(B_abs,1,1)
Bargufunc = np.frompyfunc(B_arg,1,1)
B_abs_num = Babsufunc(radii)
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

fig2   = figure(2)
axes21 = fig2.add_subplot(211)
axes21.plot(radii,B_abs_num,color=plot_color_fit,label=plot_label_fit)
axes21.scatter(radii_measured,
        voltages,
        color=plot_color_measurements,
        s=plot_size_measurements,
        label=plot_label_measurements
        )
axes21.set_xlim([rmin-5,rmax*1.1])
axes21.set_xscale(plot_scale_x)
axes21.set_xlabel(plot_label_x,fontdict=font)
axes21.set_ylabel(plot_21_label_y,fontdict=font)
axes21.set_title(plot_21_title,fontdict=titlefont)
axes21.legend(fontsize=plot_legend_fontsize,loc='upper left')
axes21.xaxis.set_major_locator(loc2)
axes21.tick_params(labelsize=9)

axes22 = fig2.add_subplot(212)
axes22.plot(radii,B_arg_num,color=plot_color_fit,label=plot_label_fit)
axes22.scatter(radii_measured,
        -phases_degrees,
        color=plot_color_measurements,
        s=plot_size_measurements,
        label=plot_label_measurements
        )
axes22.set_xlim([rmin-5,rmax*1.1])
axes22.set_xscale(plot_scale_x)
axes22.set_xlabel(plot_label_x,fontdict=font)
axes22.set_ylabel(plot_12_label_y,fontdict=font)
axes22.set_title(plot_22_title,fontdict=titlefont)
axes22.legend(fontsize=plot_legend_fontsize,loc='upper left')
axes22.xaxis.set_major_locator(loc2)
axes22.tick_params(labelsize=9)


fig2.subplots_adjust(bottom=0.1,left=0.1,right=0.9,top=0.95,hspace=0.5)

fig2.savefig('plots-pgf/massive--alu--high-freq-approx--low.pgf')
fig2.savefig('plots-pdf/massive--alu--high-freq-approx--low.pdf')


# ---------------------------------------------------------#
# Save listing to file                                     #
# ---------------------------------------------------------#
dumpfile = open('listings/massive--alu--high-freq-approx--low.tex', 'w')

table_opening = r"""
{%
    \begin{center}
    \captionof{table}{%
        Paramaterwerte       f\"ur       Fit-Funktion       aus       Abbildung
        \ref{fig:alu:rad:approx:low}
    }
    \label{tab:fitparams:alu:rad:approx:low}
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
