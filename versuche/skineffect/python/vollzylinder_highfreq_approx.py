#!/usr/bin/env python3

from sympy import *
from mpmath import *
from matplotlib.pyplot import *
import matplotlib.ticker as plticker


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
rho_kuchling   = 0.027e-6  # resistivity Kuchling 17th edition, p.649, tab. 45
sigma_kuchling = 1/rho_kuchling
#sigma = 37.7e6                 # conductivity of aluminium (de.wikipedia.org)
sigma_abs = 19e6                                                   # affects phase
sigma_arg = 18e6                                                   # affects phase
B0    = 6.0e-2                 # does not affect phase, use for scaling abs(B)
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
plot_11_label_y          = 'gemessene Spannung (mV)'
plot_21_label_y          = 'gemessene Spannung (mV)'
plot_12_label_y          = 'Phase (Grad)'
plot_11_title            = r"Hochfrequenzn\"aherung: Betrag Magnetfeld Spule mit Vollzylinder (450 Hz)"
plot_12_title            = r"Hochfrequenzn\"aherung: Phase Magnetfeld Spule mit Vollzylinder (450 Hz)"

# Set ticker intervals for plots (in millimeters)
loc1 = plticker.MultipleLocator(base=2.5)


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

s_skin_abs = sqrt(2/(2*pi*freq*mu0*sigma_abs))
s_skin_arg = sqrt(2/(2*pi*freq*mu0*sigma_arg))
x = lambda r: r0-r
B_abs = lambda r: abs(B0 * exp(-x(r)/s_skin_abs) * exp(mpc(0,-x(r)/s_skin_abs)))
B_arg = lambda r: arg(B0 * exp(-x(r)/s_skin_arg) * exp(mpc(0,-x(r)/s_skin_arg)))


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
axes11.set_xlim([rmin*0.9,rmax*1.1+2.5])
axes11.set_xscale(plot_scale_x)
axes11.set_xlabel(plot_label_x,fontdict=font)
axes11.set_ylabel(plot_11_label_y,fontdict=font)
axes11.set_title(plot_11_title,fontdict=titlefont)
axes11.legend(fontsize=plot_legend_fontsize,loc='upper left')
axes11.xaxis.set_major_locator(loc1)
axes11.tick_params(labelsize=9)

axes12 = fig1.add_subplot(212)
axes12.plot(radii,B_arg_num,color=plot_color_fit,label=plot_label_fit)
axes12.scatter(radii_measured,
        -phases_degrees,
        color=plot_color_measurements,
        s=plot_size_measurements,
        label=plot_label_measurements
        )
axes12.set_xlim([rmin*0.9,rmax*1.1+2.5])
axes12.set_xscale(plot_scale_x)
axes12.set_xlabel(plot_label_x,fontdict=font)
axes12.set_ylabel(plot_12_label_y,fontdict=font)
axes12.set_title(plot_12_title,fontdict=titlefont)
axes12.legend(fontsize=plot_legend_fontsize,loc='upper left')
axes12.xaxis.set_major_locator(loc1)
axes12.tick_params(labelsize=9)

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
        Parameterwerte            f\"ur             Fit-Funktion            aus
        Abbildung~\ref{fig:alu:rad:approx:high}
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
