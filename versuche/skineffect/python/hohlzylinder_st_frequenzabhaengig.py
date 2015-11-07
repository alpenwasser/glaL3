#!/usr/bin/env python3

from sympy import *
from mpmath import *
from matplotlib.pyplot import *
#init_printing()     # make things prettier when we print stuff for debugging.


# ************************************************************************** #
# Magnetic field inside copper coil with  hollow copper cylinder             #
# ************************************************************************** #

# All values are in standard SI units unless otherwise noted.


# ---------------------------------------------------------#
# Init, Define Variables and Constants                     #
# ---------------------------------------------------------#
mu0   = 4*pi*1e-7                                        # vacuum permeability
sigma = 1.25e6  # adjust as needed until function plot fits measurement values
r     = 0             # radial position of measurement probe. Centered on axis
dsp   = 98e-3                                               # diameter of coil
rsp   = dsp / 2                                               # radius of coil
r1    = 30e-3                                # inner radius of copper cylinder
r2    = 35e-3                                # outer radius of copper cylinder
B0    = 6.9e-2                             # adjust this as needed for scaling
npts  = 1e3                                  # number of points for plot curve
fmin  = 1                   # minimum frequency, adjusted to measurement range
fmax  = 7.5e3               # maximum frequency, adjusted to measurement range
fmin_opt = 8e1                # minimum frequency, adjusted to get a nice plot
fmax_opt = 5e4                # maximum frequency, adjusted to get a nice plot
    # -----------------------------------------------------#
    # Create  a list  for convenient  printing of  vars to #
    # file, add LaTeX where necessary.                     #
    # -----------------------------------------------------#
params = [
        '        ' + '$\mu_0'         + '$ & $' +  '\SI{'   + str(mu0)        + r'}{\newton\per\ampere\squared}' + r'$\\' + "\n",
        '        ' + '$\sigma'        + '$ & $' +  '\SI{'   + str(sigma)      + r'}{\ampere\per\volt\per\meter}' + r'$\\' + "\n",
        '        ' + '$r'             + '$ & $' +  '\SI{'   + str(r)          + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$r_1'           + '$ & $' +  '\SI{'   + str(r1)         + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$r_2'           + '$ & $' +  '\SI{'   + str(r2)         + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$B_0'           + '$ & $' +  '\SI{'   + str(B0)         + r'}{\tesla}'                     + r'$\\' + "\n",
        '        ' + '$NPTS'          + '$ & $' +  r'\num{' + str(npts)       + '}'                              + r'$\\' + "\n",
        '        ' + '$f_{min}'       + '$ & $' +  '\SI{'   + str(fmin)       + r'}{\hertz}'                     + r'$\\' + "\n",
        '        ' + '$f_{max}'       + '$ & $' +  '\SI{'   + str(fmax)       + r'}{\hertz}'                     + r'$\\' + "\n",
        ]
params_opt = [
        '        ' + '$\mu_0'         + '$ & $' +  '\SI{'   + str(mu0)        + r'}{\newton\per\ampere\squared}' + r'$\\' + "\n",
        '        ' + '$\sigma'        + '$ & $' +  '\SI{'   + str(sigma)      + r'}{\ampere\per\volt\per\meter}' + r'$\\' + "\n",
        '        ' + '$r'             + '$ & $' +  '\SI{'   + str(r)          + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$r_1'           + '$ & $' +  '\SI{'   + str(r1)         + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$r_2'           + '$ & $' +  '\SI{'   + str(r2)         + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$B_0'           + '$ & $' +  '\SI{'   + str(B0)         + r'}{\tesla}'                     + r'$\\' + "\n",
        '        ' + '$NPTS'          + '$ & $' +  r'\num{' + str(npts)       + '}'                              + r'$\\' + "\n",
        '        ' + '$f_{min_{opt}}' + '$ & $' +  '\SI{'   + str(fmin_opt)   + r'}{\hertz}'                     + r'$\\' + "\n",
        '        ' + '$f_{max_{opt}}' + '$ & $' +  '\SI{'   + str(fmax_opt)   + r'}{\hertz}'                     + r'$\\' + "\n",
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
plot_scale_x            = 'log'
plot_label_fit          = 'Fitfunktion'
plot_label_x            = 'Frequenz (Hz)'
plot_1_label_y          = 'gemessene Spannung (mV)'
plot_2_label_y          = 'Phase (Grad)'
plot_1_title            = r"Exakte L\"osung: Betrag Magnetfeld, Spule mit Stahlrohr"
plot_2_title            = r"Exakte L\"osung: Phase Magnetfeld, Spule mit Stahlrohr"

    # -----------------------------------------------------#
    # NOTE: According to  formula 26 on p.14,  the B-Field #
    # inside the  cylinder (r<r1) is equal  to the B-Field #
    # at  the  inner  boundary   of  the  copper  cylinder #
    # (B(r1)),  therefore  we  set  r to  r1  for  further #
    # calculations.                                        #
    # -----------------------------------------------------#
r     = 30e-3


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

k = lambda f: sqrt((2*np.pi*f*mu0*sigma)/2)*(mpc(1,-1))

enum  = lambda f:(
          besselj(0,k(f)*r)
        * bessely(2,k(f)*r1)
        - besselj(2,k(f)*r1)
        * bessely(0,k(f)*r)
    )
denom = lambda f:(
          besselj(0,k(f)*r2)
        * bessely(2,k(f)*r1)
        - besselj(2,k(f)*r1)
        * bessely(0,k(f)*r2)
    )

B = lambda f: enum(f) / denom(f) * B0

B_abs = lambda f: abs(B(f))
B_arg = lambda f: arg(B(f))


# ---------------------------------------------------------#
# Generate points for frequency axis                       #
# ---------------------------------------------------------#
n                = np.linspace(0,npts,npts)
expufunc         = np.frompyfunc(exp,1,1)
#frequency_vector = fmin*expufunc(n*log(fmax_opt-fmin)/npts)
frequency_vector     = expufunc((1-n/npts)*log(fmin)) * expufunc(n*log(fmax)/npts)
frequency_vector_opt = expufunc((1-n/npts)*log(fmin_opt)) * expufunc(n*log(fmax_opt)/npts)


# ---------------------------------------------------------#
# Numerically evaluate functions                           #
# ---------------------------------------------------------#
Babsufunc = np.frompyfunc(B_abs,1,1)
B_abs_num = Babsufunc(frequency_vector)
B_abs_num_opt = Babsufunc(frequency_vector_opt)
Bargufunc = np.frompyfunc(B_arg,1,1)
B_arg_num = Bargufunc(frequency_vector)
B_arg_num_opt = Bargufunc(frequency_vector_opt)


# ---------------------------------------------------------#
# Unfortunately, the  arg() function only  delivers values #
# between -pi and  +pi for the angle of  a complex number, #
# which,  while  correct,  is   not  suitable  for  pretty #
# plotting, so we  will shift the values  larger then zero #
# accordingly for a continuous curve.                      #
# ---------------------------------------------------------#
B_arg_num = np.unwrap(B_arg_num)
B_arg_num_opt = np.unwrap(B_arg_num_opt)


# ---------------------------------------------------------#
# Measurement Values from experiment                       #
# ---------------------------------------------------------#
frequencies_measured = np.array([       1,      10,      20,      40,      80,     120,     160,     200,     400,     600,    800,    1000,    1200,   1500,   1750,    2000,   2500,   3500,   5000,   7500])
phases_degrees       = np.array([       0,    0.45,    0.95,     1.8,     3.6,     5.4,     7.2,       9,    17.5,    25.4,   32.4,    38.4,    43.5,     50,     54,      58,     64,     71,     78,     88])
voltages             = np.array([ 6.96e-2, 6.97e-2, 6.97e-2, 6.97e-2, 6.92e-2, 6.91e-2, 6.87e-2, 6.62e-2, 6.27e-2, 6.27e-2, 5.9e-2, 5.45e-2, 5.05e-2, 4.5e-2, 4.1e-2, 3.72e-2, 3.2e-2, 2.4e-2, 1.8e-2, 1.2e-2])


# ---------------------------------------------------------#
# Scale values for improved legibility in plot             #
# ---------------------------------------------------------#
voltages  = 1e3 * voltages

B_abs_num = 1e3 * B_abs_num
B_abs_num_opt = 1e3 * B_abs_num_opt
B_arg_num = 180/np.pi*B_arg_num
B_arg_num_opt = 180/np.pi*B_arg_num_opt


# ---------------------------------------------------------#
# Plot the Things                                          #
# ---------------------------------------------------------#
matplotlib.pyplot.rc('text', usetex=True)
matplotlib.pyplot.rc('font', family='serif')

# Figure 1: Range includes all measurement values.
fig1 = figure(1)
axes11 = fig1.add_subplot(211)
axes11.plot(frequency_vector,B_abs_num,color=plot_color_fit,label=plot_label_fit)
axes11.scatter(frequencies_measured,
        voltages,
        color=plot_color_measurements,
        s=plot_size_measurements,
        label=plot_label_measurements
        )
axes11.set_xlim([0.9*fmin,1.1*fmax])
axes11.set_xscale(plot_scale_x)
axes11.set_xlabel(plot_label_x,fontdict=font)
axes11.set_ylabel(plot_1_label_y,fontdict=font)
axes11.set_title(plot_1_title,fontdict=font)
axes11.legend(fontsize=plot_legend_fontsize)

axes12 = fig1.add_subplot(212)
axes12.plot(frequency_vector,B_arg_num,color=plot_color_fit,label=plot_label_fit)
axes12.scatter(frequencies_measured,
        -phases_degrees,
        color=plot_color_measurements,
        s=plot_size_measurements,
        label=plot_label_measurements
        )
axes12.set_xlim([0.9*fmin,1.1*fmax])
axes12.set_xscale(plot_scale_x)
axes12.set_xlabel(plot_label_x,fontdict=font)
axes12.set_ylabel(plot_2_label_y,fontdict=font)
axes12.set_title(plot_2_title,fontdict=font)
axes12.legend(fontsize=plot_legend_fontsize)

# Figure 2: Optimized range for a more interesting plot.
fig2 = figure(2)
axes21 = fig2.add_subplot(2,1,1)
axes21.plot(frequency_vector_opt,B_abs_num_opt,color=plot_color_fit,label=plot_label_fit)
axes21.scatter(frequencies_measured,
        voltages,
        color=plot_color_measurements,
        s=plot_size_measurements,
        label=plot_label_measurements
        )
axes21.set_xlim([0.9*fmin_opt,fmax_opt*1.1])
axes21.set_xscale(plot_scale_x)
axes21.set_xlabel(plot_label_x,fontdict=font)
axes21.set_ylabel(plot_1_label_y,fontdict=font)
axes21.set_title(plot_1_title,fontdict=font)
axes21.legend(fontsize=plot_legend_fontsize)

axes22 = fig2.add_subplot(2,1,2)
axes22.plot(frequency_vector_opt,B_arg_num_opt,color=plot_color_fit,label=plot_label_fit)
axes22.scatter(frequencies_measured,
        -phases_degrees,
        color=plot_color_measurements,
        s=plot_size_measurements,
        label=plot_label_measurements
        )
axes22.set_xlim([0.9*fmin_opt,fmax_opt*1.1])
axes22.set_xscale(plot_scale_x)
axes22.set_xlabel(plot_label_x,fontdict=font)
axes22.set_ylabel(plot_2_label_y,fontdict=font)
axes22.set_title(plot_2_title,fontdict=font)
axes22.legend(fontsize=plot_legend_fontsize)


fig1.subplots_adjust(bottom=0.1,left=0.1,right=0.9,top=0.95,hspace=0.5)

fig1.savefig('plots-pgf/hollow--st--freq--exact.pgf')
fig1.savefig('plots-pdf/hollow--st--freq--exact.pdf')


fig2.subplots_adjust(bottom=0.1,left=0.1,right=0.9,top=0.95,hspace=0.5)

fig2.savefig('plots-pgf/hollow--st--freq--exact-opt.pgf')
fig2.savefig('plots-pdf/hollow--st--freq--exact-opt.pdf')


# ---------------------------------------------------------#
# Save listing to file                                     #
# ---------------------------------------------------------#
dumpfile = open('listings/hollow--st--freq--exact.tex', 'w')

table_opening = r"""
{%
    \begin{center}
    \captionof{table}{%
        Parameter f\"ur Fitfunktion in Abbildung~\ref{fig:st:freq:exact}
    }
    \label{tab:fitparams:cu:freq:exact}
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

dumpfile_opt = open('listings/hollow--st--freq--exact-opt.tex', 'w')

table_opening = r"""
{%
    \begin{center}
    \captionof{table}{%
        Parameter f\"ur Fitfunktion in Abbildung~\ref{fig:st:freq:exact:opt}
    }
    \label{tab:fitparams:st:freq:exact:opt}
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

dumpfile_opt.writelines(table_opening)

for line in params_opt:
    dumpfile_opt.writelines(line)

dumpfile_opt.writelines(table_closing)
dumpfile_opt.close()
