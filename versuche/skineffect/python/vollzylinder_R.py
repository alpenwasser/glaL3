#!/usr/bin/env python3

from sympy import *
from mpmath import *
from matplotlib.pyplot import *
#init_printing()     # make things prettier when we print stuff for debugging.


# ************************************************************************** #
# Self-Inductance L of copper coil with massive aluminium cylinder inserted  #
# ************************************************************************** #

# All values are in standard SI units unless otherwise noted.

# ---------------------------------------------------------#
# Define Variables and Constants                           #
# ---------------------------------------------------------#
mu0   = 4*pi*1e-7                                        # vacuum permeability
#sigma = 37.7e6                 # conductivity of aluminium (de.wikipedia.org)
sigma = 23.75e6                         # de.wikipedia.org/wiki/Kupfer: 58.1e6
r     = 0             # radial position of measurement probe. Centered on axis
dsp   = 98e-3                                               # diameter of coil
rsp   = dsp / 2                                               # radius of coil
r0    = 45e-3                                         # radius of alu cylinder
B0    = 6.9e-2                             # adjust this as needed for scaling
N0    = 574                                   # number of turns of copper coil
l     = 500e-3                                         # length of copper coil
R_0   = 6.3             # resistance of copper wire for coil, calculated value
npts  = 1e3
fmin  = 1
fmax  = 250
    # -----------------------------------------------------#
    # Create  a list  for convenient  printing of  vars to #
    # file, add LaTeX where necessary.                     #
    # -----------------------------------------------------#
params = [
        '        ' + '$\mu_0'        + '$ & $' +  '\SI{'   + str(mu0)    + r'}{\newton\per\ampere\squared}' + r'$\\' + "\n",
        '        ' + '$\sigma'       + '$ & $' +  '\SI{'   + str(sigma)  + r'}{\ampere\per\volt\per\meter}' + r'$\\' + "\n",
        '        ' + '$d_{Sp}'       + '$ & $' +  '\SI{'   + str(dsp)    + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$r_{Sp}'       + '$ & $' +  '\SI{'   + str(rsp)    + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$r'            + '$ & $' +  '\SI{'   + str(r)      + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$r_{0}'        + '$ & $' +  '\SI{'   + str(r0)     + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$B_0'          + '$ & $' +  '\SI{'   + str(B0)     + r'}{\tesla}'                     + r'$\\' + "\n",
        '        ' + '$l'            + '$ & $' +  '\SI{'   + str(l)      + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$R_{\Omega,0}' + '$ & $' +  '\SI{'   + str(R_0)    + '}{\ohm}'                        + r'$\\' + "\n",
        '        ' + '$NPTS'         + '$ & $' +  r'\num{' + str(npts)   + '}'                              + r'$\\' + "\n",
        '        ' + '$N_0'          + '$ & $' +  r'\num{' + str(N0)     + '}'                              + r'$\\' + "\n",
        '        ' + '$f_{min}'      + '$ & $' +  '\SI{'   + str(fmin)   + r'}{\hertz}'                     + r'$\\' + "\n",
        '        ' + '$f_{max}'      + '$ & $' +  '\SI{'   + str(fmax)   + r'}{\hertz}'                     + r'$\\' + "\n",
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
plot_color_fit = 'blue'
plot_linewidth = 1
plot_scale_x   = 'log'
plot_label_x   = 'Frequenz (Hz)'
plot_label_y   = 'Widerstand (Ohm)'
plot_title     = r"Ohm'scher Widerstand, Spule mit Vollzylinder"


# ---------------------------------------------------------#
# Functions                                                #
#                                                          #
# See formula 23 on p.12 of script for experiment.         #
#                                                          #
# NOTE: We use  frequency f  instead of  angular frequency #
# omega since that is what we actually set on the function #
# generator.                                               #
# ---------------------------------------------------------#
k = lambda f: sqrt((2*np.pi*f*mu0*sigma)/2)*(mpc(1,-1))

R_tot = lambda f:(
        -(f * mu0 * 2 * pi * r0 * N0**2) / l
        * im(besselj(1,k(f) * r0) / (k(f) * besselj(0,k(f) * r0)))
        + R_0
    )


# ---------------------------------------------------------#
# Generate points for frequency axis                       #
# ---------------------------------------------------------#
n                = np.linspace(1,npts,npts)
expufunc         = np.frompyfunc(exp,1,1)
frequency_vector = fmin*expufunc(n*log(fmax-fmin)/npts)


# ---------------------------------------------------------#
# Numerically evaluate function                            #
# ---------------------------------------------------------#
Rufunc = np.frompyfunc(R_tot,1,1)
R_num  = Rufunc(frequency_vector)


# ---------------------------------------------------------#
# Plot the Things                                          #
# ---------------------------------------------------------#
matplotlib.pyplot.rc('text', usetex=True)
matplotlib.pyplot.rc('font', family='serif')

figwidth = 8.27 # in inches
fig  = figure(1,figsize=(figwidth,figwidth*0.36))
axes = fig.add_subplot(111)
axes.plot(frequency_vector,R_num,linewidth=plot_linewidth,color=plot_color_fit)
axes.set_xscale(plot_scale_x)
axes.set_xlim([fmin*0.9,fmax*1.1])
axes.set_xlabel(plot_label_x,fontdict=font)
axes.set_ylabel(plot_label_y,fontdict=font)
axes.set_title(plot_title,fontdict=font)
axes.tick_params(labelsize=9)

fig.subplots_adjust(bottom=0.15,left=0.125,right=0.925,top=0.90)

fig.savefig('plots-pgf/massive--alu--R.pgf')
fig.savefig('plots-pdf/massive--alu--R.pdf')


# ---------------------------------------------------------#
# Save listing to file                                     #
# ---------------------------------------------------------#
dumpfile = open('listings/massive--alu--R.tex', 'w')

table_opening = r"""
{%
    \begin{center}
    \captionof{table}{%
        Parameterwerte  f\"ur  Fit-Funktion  in  Abbildung  \ref{fig:alu:freq:R}
    }
    \label{tab:fitparams:alu:R}
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
