#!/usr/bin/env python3

from sympy import *
from mpmath import *
from matplotlib.pyplot import *
#init_printing()     # make things prettier when we print stuff for debugging.


# ************************************************************************** #
# Magnetic Flow normed  by current, copper coil with  hollow copper cylinder #
# inserted.                                                                  #
# ************************************************************************** #

# All values are in standard SI units unless otherwise noted.


# ---------------------------------------------------------#
# Init, Define Variables and Constants                     #
# ---------------------------------------------------------#
mu0   = 4*pi*1e-7                                        # vacuum permeability
sigma = 1.25e6                               # fit parameter: adjust as needed
r     = 0             # radial position of measurement probe. Centered on axis
dsp   = 98e-3                                               # diameter of coil
rsp   = dsp / 2                                               # radius of coil
r1    = 30e-3                                # inner radius of copper cylinder
r2    = 35e-3                                # outer radius of copper cylinder
B0    = 6.9e-2                          # fit parameter: adjust this as needed
N0    = 574                                   # number of turns of copper coil
l     = 500e-3                                         # length of copper coil
npts  = 1e3
fmin  = 8e1
fmax  = 5e4
    # -----------------------------------------------------#
    # Create  a list  for convenient  printing of  vars to #
    # file, add LaTeX where necessary.                     #
    # -----------------------------------------------------#
params = [
        '        ' + '$\mu_0'        + '$ & $' +  '\SI{'   + str(mu0)    + r'}{\newton\per\ampere\squared}' + r'$\\' + "\n",
        '        ' + '$\sigma'       + '$ & $' +  '\SI{'   + str(sigma)  + r'}{\ampere\per\volt\per\meter}' + r'$\\' + "\n",
        '        ' + '$d_{Sp}'       + '$ & $' +  '\SI{'   + str(dsp)    + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$r_{Sp}'       + '$ & $' +  '\SI{'   + str(rsp)    + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$r_1'          + '$ & $' +  '\SI{'   + str(r1)     + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$r_2'          + '$ & $' +  '\SI{'   + str(r2)     + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$l'            + '$ & $' +  '\SI{'   + str(l)      + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$NPTS'         + '$ & $' +  r'\num{' + str(npts)   + '}'                              + r'$\\' + "\n",
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
plot_legend_fontsize    = 9
plot_color_fit          = 'blue'
plot_color_fit_approx   = 'magenta'
plot_scale_x            = 'log'
plot_label_fit          = 'Fit-Funktion'
plot_label_fit_approx   = r'Fit-Funktion, N\"aherungsl\"osung'
plot_label_x            = 'Frequenz (Hz)'
plot_1_label_y          = r"$\displaystyle \biggl| \frac{\Phi}{I} \biggr|$ $\biggl( \displaystyle \frac{Vs}{A} \biggr)$"
plot_2_label_y          = r"$\displaystyle arg\biggl( \frac{\Phi}{I} \biggr)$ (Grad)"
plot_1_title            = r"Betrag Magn. Fluss normiert auf Spulenstrom, Spule mit Stahlrohr"
plot_2_title            = r"Phase Magn. Fluss normiert auf Spulenstrom, Spule mit Stahlrohr"


# ---------------------------------------------------------#
# Functions                                                #
#                                                          #
# For  the  exact formulas,  see  formula  28 on  p.15  of #
# script, for the approximations see formula 30 on p.16 of #
# script.                                                  #
#                                                          #
# NOTE: We use  frequency f  instead of  angular frequency #
# omega since that is what we actually set on the function #
# generator.                                               #
# ---------------------------------------------------------#

k = lambda f: sqrt((2*np.pi*f*mu0*sigma)/2)*(mpc(1,-1))

# exact solution:
enum1 = lambda f:(
          besselj(0,k(f)*r1)
        * bessely(2,k(f)*r1)
        - besselj(2,k(f)*r1)
        * bessely(0,k(f) * r1)
    )
denom1 = lambda f: (
          besselj(0,k(f)*r2)
        * bessely(2,k(f)*r1)
        - besselj(2,k(f)*r1)
        * bessely(0,k(f) * r2)
    )
enum2 = lambda f:(
          r2 * (
            besselj(1,k(f)*r2)
            * bessely(2,k(f)*r1)
            - besselj(2,k(f)*r1)
            * bessely(1,k(f) * r2)
        )
        - r1 * (
            besselj(1,k(f)*r1)
            * bessely(2,k(f)*r1)
            - besselj(2,k(f)*r1)
            * bessely(1,k(f) * r1)
        )
    )
denom2 = lambda f: (
          besselj(0,k(f)*r2)
        * bessely(2,k(f)*r1)
        - besselj(2,k(f)*r1)
        * bessely(0,k(f) * r2)
    )
term3 = rsp ** 2 - r2**2
prefactor = mu0 * pi * N0**2 / l

phi_norm = lambda f:(
        prefactor * (
            r1**2    * enum1(f)/denom1(f)
            + 2/k(f) * enum2(f)/denom2(f)
            + term3
        )
    )

phi_norm_abs = lambda f: abs(phi_norm(f))
phi_norm_arg = lambda f: arg(phi_norm(f))

# Approx. solution:
u1 = lambda f: mpc(0,1) * k(f) * r1
u2 = lambda f: mpc(0,1) * k(f) * r2

enum_approx = lambda f: (
          (u1(f) / 2 + 1) * ((u2(f) - 1) * exp( u2(f) - u1(f)) - (u1(f) - 1))
        + (u1(f) / 2 - 1) * ((u2(f) + 1) * exp(-u2(f) + u1(f)) - (u1(f) + 1))
        )
denom_approx = lambda f: (
        (u1(f) / 2) * exp(u2(f) - u1(f)) - (u1(f) / 2 - 1) * exp( - u2(f) + u1(f))
        )

phi_norm_approx = lambda f: (
        mu0 * pi * N0**2/l * (2 * r1**2 /denom_approx(f) - 2/k(f)**2 * enum_approx(f) / denom_approx(f) + (rsp**2 - r2**2))
        )

phi_norm_approx_abs = lambda f: abs(phi_norm_approx(f))
phi_norm_approx_arg = lambda f: arg(phi_norm_approx(f))


# ---------------------------------------------------------#
# Generate points for omega axis                           #
# ---------------------------------------------------------#
# See also separate file stuetzpunkte.py
n                = np.linspace(0,npts,npts)
expufunc         = np.frompyfunc(exp,1,1)
frequency_vector = expufunc((1-n/npts)*log(fmin)) * expufunc(n*log(fmax)/npts)


# ---------------------------------------------------------#
# Numerically evaluate functions                           #
# ---------------------------------------------------------#
phi_norm_abs_ufunc = np.frompyfunc(phi_norm_abs,1,1)
phi_norm_abs_num   = phi_norm_abs_ufunc(frequency_vector)
phi_norm_arg_ufunc = np.frompyfunc(phi_norm_arg,1,1)
phi_norm_arg_num   = phi_norm_arg_ufunc(frequency_vector)

phi_norm_approx_abs_ufunc = np.frompyfunc(phi_norm_approx_abs,1,1)
phi_norm_approx_abs_num   = phi_norm_approx_abs_ufunc(frequency_vector)
phi_norm_approx_arg_ufunc = np.frompyfunc(phi_norm_approx_arg,1,1)
phi_norm_approx_arg_num   = phi_norm_approx_arg_ufunc(frequency_vector)

# ---------------------------------------------------------#
# Unfortunately, the  arg() function only  delivers values #
# between -pi and  +pi for the angle of  a complex number, #
# which,  while  correct,  is   not  suitable  for  pretty #
# plotting, so we  will shift the values  larger then zero #
# accordingly for a continuous curve.                      #
# ---------------------------------------------------------#
phi_norm_arg_num   = np.unwrap(phi_norm_arg_num)
phi_norm_approx_arg_num = np.unwrap(phi_norm_approx_arg_num)


# ---------------------------------------------------------#
# Scale values for improved legibility in plot             #
# ---------------------------------------------------------#
phi_norm_abs_num   = phi_norm_abs_ufunc(frequency_vector)
phi_norm_abs_num   = 1e3 * phi_norm_abs_num
phi_norm_arg_num   = 180 / pi * phi_norm_arg_num

phi_norm_approx_arg_num = 180 / pi * phi_norm_approx_arg_num
phi_norm_approx_abs_num = phi_norm_approx_abs_ufunc(frequency_vector)
phi_norm_approx_abs_num = 1e3 * phi_norm_approx_abs_num


# ---------------------------------------------------------#
# Plot the Things                                          #
# ---------------------------------------------------------#
matplotlib.pyplot.rc('text', usetex=True)
matplotlib.pyplot.rc('font', family='serif')

fig   = figure(1)
axes1 = fig.add_subplot(211)
axes1.plot(frequency_vector,phi_norm_abs_num,color=plot_color_fit,label=plot_label_fit)
axes1.plot(frequency_vector,phi_norm_approx_abs_num,color=plot_color_fit_approx,label=plot_label_fit_approx)
axes1.set_xlim([fmin*0.9,fmax*1.1])
axes1.set_xscale(plot_scale_x)
axes1.set_xlabel(plot_label_x,fontdict=font)
axes1.set_ylabel(plot_1_label_y,fontdict=font)
axes1.set_title(plot_1_title,fontdict=titlefont)
axes1.legend(fontsize=plot_legend_fontsize)
axes1.tick_params(labelsize=9)

axes2 = fig.add_subplot(212)
axes2.plot(frequency_vector,phi_norm_arg_num,color=plot_color_fit,label=plot_label_fit)
axes2.plot(frequency_vector,phi_norm_arg_num,color=plot_color_fit_approx,label=plot_label_fit_approx)
axes2.set_xlim([fmin*0.9,fmax*1.1])
axes2.set_xscale(plot_scale_x)
axes2.set_xlabel(plot_label_x,fontdict=font)
axes2.set_ylabel(plot_2_label_y,fontdict=font)
axes2.set_title(plot_2_title,fontdict=titlefont)
axes2.legend(fontsize=plot_legend_fontsize)
axes2.tick_params(labelsize=9)

fig.subplots_adjust(bottom=0.15,left=0.125,right=0.925,top=0.95,hspace=0.5)

fig.savefig('plots-pgf/hollow--st--freq--phi-norm.pgf')
fig.savefig('plots-pdf/hollow--st--freq--phi-norm.pdf')


# ---------------------------------------------------------#
# Save listing to file                                     #
# ---------------------------------------------------------#
dumpfile = open('listings/hollow--st--phi-norm.tex', 'w')

table_opening = r"""
{%
    \begin{center}
    \captionof{table}{%
        Parameterwerte f\"ur  Fit-Funktion in  Abbildung  \ref{fig:st:freq:phi}
    }
    \label{tab:fitparams:st:phi}
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
