#!/usr/bin/env python3

from sympy import *
from mpmath import *
from matplotlib.pyplot import *
#init_printing()     # make things prettier when we print stuff for debugging.


# ************************************************************************** #
# Self-Inductance L of copper coil with hollow copper cylinder inserted.     #
# ************************************************************************** #

# All values are in standard SI units unless otherwise noted.

# -------------------------------------------------------- #
# Default precision is insufficient, therefore we increase #
# precision.   One  can  increase the  number  of  decimal #
# places or bits, where the number of bits places is ~3.33 #
# times the number of decimal places.                      #
# -------------------------------------------------------- #
#mp.dps=25  # decimal places
mp.prec=80 # precision in bits


# ---------------------------------------------------------#
# Init, Define Variables and Constants                     #
# ---------------------------------------------------------#
mu0   = 4*pi*1e-7                                        # vacuum permeability
sigma = 53e6                            # de.wikipedia.org/wiki/Kupfer: 58.1e6
dsp   = 98e-3                                               # diameter of coil
rsp   = dsp / 2                                               # radius of coil
r1    = 30e-3                                # inner radius of copper cylinder
r2    = 35e-3                                # outer radius of copper cylinder
N0    = 574
l     = 500e-3                                         # length of copper coil
npts  = 1e3                                  # number of points for plot curve
fmin  = 1
fmax  = 2500
    # -----------------------------------------------------#
    # Create  a list  for convenient  printing of  vars to #
    # file, add LaTeX where necessary.                     #
    # -----------------------------------------------------#
params = [
        '        ' + '$\mu_0'    + '$ & $' +  '\SI{'   + str(mu0)    + r'}{\newton\per\ampere\squared}' + r'$\\' + "\n",
        '        ' + '$\sigma'   + '$ & $' +  '\SI{'   + str(sigma)  + r'}{\ampere\per\volt\per\meter}' + r'$\\' + "\n",
        '        ' + '$d_{Sp}'   + '$ & $' +  '\SI{'   + str(dsp)    + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$r_{Sp}'   + '$ & $' +  '\SI{'   + str(rsp)    + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$r_1'      + '$ & $' +  '\SI{'   + str(r1)     + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$r_2'      + '$ & $' +  '\SI{'   + str(r2)     + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$l'        + '$ & $' +  '\SI{'   + str(l)      + r'}{\meter}'                     + r'$\\' + "\n",
        '        ' + '$NPTS'     + '$ & $' +  r'\num{' + str(npts)   + '}'                              + r'$\\' + "\n",
        '        ' + '$N_0'      + '$ & $' +  r'\num{' + str(N0)     + '}'                              + r'$\\' + "\n",
        '        ' + '$f_{min}'  + '$ & $' +  '\SI{'   + str(fmin)   + r'}{\hertz}'                     + r'$\\' + "\n",
        '        ' + '$f_{max}'  + '$ & $' +  '\SI{'   + str(fmax)   + r'}{\hertz}'                     + r'$\\' + "\n",
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
plot_label_y            = 'Selbstinduktion L (mH)'
plot_title              = "Selbstinduktionskoeffizient, Spule mit Kupferrohr"


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

# Exact solution:
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

L = lambda f: re(phi_norm(f))

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

L_approx = lambda f: re(phi_norm_approx(f))


# ---------------------------------------------------------#
# Generate points for frequency axis                       #
# ---------------------------------------------------------#
n                = np.linspace(1,npts,npts)
expufunc         = np.frompyfunc(exp,1,1)
frequency_vector = fmin*expufunc(n*log(fmax-fmin)/npts)


# ---------------------------------------------------------#
# Numerically evaluate function                            #
# ---------------------------------------------------------#
L_ufunc = np.frompyfunc(L,1,1)
L_num   = L_ufunc(frequency_vector)
L_num = 1e3 * L_num                     # improve legibility

L_approx_ufunc = np.frompyfunc(L_approx,1,1)
L_approx_num   = L_approx_ufunc(frequency_vector)
L_approx_num = 1e3 * L_approx_num


# ---------------------------------------------------------#
# Plot the Things                                          #
# ---------------------------------------------------------#
matplotlib.pyplot.rc('text', usetex=True)
matplotlib.pyplot.rc('font', family='serif')

figwidth = 8.27 # in inches
fig  = figure(1,figsize=(figwidth,figwidth*0.4))
axes = fig.add_subplot(111)
axes.plot(frequency_vector,L_num,color=plot_color_fit,label=plot_label_fit)
axes.plot(frequency_vector,L_approx_num,color=plot_color_fit_approx,label=plot_label_fit_approx)
axes.set_xscale(plot_scale_x)
axes.set_xlim([fmin*0.9,fmax*1.1])
axes.set_xlabel(plot_label_x,fontdict=font)
axes.set_ylabel(plot_label_y,fontdict=font)
axes.set_title(plot_title,fontdict=titlefont)
axes.legend(fontsize=plot_legend_fontsize)
axes.tick_params(labelsize=9)

fig.subplots_adjust(bottom=0.15,left=0.125,right=0.925,top=0.90)

fig.savefig('plots-pgf/hollow--cu--L.pgf')
fig.savefig('plots-pdf/hollow--cu--L.pdf')


# ---------------------------------------------------------#
# Save listing to file                                     #
# ---------------------------------------------------------#
dumpfile = open('listings/hollow--cu--L.tex', 'w')

table_opening = r"""
{%
    \begin{center}
    \captionof{table}{%
        Parameterwerte  f\"ur  Fit-Funktion  in  Abbildung  \ref{fig:cu:freq:L}
    }
    \label{tab:fitparams:cu:L}
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
