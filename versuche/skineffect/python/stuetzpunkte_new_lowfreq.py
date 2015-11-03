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
# Define Variables and Constants                           #
# ---------------------------------------------------------#
npts              = 19 # careful: number of points is npts + 1 (starts at 0)
#fmin              = 5e-4
#fmax              = 5e-2
fmin              = 0.1
fmax              = 0.2
highest_frequency = fmin * exp(log(fmax-fmin))
freq_max_ratio    = highest_frequency / fmax
font = {
        'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 11,
        }
plot_legend_fontsize       = 11
plot_color_old             = 'magenta'
plot_color_new             = 'blue'
plot_color_common          = 'black'
plot_label_points_old      = r"St\"utzpunktformel A: $\displaystyle f_{kA} = f_{min} \cdot exp\Biggr(\frac{k}{NPTS} \cdot ln(f_{max}-f_{min})\Biggr)$"
plot_label_points_new      = r"St\"utzpunktformel B: $\displaystyle f_{kB} = exp\Biggr((1-\frac{k}{NPTS}) \cdot ln(f_{min})\Biggr) \cdot exp\Biggr(\frac{k}{NPTS} \cdot ln(f_{max})\Biggr)$"
plot_label_vertical_common = r"minimale Frequenz St\"utzpunkt: "
plot_label_vertical_old    = r"maximale Frequenz St\"utzpunkt, Methode A: "
plot_label_vertical_new    = r"maximale Frequenz St\"utzpunkt, Methode B: "
plot_added_text            = r"Verh\"altnis der maximalen Frequenzen Methode A -- Methode B: $\displaystyle \frac{f_{kA}}{f_{kB}}\Bigg|_{k=NPTS} \approx " + str(freq_max_ratio) + "$"
plot_size_measurements     = 24
plot_scale_x               = 'log'
plot_label_x               = r"Frequenz des St\"utzpunkts (Hz)"
plot_1_label_y             = 'k (siehe Formel)'
plot_1_title               = r"Vergleich St\"utzpunktformeln f\"ur den Frequenzbereich von " + str(fmin) + " Hz bis " + str(highest_frequency) + " Hz, " + str(npts+1) + " Punkte"
y_lim_low                  = -2
y_lim_high                 = npts + 2
x_lim_low                  = 0.67 * fmin
x_lim_high                 = 1.33 * fmin * fmax


# ---------------------------------------------------------#
# Generate points for frequency axis                       #
# ---------------------------------------------------------#
n                    = np.linspace(0,npts,npts)
expufunc             = np.frompyfunc(exp,1,1)
frequency_vector_old = fmin*expufunc(n*log(fmax-fmin)/npts)
frequency_vector_new = expufunc((1-n/npts)*log(fmin)) * expufunc(n*log(fmax)/npts)

plot_label_vertical_common += str(frequency_vector_old[0]) + " Hz"
plot_label_vertical_old    += str(frequency_vector_old[npts-1]) + " Hz"
plot_label_vertical_new    += str(frequency_vector_new[npts-1]) + " Hz"


# ---------------------------------------------------------#
# Plot the Things                                          #
# ---------------------------------------------------------#
matplotlib.pyplot.rc('text', usetex=True)
matplotlib.pyplot.rc('font', family='serif')

#fig1 = figure(1)
#fig1 = figure(1,figsize=(9,9))
#fig1 = figure(1,figsize=(8.26,11.7)) # A4 size in inches
fig1 = figure(1,figsize=(8.26,10.0))
axes1 = fig1.add_subplot(111)
axes1.set_position([0.1,0.1,0.5,0.8])
axes1.scatter(frequency_vector_old,
        n,
        color=plot_color_old,
        s=plot_size_measurements,
        label=plot_label_points_old
        )
axes1.scatter(frequency_vector_new,
        n,
        color=plot_color_new,
        s=plot_size_measurements,
        label=plot_label_points_new
        )
# Draw the common starting point black and a bit bigger
axes1.scatter(frequency_vector_old[0],
        n[0],
        color=plot_color_common,
        s=plot_size_measurements*1.5,
        )
axes1.plot([frequency_vector_old[0],frequency_vector_old[0]],
        [y_lim_low,y_lim_high],
        color=plot_color_common,
        label=plot_label_vertical_common
        )
axes1.plot([frequency_vector_old[npts-1],frequency_vector_old[npts-1]],
        [y_lim_low,y_lim_high],
        color=plot_color_old,
        label=plot_label_vertical_old
        )
axes1.plot([frequency_vector_new[npts-1],frequency_vector_new[npts-1]],
        [y_lim_low,y_lim_high],
        color=plot_color_new,
        label=plot_label_vertical_new
        )
axes1.set_xscale(plot_scale_x)
axes1.set_ylim([y_lim_low,y_lim_high])
#axes1.set_xlim([x_lim_low,x_lim_high])
axes1.set_xlabel(plot_label_x,fontdict=font)
axes1.set_ylabel(plot_1_label_y,fontdict=font)
axes1.set_title(plot_1_title,fontdict=font)

    # ---------------------------------------------------- #
    # Work some  magic to append  the fraction of  the two #
    # methods  to  the legend  instead  of  it being  some #
    # random piece of text on the plot.                    #
    # ---------------------------------------------------- #
rect = matplotlib.patches.Rectangle([0,0],0,0,color='white',label=plot_added_text)
handles,legends = axes1.get_legend_handles_labels()
handles.append(rect)
axes1.legend(handles=handles,fontsize=plot_legend_fontsize,loc='upper left',bbox_to_anchor=(0.0,-0.075))

    # ---------------------------------------------------- #
    # This  would be  necessary if  we wanted  to actually #
    # draw  the  patch,  leaving   this  here  for  future #
    # reference.                                           #
    # ---------------------------------------------------- #
#axes1.add_patch(rect)

fig1.subplots_adjust(bottom=0.3,left=0.1,right=0.9,top=0.95,hspace=0.5)

fig1.savefig('plots-pgf/stuetzpunkte-lowfreq.pgf')
fig1.savefig('plots-pdf/stuetzpunkte-lowfreq.pdf')

#show()
