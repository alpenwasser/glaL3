# ************************************************************************** #
# EXAMPLE FILE FOR APPENDIX OF LATEX. THIS FILE IS NOT FUNCIONAL!            #
# ************************************************************************** #

# ************************************************************************** #
# Magnetic field inside copper coil with massive alu cylinder                #
# ************************************************************************** #

# All values are in standard SI units unless otherwise noted.

# ---------------------------------------------------------#
# Define Variables and Constants                           #
# ---------------------------------------------------------#
mu0            = 4*pi*1e-7
rho_kuchling   = 0.027e-6
sigma_kuchling = 1/rho_kuchling
sigma_abs      = 24e6
sigma_arg      = 22.25e6
r              = 0
r0             = 45e-3
B0             = 6.9e-2
npts           = 1e3
fmin           = 1
fmax           = 250

#... LaTeX export stuff ...#

#... Plot Configuration ...#

# ---------------------------------------------------------#
# Functions                                                #
#                                                          #
# NOTE: We use  frequency f  instead of  angular frequency #
# omega since that is what we actually set on the function #
# generator.                                               #
# NOTE: We evaluate B_abs and B_arg based on two different #
# values for sigma, which allows to fit each of the curves #
# more accurately.                                         #
# ---------------------------------------------------------#

k_abs = lambda f: sqrt((2*np.pi*f*mu0*sigma_abs)/2)*(mpc(1,-1))
k_arg = lambda f: sqrt((2*np.pi*f*mu0*sigma_arg)/2)*(mpc(1,-1))

enum_abs  = lambda f: besselj(0,k_abs(f)*r)
denom_abs = lambda f: besselj(0,k_abs(f)*r0)
enum_arg  = lambda f: besselj(0,k_arg(f)*r)
denom_arg = lambda f: besselj(0,k_arg(f)*r0)

B_abs = lambda f: abs(enum_abs(f) / denom_abs(f) * B0)
B_arg = lambda f: arg(enum_arg(f) / denom_arg(f) * B0)

# ---------------------------------------------------------#
# Generate points for frequency axis                       #
# ---------------------------------------------------------#
n                = np.linspace(0,npts,npts)
expufunc         = np.frompyfunc(exp,1,1)
frequency_vector = fmin*expufunc(n*log(fmax-fmin)/npts)

# ---------------------------------------------------------#
# Numerically evaluate functions                           #
# ---------------------------------------------------------#
Babsufunc = np.frompyfunc(B_abs,1,1)
B_abs_num = Babsufunc(frequency_vector)
Bargufunc = np.frompyfunc(B_arg,1,1)
B_arg_num = Bargufunc(frequency_vector)

# ---------------------------------------------------------#
# The arg() function only gives values between -pi and +pi #
# for the angle of a complex number.                       #
# ---------------------------------------------------------#
B_arg_num = np.unwrap(B_arg_num)

# ---------------------------------------------------------#
# Measurement Values from the actual experiment            #
# ---------------------------------------------------------#
frequencies_measured = np.array([     1,     5,    #... etc ...# ])
phases_degrees       = np.array([   5.4,    26,    #... etc ...# ])
voltages             = np.array([6.9e-2,6.5e-2,    #... etc ...# ])

# ---------------------------------------------------------#
# Plot the Things                                          #
# ---------------------------------------------------------#
fig   = figure(1)
axes1 = fig.add_subplot(211)
axes1.plot(frequency_vector,B_abs_num,color=plot_color_fit,label=plot_label_fit)
axes1.scatter(frequencies_measured,
        voltages,
        color=plot_color_measurements,
        s=plot_size_measurements,
        label=plot_label_measurements
        )
axes1.set_xlim([fmin*0.9,fmax*1.1])
axes1.set_xscale(plot_scale_x)
axes1.set_xlabel(plot_label_x,fontdict=font)
axes1.set_ylabel(plot_1_label_y,fontdict=font)
axes1.set_title(plot_1_title,fontdict=titlefont)
axes1.legend(fontsize=plot_legend_fontsize)
axes1.tick_params(labelsize=9)

#... second plot analogous ...#

fig.savefig('plots-pgf/massive--alu--freq.pgf')
fig.savefig('plots-pdf/massive--alu--freq.pdf')

#... save listing to LaTeX file for automated inclusion ...#

#... save results to txt file for error analysis ... #
