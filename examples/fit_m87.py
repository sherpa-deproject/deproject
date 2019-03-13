"""
Fit the M87 data and create basic plots.

This assumes that Sherpa is using the Matplotlib backend; that is,
the plot_pkg line in ~/.sherpa.rc is set to pylab and not chips.

"""

import numpy
from deproject.deproject import Deproject
from sherpa.astro.ui import plot_fit, plot_fit_delchi, \
    set_method, set_stat, set_xlog
from astropy import units as u

from matplotlib import pyplot as plt

redshift = 0.004233                     # M87 redshift
angdist = 16 * u.Mpc                    # M87 distance
theta = 75 * u.deg                      # Covering angle of sectors

# Each pixel is 0.492 arcseconds
radii = numpy.arange(30., 640., 30) * 0.492 * u.arcsec
dep = Deproject(radii, theta=theta, angdist=angdist)

# Load datasets for each annulus
for annulus in range(len(radii) - 1):
    dep.load_pha('m87/r%dgrspec.pha' % (annulus + 1), annulus)

# Subtract background
dep.subtract()

# Set source model and ignore specified energy ranges
dep.set_source('xswabs*xsmekal')
dep.ignore(None, 0.5)
dep.ignore(1.8, 2.2)
dep.ignore(7, None)

# Specify Galactic absorption
dep.set_par('xswabs.nh', 0.0255)
dep.freeze('xswabs.nh')

# Initialize abundance to 0.5 and thaw
dep.set_par('xsmekal.abundanc', 0.5)
dep.thaw('xsmekal.abundanc')

# Set redshift
dep.set_par('xsmekal.redshift', redshift)

# Do the initial onion-peeling fit with
#   Levenberg-Marquardt optimization method
#   XSPEC variance with a Chi^2 fit statistic
#
set_method("levmar")
set_stat("chi2xspecvar")

dep.guess()

# Display the central annulus (using Sherpa routines)
set_xlog()
plot_fit(0)
plt.savefig('m87_ann0_guess.png')

onion = dep.fit()

# Display the central annulus (using Sherpa routines)
plot_fit_delchi(0)
plt.savefig('m87_ann0_fit.png')

# Create basic plots
dep.density_plot()
plt.savefig('m87_density.png')

dep.par_plot('xsmekal.kt')
plt.savefig('m87_temperature.png')

dep.fit_plot('rstat')
plt.savefig('m87_rstat.png')

# Error analysis
errs = dep.covar()

plt.subplot(2, 1, 1)
dep.covar_plot('xsmekal.kt', clearwindow=False)
plt.xlabel('')
plt.subplot(2, 1, 2)
dep.covar_plot('xsmekal.abundanc', clearwindow=False)
plt.subplots_adjust(hspace=0.3)
plt.savefig('m87_temperature_abundance.png')

dep.covar_plot('density', units='arcmin')
plt.savefig('m87_density_errs.png')

# Display some basic information on the configuration
#
print("z={:.5f} angdist={}".format(dep.redshift, dep.angdist))
