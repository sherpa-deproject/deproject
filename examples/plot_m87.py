"""
Plot density and temperature results for M87.  Compare values from
deproject to those obtained by independent XSPEC analysis.

The fit_m87.py script must have already been run, to set up the
``dep`` object with the deprojection results.

What are the errors in Paul's file (one sigma or 90%)?

"""

import numpy

try:
    from pycrates import read_file
except ImportError:
    try:
        from astropy.table import Table
    except ImportError:
        raise ImportError("Unable to load pycrates or astropy")

from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches


# Draw a rectangle for each error box, following
# https://matplotlib.org/2.2.3/gallery/statistics/errorbars_and_boxes.html#sphx-glr-gallery-statistics-errorbars-and-boxes-py
#
def make_error_boxes(ax, xdata, ydata, xerror, yerror, facecolor='gray',
                     edgecolor='None', alpha=0.5):

    # Create list for all the error patches
    errorboxes = []

    # Loop over data points; create box from errors at each point
    for x, y, xe, ye in zip(xdata, ydata, xerror.T, yerror.T):
        rect = Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum())
        errorboxes.append(rect)

    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=facecolor, alpha=alpha,
                         edgecolor=edgecolor)

    # Add collection to axes
    ax.add_collection(pc)

    """
    # Plot errorbars
    artists = ax.errorbar(xdata, ydata, xerr=xerror, yerr=yerror,
                          fmt='None', ecolor='k')

    return artists
    """

    # Create a fake artist for the legend following
    # https://matplotlib.org/users/legend_guide.html#creating-artists-specifically-for-adding-to-the-legend-aka-proxy-artists
    #
    return mpatches.Patch(facecolor=facecolor, edgecolor=edgecolor,
                          alpha=alpha)


# Read in results from XSPEC analysis of M87

arcsec_per_pix = 0.492
arcsec_per_rad = (180.*3600.)/numpy.pi

# Read in the values from XSPEC analysis, where the radius is in
# ACIS pixels.
#
try:
    phys = read_file('m87_phys.dat')

    r0 = phys.r0.values
    r1 = phys.r1.values
    ne = phys.ne.values
    nelow = phys.nelow.values
    nehigh = phys.nehigh.values
    kt = phys.kt.values
    ktlow = phys.ktlow.values
    kthigh = phys.kthigh.values

except NameError:

    phys = Table.read('m87_phys.dat', format='ascii')

    r0 = phys['r0']
    r1 = phys['r1']
    ne = phys['ne']
    nelow = phys['nelow']
    nehigh = phys['nehigh']
    kt = phys['kt']
    ktlow = phys['ktlow']
    kthigh = phys['kthigh']

##### Density ###########

x = (r0 + r1) / 2. * arcsec_per_pix / 60.  # arcmin
x_err = (r1 - r0) / 2. * arcsec_per_pix / 60.  # arcmin
dx = numpy.vstack((x_err, x_err))

y = ne

bad = nelow < 1e-4
y_err_low = kt - ktlow
y_err_high = kthigh - kt
y_err_low[bad] = 0.01
y_err_high[bad] = 0.01

dy = numpy.vstack((y_err_low, y_err_high))

# Plot up the Sherpa results
dep.covar_plot('density', ylog=True, units='arcmin')

# Access a Matplotlib artist from the plot for the legend
ax = plt.gca()
sherpa_artist = ax.lines[1]

# Plot up the XSPEC results
#
xspec_artist = make_error_boxes(ax, x, y, dx, dy)

plt.legend((sherpa_artist, xspec_artist), ('Sherpa', 'XSPEC'))
plt.savefig('m87_compare_density.png')


print "Calculating density for z=%.5f angdist=%.2e cm" % (dep.redshift, dep.angdist)

y = kt

bad = nelow < 1e-4
y_err_low = ne - nelow
y_err_high = nehigh - ne
y_err_low[bad] = 0.01
y_err_high[bad] = 0.01

dy = numpy.vstack((y_err_low, y_err_high))

dep.covar_plot('xsmekal.kt', units='arcmin')

ax = plt.gca()
sherpa_axis = ax.lines[1]

xspec_artist = make_error_boxes(ax, x, y, dx, dy)

plt.legend((sherpa_artist, xspec_artist), ('Sherpa', 'XSPEC'),
           loc='upper left')

plt.savefig('m87_compare_temperature.png')
