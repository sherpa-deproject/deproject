"""
Plot density and temperature results for M87.  Compare values from
deproject to those obtained by independent XSPEC analysis.
"""

import numpy
arcsec_per_pix = 0.492

phys = read_file('m87_phys.dat')        # values from XSPEC analysis

r0 = get_colvals(phys, 'r0')            # radius in ACIS pixels
r1 = get_colvals(phys, 'r1')
ne = get_colvals(phys, 'ne')
nelow = get_colvals(phys, 'nelow')
nehigh = get_colvals(phys, 'nehigh')
kt = get_colvals(phys, 'kt')
ktlow = get_colvals(phys, 'ktlow')
kthigh = get_colvals(phys, 'kthigh')

def add_set_clear_window(window_id):
    """Add a new window with given name ``window_id``, set it to be the current
    window, and clear it.  The exception handling is to allow running this
    function repeatedly within an interactive session.
    """
    try:
        add_window(['id', window_id])
    except RuntimeError, e:
        set_current_window(window_id)
    try:
        clear_plot()
    except:
        pass

##### Density ###########

add_set_clear_window('density')
x = (r0 + r1)/2. * arcsec_per_pix / 60. # arcmin
x_err = (r1 - r0)/2. * arcsec_per_pix / 60. # arcmin
y = ne

bad = nelow < 1e-4
y_err_low = ne - nelow
y_err_high = nehigh - ne
y_err_low[bad] = 0.01
y_err_high[bad] = 0.01

add_curve(x, y, [y_err_low, y_err_high, x_err, x_err])
log_scale()

print "Calculating density for z=%.5f angdist=%.2e cm" % (dep.redshift, dep.angdist)
d_ne = dep.get_density()
add_curve(x, d_ne)
set_curve(['symbol.color', 'red', 'line.color', 'red'])
set_plot_xlabel('Radial distance (arcmin)')
set_plot_ylabel('Density (cm^{-3})')
limits(X_AXIS, 0.2, 10)
# print_window('m87_density', ['format', 'png'])

##### Temperature ###########

add_set_clear_window('temperature')
y = kt

bad = nelow < 1e-4
y_err_low = kt - ktlow
y_err_high = kthigh - kt
y_err_low[bad] = 0.01
y_err_high[bad] = 0.01

add_curve(x, y, [y_err_low, y_err_high, x_err, x_err])
log_scale(X_AXIS)

d_kt = dep.get_par('xsmekal.kt')
add_curve(x, d_kt)
set_curve(['symbol.color', 'red', 'line.color', 'red'])
set_plot_xlabel('Radial distance (arcmin)')
set_plot_ylabel('Temperature (keV)')
limits(X_AXIS, 0.2, 10)
# print_window('m87_temperature', ['format', 'png'])

