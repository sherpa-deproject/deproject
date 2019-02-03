"""
Plot density and temperature results for M87.  Compare values from
deproject to those obtained by independent XSPEC analysis.
"""

import numpy
arcsec_per_pix = 0.492
arcsec_per_rad = (180.*3600.)/numpy.pi

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

##### Temperature ###########

add_set_clear_window('temperature')
x = (r0 + r1)/2. * arcsec_per_pix / 60. # arcmin
x_err = (r1 - r0)/2. * arcsec_per_pix / 60. # arcmin
y = kt

bad = nelow < 1e-4
y_err_low = kt - ktlow
y_err_high = kthigh - kt
y_err_low[bad] = 0.01
y_err_high[bad] = 0.01

add_curve(x, y, [y_err_low, y_err_high, x_err, x_err])
set_curve(['symbol.color','green','line.color','green','err.color','green','err.thickness',3])
log_scale(X_AXIS)

d_kt = dep.get_par('xsmekal.kt')
d_kt_m = []
d_kt_p = []
for n in range(len(radii)-1):
    d_kt_m.append(-conf_params[n].parmins[0])
    d_kt_p.append(conf_params[n].parmaxes[0])

add_curve(x, d_kt)
set_curve(['symbol.color', 'red', 'line.color', 'red'])
set_plot_xlabel('Radial distance (arcmin)')
set_plot_ylabel('Temperature (keV)')
add_curve(x, d_kt, [d_kt_m, d_kt_p, x_err, x_err])
set_curve(['symbol.color','red','line.color','red','err.color','red','err.thickness',2])
limits(X_AXIS, 0.2, 6)
# print_window('m87_temperature', ['format', 'png'])

##### Density ###########

r_sphere = dep.radii[-1] / arcsec_per_rad * angdist
volume = 4 * numpy.pi / 3 * r_sphere**3  # volume of sphere enclosing outer shell (cm^3)
mue=1.18

print "Calculating density for z=%.5f angdist=%.2e cm" % (dep.redshift, dep.angdist)

norm = dep.get_par('xsmekal.norm')   # returns array of kT values
norm_m = []
norm_p = []
for n in range(len(radii)-1):
    norm_m.append(conf_params[n].parmins[2])
    norm_p.append(conf_params[n].parmaxes[2])

dens = []
dens_low = []
dens_high = []
dens_const =  numpy.sqrt(4 * numpy.pi * angdist**2 * 1e14 * (1.0+redshift)**2 / volume * mue);
for shell in range(dep.nshell):
    dens.append(  dens_const*numpy.sqrt(norm[shell]))
    dens_low.append(dens_const*numpy.sqrt(norm[shell]+norm_m[shell]))
    dens_high.append(dens_const*numpy.sqrt(norm[shell]+norm_p[shell]))

dens_loerr = numpy.array(dens)-numpy.array(dens_low)
dens_hierr = numpy.array(dens_high)-numpy.array(dens)

add_set_clear_window('density')
y = ne

bad = nelow < 1e-4
y_err_low = ne - nelow
y_err_high = nehigh - ne
y_err_low[bad] = 0.01
y_err_high[bad] = 0.01

add_curve(x, y, [y_err_low, y_err_high, x_err, x_err])
set_curve(['symbol.color','green','line.color','green','err.color','green','err.thickness',3])
log_scale()

add_curve(x, dens)
set_curve(['symbol.color', 'red', 'line.color', 'red'])
set_plot_xlabel('Radial distance (arcmin)')
set_plot_ylabel('Density (cm^{-3})')
add_curve(x, dens, [dens_loerr, dens_hierr, x_err, x_err])
set_curve(['symbol.color','red','line.color','red','err.color','red','err.thickness',2])
limits(X_AXIS, 0.2, 6)
# print_window('m87_density', ['format', 'png'])
