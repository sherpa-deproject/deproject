from deproject import Deproject

redshift = 0.004233                     # M87 redshift
arcsec_per_pixel = 0.492                # ACIS plate scale
angdist = 4.9e25                        # M87 distance (cm) (16 Mpc)

radii = numpy.arange(30., 640., 30) * arcsec_per_pixel
dep = Deproject(radii, theta=75, angdist=angdist)

# Load datasets for each annulus
for annulus in range(len(radii)-1):
    dep.load_pha('m87/r%dgrspec.pha' % (annulus+1))

# Subtract background
dep.subtract()

# Set source model and ignore specified energy ranges
dep.set_source('xswabs*xsmekal')
dep.ignore(None, 0.5)
dep.ignore(1.8, 2.2)
dep.ignore(7, None)

# Specify Galactic absorption
dep.set_par('xswabs.nh', 0.0255)
dep.freeze("xswabs.nh")

# Initialize abundance to 0.5 and thaw
dep.set_par('xsmekal.abundanc', 0.5)
dep.thaw('xsmekal.abundanc')

# Set redshift
dep.set_par('xsmekal.redshift', redshift)

# Do the initial onion-peeling fit
set_method("levmar")                    # Levenberg-Marquardt optimization method
set_stat("chi2gehrels")                 # Gehrels Chi^2 fit statistic
dep.fit()

