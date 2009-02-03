from deproject import Deproject
from cosmocalc import cosmocalc

redshift = 0.004233
from math import pi

dep = Deproject(numpy.arange(30., 640., 30)*0.492)
dep.theta = 75
dep.angdist = cosmocalc(redshift)['DA_cm'] * 0.892

set_method("levmar")
set_stat("chi2gehrels")

for ann in range(dep.nshell):
    dep.load_pha('m87/r%dgrspec.pha' % (ann+1), annulus=ann)

dep.set_source('xswabs*xsmekal')
dep.ignore(None, 0.5)
dep.ignore(1.8, 2.2)
dep.ignore(7, None)

dep.set_par('xswabs.nh', 0.0255)
dep.freeze("xswabs.nh")

dep.set_par('xsmekal.abundanc', 0.5)
dep.thaw('xsmekal.abundanc')

dep.set_par('xsmekal.redshift', redshift)
dep.subtract()

dep.fit()

