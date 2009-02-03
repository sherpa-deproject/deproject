from cosmocalc import cosmocalc, arcsec_per_rad
from math import pi, sqrt

redshift = 1.06

set_method("neldermead")
set_method("levmar")
set_stat("cstat")

create_model_component('xsphabs', 'xsphabs1')
create_model_component('xsapec', 'xsapec1')

obsids = (9407, 9774, 9775, 9408)
ann = 0
radii = ('2.5', '6', '17')
rvals = [float(x) for x in radii]
for dataid, obsid in enumerate(obsids):
    load_pha(dataid+1, '3c186/%d/ellipse%s-%s.pi' % (obsid, radii[ann], radii[ann+1]))
    set_source(dataid+1, 'xsphabs1*xsapec1')

execfile("3c186/acis-s-bkg.py")
acis_s_bkg = get_bkg_source()
for dataid, obsid in enumerate(obsids):
    bkg_norm_name = 'bkg_norm_%d' % obsid
    create_model_component('xsconstant', bkg_norm_name)
    bkg_norm = eval(bkg_norm_name)
    set_bkg_model(dataid+1, bkg_norm * acis_s_bkg)

ignore(None, 0.5)
ignore(7, None)
freeze("xsphabs1.nh")

set_par('xsapec1.redshift', redshift)
set_par('xsphabs1.nh', 0.0564)

fit()

DA_cm = cosmocalc(redshift)['DA_cm']
r0 = rvals[ann] / arcsec_per_rad * DA_cm
r1 = rvals[ann+1] / arcsec_per_rad * DA_cm
volume_shell = 4*pi/3 * (r1**3 - r0**3)
volume_ann_shell = 4 * pi / 3 * (r1**2 - r0**2)**1.5

volume = volume_ann_shell

z = redshift
norm = xsapec1.norm.val

density = sqrt( norm * 1e14 * 4 * pi * (DA_cm * (1+z))**2 / (1.18 * volume) )

# norm =  1e-14  / (4 pi (D_A*(1+z))^2) n_e * n_p * volume
# norm =  1e-14  / (4 pi (D_A*(1+z))**2) n_p * 1.18 * n_p * volume
# n_p**2 = norm * 1e14 * 4 pi (D_A*(1+z))**2 / (1.18 * volume)
# n_p = sqrt( norm * 1e14 * 4 * pi * (D_A * (1+z))**2 / (1.18 * volume) )
