import deproject

radii = ('2.5', '6', '17')
dep = deproject.Deproject(radii=[float(x) for x in radii])

set_method("neldermead")
set_method("levmar")
set_stat("cstat")

obsids = (9407, 9774, 9775, 9408)
for ann in range(dep.nshell):
    for obsid in obsids:
        dep.load_pha('3c186/%d/ellipse%s-%s.pi' % (obsid, radii[ann], radii[ann+1]),
                     annulus=ann)

dep.set_source('xsphabs*xsapec')
dep.ignore(None, 0.5)
dep.ignore(7, None)
dep.freeze("xsphabs.nh")

dep.set_par('xsapec.redshift', 1.06)
dep.set_par('xsphabs.nh', 0.0564)

execfile("3c186/acis-s-bkg.py")
acis_s_bkg = get_bkg_source()

dep.set_bkg_model(acis_s_bkg)

dep.pre_fit()
