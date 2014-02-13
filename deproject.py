"""
Deproject from a set of 2-d annular spectra to the 3-d object properties.

:Copyright: Smithsonian Astrophysical Observatory (2009)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""
import re
from math import pi, sqrt
import numpy
import specstack
from cosmocalc import cosmocalc, arcsec_per_rad
import sherpa.astro.ui as SherpaUI

class Deproject(specstack.SpecStack):
    """
    Deproject from a set of 2-d annular spectra to the 3-d object properties.

    The ``radii`` parameter must be a list of values that starts with the inner
    radius of the inner annulus and includes each radius up through the outer
    radius of the outer annulus.  Thus the ``radii`` list will be one element
    longer than the number of annuli.

    :param radii: sorted list of circular annulus radii (arcsec) for extracted spectra
    :param theta: azimuthal extent of annuli (degrees) (default = 360)
    :param angdist: angular size distance (cm)
    """
    def __init__(self, radii, theta=360, angdist=None):
        if len(radii) < 2:
            raise ValueError('radii parameter must be a list with at least two values')
        self.radii = radii
        self.nshell = len(radii)-1
        self._theta = theta
        self._angdist = angdist
        self._redshift = None
        super(Deproject, self).__init__()

    def _get_redshift(self):
        if self._redshift is None:
            self._redshift = self.find_parval('redshift')
        return self._redshift

    def _set_redshift(self, redshift):
        self._redshift = redshift

    def _get_angdist(self):
        if self._angdist is None:
            cc = cosmocalc(self.redshift)
            self._angdist = cc['DA_cm']
        return self._angdist

    def _set_angdist(self, angdist):
        self._angdist = angdist

    redshift = property(_get_redshift, _set_redshift, None, "Source redshift")
    angdist = property(_get_angdist, _set_angdist, None, "Angular size distance")

    def _calc_vol_norm(self):
        """
        Calculate the normalized volumes of cylindrical annuli intersecting with
        spherical shells.  Sets ``self.vol_norm`` to volume[i,j] / v_sphere
        giving the normalized volume of i'th shell intersecting with j'th
        annulus (where indexes start at 0).  V_sphere = 4/3 pi r^3 where r is
        the outer radius of the outer shell.

        :rtype: None
        """
        r = self.radii
        theta_rad = self._theta / 180. * pi
        cv = numpy.zeros([self.nshell, self.nshell])
        v = numpy.zeros([self.nshell, self.nshell])
        for a, ra0 in enumerate(r[:-1]):  # Annulus
            ra1 = r[a+1]
            for s, rs0 in enumerate(r[:-1]):  # Spherical shell
                rs1 = r[s+1]
                if s >= a:
                    # Volume of cylindrical annulus (ra0,ra1) intersecting sphere (rs1)
                    cv[s,a] = 2.0 * theta_rad / 3 * ((rs1**2 - ra0**2)**1.5 - (rs1**2 - ra1**2)**1.5)

                    # Volume of annulus (ra0,ra1) intersecting spherical shell (rs0,rs1)
                    v[s,a] = cv[s,a]
                    if s - a > 0:
                        v[s,a] -= cv[s-1,a]
        self.vol_norm = v / (4. * pi / 3. * r[-1]**3)

    def _create_src_model_components(self):
        """
        Create source model components for each shell corresponding to the
        source model expression.
        """
        # Find the generic components in source model expression
        RE_model = re.compile(r'\b \w+ \b', re.VERBOSE)
        for match in RE_model.finditer(self.srcmodel):
            model_type = match.group()
            self.srcmodel_comps.append(dict(type=model_type,
                                            start=match.start(),
                                            end=match.end()))
        
        # For each shell create the corresponding model components so they can
        # be used later to create composite source models for each dataset
        for shell in range(self.nshell):
            for srcmodel_comp in self.srcmodel_comps:
                model_comp = {}
                model_comp['type'] = srcmodel_comp['type']
                model_comp['name'] = '%s_%d' % (model_comp['type'], shell)
                model_comp['shell'] = shell
                SherpaUI.create_model_component(model_comp['type'], model_comp['name'])
                model_comp['object'] = eval(model_comp['name'])  # Work-around in lieu of accessor
                self.model_comps.append(model_comp)

    def set_source(self, srcmodel='xsphabs*xsapec'):
        """
        Create a source model for each dataset.  A dataset is associated
        with a specific extraction annulus. 

        :param srcmodel: string expression defining source model
        :rtype: None
        """
        self.srcmodel = srcmodel
        self._calc_vol_norm()
        self._create_src_model_components()

        for dataset in self.datasets:
            dataid = dataset['id']
            annulus = dataset['annulus']
            modelexprs = []
            for shell in range(annulus, self.nshell):
                srcmodel = self.srcmodel
                for model_comp in reversed(self.srcmodel_comps):
                    i0 = model_comp['start']
                    i1 = model_comp['end']
                    model_comp_name = '%s_%d' % (model_comp['type'], shell)
                    srcmodel = srcmodel[:i0] + model_comp_name + srcmodel[i1:]
                modelexprs.append('%.5f * %s' % (self.vol_norm[shell, annulus], srcmodel))

            modelexpr = " + ".join(modelexprs)
            print 'Setting source model for dataset %d = %s' % (dataid, modelexpr)
            SherpaUI.set_source(dataid, modelexpr)
        
    def set_bkg_model(self, bkgmodel):
        """
        Create a source model for each dataset.  A dataset is associated
        with a specific extraction annulus. 

        :param bkgmodel: string expression defining background model
        :rtype: None
        """
        self.bkgmodel = bkgmodel

        bkg_norm = {}
        for obsid in self.obsids:
            bkg_norm_name = 'bkg_norm_%d' % obsid
            print 'Creating model component xsconstant.%s' % bkg_norm_name
            SherpaUI.create_model_component('xsconstant', bkg_norm_name)
            bkg_norm[obsid] = eval(bkg_norm_name)  # Uggh, don't know proper model accessor

        for dataset in self.datasets:
            print 'Setting bkg model for dataset %d to bkg_norm_%d' % (dataset['id'], dataset['obsid'])
            SherpaUI.set_bkg_model(dataset['id'], bkg_norm[dataset['obsid']] * bkgmodel)

    def fit(self):
        """
        Do a fit of the model parameters using the "onion-peeling" method:   

         - First fit the outside shell model using the outer annulus spectrum
         - Freeze the model parameters for the outside shell
         - Fit the next inward shell / annulus and freeze those parameters
         - Repeat until all datasets have been fit and all shell parameters determined.
         - Return model parameters to original thawed/frozen status

        :rtype: None
        """
        thawed = []                  # Parameter objects that are not already frozen
        for annulus in reversed(range(self.nshell)):
            dataids = [x['id'] for x in self.datasets if x['annulus'] == annulus]
            print 'Fitting', dataids
            SherpaUI.fit(*dataids)
            for model_comp in self.model_comps:
                name = model_comp['name']
                if model_comp['shell'] == annulus:
                    # Remember parameters that are currently thawed
                    for par in [SherpaUI.get_par('%s.%s'%(name, x))
                                for x in SherpaUI.get_model_pars(name)]:
                        if not par.frozen:
                            thawed.append(par)
                    print 'Freezing', model_comp['name']
                    SherpaUI.freeze(model_comp['name'])

        # Unfreeze parameters
        for par in thawed:
            print 'Thawing', par.fullname
            par.thaw()
            
    def get_density(self):
        """
        Get electron density (cm^-3) for each shell using the standard
        definition of normalization for Xspec thermal models::

         n_e = sqrt(norm * 4*pi * DA^2 * 1e14 * (1+z)^2 / volume * ne_nh_ratio))

         norm        = model normalization from sherpa fit
         DA          = angular size distance (cm)
         volume      = volume (cm^3)
         ne_nh_ratio = 1.18

        Note that the model components for each volume element (intersection of the
        annular cylinder ``a`` with the spherical shell ``s``) are multiplied by a volume
        normalization::

         vol_norm[s,a] = volume[s,a] / v_sphere
         v_sphere = volume of sphere enclosing outer annulus

        With this convention the ``volume`` used in calculating the electron density
        is simply ``v_sphere``.

        :rtype: numpy array of densities (cm^-3) corresponding to shells
        """
        
        ne_nh_ratio = 1.18                # Electron to proton ratio n_e/n_p
        DA_cm = self.angdist
        r_sphere = self.radii[-1] / arcsec_per_rad * DA_cm
        volume = 4 * pi / 3 * r_sphere**3  # volume of sphere enclosing outer shell (cm^3)
        z = self.redshift
        
        dens = []
        for shell in range(self.nshell):
            norm = self.find_norm(shell)
            dens.append(sqrt(norm * 4 * pi * DA_cm**2 * 1e14 * (1.0+z)**2 / volume * ne_nh_ratio))

        return numpy.array(dens)

    def conf(self):
        """
        Run conf on each of the model parameters using the "onion-peeling" method:   

         - First conf the outside shell model using the outer annulus spectrum
         - Freeze the model parameters for the outside shell
         - get confidences for the next inward shell / annulus and freeze those parameters
         - Repeat until all datasets have been conf()-ed and all shell-by-shell error parameters determined.
         - Return model parameters to original thawed/frozen status
         - WARNING: This ignores the correlations between parameters

        :rtype: None
        """
        thawed = []                  # Parameter objects that are not already frozen
        conf_results = []
        this_conf_result = []
        for annulus in reversed(range(self.nshell)):
            dataids = [x['id'] for x in self.datasets if x['annulus'] == annulus]
            print 'Getting shell-by-shell confidence for dataset ', dataids
            SherpaUI.conf(*dataids)
            this_conf_result = SherpaUI.get_conf_results()
            conf_results.append(this_conf_result)
            # get_conf_results().parmins,parmaxes,parnames,parvals
            for model_comp in self.model_comps:
                name = model_comp['name']
                if model_comp['shell'] == annulus:
                    # Remember parameters that are currently thawed
                    for par in [SherpaUI.get_par('%s.%s'%(name, x))
                                for x in SherpaUI.get_model_pars(name)]:
                        if not par.frozen:
                            thawed.append(par)
                    print 'Freezing', model_comp['name']
                    SherpaUI.freeze(model_comp['name'])

        # Unfreeze parameters
        for par in thawed:
            print 'Thawing', par.fullname
            par.thaw()
            
        return conf_results

