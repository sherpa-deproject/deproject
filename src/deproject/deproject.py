"""
Deproject 2-d circular annular spectra to 3-d object properties.

This module implements the "onion-skin" approach popular in X-ray
analysis of galaxy clusters and groups to estimate the three-dimensional
temperature, metallicity, and density distributions of an optically-thin
plasma from the observed (projected) two-dimensional data, arranged in
concentric circular annuli.

:Copyright: Smithsonian Astrophysical Observatory (2009, 2019)
:Author: Tom Aldcroft (taldcroft@cfa.harvard.edu), Douglas Burke (dburke@cfa.harvard.edu)
"""

from collections import defaultdict, OrderedDict
import copy
import logging
from math import pi
import re

import numpy
from astropy.table import Table
from astropy import units as u
from astropy.cosmology import Planck15

from sherpa.plot import plotter
from sherpa.astro import ui
from sherpa.astro.io import read_pha
from sherpa.models.parameter import CompositeParameter

from . import specstack
from . import simplegraph
from . import fieldstore

__all__ = ("Deproject", "deproject_from_xflt")


_sherpa_logger = logging.getLogger('sherpa')


arcsec_per_rad = (u.radian / u.arcsec).to(1)


class Deproject(specstack.SpecStack):
    """Support deprojecting a set of spectra (2-d concentric circular annuli).

    Parameters
    ----------
    radii : AstroPy Quantity representing an angle on the sky
        The edges of each annulus, which must be circular, concentric,
        in ascending order, and >= 0. If there are n annuli then there are n+1
        radii, since the start and end of the sequence must be given.
        The units are expected to be arcsec, arcminute, or degree.
    theta : AstroPy Quantity (scalar or array) representing an angle
        The "fill factor" of each annulus, given by the azimuthal coverage
        of the shell in degrees. The value can be a scalar, so the same
        value is used for all annuli, or a sequence with a length
        matching the number of annuli. Since the annulus assumes circular
        symmetry there is no need to define the starting point of the
        measurement, for cases when the value is less than 360 degrees.
    angdist : None or AstroPy.Quantity, optional
        The angular-diameter distance to the source. If not given then
        it is calculated using the source redshift along with the
        `cosmology` attribute.
    cosmology : None or astropy.cosmology object, optional
        The cosmology used to convert redshift to an angular-diameter
        distance. This is used when `angdist` is None. If `cosmology`
        is None then the `astropy.cosmology.Planck15` Cosmology
        object is used.

    Examples
    --------

    The following highly-simplified example fits a deprojected model
    to data from three annuli - ann1.pi, ann2.pi, and ann3.pi - and
    also calculates errors on the parameters using the confidence
    method::

    >>> dep = Deproject([0, 10, 40, 100] * u.arcsec)
    >>> dep.load_pha('ann1.pi', 0)
    >>> dep.load_pha('ann2.pi', 1)
    >>> dep.load_pha('ann3.pi', 2)
    >>> dep.subtract()
    >>> dep.notice(0.5, 7.0)
    >>> dep.set_source('xsphabs * xsapec')
    >>> dep.set_par('xsapec.redshift', 0.23)
    >>> dep.thaw('xsapec.abundanc')
    >>> dep.set_par('xsphabs.nh', 0.087)
    >>> dep.freeze('xsphabs.nh')
    >>> dep.fit()
    >>> dep.fit_plot('rstat')
    >>> errs = dep.conf()
    >>> dep.conf_plot('density')

    """

    @u.quantity_input(radii='angle', theta='angle', angdist='length')
    def __init__(self, radii,
                 theta=360 * u.deg,
                 angdist=None,
                 cosmology=None):

        nshell = numpy.size(radii) - 1
        if nshell < 1:
            raise ValueError('radii parameter must be a sequence ' +
                             'with at least two values')

        dr = radii[1:] - radii[:-1]
        if numpy.any(dr <= 0):
            raise ValueError('radii parameter must be in increasing order')

        # All values must be >= 0
        #
        if radii[0] < 0:
            raise ValueError('radii must be >= 0')

        self.radii = radii
        self.nshell = nshell

        ntheta = numpy.size(theta)
        if ntheta == 1:
            thetas = numpy.repeat(theta, nshell)
        elif ntheta == nshell:
            thetas = theta
        else:
            raise ValueError('theta must be a scalar or ' +
                             'match the number of annuli')

        theta_min = thetas.min()
        if theta_min <= (0.0 * u.deg):
            raise ValueError('theta must be > 0 degrees')

        theta_max = thetas.max()
        if theta_max > (360.0 * u.deg):
            raise ValueError('theta must be <= 360 degrees')

        self._theta = thetas

        if angdist is not None:
            self._set_angdist(angdist)
        else:
            self._angdist = None

        self._redshift = None
        self._fit_results = None
        self._covar_results = None
        self._conf_results = None

        self._cosmology = Planck15 if cosmology is None else cosmology

        super().__init__()

    def load_pha(self, specfile, annulus):
        if annulus < 0 or annulus >= self.nshell:
            raise ValueError("Expected 0 <= annulus < " +
                             "{} but sent {}".format(self.nshell, annulus))

        super().load_pha(specfile, annulus)

    def _get_redshift(self):
        if self._redshift is None:
            self._redshift = self.find_parval('redshift')
        return self._redshift

    def _set_redshift(self, redshift):
        self._redshift = redshift

    # Perhaps the angular-diameter distance shouldn't be cached if
    # not explicitly set. This lets the value be updated if the
    # redshift or cosmology object is updated. Alternatively, we
    # tell users they have to manually set da if these things
    # change.
    #
    def _get_angdist(self):
        if self._angdist is None:
            da = self.cosmology.angular_diameter_distance(self.redshift)
            self._angdist = da
        return self._angdist

    @u.quantity_input(angdist='length')
    def _set_angdist(self, angdist):

        if angdist <= 0:
            raise ValueError("angdist must be > 0")

        self._angdist = angdist

    redshift = property(_get_redshift, _set_redshift, None,
                        "Source redshift")
    angdist = property(_get_angdist, _set_angdist, None,
                       "Angular size distance (an AstroPy quantity)")

    @property
    def cosmology(self):
        """Return the cosmology object (only used if angdist not set)"""
        return self._cosmology

    def _calc_vol_norm(self):
        r"""Calculate the normalized volumes of the deprojected views.

        Sets the `vol_norm` field to a matrix of the normalized volumes
        of the cylindrical annuli intersecting with the spherical shell.
        The matrix is defined as $volume[i, j] / V_sphere$, where
        $i$ represents the shell and $j$ the annulus (with indexes
        starting at 0), $V_sphere = 4/3 \pi r_o^3$, and $r_o$ is the
        outermost radius of the shells.
        """

        # The units for the radii are not important here
        r = self.radii.value
        theta_rad = self._theta.to_value(u.rad)

        cv = numpy.zeros([self.nshell, self.nshell])
        v = numpy.zeros([self.nshell, self.nshell])

        for a, ra0 in enumerate(r[:-1]):
            # Annulus
            ra1 = r[a + 1]
            ra0sq = ra0**2
            ra1sq = ra1**2

            for s, rs0 in enumerate(r[:-1]):
                # Spherical shell
                rs1 = r[s + 1]
                if s >= a:
                    # Volume of cylindrical annulus (ra0,ra1) intersecting
                    # the sphere (rs1)
                    #
                    rs1sq = rs1**2
                    rterm = (rs1sq - ra0sq)**1.5 - (rs1sq - ra1sq)**1.5
                    cv[s, a] = 2.0 * theta_rad[a] / 3 * rterm

                    # Volume of annulus (ra0,ra1) intersecting the spherical
                    # shell (rs0,rs1)
                    v[s, a] = cv[s, a]
                    if s - a > 0:
                        v[s, a] -= cv[s - 1, a]

        self.vol_norm = v / (4. * pi / 3. * r[-1]**3)

    def _create_name(self, model_name, annulus):
        """Create the name used for a model component for the given annulus.

        Parameters
        ----------
        model_name : str
            The Sherpa model name (e.g. 'xsphabs').
        annulus : int
            The annulus number.

        Returns
        -------
        name : str
            The name of the model component.

        Notes
        -----
        At present there is no real need to allow the naming scheme to
        be changed (e.g. in a sub-class), but it is useful to help
        record where names are created.
        """

        # This is perhaps a bit "over engineered"
        #
        name = '{}_{}'.format(model_name, annulus)
        return name

    def _create_src_model_components(self):
        """Create the model components for each shell."""

        self._reset_model_comps()

        # Find the generic components in source model expression
        # and set up their names.
        #
        RE_model = re.compile(r'\b \w+ \b', re.VERBOSE)
        for match in RE_model.finditer(self.srcmodel):
            model_type = match.group()

            store = dict(type=model_type,
                         start=match.start(),
                         end=match.end())
            self.srcmodel_comps.append(store)

        # For each shell create the corresponding model components so they can
        # be used later to create composite source models for each dataset
        for shell in range(self.nshell):
            for srcmodel_comp in self.srcmodel_comps:
                model_comp = {}
                model_comp['type'] = srcmodel_comp['type']
                name = self._create_name(model_comp['type'], shell)
                model_comp['name'] = name
                model_comp['shell'] = shell
                comp = ui.create_model_component(model_comp['type'], name)
                model_comp['object'] = comp
                self.model_comps.append(model_comp)

    def set_source(self, srcmodel='xsphabs*xsapec'):
        """Create a source model for each annulus.

        Unlike the standard `set_source` command, this version just
        uses the <model name>, not <model name>.<username>, since
        the <username> is automatically created for users by appending
        the annulus number to <model name>.

        Parameters
        ----------
        srcmodel : str, optional
            The source model expression applied to each annulus.

        See Also
        --------
        set_bkg_model, set_par

        Notes
        -----
        The data must have been read in for all the data before calling
        this method (this matches Sherpa, where you can not call set_source
        unless you have already loaded the data to fit).

        Examples
        --------

        The following two calls have the same result: model instances
        called 'xsphabs<annulus>' and 'xsapec<annulus>' are created
        for each annulus, and the source expression for the annulus
        set to their multiplication:

        >>> dep.set_source()
        >>> dep.set_source('xsphabs * xsapec')

        Use the XSPEC vapec model rather than the apec model to
        represent the plasma emission:

        >>> dep.set_source('xsphabs * xsvapec')

        """

        # We can not check that all data has been loaded in (that is,
        # if there are multiple data sets per annulis), but we can at
        # least ensure that there is a dataset loaded
        # for each annulus.
        #
        seen = set([])
        for dataset in self.datasets:
            seen.add(dataset['annulus'])

        expected = set(range(self.nshell))
        diff = sorted(list(expected.difference(seen)))
        if len(diff) == 1:
            raise ValueError("missing data for annulus {}".format(diff[0]))
        elif len(diff) > 0:
            raise ValueError("missing data for annuli {}".format(diff))

        self.srcmodel = srcmodel
        self._calc_vol_norm()
        self._create_src_model_components()

        # TODO: isn't it better to loop over annuli, out to in, to
        #       avoid some repeated work?
        #
        for dataset in self.datasets:
            dataid = dataset['id']
            annulus = dataset['annulus']
            modelexprs = []
            for shell in range(annulus, self.nshell):
                srcmodel = self.srcmodel
                for model_comp in reversed(self.srcmodel_comps):
                    i0 = model_comp['start']
                    i1 = model_comp['end']
                    model_comp_name = self._create_name(model_comp['type'],
                                                        shell)
                    srcmodel = srcmodel[:i0] + model_comp_name + srcmodel[i1:]

                f = self.vol_norm[shell, annulus]
                modelexprs.append('{} * {}'.format(f, srcmodel))

            modelexpr = " + ".join(modelexprs)
            print('Setting source model for dataset %d = %s' % (dataid, modelexpr))
            ui.set_source(dataid, modelexpr)

    def set_bkg_model(self, bkgmodel):
        """Create a background model for each annulus.

        The background model is the same between the annuli, except that
        a scaling factor is added for each annulus (to allow for
        normalization uncertainities). The scaling factors are labelled
        'bkg_norm_<obsid>', and at least one of these must be frozen
        (otherwise it is likely to be degenerate with the background
        normalization, causing difficulties for the optimiser).

        Parameters
        ----------
        bkgmodel : model instance
            The background model expression applied to each annulus.
            Unlike set_source this should be the actual model instance,
            and not a string.

        See Also
        --------
        set_source, set_par

        Examples
        --------

        Model the background with a single power-law component:

        >>> dep.set_bkg_model(xspowerlaw.bpl)

        """

        self.bkgmodel = bkgmodel

        # TODO:
        #   - record the background components in the same way the source
        #     is done
        #   - should the background be allowed to have different components
        #     per annulus?
        #
        bkg_norm = {}
        for obsid in self.obsids:
            bkg_norm_name = 'bkg_norm_%d' % obsid
            print('Creating model component xsconstant.%s' % bkg_norm_name)
            bcomp = ui.create_model_component('xsconstant', bkg_norm_name)
            bkg_norm[obsid] = bcomp

        for dataset in self.datasets:
            print('Setting bkg model for dataset %d to bkg_norm_%d' % (dataset['id'], dataset['obsid']))
            ui.set_bkg_model(dataset['id'],
                             bkg_norm[dataset['obsid']] * bkgmodel)

    def get_shells(self):
        """How are the annuli grouped?

        An annulus may have multiple data sets associated with it, but
        it may also be linked to other annuli due to tied parameters.
        The return value is per group, in the ordering needed for
        the outside-to-inside onion skin fit, where the keys for
        the dictionary are 'annuli' and 'dataids'.

        Returns
        -------
        groups : list of dicts
            Each dictionary has the keys 'annuli' and 'dataids', and
            lists the annuli and data identifiers that are fit together.
            The ordering matches that of the onion-skin approach, so
            the outermost group first.

        See Also
        --------
        get_radii, tie_par

        Examples
        --------

        For a 3-annulus deprojection where there are no parameter ties
        to combine annului:

        >>> dep.get_shells()
        [{'annuli': [2], 'dataids': [2]},
         {'annuli': [1], 'dataids': [1]},
         {'annuli': [0], 'dataids': [0]}]

        After tie-ing the abundance parameter for the outer two shells,
        there are now two groups of annuli:

        >>> dep.tie_par('xsapec.abundanc', 1, 2)
        Tying xsapec_2.Abundanc to xsapec_1.Abundanc
        >>> dep.get_shells()
        [{'annuli': [1, 2], 'dataids': [1, 2]},
         {'annuli': [0], 'dataids': [0]}]

        """

        # Map from model component name (e.g. 'xsapec_2') to shell
        # number.
        #
        cpt_map = {}
        for model_comp in self.model_comps:
            mname = model_comp['name']
            assert mname not in cpt_map
            cpt_map[mname] = model_comp['shell']

        # Find the connected shells/annuli.
        #
        graph = simplegraph.SimpleGraph()
        for shell in range(self.nshell):
            graph.add_link(shell, shell)

        for model_comp in self.model_comps:
            shell = model_comp['shell']
            for par in model_comp['object'].pars:
                if par.link is None:
                    continue

                # For the moment only support tied parameters (i.e.
                # they are set equal). It should be possible to just
                # iterate through the parts of the composite parameter
                # and extract the links (since there could be more than
                # one), but do not try this yet.
                #
                if isinstance(par.link, CompositeParameter):
                    raise ValueError("Parameter link for " +
                                     "{} is not simple".format(par.fullname))

                # If the link is to a "unknown" component then we could
                # iterate through all the source expressions to find
                # the relevant data sets, and hence annuli, but leave
                # that for a later revision since the current assumption
                # is that all source components are handled by deproject.
                #
                linkname = par.link.modelname
                try:
                    lshell = cpt_map[linkname]
                except KeyError:
                    raise RuntimeError("Model component " +
                                       "{} is not ".format(linkname) +
                                       "managed by deproject")

                graph.add_link(shell, lshell)

        # Rely on the shell numbering to be numeric and in ascending
        # order to get the list of shells that must be fit together.
        #
        fit_groups = sorted([sorted(grp) for grp in graph.get_groups()])

        # It is possible for the groups to be invalid here, in that
        # the user could have tied together annuli 2 and 4, but not
        # 3, which breaks the onion-skin approach.
        #
        for grp in fit_groups:
            # ensure that the membership is n, n+1, ..., m with no
            # gaps.
            if len(grp) == 1:
                continue

            grp = numpy.asarray(grp)
            dg = grp[1:] - grp[:-1]
            if numpy.any(dg != 1):
                raise ValueError("Non-consecutive annuli are " +
                                 "tied together: {}".format(grp))

        # What datasets are used for each group?
        #
        out = []
        for anns in fit_groups:

            # What dataset ids are associated with these annuli
            dataids = [x['id'] for x in self.datasets
                       if x['annulus'] in anns]

            out.append({'annuli': anns, 'dataids': dataids})

        return list(reversed(out))

    def get_radii(self, units='arcsec'):
        """What are the radii of the shells?

        Return the inner and outer edge of each annulus, in the given
        units. Physical units (e.g. 'kpc') can only be used if a redshift or
        angular-diameter distance has been set. This does not apply the
        grouping that `get_shells` does.

        Parameters
        ----------
        units : str or astropy.units.Unit, optional
            The name of the units to use for the returned radii. They must
            be an angle - such as 'arcsec' - or a length - such as 'kpc'
            or 'Mpc' (case is important).

        See Also
        --------
        get_shells

        Returns
        -------
        rlo, rhi : astropy.units.Quantity, astropy.units.Quantity
            The inner and outer radius for each annulus.
        """

        # Do we know about this unit? Give a slightly-more helpful
        # message than the default from the AstroPy parser.
        #
        try:
            unit = u.Unit(units)
        except ValueError:
            raise ValueError("Invalid unit: expected a value like " +
                             "'arcsec' or 'kpc'.")

        radii = self.radii.copy()

        if unit.physical_type == 'angle':
            radii = radii.to(unit)

        elif unit.physical_type == 'length':
            # Treat the angular distance value as having length / radian
            rscale = self.angdist / (1 * u.radian)

            # This would convert to m
            # radii = (radii * rscale).decompose()

            radii = (radii * rscale).to(unit)

        else:
            raise u.UnitConversionError("Must be given an angle or length")

        return radii[:-1], radii[1:]

    def guess(self):
        """Guess the starting point by fitting the projected data.

        Use a fitting scheme - based on the suggestion in the XSPEC projct
        documention - to estimate the starting position of the fit (the
        initial fit parameters). This can be useful since it can reduce
        the time taken to fit the deprojected data and help avoid
        the deprojection from getting stuck in a local minimum.

        See Also
        --------
        fit

        Notes
        -----

        Each annulus, from outer to inner, is fit individually, ignoring
        the contribution from any outer annulus. After the fit, the
        model normalisation is corrected for the volume-filling factor of
        the annulus. If there are any tied parameters between annuli then
        these annuli are combined together (fit simultaneously).

        Unlike the Sherpa guess function, this does *not* change the
        limits of any parameter.

        Possible improvements include:
          - re-normalize each spectrum before fitting.
          - transfer the model parameters of the inner-most shell in a
            group to the next set of shells to fit.

        """

        groups = self.get_shells()
        ngroups = len(groups)
        assert (ngroups > 0) & (ngroups <= self.nshell)

        if ngroups != self.nshell:
            print("Note: annuli have been tied together")

        for group in groups:

            annuli = group['annuli']
            nannuli = len(annuli)
            assert nannuli > 0
            dataids = group['dataids']

            msg = 'Projected fit to '
            if len(annuli) == 1:
                msg += 'annulus {} '.format(annuli[0])
            else:
                msg += 'annuli {} '.format(annuli)

            if len(dataids) == 1:
                msg += 'dataset: {}'.format(dataids[0])
            else:
                msg += ' datasets: {}'.format(dataids)

            print(msg)

            orig_models = [(did, ui.get_source(did)) for did in dataids]

            # perhaps this logic should be packaged up
            shells = dict([(x['id'], x['annulus']) for x in self.datasets])

            try:
                # Is there a better way to re-create the "base" model?
                #
                for did in dataids:
                    srcmodel = self.srcmodel
                    shell = shells[did]

                    for model_comp in reversed(self.srcmodel_comps):
                        i0 = model_comp['start']
                        i1 = model_comp['end']
                        model_comp_name = self._create_name(model_comp['type'],
                                                            shell)
                        srcmodel = srcmodel[:i0] + model_comp_name + srcmodel[i1:]

                    ui.set_source(did, srcmodel)

                # TODO: run renormalize on each dataset before the fit

                ui.fit(*dataids)

            finally:
                for did, smdl in orig_models:
                    ui.set_source(did, smdl)

            # Correct the normalization
            #
            fs = self.vol_norm.diagonal()
            for did in dataids:
                shell = shells[did]
                f = fs[shell]

                for mdl in self.model_comps:
                    if mdl['shell'] != shell:
                        continue

                    for par in mdl['object'].pars:
                        if par.name != 'norm':
                            continue

                        par.val /= f

    def _freeze_model_pars(self):
        """Freeze, and return, all thawed parameters in the fit

        Returns
        -------
        pars : list of sherpa.parameter.Parameter instances
            The parameters that have been frozen.

        See Also
        --------
        _thaw_model_pars

        """

        out = []
        for model_comp in self.model_comps:
            out.extend([p for p in model_comp['object'].pars
                        if not p.frozen])

        return out

    def _thaw_model_pars(self, pars, message=False):
        """Thaw the parameters, with optional screen message.

        Parameters
        ----------
        pars : list of sherpa.parameter.Parameter instances
            The parameters to thaw. Note that thaw is called whatever the
            state of a parameter.
        message : bool, optional
            If True then a screen message is displayed for each parameter.

        See Also
        --------
        _freeze_model_pars

        """

        for par in pars:
            if message:
                print('Thawing {}'.format(par.fullname))

            par.thaw()

    def _apply_per_group(self, verb, store):
        """Apply a procedure per onion-skin group (from outer to inner).

        For each shell - run from outer to inner - apply the onion-skin
        approach (so free up the shell but freeze the contribution from
        outer shells) to runfunc - which is given the data ids to use,
        and then store the results of getfunc. Once all the shells have
        been processed return a structure containing the results.

        Parameters
        ----------
        verb : str
            This is the first part of the message displayed to users,
            per annulus.
        store : fieldstore.FieldStore instance
            This object will run the function and then parse the output
            into the necessary form.

        Returns
        -------
        rvals : astropy.table.Table instance
            The data, as a set of columns. The choice of columns is
            controlled by the `store` object. Additional columns include
            'annulus', 'rlo_ang', 'rhi_ang', 'rlo_phys', 'rhi_phys',
            'density', and optionally 'density_lo' and 'density_hi'.

        Notes
        -----
        Any parameter links result in annuli being grouped together
        for a fit: that is, the fit will be a simultaneous fit to all
        the datasets associated with the set of annuli.

        It would be useful to add in extra metadata to the output, and
        take advantage of the units support in AstroPy where possible.
        """

        groups = self.get_shells()
        ngroups = len(groups)
        assert (ngroups > 0) & (ngroups <= self.nshell)

        if ngroups != self.nshell:
            print("Note: annuli have been tied together")

        # Find all the thawed parameters so that they can be
        # restored at the end of the fit, or in case of an error.
        #
        thawed = self._freeze_model_pars()

        # We need to be able to map from dataset id to annulus, and
        # from annulus to model components.
        #
        annulusmap = {x['id']: x['annulus'] for x in self.datasets}

        # Get a list of all the annuli, in ascending order.
        all_annuli = sorted(list({x['annulus'] for x in self.datasets}))

        componentmap = defaultdict(list)
        for mcomp in self.model_comps:
            val = (mcomp['object'], mcomp['type'])
            componentmap[mcomp['shell']].append(val)

        # Store the data as a dictionary of arrays. It would be useful
        # if the FieldStore instance could handle this - since it
        # "knows" what the columns are going to be - but the current
        # implementation calculates these columns when the run method
        # is called (i.e. at run time), rather than before the
        # onion-peel approach is called. It is possible that a
        # re-design would be helpful (such as calculate the parameter
        # names at the start, which would have other useful consequences
        # once more-complicated model expressions are sorted), and then
        # pass that information along, but for now let's see how this
        # works.
        #
        out = OrderedDict()
        try:
            for group in groups:

                annuli = group['annuli']
                nannuli = len(annuli)
                assert nannuli > 0
                dataids = group['dataids']

                msg = '{} '.format(verb)
                if nannuli == 1:
                    msg += 'annulus {} '.format(annuli[0])
                else:
                    msg += 'annuli {} '.format(annuli)

                if len(dataids) == 1:
                    msg += ' dataset: {}'.format(dataids[0])
                else:
                    msg += ' datasets: {}'.format(dataids)

                print(msg)

                res = store.run(annulusmap, componentmap, *dataids)
                assert len(res) == nannuli

                # Extract out the per-shell results and add to the
                # output.
                #
                for shell in annuli:

                    if len(out) == 0:
                        # Note, add in extra fields to the stored fields
                        #
                        out['annulus'] = all_annuli

                        rlo_ang, rhi_ang = self.get_radii(units='arcsec')
                        out['rlo_ang'] = rlo_ang
                        out['rhi_ang'] = rhi_ang

                        rlo_phys, rhi_phys = self.get_radii(units='kpc')
                        out['rlo_phys'] = rlo_phys
                        out['rhi_phys'] = rhi_phys

                        for field in res[shell]:
                            out[field] = [None] * self.nshell

                    for field, value in res[shell].items():
                        assert out[field][shell] is None, shell
                        out[field][shell] = value

                for model_comp in self.model_comps:
                    # Freeze the current annulus
                    if model_comp['shell'] in annuli:
                        print('Freezing {}'.format(model_comp['name']))
                        ui.freeze(model_comp['object'])

        finally:
            self._thaw_model_pars(thawed, message=True)

        for k, vs in out.items():
            assert vs is not None, k

            # It's useful to have NumPy arrays for some of the following
            # calculations, but do not convert those that already have
            # units attached.
            #
            # Convert from None to numpy.NaN (which is assumed to
            # only occur in an error colum (k ends in _lo or _hi)
            # but this restriction is not checked
            #
            if not isinstance(vs, numpy.ndarray):
                out[k] = numpy.asarray([numpy.nan if v is None else v
                                        for v in vs])

        # Add in the density calculation.
        #
        # The norm field should be identified at 'set_source' time
        # so that it doesn't have to be re-discovered each time.
        #
        # For now look for .norm values in the returned structure,
        # but could do this a number of ways.
        #
        normpars = [n for n in out.keys() if n.endswith('.norm')]
        if len(normpars) == 0:
            raise RuntimeError("Unable to find norm parameter!")
        elif len(normpars) > 1:
            raise RuntimeError("Multiple norm parameters found!")

        normpar = normpars[0]
        norms = out[normpar]

        out['density'] = self._calc_density(norms)

        normparlo = '{}_lo'.format(normpar)
        normparhi = '{}_hi'.format(normpar)
        try:
            normlos = out[normparlo]
            normhis = out[normparhi]

            out['density_lo'] = self._calc_density(norms + normlos) - \
                out['density']
            out['density_hi'] = self._calc_density(norms + normhis) - \
                out['density']

        except KeyError:
            pass

        return Table(out)

    def fit(self):
        """Fit the data using the "onion-peeling" method.

        Unlike the normal Sherpa fit, this does not fit all the data
        simultaneously, but instead fits the outermost annulus first,
        then freezes its parameters and fits the annulus inside it,
        repeating this until all annuli have been fit. At the end of
        the fit all the parameters that were frozen are freed. The
        results can also be retrieved with ``get_fit_results``.

        Returns
        -------
        fits : astropy.table.Table instance
            This records per-annulus data, such as the inner and outer
            radius (`rlo_ang`, `rhi_ang`, `rlo_phys`, `rhi_phys`), the
            final fit statistic and change in fit statistic (`statval` and
            `dstatval`), the reduced statistic and q value (as `rstat`
            and `qval`) if appropriate, and the thawed parameter values
            (accessed using <model name>.<par name> syntax, where the
            match is case sensitive).

        See Also
        --------
        conf, covar, get_fit_results, guess, fit_plot

        Notes
        -----

        If there are any tied parameters between annuli then these annuli
        are combined together (fit simultaneously). The results from the
        fits to each annulus can be retrieved after ``fit`` has been called
        with the ``get_fit_results`` method.

        The results have been separated out per annulus, even if several
        annuli were combined in a fit due to tied parameters, and there is
        no information in the returned structure to note this.

        Examples
        --------

        Fit the annuli using the onion-peeling approach, and then plot
        up the reduced statistic for each dataset:

        >>> res = dep.fit()
        >>> plt.clf()
        >>> rmid = 0.5 * (res['rlo_phys'] + res['rhi_phys'])
        >>> plt.plot(rmid, res['rstat'])

        Plot the temperature-abundance values per shell, color-coded
        by annulus:

        >>> plt.clf()
        >>> plt.plot(res['xsapec.kT'], res['xsapec.Abundanc'],
        ...          c=res['annulus'])
        >>> plt.colorbar()
        >>> plt.xlabel('kT')
        >>> plt.ylabel('Abundance')

        Plot up the temperature distibution as a function of radius
        from the fit::

        >>> dep.fit()
        >>> dep.fit_plot('xsmekal.kt')

        """

        # For now return nothing
        self._fit_results = None
        self._fit_results = self._apply_per_group('Fitting',
                                                  fieldstore.FitStore())
        return self._fit_results

    def covar(self):
        """Estimate errors using covariance, using the "onion-peeling" method.

        It is assumed that ``fit`` has been called. The results can also be
        retrieved with ``get_covar_results``.

        Returns
        -------
        errors : astropy.table.Table instance
            This records per-annulus data, such as the inner and outer
            radius (`rlo_ang`, `rhi_ang`, `rlo_phys`, `rhi_phys`), the
            sigma and percent values, and parameter results (accessed
            using <model name>.<par name>, <model name>.<par name>_lo,
            and <model name>.<par name>_hi syntax, where the match is
            case sensitive). The _lo and _hi values are symmetric for
            covar, that is the _lo value will be the negative of the
            _hi value.

        See Also
        --------
        conf, fit, get_covar_results, covar_plot

        Examples
        --------

        Run a fit and then error analysis, then plot up the abundance
        against temperature values including the error bars. Since
        the covariance routine returns symmetric error bars, the
        <param>_hi values are used in the plot::

        >>> dep.fit()
        >>> errs = dep.covar()
        >>> kt, abund = errs['xsapec.kT'], errs['xsapec.Abundanc']
        >>> dkt = errs['xsapec.kT_hi']
        >>> dabund = errs['xsapec.Abundanc_hi']
        >>> plt.clf()
        >>> plt.errorbar(kt, abund, xerr=dkt, yerr=dabund, fmt='.')

        Plot up the temperature distibution as a function of radius,
        including the error bars calculated by the covar routine::

        >>> dep.fit()
        >>> dep.covar()
        >>> dep.covar_plot('xsmekal.kt')

        """

        self._covar_results = None
        self._covar_results = self._apply_per_group('Covariance for',
                                                    fieldstore.CovarStore())
        return self._covar_results

    def conf(self):
        """Estimate errors using confidence, using the "onion-peeling" method.

        It is assumed that ``fit`` has been called. The results can also be
        retrieved with ``get_conf_results``.

        Returns
        -------
        errors : astropy.table.Table instance
            This records per-annulus data, such as the inner and outer
            radius (`rlo_ang`, `rhi_ang`, `rlo_phys`, `rhi_phys`), the
            sigma and percent values, and parameter results (accessed
            using <model name>.<par name>, <model name>.<par name>_lo,
            and <model name>.<par name>_hi syntax, where the match is
            case sensitive).

        See Also
        --------
        covar, fit, get_conf_results, conf_plot

        Examples
        --------

        Run a fit and then error analysis, then plot up the abundance
        against temperature values including the error bars. Note that
        the Matplotlib `errorbar` routine requires "positive" error values
        whereas the <param>_lo values are negative, hence they are
        negated in the creation of ``dkt`` and ``dabund``::

        >>> dep.fit()
        >>> errs = dep.conf()
        >>> kt, abund = errs['xsapec.kT'], errs['xsapec.Abundanc']
        >>> ktlo, kthi = errs['xsapec.kT_lo'], errs['xsapec.kT_hi']
        >>> ablo, abhi = errs['xsapec.Abundanc_lo'], errs['xsapec.Abundanc_hi']
        >>> dkt = np.vstack((-ktlo, kthi))
        >>> dabund = np.vstack((-ablo, abhi))
        >>> plt.clf()
        >>> plt.errorbar(kt, abund, xerr=dkt, yerr=dabund, fmt='.')

        Plot up the temperature distibution as a function of radius,
        including the error bars calculated by the conf routine::

        >>> dep.fit()
        >>> dep.conf()
        >>> dep.conf_plot('xsmekal.kt')

        """

        self._conf_results = None
        self._conf_results = self._apply_per_group('Confidence for',
                                                   fieldstore.ConfStore())
        return self._conf_results

    def get_fit_results(self):
        """What are the fit results, per annulus?

        This returns the fit result for each annulus from the last time
        that the ``fit`` method was called. It *does not* check to see
        if anything has changed since the last ``fit`` call (e.g.
        parameters being tied together or untied, or a manual fit
        to a shell).  Note that ``get_shells`` should be used to find out
        if the shells were grouped together.

        Returns
        -------
        fits : astropy.table.Table instance
            This records per-annulus data, such as the inner and outer
            radius (`rlo_ang`, `rhi_ang`, `rlo_phys`, `rhi_phys`), the
            final fit statistic and change in fit statistic (`statval` and
            `dstatval`), the reduced statistic and q value (as `rstat`
            and `qval`) if appropriate, and the thawed parameter values
            (accessed using <model name>.<par name> syntax, where the
            match is case sensitive).

        See Also
        --------
        fit, get_conf_results, get_covar_results, get_radii, get_shells,
        fit_plot

        """

        if self._fit_results is None:
            raise ValueError("The fit method has not been called")

        return copy.deepcopy(self._fit_results)

    def get_covar_results(self):
        """What are the covar results, per annulus?

        This returns the fit result for each annulus from the last time
        that the ``covar`` method was called. It *does not* check to see
        if anything has changed since the last ``covar`` call (e.g.
        parameters being tied together or untied, or a manual fit
        to a shell). Note that ``get_shells`` should be used to find out
        if the shells were grouped together.

        Returns
        -------
        errors : astropy.table.Table instance
            This records per-annulus data, such as the inner and outer
            radius (`rlo_ang`, `rhi_ang`, `rlo_phys`, `rhi_phys`), the
            sigma and percent values, and parameter results (accessed
            using <model name>.<par name>, <model name>.<par name>_lo,
            and <model name>.<par name>_hi syntax, where the match is
            case sensitive).

        See Also
        --------
        fit, get_conf_results, get_fit_results, get_radii, get_shells,
        covar_plot

        """

        if self._covar_results is None:
            raise ValueError("The covar method has not been called")

        return copy.deepcopy(self._covar_results)

    def get_conf_results(self):
        """What are the conf results, per annulus?

        This returns the fit result for each annulus from the last time
        that the ``conf`` method was called. It *does not* check to see
        if anything has changed since the last ``conf`` call (e.g.
        parameters being tied together or untied, or a manual fit
        to a shell). Note that ``get_shells`` should be used to find out
        if the shells were grouped together (although this can be
        reconstructed from the `datasets` field of each `ErrorEstResults`
        instance).

        Returns
        -------
        errors : astropy.table.Table instance
            This records per-annulus data, such as the inner and outer
            radius (`rlo_ang`, `rhi_ang`, `rlo_phys`, `rhi_phys`), the
            sigma and percent values, and parameter results (accessed
            using <model name>.<par name>, <model name>.<par name>_lo,
            and <model name>.<par name>_hi syntax, where the match is
            case sensitive).

        See Also
        --------
        fit, get_covar_results, get_fit_results, get_radii, get_shells,
        conf_plot

        """

        if self._conf_results is None:
            raise ValueError("The conf method has not been called")

        return copy.deepcopy(self._conf_results)

    def _calc_density(self, norms, ne_nh_ratio=1.18):
        """Calculate the electron density for each shell.

        This performs the calculation described in `get_density`.

        Parameters
        ----------
        norms : sequence of float
            The normalization values, in annulus order.
        ne_hh_ratio : float, optional
            The n_e to n_h ratio (default 1.18).

        Returns
        -------
        dens : astropy.units.quantity.Quantity instance
            The densities calculated for each shell, in units of cm^-3.

        """

        if len(norms) != self.nshell:
            raise ValueError("norms has wrong length")

        # Manual descontruction/reconstruction of units
        #
        DA_cm = self.angdist.to_value(u.cm)
        rmax_rad = self.radii[-1].to_value(u.rad)
        z = self.redshift

        r_sphere = rmax_rad * DA_cm

        # volume of sphere enclosing outer shell (cm^3)
        #
        # volume = 4 * pi / 3 * r_sphere**3
        # factor = 4 * pi * DA_cm**2 * 1e14 * (1.0 + z)**2 / volume * ne_nh_ratio
        #
        # and after manual cancellation of the 4 pi terms
        #
        vterm = 1.0 / 3 * r_sphere**3
        factor = DA_cm**2 * 1e14 * (1.0 + z)**2 / vterm * ne_nh_ratio

        return numpy.sqrt(factor * numpy.asarray(norms)) * u.cm**(-3)

    def get_density(self):
        """Calculate the electron density for each shell.

        Convert the model normalzations (assumed to match the standard
        definition for XSPEC thermal-plasma models) for each shell.

        Returns
        -------
        dens : astropy.units.quantity.Quantity instance
            The densities calculated for each shell, in units of cm^-3.

        See Also
        --------
        find_norm

        Notes
        -----
        The electron density is taken to be::

          n_e^2 = norm * 4*pi * DA^2 * 1e14 * (1+z)^2 / volume * ne_nh_ratio

        where::

          norm        = model normalization from sherpa fit
          DA          = angular size distance (cm)
          volume      = volume (cm^3)
          ne_nh_ratio = 1.18

        The model components for each volume element (the intersection of the
        annular cylinder ``a`` with the spherical shell ``s``) are multiplied
        by a volume normalization::

          vol_norm[s,a] = volume[s,a] / v_sphere
          v_sphere = volume of sphere enclosing outer annulus

        With this convention the ``volume`` used in calculating the electron
        density is simply ``v_sphere``.

        """

        norms = [self.find_norm(s) for s in range(self.nshell)]
        return self._calc_density(norms)

    def _radial_plot(self, plottitle, xunits, ys, ylabel,
                     dys=None,
                     xlog=True, ylog=False,
                     overplot=False, clearwindow=True):
        """Create a plot of the data versus radius (of the annuli).

        Parameters
        ----------
        plottitle : str
            The title for the plot
        xunits : str or astropy.units.Unit
            The X-axis units (a length or angle, such as 'Mpc' or
            'arcsec', where the case is important).
        ys : sequence of float
            The Y values to plot (must be in annuli order).
        ylabel : str
            The label for the Y axis
        dys : None or ndarray, optional
            The error bars on the y axis. This can be None or a ndarray
            of one or two dimensions (N points or N by 2).
        xlog : bool, optional
            Should the x axis be drawn with a log scale (default True)?
        ylog : bool, optional
            Should the y axis be drawn with a log scale (default False)?
        overplot : bool, optional
            Clear the plot or add to existing plot?
        clearwindow : bool, optional
            How does this interact with overplot?

        """

        rlo, rhi = self.get_radii(units=xunits)

        # drop units support immediately as ChIPS doesn't recognize
        # this (can support them in matplotlib, but given the
        # Sherpa plotting API it isn't clear how well supported it
        # would be)
        #
        rmid = (rlo.value + rhi.value) / 2
        dr = rhi.value - rlo.value

        # Attempt to handle LaTeX differences between the backends,
        # but the support is *very* limited so may not work here.
        #
        # The aim is to support the matplotlib backend, with minimal
        # support for ChIPS.
        #
        xlabel = _add_unit('Radius', rlo)
        if ylabel.find('_') > -1 or ylabel.find('^') > -1:
            # Unfortunately the matplotlib version is a "global"
            # check, so doesn't check if parts of the term are
            # already enclosed in '$'. This is a problem for those
            # labels that have AstroPy unit strings, since they
            # have already been protected.
            #
            if not (plotter.name == 'pylab' and ylabel.find('$') > -1):
                ylabel = plotter.get_latex_for_string(ylabel)

        prefs = plotter.get_data_plot_defaults()
        prefs['xerrorbars'] = True

        # We handle error bars manually for ChIPS (it has to be done
        # for asymmetric Y errors, but there also seems to be issues
        # with the X axis errors not being drawn which I do not
        # want to investigate too much just right now).
        #
        manual_errors = plotter.name == 'chips' and dys is not None

        prefs['yerrorbars'] = dys is not None
        if manual_errors:
            prefs['yerrorbars'] = False

        prefs['xlog'] = xlog
        prefs['ylog'] = ylog

        # Access the underlying plot machinery directly, rather than
        # use the sherpa.plot.DataPlot object, since Sherpa does not
        # support asymmetric errors but plotter.plot does, at least
        # for the pylab backend.
        #
        try:
            plotter.begin()
            plotter.plot(rmid, ys, dys, dr, plottitle, xlabel, ylabel,
                         overplot, clearwindow, **prefs)

            # For some reason the X error bar isn't being drawn with ChIPS
            # so force it.
            #
            if plotter.name == 'chips':
                import pychips

                # Assume the current curve is the data we have just plotted
                # and we do not want to replot the symbol.
                #
                crv = pychips.get_curve()
                crv.symbol.style = 'none'
                crv.err.up = True
                crv.err.down = True
                crv.err.left = True
                crv.err.right = True

                ndim = numpy.asarray(dys).ndim
                if ndim == 2:
                    dylo = dys[0]
                    dyhi = dys[1]
                elif ndim == 1:
                    dylo = dys
                    dyhi = dys
                else:
                    dylo = None
                    dyhi = None

                errs = [dylo, dyhi, dr / 2, dr / 2]
                pychips.add_curve(rmid, ys, errs, crv)

        except BaseException as exc:
            plotter.exceptions()
            raise exc

        else:
            plotter.end()

    def par_plot(self, par, units='kpc',
                 xlog=True, ylog=False,
                 overplot=False, clearwindow=True):
        """Plot up the parameter as a function of radius.

        This plots up the current parameter values. The ``fit_plot``,
        ``conf_plot``, and ``covar_plot`` routines display the fit
        and error results for these parameters.

        Parameters
        ----------
        par : str
            The parameter name, specified as <model_type>.<par_name>.
        units : str or astropy.units.Unit, optional
            The X-axis units (a length or angle, such as 'Mpc' or
            'arcsec', where the case is important).
        xlog : bool, optional
            Should the x axis be drawn with a log scale (default True)?
        ylog : bool, optional
            Should the y axis be drawn with a log scale (default False)?
        overplot : bool, optional
            Clear the plot or add to existing plot?
        clearwindow : bool, optional
            How does this interact with overplot?

        See Also
        --------
        conf_plot, covar_plot, density_plot, fit_plot

        Examples
        --------

        Plot the temperature as a function of radius.

        >>> dep.par_plot('xsapec.kt')

        Label the radii with units of arcminutes for the abundanc
        parameter of the xsapec model:

        >>> dep.par_plot('xsapec.abundanc', units='arcmin')

        """

        # Assume par is "model_name.par_name" and we do not have to
        # worry about case for model_name, but may have to for par_name
        #
        mname, pname = self._split_parname(par)
        pname = pname.lower()
        cpts = [cpt['object'] for cpt in self.model_comps
                if cpt['type'] == mname]

        # Probably can not get here and this happen (thanks to the get_par)
        # call, but just in case
        if len(cpts) == 0:
            raise ValueError("No matching model {} for par={}".format(mname,
                                                                      par))

        # Assume they are all the same (they better be)
        #
        yunits = None
        for p in cpts[0].pars:
            if p.name.lower() != pname:
                continue

            yunits = p.units
            break

        # Also should not happen, so report if it does but do not
        # error out
        if yunits is None:
            print("WARNING: unable to find match for parameter {}".format(par))
            yunits = ''

        ylabel = par
        if yunits.strip() != '':
            ylabel += " ({})".format(yunits)

        pvals = self.get_par(par)
        self._radial_plot(par, units, pvals, ylabel,
                          xlog=xlog, ylog=ylog,
                          overplot=overplot, clearwindow=clearwindow)

    def density_plot(self, units='kpc',
                     xlog=True, ylog=True,
                     overplot=False, clearwindow=True):
        """Plot up the electron density as a function of radius.

        The density is displayed with units of cm^-3. This plots up the
        density calculated using the current normalization parameter
        values.  The ``fit_plot``, ``conf_plot``, and ``covar_plot``
        routines display the fit and error results for these parameters.

        Parameters
        ----------
        units : str or astropy.units.Unit, optional
            The X-axis units (a length or angle, such as 'Mpc' or
            'arcsec', where the case is important).
        xlog : bool, optional
            Should the x axis be drawn with a log scale (default True)?
        ylog : bool, optional
            Should the y axis be drawn with a log scale (default False)?
        overplot : bool, optional
            Clear the plot or add to existing plot?
        clearwindow : bool, optional
            How does this interact with overplot?

        See Also
        --------
        conf_plot, covar_plot, fit_plot, par_plot

        Examples
        --------

        Plot the density as a function of radius.

        >>> dep.density_plot()

        Label the radii with units of arcminutes:

        >>> dep.density_plot(units='arcmin')

        """

        nes = self.get_density().value

        # Unfortunately the LaTeX emulation in the two backends is not
        # comparable, which limits the fidelity of the label.
        #
        # ylabel = 'n$_e$ (cm$^{-3}$)'
        ylabel = r'n_e\ (\mathrm{cm^{-3}})'

        self._radial_plot('density', units, nes, ylabel,
                          xlog=xlog, ylog=ylog,
                          overplot=overplot, clearwindow=clearwindow)

    def fit_plot(self, field, results=None,
                 units='kpc',
                 xlog=True, ylog=False,
                 overplot=False, clearwindow=True):
        """Plot up the fit results as a function of radius.

        This method can be used to plot up the last fit results or
        a previously-stored set. To include error bars on the
        dependent values use the `conf_plot` or `covar_plot` methods.

        Parameters
        ----------
        field : str
            The column to plot from the fit results (the match is case
            insensitive).
        results : None or astropy.table.Table instance
            The return value from the ``fit`` or ``get_fit_results``
            methods.
        units : str or astropy.units.Unit, optional
            The X-axis units (a length or angle, such as 'Mpc' or
            'arcsec', where the case is important).
        xlog : bool, optional
            Should the x axis be drawn with a log scale (default True)?
        ylog : bool, optional
            Should the y axis be drawn with a log scale (default False)?
        overplot : bool, optional
            Clear the plot or add to existing plot?
        clearwindow : bool, optional
            How does this interact with overplot?

        See Also
        --------
        fit, get_fit_results, conf_plot, covar_plot, density_plot, par_plot

        Examples
        --------

        Plot the temperature as a function of radius from the last
        fit:

        >>> dep.fit_plot('xsapec.kt')

        Plot the reduced fit statistic from the last fit:

        >>> dep.fit_plot('rstat')

        Plot the density with the radii labelled in arcminutes and the
        density shown on a log scale:

        >>> dep.fit_plot('density', units='arcmin', ylog=True)

        Overplot the current fit results on those from a previous fit,
        where ``fit1`` was returned from the ``fit`` or ``get_fit_results``
        methods:

        >>> dep.fit_plot('xsapec.abundanc', results=fit1)
        >>> dep.fit_plot('xsapec.abundanc', overplot=True)

        """

        if results is None:
            plotdata = self.get_fit_results()
        else:
            plotdata = results

        try:
            ys = plotdata[field]
        except KeyError:
            flower = field.lower()
            names = [n for n in plotdata.keys() if n.lower() == flower]
            if len(names) == 0:
                raise ValueError("Unrecognized field {}".format(field))
            elif len(names) > 1:
                raise RuntimeError("Multiple fields match {}".format(field))

            field = names[0]
            ys = plotdata[field]

        ylabel = _add_unit(field, ys)
        self._radial_plot(field, units, ys, ylabel,
                          xlog=xlog, ylog=ylog,
                          overplot=overplot, clearwindow=clearwindow)

    def conf_plot(self, field, results=None,
                  units='kpc',
                  xlog=True, ylog=False,
                  overplot=False, clearwindow=True):
        """Plot up the confidence errors as a function of radius.

        This method can be used to plot up the last conf results or
        a previously-stored set. Any error bars are shown at the
        scale they were calculated (as given by the ``sigma`` and
        ``percent`` columns of the results).

        Parameters
        ----------
        field : str
            The column to plot from the fit results (the match is case
            insensitive).
        results : None or astropy.table.Table instance
            The return value from the ``conf`` or ``get_conf_results``
            methods.
        units : str or astropy.units.Unit, optional
            The X-axis units (a length or angle, such as 'Mpc' or
            'arcsec', where the case is important).
        xlog : bool, optional
            Should the x axis be drawn with a log scale (default True)?
        ylog : bool, optional
            Should the y axis be drawn with a log scale (default False)?
        overplot : bool, optional
            Clear the plot or add to existing plot?
        clearwindow : bool, optional
            How does this interact with overplot?

        See Also
        --------
        fit, get_conf_results, fit_plot, covar_plot, density_plot, par_plot

        Notes
        -----
        Error bars are included on the dependent axis if the results
        contain columns that match the requested field with suffixes
        of '_lo' and '_hi'. These error bars are asymmetric, which is
        different to ``covar_plot``.

        If a limit is missing (i.e. it is a NaN) then no error bar is
        drawn. This can make it look like the error is very small.

        Examples
        --------

        Plot the temperature as a function of radius from the last
        fit, including error bars:

        >>> dep.conf_plot('xsapec.kt')

        Plot the density with the radii labelled in arcminutes and the
        density shown on a log scale:

        >>> dep.conf_plot('density', units='arcmin', ylog=True)

        Overplot the current conf results on those from a previous fit,
        where ``conf1`` was returned from the ``conf`` or ``get_conf_results``
        methods:

        >>> dep.conf_plot('xsapec.abundanc', results=conf1)
        >>> dep.conf_plot('xsapec.abundanc', overplot=True)

        """

        if results is None:
            plotdata = self.get_conf_results()
        else:
            plotdata = results

        try:
            ys = plotdata[field]
        except KeyError:
            flower = field.lower()
            names = [n for n in plotdata.keys() if n.lower() == flower]
            if len(names) == 0:
                raise ValueError("Unrecognized field {}".format(field))
            elif len(names) > 1:
                raise RuntimeError("Multiple fields match {}".format(field))

            field = names[0]
            ys = plotdata[field]

        try:
            flo = '{}_lo'.format(field)
            fhi = '{}_hi'.format(field)
            dys = numpy.vstack((-plotdata[flo], plotdata[fhi]))
        except KeyError:
            dys = None

        ylabel = _add_unit(field, ys)
        self._radial_plot(field, units, ys, ylabel, dys=dys,
                          xlog=xlog, ylog=ylog,
                          overplot=overplot, clearwindow=clearwindow)

    def covar_plot(self, field, results=None,
                   units='kpc',
                   xlog=True, ylog=False,
                   overplot=False, clearwindow=True):
        """Plot up the covariance errors as a function of radius.

        This method can be used to plot up the last covar results or
        a previously-stored set. Any error bars are shown at the
        scale they were calculated (as given by the ``sigma`` and
        ``percent`` columns of the results).

        Parameters
        ----------
        field : str
            The column to plot from the fit results (the match is case
            insensitive).
        results : None or astropy.table.Table instance
            The return value from the ``covar`` or ``get_covar_results``
            methods.
        units : str or astropy.units.Unit, optional
            The X-axis units (a length or angle, such as 'Mpc' or
            'arcsec', where the case is important).
        xlog : bool, optional
            Should the x axis be drawn with a log scale (default True)?
        ylog : bool, optional
            Should the y axis be drawn with a log scale (default False)?
        overplot : bool, optional
            Clear the plot or add to existing plot?
        clearwindow : bool, optional
            How does this interact with overplot?

        See Also
        --------
        fit, get_covar_results, fit_plot, conf_plot, density_plot, par_plot

        Notes
        -----
        Error bars are included on the dependent axis if the results
        contain columns that match the requested field with the suffix
        '_hi'. The error bars are therefore symmetric, which is
        different to ``conf_plot``.

        If a limit is missing (i.e. it is a NaN) then no error bar is
        drawn. This can make it look like the error is very small.

        Examples
        --------

        Plot the temperature as a function of radius from the last
        fit, including error bars:

        >>> dep.covar_plot('xsapec.kt')

        Plot the density with the radii labelled in arcminutes and the
        density shown on a log scale:

        >>> dep.covar_plot('density', units='arcmin', ylog=True)

        Overplot the current covar results on those from a previous fit,
        where ``covar1`` was returned from the ``covar`` or
        ``get_covar_results`` methods:

        >>> dep.covar_plot('xsapec.abundanc', results=covar1)
        >>> dep.covar_plot('xsapec.abundanc', overplot=True)

        """

        if results is None:
            plotdata = self.get_covar_results()
        else:
            plotdata = results

        try:
            ys = plotdata[field]
        except KeyError:
            flower = field.lower()
            names = [n for n in plotdata.keys() if n.lower() == flower]
            if len(names) == 0:
                raise ValueError("Unrecognized field {}".format(field))
            elif len(names) > 1:
                raise RuntimeError("Multiple fields match {}".format(field))

            field = names[0]
            ys = plotdata[field]

        try:
            fhi = '{}_hi'.format(field)
            dys = plotdata[fhi]
        except KeyError:
            dys = None

        ylabel = _add_unit(field, ys)
        self._radial_plot(field, units, ys, ylabel, dys=dys,
                          xlog=xlog, ylog=ylog,
                          overplot=overplot, clearwindow=clearwindow)


def _add_unit(label, xs):
    """Add a unit string to the label if available.

    This also converts a label of density to n_e.

    Parameters
    ----------
    label : str
        The label
    xs : values
        The numeric values

    Returns
    -------
    label : str
        If xs is an AstroPy Quantity then " (<unit>)" is added to the
        label, using LaTeX (inline) formatting. If the input label is
        'density' then it is replaced by '$n_e$'.

    """

    if label == 'density':
        label = '$n_e$'

    try:
        unit = xs.unit.to_string(format='latex_inline')
        return r"{} ({})".format(label, unit)
    except AttributeError:
        return label


def _get_xflt_keys(infile):
    """Return the XFLT0001 to XFLT0005 keys in infile.

    Parameters
    ----------
    infile : str
        The name of the PHA file with the XFLT000n keywords.

    Returns
    -------
    xflts : list of float
        The XFLT0001 to XFLT0005 numbers, in order.

    Raises
    ------
    IOError
        If the file does not exist or a keyword is missing.

    Notes
    -----
    The Sherpa logging level is temporarily changed to
    the WARN level (if it is lower than this) when reading in
    the file, to avoid excessive screen output.

    """

    # Hide sherpa logging here
    olvl = _sherpa_logger.getEffectiveLevel()
    if olvl < logging.WARN:
        _sherpa_logger.setLevel(logging.WARN)

    try:
        pha = read_pha(infile)
    finally:
        _sherpa_logger.setLevel(olvl)

    hdr = pha.header
    try:
        out = [hdr['XFLT000{}'.format(i)] for i in range(1, 6)]
    except KeyError as e:
        raise IOError("PHA file {} is missing keyword {}".format(infile,
                                                                 e))

    return out


def _find_xflt_files(pat):
    """Return the XFLT keywords in the files matching the input pattern.

    This is only for circular annuli.

    Parameters
    ----------
    pat : str
        The pattern representing the files to use (e.g. 'ann*.pi'). If
        the stk module (provided by CIAO) is available then CIAO stack
        syntax (e.g. "@files.lis" to read the names from files.lis) can
        be used. The files do not need to be in order (of increasing
        radii).

    Returns
    -------
    vals : list of (str, list)
        The file name and corresponding five XFLT values. The list is
        ordered by increasing radius of the annuli.

    Raises
    ------
    IOError
        When a non-circular annulus is encountered.

    Notes
    -----
    There is no attempt to check that the annuli are consecutive and
    not overlapping.
    """

    try:
        import stk
        infiles = stk.build(pat)
    except ImportError:
        # Assume not in a CIAO installation, so only support a glob
        import glob
        infiles = glob.glob(pat)

    if len(infiles) == 0:
        raise IOError("No matches to {}".format(pat))

    out = []
    for infile in infiles:
        xflt = _get_xflt_keys(infile)

        # numeric comparison not great, do we need some sort of
        #
        #
        if (xflt[0] != xflt[1]) | (xflt[2] != 0):
            raise IOError("Only circular annulli allowed: {}".format(infile))

        out.append((infile, xflt))

    # Sort on radius
    #
    return sorted(out, key=lambda x: x[1][0])


def deproject_from_xflt(pat, rscale,
                        rinner=0,
                        angdist=None,
                        cosmology=None):
    """Set up the projection object from XFLT keywords in the PHA files.

    When using the XSPEC projct model, values are read from XFLT keywords
    (as used by the XSPEC deprojection code [1]_) rather than being specified
    manually. This function creates a Deproject object and loads in a set
    of PHA files matching a pattern, using the XFLT keywords to set radii
    and theta values. The annuli *must* be circular.

    Parameters
    ----------
    pat : str
        The pattern representing the files to read in. If the stk module,
        provided by CIAO, is available then CIAO stack syntax [2]_ can be
        used. The order of the files does not matter, but it is currently
        assumed that there is only one file per annulus.
    rscale : AstroPy quantity
        The scaling factor used to convert the XFLT radii (XFLT001 and
        XFLT002 keywords) to an angle. If the values are in arcseconds
        then ``rscale`` would be set to ``1 * u.arcsec``.
    rinner : float, optional
        The inner radius of the central annulus, in the same system as
        the XFLT0001 and XFLT002 keyword values (this is a unitless
        value).
    angdist : None or AstroPy.Quantity, optional
        The angular-diameter distance to the source. If not given then
        it is calculated using the source redshift along with the
        `cosmology` attribute.
    cosmology : None or astropy.cosmology object, optional
        The cosmology used to convert redshift to an angular-diameter
        distance. This is used when `angdist` is None. If `cosmology`
        is None then the `astropy.cosmology.Planck15` Cosmology
        object is used.

    Returns
    -------
    dep : Deproject instance
        The deproject instance with the files loaded and associated
        with the correct annuli.

    Notes
    -----
    This currently is *not* guaranteed to support multiple data sets in
    the same annulus. There is no check that the annuli are touching and
    do not overlap.

    References
    ----------

    .. [1] https://asd.gsfc.nasa.gov/XSPECwiki/projct_model

    .. [2] http://cxc.harvard.edu/ciao/ahelp/stack.html

    Examples
    --------

    Create a Deproject instance from the files matching the pattern
    "src*.pi", whose XFLT radii are in ACIS pixels:

    >>> dep = deproject_from_xflt('src*.pi', 0.492 * u.arcsec)

    When used in CIAO, the stack syntax can be used to specify the
    files, so if the file clus.stk contains the file names, one
    per line, then the following will read them in and create a
    Deproject instance. In this case the XFLT radii are in arcminutes:

    >>> dep = deproject_from_xflt('@clus.stk', 1 * u.arcmin)

    """

    infiles = _find_xflt_files(pat)
    print("# Found {} annuli".format(len(infiles)))

    # TODO: check that annuli are touching (ie that outer[i] = inner[i+1])
    xs = [rinner] + [x[1][0] for x in infiles]
    radii = numpy.asarray(xs) * rscale
    print("# Inner radius: {}".format(radii[0]))
    print("# Outer radius: {}".format(radii[-1]))

    # assume that theta_min < theta_max so this is valid
    #
    thetas = numpy.asarray([x[1][4] - x[1][3] for x in infiles])

    tmin = thetas.min()
    tmax = thetas.max()
    if tmin <= 0:
        raise ValueError("Found theta range of {}; ".format(tmin) +
                         "XFLT0004/5 keywords assumed to be in " +
                         "clockwise order")
    if tmax > 360:
        raise ValueError("Found theta range of {}; ".format(tmax) +
                         "XFLT0004/5 keywords assumed to be in " +
                         "clockwise order")

    if tmin == tmax:
        print("# Theta = {} degrees".format(tmin))
    else:
        print("# Theta ranges from {} to {} degrees".format(tmin, tmax))

    print("#")

    thetas = thetas * u.deg
    dep = Deproject(radii,
                    theta=thetas,
                    angdist=angdist,
                    cosmology=cosmology)

    for i, x in enumerate(infiles):
        dep.load_pha(x[0], annulus=i)

    return dep
