"""
Manipulate a stack of spectra in Sherpa.

The methods in the SpecStack class provide a way to automatically apply
familiar Sherpa commands such as `set_par`_ or `freeze`_ or `plot_fit`_
to a stack of PHA spectra.  This simplifies simultaneous fitting of
multiple spectra.

Note that the :mod:`specstack` module is currently distributed within with the
:mod:`deproject` package.  `Specstack` is not yet fully documented or tested
outside the context of `deproject`.

:Copyright: Smithsonian Astrophysical Observatory (2009, 2019)
:Author: Tom Aldcroft (taldcroft@cfa.harvard.edu), Douglas Burke (dburke@cfa.harvard.edu)
"""

import re
import numpy
from sherpa.astro import ui

try:
    from sherpa.astro.plot import backend
    backend_name = backend.name
except ImportError:
    backend_name = 'none'

# Trying to move away from direct access to I/O and plotting code.
#
if backend_name == 'pychips':
    try:
        import pychips
    except ImportError:
        pass

elif backend_name == 'pylab':
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        pass


def _sherpa_plot_func(func):
    "Wrap functions which are called once per dataset per annulus"
    def _sherpa_plot_func_(self, *args, **kwargs):
        self._sherpa_plot(func, *args, **kwargs)
    return _sherpa_plot_func_


def _sherpa_plot_func_single(func):
    "Wrap functions which are called once per annulus"
    def _sherpa_plot_func_(self, *args, **kwargs):
        kwargs['single_call_per_shell'] = True
        self._sherpa_plot(func, *args, **kwargs)
    return _sherpa_plot_func_


class SpecStack:
    """Manipulate a stack of spectra in Sherpa.

    This `SpecStack` class provides a number of wrappers around Sherpa
    routines that handle loading data, setting the source model,
    setting up the data to fit, such as: the noticed energy range, how to
    handle the background, extracting parameter values, and plotting
    data.

    """

    datasets = None
    """Information (dataset identifier, annulus, name) about the loaded data."""

    def __init__(self):
        self.datasets = []
        self.obsids = set()
        self._reset_model_comps()

    def _reset_model_comps(self):
        """Clean out the model component information.

        Notes
        -----
        This does not delete the previous model components from
        Sherpa, but does remove them from the object.
        """

        # Generic model components in source model expression
        self.srcmodel_comps = []

        # All instantiated model components for shells
        self.model_comps = []

    def load_pha(self, specfile, annulus):
        """Load a pha file and add to the datasets for stacked analysis.

        It is required that datasets for all annuli are loaded before
        the source model is created (to ensure that components are
        created for each annulus).

        Parameters
        ----------
        specfile : str or sherpa.astro.data.DataPHA object
            If a string, the name of the file containing the source spectrum,
            which must be in PHA format (the data is expected to be extracted
            on the PI column). If a DataPHA object, then this is used (and
            is assumed to contain any needed background data).
        annulus : int
            The annulus number for the data.

        Returns
        -------
        dataid : int
            The Sherpa dataset identifier used for this spectrum.

        Examples
        --------

        Load the data for four annuli from the files 'ann1.pi' to 'ann4.pi'.

        >>> dep.load_pha('ann1.pi', 0)
        >>> dep.load_pha('ann2.pi', 1)
        >>> dep.load_pha('ann3.pi', 2)
        >>> dep.load_pha('ann4.pi', 3)

        Load in the PHA files into Sherpa DataPHA objects, and then use
        these objects:

        >>> s1 = ui.unpack_pha('src1.pi')
        >>> s2 = ui.unpack_pha('src2.pi')
        >>> s3 = ui.unpack_pha('src3.pi')
        >>> dep.load_pha(s1, 0)
        >>> dep.load_pha(s2, 1)
        >>> dep.load_pha(s3, 2)

        """

        dataid = len(self.datasets)

        # If the input has a counts attribute then assume a DataPHA
        # style object.
        #
        if hasattr(specfile, 'counts'):
            print('Using spectrum {} '.format(specfile.name) +
                  ' as dataset id {}'.format(dataid))
            ui.set_data(dataid, specfile)

        else:
            print('Loading spectrum file {} '.format(specfile) +
                  ' as dataset id {}'.format(dataid))
            ui.load_pha(dataid, specfile)

        data = ui.get_data(dataid)
        try:
            obsid = int(data.header['OBS_ID'])
        except (KeyError, TypeError, ValueError):
            obsid = 0

        dataset = dict(file=specfile,
                       obsid=obsid,
                       id=dataid,
                       annulus=annulus
                       )
        self.datasets.append(dataset)
        self.obsids.add(obsid)
        return dataid

    def _split_parname(self, par):
        """Return the model type and parameter name.

        The input value can either be a "generic" specifier, such as
        'xsphabs.nh' or specific to an annulus ('xspahbs_2.nh')

        Parameters
        ----------
        par : str
            The parameter name (e.g. 'xsapec.kt').

        Returns
        -------
        rtype : str, str
            The model and parameter name.

        Examples
        --------

        The foillowing will set `mname` to "xsapec" and `pname` to "kt":

        >>> mname, pname = dep._split_parname('xsapec.kt')

        """

        pos = par.rfind('.')
        if pos == -1:
            raise ValueError("Parameter name must match 'model_type.par_name'")
        return par[:pos], par[pos + 1:]

    def find_parval(self, parname):
        """Return the value of the first parameter matching the name.

        Parameters
        ----------
        parname : str
            The parameter name. The case is ignored in the match, and the
            first match is returned.

        Returns
        -------
        parval : float
            The parameter value

        Raises
        ------
        ValueError
            There is no match for the parameter.

        See Also
        --------
        find_norm, set_par

        Examples
        --------

        >>> kt = dep.find_parval('kt')

        """
        RE_parname = re.compile(parname + '$', re.IGNORECASE)
        for model_comp in self.model_comps:
            mc = model_comp['object']   # Get the model component object
            for par in mc.pars:
                if RE_parname.match(par.name):
                    return par.val

        raise ValueError('Parameter %s not found in any model component' % parname)

    def find_norm(self, shell):
        """Return the normalization value for the given shell.

        This is limited to XSPEC-style models, where the parameter is called
        "norm".

        Parameters
        ----------
        shell : int
            The shell number.

        Returns
        -------
        norm : float
            The normalization of the shell.

        Raises
        ------
        ValueError
            If there is not one `norm` parameter for the shell.

        See Also
        --------
        find_parval, set_par

        """
        norms = []
        for model_comp in self.model_comps:
            if model_comp['shell'] == shell:
                mc = model_comp['object']   # Get the model component object
                if hasattr(mc, 'norm'):      # Special attr of model component
                    norms.append(mc.norm.val)

        if len(norms) == 0:
            raise ValueError('Model for shell %d has no norm' % shell)
        elif len(norms) > 1:
            raise ValueError('Model for shell %d has multiple norms' % shell)

        return norms[0]

    def _get_n_datasets(self):
        """How many datasets are registered?

        This is not the same as the number of annuli.

        Returns
        -------
        ndata : int
            The number of datasets.
        """
        return len(self.datasets)

    n_datasets = property(_get_n_datasets)

    def _sherpa_cmd(self, func, *args):
        """Apply a function to each dataset.

        Parameters
        ----------
        func : function reference
            This function is called with each dataset identifier as
            the first argument, with the remaining arguments next.
        *args
            The additional arguments for the function.

        """
        for dataset in self.datasets:
            func(dataset['id'], *args)

    def subtract(self, *args):
        """Subtract the background from each dataset.

        See Also
        --------
        unsubtract

        Examples
        --------

        >>> dep.substract()

        """
        self._sherpa_cmd(ui.subtract, *args)

    def unsubtract(self, *args):
        """Un-subtract the background from each dataset.

        This can be useful when you want to compare the results to
        the "wstat" stat (a Poisson-based stat which includes the
        background data as a component and provides a goodness-of-fit
        estimate).

        See Also
        --------
        subtract

        Examples
        --------

        >>> dep.unsubstract()

        """
        self._sherpa_cmd(ui.unsubtract, *args)

    def group(self, *args):
        """Apply the grouping for each data set.

        See Also
        --------
        ungroup

        Examples
        --------

        >>> dep.group()

        """
        self._sherpa_cmd(ui.group, *args)

    def ungroup(self, *args):
        """Turn off the grouping for each data set.

        See Also
        --------
        group

        Examples
        --------

        >>> dep.ungroup()

        """
        self._sherpa_cmd(ui.ungroup, *args)

    def notice(self, *args):
        """Apply Sherpa notice command to each dataset.

        The filter is applied to each data set separately.

        See Also
        --------
        ignore

        Examples
        --------

        Restrict the analysis to those bins which fall in the range
        0.5 to 7.0 keV, where the limits are included in the noticed
        range. The first call to `notice` is used to clear any
        existing filter.

        >>> dep.notice(None, None)
        >>> dep.notice(0.5, 7.0)

        """
        self._sherpa_cmd(ui.notice_id, *args)

    def ignore(self, *args):
        """Apply Sherpa ignore command to each dataset.

        The filter is applied to each data set separately.

        See Also
        --------
        notice

        Examples
        --------

        Restrict the analysis to those bins which fall in the range
        0.5 to 7.0 keV, where the limits are not included in the
        noticed range. The call to `notice` is used to clear any
        existing filter.

        >>> dep.notice(None, None)
        >>> dep.ignore(None, 0.5)
        >>> dep.ignore(7.0, None)

        """
        self._sherpa_cmd(ui.ignore_id, *args)

    def _sherpa_par(self, func, par, msg, *args):
        """Apply the function to the given parameter.

        See thaw(), freeze(), set_par() and get_par() for examples.

        Parameters
        ----------
        func : function reference
            The function to call. The first argument is the name of the
            given parameter (including the model name), and then the
            optional arguments (`args`).
        par : str
            The parameter name to apply the function to (per shell). It
            must contain the model name (e.g. `xsapec.kt`).
        msg : str or None
            If not None, a format string that is printed each time `func`
            is called. It must contain one "%s" token, which is sent the
            full name of the parameter being processed.
        *args
            Additional arguments for `func`.

        Returns
        -------
        rvals : ndarray
            The return value from `func` for each parameter, ordered by
            shell.

        Raises
        ------
        ValueError
            If the parameter is not found.

        """

        model_type, parname = self._split_parname(par)

        vals = []                       # return values
        for model_comp in self.model_comps:
            if model_comp['type'] == model_type:
                fullparname = '%s.%s' % (model_comp['name'], parname)
                if msg is not None:
                    print(msg % fullparname)
                vals.append(func(fullparname, *args))

        if len(vals) == 0:
            raise ValueError("No parameter found matching {}".format(par))

        return vals

    def thaw(self, par):
        """Thaw the given parameter in each shell.

        Parameters
        ----------
        par : str
            The parameter name, specified as <model_type>.<par_name>

        See Also
        --------
        freeze, tie_par, untie_par

        Examples
        --------

        >>> dep.thaw('clus.abundanc')

        """
        self._sherpa_par(ui.thaw, par, 'Thawing %s')

    def freeze(self, par):
        """Freeze the given parameter in each shell.

        Parameters
        ----------
        par : str
            The parameter name, specified as <model_type>.<par_name>

        See Also
        --------
        thaw, tie_par, untie_par

        Examples
        --------

        >>> dep.freeze('clus.abundanc')

        """
        self._sherpa_par(ui.freeze, par, 'Freezing %s')

    def set_par(self, par, val):
        """Set the parameter value in each shell.

        Parameters
        ----------
        par : str
            The parameter name, specified as <model_type>.<par_name>
        val : float
            The parameter value.

        See Also
        --------
        get_par, tie_par

        Examples
        --------

        >>> dep.set_par('xsapec.abundanc', 0.25)

        """
        self._sherpa_par(ui.set_par, par, 'Setting %%s=%s' % str(val), val)

    def get_par(self, par):
        """Return the parameter value for each shell.

        Parameters
        ----------
        par : str
            The parameter name, specified as <model_type>.<par_name>

        Returns
        -------
        vals : ndarray
            The parameter values, in shell order.

        See Also
        --------
        find_parval, find_norm, set_par

        Examples
        --------

        >>> kts = dep.get_par('xsapec.kt')

        """

        pars = self._sherpa_par(ui.get_par, par, 'Getting %s')
        return numpy.array([x.val for x in pars])

    def _find_shell_par(self, par, shell):
        """Return the shell parameter.

        Parameters
        ----------
        par : str
            The parameter name, in the form `model name.par name`.
        shell : int
            The shell number

        Returns
        -------
        par : sherpa.models.parameter.Parameter

        """

        model_type, parname = self._split_parname(par)
        models = []
        for model_comp in self.model_comps:
            if model_comp['shell'] != shell or \
               model_comp['type'] != model_type:
                continue

            models.append(model_comp['object'])

        if len(models) == 0:
            raise ValueError('No parameter found matching ' +
                             '{} for shell {}'.format(par, shell))
        elif len(models) > 1:
            raise ValueError('Multiple parameters found matching ' +
                             '{} for shell {}'.format(par, shell))

        return getattr(models[0], parname)

    def tie_par(self, par, base, *others):
        """Tie parameters in one or more shells to the base shell.

        This is a limited form of the Sherpa ability to link parameters,
        since it sets the parameter in the other shells to the same
        value as the parameter in the base shell. More complex
        situations will require direct calls to `sherpa.astro.ui.link`.

        Parameters
        ----------
        par : str
            The parameter specifier, as <model_type>.<par_name>.
        base : int
            The base shell number.
        *others : scalar
            The shell, or shells to link to the base shell.

        See Also
        --------
        set_par, untie_par

        Examples
        --------

        Tie the temperature and abundance parameters in shell 9 to that
        in shell 8, so that any fits will set the shell 9 values to those
        used in shell 8 (so reducing the number of free parameters in the
        fit).

        >>> dep.tie_par('xsapec.kt', 8, 9)
        Tying xsapec_9.kT to xsapec_8.kT
        >>> dep.tie_par('xsapec.abundanc', 8, 9)
        Tying xsapec_9.Abundanc to xsapec_8.Abundanc

        Tie three annuli together:

        >>> dep.tie_par('xsapec.kt', 12, 13, 14)
        Tying xsapec_13.kT to xsapec_12.kT
        Tying xsapec_14.kT to xsapec_12.kT

        """

        bpar = self._find_shell_par(par, base)

        for other in others:
            opar = self._find_shell_par(par, other)
            print('Tying {} to {}'.format(opar.fullname, bpar.fullname))
            ui.link(opar, bpar)

    def untie_par(self, par, *others):
        """Remove the parameter tie/link in the shell.

        This is intended to remove links between shells created by `tie_par`,
        but will remove any links created by `sherpa.astro.ui.link`.

        Parameters
        ----------
        par : str
            The parameter specifier, as <model_type>.<par_name>.
        *others : scalar
            The shell, or shells to un-tie/unlink.

        See Also
        --------
        tie_par

        Notes
        -----
        It is safe to call on a parameter that is not tied or linked
        to another parameter.

        Examples
        --------

        Untie the abundance parameter in shell 9; that is, it is now free
        to vary independently in a fit.

        >>> dep.untie_par('xsapec.abundanc', 9)
        Untying xsapec_9.Abundanc

        Untie multiple annuli:

        >>> dep.untie_par('xsmekal.kt', 13, 14)
        Untying xsmekal_13.kT
        Untying xsmekal_14.kT

        """

        for other in others:
            opar = self._find_shell_par(par, other)
            if opar.link is None:
                continue

            print('Untying {}'.format(opar.fullname))
            ui.unlink(opar)

    def _sherpa_plot(self, func, *args, **kwargs):
        """Call the Sherpa plot ``func`` for each shell.

        Parameters
        ----------
        func : function reference
            The Sherpa plot function
        args
            The arguments for `func`
        kwargs
            Any keyword arguments for `func`. There are special
            keywords which are not passed on: `single_call_per_shell` is
            a boolean which indicates that the function is only
            called once per shell, and `add_shell_value` which indicates
            that the first argument is formatted to accept the shell
            value.

        Notes
        -----
        This method attempts to handle the differences when using
        ChIPS or Matplotlib as the Sherpa plotting backend, but has
        not been properly tested, so there may be issues.

        It is known not to work when the input function is plot_fit_delchi
        or plot_fit_resid, at least with the Matplotlib backend. It is
        unclear whether this is a problem here or the Sherpa
        matplotlib backend.

        """

        # Attempt to support multiple datasets per annulus.
        #
        dmap = [[] for _ in range(self.nshell)]
        for d in self.datasets:
            dmap[d['annulus']].append(d['id'])

        # Extract the arguments used here
        #
        special = {}
        for key in ['add_shell_value', 'single_call_per_shell']:
            val = False
            if key in kwargs:
                val = kwargs[key]
                del kwargs[key]

            special[key] = val

        nargs = len(args)

        for shell in range(self.nshell):
            if backend_name == 'pychips':
                window_id = 'Shell%d' % shell
                try:
                    pychips.add_window(['id', window_id])
                except RuntimeError:
                    pychips.set_current_window(window_id)

            elif backend_name == 'pylab':
                plt.figure(shell)

            new_args = args
            if nargs > 0:
                if special['add_shell_value']:
                    new_args = list(args)[:]
                    new_args[0] = new_args[0].format(shell)

            # For the moment assume that if the user supplied an
            # argument then we use it, whatever the various
            # settings are. This catches errors like
            #   dep.plot_fit('rstat')
            #
            if special['single_call_per_shell'] or nargs > 0:
                func(*new_args, **kwargs)
            else:
                # call it once per data set, assuming that the
                # overplot keyword is supported by this function
                # (could check this condition but leave that for now)
                #
                overplot = False
                for did in dmap[shell]:
                    kwargs['overplot'] = overplot
                    func(did, **kwargs)
                    overplot = True

    def print_window(self, *args, **kwargs):
        """Create a hardcopy version of each plot window.

        Parameters
        ----------
        args
            The arguments to be sent to the "create a hardcopy" routine
            (print_window for ChIPS and savefig for Matplotlib).
            The first argument, if given, is assumed to be the file name
            and so will have the shell number added to it.
        kwargs
            Keyword arguments for the call.

        Notes
        -----
        This is not guaranteed to work properly for Matplotlib.

        Examples
        --------

        Create hardcopy versions of the data plots, called "data0",
        "data1", ...

        >>> dep.plot_data()
        >>> dep.print_window('data')

        """

        if len(args) > 0:
            args = list(args)[:]
            try:
                args[0] += '{}'
            except TypeError:
                raise TypeError("First argument must be a string")

            kwargs['add_shell_value'] = True

        kwargs['single_call_per_shell'] = True
        if backend_name == 'pychips':
            plotfn = pychips.plot_window
        elif backend_name == 'pylab':
            plotfn = plt.savefig
        else:
            raise RuntimeError("Unrecognized plotting backend: {}".format(backend_name))

        self._sherpa_plot(plotfn, *args, **kwargs)

    def dummyfunc(self, *args, **kwargs):
        pass

    try:
        if backend_name == 'pychips':
            log_scale = _sherpa_plot_func_single(pychips.log_scale)
            linear_scale = _sherpa_plot_func_single(pychips.linear_scale)
        elif backend_name == 'pylab':
            # should be able to handle this
            log_scale = dummyfunc
            linear_scale = dummyfunc

    except NameError:
        # Allow doc generation to work without pychips
        log_scale = dummyfunc
        linear_scale = dummyfunc

    plot_fit = _sherpa_plot_func(ui.plot_fit)
    plot_arf = _sherpa_plot_func(ui.plot_arf)
    plot_bkg_fit = _sherpa_plot_func(ui.plot_bkg_fit)
    plot_bkg_ratio = _sherpa_plot_func(ui.plot_bkg_ratio)
    plot_chisqr = _sherpa_plot_func(ui.plot_chisqr)
    plot_fit_delchi = _sherpa_plot_func(ui.plot_fit_delchi)
    plot_psf = _sherpa_plot_func(ui.plot_psf)
    plot_bkg = _sherpa_plot_func(ui.plot_bkg)
    plot_bkg_fit_delchi = _sherpa_plot_func(ui.plot_bkg_fit_delchi)
    plot_bkg_resid = _sherpa_plot_func(ui.plot_bkg_resid)
    plot_data = _sherpa_plot_func(ui.plot_data)
    plot_fit_resid = _sherpa_plot_func(ui.plot_fit_resid)
    plot_ratio = _sherpa_plot_func(ui.plot_ratio)
    plot_bkg_chisqr = _sherpa_plot_func(ui.plot_bkg_chisqr)
    plot_bkg_fit_resid = _sherpa_plot_func(ui.plot_bkg_fit_resid)
    plot_bkg_source = _sherpa_plot_func(ui.plot_bkg_source)
    plot_delchi = _sherpa_plot_func(ui.plot_delchi)
    plot_model = _sherpa_plot_func(ui.plot_model)
    plot_resid = _sherpa_plot_func(ui.plot_resid)
    plot_bkg_delchi = _sherpa_plot_func(ui.plot_bkg_delchi)
    plot_bkg_model = _sherpa_plot_func(ui.plot_bkg_model)
    plot_bkg_unconvolved = _sherpa_plot_func(ui.plot_bkg_source)
    plot_bkg_source = _sherpa_plot_func(ui.plot_bkg_source)
    plot_fit = _sherpa_plot_func(ui.plot_fit)
    plot_order = _sherpa_plot_func(ui.plot_order)
    plot_source = _sherpa_plot_func(ui.plot_source)
