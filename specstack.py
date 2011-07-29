"""
Manipulate a stack of spectra in Sherpa.

The methods in the SpecStack class provide a way to automatically apply
familiar Sherpa commands such as `set_par`_ or `freeze`_ or `plot_fit`_
to a stack of PHA spectra.  This simplifies simultaneous fitting of
multiple spectra.

Note that the :mod:`specstack` module is currently distributed within with the
:mod:`deproject` package.  `Specstack` is not yet fully documented or tested
outside the context of `deproject`.  

:Copyright: Smithsonian Astrophysical Observatory (2009)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""
import re
import numpy
import sherpa.astro.ui as SherpaUI
try:
    import pycrates
    import pychips
except ImportError:
    # Allow doc generation to work without pychips and pycrates
    print "ERROR: could not import pycrates and pychips"

def _sherpa_plot_func(func):
    def _sherpa_plot_func_(self, *args, **kwargs):
        self._sherpa_plot(func, *args, **kwargs)
    return _sherpa_plot_func_

class SpecStack(object):
    """
    Manipulate a stack of spectra in Sherpa.
    """
    def __init__(self):
        self.datasets = []
        self.model_comps = []           # Model components
        self.obsids = set()
        self.srcmodel_comps = []        # Generic model components in source model expression
        self.model_comps = []           # All instantiated model components for shells

    def load_pha(self, specfile, annulus):
        """
        Load a pha file and add to the datasets for stacked analysis.

        :param specfile: extracted source PHA/PI spectrum file
        :param annulus: annulus for spectrum file
        """
        dataid = len(self.datasets)
        print 'Loading spectrum file %s as dataset id %d' % (specfile, dataid)
        SherpaUI.load_pha(dataid, specfile)

        try:
            obsid = int(pycrates.read_file(specfile).get_key_value('OBS_ID'))
        except (TypeError, ValueError):
            obsid = 0
        dataset = dict(file=specfile,
                       obsid=obsid,
                       id=dataid,
                       annulus=annulus
                       )
        self.datasets.append(dataset)
        self.obsids.add(obsid)

    def find_parval(self, parname):
        """
        Find the value of the first parameter with the given ``parname``.  Ignore
        case when matching.

        :param parname: parameter name
        :rtype: parameter value 
        """
        RE_parname = re.compile(parname + '$', re.IGNORECASE)
        for model_comp in self.model_comps:
            mc = model_comp['object']   # Get the model component object
            for par in mc.pars:
                if RE_parname.match(par.name):
                    return par.val

        raise ValueError('Parameter %s not found in any model component' % shell)

    def find_norm(self, shell):
        """
        Find the normalization value for the ``shell`` model.  Check for multiple
        or missing norms.

        :param shell: shell index
        :rtype: norm value
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
        """
        Return the number of datasets
        :rtype: int"""
        return len(self.datasets)

    n_datasets = property(_get_n_datasets)

    def _sherpa_cmd(self, func, *args):
        """
        Apply an arbitrary Sherpa function to each of the datasets.
        :rtype: None
        """
        for dataset in self.datasets:
            func(dataset['id'], *args)

    def subtract(self, *args):
        """Subtract background from each dataset"""
        self._sherpa_cmd(SherpaUI.subtract, *args)

    def notice(self, *args):
        """Apply Sherpa notice command to each dataset."""
        self._sherpa_cmd(SherpaUI.notice_id, *args)

    def ignore(self, *args):
        """Apply Sherpa ignore command to each dataset."""
        self._sherpa_cmd(SherpaUI.ignore_id, *args)

    def _sherpa_par(self, func, par, msg, *args):
        """Apply ``func(*args)`` to all shell model component parameters named ``par``.

        See thaw(), freeze(), set_par() and get_par() for examples.

        :param func: Sherpa function that takes a full parameter name specification and
                     optional args, e.g. set_par() used as set_par('xsmekal_7.kt', 2.0)
        :param par: Model type and param name as in source model expression e.g. 'xsmekal.kt'
        :param msg: Format string to indicate action.
        :param *args: Optional function arguments

        :rtype: numpy array of function return values ordered by shell
        """
        model_type, parname = par.split('.')

        vals = []                       # return values
        for model_comp in self.model_comps:
            if model_comp['type'] == model_type:
                fullparname = '%s.%s' % (model_comp['name'], parname)
                if msg is not None:
                    print msg % fullparname
                vals.append(func(fullparname, *args))

        return vals

    def thaw(self, par):
        """Apply thaw command to specified parameter for each dataset.

        :param par: parameter specifier in format <model_type>.<par_name>
        :rtype: None
        """
        self._sherpa_par(SherpaUI.thaw, par, 'Thawing %s')

    def freeze(self, par):
        """Apply freeze command to specified parameter for each dataset.

        :param par: parameter specifier in format <model_type>.<par_name>
        :rtype: None
        """
        self._sherpa_par(SherpaUI.freeze, par, 'Freezing %s')

    def set_par(self, par, val):
        """Set parameter value for each dataset.

        :param par: parameter specifier in format <model_type>.<par_name>
        :param val: parameter value
        :rtype: None
        """
        self._sherpa_par(SherpaUI.set_par, par, 'Setting %%s=%s' % str(val), val)

    def get_par(self, par):
        """Get array of parameter values for datasets.

        :param par: parameter specifier in format <model_type>.<par_name>
        :param val: parameter value
        :rtype: numpy array of parameter value ordered by dataset
        """
        pars = self._sherpa_par(SherpaUI.get_par, par, 'Getting %s')
        return numpy.array([x.val for x in pars])

    def _sherpa_plot(self, func, *args, **kwargs):
        """Call Sherpa plot ``func`` for each dataset.

        :param func: Sherpa plot function
        :param args: plot function list arguments
        :param kwargs: plot function named (keyword) arguments
        :rtype: None
        """
        for shell in range(self.nshell):
            window_id = 'Shell%d' % shell
            try:
                pychips.add_window(['id', window_id])
            except RuntimeError:
                pass  # already exists

            new_args = args
            if len(args) > 0:
                # Try to format first arg
                try:
                    new_args = tuple([args[0] % shell]) + args[1:]
                except TypeError:
                    pass

            pychips.set_current_window(window_id)
            func(*new_args, **kwargs)

    def print_window(self, *args, **kwargs):
        """Print window for each dataset.

        :param args: list arguments to pass to print_window
        :param kwargs: named (keyword) arguments to pass to print_window
        :rtype: None
        """
        if len(args) > 0:
            args = tuple([args[0] + '%d']) + args[1:]
        self._sherpa_plot(pychips.plot_window, *args, **kwargs)

    def dummyfunc(self, *args, **kwargs):
        pass

    try:
        log_scale = _sherpa_plot_func(pychips.log_scale)
        linear_scale = _sherpa_plot_func(pychips.linear_scale)
    except NameError:
        # Allow doc generation to work without pychips
        log_scale = dummyfunc
        linear_scale = dummyfunc

    plot_fit = _sherpa_plot_func(SherpaUI.plot_fit)
    plot_arf = _sherpa_plot_func(SherpaUI.plot_arf)
    plot_bkg_fit = _sherpa_plot_func(SherpaUI.plot_bkg_fit)
    plot_bkg_ratio = _sherpa_plot_func(SherpaUI.plot_bkg_ratio)
    plot_chisqr = _sherpa_plot_func(SherpaUI.plot_chisqr)
    plot_fit_delchi = _sherpa_plot_func(SherpaUI.plot_fit_delchi)
    plot_psf = _sherpa_plot_func(SherpaUI.plot_psf)
    plot_bkg = _sherpa_plot_func(SherpaUI.plot_bkg)
    plot_bkg_fit_delchi = _sherpa_plot_func(SherpaUI.plot_bkg_fit_delchi)
    plot_bkg_resid = _sherpa_plot_func(SherpaUI.plot_bkg_resid)
    plot_data = _sherpa_plot_func(SherpaUI.plot_data)
    plot_fit_resid = _sherpa_plot_func(SherpaUI.plot_fit_resid)
    plot_ratio = _sherpa_plot_func(SherpaUI.plot_ratio)
    plot_bkg_chisqr = _sherpa_plot_func(SherpaUI.plot_bkg_chisqr)
    plot_bkg_fit_resid = _sherpa_plot_func(SherpaUI.plot_bkg_fit_resid)
    plot_bkg_source = _sherpa_plot_func(SherpaUI.plot_bkg_source)
    plot_delchi = _sherpa_plot_func(SherpaUI.plot_delchi)
    plot_model = _sherpa_plot_func(SherpaUI.plot_model)
    plot_resid = _sherpa_plot_func(SherpaUI.plot_resid)
    plot_bkg_delchi = _sherpa_plot_func(SherpaUI.plot_bkg_delchi)
    plot_bkg_model = _sherpa_plot_func(SherpaUI.plot_bkg_model)
    plot_bkg_unconvolved = _sherpa_plot_func(SherpaUI.plot_bkg_source)
    plot_bkg_source = _sherpa_plot_func(SherpaUI.plot_bkg_source)
    plot_fit = _sherpa_plot_func(SherpaUI.plot_fit)
    plot_order = _sherpa_plot_func(SherpaUI.plot_order)
    plot_source = _sherpa_plot_func(SherpaUI.plot_source)

