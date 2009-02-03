import re
import numpy
import sherpa.astro.ui as SherpaUI
import pycrates
import pychips

def sherpa_plot_func(func):
    def _sherpa_plot_func(self, *args, **kwargs):
        self._sherpa_plot(func, *args, **kwargs)
    return _sherpa_plot_func

class SpecStack(object):
    """
    Manipulate a stack of spectra in sherpa.
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
        :param annulus: annulus number corresponding to the ``radii`` list (starts at 0)
        """
        if annulus < 0 or annulus >= self.nshell:
            raise ValueError('annulus=%d must be between 0 to %d inclusive' %
                             (annulus, self.nshell-1))

        dataid = len(self.datasets)
        print 'Loading spectrum file %s as dataset id %d' % (specfile, dataid)
        SherpaUI.load_pha(dataid, specfile)

        try:
            obsid = int(pycrates.read_file(specfile).get_key_value('OBS_ID'))
        except TypeError:
            obsid = 0
        dataset = dict(file=specfile,
                       annulus=annulus,
                       obsid=obsid,
                       id=dataid,
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

    def _sherpa_cmd(self, func, *args):
        for dataset in self.datasets:
            func(dataset['id'], *args)

    def subtract(self, *args):
        """Subtract background from each dataset"""
        self._sherpa_cmd(SherpaUI.subtract, *args)

    def notice(self, *args):
        """Wrapper around sherpa notice command."""
        self._sherpa_cmd(SherpaUI.notice_id, *args)

    def ignore(self, *args):
        """Wrapper around sherpa ignore command."""
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
        self._sherpa_par(SherpaUI.thaw, par, 'Thawing %s')

    def freeze(self, par):
        self._sherpa_par(SherpaUI.freeze, par, 'Freezing %s')

    def set_par(self, par, val):
        self._sherpa_par(SherpaUI.set_par, par, 'Setting %%s=%s' % str(val), val)

    def get_par(self, par):
        pars = self._sherpa_par(SherpaUI.get_par, par, 'Getting %s')
        return numpy.array([x.val for x in pars])

    def _sherpa_plot(self, func, *args, **kwargs):
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

    log_scale = sherpa_plot_func(pychips.log_scale)
    linear_scale = sherpa_plot_func(pychips.linear_scale)
    def print_window(self, *args, **kwargs):
        if len(args) > 0:
            args = tuple([args[0] + '%d']) + args[1:]
        self._sherpa_plot(pychips.plot_window, *args, **kwargs)

    plot_fit = sherpa_plot_func(SherpaUI.plot_fit)
    plot_arf = sherpa_plot_func(SherpaUI.plot_arf)
    plot_bkg_fit = sherpa_plot_func(SherpaUI.plot_bkg_fit)
    plot_bkg_ratio = sherpa_plot_func(SherpaUI.plot_bkg_ratio)
    plot_chisqr = sherpa_plot_func(SherpaUI.plot_chisqr)
    plot_fit_delchi = sherpa_plot_func(SherpaUI.plot_fit_delchi)
    plot_psf = sherpa_plot_func(SherpaUI.plot_psf)
    plot_bkg = sherpa_plot_func(SherpaUI.plot_bkg)
    plot_bkg_fit_delchi = sherpa_plot_func(SherpaUI.plot_bkg_fit_delchi)
    plot_bkg_resid = sherpa_plot_func(SherpaUI.plot_bkg_resid)
    plot_data = sherpa_plot_func(SherpaUI.plot_data)
    plot_fit_resid = sherpa_plot_func(SherpaUI.plot_fit_resid)
    plot_ratio = sherpa_plot_func(SherpaUI.plot_ratio)
    plot_bkg_chisqr = sherpa_plot_func(SherpaUI.plot_bkg_chisqr)
    plot_bkg_fit_resid = sherpa_plot_func(SherpaUI.plot_bkg_fit_resid)
    plot_bkg_source = sherpa_plot_func(SherpaUI.plot_bkg_source)
    plot_delchi = sherpa_plot_func(SherpaUI.plot_delchi)
    plot_model = sherpa_plot_func(SherpaUI.plot_model)
    plot_resid = sherpa_plot_func(SherpaUI.plot_resid)
    plot_bkg_delchi = sherpa_plot_func(SherpaUI.plot_bkg_delchi)
    plot_bkg_model = sherpa_plot_func(SherpaUI.plot_bkg_model)
    plot_bkg_unconvolved = sherpa_plot_func(SherpaUI.plot_bkg_unconvolved)
    plot_fit = sherpa_plot_func(SherpaUI.plot_fit)
    plot_order = sherpa_plot_func(SherpaUI.plot_order)
    plot_source = sherpa_plot_func(SherpaUI.plot_source)
