"""
Convert fields from Sherpa structures into arrays.

This is highly-specialised for the deproject code base. It is probably
excessive, but the overall separation was useful in trying to capture
the necessary functionality so it has been kept for now. That is, this
is an internal module and should not be used by user code.

Copyright 2018, Center for Astrophysics | Harvard & Smithsonian
Douglas Burke (dburke@cfa.harvard.edu)
"""

from collections import defaultdict, OrderedDict

import sherpa.astro.ui


class FieldStore:
    """Run and convert an object with keys to a set of values.

    Parameters
    ----------
    runfunc : function
        The function to call, which takes as argument the dataset ids.
    getfunc : function
        The function to call to return the structure. It takes no
        arguments.
    fields : sequence of str
        The field names to report, and the order to use. These fields
        must exist in the object returned by the getfunc.
    """

    def __init__(self, runfunc, getfunc, fields):

        if not callable(runfunc):
            raise ValueError("runfunc must be callable")

        if not callable(getfunc):
            raise ValueError("getfunc must be callable")

        self._run = runfunc
        self._get = getfunc

        # Copy the fields and ensure it is a list; I'm not sure why I
        # want a list, but a list we get.
        #
        self._fields = [f for f in fields]

    def run(self, annuli, components, *ids):
        """Run the function and return the values.

        Parameters
        ----------
        annuli : dict
            The key is the dataset identifier and the value the annulus
            number.
        components : dict
            The keys are the annulus number and the values are a list of
            tuples of (model components, parent model name).
        *ids : dataset identifiers
            The dataset identifiers to use. They must all exist as keys
            in the `annuli` map.

        Returns
        -------
        data : dict of OrderedDict values
            The keys are the annulus value, and the values are the
            shell results. The order of the keys in the results dict is
            defined by the fields list given to the object constructor.

        Notes
        -----
        There is only support for simple parameter links, where one parameter
        is set equal to another, and both model components are associated
        with the given set of dataset identifiers.

        """

        anns = sorted(list(set([annuli[i] for i in ids])))

        datamap = defaultdict(list)
        for i, ann in annuli.items():
            if ann not in anns:
                continue
            datamap[ann].append(i)

        self._run(*ids)
        indata = self._get()

        # Group the data by annulus. The assumption here is that each
        # dataset in the same annulus has the same set of model components.
        #
        # All the fields are copied over to each annulus except for
        #    datasets
        #    parameter values
        # which are restricted to the given annulus.
        #
        anndata = {}
        for ann in anns:
            dids = datamap[ann]

            # Sanity check
            for did in dids:
                if did not in indata.datasets:
                    raise RuntimeError("Dataset {} not found in the ouput!".format(did))

            # Create the non-parameter columns
            #
            store = OrderedDict()
            store['datasets'] = dids
            for field in self._fields:
                if field == 'datasets':
                    continue

                store[field] = getattr(indata, field)

            for mcpt, pname in components[ann]:
                for par in mcpt.pars:

                    # It would be nice to report these values too,
                    # for completeness, but leave for now.
                    #
                    # Note: linked parameters appear to be marked as
                    #       frozen...
                    #
                    # I think a better check may be needed (or send in
                    # the actual pars to report).
                    #
                    if par.frozen and par.link is None:
                        continue

                    # It is important to use the "base" name for the
                    # parameter - e.g 'xsapec.kT' - rather than the
                    # name for this shell (e.g. 'xsapec_0.kT').
                    #
                    name = "{}.{}".format(pname, par.name)
                    self._extract(indata, par, name, store)

            anndata[ann] = store

        return anndata

    def _extract(self, datavals, par, field, store):
        """Extract the parameter value from the results structure.

        Parameters
        ----------
        datavals
            The results structure (FitResults, ErrorEstResults).
        par : sherpa.models.parameter.Parameter instance
            The model parameter of interest. It must be related to the
            datasets in datavals.
        field : str
            The name of the item of interest.
        store : dict
            The results are added to this structure, using the field key
            (multiple values can be added, in which case the field is
            expected to be a suffix).

        """

        name = par.fullname

        # Just want the value. Could just take the value from the
        # parameter input, but use the value in datavals.
        #
        if self._store(datavals, name, field, store):
            return

        # Parameter not found, which should mean it is a link
        #
        if par.link is None:
            raise ValueError("Parameter {} not found in the fit results".format(name))

        # For now only support simple links - equality and not a more-complex
        # functional relationship. Assume that checking whether modelname
        # is not empty is a sufficient check
        #
        if par.link.modelname == '':
            raise ValueError("Unable to handle complex parameter links: {}".format(name))

        if self._store(datavals, par.link.fullname, field, store):
            return

        raise ValueError("Unable to find linked value for {}".format(name))

    def _store(self, datavals, name, field, store):
        """Add the values to the output structure.

        This functionality is provided by the sub-classes.

        Parameters
        ----------
        datavals
            The results structure (FitResults, ErrorEstResults).
        name : str
            The parameter name.
        field : str
            The name of the item of interest.
        store : dict
            The results are added to this structure, using the field key
            (multiple values can be added, in which case the field is
            expected to be a suffix).

        Returns
        -------
        flag : bool
            True if the parameter was found and added to store, False
            otherwise.

        """

        raise NotImplementedError()


class FitStore(FieldStore):
    """Run and extract the fields from a fit.

    The constructor has no arguments.
    """

    def __init__(self):
        super().__init__(sherpa.astro.ui.fit,
                         sherpa.astro.ui.get_fit_results,
                         ['datasets', 'succeeded',
                          'statval', 'dstatval',
                          'numpoints', 'dof', 'qval', 'rstat',
                          'nfev'])

    def _store(self, datavals, name, field, store):

        # Assume that par.fullname can only occur once in parnames.
        #
        for pname, pval in zip(datavals.parnames, datavals.parvals):
            if pname == name:
                store[field] = pval
                return True

        return False


class ErrorStore(FieldStore):
    """Run and extract the fields from a sherpa.fit.ErrorEstResults instance.

    Parameters
    ----------
    runfunc : function
        The function to call, which takes as argument the dataset ids.
    getfunc : function
        The function to call to return the structure. It takes no
        arguments.

    """

    def __init__(self, runfunc, getfunc):
        super().__init__(runfunc, getfunc,
                         ['datasets', 'sigma', 'percent', 'nfits'])

    def _store(self, datavals, name, field, store):

        # Assume that par.fullname can only occur once in parnames.
        #
        # QUS: do we need to handle limits or cases where the limit
        #      could not be accurately determined.
        #
        for pname, pval, plo, phi in zip(datavals.parnames,
                                         datavals.parvals,
                                         datavals.parmins,
                                         datavals.parmaxes):
            if pname == name:
                store[field] = pval
                store[field + "_lo"] = plo
                store[field + "_hi"] = phi
                return True

        return False


class CovarStore(ErrorStore):
    """Run and extract the covariance results.

    The constructor takes no arguments.
    """

    def __init__(self):
        super().__init__(sherpa.astro.ui.covar,
                         sherpa.astro.ui.get_covar_results)


class ConfStore(ErrorStore):
    """Run and extract the confidence results.

    The constructor takes no arguments.
    """

    def __init__(self):
        super().__init__(sherpa.astro.ui.conf,
                         sherpa.astro.ui.get_conf_results)
