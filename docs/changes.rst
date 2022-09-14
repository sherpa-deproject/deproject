
.. include:: references.rst

*******	     
Changes
*******

.. _changes_021:

Version 0.2.1
=============

The build system has been updated to use ``pip`` rather than
``setup.py``. Please report any problems on the
`issues page <https://github.com/sherpa-deproject/deproject/issues/>`_.

Some plotting code related to ChIPS has been removed, and code has
been updated to reflect recent changes in Sherpa. The minimum
supported version of Sherpa has been bumped to 4.14.0.

.. _changes_020:

Version 0.2.0
=============

Overview
--------

The code has been updated to run with Python 3 and can now be
installed from `PyPI`_. Documentation has been moved to
`Read The Docs <https://deproject.readthedocs.io/>`_.

The deproject module now requires Astropy_, which can be
`installed within CIAO 4.11 <http://cxc.harvard.edu/ciao/scripting/index.html#install>`_. The three main areas where Astropy functionality is
used are:

- the use of `Astropy Quantity <http://docs.astropy.org/en/stable/units/>`_
  values (both for arguments to methods and returned values);
- `Astropy Data Tables <http://docs.astropy.org/en/stable/table/>`_
  are used to return tabular data;
- and Cosmology calculations now use the `Astropy cosmology
  module <http://docs.astropy.org/en/stable/cosmology/>`_ rather than
  the `cosmocalc` module.

The :py:func:`~deproject.deproject.deproject_from_xflt` helper function
has been introduced, which uses the ``XFLT0001`` to ``XFLT0005``
keywords in the input files to determine the annulus parameters (radii
and covering angle). The covering angle (:math:`\theta`) can now vary per
annulus.

Error values can now be generated using the onion-peeling approach,
for the confidence and covariance methods, and the values are returned
as an Astropy Table. Parameter values can now be tied together (to
combine annuli to try and avoid "ringing"). There is improved support
for accessing and plotting values.

Details
-------

The code has been re-arranged into the ``deproject`` package, which
means that you really should say ``from deproject.deproject import
Deproject``, but the ``deproject`` module re-exports
:py:mod:`deproject.deproject` so that existing scripts still work, and
you do not not have to type in the same word multiple times! The
package has been updated so that it is available on `PyPI`_.

The scaling between shells (calculated from the intersection between
spheres and cylinders) was limited to 5 decimal places, which could
cause problems with certain choices of annuli (such as an annulus
making no contribution to interior annuli). This restriction has been
removed.

Added support for per-annulus ``theta`` values (that is, each annulus
can have a different opening angle). The ``radii``, ``theta``, and
``angdist`` parameters to :py:class:`~deproject.deproject.Deproject`
all now require values that is an
`Astropy quantity <http://docs.astropy.org/en/stable/units/>`_ rather
than a dimensionless value.

Added the :py:func:`~deproject.deproject.deproject_from_xflt` helper
function, which creates a :py:class:`~deproject.deproject.Deproject`
instance from PHA files which contain the XSPEC XFLT0001 to XFLT0005
keywords (as used by the `projct`_ model), rather than specifying the
values from the command line. The routine will error out if the
keywords indicate elliptical annuli, and the default is to assume the
radii are in arcseconds, but a scaling factor can be given if the
radii are in some other units (such as pixels).

Added the :py:meth:`~deproject.deproject.Deproject.guess` method to do
an initial fit to each annulus, following the approach suggested in
the `XSPEC`_ documentation for `projct`_, by just fitting the
individual (not de-projected) models to each annulus. This can help
speed up the deproject fit -
:py:meth:`~deproject.deproject.Deproject.fit` - as well as help avoid
the fit getting stuck in a local minimum.

Added :py:meth:`~deproject.deproject.Deproject.covar` and
:py:meth:`~deproject.deproject.Deproject.conf` methods that estimate
errors - using the covariance and confidence methods respectively -
using the onion-skin model (i.e. the errors on the outer annuli are
evaluated, then this component is frozen and the errors on the next
annulus are evaluated).

The :py:meth:`~deproject.deproject.Deproject.fit`,
:py:meth:`~deproject.deproject.Deproject.conf`, and
:py:meth:`~deproject.deproject.Deproject.covar` methods now all return
Astropy Tables containing the results per annulus. These values can
also be retrieved with the
:py:meth:`~deproject.deproject.Deproject.get_fit_results`,
:py:meth:`~deproject.deproject.Deproject.get_conf_results`, or
:py:meth:`~deproject.deproject.Deproject.get_covar_results` methods. A
number of columns (radii and density) are returned as Astropy
quantities.

The ``cosmocalc`` module has been removed and the `Astropy cosmology
module <http://docs.astropy.org/en/stable/cosmology/>`_ is used
instead. This is only used if the angular-diameter distance to the
source is calculated rather than explicitly given. The default
cosmology is now set to `Planck15
<http://docs.astropy.org/en/stable/cosmology/index.html#built-in-cosmologies>`_.

Values, as a function of radius, can be plotted with a number of new
methods: :py:meth:`~deproject.deproject.Deproject.fit_plot`,
:py:meth:`~deproject.deproject.Deproject.conf_plot`, and
:py:meth:`~deproject.deproject.Deproject.covar_plot` display the last
fit results (with the last two including error estimates), and the
:py:meth:`~deproject.deproject.Deproject.par_plot` and
:py:meth:`~deproject.deproject.Deproject.density_plot` methods show
the current values. These support a number of options, including
switching between angular and physical distances for the radii.

The :py:meth:`~deproject.deproject.Deproject.get_shells` method has
been added to make it easy to see which annuli are combined together,
and the :py:meth:`~deproject.deproject.Deproject.get_radii` method to
find the radii of the annuli (in a range of units).

Added the :py:meth:`~deproject.deproject.SpecStack.tie_par` and
:py:meth:`~deproject.deproject.SpecStack.untie_par` methods to make it
easy to tie (or untie) parameters in neighbouring annuli. The
onion-skin approach - used when fitting or running an error analysis -
recognizes annuli that are tied together and fits these
simultaneously, rather than individually.

The :py:meth:`~deproject.deproject.Deproject.set_source` method can
now be called multiple times (previously it would lead to an error).

Added error checking for several routines, such as
:py:meth:`~deproject.deproject.SpecStack.thaw` when given an
unknown parameter name.

Updated to support Python 3.5 and to have better support when the
``pylab`` backend is selected. Support for the ChIPS backend is
limited. A basic test suite has been added.

Version 0.1.0
=============

Initial version.
