Deproject
=========

``Deproject`` is a `Sherpa <https://sherpa.readthedocs.io/>`_ extension package
to facilitate deprojection of two-dimensional annular X-ray spectra to recover
the three-dimensional source properties.  For typical thermal models this would
include the radial temperature and density profiles. This basic method has been
used extensively for X-ray cluster analysis and is the basis for the XSPEC
model ``projct``.  The ``deproject`` module brings this functionality to
*Sherpa* as a Python module that is straightforward to use and understand.

The basic physical assumption of ``deproject`` is that the extended source
emissivity is constant and optically thin within spherical shells whose radii
correspond to the annuli used to extract the specta.  Given this assumption one
constructs a model for each annular spectrum that is a linear volume-weighted
combination of shell models.

Version 0.2 of ``deproject`` is limited to circular annuli.

Further documentation is available at https://deproject.readthedocs.io/

License
-------

The ``deproject`` module is released under the
`BSD 2-Clause license <https://choosealicense.com/licenses/bsd-2-clause/>`_,
available as the file ``LICENSE`` in the source distribution.

Requirements
------------

The installation assumes that you are installing ``deproject`` into
the `CIAO environment <http://cxc.harvard.edu/ciao/>`_ (CIAO 4.14 or
later), since this is the easiest way to get the XSPEC models along
with Sherpa. The `standalone Sherpa <https://sherpa.readthedocs.io/>`_
version can be used, but in this case you will need to `build Sherpa
with XSPEC support
<https://sherpa.readthedocs.io/en/latest/install.html#xspec>`_.

The following Python packages are required:

- sherpa >= 4.14.0
- `Astropy <http://www.astropy.org/>`_
- `SciPy <https://www.scipy.org/scipylib/>`_.

Installation
------------

For Standalone Sherpa or CIAO versions newer than CIAO 4.14,
you should be able to install ``deproject`` with the command::

  pip3 install deproject

The installation requires pip version 19 or higher, but as that
was released in early 2019 it should be available.

The `installation documentation
<https://deproject.readthedocs.io/installation.html>`_ describes how
to build a development version from the `GitHub repository
<https://github.com/sherpa-deprojcet/deproject>`_.

Example
-------

If you have a set of X-ray PHA spectra called src<n>.pi, where <n> is
an integer representing the annulus number, and the files contain the
``XFLT0001`` to ``XFLT0005`` header keywords used by the
`XSPEC projct model <https://asd.gsfc.nasa.gov/XSPECwiki/projct_model>`_,
then a
`Deproject object <https://deproject-test.readthedocs.io/en/latest/modules/api/deproject.deproject.Deproject.html#deproject.deproject.Deproject>`_
can be created using the
`deproject_from_xflt <https://deproject-test.readthedocs.io/en/ciao-411/modules/api/deproject.deproject.deproject_from_xflt.html>`_
helper routine with the commands::

  >>> from deproject import deproject_from_xflt
  >>> from astropy import units as u
  >>> dep = deproject_from_xflt('src*.pi', 0.492 * u.arcsec)

where, in this example, the ``XFLT0001`` and ``XFLT0002`` keywords,
which specify the inner and outer radii of the annulus, are in
ACIS pixels, and so need to be multiplied by 0.492 arcseconds to
convert to an angle (the second parameter).

This will automatically load the spectra into separate Sherpa datasets,
which *can* be fitted individually, but it is generally easier to use
the object returned by ``deproject_from_xflt``. For instance, the
following will set the data range to be fit for *each* spectra and ensure
that the background is subtracted before fitting::

  >>> dep.ignore(None, 0.5)
  >>> dep.ignore(7.0, None)
  >>> dep.subtract()

Sherpa functions are used to change the statistic and optimiser::

  >>> from sherpa.astro import ui
  >>> ui.set_stat('chi2xspecvar')
  >>> ui.set_method('levmar')

The data can be fit, and errors estimated for all the parameter, using
the onion-skin deprojection approach, with the following commands::

  >>> onion = dep.fit()
  >>> errs = dep.conf()

The return value includes the density (and errors, if appropriate), as
an `Astropy Quantity <http://docs.astropy.org/en/stable/units/>`_.

::

  >>> print(onion['density'])
  print(onion['density'])
        density
        1 / cm3
  --------------------
    0.1100953546292787
   0.07736622021374819
   0.04164827967805805
   0.03630168106524076
  0.025221797991301052
  0.021845331641349316
                   ...
  0.012396857131392835
   0.01336640115325031
  0.012303975980575187
  0.013631563529090736
  0.013996131292837352
  0.010843683594144967
  0.023067220584935984
  Length = 20 rows

The `on-line documentation <https://deproject.readthedocs.io/>`_
contains more information, including creating the ``Deproject`` object
directly (without the need for the ``XFLTxxxx`` keywords).
