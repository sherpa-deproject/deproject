Deproject
=========

``Deproject`` is a CIAO Sherpa extension package to facilitate
deprojection of two-dimensional annular X-ray spectra to recover the
three-dimensional source properties.  For typical thermal models this would
include the radial temperature and density profiles. This basic method 
has been used extensively for X-ray cluster analysis and is the basis for the
XSPEC model ``projct``.  The ``deproject`` module brings this
functionality to *Sherpa* as a Python module that is straightforward to use and
understand.

The basic physical assumption of ``deproject`` is that the extended source
emissivity is constant and optically thin within spherical shells whose radii
correspond to the annuli used to extract the specta.  Given this assumption one
constructs a model for each annular spectrum that is a linear volume-weighted
combination of shell models.

License
-------

The ``deproject`` module is released under the
`BSD 2-Clause license <https://choosealicense.com/licenses/bsd-2-clause/>`_,
available as the file ``LICENSE`` in the source distribution.

