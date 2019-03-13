
.. include:: references.rst

Deproject
=========

:mod:`Deproject` is a `CIAO`_ `Sherpa`_ extension package to facilitate
deprojection of two-dimensional circular annular X-ray spectra to recover the
three-dimensional source properties.  For typical thermal models this would
include the radial temperature and density profiles. This basic method 
has been used extensively for X-ray cluster analysis and is the basis for the
`XSPEC`_ model `projct`_.  The :mod:`deproject` module brings this
functionality to *Sherpa* as a Python module that is straightforward to use and
understand.

The module can also be used with the `standalone Sherpa`_ release, but as
it requires support for `XSPEC`_ models - for the thermal plasma emission
codes - the documentation will focus on using Sherpa with a `CIAO`_
environment.

.. _`standalone Sherpa`: https://sherpa.readthedocs.io/

.. toctree::
   :maxdepth: 2
   :caption: Introduction

   overview
   installation
   changes
   todo
   
.. toctree::
   :maxdepth: 2
   :caption: Examples

   examples/m87
   examples/tie
   examples/3c186
   
.. toctree::
   :maxdepth: 2
   :caption: Module documentation
      
   modules/deproject
   modules/specstack

