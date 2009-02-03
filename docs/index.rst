.. deproject documentation master file, created by sphinx-quickstart on Sat Jan 31 15:06:12 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _`CIAO`: http://cxc.harvard.edu/ciao/
.. _`sherpa`: http://cxc.harvard.edu/sherpa/
.. _`XSPEC`: http://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/
.. _`projct`: https://astrophysics.gsfc.nasa.gov/XSPECwiki/projct_model
.. _`freeze`: http://cxc.harvard.edu/sherpa/ahelp/freeze.py.html
.. _`set_par`: http://cxc.harvard.edu/sherpa/ahelp/set_par.py.html
.. _`set_source`: http://cxc.harvard.edu/sherpa/ahelp/set_source.py.html
.. _`ignore`: http://cxc.harvard.edu/sherpa/ahelp/ignore.py.html


deproject
====================

Overview
----------

:mod:`Deproject` is a `CIAO`_ `Sherpa`_ extension package to facilitate
deprojection of two-dimensional annular X-ray spectra to recover the
three-dimensional source properties.  This basic method (refs) has been
used extensively for X-ray cluster analysis and is the basis for the `XSPEC`_
model `projct`_.  The :mod:`deproject` module brings this functionality to
*Sherpa* as a Python module that is straightforward to use and understand.

The :mod:`deproject` module uses :mod:`specstack` to allow for manipulation of
a stack of related input datasets and their models.  Most of the functions are
just like ordinary *Sherpa* commands (e.g. `set_par`_, `set_source`_,
`ignore`_) but operate on a stack of spectra.


Download
---------
Here is how to download

Installation
-------------
Here is how to install

Concepts
-----------

The basic physical assumption of :mod:`deproject` is that the extended source
emissivity is constant and optically thin within spherical shells whose radii
correspond to the annuli used to extract the specta.  Given this assumption one
constructs a model for each annular spectrum that is a linear volume-weighted
combination of shell models.  The geometry is illustrated in the figure below
(which would be rotated about the line to the observer in three-dimensions):

.. image:: geometry.png

Model creation
^^^^^^^^^^^^^^^
It is assumed that prior to starting :mod:`deproject` the user has extracted
source and background spectra for each annulus.  By convention the annulus
numbering starts from the inner radius at 0 and corresponds to the dataset
``id`` used within *Sherpa*.  It is not required that the annuli include the
center but they must be contiguous between the inner and outer radii.

Given a spectral model ``M[s]`` for each shell ``s``, the source model for
dataset ``a`` (i.e. annulus ``a``) is given by the sum over ``s >= a`` of
``vol_norm[s,a] * M[s]`` (normalized volume * shell model).  The image above
shows shell 3 in blue and annulus 2 in red.  The intersection of (purple) has a
physical volume defined as ``vol_norm[3,2] * v_sphere`` where ``v_sphere`` is the
volume of the sphere enclosing the outer shell.

The bookkeeping required to create all the source models is handled by the
:mod:`deproject` module.

Fitting
^^^^^^^^^^^
Once the composite source models for each dataset are created the fit analysis
can begin.  Since the parameter space is typically large the usual procedure is
to initally fit using the "onion-peeling" method:

 - First fit the outside shell model using the outer annulus spectrum
 - Freeze the model parameters for the outside shell
 - Fit the next inward annulus / shell and freeze those parameters
 - Repeat until all datasets have been fit and all shell parameters determined.

From this point the user may choose to do a simultanenous fit of the shell
models, possibly freezing some parameters as needed.  This process is made
manageable with the :mod:`specstack` methods that apply normal *Sherpa*
commands like `freeze`_ or `set_par`_ to a stack of spectral datasets.

Densities
^^^^^^^^^^^

Physical densities (cm^-3) for each shell can be calculated with
:mod:`deproject` assuming the source model is based on a thermal model with the
"standard" normalization (from the `XSPEC`_ documentation):

.. image:: thermal_norm.png

Inverting this equation and assuming a constant ratio of N_H to electrons::

 n_e = sqrt(norm * 4*pi * DA^2 * 1e14 * (1+z)^2 / volume * ne_nh_ratio))

 norm        = model normalization from Sherpa fit
 DA          = angular size distance (cm)
 volume      = volume (cm^3)
 ne_nh_ratio = 1.18

Recall that the model components for each volume element (intersection of the
annular cylinder ``a`` with the spherical shell ``s``) are multiplied by a volume
normalization::

 vol_norm[s,a] = v[s,a] / v_sphere
 v_sphere = volume of sphere enclosing outer annulus

With this convention the ``volume`` used above in calculating the electron
density for each shell is always ``v_sphere``.

Example: M87
------------------
The first thing we need to do is tell *Sherpa* about the Deproject class and
set a couple of constants::

  from deproject import Deproject

  redshift = 0.004233                     # M87 redshift
  arcsec_per_pixel = 0.492                # ACIS plate scale

Now we actually initialize the :class:`Deproject` object ``dep`` by passing a
`numpy`_ array of the the annular radii in arcsec.  The `numpy.arange`_ method
here returns an array from 30 to 640 in steps of 30.  These values were in
pixels in the original spectral extraction so we convert to arcsec.  (Note the
convenient vector multiplication that is possible with `numpy`_.)  

In this particular analysis the spectra were extracted from a 75 degree sector
of the annuli, hence ``theta=75`` in  the object initialization.  For the
default case of full 360 degree annuli this is not needed.::

  radii = numpy.arange(30., 640., 30) * arcsec_per_pixel
  dep = Deproject(radii, theta=75)

  # Load datasets for each annulus
  for ann in range(dep.nshell):
      dep.load_pha('m87/r%dgrspec.pha' % (ann+1), annulus=ann)

  # Subtract background
  dep.subtract()

  # Set source model and ignore specified energy ranges
  dep.set_source('xswabs*xsmekal')
  dep.ignore(None, 0.5)
  dep.ignore(1.8, 2.2)
  dep.ignore(7, None)

  # Specify Galactic absorption
  dep.set_par('xswabs.nh', 0.0255)
  dep.freeze("xswabs.nh")

  # Initialize abundance to 0.5 and thaw
  dep.set_par('xsmekal.abundanc', 0.5)
  dep.thaw('xsmekal.abundanc')

  # Set redshift
  dep.set_par('xsmekal.redshift', redshift)

  # Do the initial onion-peeling fit
  set_method("levmar")                    # Levenberg-Marquardt optimization method
  set_stat("chi2gehrels")                 # Gehrels Chi^2 fit statistic
  dep.fit()


.. image:: m87_density.png

Module documentation
====================

.. toctree::
   :maxdepth: 2

   deproject
   specstack
   cosmocalc
