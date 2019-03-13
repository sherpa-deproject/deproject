
.. include:: references.rst

Overview
========

The :mod:`deproject` module uses :mod:`specstack` to allow for manipulation of
a stack of related input datasets and their models.  Most of the functions
resemble ordinary *Sherpa* commands (e.g. `set_par`_, `set_source`_, `ignore`_)
but operate on a stack of spectra.

The basic physical assumption of :mod:`deproject` is that the extended source
emissivity is constant and optically thin within spherical shells whose radii
correspond to the annuli used to extract the specta.  Given this assumption one
constructs a model for each annular spectrum that is a linear volume-weighted
combination of shell models.  The geometry is illustrated in the figure below
(which would be rotated about the line to the observer in three-dimensions):

.. image:: geometry.png

Model creation
--------------
It is assumed that prior to starting :mod:`deproject` the user has extracted
source and background spectra for each annulus.  By convention the annulus
numbering starts from the inner radius at 0 and corresponds to the dataset
``id`` used within *Sherpa*.  It is not required that the annuli include the
center but they must be contiguous between the inner and outer radii.

Given a spectral model ``M[s]`` for each shell ``s``, the source model for
dataset ``a`` (i.e. annulus ``a``) is given by the sum over ``s >= a`` of
``vol_norm[s,a] * M[s]`` (normalized volume * shell model).  The image above
shows shell 3 in blue and annulus 2 in red.  The intersection of (purple) has a
physical volume defined as ``vol_norm[3,2] * v_sphere`` where ``v_sphere``
is the volume of the sphere enclosing the outer shell (as
:ref:`discussed below <densities>`).

The bookkeeping required to create all the source models is handled by the
:mod:`deproject` module.

Fitting
-------
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

.. _densities:
   
Densities
---------
Physical densities (:math:`{\rm cm}^{-3}`) for each shell can be calculated with
:mod:`deproject` assuming the source model is based on a thermal model with the
"standard" normalization (from the `XSPEC`_ documentation):

.. math::

   {\rm norm} = \frac{10^{-14}}{4\pi [D_A (1+z)]^2} \int n_e n_H dV

where :math:`D_A` is the angular size distance to the source (in cm),
:math:`z` is the source redshift, and :math:`n_e` and :math:`n_H` are the
electron and Hydrogen densities (in :math:`{\rm cm}^{-3}`).

Inverting this equation and assuming constant values for :math:`n_e`
and :math:`n_H` gives

.. math::
   
   n_e = \sqrt{\frac{4\pi \, {\rm norm} \, [D_A (1+z)]^2 \, 10^{14}}{\alpha V}}

where :math:`V` is the volume, :math:`n_H = \alpha n_e`, and
:math:`\alpha` is taken to be 1.18 (and this ratio is
constant within the source).

Recall that the model components for each volume element (intersection of the
annular cylinder ``a`` with the spherical shell ``s``) are multiplied by a
volume normalization:

.. math::

   V_{\rm norm}[s, a] = V[s, a] / V_{\rm sphere}

where :math:`V_{\rm sphere}` is the volumne of the sphere enclosing the
outermost annulus.

With this convention the volume (:math:`V`) used above in calculating
the electron density for each shell is always :math:`V_{\rm sphere}`.

