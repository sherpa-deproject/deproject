
.. include:: references.rst

Installation
============

As of :ref:`version 0.2.0 <changes_020>`, the :mod:`deproject` package
can be installed directly from `PyPI`_. This requires `CIAO 4.11
<http://cxc.harvard.edu/ciao/>`_ or later (in general the latest
release of CIAO should be used). The module *can* be used with the
`standalone release of Sherpa
<https://sherpa.readthedocs.io/en/latest/ciao.html>`_, but it is only
useful if Sherpa has been built with `XSPEC support
<https://sherpa.readthedocs.io/en/latest/install.html#xspec>`_ which
is trickier to achieve than we would like.

Requirements
------------

The package uses `Astropy`_ and `SciPy`_, for units support and
cosmological-distance calculations. It is assumed that
`Matplotlib`_ is available for plotting.

It is **strongly** advised to use the latest CIAO (at the time
of writing this is CIAO 4.14) installed with conda. Support
for CIAO installed via the ciao-install script is limited.

Using pip
---------

It should be as simple as starting the CIAO environment - this
depends on whether CIAO was installed via `ciao-install` or
`conda` - and then saying::

  pip install deproject

This approach should also work if you are using the standalone
version of Sherpa.

Manual installation
-------------------

The source is available on github at
`<https://github.com/sherpa-deproject/deproject>`_, with releases available
at `<https://github.com/sherpa-deproject/deproject/releases>`_.

After downloading the source code (whether from a release or
by cloning the repository) and moving into the directory
(`deproject-<version>` or `deproject`), installation just
requires::

  pip install .

.. note::
   This command should only be run after setting up CIAO or whatever
   Python environment contains your Standalone Sherpa installation.

Test
----

The source installation includes a basic test suite, which can be
run with

::
   
  pytest

Example data
------------

The example data can be download from either
http://cxc.cfa.harvard.edu/contrib/deproject/downloads/m87.tar.gz
or from `GitHub <https://github.com/sherpa-deproject/deproject-data>`_.
The source distribution includes scripts - in the
`examples directory <https://github.com/sherpa-deproject/deproject/tree/master/examples>`_
- that can
be used to replicate both the :doc:`basic M87 example <examples/m87>`
and the follow-on example :doc:`combining annuli <examples/tie>`.

As an example (from within an IPython session, such as the
Sherpa shell in CIAO)::

  >>> %run fit_m87.py
  ...
  ... a lot of screen output will whizz by
  ...

The current density estimates can then be displayed with::

  >>> dep.density_plot()

.. image:: examples/m87_density.png
  
the reduced-statistic for the fit to each shell with::

  >>> dep.fit_plot('rstat')

.. image:: examples/m87_rstat.png
  
and the fit results for the first annulus can be displayed using
the Sherpa functions::

  >>> set_xlog()
  >>> plot_fit_delchi(0)

.. image:: examples/m87_ann0_fit.png
