"""
Limited testing of deproject.deproject

:Copyright: Smithsonian Astrophysical Observatory (2019)
:Author: Douglas Burke (dburke@cfa.harvard.edu)
"""

# requies pytest 3.1 or greater for match argument to pytest.raises
import pytest

import numpy as np

from deproject.deproject import Deproject

from sherpa.astro import ui
from sherpa.astro.data import DataPHA

from astropy import units as u

# We need models, but do not want to depend on XSPEC availability (since
# we do not actually evaluate the models), so create "fake" versions
# of several XSPEC models. I was going to make this conditional on
# whether XSPEC is available, but it seems that over-writing things
# doesn't cause problems for us, so avoid that complexity.
#
# I use the actual parameter settings from Sherpa, but this isn't
# really needed here (could get away with only setting up a few
# parameters).
#
# from sherpa.astro.xspec import XSphabs, XSwabs, XSapec, XSmekal

from sherpa.models import Parameter, RegriddableModel1D
from sherpa.models.parameter import hugeval


class FakeXSPECModel(RegriddableModel1D):

    def calc(self, *args, **kwargs):
        "Always return 1 for each bin"

        if len(args) == 0:
            raise ValueError("No arguments")

        # When called with a single array (so XSPEC style) we normally
        # fill in the last value with 0 since it's invalid, but that
        # isn't needed here.
        #
        nbins = len(args[0])
        return np.ones(nbins)


class XSphabs(FakeXSPECModel):

    def __init__(self, name='phabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        FakeXSPECModel.__init__(self, name, (self.nH,))


class XSwabs(FakeXSPECModel):

    def __init__(self, name='wabs'):
        self.nH = Parameter(name, 'nH', 1., 0.0, 1.e5, 0.0, hugeval, '10^22 atoms / cm^2')
        FakeXSPECModel.__init__(self, name, (self.nH,))


class XSapec(FakeXSPECModel):

    def __init__(self, name='apec'):
        self.kT = Parameter(name, 'kT', 1., 0.008, 64.0, 0.0, hugeval, 'keV')
        self.Abundanc = Parameter(name, 'Abundanc', 1., 0., 5., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        FakeXSPECModel.__init__(self, name, (self.kT, self.Abundanc, self.redshift, self.norm))


class XSvapec(FakeXSPECModel):

    def __init__(self, name='vapec'):
        self.kT = Parameter(name, 'kT', 6.5, 0.0808, 68.447, 0.0, hugeval, 'keV')
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        FakeXSPECModel.__init__(self, name, (self.kT, self.He, self.C, self.N, self.O, self.Ne, self.Mg, self.Al, self.Si, self.S, self.Ar, self.Ca, self.Fe, self.Ni, self.redshift, self.norm))


class XSvvapec(FakeXSPECModel):

    def __init__(self, name='vvapec'):
        self.kT = Parameter(name, 'kT', 6.5, 0.0808, 68.447, 0.0, hugeval, 'keV')
        self.H = Parameter(name, 'H', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.He = Parameter(name, 'He', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Li = Parameter(name, 'Li', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Be = Parameter(name, 'Be', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.B = Parameter(name, 'B', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.C = Parameter(name, 'C', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.N = Parameter(name, 'N', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.O = Parameter(name, 'O', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.F = Parameter(name, 'F', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ne = Parameter(name, 'Ne', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Na = Parameter(name, 'Na', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Mg = Parameter(name, 'Mg', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Al = Parameter(name, 'Al', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Si = Parameter(name, 'Si', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.P = Parameter(name, 'P', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.S = Parameter(name, 'S', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Cl = Parameter(name, 'Cl', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ar = Parameter(name, 'Ar', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.K = Parameter(name, 'K', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ca = Parameter(name, 'Ca', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Sc = Parameter(name, 'Sc', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ti = Parameter(name, 'Ti', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.V = Parameter(name, 'V', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Cr = Parameter(name, 'Cr', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Mn = Parameter(name, 'Mn', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Fe = Parameter(name, 'Fe', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Co = Parameter(name, 'Co', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Ni = Parameter(name, 'Ni', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Cu = Parameter(name, 'Cu', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Zn = Parameter(name, 'Zn', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.Redshift = Parameter(name, 'Redshift', 0., -0.999, 10., -hugeval, hugeval, frozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)
        FakeXSPECModel.__init__(self, name, (self.kT, self.H, self.He, self.Li, self.Be, self.B, self.C, self.N, self.O, self.F, self.Ne, self.Na, self.Mg, self.Al, self.Si, self.P, self.S, self.Cl, self.Ar, self.K, self.Ca, self.Sc, self.Ti, self.V, self.Cr, self.Mn, self.Fe, self.Co, self.Ni, self.Cu, self.Zn, self.Redshift, self.norm))


class XSmekal(FakeXSPECModel):

    def __init__(self, name='mekal'):
        self.kT = Parameter(name, 'kT', 1., 0.0808, 79.9, 0.0, hugeval, 'keV')
        self.nH = Parameter(name, 'nH', 1., 1.e-5, 1.e19, 0.0, hugeval, 'cm-3', True)
        self.Abundanc = Parameter(name, 'Abundanc', 1., 0., 1000., 0.0, hugeval, frozen=True)
        self.redshift = Parameter(name, 'redshift', 0., -0.999, 10., -0.999, hugeval, frozen=True)
        self.switch = Parameter(name, 'switch', 1, 0, 2, 0, 2, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1.0, 0.0, 1.0e24, 0.0, hugeval)

        FakeXSPECModel.__init__(self, name, (self.kT, self.nH, self.Abundanc, self.redshift, self.switch, self.norm))


# Ensure the "fake" models are registered with the Sherpa session.
# I was expecting this to complain about over-writing existing classes -
# that is, if XSPEC support is available - but it doesn't seem to.
#
ui.add_model(XSphabs)
ui.add_model(XSwabs)
ui.add_model(XSapec)
ui.add_model(XSvapec)
ui.add_model(XSvvapec)
ui.add_model(XSmekal)


@pytest.fixture(autouse=True)
def cleanup_sherpa():
    """Ensure that the Sherpa state is cleaned before and after each run."""

    ui.clean()
    yield
    ui.clean()


def _test_data(annulus=0):
    """Return test data sets (very basic).

    At present there are no tests that make use of the data,
    so we could just use the same data set for each annulus.
    """

    n = 5 + annulus
    chans = np.arange(n)
    counts = np.repeat(1, n)

    return DataPHA('ann{}'.format(annulus), chans, counts)


@pytest.mark.parametrize("rad", [None, "radii", [], [10]])
def test_fail_radii_type(rad):
    """Need units for radii"""

    estr = "Argument 'radii' to function '__init__' " + \
           "has no 'unit' attribute. You may want to pass " + \
           "in an astropy Quantity instead."
    with pytest.raises(TypeError, match=estr):
        Deproject(rad)


# Pick a length since this is likely to be a common user error
# with the API. There is also the option that a future change
# would support lengths, so this acts as a check of the API (and
# is why a non-length is also included in the test).
#
@pytest.mark.parametrize("rad", [[0, 0.1, 0.2] * u.Mpc,
                                 [0.1, 0.2] * u.kg])
def test_fail_radii_convert(rad):
    """Not an angle."""

    estr = "Argument 'radii' to function '__init__' " + \
           "must be in units convertible to 'angle'."
    with pytest.raises(u.UnitsError, match=estr):
        Deproject(rad)


@pytest.mark.parametrize("rad", [[], [10]])
def test_fail_radii_num(rad):
    """Need at least 2 radii"""

    estr = 'radii parameter must be a sequence with at least two values'
    with pytest.raises(ValueError, match=estr):
        Deproject(rad * u.arcmin)


@pytest.mark.parametrize("rad", [[0, 0], [0, -10], [0, 10, 10, 20]])
def test_fail_radii_order(rad):
    """Fail if the radii are not in ascending order"""

    estr = 'radii parameter must be in increasing order'
    with pytest.raises(ValueError, match=estr):
        Deproject(rad * u.arcsec)


def test_fail_negative_radii():
    """Fail if radius is negative"""

    estr = 'radii must be >= 0'
    with pytest.raises(ValueError, match=estr):
        Deproject([-2, 0, 2] * u.rad)


def test_fail_extra_annulus():
    """We fail if more data than annuli"""

    dep = Deproject([0, 20] * u.arcsec)
    d0 = _test_data()
    dep.load_pha(d0, 0)

    d1 = _test_data(1)

    estr = 'Expected 0 <= annulus < 1 but sent 1'
    with pytest.raises(ValueError, match=estr):
        dep.load_pha(d1, 1)


def test_vnorm_1bin():
    """Check volume normalization when there's one bin"""

    dep = Deproject([0, 20] * u.arcmin)
    d0 = _test_data()
    dep.load_pha(d0, 0)

    assert not hasattr(dep, 'vol_norm')
    dep.set_source()
    assert hasattr(dep, 'vol_norm')

    vnorm = dep.vol_norm
    assert vnorm.shape == (1, 1)
    assert dep.vol_norm[0, 0] == pytest.approx(1.0)


@pytest.mark.parametrize("radii", [(0, 10, 20), (10, 30)])
def test_vnorm_unit_independence(radii):
    """The normalization terms are independent of the radius units"""

    dep1 = Deproject(radii * u.arcsec)
    dep2 = Deproject(radii * u.rad)

    dep1._calc_vol_norm()
    dep2._calc_vol_norm()

    for v1, v2 in zip(dep1.vol_norm.flatten(),
                      dep2.vol_norm.flatten()):
        diff = v1 - v2
        assert diff == pytest.approx(0.0)


def test_fail_source_missing_annulus1():
    """Does set_source fail if missing an annulus?"""

    dep = Deproject([0, 20] * u.arcsec)

    estr = 'missing data for annulus 0'
    with pytest.raises(ValueError, match=estr):
        dep.set_source()


def test_fail_source_missing_annulus2():
    """Does set_source fail if missing multiple annuli?"""

    dep = Deproject([0, 20, 40] * u.arcsec)

    # We do not want this treated as a regexp
    estr = r'missing data for annuli \[0, 1\]'
    with pytest.raises(ValueError, match=estr):
        dep.set_source()


@pytest.mark.parametrize("n,m", [(0, 1), (1, 0)])
def test_fail_source_missing_annulus3(n, m):
    """Does set_source fail if missing some annuliu"""

    dep = Deproject([0, 20, 40] * u.arcsec)
    d0 = _test_data()
    dep.load_pha(d0, n)

    estr = 'missing data for annulus {}'.format(m)
    with pytest.raises(ValueError, match=estr):
        dep.set_source()


@pytest.mark.parametrize("n,m", [(0, r"\[1, 2\]"),
                                 (1, r"\[0, 2\]"),
                                 (2, r"\[0, 1\]")])
def test_fail_source_missing_annulus4(n, m):
    """Does set_source fail if missing some annuli.

    Note the need to protect the expected string so that match
    does not take it as a regexp
    """

    dep = Deproject([0, 20, 40, 60] * u.arcsec)
    d0 = _test_data()
    dep.load_pha(d0, n)

    estr = 'missing data for annuli {}'.format(m)
    with pytest.raises(ValueError, match=estr):
        dep.set_source()


@pytest.mark.parametrize("srcmodel,comps,objs,mcomps",
                         [(None,
                           [{'type': 'xsphabs', 'start': 0, 'end': 7},
                            {'type': 'xsapec', 'start': 8, 'end': 14}],
                           [XSphabs('xsphabs.xsphabs_0'),
                            XSapec('xsapec.xsapec_0')],
                           [{'type': 'xsphabs',
                             'name': 'xsphabs_0',
                             'shell': 0},
                            {'type': 'xsapec',
                             'name': 'xsapec_0',
                             'shell': 0}]),

                          ('xswabs * xsmekal',
                          [{'type': 'xswabs', 'start': 0, 'end': 6},
                           {'type': 'xsmekal', 'start': 9, 'end': 16}],
                           [XSwabs('xswabs.xswabs_0'),
                            XSmekal('xsmekal.xsmekal_0')],
                           [{'type': 'xswabs',
                             'name': 'xswabs_0',
                             'shell': 0},
                            {'type': 'xsmekal',
                             'name': 'xsmekal_0',
                             'shell': 0}])])
def test_set_source(srcmodel, comps, objs, mcomps):
    """Check the source model can be set, with 1 annulus.

    Basic checks of the structures set in place by a set_source
    call.
    """

    dep = Deproject([0, 20] * u.arcsec)
    d0 = _test_data()
    dep.load_pha(d0, 0)

    if srcmodel is None:
        dep.set_source()
        assert dep.srcmodel == 'xsphabs*xsapec'

    else:
        dep.set_source(srcmodel)
        assert dep.srcmodel == srcmodel

    # can do a direct comparison of srcmodel_comps
    assert dep.srcmodel_comps == comps

    # for model_comps, want to remove the object key from the comparison
    #
    to_compare = []
    to_objs = []
    for mc in dep.model_comps:
        to_objs.append(mc['object'])
        mc = mc.copy()
        del mc['object']
        to_compare.append(mc)

    assert to_compare == mcomps
    for a, b in zip(objs, to_objs):
        assert a.type == b.type
        assert a.name == b.name


@pytest.mark.parametrize("nmodel", [None,
                                    "xsphabs * xsapec",
                                    "xswabs * xsmekal"])
def test_reset_source(nmodel):
    """Check the source model can be reset.

    Note there's no check on the downstream implications of this
    change, as this is just to check the initial bug.
    """

    dep = Deproject([0, 20] * u.arcsec)
    d0 = _test_data()
    dep.load_pha(d0, 0)

    dep.set_source()
    if nmodel is None:
        dep.set_source()
    else:
        dep.set_source(nmodel)


def test_reset_source_explicit():
    """Check the before and after settings of set_source.

    This overlaps previous checks, but uses simpler models
    (single component).
    """

    dep = Deproject([0, 20] * u.arcsec)
    d0 = _test_data()
    dep.load_pha(d0, 0)

    dep.set_source('xsvapec')
    assert dep.srcmodel_comps == [{'type': 'xsvapec', 'start': 0, 'end': 7}]
    assert len(dep.model_comps) == 1
    mcomps = dep.model_comps[0]
    assert len(mcomps) == 4
    assert mcomps['type'] == 'xsvapec'
    assert mcomps['name'] == 'xsvapec_0'
    assert mcomps['shell'] == 0
    assert mcomps['object'].type == 'xsvapec'
    assert mcomps['object'].name == 'xsvapec.xsvapec_0'

    dep.set_source('xsvvapec')
    assert dep.srcmodel_comps == [{'type': 'xsvvapec', 'start': 0, 'end': 8}]
    assert len(dep.model_comps) == 1
    mcomps = dep.model_comps[0]
    assert len(mcomps) == 4
    assert mcomps['type'] == 'xsvvapec'
    assert mcomps['name'] == 'xsvvapec_0'
    assert mcomps['shell'] == 0
    assert mcomps['object'].type == 'xsvvapec'
    assert mcomps['object'].name == 'xsvvapec.xsvvapec_0'


def test_set_redshift_no_model_get():
    """redshift attribute."""

    dep = Deproject([0, 20] * u.arcsec)

    estr = 'Parameter redshift not found in any model component'
    with pytest.raises(ValueError, match=estr):
        dep.redshift


def test_set_redshift_no_model_set():
    """redshift attribute."""

    dep = Deproject([0, 20] * u.arcsec)
    dep.redshift = 0.12
    assert dep.redshift == pytest.approx(0.12)


def test_set_redshift_model_get():
    """redshift attribute."""

    dep = Deproject([0, 20, 40] * u.arcsec)
    for ann in [0, 1]:
        d = _test_data(ann)
        dep.load_pha(d, ann)

    # Note: does not fail if the redshifts are different
    dep.set_source()
    ui.get_model_component('xsapec_0').redshift = 0.12
    assert dep.redshift == pytest.approx(0.12)


def test_set_redshift_model_change_get():
    """redshift attribute"""

    dep = Deproject([0, 20] * u.arcsec)
    d0 = _test_data()
    dep.load_pha(d0, 0)

    dep.set_source()
    ui.get_model_component('xsapec_0').redshift = 0.12
    ui.get_model_component('xsapec_0').redshift = 0.21
    assert dep.redshift == pytest.approx(0.21)


def test_set_redshift_model_get_change_get():
    """redshift attribute

    At present the redshift, once set, is not changed.
    """

    dep = Deproject([0, 20] * u.arcsec)
    d0 = _test_data()
    dep.load_pha(d0, 0)

    dep.set_source()
    ui.get_model_component('xsapec_0').redshift = 0.12

    # ask for the redshift, but ignore the result
    dep.redshift

    # change the redshift
    ui.get_model_component('xsapec_0').redshift = 0.21

    # redshift has not changed
    assert dep.redshift == pytest.approx(0.12)


def test_init_angdist_type_fail():
    """Does setting angdist to an invalid value (no type).
    """

    estr = "Argument 'angdist' to function '__init__' " + \
           "has no 'unit' attribute. You may want to pass " + \
           "in an astropy Quantity instead."
    with pytest.raises(TypeError, match=estr):
        Deproject([0, 10, 20] * u.arcsec, angdist=123.4)


def test_set_angdist_type_fail():
    """Does setting the angdist attribute with no type fail?"""

    dep = Deproject([0, 10, 20] * u.arcsec)

    estr = "Argument 'angdist' to function '_set_angdist' " + \
           "has no 'unit' attribute. You may want to pass " + \
           "in an astropy Quantity instead."
    with pytest.raises(TypeError, match=estr):
        dep.angdist = 123.4


def test_init_angdist_convert_fail():
    """Does setting angdist to an invalid type fail.
    """

    estr = "Argument 'angdist' to function '__init__' " + \
           "must be in units convertible to 'length'."
    with pytest.raises(u.UnitsError, match=estr):
        Deproject([0, 10, 20] * u.arcsec, angdist=123.4 * u.J)


def test_set_angdist_convert_fail():
    """Does setting angdist to an invalid type fail.
    """

    dep = Deproject([0, 10, 20] * u.arcsec)

    estr = "Argument 'angdist' to function '_set_angdist' " + \
           "must be in units convertible to 'length'."
    with pytest.raises(u.UnitsError, match=estr):
        dep.angdist = 123.4 * u.J


@pytest.mark.parametrize("angdist", [0.0 * u.Mpc, -1.2e24 * u.cm])
def test_init_angdist_invalid_fail(angdist):
    """Does setting angdist to an invalid quantity fail.
    """

    with pytest.raises(ValueError):
        Deproject([0, 10, 20] * u.arcsec, angdist=angdist)


@pytest.mark.parametrize("angdist", [0.0 * u.Mpc, -1.2e24 * u.cm])
def test_set_angdist_invalid_fail(angdist):
    """Does setting angdist to an invalid quantity fail.
    """

    dep = Deproject([0, 10, 20] * u.arcsec)

    with pytest.raises(ValueError):
        dep.angdist = angdist


@pytest.mark.parametrize("theta", [None, "theta", 0, 360, [10]])
def test_fail_theta_type(theta):
    """Need units for theta"""

    estr = "Argument 'theta' to function '__init__' " + \
           "has no 'unit' attribute. You may want to pass " + \
           "in an astropy Quantity instead."
    with pytest.raises(TypeError, match=estr):
        Deproject([0, 10, 20] * u.arcsec, theta=theta)


@pytest.mark.parametrize("theta", [10 * u.Mpc, [0.1, 0.2] * u.kg])
def test_fail_theta_convert(theta):
    """Not an angle."""

    estr = "Argument 'theta' to function '__init__' " + \
           "must be in units convertible to 'angle'."
    with pytest.raises(u.UnitsError, match=estr):
        Deproject([0, 10, 20] * u.arcmin, theta=theta)


@pytest.mark.parametrize("theta", [0 * u.deg, 0 * u.rad,
                                   -10 * u.deg, -1 * u.rad])
def test_theta_out_of_band_small_fail(theta):
    """These theta values are invalid."""

    emsg = 'theta must be > 0 degrees'
    with pytest.raises(ValueError, match=emsg):
        Deproject([0, 10, 20] * u.arcsec, theta=theta)


@pytest.mark.parametrize("theta", [361 * u.deg, 2.1 * np.pi * u.rad])
def test_theta_out_of_band_large_fail(theta):
    """These theta values are invalid."""

    emsg = 'theta must be <= 360 degrees'
    with pytest.raises(ValueError, match=emsg):
        Deproject([0, 10, 20] * u.arcsec, theta=theta)


@pytest.mark.parametrize("theta", [(10, 20), (30, 40, 350, 360)])
def test_theta_nwrong_fail(theta):
    """The number of theta are wrong."""

    emsg = 'theta must be a scalar or match the number of annuli'
    with pytest.raises(ValueError, match=emsg):
        Deproject([0, 1, 2, 3] * u.arcmin, theta=theta * u.deg)


@pytest.mark.parametrize("full", [360 * u.deg, 2 * np.pi * u.rad])
@pytest.mark.parametrize("scale", [1, 0.5, 0.1])
@pytest.mark.parametrize("radii", [(0, 10), (0, 10, 20), (0, 10, 40),
                                   (10, 20), (10, 20, 40)])
def test_theta_subset_norm(full, scale, radii):
    """Does changing theta scale the normalization?
    """

    radii = radii * u.arcmin
    nann = len(radii) - 1
    dep_full = Deproject(radii, theta=full)
    dep_scale = Deproject(radii, theta=full * scale)

    for i in range(nann):
        d0 = _test_data()
        dep_full.load_pha(d0, i)
        dep_scale.load_pha(d0, i)

    dep_full.set_source()
    dep_scale.set_source()

    v_full = dep_full.vol_norm
    v_scale = dep_scale.vol_norm
    assert v_full.shape == (nann, nann)
    assert v_scale.shape == v_full.shape

    # the following would be slightly shorter with NumPy-array-aware
    # comparison routines
    #
    # the matrix is in lower-triangular form, so the upper triangular
    # area (above the diagonal) is 0.
    #
    if nann > 1:
        zeros = np.triu_indices(nann, k=1)
        for z in v_full[zeros]:
            assert z == pytest.approx(0)

        for z in v_scale[zeros]:
            assert z == pytest.approx(0)

    # the lower-triangular area (including the diagonal) is positive
    # and scales with theta
    #
    idx = np.tril_indices(nann, k=0)
    z_full = v_full[idx]
    z_scale = v_scale[idx]
    assert np.all(z_full > 0.0)
    assert np.all(z_scale > 0.0)

    for r in z_scale / z_full:
        assert r == pytest.approx(scale)


@pytest.mark.parametrize("full", [360 * u.deg, 2 * np.pi * u.rad])
@pytest.mark.parametrize("scales,radii", [((1, 1), (0, 10, 20)),
                                          ((1, 1), (10, 20, 40)),
                                          ((0.4, 0.8), (0, 10, 20)),
                                          ((0.8, 0.4), (10, 30, 40)),
                                          ((0.5, 1.0, 0.7), (0, 2, 5, 8)),
                                          ((0.5, 1.0, 0.7), (1, 2, 5, 8))])
def test_theta_per_annuli_norm(full, scales, radii):
    """Does changing theta (per annulus) scale the normalization?

    This complements test_theta_per_annuli_norm by having the
    theta value vary per annulus.
    """

    radii = radii * u.arcmin
    nann = len(radii) - 1
    scales = np.asarray(scales)
    thetas = scales * full

    dep_full = Deproject(radii, theta=full)
    dep_scale = Deproject(radii, theta=thetas)

    for i in range(nann):
        d0 = _test_data()
        dep_full.load_pha(d0, i)
        dep_scale.load_pha(d0, i)

    dep_full.set_source()
    dep_scale.set_source()

    v_full = dep_full.vol_norm
    v_scale = dep_scale.vol_norm
    assert v_full.shape == (nann, nann)
    assert v_scale.shape == v_full.shape

    # the following would be slightly shorter with NumPy-array-aware
    # comparison routines
    #
    # the matrix is in lower-triangular form, so the upper triangular
    # area (above the diagonal) is 0.
    #
    if nann > 1:
        zeros = np.triu_indices(nann, k=1)
        for z in v_full[zeros]:
            assert z == pytest.approx(0)

        for z in v_scale[zeros]:
            assert z == pytest.approx(0)

    # the lower-triangular area (including the diagonal) is positive
    # and scales with theta
    #
    idx = np.tril_indices(nann, k=0)
    z_full = v_full[idx]
    z_scale = v_scale[idx]
    assert np.all(z_full > 0.0)
    assert np.all(z_scale > 0.0)

    # Stack the scales array up to create a n by n matrix, so that
    # we can apply the lower-triangle index and get the expected
    # scale factor
    scales_mat = np.resize(scales, (nann, nann))[idx]

    for r, scale in zip(z_scale / z_full, scales_mat):
        assert r == pytest.approx(scale)


@pytest.mark.parametrize("units", ["arcsec", "deg", "rad",
                                   "Mpc", "kpc", "cm",
                                   u.Unit("arcmin"), u.Unit("pc")])
def test_get_radii_accepts(units):
    """Can we ask for a variety of units.

    There is limited checking of the output.
    """

    dep = Deproject([0, 10, 15] * u.arcsec,
                    angdist=127.4 * u.Mpc)

    rlo, rhi = dep.get_radii(units)

    assert len(rlo) == len(rhi)
    assert len(rlo) == 2

    # Ensure the output is the correct type (length or angle)
    #
    rlo.to(u.Unit(units))

    for r1, r2 in zip(rlo, rhi):
        assert r2 > r1

    # Relies on radii[2] = 1.5 * radii[1] in the input to Deproject.
    # Apparently need to extract the value from the quantity otherwise
    # the test fails due to unbounded recursion.
    #
    delta = rhi[1] - 1.5 * rhi[0]
    assert delta.value == pytest.approx(0)


def test_get_radii_one_annulus():
    """Nothing breaks down when there's only one annulus"""

    dep = Deproject([4.8e-5, 7.2e-5] * u.rad)
    dep.angdist = 20.0 * u.Mpc
    rlo, rhi = dep.get_radii('kpc')

    assert len(rlo) == 1
    assert len(rhi) == 1

    assert rlo.unit.physical_type == 'length'
    assert rhi.unit.physical_type == 'length'

    assert rlo.unit == u.Unit('kpc')
    assert rhi.unit == u.Unit('kpc')

    # Using an absolute tolerance just as a regression test. Since
    # the angular diameter distance explicitly has been given then
    # there is no issue about changing the default Cosmology.
    #
    assert rlo[0].value == pytest.approx(0.960, abs=0.001)
    assert rhi[0].value == pytest.approx(1.440, abs=0.001)


@pytest.mark.parametrize("radii", [(0, 10), (0, 10, 20), (10, 30, 50, 90)])
def test_get_shells_no_data(radii):
    """No data."""

    dep = Deproject(radii * u.arcsec)
    shells = dep.get_shells()
    assert len(shells) == dep.nshell

    # shells is in outer-to-inner order
    shells.reverse()
    for i, shell in enumerate(shells):
        assert shell['annuli'] == [i]
        assert shell['dataids'] == []


@pytest.mark.parametrize("radii", [(0, 10), (0, 10, 20), (10, 30, 50, 90)])
def test_get_shells_no_ties(radii):
    """No ties between the shells"""

    dep = Deproject(radii * u.arcsec)
    for i in range(dep.nshell):
        d = _test_data(i)
        dep.load_pha(d, i)

    dep.set_source()
    shells = dep.get_shells()
    assert len(shells) == dep.nshell

    # shells is in outer-to-inner order
    shells.reverse()
    for i, shell in enumerate(shells):
        assert shell['annuli'] == [i]
        assert shell['dataids'] == [i]


@pytest.mark.parametrize("inner,outer", [(0, 3), (3, 0), (1, 3), (3, 1)])
def test_get_shells_not_consecutive_fail(inner, outer):
    """Can not tie parameters across groups.

    This is only checked when the annuli grouping is required,
    which is simulated with a call to get_radii.
    """

    dep = Deproject([1, 2, 3, 4, 5] * u.arcsec)
    for i in range(dep.nshell):
        d = _test_data(i)
        dep.load_pha(d, i)

    dep.set_source('xsvapec')
    dep.tie_par('xsvapec.ne', inner, outer)

    anns = [str(s) for s in sorted([inner, outer])]
    estr = r'Non-consecutive annuli are tied together: \[' + \
           r'{}\]'.format(" ".join(anns))
    with pytest.raises(ValueError, match=estr):
        dep.get_shells()


@pytest.mark.parametrize("outer,annuli,dataids",
                         [(3, [[2, 3], [1], [0]], [[0, 1], [2], [3]]),
                          (2, [[3], [1, 2], [0]], [[0], [1, 2], [3]]),
                          (1, [[3], [2], [0, 1]], [[0], [1], [2, 3]])])
def test_tie_par1(outer, annuli, dataids):
    """One pair of shells is tied."""

    dep = Deproject([1, 2, 3, 4, 5] * u.arcsec)

    # Reverse the order so that annulus != dataid
    for i in [3, 2, 1, 0]:
        d = _test_data(i)
        dep.load_pha(d, i)

    dep.set_source()
    inner = outer - 1
    dep.tie_par('xsapec.kt', inner, outer)

    shells = dep.get_shells()
    assert len(shells) == 3

    for shell, ann, dids in zip(shells, annuli, dataids):
        assert shell['annuli'] == ann
        assert shell['dataids'] == dids


@pytest.mark.parametrize("outer,annuli,dataids",
                         [(3, [[1, 2, 3], [0]], [[0, 1, 2], [3]]),
                          (2, [[3], [0, 1, 2]], [[0], [1, 2, 3]])])
def test_tie_par2(outer, annuli, dataids):
    """Two pair of shells is tied."""

    dep = Deproject([1, 2, 3, 4, 5] * u.arcsec)

    # Reverse the order so that annulus != dataid
    for i in [3, 2, 1, 0]:
        d = _test_data(i)
        dep.load_pha(d, i)

    dep.set_source()
    inner = outer - 2
    dep.tie_par('xsapec.kt', inner, outer-1, outer)

    shells = dep.get_shells()
    assert len(shells) == 2

    for shell, ann, dids in zip(shells, annuli, dataids):
        assert shell['annuli'] == ann
        assert shell['dataids'] == dids


@pytest.mark.parametrize("outer", [3, 2, 1])
def test_untie_par1(outer):
    """Can untie a parameter."""

    nshell = 4
    dep = Deproject([1, 2, 3, 4, 5] * u.arcsec)

    # Reverse the order so that annulus != dataid
    for i in [3, 2, 1, 0]:
        d = _test_data(i)
        dep.load_pha(d, i)

    dep.set_source()
    inner = outer - 1

    # Rely on previous tests to check that tie_par works
    dep.tie_par('xsapec.kt', inner, outer)

    # untie the parameter
    dep.untie_par('xsapec.kt', outer)

    shells = dep.get_shells()
    assert len(shells) == nshell

    shells.reverse()
    for i, shell in enumerate(shells):
        assert shell['annuli'] == [i]
        assert shell['dataids'] == [nshell - 1 - i]


@pytest.mark.parametrize("outer", [3, 2])
def test_untie_par2(outer):
    """Can untie two parameters."""

    nshell = 4
    dep = Deproject([1, 2, 3, 4, 5] * u.arcsec)

    # Reverse the order so that annulus != dataid
    for i in [3, 2, 1, 0]:
        d = _test_data(i)
        dep.load_pha(d, i)

    dep.set_source()
    inner = outer - 2

    # Rely on previous tests to check that tie_par works
    dep.tie_par('xsapec.kt', inner, outer - 1, outer)

    # untie the parameter
    dep.untie_par('xsapec.kt', outer - 1, outer)

    shells = dep.get_shells()
    assert len(shells) == nshell

    shells.reverse()
    for i, shell in enumerate(shells):
        assert shell['annuli'] == [i]
        assert shell['dataids'] == [nshell - 1 - i]
