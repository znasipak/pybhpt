"""Fluxes of weak-field orbits.

Every mode with |omega| <= 1e-8 used to be treated as zero-frequency and given
zero flux, so a circular orbit beyond r ~ 3.4e5 M (where the dominant
omega = 2 r^{-3/2} drops below 1e-8) radiated nothing at all. Far out, the
l = |m| = 2 modes carry the quadrupole flux dE/dt = (32/5) r^{-5}, up to a
1PN correction of relative size ~4/r.
"""
import pytest

from pybhpt.geo import KerrGeodesic
from pybhpt.teuk import TeukolskyMode
from pybhpt.flux import FluxMode


def _l2m2_energy_flux_to_infinity(a, r):
    geo = KerrGeodesic(a, r, 0., 1., 2**4)
    flux = FluxMode(geo, TeukolskyMode(-2, 2, 2, 0, 0, geo))
    return 2.*flux.infinityfluxes[0]  # m = +2 and m = -2


@pytest.mark.parametrize("r", [1e5, 4e5, 1e6, 2.6e7])
def test_weak_field_circular_flux_is_quadrupolar(r):
    edot = _l2m2_energy_flux_to_infinity(0.0029, r)
    assert edot == pytest.approx(32./5.*r**-5, rel=10./r, abs=0.)


@pytest.mark.parametrize("r", [4e5, 2.6e7])
def test_modes_below_the_old_cutoff_carry_flux(r):
    geo = KerrGeodesic(0.0029, r, 0., 1., 2**4)
    teuk = TeukolskyMode(-2, 2, 2, 0, 0, geo)
    flux = FluxMode(geo, teuk)
    assert 0 < abs(teuk.frequency) < 1e-8
    assert flux.infinityfluxes[0] > 0
    assert flux.infinityfluxes[1] > 0  # the Lz flux prefactor had its own 1e-8 cutoff


def test_a_zero_frequency_mode_still_carries_no_flux():
    """(l, m, k, n) = (2, 0, 0, 0) has omega = 0 exactly."""
    geo = KerrGeodesic(0.5, 10., 0., 1., 2**4)
    flux = FluxMode(geo, TeukolskyMode(-2, 2, 0, 0, 0, geo))
    assert flux.infinityfluxes[0] == 0.
    assert flux.horizonfluxes[0] == 0.
