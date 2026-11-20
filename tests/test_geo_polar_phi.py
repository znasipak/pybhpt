"""The polar part of the azimuthal accumulation, Delta phi_theta, on nearly polar orbits.

d phi_theta / d lambda = Lz/(1 - z^2) - a E with z = z2 cos(chi) and 1 - z2^2 = x^2, so on
a nearly polar orbit the integrand is a spike of height ~1/x^2 and width ~x at each pole
crossing, and phi_theta jumps by ~pi there.  Its Fourier series in the Mino phase decays
only like (1 - x)^n.  It used to be truncated at 2^15 coefficients with the tail still at
~7e-7 for x = 4.76e-4, which put a phase error into every m != 0 Teukolsky mode and an
8.9e-4 relative error into dL_z/dt.
"""
import numpy as np
import pytest
from scipy.integrate import quad

from pybhpt.geo import KerrGeodesic
from pybhpt.teuk import TeukolskyMode
from pybhpt.flux import FluxMode


def phi_theta_reference(geo, chi):
    """Delta phi_theta at polar anomaly chi, by adaptive quadrature of the raw integrand."""
    a = geo.a
    En, Lz, _ = geo.orbitalconstants
    z1, z2 = geo.polarroots

    x = abs(geo.x)

    def dphi_dchi(c):
        z = z2*np.cos(c)
        dmino = 1./(a*np.sqrt(1. - En*En)*np.sqrt(z1*z1 - z*z))
        # 1 - z^2 = x^2 cos^2 + sin^2 exactly (1 - z2^2 = x^2): no cancellation at the pole
        return (Lz/((x*np.cos(c))**2 + np.sin(c)**2) - a*En)*dmino
    # break the interval at the spikes so quad resolves them
    brk = [b for b in (x, 10*x, np.pi - 10*x, np.pi - x) if 0 < b < np.pi]

    def cumulative(c):
        pts = [b for b in brk if b < c]
        return quad(dphi_dchi, 0., c, points=pts or None, limit=2000,
                    epsabs=1e-15, epsrel=1e-14)[0]

    total = cumulative(np.pi)
    # chi = pi is Mino phase pi; the secular part is linear in the Mino phase
    return np.array([cumulative(c) for c in chi]), total


ORBITS = [
    (0.0029461646102816605, 112.10133694191467, 0.00047588920907111877),   # GSN near-polar row
    (0.0029461646102816605, 112.10133694191467, 0.0042776),
    (0.9, 6.0, 0.02),
    (0.9, 8.0, -0.3),
    (0.5, 10.0, 0.7),
]


@pytest.mark.parametrize("a, p, x", ORBITS + [(1e-6, 112.10133694191467, 0.00047588920907111877)])
def test_azimuthal_mino_frequency_matches_quadrature(a, p, x):
    """Upsilon_phi on circular spherical orbits.  The polar part is the mean of the spiked
    integrand; a uniform trapezoid rule capped at 2^15 points only reaches e^(-N x), so it
    was 3.4e-7 off at x = 4.76e-4, which shifted every mode frequency by m Omega 3.4e-7
    and split the dEdt of +m and -m modes by ~1e-6."""
    geo = KerrGeodesic(a, p, 0., x, 2**4)
    En, Lz, _ = geo.orbitalconstants
    z1, z2 = geo.polarroots

    def dmino(c):
        return 1./(a*np.sqrt(1. - En*En)*np.sqrt(z1*z1 - (z2*np.cos(c))**2))

    ax = abs(x)
    pts = [b for b in (ax, 10*ax, np.pi - 10*ax, np.pi - ax) if 0 < b < np.pi]
    lam = quad(dmino, 0., np.pi, points=pts, limit=2000, epsabs=0., epsrel=2e-14)[0]
    phith = quad(lambda c: (Lz/((ax*np.cos(c))**2 + np.sin(c)**2) - a*En)*dmino(c), 0., np.pi,
                 points=pts, limit=2000, epsabs=0., epsrel=2e-14)[0]
    Delta = p*p - 2*p + a*a
    upsilon_phi = phith/lam + a*(En*(p*p + a*a) - a*Lz)/Delta
    assert abs(geo.minofrequencies[3]/upsilon_phi - 1.) < 1e-12


@pytest.mark.parametrize("a, p, x", ORBITS)
def test_polar_phi_samples_match_quadrature(a, p, x):
    N = 2**8
    geo = KerrGeodesic(a, p, 0., x, N, parametrization="darwin")
    q = np.asarray(geo.polarphase)
    chi = np.linspace(0., np.pi, q.size)    # Darwin samples are uniform in chi
    # the grid really is uniform in chi: theta = arccos(z2 cos chi)
    assert np.allclose(np.cos(np.asarray(geo.polarpoints)), geo.polarroots[1]*np.cos(chi), atol=1e-13)
    cum, total = phi_theta_reference(geo, chi)
    ref = cum - total*q/np.pi
    got = np.asarray(geo.azimuthalpolar)
    assert np.max(np.abs(got - ref)) < 1e-11


def test_near_polar_l2_modes_distribute_lz_like_a_rotated_circular_orbit():
    """In the Schwarzschild limit a spherical orbit is a circular orbit tilted by
    iota = arccos(x).  Its l = 2, omega = 2 Omega radiation splits over m as
    |d^2_{m,2}(iota)|^2, and sum_m m |d^2_{m,2}|^2 = 2 cos(iota) = 2x, so
    sum_m m Edot_{2,m,2-m} / sum_m Edot_{2,m,2-m} = 2x.  The m != 0 terms are the ones
    that see phi_theta.  Spin adds a term ~2e-3 a/x to this ratio (about 5e-6 at a = 1e-6, as
    GSN's data also show), so a is taken small enough for that to be negligible."""
    a, p, x = 1e-9, 112.10133694191467, 0.00047588920907111877
    geo = KerrGeodesic(a, p, 0., x, 2**12)
    E = {}
    for m in range(-2, 3):
        mode = TeukolskyMode(-2, 2, m, 2 - m, 0, geo)
        E[m] = FluxMode(geo, mode).infinityfluxes[0]
    ratio = sum(m*E[m] for m in E)/sum(E.values())
    assert abs(ratio/(2*x) - 1.) < 1e-6
