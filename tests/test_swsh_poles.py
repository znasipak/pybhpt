"""Spin-weighted harmonics near the poles.

sY_lm vanishes like sin^|m+s|(theta/2) at theta = 0 and like cos^|m-s|(theta/2) at
theta = pi.  Building it as (sum_j A_j Y_jm)/sin^|s| (scalar harmonics divided by a
power of sin) makes the numerator cancel down to that small value, so the relative
error grows like eps/theta^(|s| + |m+s|): about 4e-3 for s = -2, m = 0 at theta = 5e-4.
Nearly polar orbits sample exactly these angles, and the noise leaked into every
high-k Teukolsky mode, so the k-sum never converged.

The reference below is the Wigner-d / Jacobi form, which contains no such cancellation:

    sY_lm(theta) = (-1)^s sqrt((2l+1)/(4 pi)) d^l_{m,-s}(theta),
    d^j_{m'm} = (-1)^lam sqrt(C(2j-k, k+a)/C(k+b, b)) sin^a(theta/2) cos^b(theta/2) P_k^(a,b)(cos theta).

It is evaluated with scipy, independently of pybhpt's C++.
"""
from math import lgamma, exp, sqrt, pi

import numpy as np
import pytest
from scipy.special import eval_jacobi

from pybhpt.geo import KerrGeodesic
from pybhpt.teuk import TeukolskyMode
from pybhpt.flux import FluxMode
from pybhpt.swsh import Yslm, Yslm_derivative


def _jacobi_params(s, l, m):
    # d^l_{m', mm} with m' = m, mm = -s
    mp, mm = m, -s
    k = min(l + mm, l - mm, l + mp, l - mp)
    if k == l + mm:
        a, lam = mp - mm, mp - mm
    elif k == l - mm:
        a, lam = mm - mp, 0
    elif k == l + mp:
        a, lam = mm - mp, 0
    else:
        a, lam = mp - mm, mp - mm
    b = 2*l - 2*k - a

    def lbin(n, r):
        return lgamma(n + 1) - lgamma(r + 1) - lgamma(n - r + 1)

    norm = (-1)**(lam + s)*sqrt((2*l + 1)/(4*pi))*exp(0.5*(lbin(2*l - k, k + a) - lbin(k + b, b)))
    return k, a, b, norm


def yslm_ref(s, l, m, th):
    k, a, b, norm = _jacobi_params(s, l, m)
    th = np.asarray(th, dtype=float)
    sh, ch = np.sin(th/2), np.cos(th/2)
    return norm*sh**a*ch**b*eval_jacobi(k, a, b, np.cos(th))


def yslm_derivative_ref(s, l, m, th):
    # d/dth [sh^a ch^b P(x)] = (a/2) sh^(a-1) ch^(b+1) P - (b/2) sh^(a+1) ch^(b-1) P
    #                          - 2 sh^(a+1) ch^(b+1) P'(x),  P' = (k+a+b+1)/2 P_(k-1)^(a+1,b+1)
    k, a, b, norm = _jacobi_params(s, l, m)
    th = np.asarray(th, dtype=float)
    sh, ch = np.sin(th/2), np.cos(th/2)
    x = np.cos(th)
    P = eval_jacobi(k, a, b, x)
    dP = 0.5*(k + a + b + 1)*eval_jacobi(k - 1, a + 1, b + 1, x) if k > 0 else 0.*x
    out = -2.*sh**(a + 1)*ch**(b + 1)*dP
    if a > 0:
        out = out + 0.5*a*sh**(a - 1)*ch**(b + 1)*P
    if b > 0:
        out = out - 0.5*b*sh**(a + 1)*ch**(b - 1)*P
    return norm*out


# Angles within 1e-2 of either pole, down to 1e-6.  All lie well inside the first node
# of every harmonic tested (about pi/l), so a pointwise relative error is meaningful.
POLE_ANGLES = np.array([1e-6, 1e-5, 1e-4, 4.7588920907e-4, 1e-3, 1e-2,
                        np.pi - 1e-2, np.pi - 1e-3, np.pi - 4.7588920907e-4, np.pi - 1e-4])
NEAR_POLE = 0.05
HARMONICS = [(s, l, m) for s in (-2, -1, 1, 2) for l in range(abs(s), 9) for m in range(-l, l + 1)]


def _scaled_err(approx, ref):
    """Max abs error over the max of |ref|: robust where ref crosses zero."""
    approx, ref = np.asarray(approx, dtype=float), np.asarray(ref, dtype=float)
    return np.max(np.abs(approx - ref))/np.max(np.abs(ref))


def _pointwise_relerr(approx, ref):
    approx, ref = np.asarray(approx, dtype=float), np.asarray(ref, dtype=float)
    return np.max(np.abs(approx - ref)/np.abs(ref))


def test_reference_matches_pybhpt_convention_in_the_bulk():
    """The scipy reference uses pybhpt's sign and normalisation, checked away from the poles."""
    th = np.linspace(0.3, 2.8, 7)
    worst = max(_scaled_err([Yslm(s, l, m, t) for t in th], yslm_ref(s, l, m, th))
                for s, l, m in HARMONICS)
    assert worst < 1e-13


@pytest.mark.parametrize("s, l, m", HARMONICS)
def test_yslm_accurate_near_poles(s, l, m):
    ref = yslm_ref(s, l, m, POLE_ANGLES)
    got = [Yslm(s, l, m, t) for t in POLE_ANGLES]
    keep = np.abs(ref) > 1e-280   # very high powers of sin(theta/2) underflow at 1e-6
    assert _pointwise_relerr(np.asarray(got)[keep], ref[keep]) < 1e-11


@pytest.mark.parametrize("s, l, m", HARMONICS)
def test_yslm_derivative_accurate_near_poles(s, l, m):
    ref = yslm_derivative_ref(s, l, m, POLE_ANGLES)
    got = [Yslm_derivative(s, l, m, t) for t in POLE_ANGLES]
    keep = np.abs(ref) > 1e-280
    assert _pointwise_relerr(np.asarray(got)[keep], ref[keep]) < 1e-10


# The p = 112, x = 4.76e-4 orbit from the GSN spherical dataset whose k-sum never
# converged.  a is tiny so the spheroidal correction (a*omega ~ 3e-8) does not mask
# the comparison.
NEAR_POLAR = dict(a=1e-6, p=112.10133694191467, e=0.0, x=4.7588920907111877e-4)


@pytest.mark.parametrize("l, m, k", [(2, 0, 2), (4, 0, 40), (3, 1, 39), (6, -2, 5), (5, 2, 1)])
def test_grid_harmonic_accurate_on_near_polar_orbit(l, m, k):
    """S and S' on the polar grid: pointwise near the poles, scaled everywhere."""
    orbit = KerrGeodesic(NEAR_POLAR["a"], NEAR_POLAR["p"], NEAR_POLAR["e"], NEAR_POLAR["x"], 2**10)
    mode = TeukolskyMode(-2, l, m, k, 0, orbit)
    mode.solve(orbit)
    th = np.asarray(mode.polarpoints)
    b = np.asarray(mode.couplingcoefficients)
    lmin = mode.mincouplingmode
    S_ref = sum(bi*yslm_ref(-2, lmin + i, m, th) for i, bi in enumerate(b) if bi != 0.)
    Sp_ref = sum(bi*yslm_derivative_ref(-2, lmin + i, m, th) for i, bi in enumerate(b) if bi != 0.)
    S, Sp = np.asarray(mode.polarsolutions), np.asarray(mode.polarderivatives)
    pole = (th < NEAR_POLE) | (th > np.pi - NEAR_POLE)
    assert pole.sum() >= 8   # the grid really does reach the poles
    assert _pointwise_relerr(S[pole], S_ref[pole]) < 1e-10
    assert _pointwise_relerr(Sp[pole], Sp_ref[pole]) < 1e-9
    assert _scaled_err(S, S_ref) < 1e-12
    assert _scaled_err(Sp, Sp_ref) < 1e-12


def test_near_polar_high_k_mode_is_not_noise():
    """In the Schwarzschild limit a polar circular orbit is a rotated equatorial one and only
    radiates at omega = m' Omega with |m'| <= l, so (l, m, k) = (2, 0, 40) must vanish up to
    O(a) corrections.  The pole error used to give it 6.5e-15, a tenth of the (2, 0, 2) flux
    times 1e-4, and it is what stopped the k-sum converging."""
    orbit = KerrGeodesic(NEAR_POLAR["a"], NEAR_POLAR["p"], NEAR_POLAR["e"], NEAR_POLAR["x"], 2**15)
    big = FluxMode(orbit, TeukolskyMode(-2, 2, 0, 2, 0, orbit)).infinityfluxes[0]
    tiny = FluxMode(orbit, TeukolskyMode(-2, 2, 0, 40, 0, orbit)).infinityfluxes[0]
    assert big > 1e-11
    assert abs(tiny) < 1e-16*big
