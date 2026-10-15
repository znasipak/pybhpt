"""Spin-weighted spheroidal harmonics at complex spheroidicity (GitHub issue #28).

For complex gamma the angular equation is complex-symmetric rather than Hermitian, so the
harmonics are normalized with the bilinear pairing  int S^2 sin(th) dth = 1/(2 pi)  (not
int |S|^2), which keeps S analytic in gamma and makes distinct-j modes orthogonal.
"""
import numpy as np
import pytest
from scipy.integrate import quad

from pybhpt.swsh import SpinWeightedSpheroidalHarmonic, _bilinear_normalize

CASES = [
    (-2, 2, 2, 0.6 - 0.4j),   # example from issue #28
    (-2, 3, 1, 1.5 - 1.0j),
    (2, 2, -2, 0.36 - 0.225j),
    (0, 3, 2, 0.8 + 0.5j),
]


def bilinear_integral(f, g):
    re = quad(lambda th: (f(th)*g(th)*np.sin(th)).real, 0, np.pi, epsabs=1e-13, epsrel=1e-12)[0]
    im = quad(lambda th: (f(th)*g(th)*np.sin(th)).imag, 0, np.pi, epsabs=1e-13, epsrel=1e-12)[0]
    return re + 1j*im


@pytest.mark.parametrize("s, j, m, g", CASES)
def test_bilinear_normalization(s, j, m, g):
    S = SpinWeightedSpheroidalHarmonic(s, j, m, g)
    assert abs(np.sum(S.couplingcoefficients**2) - 1.) < 1e-12
    assert abs(bilinear_integral(S.Sslm, S.Sslm) - 1./(2.*np.pi)) < 1e-10


def test_bilinear_orthogonality_distinct_j():
    s, m, g = -2, 2, 0.6 - 0.4j
    S2 = SpinWeightedSpheroidalHarmonic(s, 2, m, g)
    S3 = SpinWeightedSpheroidalHarmonic(s, 3, m, g)
    assert abs(bilinear_integral(S2.Sslm, S3.Sslm)) < 1e-10


@pytest.mark.parametrize("s, j, m, g", CASES)
def test_coefficients_analytic_in_gamma(s, j, m, g):
    # Cauchy-Riemann: the derivative along real and imaginary directions must agree.
    # Hermitian normalization breaks this (|.| is not analytic).
    h = 1e-5
    def b(gg):
        return SpinWeightedSpheroidalHarmonic(s, j, m, gg).couplingcoefficients
    n = min(len(b(g + h)), len(b(g - h)), len(b(g + 1j*h)), len(b(g - 1j*h)))
    d_re = (b(g + h)[:n] - b(g - h)[:n])/(2.*h)
    d_im = (b(g + 1j*h)[:n] - b(g - 1j*h)[:n])/(2j*h)
    assert np.max(np.abs(d_re - d_im)) < 1e-6*np.max(np.abs(d_re))


def test_continuous_with_real_gamma():
    s, j, m, g = -2, 2, 2, 0.5
    real = SpinWeightedSpheroidalHarmonic(s, j, m, g).couplingcoefficients
    near = SpinWeightedSpheroidalHarmonic(s, j, m, g + 1e-10j).couplingcoefficients
    n = min(len(real), len(near))
    np.testing.assert_allclose(near[:n], real[:n], atol=1e-8)


def test_bilinear_normalize_convention():
    b = np.array([0.2 + 0.1j, 3.0 - 1.0j, 0.05j])
    out = _bilinear_normalize(b, 1)
    assert abs(np.sum(out**2) - 1.) < 1e-14
    assert out[1].real > 0.
    np.testing.assert_allclose(out/out[1], b/b[1], rtol=1e-14)


def test_bilinear_normalize_self_orthogonal_raises():
    # sum b^2 = 1 + (1j)^2 = 0: no bilinear normalization exists
    with pytest.raises(ValueError):
        _bilinear_normalize(np.array([1.0, 1.0j]), 0)
