"""Direct checks of the spin-weighted spheroidal harmonic grid evaluation.

The SpinWeightedHarmonic grid path (tabulated Ylm + collapsed coupling coefficients) is
exercised through TeukolskyMode, which exposes S, S', S'' on the polar grid via
polarsolutions / polarderivatives / polarderivatives2. We compare it against the
independent series-based SpinWeightedSpheroidalHarmonic (same spectral coupling, but the
straightforward per-point Sslm/Sslm_derivative sums), so the grid refactor is validated
against a different code path across spins, mode numbers, negative m, and orbit
geometries (inclined, equatorial, near-extremal).
"""
import numpy as np
import pytest

from pybhpt.geo import KerrGeodesic
from pybhpt.teuk import TeukolskyMode
from pybhpt.swsh import SpinWeightedSpheroidalHarmonic, SpinWeightedHarmonic

# (s, l, m, k, n, a, p, e, x)
CASES = [
    (-2, 5,  2, 2, 3, 0.9,  7.0, 0.6, 0.1),   # generic inclined, s = -2
    ( 2, 5,  2, 2, 3, 0.9,  7.0, 0.6, 0.1),   # generic inclined, s = +2
    ( 0, 4,  2, 2, 2, 0.9,  7.0, 0.6, 0.1),   # scalar, s = 0
    (-2, 6, -3, 1, 0, 0.9,  6.5, 0.3, 0.5),   # negative m, inclined
    (-2, 3,  3, 0, 5, 0.9,  5.0, 0.5, 1.0),   # equatorial (theta = pi/2 constant)
    (-2, 8,  1, -2, 2, 0.99, 6.5, 0.85, 0.05),  # near-extremal, strong field
]


def _relerr(approx, ref):
    """Max abs error normalised by the overall scale (robust where ref crosses zero)."""
    scale = np.max(np.abs(ref))
    if scale == 0.0:
        return np.max(np.abs(approx))
    return np.max(np.abs(approx - ref)) / scale


@pytest.mark.parametrize("s, l, m, k, n, a, p, e, x", CASES)
def test_swsh_grid_matches_series(s, l, m, k, n, a, p, e, x):
    orbit = KerrGeodesic(a, p, e, x, 2**7)
    mode = TeukolskyMode(s, l, m, k, n, orbit)
    mode.solve(orbit)

    th = np.asarray(orbit.polarpoints)
    gamma = a * orbit.mode_frequency(m, k, n)
    series = SpinWeightedSpheroidalHarmonic(s, l, m, gamma)

    S_grid = np.asarray(mode.polarsolutions)
    assert _relerr(S_grid, np.asarray(series.Sslm(th))) < 1e-8, "S(theta) grid vs series mismatch"

    # TeukolskyMode only populates the polar derivatives for |s| = 2 (the scalar source
    # integral does not need S', S''); check them where they are computed.
    if s != 0:
        SP_grid = np.asarray(mode.polarderivatives)
        SPP_grid = np.asarray(mode.polarderivatives2)
        assert _relerr(SP_grid, np.asarray(series.Sslm_derivative(th))) < 1e-8, "S'(theta) grid vs series mismatch"
        assert _relerr(SPP_grid, np.asarray(series.Sslm_derivative2(th))) < 1e-7, "S''(theta) grid vs series mismatch"


# (s, l, m, gamma) for the standalone grid class (independent of an orbit)
GRID_CASES = [
    (-2, 2,  1, 0.5),
    ( 0, 3,  2, 1.3),
    ( 2, 5, -3, 2.0),
    (-2, 4,  0, -0.8),
    ( 2, 2,  0, 0.0),   # spheroidicity 0 -> reduces to spherical harmonic
    ( 0, 2,  0, 0.0),
]


@pytest.mark.parametrize("s, l, m, gamma", GRID_CASES)
def test_swsh_grid_class_matches_series(s, l, m, gamma):
    """The grid-backed SpinWeightedHarmonic must agree with the series class on its stored
    grid, and its arbitrary-theta eval must reproduce the stored grid values."""
    th = np.linspace(0.05, np.pi - 0.05, 40)
    grid = SpinWeightedHarmonic(s, l, m, gamma, th)
    series = SpinWeightedSpheroidalHarmonic(s, l, m, gamma)

    assert abs(grid.eigenvalue - series.eigenvalue) < 1e-10
    assert _relerr(grid.solutions, np.asarray(series.Sslm(th))) < 1e-8, "grid vs series"
    # arbitrary-theta eval reuses the stored coupling -> must match the stored grid
    assert _relerr(grid.eval(th, deriv=0), grid.solutions) < 1e-10, "eval vs stored grid"

    assert _relerr(grid.derivatives, np.asarray(series.Sslm_derivative(th))) < 1e-8, "S' grid vs series"
    assert _relerr(grid.secondderivatives, np.asarray(series.Sslm_derivative2(th))) < 1e-7, "S'' grid vs series"
    assert _relerr(grid.eval(th, deriv=1), grid.derivatives) < 1e-10, "eval S' vs grid"
    assert _relerr(grid.eval(th, deriv=2), grid.secondderivatives) < 1e-8, "eval S'' vs grid"

    # scalar in -> scalar out
    assert np.ndim(grid.eval(0.7)) == 0


def test_swsh_series_array_matches_scalar():
    """The series eval must be consistent whether called on an array or point-by-point."""
    series = SpinWeightedSpheroidalHarmonic(-2, 6, 2, 0.9 * 0.3)
    th = np.linspace(0.15, np.pi - 0.15, 17)
    arr = np.asarray(series.Sslm(th))
    pt = np.array([series.Sslm(float(t)) for t in th])
    assert np.allclose(arr, pt, rtol=1e-12, atol=0.0)
