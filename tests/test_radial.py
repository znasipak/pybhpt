from pybhpt.radial import renormalized_angular_momentum, RadialTeukolsky
import pytest
import numpy as np

test_data_renormalized_angular_momentum = [
    ([-2, 5, 2, 0.99, 2.1], -2.45191080799969j),
    ([-2, 15, 2, 0.99, 2.1], 14.855326709178254),
    ([-2, 20, 2, 0.99, 6.5], -5.953247813925858j),
    ([-2, 30, 2, 0.99, 6.5], 29.78319684515244)
]

@pytest.mark.parametrize("params, expected", test_data_renormalized_angular_momentum)
def test_renormalized_angular_momentum(params, expected):
    nu = renormalized_angular_momentum(*params)
    assert np.isclose(np.cos(2.*np.pi*nu), np.cos(2.*np.pi*expected), rtol=1e-9, atol = 0), f"Expected {expected}, got {nu} for parameters {params}"

def test_radial_teukolsky_init():
    with pytest.raises(ValueError):
        RadialTeukolsky(0, 10, 2, 1.1, 0.1, [5.])  # invalid a
    with pytest.raises(ValueError):
        RadialTeukolsky(-1, 1, 2, 0.5, 0.1, [5.])  # invalid j
    with pytest.raises(ValueError):
        RadialTeukolsky(2, 2, 2, 0.5, 0.1, [1.])  # invalid r
    with pytest.raises(AttributeError):
        RadialTeukolsky(0.5, 2, 2, 0.5, 0.1, 10.) # invalid r type

def test_radial_teukolsky_solve():
    with pytest.raises(ValueError):
        Rt = RadialTeukolsky(0, 10, 2, 0.5, 0.1, [5.])
        Rt.solve("HYPERBOLIC", "None")  # invalid method
    with pytest.raises(ValueError):
        Rt = RadialTeukolsky(0.5, 2, 2, 0.5, 0.1, [5.])
        Rt.solve("AUTO", "Down") # invalid boundary condition

@pytest.mark.parametrize("s", [0, 1, -1])
def test_radial_teukolsky_gsn_unsupported_spin(s):
    Rt = RadialTeukolsky(s, 2, 2, 0.5, 2.0, np.linspace(3, 30, 8))
    with pytest.raises(ValueError):
        Rt.solve("GSN", "Up")

@pytest.mark.parametrize("s, j, m, a, omega", [
    (2, 2, 2, 0.99, 2.0),
    (2, 3, 3, 0.99, 2.0),
    (2, 2, 2, 0.9, 2.0),
    (2, 2, 0, 0.5, 2.0),
    (-2, 2, 2, 0.99, 2.0),
])
@pytest.mark.parametrize("bc", ["In", "Up"])
def test_radial_teukolsky_gsn_matches_auto(s, j, m, a, omega, bc):
    r = np.linspace(3, 30, 8)
    ref = RadialTeukolsky(s, j, m, a, omega, r)
    ref.solve("AUTO", bc)
    gsn = RadialTeukolsky(s, j, m, a, omega, r)
    gsn.solve("GSN", bc)
    np.testing.assert_allclose(gsn.radialsolutions(bc), ref.radialsolutions(bc), rtol=1e-7, atol=0)
    np.testing.assert_allclose(gsn.radialderivatives(bc), ref.radialderivatives(bc), rtol=1e-7, atol=0)

@pytest.mark.parametrize("method", ["AUTO", "MST", "ASYMP", "HBL", "GSN", "TEUK"])
def test_radial_teukolsky_solve_auto(method):
    r = np.linspace(2.1, 100, 100)
    Rt = RadialTeukolsky(-2, 2, 2, 0.5, 0.1, r)
    Rt.solve(method, "None")
    assert Rt.radialsolutions('Up').shape == r.shape, f"Radial Up solutions shape mismatch for {method} method"
    assert Rt.radialsolutions('In').shape == r.shape, f"Radial In solutions shape mismatch for {method} method"