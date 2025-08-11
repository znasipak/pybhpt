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

@pytest.mark.parametrize("method", ["AUTO", "MST", "ASYMP", "HBL", "GSN", "TEUK"])
def test_radial_teukolsky_solve_auto(method):
    r = np.linspace(2.1, 100, 100)
    Rt = RadialTeukolsky(-2, 2, 2, 0.5, 0.1, r)
    Rt.solve(method, "None")
    assert Rt.radialsolutions('Up').shape == r.shape, f"Radial Up solutions shape mismatch for {method} method"
    assert Rt.radialsolutions('In').shape == r.shape, f"Radial In solutions shape mismatch for {method} method"