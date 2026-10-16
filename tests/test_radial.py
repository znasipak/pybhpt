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

# Reference boundary values for a = 0, s = -2, l = 2, m = 2, omega = 0.3, computed with mpmath
# (60 digits) by optimally truncating the same series. The Up asymptotic series has c_5 = 0
# exactly for this mode, which previously stopped the summation early and pushed the Up
# boundary from r = 20/omega out to r = 533.
test_data_asymptotic_boundary = [
    ("Up", 20./0.3, -289930.7941591196-57469.09621284917j, 4617.701325123715-91814.33438341816j),
    ("In", 2.2, 0.11404129915797666+0.14104370771576444j, 1.556981634585168+1.2136522927301545j),
]

@pytest.mark.parametrize("bc, rb, R_ref, Rp_ref", test_data_asymptotic_boundary)
def test_radial_teukolsky_asymptotic_boundary(bc, rb, R_ref, Rp_ref):
    Rt = RadialTeukolsky(-2, 2, 2, 0., 0.3, np.array([4., 10.]))
    Rt.solveboundarycondition("ASYMP")
    assert np.isclose(Rt.boundarypoint(bc), rb, rtol=1e-14, atol=0), f"{bc} boundary moved to {Rt.boundarypoint(bc)}"
    assert np.isclose(Rt.boundarysolution(bc), R_ref, rtol=1e-12, atol=0)
    assert np.isclose(Rt.boundaryderivative(bc), Rp_ref, rtol=1e-12, atol=0)

def test_radial_teukolsky_auto_far_zone():
    # For omega*r > 100 AUTO tries the asymptotic series for both solutions. The horizon (In)
    # series diverges this far out, so AUTO must reject it and fall back to integration
    # rather than return a truncated divergent sum.
    r = np.linspace(500., 5000., 50)
    Rt = RadialTeukolsky(-2, 2, 2, 0.9, 0.3, r)
    Rt.solve("AUTO", "None")
    Rhbl = RadialTeukolsky(-2, 2, 2, 0.9, 0.3, r)
    Rhbl.solve("HBL", "None")
    for bc in ["In", "Up"]:
        assert np.allclose(Rt.radialsolutions(bc), Rhbl.radialsolutions(bc), rtol=1e-8, atol=0), f"AUTO {bc} solution disagrees with HBL"
        assert np.allclose(Rt.radialderivatives(bc), Rhbl.radialderivatives(bc), rtol=1e-8, atol=0), f"AUTO {bc} derivative disagrees with HBL"

def test_radial_teukolsky_near_extremal_horizon_boundary():
    # For a = 0.999 the default In boundary 1.1 r_+ lies outside the horizon series' radius of
    # convergence, |r - r_+| < 2 kappa. The series must reject it so the boundary moves inward.
    # Reference solutions were computed on main, where the boundary stayed at 1.1 r_+ and the
    # truncated sum happened to be accurate, and agree with TEUK integration to ~1e-13.
    a = 0.999
    kappa = np.sqrt(1. - a*a)
    rplus = 1. + kappa
    r = np.array([3., 10.])
    Rt = RadialTeukolsky(-2, 5, 2, a, 0.02, r)
    Rt.solve("AUTO", "In")
    assert rplus < Rt.boundarypoint("In") < rplus + 2.*kappa
    R_ref = np.array([-1591.3774212456742+503.99857568443264j, -59998484.632392526-2851215.822590392j])
    Rp_ref = np.array([-5610.473993908352+1414.5623325168212j, -46577097.17367911-3197014.1241712696j])
    assert np.allclose(Rt.radialsolutions("In"), R_ref, rtol=1e-10, atol=0)
    assert np.allclose(Rt.radialderivatives("In"), Rp_ref, rtol=1e-10, atol=0)
