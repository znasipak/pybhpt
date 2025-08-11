import pytest
from pybhpt.geo import KerrGeodesic
import numpy as np

# values generated with BlackHolePerturbationToolkit Mathematica module
test_data_frequencies = [
    ([0.8, 12.0, 0.0, 1.0], [0, 0, 0.023602039749475015]),
    ([0.8, 12.0, 0.0, -1.0], [0, 0, -0.024528308737989217]),
    ([0.99, 12.0, 0.5, 1.0], [0.01319448012413408, 0.015513550082646335, 0.01614608323014866]),
    ([0.01, 12.0, 0.5, -1.0], [0.012000588609775374, 0.01710391897431342, -0.01709526006998466]),
    ([0.1, 6.0, 0.0, 0.1], [0.005454600444974914, 0.06788263510563093, 0.06880530205262363]),
    ([0.95, 6.0, 0.9, 0.3], [0.004800849870503979, 0.010875592729535219, 0.01278659649394388]),
]

test_data_minotime = [
    ([0.9, 15.0, 0.0, 1.0], 12.4, np.array([3108.0241527325566, 15., 1.5707963267948966, 52.68306319729915])), # equatorial circular orbit
    ([0.95, 2.0, 0.0, 1.0], -12.32, np.array([-158.84975882874983, 2., 1.5707963267948966, -42.04123927345041])), # equatorial circular orbit
    ([0.95, 2.72, 0.77, 1.0], 20.9, np.array([555.748324063948, 1.5991882905840717, 1.5707963267948966, 94.7718354186187])), # equatorial eccentric orbit
    ([0.999, 6.2, 0.95, 0.1], 212.21, np.array([99593.50088382357, 3.2671554349116683, 1.6973651216220647, 973.2484697350069])), # generic
    ([0.87, 106.2, 0.99, -0.5], -2.32, np.array([-9.803436708010195e6, 96.84937691609196, 0.8642840935869669, 24.03833354082708])), # generic
    ([0.87, 10.2, 0.99, -0.5], -12.32, np.array([-294656.1821307974, 6.6148714120513326, 2.5723906230582694, 49.52880519014693])), # generic
    ([0.87, 10.2, 0.99, -0.5], -12.32, np.array([-294656.1821307974, 6.6148714120513326, 2.5723906230582694, 49.52880519014693])), # generic
    ([0.87, 8., 0., -0.5], 2, np.array([164.20161867254475, 8., 1.382158827830281, -7.1309347670154])), # spherical
]

test_parameters = [
    (-0.1, 10.0, 0.5, 0.5),  # invalid a
    (1.1, 10.0, 0.5, 0.5),  # invalid a
    (0.5, 10.0, 0.5, 0.5, 100),  # nsamples not a power of two
    (0.87, 6.4, 0.99, -0.5),  # plunging geodesic
]

@pytest.mark.parametrize("apex", test_parameters)
def test_kerr_geodesic_init(apex):
    with pytest.raises(ValueError):
        KerrGeodesic(*apex)

@pytest.mark.parametrize("apex, omega", test_data_frequencies)
def test_kerr_geodesic_frequencies(apex, omega):
    a, p, e, x, = apex
    orbit = KerrGeodesic(a, p, e, x)
    omega_eval = orbit.frequencies
    assert np.allclose(omega, omega_eval, rtol = 1e-9, atol = 0), f"Failed for apex {apex}. Expected {omega}, got {omega_eval} with rel error of {np.abs(1-omega/omega_eval)}"

@pytest.mark.parametrize("apex, la, xp", test_data_minotime)
def test_kerr_geodesic_minotime(apex, la, xp):
    a, p, e, x, = apex
    orbit = KerrGeodesic(a, p, e, x)
    xpeval = orbit(la)
    assert np.allclose(xpeval, xp, rtol = 1e-9, atol = 0), f"Failed for apex {apex} with la {la}. Expected {xp}, got {xpeval} with rel error of {np.abs(1-xp/xpeval)}"

