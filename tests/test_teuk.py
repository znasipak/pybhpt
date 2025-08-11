from pybhpt.geo import KerrGeodesic
from pybhpt.teuk import TeukolskyMode
import pytest
import numpy as np

test_data_teuk_amplitudes = [
    # circular orbits
    ([0.9, 10., 0., 1., 2**2], [-2, 5, 2, 0, 0], 
     [complex(2.270392541324341e-7,2.7044100708956895e-7), 
      complex(-4.631392141748372e-8,2.133802344473894e-8)]
    ),
    ([0.9, 5, 0, 1, 2**2], [0, 6, 2, 0, 0], 
     [complex(9.801251306590833e-9,-8.378596245546433e-8), 
      complex(-1.7883676441395782e-6,2.247615864437883e-6)]
    ),
    ([0.9, 5, 0, 1, 2**2], [2, 3, 3, 0, 0], 
     [complex(0.02110394906733878,0.03593604102282673), 
      complex(-36.11531169953736,-43.391859358161284)]
    ),
    # eccentric equatorial orbits
    ([0.9, 5., 0.5, 1., 2**9], [-2, 3, 3, 0, 1], 
     [complex(-0.00005383220891991058,-0.0006147497965404015), 
      complex(0.0009114873055821887,0.0011872334313407199)]
    ),
    ([0.9, 5, 0.5, 1, 2**9], [0, 3, 1, 0, 5], 
     [complex(0.000025591092570891442,0.00002724728967365839), 
      complex(0.00002860344315787682,0.000029752874170139442)]
    ),
    ([0.65, 5.3, 0.7, 1, 2**9], [2, 13, 10, 0, 20], 
     [complex(-5.849228222663369e-12,-1.1285520078506349e-12), 
      complex(0.00036841567195881875,0.0001068390464888473)]
    ),
    # spherical orbits
    ([0.9, 7., 0., 0.1, 2**9], [-2, 3, 3, 4, 0], 
     [complex(1.1271884617178083e-9,9.482392163824508e-8), 
      complex(-6.696690048496821e-7,-1.0029341164657833e-6)]
    ),
    ([0.9, 6.5, 0., 0.5, 2**9], [0, 13, 2, -5, 0], 
     [complex(4.9032544335101e-18,9.394081304920472e-19), 
      complex(5.891851988551972e-16,-3.029620642780757e-16)]
    ),
    ([0.9, 5.3, 0., 0.05, 2**9], [2, 13, 10, -3, 0], 
     [complex(4.3363673814080184e-10,4.1293933462392084e-10), 
      complex(0.0003094972762916496,-0.0005739959201115489)]
    ),
    # generic orbits
    ([0.9, 7., 0.6, 0.1, 2**10], [-2, 3, 3, 4, -5], 
     [complex(-2.3111401512947867e-11,-2.941719101203449e-10), 
      complex(2.125582269356303e-9,2.984152461796948e-9)]
    ),
    ([0.9, 6.5, 0.3, 0.5, 2**10], [0, 4, 2, -4, 2], 
     [complex(1.1218341814508144e-7,3.466068165183326e-7), 
      complex(-1.476696172619881e-6,-3.7630253125367686e-7)]
    ),
    ([0.99, 6.5, 0.85, 0.05, 2**12], [2, 13, 10, -3, 20], 
     [complex(2.5098536599898154e-13,-5.913850827127284e-14), 
      complex(-5.50126823211977e-6,1.0551842131079352e-6)]
    ),
]

@pytest.mark.parametrize("apex, mode_params, expected_amplitude", test_data_teuk_amplitudes)
def test_teukolsky_mode_in_amplitude(apex, mode_params, expected_amplitude):
    orbit = KerrGeodesic(*apex)
    mode = TeukolskyMode(*mode_params, orbit)
    mode.solve(orbit)
    rtol = mode.precision('In')
    # TODO: Adjust the rtol based on the expected precision of the amplitude, but still unsure of accuracy of toolkit data
    # For now, we set a minimum rtol to avoid too strict comparisons
    if rtol < 1e-5:
        rtol = 1e-5
    assert np.isclose(mode.amplitude('In'), expected_amplitude[0], rtol=rtol, atol=0), f"In amplitude of {mode.amplitude('In')} with error of {np.abs(1-mode.amplitude('In')/expected_amplitude[0])} for orbit {apex} and parameters {mode_params} with mode precision {rtol}"

@pytest.mark.parametrize("apex, mode_params, expected_amplitude", test_data_teuk_amplitudes)
def test_teukolsky_mode_up_amplitude(apex, mode_params, expected_amplitude):
    orbit = KerrGeodesic(*apex)
    mode = TeukolskyMode(*mode_params, orbit)
    mode.solve(orbit)
    rtol = mode.precision('Up')
    if rtol < 1e-5:
        rtol = 1e-5
    assert np.isclose(mode.amplitude('Up'), expected_amplitude[1], rtol=rtol, atol=0), f"Up amplitude of {mode.amplitude('Up')} with error of {np.abs(1-mode.amplitude('Up')/expected_amplitude[1])} for orbit {apex} and parameters {mode_params} with mode precision {rtol}"
