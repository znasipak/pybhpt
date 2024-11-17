from cybhpt_full import KerrGeodesic as KerrGeodesicCython
from cybhpt_full import kerr_geo_V01, kerr_geo_V02, kerr_geo_V11, kerr_geo_V22, kerr_geo_V31, kerr_geo_V32
from cybhpt_full import kerr_mino_frequencies_wrapper, kerr_orbital_constants_wrapper
import numpy as np

def kerrgeo_Vt_radial(a, En, Lz, Q, r):
    return kerr_geo_V01(a, En, Lz, Q, r)

def kerrgeo_Vt_polar(a, En, Lz, Q, theta):
    return kerr_geo_V02(a, En, Lz, Q, theta)

def kerrgeo_Vr(a, En, Lz, Q, r):
    return kerr_geo_V11(a, En, Lz, Q, r)

def kerrgeo_Vtheta(a, En, Lz, Q, theta):
    return kerr_geo_V22(a, En, Lz, Q, theta)

def kerrgeo_Vphi_radial(a, En, Lz, Q, r):
    return kerr_geo_V31(a, En, Lz, Q, r)

def kerrgeo_Vphi_polar(a, En, Lz, Q, theta):
    return kerr_geo_V32(a, En, Lz, Q, theta)

def kerrgeo_mino_frequencies(a, p, e, x):
    return kerr_mino_frequencies_wrapper(a, p, e, x)

def kerr_orbital_constants(a, p, e, x):
    kerr_orbital_constants_wrapper(a, p, e, x)

class KerrGeodesic:
    def __init__(self, a, p, e, x, nsamples = 2**8):
        self.base = KerrGeodesicCython(a, p, e, x, nsamples)
        self.timeradial = self.base.get_time_accumulation(1)
        self.timepolar = self.base.get_time_accumulation(2)
        self.radialpoints = self.base.get_radial_points()
        self.polarpoints = self.base.get_polar_points()
        self.azimuthalradial = self.base.get_azimuthal_accumulation(1)
        self.azimuthalpolar = self.base.get_azimuthal_accumulation(2)
        self.nsamples = self.base.nsamples

    @property
    def blackholespin(self):
        return self.base.blackholespin
    
    @property
    def semilatusrectum(self):
        return self.base.semilatusrectum
    
    @property
    def eccentricity(self):
        return self.base.eccentricity

    @property
    def inclination(self):
        return self.base.inclination

    @property
    def orbitalenergy(self):
        return self.base.orbitalenergy

    @property
    def orbitalangularmomentum(self):
        return self.base.orbitalangularmomentum

    @property
    def carterconstant(self):
        return self.base.carterconstant
    
    @property
    def orbitalconstants(self):
        return np.array([self.orbitalenergy, self.orbitalangularmomentum, self.carterconstant])

    @property
    def radialroots(self):
        return self.base.radialroots

    @property
    def polarroots(self):
        return self.base.polarroots

    @property
    def minofrequencies(self):
        return self.base.minofrequencies

    @property
    def timefrequencies(self):
        return self.base.timefrequencies

    @property
    def frequencies(self):
        return self.timefrequencies

    @property
    def carterfrequencies(self):
        return self.base.carterfrequencies
    
    @property
    def timeradialfourier(self):
        return self.base.get_time_coefficients(1)
    
    @property
    def timepolarfourier(self):
        return self.base.get_time_coefficients(2)
    
    @property
    def radialfourier(self):
        return self.base.get_radial_coefficients()
    
    @property
    def polarfourier(self):
        return self.base.get_polar_coefficients()
    
    @property
    def azimuthalradialfourier(self):
        return self.base.get_azimuthal_coefficients(1)
    
    @property
    def azimuthalpolarfourier(self):
        return self.base.get_azimuthal_coefficients(2)
    
    def mode_frequency(self, m, k, n):
        return np.dot(np.array([n, k, m]), (self.frequencies))
    
    def minotime(self, t):
        return self.base.mino_time(t)
    
    def __call__(self, la):
        if isinstance(la, np.ndarray) or isinstance(la, list):
            return self.base.position_vec(np.array(la))
        else:
            return self.base.position(la)