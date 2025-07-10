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
    return kerr_orbital_constants_wrapper(a, p, e, x)

class KerrGeodesic:
    """
    Class that produces a Kerr geodesic given the parameters of the orbit.
    This class is a wrapper around the Cython implementation of the Kerr geodesic
    and provides a Python interface to the underlying C++ code.
    Parameters
    ----------
    a : float
        The black hole spin parameter.
    p : float
        The semilatus rectum of the orbit.
    e : float
        The eccentricity of the orbit.
    x : float
        The inclination of the orbit.
    nsamples : int
        The number of samples to use for the geodesic. Default is 256.
    
    Attributes
    ----------
    blackholespin : float
        The black hole spin parameter.
    semilatusrectum : float
        The semilatus rectum of the orbit.
    eccentricity : float
        The eccentricity of the orbit.
    inclination : float
        The inclination of the orbit.
    orbitalenergy : float
        The orbital energy En of the orbit.
    orbitalangularmomentum : float
        Th z-component of the orbital angular momentum Lz of the orbit.
    carterconstant : float
        The Carter constant Qc of the orbit.
    orbitalconstants : numpy.ndarray
        The orbital constants (En, Lz, Qc) of the orbit.
    radialroots : numpy.ndarray
        The roots of the radial equation.
    polarroots : numpy.ndarray
        The roots of the polar equation.
    minofrequencies : numpy.ndarray
        The orbital frequencies with respect to Mino time.
    timefrequencies : numpy.ndarray
        The orbital frequencies with respect to the time coordinate.
    frequencies : numpy.ndarray
        The (coordinate time) frequencies of the orbit.
    carterfrequencies : numpy.ndarray
        The frequencies for computing Carter constant fluxes.
    timeradialfourier : numpy.ndarray
        The Fourier coefficients of coordinate time with respect to the radial Mino phase.
    timepolarfourier : numpy.ndarray
        The Fourier coefficients of coordinate time with respect to the polar Mino phase.
    radialfourier : numpy.ndarray
        The Fourier coefficients of radial position with respect to the radial Mino phase.
    polarfourier : numpy.ndarray
        The Fourier coefficients of polar position with respect to the polar Mino phase.
    azimuthalradialfourier : numpy.ndarray
        The Fourier coefficients of azimuthal position with respect to the radial Mino phase.
    azimuthalpolarfourier : numpy.ndarray
        The Fourier coefficients of azimuthal position with respect to the polar Mino phase.  
    """
    def __init__(self, a, p, e, x, nsamples = 2**8):
        self.base = KerrGeodesicCython(a, p, e, x, nsamples)
        """The base class that contains the Cython implementation of the Kerr geodesic."""
        self.timeradial = self.base.get_time_accumulation(1)
        self.timepolar = self.base.get_time_accumulation(2)
        self.radialpoints = self.base.get_radial_points()
        self.polarpoints = self.base.get_polar_points()
        self.azimuthalradial = self.base.get_azimuthal_accumulation(1)
        self.azimuthalpolar = self.base.get_azimuthal_accumulation(2)
        self.nsamples = self.base.nsamples

    @property
    def blackholespin(self):
        """
        The black hole spin parameter.
        """
        return self.base.blackholespin
    
    @property
    def semilatusrectum(self):
        """
        The semilatus rectum of the orbit.
        """
        return self.base.semilatusrectum
    
    @property
    def eccentricity(self):
        """
        The eccentricity of the orbit.
        """
        return self.base.eccentricity

    @property
    def inclination(self):
        """
        The cosine of the inclination angle of the orbit with respect to the equatorial plane.
        """
        return self.base.inclination

    @property
    def orbitalenergy(self):
        """
        The orbital energy En of the orbit.
        """
        return self.base.orbitalenergy

    @property
    def orbitalangularmomentum(self):
        """
        The z-component of the orbital angular momentum Lz of the orbit.
        """
        return self.base.orbitalangularmomentum

    @property
    def carterconstant(self):
        """
        The Carter constant Qc of the orbit.
        """
        return self.base.carterconstant
    
    @property
    def orbitalconstants(self):
        """
        The orbital constants (En, Lz, Qc) of the orbit.
        """
        return np.array([self.orbitalenergy, self.orbitalangularmomentum, self.carterconstant])

    @property
    def radialroots(self):
        """
        The roots of the radial equation.
        """
        return self.base.radialroots

    @property
    def polarroots(self):
        """
        The roots of the polar equation.
        """
        return self.base.polarroots

    @property
    def minofrequencies(self):
        """
        The orbital frequencies with respect to Mino time.
        """
        return self.base.minofrequencies

    @property
    def timefrequencies(self):
        """
        The orbital frequencies with respect to the time coordinate.
        """
        return self.base.timefrequencies

    @property
    def frequencies(self):
        """
        The (coordinate time) frequencies of the orbit.
        """
        return self.timefrequencies

    @property
    def carterfrequencies(self):
        """
        The frequencies for computing Carter constant fluxes.
        """
        return self.base.carterfrequencies
    
    @property
    def timeradialfourier(self):
        """
        The Fourier coefficients of coordinate time with respect to the radial Mino phase.
        """
        return self.base.get_time_coefficients(1)
    
    @property
    def timepolarfourier(self):
        """
        The Fourier coefficients of coordinate time with respect to the polar Mino phase.
        """
        return self.base.get_time_coefficients(2)
    
    @property
    def radialfourier(self):
        """
        The Fourier coefficients of radial position with respect to the radial Mino phase.
        """
        return self.base.get_radial_coefficients()
    
    @property
    def polarfourier(self):
        """
        The Fourier coefficients of polar position with respect to the polar Mino phase.
        """
        return self.base.get_polar_coefficients()
    
    @property
    def azimuthalradialfourier(self):
        """
        The Fourier coefficients of azimuthal position with respect to the radial Mino phase.
        """
        return self.base.get_azimuthal_coefficients(1)
    
    @property
    def azimuthalpolarfourier(self):
        """
        The Fourier coefficients of azimuthal position with respect to the polar Mino phase.
        """
        return self.base.get_azimuthal_coefficients(2)
    
    def mode_frequency(self, m, k, n):
        """
        Returns the frequency of the mode with azimuthal number m, polar number k, and radial number n.
        
        Parameters
        ----------
        m : int
            The azimuthal number of the mode.
        k : int
            The polar number of the mode.
        n : int
            The radial number of the mode.
        
        Returns
        -------
        float
            The frequency of the mode with azimuthal number m, polar number k, and radial number n.
        """
        return np.dot(np.array([n, k, m]), (self.frequencies))
    
    def minotime(self, t):
        """
        Function that returns the Mino time for a given Boyer-Lindquist time.

        Parameters
        ----------
        t : float
            The Boyer-Lindquist time.
        
        Returns
        -------
        float
            The Mino time for the given Boyer-Lindquist time.
        """
        return self.base.mino_time(t)
    
    def __call__(self, la):
        """
        Function that returns the position vector for a given Mino time value.
        Parameters
        ----------
        la : float or numpy.ndarray
            The Mino time value(s). If a numpy array is provided, the function will return a numpy array of the same shape.
        Returns 
        -------
        numpy.ndarray
            The position vector for the given Mino time value(s).
        """
        if isinstance(la, np.ndarray) or isinstance(la, list):
            return self.base.position_vec(np.array(la))
        else:
            return self.base.position(la)