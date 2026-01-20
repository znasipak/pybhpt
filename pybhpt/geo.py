from cybhpt_full import KerrGeodesic as _KerrGeodesicCython
from cybhpt_full import _kerr_geo_V01, _kerr_geo_V02, _kerr_geo_V11, _kerr_geo_V22, _kerr_geo_V31, _kerr_geo_V32
from cybhpt_full import _kerr_mino_frequencies_wrapper, _kerr_orbital_constants_wrapper, _kerr_orbital_constants_array_wrapper, _kerr_kepler_parameters_wrapper, _kerr_kepler_parameters_array_wrapper
from cybhpt_full import _jacobian_ELQ_to_pex_wrapper, _jacobian_ELQ_to_pex_array_wrapper, _jacobian_pex_to_ELQ_wrapper, _jacobian_pex_to_ELQ_array_wrapper
import numpy as np
import numpy.typing as npt
from typing import Union

# Define a type alias for convenience
NumericArray = npt.NDArray[np.float64]
ScalarOrArray = Union[float, NumericArray]

def kerrgeo_Vt_radial(a: float, En: float, Lz: float, Q: float, r: float) -> float:
    """
    The radial part of the potential Vt for the geodesic evolution of $t_p$.

    Parameters
    ----------
    a : float
        The black hole spin parameter.
    En : float
        The orbital energy of the particle.
    Lz : float
        The z-component of the orbital angular momentum of the particle.
    Q : float
        The Carter constant of the particle.
    r : float
        The radial coordinate.
    
    Returns
    -------
    float
        The value of the radial potential Vt at the given parameters.
    """
    return _kerr_geo_V01(a, En, Lz, Q, r)

def kerrgeo_Vt_polar(a: float, En: float, Lz: float, Q: float, theta: float) -> float:
    """
    The polar part of the potential Vt for the geodesic evolution of $t_p$.

    Parameters
    ----------
    a : float
        The black hole spin parameter.
    En : float
        The orbital energy of the particle.
    Lz : float
        The z-component of the orbital angular momentum of the particle.
    Q : float
        The Carter constant of the particle.
    theta : float
        The polar angle coordinate.

    Returns
    -------
    float
        The value of the polar potential Vt at the given parameters.
    """
    return _kerr_geo_V02(a, En, Lz, Q, theta)

def kerrgeo_Vr(a: float, En: float, Lz: float, Q: float, r: float) -> float:
    """
    The (squared) radial potential Vr for the geodesic evolution of $r_p$.

    Parameters
    ----------
    a : float
        The black hole spin parameter.
    En : float
        The orbital energy of the particle.
    Lz : float
        The z-component of the orbital angular momentum of the particle.
    Q : float
        The Carter constant of the particle.
    r : float
        The radial coordinate.

    Returns
    -------
    float
        The value of the radial potential Vr at the given parameters.
    """
    return _kerr_geo_V11(a, En, Lz, Q, r)

def kerrgeo_Vtheta(a: float, En: float, Lz: float, Q: float, theta: float) -> float:
    """
    The (squared) polar potential Vtheta for the geodesic evolution of $\theta_p$.

    Parameters
    ----------
    a : float
        The black hole spin parameter.
    En : float
        The orbital energy of the particle.
    Lz : float
        The z-component of the orbital angular momentum of the particle.
    Q : float
        The Carter constant of the particle.
    theta : float
        The polar angle coordinate.

    Returns
    -------
    float
        The value of the polar potential Vtheta at the given parameters.
    """
    return _kerr_geo_V22(a, En, Lz, Q, theta)

def kerrgeo_Vphi_radial(a: float, En: float, Lz: float, Q: float, r: float) -> float:
    """
    The (squared) radial potential Vphi for the geodesic evolution of $r_p$.

    Parameters
    ----------
    a : float
        The black hole spin parameter.
    En : float
        The orbital energy of the particle.
    Lz : float
        The z-component of the orbital angular momentum of the particle.
    Q : float
        The Carter constant of the particle.
    r : float
        The radial coordinate.

    Returns
    -------
    float
        The value of the radial potential Vphi at the given parameters.
    """
    return _kerr_geo_V31(a, En, Lz, Q, r)

def kerrgeo_Vphi_polar(a: float, En: float, Lz: float, Q: float, theta: float) -> float:
    """
    The (squared) polar potential Vphi for the geodesic evolution of $r_p$.

    Parameters
    ----------
    a : float
        The black hole spin parameter.
    En : float
        The orbital energy of the particle.
    Lz : float
        The z-component of the orbital angular momentum of the particle.
    Q : float
        The Carter constant of the particle.
    theta : float
        The polar angle coordinate.

    Returns
    -------
    float
        The value of the polar potential Vphi at the given parameters.
    """
    return _kerr_geo_V32(a, En, Lz, Q, theta)

def kerr_mino_frequencies(a: float, p: float, e: float, x: float) -> np.ndarray:
    """
    Returns the Mino frequencies of a Kerr geodesic.

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
    
    Returns
    -------
    numpy.ndarray
        The Mino time frequencies of the orbit.
    """
    return _kerr_mino_frequencies_wrapper(a, p, e, x)

def kerr_fundamental_frequencies(a: float, p: float, e: float, x: float) -> NumericArray:
    """
    Returns the fundamental (time) frequencies of a Kerr geodesic.
    
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
    
    Returns
    -------
    numpy.ndarray
        The fundamental frequencies of the orbit.
    """
    Ups = _kerr_mino_frequencies_wrapper(a, p, e, x)
    return Ups[1:]/Ups[0]

def kerr_orbital_constants(a: ScalarOrArray, 
                           p: ScalarOrArray, 
                           e: ScalarOrArray, 
                           x: ScalarOrArray) -> NumericArray:
    """
    Returns the orbital constants of a Kerr geodesic (En, Lz, Qc).

    Parameters
    ----------
    a : float or numpy.ndarray
        The black hole spin parameter.
    p : float or numpy.ndarray
        The semilatus rectum of the orbit.
    e : float or numpy.ndarray
        The eccentricity of the orbit.
    x : float or numpy.ndarray
        The inclination of the orbit. 

    Returns
    -------
    numpy.ndarray
        The orbital constants (En, Lz, Qc) of the orbit.
    """
    if isinstance(p, (np.ndarray, list)):
        assert isinstance(e, (np.ndarray, list)) and isinstance(x, (np.ndarray, list)), "If p is a numpy array, e and x must also be numpy arrays."
        p = np.asarray(p, dtype=np.float64)
        e = np.asarray(e, dtype=np.float64)
        x = np.asarray(x, dtype=np.float64)
        assert p.shape == e.shape == x.shape, "If p, e, and x are numpy arrays, they must have the same shape."
        if isinstance(a, (np.ndarray, list)):
            a = np.asarray(a, dtype=np.float64)
            assert a.shape == p.shape, "If a is a numpy array, it must have the same shape as p, e, and x."
        else:
            a = np.full(p.shape, a)
        return _kerr_orbital_constants_array_wrapper(a, p, e, x)
    else:
        return _kerr_orbital_constants_wrapper(a, p, e, x)

def kerr_kepler_parameters(a: ScalarOrArray, 
                           En: ScalarOrArray, 
                           Lz: ScalarOrArray, 
                           Qc: ScalarOrArray) -> NumericArray:
    """
    Returns the Keplerian parameters of a Kerr geodesic (p, e, x) given the orbital constants (a, En, Lz, Qc).

    Parameters
    ----------
    a : float or numpy.ndarray
        The black hole spin parameter.
    En : float or numpy.ndarray
        The orbital energy of the particle.
    Lz : float or numpy.ndarray
        The z-component of the orbital angular momentum of the particle.
    Qc : float or numpy.ndarray
        The Carter constant of the particle.

    Returns
    -------
    numpy.ndarray
        The Keplerian parameters (p, e, x) of the orbit.
    """
    if isinstance(En, (np.ndarray, list)):
        assert isinstance(Lz, (np.ndarray, list)) and isinstance(Qc, (np.ndarray, list)), "If En is a numpy array, Lz and Qc must also be numpy arrays."
        En = np.asarray(En, dtype=np.float64)
        Lz = np.asarray(Lz, dtype=np.float64)
        Qc = np.asarray(Qc, dtype=np.float64)
        assert En.shape == Lz.shape == Qc.shape, "If En, Lz, and Qc are numpy arrays, they must have the same shape."
        if isinstance(a, (np.ndarray, list)):
            a = np.asarray(a, dtype=np.float64)
            assert a.shape == En.shape, "If a is a numpy array, it must have the same shape as En, Lz, and Qc."
        else:
            a = np.full(En.shape, a)
        return _kerr_kepler_parameters_array_wrapper(a, En, Lz, Qc)
    else:
        return _kerr_kepler_parameters_wrapper(a, En, Lz, Qc)
    
def jacobian_ELQ_to_pex(a: ScalarOrArray,
                        p: ScalarOrArray,
                        e: ScalarOrArray,
                        x: ScalarOrArray) -> NumericArray:
    """
    Returns the Jacobian matrix of the transformation from (E, Lz, Q) to (p, e, x).
    |   dp/dE  dp/dLz  dp/dQ    |
    |   de/dE  de/dLz  de/dQ    |
    |   dx/dE  dx/dLz  dx/dQ    |

    Parameters
    ----------
    a : float or numpy.ndarray
        The black hole spin parameter.
    p : float or numpy.ndarray
        The semilatus rectum of the orbit.
    e : float or numpy.ndarray
        The eccentricity of the orbit.
    x : float or numpy.ndarray
        The inclination of the orbit.

    Returns
    -------
    numpy.ndarray
        The Jacobian matrix of the transformation from (E, Lz, Q) to (p, e, x).
    """
    if isinstance(p, (np.ndarray, list)):
        assert isinstance(e, (np.ndarray, list)) and isinstance(x, (np.ndarray, list)), "If p is a numpy array, e and x must also be numpy arrays."
        p = np.asarray(p, dtype=np.float64)
        e = np.asarray(e, dtype=np.float64)
        x = np.asarray(x, dtype=np.float64)
        assert p.shape == e.shape == x.shape, "If p, e, and x are numpy arrays, they must have the same shape."
        if isinstance(a, (np.ndarray, list)):
            a = np.asarray(a, dtype=np.float64)
            assert a.shape == p.shape, "If a is a numpy array, it must have the same shape as p, e, and x."
        else:
            a = np.full(p.shape, a)
        return _jacobian_ELQ_to_pex_array_wrapper(a, p, e, x)
    else:
        return _jacobian_ELQ_to_pex_wrapper(a, p, e, x)
    
def jacobian_pex_to_ELQ(a: ScalarOrArray,
                        p: ScalarOrArray,
                        e: ScalarOrArray,
                        x: ScalarOrArray) -> NumericArray:
    """
    Returns the Jacobian matrix of the transformation from (p, e, x) to (E, Lz, Q):
    |   dE/dp  dE/de  dE/dx    |
    |   dLz/dp dLz/de dLz/dx   |
    |   dQ/dp  dQ/de  dQ/dx    |

    Parameters
    ----------
    a : float or numpy.ndarray
        The black hole spin parameter.
    p : float or numpy.ndarray
        The semilatus rectum of the orbit.
    e : float or numpy.ndarray
        The eccentricity of the orbit.
    x : float or numpy.ndarray
        The inclination of the orbit.

    Returns
    -------
    numpy.ndarray
        The Jacobian matrix of the transformation from (p, e, x) to (E, Lz, Q).
    """
    if isinstance(p, (np.ndarray, list)):
        assert isinstance(e, (np.ndarray, list)) and isinstance(x, (np.ndarray, list)), "If p is a numpy array, e and x must also be numpy arrays."
        p = np.asarray(p, dtype=np.float64)
        e = np.asarray(e, dtype=np.float64)
        x = np.asarray(x, dtype=np.float64)
        assert p.shape == e.shape == x.shape, "If p, e, and x are numpy arrays, they must have the same shape."
        if isinstance(a, (np.ndarray, list)):
            a = np.asarray(a, dtype=np.float64)
            assert a.shape == p.shape, "If a is a numpy array, it must have the same shape as p, e, and x."
        else:
            a = np.full(p.shape, a)
        return _jacobian_pex_to_ELQ_array_wrapper(a, p, e, x)
    else:
        return _jacobian_pex_to_ELQ_wrapper(a, p, e, x)

def is_power_of_two(n: int) -> bool:
    return n > 0 and (n & (n - 1)) == 0

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
        The number of samples to use for the geodesic. Must be a power of two. Default is 256.

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
        if a < 0 or a > 1:
            raise ValueError(f"Black hole spin parameter {a} must be in the range [0, 1].")
        if not is_power_of_two(nsamples):
            raise ValueError(f"Number of samples {nsamples} must be a power of 2.")

        self.base = _KerrGeodesicCython(a, p, e, x, nsamples)
        """The base class that contains the Cython implementation of the Kerr geodesic."""
        self.timeradial = self.base.get_time_accumulation(1)
        self.timepolar = self.base.get_time_accumulation(2)
        self.radialpoints = self.base.get_radial_points()
        self.polarpoints = self.base.get_polar_points()
        self.azimuthalradial = self.base.get_azimuthal_accumulation(1)
        self.azimuthalpolar = self.base.get_azimuthal_accumulation(2)
        self.nsamples = self.base.nsamples

        if np.isnan(self.frequencies).any():
            raise ValueError(f"Orbital parameters (a, p, e, x) = {self.apex} do not represent a valid bound non-plunging orbit.")

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
    def a(self):
        """
        The black hole spin parameter.
        """
        return self.blackholespin
    
    @property
    def p(self):
        """
        The semilatus rectum of the orbit.
        """
        return self.semilatusrectum
    
    @property
    def e(self):
        """
        The eccentricity of the orbit.
        """
        return self.eccentricity
    
    @property
    def x(self):
        """
        The inclination of the orbit.
        """
        return self.inclination
    
    @property
    def apex(self):
        """
        The parameters of the orbit (a, p, e, x).
        """
        return np.array([self.blackholespin, self.semilatusrectum, self.eccentricity, self.inclination])
    
    @property
    def orbitalparameters(self):
        """
        The parameters of the orbit (a, p, e, x).
        """
        return self.apex

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
    
    def psi_radial(self, la):
        """
        Function that returns the radial true anomaly for a given Mino time value.

        Parameters
        ----------
        la : float or numpy.ndarray
            The Mino time value(s). If a numpy array is provided, the function will return a numpy array of the same shape.
        
        Returns 
        -------
        numpy.ndarray
            The radial true anomaly for the given Mino time value(s).
        """
        if isinstance(la, np.ndarray) or isinstance(la, list):
            return self.base.psi_radial_vec(np.array(la))
        else:
            return self.base.psi_radial(la)
        
    def psi_polar(self, la):
        """
        Function that returns the polar Darwin phase for a given Mino time value.

        Parameters
        ----------
        la : float or numpy.ndarray
            The Mino time value(s). If a numpy array is provided, the function will return a numpy array of the same shape.
        
        Returns 
        -------
        numpy.ndarray
            The polar Darwin phase for the given Mino time value(s).
        """
        if isinstance(la, np.ndarray) or isinstance(la, list):
            return self.base.psi_polar_vec(np.array(la))
        else:
            return self.base.psi_polar(la)
        
    def psi_radial_of_t(self, t):
        """
        Function that returns the radial true anomaly for a given time value.

        Parameters
        ----------
        t : float or numpy.ndarray
            The time value(s). If a numpy array is provided, the function will return a numpy array of the same shape.
        
        Returns 
        -------
        numpy.ndarray
            The radial true anomaly for the given time value(s).
        """
        if isinstance(t, np.ndarray) or isinstance(t, list):
            return self.base.psi_radial_time_vec(np.array(t))
        else:
            return self.base.psi_radial_time(t)

    def psi_polar_of_t(self, t):
        """
        Function that returns the polar Darwin phase for a given time value.

        Parameters
        ----------
        t : float or numpy.ndarray
            The time value(s). If a numpy array is provided, the function will return a numpy array of the same shape.
        
        Returns 
        -------
        numpy.ndarray
            The polar Darwin phase for the given time value(s).
        """
        if isinstance(t, np.ndarray) or isinstance(t, list):
            return self.base.psi_polar_time_vec(np.array(t))
        else:
            return self.base.psi_polar_time(t)

    def time_position(self, la):
        """
        Function that returns the time position for a given Mino time value.

        Parameters
        ----------
        la : float or numpy.ndarray
            The Mino time value(s). If a numpy array is provided, the function will return a numpy array of the same shape.
        
        Returns 
        -------
        numpy.ndarray
            The time position for the given Mino time value(s).
        """
        if isinstance(la, np.ndarray) or isinstance(la, list):
            return self.base.time_position_vec(np.array(la))
        else:
            return self.base.time_position(la)

    def radial_position(self, la):
        """
        Function that returns the radial position for a given Mino time value.

        Parameters
        ----------
        la : float or numpy.ndarray
            The Mino time value(s). If a numpy array is provided, the function will return a numpy array of the same shape.
        
        Returns 
        -------
        numpy.ndarray
            The radial position for the given Mino time value(s).
        """
        if isinstance(la, np.ndarray) or isinstance(la, list):
            return self.base.radial_position_vec(np.array(la))
        else:
            return self.base.radial_position(la)
        
    def polar_position(self, la):
        """
        Function that returns the polar position for a given Mino time value.

        Parameters
        ----------
        la : float or numpy.ndarray
            The Mino time value(s). If a numpy array is provided, the function will return a numpy array of the same shape.
        
        Returns 
        -------
        numpy.ndarray
            The polar position for the given Mino time value(s).
        """
        if isinstance(la, np.ndarray) or isinstance(la, list):
            return self.base.polar_position_vec(np.array(la))
        else:
            return self.base.polar_position(la)
    
    def azimuthal_position(self, la):
        """
        Function that returns the azimuthal position for a given Mino time value.

        Parameters
        ----------
        la : float or numpy.ndarray
            The Mino time value(s). If a numpy array is provided, the function will return a numpy array of the same shape.
        
        Returns 
        -------
        numpy.ndarray
            The azimuthal position for the given Mino time value(s).
        """
        if isinstance(la, np.ndarray) or isinstance(la, list):
            return self.base.azimuthal_position_vec(np.array(la))
        else:
            return self.base.azimuthal_position(la)
    
    def radial_position_of_t(self, t):
        """
        Function that returns the radial position for a given Mino time value.

        Parameters
        ----------
        t : float or numpy.ndarray
            The Mino time value(s). If a numpy array is provided, the function will return a numpy array of the same shape.
        
        Returns 
        -------
        numpy.ndarray
            The radial position for the given Mino time value(s).
        """
        if isinstance(t, np.ndarray) or isinstance(t, list):
            return self.base.radial_position_time_vec(np.array(t))
        else:
            return self.base.radial_position_time(t)
        
    def polar_position_of_t(self, t):
        """
        Function that returns the polar position for a given time value.

        Parameters
        ----------
        t : float or numpy.ndarray
            The time value(s). If a numpy array is provided, the function will return a numpy array of the same shape.
        
        Returns 
        -------
        numpy.ndarray
            The polar position for the given time value(s).
        """
        if isinstance(t, np.ndarray) or isinstance(t, list):
            return self.base.polar_position_time_vec(np.array(t))
        else:
            return self.base.polar_position_time(t)
    
    def azimuthal_position_of_t(self, t):
        """
        Function that returns the azimuthal position for a given time value.

        Parameters
        ----------
        t : float or numpy.ndarray
            The time value(s). If a numpy array is provided, the function will return a numpy array of the same shape.
        
        Returns 
        -------
        numpy.ndarray
            The azimuthal position for the given time value(s).
        """
        if isinstance(t, np.ndarray) or isinstance(t, list):
            return self.base.azimuthal_position_time_vec(np.array(t))
        else:
            return self.base.azimuthal_position_time(t)
    
    def mino_of_t(self, t):
        """
        Function that returns the Mino time for a given Boyer-Lindquist time.

        Parameters
        ----------
        t : float or numpy.ndarray
            The Boyer-Lindquist time value(s). If a numpy array is provided, the function will return a numpy array of the same shape.
        
        Returns 
        -------
        numpy.ndarray
            The Mino time for the given Boyer-Lindquist time value(s).
        """
        if isinstance(t, np.ndarray) or isinstance(t, list):
            return self.base.mino_time_vec(np.array(t))
        else:
            return self.base.mino_time(t)
    
    def position_of_t(self, t):
        """
        Function that returns the position vector (r, theta, phi) for a given Boyer-Lindquist time value.
        
        Parameters
        ----------
        t : float or numpy.ndarray
            The Boyer-Lindquist time value(s). If a numpy array is provided, the function will return a numpy array of the same shape.
        
        Returns 
        -------
        numpy.ndarray
            The position vector (r, theta, phi) for the given Boyer-Lindquist time value(s).
        """
        if isinstance(t, np.ndarray) or isinstance(t, list):
            return self.base.position_time_vec(np.array(t))
        else:
            return self.base.position_time(t)
        
    def position(self, la):
        """
        Function that returns the position vector (t, r, theta, phi) for a given Mino time value.
        Parameters
        ----------
        la : float or numpy.ndarray
            The Mino time value(s). If a numpy array is provided, the function will return a numpy array of the same shape.
        Returns 
        -------
        numpy.ndarray
            The position vector (t, r, theta, phi) for the given Mino time value(s).
        """
        if isinstance(la, np.ndarray) or isinstance(la, list):
            return self.base.position_vec(np.array(la))
        else:
            return self.base.position(la)
    
    def __call__(self, la):
        return self.position(la)