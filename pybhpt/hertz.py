from cybhpt_full import HertzMode as HertzModeCython
from cybhpt_full import test_hertz_mode_cython
from cybhpt_full import teuk_to_hertz_ORG, teuk_to_hertz_IRG, teuk_to_hertz_SRG, teuk_to_hertz_ARG
from pybhpt.radial import RadialTeukolsky
import numpy as np

available_gauges = [
    "IRG",
    "ORG",
    "SRG0", 
    "SRG4", 
    "ARG0", 
    "ARG4",
]

def hertz_IRG(Zin, Zup, j, m, k, a, omega, lambdaCH):
    """
    Convert Teukolsky amplitudes to Hertz potential in the IRG gauge.
    """
    return teuk_to_hertz_IRG(Zin, Zup, j, m, k, a, omega, lambdaCH)

def hertz_ORG(Zin, Zup, j, m, k, a, omega, lambdaCH):
    """Convert Teukolsky amplitudes to Hertz potential in the ORG gauge.
    """
    return teuk_to_hertz_ORG(Zin, Zup, j, m, k, a, omega, lambdaCH)

def hertz_SRG(Zin, Zup, j, m, k, a, omega, lambdaCH):
    """Convert Teukolsky amplitudes to Hertz potential in the SRG gauge.
    """
    return teuk_to_hertz_SRG(Zin, Zup, j, m, k, a, omega, lambdaCH)

def hertz_ARG(Zin, Zup, j, m, k, a, omega, lambdaCH):
    """Convert Teukolsky amplitudes to Hertz potential in the ARG gauge.
    """
    return teuk_to_hertz_ARG(Zin, Zup, j, m, k, a, omega, lambdaCH)

def teuk_to_hertz_amplitude(gauge, Zin, Zup, j, m, k, a, omega, lambdaCH):
    if gauge == "IRG":
        return hertz_IRG(Zin, Zup, j, m, k, a, omega, lambdaCH)
    elif gauge == "ORG":
        return hertz_ORG(Zin, Zup, j, m, k, a, omega, lambdaCH)
    elif gauge == "SRG0" or gauge == "SRG4":
        return hertz_SRG(Zin, Zup, j, m, k, a, omega, lambdaCH)
    elif gauge == "ARG0" or gauge == "ARG4":
        return hertz_ARG(Zin, Zup, j, m, k, a, omega, lambdaCH)
    else:
        return (0.j, 0.j)

def test_hertz_mode(j, m, k, n, geo):
    test_hertz_mode_cython(j, m, k, n, geo.base)

def gauge_check(gauge):
    """ Check if the provided gauge is supported. 
    Raises a TypeError if the gauge is not supported.
    """
    if gauge not in available_gauges:
        TypeError("{} is not a supported gauge.".format(gauge))


class HertzMode:
    """
    Class that produces a Hertz potential mode given a Teukolsky object and a gauge.
    This class is a wrapper around the Cython implementation of the Hertz potential
    and provides a Python interface to the underlying C++ code.

    Parameters
    ----------
    teuk : Teukolsky
        The Teukolsky object to be used for the Hertz potential.
    gauge : str
        The gauge to be used for the Hertz potential. Must be one of the following:
        - "IRG"
        - "ORG"
        - "SRG0"
        - "SRG4"
        - "ARG0"
        - "ARG4"

    Attributes
    ----------
    base : HertzModeCython
        The underlying Cython implementation of the Hertz potential mode.
    gauge : str
        The gauge used for the Hertz potential.
    sampleR : int
        The number of radial samples used in the Hertz potential mode solutions
    sampleTh : int
        The number of polar samples used in the Hertz potential mode solutions
    spinweight : int
        The spin weight of the Hertz potential mode.
    spheroidalmode : int
        The spheroidal mode number of the Hertz potential mode.
    azimuthalmode : int
        The azimuthal mode number of the Hertz potential mode.
    radialmode : int
        The radial mode number of the Hertz potential mode.
    polarmode : int
        The polar mode number of the Hertz potential mode.
    blackholespin : float
        The spin of the black hole associated with the background spacetime.
    frequency : float
        The frequency of the Hertz potential mode.
    horizonfrequency : float
        The frequency of the Hertz potential mode at the horizon.
    eigenvalue : complex
        The spheroidal eigenvalue associated with the Hertz potential mode.
    mincouplingmode : int
        The minimum l-mode used for coupling the spherical and spheroidal harmonics
    maxcouplingmode : int
        The maximum l-mode used for coupling the spherical and spheroidal harmonics
    minscalarcouplingmode : int
        The minimum l-mode used for coupling the scalar harmonics
    maxscalarcouplingmode : int
        The maximum l-mode used for coupling the scalar harmonics
    j : int
        Alias for spheroidalmode.
    m : int
        Alias for azimuthalmode.
    k : int
        Alias for polarmode.
    n : int
        Alias for radialmode.
    omega : float
        Alias for frequency.
    a : float
        Alias for blackholespin.
    
    Properties
    ----------
    couplingcoefficients : np.ndarray
        The coupling coefficients for the Hertz potential mode.
    scalarcouplingcoefficients : np.ndarray
        The scalar coupling coefficients for the Hertz potential mode.
    polarpoints : np.ndarray
        The polar points used in the Hertz potential mode solutions.
    polarsolutions : np.ndarray
        The polar mode solutions of the Hertz potential.
    polarderivatives : np.ndarray
        Derivatives of the polar mode solutions of the Hertz potential.
    polarderivatives2 : np.ndarray
        Second derivatives of the polar mode solutions of the Hertz potential.
    radialpoints : np.ndarray
        The radial points used in the Hertz potential mode solutions.
    radialsolutions : np.ndarray
        The radial mode solutions of the Hertz potential.
    radialderivatives : np.ndarray
        Derivatives of the radial mode solutions of the Hertz potential.
    radialderivatives2 : np.ndarray
        Second derivatives of the radial mode solutions of the Hertz potential.

    Methods
    -------
    solve()
        Solve the Hertz potential mode equations.
    couplingcoefficient(l)
        Returns the coupling coefficient for the given l-mode.
    scalarcouplingcoefficient(l)
        Returns the scalar coupling coefficient for the given l-mode.
    radialpoint(pos)
        Returns the radial point corresponding to the given position.
    radialsolution(bc, pos)
        Returns the radial solution for the given boundary condition and position.
    radialderivative(bc, pos)
        Returns the radial derivative for the given boundary condition and position.
    radialderivative2(bc, pos)
        Returns the second radial derivative for the given boundary condition and position.
    homogeneousradialsolution(bc, pos)
        Returns the homogeneous radial solution for the given boundary condition and position.
    homogeneousradialderivative(bc, pos)
        Returns the homogeneous radial derivative for the given boundary condition and position.
    homogeneousradialderivative2(bc, pos)
        Returns the second homogeneous radial derivative for the given boundary condition and position.
    polarpoint(pos)
        Returns the polar point corresponding to the given position.
    polarsolution(pos)
        Returns the polar solution for the given position.
    polarderivative(pos)
        Returns the polar derivative for the given position.
    polarderivative2(pos)
        Returns the second polar derivative for the given position.
    amplitude(bc)
        Returns the Hertz amplitude for the given boundary condition.
    __call__(r, deriv=0)
        Returns the radial Hertz potential mode evaluated at the given radial values `r`.
        Alternatively, it can return the radial derivative if `deriv` is set to 1 or 2.
        The radial values `r` must lie outside the source region defined by the radial points.
        If r contains values inside the source region, a ValueError is raised.
    """
    def __init__(self, teuk, gauge):
        self.base = HertzModeCython(teuk.base, gauge)
        self.gauge = gauge
        self.sampleR = self.base.radialsamplenumber
        self.sampleTh = self.base.polarsamplenumber

    @property
    def spinweight(self):
        return self.base.spinweight

    @property
    def spheroidalmode(self):
        return self.base.spheroidalmode
    
    @property
    def azimuthalmode(self):
        return self.base.azimuthalmode

    @property
    def radialmode(self):
        return self.base.radialmode

    @property
    def polarmode(self):
        return self.base.spipolarmodenweight

    @property
    def blackholespin(self):
        return self.base.blackholespin
    
    @property
    def frequency(self):
        return self.base.frequency

    @property
    def horizonfrequency(self):
        return self.base.horizonfrequency

    @property
    def eigenvalue(self):
        return self.base.eigenvalue

    @property
    def mincouplingmode(self):
        return self.base.mincouplingmode
    
    @property
    def maxcouplingmode(self):
        return self.base.maxcouplingmode
    
    @property
    def minscalarcouplingmode(self):
        return self.base.minscalarcouplingmode
    
    @property
    def maxscalarcouplingmode(self):
        return self.base.maxscalarcouplingmode

    # some useful aliases
    @property
    def j(self):
        return self.spheroidalmode
    @property
    def m(self):
        return self.azimuthalmode
    @property
    def k(self):
        return self.polarmode
    @property
    def n(self):
        return self.radialmode
    @property
    def omega(self):
        return self.frequency
    @property
    def a(self):
        return self.blackholespin
    
    @property
    def couplingcoefficients(self):
        return self.base.couplingcoefficients
    
    @property
    def scalarcouplingcoefficients(self):
        return self.base.scalarcouplingcoefficients
    
    @property
    def polarpoints(self):
        return self.base.polarpoints
    
    @property
    def polarsolutions(self):
        return self.base.polarsolutions
    
    @property
    def polarderivatives(self):
        return self.base.polarderivatives
    
    @property
    def polarderivatives2(self):
        return self.base.polarderivatives2
    
    @property
    def radialpoints(self):
        return self.base.radialpoints
    
    @property
    def radialsolutions(self):
        return self.base.radialsolutions
    
    @property
    def radialderivatives(self):
        return self.base.radialderivatives
    
    @property
    def radialderivatives2(self):
        return self.base.radialderivatives2
    
    @property
    def amplitudes(self):
        return {"In": self.amplitude('In'), "Up": self.amplitude('Up')}
    
    def solve(self):
        self.base.solve()

    def couplingcoefficient(self, l):
        return self.base.couplingcoefficient(l)
    
    def scalarcouplingcoefficient(self, l):
        return self.base.scalarcouplingcoefficient(l)

    def radialpoint(self, pos):
        return self.base.radialpoint(pos)
    
    def radialsolution(self, bc, pos):
        return self.base.radialsolution(bc, pos)
    
    def radialderivative(self, bc, pos):
        return self.base.radialderivative(bc, pos)
    
    def radialderivative2(self, bc, pos):
        return self.base.radialderivative2(bc, pos)
    
    def homogeneousradialsolution(self, bc, pos):
        return self.base.homogeneousradialsolution(bc, pos)
    
    def homogeneousradialderivative(self, bc, pos):
        return self.base.homogeneousradialderivative(bc, pos)
    
    def homogeneousradialderivative2(self, bc, pos):
        return self.base.homogeneousradialderivative2(bc, pos)
    
    def polarpoint(self, pos):
        return self.base.polarpoint(pos)
    
    def polarsolution(self, pos):
        return self.base.polarsolution(pos)
    
    def polarderivative(self, pos):
        return self.base.polarderivative(pos)
    
    def polarderivative2(self, pos):
        return self.base.polarderivative2(pos)
    
    def amplitude(self, bc):
        return self.base.hertz_amplitude(bc)
    
    def __call__(self, r, deriv = 0):
        """
        Evaluate the radial Hertz potential mode at the given radial values `r`.
        If `deriv` is set to 0, it returns the Hertz potential mode.
        If `deriv` is set to 1, it returns the first radial derivative.
        If `deriv` is set to 2, it returns the second radial derivative.
        
        Parameters
        ----------
        r : array-like
            The radial values at which to evaluate the Hertz potential mode.
        deriv : int, optional
            The order of the radial derivative to compute. Default is 0 (no derivative). 
        Returns
        -------
        numpy.ndarray
            The evaluated Hertz potential mode or its radial derivative at the given radial values `r`.
        Raises
        ------
        ValueError
            If any of the radial values `r` lie within the source region defined by the radialpoints.
        """
        rmin = self.radialpoints[0]
        rmax = self.radialpoints[-1]
        rinner = r[r < rmin]
        router = r[r > rmax]
        if np.any((r >= rmin) & (r <= rmax)):
            raise ValueError("Radial points must lie outside the source region")

        if rinner.shape[0] > 0:
            Rt = RadialTeukolsky(self.spinweight, self.spheroidalmode, self.azimuthalmode, self.blackholespin, self.frequency, rinner)
            Rt.solve(bc='In')
            Rin = Rt('In', deriv=deriv)
        else:
            Rin = np.array([])

        if router.shape[0] > 0:
            Rt = RadialTeukolsky(self.spinweight, self.spheroidalmode, self.azimuthalmode, self.blackholespin, self.frequency, router)
            Rt.solve(bc='Up')
            Rup = Rt('Up', deriv=deriv)
        else:
            Rup = np.array([])

        PsiIn = self.amplitude('In')
        PsiUp = self.amplitude('Up')

        return np.concatenate((PsiIn*Rin, PsiUp*Rup))