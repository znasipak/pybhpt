from cybhpt_full import TeukolskyMode as TeukolskyModeCython

class TeukolskyMode:
    """A class for computing Teukolsky modes sourced by a point-particle orbiting in a Kerr background.
    
    Parameters
    ----------
    s : int
        The spin weight of the Teukolsky mode.
    j : int
        The spheroidal harmonic mode number.
    m : int
        The azimuthal harmonic mode number.
    k : int
        The polar harmonic mode number.
    n : int
        The radial harmonic mode number.
    geo : KerrGeodesic class instance
        KerrGeodesic object containing the background motion of the point-particle source.
    auto_solve : bool, optional
        If True, the Teukolsky equation is automatically solved upon initialization. Default is False

    Attributes
    ----------
    spinweight : int
        The spin weight of the Teukolsky mode.
    spheroidalmode : int
        The spheroidal harmonic mode number.
    azimuthalmode : int
        The azimuthal harmonic mode number.
    radialmode : int
        The radial harmonic mode number.
    polarmode : int
        The polar harmonic mode number.
    blackholespin : float
        The spin of the black hole in the Kerr background.
    frequency : float
        The frequency of the Teukolsky mode.
    horizonfrequency : float
        The frequency of the mode at the horizon.
    eigenvalue : float
        The spheroidal eigenvalue of the Teukolsky mode.
    mincouplingmode : int
        The minimum l-mode used for coupling the spherical and spheroidal harmonics
    maxcouplingmode : int
        The maximum l-mode used for coupling the spherical and spheroidal harmonics
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
    
    Methods
    -------
    solve(geo, method = "AUTO", nsamples = 256, teuk = None, swsh = None)
        Solve the Teukolsky equation for the given mode and geodesic.
    flipspinweight()
        Flip the spin weight of the Teukolsky mode.
    flipspinweightandfrequency()
        Flip the spin weight and frequency of the Teukolsky mode.
    couplingcoefficient(l)
        Compute the coupling coefficient for the given l-mode.
    radialpoint(pos)
        Compute the radial point for the given position.
    radialsolution(bc, pos)
        Compute the radial solution for the given boundary condition and position.
    radialderivative(bc, pos)
        Compute the radial derivative for the given boundary condition and position.
    radialderivative2(bc, pos)
        Compute the second radial derivative for the given boundary condition and position.
    homogeneousradialsolution(bc, pos)
        Compute the homogeneous radial solution for the given boundary condition and position.
    homogeneousradialderivative(bc, pos)
        Compute the homogeneous radial derivative for the given boundary condition and position.
    homogeneousradialderivative2(bc, pos)
        Compute the second homogeneous radial derivative for the given boundary condition and position.
    polarpoint(pos)
        Compute the polar point for the given position.
    polarsolution(pos)
        Compute the polar solution for the given position.
    polarderivative(pos)
        Compute the polar derivative for the given position.
    polarderivative2(pos)
        Compute the second polar derivative for the given position.
    amplitude(bc)
        Compute the Teukolsky amplitude for the given boundary condition.
    precision(bc)
        Compute the precision of the Teukolsky amplitude for the given boundary condition.
    """
    def __init__(self, s, j, m, k, n, geo, auto_solve = False):
        self.base = TeukolskyModeCython(s, j, m, k, n, geo.base)
        if auto_solve:
            self.solve(geo)

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
        return self.base.polarmode

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
    
    def solve(self, geo, method = "AUTO", nsamples = 256, teuk = None, swsh = None):
        """Solve the Teukolsky equation for the given mode and geodesic.
        Parameters
        ----------
        geo : KerrGeodesic class instance
            KerrGeodesic object containing the background motion of the point-particle source.
        method : str, optional
            The method to use for solving the Teukolsky equation. Default is "AUTO".
        nsamples : int, optional
            The number of samples to use for the solution. Default is 256.
        teuk : RadialTeukolsky, optional
            RadialTeukolsky object to use for constructing the radial Green function. Default is None.
        swsh : SpheroidalHarmonicMode, optional
            SpheroidalHarmonic object to use for coupling with spheroidal harmonics. Default is None.

        """
        if teuk is None or swsh is None:
            self.base.solve(geo.base, method, nsamples)
        else:
            self.base.solve(geo.base, method, nsamples, teuk.base, swsh.base)

    """
    Flips the spin-weight of the Teukolsky solutions from :math:`s \rightarrow -s`
    """
    def flipspinweight(self):
        self.base.flip_spinweight()

    """
    Flips the spin-weight and frequency of the Teukolsky solutions from :math:`s \rightarrow -s` and :math:`\omega \rightarrow -\omega`
    """
    def flipspinweightandfrequency(self):
        self.base.flip_spinweight_frequency()

    """
    Spherical-spheroidal mixing coefficient between a spherical harmonic :math:`l` mode with a spheroidal :math:`j` mode

    :param l: spherical harmonic mode
    :type l: int
    """
    def couplingcoefficient(self, l):
        return self.base.couplingcoefficient(l)

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
        return self.base.teukolsky_amplitude(bc)
    
    def precision(self, bc):
        return self.base.teukolsky_amplitude_precision(bc)