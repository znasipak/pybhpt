from cybhpt_full import TeukolskyMode as TeukolskyModeCython

class TeukolskyMode:
    """
    Class for calculating the :math:`(j,m,k,n)`-mode solution to the inhomogeneous Teukolsky equation for a point-particle
    on a bound geodesic with frequency :math:`\Omega_{mkn} = \Omega_\phi + k \Omega_\theta + n \Omega_r`

    :param s: spin-weight of the perturbation
    :type s: int
    :param j: spheroidal mode of the perturbation
    :type j: int
    :param m: azimuthal mode of the perturbation
    :type m:
    :param k: polar mode of the perturbation
    :type k:
    :param n: radial mode of the perturbation
    :type n:
    :param geo: geodesic solution for the motion of the perturbing point-particle
    :type geo: KerrGeodesic
    :param auto_solve: option to solve the Teukolsky equation upon instantiation
    :type auto_solve: bool, optional. Default is False
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
    
    """
    Solve the Teukolsky equation given the geodesic of the perturbing body

    :param geo: geodesic solution for the motion of the perturbing point-particle
    :type geo: KerrGeodesic 
    :param method: method for solving the radial Teukolsky equation. See RadialTeukolsky for available methods.
    :type method: str, optional. Default is "AUTO"
    :param nsamples:
    :type nsamples:
    :param teuk: Pre-computed radial Teukolsky solution. If None, then the solver will compute the solutions.
    :type teuk: RadialTeukolsky, optional. Default is None.
    :param swsh: Pre-computed spin-weighted spheroidal harmonic.  If None, then the solver will compute the solutions.
    :type swsh: SpinWeightedSpheroidalHarmonics, optional. Default is None.
    """
    def solve(self, geo, method = "AUTO", nsamples = 256, teuk = None, swsh = None):
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