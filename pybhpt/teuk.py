from cybhpt_full import TeukolskyMode as TeukolskyModeCython

class TeukolskyMode:

    def __init__(self, s, j, m, k, n, geo):
        self.base = TeukolskyModeCython(s, j, m, k, n, geo.base)

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
        if teuk is None or swsh is None:
            self.base.solve(geo.base, method, nsamples)
        else:
            self.base.solve(geo.base, method, nsamples, teuk.base, swsh.base)

    def flipspinweight(self):
        self.base.flip_spinweight()

    def flipspinweightandfrequency(self):
        self.base.flip_spinweight_frequency()

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