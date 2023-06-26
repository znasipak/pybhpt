from cybhpt_full import HertzMode as HertzModeCython
from cybhpt_full import test_hertz_mode_cython
from cybhpt_full import teuk_to_hertz_ORG, teuk_to_hertz_IRG

available_gauges = [
    "IRG",
    "ORG",
    "SAAB0", 
    "SAAB4", 
    "ASAAB0", 
    "ASAAB4",
]

def hertz_IRG(Zin, Zup, j, m, k, a, omega, lambdaCH):
    return teuk_to_hertz_IRG(Zin, Zup, j, m, k, a, omega, lambdaCH)

def test_hertz_mode(j, m, k, n, geo):
    test_hertz_mode_cython(j, m, k, n, geo.base)

def gauge_check(gauge):
    if gauge not in available_gauges:
        TypeError("{} is not a supported gauge.".format(gauge))


class HertzMode:

    def __init__(self, teuk, gauge):
        self.base = HertzModeCython(teuk.base, gauge)
        self.gauge = gauge
        

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