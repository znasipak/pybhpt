from cybhpt_full import HertzMode as HertzModeCython
from cybhpt_full import test_hertz_mode_cython
from cybhpt_full import teuk_to_hertz_ORG, teuk_to_hertz_IRG, teuk_to_hertz_SAAB, teuk_to_hertz_ASAAB
from pybhpt.radial import RadialTeukolsky
import numpy as np

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

def hertz_ORG(Zin, Zup, j, m, k, a, omega, lambdaCH):
    return teuk_to_hertz_ORG(Zin, Zup, j, m, k, a, omega, lambdaCH)

def hertz_SAAB(Zin, Zup, j, m, k, a, omega, lambdaCH):
    return teuk_to_hertz_SAAB(Zin, Zup, j, m, k, a, omega, lambdaCH)

def hertz_ASAAB(Zin, Zup, j, m, k, a, omega, lambdaCH):
    return teuk_to_hertz_ASAAB(Zin, Zup, j, m, k, a, omega, lambdaCH)

def teuk_to_hertz_amplitude(gauge, Zin, Zup, j, m, k, a, omega, lambdaCH):
    if gauge == "IRG":
        return hertz_IRG(Zin, Zup, j, m, k, a, omega, lambdaCH)
    elif gauge == "ORG":
        return hertz_ORG(Zin, Zup, j, m, k, a, omega, lambdaCH)
    elif gauge == "SAAB0" or gauge == "SAAB4":
        return hertz_SAAB(Zin, Zup, j, m, k, a, omega, lambdaCH)
    elif gauge == "ASAAB0" or gauge == "ASAAB4":
        return hertz_ASAAB(Zin, Zup, j, m, k, a, omega, lambdaCH)
    else:
        return (0.j, 0.j)

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
    
    def __call__(self, r, deriv = 0):
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