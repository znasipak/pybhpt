from cybhpt_full import RedshiftCoefficients as _RedshiftCoefficientsCython
from cybhpt_full import SphericalHarmonicCoupling as _SphericalHarmonicCouplingCython

class RedshiftCoefficients:
    def __init__(self, gauge, geo):
        self.base = _RedshiftCoefficientsCython(gauge, geo.base)
    def __call__(self, Ni, ai, bi, ci, di, jr = 0, jz = 0):
        return self.base(Ni, ai, bi, ci, di, jr, jz)
    
class SphericalHarmonicCoupling:
    def __init__(self, lmax, m):
        self.base = _SphericalHarmonicCouplingCython(lmax, m)
        self.azimuthalmode = self.base.azimuthalmode

    def zcouplingcoefficient(self, n, i, l):
        return self.base.zcouplingcoefficient(n, i, l)

    def dzcouplingcoefficient(self, n, i, l):
        return self.base.dzcouplingcoefficient(n, i, l)