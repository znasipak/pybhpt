from cybhpt_full import circular_redshift, run_tests, metric_coefficients_cython_ORG, metric_coefficients_cython_IRG
from cybhpt_full import RedshiftCoefficients as RedshiftCoefficientsCython
from cybhpt_full import SphericalHarmonicCoupling as SphericalHarmonicCouplingCython

# def metric_coefficients_ORG(ai, bi, nt, nr, nz, np, a, r, z):
#     return metric_coefficients_cython_ORG(ai, bi, nt, nr, nz, np, a, r, z)

# def metric_coefficients_IRG(ai, bi, nt, nr, nz, np, a, r, z):
#     return metric_coefficients_cython_IRG(ai, bi, nt, nr, nz, np, a, r, z)

class RedshiftCoefficients:
    def __init__(self, gauge, geo):
        self.base = RedshiftCoefficientsCython(gauge, geo.base)

    def __call__(self, Ni, ai, bi, ci, di, jr = 0, jz = 0):
        return self.base(Ni, ai, bi, ci, di, jr, jz)
    
class SphericalHarmonicCoupling:
    def __init__(self, lmax, m):
        self.base = SphericalHarmonicCouplingCython(lmax, m)
        self.azimuthalmode = self.base.azimuthalmode

    def zcouplingcoefficient(self, n, i, l):
        return self.base.zcouplingcoefficient(n, i, l)

    def dzcouplingcoefficient(self, n, i, l):
        return self.base.dzcouplingcoefficient(n, i, l)

def circular_redshift_export(filename, gauge, lmax, geo):
    return circular_redshift(filename, gauge, lmax, geo.base)

def run_unit_tests():
    run_tests()