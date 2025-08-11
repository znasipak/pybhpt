from libcpp.vector cimport vector
from libcpp.complex cimport complex as cpp_complex
from cython.operator cimport dereference
from libcpp.string cimport string

include "flux_wrap.pyx"

cdef extern from "redshift.hpp":
    cdef cppclass RedshiftCoefficientsCPP "RedshiftCoefficients":
        RedshiftCoefficientsCPP(Gauge gauge, GeodesicSource &geo)
        cpp_complex[double] getComponent(int Ni, int ai, int bi, int ci, int di, int jru, int jzu)

cdef extern from "metriccoeffs.hpp":
    cpp_complex[double] metric_coefficient_ORG(int ai, int bi, int nt, int nr, int nz, int nphi, double a, double r, double z)
    cpp_complex[double] metric_coefficient_IRG(int ai, int bi, int nt, int nr, int nz, int nphi, double a, double r, double z)

cdef extern from "metric.hpp":
    cdef cppclass SphericalHarmonicCouplingCPP "SphericalHarmonicCoupling":
        SphericalHarmonicCouplingCPP(int lmax, int m)
        void generateCouplings()
        int getAzimuthalModeNumber()

        double getZCouplingCoefficient(int n, int i, int l)
        double getZCouplingCoefficientNoCheck(int n, int i, int l)
        int getMinZCouplingModeNumber()
        int getMaxZCouplingModeNumber()

        double getDerivativeCouplingCoefficient(int n, int i, int l)
        double getDerivativeCouplingCoefficientNoCheck(int n, int i, int l)
        int getMinDerivativeCouplingModeNumber()
        int getMaxDerivativeCouplingModeNumber()


def metric_coefficients_cython_S4dagger(int ai, int bi, int nt, int nr, int nz, int nphi, double a, double r, double z):
    return metric_coefficient_ORG(ai, bi, nt, nr, nz, nphi, a, r, z)

def metric_coefficients_cython_S0dagger(int ai, int bi, int nt, int nr, int nz, int nphi, double a, double r, double z):
    return metric_coefficient_IRG(ai, bi, nt, nr, nz, nphi, a, r, z)

cdef class RedshiftCoefficients:
    cdef RedshiftCoefficientsCPP *huucpp

    def __cinit__(self, unicode gauge, KerrGeodesic geo):
        self.huucpp = new RedshiftCoefficientsCPP(str_to_gauge(gauge), dereference(geo.geocpp))

    def __dealloc__(self):
        del self.huucpp

    def __call__(self, int Ni, int ai, int bi, int ci, int di, int jr, int jz):
        return self.huucpp.getComponent(Ni, ai, bi, ci, di, jr, jz)

cdef class SphericalHarmonicCoupling:
    cdef SphericalHarmonicCouplingCPP *cpp

    def __cinit__(self, int lmax, int m):
        self.cpp = new SphericalHarmonicCouplingCPP(lmax, m)
        self.cpp.generateCouplings()

    def __dealloc__(self):
        del self.cpp

    @property
    def azimuthalmode(self):
        return self.cpp.getAzimuthalModeNumber()

    def zcouplingcoefficient(self, int n, int i, int l):
        return self.cpp.getZCouplingCoefficient(n, i, l)

    def dzcouplingcoefficient(self, int n, int i, int l):
        return self.cpp.getDerivativeCouplingCoefficient(n, i, l)

