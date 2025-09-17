from libcpp.vector cimport vector
from libcpp.complex cimport complex as cpp_complex
from cython.operator cimport dereference
import numpy as np
cimport numpy as np

# from geo_wrap cimport GeodesicSource
include "geo_wrap.pyx"
include "radialsolver_wrap.pyx"
include "swsh_wrap.pyx"

cdef extern from "teukolsky.hpp":
    cdef cppclass SpinWeightedHarmonic:
        pass

    cdef cppclass TeukolskyModeCPP "TeukolskyMode":
        # TeukolskyModeCPP(int L, int m, int k, int n, GeodesicSource& geo)
        TeukolskyModeCPP(int s, int L, int m, int k, int n, GeodesicSource& geo)
        TeukolskyModeCPP(int s, int L, int m, int k, int n, double a, vector[double] theta, vector[double] r)

        int generateSolutions(GeodesicSource& geo, SolutionMethod method, int samplesize)
        int generateSolutions(double omega, GeodesicTrajectory &traj, GeodesicConstants &geoConst, vector[double] r, vector[double] theta, SolutionMethod method, int samplesize)
        int generateSolutions(SpinWeightedHarmonic swsh, RadialTeukolsky teuk, GeodesicTrajectory& traj, GeodesicConstants &geoConst)

        int flipSpinWeightAndFrequency()
        int flipSpinWeight()

        int getSpinWeight()
        int getSpheroidalModeNumber()
        int getAzimuthalModeNumber()
        int getPolarModeNumber()
        int getRadialModeNumber()
        int getSampleSize()
        double getBlackHoleSpin()
        double getFrequency()
        double getHorizonFrequency()
        double getEigenvalue()
        vector[double] getCouplingCoefficient()
        double getCouplingCoefficient(int l)
        int getMinCouplingModeNumber()
        int getMaxCouplingModeNumber()

        vector[double] getRadialPoints()
        cpp_complex[double] getHomogeneousRadialSolution(BoundaryCondition bc)
        cpp_complex[double] getHomogeneousRadialDerivative(BoundaryCondition bc)
        cpp_complex[double] getHomogeneousSecondRadialDerivative(BoundaryCondition bc)
        cpp_complex[double] getTeukolskyAmplitude(BoundaryCondition bc)
        double getTeukolskyAmplitudePrecision(BoundaryCondition bc)
        cpp_complex[double] getRadialSolution(BoundaryCondition bc)
        cpp_complex[double] getRadialDerivative(BoundaryCondition bc)

        vector[double] getPolarPoints()
        vector[double] getPolarSolution()
        vector[double] getPolarDerivative()

        int getRadialSampleNumber()
        double getRadialPoints(int pos)
        cpp_complex[double] getHomogeneousRadialSolution(BoundaryCondition bc, int pos)
        cpp_complex[double] getHomogeneousRadialDerivative(BoundaryCondition bc, int pos)
        cpp_complex[double] getHomogeneousSecondRadialDerivative(BoundaryCondition bc, int pos)
        cpp_complex[double] getRadialSolution(BoundaryCondition bc, int pos)
        cpp_complex[double] getRadialDerivative(BoundaryCondition bc, int pos)

        int getPolarSampleNumber()
        double getPolarPoints(int pos)
        double getPolarSolution(int pos)
        double getPolarDerivative(int pos)
        double getPolarSecondDerivative(int pos)

cdef extern from "hertz.hpp":
    cdef enum Gauge:
        ORG, IRG, SAAB0, SAAB4, ASAAB0, ASAAB4, SAAB, ASAAB

    cdef cppclass HertzModeCPP "HertzMode":
        HertzModeCPP(TeukolskyModeCPP& teuk, Gauge gauge)

        int generateSolutions()

        Gauge getGauge()
        int getSpinWeight()
        int getSpheroidalModeNumber()
        int getAzimuthalModeNumber()
        int getPolarModeNumber()
        int getRadialModeNumber()
        int getSampleSize()
        double getBlackHoleSpin()
        double getFrequency()
        double getHorizonFrequency()
        double getEigenvalue()

        # Vector getCouplingCoefficient()
        double getCouplingCoefficient(int l)
        int getMinCouplingModeNumber()
        int getMaxCouplingModeNumber()

        # Vector getScalarCouplingCoefficient()
        double getScalarCouplingCoefficient(int l)
        int getMinScalarCouplingModeNumber()
        int getMaxScalarCouplingModeNumber()

        # Vector getRadialPoints()
        # ComplexVector getHomogeneousRadialSolution(BoundaryCondition bc)
        # ComplexVector getHomogeneousRadialDerivative(BoundaryCondition bc)
        cpp_complex[double] getHertzAmplitude(BoundaryCondition bc)
        # ComplexVector getRadialSolution(BoundaryCondition bc)
        # ComplexVector getRadialDerivative(BoundaryCondition bc)

        # Vector getPolarPoints()
        # Vector getPolarSolution()
        # Vector getPolarDerivative()

        double getRadialPoints(int pos)
        cpp_complex[double] getHomogeneousRadialSolution(BoundaryCondition bc, int pos)
        cpp_complex[double] getHomogeneousRadialDerivative(BoundaryCondition bc, int pos)
        cpp_complex[double] getHomogeneousRadialSecondDerivative(BoundaryCondition bc, int pos)
        cpp_complex[double] getHomogeneousRadialSolution(BoundaryCondition bc, int dr, int pos)
        cpp_complex[double] getRadialSolution(BoundaryCondition bc, int pos)
        cpp_complex[double] getRadialDerivative(BoundaryCondition bc, int pos)
        cpp_complex[double] getRadialSecondDerivative(BoundaryCondition bc, int pos)
        cpp_complex[double] getRadialSolution(BoundaryCondition bc, int dr, int pos)

        double getPolarPoints(int pos)
        double getPolarSolution(int pos)
        double getPolarDerivative(int pos)
        double getPolarSecondDerivative(int pos)
    
    void test_hertz_mode(int j, int m, int k, int n, GeodesicSource& geo)
    void teukolsky_to_hertz_ORG(cpp_complex[double] &Psi, cpp_complex[double] Zteuk, int L, int m, int k, double a, double omega, double lambdaCH)
    void teukolsky_to_hertz_IRG(cpp_complex[double] &Psi, cpp_complex[double] Zteuk, int L, int m, int k, double a, double omega, double lambdaCH)
    void teukolsky_to_hertz_SAAB(cpp_complex[double] &Psi, cpp_complex[double] Zteuk, int L, int m, int k, double a, double omega, double lambdaCH)
    void teukolsky_to_hertz_ASAAB(cpp_complex[double] &Psi, cpp_complex[double] Zteuk, int L, int m, double a, double omega, double lambdaCH)

    void teukolsky_to_hertz_ORG(cpp_complex[double] &PsiIn, cpp_complex[double] &PsiUp, cpp_complex[double] ZteukIn, cpp_complex[double] ZteukUp, int L, int m, int k, double a, double omega, double lambdaCH)
    void teukolsky_to_hertz_IRG(cpp_complex[double] &PsiIn, cpp_complex[double] &PsiUp, cpp_complex[double] ZteukIn, cpp_complex[double] ZteukUp, int L, int m, int k, double a, double omega, double lambdaCH)
    void teukolsky_to_hertz_SAAB(cpp_complex[double] &PsiIn, cpp_complex[double] &PsiUp, cpp_complex[double] ZteukIn, cpp_complex[double] ZteukUp, int L, int m, int k, double a, double omega, double lambdaCH)
    void teukolsky_to_hertz_ASAAB(cpp_complex[double] &PsiIn, cpp_complex[double] &PsiUp, cpp_complex[double] ZteukIn, cpp_complex[double] ZteukUp, int L, int m, double a, double omega, double lambdaCH)

cdef extern from "metriccoeffs.hpp":
    cpp_complex[double] metric_coefficient_ORG(int ai, int bi, int nt, int nr, int nz, int nphi, double a, double r, double z)
    cpp_complex[double] metric_coefficient_IRG(int ai, int bi, int nt, int nr, int nz, int nphi, double a, double r, double z)
    vector[vector[vector[cpp_complex[double]]]] metric_coefficients_ORG_11(double a, vector[double] r, vector[double] z)

cdef dict gauge_dict = {
    "ORG" : Gauge.ORG, 
    "IRG" : Gauge.IRG, 
    "SRG0" : Gauge.SAAB0, 
    "SRG4" : Gauge.SAAB4, 
    "ARG0" : Gauge.ASAAB0, 
    "ARG4" : Gauge.ASAAB4,
    # "SAAB" : Gauge.SAAB, 
    # "ASAAB" : Gauge.ASAAB,
}

cdef Gauge str_to_gauge(unicode gauge_str) except *:
    if gauge_str in gauge_dict.keys():
        return gauge_dict[gauge_str]
    else:
        print("Error")
        TypeError("{} is not a supported gauge.".format(gauge_str))

cdef class TeukolskyMode:
    cdef TeukolskyModeCPP *teukcpp
    cdef int sampleR
    cdef int sampleTh

    def __cinit__(self, int s, int j, int m, int k, int n, KerrGeodesic geo):
        self.teukcpp = new TeukolskyModeCPP(s, j, m, k, n, dereference(geo.geocpp))
        self.sampleR = 1
        self.sampleTh = 1

    def __dealloc__(self):
        del self.teukcpp

    @property
    def radialsamplenumber(self):
        return self.sampleR
    
    @property
    def polarsamplenumber(self):
        return self.sampleTh

    @property
    def spinweight(self):
        return self.teukcpp.getSpinWeight()

    @property
    def spheroidalmode(self):
        return self.teukcpp.getSpheroidalModeNumber()
    
    @property
    def azimuthalmode(self):
        return self.teukcpp.getAzimuthalModeNumber()

    @property
    def radialmode(self):
        return self.teukcpp.getRadialModeNumber()

    @property
    def polarmode(self):
        return self.teukcpp.getPolarModeNumber()

    @property
    def blackholespin(self):
        return self.teukcpp.getBlackHoleSpin()
    
    @property
    def frequency(self):
        return self.teukcpp.getFrequency()

    @property
    def horizonfrequency(self):
        return self.teukcpp.getHorizonFrequency()

    @property
    def eigenvalue(self):
        return self.teukcpp.getEigenvalue()

    @property
    def mincouplingmode(self):
        return self.teukcpp.getMinCouplingModeNumber()
    
    @property
    def maxcouplingmode(self):
        return self.teukcpp.getMaxCouplingModeNumber()

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
    def radialpoints(self):
        return np.array([self.teukcpp.getRadialPoints(i) for i in range(self.sampleR)])
    @property
    def polarpoints(self):
        return np.array([self.teukcpp.getPolarPoints(i) for i in range(self.sampleTh)])
    @property
    def radialsolutions(self):
        return {
                "In":np.array([self.teukcpp.getHomogeneousRadialSolution(BoundaryCondition.In, i) for i in range(self.sampleR)]),
                "Up":np.array([self.teukcpp.getHomogeneousRadialSolution(BoundaryCondition.Up, i) for i in range(self.sampleR)])
        }

    def teukolsky_amplitude(self, unicode bc):
        return self.teukcpp.getTeukolskyAmplitude(str_to_bc(bc))

    def teukolsky_amplitude_precision(self, unicode bc):
        return self.teukcpp.getTeukolskyAmplitudePrecision(str_to_bc(bc))

    def couplingcoefficient(self, int l):
        return self.teukcpp.getCouplingCoefficient(l)
    def radialpoint(self, int i):
        return self.teukcpp.getRadialPoints(i)
    def radialsolution(self, unicode bc, int i):
        return self.teukcpp.getRadialSolution(str_to_bc(bc), i)
    def radialderivative(self, unicode bc, int i):
        return self.teukcpp.getRadialDerivative(str_to_bc(bc), i)
    def radialderivative2(self, unicode bc, int i):
        return self.teukolsky_amplitude(bc)*self.homogeneousradialderivative2(bc, i)
    def homogeneousradialsolution(self, unicode bc, int i):
        return self.teukcpp.getHomogeneousRadialSolution(str_to_bc(bc), i)
    def homogeneousradialderivative(self, unicode bc, int i):
        return self.teukcpp.getHomogeneousRadialDerivative(str_to_bc(bc), i)
    def homogeneousradialderivative2(self, unicode bc, int i):
        return self.teukcpp.getHomogeneousSecondRadialDerivative(str_to_bc(bc), i)
    def polarpoint(self, int i):
        return self.teukcpp.getPolarPoints(i)
    def polarsolution(self, int i):
        return self.teukcpp.getPolarSolution(i)
    def polarderivative(self, int i):
        return self.teukcpp.getPolarDerivative(i)
    def polarderivative2(self, int i):
        return self.teukcpp.getPolarSecondDerivative(i)

    def solve(self, KerrGeodesic geo, unicode method = "AUTO", int nsample = 256, teuk=None, swsh=None):
        self.teukcpp.generateSolutions(dereference(geo.geocpp), str_to_method(method), nsample)
        self.sampleR = self.teukcpp.getRadialSampleNumber()
        self.sampleTh = self.teukcpp.getPolarSampleNumber()

    def flip_spinweight_frequency(self):
        self.teukcpp.flipSpinWeightAndFrequency()
    
    def flip_spinweight(self):
        self.teukcpp.flipSpinWeight()

cdef class HertzMode:
    cdef HertzModeCPP *hertzcpp
    cdef unicode gauge_str
    cdef Gauge gauge_cpp
    cdef int sampleR
    cdef int sampleTh

    def __init__(self, TeukolskyMode teuk, unicode gauge):
        if np.abs(teuk.spinweight) != 2:
            raise ValueError("Hertz mode only accepts Teukolsky solutions with spin-weight -2,+2.")
        self.gauge_cpp = str_to_gauge(gauge)
        self.gauge_str = gauge
        self.hertzcpp = new HertzModeCPP(dereference(teuk.teukcpp), self.gauge_cpp)
        self.sampleR = teuk.sampleR
        self.sampleTh = teuk.sampleTh

    def solve(self):
        self.hertzcpp.generateSolutions()
        
    def __dealloc__(self):
        del self.hertzcpp

    @property
    def radialsamplenumber(self):
        return self.sampleR
    
    @property
    def polarsamplenumber(self):
        return self.sampleTh

    @property
    def gauge(self):
        return self.gauge_str

    @property
    def spinweight(self):
        return self.hertzcpp.getSpinWeight()

    @property
    def spheroidalmode(self):
        return self.hertzcpp.getSpheroidalModeNumber()
    
    @property
    def azimuthalmode(self):
        return self.hertzcpp.getAzimuthalModeNumber()

    @property
    def radialmode(self):
        return self.hertzcpp.getRadialModeNumber()

    @property
    def polarmode(self):
        return self.hertzcpp.getPolarModeNumber()

    @property
    def blackholespin(self):
        return self.hertzcpp.getBlackHoleSpin()
    
    @property
    def frequency(self):
        return self.hertzcpp.getFrequency()

    @property
    def horizonfrequency(self):
        return self.hertzcpp.getHorizonFrequency()

    @property
    def eigenvalue(self):
        return self.hertzcpp.getEigenvalue()

    @property
    def mincouplingmode(self):
        return self.hertzcpp.getMinCouplingModeNumber()
    
    @property
    def maxcouplingmode(self):
        return self.hertzcpp.getMaxCouplingModeNumber()

    @property
    def minscalarcouplingmode(self):
        return self.hertzcpp.getMinScalarCouplingModeNumber()
    
    @property
    def maxscalarcouplingmode(self):
        return self.hertzcpp.getMaxScalarCouplingModeNumber()

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

    def hertz_amplitude(self, unicode bc):
        return self.hertzcpp.getHertzAmplitude(str_to_bc(bc))

    def couplingcoefficient(self, int l):
        return self.hertzcpp.getCouplingCoefficient(l)
    def scalarcouplingcoefficient(self, int l):
        return self.hertzcpp.getScalarCouplingCoefficient(l)
    def radialpoint(self, int i):
        return self.hertzcpp.getRadialPoints(i)
    def radialsolution(self, unicode bc, int i):
        return self.hertzcpp.getRadialSolution(str_to_bc(bc), i)
    def radialderivative(self, unicode bc, int i):
        return self.hertzcpp.getRadialDerivative(str_to_bc(bc), i)
    def radialderivative2(self, unicode bc, int i):
        return self.hertz_amplitude(bc)*self.homogeneousradialderivative2(bc, i)
    def homogeneousradialsolution(self, unicode bc, int i):
        return self.hertzcpp.getHomogeneousRadialSolution(str_to_bc(bc), i)
    def homogeneousradialderivative(self, unicode bc, int i):
        return self.hertzcpp.getHomogeneousRadialDerivative(str_to_bc(bc), i)
    def homogeneousradialderivative2(self, unicode bc, int i):
        return self.hertzcpp.getHomogeneousRadialSecondDerivative(str_to_bc(bc), i)
    def polarpoint(self, int i):
        return self.hertzcpp.getPolarPoints(i)
    def polarsolution(self, int i):
        return self.hertzcpp.getPolarSolution(i)
    def polarderivative(self, int i):
        return self.hertzcpp.getPolarDerivative(i)
    def polarderivative2(self, int i):
        return self.hertzcpp.getPolarSecondDerivative(i)

    @property
    def couplingcoefficients(self):
        return np.array([self.couplingcoefficient(i) for i in range(self.mincouplingmode, self.maxcouplingmode+1)])

    @property
    def scalarcouplingcoefficients(self):
        return np.array([self.scalarcouplingcoefficient(i) for i in range(self.minscalarcouplingmode, self.maxscalarcouplingmode+1)])

    @property
    def polarpoints(self):
        return np.array([self.polarpoint(i) for i in range(self.sampleTh)])
    
    @property
    def polarsolutions(self):
        return np.array([self.polarsolution(i) for i in range(self.sampleTh)])
        
    @property
    def polarderivatives(self):
        return np.array([self.polarderivative(i) for i in range(self.sampleTh)])

    @property
    def polarderivatives2(self):
        return np.array([self.polarderivative2(i) for i in range(self.sampleTh)])
    
    @property
    def radialpoints(self):
        return np.array([self.radialpoint(i) for i in range(self.sampleR)])
    @property
    def radialsolutions(self):
        return {
                "In":np.array([self.hertzcpp.getRadialSolution(BoundaryCondition.In, i) for i in range(self.sampleR)]),
                "Up":np.array([self.hertzcpp.getRadialSolution(BoundaryCondition.Up, i) for i in range(self.sampleR)])
        }
    @property
    def radialderivatives(self):
        return {
                "In":np.array([self.hertzcpp.getRadialDerivative(BoundaryCondition.In, i) for i in range(self.sampleR)]),
                "Up":np.array([self.hertzcpp.getRadialDerivative(BoundaryCondition.Up, i) for i in range(self.sampleR)])
        }
    @property
    def radialderivatives2(self):
        return {
                "In":np.array([self.radialderivative2('In', i) for i in range(self.sampleR)]),
                "Up":np.array([self.radialderivative2('Up', i) for i in range(self.sampleR)])
        }
        
def test_hertz_mode_cython(int j, int m, int k, int n, KerrGeodesic geo):
    test_hertz_mode(j, m, k, n, dereference(geo.geocpp))

cdef dict basis_dict = {
    "tetrad": None,
    "coordinate": None
}

def teuk_to_hertz_ORG(cpp_complex[double] ZIn, cpp_complex[double] ZUp, int j, int m, int k, double a, double omega, double lambdaCH):
    cdef cpp_complex[double] PsiIn
    cdef cpp_complex[double] PsiUp
    teukolsky_to_hertz_ORG(PsiIn, PsiUp, ZIn, ZUp, j, m, k, a, omega, lambdaCH)
    return (PsiIn, PsiUp)

def teuk_to_hertz_IRG(cpp_complex[double] ZIn, cpp_complex[double] ZUp, int j, int m, int k, double a, double omega, double lambdaCH):
    cdef cpp_complex[double] PsiIn
    cdef cpp_complex[double] PsiUp
    teukolsky_to_hertz_IRG(PsiIn, PsiUp, ZIn, ZUp, j, m, k, a, omega, lambdaCH)
    return (PsiIn, PsiUp)

def teuk_to_hertz_SRG(cpp_complex[double] ZIn, cpp_complex[double] ZUp, int j, int m, int k, double a, double omega, double lambdaCH):
    cdef cpp_complex[double] PsiIn
    cdef cpp_complex[double] PsiUp
    teukolsky_to_hertz_SAAB(PsiIn, PsiUp, ZIn, ZUp, j, m, k, a, omega, lambdaCH)
    return (PsiIn, PsiUp)

def teuk_to_hertz_ARG(cpp_complex[double] ZIn, cpp_complex[double] ZUp, int j, int m, int k, double a, double omega, double lambdaCH):
    cdef cpp_complex[double] PsiIn
    cdef cpp_complex[double] PsiUp
    teukolsky_to_hertz_ASAAB(PsiIn, PsiUp, ZIn, ZUp, j, m, a, omega, lambdaCH)
    return (PsiIn, PsiUp)

cdef dict metric_component_gauge_dict = {
    "ORG" : {
        (1, 1): None, 
        (1, 3): None,
        (1, 4): None,
        (3, 3): None,
        (4, 4): None}, 
    "IRG" : {
        (2, 2): None,
        (2, 3): None, 
        (2, 4): None,
        (3, 3): None,
        (4, 4): None}, 
    "SRG0" : {
        (2, 2): None, 
        (2, 3): None, 
        (2, 4): None,
        (3, 3): None,
        (4, 4): None}, 
    "SRG4" : {
        (1, 1): None, 
        (1, 3): None, 
        (1, 4): None,
        (3, 3): None,
        (4, 4): None}, 
    "ARG0" : {
        (2, 2): None, 
        (2, 3): None,
        (2, 4): None,
        (3, 3): None,
        (4, 4): None}, 
    "ARG4" : {
        (1, 1): None, 
        (1, 3): None, 
        (1, 4): None,
        (3, 3): None,
        (4, 4): None},
    # "SAAB" : {
    #     (1, 1): None,
    #     (1, 3): None,
    #     (1, 4): None,
    #     (2, 2): None, 
    #     (2, 3): None, 
    #     (2, 4): None,
    #     (3, 3): None,
    #     (4, 4): None}, 
    # "ASAAB" : {
    #     (1, 1): None,
    #     (1, 3): None,
    #     (1, 4): None,
    #     (2, 2): None, 
    #     (2, 3): None, 
    #     (2, 4): None,
    #     (3, 3): None,
    #     (4, 4): None},
}

def metric_11(double a, double r, double z):
    cdef vector[double] rvec = vector[double](1)
    cdef vector[double] zvec = vector[double](1)
    rvec[0] = r
    zvec[0] = z
    cdef cpp_complex[double] temp

    cdef vector[vector[vector[cpp_complex[double]]]] coeffs = metric_coefficients_ORG_11(a, rvec, zvec)
    
    return np.array(coeffs).squeeze()

def metric_coefficient_S4(int alpha, int beta, int nt, int nr, int nz, int np, double a, double r, double z):
    return metric_coefficient_ORG(alpha, beta, nt, nr, nz, np, a, r, z)

def metric_coefficient_S0(int alpha, int beta, int nt, int nr, int nz, int np, double a, double r, double z):
    return metric_coefficient_IRG(alpha, beta, nt, nr, nz, np, a, r, z)

cdef class MetricModeGenerator:
    cdef unicode gauge_str
    cdef Gauge gauge_cpp
    cdef unicode basis

    def __init__(self, unicode gauge, unicode basis="tetrad"):
        self.gauge_str = gauge
        self.gauge_cpp = str_to_gauge(gauge)
        if basis in basis_dict.keys():
            self.basis = basis
        else:
            raise TypeError('{} is not a valid basis. Must be tetrad or coordinate.'.format(basis))
    
    def __call__(self, HertzMode hertz, int ai, int bi):
        if self.basis == "tetrad":
            return self.tetradcomponent(hertz, ai, bi)
        else:
            return self.tetradcomponent(hertz, ai, bi)
    
    def tetradcomponent(self, HertzMode hertz, int ai, int bi):
        if hertz.gauge is not self.gauge_str:
            raise TypeError("Hertz potential in {} gauge. Must be in {} gauge".format(hertz.gauge, self.gauge_str))
        cdef cpp_complex[double] habIn = cpp_complex[double](0., 0.)
        cdef cpp_complex[double] habUp = cpp_complex[double](0., 0.)
        cdef cpp_complex[double] habbase, dPsiIn, dPsiUp, dS
        cdef int atemp
        cdef double a, r, z
        cdef cpp_complex[double] im = cpp_complex[double](0., hertz.azimuthalmode)
        cdef cpp_complex[double] iomega = cpp_complex[double](0., hertz.frequency)
        a = hertz.blackholespin
        r = hertz.radialpoint(0)
        z = np.cos(hertz.polarpoint(0))
        if(ai < bi):
            atemp = ai
            ai = bi
            bi = atemp
        if (ai, bi) in metric_component_gauge_dict[self.gauge_str].keys():
            for nt in range(3):
                for nr in range(3):
                    if nr == 0:
                        dPsiIn = hertz.radialsolution("In", 0)
                        dPsiUp = hertz.radialsolution("Up", 0)
                    elif nr == 1:
                        dPsiIn = hertz.radialderivative("In", 0)
                        dPsiUp = hertz.radialderivative("Up", 0)
                    elif nr == 2:
                        dPsiIn = hertz.radialderivative2("In", 0)
                        dPsiUp = hertz.radialderivative2("Up", 0)
                    else:
                        dPsiIn = cpp_complex[double](0., 0.)
                        dPsiUp = cpp_complex[double](0., 0.)
                    for nz in range(3):
                        if nz == 0:
                            dS = hertz.polarsolution(0)
                        elif nz == 1: # SWSH are functions of theta but the metric coefficients are organized wrt to d/dz
                            dS = -hertz.polarderivative(0)/np.sqrt(1. - z**2)
                        elif nz == 2:
                            dS = hertz.polarderivative2(0)/(1. - z**2) - z*hertz.polarderivative(0)/(1. - z**2)**(1.5)
                        else:
                            dS = cpp_complex[double](0., 0.)
                        for nph in range(3):
                            if nt + nr + nz + nph <= 2:
                                habbase = metric_coefficient_S0(ai, bi, nt, nr, nz, nph, a, r, z)*dS*(im)**nph*(-1.*iomega)**nt
                                habIn += habbase*dPsiIn
                                habUp += habbase*dPsiUp
            return {"In": habIn, "Up": habUp}
        else:
            return {"In": habIn, "Up": habUp}
