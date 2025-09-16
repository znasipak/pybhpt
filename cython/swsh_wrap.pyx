from libcpp.vector cimport vector
from libcpp.complex cimport complex as cpp_complex
import numpy as np
cimport numpy as np

cdef extern from "gsl/gsl_errno.h":
    void gsl_set_error_handler_off()

# call at import-time
gsl_set_error_handler_off()

cdef extern from "swsh.hpp":
    cdef cppclass SpinWeightedHarmonic:
        SpinWeightedHarmonic(int s, int L, int m, double gamma, const vector[double]& theta)

        int getSpinWeight()
        int getSpheroidalModeNumber()
        int getAzimuthalModeNumber()
        double getSpheroidicity()
        double getEigenvalue()
        vector[double] getCouplingCoefficient()
        double getCouplingCoefficient(int l)
        int getMinCouplingModeNumber()
        int getMaxCouplingModeNumber()

    	int generateSolutionsAndDerivatives()
        int generateCouplingCoefficients()
        int generateSolutions()
        int generateDerivatives()

        vector[double] getArguments()
        vector[double] getSolution()
        vector[double] getDerivative()
        vector[double] getSecondDerivative()

        double getArguments(int pos)
        double getSolution(int pos)
        double getDerivative(int pos)
        double getSecondDerivative(int pos)

    double Asljm(const int &s, const int &l, const int &j, const int &m)
    double dAsljm(const int &s, const int &l, const int &j, const int &m)

    double clebsch(const int &j1, const int &j2, const int &j, const int &m1, const int &m2, const int &m)
    double w3j(const int &j1, const int &j2, const int &j, const int &m1, const int &m2, const int &m);

    cpp_complex[double] Sslm(const int &s, const int &l, const int &m, const double &g, const double &th, const double &ph)
    double Sslm(const int &s, const int &l, const int &m, const double &g, const double &th)

    double Sslm(const int &s, const int &l, const int &m, const double &g, const vector[double]& bvec, const double &th)
    double Sslm_derivative(const int &s, const int &l, const int &m, const double &g, const vector[double]& bvec, const double &th)
    double Sslm_secondDerivative(const int &s, const int &l, const int &m, const double &g, const double& la, const double &th, const double &Slm, const double &SlmP)

    double swsh_eigenvalue(const int &s, const int &l, const int &m, const double &g)

    cpp_complex[double] Yslm(const int &s, const int &l, const int &m, const double &th, const double &ph)
    double Yslm(const int &s, const int &l, const int &m, const double &th)
    double Yslm_derivative(const int &s, const int &l, const int &m, const double &th)

def Yslm(int s, int l, int m, double theta):
    return Yslm(s, l, m, theta)

def Yslm_derivative(int s, int l, int m, double theta):
    return Yslm_derivative(s, l, m, theta)

def Sslm(int s, int l, int m, double gamma, double theta):
    return Sslm(s, l, m, gamma, theta)
def Sslm_derivative(int s, int l, int m, double gamma, double theta):
    return Sslm_derivative(s, l, m, gamma, theta)
def Sslm_secondDerivative(int s, int l, int m, double gamma, double theta):
    cdef double la = swsh_eigenvalue(s, l, m, gamma)
    cdef double Slm = Sslm(s, l, m, gamma, theta)
    cdef double SlmP = Sslm_derivative(s, l, m, gamma, theta)
    return Sslm_secondDerivative(s, l, m, gamma, la, theta, Slm, SlmP)

cdef class RadialTeukolsky:
    cdef RadialTeukolskyCPP *teukcpp

    def __cinit__(self, double a, int s, int l, int m, double omega, np.ndarray[ndim=1, dtype=np.float64_t] r not None):
        cdef int length_r = r.shape[0]
        cdef vector[double] rvec
        rvec = vector[double](length_r)
        rvec.assign(&r[0], &r[0] + length_r)
        self.teukcpp = new RadialTeukolskyCPP(a, s, l, m, omega, rvec)
        if self.teukcpp == NULL:
            raise MemoryError('Not enough memory.')
    
    def __dealloc__(self):
        del self.teukcpp
    
    @property
    def blackholespin(self):
        return self.teukcpp.getBlackHoleSpin()
    
    @property
    def spinweight(self):
        return self.teukcpp.getSpinWeight()
    @property
    def s(self):
        return self.spinweight

    @property
    def spheroidalmode(self):
        return self.teukcpp.getSpheroidalModeNumber()

    @property
    def j(self):
        return self.spheroidalmode

    @property
    def azimuthalmode(self):
        return self.teukcpp.getAzimuthalModeNumber()
    
    @property
    def m(self):
        return self.azimuthalmode

    @property
    def frequency(self):
        return self.teukcpp.getModeFrequency()
    
    @property
    def mode_frequency(self):
        return self.frequency

    @property
    def omega(self):
        return self.frequency

    @property
    def eigenvalue(self):
        return self.teukcpp.getSpinWeightedSpheroidalEigenvalue()
    
    def solve_bc(self, unicode method):
        self.teukcpp.generateRetardedBoundaryConditions(str_to_method(method))

    def set_bc(self, unicode bc, cpp_complex[double] R, cpp_complex[double] Rp, double r):
        self.teukcpp.setBoundaryConditions(str_to_bc(bc), R, Rp, r)

    def solve(self, unicode method="AUTO", unicode bc="None"):
        if bc == "None":
            self.teukcpp.generateSolutions(str_to_method(method))
        else:
            self.teukcpp.generateSolutions(str_to_bc(bc), str_to_method(method))
    
    def flip_spinweight(self):
        self.teukcpp.flipSpinWeight()

    def radialpoint(self, int pos):
        return self.teukcpp.getRadialPoints(pos)

    def boundarypoint(self, unicode bc):
        return self.teukcpp.getBoundaryPoint(str_to_bc(bc))

    def boundarysolution(self, unicode bc):
        return self.teukcpp.getBoundarySolution(str_to_bc(bc)).getValue()
    
    def boundaryderivative(self, unicode bc):
        return self.teukcpp.getBoundaryDerivative(str_to_bc(bc)).getValue()

    def solution(self, unicode bc, int pos):
        return self.teukcpp.getSolution(str_to_bc(bc), pos)

    def derivative(self, unicode bc, int pos):
        return self.teukcpp.getDerivative(str_to_bc(bc), pos)

    def second_derivative(self, unicode bc, int pos):
        return self.teukcpp.getSecondDerivative(str_to_bc(bc), pos)
    
    def derivative2(self, unicode bc, int pos):
        return self.teukcpp.getSecondDerivative(str_to_bc(bc), pos)