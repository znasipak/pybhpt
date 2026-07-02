from libcpp.vector cimport vector
from libcpp.complex cimport complex as cpp_complex
import numpy as np
cimport numpy as np

cdef extern from "gsl/gsl_errno.h":
  void gsl_set_error_handler_off()

# If you need to disable GSL error handling, do so in a targeted way within specific functions.
# gsl_set_error_handler_off()  # Removed global call to avoid masking numerical errors.

cdef extern from "swsh.hpp":
    cdef cppclass SpinWeightedHarmonic:
        SpinWeightedHarmonic(int s, int L, int m, double gamma, vector[double]& theta)

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

    double Asljm(int &s, int &l, int &j, int &m)
    double dAsljm(int &s, int &l, int &j, int &m)

    double clebsch(int &j1, int &j2, int &j, int &m1, int &m2, int &m)
    double w3j(int &j1, int &j2, int &j, int &m1, int &m2, int &m)

    cpp_complex[double] Sslm(int &s, int &l, int &m, double &g, double &th, double &ph)
    double Sslm(int &s, int &l, int &m, double &g, double &th)

    double Sslm(int &s, int &l, int &m, double &g, vector[double]& bvec, double &th)
    double Sslm_derivative(int &s, int &l, int &m, double &g, vector[double]& bvec, double &th)
    double Sslm_secondDerivative(int &s, int &l, int &m, double &g, double& la, double &th, double &Slm, double &SlmP)

    double swsh_eigenvalue(int &s, int &l, int &m, double &g)

    cpp_complex[double] Yslm(int &s, int &l, int &m, double &th, double &ph)
    double Yslm(int &s, int &l, int &m, double &th) except +
    double Yslm_derivative(int &s, int &l, int &m, double &th)
    double Yslm_derivative2(int &s, int &l, int &m, double &th)

def _YslmCy(int s, int l, int m, double theta):
  return Yslm(s, l, m, theta)

def _YslmCy_derivative(int s, int l, int m, double theta):
  return Yslm_derivative(s, l, m, theta)

def _YslmCy_derivative2(int s, int l, int m, double theta):
  return Yslm_derivative2(s, l, m, theta)

def _clebschCy(int j1, int j2, int j, int m1, int m2, int m):
  return clebsch(j1, j2, j, m1, m2, m)

def _w3jCy(int j1, int j2, int j, int m1, int m2, int m):
  return w3j(j1, j2, j, m1, m2, m)

def _swsh_eigenvalueCy(int s, int l, int m, double g):
  return swsh_eigenvalue(s, l, m, g)

# arbitrary-theta evaluation from a precomputed coupling vector (no re-solve)
def _SslmCy_bvec(int s, int l, int m, double g, vector[double] bvec, double th):
  return Sslm(s, l, m, g, bvec, th)

def _SslmCy_derivative_bvec(int s, int l, int m, double g, vector[double] bvec, double th):
  return Sslm_derivative(s, l, m, g, bvec, th)

def _SslmCy_secondDerivative(int s, int l, int m, double g, double la, double th, double Slm, double SlmP):
  return Sslm_secondDerivative(s, l, m, g, la, th, Slm, SlmP)


cdef class _SpinWeightedHarmonic:
    """Cython wrapper over the C++ SpinWeightedHarmonic (grid-based spheroidal harmonic).
    Precomputes S, S', S'' on the supplied theta grid at construction."""
    cdef SpinWeightedHarmonic *swshcpp

    def __cinit__(self, int s, int l, int m, double gamma,
                  np.ndarray[ndim=1, dtype=np.float64_t] theta not None):
        cdef int n = theta.shape[0]
        cdef vector[double] thvec = vector[double](n)
        thvec.assign(&theta[0], &theta[0] + n)
        self.swshcpp = new SpinWeightedHarmonic(s, l, m, gamma, thvec)
        if self.swshcpp == NULL:
            raise MemoryError('Not enough memory.')
        self.swshcpp.generateSolutionsAndDerivatives()

    def __dealloc__(self):
        del self.swshcpp

    @property
    def spinweight(self): return self.swshcpp.getSpinWeight()
    @property
    def spheroidalmode(self): return self.swshcpp.getSpheroidalModeNumber()
    @property
    def azimuthalmode(self): return self.swshcpp.getAzimuthalModeNumber()
    @property
    def spheroidicity(self): return self.swshcpp.getSpheroidicity()
    @property
    def eigenvalue(self): return self.swshcpp.getEigenvalue()
    @property
    def mincouplingmode(self): return self.swshcpp.getMinCouplingModeNumber()
    @property
    def maxcouplingmode(self): return self.swshcpp.getMaxCouplingModeNumber()
    def couplingcoefficient(self, int l): return self.swshcpp.getCouplingCoefficient(l)
    @property
    def couplingcoefficients(self): return np.array(self.swshcpp.getCouplingCoefficient())
    @property
    def arguments(self): return np.array(self.swshcpp.getArguments())
    @property
    def solutions(self): return np.array(self.swshcpp.getSolution())
    @property
    def derivatives(self): return np.array(self.swshcpp.getDerivative())
    @property
    def secondderivatives(self): return np.array(self.swshcpp.getSecondDerivative())