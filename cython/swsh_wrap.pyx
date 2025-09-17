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

def YslmCy(int s, int l, int m, double theta):
  return Yslm(s, l, m, theta)

def YslmCy_derivative(int s, int l, int m, double theta):
  return Yslm_derivative(s, l, m, theta)

def clebschCy(int j1, int j2, int j, int m1, int m2, int m):
  return clebsch(j1, j2, j, m1, m2, m)

def w3jCy(int j1, int j2, int j, int m1, int m2, int m):
  return w3j(j1, j2, j, m1, m2, m)