from libcpp.vector cimport vector
from libcpp.complex cimport complex as cpp_complex
import numpy as np
cimport numpy as np

cdef extern from "gsl/gsl_errno.h":
    void gsl_set_error_handler_off()

# call at import-time
gsl_set_error_handler_off()

cdef extern from "radialsolver.hpp":
    cdef enum SolutionMethod:
        AUTO, MST, ASYM, HBL, GSN, TEUK

    cdef enum BoundaryCondition:
        In, Up

    cdef cppclass Result:
        Result(cpp_complex[double], cpp_complex[double])
        cpp_complex[double] getValue()
        cpp_complex[double] getAccuracy()


    cdef cppclass RadialTeukolskyCPP "RadialTeukolsky":
        RadialTeukolskyCPP(double a, int s, int l, int m, double omega, vector[double] r) except +
        # RadialTeukolsky(double a, int s, int l, int m, double omega, double la, vector[double] r)  except +

        double getBlackHoleSpin()
        int getSpinWeight()
        int getSpheroidalModeNumber()
        int getAzimuthalModeNumber()
        double getModeFrequency()
        double getSpinWeightedSpheroidalEigenvalue()

        void generateRetardedBoundaryConditions(SolutionMethod method)
        void setBoundaryConditions(BoundaryCondition bc, cpp_complex[double] R, cpp_complex[double] Rp, double r)
        void generateSolutions(SolutionMethod method)
        void generateSolutions(BoundaryCondition bc, SolutionMethod method)
        int resampleSolutions(vector[double] radialSamples)
        void flipSpinWeight()

        vector[double] getRadialPoints()
        double getRadialPoints(int pos)
        double getBoundaryPoint(BoundaryCondition bc)
        Result getBoundarySolution(BoundaryCondition bc)
        Result getBoundaryDerivative(BoundaryCondition bc)

        vector[cpp_complex[double]] getSolution(BoundaryCondition bc)
        vector[cpp_complex[double]] getDerivative(BoundaryCondition bc)
        vector[cpp_complex[double]] getSecondDerivative(BoundaryCondition bc)
        cpp_complex[double] getSolution(BoundaryCondition bc, int pos)
        cpp_complex[double] getDerivative(BoundaryCondition bc, int pos)
        cpp_complex[double] getSecondDerivative(BoundaryCondition bc, int pos)

    cdef cppclass InterpolatedRadialTeukolskyCPP "InterpolatedRadialTeukolsky":
        InterpolatedRadialTeukolskyCPP(double a, int s, int l, int m, double omega,
                                       double rmin, double rmax,
                                       SolutionMethod method,
                                       bint solveIn, bint solveUp,
                                       double rtol) except +

        cpp_complex[double] evaluateSolution  (BoundaryCondition bc, double r)
        cpp_complex[double] evaluateDerivative(BoundaryCondition bc, double r)

        vector[double] getRadialPoints(BoundaryCondition bc)
        int getSampleCount(BoundaryCondition bc)

        double getBlackHoleSpin()
        int getSpinWeight()
        int getSpheroidalModeNumber()
        int getAzimuthalModeNumber()
        double getModeFrequency()
        double getSpinWeightedSpheroidalEigenvalue()

    void flip_spin_of_radial_teukolsky_TS(cpp_complex[double] &RinFlip, cpp_complex[double] &RinPFlip, BoundaryCondition bc, int s, int m, double a, double omega, double la, double r, cpp_complex[double] Rin, cpp_complex[double] RinP)
    double teukolsky_starobinsky_constant(int s, int m, double a, double omega, double lambdaCH)
    double teukolsky_starobinsky_constant_D(int m, double a, double omega, double lambdaCH)
    cpp_complex[double] teukolsky_starobinsky_complex_constant(int j, int m, double a, double omega, double lambdaCH)
    cpp_complex[double] teukolsky_starobinsky_amplitude(BoundaryCondition bc, int s, int m, double a, double omega, double lambdaCH)

cdef extern from "monodromy.hpp":
    cpp_complex[double] nu_solver_monodromy(int s, int l, int m, double q, double eps, double la)

cdef extern from "nusolver.hpp":
    cpp_complex[double] nu_solver(double q, int s, int l, int m, double epsilon)

cdef extern from "hypergeo_f.hpp":
    cpp_complex[double] hypergeo_2F1_complex(cpp_complex[double] a, cpp_complex[double] b, cpp_complex[double] c, cpp_complex[double] x)

def hypergeo_2F1(cpp_complex[double] a, cpp_complex[double] b, cpp_complex[double] c, cpp_complex[double] x):
    return hypergeo_2F1_complex(a, b, c, x)

def renormalized_angular_momentum(int s, int l, int m, double a, double omega):
    return nu_solver(a, s, l, m, 2.*omega)

def renormalized_angular_momentum_monodromy(int s, int l, int m, double a, double omega, double la):
    return nu_solver_monodromy(s, l, m, a, 2.*omega, la)

cdef dict bc_dict = {
    "In" : BoundaryCondition.In,
    "Up" : BoundaryCondition.Up
}

cdef dict method_dict = {
    "AUTO" : SolutionMethod.AUTO, 
    "MST" : SolutionMethod.MST, 
    "ASYMP" : SolutionMethod.ASYM, 
    "HBL" : SolutionMethod.HBL, 
    "GSN" : SolutionMethod.GSN, 
    "TEUK" : SolutionMethod.TEUK
}

def available_methods():
    """
    Returns a list of available solution methods.
    """
    return list(method_dict.keys())

def flip_spin_of_solutions(unicode bc, int s, double a, int m, double omega, double la, double r, cpp_complex[double] R, cpp_complex[double] Rp):
    cdef cpp_complex[double] R0, RP0
    flip_spin_of_radial_teukolsky_TS(R0, RP0, str_to_bc(bc), s, m, a, omega, la, r, R, Rp)
    return (R0, RP0)

cdef BoundaryCondition str_to_bc(unicode bc_str):
    if bc_str in bc_dict.keys():
        return bc_dict[bc_str]
    else:
        raise ValueError("{} is not a supported boundary condition.".format(bc_str))

cdef SolutionMethod str_to_method(unicode method_str):
    if method_str in method_dict.keys():
        return method_dict[method_str]
    else:
        raise ValueError("{} is not a supported solution method.".format(method_str))

def teukolsky_starobinsky_transformation_amplitude(unicode bc, int s, int m, double a, double omega, double lambdaCH):
    return teukolsky_starobinsky_amplitude(str_to_bc(bc), s, m, a, omega, lambdaCH)

def teukolsky_starobinsky_const(int j, int m, double a, double omega, double lambdaCH):
    return teukolsky_starobinsky_complex_constant(j, m, a, omega, lambdaCH)

def teukolsky_starobinsky_const_squared(int s, int m, double a, double omega, double lambdaCH):
    return teukolsky_starobinsky_constant(s, m, a, omega, lambdaCH)

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


cdef class InterpolatedRadialTeukolsky:
    cdef InterpolatedRadialTeukolskyCPP *interpcpp

    def __cinit__(self, double a, int s, int l, int m, double omega,
                  double rmin, double rmax,
                  unicode method=u"AUTO",
                  bint solve_in=True, bint solve_up=True,
                  double rtol=1e-10):
        self.interpcpp = new InterpolatedRadialTeukolskyCPP(
            a, s, l, m, omega, rmin, rmax,
            str_to_method(method), solve_in, solve_up, rtol)
        if self.interpcpp == NULL:
            raise MemoryError('Not enough memory.')

    def __dealloc__(self):
        del self.interpcpp

    @property
    def blackholespin(self):
        return self.interpcpp.getBlackHoleSpin()

    @property
    def spinweight(self):
        return self.interpcpp.getSpinWeight()

    @property
    def s(self):
        return self.spinweight

    @property
    def spheroidalmode(self):
        return self.interpcpp.getSpheroidalModeNumber()

    @property
    def j(self):
        return self.spheroidalmode

    @property
    def azimuthalmode(self):
        return self.interpcpp.getAzimuthalModeNumber()

    @property
    def m(self):
        return self.azimuthalmode

    @property
    def frequency(self):
        return self.interpcpp.getModeFrequency()

    @property
    def mode_frequency(self):
        return self.frequency

    @property
    def omega(self):
        return self.frequency

    @property
    def eigenvalue(self):
        return self.interpcpp.getSpinWeightedSpheroidalEigenvalue()

    def nsamples(self, unicode bc):
        """Return the number of interpolation nodes for the given boundary condition."""
        return self.interpcpp.getSampleCount(str_to_bc(bc))

    def radialpoints(self, unicode bc):
        """Return the interpolation nodes for the given boundary condition."""
        cdef vector[double] r = self.interpcpp.getRadialPoints(str_to_bc(bc))
        return np.array(r)

    def solution(self, unicode bc, double r):
        """Evaluate R(r) for the given boundary condition."""
        return self.interpcpp.evaluateSolution(str_to_bc(bc), r)

    def derivative(self, unicode bc, double r):
        """Evaluate R'(r) for the given boundary condition."""
        return self.interpcpp.evaluateDerivative(str_to_bc(bc), r)

    def solutions(self, unicode bc, np.ndarray[ndim=1, dtype=np.float64_t] r not None):
        """Evaluate R(r) at each point in the array r."""
        cdef int n = r.shape[0]
        cdef np.ndarray[ndim=1, dtype=np.complex128_t] result = np.zeros(n, dtype=np.complex128)
        cdef BoundaryCondition bce = str_to_bc(bc)
        cdef cpp_complex[double] val
        cdef int i
        for i in range(n):
            val = self.interpcpp.evaluateSolution(bce, r[i])
            result[i] = complex(val.real(), val.imag())
        return result

    def derivatives(self, unicode bc, np.ndarray[ndim=1, dtype=np.float64_t] r not None):
        """Evaluate R'(r) at each point in the array r."""
        cdef int n = r.shape[0]
        cdef np.ndarray[ndim=1, dtype=np.complex128_t] result = np.zeros(n, dtype=np.complex128)
        cdef BoundaryCondition bce = str_to_bc(bc)
        cdef cpp_complex[double] val
        cdef int i
        for i in range(n):
            val = self.interpcpp.evaluateDerivative(bce, r[i])
            result[i] = complex(val.real(), val.imag())
        return result