from libcpp.vector cimport vector
from libcpp.complex cimport complex as cpp_complex
from cython.operator cimport dereference
from libcpp.string cimport string as cpp_string
import numpy as np
cimport numpy as np

include "teukolsky_wrap.pyx"

cdef extern from "fluxes.hpp": 
    cdef cppclass FluxesCPP "Fluxes":
        FluxesCPP()
        FluxesCPP(double inf, double hor)
        double infinity
        double horizon

    cdef cppclass FluxListCPP "FluxList":
        FluxListCPP()
        FluxesCPP Edot
        FluxesCPP Ldot
        FluxesCPP Qdot

    FluxListCPP flux_mode(int s, GeodesicSource& geo, TeukolskyModeCPP& teukMode, int include_minus_m)

    FluxListCPP fluxes(GeodesicSource& geo)
    FluxListCPP fluxes(int s, GeodesicSource& geo)
    FluxListCPP flux_l(int s, int L, GeodesicSource& geo)
    FluxListCPP flux_lm(int s, int L, int m, GeodesicSource& geo)
    FluxListCPP flux_lm_sum(int s, int L, int m, GeodesicSource& geo)
    FluxListCPP flux_lmk(int s, int L, int m, int k, GeodesicSource& geo)
    FluxListCPP flux_lmk_sum(int s, int L, int m, int k, GeodesicSource& geo)

    void full_flux_parallel_l(int s, GeodesicSource geo, int modeMax, cpp_string dir)
    void full_flux_parallel_lm(GeodesicSource geo, int lMax, cpp_string dir)

    void ELQdot_to_pexdot(double &pdot, double &edot, double &xdot, double a, double p, double e, double x, double Edot, double Lzdot, double Qdot)
    void pexdot_to_ELQ_dot(double &Edot, double &Lzdot, double &Qdot, double a, double p, double e, double x, double pdot, double edot, double xdot)
    void ELQdot_to_pexdot(int n, double* pdot, double* edot, double* xdot, const double* a, const double* p, const double* e, const double* x, const double* Edot, const double* Lzdot, const double* Qdot)
    void pexdot_to_ELQ_dot(int n, double* Edot, double* Lzdot, double* Qdot, const double* a, const double* p, const double* e, const double* x, const double* pdot, const double* edot, const double* xdot)

cdef class FluxList:
    cdef FluxListCPP *fluxlistcpp

    def __cinit__(self):
        self.fluxlistcpp = new FluxListCPP()

    def __dealloc__(self):
        del self.fluxlistcpp

    def zero_fluxes(self):
        self.fluxlistcpp.Edot.infinity = 0.
        self.fluxlistcpp.Ldot.infinity = 0.
        self.fluxlistcpp.Qdot.infinity = 0.

        self.fluxlistcpp.Edot.horizon = 0.
        self.fluxlistcpp.Ldot.horizon = 0.
        self.fluxlistcpp.Qdot.horizon = 0.

    cdef set_fluxes(self, FluxListCPP fluxes):
        self.fluxlistcpp.Edot.infinity = fluxes.Edot.infinity
        self.fluxlistcpp.Ldot.infinity = fluxes.Ldot.infinity
        self.fluxlistcpp.Qdot.infinity = fluxes.Qdot.infinity

        self.fluxlistcpp.Edot.horizon = fluxes.Edot.horizon
        self.fluxlistcpp.Ldot.horizon = fluxes.Ldot.horizon
        self.fluxlistcpp.Qdot.horizon = fluxes.Qdot.horizon

    def set_infinity_fluxes(self, double Edot, double Ldot, double Qdot):
        self.fluxlistcpp.Edot.infinity = Edot
        self.fluxlistcpp.Ldot.infinity = Ldot
        self.fluxlistcpp.Qdot.infinity = Qdot

    def set_horizon_fluxes(self, double Edot, double Ldot, double Qdot):
        self.fluxlistcpp.Edot.horizon = Edot
        self.fluxlistcpp.Ldot.horizon = Ldot
        self.fluxlistcpp.Qdot.horizon = Qdot

    def add_infinity_fluxes(self, double Edot, double Ldot, double Qdot):
        self.fluxlistcpp.Edot.infinity += Edot
        self.fluxlistcpp.Ldot.infinity += Ldot
        self.fluxlistcpp.Qdot.infinity += Qdot

    def add_horizon_fluxes(self, double Edot, double Ldot, double Qdot):
        self.fluxlistcpp.Edot.horizon += Edot
        self.fluxlistcpp.Ldot.horizon += Ldot
        self.fluxlistcpp.Qdot.horizon += Qdot

    def set_infinity_fluxes(self, double Edot, double Ldot, double Qdot):
        self.fluxlistcpp.Edot.infinity = Edot
        self.fluxlistcpp.Ldot.infinity = Ldot
        self.fluxlistcpp.Qdot.infinity = Qdot

    def add_fluxes(self, double EdotH, double LdotH, double QdotH, double EdotI, double LdotI, double QdotI):
        self.fluxlistcpp.Edot.horizon += EdotH
        self.fluxlistcpp.Ldot.horizon += LdotH
        self.fluxlistcpp.Qdot.horizon += QdotH

        self.fluxlistcpp.Edot.infinity += EdotI
        self.fluxlistcpp.Ldot.infinity += LdotI
        self.fluxlistcpp.Qdot.infinity += QdotI

    @property
    def energy(self):
        return {
            "I": self.fluxlistcpp.Edot.infinity,
            "H": self.fluxlistcpp.Edot.horizon
        }

    @property
    def angularmomentum(self):
        return {
            "I": self.fluxlistcpp.Ldot.infinity,
            "H": self.fluxlistcpp.Ldot.horizon
        }

    @property
    def carterconstant(self):
        return {
            "I": self.fluxlistcpp.Qdot.infinity,
            "H": self.fluxlistcpp.Qdot.horizon
        }


def flux(int s, KerrGeodesic geo, _TeukolskyMode teuk):
    cdef FluxListCPP fluxescpp = flux_mode(s, dereference(geo.geocpp), dereference(teuk.teukcpp), include_minus_m = 0)
    fluxes = FluxList()
    fluxes.set_fluxes(fluxescpp)
    return fluxes

def _ELQdot_to_pexdot_wrapper(double a, double p, double e, double x, double Edot, double Lzdot, double Qdot):
    cdef double pdot, edot, xdot
    pdot = 0.
    edot = 0.
    xdot = 0.
    ELQdot_to_pexdot(pdot, edot, xdot, a, p, e, x, Edot, Lzdot, Qdot)
    return np.array([pdot, edot, xdot])

def _ELQdot_to_pexdot_array_wrapper(double[:] a, double[:] p, double[:] e, double[:] x, double[:] Edot, double[:] Lzdot, double[:] Qdot):
    """
    Zero-copy interface using NumPy memoryviews.
    """
    cdef int n = a.shape[0] 

    # Pre-allocate output NumPy arrays
    pdot_out = np.empty(n, dtype=np.float64)
    edot_out = np.empty(n, dtype=np.float64)
    xdot_out = np.empty(n, dtype=np.float64) 

    # Cast to memoryviews to get raw pointers
    cdef double[:] pdot_mv = pdot_out
    cdef double[:] edot_mv = edot_out
    cdef double[:] xdot_mv = xdot_out

    # Execute C++ loop
    ELQdot_to_pexdot(n, &pdot_mv[0], &edot_mv[0], &xdot_mv[0], &a[0], &p[0], &e[0], &x[0], &Edot[0], &Lzdot[0], &Qdot[0])

    return np.array([pdot_out, edot_out, xdot_out])

def _pexdot_to_ELQ_dot_wrapper(double a, double p, double e, double x, double pdot, double edot, double xdot):
    cdef double Edot, Lzdot, Qdot
    Edot = 0.
    Lzdot = 0.
    Qdot = 0.
    pexdot_to_ELQ_dot(Edot, Lzdot, Qdot, a, p, e, x, pdot, edot, xdot)
    return np.array([Edot, Lzdot, Qdot])

def _pexdot_to_ELQ_dot_array_wrapper(double[:] a, double[:] p, double[:] e, double[:] x, double[:] pdot, double[:] edot, double[:] xdot):
    """
    Zero-copy interface using NumPy memoryviews.
    """
    cdef int n = a.shape[0] 

    # Pre-allocate output NumPy arrays
    Edot_out = np.empty(n, dtype=np.float64)
    Lzdot_out = np.empty(n, dtype=np.float64)
    Qdot_out = np.empty(n, dtype=np.float64)  

    # Cast to memoryviews to get raw pointers
    cdef double[:] Edot_mv = Edot_out
    cdef double[:] Lzdot_mv = Lzdot_out
    cdef double[:] Qdot_mv = Qdot_out

    # Execute C++ loop
    pexdot_to_ELQ_dot(n, &Edot_mv[0], &Lzdot_mv[0], &Qdot_mv[0], &a[0], &p[0], &e[0], &x[0], &pdot[0], &edot[0], &xdot[0])

    return np.array([Edot_out, Lzdot_out, Qdot_out])