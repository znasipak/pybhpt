from libcpp.vector cimport vector
from libcpp.complex cimport complex as cpp_complex
from cython.operator cimport dereference
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
        self.fluxlistcpp.Edot.infinity += Edot
        self.fluxlistcpp.Ldot.infinity += Ldot
        self.fluxlistcpp.Qdot.infinity += Qdot

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


def flux(int s, KerrGeodesic geo, TeukolskyMode teuk):
    cdef FluxListCPP fluxescpp = flux_mode(s, dereference(geo.geocpp), dereference(teuk.teukcpp), include_minus_m = 0)
    fluxes = FluxList()
    fluxes.set_fluxes(fluxescpp)
    return fluxes