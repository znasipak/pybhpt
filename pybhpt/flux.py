"""
Docstring for pybhpt.flux

Module containing classes and functions for computing fluxes of energy, angular momentum, and Carter constant in Kerr spacetime.
"""

from cybhpt_full import flux as _FluxCython
from cybhpt_full import _ELQdot_to_pexdot_wrapper, _ELQdot_to_pexdot_array_wrapper
import numpy as np

class FluxList:
    """A class for storing fluxes of energy, angular momentum, and Carter constant.

    Parameters
    ----------
    fluxes : list of dicts, optional
        A list containing three dictionaries, each with keys "I" and "H" representing the fluxes at infinity and on the horizon, respectively.
        If not provided, initializes with zero fluxes for all three components.

    Attributes
    ----------
    fluxes : list of dicts
        A list containing three dictionaries, each with keys "I" and "H" representing the fluxes at infinity and on the horizon, respectively.
        The first dictionary corresponds to energy flux,
        the second to angular momentum flux, and the third to Carter constant flux.

    Methods
    -------
    __call__():
        Returns the list of fluxes.
    """
    def __init__(self, fluxes=None):
        if fluxes is None:
            self.fluxes = [{"I": 0., "H": 0.}, {"I": 0., "H": 0.}, {"I": 0., "H": 0.}]
        elif len(fluxes) != 3:
            self.fluxes = [{"I": 0., "H": 0.}, {"I": 0., "H": 0.}, {"I": 0., "H": 0.}]
        else:
            self.fluxes = fluxes
    
    @property
    def energy(self):
        return self.fluxes[0]

    @property
    def angularmomentum(self):
        return self.fluxes[1]

    @property
    def carterconstant(self):
        return self.fluxes[2]
    
    def __call__(self):
        return self.fluxes

    
class FluxMode:
    """A class for computing fluxes of energy, angular momentum, and Carter constant for a given Teukolsky mode.

    Parameters
    ----------
    geo : KerrGeodesic class instance
        KerrGeodesic object containing the background motion of the point-particle source.

    teuk : TeukolskyMode object
        TeukolskyMode object containing mode solutions to the point-particle-sourced Teukolsky equation.

    Attributes
    ----------
    energy : dict
        A dictionary containing the energy flux at infinity and on the horizon.
        The keys are "I" for infinity and "H" for the horizon.
    angularmomentum : dict
        A dictionary containing the angular momentum flux at infinity and on the horizon.
        The keys are "I" for infinity and "H" for the horizon.
    carterconstant : dict
        A dictionary containing the Carter constant flux at infinity and on the horizon.
        The keys are "I" for infinity and "H" for the horizon.
    fluxes : FluxList
        A FluxList object containing the energy, angular momentum, and Carter constant fluxes.
    horizonfluxes : list
        A list containing the energy, angular momentum, and Carter constant fluxes on the horizon.
    infinityfluxes : list
        A list containing the energy, angular momentum, and Carter constant fluxes at infinity.
    totalfluxes : list
        A list containing the total fluxes, which are the sum of the horizon and infinity fluxes.
        """
    def __init__(self, geo, teuk):
        self.base = _FluxCython(teuk.spinweight, geo.base, teuk.base)
        self.apex = np.array([geo.a, geo.p, geo.e, geo.x])

    @property
    def energy(self):
        return self.base.energy

    @property
    def angularmomentum(self):
        return self.base.angularmomentum

    @property
    def carterconstant(self):
        return self.base.carterconstant
    
    @property
    def Edot(self):
        return self.energy
    @property
    def Ldot(self):
        return self.angularmomentum
    @property
    def Qdot(self):
        return self.carterconstant
    
    @property
    def fluxes(self):
        return FluxList([self.energy, self.angularmomentum, self.carterconstant])

    @property
    def horizonfluxes(self):
        return [self.energy["H"], self.angularmomentum["H"], self.carterconstant["H"]]
    
    @property
    def infinityfluxes(self):
        return [self.energy["I"], self.angularmomentum["I"], self.carterconstant["I"]]
    
    @property
    def totalfluxes(self):
        return [self.energy["I"] + self.energy["H"], self.angularmomentum["I"] + self.angularmomentum["H"], self.carterconstant["I"] + self.carterconstant["H"]]

def transform_ELQ_fluxes_to_pex(a, p, e, x, Edot, Lzdot, Qdot):
    """Transform fluxes of energy, angular momentum, and Carter constant to fluxes of semi-latus rectum, eccentricity, and inclination.

    Parameters
    ----------
    a : float or np.ndarray
        Spin parameter of the Kerr black hole.
    p : float or np.ndarray
        Semi-latus rectum of the orbit.
    e : float or np.ndarray
        Eccentricity of the orbit.
    x : float or np.ndarray
        Cosine of the inclination angle of the orbit.
    Edot : float or np.ndarray
        Flux of energy.
    Lzdot : float or np.ndarray
        Flux of angular momentum.
    Qdot : float or np.ndarray
        Flux of Carter constant.

    Returns
    -------
    pdot : float or np.ndarray
        Flux of semi-latus rectum.
    edot : float or np.ndarray
        Flux of eccentricity.
    xdot : float or np.ndarray
        Flux of inclination.
    """
    if isinstance(p, np.ndarray):
        assert isinstance(e, np.ndarray) and isinstance(x, np.ndarray) and isinstance(Edot, np.ndarray) and isinstance(Lzdot, np.ndarray) and isinstance(Qdot, np.ndarray), "If p is an array, e, x, Edot, Lzdot, and Qdot must also be arrays of the same shape."
        assert p.shape == e.shape == x.shape == Edot.shape == Lzdot.shape == Qdot.shape, "p, e, x, Edot, Lzdot, and Qdot must have the same shape."
        if not isinstance(a, np.ndarray):
            a = np.full_like(p, a)
        return _ELQdot_to_pexdot_array_wrapper(a, p, e, x, Edot, Lzdot, Qdot)
    else:  
        return _ELQdot_to_pexdot_wrapper(a, p, e, x, Edot, Lzdot, Qdot)