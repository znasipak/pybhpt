from cybhpt_full import FluxList as FluxListCython
from cybhpt_full import flux as FluxCython

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

    # def __add__(self, fluxlist2):
    #     for i in range(3):
    #         for bc in ["H", "I"]:
    #             self.fluxes[i][bc] += fluxlist2.fluxes[i][bc]
    
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
        self.base = FluxCython(teuk.spinweight, geo.base, teuk.base)

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
