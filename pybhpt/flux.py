from cybhpt_full import FluxList as FluxListCython
from cybhpt_full import flux as FluxCython

class FluxList:
<<<<<<< HEAD
    def __init__(self, fluxlist=None):
        if fluxlist is None:
            self.fluxlist = [{"I": 0., "H": 0}, {"I": 0., "H": 0}, {"I": 0., "H": 0}]
        else:
            self.fluxlist = fluxlist

    def __add__(self, fluxlist2):
        for i in range(3):
            for bc in ["H", "I"]:
                self.fluxlist[0][bc] += fluxlist2[0][bc]
    
    @property
    def energy(self):
        return self.fluxlist.energy

    @property
    def angularmomentum(self):
        return self.fluxlist.angularmomentum

    @property
    def carterconstant(self):
        return self.fluxlist.carterconstant
=======
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
>>>>>>> e2354353061328e127ed9ca05c77018a631e6945

    
class FluxMode:
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
<<<<<<< HEAD
        return [self.energy, self.angularmomentum, self.carterconstant]
=======
        return FluxList([self.energy, self.angularmomentum, self.carterconstant])
>>>>>>> e2354353061328e127ed9ca05c77018a631e6945

    @property
    def horizonfluxes(self):
        return [self.energy["H"], self.angularmomentum["H"], self.carterconstant["H"]]
    
    @property
    def infinityfluxes(self):
        return [self.energy["I"], self.angularmomentum["I"], self.carterconstant["I"]]
    
    @property
    def totalfluxes(self):
        return [self.energy["I"] + self.energy["H"], self.angularmomentum["I"] + self.angularmomentum["H"], self.carterconstant["I"] + self.carterconstant["H"]]