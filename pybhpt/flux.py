from cybhpt_full import flux

class Fluxes:
    def __init__(self, fluxlist):
        self.fluxlist = fluxlist

    @property
    def energyflux(self):
        return self.fluxlist.energyflux

    @property
    def angularmomentumflux(self):
        return self.fluxlist.angularmomentumflux

    @property
    def carterflux(self):
        return self.fluxlist.carterflux

def fluxes(s, geo, teuk):
    fluxlist = flux(s, geo.base, teuk.base)
    return Fluxes(fluxlist)