from cybhpt_full import RadialTeukolsky as RadialTeukolskyCython
import numpy as np

class RadialTeukolsky:

    def __init__(self, s, l, m, a, omega, r):
        self.radialpoints = r
        self.base = RadialTeukolskyCython(a, s, l, m, omega, self.radialpoints)

    @property
    def blackholespin(self):
        return self.base.blackholespin
    
    @property
    def spinweight(self):
        return self.base.spinweight
    @property
    def s(self):
        return self.spinweight

    @property
    def spheroidalmode(self):
        return self.base.spheroidalmode

    @property
    def j(self):
        return self.spheroidalmode

    @property
    def azimuthalmode(self):
        return self.base.azimuthalmode
    
    @property
    def m(self):
        return self.azimuthalmode

    @property
    def frequency(self):
        return self.base.frequency
    
    @property
    def mode_frequency(self):
        return self.frequency

    @property
    def omega(self):
        return self.frequency

    @property
    def eigenvalue(self):
        return self.base.eigenvalue
    
    def solveboundarycondition(self, method):
        self.base.solve_bc(method)

    def setboundarycondition(self, bc, R, Rp, r):
        self.base.set_bc(bc, R, Rp, r)

    def solve(self, method = "AUTO", bc = None):
        if bc is None:
            self.base.solve(method, "None")
        else:
            self.base.solve(method, bc)

    def flipspinweight(self):
        self.base.flip_spinweight()

    def radialpoint(self, pos):
        return self.base.radialpoint(pos)
    
    def boundarypoint(self, bc):
        return self.base.boundarypoint(bc)
    
    def boundarysolution(self, bc):
        return self.base.boundarysolution(bc)
    
    def boundaryderivative(self, bc):
        return self.base.boundaryderivative(bc)
    
    def radialsolution(self, bc, pos):
        return self.base.solution(bc, pos)
    
    def radialderivative(self, bc, pos):
        return self.base.derivative(bc, pos)
    
    def radialderivative2(self, bc, pos):
        return self.base.derivative2(bc, pos)