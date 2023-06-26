from cybhpt_full import circular_redshift, run_tests, metric_coefficients_cython_ORG, metric_coefficients_cython_IRG
from cybhpt_full import RedshiftCoefficients as RedshiftCoefficientsCython
from cybhpt_full import SphericalHarmonicCoupling as SphericalHarmonicCouplingCython

def metric_coefficients_ORG(ai, bi, nt, nr, nz, np, a, r, z):
    return metric_coefficients_cython_ORG(ai, bi, nt, nr, nz, np, a, r, z)

def metric_coefficients_IRG(ai, bi, nt, nr, nz, np, a, r, z):
    return metric_coefficients_cython_IRG(ai, bi, nt, nr, nz, np, a, r, z)

class RedshiftCoefficients:
    def __init__(self, gauge, geo):
        self.base = RedshiftCoefficientsCython(gauge, geo.base)

    def __call__(self, Ni, ai, bi, ci, di, jr = 0, jz = 0):
        return self.base(Ni, ai, bi, ci, di, jr, jz)
    
class SphericalHarmonicCoupling:
    def __init__(self, lmax, m):
        self.base = SphericalHarmonicCouplingCython(lmax, m)
        self.azimuthalmode = self.base.azimuthalmode

    def zcouplingcoefficient(self, n, i, l):
        return self.base.zcouplingcoefficient(n, i, l)

    def dzcouplingcoefficient(self, n, i, l):
        return self.base.dzcouplingcoefficient(n, i, l)

def circular_redshift_export(filename, gauge, lmax, geo):
    return circular_redshift(filename, gauge, lmax, geo.base)

def run_unit_tests():
    run_tests()

from scipy.special import sph_harm
from scipy.special import binom
from scipy.special import factorial
import numpy as np

def fac(n):
    if n < 0:
        return 0
    return float(np.math.factorial(n))

def Yslm(s, l, m, th):
    if s == 0:
        return np.real(sph_harm(m, l, 0., th))
    elif s + m < 0:
        return (-1.)**(s+m)*YslmBase(-s, l, -m, np.cos(th))
    else:
        return YslmBase(s, l, m, np.cos(th))
    
def YslmBase(s, l, m, z):
    rmax = l - s
    pref = (0.5)**(l)*(-1.)**(m)*np.sqrt(factorial(l+m)/factorial(l+s)*factorial(l-m)/factorial(l-s)*(2*l+1)/(4.*np.pi))*np.sqrt(1. - z)**(s + m)*np.sqrt(1. + z)**(s - m)
    
    yslm = z*0.
    for r in range(0, rmax + 1):
        yslm += binom(l - s, r)*binom(l + s, r + s - m)*(z - 1.)**(rmax - r)*(z + 1.)**(r)
    
    return pref*yslm

def clebsch(l1, l2, l3, m1, m2, m3):
    return (-1.)**abs(l1 - l2 + m3)*np.sqrt(2*l3 + 1)*w3j(l1, l2, l3, m1, m2, -m3)

def w3j(l1, l2, l3, m1, m2, m3):
    if m1 + m2 + m3 != 0:
        return 0
    elif abs(l1 - l2) > l3:
        return 0
    elif l1 + l2 < l3:
        return 0
    
    if abs(m1) > l1:
        return 0
    elif abs(m2) > l2:
        return 0
    elif abs(m3) > l3:
        return 0
    
    sumTerm = w3j_tsum(l1, l2, l3, m1, m2, m3)
    if sumTerm == 0:
        return 0
    sumSign = np.sign(sumTerm)
    tempLog = 0.5*(np.log(fac(l1 + m1)) + np.log(fac(l2 + m2)) + np.log(fac(l3 + m3)))
    tempLog += 0.5*(np.log(fac(l1 - m1)) + np.log(fac(l2 - m2)) + np.log(fac(l3 - m3)))
    tempLog += np.log(triangle_coeff(l1, l2, l3))
    tempLog += np.log(abs(sumTerm))
    
    temp = sumSign*np.exp(tempLog)
    temp *= (-1.)**abs(l1-l2-m3)
    
    return temp
    
def w3j_tsum(l1, l2, l3, m1, m2, m3):
    t_min_num = w3j_t_min(l1, l2, l3, m1, m2, m3)
    t_max_num = w3j_t_max(l1, l2, l3, m1, m2, m3)
    x = 0
    if t_max_num < t_min_num:
        t_max_num = t_min_num

    for t in range(t_min_num - 1, t_max_num + 2):
        term = (fac(t)*fac(l3 - l2 + m1 + t)*fac(l3 - l1 - m2 + t)
            *fac(l1 + l2 - l3 - t)*fac(l1 - t - m1)*fac(l2 - t + m2))
        if term > 0:
            x += (-1.)**t/term
    
    return x

def w3j_t_min(l1, l2, l3, m1, m2, m3):
    temp = 0
    
    comp = l3 - l2 + m1
    if temp + comp < 0:
        temp = -comp
    comp = l3 - l1 - m2
    if temp + comp < 0:
        temp = -comp
        
    return temp

def w3j_t_max(l1, l2, l3, m1, m2, m3):
    temp = 1
    comp = l1 + l2 - l3
    if comp - temp > 0:
        temp = comp
    comp = l1 - m1
    if comp - temp > 0:
        temp = comp
    comp = l2 + m2
    if comp - temp > 0:
        temp = comp
        
    return temp

def triangle_coeff(l1, l2, l3):
    return np.sqrt(fac(l1 + l2 - l3)*fac(l3 + l1 - l2)*fac(l2 + l3 - l1)/fac(l1 + l2 + l3 + 1))