from cybhpt_full import metric_coefficients_cython_ORG, metric_coefficients_cython_IRG
import numpy as np
from pybhpt.hertz import available_gauges
from pybhpt.swsh import Yslm, spin_operator_normalization

def gauge_check(gauge):
    if gauge not in available_gauges:
        TypeError("{} is not a supported gauge.".format(gauge))

S0_gauges = ["IRG", "ASAAB0", "SAAB0"]
S4_gauges = ["ORG", "ASAAB4", "SAAB4"]

def metric_coefficients_S4_ab(ai, bi, nt, nr, nz, nph, a, r, z):
    return metric_coefficients_cython_ORG(ai, bi, nt, nr, nz, nph, a, r, z)

def metric_coefficients_S0_ab(ai, bi, nt, nr, nz, nph, a, r, z):
    return metric_coefficients_cython_IRG(ai, bi, nt, nr, nz, nph, a, r, z)

def metric_coefficients_S0(a, b, c, d, q, rvals, zvals):
    h22 = np.array([[metric_coefficients_S0_ab(2, 2, a, b, c, d, q, r, z) for z in zvals] for r in rvals])
    h24 = np.array([[metric_coefficients_S0_ab(2, 4, a, b, c, d, q, r, z) for z in zvals] for r in rvals])
    h44 = np.array([[metric_coefficients_S0_ab(4, 4, a, b, c, d, q, r, z) for z in zvals] for r in rvals])
    return np.array([2.*h22, h24, h44])

def metric_coefficients_S4(a, b, c, d, q, rvals, zvals):
    h11 = np.array([[metric_coefficients_S4_ab(1, 1, a, b, c, d, q, r, z) for z in zvals] for r in rvals])
    h13 = np.array([[metric_coefficients_S4_ab(1, 3, a, b, c, d, q, r, z) for z in zvals] for r in rvals])
    h33 = np.array([[metric_coefficients_S4_ab(3, 3, a, b, c, d, q, r, z) for z in zvals] for r in rvals])
    return np.array([2.*h11, h13, h33])

def metric_coefficients(gauge, a, b, c, d, q, rvals, zvals):
    gauge_check(gauge)
    if gauge in S0_gauges:
        return metric_coefficients_S0(a, b, c, d, q, rvals, zvals)
    else:
        return metric_coefficients_S4(a, b, c, d, q, rvals, zvals)

class MetricCoefficients:
    def __init__(self, gauge, q, r, th):
        gauge_check(gauge)
        self.gauge = gauge
        rsamples = r.shape[0]
        zsamples = th.shape[0]
        self.blackholespin = q
        self.radialpoints = r
        self.polarpoints = th
        if gauge in S0_gauges:
            self.storedcomponents = {(2, 2): 0, (2, 4): 1, (4, 4): 2}
            self.conjugatecomponents = {(2, 3): 1, (3, 3): 2}
        else:
            self.storedcomponents = {(1, 1): 0, (1, 3): 1, (3, 3): 2}
            self.conjugatecomponents = {(1, 4): 1, (4, 4): 2}
        z = np.cos(th)
        z[np.abs(z) < 1.e-15] = 0.
        self.coeffs = np.zeros((3, 3, 3, 3, 3, rsamples, zsamples), dtype=np.complex128)
        self.zeros = np.zeros((rsamples, zsamples), dtype=np.complex128)
        for ai in range(3):
            for bi in range(3):
                for ci in range(3):
                    for di in range(3):
                        self.coeffs[ai, bi, ci, di] = metric_coefficients(self.gauge, ai, bi, ci, di, q, r, z)

    def hab(self, a, b, nt, nr, ns, nphi):
        if b < a:
            atemp = a
            a = b
            b = atemp
        if (a, b) in self.storedcomponents.keys():
            return self.coeffs[nt, nr, ns, nphi][self.storedcomponents[(a, b)]]
        elif (a, b) in self.conjugatecomponents.keys():
            return np.conj(self.coeffs[nt, nr, ns, nphi][self.conjugatecomponents[(a, b)]])
        else:
            return self.zeros

    def __call__(self, a, b, nt, nr, ns, nphi):
        return self.hab(a, b, nt, nr, ns, nphi)
    
def tetrad_project_l(a, r, z, mu):
    if mu == 0:
        return -1
    elif mu == 1:
        sigma = r**2 + a**2*z**2
        delta = r**2 - 2.*r + a**2
        return sigma/delta
    elif mu == 2:
        return 0.
    elif mu == 3:
        return a*(1. - z**2)
    else:
        return 0.
    
def tetrad_project_n(a, r, z, mu):
    sigma = r**2 + a**2*z**2
    delta = r**2 - 2.*r + a**2
    if mu == 0:
        return -0.5*delta/sigma
    elif mu == 1:
        return -0.5
    elif mu == 2:
        return 0.
    elif mu == 3:
        return 0.5*delta/sigma*a*(1. - z**2)
    else:
        return 0.
    
def tetrad_project_m(a, r, z, mu):
    rhobar = -1./(r + 1j*a*z)
    sigma = r**2 + a**2*z**2
    pref = - rhobar*np.sqrt(0.5*(1. - z**2))
    if mu == 0:
        return -1j*pref*a
    elif mu == 1:
        return 0.
    elif mu == 2:
        return -pref*sigma/(1. - z**2)
    elif mu == 3:
        return 1j*pref*(r**2 + a**2)
    else:
        return 0.
    
def kinnersley_tetrad_covector(b, q, r, z, mu):
    if b == 1:
        return tetrad_project_l(q, r, z, mu)
    elif b == 2:
        return tetrad_project_n(q, r, z, mu)
    elif b == 3:
        return tetrad_project_m(q, r, z, mu)
    elif b == 4:
        return np.conj(tetrad_project_m(q, r, z, mu))
    else:
        return 0.
    
def hmunu_BL(gauge, mu, nu, q, r, hab):
    gauge_check(gauge)
    if gauge in S0_gauges:
        e1mu = kinnersley_tetrad_covector(1, q, r, 0, mu)
        e1nu = kinnersley_tetrad_covector(1, q, r, 0, nu)
        e3mu = kinnersley_tetrad_covector(3, q, r, 0, mu)
        e3nu = kinnersley_tetrad_covector(3, q, r, 0, nu)
        h22 = e1mu*e1nu*hab[0]
        h24 = -e1mu*e3nu*hab[1] - e1nu*e3mu*hab[1]
        h44 = e3mu*e3nu*hab[2]
        return (h22.real + 2.*h24.real + 2.*h44.real)
    else:
        e1mu = kinnersley_tetrad_covector(2, q, r, 0, mu)
        e1nu = kinnersley_tetrad_covector(2, q, r, 0, nu)
        e3mu = kinnersley_tetrad_covector(4, q, r, 0, mu)
        e3nu = kinnersley_tetrad_covector(4, q, r, 0, nu)
        h11 = e1mu*e1nu*hab[0]
        h13 = -e1mu*e3nu*hab[1] - e1nu*e3mu*hab[1]
        h33 = e3mu*e3nu*hab[2]
        return (h11.real + 2.*h13.real + 2.*h33.real)

# class MetricMode:
#     def __init__(self, hertz):
#         self.phi = hertz

#     def construct_coefficients(self, r, th):
#         coeffs = MetricCoefficients(self.phi.gauge, self.phi.blackholespin, r, th)
#         return coeffs

#     def solve(self, r, th):
#         self.coeffs = self.construct_coefficients(r, th)
#         q = self.phi.blackholespin
#         omega = self.phi.frequency
#         m = self.phi.azimuthalmode
#         gauge = self.phi.gauge
#         s_sgn = 1
#         if gauge in S0_gauges:
#             s_sgn = -1
#         habl = []
#         lmax = self.phi.maxcouplingmode
#         lmin = self.phi.mincouplingmode
#         for l in range(lmin, lmax + 1):
#             hab = np.zeros((3, r.shape[0], th.shape[0]), dtype=np.complex128)
#             for ns in range(0, 3):
#                 Jterm = spin_operator_normalization(s_sgn*2, ns, l)
#                 yslm = Jterm*self.phi.couplingcoefficient(l)*np.real(Yslm(s_sgn*(2-ns), l, m, th))
#                 for nt in range(0, 3):
#                     for nph in range(0, 3):
#                         pref = (1.j*m)**nph*(-1.j*omega)**nt
#                         for nr in range(0, 3):
#                             if nt + nr + ns + nph <= 2:
#                                 hNabcds = self.coeffs(nt, nr, ns, nph)
#                                 habTerm = []
#                                 phiTimesY = np.outer(yslm, )
#                                 for hNabcd in hNabcds:
#                                     habTerm.append(-pref*hNabcd[0]*yslm*hertzMode(l, m, deriv=nr)[0])
#                                 hab += np.array(habTerm)
#             habl.append(hab)

#         self.hab = habl
