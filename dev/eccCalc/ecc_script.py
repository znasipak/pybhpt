from pybhpt.geo import KerrGeodesic
from pybhpt.radial import RadialTeukolsky
from pybhpt.teuk import TeukolskyMode
from pybhpt.hertz import HertzMode
from pybhpt.redshift import RedshiftCoefficients, Yslm
import numpy as np
import matplotlib.pyplot as plt

from scipy.special import factorial as fac
from pybhpt.redshift import w3j
from scipy.special import sph_harm

def muCoupling(s, l):
    """
    Eigenvalue for the spin-weighted spherical harmonic lowering operator
    Setting s -> -s gives the negative of the eigenvalue for the raising operator
    """
    if l + s < 0 or l - s + 1. < 0:
        return 0
    return np.sqrt((l - s + 1.)*(l + s))

def Asjlm(s, j, l, m):
    if s >= 0:
        return (-1.)**(m + s)*np.sqrt(4**s*fac(s)**2*(2*l + 1)*(2*j + 1)/fac(2*s))*w3j(s, l, j, 0, m, -m)*w3j(s, l, j, s, -s, 0)
    else:
        return (-1.)**(m)*np.sqrt(4**(-s)*fac(-s)**2*(2*l + 1)*(2*j + 1)/fac(-2*s))*w3j(-s, l, j, 0, m, -m)*w3j(-s, l, j, s, -s, 0)

def spin_operator_normalization(s, ns, l):
    s_sgn = np.sign(s)
    nmax1 = np.abs(s) + 1
    Jterm = 1.
    for ni in range(1, ns + 1):
        Jterm *= -s_sgn*muCoupling((nmax1-ni), l)
    return Jterm

def full_libration(xp):
    return np.concatenate((xp, np.flip(xp[1:-1])))

def construct_hertz_radial_components(hertz, samples = 2**6):
    samples_half = int(samples/2 + 1)
    hertzR0In = np.array([hertz.radialsolution('In', i) for i in range(samples_half)])
    hertzR0Up = np.array([hertz.radialsolution('Up', i) for i in range(samples_half)])
    hertzR1In = np.array([hertz.radialderivative('In', i) for i in range(samples_half)])
    hertzR1Up = np.array([hertz.radialderivative('Up', i) for i in range(samples_half)])
    hertzR2In = np.array([hertz.radialderivative2('In', i) for i in range(samples_half)])
    hertzR2Up = np.array([hertz.radialderivative2('Up', i) for i in range(samples_half)])
    hertzR = np.array(
        [[np.concatenate((hertzR0In, np.flip(hertzR0In)[1:-1])), np.concatenate((hertzR0Up, np.flip(hertzR0Up)[1:-1]))], 
        [np.concatenate((hertzR1In, np.flip(hertzR1In)[1:-1])), np.concatenate((hertzR1Up, np.flip(hertzR1Up)[1:-1]))], 
        [np.concatenate((hertzR2In, np.flip(hertzR2In)[1:-1])), np.concatenate((hertzR2Up, np.flip(hertzR2Up)[1:-1]))]])
    return hertzR

def construct_hertz_radial_components_circ(hertz):
    i = 0
    hertzR0In = hertz.radialsolution('In', i)
    hertzR0Up = hertz.radialsolution('Up', i)
    hertzR1In = hertz.radialderivative('In', i)
    hertzR1Up = hertz.radialderivative('Up', i)
    hertzR2In = hertz.radialderivative2('In', i)
    hertzR2Up = hertz.radialderivative2('Up', i)
    hertzR = np.array(
        [[hertzR0In, hertzR0Up],
        [hertzR1In, hertzR1Up],
        [hertzR2In, hertzR2Up]]
    )
    return hertzR

def construct_hertz_polar_spin_weighting(hertz, samples = 2**6):
    samples_half = int(samples/2 + 1)
    th = np.array([hertz.polarpoint(i) for i in range(samples_half)])
    thp = full_libration(th)
    sth = np.sin(thp)
    spinWeightingPolar = np.array(
        [pow(sth, -2),
        pow(sth, -1),
        pow(sth, 0)]
    )
    return spinWeightingPolar

def redshift_mode_inc(coeffs, hertz, l, samples = 2**6):
    m = hertz.azimuthalmode
    omega = hertz.frequency
    s_sgn = np.sign(hertz.spinweight)
    huu = np.zeros((samples, 2), dtype=np.complex128)
    hertzR = construct_hertz_radial_components_circ(hertz)
    spinWeightingPolar = construct_hertz_polar_spin_weighting(hertz, samples=samples)
    for ns in range(0, 3):
        scalar_coupling = 0.
        llmin = np.max([np.abs(m), 2, l - 3 + ns])
        llmax = np.max([llmin, l + 3 - ns])
        for ll in range(llmin, llmax + 1):
            Jterm = spin_operator_normalization(s_sgn*2, ns, ll)
            scalar_coupling += Jterm*hertz.couplingcoefficient(ll)*Asjlm(s_sgn*(2-ns), l, ll, m)
        for nt in range(0, 3):
            for nph in range(0, 3):
                pref = (1.j*m)**nph*(-1.j*omega)**nt
                for nr in range(0, 3):
                    if nt + nr + ns + nph <= 2:
                        hNabcd0 = np.array([coeffs(0, nt, nr, ns, nph, 0, i) for i in range(samples)])
                        hNabcd0 *= spinWeightingPolar[ns]
                        hNabcd = np.array([hNabcd0, hNabcd0]).T
                        huuTerm = -pref*hNabcd*scalar_coupling*hertzR[nr]
                        # if np.sum(huuTerm) == 0 and np.sum(hNabcd0) != 0:
                        #     print(nt, nr, ns, nph)
                        #     print(pref, spinWeightingPolar[ns], scalar_coupling, hertzR[nr])
                        huu += huuTerm
    return huu

def redshift_mode_ecc(coeffs, hertz, l, samples = 2**6):
    m = hertz.azimuthalmode
    omega = hertz.frequency
    s_sgn = np.sign(hertz.spinweight)
    ylm = np.real(sph_harm(m, l, 0., 0.5*np.pi))
    huu = np.zeros((2, samples), dtype=np.complex128)
    if np.abs(ylm) == 0.:
        return huu
    hertzR = construct_hertz_radial_components(hertz, samples)
    for ns in range(0, 3):
        scalar_coupling = 0.
        llmin = np.max([np.abs(m), 2, l - 3 + ns])
        llmax = np.max([llmin, l + 3 - ns])
        for ll in range(llmin, llmax + 1):
            Jterm = spin_operator_normalization(s_sgn*2, ns, ll)
            scalar_coupling += Jterm*hertz.couplingcoefficient(ll)*Asjlm(s_sgn*(2-ns), l, ll, m)
        for nt in range(0, 3):
            for nph in range(0, 3):
                pref = (1.j*m)**nph*(-1.j*omega)**nt
                for nr in range(0, 3):
                    if nt + nr + ns + nph <= 2:
                        hNabcd0 = [coeffs(0, nt, nr, ns, nph, i, 0) for i in range(samples)]
                        hNabcd = np.array([hNabcd0, hNabcd0])
                        huuTerm = -pref*hNabcd*scalar_coupling*hertzR[nr]*ylm
                        huu += huuTerm
    return huu

def redshift_mode_calc(coeffs, hertz, l):
    m = hertz.azimuthalmode
    omega = hertz.frequency
    s_sgn = 1
    if np.sign(hertz.spinweight) < 0:
        s_sgn = -1
    ylm = np.real(sph_harm(m, l, 0., 0.5*np.pi))
    huu = np.array([0., 0.], dtype=np.complex128)
    if np.abs(ylm) == 0.:
        return huu
    hertzRIn = [hertz.radialsolution('In', 0), hertz.radialderivative('In', 0), hertz.radialderivative2('In', 0)]
    hertzRUp = [hertz.radialsolution('Up', 0), hertz.radialderivative('Up', 0), hertz.radialderivative2('Up', 0)]
    for ns in range(0, 3):
        scalar_coupling = 0.
        llmin = np.max([np.abs(m), 2, l - 3 + ns])
        llmax = l + 3 - ns
        for ll in range(llmin, llmax + 1):
            Jterm = spin_operator_normalization(s_sgn*2, ns, ll)
            scalar_coupling += Jterm*hertz.couplingcoefficient(ll)*Asjlm(s_sgn*(2-ns), l, ll, m)
        # print(scalar_coupling, ns, l, hertz.spheroidalmode)
        for nt in range(0, 3):
            for nph in range(0, 3):
                pref = (1.j*m)**nph*(-1.j*omega)**nt
                for nr in range(0, 3):
                    if nt + nr + ns + nph <= 2:
                        hNabcd = coeffs(0, nt, nr, ns, nph)
                        if(abs(hNabcd) > 0.):
                            huuTerm = -pref*hNabcd*scalar_coupling*np.array([hertzRIn[nr], hertzRUp[nr]])*ylm
                            # if np.abs(huuTerm[0]) + np.abs(huuTerm[1]) > 0:
                            #     print(nt, nr, ns, nph, huuTerm)
                            huu += huuTerm
    return huu

from tqdm import tqdm

def redshift_calc(gauge, lmax, geo):
    if isinstance(gauge, str):
        if gauge == "SAAB":
            return redshift_calc(["SAAB0", "SAAB4"], lmax, geo)
        elif gauge == "ASAAB":
            return redshift_calc(["ASAAB0", "ASAAB4"], lmax, geo)
        gauges = [gauge]
    else:
        gauges = gauge
    coeffs = []
    for gauge in gauges:
        coeffs.append(RedshiftCoefficients(gauge, geo))
    huuYlm = np.zeros((lmax + 1, lmax + 1, 2), dtype=np.complex128)
    s = choose_spin_from_gauge(gauges[0])
    for m in tqdm(range(0, lmax + 1)):
        jmin = np.max([2, abs(m)])
        jmax = lmax + 4 + int(20*geo.blackholespin)
        for j in range(jmin, jmax + 1):
            teuk = TeukolskyMode(s, j, m, 0, 0, geo)
            teuk.solve(geo)
            psis = []
            for i, gauge in enumerate(gauges):
                psis.append(HertzMode(teuk, gauge))
                psis[i].solve()

            for l in range(abs(m), lmax + 1):
                # huuMode = redshift_mode_calc(coeffs, hertz, l)
                # huuYlm[l, m] += huuMode
                for coeff, psi in zip(coeffs, psis):
                        huuMode = redshift_mode_calc(coeff, psi, l)
                        if np.isnan(np.sum(huuMode)):
                            print("Error ", j, l, m)
                        else:
                            huuTemp = huuMode
                            # if l == j:
                            #     print(l, m, k)
                            #     print(orbit_average_polar(huuTemp, geo, axis = 0))
                            huuYlm[l, m] += huuTemp
    
    huuYl = np.zeros((lmax + 1, 2))
    for l in range(0, lmax + 1):
        huuYl[l] += 2.*np.real(huuYlm[l, 0])
        for m in range(1, l + 1):
            huuYl[l] +=  4.*np.real(huuYlm[l, m])

    return huuYl, huuYlm

def redshift_calc_ecc(gauge, lmax, geo, nrange = [-10, 10]):
    if isinstance(gauge, str):
        if gauge == "SAAB":
            return redshift_calc_ecc(["SAAB0", "SAAB4"], lmax, geo, nrange=nrange)
        elif gauge == "ASAAB":
            return redshift_calc_ecc(["ASAAB0", "ASAAB4"], lmax, geo, nrange=nrange)
        gauges = [gauge]
    else:
        gauges = gauge
    coeffs = []
    for gauge in gauges:
        coeffs.append(RedshiftCoefficients(gauge, geo))
    samples = 2*(geo.radialpoints.shape[0] - 1)
    huuYlm = np.zeros((lmax + 1, lmax + 1, 2, samples), dtype=np.complex128)
    s = choose_spin_from_gauge(gauges[0])
    for m in tqdm(range(0, lmax + 1)):
        jmin = np.max([abs(s), abs(m)])
        jmax = lmax + 4 + int(20)
        for j in range(jmin, jmax + 1):
            nmin = nrange[0]
            nmax = nrange[1]
            for n in range(nmin, nmax + 1):
                teuk = TeukolskyMode(s, j, m, 0, n, geo)
                teuk.solve(geo)
                psis = []
                for i, gauge in enumerate(gauges):
                    psis.append(HertzMode(teuk, gauge))
                    psis[i].solve()
                phaseUp = m*geo.azimuthalradial - teuk.frequency*geo.timeradial
                phaseDown = np.flip(-phaseUp)[1:-1]
                phase = np.concatenate((phaseUp, phaseDown)) - 2.*np.pi*n*np.arange(samples)/samples
                for l in range(abs(m), lmax + 1):
                    # huuMode = redshift_mode_ecc(coeffs, hertz, l, samples=samples)
                    # # if l == j and l < 5:
                    # #     print(l, m, n, np.mean(huuMode, axis = 1))
                    # if np.isnan(np.sum(huuMode)):
                    #     print(j, l, m, n)
                    # else:
                    #     huuYlm[l, m] += huuMode*np.exp(1.j*np.array([phase, phase]))
                    for coeff, psi in zip(coeffs, psis):
                        huuMode = redshift_mode_ecc(coeff, psi, l, samples=samples)
                        if np.isnan(np.sum(huuMode)):
                            print("Error ", j, l, m, n)
                        else:
                            huuTemp = huuMode*np.exp(1.j*np.array([phase, phase]))
                            # if l == j:
                            #     print(l, m, k)
                            #     print(orbit_average_polar(huuTemp, geo, axis = 0))
                            huuYlm[l, m] += huuTemp
    
    huuYl = np.zeros((lmax + 1, 2, samples))
    for l in range(0, lmax + 1):
        huuYl[l] += 2.*np.real(huuYlm[l, 0])
        for m in range(1, l + 1):
            huuYl[l] += 4.*np.real(huuYlm[l, m])

    return huuYl, huuYlm

def choose_spin_from_gauge(gauge):
    if gauge == "ORG" or gauge == "SAAB4" or gauge == "ASAAB4":
        s = 2
    else:
        s = -2
    return s

def redshift_calc_inc(gauge, lmax, geo, krange = [-10, 10]):
    if isinstance(gauge, str):
        if gauge == "SAAB":
            return redshift_calc_inc(["SAAB0", "SAAB4"], lmax, geo, krange=krange)
        elif gauge == "ASAAB":
            return redshift_calc_inc(["ASAAB0", "ASAAB4"], lmax, geo, krange=krange)
        gauges = [gauge]
    else:
        gauges = gauge
    coeffs = []
    for gauge in gauges:
        coeffs.append(RedshiftCoefficients(gauge, geo))
    samples = 2*(geo.polarpoints.shape[0] - 1)
    huuYlm = np.zeros((lmax + 1, lmax + 1, samples, 2), dtype=np.complex128)
    Ylms = np.zeros((lmax + 1, lmax + 1, samples, 2))
    th = geo.polarpoints
    thp = full_libration(th)
    for m in range(-lmax, lmax + 1):
        for l in range(m, lmax + 1):
            ylm0 = np.array([np.real(sph_harm(m, l, 0., th)) for th in thp])
            ylm = np.array([ylm0, ylm0]).T
            Ylms[l, m] = ylm
    s = choose_spin_from_gauge(gauges[0])
    for m in tqdm(range(0, lmax + 1)):
        jmin = np.max([abs(s), abs(m)])
        jmax = lmax + 4 + int(20*geo.blackholespin)
        for j in range(jmin, jmax + 1):
            kmin = krange[0]
            kmax = krange[1]
            for k in range(kmin, kmax + 1):
                teuk = TeukolskyMode(s, j, m, k, 0, geo)
                teuk.solve(geo)
                psis = []
                for i, gauge in enumerate(gauges):
                    psis.append(HertzMode(teuk, gauge))
                    psis[i].solve()
                phaseUp = m*geo.azimuthalpolar - teuk.frequency*geo.timepolar
                phaseDown = np.flip(-phaseUp)[1:-1]
                phase = np.concatenate((phaseUp, phaseDown)) - 2.*np.pi*k*np.arange(samples)/samples
                phaseArray = np.array([phase, phase]).T
                for l in range(abs(m), lmax + 1):
                    for coeff, psi in zip(coeffs, psis):
                        huuMode = redshift_mode_inc(coeff, psi, l, samples=samples)
                        if np.isnan(np.sum(huuMode)):
                            print("Error ", j, l, m, k)
                        else:
                            huuTemp = huuMode*np.exp(1.j*phaseArray)*Ylms[l, m]
                            # if l == j:
                            #     print(l, m, k)
                            #     print(orbit_average_polar(huuTemp, geo, axis = 0))
                            huuYlm[l, m] += huuTemp
    
    huuYl = np.zeros((lmax + 1, samples, 2))
    for l in range(0, lmax + 1):
        huuYl[l] += 2.*np.real(huuYlm[l, 0])
        for m in range(1, l + 1):
            huuYl[l] += 2.*np.real(huuYlm[l, m])
            huuYl[l] += 2.*np.real(np.roll(huuYlm[l, m], int(samples/2)))

    return huuYl, huuYlm

from scipy.special import ellipk
from pybhpt.geo import kerrgeo_Vt_radial, kerrgeo_Vt_polar, kerrgeo_Vr, kerrgeo_Vphi_radial, kerrgeo_Vphi_polar

def huu_reg_ecc(geo):
    L = geo.orbitalangularmomentum
    r0 = np.concatenate((geo.radialpoints, np.flip(geo.radialpoints)[1:-1]))
    q = geo.blackholespin
    kappa = L**2 + q**2 + 2*q**2/r0
    k = kappa/(r0**2 + kappa)
    K = ellipk(k)
    return 4.*K/(np.pi*np.sqrt(kappa + r0**2))

def huu_reg_inc(geo):
    En = geo.orbitalenergy
    Lz = geo.orbitalangularmomentum
    Qc = geo.carterconstant
    r = geo.radialpoints[0]
    th0 = np.concatenate((geo.polarpoints, np.flip(geo.polarpoints)[1:-1]))
    q = geo.blackholespin
    return huu_reg_gen(q, En, Lz, Qc, r, th0)

def huu_reg_gen(q, En, Lz, Qc, r, th):
    sig = r**2 + q**2*np.cos(th)**2
    eta = pow(sig,-1)*pow(2.*(2.*r + sig)*pow(q,2)*(sig*(-1.*Qc + pow(Lz,2)) + (2.*r + sig*pow(En,2))*(-1.*sig + pow(r,2))) + pow(q,4)*pow(2.*r + sig,2) + pow(sig*(Qc + pow(Lz,2)) - (2.*r + sig*pow(En,2))*(-1.*sig + pow(r,2)),2),0.5)
    zeta = Qc - 1.*r*(2. + r*(-2. + pow(En,2))) + sig*pow(En,2) + pow(Lz,2) + pow(q,2) + 2.*r*(pow(q,2) + pow(r,2))*pow(sig,-1)
    k = 2*eta/(eta + zeta)
    K = ellipk(k)
    return 4.*K/(np.pi*np.sqrt(eta/k))

def ut_ecc(q, En, Lz, r):
    Sigma = r**2
    theta = 0.5*np.pi
    Qc = 0.
    return (np.array([kerrgeo_Vt_radial(q, En, Lz, Qc, ri) for ri in r]) + kerrgeo_Vt_polar(q, En, Lz, Qc, theta))/Sigma

def uphi_ecc(q, En, Lz, r):
    Sigma = r**2
    theta = 0.5*np.pi
    Qc = 0.
    return (np.array([kerrgeo_Vphi_radial(q, En, Lz, Qc, ri) for ri in r]) + kerrgeo_Vphi_polar(q, En, Lz, Qc, theta))/Sigma

def ur_ecc(q, En, Lz, r):
    r0 = np.unique(r)
    Sigma = r**2
    Qc = 0.
    urUp = np.sqrt(np.abs(np.array([kerrgeo_Vr(q, En, Lz, Qc, ri) for ri in r0])))
    urUp[0] = 0.
    urUp[-1] = 0.
    urDown = -np.flip(urUp)[1:-1]
    ur = np.concatenate((urUp, urDown))/Sigma
    return ur

def ut_ecc_geo(geo):
    dM = geo.orbitalenergy
    dJ = geo.orbitalangularmomentum
    r = np.concatenate((geo.radialpoints, np.flip(geo.radialpoints)[1:-1]))
    q = geo.blackholespin
    return ut_ecc(q, dM, dJ, r)

def ur_ecc_geo(geo):
    dM = geo.orbitalenergy
    dJ = geo.orbitalangularmomentum
    r = np.concatenate((geo.radialpoints, np.flip(geo.radialpoints)[1:-1]))
    q = geo.blackholespin
    return ur_ecc(q, dM, dJ, r)

def uphi_ecc_geo(geo):
    dM = geo.orbitalenergy
    dJ = geo.orbitalangularmomentum
    r = np.concatenate((geo.radialpoints, np.flip(geo.radialpoints)[1:-1]))
    q = geo.blackholespin
    return uphi_ecc(q, dM, dJ, r)

def completion_ecc(geo):
    dM = geo.orbitalenergy
    dJ = geo.orbitalangularmomentum
    r = np.concatenate((geo.radialpoints, np.flip(geo.radialpoints)[1:-1]))
    q = geo.blackholespin
    htt = 2*dM/r
    hrr = 2.*r**2*((r + q**2)*dM - q*dJ)/(r**2 - 2.*r + q**2)**2
    htphi = -2*dJ/r
    hphiphi = 2.*q*((r + 2.)*dJ - (r + 1)*q*dM)/r
    ut = ut_ecc(q, dM, dJ, r)
    uphi = uphi_ecc(q, dM, dJ, r)
    ur = ur_ecc(q, dM, dJ, r)
    # return np.mean(r**2*(htt*ut*ut + 2*htphi*ut*uphi + hrr*ur*ur + hphiphi*uphi*uphi))/np.mean(r**2)
    return (htt*ut*ut + 2*htphi*ut*uphi + hrr*ur*ur + hphiphi*uphi*uphi)

def orbit_average_ecc(func, geo, axis = 0):
    r = np.concatenate((geo.radialpoints, np.flip(geo.radialpoints)[1:-1]))
    return np.mean(func*r**2, axis = axis)/np.mean(r**2)

def orbit_average_inc(func, geo, axis = 0):
    th = np.concatenate((geo.polarpoints, np.flip(geo.polarpoints)[1:-1])).T
    r = geo.radialpoints[0]
    q = geo.blackholespin
    sig = r**2 + q**2*np.cos(th)**2
    if len(func.shape) <= 1:
        sigA = np.array(sig)
    else:
        sigA = np.array([sig, sig]).T
    return np.mean(func*sigA, axis = axis)/np.mean(sig)

def Tr_ecc_geo_ratio(geo):
    dM = geo.orbitalenergy
    dJ = geo.orbitalangularmomentum
    r = np.concatenate((geo.radialpoints, np.flip(geo.radialpoints)[1:-1]))
    q = geo.blackholespin
    return np.mean(ut_ecc(q, dM, dJ, r)*r**2)/np.mean(r**2)

nmin = -20
nmax = 20
a = 0.9
p = 12.405076555398
e = 0.1
samples = 2**8
geoEcc = KerrGeodesic(a, p, e, 1.0, samples)
gauges = ["IRG", "ORG", "SAAB"]

for gauge in gauges:
    huuData0 = redshift_calc_ecc(gauge, 30, geoEcc, nrange = [nmin, nmax])
    fn = "redshift_a{}_p{}_e{}_".format(int(100*a), int(p), int(100*e)+gauge+".npy")
    np.save(fn, huuData0[0])