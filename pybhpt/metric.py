from cybhpt_full import metric_coefficients_cython_S4dagger, metric_coefficients_cython_S0dagger
import numpy as np
from pybhpt.hertz import available_gauges
from pybhpt.swsh import Yslm, spin_operator_normalization

def gauge_check(gauge):
    if gauge not in available_gauges:
        TypeError("{} is not a supported gauge.".format(gauge))

S0_gauges = ["IRG", "ARG0", "SRG0"]
S4_gauges = ["ORG", "ARG4", "SRG4"]

def metric_coefficients_S4dagger_ab(ai, bi, nt, nr, nz, nph, a, r, z):
    """Compute the metric coefficients for the reconstructed perturbation associated with the
    S4dagger reconstruction operator.
    
    Parameters
    ----------
    ai, bi : int
        The indices of the tetrad-projected components of the metric perturbation.
    nt, nr, nz, nph : int
        The indices associated with the time, radial, spin, and azimuthal derivative operators.
    a : float
        The black hole spin parameter.
    r : float
        The radial coordinate.
    z : float
        The polar coordinate (cosine of the polar angle).
    
    Returns
    -------
    float
        The computed metric coefficient."""
    return metric_coefficients_cython_S4dagger(ai, bi, nt, nr, nz, nph, a, r, z)

def metric_coefficients_S0dagger_ab(ai, bi, nt, nr, nz, nph, a, r, z):
    """Compute the metric coefficients for the reconstructed perturbation associated with the
    S0dagger reconstruction operator.
    
    Parameters
    ----------
    ai, bi : int
        The indices of the tetrad-projected components of the metric perturbation.
    nt, nr, nz, nph : int
        The indices associated with the time, radial, spin, and azimuthal derivative operators.
    a : float
        The black hole spin parameter.
    r : float
        The radial coordinate.
    z : float
        The polar coordinate (cosine of the polar angle).
    
    Returns
    -------
    float
        The computed metric coefficient."""
    return metric_coefficients_cython_S0dagger(ai, bi, nt, nr, nz, nph, a, r, z)

def metric_coefficients_S0dagger(a, b, c, d, q, rvals, zvals):
    """Compute the metric coefficients for the reconstructed perturbation associated with the
    S0dagger reconstruction operator.
    
    Parameters
    ----------
    a, b, c, d : int
        The indices associated with the time, radial, spin, and azimuthal derivative operators.
    q : float
        The black hole spin parameter.
    rvals : numpy.ndarray
        The radial coordinate.
    zvals : np.ndarray
        The polar coordinate (cosine of the polar angle).
    
    Returns
    -------
    numpy.ndarray
        A 3D array containing the computed metric coefficients for the S0dagger operator."""
    h22 = np.array([[metric_coefficients_S0dagger_ab(2, 2, a, b, c, d, q, r, z) for z in zvals] for r in rvals])
    h24 = np.array([[metric_coefficients_S0dagger_ab(2, 4, a, b, c, d, q, r, z) for z in zvals] for r in rvals])
    h44 = np.array([[metric_coefficients_S0dagger_ab(4, 4, a, b, c, d, q, r, z) for z in zvals] for r in rvals])
    return np.array([2.*h22, h24, h44])

def metric_coefficients_S4dagger(a, b, c, d, q, rvals, zvals):
    """Compute the metric coefficients for the reconstructed perturbation associated with the
    S4dagger reconstruction operator.
    
    Parameters
    ----------
    a, b, c, d : int
        The indices associated with the time, radial, spin, and azimuthal derivative operators.
    q : float
        The black hole spin parameter.
    rvals : numpy.ndarray
        The radial coordinate.
    zvals : np.ndarray
        The polar coordinate (cosine of the polar angle).
    
    Returns
    -------
    numpy.ndarray
        A 3D array containing the computed metric coefficients for the S0dagger operator."""
    h11 = np.array([[metric_coefficients_S4dagger_ab(1, 1, a, b, c, d, q, r, z) for z in zvals] for r in rvals])
    h13 = np.array([[metric_coefficients_S4dagger_ab(1, 3, a, b, c, d, q, r, z) for z in zvals] for r in rvals])
    h33 = np.array([[metric_coefficients_S4dagger_ab(3, 3, a, b, c, d, q, r, z) for z in zvals] for r in rvals])
    return np.array([2.*h11, h13, h33])

def metric_coefficients(gauge, a, b, c, d, q, rvals, zvals):
    """Compute the metric coefficients for the reconstructed perturbation associated with the
    specified gauge.

    Parameters
    ----------
    gauge : str
        The gauge to use for the reconstruction. Must be one of the available gauges.
    a, b, c, d : int
        The indices associated with the time, radial, spin, and azimuthal derivative operators.
    q : float
        The black hole spin parameter.
    rvals : numpy.ndarray
        The radial coordinate.
    zvals : np.ndarray
        The polar coordinate (cosine of the polar angle).

    Returns
    -------
    numpy.ndarray
        A 3D array containing the computed metric coefficients for the specified gauge.
    """
    gauge_check(gauge)
    if gauge in S0_gauges:
        return metric_coefficients_S0dagger(a, b, c, d, q, rvals, zvals)
    else:
        return metric_coefficients_S4dagger(a, b, c, d, q, rvals, zvals)

class MetricCoefficients:
    """
    A class for computing the metric coefficients of the reconstructed perturbation
    associated with the S0dagger or S4dagger reconstruction operator based on the specified gauge.

    Parameters
    ----------
    gauge : str
        The gauge to use for the reconstruction. Must be one of the available gauges.
    q : float
        The black hole spin parameter.
    r : numpy.ndarray
        The radial coordinate values.
    th : numpy.ndarray
        The polar coordinate values (cosine of the polar angle).

    Attributes
    ----------
    gauge : str
        The gauge used for the reconstruction.
    blackholespin : float
        The black hole spin parameter.
    radialpoints : numpy.ndarray
        The radial coordinate values.
    polarpoints : numpy.ndarray
        The polar coordinate values (cosine of the polar angle).
    storedcomponents : dict
        A dictionary mapping pairs of indices (a, b) to the index of the stored component
        in the coefficients array.
    conjugatecomponents : dict
        A dictionary mapping pairs of indices (a, b) to the index of the conjugate component
        in the coefficients array.
    coeffs : numpy.ndarray
        A 7D array containing the computed metric coefficients for the specified gauge.
    zeros : numpy.ndarray
        A 2D array of zeros with the same shape as the radial and polar coordinate arrays.

    Methods
    -------
    hab(a, b, nt, nr, ns, nphi):
        Returns the metric coefficient for the specified indices (a, b) and derivative orders
        (nt, nr, ns, nphi). If the indices are not found in the stored or conjugate components,
        it returns a zero array.
    __call__(a, b, nt, nr, ns, nphi):
        Calls the `hab` method to retrieve the metric coefficient for the specified indices
        and derivative orders.
    """
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
        """
        Returns the metric coefficient for the specified indices (a, b) and derivative orders
        (nt, nr, ns, nphi). If the indices are not found in the stored or conjugate components,
        it returns a zero array.

        Parameters
        ----------
        a, b : int
            The indices of the tetrad-projected components of the metric perturbation.
        nt, nr, ns, nphi : int
            The indices associated with the time, radial, spin, and azimuthal derivative operators.

        Returns
        -------
        numpy.ndarray
            The computed metric coefficient for the specified indices and derivative orders.
        """
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
        """ Calls the `hab` method to retrieve the metric coefficient for the specified indices
        and derivative orders.

        Parameters
        ----------
        a, b : int
            The indices of the tetrad-projected components of the metric perturbation.
        nt, nr, ns, nphi : int
            The indices associated with the time, radial, spin, and azimuthal derivative operators.

        Returns
        -------
        numpy.ndarray
            The computed metric coefficient for the specified indices and derivative orders.
        """
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
