import scipy.sparse
import scipy.sparse.linalg
from scipy.special import sph_harm
from scipy.special import binom
from scipy.special import factorial
import numpy as np
from cybhpt_full import YslmCy, clebschCy, w3jCy

"""
Wigner 3j-symbol and Clebsch-Gordon coefficients
"""

def fac(n):
    """
    Computes the factorial of a non-negative integer n.

    Parameters
    ----------
    n : int
        A non-negative integer.

    Returns
    -------
    float
        The factorial of n.
    """
    if n < 0:
        return 0
    return float(np.math.factorial(n))

def Yslm(s, l, m, th, ph = None):
    """
    Evaluate the spin-weighted spherical harmonic $Y_{s}^{lm}$ at a given angle theta.

    Parameters
    ----------
    s : int
        The spin weight of the harmonic.
    l : int
        The angular number of the spherical harmonic.
    m : int
        The azimuthal number of the spherical harmonic.
    th : array_like
        The polar angle(s) at which to evaluate the spherical harmonic.

    Returns
    -------
    array_like
        The values of the spherical harmonic at the specified angles.
    """
    assert isinstance(s, int) and isinstance(l, int) and isinstance(m, int), "s, l, and m must be integers"
    if ph is not None:
        return Yslm(s, l, m, th)*np.exp(1.j*m*ph)
    if np.abs(s) > l:
        return 0.*th
    if s == 0:
        return np.real(sph_harm(m, l, 0., th))
    elif s + m < 0:
        return (-1.)**(s+m)*YslmBase(-s, l, -m, th)
    else:
        return YslmBase(s, l, m, th)

def YslmBase(s, l, m, th):
    assert isinstance(s, int) and isinstance(l, int) and isinstance(m, int), "s, l, and m must be integers"
    if not isinstance(th, (int, float)):
        b = np.broadcast(th)
        out = np.empty(b.shape)
        out.flat = [YslmCy(s, l, m, thi) for (thi,) in b]
    else:
        out = YslmCy(s, l, m, th)
    return out

def Yslm_eigenvalue(s, l, *args):
    return l*(l + 1.) - s*(s + 1.)

def clebsch(l1, l2, l3, m1, m2, m3):
    """
    Compute the Clebsch-Gordon coefficient <l1,m1,l2,m2|l3,m3>.

    Parameters
    ----------
    l1 : int
        The angular number of the first state.
    l2 : int
        The angular number of the second state.
    l3 : int
        The angular number of the combined state.
    m1 : int
        The azimuthal number of the first state.
    m2 : int
        The azimuthal number of the second state.
    m3 : int
        The azimuthal number of the combined state.

    Returns
    -------
    float
        The Clebsch-Gordon coefficient <l1,m1,l2,m2|l3,m3>.
    """
    assert isinstance(l1, int) and isinstance(l2, int) and isinstance(l3, int), "l1, l2, and l3 must be integers"
    assert isinstance(m1, int) and isinstance(m2, int) and isinstance(m3, int), "m1, m2, and m3 must be integers"
    return clebschCy(l1, l2, l3, m1, m2, m3)

def w3j(l1, l2, l3, m1, m2, m3):
    """
    Compute the Wigner 3j-symbol
        | l1  l2  l3 |
        | m1  m2  m3 |

    Parameters
    ----------
    l1 : int
        The angular number of the first state. 
    l2 : int
        The angular number of the second state.
    l3 : int
        The angular number of the combined state.
    m1 : int
        The azimuthal number of the first state.
    m2 : int
        The azimuthal number of the second state.
    m3 : int
        The azimuthal number of the combined state.

    Returns
    -------
    float
        The Wigner 3j-symbol $ \begin{pmatrix} l1 & l2 & l3 \\ m1 & m2 & m3 \end{pmatrix} $
    """
    assert isinstance(l1, int) and isinstance(l2, int) and isinstance(l3, int), "l1, l2, and l3 must be integers"
    assert isinstance(m1, int) and isinstance(m2, int) and isinstance(m3, int), "m1, m2, and m3 must be integers"
    return w3jCy(l1, l2, l3, m1, m2, m3)
    
"""
SWSH Eigenvalue Functions
"""

def k1(s, l, j, m):
    return np.sqrt((2*l + 1)/(2*j + 1))*clebsch(l, 1, j, m, 0, m)*clebsch(l, 1, j, -s, 0, -s);

def k2(s, l, j, m):
    ktemp = 2./3.*np.sqrt((2*l + 1)/(2*j + 1))*clebsch(l, 2, j, m, 0, m)*clebsch(l, 2, j, -s, 0, -s);
    if j == l:
        ktemp += 1/3.
    return ktemp

def k2m2(s, l, m):
    temp = (l - m - 1.)/(l - 1.)
    temp *= (l + m - 1.)/(l - 1.)
    temp *= np.float64(l - m)/l
    temp *= np.float64(l + m)/l
    temp *= (l - s)/(2.*l - 1.)
    temp *= (l + s)/(2.*l - 1.)
    temp *= (l - s - 1.)/(2.*l + 1.)
    temp *= (l + s - 1.)/(2.*l - 3.)
    return np.sqrt(temp)

def k2m1(s, l, m):
    temp = np.float64(l - m)*np.float64(l + m)/(2.*l - 1.)
    temp *= np.float64(l - s)*np.float64(l + s)/(2.*l + 1.)
    return -2.*m*s*np.sqrt(temp)/l/(l - 1.)/(l + 1.)

def k2p0(s, l, m):
    temp = np.float64(l*(l + 1.) - 3.*m*m)/(2.*l - 1.)/l
    temp *= np.float64(l*(l + 1.) - 3.*s*s)/(2.*l + 3.)/(l + 1.)
    return 1./3.*(1. + 2.*temp)

def k2p1(s, l, m):
    temp = np.float64(l - m + 1.)*np.float64(l + m + 1.)/(2.*l + 1.)
    temp *= np.float64(l - s + 1.)*np.float64(l + s + 1.)/(2.*l + 3.)
    return -2.*m*s*np.sqrt(temp)/l/(l + 1.)/(l + 2.)

def k2p2(s, l, m):
    temp = (l - m + 1.)/(l + 1.)
    temp *= (l + m + 1.)/(l + 1.)
    temp *= (l - m + 2.)/(l + 2.)
    temp *= (l + m + 2.)/(l + 2.)
    temp *= (l - s + 2.)/(2.*l + 3.)
    temp *= (l + s + 2.)/(2.*l + 3.)
    temp *= (l - s + 1.)/(2.*l + 1.)
    temp *= (l + s + 1.)/(2.*l + 5.)
    return np.sqrt(temp)

def k1m1(s, l, m):
    temp = np.float64(l - m)*np.float64(l + m)/(2.*l - 1.)
    temp *= np.float64(l - s)*np.float64(l + s)/(2.*l + 1.)
    return np.sqrt(temp)/l

def k1p0(s, l, m):
    return -np.float64(m*s)/l/(l + 1.)

def k1p1(s, l, m):
    temp = np.float64(l - m + 1.)*np.float64(l + m + 1.)/(2.*l + 3.)
    temp *= np.float64(l - s + 1.)*np.float64(l + s + 1.)/(2.*l + 1.)
    return np.sqrt(temp)/(l + 1.)

def akm2(s, l, m, g):
    return -g*g*k2m2(s, l, m)

def akm1(s, l, m, g):
    return -g*g*k2m1(s, l, m) + 2.*s*g*k1m1(s, l, m)

def akp0(s, l, m, g):
    return -g*g*k2p0(s, l, m) + 2.*s*g*k1p0(s, l, m) + l*(l + 1.) - s*(s + 1.) - 2.*m*g + g*g

def akp1(s, l, m, g):
    return -g*g*k2p1(s, l, m) + 2.*s*g*k1p1(s, l, m)

def akp2(s, l, m, g):
    return -g*g*k2p2(s, l, m)

def spectral_sparse_matrix(s, m, g, nmax):
    lmin = max(abs(s), abs(m))
    larray = np.arange(lmin, lmin + nmax)
    return scipy.sparse.diags([akm2(s, larray[2:], m, g), akm1(s, larray[1:], m, g), akp0(s, larray, m, g), akp1(s, larray[:-1], m, g), akp2(s, larray[:-2], m, g)], [-2, -1, 0, 1, 2])

def swsh_eigs(s, l, m, g, nmax=None, return_eigenvectors=True):
    lmin = max(abs(s), abs(m))
    kval = l - lmin
    
    if nmax is None:
        buffer = round(20 + 2*np.abs(g))
        Nmax = kval + buffer + 2
    else:
        if nmax < kval:
            Nmax = kval + 5
        else:
            Nmax = nmax
    
    mat = spectral_sparse_matrix(s, m, g, Nmax)
    out = scipy.sparse.linalg.eigs(mat, k=Nmax-2, which='SM', return_eigenvectors=return_eigenvectors)
    
    return out

def swsh_coeffs(s, l, m, g, th):
    if g == 0.:
        return Yslm(s, l, m, th)
    
    _, eig = swsh_eigs(s, l, m, g, nmax=None, return_eigenvectors=True)
    if g.imag == 0.:
        coeffs = np.real(eig[l - max(abs(s), abs(m))])
    else:
        coeffs = eig[l - max(abs(s), abs(m))]
    return coeffs

def swsh_eigenvalue(s, l, m, g, nmax=None):
    """
    Compute the eigenvalue of the spin-weighted spheroidal harmonic.
    
    Parameters
    ----------
    s : int
        The spin weight of the harmonic.
    l : int
        The angular number of the harmonic.
    m : int
        The azimuthal number of the harmonic.
    g : float or complex
        The spheroidicity parameter.
    nmax : int, optional
        The maximum number of basis functions to use in the computation. If None, a default value is chosen.

    Returns
    -------
    float or complex
        The eigenvalue of the spin-weighted spheroidal harmonic.
    """
    if g == 0.:
        return Yslm_eigenvalue(s, l)
    
    las = swsh_eigs(s, l, m, g, nmax=nmax, return_eigenvectors=False)
    
    idx = l - max(abs(s), abs(m))
    sorted_indices = np.argsort(np.real(las))
    if idx < 0 or idx >= len(sorted_indices):
        raise IndexError(f"Index {idx} out of bounds for eigenvalue array of length {len(sorted_indices)}")
    if g.imag == 0.:
        eigen = np.real(las)[sorted_indices[idx]]
    else:
        eigen = las[sorted_indices[idx]]
    return eigen
class SWSHBase:
    def __init__(self, *args):
        arg_num = np.array(args).shape[0]
        if arg_num < 3:
            print('Error. Not enough arguments to create class')
            pass
        
        self.s = args[0]
        self.l = args[1]
        self.m = args[2]
        self.lmin = max(abs(self.s), abs(self.m))

        if arg_num > 3:
            self.spheroidicity = args[3]
        if self.spheroidicity.imag == 0:
            self.spheroidicity = np.real(self.spheroidicity)

class SWSHSeriesBase(SWSHBase):
    def __init__(self, s, l, m, g):
        SWSHBase.__init__(self, s, l, m, g)
        
    def sparse_matrix(self, nmax):
        return spectral_sparse_matrix(self.s, self.m, self.spheroidicity, nmax)
        
    def eigs(self, nmax = None, **kwargs):
        kval = self.l - self.lmin
    
        if nmax is None:
            buffer = round(20 + np.abs(2*self.spheroidicity))
            Nmax = kval + buffer + 2
        else:
            if nmax < kval:
                Nmax = kval + 5
            else:
                Nmax = nmax
                
        if "k" not in kwargs.keys():
            kwargs["k"] = Nmax - 2
                
        if "which" not in kwargs.keys():
            kwargs["which"] = 'SM'
            
#         if "sigma" not in kwargs.keys():
#             kwargs["sigma"] = (self.l + Nmax)*(self.l + Nmax + 1) - self.s*(self.s + 1)
            
        if "return_eigenvectors" not in kwargs.keys():
            kwargs["return_eigenvectors"] = True
            
        mat = self.sparse_matrix(Nmax)
        return scipy.sparse.linalg.eigs(mat, **kwargs)
    
    def generate_eigenvalue(self):
        if self.spheroidicity.imag == 0.:
            las = np.real(self.eigs(return_eigenvectors=False))
        else:
            las = self.eigs(return_eigenvectors=False)
        pos = np.argsort(np.real(las))[self.l - self.lmin]
        return las[pos]

    def generate_eigs(self):
        las, eigs = self.eigs()
        pos_vec = np.argsort(np.real(las))
        pos = pos_vec[self.l - self.lmin]
        if self.spheroidicity.imag == 0.:
            eigs_temp = np.real(eigs[:, pos])
            eigs_return = np.sign(eigs_temp[self.l - self.lmin])*eigs_temp
            eig = np.real(las[pos])
        else:
            eigs_temp = eigs[:, pos]
            ref = eigs_temp[self.l - self.lmin]
            eigs_temp = eigs_temp/ref
            eigs_norm = np.linalg.norm(eigs_temp)
            eigs_return = eigs_temp/eigs_norm
            eig = las[pos]
        return (eig, eigs_return)
class SpinWeightedSpheroidalHarmonic(SWSHSeriesBase):
    """
    A class for generating a spin-weighted spheroidal harmonic.

    Parameters
    ----------
    s : int
        The spin weight of the harmonic.
    l : int
        The angular number of the harmonic.
    m : int
        The azimuthal number of the harmonic.
    g : float or complex
        The spheroidicity parameter.

    Attributes
    ----------
    couplingcoefficients : array_like
        The coupling coefficients between the spin-weighted spheroidal harmonic and the spin-weighted spherical  
    
    """
    def __init__(self, s, l, m, g):
        SWSHSeriesBase.__init__(self, s, l, m, g)
        if self.spheroidicity == 0.:
            self.eval = self.Yslm
            self.eigenvalue = Yslm_eigenvalue(self.s, self.l)
            self.coeffs = np.zeros(self.l - self.lmin + 1)
            self.coeffs[-1] = 1.
        else:
            self.eval = self.Sslm
            self.eigenvalue, self.coeffs = self.generate_eigs()

    @property
    def couplingcoefficients(self):
        return self.coeffs

    def Yslm(self, l, th):
        """
        Evaluate the spin-weighted spherical harmonic $Y_{s}^{lm}(theta)$ at a given angle theta.

        Parameters
        ----------
        l : int
            The angular number of the spherical harmonic.
        th : array_like
            The polar angle(s) at which to evaluate the spherical harmonic.

        Returns
        -------
        array_like
            The values of the spherical harmonic at the specified angles.
        """
        return Yslm(self.s, l, self.m, th)
    
    def Sslm(self, *args):
        """
        Evaluate the spin-weighted spheroidal harmonic $S_{s}^{lm}(theta)$ at a given angle theta.

        Parameters
        ----------
        th : array_like
            The polar angle(s) at which to evaluate the spheroidal harmonic.
        
        Returns
        -------
        array_like
            The values of the spheroidal harmonic at the specified angles.
        """
        th = args[-1]
        term_num = self.coeffs.shape[0]
        if isinstance(th, (int, float)):
            Yslm_array = np.empty(term_num)
        else:
            pts_num = th.shape[0]
            Yslm_array = np.empty((term_num, pts_num))
        for i in range(term_num):
            Yslm_array[i] = self.Yslm(self.lmin + i, th)
            
        return np.dot(self.coeffs, Yslm_array)
            
    def __call__(self, th, ph = None):
        out = self.eval(self.l, th)
        if ph is not None:
            out = out*np.exp(1.j*self.m*ph)
        return out

def muCoupling(s, l):
    """
    Eigenvalue for the spin-weighted spherical harmonic lowering operator
    Setting s -> -s gives the negative of the eigenvalue for the raising operator
    """
    if l + s < 0 or l - s + 1. < 0:
        return 0
    return np.sqrt((l - s + 1.)*(l + s))

def Asjlm(s, j, l, m):
    """
    Coupling coefficient between scalar and spin-weighted spherical harmonics
    
    Parameters
    ----------
    s : int
        The spin weight of the harmonic.
    j : int
        The angular number of the scalar harmonic.
    l : int
        The angular number of the spin-weighted harmonic.
    m : int
        The azimuthal number of the harmonics. 

    Returns
    -------
    float
        The coupling coefficient $A_{s}^{jlm}$
    """
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
        