# from geo_wrap cimport GeodesicSource, kerr_geo_orbit, kerr_geo_orbital_constants, kerr_geo_mino_frequencies

from libcpp.vector cimport vector
import numpy as np
cimport numpy as np

from libcpp.vector cimport vector

cdef extern from "geo.hpp":
    cdef cppclass GeodesicTrajectory:
        GeodesicTrajectory()
        GeodesicTrajectory(vector[double] tR, vector[double] tTheta, vector[double] r, vector[double] theta, vector[double] phiR, vector[double] phiTheta)
        vector[double] tR, tTheta, r, theta, phiR, phiTheta

    cdef cppclass GeodesicConstants:
        GeodesicConstants()
        GeodesicConstants(double a, double p, double e, double x, double En, double Lz, double Q)
        GeodesicConstants(double a, double p, double e, double x, double En, double Lz, double Q, double, double, double, double)
        GeodesicConstants(double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double)
        double a, p, e, x
        double En, Lz, Q
        double r1, r2, r3, r4, z1, z2
        double upsilonT, upsilonR, upsilonTheta, upsilonPhi
        double carterR, carterTheta, carterPhi

    cpdef cppclass GeodesicSource:
        # GeodesicSource() except +
        GeodesicSource(double a, double p, double e, double x, int Nsample) except +
        # ~GeodesicSource()

        int getOrbitalSampleNumber()
        double getBlackHoleSpin()
        double getSemiLatusRectum()
        double getEccentricity()
        double getInclination()
        double getOrbitalEnergy()
        double getOrbitalAngularMomentum()
        double getCarterConstant()
        double getRadialRoot(int i)
        double getPolarRoot(int i)
        double getMinoFrequency(int mu)
        double getTimeFrequency(int i)
        double getTimeFrequency(int m, int k, int n)
        double getCarterFrequency(int i)
        double getCarterFrequency(int m, int k, int n)

        vector[double] getTimeAccumulation(int j)
        vector[double] getRadialPosition()
        vector[double] getPolarPosition()
        vector[double] getAzimuthalAccumulation(int j)

        double getTimeAccumulation(int j, int pos)
        double getRadialPosition(int pos)
        double getPolarPosition(int pos)
        double getAzimuthalAccumulation(int j, int pos)

        double getPsiRadialOfMinoTime(double la)
        double getPsiPolarOfMinoTime(double la)

        double getTimePositionOfMinoTime(double la)
        double getRadialPositionOfMinoTime(double la)
        double getPolarPositionOfMinoTime(double la)
        double getAzimuthalPositionOfMinoTime(double la)

        vector[double] getPsiRadialOfMinoTime(vector[double] la)
        vector[double] getPsiPolarOfMinoTime(vector[double] la)

        vector[double] getTimePositionOfMinoTime(vector[double] la)
        vector[double] getRadialPositionOfMinoTime(vector[double] la)
        vector[double] getPolarPositionOfMinoTime(vector[double] la)
        vector[double] getAzimuthalPositionOfMinoTime(vector[double] la)
        vector[double] getPositionOfMinoTime(vector[double] la)

        double getMinoTimeOfTime(double t)

        vector[double] getTimeCoefficients(int j)
        vector[double] getRadialCoefficients()
        vector[double] getPolarCoefficients()
        vector[double] getAzimuthalCoefficients(int j)

        GeodesicConstants getConstants()
        GeodesicTrajectory getTrajectory()
        GeodesicTrajectory getCoefficients()

    GeodesicSource kerr_geo_orbit(double a, double p, double e, double x, int n)
    void kerr_geo_kepler_parameters(double &p, double &e, double &x, double &a, double &En, double &Lz, double &Qc)
    void kerr_geo_kepler_parameters(int n, double* p, double* e, double* x, const double* a, const double* En, const double* Lz, const double* Qc)
    void kerr_geo_orbital_constants(double &En, double &Lz, double &Qc, double &a, double &p, double &e, double &x)
    void kerr_geo_orbital_constants(int n, double* En, double* Lz, double* Qc, const double* a, const double* p, const double* e, const double* x)

    void kerr_geo_radial_roots(double &r1, double &r2, double &r3, double &r4, double &a, double &p, double &e, double &En, double &Lz, double &Qc)
    void kerr_geo_polar_roots(double &z1, double &z2, double &a, double &x, double &En, double &Lz, double &Qc)
    void kerr_geo_mino_frequencies(double &upT, double &upR, double &upTh, double &upPh, double &a, double &p, double &e, double &x)

    void jacobian_ELQ_to_pex(double &dpdE, double &dedE, double &dxdE, double &dpdLz, double &dedLz, double &dxdLz, double &dpdQ, double &dedQ, double &dxdQ, double a, double p, double e, double x)
    void jacobian_ELQ_to_pex(int n, double* dpdE, double* dedE, double* dxdE, double* dpdLz, double* dedLz, double* dxdLz, double* dpdQ, double* dedQ, double* dxdQ, const double* a, const double* p, const double* e, const double* x)
    void jacobian_pex_to_ELQ(double &dEdp, double &dEde, double &dEdx, double &dLdp, double &dLde, double &dLdx, double &dQdp, double &dQde, double &dQdx, double a, double p, double e, double x)
    void jacobian_pex_to_ELQ(int n, double* dEdp, double* dEde, double* dEdx, double* dLdp, double* dLde, double* dLdx, double* dQdp, double* dQde, double* dQdx, const double* a, const double* p, const double* e, const double* x)
    void jacobian_ELQ_to_pex_spherical(double &dpdE, double &dedE, double &dxdE, double &dpdLz, double &dedLz, double &dxdLz, double &dpdQ, double &dedQ, double &dxdQ, double a, double p, double e, double x)
    void jacobian_ELQ_to_pex_spherical(int n, double* dpdE, double* dedE, double* dxdE, double* dpdLz, double* dedLz, double* dxdLz, double* dpdQ, double* dedQ, double* dxdQ, const double* a, const double* p, const double* e, const double* x)
    void jacobian_pex_to_ELQ_spherical(double &dEdp, double &dEde, double &dEdx, double &dLdp, double &dLde, double &dLdx, double &dQdp, double &dQde, double &dQdx, double a, double p, double e, double x)
    void jacobian_pex_to_ELQ_spherical(int n, double* dEdp, double* dEde, double* dEdx, double* dLdp, double* dLde, double* dLdx, double* dQdp, double* dQde, double* dQdx, const double* a, const double* p, const double* e, const double* x)

    double kerr_geo_VtR(const double & a, const double & En, const double & Lz, const double & Q, const double & r)
    double kerr_geo_VtTheta(const double & a, const double & En, const double & Lz, const double & Q, const double & theta)
    double kerr_geo_Vr(const double & a, const double & En, const double & Lz, const double & Q, const double & r)
    double kerr_geo_Vtheta(const double & a, const double & En, const double & Lz, const double & Q, const double & theta)
    double kerr_geo_VphiR(const double & a, const double & En, const double & Lz, const double & Q, const double & r)
    double kerr_geo_VphiTheta(const double & a, const double & En, const double & Lz, const double & Q, const double & theta)
    double kerr_geo_Vz(const double & a, const double & En, const double & Lz, const double & Q, const double & z)
    double kerr_geo_Vz_dz(const double & a, const double & En, const double & Lz, const double & Q, const double & z)
    double kerr_geo_Vz_dz2(const double & a, const double & En, const double & Lz, const double & Q, const double & z)
    double kerr_isco(double a, int sgnX)
    double kerr_isco_frequency(double a)

def _kerr_geo_V01(double a, double En, double Lz, double Q, double r):
    return kerr_geo_VtR(a, En, Lz, Q, r)

def _kerr_geo_V02(double a, double En, double Lz, double Q, double theta):
    return kerr_geo_VtTheta(a, En, Lz, Q, theta)

def _kerr_geo_V11(double a, double En, double Lz, double Q, double r):
    return kerr_geo_Vr(a, En, Lz, Q, r)

def _kerr_geo_V22(double a, double En, double Lz, double Q, double theta):
    return kerr_geo_Vtheta(a, En, Lz, Q, theta)

def _kerr_geo_V31(double a, double En, double Lz, double Q, double r):
    return kerr_geo_VphiR(a, En, Lz, Q, r)

def _kerr_geo_V32(double a, double En, double Lz, double Q, double theta):
    return kerr_geo_VphiTheta(a, En, Lz, Q, theta)

cdef class KerrGeodesic:
    cdef GeodesicSource *geocpp
    cdef int nsamplescpp

    def __init__(self, double a, double p, double e, double x, int nsamples = 2**8):
        self.geocpp = new GeodesicSource(a, p, e, x, nsamples)
        self.nsamplescpp = nsamples

    def __dealloc__(self):
        del self.geocpp

    @property
    def nsamples(self):
        return self.nsamplescpp

    @property
    def blackholespin(self):
        return self.geocpp.getBlackHoleSpin()
    
    @property
    def semilatusrectum(self):
        return self.geocpp.getSemiLatusRectum()
    
    @property
    def eccentricity(self):
        return self.geocpp.getEccentricity()

    @property
    def inclination(self):
        return self.geocpp.getInclination()

    @property
    def orbitalenergy(self):
        return self.geocpp.getOrbitalEnergy()

    @property
    def orbitalangularmomentum(self):
        return self.geocpp.getOrbitalAngularMomentum()

    @property
    def carterconstant(self):
        return self.geocpp.getCarterConstant()

    @property
    def radialroots(self):
        return np.array([self.geocpp.getRadialRoot(1+i) for i in range(4)])

    @property
    def polarroots(self):
        return np.array([self.geocpp.getPolarRoot(1+i) for i in range(2)])

    @property
    def minofrequencies(self):
        return np.array([self.geocpp.getMinoFrequency(i) for i in range(4)])

    @property
    def timefrequencies(self):
        return np.array([self.geocpp.getTimeFrequency(i) for i in range(1, 4)])

    @property
    def frequencies(self):
        return self.timefrequencies

    @property
    def carterfrequencies(self):
        return np.array([self.geocpp.getCarterFrequency(i) for i in range(1, 4)]) 

    def mode_time_frequency(self, np.ndarray[ndim=1, dtype=np.int64_t] kvec):
        return np.dot(kvec, (self.frequencies))
    
    mode_frequency = mode_time_frequency

    def mode_carter_frequency(self, np.ndarray[ndim=1, dtype=np.int64_t] kvec):
        return np.dot(kvec,(self.carterfrequencies))

    cdef void getPsiRadialOfMinoTimeArray(self, np.float64_t *psi, np.float64_t *la, int n):
        for i in range(n):
            psi[i] = self.geocpp.getPsiRadialOfMinoTime(la[i])
    cdef void getPsiPolarOfMinoTimeArray(self, np.float64_t *psi, np.float64_t *la, int n):
        for i in range(n):
            psi[i] = self.geocpp.getPsiPolarOfMinoTime(la[i])

    cdef void getPsiRadialOfTimeArray(self, np.float64_t *psi, np.float64_t *t, int n):
        cdef double tmp
        for i in range(n):
            tmp = self.geocpp.getMinoTimeOfTime(t[i])
            psi[i] = self.geocpp.getPsiRadialOfMinoTime(tmp)
    cdef void getPsiPolarOfTimeArray(self, np.float64_t *psi, np.float64_t *t, int n):
        cdef double tmp
        for i in range(n):
            tmp = self.geocpp.getMinoTimeOfTime(t[i])
            psi[i] = self.geocpp.getPsiPolarOfMinoTime(tmp)

    cdef void getTimePositionOfMinoTimeArray(self, np.float64_t *t, np.float64_t *la, int n):
        for i in range(n):
            t[i] = self.geocpp.getTimePositionOfMinoTime(la[i])
    cdef void getRadialPositionOfMinoTimeArray(self, np.float64_t *t, np.float64_t *la, int n):
        for i in range(n):
            t[i] = self.geocpp.getRadialPositionOfMinoTime(la[i])
    cdef void getPolarPositionOfMinoTimeArray(self, np.float64_t *t, np.float64_t *la, int n):
        for i in range(n):
            t[i] = self.geocpp.getPolarPositionOfMinoTime(la[i])
    cdef void getAzimuthalPositionOfMinoTimeArray(self, np.float64_t *t, np.float64_t *la, int n):
        for i in range(n):
            t[i] = self.geocpp.getAzimuthalPositionOfMinoTime(la[i])

    cdef void getRadialPositionOfTimeArray(self, np.float64_t *r, np.float64_t *t, int n):
        cdef double tmp
        for i in range(n):
            tmp = self.geocpp.getMinoTimeOfTime(t[i])
            r[i] = self.geocpp.getRadialPositionOfMinoTime(tmp)
    cdef void getPolarPositionOfTimeArray(self, np.float64_t *th, np.float64_t *t, int n):
        cdef double tmp
        for i in range(n):
            tmp = self.geocpp.getMinoTimeOfTime(t[i])
            th[i] = self.geocpp.getPolarPositionOfMinoTime(tmp)
    cdef void getAzimuthalPositionOfTimeArray(self, np.float64_t *ph, np.float64_t *t, int n):
        cdef double tmp
        for i in range(n):
            tmp = self.geocpp.getMinoTimeOfTime(t[i])
            ph[i] = self.geocpp.getAzimuthalPositionOfMinoTime(tmp)

    def psi_radial(self, double la):
        return self.geocpp.getPsiRadialOfMinoTime(la)

    def psi_polar(self, double la):
        return self.geocpp.getPsiPolarOfMinoTime(la)

    def psi_radial_time(self, double t):
        cdef double tmp = self.geocpp.getMinoTimeOfTime(t)
        return self.geocpp.getPsiRadialOfMinoTime(tmp)

    def psi_polar_time(self, double t):
        cdef double tmp = self.geocpp.getMinoTimeOfTime(t)
        return self.geocpp.getPsiPolarOfMinoTime(tmp)

    def psi_radial_vec(self, np.ndarray[ndim=1, dtype=np.float64_t] la):
        cdef int n = la.shape[0]
        cdef np.ndarray[ndim=1, dtype=np.float64_t] psi = np.empty(n, dtype = np.float64)
        self.getPsiRadialOfMinoTimeArray(&psi[0], &la[0], n)
        return psi

    def psi_polar_vec(self, np.ndarray[ndim=1, dtype=np.float64_t] la):
        cdef int n = la.shape[0]
        cdef np.ndarray[ndim=1, dtype=np.float64_t] psi = np.empty(n, dtype = np.float64)
        self.getPsiPolarOfMinoTimeArray(&psi[0], &la[0], n)
        return psi

    def psi_radial_time_vec(self, np.ndarray[ndim=1, dtype=np.float64_t] t):
        cdef int n = t.shape[0]
        cdef np.ndarray[ndim=1, dtype=np.float64_t] psi = np.empty(n, dtype = np.float64)
        self.getPsiRadialOfTimeArray(&psi[0], &t[0], n)
        return psi

    def psi_polar_time_vec(self, np.ndarray[ndim=1, dtype=np.float64_t] t):
        cdef int n = t.shape[0]
        cdef np.ndarray[ndim=1, dtype=np.float64_t] psi = np.empty(n, dtype = np.float64)
        self.getPsiPolarOfTimeArray(&psi[0], &t[0], n)
        return psi

    def time_position(self, double la):
        return self.geocpp.getTimePositionOfMinoTime(la)

    def radial_position(self, double la):
        return self.geocpp.getRadialPositionOfMinoTime(la)
    
    def polar_position(self, double la):
        return self.geocpp.getPolarPositionOfMinoTime(la)
    
    def azimuthal_position(self, double la):
        return self.geocpp.getAzimuthalPositionOfMinoTime(la)

    def radial_position_time(self, double t):
        cdef double tmp = self.geocpp.getMinoTimeOfTime(t)
        return self.geocpp.getRadialPositionOfMinoTime(tmp)

    def polar_position_time(self, double t):
        cdef double tmp = self.geocpp.getMinoTimeOfTime(t)
        return self.geocpp.getPolarPositionOfMinoTime(tmp)

    def azimuthal_position_time(self, double t):
        cdef double tmp = self.geocpp.getMinoTimeOfTime(t)
        return self.geocpp.getAzimuthalPositionOfMinoTime(tmp)

    def time_position_vec(self, np.ndarray[ndim=1, dtype=np.float64_t] la):
        cdef int n = la.shape[0]
        cdef np.ndarray[ndim=1, dtype=np.float64_t] t = np.empty(n, dtype = np.float64)
        self.getTimePositionOfMinoTimeArray(&t[0], &la[0], n)
        # for some reason this is marginally slower than the line above with the cdef function
        # for i in range(n):
        #     t[i] = self.geocpp.getTimePositionOfMinoTime(la[i])
        return t
    
    def radial_position_vec(self, np.ndarray[ndim=1, dtype=np.float64_t] la):
        cdef int n = la.shape[0]
        cdef np.ndarray[ndim=1, dtype=np.float64_t] x = np.empty(n, dtype = np.float64)
        self.getRadialPositionOfMinoTimeArray(&x[0], &la[0], n)
        return x

    def polar_position_vec(self, np.ndarray[ndim=1, dtype=np.float64_t] la):
        cdef int n = la.shape[0]
        cdef np.ndarray[ndim=1, dtype=np.float64_t] x = np.empty(n, dtype = np.float64)
        self.getPolarPositionOfMinoTimeArray(&x[0], &la[0], n)
        return x

    def azimuthal_position_vec(self, np.ndarray[ndim=1, dtype=np.float64_t] la):
        cdef int n = la.shape[0]
        cdef np.ndarray[ndim=1, dtype=np.float64_t] x = np.empty(n, dtype = np.float64)
        self.getAzimuthalPositionOfMinoTimeArray(&x[0], &la[0], n)
        return x

    def radial_position_time_vec(self, np.ndarray[ndim=1, dtype=np.float64_t] t):
        cdef int n = t.shape[0]
        cdef np.ndarray[ndim=1, dtype=np.float64_t] x = np.empty(n, dtype = np.float64)
        self.getRadialPositionOfTimeArray(&x[0], &t[0], n)
        return x

    def polar_position_time_vec(self, np.ndarray[ndim=1, dtype=np.float64_t] t):
        cdef int n = t.shape[0]
        cdef np.ndarray[ndim=1, dtype=np.float64_t] x = np.empty(n, dtype = np.float64)
        self.getPolarPositionOfTimeArray(&x[0], &t[0], n)
        return x

    def azimuthal_position_time_vec(self, np.ndarray[ndim=1, dtype=np.float64_t] t):
        cdef int n = t.shape[0]
        cdef np.ndarray[ndim=1, dtype=np.float64_t] x = np.empty(n, dtype = np.float64)
        self.getAzimuthalPositionOfTimeArray(&x[0], &t[0], n)
        return x

    def position(self, double la):
        return np.array([self.geocpp.getTimePositionOfMinoTime(la), self.geocpp.getRadialPositionOfMinoTime(la), self.geocpp.getPolarPositionOfMinoTime(la), self.geocpp.getAzimuthalPositionOfMinoTime(la)])

    def position_time(self, double t):
        cdef double tmp = self.geocpp.getMinoTimeOfTime(t)
        return np.array([self.geocpp.getRadialPositionOfMinoTime(tmp), self.geocpp.getPolarPositionOfMinoTime(tmp), self.geocpp.getAzimuthalPositionOfMinoTime(tmp)])

    def position_vec(self, np.ndarray[ndim=1, dtype=np.float64_t] la):
        cdef np.ndarray[ndim=2, dtype=np.float64_t] xp = np.empty((la.shape[0], 4), dtype=np.float64)
        for i in range(la.shape[0]):
            xp[i] = self.position(la[i])
        return xp.T

    def position_time_vec(self, np.ndarray[ndim=1, dtype=np.float64_t] t):
        cdef np.ndarray[ndim=2, dtype=np.float64_t] xp = np.empty((t.shape[0], 3), dtype=np.float64)
        for i in range(t.shape[0]):
            xp[i] = self.position_time(t[i])
        return xp.T

    def mino_time(self, double t):
        return self.geocpp.getMinoTimeOfTime(t)

    def mino_time_vec(self, np.ndarray[ndim=1, dtype=np.float64_t] t):
        cdef int n = t.shape[0]
        cdef np.ndarray[ndim=1, dtype=np.float64_t] la = np.empty(n, dtype = np.float64)
        for i in range(n):
            la[i] = self.mino_time(t[i])
        return la
    
    def get_time_accumulation(self, int j):
        cdef vector[double] deltaX_cpp = self.geocpp.getTimeAccumulation(j)
        cdef int n = deltaX_cpp.size()
        cdef np.ndarray[ndim=1, dtype=np.float64_t] deltaX = np.empty(n, dtype = np.float64)
        for i in range(n):
            deltaX[i] = deltaX_cpp[i]
        return deltaX

    def get_radial_points(self):
        cdef vector[double] deltaX_cpp = self.geocpp.getRadialPosition()
        cdef int n = deltaX_cpp.size()
        cdef np.ndarray[ndim=1, dtype=np.float64_t] deltaX = np.empty(n, dtype = np.float64)
        for i in range(n):
            deltaX[i] = deltaX_cpp[i]
        return deltaX

    def get_polar_points(self):
        cdef vector[double] deltaX_cpp = self.geocpp.getPolarPosition()
        cdef int n = deltaX_cpp.size()
        cdef np.ndarray[ndim=1, dtype=np.float64_t] deltaX = np.empty(n, dtype = np.float64)
        for i in range(n):
            deltaX[i] = deltaX_cpp[i]
        return deltaX
    
    def get_azimuthal_accumulation(self, int j):
        cdef vector[double] deltaX_cpp = self.geocpp.getAzimuthalAccumulation(j)
        cdef int n = deltaX_cpp.size()
        cdef np.ndarray[ndim=1, dtype=np.float64_t] deltaX = np.empty(n, dtype = np.float64)
        for i in range(n):
            deltaX[i] = deltaX_cpp[i]
        return deltaX

    def get_time_coefficients(self, int j):
        cdef vector[double] deltaX_cpp = self.geocpp.getTimeCoefficients(j)
        cdef int n = deltaX_cpp.size()
        cdef np.ndarray[ndim=1, dtype=np.float64_t] deltaX = np.empty(n, dtype = np.float64)
        for i in range(n):
            deltaX[i] = deltaX_cpp[i]
        return deltaX

    def get_radial_coefficients(self):
        cdef vector[double] deltaX_cpp = self.geocpp.getRadialCoefficients()
        cdef int n = deltaX_cpp.size()
        cdef np.ndarray[ndim=1, dtype=np.float64_t] deltaX = np.empty(n, dtype = np.float64)
        for i in range(n):
            deltaX[i] = deltaX_cpp[i]
        return deltaX

    def get_polar_coefficients(self):
        cdef vector[double] deltaX_cpp = self.geocpp.getPolarCoefficients()
        cdef int n = deltaX_cpp.size()
        cdef np.ndarray[ndim=1, dtype=np.float64_t] deltaX = np.empty(n, dtype = np.float64)
        for i in range(n):
            deltaX[i] = deltaX_cpp[i]
        return deltaX
    
    def get_azimuthal_coefficients(self, int j):
        cdef vector[double] deltaX_cpp = self.geocpp.getAzimuthalCoefficients(j)
        cdef int n = deltaX_cpp.size()
        cdef np.ndarray[ndim=1, dtype=np.float64_t] deltaX = np.empty(n, dtype = np.float64)
        for i in range(n):
            deltaX[i] = deltaX_cpp[i]
        return deltaX

def _kerr_kepler_parameters_wrapper(double a, double En, double Lz, double Qc):
    cdef double p, e, x
    p = 0.
    e = 0.
    x = 0.
    kerr_geo_kepler_parameters(p, e, x, a, En, Lz, Qc)
    return np.array([p, e, x])

def _kerr_kepler_parameters_array_wrapper(double[:] a, double[:] E, double[:] Lz, double[:] Q):
    """
    Zero-copy interface using NumPy memoryviews.
    """
    cdef int n = E.shape[0]
    
    # Pre-allocate output NumPy arrays
    p_out = np.empty(n, dtype=np.float64)
    e_out = np.empty(n, dtype=np.float64)
    xI_out = np.empty(n, dtype=np.float64)
    
    # Cast to memoryviews to get raw pointers
    cdef double[:] p_mv = p_out
    cdef double[:] e_mv = e_out
    cdef double[:] xI_mv = xI_out
    
    # Execute C++ loop
    kerr_geo_kepler_parameters(n, &p_mv[0], &e_mv[0], &xI_mv[0], &a[0], &E[0], &Lz[0], &Q[0])
    
    return np.array([p_out, e_out, xI_out])

def _kerr_orbital_constants_wrapper(double a, double p, double e, double x):
    cdef double En, Lz, Qc
    En = 0.
    Lz = 0.
    Qc = 0.
    kerr_geo_orbital_constants(En, Lz, Qc, a, p, e, x)
    return np.array([En, Lz, Qc])

def _kerr_orbital_constants_array_wrapper(double[:] a, double[:] p, double[:] e, double[:] x):
    """
    Zero-copy interface using NumPy memoryviews.
    """
    cdef int n = a.shape[0]

    # Pre-allocate output NumPy arrays
    En_out = np.empty(n, dtype=np.float64)
    Lz_out = np.empty(n, dtype=np.float64)
    Qc_out = np.empty(n, dtype=np.float64)
    
    # Cast to memoryviews to get raw pointers
    cdef double[:] En_mv = En_out
    cdef double[:] Lz_mv = Lz_out
    cdef double[:] Qc_mv = Qc_out
    
    # Execute C++ loop
    kerr_geo_orbital_constants(n, &En_mv[0], &Lz_mv[0], &Qc_mv[0], &a[0], &p[0], &e[0], &x[0])
    
    return np.array([En_out, Lz_out, Qc_out])

def _kerr_mino_frequencies_wrapper(double a, double p, double e, double x):
    cdef double upT, upR, upTh, upPhi
    upT = 0.
    upR = 0.
    upTh = 0.
    upPhi = 0.
    kerr_geo_mino_frequencies(upT, upR, upTh, upPhi, a, p, e, x)
    return np.array([upT, upR, upTh, upPhi])

def _kerr_radial_roots_wrapper(double a, double p, double e, double En, double Lz, double Qc):
    cdef double r1, r2, r3, r4
    r1 = 0.
    r2 = 0.
    r3 = 0.
    r4 = 0.
    kerr_geo_radial_roots(r1, r2, r3, r4, a, p, e, En, Lz, Qc)
    return np.array([r1, r2, r3, r4])

def _kerr_polar_roots_wrapper(double a, double x, double En, double Lz, double Qc):
    cdef double z1, z2
    z1 = 0.
    z2 = 0.
    kerr_geo_polar_roots(z1, z2, a, x, En, Lz, Qc)
    return np.array([z1, z2])

def _kerr_isco_wrapper(double a, int sgnX):
    return kerr_isco(a, sgnX)

def _kerr_isco_frequency_wrapper(double a):
    return kerr_isco_frequency(a)

def _jacobian_ELQ_to_pex_wrapper(double a, double p, double e, double x):
    cdef double dpdE, dedE, dxdE
    cdef double dpdLz, dedLz, dxdLz
    cdef double dpdQ, dedQ, dxdQ
    dpdE = 0.
    dedE = 0.
    dxdE = 0.
    dpdLz = 0.
    dedLz = 0.
    dxdLz = 0.
    dpdQ = 0.
    dedQ = 0.
    dxdQ = 0.
    jacobian_ELQ_to_pex(dpdE, dedE, dxdE, dpdLz, dedLz, dxdLz, dpdQ, dedQ, dxdQ, a, p, e, x)
    return np.array([[dpdE, dpdLz, dpdQ],
                     [dedE, dedLz, dedQ],
                     [dxdE, dxdLz, dxdQ]])

def _jacobian_pex_to_ELQ_wrapper(double a, double p, double e, double x):
    cdef double dEdp, dEde, dEdx
    cdef double dLdp, dLde, dLdx
    cdef double dQdp, dQde, dQdx
    dEdp = 0.
    dEde = 0.
    dEdx = 0.
    dLdp = 0.
    dLde = 0.
    dLdx = 0.
    dQdp = 0.
    dQde = 0.
    dQdx = 0.
    jacobian_pex_to_ELQ(dEdp, dEde, dEdx, dLdp, dLde, dLdx, dQdp, dQde, dQdx, a, p, e, x)
    return np.array([[dEdp, dEde, dEdx],
                     [dLdp, dLde, dLdx],
                     [dQdp, dQde, dQdx]])

def _jacobian_ELQ_to_pex_array_wrapper(double[:] a, double[:] p, double[:] e, double[:] x):
    """
    Zero-copy interface using NumPy memoryviews.
    """
    cdef int n = a.shape[0]
    
    # Pre-allocate output NumPy arrays
    dpdE_out = np.empty(n, dtype=np.float64)
    dedE_out = np.empty(n, dtype=np.float64)
    dxdE_out = np.empty(n, dtype=np.float64)
    dpdLz_out = np.empty(n, dtype=np.float64)
    dedLz_out = np.empty(n, dtype=np.float64)
    dxdLz_out = np.empty(n, dtype=np.float64)
    dpdQ_out = np.empty(n, dtype=np.float64)
    dedQ_out = np.empty(n, dtype=np.float64)
    dxdQ_out = np.empty(n, dtype=np.float64)
    
    # Cast to memoryviews to get raw pointers
    cdef double[:] dpdE_mv = dpdE_out
    cdef double[:] dedE_mv = dedE_out
    cdef double[:] dxdE_mv = dxdE_out
    cdef double[:] dpdLz_mv = dpdLz_out
    cdef double[:] dedLz_mv = dedLz_out
    cdef double[:] dxdLz_mv = dxdLz_out
    cdef double[:] dpdQ_mv = dpdQ_out
    cdef double[:] dedQ_mv = dedQ_out
    cdef double[:] dxdQ_mv = dxdQ_out
    
    # Execute C++ loop
    jacobian_ELQ_to_pex(n, &dpdE_mv[0], &dedE_mv[0], &dxdE_mv[0],
                        &dpdLz_mv[0], &dedLz_mv[0], &dxdLz_mv[0],
                        &dpdQ_mv[0], &dedQ_mv[0], &dxdQ_mv[0],
                        &a[0], &p[0], &e[0], &x[0])
    
    return np.array([[dpdE_out, dpdLz_out, dpdQ_out],
                     [dedE_out, dedLz_out, dedQ_out],
                     [dxdE_out, dxdLz_out, dxdQ_out]])

def _jacobian_pex_to_ELQ_array_wrapper(double[:] a, double[:] p, double[:] e, double[:] x):
    """
    Zero-copy interface using NumPy memoryviews.
    """
    cdef int n = a.shape[0]

    # Pre-allocate output NumPy arrays
    dEdp_out = np.empty(n, dtype=np.float64)
    dEde_out = np.empty(n, dtype=np.float64)
    dEdx_out = np.empty(n, dtype=np.float64)
    dLdp_out = np.empty(n, dtype=np.float64)
    dLde_out = np.empty(n, dtype=np.float64)
    dLdx_out = np.empty(n, dtype=np.float64)
    dQdp_out = np.empty(n, dtype=np.float64)
    dQde_out = np.empty(n, dtype=np.float64)
    dQdx_out = np.empty(n, dtype=np.float64)

    # Cast to memoryviews to get raw pointers
    cdef double[:] dEdp_mv = dEdp_out
    cdef double[:] dEde_mv = dEde_out
    cdef double[:] dEdx_mv = dEdx_out
    cdef double[:] dLdp_mv = dLdp_out
    cdef double[:] dLde_mv = dLde_out
    cdef double[:] dLdx_mv = dLdx_out
    cdef double[:] dQdp_mv = dQdp_out
    cdef double[:] dQde_mv = dQde_out
    cdef double[:] dQdx_mv = dQdx_out

    # Execute C++ loop
    jacobian_pex_to_ELQ(n, &dEdp_mv[0], &dEde_mv[0], &dEdx_mv[0],
                        &dLdp_mv[0], &dLde_mv[0], &dLdx_mv[0],
                        &dQdp_mv[0], &dQde_mv[0], &dQdx_mv[0],
                        &a[0], &p[0], &e[0], &x[0])
    
    return np.array([[dEdp_out, dEde_out, dEdx_out],
                     [dLdp_out, dLde_out, dLdx_out],
                     [dQdp_out, dQde_out, dQdx_out]])