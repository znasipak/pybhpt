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
    void kerr_geo_orbital_constants(double &En, double &Lz, double &Qc, double &a, double &p, double &e, double &x)
    void kerr_geo_radial_roots(double &r1, double &r2, double &r3, double &r4, double &a, double &p, double &e, double &En, double &Lz, double &Qc)
    void kerr_geo_polar_roots(double &z1, double &z2, double &a, double &x, double &En, double &Lz, double &Qc)
    void kerr_geo_mino_frequencies(double &upT, double &upR, double &upTh, double &upPh, double &a, double &p, double &e, double &x)

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

def kerr_geo_V01(double a, double En, double Lz, double Q, double r):
    return kerr_geo_VtR(a, En, Lz, Q, r)

def kerr_geo_V02(double a, double En, double Lz, double Q, double theta):
    return kerr_geo_VtTheta(a, En, Lz, Q, theta)

def kerr_geo_V11(double a, double En, double Lz, double Q, double r):
    return kerr_geo_Vr(a, En, Lz, Q, r)

def kerr_geo_V22(double a, double En, double Lz, double Q, double theta):
    return kerr_geo_Vtheta(a, En, Lz, Q, theta)

def kerr_geo_V31(double a, double En, double Lz, double Q, double r):
    return kerr_geo_VphiR(a, En, Lz, Q, r)

def kerr_geo_V32(double a, double En, double Lz, double Q, double theta):
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


def kerr_orbital_constants_wrapper(double a, double p, double e, double x):
    cdef double En, Lz, Qc
    En = 0.
    Lz = 0.
    Qc = 0.
    kerr_geo_orbital_constants(En, Lz, Qc, a, p, e, x)
    return np.array([En, Lz, Qc])

def kerr_mino_frequencies_wrapper(double a, double p, double e, double x):
    cdef double upT, upR, upTh, upPhi
    upT = 0.
    upR = 0.
    upTh = 0.
    upPhi = 0.
    kerr_geo_mino_frequencies(upT, upR, upTh, upPhi, a, p, e, x)
    return np.array([upT, upR, upTh, upPhi])
        