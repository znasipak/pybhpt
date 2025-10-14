// geo.hpp

#ifndef GEO_HPP
#define GEO_HPP

#include "utils.hpp"
#include <chrono>
#include "specialfunc.hpp"
#include <gsl/gsl_roots.h>

class GeodesicTrajectory{
public:
	GeodesicTrajectory() {};
	GeodesicTrajectory(Vector tR, Vector tTheta, Vector r, Vector theta, Vector phiR, Vector phiTheta);
	~GeodesicTrajectory() {};
	Vector tR;
	Vector tTheta;
	Vector r;
	Vector theta;
	Vector phiR;
	Vector phiTheta;

	double getTimeAccumulationRadial(int pos);
	double getTimeAccumulationPolar(int pos);
	double getTimeAccumulation(int j, int pos);
	double getRadialPosition(int pos);
	double getPolarPosition(int pos);
	double getAzimuthalAccumulationRadial(int pos);
	double getAzimuthalAccumulationPolar(int pos);
	double getAzimuthalAccumulation(int j, int pos);
};

class GeodesicConstants{
public:
	GeodesicConstants(){};
	GeodesicConstants(double a, double p, double e, double x, double En, double Lz, double Q);
	GeodesicConstants(double a, double p, double e, double x, double En, double Lz, double Q, double upsilonT, double upsilonR, double upsilonTheta, double upsilonPhi);
	GeodesicConstants(double a, double p, double e, double x, double En, double Lz, double Q, double r1, double r2, double r3, double r4, double z1, double z2, double upsilonT, double upsilonR, double upsilonTheta, double upsilonPhi, double carterR, double carterTheta, double carterPhi);
	~GeodesicConstants(){};
	double a;
	double p;
	double e;
	double x;

	double En;
	double Lz;
	double Q;

	double r1;
	double r2;
	double r3;
	double r4;
	double z1;
	double z2;

	double upsilonT;
	double upsilonR;
	double upsilonTheta;
	double upsilonPhi;

	double carterR;
	double carterTheta;
	double carterPhi;

	double getTimeFrequency(int m, int k, int n);
};

class GeodesicSource{
public:
	GeodesicSource(){};
	GeodesicSource(double a, double p, double e, double x, int Nsample);
	~GeodesicSource(){};

	int getOrbitalSampleNumber();
	double getBlackHoleSpin();
	double getSemiLatusRectum();
	double getEccentricity();
	double getInclination();
	double getOrbitalEnergy();
	double getOrbitalAngularMomentum();
	double getCarterConstant();
	double getRadialRoot(int i);
	double getPolarRoot(int i);
	double getMinoFrequency(int mu);
	double getTimeFrequency(int i);
	double getTimeFrequency(int m, int k, int n);
	double getCarterFrequency(int i);
	double getCarterFrequency(int m, int k, int n);

	Vector getTimeAccumulation(int j);
	Vector getRadialPosition();
	Vector getPolarPosition();
	Vector getAzimuthalAccumulation(int j);

	double getTimeAccumulation(int j, int pos);
	double getRadialPosition(int pos);
	double getPolarPosition(int pos);
	double getAzimuthalAccumulation(int j, int pos);

	double getPsiRadialOfMinoTime(double lambda);
	double getPsiPolarOfMinoTime(double lambda);

	double getTimePositionOfMinoTime(double lambda);
	double getRadialPositionOfMinoTime(double lambda);
	double getPolarPositionOfMinoTime(double lambda);
	double getAzimuthalPositionOfMinoTime(double lambda);

	Vector getPsiRadialOfMinoTime(Vector lambda);
	Vector getPsiPolarOfMinoTime(Vector lambda);

	Vector getTimePositionOfMinoTime(Vector lambda);
	Vector getRadialPositionOfMinoTime(Vector lambda);
	Vector getPolarPositionOfMinoTime(Vector lambda);
	Vector getAzimuthalPositionOfMinoTime(Vector lambda);
	Vector getPositionOfMinoTime(double lambda);

	double getMinoTimeOfTime(double t);

	Vector getTimeCoefficients(int j);
	Vector getRadialCoefficients();
	Vector getPolarCoefficients();
	Vector getAzimuthalCoefficients(int j);

	GeodesicConstants getConstants();
	GeodesicTrajectory getTrajectory();
	GeodesicTrajectory getCoefficients();
	
	GeodesicConstants& getConstantsRef();
	GeodesicTrajectory& getTrajectoryRef();
	GeodesicTrajectory& getCoefficientsRef();

	void setConstants(double a, double p, double e, double x,
		double En, double Lz, double Qc, double r1, double r2, double r3, double r4, double z1, double z2, double upT, double upR, double upTh, double upPhi, double cR, double cTh, double cPh);
	void setTrajectory(Vector tR, Vector tTheta, Vector r, Vector theta, Vector phiR, Vector phiTheta);
	void setCoefficients(Vector tR, Vector tTheta, Vector r, Vector theta, Vector phiR, Vector phiTheta);
	// void save(std::string dir);

private:
	GeodesicConstants _geoConstants;
	GeodesicTrajectory _geoTrajectory;
	GeodesicTrajectory _geoCoefficients;
};

// void output_geodesic(GeodesicSource geo, std::string dir);
// GeodesicSource load_geodesic(std::string dir);

double kerr_geo_energy(double a, double p, double e, double x);
double kerr_geo_momentum(double a, double p, double e, double x);
double kerr_geo_carter(double a, double p, double e, double x);
double kerr_geo_momentum(double En, double a, double p, double e, double x);
double kerr_geo_carter(double En, double Lz, double a, double p, double e, double x);

void kerr_geo_orbital_constants(double &En, double &Lz, double &Qc, double a, double p, double e, double x);
void kerr_geo_radial_roots(double &r1, double &r2, double &r3, double &r4, double a, double p, double e, double En, double Lz, double Qc);
void kerr_geo_polar_roots(double &z1, double &z2, double a, double x, double En, double Lz, double Qc);

Vector kerr_delta(Vector &r, double a);
void kerr_delta(Vector &delta, Vector &r, double a);
double kerr_delta(double r, double a);

Vector kerr_varpi2(Vector &r, double a);
void kerr_varpi2(Vector &varpi2, Vector &r, double a);
double kerr_varpi2(double r, double a);

Vector rp_psi(Vector &psi, double p, double e);
void rp_psi(Vector &rp, Vector &psi, double p, double e);
Vector zp_chi(Vector &chi, double x);
void zp_chi(Vector &zp, Vector &chi, double x);

double rp_psi(double psi, double p, double e);
double zp_chi(double chi, double x);

double dmino_dpsi(double psi, double a, double p, double e, double En, double r3, double const&r4);
double dmino_dchi(double chi, double a, double En, double z1, double z2);
double dmino_dchi_schw(double Lz, double Qc);
double dVtrdMino(double rp, double a, double En, double Lz, double Qc);
double dVtzdMino(double zp, double a, double En, double Lz, double Qc);
double dVphirdMino(double rp, double a, double En, double Lz, double Qc);
double dVphizdMino(double zp, double a, double En, double Lz, double Qc);

void kerr_geo_radial_mino_period(double &LaR, double a, double p, double e, double En, double r3, double r4);
void kerr_geo_radial_mino_frequency(double &upR, double a, double p, double e, double En, double r3, double r4);
void kerr_geo_polar_mino_period(double &LaTh, double a, double En, double z1, double z2);
void kerr_geo_polar_mino_frequency(double &upTh, double a, double En, double z1, double z2);
void kerr_geo_polar_mino_period_schw(double &LaTh, double Lz, double Qc);
void kerr_geo_polar_mino_frequency_schw(double &LaTh, double Lz, double Qc);

void kerr_geo_time_mino_frequency_radial(double &upT, double a, double p, double e, double En, double Lz, double Qc, double r3, double r4);
void kerr_geo_time_mino_frequency_polar(double &upT, double a, double En, double Lz, double Qc, double z1, double z2);
void kerr_geo_time_mino_frequency_polar_schw(double &upT, double En, double Lz, double Qc, double z1, double z2);
void kerr_geo_time_mino_frequency(double &upT, double a, double p, double e, double En, double Lz, double Qc, double z1, double z2, double r3, double r4, double upR, double upTh);
void kerr_geo_azimuthal_mino_frequency_radial(double &upPh, double a, double p, double e, double En, double Lz, double Qc, double r3, double r4);
void kerr_geo_azimuthal_mino_frequency_polar(double &upPh, double a, double En, double Lz, double Qc, double z1, double z2);
void kerr_geo_azimuthal_mino_frequency_polar_schw(double &upPh, double En, double Lz, double Qc, double z1, double z2);
void kerr_geo_azimuthal_mino_frequency(double &upPh, double a, double p, double e, double En, double Lz, double Qc, double z1, double z2, double r3, double r4, double upR, double upTh);

void kerr_geo_mino_frequencies(double &upT, double &upR, double &upTh, double &upPh, double a, double p, double e, double x);
void kerr_geo_mino_frequencies(double &upT, double &upR, double &upTh, double &upPh, double a, double p, double e, double x, double En, double Lz, double Qc, double r1, double r2, double r3, double r4, double z1, double z2);
void kerr_geo_carter_frequencies(double &cR, double &cTh, double &cPh, double upT, double upR, double upTh, double upPh, double a, double En, double Lz, double Qc, double z1, double z2);

// void kerr_geo_time(Vector &tR, Vector &tTh, double upR, double upTh, double a, double p, double e, double x, double En, double Lz, double Qc, double r1, double r2, double r3, double r4, double z1, double z2);
// void kerr_geo_radial_time(Vector &tR, double upR, double a, double p, double e, double En, double Lz, double r1, double r2, double r3, double r4);
// void kerr_geo_polar_time(Vector &tTh, double upTh, double a, double x, double En, double Lz, double Qc, double z1, double z2);

Vector angle_vector(int Nsample);
void angle_vector(Vector &angle);
Vector angle_vector_half(int Nsample);
void angle_vector_half(Vector &angle);

Vector mino_of_psi_fourier(double a, double p, double e, double En, double r3, double r4);
Vector mino_of_chi_fourier(double a, double En, double z1, double z2);
Vector mino_of_chi_schw_fourier(double Lz, double Qc);
double integrated_dct_sum(double angle, Vector &fourier);
Vector kepler_phase_of_angle_fourier(Vector &fourier);
double mino_of_kepler_phase(double angle, Vector &fourier);
double kepler_phase_of_angle(double angle, Vector &fourier);

double rp_of_angle(double angle, double p, double e, Vector &fourier);
double zp_of_angle(double angle, double x, Vector &fourier);
Vector tp_polar_of_angle_fourier(double a, double x, double En, double Lz, double Qc, double upR, Vector &fourier);
Vector tp_radial_of_angle_fourier(double a, double p, double e, double En, double Lz, double Qc, double upTh, Vector &fourier);
Vector phip_polar_of_angle_fourier(double a, double x, double En, double Lz, double Qc, double upR, Vector &fourier);
Vector phip_radial_of_angle_fourier(double a, double p, double e, double En, double Lz, double Qc, double upTh, Vector &fourier);
double tp_of_angle(double angle, Vector &fourier);
double phip_of_angle(double angle, Vector &fourier);

// double lambda_of_time(double t, GeodesicSource & geo);

void kerr_trajectory(Vector& tR, Vector& tTh, Vector& rp, Vector& thetap, Vector& phiR, Vector& phiTh, double p, double e, double x, Vector fourier_tr, Vector fourier_tz, Vector fourier_psi, Vector fourier_chi, Vector fourier_phir, Vector fourier_phiz);

GeodesicSource kerr_geo_orbit(double a, double p, double e, double x, int Nsample);
// void mino_of_psi_test();


///////////////////////
// Special functions //
///////////////////////

double kerr_geo_radial_mino_frequency_sf(double En, double r1, double r2, double r3, double r4);
double kerr_geo_polar_mino_frequency_sf(double a, double En, double z1, double z2);
void kerr_geo_radial_mino_frequency_sf(double &upR, double En, double r1, double r2, double r3, double r4);
void kerr_geo_polar_mino_frequency_sf(double &upTh, double a, double En, double z1, double z2);
double rp_of_angle_sf(double qr, double r1, double r2, double r3, double r4);
double zp_of_angle_sf(double qth, double z1, double z2);

/////////////////////
// Circular Orbits //
/////////////////////

GeodesicSource kerr_geo_circ(double a, double r, int sgnX);
double kerr_geo_energy_circ(double a, double r);
double kerr_geo_momentum_circ(double a, double r);
double kerr_geo_time_frequency_circ(double a, double r);
double kerr_geo_azimuthal_frequency_circ(double a, double r);

double kerr_geo_azimuthal_frequency_circ_time(double a, double r, int sgnX);
double kerr_geo_azimuthal_frequency_circ_time(double a, double r);
double kerr_geo_radius_circ(double a, double Omega);

double kerr_geo_denergy_domega_circ(double a, double om);
double kerr_geo_dmomentum_domega_circ(double a, double om);
double dr_domega(double a, double om);
double denergy_dr(double a, double r);
double dmomentum_dr(double a, double r);

double kerr_geo_VtR(double a, double En, double Lz, double Q, double r);
double kerr_geo_VtTheta(double a, double En, double Lz, double Q, double theta);
double kerr_geo_Vr(double a, double En, double Lz, double Q, double r);
double kerr_geo_Vtheta(double a, double En, double Lz, double Q, double theta);
double kerr_geo_VphiR(double a, double En, double Lz, double Q, double r);
double kerr_geo_VphiTheta(double a, double En, double Lz, double Q, double theta);
double kerr_geo_Vz(double a, double En, double Lz, double Q, double z);
double kerr_geo_Vz_dz(double a, double En, double Lz, double Q, double z);
double kerr_geo_Vz_dz2(double a, double En, double Lz, double Q, double z);
double kerr_isco(double a, int sgnX);
double kerr_isco_frequency(double a);

#endif
