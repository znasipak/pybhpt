// sourceintegration.hpp

#ifndef SOURCEINT_HPP
#define SOURCEINT_HPP

#include "geo.hpp"
#include "specialfunc.hpp"
#include "radialsolver.hpp"
#include "swsh.hpp"
#include <functional>

typedef struct TeukolskyAmplitudesStruct{
	Complex in;
	Complex up;
	double inPrecision;
	double upPrecision;
} TeukolskyAmplitudes;

typedef struct DerivativesMatrixStruct{
	Vector solution;
	Vector derivative;
	Vector secondDerivative;
} DerivativesMatrix;

typedef struct ComplexDerivativesMatrixStruct{
	ComplexVector solution;
	ComplexVector derivative;
	ComplexVector secondDerivative;
} ComplexDerivativesMatrix;

class SummationHelper{
public:
	SummationHelper();
	~SummationHelper();

	Complex getSum();
	double getError();
	double getPrecision();
	double getMaxTerm();

	void setBasePrecision(double val);
	void add(Complex val);

private:
	Complex _sum;
	Complex _previousSum;
	double _maxTerm;
	double _error;
	double _basePrecision;
};

TeukolskyAmplitudes field_amplitude_circeq(int s, int L, int m, GeodesicTrajectory& traj, GeodesicConstants &geoConst, RadialTeukolsky &teuk, SpinWeightedHarmonic &swsh);
TeukolskyAmplitudes field_amplitude_ecceq(int s, int L, int m, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConst, RadialTeukolsky &teuk, SpinWeightedHarmonic &swsh);
TeukolskyAmplitudes field_amplitude_sphinc(int s, int L, int m, int k, GeodesicTrajectory& traj, GeodesicConstants &geoConst, RadialTeukolsky &teuk, SpinWeightedHarmonic &swsh);
TeukolskyAmplitudes field_amplitude(int s, int L, int m, int k, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConst, RadialTeukolsky &teuk, SpinWeightedHarmonic &swsh);

// TeukolskyAmplitudes field_amplitude_circeq(int s, int L, int m, GeodesicTrajectory& traj, GeodesicConstants &geoConst, ComplexDerivativesMatrixStruct Rin, ComplexDerivativesMatrixStruct Rup, DerivativesMatrix Slm);
// TeukolskyAmplitudes field_amplitude_ecceq(int s, int L, int m, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConst, ComplexDerivativesMatrixStruct Rin, ComplexDerivativesMatrixStruct Rup, DerivativesMatrix Slm);
// TeukolskyAmplitudes field_amplitude_sphinc(int s, int L, int m, int k, GeodesicTrajectory& traj, GeodesicConstants &geoConst, ComplexDerivativesMatrixStruct Rin, ComplexDerivativesMatrixStruct Rup, DerivativesMatrix Slm);
// TeukolskyAmplitudes field_amplitude(int s, int L, int m, int k, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConst, ComplexDerivativesMatrixStruct Rin, ComplexDerivativesMatrixStruct Rup, DerivativesMatrix Slm);

TeukolskyAmplitudes teukolsky_amplitude_circeq(int s, int L, int m, GeodesicTrajectory& traj, GeodesicConstants &geoConst, ComplexDerivativesMatrixStruct Rin, ComplexDerivativesMatrixStruct Rup, DerivativesMatrix Slm);
// TeukolskyAmplitudes teukolsky_amplitude_circeq_plus_2(int L, int m, GeodesicTrajectory& traj, GeodesicConstants &geoConst, ComplexDerivativesMatrixStruct Rin, ComplexDerivativesMatrixStruct Rup, DerivativesMatrix Slm);
TeukolskyAmplitudes teukolsky_amplitude_ecceq(int s, int L, int m, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConst, ComplexDerivativesMatrixStruct Rin, ComplexDerivativesMatrixStruct Rup, DerivativesMatrix Slm);
TeukolskyAmplitudes teukolsky_amplitude_sphinc(int s, int L, int m, int k, GeodesicTrajectory& traj, GeodesicConstants &geoConst, ComplexDerivativesMatrixStruct Rin, ComplexDerivativesMatrixStruct Rup, DerivativesMatrix Slm);
TeukolskyAmplitudes teukolsky_amplitude(int s, int L, int m, int k, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConst, ComplexDerivativesMatrixStruct Rin, ComplexDerivativesMatrixStruct Rup, DerivativesMatrix Slm);
Complex teukolskyIntegrand(int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &tR, double const &tTh, double const &rp, double const &thp, double const &phiR, double const &phiTh, double const &qr, double const &qth, Complex const &Rt, Complex const &RtP, Complex const &RtPP, double const &St, double const &StP, double const &StPP);
Complex teukolskyIntegrandPlus(int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &tR, double const &tTh, double const &rp, double const &thp, double const &phiR, double const &phiTh, double const &qr, double const &qth, Complex const &Rt, Complex const &RtP, Complex const &RtPP, double const &St, double const &StP, double const &StPP);
void teukolskyIntegrandMinus2(Complex &integrandIn, Complex &integrandUp, int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &tR, double const &tTh, double const &rp, double const &thp, double const &phiR, double const &phiTh, double const &qr, double const &qth, Complex const &Rin, Complex const &RinP, Complex const &RinPP,  Complex const &Rup, Complex const &RupP, Complex const &RupPP, double const &St, double const &StP, double const &StPP);
void teukolskyIntegrandPlus2(Complex &integrandIn, Complex &integrandUp, int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &tR, double const &tTh, double const &rp, double const &thp, double const &phiR, double const &phiTh, double const &qr, double const &qth, Complex const &Rin, Complex const &RinP, Complex const &RinPP,  Complex const &Rup, Complex const &RupP, Complex const &RupPP, double const &St, double const &StP, double const &StPP);
void teukolskyIntegrandMinus2RadialTurningPoint(Complex &integrandIn, Complex &integrandUp, int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &tR, double const &tTh, double const &rp, double const &thp, double const &phiR, double const &phiTh, double const &qr, double const &qth, Complex const &Rin, Complex const &RinP, Complex const &RinPP,  Complex const &Rup, Complex const &RupP, Complex const &RupPP, double const &St, double const &StP, double const &StPP);
void teukolskyIntegrandPlus2RadialTurningPoint(Complex &integrandIn, Complex &integrandUp, int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &tR, double const &tTh, double const &rp, double const &thp, double const &phiR, double const &phiTh, double const &qr, double const &qth, Complex const &Rin, Complex const &RinP, Complex const &RinPP,  Complex const &Rup, Complex const &RupP, Complex const &RupPP, double const &St, double const &StP, double const &StPP);
void teukolskyIntegrandMinus2PolarTurningPoint(Complex &integrandIn, Complex &integrandUp, int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &tR, double const &tTh, double const &rp, double const &thp, double const &phiR, double const &phiTh, double const &qr, double const &qth, Complex const &Rin, Complex const &RinP, Complex const &RinPP,  Complex const &Rup, Complex const &RupP, Complex const &RupPP, double const &St, double const &StP, double const &StPP);
void teukolskyIntegrandPlus2PolarTurningPoint(Complex &integrandIn, Complex &integrandUp, int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &tR, double const &tTh, double const &rp, double const &thp, double const &phiR, double const &phiTh, double const &qr, double const &qth, Complex const &Rin, Complex const &RinP, Complex const &RinPP,  Complex const &Rup, Complex const &RupP, Complex const &RupPP, double const &St, double const &StP, double const &StPP);
void teukolskyIntegrandMinus2RadialPolarTurningPoint(Complex &integrandIn, Complex &integrandUp, int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &tR, double const &tTh, double const &rp, double const &thp, double const &phiR, double const &phiTh, double const &qr, double const &qth, Complex const &Rin, Complex const &RinP, Complex const &RinPP,  Complex const &Rup, Complex const &RupP, Complex const &RupPP, double const &St, double const &StP, double const &StPP);
void teukolskyIntegrandPlus2RadialPolarTurningPoint(Complex &integrandIn, Complex &integrandUp, int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &tR, double const &tTh, double const &rp, double const &thp, double const &phiR, double const &phiTh, double const &qr, double const &qth, Complex const &Rin, Complex const &RinP, Complex const &RinPP,  Complex const &Rup, Complex const &RupP, Complex const &RupPP, double const &St, double const &StP, double const &StPP);

int scalar_integrand_I1(Complex &integrand, int m, int n, double freq, double tR, double rp, double phiR, double qr, Complex Rt);
int scalar_integrand_I2(Complex &integrand, int m, int k, double freq, double tTh, double, double phiTh, double qth, double St);
int scalar_integrand_I3(Complex &integrand, int m, int n, double freq, double tR, double, double phiR, double qr, Complex Rt);
int scalar_integrand_I4(Complex &integrand, int m, int k, double freq, double tTh, double aCosThP, double phiTh, double qth, double St);
TeukolskyAmplitudes scalar_amplitude(int L, int m, int k, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConst, RadialTeukolsky &teuk, SpinWeightedHarmonic &swsh);
TeukolskyAmplitudes scalar_amplitude_equatorial(int L, int m, int k, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConst, RadialTeukolsky &teuk, SpinWeightedHarmonic &swsh);
TeukolskyAmplitudes scalar_amplitude_spherical(int L, int m, int k, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConst, RadialTeukolsky &teuk, SpinWeightedHarmonic &swsh);
TeukolskyAmplitudes scalar_amplitude_circular(int L, int m, int k, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConst, RadialTeukolsky &teuk, SpinWeightedHarmonic &swsh);
TeukolskyAmplitudes scalar_amplitude_generic(int L, int m, int k, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConst, RadialTeukolsky &teuk, SpinWeightedHarmonic &swsh);

// TeukolskyAmplitudes scalar_amplitude_circeq(int L, int m, GeodesicTrajectory& traj, GeodesicConstants &geoConst, ComplexDerivativesMatrixStruct Rin, ComplexDerivativesMatrixStruct Rup, DerivativesMatrix Slm);
// TeukolskyAmplitudes scalar_amplitude_ecceq(int L, int m, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConst, ComplexDerivativesMatrixStruct Rin, ComplexDerivativesMatrixStruct Rup, DerivativesMatrix Slm);
// TeukolskyAmplitudes scalar_amplitude_sphinc(int L, int m, int k, GeodesicTrajectory& traj, GeodesicConstants &geoConst, ComplexDerivativesMatrixStruct Rin, ComplexDerivativesMatrixStruct Rup, DerivativesMatrix Slm);
// TeukolskyAmplitudes scalar_amplitude(int L, int m, int k, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConst, ComplexDerivativesMatrixStruct Rin, ComplexDerivativesMatrixStruct Rup, DerivativesMatrix Slm);
// Complex scalar_integrand_ecceq(int L, int m, int n, GeodesicConstants &geoConstants, double tR, double rp, double phiR, double qr, Complex Rt, double St);
// Complex scalar_integrand_sphinc(int L, int m, int k, int n, GeodesicConstants &geoConstants, double tR, double tTh, double rp, double thp, double phiR, double phiTh, double qr, double qth, Complex Rt, double St);
Complex scalar_integrand_1(int L, int m, int k, int n, GeodesicConstants &geoConstants, double tR, double rp, double phiR, double qr, Complex Rt);
Complex scalar_integrand_2(int L, int m, int k, int n, GeodesicConstants &geoConstants, double tTh, double thp, double phiTh, double qth, double St);
Complex scalar_integrand_3(int L, int m, int k, int n, GeodesicConstants &geoConstants, double tR, double rp, double phiR, double qr, Complex Rt);
Complex scalar_integrand_4(int L, int m, int k, int n, GeodesicConstants &geoConstants, double tTh, double thp, double phiTh, double qth, double St);

Complex wronskian(double a, double rp, Complex Rin, Complex RinP, Complex Rup, Complex RupP);
Complex wronskian(int s, double a, double rp, Complex Rin, Complex RinP, Complex Rup, Complex RupP);
Complex scalar_wronskian(double a, double rp, Complex Rin, Complex RinP, Complex Rup, Complex RupP);

Complex A_nn_0(int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geo, double const &rp, double const &thp, double const &Slm, double const &SlmP, double const &SlmPP);
Complex A_nmbar_0(int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geo, double const &rp, double const &thp, double const &Slm, double const &SlmP, double const &SlmPP);
Complex A_nmbar_1(int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geo, double const &rp, double const &thp, double const &Slm, double const &SlmP, double const &SlmPP);
Complex A_mbarmbar_0(int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geo, double const &rp, double const &thp, double const &Slm, double const &SlmP, double const &SlmPP);
Complex A_mbarmbar_1(int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geo, double const &rp, double const &thp, double const &Slm, double const &SlmP, double const &SlmPP);
Complex A_mbarmbar_2(int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geo, double const &rp, double const &thp, double const &Slm, double const &SlmP, double const &SlmPP);
void A_coeffs(Complex &Ann0, Complex &Anmbar0, Complex &Ambarmbar0, Complex &Anmbar1, Complex &Ambarmbar1, Complex &Ambarmbar2, int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geo, double const &rp, double const &thp, double const &Slm, double const &SlmP, double const &SlmPP);
void A13_coeffs(Complex &All0, Complex &Alm0, Complex &Amm0, Complex &Alm1, Complex &Amm1, Complex &Amm2, int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geo, double const &rp, double const &thp, double const &Slm, double const &SlmP, double const &SlmPP);

void u_24_coeffs(Complex &u2m, Complex &u2p, Complex &u4m, Complex &u4p, GeodesicConstants &geoConstants, double const &rp, double const &thp);
void u_13_coeffs(Complex &u1m, Complex &u1p, Complex &u3m, Complex &u3p, GeodesicConstants &geoConstants, double const &rp, double const &thp);
void u_24_coeffs_RadialTurningPoint(Complex &u2, Complex &u4m, Complex &u4p, GeodesicConstants &geoConstants, double const &rp, double const &thp);
void u_13_coeffs_RadialTurningPoint(Complex &u1, Complex &u3m, Complex &u3p, GeodesicConstants &geoConstants, double const &rp, double const &thp);
void u_24_coeffs_PolarTurningPoint(Complex &u2m, Complex &u2p, Complex &u4, GeodesicConstants &geoConstants, double const &rp, double const &thp);
void u_13_coeffs_PolarTurningPoint(Complex &u1m, Complex &u1p, Complex &u3, GeodesicConstants &geoConstants, double const &rp, double const &thp);
void u_24_coeffs_RadialPolarTurningPoint(Complex &u2, Complex &u4, GeodesicConstants &geoConstants, double const &rp, double const &thp);
void u_13_coeffs_RadialPolarTurningPoint(Complex &u1, Complex &u3, GeodesicConstants &geoConstants, double const &rp, double const &thp);
double u_1(GeodesicConstants &geo, double const &rp, int const &sgnUr, double const &thp, int const &sgnUth);
Complex u_3(GeodesicConstants &geo, double const &rp, int const &sgnUr, double const &thp, int const &sgnUth);
double u_2(GeodesicConstants &geo, double const &rp, int const &sgnUr, double const &thp, int const &sgnUth);
Complex u_4(GeodesicConstants &geo, double const &rp, int const &sgnUr, double const &thp, int const &sgnUth);

double C_ll(GeodesicConstants &geo, double const &rp, int const &sgnUr, double const &thp, int const &sgnUth);
Complex C_lm(GeodesicConstants &geo, double const &rp, int const &sgnUr, double const &thp, int const &sgnUth);
Complex C_mm(GeodesicConstants &geo, double const &rp, int const &sgnUr, double const &thp, int const &sgnUth);
double C_nn(GeodesicConstants &geo, double const &rp, int const &sgnUr, double const &thp, int const &sgnUth);
Complex C_nmbar(GeodesicConstants &geo, double const &rp, int const &sgnUr, double const &thp, int const &sgnUth);
Complex C_mbarmbar(GeodesicConstants &geo, double const &rp, int const &sgnUr, double const &thp, int const &sgnUth);
void C_coeffs(Complex &Cnn, Complex &Cnmbar, Complex &Cmbarmbar, GeodesicConstants &geo, double const &rp, int const &sgnUr, double const &thp, int const &sgnUth);
void C13_coeffs(Complex &Cll, Complex &Clm, Complex &Cmm, GeodesicConstants &geo, double const &rp, int const &sgnUr, double const &thp, int const &sgnUth);


#endif
