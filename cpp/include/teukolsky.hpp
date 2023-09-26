// teukolsky.hpp

#ifndef TEUKOLSKY_HPP
#define TEUKOLSKY_HPP

#include "radialsolver.hpp"
#include "sourceintegration.hpp"
#include "geo.hpp"

#define SAMPLE_SIZE_INIT 256

class TeukolskyMode{
public:
	TeukolskyMode(int L, int m, int k, int n, GeodesicSource& geo);
	TeukolskyMode(int s, int L, int m, int k, int n, GeodesicSource& geo);
	TeukolskyMode(int s, int L, int m, int k, int n, double a, double omega, Vector theta, Vector r);
	~TeukolskyMode();

	int generateSolutions(GeodesicSource& geo, SolutionMethod method = AUTO, int samplesize = SAMPLE_SIZE_INIT);
	int generateSolutions(double omega, GeodesicTrajectory &traj, GeodesicConstants &geoConst, Vector r, Vector theta, SolutionMethod method, int samplesize);
	int generateSolutions(SpinWeightedHarmonic& swsh, RadialTeukolsky& teuk, GeodesicTrajectory& traj, GeodesicConstants &geoConst);

	// int extendSolutions(Vector theta, Vector r);
	int flipSpinWeightAndFrequency();
	int flipSpinWeight();

	int getSpinWeight();
	int getSpheroidalModeNumber();
	int getAzimuthalModeNumber();
	int getPolarModeNumber();
	int getRadialModeNumber();
	int getSampleSize();
	double getBlackHoleSpin();
	double getFrequency();
	double getHorizonFrequency();
	double getEigenvalue();
	Vector getCouplingCoefficient();
	double getCouplingCoefficient(int l);
	int getMinCouplingModeNumber();
	int getMaxCouplingModeNumber();

	Vector getRadialPoints();
	ComplexVector getHomogeneousRadialSolution(BoundaryCondition bc);
	ComplexVector getHomogeneousRadialDerivative(BoundaryCondition bc);
	ComplexVector getHomogeneousSecondRadialDerivative(BoundaryCondition bc);
	Complex getTeukolskyAmplitude(BoundaryCondition bc);
	double getTeukolskyAmplitudePrecision(BoundaryCondition bc);
	ComplexVector getRadialSolution(BoundaryCondition bc);
	ComplexVector getRadialDerivative(BoundaryCondition bc);

	Vector getPolarPoints();
	Vector getPolarSolution();
	Vector getPolarDerivative();

	int getRadialSampleNumber();
	double getRadialPoints(int pos);
	Complex getHomogeneousRadialSolution(BoundaryCondition bc, int pos);
	Complex getHomogeneousRadialDerivative(BoundaryCondition bc, int pos);
	Complex getHomogeneousSecondRadialDerivative(BoundaryCondition bc, int pos);
	Complex getRadialSolution(BoundaryCondition bc, int pos);
	Complex getRadialDerivative(BoundaryCondition bc, int pos);

	int getPolarSampleNumber();
	double getPolarPoints(int pos);
	double getPolarSolution(int pos);
	double getPolarDerivative(int pos);
	double getPolarSecondDerivative(int pos);

private:
	int _s;
	int _L;
	int _m;
	int _k;
	int _n;
	int _sampleSize;

	double _a;
	double _omega;
	double _lambda;
	Vector _coupling;

	Vector _theta;
	Vector _Slm;
	Vector _SlmP;

	Vector _r;
	ComplexVector _Rin;
	ComplexVector _RinP;
	ComplexVector _Rup;
	ComplexVector _RupP;

	Complex _ZlmIn;
	Complex _ZlmUp;
	double _ZlmInPrecision;
	double _ZlmUpPrecision;
};

void flip_spin_of_coupling_coefficients(Vector &bslmo, int L, int m);
void flip_spin_of_spheroidal_harmonic(Vector &SlmFlip, Vector &SlmPFlip, int l, int m);
void flip_spin_of_radial_teukolsky_amplitude_TS(Complex &Zlm, BoundaryCondition bc, int s, int j, int m, int k, double a, double omega, double lambdaCH);
double test_teukolsky_solutions(int spin, Complex R0, Complex R1, Complex R2, double r, double a, double m, double omega, double lambda);
int minimum_radial_harmonic(int m, int k, GeodesicSource& geo);
int minimum_polar_harmonic_circ(int m, GeodesicSource& geo);

// class SelfForceModeCoupling{
// public:
// 	SelfForceModeCoupling(int n, int lmax, int jmax, int m);
// 	~SelfForceModeCoupling();
//
// 	void generateCoupling();
// 	double getModeCoupling(int l, int j, int jth);
//
// private:
// 	int _n;
// 	int _lmax;
// 	int _jmax;
// 	int _m;
//
// 	Vector _thp;
// 	RealTensor _couplingCoeffs;
// }

#endif
