// hertz.hpp

#ifndef HERTZ_HPP
#define HERTZ_HPP

#include "teukolsky.hpp"

enum Gauge {ORG, IRG, SAAB0, SAAB4, ASAAB0, ASAAB4, SAAB, ASAAB};

class HertzMode{
public:
	HertzMode(TeukolskyMode& teuk, Gauge gauge = ORG);
	~HertzMode();

	int generateSolutions();

	Gauge getGauge();
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

	Vector getScalarCouplingCoefficient();
	double getScalarCouplingCoefficient(int l);
	int getMinScalarCouplingModeNumber();
	int getMaxScalarCouplingModeNumber();

	Vector getRadialPoints();
	ComplexVector getHomogeneousRadialSolution(BoundaryCondition bc);
	ComplexVector getHomogeneousRadialDerivative(BoundaryCondition bc);
	Complex getHertzAmplitude(BoundaryCondition bc);
	ComplexVector getRadialSolution(BoundaryCondition bc);
	ComplexVector getRadialDerivative(BoundaryCondition bc);

	Vector getPolarPoints();
	Vector getPolarSolution();
	Vector getPolarDerivative();

	double getRadialPoints(int pos);
	Complex getHomogeneousRadialSolution(BoundaryCondition bc, int pos);
	Complex getHomogeneousRadialDerivative(BoundaryCondition bc, int pos);
	Complex getHomogeneousRadialSecondDerivative(BoundaryCondition bc, int pos);
	Complex getHomogeneousRadialSolution(BoundaryCondition bc, int dr, int pos);
	Complex getRadialSolution(BoundaryCondition bc, int pos);
	Complex getRadialDerivative(BoundaryCondition bc, int pos);
	Complex getRadialSecondDerivative(BoundaryCondition bc, int pos);
	Complex getRadialSolution(BoundaryCondition bc, int dr, int pos);

	double getPolarPoints(int pos);
	double getPolarSolution(int pos);
	double getPolarDerivative(int pos);
	double getPolarSecondDerivative(int pos);

private:
	Gauge _gauge;

	int _s;
	int _L;
	int _m;
	int _k;
	int _n;
	int _sampleSize;

	double _a;
	double _omega;
	double _lambda;

	double _lmin;
	Vector _coupling;
	Vector _scalarCoupling;

	Vector _theta;
	Vector _Slm;
	Vector _SlmP;
	Vector _SlmPP;

	Vector _r;
	ComplexVector _Rin;
	ComplexVector _RinP;
	ComplexVector _RinPP;
	ComplexVector _Rup;
	ComplexVector _RupP;
	ComplexVector _RupPP;

	Complex _PsilmIn;
	Complex _PsilmUp;

	void flipSpinWeight(); //helper function taken from TeukolskyMode
};

void test_hertz_mode(int j, int m, int k, int n, GeodesicSource& geo);

void teukolsky_to_hertz_ORG(Complex &Psi, Complex Zteuk, int L, int m, int k, double a, double omega, double lambdaCH);
void teukolsky_to_hertz_IRG(Complex &Psi, Complex Zteuk, int L, int m, int k, double a, double omega, double lambdaCH);
void teukolsky_to_hertz_SAAB(Complex &Psi, Complex Zteuk, int L, int m, int k, double a, double omega, double lambdaCH);
void teukolsky_to_hertz_ASAAB(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambdaCH);

void teukolsky_to_hertz_ORG(Complex &PsiIn, Complex &PsiUp, Complex ZteukIn, Complex ZteukUp, int L, int m, int k, double a, double omega, double lambdaCH);
void teukolsky_to_hertz_IRG(Complex &PsiIn, Complex &PsiUp, Complex ZteukIn, Complex ZteukUp, int L, int m, int k, double a, double omega, double lambdaCH);
void teukolsky_to_hertz_SAAB(Complex &PsiIn, Complex &PsiUp, Complex ZteukIn, Complex ZteukUp, int L, int m, int k, double a, double omega, double lambdaCH);
void teukolsky_to_hertz_ASAAB(Complex &PsiIn, Complex &PsiUp, Complex ZteukIn, Complex ZteukUp, int L, int m, double a, double omega, double lambdaCH);

// void teukolsky_to_hertz_amplitude_in(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda);
// void teukolsky_to_hertz_amplitude_up(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda);
// void teukolsky_to_hertz_amplitude_IRG_in(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda);
// void teukolsky_to_hertz_amplitude_IRG_up(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda);
// void teukolsky_to_hertz_amplitude_SAAB0_in(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda);
// void teukolsky_to_hertz_amplitude_SAAB0_up(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda);
// void teukolsky_to_hertz_amplitude_SAAB4_in(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda);
// void teukolsky_to_hertz_amplitude_SAAB4_up(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda);
// void teukolsky_to_hertz_amplitude_ASAAB0_in(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda);
// void teukolsky_to_hertz_amplitude_ASAAB0_up(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda);
// void teukolsky_to_hertz_amplitude_ASAAB4_in(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda);
// void teukolsky_to_hertz_amplitude_ASAAB4_up(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda);

void flip_spin_of_spheroidal_eigenvalue(double &lambdaFlip, int s, double lambda);

void generate_scalar_spherical_spheroidal_coupling(Vector &Bljm, int s, int lmin, int m, Vector bljm);
void generate_radial_second_derivative(ComplexVector &R2, int s, int m, double a, double omega, double lambda, Vector r, ComplexVector R0, ComplexVector R1);

#endif
