// swsh.hpp

#ifndef SWSH_HPP
#define SWSH_HPP

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_legendre.h>
#include "specialfunc.hpp"

// SWSH

class SpinWeightedHarmonic{
public:
	SpinWeightedHarmonic(int s, int L, int m, double gamma, const Vector& theta);
	~SpinWeightedHarmonic();

	int getSpinWeight();
	int getSpheroidalModeNumber();
	int getAzimuthalModeNumber();
	double getSpheroidicity();
	double getEigenvalue();
	Vector getCouplingCoefficient();
	double getCouplingCoefficient(int l);
	int getMinCouplingModeNumber();
	int getMaxCouplingModeNumber();

	int generateSolutionsAndDerivatives();
	int generateCouplingCoefficients();
	int generateSolutions();
	int generateDerivatives();

	Vector getArguments();
	Vector getSolution();
	Vector getDerivative();
	Vector getSecondDerivative();

	double getArguments(int pos);
	double getSolution(int pos);
	double getDerivative(int pos);
	double getSecondDerivative(int pos);

private:
	int _s;
	int _L;
	int _m;
	double _gamma;
	double _lambda;
	Vector _bcoupling;

	Vector _theta;
	Vector _Slm;
	Vector _SlmP;
};

typedef struct coupling_convergence_struct{
	int jmax;
	int lmax;
	int testIndex;
	double jmax_err;
	double lmax_err;
	double test_err;
} coupling_converge;

enum coupling_test {SUCCESS, FAIL, STALL};

/////////////////////
// Spectral matrix //
/////////////////////

// Matrix coefficients
double k1(const int &s, const int &l, const int &j, const int &m);
double k2(const int &s, const int &l, const int &j, const int &m);

double akm2(const int &s, const int &l, const int &m, const double &g);
double akm1(const int &s, const int &l, const int &m, const double &g);
double akp0(const int &s, const int &l, const int &m, const double &g);
double akp1(const int &s, const int &l, const int &m, const double &g);
double akp2(const int &s, const int &l, const int &m, const double &g);

double spectral_weight(const double &g);

// Spectral matrix
int spectral_matrix(const int &s, const int &lmin, const int &m, const double &g, gsl_matrix* mat);
int spectral_matrix_sparse_init(const int &s, const int &lmin, const int &m, const double &g, gsl_spmatrix* mat);
int spectral_matrix_sparse(const int &s, const int &lmin, const int &m, const double &g, gsl_spmatrix* mat, const size_t &nmax);

// Solve eigensystem of spectral matrix
int spectral_solver(const int &s, const int &l, const int &m, const double &g, double& la, Vector& bvec);

double spectral_solver_n(const int &s, const int &l, const int &m, const double &g, const unsigned int &nmax);
int spectral_solver_n(const int &s, const int &l, const int &m, const double &g, gsl_vector* la);
int spectral_solver_n(const int &s, const int &l, const int &m, const double &g, gsl_vector* la, gsl_matrix* bmat);
int spectral_solver_n(const int &s, const int &l, const int &m, const double &g, gsl_vector* la, gsl_matrix* bmat, gsl_spmatrix* mat);

coupling_test spherical_spheroidal_coupling_convergence_test(const int &s, const int &l, const int &m, const double &g, gsl_matrix* bmat, gsl_matrix* bmat2, coupling_converge &b_data);

////////////////////
// Spin Couplings //
////////////////////

// Coupling between scalar and spin-weighted harmonics
double Asljm(const int &s, const int &l, const int &j, const int &m);
double dAsljm(const int &s, const int &l, const int &j, const int &m);

// Clebsch-Gordan coefficients
double clebsch(const int &j1, const int &j2, const int &j, const int &m1, const int &m2, const int &m);

// Wigner 3J symbol
double w3j(const int &j1, const int &j2, const int &j, const int &m1, const int &m2, const int &m);

//////////////////////////////////////////////
// Spin-Weighted Spherical Harmonics (SWSH) //
//////////////////////////////////////////////

// SWSHs
Complex Sslm(const int &s, const int &l, const int &m, const double &g, const double &th, const double &ph);
double Sslm(const int &s, const int &l, const int &m, const double &g, const double &th);

double Sslm(const int &s, const int &l, const int &m, const double &g, const Vector& bvec, const double &th);
double Sslm_derivative(const int &s, const int &l, const int &m, const double &g, const Vector& bvec, const double &th);
double Sslm_secondDerivative(const int &s, const int &l, const int &m, const double &g, const double& lambda, const double &th, const double &Slm, const double &SlmP);

// SWSH Eigenvalues
double swsh_eigenvalue(const int &s, const int &l, const int &m, const double &g);

///////////////////////////////////////
// Spin-Weighted Spherical Harmonics //
///////////////////////////////////////

Complex Yslm(const int &s, const int &l, const int &m, const double &th, const double &ph);
double Yslm(const int &s, const int &l, const int &m, const double &th);
double Yslm_derivative(const int &s, const int &l, const int &m, const double &th);

////////////////////
// Test functions //
////////////////////

void test_swsh_eigenvalue(void);
void test_swsh_class();

#endif
