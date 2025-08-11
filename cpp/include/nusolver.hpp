//////////////////////////////////////////////////////////////
//
//	nusolver.hpp
//
//////////////////////////////////////////////////////////////

#ifndef NUSOLVER_HPP
#define NUSOLVER_HPP

#include <cstdlib> // for std::exit
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include "monodromy.hpp"
#include "cf.hpp"
#include "swsh.hpp"

/////////////
// Classes //
/////////////

// SpheroidalModeParameters class

class SpheroidalModeParameters{
public:
	SpheroidalModeParameters(double a, int s, int L, int m, double omega);
	SpheroidalModeParameters(double a, int s, int L, int m, double omega, double lambda);
	SpheroidalModeParameters(double a, int s, int L, int m, int k, int n, double omega);
	SpheroidalModeParameters(double a, int s, int L, int m, int k, int n, double omega, double lambda);
	~SpheroidalModeParameters();
	
	double getBlackHoleSpin() const;
	int getSpinWeight() const;
	int getSpinWeightedSpheroidalModeNumber() const;
	int getAzimuthalModeNumber() const;
	int getPolarModeNumber() const;
	int getRadialModeNumber() const;
	double getModeFrequency() const;
	double getSpinWeightedSpheroidalEigenvalue() const;
	
	void setBlackHoleSpin(double a);
	void setSpinWeight(int s);
	void setModeNumbers(int L, int m);
	void setModeNumbers(int L, int m, int k, int n);
	void setSpinWeightedSpheroidalModeNumber(int L);
	void setAzimuthalModeNumber(int m);
	void setPolarModeNumber(int k);
	void setRadialModeNumber(int n);
	void setModeFrequency(double omega);
	void setSpinWeightedSpheroidalEigenvalue(double lambda);

protected:
	double _a;
	int _s;
	int _L;
	int _m;
	int _k;
	int _n;
	double _omega;
	double _lambda;
};

// MstParameters class

class MstParameters: public SpheroidalModeParameters{
public:
	MstParameters(double q, int s, int L, int m, double epsilon);
	MstParameters(double q, int s, int L, int m, double epsilon, double lambda);
	MstParameters(double q, int s, int L, int m, double epsilon, double lambda, Complex nu);
	~MstParameters();
	
	double getMstQ() const;
	double getMstEpsilon() const;
	double getMstKappa() const;
	double getMstTau() const;
	double getMstX(double r) const;
	double getMstDXDR() const;
	double getMstDXDR(double r) const;
	Complex getRenormalizedAngularMomentum() const;
	
	void setMstQ(double q);
	void setMstEpsilon(double epsilon);
	void setRenormalizedAngularMomentum(Complex nu);
	
protected:
	Complex _nu;
};

typedef struct LR_parameters_struct{
	MstParameters mstParams;
	int n0;
} LR_parameters;

//////////////////////
// MST coefficients //
//////////////////////

double epsMST(double om);
double kappaMST(double a);
double tauMST(int m, double a, double om);
double xMST(double a, double r);
double xPMST(double a, double r);

Complex alphaMST(int n, MstParameters params);
Complex betaMST(int n, MstParameters params);
Complex gammaMST(int n, MstParameters params);

////////////////////////////////////////////
// Rn and Ln continued fraction functions //
////////////////////////////////////////////

Complex Rn_a_coeff(int n, void* p);
Complex Rn_b_coeff(int n, void* p);

Complex alphaRn_cf(int n, const MstParameters &params);
Complex Rn_cf(int n, const MstParameters &params);

Complex Ln_a_coeff(int n, void* p);
Complex Ln_b_coeff(int n, void* p);

Complex gammaLn_cf(int n, const MstParameters &params);
Complex gammaLn_recursion(int n, const MstParameters &params);
Complex Ln_cf(int n, const MstParameters &params);

//////////////////////////////
// Characteristic equations //
//////////////////////////////

Complex nu_eqn_136(const MstParameters& params);

// @brief Characteristic equation when renormalized angular momentum $\nu \in \mathbb{Z}$.
// @parameter nu The real part of $\nu$.
// @parameter p Pointer to additional parameters in the characteristic equation.
double nu_eqn_136_realnu(double nu, void* p);

// @brief Characteristic equation when renormalized angular momentum $\nu = -0.5 + \nuI*i$.
// @parameter nuI The imaginary part of $\nu$.
// @parameter p Pointer to additional parameters in the characteristic equation.
// @return The complex residual of the characteristic equation.
Complex nu_eqn_136_complexnu_half_full(double nuI, void* p);

// @brief Characteristic equation when renormalized angular momentum $\nu = -0.5 + \nuI*i$
// @parameter nuI The imaginary part of $\nu$.
// @parameter p Pointer to additional parameters in the characteristic equation.
// @return The real part of the residual of the characteristic equation.
double nu_eqn_136_complexnu_half(double nuI, void* p);

// @brief Characteristic equation when renormalized angular momentum $\nu = -0.5 + \nuI*i$.
// @parameter nuI The imaginary part of $\nu$.
// @parameter p Pointer to additional parameters in the characteristic equation.
// @return The imaginary part of the residual of the characteristic equation.
double nu_eqn_136_complexnu_half_imag(double nuI, void* p);

Complex nu_eqn_136_complexnu_full(double nu, void* p);

// @brief Characteristic equation when renormalized angular momentum $\nu = -0.5 + \nuI*i$
// @parameter nuI The imaginary part of $\nu$.
// @parameter p Pointer to additional parameters in the characteristic equation.
// @return The real part of the residual of the characteristic equation.
double nu_eqn_136_complexnu(double nuI, void* p);

// @brief Characteristic equation when renormalized angular momentum $\nu = 1.0 + \nuI*i$.
// @parameter nuI The imaginary part of $\nu$.
// @parameter p Pointer to additional parameters in the characteristic equation.
// @return The imaginary part of the residual of the characteristic equation.
double nu_eqn_136_complexnu_imag(double nuI, void* p);

double nu_eqn_133_realnu(double nu, void* p);
double nu_eqn_133_complexnu_half(double nu, void* p);
double nu_eqn_133_complexnu(double nu, void* p);

///////////////////////////
// Root finder functions //
///////////////////////////

Complex nu_solver(double q, int s, int l, int m, double epsilon);
int nu_solver(MstParameters &params);
int nu_solver_guess(MstParameters &params);

int nu_solver_noguess_rootfinder(gsl_function F, const double &x_lo, const double &x_hi, MstParameters &params);
int nu_solver_real_noguess(MstParameters &params);
int nu_solver_real_noguess(int stepNum, MstParameters &params);
int realnu_interval_test(double x_lo, double x_hi, MstParameters& params);
int realnu_interval_search(int stepNum, double& x_lo, double& x_hi, MstParameters& params);
int realnu_interval_search(double& x_lo, double& x_hi, MstParameters& params);
int realnu_root_test(Complex& nu, double& nuError, MstParameters& params);

int nu_solver_complex_noguess(MstParameters &params);
int nu_solver_noguess(MstParameters &params);

//*************************************************************************************
// Low frequency expansions
//*************************************************************************************

double nu_solver_low_freq(double q, int s, int l, int m, double eps);

double low_frequency_max_epsilon(double q, int l, int m);
double max_epsilon_power_law(double q, int l, int m);

double nu_solver_low_freq_s2(double q, int l, int m, double eps);
double nu_solver_low_freq_s0(double q, int l, int m, double eps);
double nu_solver_low_freq_s0_l0(double q, double eps);
double nu_solver_low_freq_s0_l1(double q, int m, double eps);
double nu_solver_low_freq_s0_l2(double q, int m, double eps);
double nu_solver_low_freq_s0_l3(double q, int m, double eps);
double nu_solver_low_freq_s0_l4(double q, int m, double eps);
double nu_solver_low_freq_s0_l5(double q, int m, double eps);
double nu_solver_low_freq_s0_l6(double q, int m, double eps);
double nu_solver_low_freq_s0_l(double q, int l, int m, double eps);

double nu_solver_low_freq_s2_l2(double q, int m, double eps);
double nu_solver_low_freq_s2_l3(double q, int m, double eps);
double nu_solver_low_freq_s2_l4(double q, int m, double eps);
double nu_solver_low_freq_s2_l5(double q, int m, double eps);
double nu_solver_low_freq_s2_l6(double q, int m, double eps);
double nu_solver_low_freq_s2_l(double q, int l, int m, double eps);
double gen_l_last_term_s2(double q, int l, int m);

#endif
