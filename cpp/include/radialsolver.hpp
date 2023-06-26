// radialsolver.hpp

#ifndef RADIAL_HPP
#define RADIAL_HPP

#include <gsl/gsl_odeiv2.h>
#include <boost/numeric/odeint.hpp>
#include "mst.hpp"
#include "gsn_asymp.hpp"

enum SolutionMethod {AUTO, MST, ASYM, HBL, GSN, TEUK};

class RadialTeukolsky{
public:
	RadialTeukolsky(const double a, const int s, const int L, const int m, const double om, const Vector r);
	RadialTeukolsky(const double a, const int s, const int L, const int m, const double omega, const double lambda, const Vector r);
	~RadialTeukolsky();

	double getBlackHoleSpin();
	int getSpinWeight();
	int getSpheroidalModeNumber();
	int getAzimuthalModeNumber();
	double getModeFrequency();
	double getSpinWeightedSpheroidalEigenvalue();

	void generateRetardedBoundaryConditions(SolutionMethod method = AUTO);
	void generateRetardedBoundaryCondition(BoundaryCondition bc, SolutionMethod method = AUTO);
	void setBoundaryConditions(BoundaryCondition bc, Complex R, Complex Rp, double r);
	void generateSolutions(SolutionMethod method = AUTO, bool make_stable=true);
	void generateSolutions(BoundaryCondition bc, SolutionMethod method = AUTO, bool make_stable=true);
	int resampleSolutions(Vector radialSamples);

	void flipSpinWeight();

	Vector getRadialPoints();
	double getRadialPoints(int pos);
	double getBoundaryPoint(BoundaryCondition bc);
	Result getBoundarySolution(BoundaryCondition bc);
	Result getBoundaryDerivative(BoundaryCondition bc);

	ComplexVector getSolution(BoundaryCondition bc);
	ComplexVector getDerivative(BoundaryCondition bc);
	ComplexVector getSecondDerivative(BoundaryCondition bc);
	Complex getSolution(BoundaryCondition bc, int pos);
	Complex getDerivative(BoundaryCondition bc, int pos);
	Complex getSecondDerivative(BoundaryCondition bc, int pos);

protected:
	bool failCheck();
	bool failCheck(BoundaryCondition bc);
	bool validMethodDomain(SolutionMethod method);
	int ASYMsolveBoundary(BoundaryCondition bc);
	int MSTsolveBoundary(BoundaryCondition bc);
	void solveBoundaryPoint(BoundaryCondition bc);
	void AUTOsolve();
	void ASYMsolve();
	void MSTsolve();
	void STATICsolve();
	void AUTOsolve(BoundaryCondition bc);
	void ASYMsolve(BoundaryCondition bc);
	void MSTsolve(BoundaryCondition bc);
	void STATICsolve(BoundaryCondition bc);
	void HBLsolve(BoundaryCondition bc, bool make_stable);
	void TEUKsolve(BoundaryCondition bc, bool make_stable);
	void GSNsolve(BoundaryCondition bc);
	void SpinFlipsolve(BoundaryCondition bc, int (*func)(ComplexVector &, ComplexVector &, RadialTeukolsky&, const Vector &));

	double _a;
	int _s;
	int _L;
	int _m;
	double _omega;
	double _lambda;

	Vector _radialPoints;
	double _horizonBoundary;
	double _infinityBoundary;

	Result _horizonBoundarySolution;
	Result _horizonBoundaryDerivative;
	Result _infinityBoundarySolution;
	Result _infinityBoundaryDerivative;

	Complex _inTransmissionAmplitude;
	Complex _inIncidenceAmplitude;
	Complex _inReflectionAmplitude;
	Complex _upTransmissionAmplitude;
	Complex _upIncidenceAmplitude;
	Complex _upReflectionAmplitude;

	ComplexVector _inSolution;
	ComplexVector _inDerivative;
	ComplexVector _upSolution;
	ComplexVector _upDerivative;
};

typedef struct hbl_parameters_struct{
	double a;
	int s;
	int m;
	double om;
	double la;
	int H;
} hbl_parameters;

typedef struct rstar_parameters_struct{
	double rstar;
	double a;
} rstar_parameters;

typedef struct teukolsky_parameters_struct{
	int s;
	int l;
	int m;
	double a;
	double om;
	double la;
	Complex nu;
} teukolsky_parameters;

// teukolsky_parameters generate_teuk_parameters(int s, int l, int m, double a, double om);
// MstParameters teuk_to_mst_parameters(teukolsky_parameters params);
// hbl_parameters teuk_to_hbl_parameters(teukolsky_parameters params);

/////////////////////////
// Teukolsky functions //
/////////////////////////

//*************************************************************//
//       Teukolsky-Starobinsky Identities for s = \pm 2        //
//*************************************************************//

int flip_spin(int s);
double flip_eigenvalue(int s, double lambda);

double teukolsky_starobinsky_constant(int s, int m, double a, double omega, double lambda);
double teukolsky_starobinsky_constant(int m, double a, double omega, double lambda);
double teukolsky_starobinsky_constant_D(int m, double a, double omega, double lambda);
Complex teukolsky_starobinsky_complex_constant(int j, int m, double a, double omega, double lambda);

Complex teukolsky_starobinsky_amplitude(BoundaryCondition bc, int s, int m, double a, double omega, double lambda);
Complex teukolsky_starobinsky_minus_1_in(int m, double a, double omega, double);
Complex teukolsky_starobinsky_plus_1_in(int m, double a, double omega, double lambda);
Complex teukolsky_starobinsky_minus_1_up(int m, double a, double omega, double lambda);
Complex teukolsky_starobinsky_plus_1_up(int, double, double omega, double);
Complex teukolsky_starobinsky_minus_2_in(int m, double a, double omega, double);
Complex teukolsky_starobinsky_plus_2_in(int m, double a, double omega, double lambda);
Complex teukolsky_starobinsky_minus_2_up(int m, double a, double omega, double lambda);
Complex teukolsky_starobinsky_plus_2_up(int, double, double omega, double);

// void flip_spin_of_radial_teukolsky_TS(Complex &RinFlip, Complex &RupFlip, int m, double a, double omega, double lambda, double r, Complex Rin, Complex RinP, Complex Rup, Complex RupP);
void flip_spin_of_radial_teukolsky_TS(Complex &RinFlip, Complex &RinPFlip, BoundaryCondition bc, int s, int m, double a, double omega, double lambda, double r, Complex Rin, Complex RinP);
void flip_spin_of_radial_teukolsky_TS(ComplexVector &RinFlip, ComplexVector &RinPFlip, BoundaryCondition bc, int s, int m, double a, double omega, double lambda, Vector r, ComplexVector Rin, ComplexVector RinP);
// void flip_spin_of_radial_teukolsky_TS(Complex &RinFlip, Complex &RinPFlip, Complex &RupFlip, Complex &RupPFlip, int m, double a, double omega, double lambda, double r, Complex Rin, Complex RinP, Complex Rup, Complex RupP);
// void flip_spin_of_radial_teukolsky(ComplexVector &RinFlip, ComplexVector &RinPFlip, ComplexVector &RupFlip, ComplexVector &RupPFlip, int m, double a, double omega, double lambda, Vector rVec, ComplexVector RinVec, ComplexVector RinPVec, ComplexVector RupVec, ComplexVector RupPVec);
// void flip_spin_of_radial_teukolsky_static(ComplexVector &RinFlip, ComplexVector &RinPFlip, ComplexVector &RupFlip, ComplexVector &RupPFlip, double a, Vector rVec, ComplexVector RinVec, ComplexVector RinPVec, ComplexVector RupVec, ComplexVector RupPVec);


// int teuk_Rin(Complex Rin[], double r[], teukolsky_parameters params);
// int teuk_RinP(Complex RinP[], Complex Rin[], double r[], teukolsky_parameters params);

// Complex teuk_Rin(const double &r, teukolsky_parameters &teuk_params);
// Complex teuk_RinP(const Complex &Rin, const double &r, teukolsky_parameters &teuk_params);
// void test_teuk_Rin(void);
// void test_teuk_ode(void);
// void test_teuk_ode_2(void);

//*************************************************************//
// Different methods for solving the radial Teukolsky equation //
//*************************************************************//

// Use the semi-analytic expansions of Mano, Suzuki, and Takasugi to calculate radial Teukolsky solutions
int teuk_MST_series_boundary(Result &Rin, Result &RinP, Result &Rup, Result &RupP, RadialTeukolsky &teuk, const double &rIn, const double &rUp);
int teuk_in_MST_series_boundary(Result &R, Result &Rp, RadialTeukolsky &teuk, const double &r);
int teuk_up_MST_series_boundary(Result &R, Result &Rp, RadialTeukolsky &teuk, const double &r);

int teuk_MST_series(ComplexVector &Rin, ComplexVector &RinP, ComplexVector &Rup, ComplexVector &RupP, RadialTeukolsky &teuk, const Vector &r);
int teuk_in_MST_series(ComplexVector &R, ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r);
int teuk_in_MST_series(ComplexVector &R, RadialTeukolsky &teuk, const Vector &r);
int teuk_in_derivative_MST_series(ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r);
int teuk_up_MST_series(ComplexVector &R, ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r);
int teuk_up_MST_series(ComplexVector &R, RadialTeukolsky &teuk, const Vector &r);
int teuk_up_derivative_MST_series(ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r);

// Use Frobenius and asymptotic expansions near the singular points to calculate radial Teukolsky solutions
int teuk_ASYM_series(ComplexVector &Rin, ComplexVector &RinP, ComplexVector &Rup, ComplexVector &RupP, RadialTeukolsky &teuk, Vector &r);
int teuk_in_ASYM_series(ComplexVector &R, RadialTeukolsky &teuk, const Vector &r);
Result teuk_in_ASYM_series(RadialTeukolsky &teuk, const double &r);
int teuk_in_ASYM_series(Result &R, RadialTeukolsky &teuk, const double &r);
int teuk_in_derivative_ASYM_series(ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r);
Result teuk_in_derivative_ASYM_series(RadialTeukolsky &teuk, const double &r);
int teuk_in_derivative_ASYM_series(Result &R, RadialTeukolsky &teuk, const double &r);
int teuk_up_ASYM_series(ComplexVector &R, RadialTeukolsky &teuk, const Vector &r);
Result teuk_up_ASYM_series(RadialTeukolsky &teuk, const double &r);
int teuk_up_ASYM_series(Result &R, RadialTeukolsky &teuk, const double &r);
int teuk_up_derivative_ASYM_series(ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r);
Result teuk_up_derivative_ASYM_series(RadialTeukolsky &teuk, const double &r);
int teuk_up_derivative_ASYM_series(Result &R, RadialTeukolsky &teuk, const double &r);

// static solutions
int teuk_static(ComplexVector &Rin, ComplexVector &RinP, ComplexVector &Rup, ComplexVector &RupP, RadialTeukolsky &teuk, Vector &r);
Result teuk_in_static(RadialTeukolsky &teuk, const double &r);
int teuk_in_static(Result &R, RadialTeukolsky &teuk, const double &r);
int teuk_in_static(ComplexVector &R, RadialTeukolsky &teuk, const Vector &r);
Result teuk_in_derivative_static(RadialTeukolsky &teuk, const double &r);
int teuk_in_derivative_static(Result &R, RadialTeukolsky &teuk, const double &r);
int teuk_in_derivative_static(ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r);
Result teuk_up_static(RadialTeukolsky &teuk, const double &r);
int teuk_up_static(Result &R, RadialTeukolsky &teuk, const double &r);
int teuk_up_static(ComplexVector &R, RadialTeukolsky &teuk, const Vector &r);
Result teuk_up_derivative_static(RadialTeukolsky &teuk, const double &r);
int teuk_up_derivative_static(Result &R, RadialTeukolsky &teuk, const double &r);
int teuk_up_derivative_static(ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r);

Result teuk_in_static_solution(const double &a, const int &s, const int &l, const int &m, const double &r);
Result teuk_in_static_derivative(const double &a, const int &s, const int &l, const int &m, const double &r);
Result teuk_up_static_solution(const double &a, const int &s, const int &l, const int &m, const double &r);
Result teuk_up_static_derivative(const double &a, const int &s, const int &l, const int &m, const double &r);

// Get the second derivative of a Teukolsky solution using the Teukolsky equation
Complex teuk_secondDerivative(double a, int s, int m, double omega, double lambda, double r, Complex R, Complex Rp);

//*************************************************************//
// 						Teukolsky ODE solvers 				   //
//*************************************************************//
static const size_t state_dim = 4;
typedef std::array<double, state_dim> state_type;
typedef boost::numeric::ublas::vector<double> vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;
template <typename ODE_FUNC>
int teuk_integrate_boost(ComplexVector &psi, ComplexVector &dpsidr, ODE_FUNC sys, state_type psi0, const double r0, const Vector &r);
template <typename ODE_FUNC>
int teuk_integrate_boost(ComplexVector &psi, ComplexVector &dpsidr, ODE_FUNC sys, ODE_FUNC jac, state_type psi0, const double r0, const Vector &r);

int teuk_integrate_gsl(ComplexVector &psi, ComplexVector &dpsidr, int (*sys)(double, const double*, double*, void*), state_type psi0, const double r0, const Vector &r, void *params);
int teuk_integrate_gsl(ComplexVector &psi, ComplexVector &dpsidr, int (*sys)(double, const double*, double*, void*), int (*jac)(double, const double*, double*, double*, void*), state_type psi0, const double r0, const Vector &r, void *params);
int teuk_jac_null_gsl(double r, const double y[], double f[], void* params);

//*************************************************************//
// 				Teukolsky equation in BoyerLindquist 		   //
//*************************************************************//
// Obtain solutions by directly integrating the radial Teukolsky equation for the original Teukolsky function
class teuk_blc{
	hbl_parameters _params;

public:
    teuk_blc(hbl_parameters params);
    void operator()(const state_type &psi, state_type &dpsidr, const double r) const;
};

int teuk_blc_gsl(double r, const double y[], double dydr[], void* parameters);
double teuk_potential_FR(const double r, hbl_parameters params);
double teuk_potential_GR(const double r, hbl_parameters params);
double teuk_potential_GI(const double r, hbl_parameters params);

int teuk_in_TEUK_integrate(ComplexVector &R, ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r);
int teuk_up_TEUK_integrate(ComplexVector &R, ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r);

double teuk_delta(const double &r, const hbl_parameters &params);
double teuk_deltaP(const double &r, const hbl_parameters &params);
double teuk_deltaPP(const double &r, const hbl_parameters &params);
double teuk_K(const double &r, const hbl_parameters &params);
double teuk_KP(const double &r, const hbl_parameters &params);
double teuk_KPP(const double &r, const hbl_parameters &params);
double teuk_varpi(const double &r, const hbl_parameters &params);
double teuk_varpiP(const double &r, const hbl_parameters &params);

//*************************************************************//
// 			Generalized Sasaki-Nakamura transformation 		   //
//*************************************************************//

// Use the generalized Sasaki and Nakamura transformation to numerically integrate the radial Teukolsky ODE & obtain solution
class teuk_GSN{
	hbl_parameters _params;

public:
    teuk_GSN(hbl_parameters params);
    void operator()(const state_type &psi, state_type &dpsidr, const double r) const;
};
int teuk_GSN_gsl(double r, const double y[], double dydr[], void* parameters);

int teuk_in_GSN_integrate(ComplexVector &R, ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r);
int teuk_up_GSN_integrate(ComplexVector &R, ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r);

ComplexVector GSN_X_to_teuk_R(const ComplexVector X, const ComplexVector dX, const Vector &r, const hbl_parameters &params);

Complex teuk_R_to_GSN_X(const Complex R, const Complex Rp, const double &r, const hbl_parameters &params);
Complex teuk_RP_to_GSN_dX(const Complex R, const Complex Rp, const double &r, const hbl_parameters &params);
Complex GSN_X_to_teuk_R(const Complex X, const Complex dX, const double &r, const hbl_parameters &params);
Complex GSN_dX_to_teuk_RP(const Complex X, const Complex dX, const double &r, const hbl_parameters &params);

Complex teuk_R_to_GSN_chi(const Complex R, const Complex Rp, const double &r, const hbl_parameters &params);
Complex teuk_RP_to_GSN_dchi(const Complex R, const Complex Rp, const double &r, const hbl_parameters &params);
Complex GSN_chi_to_teuk_R(const Complex chi, const Complex dchi, const double &r, const hbl_parameters &params);
Complex GSN_dchi_to_teuk_RP(const Complex chi, const Complex dchi, const double &r, const hbl_parameters &params);

Complex GSN_X_to_GSN_chi(const Complex X, const Complex dX, const double &r, const hbl_parameters &params);
Complex GSN_dX_to_GSN_dchi(const Complex X, const Complex dX, const double &r, const hbl_parameters &params);
Complex GSN_chi_to_GSN_X(const Complex chi, const Complex dchi, const double &r, const hbl_parameters &params);
Complex GSN_chi_to_GSN_X(const Complex chi, const double &r, const hbl_parameters &params);
Complex GSN_dchi_to_GSN_dX(const Complex chi, const Complex dchi, const double &r, const hbl_parameters &params);
Result GSN_chi_to_GSN_X(const Result chi, const double &r, const hbl_parameters &params);
Result GSN_chi_to_GSN_X(const Result chi, const Result dchi, const double &r, const hbl_parameters &params);
Result GSN_dchi_to_GSN_dX(const Result chi, const Result dchi, const double &r, const hbl_parameters &params);

Complex GSN_eta(const double &r, const hbl_parameters &params);
Complex GSN_etaP(const double &r, const hbl_parameters &params);

Complex GSN_F_potential(const double &r, const hbl_parameters &params);
Complex GSN_U_potential(const double &r, const hbl_parameters &params);
Complex GSN_U_potential(const double &a, const int &s, const int &m, const double &omega, const double &lambda, const double &r);
Complex GSN_F_potential(const double &a, const int &s, const int &m, const double &omega, const double &lambda, const double &r);
Complex GSN_F_potential_func(const double &r, const hbl_parameters &params);
Complex GSN_U_potential_func(const double &r, const hbl_parameters &params);

Complex GSN_F_star_potential(const double &r, const double &delta, const double &deltaP, const hbl_parameters &params);
Complex GSN_U_star_potential(const double &r, const double &delta, const double &deltaP, const Complex &beta, const Complex &betaP, const Complex &betaPP, const hbl_parameters &params);
Complex GSN_G(const double &r, const double &delta, const double &deltaP, const hbl_parameters &params);
Complex GSN_GP(const double &r, const double &delta, const double &deltaP, const hbl_parameters &params);

Complex GSN_V(const double &r, const hbl_parameters &params);

Complex GSN_U1(const double &r, const double &delta, const double &deltaP, const Complex &beta, const Complex &betaP, const Complex &betaPP, const hbl_parameters &params);
Complex GSN_alpha(const double &r, const double &delta, const double &deltaP, const Complex &beta, const hbl_parameters &params);
Complex GSN_alphaP(const double &r, const double &delta, const double &deltaP, const Complex &beta, const Complex &betaP, const hbl_parameters &params);
Complex GSN_beta(const double &r, const double &delta, const double &deltaP, const hbl_parameters &params);
Complex GSN_betaP(const double &r, const double &delta, const double &deltaP, const hbl_parameters &params);
Complex GSN_betaPP(const double &r, const double &delta, const double &deltaP, const hbl_parameters &params);

Complex GSN_alpha(const double &r, const hbl_parameters &params);
Complex GSN_alphaP(const double &r, const hbl_parameters &params);
Complex GSN_beta(const double &r, const hbl_parameters &params);
Complex GSN_betaP(const double &r, const hbl_parameters &params);
Complex GSN_betaPP(const double &r, const hbl_parameters &params);

//*************************************************************//
// 				Hyperboloidal slicing transformation 		   //
//*************************************************************//

// Use hyperboloidal slicing to numerically integrate the radial Teukolsky ODE & obtain solutions

class teuk_hbl{
	hbl_parameters _params;

public:
    teuk_hbl(hbl_parameters params);
    void operator()(const state_type &psi, state_type &dpsidr, const double r) const;
};

struct teuk_hbl_implicit{
	hbl_parameters _params;

    teuk_hbl_implicit(hbl_parameters params);
    void operator()(const vector_type &psi, vector_type &dpsidr, const double r) const;
};

struct jacobi_hbl_implicit{
	hbl_parameters _params;

	jacobi_hbl_implicit(hbl_parameters params);
    void operator()( const vector_type &psi, matrix_type &J, const double r, vector_type &dfdr) const;
};

int teuk_in_HBL_integrate(ComplexVector &R, ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r);
int teuk_up_HBL_integrate(ComplexVector &R, ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r);

double PhiHBL(double r, double a);
double hHBL(double r, double a);

Complex teuk_R_to_hbl_Psi(Complex Rin, double r, hbl_parameters params);
Complex teuk_RP_to_hbl_dPsi(Complex RinP, Complex Rin, double r, hbl_parameters params);
Complex hbl_Psi_to_teuk_R(Complex Psi, double r, hbl_parameters params);
Complex hbl_dPsi_to_teuk_RP(Complex dPsi, Complex Psi, double r, hbl_parameters params);

int teuk_hbl_gsl(double r, const double y[], double f[], void* params);
int jacobian(double r, const double y[], double *dfdy, double dfdr[], void* params);

double potential_GR(double r, hbl_parameters params);
double potential_GI(double r, hbl_parameters params);
double potential_FR(double r, hbl_parameters params);
double potential_FI(double r, hbl_parameters params);
double potential_UR(double r, hbl_parameters params);
double potential_UI(double r, hbl_parameters params);

double potential_dGR(double r, hbl_parameters params);
double potential_dGI(double r, hbl_parameters params);
double potential_dUR(double r, hbl_parameters params);
double potential_dUI(double r, hbl_parameters params);

double rstar_to_r(double rstar, double a);
double r_to_rstar(double r, double a);
double r_to_rstar_root(double r, void* params);

//*************************************************************//
// 			 Series expansions near singular points 		   //
//*************************************************************//

Result teuk_up_asymptotic_infinity(const double &a, const int &s, const int &L, const int &m, const double &omega, const double &r);
Result teuk_up_asymptotic_infinity(const double &a, const int &s, const int &L, const int &m, const double &omega, const double &lambda, const double &r);
Result teuk_up_derivative_asymptotic_infinity(const double &a, const int &s, const int &L, const int &m, const double &omega, const double &r);
Result teuk_up_derivative_asymptotic_infinity(const double &a, const int &s, const int &L, const int &m, const double &omega, const double &lambda, const double &r);

Result teuk_in_asymptotic_horizon(const double &a, const int &s, const int &L, const int &m, const double &omega, const double &r);
Result teuk_in_asymptotic_horizon(const double &a, const int &s, const int &L, const int &m, const double &omega, const double &lambda, const double &r);
Result teuk_in_derivative_asymptotic_horizon(const double &a, const int &s, const int &L, const int &m, const double &omega, const double &r);
Result teuk_in_derivative_asymptotic_horizon(const double &a, const int &s, const int &L, const int &m, const double &omega, const double &lambda, const double &r);

Complex gsn_asymptotic_initial_sum(const double &a, const int &s, const int &m, const double &omega, const double &lambda, const double &r);
Complex gsn_asymptotic_initial_sum(const double &r, hbl_parameters params);
Result gsn_up_asymptotic_infinity(const double &a, const int &s, const int &L, const int &m, const double &omega, const double &r);
Result gsn_up_asymptotic_infinity(const double &r, hbl_parameters params);
Result gsn_up_asymptotic_infinity(const double &a, const int &s, const int &L, const int &m, const double &omega, const double &lambda, const double &r);
Result gsn_up_asymptotic_infinity(const double &a, const int &s, const int &m, const double &omega, const double &lambda, const double &r);
Result gsn_up_asymptotic_infinity_chi_series(const double &a, const int &s, const int &m, const double &omega, const double &lambda, const double &r);

Complex gsn_asymptotic_derivative_initial_sum(const double &r, hbl_parameters params);
Complex gsn_asymptotic_derivative_initial_sum(const double &a, const int &s, const int &m, const double &omega, const double &lambda, const double &r);
Complex gsn_asymptotic_derivative_initial_sum_series(const double &r, hbl_parameters params);
Complex gsn_asymptotic_derivative_initial_sum_series(const double &a, const int &s, const int &m, const double &omega, const double &lambda, const double &r);
Result gsn_up_derivative_asymptotic_infinity(const double &a, const int &s, const int &L, const int &m, const double &omega, const double &r);
Result gsn_up_derivative_asymptotic_infinity(const double &r, hbl_parameters params);
Result gsn_up_derivative_asymptotic_infinity(const double &a, const int &s, const int &m, const double &omega, const double &lambda, const double &r);
Result gsn_up_derivative_asymptotic_infinity_chi_series(const double &a, const int &s, const int &m, const double &omega, const double &lambda, const double &r);

#endif
