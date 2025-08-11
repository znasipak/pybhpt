// nusolver.c

#include "nusolver.hpp"

#define eqn_136_EPS 1.e-13
#define CF_ITER_MAX 200

/////////////
// Classes //
/////////////

// SpheroidalModeParameters class

SpheroidalModeParameters::SpheroidalModeParameters(double a, int s, int L, int m, double omega):
	_a(a), _s(s), _L(L), _m(m), _k(0), _n(0), _omega(omega)
{
	_lambda = swsh_eigenvalue(s, L, m, a*omega);
}
SpheroidalModeParameters::SpheroidalModeParameters(double a, int s, int L, int m, double omega, double lambda):
	_a(a), _s(s), _L(L), _m(m), _k(0), _n(0), _omega(omega),  _lambda(lambda) {}
SpheroidalModeParameters::SpheroidalModeParameters(double a, int s, int L, int m, int k, int n, double omega):
	_a(a), _s(s), _L(L), _m(m), _k(k), _n(n), _omega(omega)
{
	_lambda = swsh_eigenvalue(s, L, m, a*omega);
}
SpheroidalModeParameters::SpheroidalModeParameters(double a, int s, int L, int m, int k, int n, double omega, double lambda):
	_a(a), _s(s), _L(L), _m(m), _k(k), _n(n), _omega(omega),  _lambda(lambda) {}
SpheroidalModeParameters::~SpheroidalModeParameters(){}

double SpheroidalModeParameters::getBlackHoleSpin() const { return _a; }
int SpheroidalModeParameters::getSpinWeight() const { return _s; }
int SpheroidalModeParameters::getSpinWeightedSpheroidalModeNumber() const { return _L; }
int SpheroidalModeParameters::getAzimuthalModeNumber() const { return _m; }
int SpheroidalModeParameters::getPolarModeNumber() const { return _k; }
int SpheroidalModeParameters::getRadialModeNumber() const { return _n; }
double SpheroidalModeParameters::getModeFrequency() const { return _omega; }
double SpheroidalModeParameters::getSpinWeightedSpheroidalEigenvalue() const { return _lambda; }

void SpheroidalModeParameters::setBlackHoleSpin(double a){ _a = a; }
void SpheroidalModeParameters::setSpinWeight(int s){ _s = s; }
void SpheroidalModeParameters::setModeNumbers(int L, int m){ _L = L; _m = m; }
void SpheroidalModeParameters::setModeNumbers(int L, int m, int k, int n){ _L = L; _m = m; _k = k; _n = n; }
void SpheroidalModeParameters::setSpinWeightedSpheroidalModeNumber(int L){ _L = L; }
void SpheroidalModeParameters::setAzimuthalModeNumber(int m){ _m = m; }
void SpheroidalModeParameters::setPolarModeNumber(int k){ _k = k; }
void SpheroidalModeParameters::setRadialModeNumber(int n){ _n = n; }
void SpheroidalModeParameters::setModeFrequency(double omega){ _omega = omega; }
void SpheroidalModeParameters::setSpinWeightedSpheroidalEigenvalue(double lambda){ _lambda = lambda; }

// MstParameters class

MstParameters::MstParameters(double q, int s, int L, int m, double epsilon):
	SpheroidalModeParameters(q, s, L, m, epsilon/2.), _nu(0.) {}
MstParameters::MstParameters(double q, int s, int L, int m, double epsilon, double lambda):
	SpheroidalModeParameters(q, s, L, m, epsilon/2., lambda), _nu(0.) {}
MstParameters::MstParameters(double q, int s, int L, int m, double epsilon, double lambda, Complex nu):
	SpheroidalModeParameters(q, s, L, m, epsilon/2., lambda), _nu(nu) {}
MstParameters::~MstParameters() {}

double MstParameters::getMstQ() const { return _a; }
double MstParameters::getMstEpsilon() const { return 2.*_omega; }
double MstParameters::getMstKappa() const { return sqrt(1. - _a*_a); }
double MstParameters::getMstTau() const { return (2.*_omega - _m*_a)/sqrt(1. - _a*_a); }
double MstParameters::getMstX(double r) const {
	return (1. + getMstKappa() - r)/(2.*getMstKappa());
}
double MstParameters::getMstDXDR(double) const {
	return getMstDXDR();
}
double MstParameters::getMstDXDR() const {
	return -1./(2.*getMstKappa());
}
Complex MstParameters::getRenormalizedAngularMomentum() const { return _nu; }

void MstParameters::setMstQ(double q){ _a = q; }
void MstParameters::setMstEpsilon(double epsilon){ _omega = epsilon/2.; }
void MstParameters::setRenormalizedAngularMomentum(Complex nu){ _nu = nu; }

//////////////////////
// MST coefficients //
//////////////////////

double epsMST(double om){
	return 2*om;
}

double kappaMST(double a){
	return sqrt(1-a*a);
}

double tauMST(int m, double a, double om){
	return (epsMST(om) - m*a)/kappaMST(a);
}

double xMST(double a, double r){
	double kappa = kappaMST(a);
	double rp = 1. + kappa;

	return (rp - r)/(2.*kappa);
}

double xPMST(double a, double){
	double kappa = kappaMST(a);

	return  - 1./(2.*kappa);
}

Complex alphaMST(int n, MstParameters params){
	Complex s = Complex(params.getSpinWeight());
	Complex eps = params.getMstEpsilon(), kappa = params.getMstKappa(), tau = params.getMstTau();
	Complex npnu = Complex(n) + params.getRenormalizedAngularMomentum();

	return I*eps*kappa*(npnu + 1. + s + I*eps)*(npnu + 1. + s - I*eps)
		*(npnu + 1. + I*tau)/(npnu + 1.)/(2.*npnu + 3.);
}

Complex betaMST(int n, MstParameters params){
	Complex s = Complex(params.getSpinWeight());
	Complex eps = params.getMstEpsilon(), kappa = params.getMstKappa(),
		tau = params.getMstTau(), la = params.getSpinWeightedSpheroidalEigenvalue();
	Complex npnu = Complex(n) + params.getRenormalizedAngularMomentum();

	return (-la - s*(s + 1.) + npnu*(npnu + 1.) + eps*eps + eps*kappa*tau +
		eps*kappa*tau*(s*s + eps*eps)/npnu/(npnu + 1.));
}

Complex gammaMST(int n, MstParameters params){
	Complex s = Complex(params.getSpinWeight());
	Complex eps = params.getMstEpsilon(), kappa = params.getMstKappa(), tau = params.getMstTau();
	Complex npnu = Complex(n) + params.getRenormalizedAngularMomentum();

	return -I*eps*kappa*(npnu - s + I*eps)*(npnu - s - I*eps)
		*(npnu - I*tau)/npnu/(2.*npnu - 1.);
}

////////////////////////////////////////////
// Rn and Ln continued fraction functions //
////////////////////////////////////////////

Complex Rn_a_coeff(int n, void* p){
	LR_parameters *params = (LR_parameters *)p;
	int n0 = params->n0;

	return -alphaMST(n + n0 - 1, params->mstParams)*gammaMST(n + n0, params->mstParams);
}

Complex Rn_b_coeff(int n, void* p){
	LR_parameters *params = (LR_parameters *)p;
	int n0 = params->n0;

	return betaMST(n + n0, params->mstParams);
}

Complex alphaRn_cf(int n, const MstParameters &params){
	cf_coeffs CF;
	bool success = 0;
	Complex Rn;

	LR_parameters LR_params = {.mstParams = params, .n0 = n};

	CF.a_coeff = &Rn_a_coeff;
	CF.b_coeff = &Rn_b_coeff;
	CF.params = &LR_params;

	cf_solver* cf = cf_solver_alloc(CF);

	int iter = 0;
	while(iter <= CF_ITER_MAX && success == 0){
		iter++;
		success = cf_lentz_iterate(cf, DBL_EPSILON);
	}

	if(iter > CF_ITER_MAX){
		std::cout << "ERROR: cf_solver for Rn did not converge within tolerance "
			<< DBL_EPSILON << " after " << iter << " iterations\n";
	}

	Rn = cf_solver_convergent(cf);
	cf_solver_free(cf);

	return Rn;
}

Complex Rn_cf(int n, const MstParameters &params){
	return alphaRn_cf(n, params)/alphaMST(n-1, params);
}

Complex Ln_a_coeff(int n, void* p){
	LR_parameters *params = (LR_parameters *)p;
	int n0 = params->n0;

	return -alphaMST(-n + n0, params->mstParams)
		*gammaMST(-n + n0 + 1, params->mstParams);
}

Complex Ln_b_coeff(int n, void* p){
	LR_parameters *params = (LR_parameters *)p;
	int n0 = params->n0;

	return betaMST(-n + n0, params->mstParams);
}

Complex gammaLn_cf(int n, const MstParameters &params){
	cf_coeffs CF;
	bool success = 0;
	Complex Ln;

	LR_parameters LR_params = {.mstParams = params, .n0 = n};

	CF.a_coeff = &Ln_a_coeff;
	CF.b_coeff = &Ln_b_coeff;
	CF.params = &LR_params;

	cf_solver* cf = cf_solver_alloc(CF);

	int iter = 0;
	while(iter <= CF_ITER_MAX && success == 0){
		iter++;
		success = cf_lentz_iterate(cf, DBL_EPSILON);
	}

	if(iter > CF_ITER_MAX){
		std::cout << "ERROR: cf_solver for Rn did not converge within tolerance "
			<< DBL_EPSILON << " after " << iter << " iterations\n";
	}

	Ln = cf_solver_convergent(cf);

	cf_solver_free(cf);

	return Ln;
}

Complex gammaLn_recursion(int n, const MstParameters &params){
	int n0 = n;
	double absnu = std::abs(params.getRenormalizedAngularMomentum());

	while(-n0 < absnu){
		n0 --;
	}

	Complex gammaLn = gammaLn_cf(n0, params);

	while(n0 < n){
		n0 ++;
		gammaLn = - gammaMST(n0 + 1, params)
			*alphaMST(n0, params)/(betaMST(n0, params) + gammaLn);
	}

	return gammaLn;
}

Complex Ln_cf(int n, const MstParameters &params){
	return gammaLn_cf(n, params)/gammaMST(n + 1, params);
}

//////////////////////////////
// Characteristic equations //
//////////////////////////////

Complex nu_eqn_136(const MstParameters& params){
	Complex beta = betaMST(0, params);
	return 1. + alphaRn_cf(1, params)/beta + gammaLn_cf(-1, params)/beta;
}

double nu_eqn_136_realnu(double nu, void* p){
	MstParameters *params = (MstParameters *) p;
	(*params).setRenormalizedAngularMomentum(nu);
	double norm, beta = std::real(betaMST(0, *params));
	if(beta == 0.){
		norm = 1.;
	}else{
		norm = beta;
	}
	double eqn136 = beta + std::real(alphaRn_cf(1, *params))
		+ std::real(gammaLn_cf(-1, *params));
	if(isinf(std::abs(eqn136)) || isnan(std::abs(eqn136))){
		eqn136 = norm*(1.e+10);
	}

	return eqn136/norm;
}

double nu_eqn_136_realnu_absnorm(double nu, void* p){
	MstParameters *params = (MstParameters *) p;
	(*params).setRenormalizedAngularMomentum(nu);
	double norm, beta = std::real(betaMST(0, *params));
	if(beta == 0.){
		norm = 1.;
	}else{
		norm = std::abs(beta);
	}
	double eqn136 = beta + std::real(alphaRn_cf(1, *params))
		+ std::real(gammaLn_cf(-1, *params));
	if(isinf(std::abs(eqn136)) || isnan(std::abs(eqn136))){
		eqn136 = norm*(1.e+10);
	}

	return eqn136/norm;
}

Complex nu_eqn_136_complexnu_norm_full(double nu, void* p){
	MstParameters *params = (MstParameters *) p;
	(*params).setRenormalizedAngularMomentum(-I*nu);
	Complex norm, beta = betaMST(0, *params);
	if(std::abs(beta) == 0.){
		norm = 1.;
	}else{
		norm = std::abs(beta);
	}
	Complex eqn136 = beta + alphaRn_cf(1, *params)
		+ gammaLn_cf(-1, *params);

	return eqn136/norm;
}

Complex nu_eqn_136_complexnu_norm_n10(double nu, void* p){
	MstParameters *params = (MstParameters *) p;
	(*params).setRenormalizedAngularMomentum(-I*nu);
	int n = (*params).getSpinWeightedSpheroidalModeNumber();
	Complex norm, beta = betaMST(n, *params);
	if(std::abs(beta) == 0.){
		norm = 1.;
	}else{
		norm = std::abs(beta);
	}
	Complex eqn136 = beta + alphaRn_cf(n + 1, *params)
		+ gammaLn_cf(n - 1, *params);

	return eqn136/norm;
}

double nu_eqn_136_complexnu_norm(double nu, void* p){
	return std::real(nu_eqn_136_complexnu_norm_full(nu, p));
}

double nu_eqn_136_complexnu_norm_imag(double nu, void* p){
	return std::imag(nu_eqn_136_complexnu_norm_full(nu, p));
}

Complex nu_eqn_136_complexnu_half_norm_full(double nu, void* p){
	MstParameters *params = (MstParameters *) p;
	(*params).setRenormalizedAngularMomentum(- 0.5 - I*nu);
	Complex norm, beta = betaMST(0, *params);
	if(std::abs(beta) == 0.){
		norm = 1.;
	}else{
		norm = std::abs(beta);
	}
	Complex eqn136 = beta + alphaRn_cf(1, *params)
		+ gammaLn_cf(-1, *params);

	return eqn136/norm;
}

Complex nu_eqn_136_complexnu_half_norm_n10(double nu, void* p){
	MstParameters *params = (MstParameters *) p;
	(*params).setRenormalizedAngularMomentum(- 0.5 - I*nu);
	int n = (*params).getSpinWeightedSpheroidalModeNumber();
	Complex norm, beta = betaMST(n, *params);
	if(std::abs(beta) == 0.){
		norm = 1.;
	}else{
		norm = std::abs(beta);
	}
	Complex eqn136 = beta + alphaRn_cf(n + 1, *params)
		+ gammaLn_cf(n - 1, *params);

	return eqn136/norm;
}

double nu_eqn_136_complexnu_half_norm(double nu, void* p){
	return std::real(nu_eqn_136_complexnu_half_norm_full(nu, p));
}

double nu_eqn_136_complexnu_half_norm_imag(double nu, void* p){
	return std::imag(nu_eqn_136_complexnu_half_norm_full(nu, p));
}

Complex nu_eqn_136_complexnu_full(double nu, void* p){
	MstParameters *params = (MstParameters *) p;
	(*params).setRenormalizedAngularMomentum(-I*nu);
	Complex eqn136 = betaMST(0, *params) + alphaRn_cf(1, *params)
		+ gammaLn_cf(-1, *params);

	return eqn136;
}

double nu_eqn_136_complexnu(double nu, void* p){
	return std::real(nu_eqn_136_complexnu_full(nu, p));
}

double nu_eqn_136_complexnu_imag(double nu, void* p){
	return std::imag(nu_eqn_136_complexnu_full(nu, p));
}

Complex nu_eqn_136_complexnu_half_full(double nu, void* p){
	MstParameters *params = (MstParameters *) p;
	(*params).setRenormalizedAngularMomentum(- 0.5 - I*nu);
	Complex eqn136 = betaMST(0, *params) + alphaRn_cf(1, *params)
		+ gammaLn_cf(-1, *params);

	return eqn136;
}

double nu_eqn_136_complexnu_half(double nu, void* p){
	return std::real(nu_eqn_136_complexnu_half_full(nu, p));
}

double nu_eqn_136_complexnu_half_imag(double nu, void* p){
	return std::imag(nu_eqn_136_complexnu_half_full(nu, p));
}

double nu_eqn_133_realnu(double nu, void* p){
	MstParameters *params = (MstParameters *) p;
	(*params).setRenormalizedAngularMomentum(nu);
	double alphaRn0 = std::real(Rn_cf(0, *params));
	double gammaLnm1 = std::real(Ln_cf(-1, *params));
	double eqn133 = 1. - alphaRn0*gammaLnm1;

	return eqn133;
}

double nu_eqn_133_complexnu_half(double nu, void* p){
	MstParameters *params = (MstParameters *) p;
	(*params).setRenormalizedAngularMomentum(-0.5 - I*nu);
	int n = 0;
	Complex eqn133 = 1. - Rn_cf(n, *params)*Ln_cf(n - 1, *params);

	return std::real(eqn133);
}

double nu_eqn_133_complexnu(double nu, void* p){
	MstParameters *params = (MstParameters *) p;
	(*params).setRenormalizedAngularMomentum(-I*nu);
	int n = 0;
	Complex eqn133 = 1. - Rn_cf(n, *params)*Ln_cf(n - 1, *params);

	return std::real(eqn133);
}

///////////////////////////
// Root finder functions //
///////////////////////////

Complex nu_solver(double q, int s, int l, int m, double epsilon){
	MstParameters params(q, s, l, m, epsilon);
	nu_solver(params);
	return params.getRenormalizedAngularMomentum();
}

int nu_solver(MstParameters &params){
	int s = params.getSpinWeight(), l = params.getSpinWeightedSpheroidalModeNumber(), m = params.getAzimuthalModeNumber();
	double q = params.getMstQ(), eps = params.getMstEpsilon(), la = params.getSpinWeightedSpheroidalEigenvalue();

	// for low frequencies we can just make use of analytic low-frequency expansions
	if(eps < low_frequency_max_epsilon(q, l, m)){
	  	params.setRenormalizedAngularMomentum(nu_solver_low_freq(q, s, l, m, eps));
	}else{
		// otherwise try to make use of monodromy techniques to get an initial guess for nu.
		params.setRenormalizedAngularMomentum(nu_solver_monodromy(s, l, m, q, eps, la));
		// the monodromy calculation may return 0. if there is catastrophic cancellation that
		// leads to an inaccurate numerical result.
	  	if(std::abs(params.getRenormalizedAngularMomentum()) != 0.){
			nu_solver_guess(params);
	  	}else if(eps < 0.5){
	  		params.setRenormalizedAngularMomentum(nu_solver_low_freq(q, s, l, m, eps));
			nu_solver_guess(params);
	  	}else{
	    	nu_solver_noguess(params);
	  	}
	}
  	// std::cout << std::setprecision(15);
  	// std::cout << "NUSOLVER: nu = " << params.getRenormalizedAngularMomentum() << "\n";
  	// std::cout << std::setprecision(6);

	return 0;
}

int nu_solver_guess(MstParameters &params){
	Complex nu0 = params.getRenormalizedAngularMomentum();
	double nuTest = std::abs(nu_eqn_136(params)/nu0);
	if(params.getMstEpsilon() < 1.e-2){
		nuTest *= params.getMstEpsilon()/2.;
	}
	Complex cl = Complex(params.getSpinWeightedSpheroidalModeNumber());
	if( nuTest > 1.e-4 ){
		if( std::abs(cos(2*M_PI*nu0)) > 1 ){
			params.setRenormalizedAngularMomentum(std::real(nu0) - I*std::imag(nu0));
		}else{
			params.setRenormalizedAngularMomentum(2.*cl - nu0);
		}
		nuTest = std::abs(nu_eqn_136(params)/nu0);
		if( nuTest > 1.e-6 ){
			return nu_solver_noguess(params);
		}
	}
	double deltaNu = std::abs(cl - params.getRenormalizedAngularMomentum());
	if(std::abs(deltaNu) < 1.){
		nuTest *= std::abs(deltaNu);
	}
	if(nuTest < eqn_136_EPS){
		// std::cout << "NUSOLVER: Using monodromy eigenvalue for renormalized angular momentum. \n";
		return 0;
	};
	// std::cout << std::setprecision(15);
	// std::cout << "NUSOLVER: Using monodromy eigenvalue " << params.getRenormalizedAngularMomentum() << " for initial guess value of renormalized angular momentum. \n";
	// std::cout << std::setprecision(6);

	double deltaX = sqrt(nuTest);
	deltaX = deltaX < 0.01 ? deltaX : 0.01;

	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	gsl_function F;
	nu0 = params.getRenormalizedAngularMomentum();

	int status;
	int iter = 0, max_iter = 100;
	double r;

	if(std::abs(cos(2.*M_PI*nu0)) < 1.){
		r = std::real(nu0);
		F.function = &nu_eqn_136_realnu_absnorm;
	}else if(std::real(cos(2.*M_PI*nu0)) < -1.){
		r = std::imag(nu0);
		F.function = &nu_eqn_136_complexnu_half_norm;
	}else {
		r = std::imag(nu0);
		F.function = &nu_eqn_136_complexnu_norm;
	}

	F.params = &params;

	double x_lo = r - deltaX;
	double x_hi = r + deltaX;

	double x_lo_test = (*F.function)(x_lo, &params);
	double x_hi_test = (*F.function)(x_hi, &params);

	if((x_lo_test < 0. && x_hi_test < 0.) || (x_lo_test > 0. && x_hi_test > 0.)){
		r = 2.*params.getSpinWeightedSpheroidalModeNumber() - r;
		x_lo = r - deltaX;
		x_hi = r + deltaX;
	}

	if (!std::isfinite(x_lo_test) || !std::isfinite(x_hi_test)) {
		gsl_root_fsolver_free(s);
		// handle failure gracefully
		std::cerr << "NUSOLVER ERROR: Root solver bounds invalid, exiting." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	if (x_lo_test*x_hi_test > 0) { 
		std::cerr << "NUSOLVER ERROR: Root solver bounds invalid, exiting." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &F, x_lo, x_hi);

	do{
		iter++;
	  status = gsl_root_fsolver_iterate(s);
	  r = gsl_root_fsolver_root(s);
	  x_lo = gsl_root_fsolver_x_lower(s);
	  x_hi = gsl_root_fsolver_x_upper(s);
	  status = gsl_root_test_interval(x_lo, x_hi, 0, DBL_EPSILON);
	}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free(s);

	if(std::abs(cos(2.*M_PI*nu0)) < 1.){
		params.setRenormalizedAngularMomentum(r);
	}else if(std::real(cos(2.*M_PI*nu0)) < -1.){
		params.setRenormalizedAngularMomentum(- 0.5 - I*r);
	}else if(std::real(cos(2.*M_PI*nu0)) > 1.){
		params.setRenormalizedAngularMomentum(- I*r);
	}

	return 0;
};

int nu_solver_noguess_rootfinder(gsl_function F, const double &x_lo, const double &x_hi, MstParameters &params){
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;

	int status;
	int iter = 0, max_iter = 100;
	double r;
	double xlo = x_lo, xhi = x_hi;

	double x_lo_test = (*F.function)(x_lo, &params);
	double x_hi_test = (*F.function)(x_hi, &params);

	if (!std::isfinite(x_lo_test) || !std::isfinite(x_hi_test)) {
		gsl_root_fsolver_free(s);
		// handle failure gracefully
		std::cerr << "NUSOLVER ERROR: Root solver bounds invalid, exiting." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	if (x_lo_test*x_hi_test > 0) {
		std::cerr << "NUSOLVER ERROR: Root solver bounds invalid, exiting." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &F, xlo, xhi);

	do{
		iter++;
	  status = gsl_root_fsolver_iterate(s);
	  r = gsl_root_fsolver_root(s);
	  xlo = gsl_root_fsolver_x_lower(s);
	  xhi = gsl_root_fsolver_x_upper(s);
	  status = gsl_root_test_interval(xlo, xhi, 0, DBL_EPSILON);
	}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free(s);
	params.setRenormalizedAngularMomentum(r);

	return 0;
}

//uses naive secant like method
int realnu_interval_search(double& x_lo, double& x_hi, MstParameters& params){
	return realnu_interval_search(30, x_lo, x_hi, params);
}


int realnu_interval_search(int stepNum, double& x_lo, double& x_hi, MstParameters& params){
	int check, success = 0;
	double stepSize = (x_hi - x_lo)/stepNum;

	int i = 0;
	double x0 = x_lo;
	double x1 = x0 + stepSize;
	double x0_test = nu_eqn_136_realnu_absnorm(x0, &params);
	double x1_test = nu_eqn_136_realnu_absnorm(x1, &params);

	double derivative = (x1_test - x0_test)/stepSize;
	double previousDerivativeNorm;
	i++;
	while(i < stepNum && x0_test/x1_test >= 0){
		previousDerivativeNorm = derivative/x1_test;
		x0 = x1;
		x1 += stepSize;
		x0_test = x1_test;
		x1_test = nu_eqn_136_realnu_absnorm(x1, &params);
		derivative = (x1_test - x0_test)/stepSize;
		if( derivative/x1_test < 0 && previousDerivativeNorm > 0  && stepSize > 1.e-5 && x0_test/x1_test > 0 ){
			check = realnu_interval_search(stepNum, x0, x1, params);
			if(check > 0){
				x_lo = x0;
				x_hi = x1;
				return check;
			}
		}
		i++;
	}

	double x2 = 2.*x1 - x0;
	if( x2 - x1 < (1.e-5) ){ x2 = x1 + (1.e-5); }
	double x2_test = nu_eqn_136_realnu_absnorm(x2, &params);

	if( x1_test/x0_test < 0 && x2_test/x1_test < 0 ){
		x_lo = x0;
		x_hi = x1;
		success = 2;
	}else if( x1_test/x0_test < 0 ){
		x_lo = x0;
		x_hi = x1;
		success = 1;
	}

	return success;
}

int realnu_root_test(Complex& nu, double& nuError, MstParameters& params){
	double realNu = std::real(params.getRenormalizedAngularMomentum());
	double nu_test = nu_eqn_136_realnu(realNu, &params);
	int l = params.getSpinWeightedSpheroidalModeNumber();

	// weight points by how close they are to the singular points of the characteristic equation
	double deltaNu = std::abs(l - realNu);
	if( deltaNu > std::abs(realNu - l - 0.5) ){
		deltaNu = std::abs(realNu - l - 0.5);
	}else if( deltaNu > std::abs(realNu - l + 0.5) ){
		deltaNu = std::abs(realNu - l + 0.5);
	}
	double weight = 10./deltaNu; // penalize any roots that are found very close to integer and half integer values

	if(std::abs(nu_test) == 0. || weight*std::abs(nu_test/realNu) < eqn_136_EPS){
		nu = Complex(realNu);
		nuError = std::abs(nu_test/realNu);
		return 1;
	}else{
		if(std::abs(nu_test/realNu) < nuError){
			nu = Complex(realNu);
			nuError = std::abs(nu_test/realNu);
		}
	}

	return 0;
}

int realnu_interval_search_and_test(int stepNum, Complex& nu, double& nuError, double x_lo, double x_hi, MstParameters& params){
	double x_lo_ref = x_lo;
	double x_hi_ref = x_hi;
	int testSearch, testRoot = 0;
	int rootFlag = 0;
	gsl_function F;
	F.function = &nu_eqn_136_realnu_absnorm;
	F.params = &params;

	testSearch = realnu_interval_search(stepNum, x_lo, x_hi, params);
	if(testSearch > 0){
		// std::cout << std::setprecision(8);
		// std::cout << "NUSOLVER: Evidence of zero crossing in region Re(nu) = ["<< x_lo << ", "<< x_hi << "] \n";
		// std::cout << std::setprecision(6);
		rootFlag = 1;
		nu_solver_noguess_rootfinder(F, x_lo, x_hi, params);
		testRoot = realnu_root_test(nu, nuError, params);
		if( testRoot == 0 && testSearch == 2 ){
			x_lo_ref = x_lo;
			x_lo = x_hi;
			x_hi = 2.*x_lo - x_lo_ref;
			if( x_hi - x_lo < (1.e-5) ){ x_hi = x_lo + (1.e-5); }
			if( x_hi < x_hi_ref ){
				// std::cout << std::setprecision(8);
				// std::cout << "NUSOLVER: Evidence of zero crossing in region Re(nu) = ["<< x_lo << ", "<< x_hi << "] \n";
				// std::cout << std::setprecision(6);
				nu_solver_noguess_rootfinder(F, x_lo, x_hi, params);
				testRoot = realnu_root_test(nu, nuError, params);
			}else{
				x_hi = x_lo;
			}
		}
	}

	while( testRoot == 0 && std::abs(x_hi - x_lo) > (1.e-5) ){
		x_lo = x_hi;
		x_hi = x_hi_ref;
		testSearch = realnu_interval_search(stepNum, x_lo, x_hi, params);
		if(testSearch > 0){
			// std::cout << std::setprecision(8);
			// std::cout << "NUSOLVER: Evidence of zero crossing in region Re(nu) = ["<< x_lo << ", "<< x_hi << "] \n";
			// std::cout << std::setprecision(6);
			rootFlag = 1;
			nu_solver_noguess_rootfinder(F, x_lo, x_hi, params);
			testRoot = realnu_root_test(nu, nuError, params);
			if( testRoot == 0 && testSearch == 2 ){
				x_lo_ref = x_lo;
				x_lo = x_hi;
				if( x_hi - x_lo < (1.e-5) ){ x_hi = x_lo + (1.e-5); }
				if( x_hi < x_hi_ref ){
					// std::cout << std::setprecision(8);
					// std::cout << "NUSOLVER: Evidence of zero crossing in region Re(nu) = ["<< x_lo << ", "<< x_hi << "] \n";
					// std::cout << std::setprecision(6);
					nu_solver_noguess_rootfinder(F, x_lo, x_hi, params);
					testRoot = realnu_root_test(nu, nuError, params);
				}else{
					x_hi = x_lo;
				}
			}
		}
	}

	if(rootFlag > 0 && testRoot > 0){
		return 1;
	}else if(rootFlag > 0){
		return -1;
	}else{
		return 0;
	}
}

int nu_solver_real_noguess(MstParameters &params){
	double freqWeight = params.getMstEpsilon();
	int stepNum = (10 + round(15*freqWeight));

	return nu_solver_real_noguess(stepNum, params);
}

int nu_solver_real_noguess(int stepNum, MstParameters &params){
	Complex nu = 0.;
	double nuError = 1.;
	int l = params.getSpinWeightedSpheroidalModeNumber();
	int testSearch;
	int stepNumHighRes = stepNum/8;

	double x_lo, x_hi;
	double singularWidth = 0.1;

	x_lo = l - 0.5 + (1.e-6);
	x_hi = l - 0.5 + singularWidth;
	testSearch = realnu_interval_search_and_test(stepNumHighRes, nu, nuError, x_lo, x_hi, params);
	if( testSearch == 1 ){
		return 0;
	}

	x_lo = x_hi;
	x_hi = l - singularWidth;
	testSearch = realnu_interval_search_and_test(stepNum, nu, nuError, x_lo, x_hi, params);
	if( testSearch == 1 ){
		return 0;
	}

	x_lo = x_hi;
	x_hi = l - (1.e-6);
	testSearch = realnu_interval_search_and_test(stepNumHighRes, nu, nuError, x_lo, x_hi, params);
	if( testSearch == 1 ){
		return 0;
	}

	x_hi = l + 0.5 - (1.e-6);
	x_lo = l + 0.5 - singularWidth;
	testSearch = realnu_interval_search_and_test(stepNumHighRes, nu, nuError, x_lo, x_hi, params);
	if( testSearch == 1 ){
		return 0;
	}

	x_hi = x_lo;
	x_lo = l + singularWidth;
	testSearch = realnu_interval_search_and_test(stepNum, nu, nuError, x_lo, x_hi, params);
	if( testSearch == 1 ){
		return 0;
	}

	x_hi = x_lo;
	x_lo = l + (1.e-6);
	testSearch = realnu_interval_search_and_test(stepNumHighRes, nu, nuError, x_lo, x_hi, params);
	if( testSearch == 1 ){
		return 0;
	}

	double freq = std::abs(params.getMstEpsilon())/2.;
	double detectionThreshold = 1.e-10;
	if( freq > 6.){
		detectionThreshold *= 10;
	}
	if( freq > 7.){
		detectionThreshold *= 10;
	}
	if( freq > 8.){
		detectionThreshold *= 10;
	}
	if( freq > 9.){
		detectionThreshold *= 10;
	}

	if(nuError < detectionThreshold){
		// std::cout << "NUSOLVER: Using less accurate root at nu = " << nu << " with precision of " << nuError << ". \n";
		params.setRenormalizedAngularMomentum(nu);
	}else{
		// std::cout << "NUSOLVER: Not using less accurate root at nu = " << nu << " with precision of " << nuError << ". \n";
		params.setRenormalizedAngularMomentum(0.);
	}

	// if(rootFlag == 0){
	// 	std::cout << "NUSOLVER: Unable to identify suitable search region to apply root-finding. \n";
	// }

	return 0;
}

int nu_solver_complex_noguess(MstParameters &params){
	double r, r_test;
	gsl_function Fhalf, Fone, Fhalf2, Fone2;

	// blindly stepping through the complex plane is tricky, so need to evaluate using four different characteristic equations to help with root-finding,
	// two help search along Re(nu) = -0.5 and two help search along Re(nu) = 0.
	Fhalf.function = &nu_eqn_133_complexnu_half;
	Fhalf.params = &params;

	Fone.function = &nu_eqn_133_complexnu;
	Fone.params = &params;

	Fhalf2.function = &nu_eqn_136_complexnu_half;
	Fhalf2.params = &params;

	Fone2.function = &nu_eqn_136_complexnu;
	Fone2.params = &params;

	double x_min = (1.e+2)*DBL_EPSILON;
	double x_max = params.getMstEpsilon();
	int x_range_inc = 100, x_range_increase = 100;
	double x_inc = (x_max-x_min)/(x_range_inc-1);
	Complex potentialRoot = 0.;
	double potentialRootError = 1.;

	double x_range_val[x_range_inc];
	double x_lo, x_hi, x_lo_test, x_hi_test;
	double endpoint_min = 1.e-6;

	for(int i = 0; i < x_range_inc; i++){
		x_range_val[i] = x_min + x_inc*i;
	}

	for(int i = 0; i < x_range_inc; i++){
		x_lo = x_range_val[x_range_inc-i-2];
		x_hi = x_range_val[x_range_inc-i-1];

		// Test Eqn 133 assuming Re(nu) = -0.5
		x_lo_test = nu_eqn_133_complexnu_half(x_lo, &params);
		x_hi_test = nu_eqn_133_complexnu_half(x_hi, &params);

		if(((x_lo_test > 0. && x_hi_test < 0.) || (x_lo_test < 0. && x_hi_test > 0.)) && std::abs(x_lo_test) > endpoint_min && std::abs(x_hi_test) > endpoint_min){
			nu_solver_noguess_rootfinder(Fhalf, x_lo, x_hi, params);
			r = std::real(params.getRenormalizedAngularMomentum()); // I call real part because nu_solver_guess_rootfinder only works with real variables
			// therefore, even though I'm searching for the imaginary part of nu, it gets temporarily stored as the real part of nu by the root finder
			r_test = std::abs(nu_eqn_136_complexnu_half_norm_full(r, &params)); // Complex form of nu is then properly stored by call to nu_eqn_136_complexnu_half_full
			if(std::abs(r_test) == 0 || std::abs(r_test/r) < eqn_136_EPS){
				return 0;
			}else{
				//std::cout << "NUSOLVER: Root found at nu = "<< -0.5 - I*r <<", but did not meet initial precision requirements \n";
				if(std::abs(r_test/r) < potentialRootError){
					potentialRoot = -0.5 - I*r;
					potentialRootError = std::abs(r_test/r);
				}
			}
		}

		if((std::abs(x_lo_test) > 1. && std::abs(x_hi_test) < 1.) || (std::abs(x_lo_test) < 1. && std::abs(x_hi_test) > 1.)){
			// std::cout << "NUSOLVER: Evidence of zero crossing in region Im(nu) = ["<< x_lo << ", "<< x_hi << "] \n";
			for(int j = 0; j < x_range_increase + 1; j++){
				x_lo = x_range_val[x_range_inc-i-2] + j*x_inc/x_range_increase;
				x_hi = x_range_val[x_range_inc-i-2] + (j+1)*x_inc/x_range_increase;

				x_lo_test = nu_eqn_133_complexnu_half(x_lo, &params);
				x_hi_test = nu_eqn_133_complexnu_half(x_hi, &params);

				if(((x_lo_test > 0. && x_hi_test < 0.) || (x_lo_test < 0. && x_hi_test > 0.)) && std::abs(x_lo_test) > endpoint_min && std::abs(x_hi_test) > endpoint_min)
				{
					nu_solver_noguess_rootfinder(Fhalf, x_lo, x_hi, params);
					r = std::real(params.getRenormalizedAngularMomentum());
					r_test = std::abs(nu_eqn_136_complexnu_half_norm_full(r, &params));
					if(std::abs(r_test) == 0. || std::abs(r_test)/std::abs(r) < eqn_136_EPS){
						return 0;
					}else{
						//std::cout << "NUSOLVER: Root found at nu = "<< -0.5 - I*r <<", but did not meet initial precision requirements \n";
						if(std::abs(r_test/r) < potentialRootError){
							potentialRoot = -0.5 - I*r;
							potentialRootError = std::abs(r_test/r);
						}
					}
				}
			}
		}

		// Test Eqn 136 assuming Re(nu) = -0.5
		x_lo_test = nu_eqn_136_complexnu_half(x_lo, &params);
		x_hi_test = nu_eqn_136_complexnu_half(x_hi, &params);

		if(((x_lo_test > 0. && x_hi_test < 0.) || (x_lo_test < 0. && x_hi_test > 0.)) && std::abs(x_lo_test) > endpoint_min && std::abs(x_hi_test) > endpoint_min){
			nu_solver_noguess_rootfinder(Fhalf2, x_lo, x_hi, params);
			r = std::real(params.getRenormalizedAngularMomentum()); // I call real part because nu_solver_guess_rootfinder only works with real variables
			// therefore, even though I'm searching for the imaginary part of nu, it gets temporarily stored as the real part of nu by the root finder
			r_test = std::abs(nu_eqn_136_complexnu_half_norm_full(r, &params)); // Complex form of nu is then properly stored by call to nu_eqn_136_complexnu_half_full
			if(std::abs(r_test) == 0 || std::abs(r_test/r) < eqn_136_EPS){
				return 0;
			}else{
				//std::cout << "NUSOLVER: Root found at nu = "<< -0.5 - I*r <<", but did not meet initial precision requirements \n";
				if(std::abs(r_test/r) < potentialRootError){
					potentialRoot = -0.5 - I*r;
					potentialRootError = std::abs(r_test/r);
				}
			}
		}

		if(((x_lo_test > 0. && x_hi_test < 0.) || (x_lo_test < 0. && x_hi_test > 0.)) && std::abs(x_lo_test) > endpoint_min && std::abs(x_hi_test) > endpoint_min){
			// std::cout << "NUSOLVER: Evidence of zero crossing in region Im(nu) = ["<< x_lo << ", "<< x_hi << "] \n";
			for(int j = 0; j < x_range_increase + 1; j++){
				x_lo = x_range_val[x_range_inc-i-2] + j*x_inc/x_range_increase;
				x_hi = x_range_val[x_range_inc-i-2] + (j+1)*x_inc/x_range_increase;

				x_lo_test = nu_eqn_136_complexnu_half(x_lo, &params);
				x_hi_test = nu_eqn_136_complexnu_half(x_hi, &params);

				if((x_lo_test > 0. && x_hi_test < 0.) || (x_lo_test < 0. && x_hi_test > 0.))
				{
					nu_solver_noguess_rootfinder(Fhalf2, x_lo, x_hi, params);
					r = std::real(params.getRenormalizedAngularMomentum());
					r_test = std::abs(nu_eqn_136_complexnu_half_norm_full(r, &params));
					if(std::abs(r_test) == 0. || std::abs(r_test)/std::abs(r) < eqn_136_EPS){
						return 0;
					}else{
						//std::cout << "NUSOLVER: Root found at nu = "<< -0.5 - I*r <<", but did not meet initial precision requirements \n";
						if(std::abs(r_test/r) < potentialRootError){
							potentialRoot = -0.5 - I*r;
							potentialRootError = std::abs(r_test/r);
						}
					}
				}
			}
		}

		// Test Eqn 133 assuming Re(nu) = 0
		x_lo_test = nu_eqn_133_complexnu(x_lo, &params);
		x_hi_test = nu_eqn_133_complexnu(x_hi, &params);

		if(((x_lo_test > 0. && x_hi_test < 0.) || (x_lo_test < 0. && x_hi_test > 0.)) && std::abs(x_lo_test) > endpoint_min && std::abs(x_hi_test) > endpoint_min){
			nu_solver_noguess_rootfinder(Fone, x_lo, x_hi, params);
			r = std::real(params.getRenormalizedAngularMomentum());
			r_test = std::abs(nu_eqn_136_complexnu_norm_full(r, &params));
			if(std::abs(r_test) == 0 || std::abs(r_test)/std::abs(r) < eqn_136_EPS){
				return 0;
			}else{
				//std::cout << "NUSOLVER: Root found at nu = "<< -I*r <<", but did not meet initial precision requirements \n";
				if(std::abs(r_test/r) < potentialRootError){
					potentialRoot = - I*r;
					potentialRootError = std::abs(r_test/r);
				}
			}
		}

		if(((x_lo_test > 0. && x_hi_test < 0.) || (x_lo_test < 0. && x_hi_test > 0.)) && std::abs(x_lo_test) > endpoint_min && std::abs(x_hi_test) > endpoint_min){
			// std::cout << "NUSOLVER: Evidence of zero crossing in region Im(nu) = ["<< x_lo << ", "<< x_hi << "] \n";
			for(int j = 0; j < x_range_increase + 1; j++){
				x_lo = x_range_val[x_range_inc-i-2] + j*x_inc/x_range_increase;
				x_hi = x_range_val[x_range_inc-i-2] + (j+1)*x_inc/x_range_increase;

				x_lo_test = nu_eqn_133_complexnu(x_lo, &params);
				x_hi_test = nu_eqn_133_complexnu(x_hi, &params);

				if((x_lo_test > 0. && x_hi_test < 0.) || (x_lo_test < 0. && x_hi_test > 0.))
				{
					nu_solver_noguess_rootfinder(Fone, x_lo, x_hi, params);
					r = std::real(params.getRenormalizedAngularMomentum());
					r_test = std::abs(nu_eqn_136_complexnu_norm_full(r, &params));
					if(std::abs(r_test) == 0. || std::abs(r_test/r) < eqn_136_EPS){
						return 0;
					}else{
						//std::cout << "NUSOLVER: Root found at nu = "<< I*r <<", but did not meet initial precision requirements \n";
						if(std::abs(r_test/r) < potentialRootError){
							potentialRoot = - I*r;
							potentialRootError = std::abs(r_test/r);
						}
					}
				}
			}
		}

		// Test Eqn 136 assuming Re(nu) = 0
		x_lo_test = nu_eqn_136_complexnu(x_lo, &params);
		x_hi_test = nu_eqn_136_complexnu(x_hi, &params);

		if(((x_lo_test > 0. && x_hi_test < 0.) || (x_lo_test < 0. && x_hi_test > 0.)) && std::abs(x_lo_test) > endpoint_min && std::abs(x_hi_test) > endpoint_min){
			nu_solver_noguess_rootfinder(Fone2, x_lo, x_hi, params);
			r = std::real(params.getRenormalizedAngularMomentum());
			r_test = std::abs(nu_eqn_136_complexnu_norm_full(r, &params));
			if(std::abs(r_test) == 0 || std::abs(r_test)/std::abs(r) < eqn_136_EPS){
				// std::cout << "Searching region endpoints x_lo_test = " << x_lo_test << ", x_hi_test = " << x_hi_test << "\n";
				return 0;
			}else{
				//std::cout << "NUSOLVER: Root found at nu = "<< -I*r <<", but did not meet initial precision requirements \n";
				if(std::abs(r_test/r) < potentialRootError){
					potentialRoot = - I*r;
					potentialRootError = std::abs(r_test/r);
				}
			}
		}

		if((std::abs(x_lo_test) > 1. && std::abs(x_hi_test) < 1.) || (std::abs(x_lo_test) < 1. && std::abs(x_hi_test) > 1.)){
			// std::cout << "NUSOLVER: Evidence of zero crossing in region Im(nu) = ["<< x_lo << ", "<< x_hi << "] \n";
			for(int j = 0; j < x_range_increase + 1; j++){
				x_lo = x_range_val[x_range_inc-i-2] + j*x_inc/x_range_increase;
				x_hi = x_range_val[x_range_inc-i-2] + (j+1)*x_inc/x_range_increase;

				x_lo_test = nu_eqn_136_complexnu(x_lo, &params);
				x_hi_test = nu_eqn_136_complexnu(x_hi, &params);

				if(((x_lo_test > 0. && x_hi_test < 0.) || (x_lo_test < 0. && x_hi_test > 0.)) && std::abs(x_lo_test) > endpoint_min && std::abs(x_hi_test) > endpoint_min)
				{
					nu_solver_noguess_rootfinder(Fone2, x_lo, x_hi, params);
					r = std::real(params.getRenormalizedAngularMomentum());
					r_test = std::abs(nu_eqn_136_complexnu_norm_full(r, &params));
					if(std::abs(r_test) == 0. || std::abs(r_test/r) < eqn_136_EPS){
						return 0;
					}else{
						//std::cout << "NUSOLVER: Root found at nu = "<< I*r <<", but did not meet initial precision requirements \n";
						if(std::abs(r_test/r) < potentialRootError){
							potentialRoot = - I*r;
							potentialRootError = std::abs(r_test/r);
						}
					}
				}
			}
		}
	}

	if(potentialRootError < std::abs(x_max*(1.e-10))){
		// std::cout << "NUSOLVER: Using less accurate root at nu = " << potentialRoot << " with precision of " << potentialRootError << ". \n";
		params.setRenormalizedAngularMomentum(potentialRoot);
	}else{
		params.setRenormalizedAngularMomentum(0.);
	}

	return 0;
};

int nu_solver_noguess(MstParameters &params){
	double freqWeight = params.getMstEpsilon()/12.;
	int stepNum = round(12*exp(freqWeight));

	// std::cout << "NUSOLVER: Searching for renormalized angular momentum along the imaginary line. \n";
	nu_solver_complex_noguess(params);

	if( std::abs(params.getRenormalizedAngularMomentum()) == 0.){
		// std::cout << "NUSOLVER: Searching for renormalized angular momentum along the real line with "<< stepNum <<" steps. \n";
		nu_solver_real_noguess(stepNum, params);
		if( std::abs(params.getRenormalizedAngularMomentum()) == 0. ){
			int resIncreaseMax = 0;
			while( std::abs(params.getRenormalizedAngularMomentum()) == 0. && resIncreaseMax < 2){
				stepNum *= 4;
				// std::cout << "NUSOLVER: Searching for renormalized angular momentum along the real line with "<< stepNum <<" steps. \n";
				resIncreaseMax ++;
				nu_solver_real_noguess(stepNum, params);
			}
			if( std::abs(params.getRenormalizedAngularMomentum()) == 0. ){
				std::cout << "NUSOLVER: Root-finding methods failed for computing renormalized angular momentum. \n";
			}
		}
	}

	return 0;
};

int nu_solver_noguess2(MstParameters &params){
	double freqWeight = params.getMstEpsilon()/12.;
	int stepNum = round(20*exp(freqWeight));
	// std::cout << "NUSOLVER: Searching for renormalized angular momentum along the real line with "<< stepNum <<" steps. \n";
	nu_solver_real_noguess(stepNum, params);
	Complex nuTemp = params.getRenormalizedAngularMomentum();
	Complex l = Complex(params.getSpinWeightedSpheroidalModeNumber());
	int nuFlag = 0;
	if(std::abs(l + 0.5 - nuTemp) < 1.e-1 || std::abs(l - 0.5 - nuTemp) < 1.e-1 || std::abs(l - nuTemp) < 1.e-1 ){
		nuFlag = 1;
	}
	//////////////////////////////
	// FINISH SETTING THIS UP!!!! You left to go make dinner before you could finish this
	//////////////////////////////
	if( std::abs(params.getRenormalizedAngularMomentum()) == 0. || nuFlag == 1){
		// std::cout << "NUSOLVER: Initial search for real nu inconclusive, searching for renormalized angular momentum along the imaginary line. \n";
		nu_solver_complex_noguess(params);
		if( std::abs(params.getRenormalizedAngularMomentum()) == 0. ){
			if( nuFlag == 1 ){
				params.setRenormalizedAngularMomentum(nuTemp);
			}else{
				int resIncreaseMax = 0;
				while( std::abs(params.getRenormalizedAngularMomentum()) == 0. && resIncreaseMax < 2){
					stepNum *= 4;
					// std::cout << "NUSOLVER: Searching for renormalized angular momentum along the real line with "<< stepNum <<" steps. \n";
					resIncreaseMax ++;
					nu_solver_real_noguess(stepNum, params);
				}
				if( std::abs(params.getRenormalizedAngularMomentum()) == 0. ){
					std::cout << "NUSOLVER: Root-finding methods failed for computing renormalized angular momentum. \n";
				}
			}
		}
	}

	return 0;
};

//*************************************************************************************
// Low frequency expansions
//*************************************************************************************

double nu_solver_low_freq(double q, int s, int l, int m, double eps){
	if( std::abs(s) == 2 ){
		return nu_solver_low_freq_s2(q, l, m, eps);
	}else if( s == 0 ){
		return nu_solver_low_freq_s0(q, l, m, eps);
	}else{
		std::cout << "NUSOLVER: No low-frequency expansions for s = " << s << "\n";
		return 0.;
	}
}

double low_frequency_max_epsilon(double q, int l, int m){
	if(l == 0){
		return 5.8e-2; // based on s = 0, l = 0 expansions truncated at order O[eps]^12
	}else if(l == 1){
		return 6.7e-2; // based on s = 0, l = 1 expansions truncated at order O[eps]^12
	}else if(l == 2){
		return 6.1e-2; // based on s = -2, l = 2 expansions truncated at order O[eps]^12
	}else if(l == 3){
		return 8.6e-2; // based on s = -2, l = 3 expansions truncated at order O[eps]^12
	}else if(l == 4){
		return 0.13; // based on s = -2, l = 4 expansions truncated at order O[eps]^12
	}else if(l == 5){
		return 0.17; // based on s = -2, l = 5 expansions truncated at order O[eps]^12
	}else if(l == 6){
		return 0.21; // based on s = -2, l = 6 expansions truncated at order O[eps]^12
	}else{
		return max_epsilon_power_law(q, l, m); // based on general l expansions truncated at order O[eps]^6
	}
}

double max_epsilon_power_law(double q, int l, int m){
	double alpha = -0.142;
	double beta = -2.094;
	double x = log10(std::abs(gen_l_last_term_s2(q, l, m)));

	return 1.2*pow(10., alpha*x + beta);
}

double nu_solver_low_freq_s2(double q, int l, int m, double eps){
	if( l == 2 ){
		return nu_solver_low_freq_s2_l2(q, m, eps);
	}else if( l == 3 ){
		return nu_solver_low_freq_s2_l3(q, m, eps);
	}else if( l == 4 ){
		return nu_solver_low_freq_s2_l4(q, m, eps);
	}else if( l == 5 ){
		return nu_solver_low_freq_s2_l5(q, m, eps);
	}else if( l == 6 ){
		return nu_solver_low_freq_s2_l6(q, m, eps);
	}else{
		return nu_solver_low_freq_s2_l(q, l, m, eps);
	}
}

double nu_solver_low_freq_s0(double q, int l, int m, double eps){
	if( l == 2 ){
		return nu_solver_low_freq_s0_l0(q, eps);
	}else if( l == 3 ){
		return nu_solver_low_freq_s0_l1(q, m, eps);
	}else if( l == 4 ){
		return nu_solver_low_freq_s0_l2(q, m, eps);
	}else if( l == 5 ){
		return nu_solver_low_freq_s0_l5(q, m, eps);
	}else if( l == 6 ){
		return nu_solver_low_freq_s0_l6(q, m, eps);
	}else{
		return nu_solver_low_freq_s0_l(q, l, m, eps);
	}
}

double nu_solver_low_freq_s0_l0(double q, double eps){
	double nu2 = -(7./6);
	double nu3 = 0.;
	double nu4 = (-9449 + 648*pow(q,2))/7560.;
	double nu5 = 0.;
	double nu6 = (-3.0675478110115297 + 0.1049757045675413*pow(q,2) - 0.0007968901846452867*pow(q,4));
	double nu7 = 0.;
	double nu8 = (-8.482945357396959 + 0.6114670048561482*pow(q,2) - 0.0000378171352765272*pow(q,4) +
   	 0.00005889888619306917*pow(q,6));
	double nu9 = 0.;
	double nu10 = (-27.852492619455777 + 1.946571784834161*pow(q,2) - 0.04771972578996881*pow(q,4) +
   	 0.00009723309320360614*pow(q,6) + 2.1461504040516796e-6*pow(q,8));
	double nu11 = 0.;
	double nu12 = (-94.9403697061614 + 8.661870721237182*pow(q,2) - 0.1606661480625339*pow(q,4) +
   	 0.0022603426139676825*pow(q,6) + 5.8686933482897e-6*pow(q,8) +
   	  1.5647496358048979e-7*pow(q,10));

	return 0. + (nu2*pow(eps, 2) + nu3*pow(eps, 3) + nu4*pow(eps, 4) + nu5*pow(eps, 5) + nu6*pow(eps, 6)
		+ nu7*pow(eps, 7) + nu8*pow(eps, 8) + nu9*pow(eps, 9) + nu10*pow(eps, 10) + nu11*pow(eps, 11) + nu12*pow(eps, 12));
}

double nu_solver_low_freq_s0_l1(double q, int m, double eps){
	double nu2 = -(19./30);
	double nu3 = (432*m*q + 100*m*pow(q,3) - 100*pow(m,3)*pow(q,3))/2280.;
	double nu4 = (-1913593132 + 182463840*pow(q,2) - 156620160*pow(m,2)*pow(q,2) +
     50035125*pow(m,2)*pow(q,4) - 50035125*pow(m,4)*pow(q,4) +
     7875000*pow(m,2)*pow(q,6) - 15750000*pow(m,4)*pow(q,6) + 7875000*pow(m,6)*pow(q,6))/5.185404e9;
	double nu5 = (0.2934258915988031*m*q + 0.1294277167217829*m*pow(q,3) -
   	 0.14175873523883237*pow(m,3)*pow(q,3) + 0.0037526294961781185*m*pow(q,5) -
   	  0.006941491811328902*pow(m,3)*pow(q,5) + 0.0031888623151507833*pow(m,5)*pow(q,5) +
   	   0.0011225725200989317*pow(m,3)*pow(q,7) - 0.0022451450401978633*pow(m,5)*pow(q,7) +
   		0.0011225725200989317*pow(m,7)*pow(q,7) + 0.00010517215453286265*pow(m,3)*pow(q,9) -
   		 0.00031551646359858796*pow(m,5)*pow(q,9) + 0.00031551646359858796*pow(m,7)*pow(q,9) -
   		  0.00010517215453286265*pow(m,9)*pow(q,9));
	double nu6 = (-0.4595887076118041 + 0.053250988174829636*pow(q,2) - 0.04707682858820499*pow(m,2)*pow(q,2) +
   	 0.0002402492182176851*pow(q,4) + 0.0407393674870512*pow(m,2)*pow(q,4) -
   	  0.04130420255964668*pow(m,4)*pow(q,4) + 0.008905757076141217*pow(m,2)*pow(q,6) -
   	   0.01865530820341448*pow(m,4)*pow(q,6) + 0.00974955112727326*pow(m,6)*pow(q,6) +
   		0.0003442555060543679*pow(m,2)*pow(q,8) - 0.0005819467256328879*pow(m,4)*pow(q,8) +
   		 0.00013112693310267216*pow(m,6)*pow(q,8) + 0.00010656428647584788*pow(m,8)*pow(q,8) +
   		  0.0001323428416855527*pow(m,4)*pow(q,10) - 0.00039702852505665804*pow(m,6)*pow(q,10) +
   		   0.00039702852505665804*pow(m,8)*pow(q,10) - 0.0001323428416855527*pow(m,10)*pow(q,10) +
   			9.10423775388354e-6*pow(m,4)*pow(q,12) - 0.00003641695101553416*pow(m,6)*pow(q,12) +
   			 0.000054625426523301236*pow(m,8)*pow(q,12) -
   			  0.00003641695101553416*pow(m,10)*pow(q,12) + 9.10423775388354e-6*pow(m,12)*pow(q,12));
	double nu7 = (0.5291169514955639*m*q + 0.2123086914772289*m*pow(q,3) -
   	 0.23085536762225395*pow(m,3)*pow(q,3) + 0.01189385055743061*m*pow(q,5) -
   	  0.015224349506993208*pow(m,3)*pow(q,5) + 0.003402073575104707*pow(m,5)*pow(q,5) +
   	   0.00024715300814972683*m*pow(q,7) + 0.006941125522021893*pow(m,3)*pow(q,7) -
   		0.014961907932857002*pow(m,5)*pow(q,7) + 0.007773629402685382*pow(m,7)*pow(q,7) +
   		 0.0010883520689754704*pow(m,3)*pow(q,9) - 0.0032007649938002453*pow(m,5)*pow(q,9) +
   		  0.003136473780674079*pow(m,7)*pow(q,9) - 0.0010240608558493042*pow(m,9)*pow(q,9) +
   		   0.000038682293465563326*pow(m,3)*pow(q,11) - 0.00006463350040419206*pow(m,5)*pow(q,11) -
   			0.000038193259580803784*pow(m,7)*pow(q,11) + 0.00011555784651193042*pow(m,9)*pow(q,11) -
   			 0.000051413379992497914*pow(m,11)*pow(q,11) +
   			  0.000016182929721249774*pow(m,5)*pow(q,13) - 0.0000647317188849991*pow(m,7)*pow(q,13) +
   			   0.00009709757832749865*pow(m,9)*pow(q,13) - 0.0000647317188849991*pow(m,11)*pow(q,13) +
   				0.000016182929721249774*pow(m,13)*pow(q,13) + 8.826823307089304e-7*pow(m,5)*pow(q,15) -
   				 4.413411653544652e-6*pow(m,7)*pow(q,15) + 8.826823307089304e-6*pow(m,9)*pow(q,15) -
   				  8.826823307089304e-6*pow(m,11)*pow(q,15) + 4.413411653544652e-6*pow(m,13)*pow(q,15) -
   				   8.826823307089304e-7*pow(m,15)*pow(q,15));
	double nu8 = (-0.7178540903925492 + 0.09894096923257033*pow(q,2) - 0.182676854160643*pow(m,2)*pow(q,2) +
   	 0.0012983910802961077*pow(q,4) + 0.031234119008236815*pow(m,2)*pow(q,4) -
   	  0.035404709581739253*pow(m,4)*pow(q,4) + 0.00002488560557476383*pow(q,6) +
   	   0.01964916397010983*pow(m,2)*pow(q,6) - 0.039945201321303254*pow(m,4)*pow(q,6) +
   		0.020289961390480338*pow(m,6)*pow(q,6) + 0.0020520835169713864*pow(m,2)*pow(q,8) -
   		 0.0023130225794602056*pow(m,4)*pow(q,8) - 0.001641580555309943*pow(m,6)*pow(q,8) +
   		  0.001902519617798762*pow(m,8)*pow(q,8) + 0.00004793632702314239*pow(m,2)*pow(q,10) +
   		   0.0011322816439502386*pow(m,4)*pow(q,10) - 0.0037195804010853252*pow(m,6)*pow(q,10) +
   			0.003850570562227365*pow(m,8)*pow(q,10) - 0.0013112081321154212*pow(m,10)*pow(q,10) +
   			 0.00014373592300348628*pow(m,4)*pow(q,12) - 0.0005386993827058325*pow(m,6)*pow(q,12) +
   			  0.0007536826100965799*pow(m,8)*pow(q,12) - 0.00046621076408960723*pow(m,10)*pow(q,12) +
   			   0.00010749161369537366*pow(m,12)*pow(q,12) + 4.633325651366992e-6*pow(m,4)*pow(q,14) -
   				7.828393501040643e-6*pow(m,6)*pow(q,14) - 0.000015019682509507349*pow(m,8)*pow(q,14) +
   				 0.000045696152021095986*pow(m,10)*pow(q,14) -
   				  0.00003818631076634231*pow(m,12)*pow(q,14) +
   				   0.000010704909104427325*pow(m,14)*pow(q,14) + 2.0272454244704127e-6*pow(m,6)*pow(q,16) -
   					0.000010136227122352061*pow(m,8)*pow(q,16) +
   					 0.000020272454244704122*pow(m,10)*pow(q,16) -
   					  0.000020272454244704122*pow(m,12)*pow(q,16) +
   					   0.000010136227122352061*pow(m,14)*pow(q,16) -
   						2.0272454244704127e-6*pow(m,16)*pow(q,16) + 9.169137784372546e-8*pow(m,6)*pow(q,18) -
   						 5.501482670623528e-7*pow(m,8)*pow(q,18) + 1.375370667655882e-6*pow(m,10)*pow(q,18) -
   						  1.8338275568745096e-6*pow(m,12)*pow(q,18) + 1.375370667655882e-6*pow(m,14)*pow(q,18) -
   						   5.501482670623528e-7*pow(m,16)*pow(q,18) + 9.169137784372546e-8*pow(m,18)*pow(q,18));
	double nu9 = (1.0672635589215158*m*q + 0.3841716788791914*m*pow(q,3) -
   	 0.40736952202561016*pow(m,3)*pow(q,3) + 0.013617626593511135*m*pow(q,5) -
   	  0.003853876229172692*pow(m,3)*pow(q,5) - 0.010460279274473626*pow(m,5)*pow(q,5) +
   	   0.001037098932290521*m*pow(q,7) + 0.01892480744897247*pow(m,3)*pow(q,7) -
   		0.04229716909657789*pow(m,5)*pow(q,7) + 0.022348314672579295*pow(m,7)*pow(q,7) +
   		 0.000017120136479663127*m*pow(q,9) + 0.004646896647648948*pow(m,3)*pow(q,9) -
   		  0.01293641762348213*pow(m,5)*pow(q,9) + 0.011842648694196798*pow(m,7)*pow(q,9) -
   		   0.003570247854843279*pow(m,9)*pow(q,9) + 0.00035905645346297367*pow(m,3)*pow(q,11) -
   			0.0004529878387759747*pow(m,5)*pow(q,11) - 0.0008281732343674925*pow(m,7)*pow(q,11) +
   			 0.0015790841712110142*pow(m,9)*pow(q,11) - 0.0006569795515305208*pow(m,11)*pow(q,11) +
   			  8.141211391826661e-6*pow(m,3)*pow(q,13) + 0.00018174293697810153*pow(m,5)*pow(q,13) -
   			   0.000802820040395433*pow(m,7)*pow(q,13) + 0.0012365903853994232*pow(m,9)*pow(q,13) -
   				0.0008323984544840867*pow(m,11)*pow(q,13) + 0.0002087439611101684*pow(m,13)*pow(q,13) +
   				 0.000019657874956243357*pow(m,5)*pow(q,15) - 0.00008893698390051544*pow(m,7)*pow(q,15) +
   				  0.0001591691860396282*pow(m,9)*pow(q,15) - 0.00014046440427822548*pow(m,11)*pow(q,15) +
   				   0.000060879811258411385*pow(m,13)*pow(q,15) -
   					0.000010305484075542006*pow(m,15)*pow(q,15) + 5.737784450095071e-7*pow(m,5)*pow(q,17) -
   					 9.86888006268614e-7*pow(m,7)*pow(q,17) - 3.6722366437995366e-6*pow(m,9)*pow(q,17) +
   					  0.000013082257737694144*pow(m,11)*pow(q,17) -
   					   0.000015951149962741678*pow(m,13)*pow(q,17) + 8.8362426488851e-6*pow(m,15)*pow(q,17) -
   						1.8820042187789215e-6*pow(m,17)*pow(q,17) + 2.58189966378374e-7*pow(m,7)*pow(q,19) -
   						 1.549139798270244e-6*pow(m,9)*pow(q,19) + 3.87284949567561e-6*pow(m,11)*pow(q,19) -
   						  5.16379932756748e-6*pow(m,13)*pow(q,19) + 3.87284949567561e-6*pow(m,15)*pow(q,19) -
   						   1.549139798270244e-6*pow(m,17)*pow(q,19) + 2.58189966378374e-7*pow(m,19)*pow(q,19) +
   							9.978286073219034e-9*pow(m,7)*pow(q,21) - 6.984800251253324e-8*pow(m,9)*pow(q,21) +
   							 2.095440075375997e-7*pow(m,11)*pow(q,21) - 3.4924001256266616e-7*pow(m,13)*pow(q,21) +
   							  3.4924001256266616e-7*pow(m,15)*pow(q,21) - 2.095440075375997e-7*pow(m,17)*pow(q,21) +
   							   6.984800251253324e-8*pow(m,19)*pow(q,21) - 9.978286073219034e-9*pow(m,21)*pow(q,21));
	double nu10 = (-1.2471407556959841 + 0.19666898195792082*pow(q,2) - 0.5659436240836003*pow(m,2)*pow(q,2) -
   	 0.001599333222875175*pow(q,4) - 0.09071321123054492*pow(m,2)*pow(q,4) +
   	  0.10415260296176723*pow(m,4)*pow(q,4) + 0.00010058664250254456*pow(q,6) +
   	   0.011108785507913294*pow(m,2)*pow(q,6) - 0.012796447213741879*pow(m,4)*pow(q,6) +
   		0.0014743655479493812*pow(m,6)*pow(q,6) + 1.284773045328378e-6*pow(q,8) +
   		 0.005136923118369371*pow(m,2)*pow(q,8) - 0.0009752447309743451*pow(m,4)*pow(q,8) -
   		  0.014230894116765334*pow(m,6)*pow(q,8) + 0.01007311288240538*pow(m,8)*pow(q,8) +
   		   0.0003186244656416432*pow(m,2)*pow(q,10) + 0.004620785151170218*pow(m,4)*pow(q,10) -
   			0.01577424580747269*pow(m,6)*pow(q,10) + 0.016412691466599565*pow(m,8)*pow(q,10) -
   			 0.005577855275938736*pow(m,10)*pow(q,10) + 5.503641182781793e-6*pow(m,2)*pow(q,12) +
   			  0.0009073197946210582*pow(m,4)*pow(q,12) - 0.003199742537207509*pow(m,6)*pow(q,12) +
   			   0.004085952109078974*pow(m,8)*pow(q,12) - 0.00222947399093356*pow(m,10)*pow(q,12) +
   				0.0004304409832582562*pow(m,12)*pow(q,12) + 0.0000606617765305937*pow(m,4)*pow(q,14) -
   				 0.0000889439437073267*pow(m,6)*pow(q,14) - 0.00025491864260893056*pow(m,8)*pow(q,14) +
   				  0.000691801824764815*pow(m,10)*pow(q,14) - 0.0005663808295264999*pow(m,12)*pow(q,14) +
   				   0.0001577798145473485*pow(m,14)*pow(q,14) + 1.3034295249181442e-6*pow(m,4)*pow(q,16) +
   					0.0000287255095173468*pow(m,6)*pow(q,16) - 0.00015960466662444264*pow(m,8)*pow(q,16) +
   					 0.00032509498082593973*pow(m,10)*pow(q,16) -
   					  0.00032446348077840346*pow(m,12)*pow(q,16) +
   					   0.00016061256084051544*pow(m,14)*pow(q,16) -
   						0.00003166833330587402*pow(m,16)*pow(q,16) + 2.7329976164226455e-6*pow(m,6)*pow(q,18) -
   						 0.00001442940389099287*pow(m,8)*pow(q,18) + 0.00003115205520862466*pow(m,10)*pow(q,18) -
   						  0.00003497413425302288*pow(m,12)*pow(q,18) +
   						   0.00002130914617090965*pow(m,14)*pow(q,18) - 6.555076660820858e-6*pow(m,16)*pow(q,18) +
   							7.644158088796423e-7*pow(m,18)*pow(q,18) + 7.254250259906868e-8*pow(m,6)*pow(q,20) -
   							 1.2736019228305327e-7*pow(m,8)*pow(q,20) - 7.592314008821227e-7*pow(m,10)*pow(q,20) +
   							  3.167572297689009e-6*pow(m,12)*pow(q,20) - 5.069758927241146e-6*pow(m,14)*pow(q,20) +
   							   4.183167334075971e-6*pow(m,16)*pow(q,20) - 1.7748264372690842e-6*pow(m,18)*pow(q,20) +
   								3.0789482331135884e-7*pow(m,20)*pow(q,20) + 3.3277883494903286e-8*pow(m,8)*pow(q,22) -
   								 2.32945184464323e-7*pow(m,10)*pow(q,22) + 6.98835553392969e-7*pow(m,12)*pow(q,22) -
   								  1.164725922321615e-6*pow(m,14)*pow(q,22) + 1.164725922321615e-6*pow(m,16)*pow(q,22) -
   								   6.98835553392969e-7*pow(m,18)*pow(q,22) + 2.32945184464323e-7*pow(m,20)*pow(q,22) -
   									3.3277883494903286e-8*pow(m,22)*pow(q,22) + 1.1229026917576822e-9*pow(m,8)*pow(q,24) -
   									 8.983221534061458e-9*pow(m,10)*pow(q,24) + 3.1441275369215106e-8*pow(m,12)*pow(q,24) -
   									  6.288255073843021e-8*pow(m,14)*pow(q,24) + 7.860318842303775e-8*pow(m,16)*pow(q,24) -
   									   6.288255073843021e-8*pow(m,18)*pow(q,24) + 3.1441275369215106e-8*pow(m,20)*pow(q,24) -
   										8.983221534061458e-9*pow(m,22)*pow(q,24) + 1.1229026917576822e-9*pow(m,24)*pow(q,24));
	double nu11 = (2.285742827552359*m*q + 0.7226382638570854*m*pow(q,3) -
   	 0.7476319175116637*pow(m,3)*pow(q,3) - 0.00219166062058118*m*pow(q,5) +
   	  0.03428619897006902*pow(m,3)*pow(q,5) - 0.03524771697697731*pow(m,5)*pow(q,5) +
   	   0.0018344802687887316*m*pow(q,7) + 0.030604698067605747*pow(m,3)*pow(q,7) -
   		0.06682450181910474*pow(m,5)*pow(q,7) + 0.03441042285434022*pow(m,7)*pow(q,7) +
   		 0.00009360184403011755*m*pow(q,9) + 0.011071904076376541*pow(m,3)*pow(q,9) -
   		  0.02747806657702974*pow(m,5)*pow(q,9) + 0.021096587482392123*pow(m,7)*pow(q,9) -
   		   0.004782909978302909*pow(m,9)*pow(q,9) + 1.2831863099243466e-6*m*pow(q,11) +
   			0.001488921228615683*pow(m,3)*pow(q,11) - 0.0010581702306936586*pow(m,5)*pow(q,11) -
   			 0.005988500755069408*pow(m,7)*pow(q,11) + 0.00919245971394552*pow(m,9)*pow(q,11) -
   			  0.0036359931431080624*pow(m,11)*pow(q,11) + 0.00007748573497747353*pow(m,3)*pow(q,13) +
   			   0.0009943001173283754*pow(m,5)*pow(q,13) - 0.00461388459790664*pow(m,7)*pow(q,13) +
   				0.007097324002161559*pow(m,9)*pow(q,13) - 0.004717623354804065*pow(m,11)*pow(q,13) +
   				 0.001162398098243298*pow(m,13)*pow(q,13) + 1.3116983763057892e-6*pow(m,3)*pow(q,15) +
   				  0.00016755004751823423*pow(m,5)*pow(q,15) - 0.0007138263179424906*pow(m,7)*pow(q,15) +
   				   0.001149564364880572*pow(m,9)*pow(q,15) - 0.000860911137813701*pow(m,11)*pow(q,15) +
   					0.000282814024906643*pow(m,13)*pow(q,15) - 0.00002650267992556324*pow(m,15)*pow(q,15) +
   					 9.99655361414996e-6*pow(m,5)*pow(q,17) - 0.000017083796878677634*pow(m,7)*pow(q,17) -
   					  0.00006416191104300251*pow(m,9)*pow(q,17) + 0.0002275545406757872*pow(m,11)*pow(q,17) -
   					   0.00027680249119481956*pow(m,13)*pow(q,17) +
   						0.00015302866724277602*pow(m,15)*pow(q,17) -
   						 0.00003253156241621345*pow(m,17)*pow(q,17) + 2.0199397993704069e-7*pow(m,5)*pow(q,19) +
   						  4.483785295955416e-6*pow(m,7)*pow(q,19) - 0.00003004967541389721*pow(m,9)*pow(q,19) +
   						   0.00007592180833235841*pow(m,11)*pow(q,19) -
   							0.00009993597440736625*pow(m,13)*pow(q,19) +
  							  0.00007327517434891132*pow(m,15)*pow(q,19) -
   							   0.000028497951370963253*pow(m,17)*pow(q,19) + 4.600839235064527e-6*pow(m,19)*pow(q,19) +
   								3.8314711583873637e-7*pow(m,7)*pow(q,21) - 2.306693885284017e-6*pow(m,9)*pow(q,21) +
   								 5.794073879090637e-6*pow(m,11)*pow(q,21) - 7.780110170548704e-6*pow(m,13)*pow(q,21) +
   								  5.903430542613016e-6*pow(m,15)*pow(q,21) - 2.416050548806396e-6*pow(m,17)*pow(q,21) +
   								   4.300142573483274e-7*pow(m,19)*pow(q,21) - 7.811190251598497e-9*pow(m,21)*pow(q,21) +
   									9.302550380417304e-9*pow(m,7)*pow(q,23) - 1.6690592987003388e-8*pow(m,9)*pow(q,23) -
   									 1.436372597426608e-7*pow(m,11)*pow(q,23) + 6.913831898796669e-7*pow(m,13)*pow(q,23) -
   									  1.3693648253425153e-6*pow(m,15)*pow(q,23) + 1.4996005306683577e-6*pow(m,17)*pow(q,23) -
   									   9.518546005313515e-7*pow(m,19)*pow(q,23) + 3.296882673510069e-7*pow(m,21)*pow(q,23) -
   										4.8427259675917734e-8*pow(m,23)*pow(q,23) + 4.327893006154886e-9*pow(m,9)*pow(q,25) -
   										 3.462314404923909e-8*pow(m,11)*pow(q,25) + 1.2118100417233681e-7*pow(m,13)*pow(q,25) -
   										  2.4236200834467363e-7*pow(m,15)*pow(q,25) + 3.02952510430842e-7*pow(m,17)*pow(q,25) -
   										   2.4236200834467363e-7*pow(m,19)*pow(q,25) + 1.2118100417233681e-7*pow(m,21)*pow(q,25) -
   											3.462314404923909e-8*pow(m,23)*pow(q,25) + 4.327893006154886e-9*pow(m,25)*pow(q,25) +
   											 1.296055738409144e-10*pow(m,9)*pow(q,27) - 1.1664501645682296e-9*pow(m,11)*pow(q,27) +
   											  4.665800658272918e-9*pow(m,13)*pow(q,27) - 1.0886868202636808e-8*pow(m,15)*pow(q,27) +
   											   1.6330302303955214e-8*pow(m,17)*pow(q,27) - 1.6330302303955214e-8*pow(m,19)*pow(q,27) +
   												1.0886868202636808e-8*pow(m,21)*pow(q,27) - 4.665800658272918e-9*pow(m,23)*pow(q,27) +
   												 1.1664501645682296e-9*pow(m,25)*pow(q,27) - 1.296055738409144e-10*pow(m,27)*pow(q,27));
	double nu12 = (-2.3214691102682137 + 0.4233332887205478*pow(q,2) - 1.6535156044744763*pow(m,2)*pow(q,2) -
   	 0.009575619252781869*pow(q,4) - 0.49908235795384864*pow(m,2)*pow(q,4) +
   	  0.5605565168253283*pow(m,4)*pow(q,4) + 0.0002470037374041981*pow(q,6) -
   	   0.055879650896440175*pow(m,2)*pow(q,6) + 0.14582215694883663*pow(m,4)*pow(q,6) -
   		0.09163431302063531*pow(m,6)*pow(q,6) + 7.160197575856309e-6*pow(q,8) +
   		 0.006708444441198345*pow(m,2)*pow(q,8) + 0.01200732683845157*pow(m,4)*pow(q,8) -
   		  0.04599222621333539*pow(m,6)*pow(q,8) + 0.027296885142868994*pow(m,8)*pow(q,8) +
   		   8.22482351741291e-8*pow(q,10) + 0.0009479904980529705*pow(m,2)*pow(q,10) +
   			0.011795594862412497*pow(m,4)*pow(q,10) - 0.03960866673099132*pow(m,6)*pow(q,10) +
   			 0.039984402815029504*pow(m,8)*pow(q,10) - 0.013119321584604212*pow(m,10)*pow(q,10) +
   			  0.000041747367026636167*pow(m,2)*pow(q,12) + 0.0032656636590789098*pow(m,4)*pow(q,12) -
   			   0.010562578289446673*pow(m,6)*pow(q,12) + 0.011529680198964425*pow(m,8)*pow(q,12) -
   				0.004641527797350074*pow(m,10)*pow(q,12) + 0.0003670148617267761*pow(m,12)*pow(q,12) +
   				 5.79102211792999e-7*pow(m,2)*pow(q,14) + 0.00035632335836380887*pow(m,4)*pow(q,14) -
   				  0.00039702105798530065*pow(m,6)*pow(q,14) - 0.001997777836385748*pow(m,8)*pow(q,14) +
   				   0.00480348924636569*pow(m,10)*pow(q,14) - 0.003810652156722736*pow(m,12)*pow(q,14) +
   					0.0010450593441524928*pow(m,14)*pow(q,14) + 0.000016619643156119418*pow(m,4)*pow(q,16) +
   					 0.00019804208146357232*pow(m,6)*pow(q,16) - 0.0011730169804747388*pow(m,8)*pow(q,16) +
   					  0.0023773097322354886*pow(m,10)*pow(q,16) - 0.002323542775835549*pow(m,12)*pow(q,16) +
   					   0.0011172960106090083*pow(m,14)*pow(q,16) - 0.00021270771115390103*pow(m,16)*pow(q,16) +
   						2.7080209022507177e-7*pow(m,4)*pow(q,18) + 0.000029823255131513415*pow(m,6)*pow(q,18) -
   						 0.00014883012158573418*pow(m,8)*pow(q,18) + 0.0002867019165005833*pow(m,10)*pow(q,18) -
   						  0.00026444962593313734*pow(m,12)*pow(q,18) +
   						   0.00010840292378582608*pow(m,14)*pow(q,18) - 6.949173266558589e-6*pow(m,16)*pow(q,18) -
   							4.969976722717802e-6*pow(m,18)*pow(q,18) + 1.6160402213534573e-6*pow(m,6)*pow(q,20) -
   							 3.1931568739504298e-6*pow(m,8)*pow(q,20) - 0.000014403807312348288*pow(m,10)*pow(q,20) +
   							  0.00006335498192362688*pow(m,12)*pow(q,20) -
   							   0.00010208012483938706*pow(m,14)*pow(q,20) +
   								0.00008410906456071661*pow(m,16)*pow(q,20) -
   								 0.00003553198604180975*pow(m,18)*pow(q,20) + 6.128988361798578e-6*pow(m,20)*pow(q,20) +
   								  3.066431538069726e-8*pow(m,6)*pow(q,22) + 6.927626521816531e-7*pow(m,8)*pow(q,22) -
   								   5.440789660149455e-6*pow(m,10)*pow(q,22) + 0.000016379520603762977*pow(m,12)*pow(q,22) -
   									0.000026678953019579707*pow(m,14)*pow(q,22) +
   									 0.000025772504756001275*pow(m,16)*pow(q,22) -
   									  0.000014833773812387754*pow(m,18)*pow(q,22) +
   									   4.7183272888552094e-6*pow(m,20)*pow(q,22) - 6.402631240648955e-7*pow(m,22)*pow(q,22) +
   										5.3951488021981545e-8*pow(m,8)*pow(q,24) - 3.6418937138651074e-7*pow(m,10)*pow(q,24) +
   										 1.038683935090092e-6*pow(m,12)*pow(q,24) - 1.605410140654793e-6*pow(m,14)*pow(q,24) +
   										  1.4168155139117526e-6*pow(m,16)*pow(q,24) - 6.61494681604011e-7*pow(m,18)*pow(q,24) +
   										   9.476847603930991e-8*pow(m,20)*pow(q,24) + 4.034582534953875e-8*pow(m,22)*pow(q,24) -
   											1.347104476736004e-8*pow(m,24)*pow(q,24) + 1.2053229494811545e-9*pow(m,8)*pow(q,26) -
   											 2.2110951847400125e-9*pow(m,10)*pow(q,26) - 2.570286470340147e-8*pow(m,12)*pow(q,26) +
   											  1.4058359034011364e-7*pow(m,14)*pow(q,26) - 3.317907445584358e-7*pow(m,16)*pow(q,26) +
   											   4.527061036067011e-7*pow(m,18)*pow(q,26) - 3.8241430843664427e-7*pow(m,20)*pow(q,26) +
   												1.9843909191520907e-7*pow(m,22)*pow(q,26) - 5.824658433939265e-8*pow(m,24)*pow(q,26) +
   												 7.431488411109225e-9*pow(m,26)*pow(q,26) + 5.668134623555751e-10*pow(m,10)*pow(q,28) -
   												  5.1013211612001755e-9*pow(m,12)*pow(q,28) + 2.0405284644800702e-8*pow(m,14)*pow(q,28) -
   												   4.761233083786831e-8*pow(m,16)*pow(q,28) + 7.141849625680247e-8*pow(m,18)*pow(q,28) -
   													7.141849625680247e-8*pow(m,20)*pow(q,28) + 4.761233083786831e-8*pow(m,22)*pow(q,28) -
   													 2.0405284644800702e-8*pow(m,24)*pow(q,28) + 5.1013211612001755e-9*pow(m,26)*pow(q,28) -
   													  5.668134623555751e-10*pow(m,28)*pow(q,28) + 1.5258273928639506e-11*pow(m,10)*pow(q,30) -
   													   1.5258273928639504e-10*pow(m,12)*pow(q,30) + 6.866223267887776e-10*pow(m,14)*pow(q,30) -
   														1.8309928714367404e-9*pow(m,16)*pow(q,30) + 3.204237525014296e-9*pow(m,18)*pow(q,30) -
   														 3.845085030017155e-9*pow(m,20)*pow(q,30) + 3.204237525014296e-9*pow(m,22)*pow(q,30) -
   														  1.8309928714367404e-9*pow(m,24)*pow(q,30) + 6.866223267887776e-10*pow(m,26)*pow(q,30) -
   														   1.5258273928639504e-10*pow(m,28)*pow(q,30) + 1.5258273928639506e-11*pow(m,30)*pow(q,30));

	return 1. + (nu2*pow(eps, 2) + nu3*pow(eps, 3) + nu4*pow(eps, 4) + nu5*pow(eps, 5) + nu6*pow(eps, 6)
		+ nu7*pow(eps, 7) + nu8*pow(eps, 8) + nu9*pow(eps, 9) + nu10*pow(eps, 10) + nu11*pow(eps, 11) + nu12*pow(eps, 12));
}

double nu_solver_low_freq_s0_l2(double q, int m, double eps){
	double nu2 = -(79./210);
	double nu3 = (3*m*q)/70.;
	double nu4 = (-1416494 + 189000*pow(q,2) - 123900*pow(m,2)*pow(q,2))/1.8522e7;
	double nu5 = (0.03337083030921832*m*q - 0.0032079709520939177*m*pow(q,3) +
   	 0.0009597720944918051*pow(m,3)*pow(q,3) - 0.0003281762775433661*m*pow(q,5) +
   	  0.0004102203469292077*pow(m,3)*pow(q,5) - 0.00008204406938584153*pow(m,5)*pow(q,5));
	double nu6 = (-0.038224220512750595 + 0.011547916289123803*pow(q,2) -
   	 0.011291447625807173*pow(m,2)*pow(q,2) - 0.00021231821338844257*pow(q,4) +
   	  0.00040651567318123946*pow(m,2)*pow(q,4) - 0.000028872969969697685*pow(m,4)*pow(q,4) -
   	   0.00003738717085937083*pow(m,2)*pow(q,6) + 0.00004673396357421353*pow(m,4)*pow(q,6) -
   		9.346792714842707e-6*pow(m,6)*pow(q,6));
	double nu7 = (0.03097310767493913*m*q - 0.005998684278133017*m*pow(q,3) +
   	 0.00212201393330343*pow(m,3)*pow(q,3) - 0.0010442008135466826*m*pow(q,5) +
   	  0.0013580816972980424*pow(m,3)*pow(q,5) - 0.00027961582644451616*pow(m,5)*pow(q,5) -
   	   9.831912215789444e-6*m*pow(q,7) + 0.0000175869755025765*pow(m,3)*pow(q,7) -
   		9.07933459499698e-6*pow(m,5)*pow(q,7) + 1.3242713082099238e-6*pow(m,7)*pow(q,7));
	double nu8 = (-0.026817650745141034 + 0.012250703102111085*pow(q,2) -
   	 0.015279329269920525*pow(m,2)*pow(q,2) - 0.00043265751190068924*pow(q,4) +
   	  0.0007511976084729975*pow(m,2)*pow(q,4) - 4.690602901717264e-6*pow(m,4)*pow(q,4) -
   	   3.2977008790858414e-6*pow(q,6) - 0.00011432952719898795*pow(m,2)*pow(q,6) +
   		0.0001515402560477707*pow(m,4)*pow(q,6) - 0.00003135066877123091*pow(m,6)*pow(q,6) +
   		 5.035681698138224e-7*pow(m,2)*pow(q,8) - 3.770653205615563e-8*pow(m,4)*pow(q,8) -
   		  6.138000578104473e-7*pow(m,6)*pow(q,8) + 1.479384200527806e-7*pow(m,8)*pow(q,8) +
   		   1.4314512987257154e-7*pow(m,2)*pow(q,10) - 3.578628246814289e-7*pow(m,4)*pow(q,10) +
   			2.952368303621788e-7*pow(m,6)*pow(q,10) - 8.946570617035722e-8*pow(m,8)*pow(q,10) +
   			 8.946570617035721e-9*pow(m,10)*pow(q,10));
	double nu9 = (0.031288716454982955*m*q - 0.008459848325776249*m*pow(q,3) +
   	 0.0036248472218314673*pow(m,3)*pow(q,3) - 0.001535991677630517*m*pow(q,5) +
   	  0.002038542353687849*pow(m,3)*pow(q,5) - 0.0004254839566205436*pow(m,5)*pow(q,5) -
   	   0.000035670447514678786*m*pow(q,7) + 0.00006407165666180965*pow(m,3)*pow(q,7) -
   		0.00003406180896228998*pow(m,5)*pow(q,7) + 5.1117746314102095e-6*pow(m,7)*pow(q,7) +
   		 1.6086573590508708e-8*m*pow(q,9) + 4.877125462701701e-8*pow(m,3)*pow(q,9) -
   		  1.351410184216712e-7*pow(m,5)*pow(q,9) + 8.35490207792348e-8*pow(m,7)*pow(q,9) -
   		   1.3265830575089315e-8*pow(m,9)*pow(q,9) + 4.892301907037255e-8*pow(m,3)*pow(q,11) -
   			1.223075476759314e-7*pow(m,5)*pow(q,11) + 1.0090372683264339e-7*pow(m,7)*pow(q,11) -
   			 3.057688691898285e-8*pow(m,9)*pow(q,11) + 3.0576886918982846e-9*pow(m,11)*pow(q,11));
	double nu10 = (-0.022443371129554832 + 0.0129516076087242*pow(q,2) -
   	 0.019312473642693243*pow(m,2)*pow(q,2) - 0.0006677049525514774*pow(q,4) +
   	  0.0015504227145146483*pow(m,2)*pow(q,4) - 0.00020241245465524776*pow(m,4)*pow(q,4) -
   	   0.000012343387381072892*pow(q,6) - 0.0001144367131740799*pow(m,2)*pow(q,6) +
   		0.00018441513538279157*pow(m,4)*pow(q,6) - 0.000041730507906934076*pow(m,6)*pow(q,6) +
   		 1.075186966462564e-9*pow(q,8) + 3.4709339137151903e-6*pow(m,2)*pow(q,8) -
   		  2.8400370886851643e-6*pow(m,4)*pow(q,8) - 1.1283310292202208e-6*pow(m,6)*pow(q,8) +
   		   4.159089991945201e-7*pow(m,8)*pow(q,8) + 9.645019464777811e-7*pow(m,2)*pow(q,10) -
   			2.3982654743099597e-6*pow(m,4)*pow(q,10) + 1.9547060711576652e-6*pow(m,6)*pow(q,10) -
   			 5.766383015810625e-7*pow(m,8)*pow(q,10) + 5.569575825557584e-8*pow(m,10)*pow(q,10) +
   			  1.2459818201974027e-8*pow(m,2)*pow(q,12) - 3.088458242008442e-8*pow(m,4)*pow(q,12) +
   			   2.503596732944482e-8*pow(m,6)*pow(q,12) - 7.240900013729314e-9*pow(m,8)*pow(q,12) +
   				6.131367095917242e-10*pow(m,10)*pow(q,12) + 1.656019280316526e-11*pow(m,12)*pow(q,12));
	double nu11 = (0.03296832711715091*m*q - 0.011062907380531905*m*pow(q,3) +
   	 0.00584499535725277*pow(m,3)*pow(q,3) - 0.0016840837158371693*m*pow(q,5) +
   	  0.0022399515259529943*pow(m,3)*pow(q,5) - 0.00046997228141925013*pow(m,5)*pow(q,5) -
   	   0.00005510682910221688*m*pow(q,7) + 0.00010352219237277097*pow(m,3)*pow(q,7) -
   		0.0000587624768938038*pow(m,5)*pow(q,7) + 9.183923080651106e-6*pow(m,7)*pow(q,7) -
   		 2.8076846369874315e-8*m*pow(q,9) + 1.0661058953559986e-6*pow(m,3)*pow(q,9) -
   		  1.8350630069025238e-6*pow(m,5)*pow(q,9) + 9.102317304537478e-7*pow(m,7)*pow(q,9) -
   		   1.2977886497116605e-7*pow(m,9)*pow(q,9) + 9.559334567845373e-9*m*pow(q,11) +
   			2.81354960445681e-7*pow(m,3)*pow(q,11) - 7.596888805900565e-7*pow(m,5)*pow(q,11) +
   			 6.485551561542219e-7*pow(m,7)*pow(q,11) - 2.0001017863382373e-7*pow(m,9)*pow(q,11) +
   			  2.0229608056131842e-8*pow(m,11)*pow(q,11) + 2.06366076220171e-9*pow(m,3)*pow(q,13) -
   			   6.676370311771543e-9*pow(m,5)*pow(q,13) + 8.049346337709196e-9*pow(m,7)*pow(q,13) -
   				4.419050939302308e-9*pow(m,9)*pow(q,13) + 1.0772403015546491e-9*pow(m,11)*pow(q,13) -
   				 9.482615039170424e-11*pow(m,13)*pow(q,13) - 1.248751333254287e-10*pow(m,3)*pow(q,15) +
   				  4.682817499703577e-10*pow(m,5)*pow(q,15) - 6.790085374570186e-10*pow(m,7)*pow(q,15) +
   				   4.780376197614069e-10*pow(m,9)*pow(q,15) - 1.6975213436425466e-10*pow(m,11)*pow(q,15) +
   					2.9267609373147355e-11*pow(m,13)*pow(q,15) - 1.9511739582098236e-12*pow(m,15)*pow(q,15));
	double nu12 = (-0.020681283934631207 + 0.013774583436434614*pow(q,2) -
   	 0.023938561866459717*pow(m,2)*pow(q,2) - 0.0009812269288954113*pow(q,4) +
   	  0.0031294830247713167*pow(m,2)*pow(q,4) - 0.0008114000313445158*pow(m,4)*pow(q,4) -
   	   0.000010589048582507119*pow(q,6) - 0.000025895191163384186*pow(m,2)*pow(q,6) +
   		0.0000708864028708062*pow(m,4)*pow(q,6) - 0.00002033471571115035*pow(m,6)*pow(q,6) -
   		 8.864506707258458e-8*pow(q,8) + 8.117694489178865e-6*pow(m,2)*pow(q,8) -
   		  7.532526134191691e-6*pow(m,4)*pow(q,8) - 1.765481385849613e-6*pow(m,6)*pow(q,8) +
   		   8.477090395428887e-7*pow(m,8)*pow(q,8) + 2.7134988081085585e-9*pow(q,10) +
   			2.8705488755060215e-6*pow(m,2)*pow(q,10) - 7.058827776925073e-6*pow(m,4)*pow(q,10) +
   			 5.671859057174495e-6*pow(m,6)*pow(q,10) - 1.638035819272373e-6*pow(m,8)*pow(q,10) +
   			  1.5390158227700343e-7*pow(m,10)*pow(q,10) + 8.706725398285992e-8*pow(m,2)*pow(q,12) -
   			   2.2052878214537159e-7*pow(m,4)*pow(q,12) + 1.8447384397230875e-7*pow(m,6)*pow(q,12) -
   				5.6523005646327777e-8*pow(m,8)*pow(q,12) + 5.4453224742197205e-9*pow(m,10)*pow(q,12) +
   				 6.536736231095242e-11*pow(m,12)*pow(q,12) + 3.716274327755884e-10*pow(m,2)*pow(q,14) -
   				  1.4194229447253746e-9*pow(m,4)*pow(q,14) + 1.9534125728839646e-9*pow(m,6)*pow(q,14) -
   				   1.1462357332774618e-9*pow(m,8)*pow(q,14) + 2.4935368079023007e-10*pow(m,10)*pow(q,14) -
   					6.300326310590689e-12*pow(m,12)*pow(q,14) - 2.4346821363559534e-12*pow(m,14)*pow(q,14) -
   					 7.113140505878851e-11*pow(m,4)*pow(q,16) + 2.6674276897045687e-10*pow(m,6)*pow(q,16) -
   					  3.867770150071625e-10*pow(m,8)*pow(q,16) + 2.7229990999067477e-10*pow(m,10)*pow(q,16) -
   					   9.669425375179062e-11*pow(m,12)*pow(q,16) + 1.6671423060653555e-11*pow(m,14)*pow(q,16) -
   						1.1114282040435705e-12*pow(m,16)*pow(q,16));

	return 2. + (nu2*pow(eps, 2) + nu3*pow(eps, 3) + nu4*pow(eps, 4) + nu5*pow(eps, 5) + nu6*pow(eps, 6)
		+ nu7*pow(eps, 7) + nu8*pow(eps, 8) + nu9*pow(eps, 9) + nu10*pow(eps, 10) + nu11*pow(eps, 11) + nu12*pow(eps, 12));
}

double nu_solver_low_freq_s0_l3(double q, int m, double eps){
	double nu2 = -(169./630);
	double nu3 = (19*m*q)/1260.;
	double nu4 = (-148760842 + 18151560*pow(q,2) - 5927040*pow(m,2)*pow(q,2))/5.501034e9;
	double nu5 = (0.005640185917774731*m*q - 0.0005043824488268933*m*pow(q,3) +
   	 0.00008751536529314307*pow(m,3)*pow(q,3));
	double nu6 = (-0.006461623912280694 + 0.001765178127355456*pow(q,2) -
   	 0.0008945126785437488*pow(m,2)*pow(q,2) - 0.000044888652520526816*pow(q,4) +
   	  0.0000653345324501329*pow(m,2)*pow(q,4) - 8.11765447006749e-6*pow(m,4)*pow(q,4));
	double nu7 = (0.002494276058961166*m*q - 0.0005586197410324905*m*pow(q,3) +
   	 0.00013035767824169948*pow(m,3)*pow(q,3) + 0.000013027942816090196*m*pow(q,5) -
   	  7.156626308736202e-6*pow(m,3)*pow(q,5) + 6.010482487023779e-7*pow(m,5)*pow(q,5) +
   	   8.453085376162299e-7*m*pow(q,7) - 1.150558842866535e-6*pow(m,3)*pow(q,7) +
   		3.287310979618672e-7*pow(m,5)*pow(q,7) - 2.3480792711561942e-8*pow(m,7)*pow(q,7));
	double nu8 = (-0.00212149594160692 + 0.0009269789023106594*pow(q,2) -
   	 0.0006418449425535479*pow(m,2)*pow(q,2) - 0.0000747390170490615*pow(q,4) +
   	  0.00012436509658936737*pow(m,2)*pow(q,4) - 0.000018526935005462478*pow(m,4)*pow(q,4) +
   	   5.223303629929774e-7*pow(q,6) - 1.3806028759546738e-6*pow(m,2)*pow(q,6) +
   		2.8425836733942186e-7*pow(m,4)*pow(q,6) - 9.648987340334575e-10*pow(m,6)*pow(q,6) +
   		 3.6888092129817866e-8*pow(m,2)*pow(q,8) - 5.020879206558542e-8*pow(m,4)*pow(q,8) +
   		  1.4345369161595836e-8*pow(m,6)*pow(q,8) - 1.024669225828274e-9*pow(m,8)*pow(q,8));
	double nu9 = (0.0012403776203171867*m*q - 0.0004736895685836144*m*pow(q,3) +
   	 0.0001381048896540659*pow(m,3)*pow(q,3) + 0.000028728249299461554*m*pow(q,5) -
   	  0.000017288641454560554*pow(m,3)*pow(q,5) + 1.5770279084053183e-6*pow(m,5)*pow(q,5) +
   	   2.7342350821210253e-6*m*pow(q,7) - 3.805234582699139e-6*pow(m,3)*pow(q,7) +
   		1.1101875794861178e-6*pow(m,5)*pow(q,7) - 8.05665732863411e-8*pow(m,7)*pow(q,7) +
   		 1.112824811808274e-8*m*pow(q,9) - 1.8792921706723137e-8*pow(m,3)*pow(q,9) +
   		  9.290453094642478e-9*pow(m,5)*pow(q,9) - 1.7270611600576227e-9*pow(m,7)*pow(q,9) +
   		   1.0128165405553904e-10*pow(m,9)*pow(q,9));
	double nu10 = (-0.0008420801443405631 + 0.0005182379760872247*pow(q,2) -
   	 0.00045312631858195843*pow(m,2)*pow(q,2) - 0.00007724793611333435*pow(q,4) +
   	  0.000148916821922478*pow(m,2)*pow(q,4) - 0.000025567502190973252*pow(m,4)*pow(q,4) +
   	   1.258023411052832e-6*pow(q,6) - 3.142051306035189e-6*pow(m,2)*pow(q,6) +
   		5.688791691529919e-7*pow(m,4)*pow(q,6) + 1.7396639765356648e-8*pow(m,6)*pow(q,6) +
   		 5.527179307795095e-9*pow(q,8) + 1.3108147337997035e-7*pow(m,2)*pow(q,8) -
   		  1.9364260725973623e-7*pow(m,4)*pow(q,8) + 5.709139069753163e-8*pow(m,6)*pow(q,8) -
   		   4.175886723732425e-9*pow(m,8)*pow(q,8) - 3.6884636610224486e-10*pow(m,2)*pow(q,10) +
   			3.4000909606116023e-10*pow(m,4)*pow(q,10) + 7.710301778091962e-11*pow(m,6)*pow(q,10) -
   			 5.2766630826878384e-11*pow(m,8)*pow(q,10) + 4.500883087043387e-12*pow(m,10)*pow(q,10));
	double nu11 = (0.0006771272333630794*m*q - 0.000370696257955747*m*pow(q,3) +
   	 0.00012799181837228467*pow(m,3)*pow(q,3) + 0.00003758917351116027*m*pow(q,5) -
   	  0.000025412249306067784*pow(m,3)*pow(q,5) + 2.6344763383580005e-6*pow(m,5)*pow(q,5) +
   	   3.8262721643824205e-6*m*pow(q,7) - 5.427184205687748e-6*pow(m,3)*pow(q,7) +
   		1.613605914885263e-6*pow(m,5)*pow(q,7) - 1.1892225414625859e-7*pow(m,7)*pow(q,7) +
   		 3.65381230595975e-8*m*pow(q,9) - 6.085865705193036e-8*pow(m,3)*pow(q,9) +
   		  3.02579678493725e-8*pow(m,5)*pow(q,9) - 5.743233430423714e-9*pow(m,7)*pow(q,9) +
   		   3.4463600415654533e-10*pow(m,9)*pow(q,9) - 3.722125196935459e-11*m*pow(q,11) +
   			9.568498345956707e-11*pow(m,3)*pow(q,11) - 7.816786135396428e-11*pow(m,5)*pow(q,11) +
   			 2.182576121454718e-11*pow(m,7)*pow(q,11) - 2.188631361039012e-12*pow(m,9)*pow(q,11) +
   			  6.700001024363536e-14*pow(m,11)*pow(q,11));
	double nu12 = (-0.0003845774226816963 + 0.0003085712877624197*pow(q,2) -
   	 0.0003222817069907867*pow(m,2)*pow(q,2) - 0.00006685653744969744*pow(q,4) +
   	  0.00014642215046358427*pow(m,2)*pow(q,4) - 0.000028326657178750664*pow(m,4)*pow(q,4) +
   	   1.8521457788424656e-6*pow(q,6) - 5.281015268956064e-6*pow(m,2)*pow(q,6) +
   		1.3881013520102658e-6*pow(m,4)*pow(q,6) - 3.874248499096748e-8*pow(m,6)*pow(q,6) +
   		 1.994174863961969e-8*pow(q,8) + 1.0008289244173893e-7*pow(m,2)*pow(q,8) -
   		  1.9901344039854652e-7*pow(m,4)*pow(q,8) + 6.404630416914848e-8*pow(m,6)*pow(q,8) -
   		   4.953290990110424e-9*pow(m,8)*pow(q,8) - 5.563863406627018e-11*pow(q,10) -
   			1.5159873903554614e-9*pow(m,2)*pow(q,10) + 1.9437433586902602e-9*pow(m,4)*pow(q,10) -
   			 2.6026554363579613e-10*pow(m,6)*pow(q,10) - 6.623086742352842e-11*pow(m,8)*pow(q,10) +
   			  8.184253699536045e-12*pow(m,10)*pow(q,10) + 5.438501323315185e-12*pow(m,2)*pow(q,12) -
   			   5.364195325903459e-12*pow(m,4)*pow(q,12) - 1.4910508469426654e-12*pow(m,6)*pow(q,12) +
   				1.7737319453083185e-12*pow(m,8)*pow(q,12) - 3.800924936297998e-13*pow(m,10)*pow(q,12) +
   				 2.3105397852419587e-14*pow(m,12)*pow(q,12) + 1.3318470709264503e-12*pow(m,2)*pow(q,14) -
   				  3.625583693077559e-12*pow(m,4)*pow(q,14) + 3.503292179620578e-12*pow(m,6)*pow(q,14) -
   				   1.4839407179149647e-12*pow(m,8)*pow(q,14) + 3.0213197442312995e-13*pow(m,10)*pow(q,14) -
   					2.8774473754583803e-14*pow(m,12)*pow(q,14) + 1.0276597769494215e-15*pow(m,14)*pow(q,14));

	return 3. + (nu2*pow(eps, 2) + nu3*pow(eps, 3) + nu4*pow(eps, 4) + nu5*pow(eps, 5) + nu6*pow(eps, 6)
		+ nu7*pow(eps, 7) + nu8*pow(eps, 8) + nu9*pow(eps, 9) + nu10*pow(eps, 10) + nu11*pow(eps, 11) + nu12*pow(eps, 12));
}

double nu_solver_low_freq_s0_l4(double q, int m, double eps){
	double nu2 = -(289./1386);
	double nu3 = (97*m*q)/13860.;
	double nu4 = (-4354548790 + 514078488*pow(q,2) - 100740024*pow(m,2)*pow(q,2))/3.4612505928e11;
	double nu5 = (0.0015465100949740168*m*q - 0.0001273638178400083*m*pow(q,3) +
   	 0.000013192275097037001*pow(m,3)*pow(q,3));
	double nu6 = (-0.001774660670779825 + 0.00046245475933980805*pow(q,2) -
   	 0.00014074974097612523*pow(m,2)*pow(q,2) - 0.000010005308831124607*pow(q,4) +
   	  8.693080285287481e-6*pow(m,2)*pow(q,4) - 6.376681311068338e-7*pow(m,4)*pow(q,4));
	double nu7 = (0.00039706691014570106*m*q - 0.00008106862327994429*m*pow(q,3) +
   	 0.00001135206886268987*pow(m,3)*pow(q,3) + 1.661676939939185e-6*m*pow(q,5) -
   	  5.889513498478437e-7*pow(m,3)*pow(q,5) + 3.3730767977967904e-8*pow(m,5)*pow(q,5));
	double nu8 = (-0.0003380130059294916 + 0.0001381283622642085*pow(q,2) -
   	 0.00005747576181182069*pow(m,2)*pow(q,2) - 9.620986987008067e-6*pow(q,4) +
   	  9.650564497264167e-6*pow(m,2)*pow(q,4) - 8.744383638207239e-7*pow(m,4)*pow(q,4) +
   	   9.180000325587916e-8*pow(q,6) - 1.898810198411448e-7*pow(m,2)*pow(q,6) +
   		4.2209261850448685e-8*pow(m,4)*pow(q,6) - 2.0139784467701803e-9*pow(m,6)*pow(q,6));
	double nu9 = (0.00011219511709157314*m*q - 0.0000387045758171058*m*pow(q,3) +
   	 6.856771050272362e-6*pow(m,3)*pow(q,3) + 2.431431161591878e-6*m*pow(q,5) -
   	  1.0059761944614068e-6*pow(m,3)*pow(q,5) + 6.752163979409028e-8*pow(m,5)*pow(q,5) -
   	   2.3224030302649908e-8*m*pow(q,7) + 1.7013489495611953e-8*pow(m,3)*pow(q,7) -
   		2.4908545665482975e-9*pow(m,5)*pow(q,7) + 8.870163090529585e-11*pow(m,7)*pow(q,7) -
   		 1.095987725933822e-9*m*pow(q,9) + 1.5602603042807881e-9*pow(m,3)*pow(q,9) -
   		  5.194525159373844e-10*pow(m,5)*pow(q,9) + 5.708269405905323e-11*pow(m,7)*pow(q,9) -
   		   1.9027564686351076e-12*pow(m,9)*pow(q,9));
	double nu10 = (-0.00007609350885288868 + 0.00004294837718441508*pow(q,2) -
   	 0.000022713962354167934*pow(m,2)*pow(q,2) - 5.723926549660549e-6*pow(q,4) +
   	  6.855555480778363e-6*pow(m,2)*pow(q,4) - 7.453380543728812e-7*pow(m,4)*pow(q,4) +
   	   1.8537867882336953e-7*pow(q,6) - 4.053697787771432e-7*pow(m,2)*pow(q,6) +
   		1.0072607558114304e-7*pow(m,4)*pow(q,6) - 5.342661361836112e-9*pow(m,6)*pow(q,6) -
   		 6.563566061875882e-10*pow(q,8) + 2.082970743960046e-9*pow(m,2)*pow(q,8) -
   		  6.415060143667802e-10*pow(m,4)*pow(q,8) + 3.396179972915292e-11*pow(m,6)*pow(q,8) +
   		   6.6337420653362e-13*pow(m,8)*pow(q,8) - 3.678574720262309e-11*pow(m,2)*pow(q,10) +
   			5.23685984481787e-11*pow(m,4)*pow(q,10) - 1.743491143457657e-11*pow(m,6)*pow(q,10) +
   			 1.9159243334699525e-12*pow(m,8)*pow(q,10) - 6.386414444899843e-14*pow(m,10)*pow(q,10));
	double nu11 = (0.00003406721356000823*m*q - 0.000017051665569887862*m*pow(q,3) +
   	 3.6599839119602858e-6*pow(m,3)*pow(q,3) + 2.0899609566013465e-6*m*pow(q,5) -
   	  9.99449062188068e-7*pow(m,3)*pow(q,5) + 7.712841807030648e-8*pow(m,5)*pow(q,5) -
   	   5.6639334242094956e-8*m*pow(q,7) + 4.391317340363255e-8*pow(m,3)*pow(q,7) -
   		6.836614874105728e-9*pow(m,5)*pow(q,7) + 2.5679513233739946e-10*pow(m,7)*pow(q,7) -
   		 3.5673471710952977e-9*m*pow(q,9) + 5.1477875092018095e-9*pow(m,3)*pow(q,9) -
   		  1.739595267449129e-9*pow(m,5)*pow(q,9) + 1.9360920945019364e-10*pow(m,7)*pow(q,9) -
   		   6.517456977236887e-12*pow(m,9)*pow(q,9) - 8.222616639407466e-12*m*pow(q,11) +
   			1.3664618533957099e-11*pow(m,3)*pow(q,11) - 6.685761534693834e-12*pow(m,5)*pow(q,11) +
   			 1.3566556648428198e-12*pow(m,7)*pow(q,11) - 1.1629673671891181e-13*pow(m,9)*pow(q,11) +
   			  3.4007120202943138e-15*pow(m,11)*pow(q,11));
	double nu12 = (-0.000019223974708540905 + 0.000014055704732401252*pow(q,2) -
   	 9.02698499525213e-6*pow(m,2)*pow(q,2) - 2.9245595553741652e-6*pow(q,4) +
   	  4.11931106599532e-6*pow(m,2)*pow(q,4) - 5.226280874236734e-7*pow(m,4)*pow(q,4) +
   	   1.9556388898926305e-7*pow(q,6) - 4.6638080530654496e-7*pow(m,2)*pow(q,6) +
   		1.2841163713979226e-7*pow(m,4)*pow(q,6) - 7.487237205340593e-9*pow(m,6)*pow(q,6) -
   		 1.7560471366298662e-9*pow(q,8) + 5.295790435014114e-9*pow(m,2)*pow(q,8) -
   		  1.5939817809091327e-9*pow(m,4)*pow(q,8) + 6.697715913198669e-11*pow(m,6)*pow(q,8) +
   		   3.346878886854118e-12*pow(m,8)*pow(q,8) - 4.698107572087482e-12*pow(q,10) -
   			9.888642283839855e-11*pow(m,2)*pow(q,10) + 1.5456555624443152e-10*pow(m,4)*pow(q,10) -
   			 5.305990934248017e-11*pow(m,6)*pow(q,10) + 5.949256896831782e-12*pow(m,8)*pow(q,10) -
   			  2.0177534761448037e-13*pow(m,10)*pow(q,10) + 9.424488124635083e-14*pow(m,2)*pow(q,12) -
   			   4.9221379639091777e-14*pow(m,4)*pow(q,12) - 7.62628913263295e-14*pow(m,6)*pow(q,12) +
   				3.535259953216403e-14*pow(m,8)*pow(q,12) - 4.260686688907054e-15*pow(m,10)*pow(q,12) +
   				 1.4747687581347302e-16*pow(m,12)*pow(q,12));

	return 4. + (nu2*pow(eps, 2) + nu3*pow(eps, 3) + nu4*pow(eps, 4) + nu5*pow(eps, 5) + nu6*pow(eps, 6)
		+ nu7*pow(eps, 7) + nu8*pow(eps, 8) + nu9*pow(eps, 9) + nu10*pow(eps, 10) + nu11*pow(eps, 11) + nu12*pow(eps, 12));
}

double nu_solver_low_freq_s0_l5(double q, int m, double eps){
	double nu2 = -(439./2574);
	double nu3 = (49*m*q)/12870.;
	double nu4 = (-8180730394 + 950331096*pow(q,2) - 124185204*pow(m,2)*pow(q,2))/1.19377826568e12;
	double nu5 = (0.0005569197021871516*m*q - 0.000044119381726219335*m*pow(q,3) +
   	 3.0398534672038946e-6*pow(m,3)*pow(q,3));
	double nu6 = (-0.0006396885737871447 + 0.00016303996921880097*pow(q,2) -
   	 0.000033116477561996*pow(m,2)*pow(q,2) - 3.2887497031771704e-6*pow(q,4) +
   	  1.9017723984030875e-6*pow(m,2)*pow(q,4) - 9.239962907837447e-8*pow(m,4)*pow(q,4));
	double nu7 = (0.00009396100679264214*m*q - 0.00001838359835766925*m*pow(q,3) +
   	 1.715245279645849e-6*pow(m,3)*pow(q,3) + 3.3609681107602414e-7*m*pow(q,5) -
   	  7.851977302416622e-8*pow(m,3)*pow(q,5) + 2.9622489343334737e-9*pow(m,5)*pow(q,5));
	double nu8 = (-0.00008004631964909296 + 0.000031738461029965736*pow(q,2) -
   	 8.814965375412861e-6*pow(m,2)*pow(q,2) - 2.0641641893717133e-6*pow(q,4) +
   	  1.374896479156856e-6*pow(m,2)*pow(q,4) - 8.276100317584947e-8*pow(m,4)*pow(q,4) +
   	   1.621886037063658e-8*pow(q,6) - 2.230862402874109e-8*pow(m,2)*pow(q,6) +
   		3.2345998709608336e-9*pow(m,4)*pow(q,6) - 1.0074570169186828e-10*pow(m,6)*pow(q,6));
	double nu9 = (0.000017294019170899545*m*q - 5.6659137480121945e-6*m*pow(q,3) +
   	 6.68947992208032e-7*pow(m,3)*pow(q,3) + 3.1822277752817874e-7*m*pow(q,5) -
   	  8.711834511226844e-8*pow(m,3)*pow(q,5) + 3.877196773586953e-9*pow(m,5)*pow(q,5) -
   	   2.6056186113933416e-9*m*pow(q,7) + 1.3304160519307854e-9*pow(m,3)*pow(q,7) -
   		1.3938878958424947e-10*pow(m,5)*pow(q,7) + 3.718373616885884e-12*pow(m,7)*pow(q,7));
	double nu10 = (-0.00001173469568476099 + 6.360479700774421e-6*pow(q,2) -
   	 2.2444410586153035e-6*pow(m,2)*pow(q,2) - 7.837458921130684e-7*pow(q,4) +
   	  6.227918703554257e-7*pow(m,2)*pow(q,4) - 4.5029502292248026e-8*pow(m,4)*pow(q,4) +
   	   2.1347915595317195e-8*pow(q,6) - 3.112882466636513e-8*pow(m,2)*pow(q,6) +
   		5.1365230054998835e-9*pow(m,4)*pow(q,6) - 1.821227951488509e-10*pow(m,6)*pow(q,6) -
   		 1.0291128576566343e-10*pow(q,8) + 2.6398789580776044e-10*pow(m,2)*pow(q,8) -
   		  7.894388740918063e-11*pow(m,4)*pow(q,8) + 6.541296855772758e-12*pow(m,6)*pow(q,8) -
   		   1.5409675392959295e-13*pow(m,8)*pow(q,8));
	double nu11 = (3.3818203761426274e-6*m*q - 1.5863923402221822e-6*m*pow(q,3) +
   	 2.269111877907427e-7*pow(m,3)*pow(q,3) + 1.728713166378909e-7*m*pow(q,5) -
   	  5.497823294461041e-8*pow(m,3)*pow(q,5) + 2.8369223019368505e-9*pow(m,5)*pow(q,5) -
   	   4.56217600876917e-9*m*pow(q,7) + 2.5584200579153822e-9*pow(m,3)*pow(q,7) -
   		2.9829638040254077e-10*pow(m,5)*pow(q,7) + 8.818942670723834e-12*pow(m,7)*pow(q,7) +
   		 2.2949135315161178e-11*m*pow(q,9) - 2.0025069237349458e-11*pow(m,3)*pow(q,9) +
   		  3.86380492980225e-12*pow(m,5)*pow(q,9) - 2.399906580518349e-13*pow(m,7)*pow(q,9) +
   		   4.418453750381894e-15*pow(m,9)*pow(q,9) + 8.544643611109884e-13*m*pow(q,11) -
   			1.2506035329704995e-12*pow(m,3)*pow(q,11) + 4.5363750282593796e-13*pow(m,5)*pow(q,11) -
   			 6.070257232059315e-14*pow(m,7)*pow(q,11) + 3.263579157021137e-15*pow(m,9)*pow(q,11) -
   			  5.933780285492975e-17*pow(m,11)*pow(q,11));
	double nu12 = (-1.908298770727694e-6 + 1.3216999950017323e-6*pow(q,2) -
   	 5.660337703238208e-7*pow(m,2)*pow(q,2) - 2.5039192341829626e-7*pow(q,4) +
   	  2.3427899992525377e-7*pow(m,2)*pow(q,4) - 1.9832696004754464e-8*pow(m,4)*pow(q,4) +
   	   1.4440177099948324e-8*pow(q,6) - 2.3258912910135888e-8*pow(m,2)*pow(q,6) +
   		4.349170982388321e-9*pow(m,4)*pow(q,6) - 1.741992203074586e-10*pow(m,6)*pow(q,6) -
   		 2.3655558577904256e-10*pow(q,8) + 6.13078545600735e-10*pow(m,2)*pow(q,8) -
   		  1.975297138162244e-10*pow(m,4)*pow(q,8) + 1.767553889306433e-11*pow(m,6)*pow(q,8) -
   		   4.467502405093565e-13*pow(m,8)*pow(q,8) + 5.01821319205909e-13*pow(q,10) -
   			1.773209652882549e-12*pow(m,2)*pow(q,10) + 6.814296134975218e-13*pow(m,4)*pow(q,10) -
   			 6.230429429957448e-14*pow(m,6)*pow(q,10) + 6.464896620162339e-16*pow(m,8)*pow(q,10) +
   			  4.980464546393227e-17*pow(m,10)*pow(q,10) + 1.480515995764092e-14*pow(m,2)*pow(q,12) -
   			   2.1668996615780555e-14*pow(m,4)*pow(q,12) + 7.86010054695589e-15*pow(m,6)*pow(q,12) -
   				1.051783238657407e-15*pow(m,8)*pow(q,12) + 5.654748594932295e-17*pow(m,10)*pow(q,12) -
   				 1.0281361081695082e-18*pow(m,12)*pow(q,12));

	return 5. + (nu2*pow(eps, 2) + nu3*pow(eps, 3) + nu4*pow(eps, 4) + nu5*pow(eps, 5) + nu6*pow(eps, 6)
		+ nu7*pow(eps, 7) + nu8*pow(eps, 8) + nu9*pow(eps, 9) + nu10*pow(eps, 10) + nu11*pow(eps, 11) + nu12*pow(eps, 12));
}

double nu_solver_low_freq_s0_l6(double q, int m, double eps){
	double nu2 = -(619./4290);
	double nu3 = (23*m*q)/10010.;
	double nu4 = (-77770592698 + 8952423480*pow(q,2) - 835769220*pow(m,2)*pow(q,2))/1.8790954182e13;
	double nu5 = (0.00023914201485180487*m*q - 0.000018540827364356776*m*pow(q,3) +
   	 9.113979702214997e-7*pow(m,3)*pow(q,3));
	double nu6 = (-0.00027484134500451227 + 0.00006919181879155261*pow(q,2) -
   	 0.000010046033065130574*pow(m,2)*pow(q,2) - 1.3441582673715962e-6*pow(q,4) +
   	  5.547164754767922e-7*pow(m,2)*pow(q,4) - 1.9183074351752726e-8*pow(m,4)*pow(q,4));
	double nu7 = (0.00002859721986200323*m*q - 5.467095970528956e-6*m*pow(q,3) +
   	 3.643320075523762e-7*pow(m,3)*pow(q,3) + 9.427858433789445e-8*m*pow(q,5) -
   	  1.5636942515842828e-8*pow(m,3)*pow(q,5) + 4.186213712256188e-10*pow(m,5)*pow(q,5));
	double nu8 = (-0.000024374238141059716 + 9.507659621449959e-6*pow(q,2) -
   	 1.8879997932486744e-6*pow(m,2)*pow(q,2) - 5.966747087806834e-7*pow(q,4) +
   	  2.83365928474378e-7*pow(m,2)*pow(q,4) - 1.2165582358257059e-8*pow(m,4)*pow(q,4) +
   	   4.280946017663346e-9*pow(q,6) - 4.2010491109427126e-9*pow(m,2)*pow(q,6) +
   		4.3032518179041774e-10*pow(m,4)*pow(q,6) - 9.4728190663287e-12*pow(m,6)*pow(q,6));
	double nu9 = (3.714779906355213e-6*m*q - 1.1845986335499942e-6*m*pow(q,3) +
   	 9.9931685528607e-8*pow(m,3)*pow(q,3) + 6.292986746532372e-8*m*pow(q,5) -
   	  1.2254909534456153e-8*pow(m,3)*pow(q,5) + 3.8821192753245914e-10*pow(m,5)*pow(q,5) -
   	   4.519582010491235e-10*m*pow(q,7) + 1.6287424638117368e-10*pow(m,3)*pow(q,7) -
   		1.1987164851078195e-11*pow(m,5)*pow(q,7) + 2.2490366327020875e-13*pow(m,7)*pow(q,7));
	double nu10 = (-2.5216426400477936e-6 + 1.33824946372184e-6*pow(q,2) -
   	 3.3765539949693756e-7*pow(m,2)*pow(q,2) - 1.5860003824560994e-7*pow(q,4) +
   	  8.983527574762913e-8*pow(m,2)*pow(q,4) - 4.635545368452535e-9*pow(m,4)*pow(q,4) +
   	   3.971716836354415e-9*pow(q,6) - 4.124002861887549e-9*pow(m,2)*pow(q,6) +
   		4.824988327324821e-10*pow(m,4)*pow(q,6) - 1.2141437955777706e-11*pow(m,6)*pow(q,6) -
   		 1.5430354684826906e-11*pow(q,8) + 2.8497171572254127e-11*pow(m,2)*pow(q,8) -
   		  5.95561400193269e-12*pow(m,4)*pow(q,8) + 3.437983368038397e-13*pow(m,6)*pow(q,8) -
   		   5.653050979421336e-15*pow(m,8)*pow(q,8));
	double nu11 = (5.101278704227192e-7*m*q - 2.3165158722224514e-7*m*pow(q,3) +
   	 2.3676014541594396e-8*pow(m,3)*pow(q,3) + 2.3789918291835723e-8*m*pow(q,5) -
   	  5.3829121996630904e-9*pow(m,3)*pow(q,5) + 1.9788149245229367e-10*pow(m,5)*pow(q,5) -
   	   5.544580627858854e-10*m*pow(q,7) + 2.202609152954293e-10*pow(m,3)*pow(q,7) -
   		1.815163581273098e-11*pow(m,5)*pow(q,7) + 3.8003160677219867e-13*pow(m,7)*pow(q,7) +
   		 2.359081686260835e-12*m*pow(q,9) - 1.5207235243953144e-12*pow(m,3)*pow(q,9) +
   		  2.191160602111118e-13*pow(m,5)*pow(q,9) - 1.044835338006979e-14*pow(m,7)*pow(q,9) +
   		   1.5361771035378103e-16*pow(m,9)*pow(q,9));
	double nu12 = (-2.8793581378728597e-7 + 1.9414190304390657e-7*pow(q,2) -
   	 5.944185472521236e-8*pow(m,2)*pow(q,2) - 3.514939240131911e-8*pow(q,4) +
   	  2.3433972846654205e-8*pow(m,2)*pow(q,4) - 1.415802581860885e-9*pow(m,4)*pow(q,4) +
   	   1.857797808056773e-9*pow(q,6) - 2.1267210485151262e-9*pow(m,2)*pow(q,6) +
   		2.821889952798674e-10*pow(m,4)*pow(q,6) - 8.033706919020842e-12*pow(m,6)*pow(q,6) -
   		 2.521386501117784e-11*pow(q,8) + 4.6748223167983136e-11*pow(m,2)*pow(q,8) -
   		  1.0659438000973468e-11*pow(m,4)*pow(q,8) + 6.757233408327805e-13*pow(m,6)*pow(q,8) -
   		   1.2158141235766626e-14*pow(m,8)*pow(q,8) + 7.306626800373558e-14*pow(q,10) -
   			2.1448094724646524e-13*pow(m,2)*pow(q,10) + 7.741188805034364e-14*pow(m,4)*pow(q,10) -
   			 8.519892038378952e-15*pow(m,6)*pow(q,10) + 3.4907886859102076e-16*pow(m,8)*pow(q,10) -
   			  4.6759913086735494e-18*pow(m,10)*pow(q,10));

	return 6. + (nu2*pow(eps, 2) + nu3*pow(eps, 3) + nu4*pow(eps, 4) + nu5*pow(eps, 5) + nu6*pow(eps, 6)
		+ nu7*pow(eps, 7) + nu8*pow(eps, 8) + nu9*pow(eps, 9) + nu10*pow(eps, 10) + nu11*pow(eps, 11) + nu12*pow(eps, 12));
}

double nu_solver_low_freq_s0_l(double q, int l, int m, double eps){
	double nu2 = (11 - 15*l - 15*pow(l,2))/(2.*(-1 + 2*l)*(1 + 2*l)*(3 + 2*l));
	double nu3 = ((-3 + 5*l + 5*pow(l,2))*m*q)/(l*(1 + l)*(-1 + 2*l)*(1 + 2*l)*(3 + 2*l));
	double nu4 = (-3240 - 8733*l + 73892*pow(l,2) + 9955*pow(l,3) - 278260*pow(l,4) -
     64365*pow(l,5) + 382305*pow(l,6) + 235200*pow(l,7) - 79800*pow(l,8) -
     92400*pow(l,9) - 18480*pow(l,10) +
     (-756*l - 1116*pow(l,2) + 7338*pow(l,3) + 13942*pow(l,4) - 13234*pow(l,5) -
        40774*pow(l,6) - 18688*pow(l,7) + 10928*pow(l,8) + 10400*pow(l,9) +
        2080*pow(l,10))*pow(q,2) +
     (3240 + 2106*l - 31068*pow(l,2) - 27420*pow(l,3) + 75450*pow(l,4) +
        84144*pow(l,5) - 10032*pow(l,6) - 32640*pow(l,7) - 8160*pow(l,8))*
      pow(m,2)*pow(q,2))/
   (8.*l*(1 + l)*(-3 + 2*l)*pow(-1 + 2*l,3)*pow(1 + 2*l,3)*pow(3 + 2*l,3)*
     (5 + 2*l));
	double nu5 = ((-1035*l - 2385*pow(l,2) + 6265*pow(l,3) + 17185*pow(l,4) - 4865*pow(l,5) -
        32795*pow(l,6) - 16640*pow(l,7) + 8440*pow(l,8) + 8400*pow(l,9) +
        1680*pow(l,10))*pow(m,3)*pow(q,3) +
     m*((-2970 - 288*l + 61605*pow(l,2) - 51024*pow(l,3) - 275448*pow(l,4) +
           143790*pow(l,5) + 558132*pow(l,6) + 13524*pow(l,7) -
           447951*pow(l,8) - 178920*pow(l,9) + 86184*pow(l,10) +
           66528*pow(l,11) + 11088*pow(l,12))*q +
        (-45*l + 387*pow(l,2) + 1874*pow(l,3) - 1201*pow(l,4) - 11470*pow(l,5) -
           7024*pow(l,6) + 17972*pow(l,7) + 24617*pow(l,8) + 4440*pow(l,9) -
           8088*pow(l,10) - 4896*pow(l,11) - 816*pow(l,12))*pow(q,3)))/
   (2.*(-1 + l)*pow(l,2)*pow(1 + l,2)*(2 + l)*(-3 + 2*l)*pow(-1 + 2*l,3)*
     pow(1 + 2*l,3)*pow(3 + 2*l,3)*(5 + 2*l));
	double nu6 = (112266000*l + 260690400*pow(l,2) - 2287534590*pow(l,3) - 8604512637*pow(l,4) +
     29309477879*pow(l,5) + 71867490821*pow(l,6) - 160551523104*pow(l,7) -
     337532467695*pow(l,8) + 413178904251*pow(l,9) + 989118681039*pow(l,10) -
     379995732626*pow(l,11) - 1674785645618*pow(l,12) - 354515937314*pow(l,13) +
     1372123946578*pow(l,14) + 890874540864*pow(l,15) - 369348082488*pow(l,16) -
     536988010944*pow(l,17) - 84800033472*pow(l,18) + 106275513344*pow(l,19) +
     51862883072*pow(l,20) + 1200815616*pow(l,21) - 4700247552*pow(l,22) -
     1254629376*pow(l,23) - 104552448*pow(l,24) +
     (32319000*pow(l,2) + 113785560*pow(l,3) - 699609744*pow(l,4) -
        2558184048*pow(l,5) + 4929087906*pow(l,6) + 22980548136*pow(l,7) -
        7380477006*pow(l,8) - 97899968196*pow(l,9) - 55043954118*pow(l,10) +
        184813002000*pow(l,11) + 236007376830*pow(l,12) - 83196181068*pow(l,13) -
        298607171604*pow(l,14) - 114255368064*pow(l,15) + 117987798096*pow(l,16) +
        111224750976*pow(l,17) + 7772312448*pow(l,18) - 26938899456*pow(l,19) -
        11561238528*pow(l,20) + 51803136*pow(l,21) + 1179859968*pow(l,22) +
        306561024*pow(l,23) + 25546752*pow(l,24))*pow(q,2) +
     (196020*pow(l,3) + 243702*pow(l,4) - 7026984*pow(l,5) - 24991539*pow(l,6) +
        16309614*pow(l,7) + 191649249*pow(l,8) + 176532192*pow(l,9) -
        498917181*pow(l,10) - 1014435954*pow(l,11) + 180622017*pow(l,12) +
        2040168216*pow(l,13) + 1600408296*pow(l,14) - 832381056*pow(l,15) -
        1928205840*pow(l,16) - 828686592*pow(l,17) + 358451712*pow(l,18) +
        473230848*pow(l,19) + 143785728*pow(l,20) - 18450432*pow(l,21) -
        22591488*pow(l,22) - 5455872*pow(l,23) - 454656*pow(l,24))*pow(q,4) +
     (-17010000*pow(l,2) - 87080670*pow(l,3) + 58122873*pow(l,4) +
        988041258*pow(l,5) + 962718684*pow(l,6) - 3795943844*pow(l,7) -
        7514516189*pow(l,8) + 3199152504*pow(l,9) + 18347160072*pow(l,10) +
        10742965728*pow(l,11) - 11670601776*pow(l,12) - 16269354752*pow(l,13) -
        3022658816*pow(l,14) + 5076289536*pow(l,15) + 3197951232*pow(l,16) +
        259184640*pow(l,17) - 329840640*pow(l,18) - 113254400*pow(l,19) -
        11325440*pow(l,20))*pow(m,4)*pow(q,4) +
     pow(m,2)*((30618000 - 157609800*l - 1355395680*pow(l,2) +
           3916724868*pow(l,3) + 19202952798*pow(l,4) - 31841345790*pow(l,5) -
           134993476350*pow(l,6) + 95986837644*pow(l,7) + 515885854074*pow(l,8) +
           2677994550*pow(l,9) - 1029256278030*pow(l,10) -
           550428150312*pow(l,11) + 896887674828*pow(l,12) +
           929499411120*pow(l,13) - 157316838960*pow(l,14) -
           535015756800*pow(l,15) - 153173865600*pow(l,16) +
           97808578560*pow(l,17) + 65812239360*pow(l,18) + 5329766400*pow(l,19) -
           5477606400*pow(l,20) - 1717309440*pow(l,21) - 156119040*pow(l,22))*
         pow(q,2) + (-3402000*pow(l,2) - 4374000*pow(l,3) + 93589020*pow(l,4) +
           238500936*pow(l,5) - 560823552*pow(l,6) - 2284756288*pow(l,7) +
           291783596*pow(l,8) + 8733820008*pow(l,9) + 7282297224*pow(l,10) -
           12566212224*pow(l,11) - 22984392768*pow(l,12) - 3121970944*pow(l,13) +
           18739513088*pow(l,14) + 14592540672*pow(l,15) - 688899072*pow(l,16) -
           5695518720*pow(l,17) - 2482452480*pow(l,18) + 21463040*pow(l,19) +
           304922624*pow(l,20) + 86507520*pow(l,21) + 7864320*pow(l,22))*pow(q,4)
        ))/(16.*(-1 + l)*pow(l,3)*pow(1 + l,3)*(2 + l)*(-5 + 2*l)*pow(-3 + 2*l,2)*
     pow(-1 + 2*l,5)*pow(1 + 2*l,5)*pow(3 + 2*l,5)*pow(5 + 2*l,2)*(7 + 2*l));

	return l + nu2*pow(eps, 2) + nu3*pow(eps, 3) + nu4*pow(eps, 4) + nu5*pow(eps, 5) + nu6*pow(eps, 6);
}

double nu_solver_low_freq_s2_l2(double q, int m, double eps){
	double nu2 = -0.5095238095238095;
	double nu3 = (m*q)/12.;
	double nu4 = (-10171398 + 2425500*pow(q,2) - 1127000*pow(m,2)*pow(q,2))/5.5566e7;
	double nu5 = (0.11154149389073732*m*q - 0.006000284329723582*m*pow(q,3) -
   	 0.004320289421102851*pow(m,3)*pow(q,3) - 0.008722741433021807*m*pow(q,5) +
   	  0.010903426791277258*pow(m,3)*pow(q,5) - 0.0021806853582554517*pow(m,5)*pow(q,5));
	double nu6 = (-0.1596012468191163 + 0.07184075359422912*pow(q,2) - 0.04804247212649977*pow(m,2)*pow(q,2) +
   	 0.000040996328374207876*pow(q,4) - 0.0025083711304996403*pow(m,2)*pow(q,4) -
   	  0.0008433593396013377*pow(m,4)*pow(q,4) - 0.006272584052302805*pow(m,2)*pow(q,6) +
   	   0.007840730065378505*pow(m,4)*pow(q,6) - 0.0015681460130757012*pow(m,6)*pow(q,6));
	double nu7 = (0.1611439159638726*m*q - 0.008224523769804466*m*pow(q,3) -
   	 0.01493226808628193*pow(m,3)*pow(q,3) - 0.02756414056326728*m*pow(q,5) +
   	  0.033160093306407*pow(m,3)*pow(q,5) - 0.0068100241054788214*pow(m,5)*pow(q,5) -
   	   0.0014642807463141117*m*pow(q,7) - 0.00015321499641052325*pow(m,3)*pow(q,7) +
   		0.0021133872250504258*pow(m,5)*pow(q,7) - 0.0004958914823257908*pow(m,7)*pow(q,7));
	double nu8 = (-0.18404997180331104 + 0.10538411700149444*pow(q,2) - 0.08365938254027563*pow(m,2)*pow(q,2) +
   	 0.0020633625611387323*pow(q,4) - 0.008975969208948446*pow(m,2)*pow(q,4) -
   	  0.0026628715936777795*pow(m,4)*pow(q,4) + 0.000034253467382628035*pow(q,6) -
   	   0.020258975041701266*pow(m,2)*pow(q,6) + 0.024802147825988694*pow(m,4)*pow(q,6) -
   		0.004929157672353353*pow(m,6)*pow(q,6) - 0.0010731579291939123*pow(m,2)*pow(q,8) +
   		 0.0011157961613449594*pow(m,4)*pow(q,8) + 0.000013774580385810716*pow(m,6)*pow(q,8) -
   		  0.00005641281253685776*pow(m,8)*pow(q,8) + 0.00007466404580628326*pow(m,2)*pow(q,10) -
   		   0.00018666011451570815*pow(m,4)*pow(q,10) + 0.00015399459447545925*pow(m,6)*pow(q,10) -
   			0.000046665028628927037*pow(m,8)*pow(q,10) + 4.666502862892704e-6*pow(m,10)*pow(q,10));
	double nu9 = (0.2425453944472263*m*q - 0.02524853584073296*m*pow(q,3) -
   	 0.01563834844922272*pow(m,3)*pow(q,3) - 0.040574631925298615*m*pow(q,5) +
   	  0.04612310130650436*pow(m,3)*pow(q,5) - 0.00966260323955135*pow(m,5)*pow(q,5) -
   	   0.005137620410099943*m*pow(q,7) + 0.00020811980534732373*pow(m,3)*pow(q,7) +
   		0.0063635274067174175*pow(m,5)*pow(q,7) - 0.0015031234490539666*pow(m,7)*pow(q,7) -
   		 0.00015691168890066365*m*pow(q,9) - 0.0000798772502072253*pow(m,3)*pow(q,9) +
   		  0.0003877051450583064*pow(m,5)*pow(q,9) - 0.00017139420360470588*pow(m,7)*pow(q,9) +
   		   0.000020477997654288433*pow(m,9)*pow(q,9) + 0.00011959427897839245*pow(m,3)*pow(q,11) -
   			0.0002989856974459811*pow(m,5)*pow(q,11) + 0.00024666320039293444*pow(m,7)*pow(q,11) -
   			 0.00007474642436149528*pow(m,9)*pow(q,11) + 7.474642436149528e-6*pow(m,11)*pow(q,11));
	double nu10 = (-0.24018689031038493 + 0.15850435014906356*pow(q,2) - 0.15339324400328044*pow(m,2)*pow(q,2) +
   	 0.0001794725228493287*pow(q,4) - 0.006949235555770597*pow(m,2)*pow(q,4) -
   	  0.004735109016885951*pow(m,4)*pow(q,4) + 0.00022929702770840485*pow(q,6) -
   	   0.029358907856902185*pow(m,2)*pow(q,6) + 0.03455413626613727*pow(m,4)*pow(q,6) -
   		0.006643325301892074*pow(m,6)*pow(q,6) + 1.928011168318525e-6*pow(q,8) -
   		 0.004069865197902102*pow(m,2)*pow(q,8) + 0.004903988300065199*pow(m,4)*pow(q,8) -
   		  0.0007944600680882113*pow(m,6)*pow(q,8) - 0.000036346027015133126*pow(m,8)*pow(q,8) +
   		   0.0003275788565305283*pow(m,2)*pow(q,10) - 0.0009164281358539712*pow(m,4)*pow(q,10) +
   			0.0008724316408274613*pow(m,6)*pow(q,10) - 0.00032279329155559747*pow(m,8)*pow(q,10) +
   			 0.000039210930051578956*pow(m,10)*pow(q,10) +
   			  0.00003146405897509786*pow(m,2)*pow(q,12) + 0.000010494916552660167*pow(m,4)*pow(q,12) -
   			   0.00015799303833987272*pow(m,6)*pow(q,12) + 0.0001642172826207738*pow(m,8)*pow(q,12) -
   				0.0000537554113080594*pow(m,10)*pow(q,12) + 5.572191499400301e-6*pow(m,12)*pow(q,12));
	double nu11 = (0.3796285607732202*m*q - 0.06989090954578227*m*pow(q,3) -
   	 0.00114330809777122*pow(m,3)*pow(q,3) - 0.05888326160337101*m*pow(q,5) +
   	  0.06814175569945309*pow(m,3)*pow(q,5) - 0.014174388027653941*pow(m,5)*pow(q,5) -
   	   0.007702526260697707*m*pow(q,7) + 0.0007898466420559558*pow(m,3)*pow(q,7) +
   		0.008646961984779952*pow(m,5)*pow(q,7) - 0.0019444383647775863*pow(m,7)*pow(q,7) -
   		 0.0007125732162696677*m*pow(q,9) - 0.000047994152427082276*pow(m,3)*pow(q,9) +
   		  0.0016880902585819394*pow(m,5)*pow(q,9) - 0.001086450167525426*pow(m,7)*pow(q,9) +
   		   0.00017118131917642752*pow(m,9)*pow(q,9) - 0.000015018199721856467*m*pow(q,11) +
   			0.0007299348901492759*pow(m,3)*pow(q,11) - 0.0017353265354060024*pow(m,5)*pow(q,11) +
   			 0.0014200222534963185*pow(m,7)*pow(q,11) - 0.00044685180378127294*pow(m,9)*pow(q,11) +
   			  0.00004723939526353739*pow(m,11)*pow(q,11) + 0.00005091045268104386*pow(m,3)*pow(q,13) -
   			   0.00008980630983867067*pow(m,5)*pow(q,13) + 0.000011328253994805542*pow(m,7)*pow(q,13) +
   				0.000045462474668721706*pow(m,9)*pow(q,13) -
   				 0.00002023673537239661*pow(m,11)*pow(q,13) + 2.3418638664961855e-6*pow(m,13)*pow(q,13) -
   				  1.2782035966487483e-6*pow(m,3)*pow(q,15) + 4.7932634874328064e-6*pow(m,5)*pow(q,15) -
   				   6.950232056777569e-6*pow(m,7)*pow(q,15) + 4.893123143420989e-6*pow(m,9)*pow(q,15) -
   					1.7375580141943922e-6*pow(m,11)*pow(q,15) + 2.995789679645504e-7*pow(m,13)*pow(q,15) -
   					 1.9971931197636692e-8*pow(m,15)*pow(q,15));
	double nu12 = (-0.33445503687254396 + 0.24569213330655992*pow(q,2) - 0.28721683994114044*pow(m,2)*pow(q,2) -
   	 0.008928600992991306*pow(q,4) + 0.009611478583084991*pow(m,2)*pow(q,4) -
   	  0.0076865938220997895*pow(m,4)*pow(q,4) + 0.0008585063488672768*pow(q,6) -
   	   0.03487494139741994*pow(m,2)*pow(q,6) + 0.038961471439744155*pow(m,4)*pow(q,6) -
   		0.007022450677183902*pow(m,6)*pow(q,6) + 0.000021827790331568266*pow(q,8) -
   		 0.007195781410130082*pow(m,2)*pow(q,8) + 0.009730180939334738*pow(m,4)*pow(q,8) -
   		  0.002722450851333896*pow(m,6)*pow(q,8) + 0.00023604883098223285*pow(m,8)*pow(q,8) +
   		   1.776711858188358e-7*pow(q,10) + 0.0004051482600611335*pow(m,2)*pow(q,10) -
   			0.001273310486599961*pow(m,4)*pow(q,10) + 0.0016154434634623345*pow(m,6)*pow(q,10) -
   			 0.0008765693276733532*pow(m,8)*pow(q,10) + 0.00013640139099122315*pow(m,10)*pow(q,10) +
   			  0.00019288908267118274*pow(m,2)*pow(q,12) + 0.00007382656456509507*pow(m,4)*pow(q,12) -
   			   0.0009534480865323845*pow(m,6)*pow(q,12) + 0.0009727627713826542*pow(m,8)*pow(q,12) -
   				0.00031953804261942894*pow(m,10)*pow(q,12) +
   				 0.000033507710532881375*pow(m,12)*pow(q,12) + 7.49180149826148e-6*pow(m,2)*pow(q,14) +
   				  0.000017843422502456556*pow(m,4)*pow(q,14) - 0.00006930215869952482*pow(m,6)*pow(q,14) +
   				   0.000054053493623847675*pow(m,8)*pow(q,14) - 8.615813879592867e-6*pow(m,10)*pow(q,14) -
   					1.8881398161096846e-6*pow(m,12)*pow(q,14) + 4.173947706616576e-7*pow(m,14)*pow(q,14) -
   					 3.1755992782784327e-6*pow(m,4)*pow(q,16) + 0.000011908497293544121*pow(m,6)*pow(q,16) -
   					  0.00001726732107563898*pow(m,8)*pow(q,16) +
   					   0.000012156590987159624*pow(m,10)*pow(q,16) - 4.316830268909745e-6*pow(m,12)*pow(q,16) +
   						7.442810808465076e-7*pow(m,14)*pow(q,16) - 4.961873872310051e-8*pow(m,16)*pow(q,16));

	return 2. + (nu2*pow(eps, 2) + nu3*pow(eps, 3) + nu4*pow(eps, 4) + nu5*pow(eps, 5) + nu6*pow(eps, 6)
		+ nu7*pow(eps, 7) + nu8*pow(eps, 8) + nu9*pow(eps, 9) + nu10*pow(eps, 10) + nu11*pow(eps, 11) + nu12*pow(eps, 12));
}

double nu_solver_low_freq_s2_l3(double q, int m, double eps){
	double nu2 = -(13./42);
	double nu3 = (7*m*q)/360.;
	double nu4 = (-3931560 + 658560*pow(q,2) - 164297*pow(m,2)*pow(q,2))/9.779616e7;
	double nu5 = (0.009574281171503394*m*q - 0.0010891320613542837*m*pow(q,3) +
   	 0.0001691490522509041*pow(m,3)*pow(q,3));
	double nu6 = (-0.01236891798047223 + 0.004384301409313434*pow(q,2) -
   	 0.0018144039509680042*pow(m,2)*pow(q,2) - 0.00018784059188099592*pow(q,4) +
   	  0.00017927465247891735*pow(m,2)*pow(q,4) - 0.000020319307969390976*pow(m,4)*pow(q,4));
	double nu7 = (0.0055301133932531285*m*q - 0.0015858594274482656*m*pow(q,3) +
   	 0.00031868280248856294*pow(m,3)*pow(q,3) + 0.00004518800555557009*m*pow(q,5) -
   	  0.000011563061102631187*pow(m,3)*pow(q,5) - 1.6661302213665483e-7*pow(m,5)*pow(q,5) +
   	   8.140008140008141e-6*m*pow(q,7) - 0.00001107945552389997*pow(m,3)*pow(q,7) +
   		3.165558721114277e-6*pow(m,5)*pow(q,7) - 2.2611133722244833e-7*pow(m,7)*pow(q,7));
	double nu8 = (-0.005264917347457815 + 0.0029297097923633862*pow(q,2) -
   	 0.0016992204777677006*pow(m,2)*pow(q,2) - 0.0003415288249901822*pow(q,4) +
   	  0.0004146475180269634*pow(m,2)*pow(q,4) - 0.00005309909319009391*pow(m,4)*pow(q,4) +
   	   3.0669075281869897e-6*pow(q,6) - 6.600092564946471e-7*pow(m,2)*pow(q,6) -
   		8.862165686137622e-7*pow(m,4)*pow(q,6) + 6.541294229087476e-8*pow(m,6)*pow(q,6) +
   		 2.0941388462755983e-6*pow(m,2)*pow(q,8) - 2.85035565187512e-6*pow(m,4)*pow(q,8) +
   		  8.143873291071772e-7*pow(m,6)*pow(q,8) - 5.8170523507655514e-8*pow(m,8)*pow(q,8));
	double nu9 = (0.0035977692421220865*m*q - 0.0017091054362959714*m*pow(q,3) +
   	 0.00042050333567595185*pow(m,3)*pow(q,3) + 0.00010583965223166894*m*pow(q,5) -
   	  0.000025379636401437505*pow(m,3)*pow(q,5) - 1.2710092243674545e-6*pow(m,5)*pow(q,5) +
   	   0.000026967815465898808*m*pow(q,7) - 0.0000368560313929517*pow(m,3)*pow(q,7) +
   		0.000010469993044216562*pow(m,5)*pow(q,7) - 7.472508292972927e-7*pow(m,7)*pow(q,7) +
   		 1.770941941882113e-7*m*pow(q,9) + 3.318458556240717e-8*pow(m,3)*pow(q,9) -
   		  3.043868019345085e-7*pow(m,5)*pow(q,9) + 1.0172550721125815e-7*pow(m,7)*pow(q,9) -
   		   7.617485027368065e-9*pow(m,9)*pow(q,9));
	double nu10 = (-0.0027281313156789395 + 0.002098061476132987*pow(q,2) -
   	 0.001545394240900196*pow(m,2)*pow(q,2) - 0.0004085285496756428*pow(q,4) +
   	  0.000587255938822716*pow(m,2)*pow(q,4) - 0.00008295349071169261*pow(m,4)*pow(q,4) +
   	   6.46280335081645e-6*pow(q,6) + 2.131940057088963e-6*pow(m,2)*pow(q,6) -
   		3.771120011361667e-6*pow(m,4)*pow(q,6) + 2.5712998241797006e-7*pow(m,6)*pow(q,6) +
   		 6.17533988361362e-8*pow(q,8) + 6.643815326529288e-6*pow(m,2)*pow(q,8) -
   		  9.276227726280598e-6*pow(m,4)*pow(q,8) + 2.6677214857547045e-6*pow(m,6)*pow(q,8) -
   		   1.9125887593064704e-7*pow(m,8)*pow(q,8) + 3.3318683329641027e-8*pow(m,2)*pow(q,10) -
   			1.802869819987563e-8*pow(m,4)*pow(q,10) - 2.423064710783057e-8*pow(m,6)*pow(q,10) +
   			 9.699598974945487e-9*pow(m,8)*pow(q,10) - 7.589369968803146e-10*pow(m,10)*pow(q,10));
	double nu11 = (0.0025552217337465725*m*q - 0.0016655360654594241*m*pow(q,3) +
   	 0.00047717576484634477*pow(m,3)*pow(q,3) + 0.00016021736937035788*m*pow(q,5) -
   	  0.00005066460120878628*pow(m,3)*pow(q,5) + 6.651218312692636e-8*pow(m,5)*pow(q,5) +
   	   0.00003899033786135922*m*pow(q,7) - 0.00005304709709613082*pow(m,3)*pow(q,7) +
   		0.0000149047165311553*pow(m,5)*pow(q,7) - 1.059982002247987e-6*pow(m,7)*pow(q,7) +
   		 6.049399305030931e-7*m*pow(q,9) + 5.068715107103255e-8*pow(m,3)*pow(q,9) -
   		  9.69008215070852e-7*pow(m,5)*pow(q,9) + 3.285584061488475e-7*pow(m,7)*pow(q,9) -
   		   2.471506685275762e-8*pow(m,9)*pow(q,9) - 1.8137706660010567e-9*m*pow(q,11) +
  			 9.59296290432776e-9*pow(m,3)*pow(q,11) - 8.481028421428965e-9*pow(m,5)*pow(q,11) +
   			  2.0597130796564514e-10*pow(m,7)*pow(q,11) + 5.492310183048747e-10*pow(m,9)*pow(q,11) -
   			   5.336614316825756e-11*pow(m,11)*pow(q,11));
	double nu12 = (-0.0016265713032847242 + 0.0015843711566805255*pow(q,2) -
   	 0.0013883996607896357*pow(m,2)*pow(q,2) - 0.0004161583109360726*pow(q,4) +
   	  0.0006848614392841574*pow(m,2)*pow(q,4) - 0.00010724648167146021*pow(m,4)*pow(q,4) +
   	   0.000011005796793159422*pow(q,6) - 4.427179775779799e-6*pow(m,2)*pow(q,6) -
   		2.9446796185692183e-6*pow(m,4)*pow(q,6) + 2.2185209369418957e-7*pow(m,6)*pow(q,6) +
   		 2.4167561411503264e-7*pow(q,8) + 8.549087112478261e-6*pow(m,2)*pow(q,8) -
   		  0.000012487789010214964*pow(m,4)*pow(q,8) + 3.631937289484359e-6*pow(m,6)*pow(q,8) -
   		   2.620365703745514e-7*pow(m,8)*pow(q,8) - 1.0440152502229138e-9*pow(q,10) +
   			1.0439488997013763e-7*pow(m,2)*pow(q,10) - 4.846764141730257e-8*pow(m,4)*pow(q,10) -
   			 8.350474495212279e-8*pow(m,6)*pow(q,10) + 3.227283713499168e-8*pow(m,8)*pow(q,10) -
   			  2.503163674918415e-9*pow(m,10)*pow(q,10) + 2.5614011525827815e-10*pow(m,2)*pow(q,12) +
   			   7.91847879411638e-10*pow(m,4)*pow(q,12) - 1.4111408947920634e-9*pow(m,6)*pow(q,12) +
   				3.7982044250549976e-10*pow(m,8)*pow(q,12) - 1.5512731463501408e-11*pow(m,10)*pow(q,12) -
   				 1.1548109198510125e-12*pow(m,12)*pow(q,12) + 1.070349525313365e-10*pow(m,2)*pow(q,14) -
   				  2.913729263353049e-10*pow(m,4)*pow(q,14) + 2.8154487128034425e-10*pow(m,6)*pow(q,14) -
   				   1.192580798265817e-10*pow(m,8)*pow(q,14) + 2.428107719460874e-11*pow(m,10)*pow(q,14) -
   					2.3124835423436896e-12*pow(m,12)*pow(q,14) + 8.258869794084606e-14*pow(m,14)*pow(q,14));

	return 3. + (nu2*pow(eps, 2) + nu3*pow(eps, 3) + nu4*pow(eps, 4) + nu5*pow(eps, 5) + nu6*pow(eps, 6)
		+ nu7*pow(eps, 7) + nu8*pow(eps, 8) + nu9*pow(eps, 9) + nu10*pow(eps, 10) + nu11*pow(eps, 11) + nu12*pow(eps, 12));
}

double nu_solver_low_freq_s2_l4(double q, int m, double eps){
	double nu2 = -(1571./6930);
	double nu3 = (53*m*q)/6600.;
	double nu4 = (-2739296734760 + 399322677600*pow(q,2) - 65253658239*pow(m,2)*pow(q,2))/1.7306252964e14;
	double nu5 = (0.002082966917812123*m*q - 0.0001951077541986633*m*pow(q,3) +
   	 0.000018835639949276312*pow(m,3)*pow(q,3));
	double nu6 = (-0.0025788774026921217 + 0.0007935572261318273*pow(q,2) -
   	 0.00020882036219149167*pow(m,2)*pow(q,2) - 0.00002550474281364407*pow(q,4) +
   	  0.000015002432823293192*pow(m,2)*pow(q,4) - 1.0115583670684378e-6*pow(m,4)*pow(q,4));
	double nu7 = (0.0006216699197513281*m*q - 0.00014703023388492917*m*pow(q,3) +
   	 0.000018667962084306105*pow(m,3)*pow(q,3) + 3.925249564187367e-6*m*pow(q,5) -
   	  1.1468752005979093e-6*pow(m,3)*pow(q,5) + 6.109849988390876e-8*pow(m,5)*pow(q,5));
	double nu8 = (-0.0005681515723797658 + 0.0002703351446715027*pow(q,2) -
   	 0.00009882784796586344*pow(m,2)*pow(q,2) - 0.00002463271670695482*pow(q,4) +
   	  0.000019329127506903887*pow(m,2)*pow(q,4) - 1.6110513371467068e-6*pow(m,4)*pow(q,4) +
   	   3.6803871323042377e-7*pow(q,6) - 4.933590882435913e-7*pow(m,2)*pow(q,6) +
   		9.617869815356124e-8*pow(m,4)*pow(q,6) - 4.346072662933072e-9*pow(m,6)*pow(q,6));
	double nu9 = (0.0002042170450555224*m*q - 0.00008157649745638222*m*pow(q,3) +
   	 0.000013087938269831425*pow(m,3)*pow(q,3) + 6.451828392645588e-6*m*pow(q,5) -
   	  2.259332601269539e-6*pow(m,3)*pow(q,5) + 1.4017788091545165e-7*pow(m,5)*pow(q,5) -
   	   7.904909733133957e-8*m*pow(q,7) + 3.824667097790258e-8*pow(m,3)*pow(q,7) -
   		3.602862286465049e-9*pow(m,5)*pow(q,7) + 4.5973545245460956e-11*pow(m,7)*pow(q,7) -
   		 6.300534150120834e-9*m*pow(q,9) + 8.96951042204702e-9*pow(m,3)*pow(q,9) -
   		  2.98619066490102e-9*pow(m,5)*pow(q,9) + 3.2815282031879344e-10*pow(m,7)*pow(q,9) -
   		   1.093842734395978e-11*pow(m,9)*pow(q,9));
	double nu10 = (-0.00014832751516800276 + 0.00009687884797851326*pow(q,2) -
   	 0.00004548989847978439*pow(m,2)*pow(q,2) - 0.000016109123139986*pow(q,4) +
   	  0.00001598605349029087*pow(m,2)*pow(q,4) - 1.5965774322538778e-6*pow(m,4)*pow(q,4) +
   	   7.242553409045931e-7*pow(q,6) - 1.171224766186511e-6*pow(m,2)*pow(q,6) +
   		2.5341496668764327e-7*pow(m,4)*pow(q,6) - 1.236954257551198e-8*pow(m,6)*pow(q,6) -
   		 3.853832136325672e-9*pow(q,8) + 4.5502360524470775e-9*pow(m,2)*pow(q,8) -
   		  3.3820888921067554e-12*pow(m,4)*pow(q,8) - 1.545977201725239e-10*pow(m,6)*pow(q,8) +
   		   8.580387264691758e-12*pow(m,8)*pow(q,8) - 7.902337676255947e-10*pow(m,2)*pow(q,10) +
   			1.1249855719669925e-9*pow(m,4)*pow(q,10) - 3.7453787944754747e-10*pow(m,6)*pow(q,10) +
   			 4.115800873049972e-11*pow(m,8)*pow(q,10) - 1.371933624349991e-12*pow(m,10)*pow(q,10));
	double nu11 = (0.00007229731995053895*m*q - 0.000041851850975868395*m*pow(q,3) +
   	 8.152682365850613e-6*pow(m,3)*pow(q,3) + 6.267275613006614e-6*m*pow(q,5) -
   	  2.578903017123856e-6*pow(m,3)*pow(q,5) + 1.8218618840577633e-7*pow(m,5)*pow(q,5) -
   	   1.9839434977029684e-7*m*pow(q,7) + 1.0174880617716764e-7*pow(m,3)*pow(q,7) -
   		9.536786829271686e-9*pow(m,5)*pow(q,7) + 8.982997701422847e-11*pow(m,7)*pow(q,7) -
   		 2.0641806335547422e-8*m*pow(q,9) + 2.989957739016302e-8*pow(m,3)*pow(q,9) -
   		  9.982326097378968e-9*pow(m,5)*pow(q,9) + 1.0953896923390268e-9*pow(m,7)*pow(q,9) -
   		   3.644068841470398e-11*pow(m,9)*pow(q,9) - 3.5054645975987895e-11*m*pow(q,11) -
   			2.674508234642293e-13*pow(m,3)*pow(q,11) + 5.481045451377056e-11*pow(m,5)*pow(q,11) -
   			 2.1953501376855717e-11*pow(m,7)*pow(q,11) + 2.5522471943618434e-12*pow(m,9)*pow(q,11) -
   			  8.710353182456075e-14*pow(m,11)*pow(q,11));
	double nu12 = (-0.000043623571802934326 + 0.00003680284225446887*pow(q,2) -
   	 0.000021155255980152927*pow(m,2)*pow(q,2) - 9.309964450656558e-6*pow(q,4) +
   	  0.00001118412056833486*pow(m,2)*pow(q,4) - 1.300008664207689e-6*pow(m,4)*pow(q,4) +
   	   7.992730146556814e-7*pow(q,6) - 1.4928178610511828e-6*pow(m,2)*pow(q,6) +
   		3.5548475137330175e-7*pow(m,4)*pow(q,6) - 1.8679663466185905e-8*pow(m,6)*pow(q,6) -
   		 9.425910186587566e-9*pow(q,8) + 1.0373327085768051e-8*pow(m,2)*pow(q,8) +
   		  1.106452960890174e-9*pow(m,4)*pow(q,8) - 5.960767360447302e-10*pow(m,6)*pow(q,8) +
   		   3.047382414772888e-11*pow(m,8)*pow(q,8) - 2.044387091905525e-11*pow(q,10) -
   			2.458401881916255e-9*pow(m,2)*pow(q,10) + 3.618318173491282e-9*pow(m,4)*pow(q,10) -
   			 1.2192399768132924e-9*pow(m,6)*pow(q,10) + 1.3475518918634421e-10*pow(m,8)*pow(q,10) -
   			  4.5072323890035306e-12*pow(m,10)*pow(q,10) + 4.13136457672929e-13*pow(m,2)*pow(q,12) -
   			   3.3648837479652292e-12*pow(m,4)*pow(q,12) + 4.148804673622602e-12*pow(m,6)*pow(q,12) -
   				1.3375756841180878e-12*pow(m,8)*pow(q,12) + 1.4533902664962118e-13*pow(m,10)*pow(q,12) -
   				 4.8207258618350046e-15*pow(m,12)*pow(q,12));

	return 4. + (nu2*pow(eps, 2) + nu3*pow(eps, 3) + nu4*pow(eps, 4) + nu5*pow(eps, 5) + nu6*pow(eps, 6)
		+ nu7*pow(eps, 7) + nu8*pow(eps, 8) + nu9*pow(eps, 9) + nu10*pow(eps, 10) + nu11*pow(eps, 11) + nu12*pow(eps, 12));
}

double nu_solver_low_freq_s2_l5(double q, int m, double eps){
	double nu2 = -(773./4290);
	double nu3 = (3743*m*q)/900900.;
	double nu4 = (-132069805110 + 17824134900*pow(q,2) - 2046900856*pow(m,2)*pow(q,2))/1.658025369e13;
	double nu5 = (0.0006752028412835749*m*q - 0.00005779703722316665*m*pow(q,3) +
   	 3.797099416961037e-6*pow(m,3)*pow(q,3));
	double nu6 = (-0.0008167184296781829 + 0.00023366834402763174*pow(q,2) -
   	 0.00004272052662039862*pow(m,2)*pow(q,2) - 6.469792835930856e-6*pow(q,4) +
   	  2.7098645446004687e-6*pow(m,2)*pow(q,4) - 1.235266790283047e-7*pow(m,4)*pow(q,4));
	double nu7 = (0.00012558908036495683*m*q - 0.00002706637496954619*m*pow(q,3) +
   	 2.3555320539280007e-6*pow(m,3)*pow(q,3) + 5.898842899406911e-7*m*pow(q,5) -
   	  1.194319156903998e-7*pow(m,3)*pow(q,5) + 4.279760435219111e-9*pow(m,5)*pow(q,5));
	double nu8 = (-0.00011229242091871227 + 0.000049443756761745385*pow(q,2) -
   	 0.000012489175218798599*pow(m,2)*pow(q,2) - 3.924673341984395e-6*pow(q,4) +
   	  2.1476295967662597e-6*pow(m,2)*pow(q,4) - 1.2150723370327796e-7*pow(m,4)*pow(q,4) +
   	   4.5437159948553737e-8*pow(q,6) - 4.0896727652535754e-8*pow(m,2)*pow(q,6) +
   		5.2990437286850295e-9*pow(m,4)*pow(q,6) - 1.5721752184479494e-10*pow(m,6)*pow(q,6));
	double nu9 = (0.00002546037742907608*m*q - 9.191279242382248e-6*m*pow(q,3) +
   	 1.0092139571666433e-6*pow(m,3)*pow(q,3) + 6.049285249243261e-7*m*pow(q,5) -
  	   1.4549298132929452e-7*pow(m,3)*pow(q,5) + 6.119490856725894e-9*pow(m,5)*pow(q,5) -
   		6.378660227599314e-9*m*pow(q,7) + 2.6085062084204085e-9*pow(m,3)*pow(q,7) -
  		  2.4916655635798147e-10*pow(m,5)*pow(q,7) + 6.3399783836769846e-12*pow(m,7)*pow(q,7));
	double nu10 = (-0.00001810219184882074 + 0.000010839940083123235*pow(q,2) -
   	 3.500576115607208e-6*pow(m,2)*pow(q,2) - 1.5659961764160183e-6*pow(q,4) +
   	  1.0708070993811506e-6*pow(m,2)*pow(q,4) - 7.267257576516977e-8*pow(m,4)*pow(q,4) +
   	   5.5674687348502194e-8*pow(q,6) - 6.208755374169536e-8*pow(m,2)*pow(q,6) +
   		9.248889566982133e-9*pow(m,4)*pow(q,6) - 3.1138797658202e-10*pow(m,6)*pow(q,6) -
  		  3.938695745905834e-10*pow(q,8) + 6.604816042283329e-10*pow(m,2)*pow(q,8) -
  			1.7023660950012862e-10*pow(m,4)*pow(q,8) + 1.3139985313629105e-11*pow(m,6)*pow(q,8) -
  			  2.9821202673243604e-13*pow(m,8)*pow(q,8));
	double nu11 = (5.4865209960662035e-6*m*q - 2.833366915784419e-6*m*pow(q,3) +
  	  3.7700606397980353e-7*pow(m,3)*pow(q,3) + 3.559696675253938e-7*m*pow(q,5) -
   	   1.0110609667859859e-7*pow(m,3)*pow(q,5) + 4.926593038286241e-9*pow(m,5)*pow(q,5) -
   		1.1629676555293055e-8*m*pow(q,7) + 5.457849681488183e-9*pow(m,3)*pow(q,7) -
  		 5.835870296014355e-10*pow(m,5)*pow(q,7) + 1.6412025667534726e-11*pow(m,7)*pow(q,7) +
		  7.402662793738869e-11*m*pow(q,9) - 4.644454073540921e-11*pow(m,3)*pow(q,9) +
  		   6.970959942516934e-12*pow(m,5)*pow(q,9) - 3.192867850226751e-13*pow(m,7)*pow(q,9) +
  			3.3283089963300053e-15*pow(m,9)*pow(q,9) + 3.566697840593494e-12*m*pow(q,11) -
   			 5.220258589468645e-12*pow(m,3)*pow(q,11) + 1.8935697910650877e-12*pow(m,5)*pow(q,11) -
   			  2.533841590921628e-13*pow(m,7)*pow(q,11) + 1.3622804252266817e-14*pow(m,9)*pow(q,11) -
  			   2.476873500412149e-16*pow(m,11)*pow(q,11));
	double nu12 = (-3.2407151335377104e-6 + 2.4727749395213963e-6*pow(q,2) -
   	 9.7360326128864e-7*pow(m,2)*pow(q,2) - 5.386704122295325e-7*pow(q,4) +
   	  4.4423247129101386e-7*pow(m,2)*pow(q,4) - 3.529729878477354e-8*pow(m,4)*pow(q,4) +
   	   3.824363562908837e-8*pow(q,6) - 5.0647882523720875e-8*pow(m,2)*pow(q,6) +
   		8.63134390388597e-9*pow(m,4)*pow(q,6) - 3.280066455794238e-10*pow(m,6)*pow(q,6) -
   		 8.346535905860814e-10*pow(q,8) + 1.6168178952807542e-9*pow(m,2)*pow(q,8) -
   		  4.546217062798678e-10*pow(m,4)*pow(q,8) + 3.758977625089462e-11*pow(m,6)*pow(q,8) -
   		   9.012500391316223e-13*pow(m,8)*pow(q,8) + 2.7032457766077796e-12*pow(q,10) -
   			4.666342640175444e-12*pow(m,2)*pow(q,10) + 8.164067490619133e-13*pow(m,4)*pow(q,10) +
   			 4.505545580036825e-14*pow(m,6)*pow(q,10) - 1.0488887363113028e-14*pow(m,8)*pow(q,10) +
   			  3.162317251044508e-16*pow(m,10)*pow(q,10) + 2.566127003018758e-13*pow(m,2)*pow(q,12) -
   			   3.75581199414051e-13*pow(m,4)*pow(q,12) + 1.3623639540332227e-13*pow(m,6)*pow(q,12) -
   				1.8230193917279094e-14*pow(m,8)*pow(q,12) + 9.801179525418868e-16*pow(m,10)*pow(q,12) -
   				 1.7820326409852486e-17*pow(m,12)*pow(q,12));

	return 5. + (nu2*pow(eps, 2) + nu3*pow(eps, 3) + nu4*pow(eps, 4) + nu5*pow(eps, 5) + nu6*pow(eps, 6)
		+ nu7*pow(eps, 7) + nu8*pow(eps, 8) + nu9*pow(eps, 9) + nu10*pow(eps, 10) + nu11*pow(eps, 11) + nu12*pow(eps, 12));
}

double nu_solver_low_freq_s2_l6(double q, int m, double eps){
	double nu2 = -(901./6006);
	double nu3 = (2053*m*q)/840840.;
	double nu4 = (-3560502411816 + 458208316566*pow(q,2) - 38859182648*pow(m,2)*pow(q,2))/7.7343567413112e14;
	double nu5 = (0.0002737793553863499*m*q - 0.000022339776246208947*m*pow(q,3) +
   	 1.0614973554348732e-6*pow(m,3)*pow(q,3));
	double nu6 = (-0.0003266209725655114 + 0.00008950830973746*pow(q,2) -
   	 0.000012012228269642315*pow(m,2)*pow(q,2) - 2.2483188451757706e-6*pow(q,4) +
   	  7.13386829982742e-7*pow(m,2)*pow(q,4) - 2.3495449459149105e-8*pow(m,4)*pow(q,4));
	double nu7 = (0.000035064729081712864*m*q - 7.179607351374349e-6*m*pow(q,3) +
   	 4.545233617590309e-7*pow(m,3)*pow(q,3) + 1.3989768202517944e-7*m*pow(q,5) -
   	  2.0864786019329748e-8*pow(m,3)*pow(q,5) + 5.379813157013021e-10*pow(m,5)*pow(q,5));
	double nu8 = (-0.00003094928739025956 + 0.000013032928856085444*pow(q,2) -
   	 2.409894279294523e-6*pow(m,2)*pow(q,2) - 9.5348760166523e-7*pow(q,4) +
   	  3.87418402249647e-7*pow(m,2)*pow(q,4) - 1.5875123139804374e-8*pow(m,4)*pow(q,4) +
   	   9.616657903311572e-9*pow(q,6) - 6.4422537276856215e-9*pow(m,2)*pow(q,6) +
   		6.0245363195442e-10*pow(m,4)*pow(q,6) - 1.2760294290650075e-11*pow(m,6)*pow(q,6));
	double nu9 = (4.873867150989921e-6*m*q - 1.6657864306072183e-6*m*pow(q,3) +
   	 1.3310064391031263e-7*pow(m,3)*pow(q,3) + 9.935858854452846e-8*m*pow(q,5) -
   	  1.749830338814879e-8*pow(m,3)*pow(q,5) + 5.310435805193492e-10*pow(m,5)*pow(q,5) -
   	   8.654271570641682e-10*m*pow(q,7) + 2.594618348437806e-10*pow(m,3)*pow(q,7) -
  		 1.7744472175943182e-11*pow(m,5)*pow(q,7) + 3.2090709878370047e-13*pow(m,7)*pow(q,7));
	double nu10 = (-3.4221052837127867e-6 + 1.952769373781023e-6*pow(q,2) -
   	 4.6074735102981277e-7*pow(m,2)*pow(q,2) - 2.6093479901672423e-7*pow(q,4) +
   	  1.3127867505767914e-7*pow(m,2)*pow(q,4) - 6.4569570634870805e-9*pow(m,4)*pow(q,4) +
   	   8.12004897653366e-9*pow(q,6) - 6.695971293736297e-9*pow(m,2)*pow(q,6) +
   		7.223399194687858e-10*pow(m,4)*pow(q,6) - 1.7457843576398517e-11*pow(m,6)*pow(q,6) -
   		 4.538736013494196e-11*pow(q,8) + 5.423337365652941e-11*pow(m,2)*pow(q,8) -
   		  9.938381208202007e-12*pow(m,4)*pow(q,8) + 5.395913042908788e-13*pow(m,6)*pow(q,8) -
   		   8.580566367862752e-15*pow(m,8)*pow(q,8));
	double nu11 = (7.161527464431974e-7*m*q - 3.4831817658559383e-7*m*pow(q,3) +
   	 3.3703686579929314e-8*pow(m,3)*pow(q,3) + 3.9681118959456256e-8*m*pow(q,5) -
   	  8.21465211415013e-9*pow(m,3)*pow(q,5) + 2.889537118776301e-10*pow(m,5)*pow(q,5) -
   	   1.0886298188642386e-9*m*pow(q,7) + 3.721642497239884e-10*pow(m,3)*pow(q,7) -
   		2.8620042817845484e-11*pow(m,5)*pow(q,7) + 5.767381770376647e-13*pow(m,7)*pow(q,7) +
   		 5.8459996744662154e-12*m*pow(q,9) - 2.9788295152783832e-12*pow(m,3)*pow(q,9) +
   		  3.873525033583375e-13*pow(m,5)*pow(q,9) - 1.7485730504684397e-14*pow(m,7)*pow(q,9) +
   		   2.488259824758018e-16*pow(m,9)*pow(q,9));
	double nu12 = (-4.1781652281424653e-7 + 3.02206154448098e-7*pow(q,2) -
   	 8.67813315675013e-8*pow(m,2)*pow(q,2) - 6.071468509643008e-8*pow(q,4) +
   	  3.6627780099644524e-8*pow(m,2)*pow(q,4) - 2.107946379321995e-9*pow(m,4)*pow(q,4) +
   	   3.7795978593450215e-9*pow(q,6) - 3.6653886874637096e-9*pow(m,2)*pow(q,6) +
   		4.515320724091072e-10*pow(m,4)*pow(q,6) - 1.2335435195470879e-11*pow(m,6)*pow(q,6) -
   		 6.613014452436287e-11*pow(q,8) + 9.294457353767022e-11*pow(m,2)*pow(q,8) -
   		  1.8909492528229306e-11*pow(m,4)*pow(q,8) + 1.1286426503240295e-12*pow(m,6)*pow(q,8) -
   		   1.959306101292229e-14*pow(m,8)*pow(q,8) + 2.648674522648088e-13*pow(q,10) -
   			5.19513960194141e-13*pow(m,2)*pow(q,10) + 1.6050298143471826e-13*pow(m,4)*pow(q,10) -
   			 1.632435261169784e-14*pow(m,6)*pow(q,10) + 6.395061906412744e-16*pow(m,8)*pow(q,10) -
   			  8.338874874295073e-18*pow(m,10)*pow(q,10));

	return 6. + (nu2*pow(eps, 2) + nu3*pow(eps, 3) + nu4*pow(eps, 4) + nu5*pow(eps, 5) + nu6*pow(eps, 6)
		+ nu7*pow(eps, 7) + nu8*pow(eps, 8) + nu9*pow(eps, 9) + nu10*pow(eps, 10) + nu11*pow(eps, 11) + nu12*pow(eps, 12));
}

double nu_solver_low_freq_s2_l(double q, int l, int m, double eps){
	double nu2 = -0.5*(24 + 13*l + 28*pow(l,2) + 30*pow(l,3) + 15*pow(l,4))/(l*(1 + l)*(-1 + 2*l)*(1 + 2*l)*(3 + 2*l));
	double nu3 = ((108 + 18*l + 17*pow(l,2) + 3*pow(l,3) + 14*pow(l,4) + 15*pow(l,5) +
       5*pow(l,6))*m*q)/((-1 + l)*pow(l,2)*pow(1 + l,2)*(2 + l)*(-1 + 2*l)*(1 + 2*l)*(3 + 2*l));
	double nu4 = (51840*l + 154656*pow(l,2) - 747792*pow(l,3) - 2705934*pow(l,4) -
     2530951*pow(l,5) + 57648*pow(l,6) + 1950653*pow(l,7) + 2614892*pow(l,8) +
     2753030*pow(l,9) + 2875973*pow(l,10) + 3198765*pow(l,11) +
     2775140*pow(l,12) + 1088535*pow(l,13) - 613935*pow(l,14) -
     1061760*pow(l,15) - 603960*pow(l,16) - 166320*pow(l,17) - 18480*pow(l,18) +
     (233280*pow(l,2) + 1045872*pow(l,3) + 602424*pow(l,4) - 3527460*pow(l,5) -
        6503040*pow(l,6) - 2944886*pow(l,7) + 1561218*pow(l,8) +
        1538990*pow(l,9) - 365930*pow(l,10) - 1087194*pow(l,11) -
        518458*pow(l,12) + 197750*pow(l,13) + 371610*pow(l,14) +
        222208*pow(l,15) + 80816*pow(l,16) + 18720*pow(l,17) + 2080*pow(l,18))*
      pow(q,2) + (-233280 - 1220832*l - 1129248*pow(l,2) + 4533264*pow(l,3) +
        10689684*pow(l,4) + 3358854*pow(l,5) - 9724656*pow(l,6) -
        8311080*pow(l,7) - 413946*pow(l,8) + 1846326*pow(l,9) +
        1147788*pow(l,10) + 335916*pow(l,11) - 230118*pow(l,12) -
        360528*pow(l,13) - 214704*pow(l,14) - 65280*pow(l,15) - 8160*pow(l,16))*
      pow(m,2)*pow(q,2))/(8.*(-1 + l)*pow(l,4)*pow(1 + l,4)*(2 + l)*(-3 + 2*l)*pow(-1 + 2*l,3)*
     pow(1 + 2*l,3)*pow(3 + 2*l,3)*(5 + 2*l));
	double nu5 = ((466560 + 3141504*l + 6737472*pow(l,2) - 1347552*pow(l,3) - 30270684*pow(l,4) -
        51686442*pow(l,5) - 24965217*pow(l,6) + 26877418*pow(l,7) +
        37523604*pow(l,8) + 3377306*pow(l,9) - 16054631*pow(l,10) -
        8016186*pow(l,11) + 430175*pow(l,12) + 1753472*pow(l,13) +
        1032036*pow(l,14) + 519880*pow(l,15) + 310805*pow(l,16) +
        182520*pow(l,17) + 73480*pow(l,18) + 16800*pow(l,19) + 1680*pow(l,20))*
      pow(m,3)*pow(q,3) + m*((-233280*pow(l,2) - 841752*pow(l,3) +
           2993004*pow(l,4) + 14331042*pow(l,5) + 15546708*pow(l,6) -
           4839555*pow(l,7) - 22645143*pow(l,8) - 21618708*pow(l,9) -
           14359584*pow(l,10) - 9231129*pow(l,11) - 4906113*pow(l,12) -
           2088807*pow(l,13) - 2007333*pow(l,14) - 2790444*pow(l,15) -
           1746996*pow(l,16) + 598185*pow(l,17) + 1831641*pow(l,18) +
           1411200*pow(l,19) + 568008*pow(l,20) + 121968*pow(l,21) +
           11088*pow(l,22))*q + (-466560*pow(l,2) - 3141504*pow(l,3) -
           6645132*pow(l,4) + 1657458*pow(l,5) + 29957973*pow(l,6) +
           49056770*pow(l,7) + 22506321*pow(l,8) - 20931941*pow(l,9) -
           26450692*pow(l,10) - 4024518*pow(l,11) + 7358203*pow(l,12) +
           4326322*pow(l,13) + 54234*pow(l,14) - 1289548*pow(l,15) -
           1031636*pow(l,16) - 561663*pow(l,17) - 289039*pow(l,18) -
           134400*pow(l,19) - 44856*pow(l,20) - 8976*pow(l,21) - 816*pow(l,22))*
         pow(q,3)))/(2.*pow(-1 + l,2)*pow(l,6)*pow(1 + l,6)*pow(2 + l,2)*(-3 + 2*l)*
     pow(-1 + 2*l,3)*pow(1 + 2*l,3)*pow(3 + 2*l,3)*(5 + 2*l));
	double nu6 = gen_l_last_term_s2(q, l, m);

	return l + nu2*pow(eps, 2) + nu3*pow(eps, 3) + nu4*pow(eps, 4) + nu5*pow(eps, 5) + nu6*pow(eps, 6);
}

double gen_l_last_term_s2(double q, int l, int m){
	return (-7838208000*pow(l,3) - 45545587200*pow(l,4) + 65591959680*pow(l,5) +
     887604625152*pow(l,6) + 395114098176*pow(l,7) - 8719848407592*pow(l,8) -
     22814668792116*pow(l,9) - 11682887193706*pow(l,10) +
     36281358487429*pow(l,11) + 67490569271634*pow(l,12) +
     33662476594919*pow(l,13) - 28227977657687*pow(l,14) -
     57571453096278*pow(l,15) - 50402266291285*pow(l,16) -
     31904262017498*pow(l,17) - 12370398215822*pow(l,18) +
     4968715105165*pow(l,19) + 13018873771642*pow(l,20) +
     12217387207795*pow(l,21) + 15110431967965*pow(l,22) +
     27544875558400*pow(l,23) + 34615024523815*pow(l,24) +
     20152464395596*pow(l,25) - 6912395124902*pow(l,26) -
     22335331124876*pow(l,27) - 17273201023078*pow(l,28) -
     4205210307144*pow(l,29) + 3278571905224*pow(l,30) + 3497847643136*pow(l,31) +
     1378635779136*pow(l,32) + 122382996736*pow(l,33) - 129089728768*pow(l,34) -
     67464004608*pow(l,35) - 15998062080*pow(l,36) - 1986496512*pow(l,37) -
     104552448*pow(l,38) + (-35271936000*pow(l,4) - 257863046400*pow(l,5) +
        34352112960*pow(l,6) + 4718880172224*pow(l,7) + 12732037653888*pow(l,8) -
        2706064682328*pow(l,9) - 62480580386256*pow(l,10) -
        90542401780824*pow(l,11) + 20900785411920*pow(l,12) +
        179951381622762*pow(l,13) + 162324254590542*pow(l,14) -
        20320687550940*pow(l,15) - 131939190628464*pow(l,16) -
        90084325409730*pow(l,17) - 16988615892198*pow(l,18) +
        11135989134072*pow(l,19) + 18412729702956*pow(l,20) +
        18737963069238*pow(l,21) + 2347552946034*pow(l,22) -
        17664798631068*pow(l,23) - 17934372243144*pow(l,24) -
        1058496916662*pow(l,25) + 11751025587174*pow(l,26) +
        10718356882968*pow(l,27) + 3244891496508*pow(l,28) -
        1598659635024*pow(l,29) - 2181581281200*pow(l,30) -
        1085013510144*pow(l,31) - 241149239424*pow(l,32) + 36435065856*pow(l,33) +
        49067784192*pow(l,34) + 18819827712*pow(l,35) + 4038773760*pow(l,36) +
        485388288*pow(l,37) + 25546752*pow(l,38))*pow(q,2) +
     (-35271936000*pow(l,4) - 323997926400*pow(l,5) - 873322871040*pow(l,6) +
        1012804552224*pow(l,7) + 10459114882128*pow(l,8) +
        18846844992000*pow(l,9) - 9345036979332*pow(l,10) -
        83986490731242*pow(l,11) - 111190565128074*pow(l,12) +
        16882236404715*pow(l,13) + 199506902378427*pow(l,14) +
        186487607475258*pow(l,15) - 38121565176798*pow(l,16) -
        187673721993369*pow(l,17) - 102886332297633*pow(l,18) +
        45203827580418*pow(l,19) + 72179065041834*pow(l,20) +
        13411959110217*pow(l,21) - 22173978422439*pow(l,22) -
        14227631756082*pow(l,23) + 942899436054*pow(l,24) +
        4488283179477*pow(l,25) + 1698721701921*pow(l,26) -
        206309255616*pow(l,27) - 263010717864*pow(l,28) + 52583887632*pow(l,29) +
        104935189488*pow(l,30) + 37438310400*pow(l,31) - 401663232*pow(l,32) -
        4927895808*pow(l,33) - 2072084736*pow(l,34) - 497295360*pow(l,35) -
        80898048*pow(l,36) - 8638464*pow(l,37) - 454656*pow(l,38))*pow(q,4) +
     (-176359680000 - 1619989632000*l - 4412173939200*pow(l,2) +
        4661140619520*pow(l,3) + 50301579674880*pow(l,4) +
        86123536505856*pow(l,5) - 69666921324288*pow(l,6) -
        447071478024864*pow(l,7) - 505034572944696*pow(l,8) +
        339582287581020*pow(l,9) + 1413607907459922*pow(l,10) +
        1136321603296335*pow(l,11) - 573542924461865*pow(l,12) -
        1690954136868115*pow(l,13) - 909494713796445*pow(l,14) +
        522299404621750*pow(l,15) + 852331803506390*pow(l,16) +
        192831370799430*pow(l,17) - 278930750183640*pow(l,18) -
        196409207714325*pow(l,19) + 4700762846115*pow(l,20) +
        54136622645265*pow(l,21) + 22347532657155*pow(l,22) +
        698041971840*pow(l,23) - 2335413621000*pow(l,24) -
        716430654576*pow(l,25) + 153983221968*pow(l,26) + 210740753664*pow(l,27) +
        74944935936*pow(l,28) + 476025600*pow(l,29) - 11443340032*pow(l,30) -
        5609553920*pow(l,31) - 1409525760*pow(l,32) - 192532480*pow(l,33) -
        11325440*pow(l,34))*pow(m,4)*pow(q,4) +
     pow(m,2)*((35271936000*pow(l,2) + 323997926400*pow(l,3) +
           325873031040*pow(l,4) - 5685312556224*pow(l,5) -
           21455158976928*pow(l,6) - 4419261402480*pow(l,7) +
           113267813176440*pow(l,8) + 199755040350960*pow(l,9) -
           71309668642380*pow(l,10) - 543094428674070*pow(l,11) -
           441388863764760*pow(l,12) + 322906463784210*pow(l,13) +
           731285798710410*pow(l,14) + 305107023107220*pow(l,15) -
           218766192285930*pow(l,16) - 264250132685280*pow(l,17) -
           79560689317020*pow(l,18) + 6539396312610*pow(l,19) +
           1910037195540*pow(l,20) + 2344750055730*pow(l,21) +
           17958354125490*pow(l,22) + 20578573406160*pow(l,23) +
           4380033978630*pow(l,24) - 11806562615796*pow(l,25) -
           13092429334932*pow(l,26) - 4186824959520*pow(l,27) +
           2753468704080*pow(l,28) + 3547810788480*pow(l,29) +
           1592055292800*pow(l,30) + 201978470400*pow(l,31) -
           137620224000*pow(l,32) - 83051566080*pow(l,33) -
           21277885440*pow(l,34) - 2810142720*pow(l,35) - 156119040*pow(l,36))*
         pow(q,2) + (211631616000*pow(l,2) + 1943987558400*pow(l,3) +
           5285496810240*pow(l,4) - 5677068207744*pow(l,5) -
           60793472438208*pow(l,6) - 104968699976736*pow(l,7) +
           79520697587520*pow(l,8) + 531621445752096*pow(l,9) +
           613496966639520*pow(l,10) - 361731014172468*pow(l,11) -
           1609065854771956*pow(l,12) - 1302975251649356*pow(l,13) +
           625923966128652*pow(l,14) + 1846744302016280*pow(l,15) +
           949777821767680*pow(l,16) - 560352693187800*pow(l,17) -
           838081861557480*pow(l,18) - 169436355628260*pow(l,19) +
           253315459570860*pow(l,20) + 170539594942980*pow(l,21) -
           621525142020*pow(l,22) - 44224946295840*pow(l,23) -
           20407159448760*pow(l,24) - 1714781080896*pow(l,25) +
           1784481271488*pow(l,26) + 448888813056*pow(l,27) -
           309487668480*pow(l,28) - 235909057536*pow(l,29) -
           52062417920*pow(l,30) + 13198327808*pow(l,31) + 13764888576*pow(l,32) +
           5160288256*pow(l,33) + 1129299968*pow(l,34) + 141557760*pow(l,35) +
           7864320*pow(l,36))*pow(q,4)))/(16.*pow(-1 + l,3)*pow(l,8)*pow(1 + l,8)*pow(2 + l,3)*(-5 + 2*l)*
     pow(-3 + 2*l,2)*pow(-1 + 2*l,5)*pow(1 + 2*l,5)*pow(3 + 2*l,5)*
     pow(5 + 2*l,2)*(7 + 2*l));
}
