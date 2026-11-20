// swsh.c

#include "swsh.hpp"
#include <map>
#include <tuple>

// LAPACK dsyev: eigenvalues (and optionally eigenvectors) of a real symmetric matrix.
// Available via -llapack (conda) or -framework Accelerate (macOS).
// Much faster than gsl_eigen_symm for small matrices (8 µs vs 600 µs for 15×15).
#ifdef PYBHPT_HAS_LAPACK
extern "C" {
    void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda,
                double* w, double* work, int* lwork, int* info);
}
#endif

#define SPECTRAL_NMAX 600
#define SPECTRAL_NMAX_INIT_ADD 15
/* The spectral basis is wide enough once every coupling coefficient beyond some index
inside it has decayed to this fraction of the peak coefficient. Demanding decay far below
double precision (the old 1.e-25) only ever probes round-off. */
#define SPECTRAL_COUPLING_TEST_EPS 1.e-14
#define SPECTRAL_COUPLING_JMAX_EPS 1.e-25
/* Tolerance on how far the coupling vector moves between successive truncations, measured
against the peak coefficient so that two round-off-level tails compare as converged. */
#define SPECTRAL_COUPLING_RESIDUAL_EPS 1.e-14
#define SPECTRAL_COUPLING_CONVERGE_EPS 1.e-25
#define COUPLING_VECTOR_MAX 50
#define ZERO_FREQ_LIMIT 1.e-11

//////////////////////////////////////////////////////
// Jacobi form of the spin-weighted harmonics        //
//////////////////////////////////////////////////////

// sY_lm(theta) = (-1)^s sqrt((2l+1)/(4 pi)) d^l_{m,-s}(theta), and the Wigner d-function is
//
//   d^j_{m'm} = (-1)^lam sqrt(C(2j-k, k+a)/C(k+b, b)) sin^a(theta/2) cos^b(theta/2) P_k^(a,b)(cos theta),
//
// with k = min(j+m, j-m, j+m', j-m') and (a, lam) fixed by which of the four is smallest,
// b = 2j - 2k - a. Here m' = m and m = -s, so a = |m+s| or |m-s| up to the branch, and
// k = l - max(|m|, |s|). The branch depends only on (s, m), so every l in a coupling sum
// shares (a, b) and one recurrence in k serves them all.
//
// This replaces sY = (sum_j A_slj Y_jm)/sin^|s|. That numerator has to cancel down to
// sY*sin^|s| ~ theta^(|s| + |m+s|) near theta = 0, so its relative error grew like
// eps/theta^(|s| + |m+s|): 4e-3 for s = -2, m = 0 at theta = 5e-4, the angles a nearly
// polar orbit samples. The Jacobi form has the vanishing factors explicit and no
// cancellation anywhere on [0, pi].
static void jacobi_parameters(int &a, int &b, int &lam, int s, int m){
	// candidates j + (-s), j - (-s), j + m, j - m, with j dropped (common to all four)
	int mp = m, mm = -s;
	int c[4] = { mm, -mm, mp, -mp };
	int which = 0;
	for(int i = 1; i < 4; i++){ if(c[i] < c[which]){ which = i; } }
	if(which == 0 || which == 3){ a = mp - mm; lam = mp - mm; }
	else{                         a = mm - mp; lam = 0; }
	// b = 2j - 2k - a with k = j + c[which], so b = -2 c[which] - a
	b = -2*c[which] - a;
}

static inline double log_binomial(int n, int r){
	return std::lgamma(n + 1.) - std::lgamma(r + 1.) - std::lgamma(n - r + 1.);
}

// The theta-independent factor of sY_lm, including both signs.
static double jacobi_norm(int a, int b, int lam, int s, int l, int m){
	int k = l - std::max(std::abs(m), std::abs(s));
	double lognorm = 0.5*(log_binomial(2*l - k, k + a) - log_binomial(k + b, b));
	double sign = ((lam + s) % 2 == 0) ? 1. : -1.;
	return sign*sqrt((2.*l + 1.)/(4.*M_PI))*exp(lognorm);
}

// P[n] = P_n^(a,b)(x) for n = 0 .. nmax-1 (three-term recurrence in n, stable on [-1, 1]).
static void jacobi_recurrence(Vector &P, int nmax, int a, int b, double x){
	if(nmax < 1){ return; }
	P[0] = 1.;
	if(nmax < 2){ return; }
	double ab = double(a + b);
	P[1] = 0.5*(double(a - b) + (ab + 2.)*x);
	for(int n = 2; n < nmax; n++){
		double nn = double(n);
		double c = 2.*nn + ab;
		double a1 = 2.*nn*(nn + ab)*(c - 2.);
		double a2 = (c - 1.)*(double(a*a - b*b));
		double a3 = (c - 1.)*c*(c - 2.);
		double a4 = 2.*(nn + a - 1.)*(nn + b - 1.)*c;
		P[n] = ((a2 + a3*x)*P[n - 1] - a4*P[n - 2])/a1;
	}
}

static inline double int_power(double v, int n){
	double out = 1.;
	for(int i = 0; i < n; i++){ out *= v; }
	return out;
}

// sin^a(theta/2) cos^b(theta/2)
static inline double jacobi_prefactor(double sh, double ch, int a, int b){
	return int_power(sh, a)*int_power(ch, b);
}

// d/dtheta [sin^a(theta/2) cos^b(theta/2)] = (a/2) sh^(a-1) ch^(b+1) - (b/2) sh^(a+1) ch^(b-1),
// written without dividing by sh or ch so it is exact at the poles.
static inline double jacobi_prefactor_derivative(double sh, double ch, int a, int b){
	double out = 0.;
	if(a > 0){ out += 0.5*a*int_power(sh, a - 1)*int_power(ch, b + 1); }
	if(b > 0){ out -= 0.5*b*int_power(sh, a + 1)*int_power(ch, b - 1); }
	return out;
}

static double yslm_jacobi(int s, int l, int m, double th){
	if(l < std::max(std::abs(m), std::abs(s))){ return 0.; }
	int a, b, lam;
	jacobi_parameters(a, b, lam, s, m);
	int k = l - std::max(std::abs(m), std::abs(s));
	Vector P(k + 1, 0.);
	jacobi_recurrence(P, k + 1, a, b, cos(th));
	return jacobi_norm(a, b, lam, s, l, m)*jacobi_prefactor(sin(0.5*th), cos(0.5*th), a, b)*P[k];
}

static double yslm_derivative_jacobi(int s, int l, int m, double th){
	if(l < std::max(std::abs(m), std::abs(s))){ return 0.; }
	int a, b, lam;
	jacobi_parameters(a, b, lam, s, m);
	int k = l - std::max(std::abs(m), std::abs(s));
	double sh = sin(0.5*th), ch = cos(0.5*th), x = cos(th);
	Vector P(k + 1, 0.);
	jacobi_recurrence(P, k + 1, a, b, x);
	double out = jacobi_prefactor_derivative(sh, ch, a, b)*P[k];
	if(k > 0){
		// d/dtheta P_k(cos theta) = -sin(theta) P_k'(x), P_k' = (k+a+b+1)/2 P_{k-1}^(a+1,b+1),
		// and sin(theta) = 2 sh ch folds into the prefactor
		Vector Q(k, 0.);
		jacobi_recurrence(Q, k, a + 1, b + 1, x);
		out -= 2.*jacobi_prefactor(sh, ch, a + 1, b + 1)*0.5*double(k + a + b + 1)*Q[k - 1];
	}
	return jacobi_norm(a, b, lam, s, l, m)*out;
}


SpinWeightedHarmonic::SpinWeightedHarmonic(int s, int L, int m, double gamma, const Vector& theta): _s(s), _L(L), _m(m), _gamma(gamma),
	_bcoupling(COUPLING_VECTOR_MAX + L - std::abs(m), 0.), _theta(theta), _Slm(_theta.size(), 0.), _SlmP(_theta.size(), 0.){}
SpinWeightedHarmonic::~SpinWeightedHarmonic(){}

int SpinWeightedHarmonic::getSpinWeight(){ return _s; }
int SpinWeightedHarmonic::getSpheroidalModeNumber(){ return _L; }
int SpinWeightedHarmonic::getAzimuthalModeNumber(){ return _m; }
double SpinWeightedHarmonic::getSpheroidicity(){ return _gamma; }
double SpinWeightedHarmonic::getEigenvalue(){	return _lambda; }
Vector SpinWeightedHarmonic::getCouplingCoefficient(){	return _bcoupling;}
double SpinWeightedHarmonic::getCouplingCoefficient(int l){	return _bcoupling[l - getMinCouplingModeNumber()]; }
int SpinWeightedHarmonic::getMinCouplingModeNumber(){
	return std::abs(_m) < std::abs(_s) ? std::abs(_s) : std::abs(_m);
}
int SpinWeightedHarmonic::getMaxCouplingModeNumber(){
	return _bcoupling.size() + getMinCouplingModeNumber() - 1;
}

int SpinWeightedHarmonic::generateSolutionsAndDerivatives(){
	generateCouplingCoefficients();
	fillYlmGrid();
	for(size_t i = 0; i < _Slm.size(); i++){
		_Slm[i] = SslmGrid(int(i));
	}
	for(size_t i = 0; i < _SlmP.size(); i++){
		_SlmP[i] = SslmDerivativeGrid(int(i));
	}

	return 0;
}

int SpinWeightedHarmonic::generateCouplingCoefficients(){
	if(_bcoupling[_L - getMinCouplingModeNumber()] == 0 && _bcoupling[_L - getMinCouplingModeNumber() + 1] == 0){
		// spectral_solver returns nonzero on success/stall and 0 when the convergence
		// test failed at the truncation limit; record the latter as a failure status.
		_couplingStatus = spectral_solver(_s, _L, _m, _gamma, _lambda, _bcoupling) == 0 ? 1 : 0;
	}
	return _couplingStatus;
}

int SpinWeightedHarmonic::getCouplingStatus(){ return _couplingStatus; }

int SpinWeightedHarmonic::generateSolutions(){
	generateCouplingCoefficients();
	fillYlmGrid();
	for(size_t i = 0; i < _Slm.size(); i++){
		_Slm[i] = SslmGrid(int(i));
	}

	return 0;
}

int SpinWeightedHarmonic::generateDerivatives(){
	generateCouplingCoefficients();
	fillYlmGrid();
	for(size_t i = 0; i < _SlmP.size(); i++){
		_SlmP[i] = SslmDerivativeGrid(int(i));
	}

	return 0;
}

Vector SpinWeightedHarmonic::getArguments(){ return _theta; }
Vector SpinWeightedHarmonic::getSolution(){	return _Slm; }
Vector SpinWeightedHarmonic::getDerivative(){ return _SlmP; }
Vector SpinWeightedHarmonic::getSecondDerivative(){
	Vector SlmPP(_theta.size(), 0.);
	for(size_t i = 0; i < SlmPP.size(); i++){
		SlmPP[i] = Sslm_secondDerivative(_s, _L, _m, _gamma, _lambda, _theta[i], _Slm[i], _SlmP[i]);
	}
	return SlmPP;
}

double SpinWeightedHarmonic::getArguments(int pos){ return _theta[pos]; }
double SpinWeightedHarmonic::getSolution(int pos){
	return _Slm[pos];
}
double SpinWeightedHarmonic::getDerivative(int pos){
	return _SlmP[pos];
}
double SpinWeightedHarmonic::getSecondDerivative(int pos){
	return Sslm_secondDerivative(_s, _L, _m, _gamma, _lambda, _theta[pos], _Slm[pos], _SlmP[pos]);
}

/////////////////////
// Spectral matrix //
/////////////////////

// Matrix coefficients
double k1(const int &s, const int &l, const int &j, const int &m){
	double k = sqrt((2*l + 1.)/(2*j + 1.));
	k *= clebsch(l, 1, j, m, 0, m);
	k *= clebsch(l, 1, j, -s, 0, -s);

	return k;
}

double k2(const int &s, const int &l, const int &j, const int &m){
	double k = 2*sqrt((2*l + 1.)/(2*j + 1.))/3.;
	k *= clebsch(l, 2, j, m, 0, m);
	k *= clebsch(l, 2, j, -s, 0, -s);
	if( j == l ){
		k += 1/3.;
	}

	return k;
}

// m(i, i-2)
double akm2(const int &s, const int &l, const int &m, const double &g){
	if( l < 0 || std::abs(m) > l || std::abs(s) > l ){
		return 0;
	}
	return -g*g*k2(s, l-2, l, m);
}

// m(i, i-1)
double akm1(const int &s, const int &l, const int &m, const double &g){
	if( l < 0 || std::abs(m) > l || std::abs(s) > l ){
		return 0;
	}
	return -g*g*k2(s, l-1, l, m) + 2*s*g*k1(s, l-1, l, m);
}

// m(i, i)
double akp0(const int &s, const int &l, const int &m, const double &g){
	if( l < 0 || std::abs(m) > l || std::abs(s) > l ){
		return 0;
	}
	return -g*g*k2(s, l, l, m) + 2*s*g*k1(s, l, l, m) + l*(l + 1) - s*(s+1) - 2*m*g + g*g;
}

// m(i, i+1)
double akp1(const int &s, const int &l, const int &m, const double &g){
	if( l < 0 || std::abs(m) > l + 1 || std::abs(s) > l + 1 ){
		return 0;
	}
	return -g*g*k2(s, l+1, l, m) + 2*s*g*k1(s, l+1, l, m);
}

// m(i, i+2)
double akp2(const int &s, const int &l, const int &m, const double &g){
	if( l < 0 || std::abs(m) > l + 2 || std::abs(s) > l + 2 ){
		return 0;
	}
	return -g*g*k2(s, l+2, l, m);
}

double spectral_weight(const double &g){
	double weight;
	if( g > 1 ){
		weight = pow(g, -1);
	}else{
		weight = 1.;
	}
	return weight;
}

// Spectral matrix
int spectral_matrix(const int &s, const int &lmin, const int &m, const double &g, gsl_matrix* mat){
	size_t nmax = mat->size1;
	if( nmax != mat->size2 ) return 1; // error, not a NxN matrix
	double weight = spectral_weight(g);

	for(size_t i = 0; i < nmax; i++){
		if(i >= 2){gsl_matrix_set(mat, i, i - 2, weight*akm2(s, lmin + i, m, g));}
		if(i >= 1){gsl_matrix_set(mat, i, i - 1, weight*akm1(s, lmin + i, m, g));}
		gsl_matrix_set(mat, i, i, weight*akp0(s, lmin + i, m, g));
		if(i < nmax - 1){gsl_matrix_set(mat, i, i + 1, weight*akp1(s, lmin + i, m, g));}
		if(i < nmax - 2){gsl_matrix_set(mat, i, i + 2, weight*akp2(s, lmin + i, m, g));}
	}
	return 0;
}

int spectral_matrix_sparse_init(const int &s, const int &lmin, const int &m, const double &g, gsl_spmatrix* mat){
	size_t nmax = mat->size1;
	if( nmax != mat->size2 ) return 1; // error, not a NxN matrix
	double weight = spectral_weight(g);

	for(size_t i = 0; i < nmax; i++){
		if(i >= 2){gsl_spmatrix_set(mat, i, i - 2, weight*akm2(s, lmin + i, m, g));}
		if(i >= 1){gsl_spmatrix_set(mat, i, i - 1, weight*akm1(s, lmin + i, m, g));}
		gsl_spmatrix_set(mat, i, i, weight*akp0(s, lmin + i, m, g));
		if(i < nmax - 1){gsl_spmatrix_set(mat, i, i + 1, weight*akp1(s, lmin + i, m, g));}
		if(i < nmax - 2){gsl_spmatrix_set(mat, i, i + 2, weight*akp2(s, lmin + i, m, g));}
	}
	return 0;
}

int spectral_matrix_sparse(const int &s, const int &lmin, const int &m, const double &g, gsl_spmatrix* mat, const size_t &nmax){
	size_t dmax = mat->size1;
	if( dmax != mat->size2 ) return 1; // error, not a NxN matrix
	double weight = spectral_weight(g);
	mat->size1 = mat->size2 = nmax;

	for(size_t i = dmax; i < nmax; i++){
		if(i >= 2){gsl_spmatrix_set(mat, i, i - 2, weight*akm2(s, lmin + i, m, g));}
		if(i >= 1){gsl_spmatrix_set(mat, i, i - 1, weight*akm1(s, lmin + i, m, g));}
		gsl_spmatrix_set(mat, i, i, weight*akp0(s, lmin + i, m, g));
		if(i < nmax - 1){gsl_spmatrix_set(mat, i, i + 1, weight*akp1(s, lmin + i, m, g));}
		if(i < nmax - 2){gsl_spmatrix_set(mat, i, i + 2, weight*akp2(s, lmin + i, m, g));}
	}

	return 0;
}

// Solve eigensystem of spectral matrix
double spectral_solver(const int &s, const int &l, const int &m, const double &g, const unsigned int &nmax){
	gsl_vector* la = gsl_vector_alloc(nmax);
	int error = spectral_solver_n(s, l, m, g, la);
	if(error) return 0;

	unsigned int lmin = std::max(std::abs(s), std::abs(m));
	double la_lm = gsl_vector_get(la, l - lmin);
	gsl_vector_free(la);

	return la_lm;
}

int spectral_solver(const int &s, const int &l, const int &m, const double &g, double& lambda, Vector& bvec){
	int lmin = std::max(std::abs(s), std::abs(m));
	int nmax = l - lmin + SPECTRAL_NMAX_INIT_ADD;
	if(g == 0 || std::abs(g) < ZERO_FREQ_LIMIT){
		lambda = l*(l + 1) - s*(s + 1);
		bvec[l - lmin] = 1.;
		return 1;
	}
	int error;
	coupling_converge bkData;
	coupling_test test;
	bkData.jmax = bkData.testIndex = l - lmin;
	bkData.jmax_err = bkData.test_err = 10.;

	// initialize spectral matrix using a sparse matrix representation
	gsl_spmatrix* mat = gsl_spmatrix_alloc_nzmax(nmax, nmax, SPECTRAL_NMAX, GSL_SPMATRIX_COO);
	spectral_matrix_sparse_init(s, lmin, m, g, mat);

	gsl_matrix* bmat = gsl_matrix_alloc(nmax, nmax);
	gsl_vector* la = gsl_vector_alloc(nmax);
	error = spectral_solver_n(s, l, m, g, la, bmat, mat);

	nmax += 2;
	spectral_matrix_sparse(s, lmin, m, g, mat, nmax); // update spectral matrix to include additional entries
	gsl_matrix* bmat2 = gsl_matrix_alloc(nmax, nmax);
	gsl_vector* la2 = gsl_vector_alloc(nmax);
	error = spectral_solver_n(s, l, m, g, la2, bmat2, mat);
	test = spherical_spheroidal_coupling_convergence_test(s, l, m, g, bmat, bmat2, bkData);

	while(nmax < SPECTRAL_NMAX - 10 && test == FAIL){
		gsl_matrix_free(bmat);
		gsl_vector_free(la);
		bmat = bmat2;
		la = la2;

		nmax += 10;
		spectral_matrix_sparse(s, lmin, m, g, mat, nmax);
		bmat2 = gsl_matrix_alloc(nmax, nmax);
		la2 = gsl_vector_alloc(nmax);

		error = spectral_solver_n(s, l, m, g, la2, bmat2, mat);
		test = spherical_spheroidal_coupling_convergence_test(s, l, m, g, bmat, bmat2, bkData);
	}
	gsl_spmatrix_free(mat);

	if(test == SUCCESS){ // if convergence test was successful, return data from the most resolved (highest nmax) spectral eigenvalue problem
		gsl_matrix_free(bmat);
		gsl_vector_free(la);
		bmat = bmat2;
		la = la2;
	}else if(test == FAIL){
		gsl_matrix_free(bmat); gsl_matrix_free(bmat2);
		gsl_vector_free(la); gsl_vector_free(la2);
		// Convergence failed at the truncation limit; the caller surfaces this to Python
		// as a warning via the coupling status rather than printing here.
		return 0;
	}else{ // if convergence test was stalled, return data from the second most resolved (second highest nmax) spectral eigenvalue problem
		// std::cout << "(SWSH) Calculation of coupling coefficients stalled for (s, l, m, gamma) = ("<<s<<", "<<l<<", "<<m<<", "<<g<<") \n";
		gsl_matrix_free(bmat2);
		gsl_vector_free(la2);
		nmax = la->size;
	}

	int rescale = 1;
	if( gsl_matrix_get(bmat, l - lmin, l - lmin) < 0 ){
		rescale = -1;
	}

	int i = 0;
	int bvecSize = bvec.size();
	bvecSize = nmax < bvecSize ? nmax : bvecSize;

	while(i < l - lmin + 2 && i < bvecSize){
		bvec[i] = rescale*gsl_matrix_get(bmat, i, l - lmin);
		i++;
	}
	double bError = std::abs(bvec[i - 2]);
	while(i < bvecSize - 1 && bError > SPECTRAL_COUPLING_CONVERGE_EPS){
		bvec[i] = rescale*gsl_matrix_get(bmat, i, l - lmin);
		bvec[i + 1] = rescale*gsl_matrix_get(bmat, i + 1, l - lmin);
		bError = std::abs(bvec[i]);
		i+=2;
	}
	lambda = gsl_vector_get(la, l - lmin);

	gsl_matrix_free(bmat);
	gsl_vector_free(la);

	return 1;
}

int spectral_solver_n(const int &s, const int &l, const int &m, const double &g, gsl_vector* la){
	int nmax = la->size;
	gsl_matrix* bmat = gsl_matrix_alloc(nmax, nmax);
	int error = spectral_solver_n(s, l, m, g, la, bmat);
	if( error ) return 1;
	gsl_matrix_free(bmat);

	return 0;
}

int spectral_solver_n(const int &s, const int &l, const int &m, const double &g, gsl_vector* la, gsl_matrix* bmat){
	int lmin = std::max(std::abs(s), std::abs(m));
	int nmax = la->size;
	if( nmax < l-lmin ){
		return 1; // error, need larger nmax tolerance
	}

	double weight = spectral_weight(g);

#ifdef PYBHPT_HAS_LAPACK
	// Use bmat directly as working buffer: zero it, fill with the spectral matrix,
	// run dsyev in-place, then transpose so eigenvectors end up as columns (GSL convention).
	// LAPACK dsyev writes eigenvectors as columns in Fortran column-major order; in the
	// C row-major layout of bmat->data this appears as eigenvectors stored as rows.
	// A subsequent in-place transpose restores the GSL convention used by the rest of the code.
	gsl_matrix_set_zero(bmat);
	int error = spectral_matrix(s, lmin, m, g, bmat);
	if( error ) return 1;

	int n = nmax;
	// Workspace query for optimal lwork
	int lwork_query = -1;
	double work_opt = 0.;
	int info = 0;
	char jobz = 'V', uplo = 'U';
	dsyev_(&jobz, &uplo, &n, bmat->data, &n, la->data, &work_opt, &lwork_query, &info);
	int lwork = std::max((int)work_opt, 3*n + 64);
	std::vector<double> work(lwork);
	dsyev_(&jobz, &uplo, &n, bmat->data, &n, la->data, work.data(), &lwork, &info);
	if( info != 0 ) return 1;
	// dsyev returns eigenvalues in ascending order; no explicit sort needed.
	// Transpose in-place to convert LAPACK's row-stored eigenvectors to GSL column convention.
	gsl_matrix_transpose(bmat);
	gsl_vector_scale(la, 1./weight);
#else
	gsl_matrix* specMat = gsl_matrix_calloc(nmax, nmax);
	int error = spectral_matrix(s, lmin, m, g, specMat);
	if( error ) return 1;

	/* while not completely obvious from the way we have expressed the coefficients
	for our sparse matrix, it is in fact real, symmetric (one can verify this through
	the symmetries of the Wigner 3j-symbol), simplifying our calculation. */
	gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(nmax);
	gsl_eigen_symmv(specMat, la, bmat, w);
	gsl_eigen_symmv_sort(la, bmat, GSL_EIGEN_SORT_VAL_ASC);
	gsl_eigen_symmv_free(w);
	gsl_matrix_free(specMat);

	gsl_vector_scale(la, 1/weight);
#endif

	return 0;
}

int spectral_solver_n(const int &s, const int &l, const int &m, const double &g, gsl_vector* la, gsl_matrix* bmat, gsl_spmatrix* mat){
	int lmin = std::max(std::abs(s), std::abs(m));
	int nmax = la->size;
	if( nmax < l-lmin ){
		return 1; // error, need larger nmax tolerance
	}

	double weight = spectral_weight(g);

#ifdef PYBHPT_HAS_LAPACK
	// Zero bmat first: gsl_spmatrix_sp2d only sets non-zero entries, leaving the
	// remaining positions untouched. The matrix must be zeroed to avoid garbage
	// in off-band positions that would corrupt the eigenvalue solve.
	gsl_matrix_set_zero(bmat);
	int error = gsl_spmatrix_sp2d(bmat, mat);
	if( error ) return 1;

	int n = nmax;
	int lwork_query = -1;
	double work_opt = 0.;
	int info = 0;
	char jobz = 'V', uplo = 'U';
	dsyev_(&jobz, &uplo, &n, bmat->data, &n, la->data, &work_opt, &lwork_query, &info);
	int lwork = std::max((int)work_opt, 3*n + 64);
	std::vector<double> work(lwork);
	dsyev_(&jobz, &uplo, &n, bmat->data, &n, la->data, work.data(), &lwork, &info);
	if( info != 0 ) return 1;
	gsl_matrix_transpose(bmat);
	gsl_vector_scale(la, 1./weight);
#else
	gsl_matrix* specMat = gsl_matrix_calloc(nmax, nmax);
	int error = gsl_spmatrix_sp2d(specMat, mat);
	if( error ) return 1;

	/* while not completely obvious from the way we have expressed the coefficients
	for our sparse matrix, it is in fact real, symmetric (one can verify this through
	the symmetries of the Wigner 3j-symbol), simplifying our calculation. */
	gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(nmax);
	gsl_eigen_symmv(specMat, la, bmat, w);
	gsl_eigen_symmv_sort(la, bmat, GSL_EIGEN_SORT_VAL_ASC);
	gsl_eigen_symmv_free(w);
	gsl_matrix_free(specMat);

	gsl_vector_scale(la, 1/weight);
#endif

	return 0;
}

// Compute eigenvalues only (no eigenvectors) from an existing sparse matrix.
// When LAPACK is available (PYBHPT_HAS_LAPACK), uses dsyev_ which is ~60-100x
// faster than gsl_eigen_symm for the small matrices we encounter here (15-35×15-35).
// Falls back to gsl_eigen_symm when LAPACK is unavailable.
// Note: gsl_eigen_symm only uses the lower triangular part, which is consistent
// with how spectral_matrix_sparse fills both triangles of the COO sparse matrix.
static int spectral_solver_eigenvalues_n(const int &s, const int &l, const int &m, const double &g, gsl_vector* la, gsl_spmatrix* mat){
	int lmin = std::max(std::abs(s), std::abs(m));
	int nmax = la->size;
	if( nmax < l - lmin ){
		return 1;
	}

	gsl_matrix* specMat = gsl_matrix_calloc(nmax, nmax);
	int error = gsl_spmatrix_sp2d(specMat, mat);
	if( error ) { gsl_matrix_free(specMat); return 1; }

	double weight = spectral_weight(g);

#ifdef PYBHPT_HAS_LAPACK
	// LAPACK dsyev: eigenvalues only ('N'), use upper triangle ('U').
	// For a symmetric matrix, row-major == col-major (A = A^T), so passing
	// the GSL row-major buffer to Fortran LAPACK is always correct here.
	// dsyev returns eigenvalues in la->data in ascending order.
	int n = nmax;
	int lwork = std::max(1, 3*n - 1);  // minimum workspace for N mode
	std::vector<double> work(lwork);
	int info = 0;
	char jobz = 'N', uplo = 'U';
	dsyev_(&jobz, &uplo, &n, specMat->data, &n, la->data, work.data(), &lwork, &info);
	if( info != 0 ){
		// dsyev did not converge. Rebuild the matrix (dsyev overwrote it) and fall back to
		// GSL: returning here would leave la unwritten, and callers read it either way.
		gsl_matrix_set_zero(specMat);
		if( gsl_spmatrix_sp2d(specMat, mat) ){ gsl_matrix_free(specMat); return 1; }
		gsl_eigen_symm_workspace* w = gsl_eigen_symm_alloc(nmax);
		int gslError = gsl_eigen_symm(specMat, la, w);
		gsl_eigen_symm_free(w);
		gsl_matrix_free(specMat);
		if( gslError ) return 1;
		std::sort(la->data, la->data + nmax);
		gsl_vector_scale(la, 1.0/weight);
		return 0;
	}
	gsl_matrix_free(specMat);
	// eigenvalues already in ascending order; just unscale
	gsl_vector_scale(la, 1.0/weight);
#else
	// Fallback: gsl_eigen_symm (~60-100x slower than LAPACK for small n)
	gsl_eigen_symm_workspace* w = gsl_eigen_symm_alloc(nmax);
	gsl_eigen_symm(specMat, la, w);
	// gsl_eigen_symm has no sort function; sort manually ascending.
	std::sort(la->data, la->data + nmax);
	gsl_eigen_symm_free(w);
	gsl_matrix_free(specMat);
	gsl_vector_scale(la, 1.0/weight);
#endif

	return 0;
}

coupling_test spherical_spheroidal_coupling_convergence_test(const int &s, const int &l, const int &m, const double &g, gsl_matrix* bmat, gsl_matrix* bmat2, coupling_converge &b_data){
	int lmin = std::max(std::abs(s), std::abs(m));
	int dim = bmat->size1;
	int dim2 = bmat2->size2;
	gsl_vector* bcol = gsl_vector_alloc(dim);
	gsl_vector* bcol2 = gsl_vector_alloc(dim2);
	double relerror = b_data.test_err;

	gsl_matrix_get_col(bcol, bmat, l - lmin);
	if( gsl_vector_get(bcol, l - lmin) < 0 ){
		gsl_vector_scale(bcol, -1.);
	}
	gsl_matrix_get_col(bcol2, bmat2, l - lmin);
	if( gsl_vector_get(bcol2, l - lmin) < 0 ){
		gsl_vector_scale(bcol2, -1.);
	}

	double bmax = gsl_vector_max(bcol);
	double bmin = gsl_vector_min(bcol);
	double norm = std::abs(bmax) > std::abs(bmin) ? std::abs(bmax) : std::abs(bmin);
	if(norm == 0.){
		norm = 1.;
	}
	while(b_data.jmax < dim - 1 && std::abs(gsl_vector_get(bcol, b_data.jmax)/norm) > SPECTRAL_COUPLING_JMAX_EPS){
		if(s == 0){
			b_data.jmax += 2;
		}else{
			b_data.jmax++;
		}
	}
	if(b_data.jmax == dim){
		if(s == 0){
			b_data.jmax -= 2;
		}else{
			b_data.jmax--;
		}
	}

	/* Walk in from the truncation edge to the first index past which every coupling
	coefficient has decayed below SPECTRAL_COUPLING_TEST_EPS of the peak. Scanning outward
	from l - lmin instead would stop at the first sign change of an oscillating tail, which
	at large spheroidicity sits far inside a band that is still growing. */
	b_data.testIndex = dim;
	while(b_data.testIndex > 0 && std::abs(gsl_vector_get(bcol, b_data.testIndex - 1)/norm) <= SPECTRAL_COUPLING_TEST_EPS){
		b_data.testIndex--;
	}
	bool tailResolved = (b_data.testIndex < dim);
	if(b_data.testIndex > dim - 1){
		b_data.testIndex = dim - 1;
	}

	double jmax_denom = gsl_vector_get(bcol2, b_data.jmax);
	if(jmax_denom == 0){
		jmax_denom = DBL_EPSILON*pow(g, 1);
	}

	double jmax_num = gsl_vector_get(bcol, b_data.jmax);
	if(jmax_num == 0){
		jmax_num = DBL_EPSILON*pow(g, 1);
	}

	b_data.jmax_err = std::abs(1 - jmax_num/jmax_denom);
	/* Measure how much widening the basis moved the coupling vector, on the scale of the
	peak coefficient. The old test took the relative change of a single coefficient chosen
	where it had already fallen to 1e-25 of the peak: a ratio of two round-off-level
	numbers, decided by the low bits of the input, so it failed at random and drove the
	truncation loop upward with nothing left to gain. */
	b_data.test_err = 0.;
	for(int i = 0; i < dim; i++){
		double diff = std::abs(gsl_vector_get(bcol2, i) - gsl_vector_get(bcol, i))/norm;
		if(diff > b_data.test_err){
			b_data.test_err = diff;
		}
	}

	double convergenceCriteria = SPECTRAL_COUPLING_RESIDUAL_EPS;
	if(g > 1.){
		convergenceCriteria *= pow(g, 2);
	}

	gsl_vector_free(bcol);
	gsl_vector_free(bcol2);

	if( tailResolved && b_data.test_err < convergenceCriteria ){
		return SUCCESS;
	}else if( tailResolved && b_data.test_err >= relerror ){
		/* Widening the basis stopped moving the coupling vector any closer: stop
		regardless of whether the test index moved, and keep the previous truncation. */
		return STALL;
	}else{
		return FAIL;
	}
}

////////////////////
// Spin Couplings //
////////////////////

// Coupling between scalar Yjm and spin-weighted harmonics sYlm
double Asljm(const int &s, const int &l, const int &j, const int &m){
	// (s,l,j,m)-only coupling (Wigner-3j) -- independent of theta, but Yslm() calls this
	// per grid point. Memoize so the 3j / factorial work runs once per distinct key.
	static thread_local std::map<std::tuple<int,int,int,int>, double> cache;
	// Bound the table: entries are never invalidated, so a long-lived process sweeping a
	// wide mode range would otherwise grow it without limit. Refilling is cheap.
	if(cache.size() > 200000) cache.clear();
	auto key = std::make_tuple(s, l, j, m);
	auto it = cache.find(key);
	if(it != cache.end()) return it->second;

	double aslmg;
	if( std::abs(l-j) > std::abs(s) || j < std::abs(m) ){
		aslmg = 0.;
	}else{
		int lmin = std::abs(m) < std::abs(s) ? std::abs(s) : std::abs(m);
		if( l < lmin ){
			aslmg = 0.;
		}else{
			aslmg = pow(-1., m + s*(1 + sgn<int>(s))/2);
			aslmg *= sqrt(pow(4, std::abs(s))*pow(factorial(std::abs(s)), 2)*(2*j + 1)*(2*l + 1)/factorial(std::abs(2*s)));
			aslmg *= w3j(std::abs(s), l, j, 0, m, -m);
			aslmg *= w3j(std::abs(s), l, j, s, -s, 0);
		}
	}
	cache[key] = aslmg;
	return aslmg;
}

double dAsljm(const int &s, const int &l, const int &j, const int &m){
	return double(j + 2 + std::abs(s))*clm(j + 1, m)*Asljm(s, l, j + 1, m)
		- double(j - 1 - std::abs(s))*clm(j, m)*Asljm(s, l, j - 1, m);
}

// Clebsch-Gordan coefficients
double clebsch(const int &j1, const int &j2, const int &j, const int &m1, const int &m2, const int &m){
	double cg = pow(-1., j1 - j2 + m);
	cg *= sqrt(2*j + 1.);
	cg *= w3j(j1, j2, j, m1, m2, -m);

	return cg;
}

// Wigner 3J symbol
double w3j(const int &j1, const int &j2, const int &j, const int &m1, const int &m2, const int &m){
	if(j1 < std::abs(m1) || j2 < std::abs(m2) || j < std::abs(m)){
		return 0.;
	}
	return gsl_sf_coupling_3j(2*j1, 2*j2, 2*j, 2*m1, 2*m2, 2*m);
}

//////////////////////////////////////////////
// Spin-Weighted Spherical Harmonics (SWSH) //
//////////////////////////////////////////////

// SWSHs
Complex Sslm(const int &s, const int &l, const int &m, const double &g, const double &th, const double &ph){
	return Sslm(s, l, m, g, th)*exp(I*Complex(m)*ph);
}

// double Sslm(const int &s, const int &l, const int &m, const double &g, const double &th){
// 	unsigned int lmin = std::max(std::abs(s), std::abs(m));
// 	unsigned int nmax = l - lmin + SPECTRAL_NMAX_INIT_ADD;
// 	int error;
// 	double swsh = 0;
// 	coupling_converge bkData;
// 	coupling_test test;
// 	bkData.jmax = bkData.testIndex = l - lmin;
// 	bkData.jmax_err = bkData.test_err = 10.;

// 	// initialize spectral matrix using a sparse matrix representation
// 	gsl_spmatrix* mat = gsl_spmatrix_alloc_nzmax(nmax, nmax, SPECTRAL_NMAX, GSL_SPMATRIX_COO);
// 	spectral_matrix_sparse_init(s, lmin, m, g, mat);

// 	gsl_matrix* bmat = gsl_matrix_alloc(nmax, nmax);
// 	gsl_vector* la = gsl_vector_alloc(nmax);
// 	error = spectral_solver_n(s, l, m, g, la, bmat, mat);

// 	nmax += 2;
// 	spectral_matrix_sparse(s, lmin, m, g, mat, nmax); // update spectral matrix to include additional entries
// 	gsl_matrix* bmat2 = gsl_matrix_alloc(nmax, nmax);
// 	gsl_vector* la2 = gsl_vector_alloc(nmax);
// 	error = spectral_solver_n(s, l, m, g, la2, bmat2, mat);
// 	test = spherical_spheroidal_coupling_convergence_test(s, l, m, g, bmat, bmat2, bkData);

// 	while(nmax < SPECTRAL_NMAX && test == FAIL){
// 		gsl_matrix_free(bmat);
// 		gsl_vector_free(la);
// 		bmat = bmat2;
// 		la = la2;

// 		nmax += 10;
// 		spectral_matrix_sparse(s, lmin, m, g, mat, nmax);
// 		bmat2 = gsl_matrix_alloc(nmax, nmax);
// 		la2 = gsl_vector_alloc(nmax);

// 		error = spectral_solver_n(s, l, m, g, la2, bmat2, mat);
// 		test = spherical_spheroidal_coupling_convergence_test(s, l, m, g, bmat, bmat2, bkData);
// 	}
// 	gsl_spmatrix_free(mat);
// 	if(nmax == SPECTRAL_NMAX){
// 		std::cout << "(SWSH) Max number of iterations executed for spectral solver. \n";
// 	}else{
// 		std::cout << "(SWSH) "<< nmax <<" iterations executed for spectral solver. \n";
// 	}

// 	if(test == SUCCESS){ // if convergence test was successful, return data from the most resolved (highest nmax) spectral eigenvalue problem
// 		gsl_matrix_free(bmat);
// 		gsl_vector_free(la);
// 		bmat = bmat2;
// 		la = la2;
// 	}else if(test == FAIL){
// 		gsl_matrix_free(bmat); gsl_matrix_free(bmat2);
// 		gsl_vector_free(la); gsl_vector_free(la2);

// 		return 0;
// 	}else{ // if convergence test was stalled, return data from the second most resolved (second highest nmax) spectral eigenvalue problem
// 		gsl_matrix_free(bmat2);
// 		gsl_vector_free(la2);
// 		nmax = la->size;
// 	}

// 	gsl_vector* bcol = gsl_vector_alloc(nmax);
// 	gsl_matrix_get_col(bcol, bmat, l - lmin);
// 	if( gsl_vector_get(bcol, l - lmin) < 0 ){
// 		gsl_vector_scale(bcol, -1.);
// 	}
// 	gsl_matrix_free(bmat);
// 	gsl_vector_free(la);

// 	double bj;
// 	for(int i = 0; i < bkData.jmax; i++){
// 		bj = gsl_vector_get(bcol, i);
// 		// std::cout << "b_"<< lmin + i << " = " << bj << "\n";
// 		swsh += bj*Yslm(s, lmin + i, m, th);
// 	}

// 	if(error) return 0.;

// 	return swsh;
// }

double Sslm(const int &s, const int &l, const int &m, const double &g, const double &th){
	Vector bvec(COUPLING_VECTOR_MAX + l - std::abs(m), 0.);
	double lambda;
	spectral_solver(s, l, m, g, lambda, bvec);
	return Sslm(s, l, m, g, bvec, th);
}

double Sslm(const int &s, const int &l, const int &m, const double &, const Vector& bvec, const double &th){
	int lmin = std::max(std::abs(s), std::abs(m));
	int i = l - lmin;

	double swsh = bvec[i]*Yslm(s, lmin + i, m, th);
	double comp = swsh;
	swsh += bvec[i + 1]*Yslm(s, lmin + i + 1, m, th);
	i--;

	double error = 1.;
	while(i > 0 && error > 1.e-14){
		comp = swsh;
		swsh += bvec[i]*Yslm(s, lmin + i, m, th);
		swsh += bvec[i - 1]*Yslm(s, lmin + i - 1, m, th);
		error = std::abs(1. - comp/swsh);
		i -= 2;
	}
	if(i == 0){
		swsh += bvec[i]*Yslm(s, lmin + i, m, th);
	}

	i = l - lmin + 2;
	error = 1.;
	int imax = bvec.size() - 1;
	while(i < imax && error > 1.e-14){
		comp = swsh;
		swsh += bvec[i]*Yslm(s, lmin + i, m, th);
		swsh += bvec[i + 1]*Yslm(s, lmin + i + 1, m, th);
		error = std::abs(1. - comp/swsh);
		i += 2;
	}
	if(i == imax){
		swsh += bvec[i]*Yslm(s, lmin + i, m, th);
	}

	return swsh;
}

double Sslm_derivative(const int &s, const int &l, const int &m, const double &, const Vector& bvec, const double &th){
	unsigned int lmin = std::max(std::abs(s), std::abs(m));
	int i = l - lmin;

	double comp = 0.;
	double swsh = bvec[i]*Yslm_derivative(s, lmin + i, m, th);
	swsh += bvec[i + 1]*Yslm_derivative(s, lmin + i + 1, m, th);
	double error = std::abs(1. - comp/swsh);
	i--;

	while(i > 0 && error > 1.e-14){
		comp = swsh;
		swsh += bvec[i]*Yslm_derivative(s, lmin + i, m, th);
		swsh += bvec[i - 1]*Yslm_derivative(s, lmin + i - 1, m, th);
		error = std::abs(1. - comp/swsh);
		i -= 2;
	}
	if(i == 0){
		swsh += bvec[i]*Yslm_derivative(s, lmin + i, m, th);
	}

	i = l - lmin + 2;
	error = 1.;
	int imax = bvec.size() - 1;
	while(i < imax && error > 1.e-14){
		comp = swsh;
		swsh += bvec[i]*Yslm_derivative(s, lmin + i, m, th);
		swsh += bvec[i + 1]*Yslm_derivative(s, lmin + i + 1, m, th);
		error = std::abs(1. - comp/swsh);
		i += 2;
	}
	if(i == imax){
		swsh += bvec[i]*Yslm_derivative(s, lmin + i, m, th);
	}

	return swsh;
}

double Sslm_secondDerivative(const int &s, const int &, const int &m, const double &g, const double& lambda, const double &th, const double &Slm, const double &SlmP){
	return -cos(th)*SlmP/sin(th) + (pow(g*sin(th), 2) + pow((m + s*cos(th))/sin(th), 2) + 2.*s*g*cos(th) - s - 2.*m*g - lambda)*Slm;
}

// SWSH Eigenvalues
// Optimized version: builds the spectral matrix once using the sparse
// incremental approach (adding only new rows on each extension) and uses
// gsl_eigen_symm (eigenvalues only, ~2x faster than gsl_eigen_symmv).
// This avoids rebuilding the full matrix from scratch on each convergence
// iteration, which was the dominant cost in the old implementation.
double swsh_eigenvalue(const int &s, const int &l, const int &m, const double &g){
	unsigned int lmin = std::max(std::abs(s), std::abs(m));

	// Spherical limit: eigenvalue is exact, no matrix needed.
	if(g == 0 || std::abs(g) < ZERO_FREQ_LIMIT){
		return l*(l + 1) - s*(s + 1);
	}

	unsigned int nmax = l - lmin + SPECTRAL_NMAX_INIT_ADD;

	// Allocate sparse matrix with enough storage for SPECTRAL_NMAX non-zeros.
	// The pentadiagonal structure means ~5*SPECTRAL_NMAX entries at most.
	gsl_spmatrix* mat = gsl_spmatrix_alloc_nzmax(nmax, nmax, 5*SPECTRAL_NMAX, GSL_SPMATRIX_COO);
	spectral_matrix_sparse_init(s, lmin, m, g, mat);

	// calloc, not alloc: if a solver path ever fails to write la, a zero eigenvalue is a
	// visible wrong answer rather than whatever happened to be on the heap.
	gsl_vector* la_prev = gsl_vector_calloc(nmax);
	if( spectral_solver_eigenvalues_n(s, l, m, g, la_prev, mat) ){
		gsl_vector_free(la_prev); gsl_spmatrix_free(mat);
		std::cout << "(SWSH) Eigenvalue solve failed for (s, l, m, gamma) = ("<<s<<", "<<l<<", "<<m<<", "<<g<<"); returning the spherical limit.\n";
		return l*(l + 1) - s*(s + 1);
	}
	double val_prev = gsl_vector_get(la_prev, l - lmin);

	nmax += 2;
	spectral_matrix_sparse(s, lmin, m, g, mat, nmax); // only 2 new rows added
	gsl_vector* la_cur = gsl_vector_calloc(nmax);
	if( spectral_solver_eigenvalues_n(s, l, m, g, la_cur, mat) ){
		// keep the coarser but valid result rather than reading an unwritten vector
		double coarse = val_prev;
		gsl_vector_free(la_cur); gsl_vector_free(la_prev); gsl_spmatrix_free(mat);
		return coarse;
	}
	double val_cur = gsl_vector_get(la_cur, l - lmin);

	double relerror = std::abs(1. - val_prev/val_cur);

	while(0.01*relerror > DBL_EPSILON && relerror > 0. && nmax < SPECTRAL_NMAX){
		gsl_vector_free(la_prev);
		la_prev = la_cur;
		val_prev = val_cur;

		nmax += 5;
		spectral_matrix_sparse(s, lmin, m, g, mat, nmax); // only 5 new rows
		la_cur = gsl_vector_calloc(nmax);
		if( spectral_solver_eigenvalues_n(s, l, m, g, la_cur, mat) ){
			// stop refining and keep the last converged value; la_cur is zeroed, not stale
			val_cur = val_prev;
			break;
		}
		val_cur = gsl_vector_get(la_cur, l - lmin);
		relerror = std::abs(1. - val_prev/val_cur);
	}

	if(nmax >= SPECTRAL_NMAX){
		std::cout << "(SWSH) Max number of iterations executed for swsh_eigenvalue. \n";
	}

	gsl_spmatrix_free(mat);
	gsl_vector_free(la_prev);
	double result = val_cur;
	gsl_vector_free(la_cur);

	return result;
}

///////////////////////////////////////
// Spin-Weighted Spherical Harmonics //
///////////////////////////////////////

Complex Yslm(const int &s, const int &l, const int &m, const double &th, const double &ph){
	return Yslm(s, l, m, th)*exp(I*Complex(m)*ph);
}

double Yslm(const int &s, const int &l, const int &m, const double &th){
	// Jacobi form for every s, accurate to the poles (see jacobi_parameters). The old sum of
	// scalar Ylm over sin^|s| lost up to all of its digits there, and even the scalar Ylm
	// (s = 0) was off by ~1e-4 relative within 1e-6 of a pole.
	return yslm_jacobi(s, l, m, th);
}

double Yslm_derivative(const int &s, const int &l, const int &m, const double &th){
	return yslm_derivative_jacobi(s, l, m, th);
}

double Yslm_derivative2(const int &s, const int &l, const int &m, const double &th){
	if( s == 0 ) return Ylm_derivative2(l, m, th);
	if( std::abs(th) < 1.e-14 ){
		if (m == -s + 2) {
			double phase = 0.25 * pow(-1., -s + 2);
			double prefactor = std::sqrt(l + s) * std::sqrt(l - s + 1.) * std::sqrt(l + s - 1.) * std::sqrt(l - s + 2.);
			return phase * prefactor * std::sqrt((2.*l + 1.)/(4.*M_PI));
		}
		else if (m == -s){
			double phase = 0.25 * pow(-1., -s + 1);
			double prefactor = std::sqrt(l + s) * std::sqrt(l - s + 1.) * std::sqrt(l + s) * std::sqrt(l - s + 1.);
			prefactor += std::sqrt(l - s) * std::sqrt(l + s + 1.) * std::sqrt(l - s) * std::sqrt(l + s + 1.);
			return phase * prefactor * std::sqrt((2.*l + 1.)/(4.*M_PI));
		}
		else if (m == -s - 2) {
			double phase = 0.25 * pow(-1., -s - 2);
			double prefactor = std::sqrt(l - s) * std::sqrt(l + s + 1.) * std::sqrt(l - s - 1.) * std::sqrt(l + s + 2.);
			return phase * prefactor * std::sqrt((2.*l + 1.)/(4.*M_PI));
		}
		else return 0.;
	}
	if( std::abs(th - M_PI) < 1.e-14 ){
		if (m == s - 2) {
			double phase = 0.25 * pow(-1., l + 2);
			double prefactor = std::sqrt(l + s) * std::sqrt(l - s + 1.) * std::sqrt(l + s - 1.) * std::sqrt(l - s + 2.);
			return phase * prefactor * std::sqrt((2.*l + 1.)/(4.*M_PI));
		}
		else if (m == s){
			double phase = 0.25 * pow(-1., l + 1);
			double prefactor = std::sqrt(l + s) * std::sqrt(l - s + 1.) * std::sqrt(l + s) * std::sqrt(l - s + 1.);
			prefactor += std::sqrt(l - s) * std::sqrt(l + s + 1.) * std::sqrt(l - s) * std::sqrt(l + s + 1.);
			return phase * prefactor * std::sqrt((2.*l + 1.)/(4.*M_PI));
		}
		else if (m == s + 2) {
			double phase = 0.25 * pow(-1., l - 2);
			double prefactor = std::sqrt(l - s) * std::sqrt(l + s + 1.) * std::sqrt(l - s - 1.) * std::sqrt(l + s + 2.);
			return phase * prefactor * std::sqrt((2.*l + 1.)/(4.*M_PI));
		}
		else return 0.;
	}
	double yslm = Yslm(s, l, m, th);
	double DyslmDtheta = Yslm_derivative(s, l, m, th);
	double eigen = double(l*(l + 1) - s*(s + 1));
	return (-s - eigen + (s / std::tan(th) + m / std::sin(th))*(s / std::tan(th) + m / std::sin(th)))*yslm - DyslmDtheta / std::tan(th);
}

// double Yslm_derivative(const int &s, const int &l, const int &m, const double &th){
// 	if( s == 0 ) return Ylm_derivative(l, m, th);
//
// 	int lmin = std::max(std::abs(m), l - std::abs(s));
// 	int lmax = l + std::abs(s);
// 	double DyslmDtheta = 0;
//
// 	for(int i = lmin; i <= lmax; i++){
// 		DyslmDtheta += Asljm(s, l, i, m)*(Ylm_derivative(i, m, th)*pow(sin(th), -std::abs(s)) - std::abs(s)*cos(th)*Ylm(i, m, th)*pow(sin(th), -std::abs(s)-1));
// 	}
// 	return DyslmDtheta*pow(sin(th), -std::abs(s));
// }

////////////////////
// Test functions //
////////////////////

void test_swsh_eigenvalue(void){
	int s = 2, l = 8, m = 2;
	double g = -4.8, th = 0.1, ph = 0.5;
	int start, stop;
	double duration;
	double la;

	start = clock();
	std::cout << "s = " << s << ", l = " << l << ", m = " << m << ", g = " << g << "\n";
	la = swsh_eigenvalue(s, l, m, g);
	std::cout << "la = " << la << "\n";
	std::cout << Sslm(s, l, m, g, th, ph) << "\n";
	stop = clock();
	duration = (stop-start)/double(CLOCKS_PER_SEC);
	std::cout << "test time 1: " << duration << " s" << std::endl;
};

void test_swsh_class(){
	int s = -2, l = 8, m = 1;
	double g = 4.8;
	Vector th = {0.1, 0.2};
	int start, stop;
	double duration;
	double la;

	start = clock();
	std::cout << "s = " << s << ", l = " << l << ", m = " << m << ", g = " << g << "\n";
	SpinWeightedHarmonic Slm(s, l, m, g, th);
	Slm.generateSolutions();
	Slm.generateDerivatives();
	la = Slm.getEigenvalue();
	std::cout << "la = " << la << "\n";
	std::cout << Slm.getSolution(0) << "\n";
	std::cout << Slm.getDerivative(0) << "\n";
	std::cout << Slm.getSecondDerivative(0) << "\n";
	stop = clock();
	duration = (stop-start)/double(CLOCKS_PER_SEC);
	std::cout << "test time 1: " << duration << " s" << std::endl;
};

/////////////////////////////////////////////////////////
// Grid-accelerated spheroidal-harmonic evaluation      //
/////////////////////////////////////////////////////////

// The spheroidal harmonic is S = sum_i b_i sY_{l_i m}, with every l_i sharing the same
// Jacobi parameters (a, b): one recurrence in n per theta gives every coupling term, and
// the derivative needs one more recurrence at (a + 1, b + 1). The old grid built sY from
// scalar Ylm divided by sin^|s|, which lost all precision near the poles (see
// jacobi_parameters).
void SpinWeightedHarmonic::fillYlmGrid(){
	if(_gridFilled){ return; }
	int lmin = getMinCouplingModeNumber();
	int imax = int(_bcoupling.size()) - 1;
	while(imax > 0 && _bcoupling[imax] == 0.){ imax--; }
	_imaxCoupling = imax;

	int a, b, lam;
	jacobi_parameters(a, b, lam, _s, _m);
	int nterm = imax + 1;
	// theta-independent weights: c_i = b_i * norm(l_i), and the P' factor (n+a+b+1)/2
	Vector coef(nterm, 0.), dcoef(nterm, 0.);
	for(int i = 0; i < nterm; i++){
		if(_bcoupling[i] == 0.){ continue; }
		coef[i] = _bcoupling[i]*jacobi_norm(a, b, lam, _s, lmin + i, _m);
		dcoef[i] = 0.5*double(i + a + b + 1)*coef[i];
	}

	int nth = int(_theta.size());
	_Sgrid.assign(nth, 0.);
	_SPgrid.assign(nth, 0.);
	Vector P(nterm, 0.), Q(nterm > 1 ? nterm - 1 : 1, 0.);
	for(int ith = 0; ith < nth; ith++){
		if(ith > 0 && _theta[ith] == _theta[ith - 1]){
			// repeated grid point (e.g. equatorial orbits): identical value
			_Sgrid[ith] = _Sgrid[ith - 1];
			_SPgrid[ith] = _SPgrid[ith - 1];
			continue;
		}
		double th = _theta[ith];
		double sh = sin(0.5*th), ch = cos(0.5*th), x = cos(th);
		jacobi_recurrence(P, nterm, a, b, x);
		double sumP = 0.;
		for(int i = 0; i < nterm; i++){ sumP += coef[i]*P[i]; }
		double sumQ = 0.;
		if(nterm > 1){
			jacobi_recurrence(Q, nterm - 1, a + 1, b + 1, x);
			for(int i = 1; i < nterm; i++){ sumQ += dcoef[i]*Q[i - 1]; }
		}
		_Sgrid[ith] = jacobi_prefactor(sh, ch, a, b)*sumP;
		_SPgrid[ith] = jacobi_prefactor_derivative(sh, ch, a, b)*sumP
			- 2.*jacobi_prefactor(sh, ch, a + 1, b + 1)*sumQ;
	}
	_gridFilled = true;
}

double SpinWeightedHarmonic::SslmGrid(int ith){
	return _Sgrid[ith];
}

double SpinWeightedHarmonic::SslmDerivativeGrid(int ith){
	return _SPgrid[ith];
}

