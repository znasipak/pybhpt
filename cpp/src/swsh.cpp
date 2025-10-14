// swsh.c

#include "swsh.hpp"

#define SPECTRAL_NMAX 300
#define SPECTRAL_NMAX_INIT_ADD 15
#define SPECTRAL_COUPLING_TEST_EPS 1.e-25
#define SPECTRAL_COUPLING_JMAX_EPS 1.e-25
#define SPECTRAL_COUPLING_CONVERGE_EPS 1.e-25
#define COUPLING_VECTOR_MAX 50
#define ZERO_FREQ_LIMIT 1.e-11

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
	for(size_t i = 0; i < _Slm.size(); i++){
		_Slm[i] = Sslm(_s, _L, _m, _gamma, _bcoupling, _theta[i]);
	}
	for(size_t i = 0; i < _SlmP.size(); i++){
		_SlmP[i] = Sslm_derivative(_s, _L, _m, _gamma, _bcoupling, _theta[i]);
	}

	return 0;
}

int SpinWeightedHarmonic::generateCouplingCoefficients(){
	if(_bcoupling[_L - getMinCouplingModeNumber()] == 0 && _bcoupling[_L - getMinCouplingModeNumber() + 1] == 0){
		spectral_solver(_s, _L, _m, _gamma, _lambda, _bcoupling);
	}
	return 0;
}

int SpinWeightedHarmonic::generateSolutions(){
	generateCouplingCoefficients();
	for(size_t i = 0; i < _Slm.size(); i++){
		_Slm[i] = Sslm(_s, _L, _m, _gamma, _bcoupling, _theta[i]);
	}

	return 0;
}

int SpinWeightedHarmonic::generateDerivatives(){
	generateCouplingCoefficients();
	for(size_t i = 0; i < _SlmP.size(); i++){
		_SlmP[i] = Sslm_derivative(_s, _L, _m, _gamma, _bcoupling, _theta[i]);
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
	if(nmax == SPECTRAL_NMAX){
		std::cout << "(SWSH) Max number of iterations executed for spectral solver. \n";
	}

	if(test == SUCCESS){ // if convergence test was successful, return data from the most resolved (highest nmax) spectral eigenvalue problem
		gsl_matrix_free(bmat);
		gsl_vector_free(la);
		bmat = bmat2;
		la = la2;
	}else if(test == FAIL){
		gsl_matrix_free(bmat); gsl_matrix_free(bmat2);
		gsl_vector_free(la); gsl_vector_free(la2);
		std::cout << "(SWSH) Error in computing coupling coefficients for (s, l, m, gamma) = ("<<s<<", "<<l<<", "<<m<<", "<<g<<") \n";
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

	double weight = spectral_weight(g);
	gsl_vector_scale(la, 1/weight);

	return 0;
}

int spectral_solver_n(const int &s, const int &l, const int &m, const double &g, gsl_vector* la, gsl_matrix* bmat, gsl_spmatrix* mat){
	int lmin = std::max(std::abs(s), std::abs(m));
	int nmax = la->size;
	if( nmax < l-lmin ){
		return 1; // error, need larger nmax tolerance
	}

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

	double weight = spectral_weight(g);
	gsl_vector_scale(la, 1/weight);

	return 0;
}

coupling_test spherical_spheroidal_coupling_convergence_test(const int &s, const int &l, const int &m, const double &g, gsl_matrix* bmat, gsl_matrix* bmat2, coupling_converge &b_data){
	int lmin = std::max(std::abs(s), std::abs(m));
	int dim = bmat->size1;
	int dim2 = bmat2->size2;
	gsl_vector* bcol = gsl_vector_alloc(dim);
	gsl_vector* bcol2 = gsl_vector_alloc(dim2);
	double relerror = b_data.test_err;
	int testIndex = b_data.testIndex;

	gsl_matrix_get_col(bcol, bmat, l - lmin);
	if( gsl_vector_get(bcol, l - lmin) < 0 ){
		gsl_vector_scale(bcol, -1.);
	}
	gsl_matrix_get_col(bcol2, bmat2, l - lmin);
	if( gsl_vector_get(bcol2, l - lmin) < 0 ){
		gsl_vector_scale(bcol2, -1.);
	}

	double norm = gsl_vector_max(bcol);
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

	while(b_data.testIndex < dim - 1 && std::abs(gsl_vector_get(bcol, b_data.testIndex)/norm) > SPECTRAL_COUPLING_TEST_EPS){
		if(s == 0){
			b_data.testIndex += 2;
		}else{
			b_data.testIndex++;
		}
	}
	if(b_data.testIndex == dim){
		if(s == 0){
			b_data.testIndex-= 2;
		}else{
			b_data.testIndex--;
		}
	}

	double jmax_denom = gsl_vector_get(bcol2, b_data.jmax);
	double test_denom = gsl_vector_get(bcol2, b_data.testIndex);
	if(jmax_denom == 0){
		jmax_denom = DBL_EPSILON*pow(g, 1);
	}
	if(test_denom == 0){
		test_denom = DBL_EPSILON*pow(g, 1);
	}

	double jmax_num = gsl_vector_get(bcol, b_data.jmax);
	double test_num = gsl_vector_get(bcol, b_data.testIndex);
	if(jmax_num == 0){
		jmax_num = DBL_EPSILON*pow(g, 1);
	}
	if(test_num == 0){
		test_num = DBL_EPSILON*pow(g, 1);
	}


	b_data.jmax_err = std::abs(1 - jmax_num/jmax_denom);
	b_data.test_err = std::abs(1 - test_num/test_denom);

	double convergenceCriteria = SPECTRAL_COUPLING_CONVERGE_EPS;
	if(g > 1.){
		convergenceCriteria *= pow(g, 2);
	}

	if( b_data.test_err < convergenceCriteria ){
		return SUCCESS;
	}else if( relerror < b_data.test_err && testIndex == b_data.testIndex ){
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
	if( std::abs(l-j) > std::abs(s) ) return 0.;
	if( j < std::abs(m) ) return 0.;
	int lmin = std::abs(m) < std::abs(s) ? std::abs(s) : std::abs(m);
	if( l < lmin ) return 0;
	double aslmg = pow(-1., m + s*(1 + sgn<int>(s))/2);
	aslmg *= sqrt(pow(4, std::abs(s))*pow(factorial(std::abs(s)), 2)*(2*j + 1)*(2*l + 1)/factorial(std::abs(2*s)));
	aslmg *= w3j(std::abs(s), l, j, 0, m, -m);
	aslmg *= w3j(std::abs(s), l, j, s, -s, 0);

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
	// for(int i = 0; i < bvec.size(); i++){
	// 	swsh += bvec[i]*Yslm(s, lmin + i, m, th);
	// }
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
	// for(int i = 0; i < bvec.size(); i++){
	// 	swsh += bvec[i]*Yslm_derivative(s, lmin + i, m, th);
	// }

	return swsh;
}

double Sslm_secondDerivative(const int &s, const int &, const int &m, const double &g, const double& lambda, const double &th, const double &Slm, const double &SlmP){
	return -cos(th)*SlmP/sin(th) + (pow(g*sin(th), 2) + pow((m + s*cos(th))/sin(th), 2) + 2.*s*g*cos(th) - s - 2.*m*g - lambda)*Slm;
}

// SWSH Eigenvalues
double swsh_eigenvalue(const int &s, const int &l, const int &m, const double &g){
	unsigned int lmin = std::max(std::abs(s), std::abs(m));
	unsigned int nmax = l - lmin + SPECTRAL_NMAX_INIT_ADD;

	double swshtest = spectral_solver(s, l, m, g, nmax);
	nmax += 2;
	double swshtest2 = spectral_solver(s, l, m, g, nmax);
	double relerror = std::abs(1 - swshtest/swshtest2);

	while(0.01*relerror > DBL_EPSILON && relerror > 0. && nmax < SPECTRAL_NMAX){
		swshtest = swshtest2;
		nmax += 5;
		swshtest2 = spectral_solver(s, l, m, g, nmax);
		relerror = std::abs(1-swshtest/swshtest2);
	}

	return swshtest2;
}

///////////////////////////////////////
// Spin-Weighted Spherical Harmonics //
///////////////////////////////////////

Complex Yslm(const int &s, const int &l, const int &m, const double &th, const double &ph){
	return Yslm(s, l, m, th)*exp(I*Complex(m)*ph);
}

double Yslm(const int &s, const int &l, const int &m, const double &th){
	if( s == 0 ) return Ylm(l, m, th);
	if( th == 0. ){
		if (m != -s) return 0.;
		else return pow(-1, s)*Yslm(0, l, 0, 0.);
	}
	int lmin = std::max(std::abs(m), l - std::abs(s));
	int lmax = l + std::abs(s);
	double yslm = 0;

	for(int i = lmin; i <= lmax; i++){
		yslm += Asljm(s, l, i, m)*Ylm(i, m, th);
	}

	return yslm/pow(sin(th), std::abs(s));
}

double Yslm_derivative(const int &s, const int &l, const int &m, const double &th){
	if( s == 0 ) return Ylm_derivative(l, m, th);

	int lmin = std::max(std::abs(m), l - std::abs(s) - 1);
	int lmax = l + std::abs(s) + 1;
	double DyslmDtheta = 0;

	for(int i = lmin; i <= lmax; i++){
		DyslmDtheta += dAsljm(s, l, i, m)*Ylm(i, m, th);
	}
	return -DyslmDtheta*pow(sin(th), -1 - std::abs(s));
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
