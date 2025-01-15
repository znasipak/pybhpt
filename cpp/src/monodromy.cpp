// monodromy.cpp

#include "monodromy.hpp"

void test_monodromy(){
	int s = 0, l = 62, m = 0;
	double q = 0., eps = 2*0.066, la = 62.*63.;
	Complex nu = nu_solver_monodromy(s, l, m, q, eps, la);

	std::cout << "nu = " << nu << "\n";
}

Complex nu_solver_monodromy(int s, int l, int m, double q, double eps, double la){
	double f, error, errorCheck, nuError;
	Complex nu, stokes, stokesTest;

	double kappa = sqrt(1. - pow(q,2));
	double tau = (eps - m*q)/kappa;
	int n, nmax;

	Complex ggCH, dCH, eCH, aCH, qCH;
	Complex m1C, m2C, mu, cosmu;

	aCH = 2.*I*eps*kappa*(1. - s + I*eps - I*tau);
	ggCH = 1. - s - I*(eps + tau);
	dCH = 1. + s + I*(eps - tau);
	eCH = 2.*I*eps*kappa;
	qCH = -(-s*(1. + s) + pow(eps, 2) + I*(2.*s - 1.)*eps*kappa - la - tau*(I + tau));

	CH_parameters params = {
		.aCH = aCH,
		.ggCH = ggCH,
		.dCH = dCH,
		.eCH = eCH,
		.qCH = qCH
	};

	m1C = aCH/eCH - (ggCH + dCH);
	m2C = -(aCH/eCH);
	mu = (m2C - m1C)/2.;
	cosmu = cos(2.*M_PI*mu);

	n = 50;
	nmax = 2000;

	// std::cout << "MONODROMY: Using monodromy methods for (s,l,m,om) = ("<<s<<","<<l<<","<<m<<","<<eps/2.<<"). \n";

	series_coeff a1 = series_coeff_init(n);
	series_coeff a2 = series_coeff_init(n);

	stokesTest = 0;
	stokes = monodromy_eigenvalue(a1, a2, params);
	error = std::abs(1. - std::real(stokesTest/stokes));
	errorCheck = std::abs(std::imag(stokes)/std::real(stokes));
	errorCheck = errorCheck > 0.1 ? 0.1 : errorCheck;

	while(error > errorCheck && errorCheck > DBL_EPSILON && n < nmax){
		stokesTest = stokes;
		n += 50;
		a1.nmax = n;
		a2.nmax = n;
		stokes = monodromy_eigenvalue(a1, a2, params);
		error = std::abs(1. - std::real(stokesTest/stokes));
		errorCheck = std::abs(std::imag(stokes)/std::real(stokes));
		errorCheck = errorCheck > 0.1 ? 0.1 : errorCheck;
	}

	if(isnan(std::real(stokes))){
		n -= 50;
		a1.nmax = n;
		a2.nmax = n;
		stokes = monodromy_eigenvalue(a1, a2, params);
	}

	f = std::real(cosmu + stokes);
	errorCheck = std::abs(std::imag(stokes)/f);
	nuError = std::abs(std::real(cosmu)/f)*DBL_EPSILON;
	errorCheck = (errorCheck < nuError)? nuError : errorCheck;
	// if(errorCheck > 1.e-3){
	// 	std::cout << "MONODROMY: Monodromy eigenvalue at infinity could not be determined due to precision loss of " << log10(errorCheck/DBL_EPSILON) << " digits \n";
	// 	return 0;
	// }

	if(f < -1.){
		nu = - 0.5 - std::abs(std::imag(cacos(f)/(2.0*M_PI)))*I;
	}else if(f < 1. && f >= -1.){
			nu = Complex(l) - acos(f)/(2.0*M_PI);
	}else{
		nu = - std::abs(std::imag(cacos(f)/(2.0*M_PI)))*I;
	}

	return nu;
}

Complex monodromy_eigenvalue(series_coeff &a1, series_coeff &a2, const CH_parameters &params){
	Complex m1C, m2C, g1, g2;
	Complex a1sum, a1sumN, a2sum, a2sumN, a12, a12recip, stokes;
	Complex aCH = params.aCH, ggCH = params.ggCH, dCH = params.dCH, eCH = params.eCH;

	m1C = aCH/eCH - (ggCH + dCH);
	m2C = -(aCH/eCH);

	g1 = cgamma(m1C - m2C);
	g2 = cgamma(m2C - m1C);

	generate_weighted_a1(a1, params);
	generate_weighted_a2(a2, params);

	int ninit = 5;
	int n = a1.nmax;
	if( n != a2.nmax ){
		std::cout << "MONODROMY: ERROR: a1 and a2 coefficient lists have unequal lengths. a1.nmax = "<< a1.nmax <<", a2.nmax = "<<a2.nmax<<". \n";
		return 0;
	}
	a1sumN = sum_weighted_a1(a1, ninit);
	a1sum = sum_weighted_a1_drop(a1, ninit);
	a2sumN = sum_weighted_a2(a2, ninit);
	a2sum = sum_weighted_a2_drop(a2, ninit);

	a12 = (a1sum*a2sum) + (a1sumN*a2sumN) + (a1sum*a2sumN) + (a1sumN*a2sum);
	if( std::abs(a1.coeffs[n]) != 1. ){
		std::cout << "MONODROMY: ERROR: a1 coefficients not properly normalized\n";
	}
	if( std::abs(a2.coeffs[n]) != 1. ){
		std::cout << "MONODROMY: ERROR: a2 coefficients not properly normalized\n";
	}
	a12recip = a1.coeffs[n]*a2.coeffs[n]/a12;
	stokes = pow(2.*M_PI, 2)*pow(-1., n-1)*a12recip/(2.*g1*g2);

	return stokes;
}

series_coeff series_coeff_init(int nmax){
	series_coeff a = {.nmax = nmax, .ndata = 1};
	a.coeffs[0] = 1.;
	return a;
}

int generate_weighted_a1(series_coeff &a, const CH_parameters &params){
	Complex m1C, m2C, prefm2, prefm1, am2, am1, am0;
	Complex aCH = params.aCH, ggCH = params.ggCH, dCH = params.dCH, eCH = params.eCH, qCH = params.qCH;
	Complex cn;

	m1C = aCH/eCH - (ggCH + dCH);
	m2C = -(aCH/eCH);

	int n = a.ndata - 1;
	int nmax = a.nmax;
	if( n == -1 ){
		std::cout << "Error: series coefficients not properly initialized \n";
		return 1;
	}else if( n == 0 ){
		am1 = 0.;
	}else{
		am1 = a.coeffs[n - 1]/(m1C - m2C); // rescale entry because we are actually storing the weighted coefficients while these recursion relations are for the unweighted coefficients
	}
	am0 = a.coeffs[n];

	while( n < nmax ){
		n++;
		am2 = am1/am0;
		am1 = 1.;
		cn = Complex(n);

		prefm2 = (aCH - eCH*(cn - 2. + ggCH + dCH))*(aCH - eCH*(cn + dCH - 1.))/(eCH*cn);
		prefm1 = 1. - dCH + (2.*aCH)/eCH + eCH - ggCH - cn + (-pow(aCH/eCH, 2.)
			+ (-1. + dCH)*eCH + (aCH*(-1. + dCH - eCH + ggCH))/eCH + qCH)/cn;

		am0 = prefm2*am2 + prefm1*am1;

		a.coeffs[n] = am0;
		a.ndata = n + 1;
		weighted_renormalize(a, cn - 1. - m2C + m1C);
	}

	return 0;
}

int generate_weighted_a2(series_coeff &a, const CH_parameters &params){
	Complex m1C, m2C, prefm2, prefm1, am2, am1, am0;
	Complex aCH = params.aCH, ggCH = params.ggCH, dCH = params.dCH, eCH = params.eCH, qCH = params.qCH;
	Complex cn;

	m1C = aCH/eCH - (ggCH + dCH);
	m2C = -(aCH/eCH);

	int n = a.ndata - 1;
	int nmax = a.nmax;
	if( n == -1 ){
		std::cout << "Error: series coefficients not properly initialized \n";
		return 1;
	}else if( n == 0 ){
		am1 = 0.;
	}else{
		am1 = a.coeffs[n - 1]/(m2C - m1C); // rescale entry because we are actually storing the weighted coefficients while these recursion relations are for the unweighted coefficients
	}
	if( std::abs(a.coeffs[n]) != 1. ){
		std::cout << "Error: series coefficients not properly normalized\n";
	}
	am0 = a.coeffs[n];

	while( n < nmax ){
		n++;
		am2 = am1/am0; // rescale coefficients by nmax coefficient
		am1 = 1.;
		cn = Complex(n);

		prefm2 = -(aCH + eCH*(cn - 2.))*(aCH + eCH*(cn - ggCH - 1.))/(eCH*cn);
		prefm1 = ((aCH + eCH*(cn - 1.))*(aCH + eCH*(cn - dCH + eCH - ggCH)) - pow(eCH, 2)*qCH)/(eCH*eCH*cn);

		am0 = prefm2*am2 + prefm1*am1;

		a.coeffs[n] = am0;
		a.ndata = n + 1;
		weighted_renormalize(a, cn - 1. + m2C - m1C);
	}

	return 0;
}

int weighted_renormalize(series_coeff &a, Complex weight){
	int nmax = a.ndata - 1;
	Complex norm, an;

	norm = a.coeffs[nmax];

	for(int n = 0; n < nmax; n++){
		an = a.coeffs[n];
		a.coeffs[n] = an*(weight - Complex(n))/norm;
	}

	a.coeffs[nmax] = 1.;

	return 0;
}

Complex sum_weighted_a1(const series_coeff &a){
	Complex aSumN, aSumNeg, aSumPos;
	aSumN = aSumNeg = aSumPos = 0.;
	int nmax = ceil(a.nmax/2);

	for(int j = 0; j<= nmax; j++){
		aSumN = a.coeffs[j];
		if( std::real(aSumN) < 0 ){
			aSumNeg += aSumN;
		}else{
			aSumPos += aSumN;
		}
	}

	return aSumNeg + aSumPos;
};

Complex sum_weighted_a1(const series_coeff &a, const int &n){
	Complex aSumN, aSumNeg, aSumPos;
	aSumN = aSumNeg = aSumPos = 0.;
	int nmax = n;
	if( nmax > a.ndata ){
		std::cout << "Error: out of bound index\n";
		return 0;
	}

	for(int j = 0; j <= nmax; j++){
		aSumN = a.coeffs[j];
		if( std::real(aSumN) < 0 ){
			aSumNeg += aSumN;
		}else{
			aSumPos += aSumN;
		}
	}

	return aSumNeg + aSumPos;
};

Complex sum_weighted_a1_drop(const series_coeff &a, const int &n){
	Complex aSumN, aSumNeg, aSumPos;
	aSumN = aSumNeg = aSumPos = 0.;
	int nmax = std::max(ceil(a.nmax/2), double(n + 2));
	if( nmax > a.ndata ){
		std::cout << "Error: out of bound index\n";
		return 0;
	}

	for(int j = n + 1; j<= nmax; j++){
		aSumN = a.coeffs[j];
		if( std::real(aSumN) < 0 ){
			aSumNeg += aSumN;
		}else{
			aSumPos += aSumN;
		}
	}

	return aSumNeg + aSumPos;
};

Complex sum_weighted_a2(const series_coeff &a){
	Complex aSumN, aSumNeg, aSumPos;
	aSumN = aSumNeg = aSumPos = 0.;
	int nmax = ceil(a.nmax/2);

	for(int j = 0; j<= nmax; j++){
		aSumN = pow(-1., j)*a.coeffs[j];
		if( std::real(aSumN) < 0 ){
			aSumNeg += aSumN;
		}else{
			aSumPos += aSumN;
		}
	}

	return aSumNeg + aSumPos;
};

Complex sum_weighted_a2(const series_coeff &a, const int &n){
	Complex aSumN, aSumNeg, aSumPos;
	aSumN = aSumNeg = aSumPos = 0.;
	int nmax = n;
	if( nmax > a.ndata ){
		std::cout << "Error: out of bound index\n";
		return 0;
	}

	for(int j = 0; j <= nmax; j++){
		aSumN = pow(-1, j)*a.coeffs[j];
		if( std::real(aSumN) < 0 ){
			aSumNeg += aSumN;
		}else{
			aSumPos += aSumN;
		}
	}

	return aSumNeg + aSumPos;
};

Complex sum_weighted_a2_drop(const series_coeff &a, const int &n){
	Complex aSumN, aSumNeg, aSumPos;
	aSumN = aSumNeg = aSumPos = 0.;
	int nmax = std::max(ceil(a.nmax/2), double(n + 2));
	if( nmax > a.ndata ){
		std::cout << "Error: out of bound index\n";
		return 0;
	}

	for(int j = n + 1; j<= nmax; j++){
		aSumN = pow(-1, j)*a.coeffs[j];
		if( std::real(aSumN) < 0 ){
			aSumNeg += aSumN;
		}else{
			aSumPos += aSumN;
		}
	}

	return aSumNeg + aSumPos;
};
