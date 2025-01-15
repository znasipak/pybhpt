// src bessel.cpp

#include "bessel.hpp"

#define BESSEL_EPSILON 1.e-15
#define NMAX 300
#define BESSEL_MAX 300


double bessel_J(int n, double x){
	if( n == 0 ){
		return gsl_sf_bessel_J0(x);
	}else if( n == 1 ){
		return gsl_sf_bessel_J1(x);
	}else if( n < 0 ){
		return pow(-1, n)*bessel_J(-n, x);
	}
	return gsl_sf_bessel_Jn(n, x);
}

double bessel_I(int n, double x){
	if( n == 0 ){
		return gsl_sf_bessel_I0(x);
	}else if( n == 1 ){
		return gsl_sf_bessel_I1(x);
	}else if( n < 0 ){
		return bessel_I(-n, x);
	}
	return gsl_sf_bessel_In(n, x);
}

// transform to real argument using DLMF (10.27.6)
Complex bessel_I(Complex nu, Complex z){
	if( std::imag(z) == 0 ){
		return bessel_I_bessel_series(nu, std::real(z));
	}else if( std::real(z) == 0){
		return pow(I, -nu)*bessel_J_bessel_series(nu, std::real(I*z));
	}else{
		std::cout << "BESSEL: ERROR: No algorithm implemented to calculate Bessel functions for arguments with real and imaginary parts. \n";
		return 0.;
	}
}

Complex bessel_I(int n, Complex z){
	if( std::imag(z) == 0 ){
		return bessel_I(n, std::real(z));
	}else if( std::real(z) == 0){
		return pow(I, -n)*bessel_J(n, std::real(I*z));
	}else{
		std::cout << "BESSEL: ERROR: No algorithm implemented to calculate Bessel functions for arguments with real and imaginary parts. \n";
		return 0.;
	}
}

// transform to real argument using DLMF (10.27.6)
Complex bessel_J(Complex nu, Complex z){
	if( std::imag(z) == 0 ){
		return bessel_J_bessel_series(nu, std::real(z));
	}else if( std::real(z) == 0 && std::imag(z) < 0 ){
		return pow(I, -nu)*bessel_I_bessel_series(nu, std::real(I*z));
	}else if( std::real(z) == 0 && std::imag(z) > 0 ){
			return pow(I, nu)*bessel_I_bessel_series(nu, std::real(-I*z));
	}else{
		std::cout << "BESSEL: ERROR: No algorithm implemented to calculate Bessel functions for arguments with real and imaginary parts. \n";
		return 0.;
	}
}

// these series expansions come from Eq (1.24) in Luke (1959) Expansion of Confluent Hypergeometric Functions and Eq (1) of Sec 5.21 in Watson (1922) A Treatise on the Theory of Bessel Functions
Complex bessel_I_bessel_series(Complex nu, double x){
	Complex term = 0., sum = 0.;
	Complex sumFactor = 1.;
	double error;
	term = sumFactor*bessel_I(0, x)/cgamma(nu + 1.);
	sum += term;
	error = std::abs(term/sum);
	
	int n = 1;
	while(error > BESSEL_EPSILON && n < NMAX){
		sumFactor *= -1.;
		term = 2.*sumFactor*bessel_I(2*n, x)*exp(lgamma(nu + 1.) - lgamma(nu + Complex(n) + 1.) - lgamma(nu - Complex(n) + 1.));
		error = std::abs(term/sum);
		sum += term;
		n += 1;
	}
	
	return pow(0.5*x, nu)*sum;
}

Complex bessel_J_bessel_series(Complex nu, double x){
	Complex term = 0., sum = 0.;
	double error;
	term = bessel_J(0, x)/cgamma(nu + 1.);
	sum += term;
	error = std::abs(term/sum);
	
	int n = 1;
	while(error > BESSEL_EPSILON && n < NMAX){
		term = 2.*bessel_J(2*n, x)*exp(lgamma(nu + 1.) - lgamma(nu + Complex(n) + 1.) - lgamma(nu - Complex(n) + 1.));
		error = std::abs(term/sum);
		sum += term;
		n += 1;
	}
	
	return pow(0.5*x, nu)*sum;
}

// DLMF 13.11.1
Complex hypergeo_1f1_bessel_series(Complex a, Complex b, Complex z){
	if( std::real(z) != 0. && std::imag(z) != 0. ){
		std::cout << "BESSEL: ERROR: No algorithm implemented to calculate the regular confluent hypergeometric function as a series of Bessel functions for arguments with real and imaginary parts. \n";
		return 0.;
	}
	
	Complex term = 0., sum = 0.;
	double error;
	int n = 0;
	term = (a - 0.5)*bessel_I(a - 0.5, 0.5*z)*exp(lgamma(2.*a - 1.) + lgamma(2.*a - b) - lgamma(b));
	sum += term;
	error = std::abs(term/sum);
	n += 1;
	
	while(error > BESSEL_EPSILON && n < NMAX){
		term = (a - 0.5 + Complex(n))*bessel_I(a - 0.5 + Complex(n), 0.5*z)*exp(lgamma(2.*a - 1. + Complex(n)) + lgamma(2.*a - b + Complex(n)) - lgamma(b + Complex(n)) - lgamma(1. + Complex(n)));
		error = std::abs(term/sum);
		sum += term;
		n += 1;
	}
		
	return exp(0.5*z + lgamma(a - 0.5) - lgamma(2.*a - 1.) - lgamma(2.*a - b) + lgamma(b))*pow(0.25*z, 0.5 - a)*sum;
}

// DLMF 13.11.2
Result hypergeo_1f1_bessel_series_2(Complex a, Complex b, Complex z){
	if( std::real(z) != 0. && std::imag(z) != 0. ){
		std::cout << "BESSEL: ERROR: No algorithm implemented to calculate the regular confluent hypergeometric function as a series of Bessel functions for arguments with real and imaginary parts. \n";
		return Result(0., 0.);
	}
	
	Complex term = 0., sum = 0.;
	double error, maxAbsTerm = 0.;
	int n = 0;
	term = (b - a - 0.5)*bessel_I(b - a - 0.5, 0.5*z)*exp(lgamma(2.*(b - a) - 1.) + lgamma(b - 2.*a) - lgamma(b));
	maxAbsTerm = std::abs(term);
	sum += term;
	error = std::abs(term/sum);
	n += 1;
	
	while(error > BESSEL_EPSILON && n < NMAX){
		term = pow(-1, n)*(b - a - 0.5 + Complex(n))*bessel_I(b - a - 0.5 + Complex(n), 0.5*z)*exp(lgamma(2.*(b - a) - 1. + Complex(n)) + lgamma(b - 2.*a + Complex(n)) - lgamma(b + Complex(n)) - lgamma(1. + Complex(n)));
		error = std::abs(term/sum);
		sum += term;
		n += 1;
		maxAbsTerm = std::abs(term) > maxAbsTerm ? std::abs(term) : maxAbsTerm;
	}
	Complex hypergeo1F1 = exp(0.5*z + lgamma(b - a - 0.5) - lgamma(2.*(b - a) - 1.) - lgamma(b - 2.*a) + lgamma(b))*pow(0.25*z, a - b + 0.5)*sum;
	error = error > std::abs(maxAbsTerm/sum*DBL_EPSILON) ? error : std::abs(maxAbsTerm/sum*DBL_EPSILON);

	return Result(hypergeo1F1, hypergeo1F1*error);
}

// Eq (1.10) Luke (1959) Expansion of Confluent Hypergeometric Functions
Complex Gk(int k, Complex a, Complex c){
	Result f21 = hypergeo_2F1(1. - 2.*a, 2.*(c - a), 1. + Complex(k) + c - 2.*a, 0.5);
	Complex prefactor = pow(-1., k)*cos(M_PI*(c - a))/pow(2., 2.*(c - a) - 1.)*exp(lgamma(Complex(k) - c + 1.) + lgamma(2.*c - 2.*a) - lgamma(1. + Complex(k) + c - 2.*a));
	return prefactor*f21.getValue();
}

Complex Rk(int k, Complex a, Complex c){
	if( a == c ){
		return pow(-1., k);
	}
	return exp(lgamma(c) - lgamma(a) - lgamma(c - a))*(Gk(k, a, c) + pow(-1., k)*Gk(k, c - a, c));
}

// Eq (1.8) Luke (1959) Expansion of Confluent Hypergeometric Functions
Complex hypergeo_1f1_bessel_series_3(Complex a, Complex b, Complex z){
	if( std::real(z) != 0. && std::imag(z) != 0. ){
		std::cout << "BESSEL: ERROR: No algorithm implemented to calculate the regular confluent hypergeometric function as a series of Bessel functions for arguments with real and imaginary parts. \n";
		return 0.;
	}
	
	Complex term = 0., sum = 0.;
	double error;
	int n = 0;
	term = Rk(n, a, b)*bessel_I(n, 0.5*z);
	sum += term;
	error = std::abs(term/sum);
	n += 1;
	
	while(error > BESSEL_EPSILON && n < NMAX){
		term = 2*pow(-1., n)*Rk(n, a, b)*bessel_I(n, 0.5*z);;
		error = std::abs(term/sum);
		sum += term;
		n += 1;
	}
		
	return exp(0.5*z)*sum;
}

Result hypergeo_u_bessel_series(Complex a, Complex b, Complex z){
	// calculate U(a, b, z) using its relation to the Bessel series expansions of the regular confluent hypergeometric functions M(a, b, z)
	// note that this only works for z is purely real or imaginary
	Result logMSeries1 = log(hypergeo_1f1_bessel_series_2(a, b, z));
	Result logMSeries2 = log(hypergeo_1f1_bessel_series_2(a - b + 1., 2. - b, z));
	Result term1 = exp(lgamma(1. - b) - lgamma(a - b + 1.) + logMSeries1);
	Result term2 = exp(lgamma(b - 1.) - lgamma(a) + (1. - b)*log(z) + logMSeries2);
	
	
	return term1 + term2;
}

// DLMF 13.11.2 + 13.2.42
Result hypergeo_u_bessel_series_2(Complex a, Complex b, Complex z){
	if( std::real(z) != 0. && std::imag(z) != 0. ){
		std::cout << "BESSEL: ERROR: No algorithm implemented to calculate the regular confluent hypergeometric function as a series of Bessel functions for arguments with real and imaginary parts. \n";
		return Result(0., 0.);
	}
	
	Complex term = 0., hypergeo1F1 = 0.;
	double error, maxAbsTerm = 0.;
	int n = 0;
	term = exp(lgamma(1. - b) - lgamma(a - b + 1.) + log(hypergeo_u_bessel_series_2_prefactor(a, b, z)) + log(hypergeo_u_bessel_series_2_term(n, a, b, z)));
	term += pow(z, 1. - b)*exp(lgamma(b - 1.) - lgamma(a) + log(hypergeo_u_bessel_series_2_prefactor(a - b + 1., 2. - b, z)) + log(hypergeo_u_bessel_series_2_term(n, a - b + 1., 2. - b, z)));
	maxAbsTerm = std::abs(term);
	hypergeo1F1 += term;
	error = std::abs(term/hypergeo1F1);
	n += 1;
	
	while(error > BESSEL_EPSILON && n < BESSEL_MAX){
		term = exp(lgamma(1. - b) - lgamma(a - b + 1.) + log(hypergeo_u_bessel_series_2_prefactor(a, b, z)) + log(hypergeo_u_bessel_series_2_term(n, a, b, z)));
		term += pow(z, 1. - b)*exp(lgamma(b - 1.) - lgamma(a) + log(hypergeo_u_bessel_series_2_prefactor(a - b + 1., 2. - b, z)) + log(hypergeo_u_bessel_series_2_term(n, a - b + 1., 2. - b, z)));
		error = std::abs(term/hypergeo1F1);
		hypergeo1F1 += term;
		n += 1;
		maxAbsTerm = std::abs(term) > maxAbsTerm ? std::abs(term) : maxAbsTerm;
	}
	error = error > std::abs(maxAbsTerm/hypergeo1F1*DBL_EPSILON) ? error : std::abs(maxAbsTerm/hypergeo1F1*DBL_EPSILON);

	return Result(hypergeo1F1, hypergeo1F1*error);
}

Complex hypergeo_u_bessel_series_2_term(int n, Complex a, Complex b, Complex z){
	return pow(-1, n)*(b - a - 0.5 + Complex(n))*bessel_I(b - a - 0.5 + Complex(n), 0.5*z)*exp(lgamma(2.*(b - a) - 1. + Complex(n)) + lgamma(b - 2.*a + Complex(n)) - lgamma(b + Complex(n)) - lgamma(1. + Complex(n)));

}

Complex hypergeo_u_bessel_series_2_prefactor(Complex a, Complex b, Complex z){
	return pow(0.25*z, a - b + 0.5)*exp(0.5*z + lgamma(b - a - 0.5) - lgamma(2.*(b - a) - 1.) - lgamma(b - 2.*a) + lgamma(b));
}