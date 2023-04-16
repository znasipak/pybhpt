// src hypergeo_u.cpp

#include "hypergeo_u.hpp"

#define HYPERGEO_U_EPSILON 2*DBL_EPSILON
#define CF_ITER_MAX 300

//**********************************************************************************
// Special functions that can defined in terms of confluent hypergeometric functions
//**********************************************************************************

Result hyper_incomplete_gamma(Complex a, Complex z){ return exp(log_hyper_incomplete_gamma(a, z)); } // incomplete \gamma function
Result hyper_incomplete_Gamma(Complex a, Complex z){ return exp(log_hyper_incomplete_Gamma(a, z)); } // incomplete \Gamma function
Result hyper_exponential_integral(Complex a, Complex z){ return exp(log_hyper_exponential_integral(a, z)); } // generalized exponential integral
Result hyper_erf(Complex z){ return exp(log_hyper_erf(z)); } // error function
Result hyper_erfc(Complex z){ return exp(log_hyper_erfc(z)); } // complementary error function
Result hyper_bessel_J(Complex nu, Complex z){ return exp(log_hyper_bessel_J(nu, z)); } // Cylindrical Bessel function J
Result hyper_bessel_I(Complex nu, Complex z){ return exp(log_hyper_bessel_I(nu, z)); } // Cylindrical Bessel function I
Result hyper_bessel_K(Complex nu, Complex z){ return exp(log_hyper_bessel_K(nu, z)); } // Cylindrical Bessel function K
Result hyper_Ai(Complex z){ return exp(log_hyper_Ai(z)); } // Airy function

Result hyper_1f1(Complex a, Complex b, Complex z){ return exp(log_hyper_1f1(a, b, z)); } // Regular confluent hypergeometric function {}_1 F_1
Result hyper_u(Complex a, Complex b, Complex z){ return exp(log_hyper_u(a, b, z)); } // Irregular confluent hypergeometric function ~ {}_2 F_0
Result hyper_0f1(Complex b, Complex z){ return exp(log_hyper_0f1(b, z)); } // Hypergeometric function {}_0 F_1

Result dhyper_1f1(Complex a, Complex b, Complex z){ return exp(log_dhyper_1f1(a, b, z)); } // Derivative of regular confluent hypergeometric function
Result dhyper_u(Complex a, Complex b, Complex z){ return exp(log_dhyper_u(a, b, z)); } // Derivative of irregular confluent hypergeometric function

//////
// Log of special functions (helps avoid overflow when evaluating for large parameters or arguments)
//////

Result log_hyper_incomplete_gamma(Complex a, Complex z){
	Result test = log_hyper_1f1(a, a + 1., -z) - log(a) + a*log(z);
	if( abs(test.getPrecision()) > 1.e-10 ){
		Result test_compare = log(cgamma(a) - hyper_incomplete_Gamma(a, z));
		if( abs(test.getPrecision()) > abs(test_compare.getPrecision()) ){
			test = test_compare;
		}
	}

	return test;
}

Result log_hyper_incomplete_Gamma(Complex a, Complex z){
	return log_hyper_u(1. - a, 1. - a, z) - z;
}

Result log_hyper_exponential_integral(Complex a, Complex z){
	return log_hyper_u(a, a, z) - z - (1. - a)*log(z);
}

Result log_hyper_erf(Complex z){
	// for some reason my 1F1 function fails for small a, but medium z
	// so once we get to larger arguments, just use erf = 1 - erfc, since
	// erfc is constructed from U, which has an asymptotic expansion for large z
	if( abs(z) > 20. ){
		return log(Result(1., DBL_EPSILON) - hyper_erfc(z));
	}
	Result test = log_hyper_1f1(1., 1.5, pow(z, 2.))  + log(2.*z) - pow(z, 2) - 0.5*log(M_PI);
	if( abs(test.getPrecision()) > 1.e-10 ){
		Result test_compare = log(Result(1., DBL_EPSILON) - hyper_erfc(z));
		if( abs(test.getPrecision()) > abs(test_compare.getPrecision()) ){
			test = test_compare;
		}
	}

	return test;
}

Result log_hyper_erfc(Complex z){
	return log_hyper_u(0.5, 0.5, pow(z, 2)) - pow(z, 2) - 0.5*log(M_PI);
}

Result log_hyper_bessel_J(Complex nu, Complex z){
	return -lgamma(nu + 1.) + nu*log(0.5*z) + log_hyper_0f1(nu + 1., -0.25*pow(z, 2));
}

Result log_hyper_bessel_I(Complex nu, Complex z){
	return -lgamma(nu + 1.) + nu*log(0.5*z) + log_hyper_0f1(nu + 1., 0.25*pow(z, 2));
}

Result log_hyper_bessel_K(Complex nu, Complex z){
	return nu*log(2.*z) - z + 0.5*log(M_PI) + log_hyper_u(nu + 0.5, 2.*nu + 1., 2.*z);
}

Result log_hyper_Ai(Complex z){
	return -0.5*log(M_PI) + log(z) + 2.*log(2.)/3. - 5.*log(3.)/6. - 2.*pow(z, 1.5)/3. + log_hyper_u(5./6., 5./3., 4.*pow(z, 1.5)/3.);
}

Result log_hyper_1f1(Complex a, Complex b, Complex z){
	return log_hypergeo_m(a, b, z);
}

Result log_hyper_u(Complex a, Complex b, Complex z){
	return log_hypergeo_u(a, b, z);
}

Result log_hyper_0f1(Complex b, Complex z){
	return -2.*pow(z, 0.5) + log_hyper_1f1(b - 0.5, 2.*b - 1., 4.*pow(z, 0.5));
}

Result log_dhyper_1f1(Complex a, Complex b, Complex z){
	return log_dhypergeo_m(a, b, z);
}

Result log_dhyper_u(Complex a, Complex b, Complex z){
	return log_dhypergeo_u(a, b, z);
}

//**********************************************************************************
// Irregular and regular confluent hypergeometric functions
//**********************************************************************************

Result log_hypergeo_u(Complex a, Complex b, Complex z){
	// test to see if any of the special cases are satisfied
	if( a == b - 1. ){
		return Result(-a*log(z), DBL_EPSILON/pow(z, -a));
	}else if( abs(a) == 0. ){
		return Result(0., 0.);
	}
	// if Re(z) < 0, use transformation to evaluate hypergeometric functions with Re(z) > 0.
	// This is helpful because some recurrence relations break down for Re(z) < 0. However,
	// if the parameter a is small, we do not need recurrence relations and therefore do
	// not need to make this transformation.
	if( std::real(z) < 0 && abs(std::imag(z)) < 2. && abs(a) > 1.){
		return log_hypergeo_u_negative_z(a, b, z);
	}

	Result hyperU(0., 1.);
	Result hyperU_compare = hyperU;

	double arg_max = -log(DBL_EPSILON);
	if( abs(z) > 100. && abs(a) < 2. && abs(b) < 2.){
		// if the argument is large and the parameters are small, go straight to the asymptotic expansion
		hyperU = hypergeo_u_asymptotic_series(a, b, z);
	}else if( abs(a*z/b) < arg_max || abs(z) < arg_max ){
		// if the series terms appear to be small or the argument is small, try Taylor series around z = 0
		hyperU = log_hypergeo_u_connection(a, b, z);
	}

	// check to see if the initial attempts at calculating U returned with enough accuracy/precision. If not
	// try the use of recurrence relations to improve calculation
	if( abs(hyperU.getPrecision()) > 100*HYPERGEO_U_EPSILON ){
		Result hyperU_compare = log_hypergeo_u_recurrence(a, b, z);
		if( abs(hyperU.getPrecision()) > abs(hyperU_compare.getPrecision()) )
			hyperU = hyperU_compare;
	}

	// if the recurrence relations failed, but the direct Taylor series calculation was not performed yet,
	// try to compute U using the Taylor series as last attempt
	if(  abs(hyperU.getPrecision()) > 1.e-10 && abs(z) >= arg_max && abs(a*z/b) >= arg_max ){
		hyperU_compare = log_hypergeo_u_connection(a, b, z);
		if( abs(hyperU.getPrecision()) > abs(hyperU_compare.getPrecision()) ){
			hyperU = hyperU_compare;
		}
	}

	// if(isnan(abs(hyperU.getValue()))){
	// 	std::cout << "a = "<<a<<", b = "<<b<<", z = "<<z<<" \n";
	// }

	return hyperU;
}

Result log_hypergeo_u_negative_z(Complex a, Complex b, Complex z){
	Result term1 = log_hypergeo_m(1. - a, 2. - b, -z);
	term1 += lgamma(b - 1.) - lgamma(a) + log(sin(M_PI*(a - b))/sin(M_PI*a) - exp(-I*M_PI*b)) + (1. - b)*log(-z) + z;
	Result term2 = log_hypergeo_u(b - a, b, -z);
	term2 += lgamma(1. - a) - lgamma(a - b + 1.) + z;

	return log(exp(term1) + exp(term2));
}

Result hypergeo_u(Complex a, Complex b, Complex z){
	return exp(log_hypergeo_u(a, b, z));
}

Result log_dhypergeo_u(Complex a, Complex b, Complex z){
	return log_hypergeo_u(a + 1., b + 1., z) + log(-a);
};

Result dhypergeo_u(Complex a, Complex b, Complex z){
	return exp(log_dhypergeo_u(a, b, z));
}

Result log_hypergeo_m(Complex a, Complex b, Complex z){
	if( a == b ){
		return Result(z, DBL_EPSILON/exp(z));
	}else if( abs(a) == 0. ){
		return Result(0., 0.);
	}
	// if the argument and parameters are small in magnitude, calculate the
	// irregular confluent hypergeometric function via a series expansion
	double arg_max = -log(DBL_EPSILON);
	Result hyperM(0., 1.);
	Result hyperM_compare = hyperM;
	if( abs(z) < arg_max || abs(a*z/b) < arg_max ){
		hyperM = log_hypergeo_m_recurrence(a, b, z);
	}

	if( abs(hyperM.getPrecision()) > 1.e-10 ){
		hyperM_compare = log_hypergeo_m_connection(a, b, z);
		if( abs(hyperM.getPrecision()) > abs(hyperM_compare.getPrecision()) ){
			hyperM = hyperM_compare;
		}
	}

	if(  abs(hyperM.getPrecision()) > 1.e-10 && abs(z) >= arg_max && abs(a*z/b) >= arg_max ){
		hyperM_compare = log_hypergeo_m_recurrence(a, b, z);
		if( abs(hyperM.getPrecision()) > abs(hyperM_compare.getPrecision()) ){
			hyperM = hyperM_compare;
		}
	}

	return hyperM;
}

Result hypergeo_m(Complex a, Complex b, Complex z){
	return exp(log_hypergeo_m(a, b, z));
}

Result log_dhypergeo_m(Complex a, Complex b, Complex z){
	return log_hypergeo_m(a + 1., b + 1., z) + log(a) - log(b);
};

Result dhypergeo_m(Complex a, Complex b, Complex z){
	return exp(log_dhypergeo_m(a, b, z));
};

//****************************************************************************************
// Connection formulae, recurrence relations, continued fractions, & asymptotic expansions
// for calculating the confluent hypergeometric functions
//****************************************************************************************

Result log_hypergeo_m_connection(Complex a, Complex b, Complex z){
	Domain region = error_bound_region(a, b, z);
	if(region == 0){
		return Result(0., 1.);
	}
	// There is a branch cut ambiguity here. Can choose exp(+I*a*M_PI) & exp(+I*(a-b)*M_PI) or exp(-I*a*M_PI) & exp(-I*(a-b)*M_PI)
	Result term1 = -I*a*M_PI - lgamma(b - a) + log(hypergeo_u_asymptotic_series(a, b, z));
	Result term2 = -I*(a - b)*M_PI - lgamma(a) + z + log(hypergeo_u_asymptotic_series(b - a, b, -z));
	if(std::real(term1.getValue() - term2.getValue()) > -log(DBL_EPSILON)){
		return term1 + lgamma(b);
	}else if(std::real(term2.getValue() - term1.getValue()) > -log(DBL_EPSILON)){
		return term2 + lgamma(b);
	}
	if( abs( 1. + (exp(term1 - term2)).getValue() ) < sqrt(DBL_EPSILON) ){
		//std::cout << "HYPERGEO_U: ERROR: Catastrophic cancellation when using the Kummer connection formula to compute M. \n";
	}
	return lgamma(b) + log(exp(term1) + exp(term2));
}

Result log_hypergeo_u_connection(Complex a, Complex b, Complex z){
	Result term1 = lgamma(1. - b) - lgamma(a - b + 1.) + a*log(z) + log(hypergeo_m_taylor_series(a, b, z));
	Result term2 = lgamma(b - 1.) - lgamma(a) + (a - b + 1.)*log(z) + log(hypergeo_m_taylor_series(a - b + 1., 2. - b, z));
	if(std::real(term1.getValue() - term2.getValue()) > -log(DBL_EPSILON)){
		return term1 - a*log(z);
	}else if(std::real(term2.getValue() - term1.getValue()) > -log(DBL_EPSILON)){
		return term2 - a*log(z);
	}
	if( abs( 1. + (exp(term1 - term2)).getValue() ) < sqrt(DBL_EPSILON) ){
		//std::cout << "HYPERGEO_U: ERROR: Catastrophic cancellation when using the Kummer connection formula to compute U. \n";
	}
	//std::cout << "log result = " << log(term1 + term2) << ", result = " << term1 + term2 << "\n";
	return log(exp(term1) + exp(term2)) - a*log(z);
}

// when the reccurence relations for U break down, use M solutions and their recurrence relations
Result log_hypergeo_u_connection_recurrence(Complex a, Complex b, Complex z){
	Result term1 = exp(lgamma(1. - b) - lgamma(a - b + 1.) + a*log(z) + log_hypergeo_m_recurrence(a, b, z));
	Result term2 = exp(lgamma(b - 1.) - lgamma(a) + (a - b + 1.)*log(z) + log_hypergeo_m_recurrence(a - b + 1., 2. - b, z));
	if( abs( 1. + (term1/term2).getValue() ) < sqrt(DBL_EPSILON) ){
		//std::cout << "HYPERGEO_U: ERROR: Catastrophic cancellation when using the Kummer connection formula to compute U. \n";
	}
	return log(term1 + term2) - a*log(z);
}

Result log_hypergeo_u_series(Complex a, Complex b, Complex z){
	Result hyperU(0., 1.);
	Domain region = error_bound_region(a, b, z);
	double tol = 1.e-12;
	double arg_max = -log(DBL_EPSILON);

	if(abs(z) < arg_max || region == 0){
		hyperU = log_hypergeo_u_connection(a, b, z);
	}

	if(abs(hyperU.getPrecision()) > tol && region > 0){
		Result hyperU_compare = log(hypergeo_u_asymptotic_series(a, b, z));
		if(abs(hyperU.getPrecision()) > abs(hyperU_compare.getPrecision())){
			hyperU = hyperU_compare;
		}
	}

	if( abs(z) > arg_max && abs(hyperU.getPrecision()) == 1. ){
		Result hyperU_compare = log_hypergeo_u_connection(a, b, z);
		if(abs(hyperU.getPrecision()) > abs(hyperU_compare.getPrecision())){
			hyperU = hyperU_compare;
		}
	}

	return hyperU;
}

Result log_hypergeo_u_recurrence(Complex aTemp, Complex bTemp, Complex z){
	// these recurrence relations break down for z is real and negative
	if(std::imag(z) == 0. && std::real(z) < 0.){
		//std::cout << "HYPERGEO_U: Recurrence relations cannot not be applied to U when z is real and < 0.\n";
		return log_hypergeo_u_connection_recurrence(aTemp, bTemp, z);
	}

	Complex aTransform, bTransform, logPrefactor;
	Result logHyperU(0., 1.);
	Result recurTest = logHyperU;
	Result logHyperU2 = logHyperU;
	Result logHyperU3 = logHyperU;
	if(std::real(aTemp)/std::real(bTemp) < 0){
		//if their real parts of the confluent parameters have opposite signs
		// transform them so they have the same sign
		aTransform = aTemp - bTemp + 1.;
		bTransform = 2. - bTemp;
		logPrefactor = (1. - bTemp)*log(z);
	}else{
		// otherwise use the original user-specified values for the parameters
		aTransform = aTemp;
		bTransform = bTemp;
		logPrefactor = 0.;
	}
	double minArgMagnitudeA = 1.;
	double minArgMagnitudeB = 1.;

	Complex a = aTransform, b = bTransform;
	// minimize the magntiude of the parameter arguments then use continued fraction ratios to recursively find the solution value for larger parameter values
	if(std::real(aTransform) > 0){
		while(std::real(a) > minArgMagnitudeA){
			a -= 1.;
			b -= 1.;
		}

		logHyperU = log_hypergeo_u_series(a, b, z);
		//std::cout << "U("<<a<<", "<<b<<", "<<z<<") = " << exp(logHyperU) << "\n";

		// if first calculation is too inaccurate, minimize argument b as well
		if( abs( logHyperU.getPrecision() ) > 1.e-5 ){
			while(std::real(b) > minArgMagnitudeB){
				b -= 1.;
			}
			while(std::real(b) < -minArgMagnitudeB){
				b += 1.;
			}
			logHyperU = log_hypergeo_u_series(a, b, z);
			//std::cout << "U("<<a<<", "<<b<<", "<<z<<") = " << exp(logHyperU) << "\n";
		}

		// make use recurrence relations or CF = U(a, b, z)/U(a + 1, b + 1, z)
		if(std::real(a) < std::real(aTransform) && std::real(b) < std::real(bTransform)){
			logHyperU2 = log(hypergo_u_ab_over_u_ap1bp1(a - 1., b - 1., z));
			// if the continued fraction method appears to be unstable (i.e. CF ~ 1) then
			// explicitly calculate U(a + 1, b + 1, z)
			if(abs(logHyperU2.getValue()) < 4.8){
				logHyperU2 = log_hypergeo_u_series(a - 1., b - 1., z);
			}else{
				logHyperU2 = logHyperU - logHyperU2;
			}
			//std::cout << "U("<<a - 1.<<", "<<b - 1.<<", "<<z<<") = " << exp(logHyperU2) << "\n";
		}
		while(std::real(a) < std::real(aTransform) && std::real(b) < std::real(bTransform)){
			logHyperU3 = logHyperU2;
			logHyperU2 = logHyperU;
			recurTest = log(hypergo_u_ap1bp1(exp(logHyperU2), exp(logHyperU3), a, b, z));
			if( abs(std::real((recurTest - logHyperU2).getValue())) < 1.2 ){
				// if the functions are of similar magnitude then these recurrence relations and
				// continued fraction methods are most likely note very stable. Calculate functions
				// directly
				logHyperU = log_hypergeo_u_series(a + 1., b + 1., z);
			}else if( std::real(recurTest.getValue()) < std::real(logHyperU.getValue()) ){
				// if functions are decreasing in magnitude make use of the continued fractions
				//std::cout << "Using CF \n";
				logHyperU -= log(hypergo_u_ab_over_u_ap1bp1(a, b, z));
			}else{
				// otherwise use value from recurrence
				//std::cout << "Using recurrence \n";
				logHyperU = recurTest;
			}
			a += 1.;
			b += 1.;
			//std::cout << "U("<<a<<", "<<b<<", "<<z<<") = " << exp(logHyperU) << "\n";
		}

		// make use of CF = U(a, b, z)/U(a + 1, b, z)
		if(std::real(a) < std::real(aTransform)){
			logHyperU2 = logHyperU;
			logHyperU2 += log(hypergo_u_a_over_u_ap1(a - 1., b, z));
		}
		while(std::real(a) < std::real(aTransform)){
			logHyperU3 = logHyperU2;
			logHyperU2 = logHyperU;
			recurTest = log(hypergo_u_ap1(exp(logHyperU2), exp(logHyperU3), a, b, z));
			if( std::real(recurTest.getValue()) < std::real(logHyperU.getValue()) ){
				// if functions are decreasing in magnitude make use of the continued fractions
				logHyperU -= log(hypergo_u_a_over_u_ap1(a, b, z));
			}else{
				// otherwise use value from recurrence
				logHyperU = recurTest;
			}
			a += 1.;
			//std::cout << "U("<<a<<", "<<b<<", "<<z<<") = " << exp(logHyperU) << "\n";
		}

		// make use of CF = U(a, b + 1, z)/U(a, b, z)
		if(std::real(b) < std::real(bTransform)){
			logHyperU2 = logHyperU;
			logHyperU2 -= log(hypergo_u_bp1_over_u_b(a, b - 1., z));
			//std::cout << "U("<<a<<", "<<b - 1.<<", "<<z<<") = " << exp(logHyperU2) << "\n";
		}
		while(std::real(b) < std::real(bTransform)){
			logHyperU3 = logHyperU2;
			logHyperU2 = logHyperU;
			recurTest = log(hypergo_u_bp1(exp(logHyperU2), exp(logHyperU3), a, b, z));
			if( std::real(recurTest.getValue()) < std::real(logHyperU.getValue()) ){
				// if functions are decreasing in magnitude make use of the continued fractions
				logHyperU += log(hypergo_u_bp1_over_u_b(a, b, z));
			}else{
				// otherwise use value from recurrence
				logHyperU = recurTest;
			}
			b += 1.;
			//std::cout << "U("<<a<<", "<<b<<", "<<z<<") = " << exp(logHyperU) << "\n";
		}
	}else{
		// first minimize the argument a
		while(std::real(a) < -minArgMagnitudeA){
			a += 1.;
			b += 1.;
		}

		logHyperU = log_hypergeo_u_series(a, b, z);
		//std::cout << "U("<<a<<", "<<b<<", "<<z<<") = " << exp(logHyperU) << "\n";

		// if first calculation is too inaccurate, minimize argument b as well
		if( abs( logHyperU.getPrecision() ) > 1.e-5 ){
			while(std::real(b) > minArgMagnitudeB){
				b -= 1.;
			}
			while(std::real(b) < -minArgMagnitudeB){
				b += 1.;
			}
			logHyperU = log_hypergeo_u_series(a, b, z);
			//std::cout << "U("<<a<<", "<<b<<", "<<z<<") = " << exp(logHyperU) << "\n";
		}

		// make use recurrence relations or CF = U(a, b, z)/U(a + 1, b + 1, z)
		if(std::real(a) > std::real(aTransform) && std::real(b) > std::real(bTransform)){
			logHyperU2 = log(hypergo_u_ab_over_u_ap1bp1(a, b, z));
			// if the continued fraction method appears to be unstable (i.e. CF ~ 1) then
			// explicitly calculate U(a + 1, b + 1, z)
			if(abs(logHyperU2.getValue()) < 4.8){
				logHyperU2 = log_hypergeo_u_series(a + 1., b + 1., z);
			}else{
				logHyperU2 = logHyperU - logHyperU2;
			}
			//std::cout << "U("<<a + 1.<<", "<<b + 1.<<", "<<z<<") = " << exp(logHyperU2) << "\n";
		}
		while(std::real(a) > std::real(aTransform) && std::real(b) > std::real(bTransform)){
			logHyperU3 = logHyperU2;
			logHyperU2 = logHyperU;
			// attempt using recurrence relation first because they are cheaper to evaluate
			recurTest = log(hypergo_u_am1bm1(exp(logHyperU2), exp(logHyperU3), a, b, z));
			if( abs(std::real((recurTest - logHyperU2).getValue())) < 1.2 ){
				// if the functions are of similar magnitude then these recurrence relations and
				// continued fraction methods are most likely not very stable. Calculate functions
				// directly
				//std::cout << "Use direct calculation \n";
				logHyperU = log_hypergeo_u_series(a - 1., b - 1., z);
			}else if( std::real(recurTest.getValue()) < std::real(logHyperU.getValue())){
				// if functions are decreasing in magnitude make use of the continued fractions
				//std::cout << "Use CF \n";
				logHyperU += log(hypergo_u_ab_over_u_ap1bp1(a - 1., b - 1., z));
			}else{
				// otherwise use value from recurrence
				//std::cout << "Use recurrence \n";
				logHyperU = recurTest;
			}
			a -= 1.;
			b -= 1.;
			//std::cout << "U("<<a<<", "<<b<<", "<<z<<") = " << exp(logHyperU) << "\n";
		}

		// make use of recurrence relations or CF = U(a, b, z)/U(a + 1, b, z)
		if(std::real(a) > std::real(aTransform)){
			logHyperU2 = logHyperU;
			logHyperU2 -= log(hypergo_u_a_over_u_ap1(a, b, z));
		}
		while(std::real(a) > std::real(aTransform)){
			logHyperU3 = logHyperU2;
			logHyperU2 = logHyperU;
			// attempt using recurrence relation first because they are cheaper to evaluate
			recurTest = log(hypergo_u_am1(exp(logHyperU2), exp(logHyperU3), a, b, z));
			if( std::real(recurTest.getValue()) < std::real(logHyperU.getValue()) ){
				// if functions are decreasing in magnitude make use of the continued fractions
				logHyperU += log(hypergo_u_a_over_u_ap1(a - 1., b, z));
			}else{
				// otherwise use value from recurrence
				logHyperU = recurTest;
			}
			a -= 1.;
			//std::cout << "U("<<a<<", "<<b<<", "<<z<<") = " << exp(logHyperU) << "\n";
		}

		// make use of recurrence relations or CF = U(a, b + 1, z)/U(a, b, z)
		if(std::real(b) > std::real(bTransform)){
			logHyperU2 = logHyperU;
			logHyperU2 += log(hypergo_u_bp1_over_u_b(a, b, z));
		}
		while(std::real(b) > std::real(bTransform)){
			logHyperU3 = logHyperU2;
			logHyperU2 = logHyperU;
			recurTest = log(hypergo_u_bm1(exp(logHyperU2), exp(logHyperU3), a, b, z));
			if( std::real(recurTest.getValue()) < std::real(logHyperU.getValue()) ){
				// if functions are decreasing in magnitude make use of the continued fractions
				logHyperU -= log(hypergo_u_bp1_over_u_b(a, b - 1., z));
			}else{
				// otherwise use value from recurrence
				logHyperU = recurTest;
			}
			b -= 1.;
			//std::cout << "U("<<a<<", "<<b<<", "<<z<<") = " << exp(logHyperU) << "\n";
		}
	}
	logHyperU += logPrefactor;

	return logHyperU;
}

Result hypergo_u_ap1bp1(Result u_ab, Result u_am1bm1, Complex a, Complex b, Complex z){
	//std::cout << "recurrence term 1 " << (z - b + 1.)*u_ab << "\n";
	//std::cout << "recurrence term 2 " <<  u_am1bm1 << "\n";
	return ((z - b + 1.)*u_ab - u_am1bm1)/(-a*z);
}

Result hypergo_u_ap1(Result u_a, Result u_am1, Complex a, Complex b, Complex z){
	return -1.*((b - 2.*a -z)*u_a + u_am1)/(a*(a - b + 1.));
}

Result hypergo_u_bp1(Result u_b, Result u_bm1, Complex a, Complex b, Complex z){
	return -1.*((1. - b - z)*u_b + (b - a - 1.)*u_bm1)/(z);
}

Result hypergo_u_am1bm1(Result u_ab, Result u_ap1bp1, Complex a, Complex b, Complex z){
	//std::cout << "recurrence term 1 " << (z - b + 1.)*u_ab << "\n";
	//std::cout << "recurrence term 2 " <<  a*z*u_ap1bp1 << "\n";
	return ((z - b + 1.)*u_ab + a*z*u_ap1bp1);
}

Result hypergo_u_am1(Result u_a, Result u_ap1, Complex a, Complex b, Complex z){
	return -1.*((b - 2.*a -z)*u_a + a*(a - b + 1.)*u_ap1);
}

Result hypergo_u_bm1(Result u_b, Result u_bp1, Complex a, Complex b, Complex z){
	return -1.*((1. - b - z)*u_b + z*u_bp1)/(b - a - 1.);
}

Result log_hypergeo_m_recurrence(Complex aTemp, Complex bTemp, Complex zTemp){
	Complex aTransform, bTransform, zTransform, logPrefactor;
	Result logHyperM(0., DBL_EPSILON);
	if(std::real(zTemp) < 0){
		//if the real part of the argument is negative
		// transform the function so that it is positive
		aTransform = bTemp - aTemp;
		bTransform = bTemp;
		zTransform = -zTemp;
		logPrefactor = zTemp;
	}else{
		// otherwise use the original user-specified values for the parameters
		aTransform = aTemp;
		bTransform = bTemp;
		zTransform = zTemp;
		logPrefactor = 0.;
	}
	double minArgMagnitudeA = 1.;
	double minArgMagnitudeB = 5.;

	Complex a = aTransform, b = bTransform;
	// minimize the magntiude of the parameter arguments then use continued fraction ratios to recursively find the solution value for larger parameter values
	if(std::real(aTransform) > 0){
		while(std::real(a) > minArgMagnitudeA){
			a -= 1.;
		}
		while(std::real(b) > minArgMagnitudeB){
			b -= 1.;
		}
		logHyperM = log(hypergeo_m_taylor_series(a, b, zTransform));
		//std::cout << "log M(a = "<<a<<", b = "<<b<<", z = "<<zTransform<<") = " << logHyperM << "\n";

		// make use of CF = M(a, b, z)/M(a + 1, b + 1, z)
		while(std::real(a) < std::real(aTransform)){
			logHyperM -= log(hypergo_m_ab_over_m_ap1bp1(a, b, zTransform));
			a += 1.;
			b += 1.;
			//std::cout << "log M(a = "<<a<<", b = "<<b<<", z = "<<zTransform<<") = " << logHyperM << "\n";
		}

		// make use of CF = M(a, b, z)/M(a, b + 1, z)
		while(std::real(b) < std::real(bTransform)){
			logHyperM -= log(hypergo_m_b_over_m_bp1(a, b, zTransform));
			b += 1.;
		}

		// make use of CF = M(a, b - 1, z)/M(a, b, z)
		while(std::real(b) > std::real(bTransform)){
			logHyperM += log(hypergo_m_b_over_m_bp1(a, b - 1., zTransform));
			b -= 1.;
		}
	}else{
		while(std::real(a) < -minArgMagnitudeA){
			a += 1.;
		}
		while(std::real(b) < -minArgMagnitudeB){
			b += 1.;
		}

		logHyperM = log(hypergeo_m_taylor_series(a, b, zTransform));
		//std::cout << "log M(a = "<<a<<", b = "<<b<<", z = "<<zTransform<<") = " << logHyperM << "\n";

		// make use of CF = M(a - 1, b - 1, z)/M(a, b, z)
		while(std::real(a) > std::real(aTransform)){
			logHyperM += log(hypergo_m_ab_over_m_ap1bp1(a - 1., b - 1., zTransform));
			a -= 1.;
			b -= 1.;
		}

		// make use of CF = M(a, b, z)/M(a, b + 1, z)
		while(std::real(b) < std::real(bTransform)){
			logHyperM -= log(hypergo_m_b_over_m_bp1(a, b, zTransform));
			b += 1.;
		}

		// make use of CF = M(a, b - 1, z)/M(a, b, z)
		while(std::real(b) > std::real(bTransform)){
			logHyperM += log(hypergo_m_b_over_m_bp1(a, b - 1., zTransform));
			b -= 1.;
		}
	}
	logHyperM += logPrefactor;

	return logHyperM;
}

Result hypergeo_m_taylor_series(Complex aTemp, Complex bTemp, Complex zTemp){
	//std::cout << "HYPEGEO: Using Taylor series expansion for M \n";
	Complex a, b, z, logPrefactor;
	if( std::real(zTemp) < 0. || abs(bTemp - aTemp) < abs(aTemp) ){
		a = bTemp - aTemp;
		b = bTemp;
		z = -zTemp;
		logPrefactor = zTemp;
	}else{
		a = aTemp;
		b = bTemp;
		z = zTemp;
		logPrefactor = 0.;
	}
	int s = 0;
	double maxRealTerm = 0, maxImaginaryTerm = 0;
	Complex previousTerm, term = 1.;
	set_max_compare(maxRealTerm, abs(std::real(term)));
	set_max_compare(maxImaginaryTerm, abs(std::imag(term)));
	Complex hypergeoM = 0;
	double realError, imaginaryError;

	while(s < 4 || s < -std::real(a) || s < -std::real(b)){
		hypergeoM += term;
		s++;
		previousTerm = term;
		term *= z*(a + Complex(s - 1))/(b + Complex(s - 1))/Complex(s);
		set_max_compare(maxRealTerm, abs(std::real(term)));
		set_max_compare(maxImaginaryTerm, abs(std::imag(term)));

	}

	if( abs(std::real(hypergeoM)) > 0. ){
		realError = abs(std::real(term)/std::real(hypergeoM));
	}else{
		realError = 0.;
	}

	if( abs(std::imag(hypergeoM)) > 0. ){
		imaginaryError = abs(std::imag(term)/std::real(hypergeoM));
	}else{
		imaginaryError = 0.;
	}

	while((realError > HYPERGEO_U_EPSILON/100 || imaginaryError > HYPERGEO_U_EPSILON/100 )  && s < 500){
		hypergeoM += term;
		s++;
		previousTerm = term;
		term *= z*(a + Complex(s - 1))/(b + Complex(s - 1))/Complex(s);
		set_max_compare(maxRealTerm, abs(std::real(term)));
		set_max_compare(maxImaginaryTerm, abs(std::imag(term)));
		if( abs(std::real(hypergeoM)) > 0. ){
			realError = abs(std::real(term)/std::real(hypergeoM));
		}else{
			realError = 0.;
		}

		if( abs(std::imag(hypergeoM)) > 0. ){
			imaginaryError = abs(std::imag(term)/std::real(hypergeoM));
		}else{
			imaginaryError = 0.;
		}
	}
	hypergeoM += term;

	if( abs(std::real(hypergeoM)) > 0. ){
		realError = abs(maxRealTerm/std::real(hypergeoM))*DBL_EPSILON;
		set_max_compare(realError, abs(std::real(term)/std::real(hypergeoM)));
	}else{
		realError = 0.;
	}

	if( abs(std::imag(hypergeoM)) > 0. ){
		imaginaryError = abs(maxImaginaryTerm/std::real(hypergeoM))*DBL_EPSILON;
		set_max_compare(imaginaryError, abs(std::imag(term)/std::imag(hypergeoM)));
	}else{
		imaginaryError = 0.;
	}

	//std::cout << "real error = " << realError << "\n";
	//std::cout << "imaginary error = " << imaginaryError << "\n";

	hypergeoM = exp(log(hypergeoM) + logPrefactor);
	//std::cout << "hypergeoM = " << hypergeoM << "\n";

	return Result(hypergeoM, std::real(hypergeoM)*realError + I*std::imag(hypergeoM)*imaginaryError);
}

// this function is not very efficient because you have to sum three series,
// but I'm not sure if there is a more straightforward way of determining the
// accuracy of the asymptotic expansions. I found that the error bounds tend to
// largely overestimate the error of the asymptotic expansion (when it does)
// "converge"
Result hypergeo_u_asymptotic_series(Complex a, Complex b, Complex z){
	Complex U, U1, U2, at, bt;
	double error;
	// see if the asymptotic series satisfies the Kummer differential equation
	at = a;
	bt = b;
	if(abs(a) > 0){
		// we shift the variables to keep the parameters as small as possible
		// when calculating U, U1, and U2 below
		at -= 2;
		bt -= 2;
	}
	U = (hypergeo_u_asymptotic_series_sum(at, bt, z));
	U1 = (hypergeo_u_asymptotic_series_sum(at + 1., bt + 1., z));
	U2 = (hypergeo_u_asymptotic_series_sum(at + 2., bt + 2., z));
	error = abs(1. + exp(log((z - bt)*U1 - U) - log((1. + at)*z*U2)));
	// error is set by how well U, U1, and U2, do not satisfy Kummer's equation
	if(abs(a) > 0){
		// we need to account for the fact that we shifted the parameters by two
		U = U2;
	}

	return Result(U, U*error);
}

Complex hypergeo_u_asymptotic_series_sum(Complex a, Complex b, Complex z){
	int s = 0;
	Complex previousTerm, term = 1.;
	Complex hypergeoU = 0;
	double errorBoundConstantFactor = error_bound_factor(a, b, z);
	double errorBound, previousBound;

	while(s < 6 || s < -std::real(a) || s < -std::real(a - b + 1.)){
		hypergeoU += term;
		s++;
		previousTerm = term;
		term *= (a + Complex(s - 1))*(a - b + 1. + Complex(s - 1))/(-z*Complex(s));
	}

	while(abs(term) > abs(previousTerm) && s < 200){
		hypergeoU += term;
		s++;
		previousTerm = term;
		term *= (a + Complex(s - 1))*(a - b + 1. + Complex(s - 1))/(-z*Complex(s));
	}

	errorBound = abs(term)*errorBoundConstantFactor*error_bound_Cn(s, a, b, z)/abs(hypergeoU);
	previousBound = abs(previousTerm)*errorBoundConstantFactor*error_bound_Cn(s - 1, a, b, z)/abs(hypergeoU);
	while(errorBound > HYPERGEO_U_EPSILON && (errorBound < previousBound || abs(term) <= abs(previousTerm))){
		// sum asymptotic expansion until the tolerence is met or the error bound and successive terms in the series start to grow
		previousBound = errorBound;
		hypergeoU += term;
		s++;
		previousTerm = term;
		term *= (a + Complex(s - 1))*(a - b + 1. + Complex(s - 1))/(-z*Complex(s));
		errorBound = abs(term)*errorBoundConstantFactor*error_bound_Cn(s, a, b, z)/abs(hypergeoU);
	}
	if(errorBound >= previousBound){
		s--;
		hypergeoU -= previousTerm;
		errorBound = previousBound;
	}

	return exp(log(hypergeoU) - a*log(z));
}

// recursion with continued fractions
// M(a, b, z)/M(a + 1, b + 1, z) = CF
Complex acoeff_mab(int n, void* p){
	U_parameters *params = (U_parameters *)p;
	if(n%2 == 1){
		return params->z*(params->a + 0.5*Complex(n) + 0.5)/(params->b + Complex(n))/(params->b + Complex(n) + 1.);
	}else{
		return params->z*(params->a - params->b - 0.5*Complex(n))/(params->b + Complex(n))/(params->b + Complex(n) + 1.);
	};
}
Complex bcoeff_mab(int, void*){
	return 1.;
}

Result hypergo_m_ab_over_m_ap1bp1(Complex a, Complex b, Complex z){
	cf_coeffs CF;
	U_parameters params = {.a = a, .b = b, .z = z};
	Complex Ma_over_Map1bp1;

	CF.a_coeff = &acoeff_mab;
	CF.b_coeff = &bcoeff_mab;
	CF.params = &params;

	cf_solver* cf = cf_solver_alloc(CF);

	int iter = 0;
	int success = 0;
	while(iter <= CF_ITER_MAX && success == 0){
		iter++;
		success = cf_lentz_iterate(cf, 10*DBL_EPSILON);
	}

	// if(iter > CF_ITER_MAX && cf->cf_error > sqrt(DBL_EPSILON)){
	// 	std::cout << "HYPERGEO_U: ERROR: cf_solver for M(a, b, z)/M(a + 1, b + 1, z) with a = "<<a<<", b = "<<b<<" only converged to "<< cf->cf_error << " after " << iter << " iterations \n";
	// }

	Ma_over_Map1bp1 = cf_solver_convergent(cf);
	cf_solver_free(cf);
	Complex error = cf->cf_error;
	if( abs(error) == 0. ){
		error += DBL_EPSILON;
	}
	error *= 1. + Ma_over_Map1bp1;

	return Result(1. + Ma_over_Map1bp1, error);
}

Result hypergo_m_b_over_m_bp1(Complex a, Complex b, Complex z){
	return hypergo_m_ab_over_m_ap1bp1(b - a, b, -z);
}

// recursion with continued fractions
// U(a, b, z)/U(a+1, b, z) = CF
Complex acoeff_ua(int n, void* p){
	U_parameters *params = (U_parameters *)p;
	return (params->a + Complex(n) + 1.)*(params->b - params->a - Complex(n) - 2.);
}
Complex bcoeff_ua(int n, void* p){
	U_parameters *params = (U_parameters *)p;
	return params->b - 2.*params->a - 2.* Complex(n) - 4. - params->z;
}

Result hypergo_u_a_over_u_ap1(Complex a, Complex b, Complex z){
	cf_coeffs CF;
	U_parameters params = {.a = a, .b = b, .z = z};
	Complex Ua_over_Uap1;

	CF.a_coeff = &acoeff_ua;
	CF.b_coeff = &bcoeff_ua;
	CF.params = &params;

	cf_solver* cf = cf_solver_alloc(CF);

	int iter = 0;
	int success = 0;
	while(iter <= CF_ITER_MAX && success == 0){
		iter++;
		success = cf_lentz_iterate(cf, 10*DBL_EPSILON);
	}

	// if(iter > CF_ITER_MAX && cf->cf_error > sqrt(DBL_EPSILON)){
	// 	std::cout << "HYPERGEO_U: ERROR: cf_solver for U(a, b, z)/U(a + 1, b, z) with a = "<<a<<", b = "<<b<<" only converged to "<< cf->cf_error << " after " << iter << " iterations \n";
	// }

	Ua_over_Uap1 = cf_solver_convergent(cf);
	cf_solver_free(cf);
	Complex error = cf->cf_error;
	if( abs(error) == 0. ){
		error += DBL_EPSILON;
	}
	error *= 2.*a - b + 2. + z - Ua_over_Uap1;

	return Result(2.*a - b + 2. + z - Ua_over_Uap1, error);
}

// recursion with continued fractions
// U(a, b + 1, z)/U(a, b, z) = CF
Complex acoeff_ub(int n, void* p){
	U_parameters *params = (U_parameters *)p;
	if(n%2 == 1){
		return (params->a - params->b + 0.5*Complex(n) + 0.5)/params->z;
	}else{
		return (params->a + 0.5*Complex(n))/params->z;
	};
}
Complex bcoeff_ub(int, void*){
	return 1.;
}

Result hypergo_u_bp1_over_u_b(Complex a, Complex b, Complex z){
	cf_coeffs CF;
	U_parameters params = {.a = a, .b = b, .z = z};
	Complex Ubp1_over_Ub;

	CF.a_coeff = &acoeff_ub;
	CF.b_coeff = &bcoeff_ub;
	CF.params = &params;

	cf_solver* cf = cf_solver_alloc(CF);

	int iter = 0;
	int success = 0;
	while(iter <= CF_ITER_MAX && success == 0){
		iter++;
		success = cf_lentz_iterate(cf, 10*DBL_EPSILON);
	}

	// if(iter > CF_ITER_MAX && cf->cf_error > sqrt(DBL_EPSILON)){
	// 	std::cout << "HYPERGEO_U: ERROR: cf_solver for U(a, b + 1, z)/U(a, b, z) with a = "<<a<<", b = "<<b<<" only converged to "<< cf->cf_error << " after " << iter << " iterations \n";
	// }

	Ubp1_over_Ub = cf_solver_convergent(cf);
	cf_solver_free(cf);
	Complex error = cf->cf_error;
	if( abs(error) == 0. ){
		error += DBL_EPSILON;
	}
	error *= 1. + Ubp1_over_Ub;

	return Result(1. + Ubp1_over_Ub, error);
}

Result hypergo_u_ab_over_u_ap1bp1(Complex a, Complex b, Complex z){
	return z*hypergo_u_bp1_over_u_b(a - b + 1., 1. - b, z);
}

// error bound functions

double error_bound_factor(int n, Complex a, Complex b, Complex z){
	return error_bound_Cn(n, a, b, z)*error_bound_factor(a, b, z);
}

double error_bound_factor(Complex a, Complex b, Complex z){
	double alpha = error_bound_alpha(a, b, z);
	double rho = error_bound_rho(a, b, z);
	double C1 = error_bound_Cn(1, a, b, z);
	return 2*alpha*exp(2*alpha*rho*C1/abs(z));
}

Domain error_bound_region(Complex a, Complex b, Complex z){
	double r = abs(b - 2.*a);
	double zR = std::real(z);
	double zI = std::imag(z);
	double zAbs = abs(z);
	if(zR >= r){
		return Region1;
	}else if(zR <= 0 && abs(zI) >= r){
		return Region2;
	}else if(zR < r && zR > 0 && zAbs > r){
		return Region2;
	}else if(zR <= 0 && zAbs > 2*r  && abs(zI) < r){
		return Region3;
	}

	return Region0;
}

double error_bound_Cn(int n, Complex a, Complex b, Complex z){
	switch(error_bound_region(a, b, z)){
		case Region1:
			return error_bound_Cn_1(n, a, b, z);
		case Region2:
			return error_bound_Cn_2(n, a, b, z);
		case Region3:
			return error_bound_Cn_3(n, a, b, z);
		default:
			//std::cout << "HYPERGEO_U: ERROR: Asymptotic expansion is being applied in an invalid domain.\n";
		return 0;
	}
}

double error_bound_Cn_1(int, Complex, Complex, Complex){ return 1.; }
double error_bound_Cn_2(int n, Complex, Complex, Complex){ return chi(n); }
double error_bound_Cn_3(int n, Complex a, Complex b, Complex z){
	double sigma = error_bound_sigma(a, b, z);
	double nu = error_bound_nu(a, b, z);
	return (chi(n) + sigma*nu*nu*n)*pow(nu, n);
}

double error_bound_alpha(Complex a, Complex b, Complex z){
	double sigma = error_bound_sigma(a, b, z);
	if(error_bound_region(a, b, z) == Region3){
		sigma *= error_bound_nu(a, b, z);
	}

	return 1.0/(1.0 - sigma);
}
double error_bound_rho(Complex a, Complex b, Complex z){
	double sigma = error_bound_sigma(a, b, z);
	if(error_bound_region(a, b, z) == Region3){
		sigma *= error_bound_nu(a, b, z);
	}

	return 0.5*abs(2.0*a*a - 2.0*a*b - b) + sigma*(1.0 + 0.25*sigma)/pow(1.0 - sigma, 2);
}
double error_bound_sigma(Complex a, Complex b, Complex z){
	return abs((b - 2.0*a)/z);
}
double error_bound_nu(Complex a, Complex b, Complex z){
	return 1.0/sqrt(0.5 + 0.5*sqrt(1.0 - 4.0*pow(error_bound_sigma(a, b, z), 2)));
}
