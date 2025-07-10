// Hypergeometric functions

#include "hypergeo_f.hpp"

///////////////////////////////////////////////////
// define constants, structures, and other types //
///////////////////////////////////////////////////

#define HYPERGEO_SERIES_MAX 800
#define HYPERGEO_EPSILON DBL_EPSILON

solution_hypergeo_2F1 solution_hypergeo_2F1_default = {
	.val = 0.,
	.err = 0.
};

test_hypergeo_2F1 test_hypergeo_2F1_default = {
	.w1 = 0,
	.w2 = 0,
	.w3 = 0,
	.w4 = 0,
	.w5 = 0,
	.w6 = 0
};

// functions

double hypergeo_2F1_real(double a, double b, double c, double x){
	Complex hypegeo;
	if(hypergeo_2F1_special_value(a, b, c, x, hypegeo)){
		return std::real(hypegeo);
	}else{
		return hypergeo_2F1(a, b, c, x).getRealValue();
	}
}

Complex hypergeo_2F1_complex(Complex a, Complex b, Complex c, Complex x){
	return hypergeo_2F1(a, b, c, x).getValue();
}

int hypergeo_2F1_special_value(Complex a, Complex b, Complex c, Complex z, Complex &hypergeo){
	if( a == c){
		hypergeo = pow(1. - z, -b);
	}else if( b == c){
		hypergeo = pow(1. - z, -a);
	}else if( a - 1. == c ){
		hypergeo = (1. - (1. - b/(a - 1.))*z)*pow(1. - z, -1. - b); // DLMF 15.4.19
	}else if( b - 1. == c ){
		hypergeo = (1. - (1. - a/(b - 1.))*z)*pow(1. - z, -1. - a); // DLMF 15.4.19
	}else if(a == 1. && b == 1. && c == 2.){
		hypergeo = -log(1. - z)/z; // DLMF 15.4.2
 	}else if(c == 1.5){
		if(std::real(z) >= 0){
			if((a == 0.5 && b == 1.) || (a == 1. && b == 0.5)){
				hypergeo = 0.5*log((1. + sqrt(z))/(1. - sqrt(z)))/sqrt(z); // DLMF 15.4.2
			}else if(a == 0.5 && b == 0.5){
				hypergeo = asin(sqrt(z))/sqrt(z); // DLMF 15.4.4
			}
		}else{
			if((a == 0.5 && b == 1.) || (a == 1. && b == 0.5)){
				hypergeo = atan(sqrt(-z))/sqrt(-z); // DLMF 15.4.3
			}else if(a == 0.5 && b == 0.5){
				hypergeo = log(sqrt(-z) + sqrt(1. - z))/sqrt(-z); // DLMF 15.4.5
			}
		}
	}else if(2.*b == c){
		hypergeo = hypergeo_2F1_group1to3(a, b, z);
	}else if(2.*a == c){
		hypergeo = hypergeo_2F1_group1to3(b, a, z);
	}else if(b == 1. - a && std::real(z) < 0.5){
		hypergeo = hypergeo_2F1_group2to3(a, c, z);
	}else{
		return 0;
	}

	return 1;
}

Complex hypergeo_2F1_group1to3(Complex a, Complex b, Complex x){
	return pow(1. - 0.5*x, -a)*hypergeo_2F1(0.5*a, 0.5*a + 0.5, b + 0.5, pow(x/(2. - x), 2)).getValue();
}

Complex hypergeo_2F1_group2to3(Complex a, Complex c, Complex x){
	return pow(1. - 2.*x, 1. - a - c)*pow(1. - x, c - 1.)*hypergeo_2F1(0.5*(a + c), 0.5*(a + c - 1.), c, (4.*x*(x - 1.))*pow((1. - 2.*x), -2)).getValue();
}

Complex hypergeo_2F1_cgamma_ratio(Complex num1, Complex num2, Complex den1, Complex den2){
	// std::cout << "num1 = "<<num1<<" num2 = "<<num2<<" num3 = "<<den1 <<" num4 = "<<den2<<" \n";
	if(!is_negative_integer(num1) && !is_negative_integer(num2) && !is_negative_integer(den1) && !is_negative_integer(den2)){
		return exp(lgamma(num1) + lgamma(num2) - lgamma(den1) - lgamma(den2));
	}else if((!is_negative_integer(num1) || !is_negative_integer(num2)) && (is_negative_integer(den1) && is_negative_integer(den2))){
		return 0.;
	}else if((!is_negative_integer(num1) && !is_negative_integer(num2)) && (is_negative_integer(den1) || is_negative_integer(den2))){
		return 0.;
	}else{
		return exp(lgamma(num1) + lgamma(num2) - lgamma(den1) - lgamma(den2));
	}
}

// Complex hypergeo_2F1_cgamma_ratio(Complex num1, Complex num2, Complex den1, Complex den2){
// 	Complex cgam[4];
// 	Complex arg[4];
// 	int intReduce[4];
// 	arg[0] = num1;
// 	arg[1] = num2;
// 	arg[2] = den1;
// 	arg[3] = den2;
// 	if(abs(num1) < 10. && abs(num2) < 10. && abs(den1) < 10. && abs(den2) < 10.){
// 		return cgamma(num1)*cgamma(num2)/cgamma(den1)/cgamma(den2);
// 	}
//
// 	for(int i = 0; i < 4; i++){
// 		cgam[i] = cgamma(arg[i]);
// 	}
//
// 	Complex cgamRatio = cgam[0]/cgam[2]*cgam[1]/cgam[3];
// 	if( abs(cgamRatio) != 0. && !isinf(abs(cgamRatio)) &&
// 	!isnan(abs(cgamRatio)) ){
// 		return cgamRatio;
// 	}
//
// 	for(int i = 0; i < 4; i++){
// 		if(abs(cgam[i]) == 0 || isinf(abs(cgam[i])) || isnan(abs(cgam[i]))){
// //			printf("Large argument encountered (|z| = %f) for gamma term %d\n",abs(arg[i]),i);
// 			intReduce[i] = round(std::real(arg[i])/1.5);
// //			printf("Reduce by integer %d\n", intReduce[i]);
// 		}else{
// 			intReduce[i] = 0;
// 		}
//
// 		if(intReduce[i] < 0){
// 			cgam[i] =  M_PI/cgamma(1. - arg[i] + Complex(intReduce[i]))
// 				/sin(M_PI*arg[i]);
// 		}else if(intReduce[i] > 0){
// 			cgam[i] = cgamma(arg[i] - Complex(intReduce[i]));
// 		}
// 	}
// 	cgamRatio = cgam[0]/cgam[2]*cgam[1]/cgam[3];
//
// 	for(int i = 0; i < 2; i++){
// 		if(intReduce[i] < 0){
// 			cgamRatio = 1./times_phammer(1./cgamRatio,
// 				1. - arg[i] + Complex(intReduce[i]), -intReduce[i]);
// 		}else if(intReduce[i] > 0){
// 			cgamRatio = times_phammer(cgamRatio, arg[i] - Complex(intReduce[i]),
// 				intReduce[i]);
// 		}
// 	}
//
// 	for(int i = 2; i < 4; i++){
// 		if(intReduce[i] < 0){
// 			cgamRatio = times_phammer(cgamRatio, 1. - arg[i] + Complex(intReduce[i]),
// 				-intReduce[i]);
// 		}else if(intReduce[i] > 0){
// 			cgamRatio = 1./times_phammer(1./cgamRatio,
// 				arg[i] - Complex(intReduce[i]), intReduce[i]);
// 		}
// 	}
//
// 	return cgamRatio;
// }

int hypergeo_2F1_series(solution_hypergeo_2F1* sol2F1, Complex a, Complex b, Complex c, Complex z){
	if( (is_negative_integer(std::real(c)) && abs(std::imag(c)/c) < 10.*DBL_EPSILON) ){
		sol2F1->val = INFINITY;
		sol2F1->err = abs(sol2F1->val)*DBL_EPSILON;
		return 0;
	}

	Complex an, anp1, ii;
	Complex f21n;
	double tol, test[3];
	int i = 0;

	if( a == c && abs(z) < 1.){
		sol2F1->val = pow(1. - z, -b);
		sol2F1->err = abs(sol2F1->val)*DBL_EPSILON;
		return 1;
	}else if( b == c && abs(z) < 1.){
		sol2F1->val = pow(1. - z, -a);
		sol2F1->err = abs(sol2F1->val)*DBL_EPSILON;
		return 1;
	}
	// else if( a == c + 1. && abs(z) < 1.){
	// 	sol2F1->val = (1. - (1. - b/a)*z)*pow(1. - z, -1. - b);
	// 	sol2F1->err = abs(sol2F1->val)*DBL_EPSILON;
	// 	return 1;
	// }

	Complex termMax;

	an = 1.;
	f21n = an;
	termMax = f21n;
	while(i < 3){
		ii = Complex(i);
		anp1 = an*(a + ii)*(b + ii)/(c + ii)*z/(ii + 1.);
		f21n += anp1;
		an = anp1;
		termMax = (abs(termMax) < abs(f21n)) ? f21n : termMax;
		test[i] = abs(an/f21n);
		i++;
	}

	tol = test[0] + test[1] + test[2];
	while(i < HYPERGEO_SERIES_MAX && tol > HYPERGEO_EPSILON){
		ii = Complex(i);
		anp1 = an*(a + ii)*(b + ii)/(c + ii)*z/(ii + 1.);
		f21n += anp1;
		an = anp1;
		termMax = (abs(termMax) < abs(f21n)) ? f21n : termMax;
		test[0] = test[1];
		test[1] = test[2];
		test[2] = abs(an)/abs(f21n);
		if(test[2] > test[1] || test[1] > test[0] || test[2] > test[0]){
			tol = 1.;
		}else{
			tol = test[0] + test[1] + test[2];
		}
		i++;
	}
	sol2F1->val = f21n;
	sol2F1->err = std::max(abs(termMax/f21n)*DBL_EPSILON, tol);
	//printf("2F1 required %d terms; exited with tolerance of %e and error of %e\n", i, tol, sol2F1->err);
	if( sol2F1->err < DBL_EPSILON ){
		sol2F1->err = DBL_EPSILON;
	}

	return 1;
}

int hypergeo_2F1_test_sign(Complex a, Complex b, Complex c, Complex z){
	if( abs(z) > 1. || is_negative_integer(c) ){
		return 0;
	}

	Complex n = 20.;
	// is a or b are a negative integer, then the series will terminate at some finite n
	if(is_negative_integer(a)){
		return 1;
	}
	if(is_negative_integer(b)){
		return 1;
	}
	Complex test = (a + n)*(b + n)/(c + n)*z/n;

	if( std::real(test) < 0. ){
		return -1;
	}else{
		return 1;
	}
}

double hypergeo_2F1_test_weighted(Complex a, Complex b, Complex c, Complex z){
	if( abs(z) > 1. || (is_negative_integer(std::real(c)) && abs(std::imag(c)/c) < 10.*DBL_EPSILON) ){
		return 0.;
	}

	Complex n = 20.;
	Complex test = (a + n)*(b + n)/(c + n)*z/n;
	// is a or b are a negative integer, then the series will terminate at some finite n
	if(is_negative_integer(a)){
		return 1.e-10;
	}
	if(is_negative_integer(b)){
		return 1.e-10;
	}

	if( std::real(test) < 0. ){
		return -abs(test);
	}else{
		return abs(test);
	}
}

//////////////////
//DLMF 15.10.11 //
//////////////////
int hypergeo_2F1_test_w1_z_1(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(a, b, c, z);
}

int hypergeo_2F1_test_w1_z_2(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(c - a, c - b, c, z);
}

int hypergeo_2F1_test_w1_zoverzm1_1(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(a, c - b, c, z/(z - 1.));
}

int hypergeo_2F1_test_w1_zoverzm1_2(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(c - a, b, c, z/(z - 1.));
}

double hypergeo_2F1_test_w1_z_1_weighted(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_weighted(a, b, c, z);
}

double hypergeo_2F1_test_w1_z_2_weighted(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_weighted(c - a, c - b, c, z);
}

double hypergeo_2F1_test_w1_zoverzm1_1_weighted(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_weighted(a, c - b, c, z/(z - 1.));
}

double hypergeo_2F1_test_w1_zoverzm1_2_weighted(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_weighted(c - a, b, c, z/(z - 1.));
}

int hypergeo_2F1_w1_z(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	if(abs(z) > 1.){
		return 0;
	}

	hypergeo_2F1_series(f21, a, b, c, z);
	return 1;
}

int hypergeo_2F1_w1_z_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	int success = hypergeo_2F1_w1_z(f21, c - a, c - b, c, z);
	f21->val *= pow(1. - z, c - a - b);
	return success;
}

int hypergeo_2F1_w1_zoverzm1_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	int success = hypergeo_2F1_w1_z(f21, a, c - b, c, z/(z - 1.));
	f21->val *= pow(1. - z, -a);
	return success;
}

int hypergeo_2F1_w1_zoverzm1_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	int success = hypergeo_2F1_w1_z(f21, c - a, b, c, z/(z - 1.));
	f21->val *= pow(1. - z, -b);
	return success;
}

Result hypergeo_2F1_w1(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z){
	int success;
	solution_hypergeo_2F1 f21 = solution_hypergeo_2F1_default;

	switch(test.w1){
		case 1:
			success = hypergeo_2F1_w1_z(&f21, a, b, c, z);
			break;

		case 2:
			success = hypergeo_2F1_w1_z_2(&f21, a, b, c, z);
			break;

		case 3:
			success = hypergeo_2F1_w1_zoverzm1_1(&f21, a, b, c, z);
			break;

		case 4:
			success = hypergeo_2F1_w1_zoverzm1_2(&f21, a, b, c, z);
			break;

		default:
			return Result(0., 0.);
	}

	if( success == 0 ){
		std::cout << "HYPERGEO_F: ERROR: Unsuccessful summation of hypergeometric series w_1 \n";
	}

	return Result(f21.val, f21.val*f21.err);
}

//////////////////
//DLMF 15.10.12 //
//////////////////
int hypergeo_2F1_test_w2_z_1(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(a - c + 1., b - c + 1., 2. - c, z);
}

int hypergeo_2F1_test_w2_z_2(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(1. - a, 1. - b, 2. - c, z);
}

int hypergeo_2F1_test_w2_zoverzm1_1(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(a - c + 1., 1. - b, 2. - c, z/(z - 1.));
}

int hypergeo_2F1_test_w2_zoverzm1_2(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(1. - a, b - c + 1., 2. - c, z/(z - 1.));
}

double hypergeo_2F1_test_w2_z_1_weighted(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_weighted(a - c + 1., b - c + 1., 2. - c, z);
}

double hypergeo_2F1_test_w2_z_2_weighted(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_weighted(1. - a, 1. - b, 2. - c, z);
}

double hypergeo_2F1_test_w2_zoverzm1_1_weighted(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_weighted(a - c + 1., 1. - b, 2. - c, z/(z - 1.));
}

double hypergeo_2F1_test_w2_zoverzm1_2_weighted(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_weighted(1. - a, b - c + 1., 2. - c, z/(z - 1.));
}

int hypergeo_2F1_w2_z_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	int success = hypergeo_2F1_w1_z(f21, a - c + 1., b - c + 1., 2. - c, z);
	f21->val *= pow(z, 1. - c);
	return success;
}

int hypergeo_2F1_w2_z_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	int success = hypergeo_2F1_w1_z(f21, 1. - a, 1. - b, 2. - c, z);
	f21->val *= pow(z, 1. - c)*pow(1. - z, c - a - b);
	return success;
}

int hypergeo_2F1_w2_zoverzm1_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	// if(abs(z) > 0.5){
	// 	return 0;
	// }

	int success = hypergeo_2F1_w1_z(f21, a - c + 1., 1. - b, 2. - c, z/(z - 1.));
	f21->val *= pow(z, 1. - c)*pow(1. - z, c - a - 1.);
	return success;
}

int hypergeo_2F1_w2_zoverzm1_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	// if(abs(z) > 0.5){
	// 	return 0;
	// }

	int success = hypergeo_2F1_w1_z(f21, 1. - a, b - c + 1., 2. - c, z/(z - 1.));
	f21->val *= pow(z, 1. - c)*pow(1. - z, c - b - 1.);
	return success;
}

Result hypergeo_2F1_w2(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z){
	int success;
	solution_hypergeo_2F1 f21 = solution_hypergeo_2F1_default;

	switch(test.w2){
		case 1:
			success = hypergeo_2F1_w2_z_1(&f21, a, b, c, z);
			break;

		case 2:
			success = hypergeo_2F1_w2_z_2(&f21, a, b, c, z);
			break;

		case 3:
			success = hypergeo_2F1_w2_zoverzm1_1(&f21, a, b, c, z);
			break;

		case 4:
			success = hypergeo_2F1_w2_zoverzm1_2(&f21, a, b, c, z);
			break;

		default:
			return Result(0., 0.);
	}

	if( success == 0 ){
		std::cout << "HYPERGEO_F: ERROR: Unsuccessful summation of hypergeometric series w_2 \n";
	}

	return Result(f21.val, f21.val*f21.err);
}

//////////////////
//DLMF 15.10.13 //
//////////////////
int hypergeo_2F1_test_w3_1mz_1(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(a, b, a + b - c + 1., 1. - z);
}

int hypergeo_2F1_test_w3_1mz_2(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(a - c + 1., b - c + 1., a + b - c + 1., 1. - z);
}

int hypergeo_2F1_test_w3_1m1overz_1(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(a, a - c + 1., a + b - c + 1., 1. - 1./z);
}

int hypergeo_2F1_test_w3_1m1overz_2(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(b, b - c + 1., a + b - c + 1., 1. - 1./z);
}

int hypergeo_2F1_w3_1mz_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	int success = hypergeo_2F1_w1_z(f21, a, b, a + b - c + 1., 1. - z);
	f21->val *= pow(z, 1. - c);
	return success;
}

int hypergeo_2F1_w3_1mz_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	int success = hypergeo_2F1_w1_z(f21, a - c + 1., b - c + 1., a + b - c + 1.,
		1. - z);
	return success;
}

int hypergeo_2F1_w3_1m1overz_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	// if(abs(1. - z) > 0.5){
	// 	return 0;
	// }

	int success = hypergeo_2F1_w1_z(f21, a, a - c + 1., a + b - c + 1.,
		1. - 1./z);
	f21->val *= pow(z, -a);
	return success;
}

int hypergeo_2F1_w3_1m1overz_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	// if(abs(1. - z) > 0.5){
	// 	return 0;
	// }

	int success = hypergeo_2F1_w1_z(f21, b, b - c + 1., a + b - c + 1.,
		1. - 1./z);
	f21->val *= pow(z, -b);
	return success;
}

Result hypergeo_2F1_w3(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z){
	int success;
	solution_hypergeo_2F1 f21 = solution_hypergeo_2F1_default;

	switch(test.w3){
		case 1:
			success = hypergeo_2F1_w3_1mz_1(&f21, a, b, c, z);
			break;

		case 2:
			success = hypergeo_2F1_w3_1mz_2(&f21, a, b, c, z);
			break;

		case 3:
			success = hypergeo_2F1_w3_1m1overz_1(&f21, a, b, c, z);
			break;

		case 4:
			success = hypergeo_2F1_w3_1m1overz_2(&f21, a, b, c, z);
			break;

		default:
			return Result(0., 0.);
	}

	if( success == 0 ){
		std::cout << "HYPERGEO_F: ERROR: Unsuccessful summation of hypergeometric series w_3 \n";
	}

	return Result(f21.val, f21.val*f21.err);
}

//////////////////
//DLMF 15.10.14 //
//////////////////
int hypergeo_2F1_test_w4_1mz_1(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(c - a, c - b, c - a - b + 1., 1. - z);
}

int hypergeo_2F1_test_w4_1mz_2(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(1. - a, 1. - b, c - a - b + 1., 1. - z);
}

int hypergeo_2F1_test_w4_1m1overz_1(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(1. - a, c - a, c - a - b + 1., 1. - 1./z);
}

int hypergeo_2F1_test_w4_1m1overz_2(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(1. - b, c - b, c - a - b + 1., 1. - 1./z);
}

int hypergeo_2F1_w4_1mz_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	int success = hypergeo_2F1_w1_z(f21, c - a, c - b, c - a - b + 1., 1. - z);
	f21->val *= pow(1. - z, c - a - b);
	return success;
}

int hypergeo_2F1_w4_1mz_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	int success = hypergeo_2F1_w1_z(f21, 1. - a, 1. - b, c - a - b + 1., 1. - z);
	f21->val *= pow(z, 1. - c)*pow(1. - z, c - a - b);
	return success;
}

int hypergeo_2F1_w4_1m1overz_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	// if(abs(1. - z) > 0.5){
	// 	return 0;
	// }

	int success = hypergeo_2F1_w1_z(f21, 1. - a, c - a, c - a - b + 1.,
		1. - 1./z);
	f21->val *= pow(z, a - c)*pow(1. - z, c - a - b);
	return success;
}

int hypergeo_2F1_w4_1m1overz_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	// if(abs(1. - z) > 0.5){
	// 	return 0;
	// }

	int success = hypergeo_2F1_w1_z(f21, 1. - b, c - b, c - a - b + 1.,
		1. - 1./z);
	f21->val *= pow(z, b - c)*pow(1. - z, c - a - b);
	return success;
}

Result hypergeo_2F1_w4(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z){
	int success;
	solution_hypergeo_2F1 f21 = solution_hypergeo_2F1_default;

	switch(test.w4){
		case 1:
			success = hypergeo_2F1_w4_1mz_1(&f21, a, b, c, z);
			break;

		case 2:
			success = hypergeo_2F1_w4_1mz_2(&f21, a, b, c, z);
			break;

		case 3:
			success = hypergeo_2F1_w4_1m1overz_1(&f21, a, b, c, z);
			break;

		case 4:
			success = hypergeo_2F1_w4_1m1overz_2(&f21, a, b, c, z);
			break;

		default:
			return Result(0., 0.);
	}

	if( success == 0 ){
		std::cout << "HYPERGEO_F: ERROR: Unsuccessful summation of hypergeometric series w_4 \n";
	}

	return Result(f21.val, f21.val*f21.err);
}

//////////////////
//DLMF 15.10.15 //
//////////////////
int hypergeo_2F1_test_w5_1overz_1(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(a, a - c + 1., a - b + 1., 1./z);
}

int hypergeo_2F1_test_w5_1overz_2(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(1. - b, c - b, a - b + 1., 1./z);
}

int hypergeo_2F1_test_w5_1over1mz_1(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(a, c - b, a - b + 1., 1./(1. - z));
}

int hypergeo_2F1_test_w5_1over1mz_2(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(1. - b, a - c + 1., a - b + 1., 1./(1. - z));
}

int hypergeo_2F1_w5_1overz_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	int success = hypergeo_2F1_w1_z(f21, a, a - c + 1., a - b + 1., 1./z);
	f21->val *= pow(-z, -a);
	return success;
}

int hypergeo_2F1_w5_1overz_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	int success = hypergeo_2F1_w1_z(f21, 1. - b, c - b, a - b + 1., 1./z);
	f21->val *= pow(-z, b - c)*pow(1. - z, c - a - b);
	return success;
}

int hypergeo_2F1_w5_1over1mz_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	// if(abs(1./z) > 0.5){
	// 	return 0;
	// }

	int success = hypergeo_2F1_w1_z(f21, a, c - b, a - b + 1., 1./(1. - z));
	f21->val *= pow(1. - z, -a);
	return success;
}

int hypergeo_2F1_w5_1over1mz_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	// if(abs(1./z) > 0.5){
	// 	return 0;
	// }

	int success = hypergeo_2F1_w1_z(f21, 1. - b, a - c + 1., a - b + 1.,
		1./(1. - z));
	f21->val *= pow(-z, 1. - c)*pow(1. - z, c - a - 1.);
	return success;
}

Result hypergeo_2F1_w5(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z){
	int success;
	solution_hypergeo_2F1 f21 = solution_hypergeo_2F1_default;

	switch(test.w5){
		case 1:
			success = hypergeo_2F1_w5_1overz_1(&f21, a, b, c, z);
			break;

		case 2:
			success = hypergeo_2F1_w5_1overz_2(&f21, a, b, c, z);
			break;

		case 3:
			success = hypergeo_2F1_w5_1over1mz_1(&f21, a, b, c, z);
			break;

		case 4:
			success = hypergeo_2F1_w5_1over1mz_2(&f21, a, b, c, z);
			break;

		default:
			return Result(0., 0.);
	}

	if( success == 0 ){
		std::cout << "HYPERGEO_F: ERROR: Unsuccessful summation of hypergeometric series w_5 \n";
	}

	return Result(f21.val, f21.val*f21.err);
}

//////////////////
//DLMF 15.10.16 //
//////////////////
int hypergeo_2F1_test_w6_1overz_1(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(b, b - c + 1., b - a + 1., 1./z);
}

int hypergeo_2F1_test_w6_1overz_2(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(1. - a, c - a, b - a + 1., 1./z);
}

int hypergeo_2F1_test_w6_1over1mz_1(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(b, c - a, b - a + 1., 1./(1. - z));
}

int hypergeo_2F1_test_w6_1over1mz_2(Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_test_sign(1. - a, b - c + 1., b - a + 1., 1./(1. - z));
}

int hypergeo_2F1_w6_1overz_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	int success = hypergeo_2F1_w1_z(f21, b, b - c + 1., b - a + 1., 1./z);
	f21->val *= pow(-z, -b);
	return success;
}

int hypergeo_2F1_w6_1overz_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	int success = hypergeo_2F1_w1_z(f21, 1. - a, c - a, b - a + 1., 1./z);
	f21->val *= pow(-z, a - c)*pow(1. - z, c - a - b);
	return success;
}

int hypergeo_2F1_w6_1over1mz_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	// if(abs(1./z) > 0.5){
	// 	return 0;
	// }

	int success = hypergeo_2F1_w1_z(f21, b, c - a, b - a + 1., 1./(1. - z));
	f21->val *= pow(1. - z, -b);
	return success;
}

int hypergeo_2F1_w6_1over1mz_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z){
	// if(abs(1./z) > 0.5){
	// 	return 0;
	// }

	int success = hypergeo_2F1_w1_z(f21, 1. - a, b - c + 1., b - a + 1.,
		1./(1. - z));
	f21->val *= pow(-z, 1. - c)*pow(1. - z, c - b - 1.);
	return success;
}

Result hypergeo_2F1_w6(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z){
	int success;
	solution_hypergeo_2F1 f21 = solution_hypergeo_2F1_default;

	switch(test.w6){
		case 1:
			success = hypergeo_2F1_w6_1overz_1(&f21, a, b, c, z);
			break;

		case 2:
			success = hypergeo_2F1_w6_1overz_2(&f21, a, b, c, z);
			break;

		case 3:
			success = hypergeo_2F1_w6_1over1mz_1(&f21, a, b, c, z);
			break;

		case 4:
			success = hypergeo_2F1_w6_1over1mz_2(&f21, a, b, c, z);
			break;

		default:
			return Result(0., 0.);
	}

	if( success == 0 ){
		std::cout << "HYPERGEO_F: ERROR: Unsuccessful summation of hypergeometric series w_6 \n";
	}

	return Result(f21.val, f21.val*f21.err);
}

// DLMF Equations for w1

Result hypergeo_2F1_w1_151011(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z){
	return hypergeo_2F1_w1(test, a, b, c, z);
}

Result hypergeo_2F1_w1_151021(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z){
	Complex gam1, gam2;

	gam1 = hypergeo_2F1_cgamma_ratio(c, c - a - b, c - a, c - b);
	gam2 = hypergeo_2F1_cgamma_ratio(c, a + b - c, a, b);

	if(abs(gam1) == 0.){
		Result w4 = hypergeo_2F1_w4(test, a, b, c, z);
		return gam2*w4;
	}else if(abs(gam2) == 0.){
		Result w3 = hypergeo_2F1_w3(test, a, b, c, z);
		return gam1*w3;
	}

	Result w3 = hypergeo_2F1_w3(test, a, b, c, z);
	Result w4 = hypergeo_2F1_w4(test, a, b, c, z);
	//std::cout << "w3 = " << w3 << "\n";
	//std::cout << "w4 = " << w4 << "\n";
	//std::cout << "w1 from combo = " << gam1*w3 + gam2*w4 << "\n";
	//std::cout << "precision from combo = " << (gam1*w3 + gam2*w4).getPrecision() << "\n";

	return gam1*w3 + gam2*w4;
}

Result hypergeo_2F1_w1_151025(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z){
	Complex gam1, gam2;
	gam1 = hypergeo_2F1_cgamma_ratio(c, b - a, b, c - a);
	gam2 = hypergeo_2F1_cgamma_ratio(c, a - b, a, c - b);

	if(abs(gam1) == 0.){
		Result w6 = hypergeo_2F1_w6(test, a, b, c, z);
		return gam2*w6;
	}else if(abs(gam2) == 0.){
		Result w5 = hypergeo_2F1_w5(test, a, b, c, z);
		return gam1*w5;
	}

	Result w5 = hypergeo_2F1_w5(test, a, b, c, z);
	Result w6 = hypergeo_2F1_w6(test, a, b, c, z);

	return gam1*w5 + gam2*w6;
}

Result hypergeo_2F1_w1_151029(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z){
	Complex gam1, gam2;

	Result w3 = hypergeo_2F1_w3(test, a, b, c, z);
	Result w5 = hypergeo_2F1_w5(test, a, b, c, z);
	gam1 = exp(M_PI*b*I)*hypergeo_2F1_cgamma_ratio(c, a - c + 1.,
		a + b - c + 1., c-a);
	gam2 = exp(M_PI*(b-c)*I)*hypergeo_2F1_cgamma_ratio(c, a - c + 1., b,
		a - b + 1.);

	return gam1*w3 + gam2*w5;
}

Result hypergeo_2F1_w1_151030(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z){
	Complex gam1, gam2;

	Result w3 = hypergeo_2F1_w3(test, a, b, c, z);
	Result w6 = hypergeo_2F1_w6(test, a, b, c, z);
	gam1 = exp(M_PI*a*I)*hypergeo_2F1_cgamma_ratio(c, b - c + 1.,
		a + b - c + 1., c - b);
	gam2 = exp(M_PI*(a-c)*I)*hypergeo_2F1_cgamma_ratio(c, b - c + 1., a,
		b - a + 1.);

	return gam1*w3 + gam2*w6;
}

Result hypergeo_2F1_w1_151033(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z){
	Complex gam1, gam2;

	Result w4 = hypergeo_2F1_w4(test, a, b, c, z);
	Result w5 = hypergeo_2F1_w5(test, a, b, c, z);
	gam1 = exp(M_PI*(c-a)*I)*hypergeo_2F1_cgamma_ratio(c, 1. - b, a,
		c - a - b + 1.);
	gam2 = exp(-M_PI*a*I)*hypergeo_2F1_cgamma_ratio(c, 1. - b, a - b + 1.,
		c - a);

	return gam1*w4 + gam2*w5;
}

Result hypergeo_2F1_w1_151034(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z){
	Complex gam1, gam2;

	Result w4 = hypergeo_2F1_w4(test, a, b, c, z);
	Result w6 = hypergeo_2F1_w6(test, a, b, c, z);
	gam1 = exp(M_PI*(c-b)*I)*hypergeo_2F1_cgamma_ratio(c, 1. - a, b,
		c - a - b + 1.);
	gam2 = exp(-M_PI*b*I)*hypergeo_2F1_cgamma_ratio(c, 1. - a, b - a + 1.,
		c - b);

	return gam1*w4 + gam2*w6;
}

// Tests

int hypergeo_2F1_arg_test(Complex z){
	if(abs(z) == 1){
		return 6;
	}

	int testNum = 0;
	double test[6];
	test[0] = abs(z);
	test[1] = abs(z/(z - 1.));
	test[2] = abs(1. - z);
	test[3] = abs(1./(1. - z));
	test[4] = abs(1./z);
	test[5] = abs(1. - 1./z);
	double testTemp;
	bool testFlag = false;
	int i;

	while(testNum < 5 && testFlag == false){
		i = 1;
		testTemp = test[0];
		testNum++;
		test[0] = test[testNum];
		test[testNum] = testTemp;
		testFlag = test[0] < test[i];
		while(i < 5 && testFlag == true){
			testFlag = test[0] < test[i];
			i++;
		}
	}

	return testNum;
}

int hypergeo_2F1_test_w1(Complex a, Complex b, Complex c, Complex z, double testLimit){
	int testNum = 0;
	double testTemp;
	double test[4];
	// test numerical stability of different series representations of w3
	// if 0 then series does not converge, if -1 then series rapidly oscillates
	// in sign, leading to large cancellations, and if +1 then the series may
	// be relatively stable.
//	test[0] =	hypergeo_2F1_test_w1_z_1(a, b, c, z)*abs(z);
//	test[1] =	hypergeo_2F1_test_w1_z_2(a, b, c, z)*abs(z);
//	test[2] =	hypergeo_2F1_test_w1_zoverzm1_1(a, b, c, z)*abs(z/(z-1));
//	test[3] =	hypergeo_2F1_test_w1_zoverzm1_2(a, b, c, z)*abs(z/(z-1));

	test[0] =	hypergeo_2F1_test_w1_z_1_weighted(a, b, c, z);
	test[1] =	hypergeo_2F1_test_w1_z_2_weighted(a, b, c, z);
	test[2] =	hypergeo_2F1_test_w1_zoverzm1_1_weighted(a, b, c, z);
	test[3] =	hypergeo_2F1_test_w1_zoverzm1_2_weighted(a, b, c, z);
	// we weight each test with the magnitude of its argument. Typically the
	// smaller the argument, the more convergent the series.
	bool testFlag = false;
	int i;
	// printf("Weighted test results for w1:\n %f, %f, %f, %f\n", test[0], test[1], test[2], test[3]);

	// if a solution seems numerically unstable (test < 0) or is not convergent
	// (test == 0), then set test to greater than or equal to 1.
	for(int j = 0; j < 4; j++){
		if( test[j] <= 0 ){
			test[j] = 1 - test[j];
		}
		// printf("Test[%d] = %.10f\n", j, test[j]);
	}

	// find series with smallest weighted test result that is less than 1
	// (the less than one condition ensures that we choose a convergent series)
	while(testNum < 4 && testFlag == false){
		i = 1;
		testTemp = test[0];
		test[0] = test[testNum];
		test[testNum] = testTemp;
		testFlag = (test[0] <= test[i]) && (test[0] < testLimit);
		while(i < 4 && testFlag == true){
			testFlag = (test[0] <= test[i]) && (test[0] < testLimit);
			i++;
		}
		testNum++;
	}
	// note that the test number is actually one greater than array number

	// if none of our tests estimate a numerically stable or convergent series,
	// return 0, which tells us that none of tests were successful
	if( testFlag == false ){
		testNum = 0;
	}

	return testNum;
}

int hypergeo_2F1_test_w3(Complex a, Complex b, Complex c, Complex z, double testLimit){
	int testNum = 0;
	double testTemp;
	double test[4];
	// test numerical stability of different series representations of w3
	// if 0 then series does not converge, if -1 then series rapidly oscillates
	// in sign, leading to large cancellations, and if +1 then the series may
	// be relatively stable.
	test[0] =	hypergeo_2F1_test_w3_1mz_1(a, b, c, z)*abs(1. - z);
	test[1] =	hypergeo_2F1_test_w3_1mz_2(a, b, c, z)*abs(1. - z);
	if(abs(1. - 1./z) > 1.){
		test[2] =	0.;
		test[3] = 0.;
	}else{
		test[2] =	hypergeo_2F1_test_w3_1m1overz_1(a, b, c, z)*abs(1. - 1./z);
		test[3] =	hypergeo_2F1_test_w3_1m1overz_2(a, b, c, z)*abs(1. - 1./z);
	}
	// we weight each test with the magnitude of its argument. Typically the
	// smaller the argument, the more convergent the series.
	bool testFlag = false;
	int i;
	// printf("Weighted test results for w3:\n %f, %f, %f, %f\n", test[0], test[1], test[2], test[3]);
	// printf("Parameters for w3:\n %f, %f, %f, %f\n", std::real(a), std::real(b), std::real(c), std::real(z));
	// find series with smallest weighted test result that is less than 1
	// (the less than one condition ensures that we choose a convergent series)
	for(int j = 0; j < 4; j++){
		if( test[j] <= 0 ){
			test[j] = 1 - test[j];
		}
	}

	// find series with smallest weighted test result that is less than 1
	// (the less than one condition ensures that we choose a convergent series)
	while(testNum < 4 && testFlag == false){
		i = 1;
		testTemp = test[0];
		test[0] = test[testNum];
		test[testNum] = testTemp;
		testFlag = (test[0] <= test[i]) && (test[0] < testLimit);
		while(i < 4 && testFlag == true){
			testFlag = (test[0] <= test[i]) && (test[0] < testLimit);
			i++;
		}
		testNum++;
	}
	// note that the test number is actually one greater than array number

	// if none of our tests estimate a numerically stable or convergent series,
	// return 0, which tells us that none of tests were successful
	if( testFlag == false ){
		testNum = 0;
	}

	return testNum;
}

int hypergeo_2F1_test_w4(Complex a, Complex b, Complex c, Complex z, double testLimit){
	int testNum = 0;
	double testTemp;
	double test[4];
	// test numerical stability of different series representations of w3
	// if 0 then series does not converge, if -1 then series rapidly oscillates
	// in sign, leading to large cancellations, and if +1 then the series may
	// be relatively stable.
	test[0] =	hypergeo_2F1_test_w4_1mz_1(a, b, c, z)*abs(1. - z);
	test[1] =	hypergeo_2F1_test_w4_1mz_2(a, b, c, z)*abs(1. - z);
	if(abs(1. - 1./z) > 1.){
		test[2] =	0.;
		test[3] = 0.;
	}else{
		test[2] =	hypergeo_2F1_test_w4_1m1overz_1(a, b, c, z)*abs(1. - 1./z);
		test[3] =	hypergeo_2F1_test_w4_1m1overz_2(a, b, c, z)*abs(1. - 1./z);
	}
	// we weight each test with the magnitude of its argument. Typically the
	// smaller the argument, the more convergent the series.
	bool testFlag = false;
	int i;
	// printf("Weighted test results for w4:\n %f, %f, %f, %f\n", test[0], test[1], test[2], test[3]);

	// find series with smallest weighted test result that is less than 1
	// (the less than one condition ensures that we choose a convergent series)
	for(int j = 0; j < 4; j++){
		if( test[j] <= 0 ){
			test[j] = 1 - test[j];
		}
	}

	// find series with smallest weighted test result that is less than 1
	// (the less than one condition ensures that we choose a convergent series)
	while(testNum < 4 && testFlag == false){
		i = 1;
		testTemp = test[0];
		test[0] = test[testNum];
		test[testNum] = testTemp;
		testFlag = (test[0] <= test[i]) && (test[0] < testLimit);
		while(i < 4 && testFlag == true){
			testFlag = (test[0] <= test[i]) && (test[0] < testLimit);
			i++;
		}
		testNum++;
	}
	// note that the test number is actually one greater than array number

	// if none of our tests estimate a numerically stable or convergent series,
	// return 0, which tells us that none of tests were successful
	if( testFlag == false ){
		testNum = 0;
	}

	return testNum;
}

int hypergeo_2F1_test_w5(Complex a, Complex b, Complex c, Complex z, double testLimit){
	int testNum = 0;
	double testTemp;
	double test[4];
	// test numerical stability of different series representations of w3
	// if 0 then series does not converge, if -1 then series rapidly oscillates
	// in sign, leading to large cancellations, and if +1 then the series may
	// be relatively stable.
	test[0] =	hypergeo_2F1_test_w5_1overz_1(a, b, c, z)*abs(1./z);
	test[1] =	hypergeo_2F1_test_w5_1overz_2(a, b, c, z)*abs(1./z);
	if(abs(1./(1. - z)) > 1.){
		test[2] =	0.;
		test[3] = 0.;
	}else{
		test[2] =	hypergeo_2F1_test_w5_1over1mz_1(a, b, c, z)*abs(1./(1. - z));
		test[3] =	hypergeo_2F1_test_w5_1over1mz_2(a, b, c, z)*abs(1./(1. - z));
	}
	// we weight each test with the magnitude of its argument. Typically the
	// smaller the argument, the more convergent the series.
	bool testFlag = true;
	int i;
	// printf("Weighted test results for w5:\n %f, %f, %f, %f\n", test[0], test[1], test[2], test[3]);
	// printf("Parameters for w5:\n %f, %f, %f, %f\n", std::real(a), std::real(b), std::real(c), std::real(z));

	// find series with smallest weighted test result that is less than 1
	// (the less than one condition ensures that we choose a convergent series)
	for(int j = 0; j < 4; j++){
		if( test[j] > 0 ){
			testFlag = false;
		}
	}

	for(int j = 0; j < 4; j++){
		if( test[j] <= 0 ){
			if(testFlag){
				test[j] = 1 - test[j];
			}else{
				test[j] = abs(test[j]);
			}
		}
	}

	testFlag = false;
	// find series with smallest weighted test result that is less than 1
	// (the less than one condition ensures that we choose a convergent series)
	while(testNum < 4 && testFlag == false){
		i = 1;
		testTemp = test[0];
		test[0] = test[testNum];
		test[testNum] = testTemp;
		testFlag = (test[0] <= test[i]) && (test[0] < testLimit);
		while(i < 4 && testFlag == true){
			testFlag = (test[0] <= test[i]) && (test[0] < testLimit);
			i++;
		}
		testNum++;
	}
	// note that the test number is actually one greater than array number

	// if none of our tests estimate a numerically stable or convergent series,
	// return 0, which tells us that none of tests were successful
	if( testFlag == false ){
		testNum = 0;
	}

	return testNum;
}

int hypergeo_2F1_test_w6(Complex a, Complex b, Complex c, Complex z, double){
	int testNum = 0;
	double test[4];
	// test numerical stability of different series representations of w3
	// if 0 then series does not converge, if -1 then series rapidly oscillates
	// in sign, leading to large cancellations, and if +1 then the series may
	// be relatively stable.
	test[0] =	hypergeo_2F1_test_w6_1overz_1(a, b, c, z)*abs(1./z);
	test[1] =	hypergeo_2F1_test_w6_1overz_2(a, b, c, z)*abs(1./z);
	if(abs(1./(1. - z)) > 1.){
		test[2] =	0.;
		test[3] = 0.;
	}else{
		test[2] =	hypergeo_2F1_test_w6_1over1mz_1(a, b, c, z)*abs(1./(1. - z));
		test[3] =	hypergeo_2F1_test_w6_1over1mz_2(a, b, c, z)*abs(1./(1. - z));
	}
	// we weight each test with the magnitude of its argument. Typically the
	// smaller the argument, the more convergent the series.
	bool testFlag = true;
	// printf("Weighted test results for w6:\n %f, %f, %f, %f\n", test[0], test[1], test[2], test[3]);
	// printf("Parameters for w6:\n %f, %f, %f, %f\n", std::real(a), std::real(b), std::real(c), std::real(z));

	// find series with smallest weighted test result that is less than 1
	// (the less than one condition ensures that we choose a convergent series)
	for(int j = 0; j < 4; j++){
		if( test[j] > 0 ){
			testFlag = false;
		}
	}

	for(int j = 0; j < 4; j++){
		if( test[j] <= 0 ){
			if(testFlag){
				test[j] = 1 - test[j];
			}else{
				test[j] = abs(test[j]);
			}
		}
	}

	testFlag = false;
	// note that the test number is actually one greater than array number

	// if none of our tests estimate a numerically stable or convergent series,
	// return 0, which tells us that none of tests were successful
	if( testFlag == false ){
		testNum = 0;
	}

	return testNum;
}

// test_hypergeo_2F1 estimates which representations of the Kummer solutions
// w1, w3, w4, w5, and w6 (see DLMF Sec. 15.10) are the most numerically
// stable. It then returns a test_hypergeo_2F1_struct with this information.
test_hypergeo_2F1 hypergeo_2F1_test_w(Complex a, Complex b, Complex c, Complex z, double testLimit){
	test_hypergeo_2F1 test = test_hypergeo_2F1_default;
	test.w1 = hypergeo_2F1_test_w1(a, b, c, z, testLimit);
	test.w3 = hypergeo_2F1_test_w3(a, b, c, z, testLimit);
	test.w4 = hypergeo_2F1_test_w4(a, b, c, z, testLimit);
	test.w5 = hypergeo_2F1_test_w5(a, b, c, z, testLimit);
	test.w6 = hypergeo_2F1_test_w6(a, b, c, z, testLimit);

	return test;
}

// hypergeo_2F1_test takes the input test_hypergeo_2F1_struct, which contains
// estimates of the most numerically stable representations of the Kummer
// solutions w1, w3, w4, w5, w6 (see DLMF Sec. 15.10) and from these
// constructs w1 numerically either from one of the 4 direct representations
// of w1 (DLMF 15.10.11) or from one of the 6 connection formulas (DLMF
// 15.10.21, 15.10.25, 15.10.29, 15.10.30, 15.10.33, and 15.10.34) which
// relate w1 to a pair of the other other Kummer solutions, w3, w4, w5, w6.

// NOTE: We ignore w2, because w1 and w2 form independent solutions of the
// hypergeometric equation near one of the regular singular points, and
// therefore cannot be represented in terms of one another.
int hypergeo_2F1_test(test_hypergeo_2F1 test){
	if(test.w1 > 0){
		return 1; // use DLMF 15.10.11
	}else if(test.w3 + test.w4 + test.w5 + test.w6 == 0){
		return 0; // all tests failed
	}else if(test.w3 > 0 && test.w4 > 0){
		return 2; // use DLMF 15.10.21
	}else if(test.w5 > 0 && test.w6 > 0){
		return 3; // use DLMF 15.10.25
	}else if(test.w3 > 0 && test.w5 > 0){
		return 4; // use DLMF 15.10.29
	}else if(test.w3 > 0 && test.w6 > 0){
		return 5; // use DLMF 15.10.30
	}else if(test.w4 > 0 && test.w5 > 0){
		return 6; // use DLMF 15.10.33
	}else if(test.w4 > 0 && test.w6 > 0){
		return 7; // use DLMF 15.10.34
	}

	// all tests failed
	return 0;
}

Result hypergeo_2F1(Complex a, Complex b, Complex c, Complex z){
	Result f21(0., DBL_EPSILON);

	if(abs(z) == 0.){
		return Result(1., 0);
	}

	Complex hypergeo;
	if(hypergeo_2F1_special_value(a, b, c, z, hypergeo)){
		return Result(hypergeo, DBL_EPSILON);
	}

	if(z == 1.){
		if(std::real(c - a - b) > 0.){
			f21 += cgamma(c)*cgamma(c - a - b)/cgamma(c - a)/cgamma(c - b);
			return f21;
		}else{
			f21 += INFINITY;
			return f21;
		}
	}

	test_hypergeo_2F1 test = hypergeo_2F1_test_w(a, b, c, z);
	// printf("Test results:\nw1 = %d, w2 = %d, w3 = %d, w4 = %d, w5 = %d, w6 = %d\n", test.w1, test.w2, test.w3, test.w4, test.w5, test.w6);
	int testNum = hypergeo_2F1_test(test);
	// if( testNum == 0 ){
	// 	test = hypergeo_2F1_test_w(0, 0, 0, z);
	// 	testNum = hypergeo_2F1_test(test);
	// }
	switch(testNum){
		case 0:
			// printf("(HYPERGEO_F) ERROR: Case 0 in Gauss hypergeometric function\n");
			f21 = hypergeo_2F1_naive(a, b, c, z);
			break;

		case 1:
			f21 = hypergeo_2F1_w1_151011(test, a, b, c, z);
			break;

		case 2:
			f21 = hypergeo_2F1_w1_151021(test, a, b, c, z);
			break;

		case 3:
			f21 = hypergeo_2F1_w1_151025(test, a, b, c, z);
			break;

		case 4:
			f21 = hypergeo_2F1_w1_151029(test, a, b, c, z);
			break;

		case 5:
			f21 = hypergeo_2F1_w1_151030(test, a, b, c, z);
			break;

		case 6:
			f21 = hypergeo_2F1_w1_151033(test, a, b, c, z);
			break;

		case 7:
			f21 = hypergeo_2F1_w1_151034(test, a, b, c, z);
			break;

		default:
			// printf("HYPERGEO_F: ERROR: Argument test failed.\n");
			testNum = 8;
	}

	if(abs(f21.getPrecision()) > 0.1 || isnan(abs(f21.getPrecision()))){
		f21 = pow(1. - z, -b)*hypergeo_2F1_naive(c - a, b, c, z/(z - 1.));
	}

	if(abs(f21.getPrecision()) > 0.1 || isnan(abs(f21.getPrecision()))){
		f21 = pow(1. - z, -a)*hypergeo_2F1_naive(a, c - b, c, z/(z - 1.));
	}

	if(abs(f21.getPrecision()) > 0.1 || isnan(abs(f21.getPrecision()))){
		f21 = hypergeo_2F1_naive(a, b, c, z);
	}

	if(abs(f21.getPrecision()) > 0.1 || isnan(abs(f21.getPrecision()))){
		f21 = pow(1. - z, c - a - b)*hypergeo_2F1_naive(c - a, c - b, c, z);
	}

	// if(abs(f21.getValue()) == 0.){
	// 	std::cout << "(HYPERGEO_F) ERROR: Failed to construct hypergeometric function for (a, b, c, z) = ("<<a<<", "<<b<<", "<<c<<", "<<z<<").\n";
	// }

	// std::cout << "F(a = "<<a<<", b = "<<b<<", c = "<<c<<", z = "<<z<<") = " << f21.getValue() << "\n";

	return f21;
}

Result hypergeo_2F1_naive(Complex a, Complex b, Complex c, Complex z){
	Result f21(0., DBL_EPSILON);

	test_hypergeo_2F1 test = hypergeo_2F1_test_w(a, b, c, z, 1.);
	// std::cout << " z = " << z << "\n";
	// printf("Test results:\n w1 = %d, w2 = %d, w3 = %d, w4 = %d, w5 = %d, w6 = %d\n", test.w1, test.w2, test.w3, test.w4, test.w5, test.w6);
	int testNum = hypergeo_2F1_test(test);
	// printf("testNum = %d\n",testNum);
	switch(testNum){
		case 1:
			f21 = hypergeo_2F1_w1_151011(test, a, b, c, z);
			break;

		case 2:
			f21 = hypergeo_2F1_w1_151021(test, a, b, c, z);
			break;

		case 3:
			f21 = hypergeo_2F1_w1_151025(test, a, b, c, z);
			break;

		case 4:
			f21 = hypergeo_2F1_w1_151029(test, a, b, c, z);
			break;

		case 5:
			f21 = hypergeo_2F1_w1_151030(test, a, b, c, z);
			break;

		case 6:
			f21 = hypergeo_2F1_w1_151033(test, a, b, c, z);
			break;

		case 7:
			f21 = hypergeo_2F1_w1_151034(test, a, b, c, z);
			break;

		default:
			// printf("(HYPERGEO_F) ERROR: Argument test failed.\n");
			testNum = 0;
	}

	return f21;
}

Result dhypergeo_2F1(Complex a, Complex b, Complex c, Complex z){
	Result dF(0., 0.);
	if(abs(c) != 0.){
		dF = a*b/c*hypergeo_2F1(a + 1., b + 1., c + 1., z); // DLMF 15.5.1
	}else{
		return Result(0., 0.);
	}

	// if calculation appears to have failed try new relation DLMF 15.5.21
	if(abs(dF.getAccuracy()) == 0.){
		dF = ((c - a)*(c - b)*hypergeo_2F1(a, b, c + 1., z) + c*(a + b - c)*hypergeo_2F1(a, b, c, z))/c/(1. - z);
	}

	// if calculation appears to have failed try new relation 15.5.20
	if(abs(dF.getAccuracy()) == 0.){
		dF = ((c - a)*hypergeo_2F1(a - 1., b, c, z) + c*(a - c + b*z)*hypergeo_2F1(a, b, c, z))/z/(1. - z);
	}

	return dF;
}

void test_hypergeo_2F1_series(){

	Complex a = 100.234 + I*2.45;
	Complex b = -100.524 + I*0.34;
	Complex c = 3.2 + 10.23*I;
	Complex z = 21.9;

	Result f21 = hypergeo_2F1(a, b, c, z)/hypergeo_2F1(a, b, c, z+1.);

	std::cout << f21 << "\n";
}
