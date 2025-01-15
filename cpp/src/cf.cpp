// cf.cpp

#include "cf.hpp"

///////////////////////////////////////////////////
// define constants, structures, and other types //
///////////////////////////////////////////////////

// define a small constant to be used in the modified Lentz continued fraction
// algorithm

#define LENTZ_SMALL 1.e-50

static cf_solver cf_solver_default = {
	.Cn = LENTZ_SMALL,
	.Dn = 0.,
	.convergent_val = LENTZ_SMALL,
	.convergent_num = -1,
	.cf_error = 0.
};

///////////////////////
// declare functions //
///////////////////////

// cf_solver_alloc allocates memory for the cf_solver structure. It then
// stores data from the cf_ab structure in a new cf_solver structure cf,
// then returns a pointer to cf
cf_solver* cf_solver_alloc(cf_coeffs cf_ab){
	cf_solver* cf;
	cf = (cf_solver *)malloc(sizeof(cf_solver_default));
	cf->a_coeff = (cf_ab.a_coeff);
	cf->b_coeff = (cf_ab.b_coeff);
	cf->params = (cf_ab.params);

	cf->Cn = cf_solver_default.Cn;
	cf->Dn = cf_solver_default.Dn;
	cf->convergent_val = cf_solver_default.convergent_val;
	cf->convergent_num = cf_solver_default.convergent_num;
	cf->cf_error = cf_solver_default.cf_error;

	return cf;
}

cf_solver_x* cf_solver_x_alloc(cf_coeffs_x cf_ab){
	cf_solver_x* cf;
	cf = (cf_solver_x *)malloc(sizeof(cf_solver_default));
	cf->a_coeff = (cf_ab.a_coeff);
	cf->b_coeff = (cf_ab.b_coeff);
	cf->params = (cf_ab.params);

	cf->Cn = cf_solver_default.Cn;
	cf->Dn = cf_solver_default.Cn;
	cf->convergent_val = cf_solver_default.convergent_val;
	cf->convergent_num = cf_solver_default.convergent_num;
	cf->cf_error = cf_solver_default.cf_error;

	return cf;
}

// cf_solver_free frees up the memory allocated to the cf_solver cf.
void cf_solver_free(cf_solver* cf){
	free(cf);
}

void cf_solver_x_free(cf_solver_x* cf){
	free(cf);
}

// cf_solver_convergent_num returns the current iteration that the cf_solver
// cf has computed in the continued fraction.
int cf_solver_convergent_num(cf_solver* cf){
	return cf->convergent_num;
}

int cf_solver_x_convergent_num(cf_solver_x* cf){
	return cf->convergent_num;
}

// cf_solver_convergent returns the convergent_num-th convergent
// that the cf_solver cf has computed in the continued fraction.
Complex cf_solver_convergent(cf_solver* cf){
	return cf->convergent_val;
}

Complex cf_solver_x_convergent(cf_solver_x* cf){
	return cf->convergent_val;
}

// cf_solver_convergent_num returns the fractional error between the
// convergent_num-th convergent and the (convergent_num-1)-th convergent
// computed by the cf_solver cf.
double cf_solver_error(cf_solver* cf){
	return cf->cf_error;
}

double cf_solver_x_error(cf_solver_x* cf){
	return cf->cf_error;
}

// cf_lentz_iterate uses the continued fraction data pointed to by cf
// to compute one additional iteration of the modified Lentz method for
// computing a continued fraction. The function returns 0 if the error
// resulting from the additional iteration is greater than tol, and returns
// 1 if the error is less than or equal to tol.
int cf_lentz_iterate(cf_solver* cf, double tol){
	Complex resn, resnm1;
	Complex Dnm1, Dn, Cnm1, Cn;
	cf_func_n* a = NULL;
	cf_func_n* b = NULL;
	void* params = cf->params;
	bool convergence_flag = 0;
	double cf_error;
	int modifyFlag = 0;

	Dnm1 = cf->Dn;
	Cnm1 = cf->Cn;

	resnm1 = cf->convergent_val;
	int j = cf->convergent_num;
	j++;

	a = &(cf->a_coeff);
	Complex an = (*a)(j, params);
	b = &(cf->b_coeff);
	Complex bn = (*b)(j, params);

	Dn = bn + an*Dnm1;
	Cn = bn + an/Cnm1;

	double precisionLoss = std::abs(1. + bn*Cnm1/an);
	precisionLoss = precisionLoss < std::abs(1. + bn/an/Dnm1) ? precisionLoss : std::abs(1. + bn/an/Dnm1);

	if(precisionLoss <= 1.e-5 && j > 2){
		//std::cout << "CF: Significant precision loss of " << precisionLoss << " encountered at n = "<< j;
		//std::cout << ". Shifting to next convergent \n";
		j++;

		Complex Cnm2 = Cnm1;
		Complex Dnm2 = Dnm1;
		Complex anm1 = an;
		Complex bnm1 = bn;
		an = (*a)(j, params);
		bn = (*b)(j, params);

		Cnm1 = (bnm1*bnm1*Cnm2*Cnm2 - anm1*anm1)/Cnm2/(bnm1*Cnm2 - anm1);
		Dnm1 = (bnm1 - anm1*Dnm2)/(bnm1*bnm1 - anm1*anm1*Dnm2*Dnm2);

		/*
		if(std::abs(1. + bn*Cnm1/an) <= 1.e-5 || std::abs(1. + bn/an/Dnm1) <= 1.e-5){
			std::cout << "CF: Significant precision loss still encountered after skipping to next convergent. \n";
		}else{
			std::cout << "CF: Significant precision avoided by skipping to next convergent. \n";
		}
		*/

		Dn = bn + an*Dnm1;
		Cn = bn + an/Cnm1;
	}

	if(std::abs(Cn) == 0. || std::abs(1. + bn*Cnm1/an) <= 1.e-14){
		Cn = LENTZ_SMALL;
		modifyFlag = 1;
	}
	if(std::abs(Dn) == 0. || std::abs(1. + bn/an/Dnm1) <= 1.e-14){
		Dn = 1./LENTZ_SMALL;
		modifyFlag = 1;
	}else{
		Dn = 1./Dn;
	}

	resn = resnm1*Dn*Cn;
	cf_error = std::abs(1. - resnm1/resn);

	cf->Cn = Cn;
	cf->Dn = Dn;
	cf->convergent_val = resn;
	cf->convergent_num = j;
	cf->cf_error = cf_error;

	if(cf_error > tol || modifyFlag == 1){
		convergence_flag = 0;
	}else{
		convergence_flag = 1;
	}

	return convergence_flag;
}

int cf_lentz_x_iterate(cf_solver_x* cf, Complex x, double tol){
	Complex resn, resnm1;
	Complex Dnm1, Dn, Cnm1, Cn;
	cf_func_n_x* a = NULL;
	cf_func_n_x* b = NULL;
	void* params = cf->params;
	bool convergence_flag = 0;
	double cf_error;
	int modifyFlag = 0;

	Dnm1 = cf->Dn;
	Cnm1 = cf->Cn;

	resnm1 = cf->convergent_val;
	int j = cf->convergent_num;
	j++;

	a = &(cf->a_coeff);
	Complex an = (*a)(j, x, params);
	b = &(cf->b_coeff);
	Complex bn = (*b)(j, x, params);

	double SMALL = std::abs(bn)*LENTZ_SMALL;
	if( SMALL == 0. ){
		SMALL = LENTZ_SMALL;
	}

	Dn = bn + an*Dnm1;
	Cn = bn + an/Cnm1;

	if(std::abs(Cn) == 0.){
		Cn = SMALL;
		modifyFlag = 1;
	}
	if(std::abs(Dn) == 0.){
		Dn = 1./SMALL;
		modifyFlag = 1;
	}else{
		Dn = 1./Dn;
	}

	resn = resnm1*Dn*Cn;

	cf_error = std::abs(1. - resnm1/resn);

	cf->Cn = Cn;
	cf->Dn = Dn;
	cf->convergent_val = resn;
	cf->convergent_num = j;
	cf->cf_error = cf_error;

	if(cf_error > tol || modifyFlag == 1){
		convergence_flag = 0;
	}else{
		convergence_flag = 1;
	}

	return convergence_flag;
}

////////////////////
// test functions //
////////////////////

Complex test_a_coeff(int n, void* p){
	struct test_params *params = (struct test_params *)p;
	double a = params->a;
	double b = params->b;

	return a*Complex(n) - b;
}

Complex test_b_coeff(int, void*){
	return 2.;
}

int test_cf_solver(){

	double x = 0.941352;
	double cfx = std::real(sqrt(x) - 1.);

	cf_coeffs CF;
	struct test_params params = {0.0, 1. - x};
	bool success = 0;

	CF.a_coeff = &test_a_coeff;
	CF.b_coeff = &test_b_coeff;
	CF.params = &params;

	cf_solver* cf = cf_solver_alloc(CF);

	int iter = 0;
	while(iter < 100 && success == 0){
		iter++;
		success = cf_lentz_iterate(cf, DBL_EPSILON);
	}

	printf("jmax = %d\n", cf_solver_convergent_num(cf));
	printf("error = %g\n", cf_solver_error(cf));
	printf("cf = %.16e\n", std::real(cf_solver_convergent(cf)));
	printf("1 - (sqrt(x)-1)/cf = %.16e\n",
		std::real(1. - cfx/cf_solver_convergent(cf)));

	cf_solver_free(cf);

	return success;
}

//

struct test_x_params{
	double a;
	double b;
};

Complex test_a_coeff_x(int, Complex x, void* p){
	struct test_x_params *params = (struct test_x_params *)p;
	double a = params->a;

	return x + a;
}

Complex test_b_coeff_x(int, Complex, void* p){
	struct test_x_params *params = (struct test_x_params *)p;
	double b = params->b;

	return b;
}

int test_cf_solver_x(){

	Complex x = 5.941352;
	double cfx = std::real(sqrt(x) - 1.);

	cf_coeffs_x CF;
	struct test_x_params params = {-1., 2.};
	bool success = 0;

	CF.a_coeff = &test_a_coeff_x;
	CF.b_coeff = &test_b_coeff_x;
	CF.params = &params;

	cf_solver_x* cf = cf_solver_x_alloc(CF);

	int iter = 0;
	while(iter < 100 && success == 0){
		iter++;
		success = cf_lentz_x_iterate(cf, x, DBL_EPSILON);
	}

	printf("jmax = %d\n", cf_solver_x_convergent_num(cf));
	printf("error = %g\n", cf_solver_x_error(cf));
	printf("cf_x = %.16e\n", std::real(cf_solver_x_convergent(cf)));
	printf("1 - (sqrt(x)-1)/cf = %.16e\n",
		std::real(1. - cfx/cf_solver_x_convergent(cf)));

	cf_solver_x_free(cf);

	return success;
}
