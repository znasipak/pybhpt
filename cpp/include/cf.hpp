// Continued fraction functions

#ifndef CF_HPP
#define CF_HPP

#include <gsl/gsl_math.h>
#include "utils.hpp"

///////////////////////////////////////////////////
// define constants, structures, and other types //
///////////////////////////////////////////////////

// typedef for a pointer cf_func_n that points to a function of type
// Complex that takes in an int and void pointer for arguments
typedef Complex (*cf_func_n)(int, void*);

// data structure cf_coeffs_struct is used to store the continued fraction
// coefficients a_coeff and b_coeff, which are functions of an
// integer n and depend on the parameters, with the void pointer params
// pointing to those additional parameters.
typedef struct cf_coeffs_struct{
	cf_func_n a_coeff;
	cf_func_n b_coeff;
	void* params;
} cf_coeffs;

// same typedef and structs as above, but now for continued fraction
// coefficients that can be a function of some complex argument
typedef Complex (*cf_func_n_x)(int, Complex, void*);

typedef struct cf_coeffs_x_struct{
	cf_func_n_x a_coeff;
	cf_func_n_x b_coeff;
	void* params;
} cf_coeffs_x;

// data structure cf_solver_struct contains all of the relevant information
// to implement the continued fraction solver
typedef struct cf_solver_struct{
	cf_func_n a_coeff;
	cf_func_n b_coeff;
	void* params;
	Complex Cn;
	Complex Dn;
	Complex convergent_val;
	int convergent_num;
	double cf_error;
} cf_solver;

typedef struct cf_solver_x_struct{
	cf_func_n_x a_coeff;
	cf_func_n_x b_coeff;
	void* params;
	Complex Cn;
	Complex Dn;
	Complex convergent_val;
	int convergent_num;
	double cf_error;
} cf_solver_x;

// a test structure for testing the continued fraction solver
struct test_params{
	double a;
	double b;
};

///////////////////////
// declare functions //
///////////////////////

// cf_solver_alloc allocates memory for the cf_solver structure. It then
// stores data from the cf_ab structure in a new cf_solver structure cf,
// then returns a pointer to cf
cf_solver* cf_solver_alloc(cf_coeffs cf_ab);

cf_solver_x* cf_solver_x_alloc(cf_coeffs_x cf_ab);

// cf_solver_free frees up the memory allocated to the cf_solver cf.
void cf_solver_free(cf_solver* cf);
void cf_solver_x_free(cf_solver_x* cf);

// cf_solver_convergent_num returns the current iteration that the cf_solver
// cf has computed in the continued fraction.
int cf_solver_convergent_num(cf_solver* cf);
int cf_solver_x_convergent_num(cf_solver_x* cf);

// cf_solver_convergent returns the convergent_num-th convergent 
// that the cf_solver cf has computed in the continued fraction.
Complex cf_solver_convergent(cf_solver* cf);
Complex cf_solver_x_convergent(cf_solver_x* cf);

// cf_solver_convergent_num returns the fractional error between the
// convergent_num-th convergent and the (convergent_num-1)-th convergent
// computed by the cf_solver cf.
double cf_solver_error(cf_solver* cf);
double cf_solver_x_error(cf_solver_x* cf);

// cf_lentz_iterate uses the continued fraction data pointed to by cf
// to compute one additional iteration of the modified Lentz method for
// computing a continued fraction. The function returns 0 if the error
// resulting from the additional iteration is greater than tol, and returns
// 1 if the error is less than or equal to tol.
int cf_lentz_iterate(cf_solver* cf, double tol);
int cf_lentz_x_iterate(cf_solver_x* cf, Complex x, double tol);

////////////////////
// test functions //
////////////////////

Complex test_a_coeff(int n, void* p);
Complex test_b_coeff(int n, void* p);
int test_cf_solver();

Complex test_a_coeff_x(int n, Complex x, void* p);
Complex test_b_coeff_x(int n, Complex x, void* p);
int test_cf_solver_x();

#endif
