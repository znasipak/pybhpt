// hypergeo_f.hpp

#ifndef HYPERGEO_HPP
#define HYPERGEO_HPP

#include "specialfunc.hpp"
#include <gsl/gsl_sf_hyperg.h>

///////////////////////////////////////////////////
// define constants, structures, and other types //
///////////////////////////////////////////////////

typedef struct solution_hypergeo_2F1_struct{
	Complex val;
	double err;
} solution_hypergeo_2F1;

typedef struct test_hypergeo_2F1_struct{
	int w1;
	int w2;
	int w3;
	int w4;
	int w5;
	int w6;
} test_hypergeo_2F1;

///////////////////////////////////////
// Gauss hypergeometric function 2F1 //
///////////////////////////////////////

double hypergeo_2F1_real(double a, double b, double c, double x);
Complex hypergeo_2F1_group1to3(Complex a, Complex b, Complex x);
Complex hypergeo_2F1_group2to3(Complex a, Complex c, Complex x);

int hypergeo_2F1_special_value(Complex a, Complex b, Complex c, Complex z, Complex &hypergeo);

Result hypergeo_2F1(Complex a, Complex b, Complex c, Complex z);
Complex hypergeo_2F1_complex(Complex a, Complex b, Complex c, Complex z);
Result dhypergeo_2F1(Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_series(solution_hypergeo_2F1* sol2F1, Complex a, Complex b, Complex c, Complex z);

//////////////////
//DLMF 15.10.11 //
//////////////////
Result hypergeo_2F1_w1(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z);

int hypergeo_2F1_w1_z(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_w1_z_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_w1_zoverzm1_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_w1_zoverzm1_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);

int hypergeo_2F1_test_w1_z_1(Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_test_w1_z_2(Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_test_w1_zoverzm1_1(Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_test_w1_zoverzm1_2(Complex a, Complex b, Complex c, Complex z);

double hypergeo_2F1_test_w1_z_1_weighted(Complex a, Complex b, Complex c, Complex z);
double hypergeo_2F1_test_w1_z_2_weighted(Complex a, Complex b, Complex c, Complex z);
double hypergeo_2F1_test_w1_zoverzm1_1_weighted(Complex a, Complex b, Complex c, Complex z);
double hypergeo_2F1_test_w1_zoverzm1_2_weighted(Complex a, Complex b, Complex c, Complex z);

//////////////////
//DLMF 15.10.12 //
//////////////////
Result hypergeo_2F1_w2(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z);

int hypergeo_2F1_w2_z_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_w2_z_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_w2_zoverzm1_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_w2_zoverzm1_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);

int hypergeo_2F1_test_w2_z_1(Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_test_w2_z_2(Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_test_w2_zoverzm1_1(Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_test_w2_zoverzm1_2(Complex a, Complex b, Complex c, Complex z);

double hypergeo_2F1_test_w2_z_1_weighted(Complex a, Complex b, Complex c, Complex z);
double hypergeo_2F1_test_w2_z_2_weighted(Complex a, Complex b, Complex c, Complex z);
double hypergeo_2F1_test_w2_zoverzm1_1_weighted(Complex a, Complex b, Complex c, Complex z);
double hypergeo_2F1_test_w2_zoverzm1_2_weighted(Complex a, Complex b, Complex c, Complex z);

//////////////////
//DLMF 15.10.13 //
//////////////////
Result hypergeo_2F1_w3(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z);

int hypergeo_2F1_w3_1mz_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_w3_1mz_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_w3_1m1overz_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_w3_1m1overz_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);

int hypergeo_2F1_test_w3_1mz_1(Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_test_w3_1mz_2(Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_test_w3_1m1overz_1(Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_test_w3_1m1overz_2(Complex a, Complex b, Complex c, Complex z);

double hypergeo_2F1_test_w3_1mz_1_weighted(Complex a, Complex b, Complex c, Complex z);
double hypergeo_2F1_test_w3_1mz_2_weighted(Complex a, Complex b, Complex c, Complex z);
double hypergeo_2F1_test_w3_1m1overz_1_weighted(Complex a, Complex b, Complex c, Complex z);
double hypergeo_2F1_test_w3_1m1overz_2_weighted(Complex a, Complex b, Complex c, Complex z);

//////////////////
//DLMF 15.10.14 //
//////////////////
Result hypergeo_2F1_w4(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z);

int hypergeo_2F1_w4_1mz_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_w4_1mz_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_w4_1m1overz_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_w4_1m1overz_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);

int hypergeo_2F1_test_w4_1mz_1(Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_test_w4_1mz_2(Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_test_w4_1m1overz_1(Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_test_w4_1m1overz_2(Complex a, Complex b, Complex c, Complex z);

double hypergeo_2F1_test_w4_1mz_1_weighted(Complex a, Complex b, Complex c, Complex z);
double hypergeo_2F1_test_w4_1mz_2_weighted(Complex a, Complex b, Complex c, Complex z);
double hypergeo_2F1_test_w4_1m1overz_1_weighted(Complex a, Complex b, Complex c, Complex z);
double hypergeo_2F1_test_w4_1m1overz_2_weighted(Complex a, Complex b, Complex c, Complex z);

//////////////////
//DLMF 15.10.15 //
//////////////////
Result hypergeo_2F1_w5(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z);

int hypergeo_2F1_w5_1overz_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_w5_1overz_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_w5_1over1mz_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_w5_1over1mz_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);

int hypergeo_2F1_test_w5_1overz_1(Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_test_w5_1overz_2(Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_test_w5_1over1mz_1(Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_test_w5_1over1mz_2(Complex a, Complex b, Complex c, Complex z);

double hypergeo_2F1_test_w5_1overz_1_weighted(Complex a, Complex b, Complex c, Complex z);
double hypergeo_2F1_test_w5_1overz_2_weighted(Complex a, Complex b, Complex c, Complex z);
double hypergeo_2F1_test_w5_1over1mz_1_weighted(Complex a, Complex b, Complex c, Complex z);
double hypergeo_2F1_test_w5_1over1mz_2_weighted(Complex a, Complex b, Complex c, Complex z);

//////////////////
//DLMF 15.10.16 //
//////////////////
Result hypergeo_2F1_w6(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z);

int hypergeo_2F1_w6_1overz_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_w6_1overz_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_w6_1over1mz_1(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_w6_1over1mz_2(solution_hypergeo_2F1* f21, Complex a, Complex b, Complex c, Complex z);

int hypergeo_2F1_test_w6_1overz_1(Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_test_w6_1overz_2(Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_test_w6_1over1mz_1(Complex a, Complex b, Complex c, Complex z);
int hypergeo_2F1_test_w6_1over1mz_2(Complex a, Complex b, Complex c, Complex z);

double hypergeo_2F1_test_w6_1overz_1_weighted(Complex a, Complex b, Complex c, Complex z);
double hypergeo_2F1_test_w6_1overz_2_weighted(Complex a, Complex b, Complex c, Complex z);
double hypergeo_2F1_test_w6_1over1mz_1_weighted(Complex a, Complex b, Complex c, Complex z);
double hypergeo_2F1_test_w6_1over1mz_2_weighted(Complex a, Complex b, Complex c, Complex z);

///////////////////////////
// DLMF Equations for w1 //
///////////////////////////

// Eqn 15.10.11
Result hypergeo_2F1_w1_151011(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z);
// Eqn 15.10.21
Result hypergeo_2F1_w1_151021(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z);
// Eqn 15.10.25
Result hypergeo_2F1_w1_151025(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z);
// Eqn 15.10.29
Result hypergeo_2F1_w1_151029(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z);
// Eqn 15.10.30
Result hypergeo_2F1_w1_151030(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z);
// Eqn 15.10.33
Result hypergeo_2F1_w1_151033(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z);
// Eqn 15.10.34
Result hypergeo_2F1_w1_151034(test_hypergeo_2F1 test, Complex a, Complex b, Complex c, Complex z);

///////////
// Tests //
///////////

int hypergeo_2F1_arg_test(Complex z);
int hypergeo_2F1_test(test_hypergeo_2F1 test);

test_hypergeo_2F1 hypergeo_2F1_test_w(Complex a, Complex b, Complex c, Complex z, double testLimit = 0.7);
int hypergeo_2F1_test_w1(Complex a, Complex b, Complex c, Complex z, double testLimit = 0.7);
int hypergeo_2F1_test_w3(Complex a, Complex b, Complex c, Complex z, double testLimit = 0.7);
int hypergeo_2F1_test_w4(Complex a, Complex b, Complex c, Complex z, double testLimit = 0.7);
int hypergeo_2F1_test_w5(Complex a, Complex b, Complex c, Complex z, double testLimit = 0.7);
int hypergeo_2F1_test_w6(Complex a, Complex b, Complex c, Complex z, double testLimit = 0.7);

int hypergeo_2F1_test_sign(Complex a, Complex b, Complex c, Complex z);
double hypergeo_2F1_test_weighted(Complex a, Complex b, Complex c, Complex z);

void test_hypergeo_2F1_series(void);
Result hypergeo_2F1_naive(Complex a, Complex b, Complex c, Complex z);

////////////////////
// Gamma function //
////////////////////

Complex hypergeo_2F1_cgamma_ratio(Complex num1, Complex num2, Complex den1, Complex den2);


#endif
