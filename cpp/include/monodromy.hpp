// monodromy.h

#ifndef MONODROMY_HPP
#define MONODROMY_HPP

#include "specialfunc.hpp"

#define SERIES_COEFF_NMAX 2000

// declare structure
typedef struct series_coeff_struct{
	int nmax;
	int ndata;
	ComplexVector coeffs;
} series_coeff;

typedef struct confluent_heun_parameters_struct{
	Complex aCH;
	Complex ggCH;
	Complex dCH;
	Complex eCH;
	Complex qCH;
} CH_parameters;

// declare functions

//////////////////////////
// Monodromy eigenvalue //
//////////////////////////
Complex nu_solver_monodromy(int s, int l, int m, double q, double eps, double la);
Complex monodromy_eigenvalue(series_coeff &a1, series_coeff &a2, const CH_parameters &params);

////////////////////////////////////
// Stokes multiplier coefficients //
////////////////////////////////////
series_coeff series_coeff_init(int nmax);

int generate_weighted_a1(series_coeff &a1, const CH_parameters &params);
int generate_weighted_a2(series_coeff &a2, const CH_parameters &params);
int weighted_renormalize(series_coeff &a, Complex weight);

Complex sum_weighted_a1(const series_coeff &a);
Complex sum_weighted_a1(const series_coeff &a, const int &n);
Complex sum_weighted_a2(const series_coeff &a);
Complex sum_weighted_a2(const series_coeff &a, const int &n);
Complex sum_weighted_a1_drop(const series_coeff &a, const int &n);
Complex sum_weighted_a2_drop(const series_coeff &a, const int &n);

//////////////////
// IO functions //
//////////////////
void cnumprint(Complex cnum);
void cnumfullprint(Complex cnum);
char* complex_to_string(char* buf, Complex cnum);

void test_monodromy(void);

#endif
