// header bessel.hpp

#ifndef BESSEL_HPP
#define BESSEL_HPP

#include "specialfunc.hpp"
#include "hypergeo_f.hpp"
#include "gsl/gsl_sf_bessel.h"

Complex bessel_I(Complex nu, Complex z);
Complex bessel_J(Complex nu, Complex z);
Complex bessel_I_bessel_series(Complex nu, double x);
Complex bessel_J_bessel_series(Complex nu, double x);

Complex hypergeo_1f1_bessel_series(Complex a, Complex b, Complex z);
Result hypergeo_1f1_bessel_series_2(Complex a, Complex b, Complex z);
Complex hypergeo_1f1_bessel_series_3(Complex a, Complex b, Complex z);

Result hypergeo_u_bessel_series_2(Complex a, Complex b, Complex z);
Complex hypergeo_u_bessel_series_2_term(int n, Complex a, Complex b, Complex z);
Complex hypergeo_u_bessel_series_2_prefactor(Complex a, Complex b, Complex z);

#endif