// hypergeo_u.hpp

#ifndef HYPERGEOU_HPP
#define HYPERGEOU_HPP

#include "specialfunc.hpp"
#include "bessel.hpp"
#include "cf.hpp"

enum Domain {Region0, Region1, Region2, Region3};

typedef struct U_parameters_struct{
	Complex a;
	Complex b;
	Complex z;
} U_parameters;

Result hyper_incomplete_gamma(Complex a, Complex z);
Result hyper_incomplete_Gamma(Complex a, Complex z);
Result hyper_exponential_integral(Complex a, Complex z);
Result hyper_erf(Complex z);
Result hyper_erfc(Complex z);
Result hyper_bessel_J(Complex nu, Complex z);
Result hyper_bessel_I(Complex nu, Complex z);
Result hyper_bessel_K(Complex nu, Complex z);
Result hyper_Ai(Complex z);

Result log_hyper_incomplete_gamma(Complex a, Complex z);
Result log_hyper_incomplete_Gamma(Complex a, Complex z);
Result log_hyper_exponential_integral(Complex a, Complex z);
Result log_hyper_erf(Complex z);
Result log_hyper_erfc(Complex z);
Result log_hyper_bessel_J(Complex nu, Complex z);
Result log_hyper_bessel_I(Complex nu, Complex z);
Result log_hyper_bessel_K(Complex nu, Complex z);
Result log_hyper_Ai(Complex z);

Result hyper_1f1(Complex a, Complex b, Complex z);
Result hyper_u(Complex a, Complex b, Complex z);
Result hyper_0f1(Complex b, Complex z);

Result dhyper_1f1(Complex a, Complex b, Complex z);
Result dhyper_u(Complex a, Complex b, Complex z);

Result log_hyper_1f1(Complex a, Complex b, Complex z);
Result log_hyper_u(Complex a, Complex b, Complex z);
Result log_hyper_0f1(Complex b, Complex z);

Result log_dhyper_1f1(Complex a, Complex b, Complex z);
Result log_dhyper_u(Complex a, Complex b, Complex z);

/** \brief Irregular confluent hypergeometric function U(a,b,z)
	* \param a Complex parameter
	* \param b Complex parameter
	* \param z Complex argument
	* \return Result for the value of the U
  */
Result hypergeo_u(Complex a, Complex b, Complex z);
Result log_hypergeo_u(Complex a, Complex b, Complex z);

/** \brief Derivative of the irregular confluent hypergeometric function U(a,b,z)
	* \param a Complex parameter
	* \param b Complex parameter
	* \param z Complex argument
	* \return Result for the value of the derivative of U
  */
Result dhypergeo_u(Complex a, Complex b, Complex z);
Result log_dhypergeo_u(Complex a, Complex b, Complex z);

/** \brief Regular confluent hypergeometric function M(a,b,z)
	* \param a Complex parameter
	* \param b Complex parameter
	* \param z Complex argument
	* \return Result for the value of the M
  */
Result hypergeo_m(Complex a, Complex b, Complex z);
Result log_hypergeo_m(Complex a, Complex b, Complex z);

/** \brief Derivative of the regular confluent hypergeometric function M(a,b,z)
	* \param a Complex parameter
	* \param b Complex parameter
	* \param z Complex argument
	* \return Result for the value of the derivative of M
  */
Result dhypergeo_m(Complex a, Complex b, Complex z);
Result log_dhypergeo_m(Complex a, Complex b, Complex z);

Result log_hypergeo_u_negative_z(Complex a, Complex b, Complex z);
Result hypergeo_m_connection(Complex a, Complex b, Complex z);
Result hypergeo_u_connection(Complex a, Complex b, Complex z);
Result log_hypergeo_m_connection(Complex a, Complex b, Complex z);
Result log_hypergeo_u_connection(Complex a, Complex b, Complex z);
Result log_hypergeo_u_connection_recurrence(Complex a, Complex b, Complex z);

Result hypergeo_u_recurrence(Complex a, Complex b, Complex z);
Result hypergeo_u_series(Complex a, Complex b, Complex z);
Result log_hypergeo_u_recurrence(Complex a, Complex b, Complex z);
Result hypergo_u_ap1bp1(Result u, Result u_m1, Complex a, Complex b, Complex z);
Result hypergo_u_ap1(Result u, Result u_m1, Complex a, Complex b, Complex z);
Result hypergo_u_bp1(Result u, Result u_m1, Complex a, Complex b, Complex z);
Result hypergo_u_am1bm1(Result u, Result u_p1, Complex a, Complex b, Complex z);
Result hypergo_u_am1(Result u, Result u_p1, Complex a, Complex b, Complex z);
Result hypergo_u_bm1(Result u, Result u_p1, Complex a, Complex b, Complex z);
Result log_hypergeo_u_series(Complex a, Complex b, Complex z);

Result log_hypergeo_m_recurrence(Complex a, Complex b, Complex z);

Result hypergeo_m_taylor_series(Complex a, Complex b, Complex z);
Result hypergeo_u_asymptotic_series(Complex a, Complex b, Complex z);
Complex hypergeo_u_asymptotic_series_sum(Complex a, Complex b, Complex z);

double error_bound_factor(int n, Complex a, Complex b, Complex z);
double error_bound_factor(Complex a, Complex b, Complex z);

Domain error_bound_region(Complex a, Complex b, Complex z);
double error_bound_Cn(int n, Complex a, Complex b, Complex z);
double error_bound_Cn_1(int n, Complex a, Complex b, Complex z);
double error_bound_Cn_2(int n, Complex a, Complex b, Complex z);
double error_bound_Cn_3(int n, Complex a, Complex b, Complex z);

double error_bound_alpha(Complex a, Complex b, Complex z);
double error_bound_rho(Complex a, Complex b, Complex z);
double error_bound_sigma(Complex a, Complex b, Complex z);
double error_bound_nu(Complex a, Complex b, Complex z);

Result hypergo_u_a_over_u_ap1(Complex a, Complex b, Complex z);
Result hypergo_u_bp1_over_u_b(Complex a, Complex b, Complex z);
Result hypergo_u_ab_over_u_ap1bp1(Complex a, Complex b, Complex z);

Result hypergo_m_b_over_m_bp1(Complex a, Complex b, Complex z);
Result hypergo_m_ab_over_m_ap1bp1(Complex a, Complex b, Complex z);


#endif