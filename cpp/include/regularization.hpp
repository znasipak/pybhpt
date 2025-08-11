// regularization.hpp

#ifndef REG_HPP
#define REG_HPP

#include "teukolsky.hpp"
#include "kerr.hpp"

Vector regularization_parameter_from_source(GeodesicSource &geo, const int &sampleSize, double (*reg_func)(const double &, const double &, const int &, const int &, GeodesicSource &));
Vector regularization_parameter_from_source_circular(GeodesicSource &geo, const int &sampleSize, double (*reg_func)(const double &, const double &, const int &, const int &, GeodesicSource &));
Vector regularization_parameter_from_source_equatorial(GeodesicSource &geo, const int &sampleSize, double (*reg_func)(const double &, const double &, const int &, const int &, GeodesicSource &));
Vector regularization_parameter_from_source_spherical(GeodesicSource &geo, const int &sampleSize, double (*reg_func)(const double &, const double &, const int &, const int &, GeodesicSource &));
Vector regularization_parameter_from_source_generic(GeodesicSource &geo, const int &sampleSize, double (*reg_func)(const double &, const double &, const int &, const int &, GeodesicSource &));

Vector regularization_parameter_At_from_source(GeodesicSource &geo, const int &sampleSize);
Vector regularization_parameter_Ar_from_source(GeodesicSource &geo, const int &sampleSize);
double regularization_parameter_At(const double &rp, const double &thp, const int &sgnUr, const int &sgnUth, GeodesicSource &geo);
double regularization_parameter_Ar(const double &rp, const double &thp, const int &sgnUr, const int &sgnUth, GeodesicSource &geo);

double regularization_parameter_At_low_level(const double &a, const double &rp, const double &thp, const double &ut, const double &ur, const double &uth, const double &uph);
double regularization_parameter_Ar_low_level(const double &a, const double &rp, const double &thp, const double &ut, const double &ur, const double &uth, const double &uph);

double regularization_parameter_IK(const double &alpha, const double &beta, const int &N);
double regularization_parameter_IE(const double &alpha, const double &beta, const int &N);
double regularization_parameter_IN_norm(int N, const double &alpha, const double &beta);
int phi_index_count(int a, int b, int c, int d);
double regularization_P_alpha_beta(int alpha, int beta, const double &a, const double &rp, const double &thp, const double &ut, const double &ur, const double &uth, const double &uph);
double regularization_P_alpha_beta_gamma(int alpha, int beta, int gamma, const double &a, const double &rp, const double &thp, const double &ut, const double &ur, const double &uth, const double &uph);
double regularization_C(int i, int j, int k, const double &thp);
double regularization_P_mu_abcd(int mu, int a, int b, int c, int d, const double &spin, const double &rp, const double &thp, const double &ut, const double &ur, const double &uth, const double &uph);
double regularization_parameter_I_abcd(int a, int b, int c, int d, const double &spin, const double &rp, const double &thp, const double &ut, const double &ur, const double &uth, const double &uph);

Vector regularization_parameter_Bt_from_source(GeodesicSource &geo, const int &sampleSize);
Vector regularization_parameter_Br_from_source(GeodesicSource &geo, const int &sampleSize);
Vector regularization_parameter_Btheta_from_source(GeodesicSource &geo, const int &sampleSize);
Vector regularization_parameter_Bphi_from_source(GeodesicSource &geo, const int &sampleSize);
double regularization_parameter_Bt(const double &rp, const double &thp, const int &sgnUr, const int &sgnUth, GeodesicSource &geo);
double regularization_parameter_Br(const double &rp, const double &thp, const int &sgnUr, const int &sgnUth, GeodesicSource &geo);
double regularization_parameter_Btheta(const double &rp, const double &thp, const int &sgnUr, const int &sgnUth, GeodesicSource &geo);
double regularization_parameter_Bphi(const double &rp, const double &thp, const int &sgnUr, const int &sgnUth, GeodesicSource &geo);


double regularization_parameter_B_alpha(int alpha, const double &rp, const double &thp, const int &sgnUr, const int &sgnUth, GeodesicSource &geo);
double regularization_parameter_B_alpha(int alpha, const double &spin, const double &rp, const double &thp, const double &ut, const double &ur, const double &uth, const double &uph);

#endif
