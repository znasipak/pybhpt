// test.hpp

#ifndef TEST_HPP
#define TEST_HPP


#include "boost/filesystem.hpp"
#include "inspiral.hpp"
#include "ssf.hpp"
#include "omp.h"

void test_special_functions();
void test_teukolsky_mode(int s, int L, int m, int k, int n, GeodesicSource geo);
void test_adiabatic_data_generation();

void full_flux_parallel_l(int s, GeodesicSource geo, int modeMax, std::string dir);
void full_flux_parallel_lm(GeodesicSource geo, int lMax, std::string dir);
void flux_parallel_lm(GeodesicSource geo, int lMax = 16);
void flux_parallel_l(int s, GeodesicSource geo, int modeMax = 16);

void ssf_parallel_m(GeodesicSource geo, int lmax);
void ssf_components_parallel_m(GeodesicSource geo, int lmax, int mmin, int mmax, int sampleSSF, const std::string& dir);
void ssf_components_parallel_m(GeodesicSource geo, int lmax, int sampleSSF, const std::string& dir);
void ssf_components_m(GeodesicSource geo, int lmax, int m, int sampleSSF, const std::string& dir);

void generate_adiabatic_inspiral_data(double a, int parallelFlag = 0);

void test_teukolsky_ZlmUp_PN(const std::string& dir);

void test_teukolsky_data_generation(const std::string& dir);
void test_teukolsky_data_mode_generation(double a, int s, int l, int m, double omega, Vector rPts, const std::string& dir);
void test_teukolsky_data_mode_generation_2(double a, int s, int l, int m, double omega, Vector rPts, const std::string& dir);

#endif
