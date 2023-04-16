//unit_test.cpp

#ifndef UNIT_HPP
#define UNIT_HPP

#include "fluxes.hpp"

void run_unit_tests();
int test_kerr_circ_eq_flux();
int test_kerr_circ_eq_flux_2();
int test_schwarzschild_circ_eq_flux();
int test_teukolsky_radial_function();
int test_teukolsky_radial_function_highModes();
int test_teukolsky_polar_function();

#endif
