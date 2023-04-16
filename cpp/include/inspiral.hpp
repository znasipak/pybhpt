// inspiral.hpp

#ifndef INSPIRAL_HPP
#define INSPIRAL_HPP

#include <chrono>
#include <thread>
#include "boost/filesystem.hpp"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include "fluxes.hpp"
#include "waveform.hpp"

typedef struct AdiabaticModeStruct{
	double spin;
	double chi;
	double freqISCO;
	Vector r0;
	Vector freq;
	Vector alpha;
	Vector A;
	Vector Phi;
	Vector EdotUp;
	Vector EdotIn;
} AdiabaticMode;

typedef struct TrajectoryDataStruct{
	Vector x;
	Vector y;
	Vector z;
} TrajectoryData;

typedef struct SplineDataStruct{
	gsl_spline* spline;
	gsl_interp_accel* acc;
	double a;
	double rescale;
} SplineData;

double alpha_of_a_omega(double a, double omega);
double chi_of_a(double a);

double omega_of_a_alpha(double a, double omega);
double radius_of_a_alpha(double a, double omega);

Vector bh_spin_sample_chi(double amax, int sampleNum);
Vector radius_sample_a_alpha(double a, double alphaMin, double alphaMax, int sampleNum);

int generate_adiabatic_circ_data_a_sample_parallel(int lMax, const std::string& dir);
int generate_adiabatic_circ_mode_data_a_sample(int L, int m, int sampleSize, const std::string& dir);
int generate_adiabatic_circ_remainder_data_a_sample(int Lmin, int sampleSize, const std::string& dir);

int generate_adiabatic_circ_data_radial_sample_parallel(double a, int Lmax, int sampleNum, const std::string& dir);
int generate_adiabatic_circ_mode_data_radial_sample(double a, int L, int m, int sampleNum, const std::string& dir);
int generate_adiabatic_circ_remainder_data_radial_sample(double a, int Lmin, int sampleNum, const std::string& dir);

int generate_adiabatic_circ_mode_data_a_r0(double a, double r0, int L, int m, std::ofstream& file);
int generate_adiabatic_circ_remainder_data_a_r0(double a, double r0, int Lmin, std::ofstream& file);

int output_adiabatic_circ_spherical_mode_data_radial_sample(double a, int lmax, int sampleNum, const std::string& dir);
int output_adiabatic_circ_spherical_mode_data_radial_sample_parallel(double a, int lmax, int sampleNum, const std::string& dir);
int output_adiabatic_circ_spherical_mode_data_radial_sample(double a, int lmax, int m, int sampleNum, const std::string& dir);
int generate_adiabatic_circ_spherical_mode_data_radial_sample(ComplexMatrix &amplitude, ComplexMatrix &flux, double a, int m);
int generate_adiabatic_circ_spherical_mode_data_a_r0(ComplexVector &amplitude, ComplexVector &flux, double a, double r0, int m);

int generate_adiabatic_circ_data(double a, int sgnX, int Lmax, const std::string& dir);
int generate_adiabatic_circ_mode_data(double a, int sgnX, int L, int m, const std::string& dir);
int generate_adiabatic_circ_mode_data_r0(double a, int sgnX, double risco, double r0, int L, int m, std::ofstream& file);
Vector radial_frequency_sample(double a, double rmin, double rmax, double sgnX, int sampleNum = 101);
Vector radial_log_sample(double risco);
double log_rescale(double r0, double r);
double exp_rescale(double r0, double x);

AdiabaticMode read_adiabatic_circ_mode(int L, int m, const std::string& dir);
// int generate_adiabatic_circ_trajectory(double a, const std::string& dir);
int generate_adiabatic_circ_trajectory(const std::string& dir, int outputSampleSize = 500);
void compile_and_output_trajectory_data(std::vector<std::string> inputDirs, const std::string& outputDir);
void compile_and_output_mode_amplitude_data(std::vector<std::string> inputDirs, const std::string& outputDir, int lmax, int radialSampleN);

TrajectoryData read_trajectory_data(const std::string& filename);

double energy_flux_of_omega_spline(double omega, SplineData spline);

void generate_adiabatic_inspiral_data_2d(std::string mainDir, int sampleSize_spin, int sampleSize_radius);
void generate_adiabatic_inspiral_data_1d(std::string dir, double spin, int sampleSize_radius);

int integrate_first_order_real_ode(Vector &psi, const Vector &x, int (*sys)(double, const double*, double*, void*), double psi0, const double x0, void *params);

#endif
