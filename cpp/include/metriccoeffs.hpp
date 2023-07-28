// metric.hpp

#include "utils.hpp"
#include "geo.hpp"

#ifndef METRIC_COEFFS_HPP
#define METRIC_COEFFS_HPP

Complex metric_coefficient_ORG(int ai, int bi, int nt, int nr, int nz, int np, double a, double r, double z);
Complex metric_coefficient_IRG(int ai, int bi, int nt, int nr, int nz, int np, double a, double r, double z);
Complex metric_coefficient_ORG_11(int nt, int nr, int nz, int np, double a, double r, double z);
Complex metric_coefficient_ORG_13(int nt, int nr, int nz, int np, double a, double r, double z);
Complex metric_coefficient_ORG_33(int nt, int nr, int nz, int np, double a, double r, double z);
Complex metric_coefficient_IRG_22(int nt, int nr, int nz, int np, double a, double r, double z);
Complex metric_coefficient_IRG_24(int nt, int nr, int nz, int np, double a, double r, double z);
Complex metric_coefficient_IRG_44(int nt, int nr, int nz, int np, double a, double r, double z);

ComplexTensor metric_coefficients_ORG_11(double a, Vector r, Vector z);
ComplexTensor metric_coefficients_ORG_13(double a, Vector r, Vector z);
ComplexTensor metric_coefficients_ORG_33(double a, Vector r, Vector z);

ComplexTensor metric_coefficients_ORG_11_dz(double a, Vector r, Vector z);
ComplexTensor metric_coefficients_ORG_13_dz(double a, Vector r, Vector z);
ComplexTensor metric_coefficients_ORG_33_dz(double a, Vector r, Vector z);

ComplexTensor metric_coefficients_ORG_11_dz2(double a, Vector r, Vector z);
ComplexTensor metric_coefficients_ORG_13_dz2(double a, Vector r, Vector z);
ComplexTensor metric_coefficients_ORG_33_dz2(double a, Vector r, Vector z);

ComplexTensor metric_coefficients_IRG_22(double a, Vector r, Vector z);
ComplexTensor metric_coefficients_IRG_24(double a, Vector r, Vector z);
ComplexTensor metric_coefficients_IRG_44(double a, Vector r, Vector z);

ComplexTensor metric_coefficients_IRG_22_dz(double a, Vector r, Vector z);
ComplexTensor metric_coefficients_IRG_24_dz(double a, Vector r, Vector z);
ComplexTensor metric_coefficients_IRG_44_dz(double a, Vector r, Vector z);

ComplexTensor metric_coefficients_IRG_22_dz2(double a, Vector r, Vector z);
ComplexTensor metric_coefficients_IRG_24_dz2(double a, Vector r, Vector z);
ComplexTensor metric_coefficients_IRG_44_dz2(double a, Vector r, Vector z);

ComplexMatrix tetrad_velocity_1(double a, double En, double Lz, double Qc, Vector r, Vector z);
ComplexMatrix tetrad_velocity_1_dz(double a, double En, double Lz, double Qc, Vector r, Vector z);
ComplexMatrix tetrad_velocity_1_dz2(double a, double En, double Lz, double Qc, Vector r, Vector z);

ComplexMatrix tetrad_velocity_2(double a, double En, double Lz, double Qc, Vector r, Vector z);
ComplexMatrix tetrad_velocity_2_dz(double a, double En, double Lz, double Qc, Vector r, Vector z);
ComplexMatrix tetrad_velocity_2_dz2(double a, double En, double Lz, double Qc, Vector r, Vector z);

ComplexMatrix tetrad_velocity_3(double a, double En, double Lz, double Qc, Vector r, Vector z);
ComplexMatrix tetrad_velocity_3_dz(double a, double En, double Lz, double Qc, Vector r, Vector z);
ComplexMatrix tetrad_velocity_3_dz2(double a, double En, double Lz, double Qc, Vector r, Vector z);

ComplexMatrix tetrad_velocity_4(double a, double En, double Lz, double Qc, Vector r, Vector z);
ComplexMatrix tetrad_velocity_4_dz(double a, double En, double Lz, double Qc, Vector r, Vector z);
ComplexMatrix tetrad_velocity_4_dz2(double a, double En, double Lz, double Qc, Vector r, Vector z);

#endif //METRIC_COEFFS_HPP