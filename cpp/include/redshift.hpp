// redshift.hpp

#include "metric.hpp"

#ifndef REDSHIFT_HPP
#define REDSHIFT_HPP

// class RedshiftData{
// public:
//     RedshiftData(Gauge gauge, int lmax, GeodesicSource &geo);

// private:
//     GeodesicSource& _geo;
//     int _lmax;
//     Gauge _gauge;
//     ComplexTensor _huuCoeffs;
//     ComplexMatrix _huuIn;
//     ComplexMatrix _huuUp;
// };

class RedshiftCoefficients{
public:
    RedshiftCoefficients(Gauge gauge, GeodesicSource &geo);

    Complex getComponent(int Ni, int ai, int bi, int ci, int di, int jru, int jzu);

private:
    ComplexTensor _coeffs;
};

void export_circular_redshift_data(std::string filename, ComplexVector HuuIn, ComplexVector HuuUp);
void export_circular_redshift_data_lm(std::string filename, ComplexMatrix HuuIn, ComplexMatrix HuuUp);
void redshift_circular(std::string filename, Gauge gauge, int lmax, GeodesicSource &geoCirc);

Complex redshift_mode(ComplexTensor &huuCoeff, SphericalHarmonicCoupling &Cjlm, HertzMode &hertz, BoundaryCondition bc, int l, int jr, int jz, int sgnUr, int sgnUz);
Complex redshift_coefficient_components(ComplexTensor &huu, int Ni, int ai, int bi, int ci, int di, int jr, int jz);
ComplexTensor redshift_coefficients_ORG(GeodesicSource geo, int sampleNum);
ComplexTensor redshift_coefficients_ORG(double a, double En, double Lz, double Qc, Vector r, Vector z);
ComplexTensor redshift_coefficients_IRG(GeodesicSource geo, int sampleNum);
ComplexTensor redshift_coefficients_IRG(double a, double En, double Lz, double Qc, Vector r, Vector z);

int redshift_mode_circular_radial(Complex &radial0, Complex &radial1, HertzMode &hertz, BoundaryCondition bc, GeodesicSource &geo);
Complex redshift_mode_circular_l(HertzMode &hertz, int l, Complex radial0, Complex radial1);
Complex redshift_mode_circular(HertzMode &hertz, BoundaryCondition bc, int l, GeodesicSource &geo);
double redshift_regularization_circular(GeodesicSource &geo);
double redshift_completion_circular(GeodesicSource &geo);

#endif REDSHIFT_HPP
