// metric.hpp

#include "utils.hpp"
#include "hertz.hpp"
#include "metriccoeffs.hpp"
#include "kerr.hpp"

#ifndef METRIC_HPP
#define METRIC_HPP

class SphericalHarmonicCoupling{
public:
	SphericalHarmonicCoupling(int lmax, int m);
	~SphericalHarmonicCoupling();

	void generateCouplings();

	int getAzimuthalModeNumber();

	Vector getZCouplingCoefficient();
	double getZCouplingCoefficient(int n, int i, int l);
  	double getZCouplingCoefficientNoCheck(int n, int i, int l);
	int getMinZCouplingModeNumber();
	int getMaxZCouplingModeNumber();

	Vector getDerivativeCouplingCoefficient();
	double getDerivativeCouplingCoefficient(int n, int i, int l);
  	double getDerivativeCouplingCoefficientNoCheck(int n, int i, int l);
	int getMinDerivativeCouplingModeNumber();
	int getMaxDerivativeCouplingModeNumber();

private:
	int _lmax;
  	int _lmin;
	int _m;

	Vector _zCoupling;
	Vector _dzCoupling;
};

class MetricPerturbation{
public:
	MetricPerturbation(HertzMode Phi, GeodesicSource &geo);
	~MetricPerturbation();

	Gauge getGauge();
	int getSpheroidalModeNumber();
	int getAzimuthalModeNumber();
	double getFrequency();
	double getBlackHoleSpin();
	double getEigenvalue();
	Complex getHertzAmplitude(BoundaryCondition bc);
	HertzMode getHertzMode();

	double getMetricComponent(int mu, int nu, double t, double r, double z, double phi);
	Complex getMetricTetradComponent(int a, int b, double t, double r, double z, double phi);

	Vector getMetricComponent(int mu, int nu, Vector t, Vector r, Vector z, Vector phi);
	ComplexVector getMetricTetradComponent(int a, int b, Vector t, Vector r, Vector z, Vector phi);

	RealMatrix getMetricComponent(Vector t, Vector r, Vector z, Vector phi);

private:
	HertzMode _Phi;
	GeodesicSource _geo;
};

int test_spherical_coupling(int l, int m, double z, int Nz, int Ndx);

int coupling_ni_to_iter(int n, int i);
void generate_derivative_scalar_spherical_coupling(Vector &Aljm, int lmin, int m);
void generate_z_scalar_spherical_coupling(Vector &Dljm, int lmin, int m);
Vector generate_dz_spherical_harmonic_coupling(int l, int m);
Vector generate_z_spherical_harmonic_coupling(int l, int m);

double dz_spherical_harmonic_coupling(int n, int i, int l, int m);
double z_spherical_harmonic_coupling(int n, int i, int l, int m);

RealMatrix worldline_grid(double t0, double rmin, double rmax, int sampleNumR, double zmin, double zmax, int sampleNumZ, double phi0, GeodesicSource &geo);
ComplexVector metric_perturbation_tetrad_circ(double a, double b, Gauge gauge, int l, int m, Vector t, Vector r, Vector z, Vector phi, GeodesicSource &geo);
RealMatrix metric_perturbation_circ(Gauge gauge, int l, int m, Vector t, Vector r, Vector z, Vector phi, GeodesicSource &geo);
Vector metric_perturbation_circ(double a, double b, Gauge gauge, int l, int m, Vector t, Vector r, Vector z, Vector phi, GeodesicSource &geo);

#endif //METRIC_HPP
