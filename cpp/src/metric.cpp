// metric.cpp

#include "metric.hpp"

#define COUPLING_COEFF_NUM 6
#define HUU_EXPANSION_ORDER 0

SphericalHarmonicCoupling::SphericalHarmonicCoupling(int lmax, int m): _lmax(lmax), _lmin(std::abs(m)), _m(m), _zCoupling(COUPLING_COEFF_NUM*(lmax - _lmin + 1)), _dzCoupling(COUPLING_COEFF_NUM*(lmax - _lmin + 1)) {}
SphericalHarmonicCoupling::~SphericalHarmonicCoupling() {}

void SphericalHarmonicCoupling::generateCouplings(){
  generate_z_scalar_spherical_coupling(_zCoupling, _lmin, _m);
  generate_derivative_scalar_spherical_coupling(_dzCoupling, _lmin, _m);
}

int SphericalHarmonicCoupling::getAzimuthalModeNumber(){ return _m; }

Vector SphericalHarmonicCoupling::getZCouplingCoefficient(){ return _zCoupling; }
double SphericalHarmonicCoupling::getZCouplingCoefficient(int n, int i, int l){
	if(l <= getMaxZCouplingModeNumber() && l >= _lmin && std::abs(i) <= n){
		return _zCoupling[6*(l - _lmin) + (n*(n+1))/2 + i];
	}else{
		return 0.;
	}
}
double SphericalHarmonicCoupling::getZCouplingCoefficientNoCheck(int n, int i, int l){
	return _zCoupling[6*(l - _lmin) + (n*(n+1))/2 + i];
}
int SphericalHarmonicCoupling::getMinZCouplingModeNumber(){
	return _lmin;
}
int SphericalHarmonicCoupling::getMaxZCouplingModeNumber(){
	return _zCoupling.size()/6 + _lmin - 1;
}

Vector SphericalHarmonicCoupling::getDerivativeCouplingCoefficient(){ return _dzCoupling; }
double SphericalHarmonicCoupling::getDerivativeCouplingCoefficient(int n, int i, int l){
	if(l <= getMaxDerivativeCouplingModeNumber() && l >= _lmin && std::abs(i) <= n){
		return _dzCoupling[6*(l - _lmin) + (n*(n+1))/2 + i];
	}else{
		return 0.;
	}
}
double SphericalHarmonicCoupling::getDerivativeCouplingCoefficientNoCheck(int n, int i, int l){
	return _dzCoupling[6*(l - _lmin) + (n*(n+1))/2 + i];
}
int SphericalHarmonicCoupling::getMinDerivativeCouplingModeNumber(){
	return _lmin;
}
int SphericalHarmonicCoupling::getMaxDerivativeCouplingModeNumber(){
	return _dzCoupling.size()/6 + _lmin - 1;
}

int test_spherical_coupling(int l, int m, double z, int Nz, int Ndx){
  if(Nz > 2 || Nz < 0 || Ndx > 2 || Ndx < 0){
    std::cout << "Invalid test parameters\n";
    return 0;
  }
  double theta = acos(z);
  SphericalHarmonicCoupling Cljm(l+3, m);
  Cljm.generateCouplings();
  double zNYlm = pow(z, Nz)*Ylm(l, m, theta);
  double dNdxYlm = 0.;
  if(Ndx == 0){
    dNdxYlm = Ylm(l, m, theta);
  }else if(Ndx == 1){
    dNdxYlm = -pow(1. - z*z, 0.5)*Ylm_derivative(l, m, theta);
  }else if(Ndx == 2){
    dNdxYlm = (m*m - (1. - z*z)*l*(l + 1.))*Ylm(l, m, theta);
  }

  double testZCoupling = 0.;
  for(int i = 0; i <= Nz; i++){
    testZCoupling += Cljm.getZCouplingCoefficient(Nz, i, l)*Ylm(l + 2*i - Nz, m, theta);
  }

  double testDXCoupling = 0.;
  for(int i = 0; i <= Ndx; i++){
    testDXCoupling += Cljm.getDerivativeCouplingCoefficient(Ndx, i, l)*Ylm(l + 2*i - Ndx, m, theta);
  }

  std::cout << "Testing z^" << Nz << " Y_{"<<l<<","<<m<<"} \n";
  std::cout << "Actual value = " << zNYlm << "\n";
  std::cout << "Test value = " << testZCoupling << "\n";
  std::cout << "Fractional error = " << std::abs(1. - testZCoupling/zNYlm) << "\n\n";

  std::cout << "Testing partial_x^" << Ndx << " Y_{"<<l<<","<<m<<"} \n";
  std::cout << "Actual value = " << dNdxYlm << "\n";
  std::cout << "Test value = " << testDXCoupling << "\n";
  std::cout << "Fractional error = " << std::abs(1. - testDXCoupling/dNdxYlm) << "\n\n";

  return 0.;
}


int coupling_ni_to_iter(int n, int i){
  return (n*(n+1) + n + i)/2;
}

// Coupling between scalar spherical harmonics, z-weighting, and their derivatives for reprojecting onto Ylm basis
void generate_derivative_scalar_spherical_coupling(Vector &Aljm, int lmin, int m){
  int lnum = Aljm.size()/COUPLING_COEFF_NUM;
  int l;
  for(int lp = 0; lp < lnum; lp++){
    l = lp + lmin;
    Aljm[6*lp + 0] = 1.;
    Aljm[6*lp + 1] = (l + 1.)*clm(l, m);
    Aljm[6*lp + 2] = -l*clm(l + 1, m);
    Aljm[6*lp + 3] = l*(l + 1.)*clm(l, m)*clm(l - 1, m);
    Aljm[6*lp + 4] = -((l - 1.)*(l + 1.)*pow(clm(l, m), 2) + l*(l + 2.)*pow(clm(l + 1, m), 2));
    Aljm[6*lp + 5] = l*(l + 1.)*clm(l + 1, m)*clm(l + 2, m);
  }
}

void generate_z_scalar_spherical_coupling(Vector &Dljm, int lmin, int m){
  int lnum = Dljm.size()/COUPLING_COEFF_NUM;
  int l;
  for(int lp = 0; lp < lnum; lp++){
    l = lp + lmin;
    Dljm[6*lp + 0] = 1.;
    Dljm[6*lp + 1] = clm(l, m);
    Dljm[6*lp + 2] = clm(l + 1, m);
    Dljm[6*lp + 3] = clm(l, m)*clm(l - 1, m);
    Dljm[6*lp + 4] = pow(clm(l, m), 2) + pow(clm(l + 1, m), 2);
    Dljm[6*lp + 5] = clm(l + 1, m)*clm(l + 2, m);
  }
}

Vector generate_dz_spherical_harmonic_coupling(int l, int m){
  Vector Alm(COUPLING_COEFF_NUM);
	// coupling with derivatives (nth derivative of [(1-z^2)\partial_z]^n Ylm, coupling to j = l+i Yjm)
	Alm[0] = 1.; // (n = 0, i = 0)
	Alm[1] = (l + 1.)*clm(l, m); // (n = 1, i = 0)
	Alm[2] = -l*clm(l + 1, m); // (n = 1, i = 1);
	Alm[3] = l*(l + 1.)*clm(l, m)*clm(l - 1, m); // (n = 2, i = 0);
	Alm[4] = -((l - 1.)*(l + 1.)*pow(clm(l, m), 2) + l*(l + 2.)*pow(clm(l + 1, m), 2)); // (n = 2, i = 1);
	Alm[5] = l*(l + 1.)*clm(l + 1, m)*clm(l + 2, m); // (n = 2, i = 2);
  return Alm;
}

Vector generate_z_spherical_harmonic_coupling(int l, int m){
  Vector Dlm(COUPLING_COEFF_NUM);
	// coupling with z-weight (z^n Ylm coupling to j = l+i Yjm)
	Dlm[0] = 1.; // (n = 0, i = 0)
	Dlm[1] = clm(l, m); // (n = 1, i = 0)
	Dlm[2] = clm(l + 1, m); // (n = 1, i = 1);
	Dlm[3] = clm(l, m)*clm(l - 1, m); // (n = 2, i = 0);
	Dlm[4] = pow(clm(l, m), 2) + pow(clm(l + 1, m), 2); // (n = 2, i = 1);
	Dlm[5] = clm(l + 1, m)*clm(l + 2, m); // (n = 2, i = 2);
  return Dlm;
}

double dz_spherical_harmonic_coupling(int n, int i, int l, int m){
  int iter = coupling_ni_to_iter(n, i);
  if(iter == 0){
    return 1.;
  }else if(iter == 1){
    return (l + 1.)*clm(l, m);
  }else if(iter == 2){
    return -l*clm(l + 1, m);
  }else if(iter == 3){
    return l*(l + 1.)*clm(l, m)*clm(l - 1, m);
  }else if(iter == 4){
    return -((l - 1.)*(l + 1.)*pow(clm(l, m), 2) + l*(l + 2.)*pow(clm(l + 1, m), 2));
  }else{
    return l*(l + 1.)*clm(l + 1, m)*clm(l + 2, m);
  }
}
double z_spherical_harmonic_coupling(int n, int i, int l, int m){
  int iter = coupling_ni_to_iter(n, i);
  if(iter == 0){
    return 1.;
  }else if(iter == 1){
    return clm(l, m);
  }else if(iter == 2){
    return clm(l + 1, m);
  }else if(iter == 3){
    return clm(l, m)*clm(l - 1, m);
  }else if(iter == 4){
    return pow(clm(l, m), 2) + pow(clm(l + 1, m), 2);
  }else{
    return clm(l + 2, m);
  }
}

Complex tetrad_covector(int b, int mu, double a, double r, double z){
  double sigma = r*r + a*a*z*z;
  if(b == 1){
    double delta = r*r - 2.*r + a*a;
    if(mu == 0){ // l_t
      return -1;
    }else if(mu == 1){ // l_r
      return sigma/delta;
    }else if(mu == 3){ // l_phi
      return a*(1. - z*z);
    }else{
      return 0.;
    }
  }else if(b == 2){
    double delta = r*r - 2.*r + a*a;
    if(mu == 0){ // n_t
      return -0.5*delta/sigma;
    }else if(mu == 1){ // n_r
      return -0.5;
    }else if(mu == 3){ // n_phi
      return 0.5*a*delta*(1. - z*z)/sigma;
    }else{
      return 0.;
    }
  }else if(b == 3){
    Complex rhobar = -1./(r + I*a*z);
    if(mu == 0){ // m_t
      return I*a*sqrt(0.5*(1. - z*z))*rhobar;
    }else if(mu == 2){ // m_z
      return sigma*sqrt(0.5/(1. - z*z))*rhobar;
    }else if(mu == 3){ // m_phi
      return -I*(r*r + a*a)*sqrt(0.5*(1. - z*z))*rhobar;
    }else{
      return 0.;
    }
  }else if(b == 4){ //mbar
    return std::conj(tetrad_covector(3, mu, a, r, z));
  }

  return 0.;
}

Complex tetrad_vector(int b, int mu, double a, double r, double z){
  if(b == 1){
    double delta = r*r - 2.*r + a*a;
    if(mu == 0){ // l^t
      return (r*r + a*a)/delta;
    }else if(mu == 1){ // l^r
      return 1.;
    }else if(mu == 3){ // l^phi
      return a/delta;
    }else{
      return 0.;
    }
  }else if(b == 2){
    double delta = r*r - 2.*r + a*a;
    double sigma = r*r + a*a*z*z;
    if(mu == 0){ // n^t
      return 0.5*(r*r + a*a)/sigma;
    }else if(mu == 1){ // n^r
      return -0.5*delta/sigma;
    }else if(mu == 3){ // n^phi
      return 0.5*a/sigma;
    }else{
      return 0.;
    }
  }else if(b == 3){
    Complex rhobar = -1./(r + I*a*z);
    if(mu == 0){ // m^t
      return -I*a*sqrt(0.5*(1. - z*z))*rhobar;
    }else if(mu == 2){ // m^z
      return sqrt(0.5*(1. - z*z))*rhobar;
    }else if(mu == 3){ // m^phi
      return -I*sqrt(0.5/(1. - z*z))*rhobar;
    }else{
      return 0.;
    }
  }else if(b == 4){ //mbar
    return std::conj(tetrad_vector(3, mu, a, r, z));
  }

  return 0.;
}

int metric_component_to_iter(int mu, int nu){
  if(mu > nu){
    return metric_component_to_iter(nu, mu);
  }
  if(mu == 0){
    return mu + nu;
  }else if(mu == 1){
    return 3 + nu;
  }else if(mu == 2){
    return 5 + nu;
  }else if(mu == 3){
    return 9;
  }

  return 0.;
}

MetricPerturbation::MetricPerturbation(HertzMode Phi, GeodesicSource &geo): _Phi(Phi), _geo(geo) {}
MetricPerturbation::~MetricPerturbation() {}

Gauge MetricPerturbation::getGauge(){ return _Phi.getGauge(); }
int MetricPerturbation::getSpheroidalModeNumber(){ return _Phi.getSpheroidalModeNumber(); }
int MetricPerturbation::getAzimuthalModeNumber(){ return _Phi.getAzimuthalModeNumber(); }
double MetricPerturbation::getFrequency(){ return _Phi.getFrequency(); }
double MetricPerturbation::getBlackHoleSpin(){ return _Phi.getBlackHoleSpin(); }
double MetricPerturbation::getEigenvalue(){ return _Phi.getEigenvalue(); }
Complex MetricPerturbation::getHertzAmplitude(BoundaryCondition bc){ return _Phi.getHertzAmplitude(bc); }
HertzMode MetricPerturbation::getHertzMode(){ return _Phi; }

Vector MetricPerturbation::getMetricComponent(int mu, int nu, Vector t, Vector r, Vector z, Vector phi){
  int component1 = 2, component2 = 4;
  if(_Phi.getGauge() == ORG || _Phi.getGauge() == SAAB0 || _Phi.getGauge() == ASAAB0){
    component1 = 1;
    component2 = 3;
  }
  ComplexVector hc1c2 = getMetricTetradComponent(component1, component2, t, r, z, phi);
  Vector hmunu(hc1c2.size(), 0.);
  double a = _geo.getBlackHoleSpin();
  for(size_t i = 0; i < hmunu.size(); i++){
    hmunu[i] += 4.*std::real(tetrad_covector(component1, mu, a, r[i], z[i])*tetrad_covector(component2, nu, a, r[i], z[i])*hc1c2[i]);
  }

  if(mu != 2 && nu != 2){
    ComplexVector hc1c1 = getMetricTetradComponent(component1, component2, t, r, z, phi);
    for(size_t i = 0; i < hmunu.size(); i++){
      hmunu[i] += 2.*std::real(tetrad_covector(component1, mu, a, r[i], z[i])*tetrad_covector(component1, nu, a, r[i], z[i])*hc1c1[i]);
    }
  }

  if(mu != 1 && nu != 1){
    ComplexVector hc2c2 = getMetricTetradComponent(component2, component2, t, r, z, phi);
    for(size_t i = 0; i < hmunu.size(); i++){
      hmunu[i] += 2.*std::real(tetrad_covector(component2, mu, a, r[i], z[i])*tetrad_covector(component2, nu, a, r[i], z[i])*hc2c2[i]);
    }
  }

  return hmunu;
}

RealMatrix MetricPerturbation::getMetricComponent(Vector t, Vector r, Vector z, Vector phi){
  int component1 = 1, component2 = 3;
  int metricComponent1 = 2, metricComponent2 = 4;
  if(_Phi.getGauge() == ORG || _Phi.getGauge() == SAAB0 || _Phi.getGauge() == ASAAB0){
    // std::cout <<"Switching to outgoing radiation gauge\n";
    component1 = 2;
    component2 = 4;
    metricComponent1 = 1;
    metricComponent2 = 3;
  }
  ComplexVector hc1c1 = getMetricTetradComponent(metricComponent1, metricComponent1, t, r, z, phi);
  ComplexVector hc1c2 = getMetricTetradComponent(metricComponent1, metricComponent2, t, r, z, phi);
  ComplexVector hc2c2 = getMetricTetradComponent(metricComponent2, metricComponent2, t, r, z, phi);
  RealMatrix hmunu(10, Vector(hc1c2.size(), 0.));
  double a = _geo.getBlackHoleSpin();
  ComplexVector e1(4);
  ComplexVector e2(4);
  for(size_t i = 0; i < hmunu[0].size(); i++){
    //check that the tetrad covectors are correct
    for(int mu = 0; mu < 4; mu++){
      Complex test = 0.;
      for(int nu = 0; nu < 4; nu++){
        test += kerr_metric_blc_z(mu, nu, a, r[i], z[i])*tetrad_vector(component2, nu, a, r[i], z[i]);
      }
      if(std::abs(tetrad_covector(component2, mu, a, r[i], z[i])-test) > 1.e-10){
        std::cout << "Tetrad vector calculation failed for mu = "<<mu<<" \n";
        std::cout << tetrad_covector(component2, mu, a, r[i], z[i]) << "\n";
        std::cout << test << "\n";
      }
      test = 0.;
      for(int nu = 0; nu < 4; nu++){
        test += kerr_metric_blc_z(mu, nu, a, r[i], z[i])*tetrad_vector(component1, nu, a, r[i], z[i]);
      }
      if(std::abs(tetrad_covector(component1, mu, a, r[i], z[i])-test) > 1.e-10){
        std::cout << "Tetrad vector calculation failed for mu = "<<mu<<" \n";
        std::cout << tetrad_covector(component1, mu, a, r[i], z[i]) << "\n";
        std::cout << test << "\n";
      }
      Complex test2 = 0.;
      for(int nu = 0; nu < 4; nu++){
        for(int ai = 1; ai <= 4; ai++){
          if(ai != metricComponent1){
            test2 += tetrad_covector(component1, nu, a, r[i], z[i])*tetrad_vector(ai, nu, a, r[i], z[i]);
          }
        }
      }
      if(std::abs(test2) > 1.e-10){
        std::cout << "Tetrad vector calculation failed for component1 = "<<component1<<" \n";
        std::cout << test2 << "\n";
      }
      test2 = 0.;
      for(int nu = 0; nu < 4; nu++){
        for(int ai = 1; ai <= 4; ai++){
          if(ai != metricComponent2){
            test2 += tetrad_covector(component2, nu, a, r[i], z[i])*tetrad_vector(ai, nu, a, r[i], z[i]);
          }
        }
      }
      if(std::abs(test2) > 1.e-10){
        std::cout << "Tetrad vector calculation failed for component2 = "<<component2<<" \n";
        std::cout << test2 << "\n";
      }
    }

    for(int mu = 0; mu < 4; mu++){
      e1[mu] = tetrad_covector(component1, mu, a, r[i], z[i]);
      e2[mu] = tetrad_covector(component2, mu, a, r[i], z[i]);
    }

    int iter = 0;
    for(int mu = 0; mu < 4; mu++){
      for(int nu = mu; nu < 4; nu++){
        hmunu[iter][i] = 0.;
        hmunu[iter][i] -= 2.*std::real(e1[mu]*e2[nu]*hc1c2[i]);
        hmunu[iter][i] -= 2.*std::real(e2[mu]*e1[nu]*hc1c2[i]);
        hmunu[iter][i] += 2.*std::real(e1[mu]*e1[nu]*hc1c1[i]);
        hmunu[iter][i] += 2.*std::real(e2[mu]*e2[nu]*hc2c2[i]);
        iter++;
      }
    }

    // Tests
    // if(i==0){
    //   // std::cout << "Test gauge conditions" << "\n";
    //   // for(int mu = 0; mu < 4; mu++){
    //   //   Complex test = 0.;
    //   //   for(int nu = 0; nu < 4; nu++){
    //   //     test += tetrad_vector(component1, nu, a, r[i], z[i])*hmunu[metric_component_to_iter(mu, nu)][i];
    //   //   }
    //   //   std::cout << test << "\n";
    //   // }
    //   // std::cout << "Test properties of null vectors\n";
    //   // Complex test = 0.;
    //   // for(int mu = 0; mu < 4; mu++){
    //   //   test += e1[mu]*tetrad_vector(component1, mu, a, r[i], z[i]);
    //   // }
    //   // std::cout << test << "\n";
    //   // test = 0.;
    //   // for(int mu = 0; mu < 4; mu++){
    //   //   test += e2[mu]*tetrad_vector(component1, mu, a, r[i], z[i]);
    //   // }
    //
    //   std::cout << "Test reprojection onto tetrad" << "\n";
    //   Complex test = 0.;
    //   for(int mu = 0; mu < 4; mu++){
    //     for(int nu = 0; nu < 4; nu ++){
    //       test += tetrad_vector(metricComponent1, mu, a, r[i], z[i])*tetrad_vector(metricComponent1, nu, a, r[i], z[i])*hmunu[metric_component_to_iter(mu, nu)][i];
    //     }
    //   }
    //   std::cout << test << "\n";
    //   std::cout << 2.*std::real(hc1c1[i]) << "\n";
    //   test = 0.;
    //   for(int mu = 0; mu < 4; mu++){
    //     for(int nu = 0; nu < 4; nu ++){
    //       test += tetrad_vector(metricComponent1, mu, a, r[i], z[i])*tetrad_vector(metricComponent2, nu, a, r[i], z[i])*hmunu[metric_component_to_iter(mu, nu)][i];
    //     }
    //   }
    //   std::cout << test << "\n";
    //   std::cout << hc1c2[i] << "\n";
    //   test = 0.;
    //   for(int mu = 0; mu < 4; mu++){
    //     for(int nu = 0; nu < 4; nu ++){
    //       test += tetrad_vector(metricComponent2, mu, a, r[i], z[i])*tetrad_vector(metricComponent2, nu, a, r[i], z[i])*hmunu[metric_component_to_iter(mu, nu)][i];
    //     }
    //   }
    //   std::cout << test << "\n";
    //   std::cout << hc2c2[i] << "\n";
    // }

    // my attempt at using gauge conditions, but I messed up
    // e1[0] = tetrad_vector(component1, 0, a, r[i], z[i]);
    // e1[3] = tetrad_vector(component1, 3, a, r[i], z[i]);
    //
    // // calculate htr
    // hmunu[1][i] += std::real(-(e1[0]*hmunu[0][i] + e1[3]*hmunu[3][i])/e1[1]);
    // // calculate hrr
    // hmunu[4][i] += std::real((e1[0]*e1[0]*hmunu[0][i] + 2.*e1[0]*e1[3]*hmunu[3][i] + e1[3]*e1[3]*hmunu[9][i])/e1[1]/e1[1]);
    // // calculate hrz
    // hmunu[5][i] += std::real(-(e1[0]*hmunu[2][i] + e1[3]*hmunu[8][i])/e1[1]);
    // // calculate hrphi
    // hmunu[6][i] += std::real(-(e1[0]*hmunu[3][i] + e1[3]*hmunu[9][i])/e1[1]);
  }

  return hmunu;
}

ComplexVector MetricPerturbation::getMetricTetradComponent(int ai, int bi, Vector t, Vector r, Vector z, Vector phi){
  double a = _Phi.getBlackHoleSpin();
  int s = _Phi.getSpinWeight();
  int j = _Phi.getSpheroidalModeNumber();
  int m = _Phi.getAzimuthalModeNumber();
  double omega = _Phi.getFrequency();
  Vector theta(z.size());
  for(size_t i = 0; i < theta.size(); i++){
    theta[i] = acos(z[i]);
  }
  RadialTeukolsky Rjmo(a, s, j, m, omega, r);
  Rjmo.generateSolutions(HBL);
  SpinWeightedHarmonic Sjmo(s, j, m, a*omega, theta);
  Sjmo.generateSolutionsAndDerivatives();
  ComplexVector hab(theta.size());
  Complex dPsijmo;
  double dSjmo;
  BoundaryCondition bc = In;
  Gauge gauge = _Phi.getGauge();
  for(size_t i = 0; i < theta.size(); i++){
    for(int nt = 0; nt <= 2; nt++){
      for(int nr = 0; nr <= 2; nr++){
        if(r[i] <= _geo.getRadialPositionOfMinoTime(_geo.getMinoTimeOfTime(t[i]))){
          bc = In;
        }else{
          bc = Up;
        }
        if(nr == 0){
          dPsijmo = _Phi.getHertzAmplitude(bc)*Rjmo.getSolution(bc, i);
        }else if(nr == 1){
          dPsijmo = _Phi.getHertzAmplitude(bc)*Rjmo.getDerivative(bc, i);
        }else{
          dPsijmo = _Phi.getHertzAmplitude(bc)*Rjmo.getSecondDerivative(bc, i);
        }
        for(int nz = 0; nz <= 2; nz++){
          if(nz == 0){
            dSjmo = Sjmo.getSolution(i);
          }else if(nz == 1){ // SWSH are functions of theta but the metric coefficients are organized wrt to d/dz
            dSjmo = -Sjmo.getDerivative(i)/pow(1. - pow(z[i], 2), 0.5);
          }else{
            dSjmo = Sjmo.getSecondDerivative(i)/(1. - pow(z[i], 2)) - z[i]*Sjmo.getDerivative(i)/pow(1. - pow(z[i], 2), 1.5);
          }
          for(int np = 0; np <= 2; np++){
            if(nt + nr + nz + np <= 2){
              if(gauge == ORG || gauge == ASAAB0 || gauge == SAAB0){
                // std::cout << "metric coeff = " << metric_coefficient_ORG(ai, bi, nt, nr, nz, np, a, r[i], z[i]) << "\n";
                // std::cout << dPsijmo << "\n";
                // std::cout << dSjmo << "\n";
                hab[i] += metric_coefficient_ORG(ai, bi, nt, nr, nz, np, a, r[i], z[i])*pow(I*double(m), np)*pow(-I*omega, nt)*dPsijmo*dSjmo*exp(I*(m*phi[i] - omega*t[i]));
              }else{
                hab[i] += metric_coefficient_IRG(ai, bi, nt, nr, nz, np, a, r[i], z[i])*pow(I*double(m), np)*pow(-I*omega, nt)*dPsijmo*dSjmo*exp(I*(m*phi[i] - omega*t[i]));
              }
            }
          }
        }
      }
    }
  }

  return hab;
}

RealMatrix metric_perturbation_circ(Gauge gauge, int l, int m, Vector t, Vector r, Vector z, Vector phi, GeodesicSource &geo){
  TeukolskyMode teuk(-2, l, m, 0, 0, geo);
  teuk.generateSolutions(geo);
  HertzMode hertz(teuk, gauge);
  hertz.generateSolutions();
  MetricPerturbation hab(hertz, geo);
  return hab.getMetricComponent(t, r, z, phi);
}

Vector metric_perturbation_circ(double a, double b, Gauge gauge, int l, int m, Vector t, Vector r, Vector z, Vector phi, GeodesicSource &geo){
  TeukolskyMode teuk(-2, l, m, 0, 0, geo);
  teuk.generateSolutions(geo);
  HertzMode hertz(teuk, gauge);
  hertz.generateSolutions();
  MetricPerturbation hab(hertz, geo);
  return hab.getMetricComponent(a, b, t, r, z, phi);
}

ComplexVector metric_perturbation_tetrad_circ(double a, double b, Gauge gauge, int l, int m, Vector t, Vector r, Vector z, Vector phi, GeodesicSource &geo){
  TeukolskyMode teuk(-2, l, m, 0, 0, geo);
  teuk.generateSolutions(geo);
  HertzMode hertz(teuk, gauge);
  hertz.generateSolutions();
  MetricPerturbation hab(hertz, geo);
  return hab.getMetricTetradComponent(a, b, t, r, z, phi);
}

RealMatrix worldline_grid(double t0, double rmin, double rmax, int sampleNumR, double zmin, double zmax, int sampleNumZ, double phi0, GeodesicSource &geo){
  double r0 = geo.getRadialPositionOfMinoTime(geo.getMinoTimeOfTime(t0));

  int sampleNumHalf = sampleNumR/2;
  RealMatrix xp(4, Vector(2*sampleNumHalf*sampleNumZ, 0.));

  double deltaR1 = (r0 - 1.e-8 - rmin);
  double deltaR2 = (r0 + 1.e-8 - rmax);
  if(sampleNumHalf > 1){
    deltaR1 /= (sampleNumHalf - 1);
    deltaR2 /= (sampleNumHalf - 1);
  }
  double deltaZ = (zmax - zmin);
  if(sampleNumZ > 1){
      deltaZ /= (sampleNumZ - 1);;
  }
  for(int jr = 0; jr < sampleNumHalf; jr++){
    for(int jz = 0; jz < sampleNumZ; jz++){
      xp[0][sampleNumZ*jr + jz] = t0;
      xp[0][sampleNumZ*(2*sampleNumHalf - 1 - jr) + jz] = t0;
      xp[1][sampleNumZ*jr + jz] = rmin + deltaR1*double(jr);
      xp[1][sampleNumZ*(2*sampleNumHalf - 1 - jr) + jz] = rmax + deltaR2*double(jr);
      xp[2][sampleNumZ*jr + jz] = zmin + deltaZ*double(jz);
      xp[2][sampleNumZ*(2*sampleNumHalf - 1 - jr) + jz] = zmin + deltaZ*double(jz);
      xp[3][sampleNumZ*jr + jz] = phi0;
      xp[3][sampleNumZ*(2*sampleNumHalf - 1 - jr) + jz] = phi0;
    }
  }

  return xp;
}
