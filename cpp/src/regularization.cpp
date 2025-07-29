// regularization.cpp

#include "regularization.hpp"

double lower_U_alpha(int alpha, const double &a, const double &rp, const double &thp, const double &ut, const double &ur, const double &uth, const double &uph){
  if(alpha == 0){
    return kerr_metric_blc(0, 0, a, rp, thp)*ut + kerr_metric_blc(0, 3, a, rp, thp)*uph;
  }else if(alpha == 1){
    return kerr_metric_blc(1, 1, a, rp, thp)*ur;
  }else if(alpha == 2){
    return kerr_metric_blc(2, 2, a, rp, thp)*uth;
  }else if(alpha == 3){
    return kerr_metric_blc(0, 3, a, rp, thp)*ut + kerr_metric_blc(3, 3, a, rp, thp)*uph;
  }else{
    return 0.;
  }
}

Vector regularization_parameter_from_source(GeodesicSource &geo, const int &sampleSize, double (*reg_func)(const double &, const double &, const int &, const int &, GeodesicSource &)){
  if(std::abs(geo.getInclination()) == 1. && geo.getEccentricity() == 0.){
    return regularization_parameter_from_source_circular(geo, sampleSize, *reg_func);
  }else if(std::abs(geo.getInclination()) == 1.){
    return regularization_parameter_from_source_equatorial(geo, sampleSize, *reg_func);
  }else if(geo.getEccentricity() == 0.){
    return regularization_parameter_from_source_spherical(geo, sampleSize, *reg_func);
  }
  return regularization_parameter_from_source_generic(geo, sampleSize, *reg_func);
}

Vector regularization_parameter_from_source_circular(GeodesicSource &geo, const int &, double (*reg_func)(const double &, const double &, const int &, const int &, GeodesicSource &)){
  double rp = geo.getRadialPosition()[0];
  double thp = geo.getPolarPosition()[0];
  Vector RP(1);

  RP[0] = reg_func(rp, thp, 1, 1, geo);

  return RP;
}

Vector regularization_parameter_from_source_spherical(GeodesicSource &geo, const int &sampleSize, double (*reg_func)(const double &, const double &, const int &, const int &, GeodesicSource &)){
  double rp = geo.getRadialPosition()[0];
  Vector thp = geo.getPolarPosition();
  Vector RP(sampleSize);

  int polarSampleSize = thp.size() - 1;
  int polarSampleRate = 2*polarSampleSize/sampleSize;

  for(int jth = 0; jth <= sampleSize/2; jth++){
    RP[jth] = reg_func(rp, thp[jth*polarSampleRate], 1, 1, geo);
    if(jth > 0 && jth < sampleSize/2){
      RP[sampleSize - jth] = reg_func(rp, thp[jth*polarSampleRate], 1, -1, geo);
    }
  }

  return RP;
}

Vector regularization_parameter_from_source_equatorial(GeodesicSource &geo, const int &sampleSize, double (*reg_func)(const double &, const double &, const int &, const int &, GeodesicSource &)){
  Vector rp = geo.getRadialPosition();
  double thp = geo.getPolarPosition()[0];
  Vector RP(sampleSize);

  int radialSampleSize = rp.size() - 1;
  int radialSampleRate = 2*radialSampleSize/sampleSize;

  for(int jr = 0; jr <= sampleSize/2; jr++){
    RP[jr] = reg_func(rp[jr*radialSampleRate], thp, 1, 1, geo);
    if(jr > 0 && jr < sampleSize/2){
      RP[(sampleSize - jr)] = reg_func(rp[jr*radialSampleRate], thp, -1, 1, geo);
    }
  }

  return RP;
}

Vector regularization_parameter_from_source_generic(GeodesicSource &geo, const int &sampleSize, double (*reg_func)(const double &, const double &, const int &, const int &, GeodesicSource &)){
  int sampleSizeSquared = sampleSize*sampleSize;
  Vector rp = geo.getRadialPosition();
  Vector thp = geo.getPolarPosition();
  Vector RP(sampleSizeSquared);

  int radialSampleSize = rp.size() - 1;
  int polarSampleSize = thp.size() - 1;
  int radialSampleRate = 2*radialSampleSize/sampleSize;
  int polarSampleRate = 2*polarSampleSize/sampleSize;

  // std::cout << "("<<radialSampleSize<<", "<<polarSampleSize<<") \n";
  // std::cout << "("<<sampleSize<<", "<<sampleSize<<") \n";
  // std::cout << "("<<radialSampleRate<<", "<<polarSampleRate<<") \n";
  for(int jr = 0; jr <= sampleSize/2; jr++){
    for(int jth = 0; jth <= sampleSize/2; jth++){
      // std::cout << "("<<jr<<", "<<jth<<") \n";
      RP[jr*sampleSize + jth] = reg_func(rp[jr*radialSampleRate], thp[jth*polarSampleRate], 1, 1, geo);
      if(jr > 0 && jth > 0 && jr < sampleSize/2 && jth < sampleSize/2){
        RP[(sampleSize - jr)*sampleSize + sampleSize - jth] = reg_func(rp[jr*radialSampleRate], thp[jth*polarSampleRate], -1, -1, geo);
        RP[(sampleSize - jr)*sampleSize + jth] = reg_func(rp[jr*radialSampleRate], thp[jth*polarSampleRate], -1, 1, geo);
        RP[jr*sampleSize + sampleSize - jth] = reg_func(rp[jr*radialSampleRate], thp[jth*polarSampleRate], 1, -1, geo);
      }else if(jr > 0 && jr < sampleSize/2){
        RP[(sampleSize - jr)*sampleSize + jth] = reg_func(rp[jr*radialSampleRate], thp[jth*polarSampleRate], -1, 1, geo);
      }else if(jth > 0 && jth < sampleSize/2){
        RP[jr*sampleSize + sampleSize - jth] = reg_func(rp[jr*radialSampleRate], thp[jth*polarSampleRate], 1, -1, geo);
      }
    }
  }

  return RP;
}

Vector regularization_parameter_At_from_source(GeodesicSource &geo, const int &sampleSize){
  return regularization_parameter_from_source(geo, sampleSize, *regularization_parameter_At);
}
Vector regularization_parameter_Ar_from_source(GeodesicSource &geo, const int &sampleSize){
  return regularization_parameter_from_source(geo, sampleSize, *regularization_parameter_Ar);
}

double regularization_parameter_At(const double &rp, const double &thp, const int &sgnUr, const int &sgnUth, GeodesicSource &geo){
  double a = geo.getBlackHoleSpin();
  double En = geo.getOrbitalEnergy();
  double Lz = geo.getOrbitalAngularMomentum();
  // std::cout << "Lz = " <<Lz << "\n";
  double Qc = geo.getCarterConstant();
  double ut = (kerr_geo_VtR(a, En, Lz, Qc, rp) + kerr_geo_VtTheta(a, En, Lz, Qc, thp))/(rp*rp + pow(a*cos(thp), 2));
  double ur = kerr_geo_Vr(a, En, Lz, Qc,  rp);
  if(ur/pow(rp, 4) < DBL_EPSILON){
    ur = 0.;
  }else{
    ur = sgnUr*sqrt(ur)/(rp*rp + pow(a*cos(thp), 2));
  }
  double uth = kerr_geo_Vtheta(a, En, Lz, Qc, thp);
  if(uth < DBL_EPSILON){
    uth = 0.;
  }else{
    uth = sgnUth*sqrt(uth)/(rp*rp + pow(a*cos(thp), 2));
  }
  double uph = (kerr_geo_VphiR(a, En, Lz, Qc, rp) + kerr_geo_VphiTheta(a, En, Lz, Qc, thp))/(rp*rp + pow(a*cos(thp), 2));

  return regularization_parameter_At_low_level(geo.getBlackHoleSpin(), rp, thp, ut, ur, uth, uph);
}

double regularization_parameter_Ar(const double &rp, const double &thp, const int &sgnUr, const int &sgnUth, GeodesicSource &geo){
  double a = geo.getBlackHoleSpin();
  double En = geo.getOrbitalEnergy();
  double Lz = geo.getOrbitalAngularMomentum();
  double Qc = geo.getCarterConstant();
  double ut = (kerr_geo_VtR(a, En, Lz, Qc, rp) + kerr_geo_VtTheta(a, En, Lz, Qc, thp))/(rp*rp + pow(a*cos(thp), 2));
  double ur = kerr_geo_Vr(a, En, Lz, Qc,  rp);
  if(ur < 0 || ur/pow(rp, 4) < DBL_EPSILON){
    // recall Vr potential goes like r^4. Rescale to see if we are at a radial turning point
    ur = 0.;
  }else{
    ur = sgnUr*sqrt(ur)/(rp*rp + pow(a*cos(thp), 2));
  }
  double uth = kerr_geo_Vtheta(a, En, Lz, Qc, thp);
  if(uth < 0 || uth < DBL_EPSILON){
    uth = 0.;
  }else{
    uth = sgnUth*sqrt(uth)/(rp*rp + pow(a*cos(thp), 2));
  }
  double uph = (kerr_geo_VphiR(a, En, Lz, Qc, rp) + kerr_geo_VphiTheta(a, En, Lz, Qc, thp))/(rp*rp + pow(a*cos(thp), 2));

  return regularization_parameter_Ar_low_level(geo.getBlackHoleSpin(), rp, thp, ut, ur, uth, uph);
}

double regularization_parameter_At_low_level(const double &a, const double &rp, const double &thp, const double &ut, const double &ur, const double &uth, const double &uph){
  return -ur/ut*regularization_parameter_Ar_low_level(a, rp, thp, ut, ur, uth, uph);
}

double regularization_parameter_Ar_low_level(const double &a, const double &rp, const double &thp, const double &ut, const double &ur, const double &uth, const double &uph){
  double grr = kerr_metric_blc(1, 1, a, rp, thp);
  double gthth = kerr_metric_blc(2, 2, a, rp, thp);
  double gphph = kerr_metric_blc(3, 3, a, rp, thp);
  double gtph = kerr_metric_blc(0, 3, a, rp, thp);

  double Ur = grr*ur;
  double Uth = gthth*uth;
  double Uph = gtph*ut + gphph*uph;
  double Vfactor = 1. + pow(Uth, 2)/gthth + pow(Uph, 2)/gphph;

  return -sqrt(pow(sin(thp), 2)*grr/gthth/gphph)*sqrt(Vfactor + pow(Ur, 2)/grr)/Vfactor;
}

double regularization_parameter_IK(const double &alpha, const double &beta, const int &N){
  if(N == 0){
    return 4.*(12.*pow(alpha, 3) + pow(alpha, 2)*(8. - 3.*pow(beta, 2)) - 4.*alpha*pow(beta, 2) + pow(beta, 2)*(pow(beta, 2) - 8.));
  }else if(N == 1){
    return 8.*beta*(9.*pow(alpha, 2) - 2*alpha*(pow(beta, 2) - 4.) + pow(beta, 2));
  }else if(N == 2){
    return -4.*(8.*pow(alpha, 3) - pow(alpha, 2)*(pow(beta, 2) - 8.) - 8.*alpha*pow(beta, 2) + pow(beta, 2)*(3.*pow(beta, 2) - 8.));
  }else if(N == 3){
    return 8.*beta*(pow(alpha, 3) - 7.*pow(alpha, 2) + alpha*(3.*pow(beta, 2) - 8.) + pow(beta, 2));
  }else{
    return -4.*(4.*pow(alpha, 4) - 4.*pow(alpha, 3) + pow(alpha, 2)*(7.*pow(beta, 2) - 8.) + 12.*alpha*pow(beta, 2) - pow(beta, 2)*(pow(beta, 2) - 8.));
  }
}

double regularization_parameter_IE(const double &alpha, const double &beta, const int &N){
  if(N == 0){
    return -16*(8.*pow(alpha, 3) + pow(alpha, 2)*(4. - 7.*pow(beta, 2)) + alpha*pow(beta, 2)*(pow(beta, 2) - 4.) - pow(beta, 2)*(pow(beta, 2) + 4.));
  }else if(N == 1){
    return -4.*beta*(12.*pow(alpha, 3) - pow(alpha, 2)*(pow(beta, 2) - 52.) - alpha*(12.*pow(beta, 2) - 32.) + pow(beta, 2)*(3.*pow(beta, 2) + 4.));
  }else if(N == 2){
    return 8.*(4.*pow(alpha, 4) + pow(alpha, 3)*(pow(beta, 2) + 12.) - 2.*pow(alpha, 2)*(pow(beta, 2) - 4.) + 3.*alpha*pow(beta, 2)*(pow(beta, 2) - 4.) +  2.*pow(beta, 2)*(3.*pow(beta, 2) - 4.));
  }else if(N == 3){
    return -4.*beta*(8.*pow(alpha, 4) - 4.*pow(alpha, 3) + pow(alpha, 2)*(15.*pow(beta, 2) - 44.) + 4.*alpha*(5.*pow(beta, 2) - 8.) + pow(beta, 2)*(3.*pow(beta, 2) + 4.));
  }else{
    return 16.*(4.*pow(alpha, 5) + 4.*pow(alpha, 4) + pow(alpha, 3)*(7.*pow(beta, 2) - 4.) + pow(alpha, 2)*(11.*pow(beta, 2) - 4.) + (2.*alpha + 1.)*pow(beta, 2)*(pow(beta, 2) + 4.));
  }
}

double regularization_parameter_IN_norm(int N, const double &alpha, const double &beta){
  double Q = alpha + 2. - sqrt(alpha*alpha + beta*beta);
  double w = 2.*sqrt(alpha*alpha + beta*beta)/(alpha + 2. + sqrt(alpha*alpha + beta*beta));

  return Q*regularization_parameter_IK(alpha, beta, N)*elliptic_k(sqrt(w)) + regularization_parameter_IE(alpha, beta, N)*elliptic_e(sqrt(w));
}

int phi_index_count(int a, int b, int c, int d){
  int N = 0;
  if(a == 3){  N += 1; }
  if(b == 3){  N += 1; }
  if(c == 3){  N += 1; }
  if(d == 3){  N += 1; }
  return N;
}

double regularization_P_alpha_beta(int alpha, int beta, const double &a, const double &rp, const double &thp, const double &ut, const double &ur, const double &uth, const double &uph){
  return kerr_metric_blc(alpha, beta, a, rp, thp) + lower_U_alpha(alpha, a, rp, thp, ut, ur, uth, uph)*lower_U_alpha(beta, a, rp, thp, ut, ur, uth, uph);
}

double regularization_P_alpha_beta_gamma(int alpha, int beta, int gamma, const double &a, const double &rp, const double &thp, const double &ut, const double &ur, const double &uth, const double &uph){
  double Pabg = 0.;
  for(int i = 0; i < 4; i++){
    Pabg += lower_U_alpha(i, a, rp, thp, ut, ur, uth, uph)*kerr_connection_blc(i, alpha, beta, a, rp, thp);
  }
  Pabg *= lower_U_alpha(gamma, a, rp, thp, ut, ur, uth, uph);

  return 0.5*partial_kerr_metric_blc(alpha, beta, gamma, a, rp, thp) + Pabg;
}

double regularization_C(int i, int j, int k, const double &thp){
  if(i == 2 && j == 3 && k == 3){
    return 0.5*sin(thp)*cos(thp);
  }else if(i == 3 && j == 2 && k == 3){
    return -0.5*cot(thp);
  }else if(i == 3 && j == 3 && k == 2){
    return -0.5*cot(thp);
  }else{
    return 0.;
  }
}

double regularization_P_mu_abcd(int mu, int a, int b, int c, int d, const double &spin, const double &rp, const double &thp, const double &ut, const double &ur, const double &uth, const double &uph){
  double Pmuabcd = 0.;
  for(int e = 2; e < 4; e++){
    Pmuabcd += (3.*regularization_P_alpha_beta(mu, a, spin, rp, thp, ut, ur, uth, uph)*regularization_P_alpha_beta(b, e, spin, rp, thp, ut, ur, uth, uph) - regularization_P_alpha_beta(mu, e, spin, rp, thp, ut, ur, uth, uph)*regularization_P_alpha_beta(a, b, spin, rp, thp, ut, ur, uth, uph))*regularization_C(e, c, d, thp);
  }

  Pmuabcd += 0.5*(3.*regularization_P_alpha_beta(mu, d, spin, rp, thp, ut, ur, uth, uph)*regularization_P_alpha_beta_gamma(a, b, c, spin, rp, thp, ut, ur, uth, uph) - (2.*regularization_P_alpha_beta_gamma(mu, a, b, spin, rp, thp, ut, ur, uth, uph) + regularization_P_alpha_beta_gamma(a, b, mu, spin, rp, thp, ut, ur, uth, uph))*regularization_P_alpha_beta(c, d, spin, rp, thp, ut, ur, uth, uph));

  return Pmuabcd;
}

double regularization_parameter_I_abcd(int a, int b, int c, int d, const double &spin, const double &rp, const double &thp, const double &ut, const double &ur, const double &uth, const double &uph){
  double alpha = pow(sin(thp), 2)*regularization_P_alpha_beta(2, 2, spin, rp, thp, ut, ur, uth, uph)/regularization_P_alpha_beta(3, 3, spin, rp, thp, ut, ur, uth, uph) - 1.;
  double beta = 2.*sin(thp)*regularization_P_alpha_beta(2, 3, spin, rp, thp, ut, ur, uth, uph)/regularization_P_alpha_beta(3, 3, spin, rp, thp, ut, ur, uth, uph);
  int N = phi_index_count(a, b, c, d);
  double dcoeff = 3.*pow(regularization_P_alpha_beta(3, 3, spin, rp, thp, ut, ur, uth, uph)*pow(sin(thp), -2), 2.5)*pow(alpha*alpha + beta*beta, 2)*pow(4.*alpha + 4. - beta*beta, 1.5)*sqrt(0.5*(alpha + 2. - sqrt(alpha*alpha + beta*beta)));
  return pow(sin(thp), -N)*regularization_parameter_IN_norm(N, alpha, beta)/dcoeff;
}

double regularization_parameter_B_alpha(int alpha, const double &spin, const double &rp, const double &thp, const double &ut, const double &ur, const double &uth, const double &uph){
  double Balpha = 0.;
  for(int a = 2; a < 4; a++){
    for(int b = 2; b < 4; b++){
      for(int c = 2; c < 4; c++){
        for(int d = 2; d < 4; d++){
          Balpha += regularization_P_mu_abcd(alpha, a, b, c, d, spin, rp, thp, ut, ur, uth, uph)*regularization_parameter_I_abcd(a, b, c, d, spin, rp, thp, ut, ur, uth, uph);
        }
      }
    }
  }

  return Balpha/(2.*M_PI);
}

Vector regularization_parameter_Bt_from_source(GeodesicSource &geo, const int &sampleSize){
  return regularization_parameter_from_source(geo, sampleSize, *regularization_parameter_Bt);
}
Vector regularization_parameter_Br_from_source(GeodesicSource &geo, const int &sampleSize){
  return regularization_parameter_from_source(geo, sampleSize, *regularization_parameter_Br);
}
Vector regularization_parameter_Btheta_from_source(GeodesicSource &geo, const int &sampleSize){
  return regularization_parameter_from_source(geo, sampleSize, *regularization_parameter_Btheta);
}
Vector regularization_parameter_Bphi_from_source(GeodesicSource &geo, const int &sampleSize){
  return regularization_parameter_from_source(geo, sampleSize, *regularization_parameter_Bphi);
}

double regularization_parameter_Bt(const double &rp, const double &thp, const int &sgnUr, const int &sgnUth, GeodesicSource &geo){
  return regularization_parameter_B_alpha(0, rp, thp, sgnUr, sgnUth, geo);
}
double regularization_parameter_Br(const double &rp, const double &thp, const int &sgnUr, const int &sgnUth, GeodesicSource &geo){
  return regularization_parameter_B_alpha(1, rp, thp, sgnUr, sgnUth, geo);
}
double regularization_parameter_Btheta(const double &rp, const double &thp, const int &sgnUr, const int &sgnUth, GeodesicSource &geo){
  return regularization_parameter_B_alpha(2, rp, thp, sgnUr, sgnUth, geo);
}
double regularization_parameter_Bphi(const double &rp, const double &thp, const int &sgnUr, const int &sgnUth, GeodesicSource &geo){
  return regularization_parameter_B_alpha(3, rp, thp, sgnUr, sgnUth, geo);
}

double regularization_parameter_B_alpha(int alpha, const double &rp, const double &thp, const int &sgnUr, const int &sgnUth, GeodesicSource &geo){
  double a = geo.getBlackHoleSpin();
  double En = geo.getOrbitalEnergy();
  double Lz = geo.getOrbitalAngularMomentum();
  // std::cout << "Lz = " <<Lz << "\n";
  double Qc = geo.getCarterConstant();

  double sig = rp*rp + pow(a*cos(thp), 2);
  double ut = (kerr_geo_VtR(a, En, Lz, Qc, rp) + kerr_geo_VtTheta(a, En, Lz, Qc, thp))/sig;
  double ur = kerr_geo_Vr(a, En, Lz, Qc, rp);
  if(ur/pow(rp, 4) < DBL_EPSILON || ur < 0.){
    ur = 0.;
  }else{
    ur = sgnUr*sqrt(ur)/sig;
  }
  double uth = kerr_geo_Vtheta(a, En, Lz, Qc, thp);
  if(std::abs(uth/Qc) < DBL_EPSILON || uth < 0.){
    uth = 0.;
  }else{
    uth = sgnUth*sqrt(uth)/sig;
  }
  double uph = (kerr_geo_VphiR(a, En, Lz, Qc, rp) + kerr_geo_VphiTheta(a, En, Lz, Qc, thp))/sig;

  return regularization_parameter_B_alpha(alpha, geo.getBlackHoleSpin(), rp, thp, ut, ur, uth, uph);
}
