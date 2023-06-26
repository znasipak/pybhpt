// redshift.cpp

#include "redshift.hpp"

#define HUU_EXPANSION_ORDER 0


RedshiftCoefficients::RedshiftCoefficients(Gauge gauge, GeodesicSource &geo){
  if(gauge == IRG || gauge == SAAB0 || gauge == ASAAB0){
    _coeffs = redshift_coefficients_IRG(geo, geo.getOrbitalSampleNumber());
  }else{
    _coeffs = redshift_coefficients_ORG(geo, geo.getOrbitalSampleNumber());
  }
}

Complex RedshiftCoefficients::getComponent(int Ni, int ai, int bi, int ci, int di, int jr, int jz){
  return redshift_coefficient_components(_coeffs, Ni, ai, bi, ci, di, jr, jz);
}

void redshift_circular_reconstruction_RG(std::string filename, Gauge gauge, int lmax, GeodesicSource &geoCirc){
  double regParam = redshift_regularization_circular(geoCirc);
  ComplexTensor huuCoeffCirc = redshift_coefficients_IRG(geoCirc, pow(2, 2));

  int s = -2;
  if(gauge == ORG){
    s = 2;
    huuCoeffCirc = redshift_coefficients_ORG(geoCirc, pow(2, 2));
  }
  int modeCoupling = 20;
  if(geoCirc.getBlackHoleSpin() == 0.){
    modeCoupling = 5;
  }

  // first we create the matrices for the spherical harmonic basis
  ComplexMatrix huuCircInYlm(lmax + 1, ComplexVector(lmax + 1, 0.));
  ComplexMatrix huuCircUpYlm(lmax + 1, ComplexVector(lmax + 1, 0.));

  // then we create the matrices for the unprojected basis
  ComplexMatrix huuCircInSlm(lmax + 1, ComplexVector(lmax + 1, 0.));
  ComplexMatrix huuCircUpSlm(lmax + 1, ComplexVector(lmax + 1, 0.));
  for(int m = 0; m <= lmax; m++){
    int lmin = abs(m) < 2 ? 2 : abs(m);
    SphericalHarmonicCoupling Cjlm(lmax + 6, m);
    Cjlm.generateCouplings();
    for(int j = lmin; j <= lmax + modeCoupling; j++){
      TeukolskyMode teuk(s, j, m, 0, 0, geoCirc);
      teuk.generateSolutions(geoCirc);
      HertzMode hertz(teuk, gauge);
      hertz.generateSolutions();

      for(int l = abs(m); l <= lmax + modeCoupling; l++){
        Complex tempMode = redshift_mode(huuCoeffCirc, Cjlm, hertz, In, l, 0, 0, 0, 0);
        std::cout << j << " " << l << " " << m <<  "\n";
        std::cout << tempMode << "\n";
        if(abs(tempMode/regParam) > DBL_EPSILON){
          if(l <= lmax){
            huuCircInYlm[l][m] += tempMode;
          }

          if(j <= lmax){
            huuCircInSlm[j][m] += tempMode;
          }
        }
        // std::cout << huuCircInYlm[l][m] << "\n";

        tempMode = redshift_mode(huuCoeffCirc, Cjlm, hertz, Up, l, 0, 0, 0, 0);
        if(abs(tempMode/regParam) > DBL_EPSILON){
          if(l <= lmax){
            huuCircUpYlm[l][m] += tempMode;
          }

          if(j <= lmax){
            huuCircUpSlm[j][m] += tempMode;
          }
        }
      }
    }
  }

  // after calculating all of the individual m-modes, we sum over them, first starting with the Ylm basis
  ComplexVector huuCircInL(lmax + 1);
  ComplexVector huuCircUpL(lmax + 1);
  for(int l = 0; l <= lmax; l++){
    huuCircInL[l] = huuCircInYlm[l][0];
    huuCircUpL[l] = huuCircUpYlm[l][0];
    for(int m = 1; m <= l; m++){
      huuCircInL[l] += std::real((1. + pow(-1., l + m))*huuCircInYlm[l][m]);
      huuCircUpL[l] += std::real((1. + pow(-1., l + m))*huuCircUpYlm[l][m]);
    }
    std::cout << "HuuIn(l = "<<l<<") = " << huuCircInL[l] << ", HuuUp(l = "<<l<<") = " << huuCircUpL[l] << "\n";
  }
  // export results to text files
  std::string file_gauge;
  if(gauge == ORG){
    file_gauge = "ORG";
  }else{
    file_gauge = "IRG";
  }
  std::string file = filename + "_" + file_gauge + "_Yl.txt";
  export_circular_redshift_data(file, huuCircInL, huuCircUpL);
  file = filename + "_" + file_gauge + "_Ylm.txt";
  export_circular_redshift_data_lm(file, huuCircInYlm, huuCircUpYlm);

  // now we repeat the summation and export in the unprojected basis
  for(int l = 0; l <= lmax; l++){
    huuCircInL[l] = huuCircInSlm[l][0];
    huuCircUpL[l] = huuCircUpSlm[l][0];
    for(int m = 1; m <= l; m++){
      huuCircInL[l] += 2.*huuCircInSlm[l][m];
      huuCircUpL[l] += 2.*huuCircUpSlm[l][m];
    }
    std::cout << "HuuIn(j = "<<l<<") = " << huuCircInL[l] << ", HuuUp(l = "<<l<<") = " << huuCircUpL[l] << "\n";
  }
  file = filename + "_" + file_gauge + "_l.txt";
  export_circular_redshift_data(file, huuCircInL, huuCircUpL);
  file = filename + "_" + file_gauge + "_lm.txt";
  export_circular_redshift_data_lm(file, huuCircInSlm, huuCircUpSlm);

  std::cout << redshift_regularization_circular(geoCirc) << "\n";
  std::cout << redshift_completion_circular(geoCirc) << "\n";
  std::cout << geoCirc.getMinoFrequency(0)/pow(geoCirc.getRadialPosition(0), 2) << "\n";
}

void redshift_circular_reconstruction_AAB(std::string filename, Gauge gauge0, Gauge gauge4, int lmax, GeodesicSource &geoCirc){
  double regParam = redshift_regularization_circular(geoCirc);
  ComplexTensor huuCoeffCircORG = redshift_coefficients_ORG(geoCirc, pow(2, 2));
  ComplexTensor huuCoeffCircIRG = redshift_coefficients_IRG(geoCirc, pow(2, 2));
  double prefactor0, prefactor4;
  if(gauge0 == SAAB0){
    prefactor0 = 4.;
    prefactor4 = 4.;
  }else{
    prefactor0 = 1./3.;
    prefactor4 = -1./3.;
  }
  int modeCoupling = 20;
  if(geoCirc.getBlackHoleSpin() == 0.){
    modeCoupling = 5;
  }

  // first we create the matrices for the spherical harmonic basis
  ComplexMatrix huuCircIn0Ylm(lmax + 1, ComplexVector(lmax + 1));
  ComplexMatrix huuCircUp0Ylm(lmax + 1, ComplexVector(lmax + 1));
  ComplexMatrix huuCircUp4Ylm(lmax + 1, ComplexVector(lmax + 1));
  ComplexMatrix huuCircIn4Ylm(lmax + 1, ComplexVector(lmax + 1));
  ComplexMatrix huuCircInYlm(lmax + 1, ComplexVector(lmax + 1));
  ComplexMatrix huuCircUpYlm(lmax + 1, ComplexVector(lmax + 1));

  // then we create the matrices for the unprojected basis
  ComplexMatrix huuCircIn0Slm(lmax + 1, ComplexVector(lmax + 1));
  ComplexMatrix huuCircIn4Slm(lmax + 1, ComplexVector(lmax + 1));
  ComplexMatrix huuCircUp0Slm(lmax + 1, ComplexVector(lmax + 1));
  ComplexMatrix huuCircUp4Slm(lmax + 1, ComplexVector(lmax + 1));
  ComplexMatrix huuCircInSlm(lmax + 1, ComplexVector(lmax + 1));
  ComplexMatrix huuCircUpSlm(lmax + 1, ComplexVector(lmax + 1));
  for(int m = 0; m <= lmax; m++){
    int lmin = abs(m) < 2 ? 2 : abs(m);
    SphericalHarmonicCoupling Cjlm(lmax + 6, m);
    Cjlm.generateCouplings();
    for(int j = lmin; j <= lmax + modeCoupling; j++){
      TeukolskyMode teuk(-2, j, m, 0, 0, geoCirc);
      teuk.generateSolutions(geoCirc);
      HertzMode hertz0(teuk, gauge0);
      hertz0.generateSolutions();
      HertzMode hertz4(teuk, gauge4);
      hertz4.generateSolutions();

      for(int l = abs(m); l <= lmax + modeCoupling; l++){
        Complex tempMode0 = redshift_mode(huuCoeffCircORG, Cjlm, hertz0, In, l, 0, 0, 0, 0);
        Complex tempMode4 = redshift_mode(huuCoeffCircIRG, Cjlm, hertz4, In, l, 0, 0, 0, 0);
        if(abs(tempMode0/regParam) > DBL_EPSILON || abs(tempMode4/regParam) > DBL_EPSILON){
          if(l <= lmax){
            huuCircIn0Ylm[l][m] += tempMode0;
            huuCircIn4Ylm[l][m] += tempMode4;
            huuCircInYlm[l][m] += prefactor0*tempMode0 + prefactor4*tempMode4;
          }

          if(j <= lmax){
            huuCircIn0Slm[j][m] += tempMode0;
            huuCircIn4Slm[j][m] += tempMode4;
            huuCircInSlm[j][m] += prefactor0*tempMode0 + prefactor4*tempMode4;
          }
        }
        tempMode0 = redshift_mode(huuCoeffCircORG, Cjlm, hertz0, Up, l, 0, 0, 0, 0);
        tempMode4 = redshift_mode(huuCoeffCircIRG, Cjlm, hertz4, Up, l, 0, 0, 0, 0);
        if(abs(tempMode0/regParam) > DBL_EPSILON || abs(tempMode4/regParam) > DBL_EPSILON){
          if(l <= lmax){
            huuCircUp0Ylm[l][m] += tempMode0;
            huuCircUp4Ylm[l][m] += tempMode4;
            huuCircUpYlm[l][m] += prefactor0*tempMode0 + prefactor4*tempMode4;
          }

          if(j <= lmax){
            huuCircUp0Slm[j][m] += tempMode0;
            huuCircUp4Slm[j][m] += tempMode4;
            huuCircUpSlm[j][m] += prefactor0*tempMode0 + prefactor4*tempMode4;
          }
        }
      }
    }
  }

  // after calculating all of the individual m-modes, we sum over them, first starting with the Ylm basis
  ComplexVector huuCircIn0L(lmax + 1);
  ComplexVector huuCircIn4L(lmax + 1);
  ComplexVector huuCircUp0L(lmax + 1);
  ComplexVector huuCircUp4L(lmax + 1);
  ComplexVector huuCircInL(lmax + 1);
  ComplexVector huuCircUpL(lmax + 1);
  for(int l = 0; l <= lmax; l++){
    huuCircIn0L[l] = huuCircIn0Ylm[l][0];
    huuCircUp0L[l] = huuCircUp0Ylm[l][0];
    huuCircIn4L[l] = huuCircIn4Ylm[l][0];
    huuCircUp4L[l] = huuCircUp4Ylm[l][0];
    huuCircInL[l] = huuCircInYlm[l][0];
    huuCircUpL[l] = huuCircUpYlm[l][0];
    for(int m = 1; m <= l; m++){
      huuCircIn0L[l] += 2.*huuCircIn0Ylm[l][m];
      huuCircUp0L[l] += 2.*huuCircUp0Ylm[l][m];
      huuCircIn4L[l] += 2.*huuCircIn4Ylm[l][m];
      huuCircUp4L[l] += 2.*huuCircUp4Ylm[l][m];
      huuCircInL[l] += 2.*huuCircInYlm[l][m];
      huuCircUpL[l] += 2.*huuCircUpYlm[l][m];
    }
    std::cout << "HuuIn(l = "<<l<<") = " << huuCircInL[l] << ", HuuUp(l = "<<l<<") = " << huuCircUpL[l] << "\n";
  }
  // export results to text files
  std::string file_gauge;
  if(gauge0 == SAAB0){
    file_gauge = "SAAB";
  }else{
    file_gauge = "ASAAB";
  }
  std::string file = filename + "_" + file_gauge + "0_Yl.txt";
  export_circular_redshift_data(file, huuCircIn0L, huuCircUp0L);
  file = filename + "_" + file_gauge + "4_Yl.txt";
  export_circular_redshift_data(file, huuCircIn4L, huuCircUp4L);
  file = filename + "_" + file_gauge + "0_Ylm.txt";
  export_circular_redshift_data_lm(file, huuCircIn0Ylm, huuCircUp0Ylm);
  file = filename + "_" + file_gauge + "4_Ylm.txt";
  export_circular_redshift_data_lm(file, huuCircIn4Ylm, huuCircUp4Ylm);
  file = filename + "_" + file_gauge + "_Yl.txt";
  export_circular_redshift_data(file, huuCircInL, huuCircUpL);
  file = filename + "_" + file_gauge + "_Ylm.txt";
  export_circular_redshift_data_lm(file, huuCircInYlm, huuCircUpYlm);

  // now we repeat the summation and export in the unprojected basis
  for(int l = 0; l <= lmax; l++){
    huuCircIn0L[l] = huuCircIn0Slm[l][0];
    huuCircUp0L[l] = huuCircUp0Slm[l][0];
    huuCircIn4L[l] = huuCircIn4Slm[l][0];
    huuCircUp4L[l] = huuCircUp4Slm[l][0];
    huuCircInL[l] = huuCircInSlm[l][0];
    huuCircUpL[l] = huuCircUpSlm[l][0];
    for(int m = 1; m <= l; m++){
      huuCircIn0L[l] += 2.*huuCircIn0Slm[l][m];
      huuCircUp0L[l] += 2.*huuCircUp0Slm[l][m];
      huuCircIn4L[l] += 2.*huuCircIn4Slm[l][m];
      huuCircUp4L[l] += 2.*huuCircUp4Slm[l][m];
      huuCircInL[l] += 2.*huuCircInSlm[l][m];
      huuCircUpL[l] += 2.*huuCircUpSlm[l][m];
    }
    std::cout << "HuuIn(j = "<<l<<") = " << huuCircInL[l] << ", HuuUp(l = "<<l<<") = " << huuCircUpL[l] << "\n";
  }
  file = filename + "_" + file_gauge + "0_l.txt";
  export_circular_redshift_data(file, huuCircIn0L, huuCircUp0L);
  file = filename + "_" + file_gauge + "4_l.txt";
  export_circular_redshift_data(file, huuCircIn4L, huuCircUp4L);
  file = filename + "_" + file_gauge + "0_lm.txt";
  export_circular_redshift_data_lm(file, huuCircIn0Slm, huuCircUp0Slm);
  file = filename + "_" + file_gauge + "4_lm.txt";
  export_circular_redshift_data_lm(file, huuCircIn4Slm, huuCircUp4Slm);
  file = filename + "_" + file_gauge + "_l.txt";
  export_circular_redshift_data(file, huuCircInL, huuCircUpL);
  file = filename + "_" + file_gauge + "_lm.txt";
  export_circular_redshift_data_lm(file, huuCircInSlm, huuCircUpSlm);

  std::cout << redshift_regularization_circular(geoCirc) << "\n";
  std::cout << redshift_completion_circular(geoCirc) << "\n";
  std::cout << geoCirc.getMinoFrequency(0)/pow(geoCirc.getRadialPosition(0), 2) << "\n";
}

void redshift_circular(std::string filename, Gauge gauge, int lmax, GeodesicSource &geoCirc){
  if(gauge == ORG || gauge == IRG){
    redshift_circular_reconstruction_RG(filename, gauge, lmax, geoCirc);
  }else{
    if(gauge == SAAB){
      redshift_circular_reconstruction_AAB(filename, SAAB0, SAAB4, lmax, geoCirc);
    }else{
      redshift_circular_reconstruction_AAB(filename, ASAAB0, ASAAB4, lmax, geoCirc);
    }
  }
}

int coefficient_abcd_to_iter(int ai, int bi, int ci, int di){
  return 3*(3*((3*ai) + bi) + ci) + di;
}

void export_circular_redshift_data(std::string filename, ComplexVector HuuIn, ComplexVector HuuUp){
	char buff[500];

	std::cout << "Saving file to " << filename << "\n";
	std::ofstream file;
	file.open(filename);

	file << "l\tRe[HuuIn]\tIm[HuuIn]\tRe[HuuUp]\tIm[HuuUp]\n";
  int lmax = HuuIn.size() - 1;

	for(int l = 0; l <= lmax; l++) {
    sprintf(buff, "%d\t%.15e\t%.15e\t%.15e\t%.15e\n", l, std::real(HuuIn[l]), std::imag(HuuIn[l]), std::real(HuuUp[l]), std::imag(HuuUp[l]));
    file << buff;
	}
	file.close();
}

void export_circular_redshift_data_lm(std::string filename, ComplexMatrix HuuIn, ComplexMatrix HuuUp){
	char buff[500];

	std::cout << "Saving file to " << filename << "\n";
	std::ofstream file;
	file.open(filename);

	file << "l\tm\tRe[HuuIn]\tIm[HuuIn]\tRe[HuuUp]\tIm[HuuUp]\n";
  int lmax = HuuIn.size() - 1;

	for(int l = 0; l <= lmax; l++) {
    for(int m = 0; m <= l; m++){
      sprintf(buff, "%d\t%d\t%.15e\t%.15e\t%.15e\t%.15e\n", l, m, std::real(HuuIn[l][m]), std::imag(HuuIn[l][m]), std::real(HuuUp[l][m]), std::imag(HuuUp[l][m]));
      file << buff;
    }
	}
	file.close();
}

Complex redshift_mode(ComplexTensor &huuCoeff, SphericalHarmonicCoupling &Cjlm, HertzMode &hertz, BoundaryCondition bc, int l, int jr, int jz, int sgnUr, int sgnUz){
  double m = hertz.getAzimuthalModeNumber();
  double omega = hertz.getFrequency();
  Complex huu = 0.;
  Complex huuTerm = 0.;
  int jru = jr;
  int jzu = jz;
  Complex hNabcd = 1.;
  if(sgnUr == -1 && jru > 0){
    jru = huuCoeff[0].size() - jr;
  }
  if(sgnUz == -1 && jzu > 0){
    jzu = huuCoeff[0].size() - jz;
  }
  double ylm = Ylm(l, m, hertz.getPolarPoints(jz));
  ComplexVector hertzR = {hertz.getRadialSolution(bc, 0, jr), hertz.getRadialSolution(bc, 1, jr), hertz.getRadialSolution(bc, 2, jr)};
  // if(isnan(abs(hertzR[0])) || isinf(abs(hertzR[0]))){
  //   std::cout << "("<<l<<", "<<hertz.getAzimuthalModeNumber()<<") \n";
  // }
  Vector freqFactor = {1., omega, omega*omega};
  Vector mFactor = {1., -m, m*m};
  ComplexVector imagFactor = {1., -I, -1., I, 1.};
  for(int ai = 0; ai <= 2; ai++){
    for(int bi = 0; bi <= 2; bi++){
      for(int ci = 0; ci <= 2; ci++){
        for(int di = 0; di <= 2; di++){
          if(ai + bi + ci + di <= 2){
            for(int Ni = 0; Ni <= HUU_EXPANSION_ORDER; Ni++){
              hNabcd = redshift_coefficient_components(huuCoeff, Ni, ai, bi, ci, di, jru, jzu);
              for(int Lc = 0; Lc <= ci; Lc++){
                for(int Np = 0; Np <= Ni; Np++){
                  huuTerm = imagFactor[ai+di]*mFactor[di]*freqFactor[ai]*hNabcd*Cjlm.getDerivativeCouplingCoefficient(ci, Lc, l - 2*Lc + ci - 2*Np + Ni)*Cjlm.getZCouplingCoefficient(Ni, Np, l - 2*Np + Ni)*hertzR[bi]*ylm;
                  huu += hertz.getScalarCouplingCoefficient(l - 2*Lc + ci - 2*Np + Ni)*huuTerm;
                }
              }
            }
          }
        }
      }
    }
  }

  return huu;
}

Complex redshift_mode_no_projection(ComplexTensor &huuCoeff, SphericalHarmonicCoupling &Cjlm, HertzMode &hertz, BoundaryCondition bc, int l, int jr, int jz, int sgnUr, int sgnUz){
  double m = hertz.getAzimuthalModeNumber();
  double omega = hertz.getFrequency();
  Complex huu = 0.;
  Complex huuTerm = 0.;
  int jru = jr;
  int jzu = jz;
  Complex hNabcd = 1.;
  if(sgnUr < 0 && jru > 0){
    jru = huuCoeff[0].size() - jr;
  }
  if(sgnUz > 0 && jzu > 0){
    jzu = huuCoeff[0].size() - jz;
  }
  double ylm = Ylm(l, hertz.getAzimuthalModeNumber(), hertz.getPolarPoints(jz));
  ComplexVector hertzR = {hertz.getRadialSolution(bc, 0, jr), hertz.getRadialSolution(bc, 1, jr), hertz.getRadialSolution(bc, 2, jr)};
  Vector freqFactor = {1., omega, omega*omega};
  Vector mFactor = {1., -m, m*m};
  ComplexVector imagFactor = {1., -I, -1., I, 1.};
  for(int ai = 0; ai <= 2; ai++){
    for(int bi = 0; bi <= 2; bi++){
      for(int ci = 0; ci <= 2; ci++){
        for(int di = 0; di <= 2; di++){
          if(ai + bi + ci + di <= 2){
            hNabcd = redshift_coefficient_components(huuCoeff, 0, ai, bi, ci, di, jru, jzu);
            for(int Lc = 0; Lc <= ci; Lc++){
              huuTerm = imagFactor[ai+di]*mFactor[di]*freqFactor[ai]*hNabcd*Cjlm.getDerivativeCouplingCoefficient(ci, Lc, l - 2*Lc + ci )*hertzR[bi]*ylm;
              huu += hertz.getScalarCouplingCoefficient(l - 2*Lc + ci)*huuTerm;
            }
          }
        }
      }
    }
  }

  return huu;
}

double u1_circ(double a, double En, double Lz, double r){
  return (En*(r*r + a*a) - a*Lz)/(2.*r*r);
}

Complex u3_circ(double a, double En, double Lz, double r){
  return I*(a*En - Lz)/(sqrt(2.)*r);
}

Complex redshift_coefficient_components_circular(int bi, int ci, int mInt, double a, double omega, double lambda, double r, double u1, Complex u3){
  double m = mInt;
  if(bi == 0 && ci == 0){
    return 1.4142135623730951*(m - 1.*a*omega)*u1*u3*(Complex(0.,-1.)*a*m*r - 2.*pow(a,2) + 2.*pow(r,2) + Complex(0.,1.)*omega*r*(pow(a,2) + pow(r,2))) -1.*r*(-2.*a*m*(Complex(0.,1.) + 2.*omega*r) + 2.*omega*(Complex(0.,1.) + omega*r)*pow(a,2) + r*(-4. - 1.*lambda + 2.*pow(m,2)))*pow(u1,2) - 0.5*pow(r,-1)*(2.*a*m*r*(Complex(0.,3.) + 2.*r*(Complex(0.,-2.) + omega*r)) + 2.*m*(Complex(0.,1.) + 2.*omega*r)*pow(a,3) - 2.*omega*(Complex(0.,1.) + omega*r)*pow(a,4) + r*(-8. + 2.*r*(-1.*lambda + Complex(0.,3.)*omega*r) + (4. + lambda - 2.*omega*r*(Complex(0.,-1.) + omega*r))*pow(r,2)) + pow(a,2)*(1.*(8. - Complex(0.,6.)*omega*r) + r*(-4. + lambda - 2.*pow(m,2) - 4.*pow(omega,2)*pow(r,2))))*pow(u3,2);
  }else if(bi == 1 && ci == 0){
    return 1.4142135623730951*(m - 1.*a*omega)*r*u1*u3*((-2. + r)*r + pow(a,2)) + ((-2. + r)*r + pow(a,2))*pow(r,-1)*(-1.*r + Complex(0.,1.)*a*m*r + pow(a,2) - Complex(0.,1.)*omega*r*(pow(a,2) + pow(r,2)))*pow(u3,2);
  }else if(bi == 1 && ci == 1){
    return -1.4142135623730951*r*u1*u3*((-2. + r)*r + pow(a,2));
  }else if(bi == 0 && ci == 1){
    return 1.4142135623730951*u1*u3*(Complex(0.,1.)*a*m*r + 2.*pow(a,2) - 2.*pow(r,2) - Complex(0.,1.)*omega*r*(pow(a,2) + pow(r,2))) + 2.*r*(m*r - 1.*a*(Complex(0.,1.) + omega*r))*pow(u1,2);
  }

  return 0.;
}

int redshift_mode_circular_radial(Complex &radial0, Complex &radial1, HertzMode &hertz, BoundaryCondition bc, GeodesicSource &geo){
  int m = hertz.getAzimuthalModeNumber();
  double a = hertz.getBlackHoleSpin();
  double omega = hertz.getFrequency();
  double lambda = hertz.getEigenvalue();

  double r = geo.getRadialPosition(0);
  double u1 = u1_circ(a, geo.getOrbitalEnergy(), geo.getOrbitalAngularMomentum(), r);
  Complex u3 = u3_circ(a, geo.getOrbitalEnergy(), geo.getOrbitalAngularMomentum(), r);

  Complex hab00 = redshift_coefficient_components_circular(0, 0, m, a, omega, lambda, r, u1, u3);
  Complex hab01 = redshift_coefficient_components_circular(0, 1, m, a, omega, lambda, r, u1, u3);
  Complex hab10 = redshift_coefficient_components_circular(1, 0, m, a, omega, lambda, r, u1, u3);
  Complex hab11 = redshift_coefficient_components_circular(1, 1, m, a, omega, lambda, r, u1, u3);
  Complex R0 = hertz.getRadialSolution(bc, 0, 0);
  Complex R1 = hertz.getRadialSolution(bc, 1, 0);
  radial0 = hab00*R0 + hab10*R1;
  radial1 = hab01*R0 + hab11*R1;

  return 0;
}

Complex redshift_mode_circular_l(HertzMode &hertz, int l, Complex radial0, Complex radial1){
  int m = hertz.getAzimuthalModeNumber();
  if( (l + m) % 2 ){ //if l+m is odd
    return 0.;
  }
  double ylm = Ylm(l, m, 0.5*M_PI);

  double Bljm0 = 0.;
  double Bljm1 = 0.;
  for(int n = - 3; n <= 3; n++){
    Bljm0 += hertz.getCouplingCoefficient(l + n)*Asljm(2, l + n, l, m);
    Bljm1 += hertz.getCouplingCoefficient(l + n)*dAsljm(2, l + n, l, m);
  }
  // std::cout << Bljm1 << "\n";

  Complex huu = (Bljm0*radial0 + Bljm1*radial1)*ylm;

  return huu;
}

Complex redshift_mode_circular(HertzMode &hertz, BoundaryCondition bc, int l,  GeodesicSource &geo){
  Complex radial0, radial1;
  redshift_mode_circular_radial(radial0, radial1, hertz, bc, geo);
  return redshift_mode_circular_l(hertz, l, radial0, radial1);
};

double redshift_regularization_circular(GeodesicSource &geo){
  double r = geo.getRadialPosition(0);
  double a = geo.getBlackHoleSpin();
  double Lz = geo.getOrbitalAngularMomentum();
  double zeta = pow(Lz, 2) + pow(r, 2) + 2.*pow(a, 2)/r + pow(a, 2);
  double k = (zeta - pow(r, 2))/zeta;
  return 2.*elliptic_k(sqrt(k))/M_PI/sqrt(zeta);
}

double redshift_completion_circular(GeodesicSource &geo){
  double rp = geo.getRadialPosition(0);
  double a = geo.getBlackHoleSpin();
  double En = geo.getOrbitalEnergy();
  double Lz = geo.getOrbitalAngularMomentum();
  return ((a*(a*En - 2.*Lz)*pow(-1.*a*En + Lz,2) + (a*En - 1.*Lz)*(4.*a*pow(En,2) - 4.*En*Lz - 1.*a*pow(Lz,2))*rp + 2.*En*(a*En - 1.*Lz)*(a*En + Lz)*pow(rp,2) + pow(En,3)*pow(rp,4)))/(rp*pow(pow(a,2) + (-2. + rp)*rp,2));
}

Complex redshift_coefficient_components(ComplexTensor &huu, int Ni, int ai, int bi, int ci, int di, int jr, int jz){
  int compABCD = coefficient_abcd_to_iter(ai, bi, ci, di);
  return huu[3*compABCD + Ni][jr][jz];
}

void huu_coeffs_subfunction_ORG(Complex &huu0, Complex &huu1, Complex &huu2, double z, Complex u1, Complex u1dz, Complex u1dz2, Complex u3, Complex u3dz, Complex u3dz2, Complex h11, Complex h11dz, Complex h11dz2, Complex h13, Complex h13dz, Complex h13dz2, Complex h33, Complex h33dz, Complex h33dz2){
  Complex h0 = u1*u1*h11 + 2.*u1*u3*h13 + u3*u3*h33;

  if(HUU_EXPANSION_ORDER > 0){
    Complex h1 = 2.*(u1*u1dz*h11 + u1*u3dz*h13 + u3*u1dz*h13 + u3*u3dz*h33) + (u1*u1*h11dz + 2.*u1*u3*h13dz + u3*u3*h33dz);
    Complex h2 = 2.*(u1dz*u1dz*h11 + 2.*u1dz*u3dz*h13 + u3dz*u3dz*h33) + 2.*(u1*u1dz2*h11 + u1*u3dz2*h13 + u3*u1dz2*h13 + u3*u3dz2*h33) + 4.*(u1*u1dz*h11dz + u1*u3dz*h13dz + u3*u1dz*h13dz + u3*u3dz*h33dz) + (u1*u1*h11dz2 + 2.*u1*u3*h13dz2 + u3*u3*h33dz2);
    huu0 = h0 - z*h1 + 0.5*z*z*h2;
    huu1 = h1 - z*h2;
    huu2 = 0.5*h2;

  }else{// zeroth order expansion
    huu0 = h0;
    huu1 = 0.;
    huu2 = 0.;
  }

  // huu0 = h33;
}

void huu_coeffs_subfunction_IRG(Complex &huu0, Complex &huu1, Complex &huu2, double z, Complex u2, Complex u2dz, Complex u2dz2, Complex u4, Complex u4dz, Complex u4dz2, Complex h22, Complex h22dz, Complex h22dz2, Complex h24, Complex h24dz, Complex h24dz2, Complex h44, Complex h44dz, Complex h44dz2){
  Complex h0 = u2*u2*h22 + 2.*u2*u4*h24 + u4*u4*h44;
  // std::cout << "u2 = " << u2 << ", u4 = " << u4 << ", h22 = " << h22 << ", h24 = " << h24 << ", h44 = " << h44 << "\n";

  if(HUU_EXPANSION_ORDER > 0){
    Complex h1 = 2.*(u2*u2dz*h22 + u2*u4dz*h24 + u4*u2dz*h24 + u4*u4dz*h44) + (u2*u2*h22dz + 2.*u2*u4*h24dz + u4*u4*h44dz);
    Complex h2 = 2.*(u2dz*u2dz*h22 + 2.*u2dz*u4dz*h24 + u4dz*u4dz*h44) + 2.*(u2*u2dz2*h22 + u2*u4dz2*h24 + u4*u2dz2*h24 + u4*u4dz2*h44) + 4.*(u2*u2dz*h22dz + u2*u4dz*h24dz + u4*u2dz*h24dz + u4*u4dz*h44dz) + (u2*u2*h22dz2 + 2.*u2*u4*h24dz2 + u4*u4*h44dz2);
    huu0 = h0 - z*h1 + 0.5*z*z*h2;
    huu1 = h1 - z*h2;
    huu2 = 0.5*h2;

  }else{// zeroth order expansion
    huu0 = h0;
    huu1 = 0.;
    huu2 = 0.;
  }

  // huu0 = h44;
}

ComplexTensor redshift_coefficients_ORG(GeodesicSource geo, int sampleNum){
  double a = geo.getBlackHoleSpin();
  double En = geo.getOrbitalEnergy();
  double Lz = geo.getOrbitalAngularMomentum();
  double Qc = geo.getCarterConstant();
  Vector r(sampleNum/2 + 1);
  Vector z(sampleNum/2 + 1);
  int geoSampleNum = geo.getRadialPosition().size() - 1;
  int sampleFactor = 2*geoSampleNum/sampleNum;

  for(size_t ii = 0; ii < r.size(); ii++){
    r[ii] = geo.getRadialPosition(sampleFactor*ii);
    z[ii] = cos(geo.getPolarPosition(sampleFactor*ii));
  }

  return redshift_coefficients_ORG(a, En, Lz, Qc, r, z);
}

ComplexTensor redshift_coefficients_IRG(GeodesicSource geo, int sampleNum){
  double a = geo.getBlackHoleSpin();
  double En = geo.getOrbitalEnergy();
  double Lz = geo.getOrbitalAngularMomentum();
  double Qc = geo.getCarterConstant();
  Vector r(sampleNum/2 + 1);
  Vector z(sampleNum/2 + 1);
  int geoSampleNum = geo.getRadialPosition().size() - 1;
  int sampleFactor = 2*geoSampleNum/sampleNum;

  for(size_t ii = 0; ii < r.size(); ii++){
    r[ii] = geo.getRadialPosition(sampleFactor*ii);
    z[ii] = cos(geo.getPolarPosition(sampleFactor*ii));
  }

  return redshift_coefficients_IRG(a, En, Lz, Qc, r, z);
}

void rescale_redshift_coefficients_ORG(ComplexTensor &huu){
  Complex rescale = 1.;
  for(int ai = 0; ai <= 2; ai++){
    rescale *= pow(-I, ai);
    for(int di = 0; di <= 2; di++){
      rescale *= pow(I, di);
      for(int bi = 0; bi <= 2; bi++){
        for(int ci = 0; ci <= 2; ci++){
          for(int Ni = 0; Ni <= 2; Ni++){
            for(size_t jr = 0; jr < huu[0].size(); jr++){
              for(size_t jz = 0; jz < huu[0].size(); jz++){
                huu[3*coefficient_abcd_to_iter(ai, bi, ci, di) + Ni][jr][jz] *= rescale;
              }
            }
          }
        }
      }
    }
  }
}


ComplexTensor redshift_coefficients_ORG(double a, double En, double Lz, double Qc, Vector r, Vector z){
  ComplexTensor h11 = metric_coefficients_ORG_11(a, r, z);
  ComplexTensor h11dz = metric_coefficients_ORG_11_dz(a, r, z);
  ComplexTensor h11dz2 = metric_coefficients_ORG_11_dz2(a, r, z);
  ComplexTensor h13 = metric_coefficients_ORG_13(a, r, z);
  ComplexTensor h13dz = metric_coefficients_ORG_13_dz(a, r, z);
  ComplexTensor h13dz2 = metric_coefficients_ORG_13_dz2(a, r, z);
  ComplexTensor h33 = metric_coefficients_ORG_33(a, r, z);
  ComplexTensor h33dz = metric_coefficients_ORG_33_dz(a, r, z);
  ComplexTensor h33dz2 = metric_coefficients_ORG_33_dz2(a, r, z);

  ComplexMatrix u1 = tetrad_velocity_1(a, En, Lz, Qc, r, z);
  ComplexMatrix u1dz = tetrad_velocity_1_dz(a, En, Lz, Qc, r, z);
  ComplexMatrix u1dz2 = tetrad_velocity_1_dz2(a, En, Lz, Qc, r, z);
  ComplexMatrix u3 = tetrad_velocity_3(a, En, Lz, Qc, r, z);
  ComplexMatrix u3dz = tetrad_velocity_3_dz(a, En, Lz, Qc, r, z);
  ComplexMatrix u3dz2 = tetrad_velocity_3_dz2(a, En, Lz, Qc, r, z);

  int compNum = h11.size();
  int rCompNum = u1.size();
  int zCompNum = u1[0].size();
  ComplexTensor huu(3*compNum, ComplexMatrix(rCompNum, ComplexVector(zCompNum)));

  Complex huu0, huu1, huu2;
  size_t jr, jz;
  for(size_t i = 0; i < h11.size(); i++){
    // first consider ur = 0 and uz = 0
    jr = 0;
    jz = 0;
    huu_coeffs_subfunction_ORG(huu0, huu1, huu2, z[jz], u1[jr][jz], u1dz[jr][jz], u1dz2[jr][jz], u3[jr][jz], u3dz[jr][jz], u3dz2[jr][jz], h11[i][jr][jz], h11dz[i][jr][jz], h11dz2[i][jr][jz], h13[i][jr][jz], h13dz[i][jr][jz], h13dz2[i][jr][jz], h33[i][jr][jz], h33dz[i][jr][jz], h33dz2[i][jr][jz]);

    huu[3*i][jr][jz] = huu0; //n = 0;
    huu[3*i + 1][jr][jz] = huu1; //n = 1;
    huu[3*i + 2][jr][jz] = huu2; //n = 2;

    jr = 0;
    jz = z.size() - 1;
    huu_coeffs_subfunction_ORG(huu0, huu1, huu2, z[jz], u1[jr][jz], u1dz[jr][jz], u1dz2[jr][jz], u3[jr][jz], u3dz[jr][jz], u3dz2[jr][jz], h11[i][jr][jz], h11dz[i][jr][jz], h11dz2[i][jr][jz], h13[i][jr][jz], h13dz[i][jr][jz], h13dz2[i][jr][jz], h33[i][jr][jz], h33dz[i][jr][jz], h33dz2[i][jr][jz]);

    huu[3*i][jr][jz] = huu0; //n = 0;
    huu[3*i + 1][jr][jz] = huu1; //n = 1;
    huu[3*i + 2][jr][jz] = huu2; //n = 2;

    jr = r.size() - 1;
    jz = 0;
    huu_coeffs_subfunction_ORG(huu0, huu1, huu2, z[jz], u1[jr][jz], u1dz[jr][jz], u1dz2[jr][jz], u3[jr][jz], u3dz[jr][jz], u3dz2[jr][jz], h11[i][jr][jz], h11dz[i][jr][jz], h11dz2[i][jr][jz], h13[i][jr][jz], h13dz[i][jr][jz], h13dz2[i][jr][jz], h33[i][jr][jz], h33dz[i][jr][jz], h33dz2[i][jr][jz]);

    huu[3*i][jr][jz] = huu0; //n = 0;
    huu[3*i + 1][jr][jz] = huu1; //n = 1;
    huu[3*i + 2][jr][jz] = huu2; //n = 2;

    jr = r.size() - 1;
    jz = z.size() - 1;
    huu_coeffs_subfunction_ORG(huu0, huu1, huu2, z[jz], u1[jr][jz], u1dz[jr][jz], u1dz2[jr][jz], u3[jr][jz], u3dz[jr][jz], u3dz2[jr][jz], h11[i][jr][jz], h11dz[i][jr][jz], h11dz2[i][jr][jz], h13[i][jr][jz], h13dz[i][jr][jz], h13dz2[i][jr][jz], h33[i][jr][jz], h33dz[i][jr][jz], h33dz2[i][jr][jz]);

    huu[3*i][jr][jz] = huu0; //n = 0;
    huu[3*i + 1][jr][jz] = huu1; //n = 1;
    huu[3*i + 2][jr][jz] = huu2; //n = 2;

    for(jz = 1; jz < z.size() - 1; jz++){
      // next consider ur = 0 and uz < 0
      jr = 0;
      huu_coeffs_subfunction_ORG(huu0, huu1, huu2, z[jz], u1[jr][jz], u1dz[jr][jz], u1dz2[jr][jz], u3[jr][jz], u3dz[jr][jz], u3dz2[jr][jz], h11[i][jr][jz], h11dz[i][jr][jz], h11dz2[i][jr][jz], h13[i][jr][jz], h13dz[i][jr][jz], h13dz2[i][jr][jz], h33[i][jr][jz], h33dz[i][jr][jz], h33dz2[i][jr][jz]);

      huu[3*i][jr][jz] = huu0; //n = 0;
      huu[3*i + 1][jr][jz] = huu1; //n = 1;
      huu[3*i + 2][jr][jz] = huu2; //n = 2;

      // next consider ur = 0 and uz < 0
      jr = r.size() - 1;
      huu_coeffs_subfunction_ORG(huu0, huu1, huu2, z[jz], u1[jr][jz], u1dz[jr][jz], u1dz2[jr][jz], u3[jr][jz], u3dz[jr][jz], u3dz2[jr][jz], h11[i][jr][jz], h11dz[i][jr][jz], h11dz2[i][jr][jz], h13[i][jr][jz], h13dz[i][jr][jz], h13dz2[i][jr][jz], h33[i][jr][jz], h33dz[i][jr][jz], h33dz2[i][jr][jz]);

      huu[3*i][jr][jz] = huu0; //n = 0;
      huu[3*i + 1][jr][jz] = huu1; //n = 1;
      huu[3*i + 2][jr][jz] = huu2; //n = 2;

      // next consider ur = 0 and uz > 0
      jr = 0;
      huu_coeffs_subfunction_ORG(huu0, huu1, huu2, z[jz], u1[jr][zCompNum - jz], u1dz[jr][zCompNum - jz], u1dz2[jr][zCompNum - jz], u3[jr][zCompNum - jz], u3dz[jr][zCompNum - jz], u3dz2[jr][zCompNum - jz], h11[i][jr][jz], h11dz[i][jr][jz], h11dz2[i][jr][jz], h13[i][jr][jz], h13dz[i][jr][jz], h13dz2[i][jr][jz], h33[i][jr][jz], h33dz[i][jr][jz], h33dz2[i][jr][jz]);

      huu[3*i][jr][zCompNum - jz] = huu0; //n = 0;
      huu[3*i + 1][jr][zCompNum - jz] = huu1; //n = 1;
      huu[3*i + 2][jr][zCompNum - jz] = huu2; //n = 2;

      // next consider ur = 0 and uz > 0
      jr = r.size() - 1;
      huu_coeffs_subfunction_ORG(huu0, huu1, huu2, z[jz], u1[jr][zCompNum - jz], u1dz[jr][zCompNum - jz], u1dz2[jr][zCompNum - jz], u3[jr][zCompNum - jz], u3dz[jr][zCompNum - jz], u3dz2[jr][zCompNum - jz], h11[i][jr][jz], h11dz[i][jr][jz], h11dz2[i][jr][jz], h13[i][jr][jz], h13dz[i][jr][jz], h13dz2[i][jr][jz], h33[i][jr][jz], h33dz[i][jr][jz], h33dz2[i][jr][jz]);

      huu[3*i][jr][zCompNum - jz] = huu0; //n = 0;
      huu[3*i + 1][jr][zCompNum - jz] = huu1; //n = 1;
      huu[3*i + 2][jr][zCompNum - jz] = huu2; //n = 2;
    }

    for(jr = 1; jr < r.size() - 1; jr++){
      // next consider ur > 0 and uz = 0
      jz = 0;
      huu_coeffs_subfunction_ORG(huu0, huu1, huu2, z[jz], u1[jr][jz], u1dz[jr][jz], u1dz2[jr][jz], u3[jr][jz], u3dz[jr][jz], u3dz2[jr][jz], h11[i][jr][jz], h11dz[i][jr][jz], h11dz2[i][jr][jz], h13[i][jr][jz], h13dz[i][jr][jz], h13dz2[i][jr][jz], h33[i][jr][jz], h33dz[i][jr][jz], h33dz2[i][jr][jz]);

      huu[3*i][jr][jz] = huu0; //n = 0;
      huu[3*i + 1][jr][jz] = huu1; //n = 1;
      huu[3*i + 2][jr][jz] = huu2; //n = 2;

      // next consider ur < 0 and uz = 0
      jz = 0;
      huu_coeffs_subfunction_ORG(huu0, huu1, huu2, z[jz], u1[rCompNum - jr][jz], u1dz[rCompNum - jr][jz], u1dz2[rCompNum - jr][jz], u3[rCompNum - jr][jz], u3dz[rCompNum - jr][jz], u3dz2[rCompNum - jr][jz], h11[i][jr][jz], h11dz[i][jr][jz], h11dz2[i][jr][jz], h13[i][jr][jz], h13dz[i][jr][jz], h13dz2[i][jr][jz], h33[i][jr][jz], h33dz[i][jr][jz], h33dz2[i][jr][jz]);

      huu[3*i][rCompNum - jr][jz] = huu0; //n = 0;
      huu[3*i + 1][rCompNum - jr][jz] = huu1; //n = 1;
      huu[3*i + 2][rCompNum - jr][jz] = huu2; //n = 2;

      // next consider ur > 0 and uz = 0
      jz = z.size() - 1;
      huu_coeffs_subfunction_ORG(huu0, huu1, huu2, z[jz], u1[jr][jz], u1dz[jr][jz], u1dz2[jr][jz], u3[jr][jz], u3dz[jr][jz], u3dz2[jr][jz], h11[i][jr][jz], h11dz[i][jr][jz], h11dz2[i][jr][jz], h13[i][jr][jz], h13dz[i][jr][jz], h13dz2[i][jr][jz], h33[i][jr][jz], h33dz[i][jr][jz], h33dz2[i][jr][jz]);

      huu[3*i][jr][jz] = huu0; //n = 0;
      huu[3*i + 1][jr][jz] = huu1; //n = 1;
      huu[3*i + 2][jr][jz] = huu2; //n = 2;

      // next consider ur < 0 and uz = 0
      jz = z.size() - 1;
      huu_coeffs_subfunction_ORG(huu0, huu1, huu2, z[jz], u1[rCompNum - jr][jz], u1dz[rCompNum - jr][jz], u1dz2[rCompNum - jr][jz], u3[rCompNum - jr][jz], u3dz[rCompNum - jr][jz], u3dz2[rCompNum - jr][jz], h11[i][jr][jz], h11dz[i][jr][jz], h11dz2[i][jr][jz], h13[i][jr][jz], h13dz[i][jr][jz], h13dz2[i][jr][jz], h33[i][jr][jz], h33dz[i][jr][jz], h33dz2[i][jr][jz]);

      huu[3*i][rCompNum - jr][jz] = huu0; //n = 0;
      huu[3*i + 1][rCompNum - jr][jz] = huu1; //n = 1;
      huu[3*i + 2][rCompNum - jr][jz] = huu2; //n = 2;
    }

    for(jr = 1; jr < r.size() - 1; jr++){
      for(jz = 1; jz < z.size() - 1; jz++){
        // next consider ur > 0 and uz < 0
        huu_coeffs_subfunction_ORG(huu0, huu1, huu2, z[jz], u1[jr][jz], u1dz[jr][jz], u1dz2[jr][jz], u3[jr][jz], u3dz[jr][jz], u3dz2[jr][jz], h11[i][jr][jz], h11dz[i][jr][jz], h11dz2[i][jr][jz], h13[i][jr][jz], h13dz[i][jr][jz], h13dz2[i][jr][jz], h33[i][jr][jz], h33dz[i][jr][jz], h33dz2[i][jr][jz]);

        huu[3*i][jr][jz] = huu0; //n = 0;
        huu[3*i + 1][jr][jz] = huu1; //n = 1;
        huu[3*i + 2][jr][jz] = huu2; //n = 2;

        // next consider ur < 0 and uz < 0
        huu_coeffs_subfunction_ORG(huu0, huu1, huu2, z[jz], u1[rCompNum - jr][jz], u1dz[rCompNum - jr][jz], u1dz2[rCompNum - jr][jz], u3[rCompNum - jr][jz], u3dz[rCompNum - jr][jz], u3dz2[rCompNum - jr][jz], h11[i][jr][jz], h11dz[i][jr][jz], h11dz2[i][jr][jz], h13[i][jr][jz], h13dz[i][jr][jz], h13dz2[i][jr][jz], h33[i][jr][jz], h33dz[i][jr][jz], h33dz2[i][jr][jz]);

        huu[3*i][rCompNum - jr][jz] = huu0; //n = 0;
        huu[3*i + 1][rCompNum - jr][jz] = huu1; //n = 1;
        huu[3*i + 2][rCompNum - jr][jz] = huu2; //n = 2;

        // next consider ur > 0 and uz > 0
        huu_coeffs_subfunction_ORG(huu0, huu1, huu2, z[jz], u1[jr][zCompNum - jz], u1dz[jr][zCompNum - jz], u1dz2[jr][zCompNum - jz], u3[jr][zCompNum - jz], u3dz[jr][zCompNum - jz], u3dz2[jr][zCompNum - jz], h11[i][jr][jz], h11dz[i][jr][jz], h11dz2[i][jr][jz], h13[i][jr][jz], h13dz[i][jr][jz], h13dz2[i][jr][jz], h33[i][jr][jz], h33dz[i][jr][jz], h33dz2[i][jr][jz]);

        huu[3*i][jr][zCompNum - jz] = huu0; //n = 0;
        huu[3*i + 1][jr][zCompNum - jz] = huu1; //n = 1;
        huu[3*i + 2][jr][zCompNum - jz] = huu2; //n = 2;

        // next consider ur < 0 and uz > 0
        huu_coeffs_subfunction_ORG(huu0, huu1, huu2, z[jz], u1[rCompNum - jr][zCompNum - jz], u1dz[rCompNum - jr][zCompNum - jz], u1dz2[rCompNum - jr][zCompNum - jz], u3[rCompNum - jr][zCompNum - jz], u3dz[rCompNum - jr][zCompNum - jz], u3dz2[rCompNum - jr][zCompNum - jz], h11[i][jr][jz], h11dz[i][jr][jz], h11dz2[i][jr][jz], h13[i][jr][jz], h13dz[i][jr][jz], h13dz2[i][jr][jz], h33[i][jr][jz], h33dz[i][jr][jz], h33dz2[i][jr][jz]);

        huu[3*i][rCompNum - jr][zCompNum - jz] = huu0; //n = 0;
        huu[3*i + 1][rCompNum - jr][zCompNum - jz] = huu1; //n = 1;
        huu[3*i + 2][rCompNum - jr][zCompNum - jz] = huu2; //n = 2;
      }
    }
  }
  // rescale_redshift_coefficients_ORG(huu);

  return huu;
}

ComplexTensor redshift_coefficients_IRG(double a, double En, double Lz, double Qc, Vector r, Vector z){
  ComplexTensor h22 = metric_coefficients_IRG_22(a, r, z);
  ComplexTensor h22dz = metric_coefficients_IRG_22_dz(a, r, z);
  ComplexTensor h22dz2 = metric_coefficients_IRG_22_dz2(a, r, z);
  ComplexTensor h24 = metric_coefficients_IRG_24(a, r, z);
  ComplexTensor h24dz = metric_coefficients_IRG_24_dz(a, r, z);
  ComplexTensor h24dz2 = metric_coefficients_IRG_24_dz2(a, r, z);
  ComplexTensor h44 = metric_coefficients_IRG_44(a, r, z);
  ComplexTensor h44dz = metric_coefficients_IRG_44_dz(a, r, z);
  ComplexTensor h44dz2 = metric_coefficients_IRG_44_dz2(a, r, z);

  ComplexMatrix u2 = tetrad_velocity_2(a, En, Lz, Qc, r, z);
  ComplexMatrix u2dz = tetrad_velocity_2_dz(a, En, Lz, Qc, r, z);
  ComplexMatrix u2dz2 = tetrad_velocity_2_dz2(a, En, Lz, Qc, r, z);
  ComplexMatrix u4 = tetrad_velocity_4(a, En, Lz, Qc, r, z);
  ComplexMatrix u4dz = tetrad_velocity_4_dz(a, En, Lz, Qc, r, z);
  ComplexMatrix u4dz2 = tetrad_velocity_4_dz2(a, En, Lz, Qc, r, z);

  int compNum = h22.size();
  int rCompNum = u2.size();
  int zCompNum = u2[0].size();
  ComplexTensor huu(3*compNum, ComplexMatrix(rCompNum, ComplexVector(zCompNum)));

  Complex huu0, huu1, huu2;
  size_t jr, jz;
  for(size_t i = 0; i < h22.size(); i++){
    // first consider ur = 0 and uz = 0
    jr = 0;
    jz = 0;
    huu_coeffs_subfunction_IRG(huu0, huu1, huu2, z[jz], u2[jr][jz], u2dz[jr][jz], u2dz2[jr][jz], u4[jr][jz], u4dz[jr][jz], u4dz2[jr][jz], h22[i][jr][jz], h22dz[i][jr][jz], h22dz2[i][jr][jz], h24[i][jr][jz], h24dz[i][jr][jz], h24dz2[i][jr][jz], h44[i][jr][jz], h44dz[i][jr][jz], h44dz2[i][jr][jz]);

    huu[3*i][jr][jz] = huu0; //n = 0;
    huu[3*i + 1][jr][jz] = huu1; //n = 1;
    huu[3*i + 2][jr][jz] = huu2; //n = 2;

    jr = 0;
    jz = z.size() - 1;
    huu_coeffs_subfunction_IRG(huu0, huu1, huu2, z[jz], u2[jr][jz], u2dz[jr][jz], u2dz2[jr][jz], u4[jr][jz], u4dz[jr][jz], u4dz2[jr][jz], h22[i][jr][jz], h22dz[i][jr][jz], h22dz2[i][jr][jz], h24[i][jr][jz], h24dz[i][jr][jz], h24dz2[i][jr][jz], h44[i][jr][jz], h44dz[i][jr][jz], h44dz2[i][jr][jz]);

    huu[3*i][jr][jz] = huu0; //n = 0;
    huu[3*i + 1][jr][jz] = huu1; //n = 1;
    huu[3*i + 2][jr][jz] = huu2; //n = 2;

    jr = r.size() - 1;
    jz = 0;
    huu_coeffs_subfunction_IRG(huu0, huu1, huu2, z[jz], u2[jr][jz], u2dz[jr][jz], u2dz2[jr][jz], u4[jr][jz], u4dz[jr][jz], u4dz2[jr][jz], h22[i][jr][jz], h22dz[i][jr][jz], h22dz2[i][jr][jz], h24[i][jr][jz], h24dz[i][jr][jz], h24dz2[i][jr][jz], h44[i][jr][jz], h44dz[i][jr][jz], h44dz2[i][jr][jz]);

    huu[3*i][jr][jz] = huu0; //n = 0;
    huu[3*i + 1][jr][jz] = huu1; //n = 1;
    huu[3*i + 2][jr][jz] = huu2; //n = 2;

    jr = r.size() - 1;
    jz = z.size() - 1;
    huu_coeffs_subfunction_IRG(huu0, huu1, huu2, z[jz], u2[jr][jz], u2dz[jr][jz], u2dz2[jr][jz], u4[jr][jz], u4dz[jr][jz], u4dz2[jr][jz], h22[i][jr][jz], h22dz[i][jr][jz], h22dz2[i][jr][jz], h24[i][jr][jz], h24dz[i][jr][jz], h24dz2[i][jr][jz], h44[i][jr][jz], h44dz[i][jr][jz], h44dz2[i][jr][jz]);

    huu[3*i][jr][jz] = huu0; //n = 0;
    huu[3*i + 1][jr][jz] = huu1; //n = 1;
    huu[3*i + 2][jr][jz] = huu2; //n = 2;

    for(jz = 1; jz < z.size() - 1; jz++){
      // next consider ur = 0 and uz < 0
      jr = 0;
      huu_coeffs_subfunction_IRG(huu0, huu1, huu2, z[jz], u2[jr][jz], u2dz[jr][jz], u2dz2[jr][jz], u4[jr][jz], u4dz[jr][jz], u4dz2[jr][jz], h22[i][jr][jz], h22dz[i][jr][jz], h22dz2[i][jr][jz], h24[i][jr][jz], h24dz[i][jr][jz], h24dz2[i][jr][jz], h44[i][jr][jz], h44dz[i][jr][jz], h44dz2[i][jr][jz]);

      huu[3*i][jr][jz] = huu0; //n = 0;
      huu[3*i + 1][jr][jz] = huu1; //n = 1;
      huu[3*i + 2][jr][jz] = huu2; //n = 2;

      // next consider ur = 0 and uz < 0
      jr = r.size() - 1;
      huu_coeffs_subfunction_IRG(huu0, huu1, huu2, z[jz], u2[jr][jz], u2dz[jr][jz], u2dz2[jr][jz], u4[jr][jz], u4dz[jr][jz], u4dz2[jr][jz], h22[i][jr][jz], h22dz[i][jr][jz], h22dz2[i][jr][jz], h24[i][jr][jz], h24dz[i][jr][jz], h24dz2[i][jr][jz], h44[i][jr][jz], h44dz[i][jr][jz], h44dz2[i][jr][jz]);

      huu[3*i][jr][jz] = huu0; //n = 0;
      huu[3*i + 1][jr][jz] = huu1; //n = 1;
      huu[3*i + 2][jr][jz] = huu2; //n = 2;

      // next consider ur = 0 and uz > 0
      jr = 0;
      huu_coeffs_subfunction_IRG(huu0, huu1, huu2, z[jz], u2[jr][zCompNum - jz], u2dz[jr][zCompNum - jz], u2dz2[jr][zCompNum - jz], u4[jr][zCompNum - jz], u4dz[jr][zCompNum - jz], u4dz2[jr][zCompNum - jz], h22[i][jr][jz], h22dz[i][jr][jz], h22dz2[i][jr][jz], h24[i][jr][jz], h24dz[i][jr][jz], h24dz2[i][jr][jz], h44[i][jr][jz], h44dz[i][jr][jz], h44dz2[i][jr][jz]);

      huu[3*i][jr][zCompNum - jz] = huu0; //n = 0;
      huu[3*i + 1][jr][zCompNum - jz] = huu1; //n = 1;
      huu[3*i + 2][jr][zCompNum - jz] = huu2; //n = 2;

      // next consider ur = 0 and uz > 0
      jr = r.size() - 1;
      huu_coeffs_subfunction_IRG(huu0, huu1, huu2, z[jz], u2[jr][zCompNum - jz], u2dz[jr][zCompNum - jz], u2dz2[jr][zCompNum - jz], u4[jr][zCompNum - jz], u4dz[jr][zCompNum - jz], u4dz2[jr][zCompNum - jz], h22[i][jr][jz], h22dz[i][jr][jz], h22dz2[i][jr][jz], h24[i][jr][jz], h24dz[i][jr][jz], h24dz2[i][jr][jz], h44[i][jr][jz], h44dz[i][jr][jz], h44dz2[i][jr][jz]);

      huu[3*i][jr][zCompNum - jz] = huu0; //n = 0;
      huu[3*i + 1][jr][zCompNum - jz] = huu1; //n = 1;
      huu[3*i + 2][jr][zCompNum - jz] = huu2; //n = 2;
    }

    for(jr = 1; jr < r.size() - 1; jr++){
      // next consider ur > 0 and uz = 0
      jz = 0;
      huu_coeffs_subfunction_IRG(huu0, huu1, huu2, z[jz], u2[jr][jz], u2dz[jr][jz], u2dz2[jr][jz], u4[jr][jz], u4dz[jr][jz], u4dz2[jr][jz], h22[i][jr][jz], h22dz[i][jr][jz], h22dz2[i][jr][jz], h24[i][jr][jz], h24dz[i][jr][jz], h24dz2[i][jr][jz], h44[i][jr][jz], h44dz[i][jr][jz], h44dz2[i][jr][jz]);

      huu[3*i][jr][jz] = huu0; //n = 0;
      huu[3*i + 1][jr][jz] = huu1; //n = 1;
      huu[3*i + 2][jr][jz] = huu2; //n = 2;

      // next consider ur < 0 and uz = 0
      jz = 0;
      huu_coeffs_subfunction_IRG(huu0, huu1, huu2, z[jz], u2[rCompNum - jr][jz], u2dz[rCompNum - jr][jz], u2dz2[rCompNum - jr][jz], u4[rCompNum - jr][jz], u4dz[rCompNum - jr][jz], u4dz2[rCompNum - jr][jz], h22[i][jr][jz], h22dz[i][jr][jz], h22dz2[i][jr][jz], h24[i][jr][jz], h24dz[i][jr][jz], h24dz2[i][jr][jz], h44[i][jr][jz], h44dz[i][jr][jz], h44dz2[i][jr][jz]);

      huu[3*i][rCompNum - jr][jz] = huu0; //n = 0;
      huu[3*i + 1][rCompNum - jr][jz] = huu1; //n = 1;
      huu[3*i + 2][rCompNum - jr][jz] = huu2; //n = 2;

      // next consider ur > 0 and uz = 0
      jz = z.size() - 1;
      huu_coeffs_subfunction_IRG(huu0, huu1, huu2, z[jz], u2[jr][jz], u2dz[jr][jz], u2dz2[jr][jz], u4[jr][jz], u4dz[jr][jz], u4dz2[jr][jz], h22[i][jr][jz], h22dz[i][jr][jz], h22dz2[i][jr][jz], h24[i][jr][jz], h24dz[i][jr][jz], h24dz2[i][jr][jz], h44[i][jr][jz], h44dz[i][jr][jz], h44dz2[i][jr][jz]);

      huu[3*i][jr][jz] = huu0; //n = 0;
      huu[3*i + 1][jr][jz] = huu1; //n = 1;
      huu[3*i + 2][jr][jz] = huu2; //n = 2;

      // next consider ur < 0 and uz = 0
      jz = z.size() - 1;
      huu_coeffs_subfunction_IRG(huu0, huu1, huu2, z[jz], u2[rCompNum - jr][jz], u2dz[rCompNum - jr][jz], u2dz2[rCompNum - jr][jz], u4[rCompNum - jr][jz], u4dz[rCompNum - jr][jz], u4dz2[rCompNum - jr][jz], h22[i][jr][jz], h22dz[i][jr][jz], h22dz2[i][jr][jz], h24[i][jr][jz], h24dz[i][jr][jz], h24dz2[i][jr][jz], h44[i][jr][jz], h44dz[i][jr][jz], h44dz2[i][jr][jz]);

      huu[3*i][rCompNum - jr][jz] = huu0; //n = 0;
      huu[3*i + 1][rCompNum - jr][jz] = huu1; //n = 1;
      huu[3*i + 2][rCompNum - jr][jz] = huu2; //n = 2;
    }

    for(jr = 1; jr < r.size() - 1; jr++){
      for(jz = 1; jz < z.size() - 1; jz++){
        // next consider ur > 0 and uz < 0
        huu_coeffs_subfunction_IRG(huu0, huu1, huu2, z[jz], u2[jr][jz], u2dz[jr][jz], u2dz2[jr][jz], u4[jr][jz], u4dz[jr][jz], u4dz2[jr][jz], h22[i][jr][jz], h22dz[i][jr][jz], h22dz2[i][jr][jz], h24[i][jr][jz], h24dz[i][jr][jz], h24dz2[i][jr][jz], h44[i][jr][jz], h44dz[i][jr][jz], h44dz2[i][jr][jz]);

        huu[3*i][jr][jz] = huu0; //n = 0;
        huu[3*i + 1][jr][jz] = huu1; //n = 1;
        huu[3*i + 2][jr][jz] = huu2; //n = 2;

        // next consider ur < 0 and uz < 0
        huu_coeffs_subfunction_IRG(huu0, huu1, huu2, z[jz], u2[rCompNum - jr][jz], u2dz[rCompNum - jr][jz], u2dz2[rCompNum - jr][jz], u4[rCompNum - jr][jz], u4dz[rCompNum - jr][jz], u4dz2[rCompNum - jr][jz], h22[i][jr][jz], h22dz[i][jr][jz], h22dz2[i][jr][jz], h24[i][jr][jz], h24dz[i][jr][jz], h24dz2[i][jr][jz], h44[i][jr][jz], h44dz[i][jr][jz], h44dz2[i][jr][jz]);

        huu[3*i][rCompNum - jr][jz] = huu0; //n = 0;
        huu[3*i + 1][rCompNum - jr][jz] = huu1; //n = 1;
        huu[3*i + 2][rCompNum - jr][jz] = huu2; //n = 2;

        // next consider ur > 0 and uz > 0
        huu_coeffs_subfunction_IRG(huu0, huu1, huu2, z[jz], u2[jr][zCompNum - jz], u2dz[jr][zCompNum - jz], u2dz2[jr][zCompNum - jz], u4[jr][zCompNum - jz], u4dz[jr][zCompNum - jz], u4dz2[jr][zCompNum - jz], h22[i][jr][jz], h22dz[i][jr][jz], h22dz2[i][jr][jz], h24[i][jr][jz], h24dz[i][jr][jz], h24dz2[i][jr][jz], h44[i][jr][jz], h44dz[i][jr][jz], h44dz2[i][jr][jz]);

        huu[3*i][jr][zCompNum - jz] = huu0; //n = 0;
        huu[3*i + 1][jr][zCompNum - jz] = huu1; //n = 1;
        huu[3*i + 2][jr][zCompNum - jz] = huu2; //n = 2;

        // next consider ur < 0 and uz > 0
        huu_coeffs_subfunction_IRG(huu0, huu1, huu2, z[jz], u2[rCompNum - jr][zCompNum - jz], u2dz[rCompNum - jr][zCompNum - jz], u2dz2[rCompNum - jr][zCompNum - jz], u4[rCompNum - jr][zCompNum - jz], u4dz[rCompNum - jr][zCompNum - jz], u4dz2[rCompNum - jr][zCompNum - jz], h22[i][jr][jz], h22dz[i][jr][jz], h22dz2[i][jr][jz], h24[i][jr][jz], h24dz[i][jr][jz], h24dz2[i][jr][jz], h44[i][jr][jz], h44dz[i][jr][jz], h44dz2[i][jr][jz]);

        huu[3*i][rCompNum - jr][zCompNum - jz] = huu0; //n = 0;
        huu[3*i + 1][rCompNum - jr][zCompNum - jz] = huu1; //n = 1;
        huu[3*i + 2][rCompNum - jr][zCompNum - jz] = huu2; //n = 2;
      }
    }
  }
  // rescale_redshift_coefficients_ORG(huu);

  return huu;
}
