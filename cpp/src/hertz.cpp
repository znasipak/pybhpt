// hertz.cpp

#include "hertz.hpp"

HertzMode::HertzMode(TeukolskyMode& teuk, Gauge gauge): _gauge(gauge), _s(teuk.getSpinWeight()), _L(teuk.getSpheroidalModeNumber()), _m(teuk.getAzimuthalModeNumber()), _k(teuk.getPolarModeNumber()), _n(teuk.getRadialModeNumber()), _sampleSize(teuk.getSampleSize()), _a(teuk.getBlackHoleSpin()), _omega(teuk.getFrequency()), _lambda(teuk.getEigenvalue()), _lmin(teuk.getMinCouplingModeNumber()), _coupling(teuk.getCouplingCoefficient()), _scalarCoupling(teuk.getCouplingCoefficient().size() - 2 + (_lmin - std::abs(_m))), _theta(teuk.getPolarPoints()), _Slm(teuk.getPolarSolution()), _SlmP(teuk.getPolarDerivative()), _SlmPP(teuk.getPolarDerivative().size()), _r(teuk.getRadialPoints()), _Rin(teuk.getHomogeneousRadialSolution(In)), _RinP(teuk.getHomogeneousRadialDerivative(In)), _RinPP(teuk.getHomogeneousRadialDerivative(In)), _Rup(teuk.getHomogeneousRadialSolution(Up)), _RupP(teuk.getHomogeneousRadialDerivative(Up)), _RupPP(teuk.getHomogeneousRadialDerivative(Up)), _PsilmIn(teuk.getTeukolskyAmplitude(In)), _PsilmUp(teuk.getTeukolskyAmplitude(Up)){ }

HertzMode::~HertzMode(){}

void HertzMode::flipSpinWeight(){
  // flip everything to the opposite signed spin-weight
	// this assumes that the polar points are distributed symmetrically on the interval 0 to pi
	// and so parity transformations are achieved by simply reversing the polar data
	double lambdaCH = _lambda + _s*(_s + 1.);

	ComplexVector inSolutionTemp = _Rin;
	ComplexVector upSolutionTemp = _Rup;
	ComplexVector inDerivativeTemp = _RinP;
	ComplexVector upDerivativeTemp = _RupP;
	flip_spin_of_radial_teukolsky_TS(_Rin, _RinP, In, _s, _m, _a, _omega, lambdaCH, _r, inSolutionTemp, inDerivativeTemp);
	flip_spin_of_radial_teukolsky_TS(_Rup, _RupP, Up, _s, _m, _a, _omega, lambdaCH, _r, upSolutionTemp, upDerivativeTemp);

	flip_spin_of_radial_teukolsky_amplitude_TS(_PsilmIn, In, _s, _L, _m, _k, _a, _omega, lambdaCH);
	flip_spin_of_radial_teukolsky_amplitude_TS(_PsilmUp, Up, _s, _L, _m, _k, _a, _omega, lambdaCH);
	
	flip_spin_of_spheroidal_harmonic(_Slm, _SlmP, _L, _m);
	flip_spin_of_coupling_coefficients(_coupling, _L, _m);

	_lambda = flip_eigenvalue(_s, _lambda);
	_s = flip_spin(_s);
}

int HertzMode::generateSolutions(){
  if(_s == -2){
    if(_gauge == ORG || _gauge == SAAB4 || _gauge == ASAAB4){
      flipSpinWeight();
    }
  }else if(_s == 2){
    if(_gauge == IRG || _gauge == SAAB0 || _gauge == ASAAB0){
      flipSpinWeight();
    }
  }else{
    std::cout << "(HERTZ) Error: HertzMode only setup to generate Hertz potential modes from |s| = 2 Teukolsky modes. Supplied Teukolsky has spin-weight s = " << _s << ".\n";
    return 0;
  }

  double lambdaCH = _lambda + _s*(_s + 1.);

  if(_gauge == ORG){
    teukolsky_to_hertz_ORG(_PsilmIn, _PsilmUp, _PsilmIn, _PsilmUp, _L, _m, _k, _a, _omega, lambdaCH);
  }else if(_gauge == IRG){
    teukolsky_to_hertz_IRG(_PsilmIn, _PsilmUp, _PsilmIn, _PsilmUp, _L, _m, _k, _a, _omega, lambdaCH);
  }else if(_gauge == SAAB0 || _gauge == SAAB4){
    teukolsky_to_hertz_SAAB(_PsilmIn, _PsilmUp, _PsilmIn, _PsilmUp, _L, _m, _k, _a, _omega, lambdaCH);
  }else if(_gauge == ASAAB0 || _gauge == ASAAB4){
    teukolsky_to_hertz_ASAAB(_PsilmIn, _PsilmUp, _PsilmIn, _PsilmUp, _L, _m, _a, _omega, lambdaCH);
  }

  generate_radial_second_derivative(_RinPP, _s, _m, _a, _omega, _lambda, _r, _Rin, _RinP);
  generate_radial_second_derivative(_RupPP, _s, _m, _a, _omega, _lambda, _r, _Rup, _RupP);
  for(size_t i = 0; i < _Slm.size(); i++){
    _SlmPP[i] = Sslm_secondDerivative(_s, _L, _m, _a*_omega, _lambda, _theta[i], _Slm[i], _SlmP[i]);
  }
  generate_scalar_spherical_spheroidal_coupling(_scalarCoupling, _s, _lmin, _m, _coupling);

  return 1;
}

// int HertzMode::generateSolutions(){
//   // std::cout << "PsiIn = " << _PsilmIn << "\n";
//   // std::cout << "PsiUp = " << _PsilmUp << "\n";

//   // check that we have spin -2 Teukolsky solutions. Code is not yet setup for spin +2
//   if(_s == -2){
//     Complex PsiTemp;
//     if(_gauge == ORG){
//       teukolsky_to_hertz_amplitude_in(PsiTemp, _PsilmIn, _L, _m, _a, _omega, _lambda);
//       _PsilmIn = PsiTemp;
//       teukolsky_to_hertz_amplitude_up(PsiTemp, _PsilmUp, _L, _m, _a, _omega, _lambda);
//       _PsilmUp = PsiTemp;
//     }else if(_gauge == IRG){
//       teukolsky_to_hertz_amplitude_IRG_in(PsiTemp, _PsilmIn, _L, _m, _a, _omega, _lambda);
//       _PsilmIn = PsiTemp;
//       teukolsky_to_hertz_amplitude_IRG_up(PsiTemp, _PsilmUp, _L, _m, _a, _omega, _lambda);
//       _PsilmUp = PsiTemp;
//     }else if(_gauge == SAAB0){
//       teukolsky_to_hertz_amplitude_SAAB0_in(PsiTemp, _PsilmIn, _L, _m, _a, _omega, _lambda);
//       _PsilmIn = PsiTemp;
//       teukolsky_to_hertz_amplitude_SAAB0_up(PsiTemp, _PsilmUp, _L, _m, _a, _omega, _lambda);
//       _PsilmUp = PsiTemp;
//     }else if(_gauge == SAAB4){
//       teukolsky_to_hertz_amplitude_SAAB4_in(PsiTemp, _PsilmIn, _L, _m, _a, _omega, _lambda);
//       _PsilmIn = PsiTemp;
//       teukolsky_to_hertz_amplitude_SAAB4_up(PsiTemp, _PsilmUp, _L, _m, _a, _omega, _lambda);
//       _PsilmUp = PsiTemp;
//     }else if(_gauge == ASAAB0){
//       teukolsky_to_hertz_amplitude_ASAAB0_in(PsiTemp, _PsilmIn, _L, _m, _a, _omega, _lambda);
//       _PsilmIn = PsiTemp;
//       teukolsky_to_hertz_amplitude_ASAAB0_up(PsiTemp, _PsilmUp, _L, _m, _a, _omega, _lambda);
//       _PsilmUp = PsiTemp;
//     }else if(_gauge == ASAAB4){
//       teukolsky_to_hertz_amplitude_ASAAB4_in(PsiTemp, _PsilmIn, _L, _m, _a, _omega, _lambda);
//       _PsilmIn = PsiTemp;
//       teukolsky_to_hertz_amplitude_ASAAB4_up(PsiTemp, _PsilmUp, _L, _m, _a, _omega, _lambda);
//       _PsilmUp = PsiTemp;
//     }
//   }else{
//     std::cout << "(HERTZ) Error: HertzMode only setup to generate Hertz potential modes from s = -2 Teukolsky modes. Supplied Teukolsky has spin-weight s = " << _s << ".\n";
//     return 0;
//   }

//   // std::cout << "PsiIn = " << _PsilmIn << "\n";
//   // std::cout << "PsiUp = " << _PsilmUp << "\n";

//   if(_gauge == ORG || _gauge == SAAB0 || _gauge == ASAAB0){
//     _s = 2;
//     ComplexVector RinTemp = _Rin;
//     ComplexVector RinPTemp = _Rin;
//     ComplexVector RupTemp = _Rin;
//     ComplexVector RupPTemp = _Rin;

//     flip_spin_of_radial_teukolsky(RinTemp, RinPTemp, RupTemp, RupPTemp, _m, _a, _omega, _lambda, _r, _Rin, _RinP, _Rup, _RupP);

//     _Rin = RinTemp;
//     _RinP = RinPTemp;
//     _Rup = RupTemp;
//     _RupP = RupPTemp;

//     flip_spin_of_spheroidal_harmonic(_Slm, _SlmP, _L, _m);
//     // note that this assumes that theta runs from theta_min to pi - theta_min and is symmetric about pi/2

//     // std::cout << "Coupling coefficients originally given by \n";
//     // for(size_t i = 0; i < _coupling.size(); i++){
//     //   std::cout << _coupling[i] << "\n";
//     // }
//     flip_spin_of_coupling_coefficients(_coupling, _L, _m);
//     // std::cout << "Coupling coefficients now given by \n";
//     // for(size_t i = 0; i < _coupling.size(); i++){
//     //   std::cout << _coupling[i] << "\n";
//     // }

//     double laTemp;
//     flip_spin_of_spheroidal_eigenvalue(laTemp, -2, _lambda);
//     _lambda = laTemp;

//     generate_radial_second_derivative(_RinPP, 2, _m, _a, _omega, _lambda, _r, _Rin, _RinP);
//     generate_radial_second_derivative(_RupPP, 2, _m, _a, _omega, _lambda, _r, _Rup, _RupP);

//     generate_scalar_spherical_spheroidal_coupling(_scalarCoupling, 2, _lmin, _m, _coupling);
//   }else{
//     _s = -2;
//     generate_radial_second_derivative(_RinPP, -2, _m, _a, _omega, _lambda, _r, _Rin, _RinP);
//     generate_radial_second_derivative(_RupPP, -2, _m, _a, _omega, _lambda, _r, _Rup, _RupP);

//     Vector bslmo = _coupling;
//     generate_scalar_spherical_spheroidal_coupling(_scalarCoupling, -2, _lmin, _m, bslmo);
//   }
//   for(size_t i = 0; i < _Slm.size(); i++){
//     _SlmPP[i] = Sslm_secondDerivative(_s, _L, _m, _a*_omega, _lambda, _theta[i], _Slm[i], _SlmP[i]);
//   }

//   return 0;
// }

Gauge HertzMode::getGauge(){ return _gauge; }
int HertzMode::getSpinWeight(){ return _s; }
int HertzMode::getSpheroidalModeNumber(){ return _L; }
int HertzMode::getAzimuthalModeNumber(){ return _m; }
int HertzMode::getPolarModeNumber(){ return _k; }
int HertzMode::getRadialModeNumber(){ return _n; }
int HertzMode::getSampleSize(){ return _sampleSize; }
double HertzMode::getBlackHoleSpin(){ return _a; }
double HertzMode::getFrequency(){ return _omega; }
double HertzMode::getHorizonFrequency(){ return _omega - 0.5*_m*_a/(1. + sqrt(1. - _a*_a)); }
double HertzMode::getEigenvalue(){ return _lambda; }

Vector HertzMode::getCouplingCoefficient(){ return _coupling; }
double HertzMode::getCouplingCoefficient(int l){
	if(l <= getMaxCouplingModeNumber() && l >= getMinCouplingModeNumber()){
		return _coupling[l - getMinCouplingModeNumber()];
	}else{
		return 0.;
	}
}
int HertzMode::getMinCouplingModeNumber(){
	return _lmin;
}
int HertzMode::getMaxCouplingModeNumber(){
	return _coupling.size() + _lmin - 1;
}

Vector HertzMode::getScalarCouplingCoefficient(){ return _scalarCoupling; }
double HertzMode::getScalarCouplingCoefficient(int l){
	if(l <= getMaxScalarCouplingModeNumber() && l >= std::abs(_m)){
		return _scalarCoupling[l - std::abs(_m)];
	}else{
    // if(l <= getMaxScalarCouplingModeNumber()){
    //   std::cout << "Requested coupling between j = "<<_L<<" and l = "<<l<<" mode too high\n";
    //   std::cout << "Max mode coupling truncated at l = "<<getMaxScalarCouplingModeNumber()<<" with coupling coefficient "<<_scalarCoupling[getMaxScalarCouplingModeNumber() - _lmin]<<"\n";
    // }
		return 0.;
	}
}
int HertzMode::getMinScalarCouplingModeNumber(){
	return std::abs(_m);
}
int HertzMode::getMaxScalarCouplingModeNumber(){
	return _scalarCoupling.size() + std::abs(_m) - 1;
}

Vector HertzMode::getRadialPoints(){ return _r; }
ComplexVector HertzMode::getHomogeneousRadialSolution(BoundaryCondition bc){
	if(bc == In){
		return _Rin;
	}else{
		return _Rup;
	}
}
ComplexVector HertzMode::getHomogeneousRadialDerivative(BoundaryCondition bc){
	if(bc == In){
		return _RinP;
	}else{
		return _RupP;
	}
}
Complex HertzMode::getHertzAmplitude(BoundaryCondition bc){
	if(bc == In){
		return _PsilmIn;
	}else{
		return _PsilmUp;
	}
}
ComplexVector HertzMode::getRadialSolution(BoundaryCondition bc){
	ComplexVector ZR(_r.size());
	for(size_t i = 0; i < ZR.size(); i++){
		ZR[i] = getRadialSolution(bc, i);
	}
	return ZR;
}
ComplexVector HertzMode::getRadialDerivative(BoundaryCondition bc){
	ComplexVector ZR(_r.size());
	for(size_t i = 0; i < ZR.size(); i++){
		ZR[i] = getRadialDerivative(bc, i);
	}
	return ZR;
}

Vector HertzMode::getPolarPoints(){ return _theta; }
Vector HertzMode::getPolarSolution(){ return _Slm; }
Vector HertzMode::getPolarDerivative(){ return _SlmP; }

double HertzMode::getRadialPoints(int pos){ return _r[pos]; }
Complex HertzMode::getHomogeneousRadialSolution(BoundaryCondition bc, int pos){
	if(bc == In){ return _Rin[pos]; }else{ return _Rup[pos]; }
}
Complex HertzMode::getHomogeneousRadialDerivative(BoundaryCondition bc, int pos){
	if(bc == In){ return _RinP[pos]; }else{ return _RupP[pos]; }
}
Complex HertzMode::getHomogeneousRadialSecondDerivative(BoundaryCondition bc, int pos){
	if(bc == In){ return _RinPP[pos]; }else{ return _RupPP[pos]; }
}
Complex HertzMode::getHomogeneousRadialSolution(BoundaryCondition bc, int dr, int pos){
  if(dr == 0){
    return getHomogeneousRadialSolution(bc, pos);
  }else if(dr == 1){
    return getHomogeneousRadialDerivative(bc, pos);
  }else{
    return getHomogeneousRadialSecondDerivative(bc, pos);
  }
}

Complex HertzMode::getRadialSolution(BoundaryCondition bc, int pos){
	if(bc == In){ return _PsilmIn*_Rin[pos]; }else{ return _PsilmUp*_Rup[pos]; }
}
Complex HertzMode::getRadialDerivative(BoundaryCondition bc, int pos){
	if(bc == In){ return _PsilmIn*_RinP[pos]; }else{ return _PsilmUp*_RupP[pos]; }
}
Complex HertzMode::getRadialSecondDerivative(BoundaryCondition bc, int pos){
	if(bc == In){ return _PsilmIn*_RinPP[pos]; }else{ return _PsilmUp*_RupPP[pos]; }
}
Complex HertzMode::getRadialSolution(BoundaryCondition bc, int dr, int pos){
  if(dr == 0){
    return getRadialSolution(bc, pos);
  }else if(dr == 1){
    return getRadialDerivative(bc, pos);
  }else{
    return getRadialSecondDerivative(bc, pos);
  }
}

double HertzMode::getPolarPoints(int pos){ return _theta[pos]; }
double HertzMode::getPolarSolution(int pos){ return _Slm[pos]; }
double HertzMode::getPolarDerivative(int pos){ return _SlmP[pos]; }
double HertzMode::getPolarSecondDerivative(int pos){ return _SlmPP[pos]; }

void test_hertz_mode(int j, int m, int k, int n, GeodesicSource& geo){
  TeukolskyMode teuk(-2, j, m, k, n, geo);
  teuk.generateSolutions(geo);
  HertzMode hertz(teuk);
  hertz.generateSolutions();
  std::cout << "(a, omega) = (" << teuk.getBlackHoleSpin() << ", "<< teuk.getFrequency() << ")\n";

  int radialPoint = 0;
  std::cout << "R^+_{-2"<<j<<""<<m<<""<<k<<""<<n<<"}(r = "<<teuk.getRadialPoints(radialPoint)<<") = " << teuk.getHomogeneousRadialSolution(Up, radialPoint) << "\n";
  // std::cout << "dR^+_{-2"<<j<<""<<m<<""<<k<<""<<n<<"}(r = "<<teuk.getRadialPoints(radialPoint)<<") = " << teuk.getHomogeneousRadialDerivative(Up, radialPoint) << "\n";
  std::cout << "R^-_{-2"<<j<<""<<m<<""<<k<<""<<n<<"}(r = "<<teuk.getRadialPoints(radialPoint)<<") = " << teuk.getHomogeneousRadialSolution(In, radialPoint) << "\n";
  // std::cout << "dR^-_{-2"<<j<<""<<m<<""<<k<<""<<n<<"}(r = "<<teuk.getRadialPoints(radialPoint)<<") = " << teuk.getHomogeneousRadialDerivative(In, radialPoint) << "\n";
  std::cout << "R^+_{+2"<<j<<""<<m<<""<<k<<""<<n<<"}(r = "<<hertz.getRadialPoints(radialPoint)<<") = " << hertz.getHomogeneousRadialSolution(Up, radialPoint) << "\n";
  // std::cout << "dR^+_{+2"<<j<<""<<m<<""<<k<<""<<n<<"}(r = "<<hertz.getRadialPoints(radialPoint)<<") = " << hertz.getHomogeneousRadialDerivative(Up, radialPoint) << "\n";
  std::cout << "R^-_{+2"<<j<<""<<m<<""<<k<<""<<n<<"}(r = "<<hertz.getRadialPoints(radialPoint)<<") = " << hertz.getHomogeneousRadialSolution(In, radialPoint) << "\n";
  // std::cout << "dR^-_{+2"<<j<<""<<m<<""<<k<<""<<n<<"}(r = "<<hertz.getRadialPoints(radialPoint)<<") = " << hertz.getHomogeneousRadialDerivative(In, radialPoint) << "\n";

  std::cout << "C^+_{-2"<<j<<""<<m<<""<<k<<""<<n<<"} = " << teuk.getTeukolskyAmplitude(Up) << "\n";
  std::cout << "C^-_{-2"<<j<<""<<m<<""<<k<<""<<n<<"} = " << teuk.getTeukolskyAmplitude(In) << "\n";
  std::cout << "Z^+_{+2"<<j<<""<<m<<""<<k<<""<<n<<"} = " << hertz.getHertzAmplitude(Up) << "\n";
  std::cout << "Z^-_{+2"<<j<<""<<m<<""<<k<<""<<n<<"} = " << hertz.getHertzAmplitude(In) << "\n";

  std::cout << "Z^+*R^+ = " << hertz.getRadialSolution(Up, radialPoint) << "\n";
  std::cout << "Z^-*R^- = " << hertz.getRadialSolution(In, radialPoint) << "\n";

  double Slm = 0.;
  double theta = 0.4*M_PI;
  double SlmRef = Sslm(2, j, m, teuk.getBlackHoleSpin()*teuk.getFrequency(), theta);
  for(int l = std::abs(m); l <= hertz.getMaxScalarCouplingModeNumber(); l++){
    // std::cout << "A^{l="<<l<<"}_{jm} = " << hertz.getScalarCouplingCoefficient(l) << "\n";
    Slm += hertz.getScalarCouplingCoefficient(l)*Ylm(l, m, theta);
  }
  Slm *= pow(sin(theta), -2);

  std::cout << "S_{+2"<<j<<""<<m<<""<<k<<""<<n<<"}(theta = "<<theta<<") reference value = " << SlmRef << "\n";
  std::cout << "S_{+2"<<j<<""<<m<<""<<k<<""<<n<<"}(theta = "<<theta<<") coupling value = " << Slm << "\n";
  std::cout << "Fractional error = " << std::abs(1. - Slm/SlmRef) << "\n";
}

void teukolsky_to_hertz_ORG(Complex &Psi, Complex Zteuk, int L, int m, int k, double a, double omega, double lambdaCH){
  Complex plm = teukolsky_starobinsky_complex_constant(L, m, a, omega, lambdaCH);
  Complex plmk = (std::real(plm) - I*pow(-1., k)*std::imag(plm));
  Psi = pow(-1., L + m + k)*8.*Zteuk/plmk;
}

void teukolsky_to_hertz_ORG(Complex &PsiIn, Complex &PsiUp, Complex ZteukIn, Complex ZteukUp, int L, int m, int k, double a, double omega, double lambdaCH){
  Complex plm = teukolsky_starobinsky_complex_constant(L, m, a, omega, lambdaCH);
  Complex plmk = (std::real(plm) - I*pow(-1., k)*std::imag(plm));
  PsiIn = pow(-1., L + m + k)*8.*ZteukIn/plmk;
  PsiUp = pow(-1., L + m + k)*8.*ZteukUp/plmk;
}

void teukolsky_to_hertz_IRG(Complex &Psi, Complex Zteuk, int L, int m, int k, double a, double omega, double lambdaCH){
  Complex plm = teukolsky_starobinsky_complex_constant(L, m, a, omega, lambdaCH);
  Complex plmk = (std::real(plm) + I*pow(-1., k)*std::imag(plm));
  Psi = pow(-1., L + m + k)*8.*Zteuk/plmk;
}

void teukolsky_to_hertz_IRG(Complex &PsiIn, Complex &PsiUp, Complex ZteukIn, Complex ZteukUp, int L, int m, int k, double a, double omega, double lambdaCH){
  Complex plm = teukolsky_starobinsky_complex_constant(L, m, a, omega, lambdaCH);
  Complex plmk = (std::real(plm) + I*pow(-1., k)*std::imag(plm));
  PsiIn = pow(-1., L + m + k)*8.*ZteukIn/plmk;
  PsiUp = pow(-1., L + m + k)*8.*ZteukUp/plmk;
}

void teukolsky_to_hertz_SAAB(Complex &Psi, Complex Zteuk, int L, int m, int k, double a, double omega, double lambdaCH){
  double D = teukolsky_starobinsky_constant_D(m, a, omega, lambdaCH);
  Psi = pow(-1., L + m + k)*Zteuk/D;
}

void teukolsky_to_hertz_SAAB(Complex &PsiIn, Complex &PsiUp, Complex ZteukIn, Complex ZteukUp, int L, int m, int k, double a, double omega, double lambdaCH){
  double D = teukolsky_starobinsky_constant_D(m, a, omega, lambdaCH);
  PsiIn = pow(-1., L + m + k)*4.*ZteukIn/D;
  PsiUp = pow(-1., L + m + k)*4.*ZteukUp/D;
}

void teukolsky_to_hertz_ASAAB(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambdaCH){
  Psi = Zteuk/(-3.*I*omega);
}

void teukolsky_to_hertz_ASAAB(Complex &PsiIn, Complex &PsiUp, Complex ZteukIn, Complex ZteukUp, int L, int m, double a, double omega, double lambdaCH){
  PsiIn = ZteukIn/(-3.*I*omega);
  PsiUp = ZteukUp/(-3.*I*omega);
}

// void teukolsky_to_hertz_amplitude_in(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda){
//   Psi = 32.*pow(-1., L + m)*Zteuk/teukolsky_starobinsky_minus_2_in(m, a, omega, lambda);
// }

// void teukolsky_to_hertz_amplitude_up(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda){
//   Psi = 32.*pow(-1., L + m)*Zteuk/teukolsky_starobinsky_minus_2_up(m, a, omega, lambda);
// }

// void teukolsky_to_hertz_amplitude_IRG_in(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda){
//   Psi = 8.*pow(-1., L + m)*Zteuk/teukolsky_starobinsky_complex_constant(L, m, a, omega, lambda);
// }
// void teukolsky_to_hertz_amplitude_IRG_up(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda){
//   Psi = 8.*pow(-1., L + m)*Zteuk/teukolsky_starobinsky_complex_constant(L, m, a, omega, lambda);
// }

// void teukolsky_to_hertz_amplitude_SAAB0_in(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda){
//   Psi = pow(-1., L + m)*4.*std::conj(teukolsky_starobinsky_complex_constant(L, m, a, omega, lambda))/teukolsky_starobinsky_amplitude(In, -2, m, a, omega, lambda)*Zteuk/teukolsky_starobinsky_constant_D(m, a, omega, lambda);
// }
// void teukolsky_to_hertz_amplitude_SAAB0_up(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda){
//   Psi = pow(-1., L+m)*4.*std::conj(teukolsky_starobinsky_complex_constant(L, m, a, omega, lambda))/teukolsky_starobinsky_amplitude(Up, -2, m, a, omega, lambda)*Zteuk/teukolsky_starobinsky_constant_D(m, a, omega, lambda);
// }
// void teukolsky_to_hertz_amplitude_SAAB4_in(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda){
//   Psi = pow(-1., L+m)*Zteuk/teukolsky_starobinsky_constant_D(m, a, omega, lambda);
// }
// void teukolsky_to_hertz_amplitude_SAAB4_up(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda){
//   Psi = pow(-1., L+m)*Zteuk/teukolsky_starobinsky_constant_D(m, a, omega, lambda);
// }
// void teukolsky_to_hertz_amplitude_ASAAB0_in(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda){
//   Psi = 4.*std::conj(teukolsky_starobinsky_complex_constant(L, m, a, omega, lambda))/teukolsky_starobinsky_amplitude(In, -2, m, a, omega, lambda)*Zteuk/(-I*omega);
// }
// void teukolsky_to_hertz_amplitude_ASAAB0_up(Complex &Psi, Complex Zteuk, int L, int m, double a, double omega, double lambda){
//   Psi = 4.*std::conj(teukolsky_starobinsky_complex_constant(L, m, a, omega, lambda))/teukolsky_starobinsky_amplitude(Up, -2, m, a, omega, lambda)*Zteuk/(-I*omega);
// }
// void teukolsky_to_hertz_amplitude_ASAAB4_in(Complex &Psi, Complex Zteuk, int L, int m, double , double omega, double ){
//   Psi = pow(-1., L+m)*Zteuk/(-I*omega);
// }
// void teukolsky_to_hertz_amplitude_ASAAB4_up(Complex &Psi, Complex Zteuk, int L, int m, double , double omega, double ){
//   Psi = pow(-1., L+m)*Zteuk/(-I*omega);
// }
// void teukolsky_to_hertz_amplitude_ASAAB4_in(Complex &Psi, Complex Zteuk, int, int, double, double omega, double){
//   Psi = Zteuk/(-I*omega);
// }
// void teukolsky_to_hertz_amplitude_ASAAB4_up(Complex &Psi, Complex Zteuk, int, int, double, double omega, double){
//   Psi = Zteuk/(-I*omega);
// }

void flip_spin_of_spheroidal_eigenvalue(double &lambdaFlip, int s, double lambda){
  lambdaFlip = lambda + 2.*s;
}

void generate_scalar_spherical_spheroidal_coupling(Vector &Bljm, int s, int lmin, int m, Vector bljm){
  int bsize = bljm.size();
  int blmax = bsize + lmin - 1;
  int lmax = Bljm.size() + std::abs(m) - 1;
  if(lmax > blmax) lmax = blmax;

  for(int lp = std::abs(m); lp <= lmax; lp++){
    for(int i = -std::abs(s); i <= std::abs(s); i++){
      if(lp + i >= lmin){
        Bljm[lp - std::abs(m)] += bljm[lp - lmin + i]*Asljm(s, lp + i, lp, m);
      }
    }
  }
}

void generate_radial_second_derivative(ComplexVector &R2, int s, int m, double a, double omega, double lambda, Vector r, ComplexVector R0, ComplexVector R1){
  for(size_t pos = 0; pos < R2.size(); pos++){
    R2[pos] = teuk_secondDerivative(a, s, m, omega, lambda, r[pos], R0[pos], R1[pos]);
  }
}
