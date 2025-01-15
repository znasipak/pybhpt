// teukolsky.cpp

#include "teukolsky.hpp"

#define ZERO_FREQ_LIMIT 1.e-10
#define MST_AUTO_FREQ 0.1
#define SAMPLE_SIZE_INIT 256

// void polar_coupling_coefficient_general(int order, Vector &couplingVector, TeukolskyMode &teuk, int l, int m);
void flip_spin_of_radial_teukolsky_amplitude_TS(Complex &Zlm, BoundaryCondition bc, int s, int j, int m, int k, double a, double omega, double lambdaCH){
	Complex Clm = teukolsky_starobinsky_amplitude(bc, s, m, a, omega, lambdaCH);
	Complex plm = teukolsky_starobinsky_complex_constant(j, m, a, omega, lambdaCH);
	// I have not seen a reference in the literature, but this appears to depend
	// on k due to some parity-transformation that also occurs
	
	if(s == 2){
		Zlm *= 0.25*(std::real(plm) + I*pow(-1, k)*std::imag(plm))/Clm;
	}else{
		Zlm *= 4.*(std::real(plm) - I*pow(-1, k)*std::imag(plm))/Clm;
	}
}

void flip_spin_of_coupling_coefficients(Vector &bslmo, int L, int m){
  int lmin = (std::abs(m) < 2) ? 2 : std::abs(m);
  for(size_t i = 0; i < bslmo.size(); i++){
    bslmo[i] *= pow(-1, lmin + i + L);
  }
}

void flip_spin_of_spheroidal_harmonic(Vector &SlmFlip, Vector &SlmPFlip, int l, int m){
  std::reverse(SlmFlip.begin(),SlmFlip.end());
  std::reverse(SlmPFlip.begin(),SlmPFlip.end());
  for(size_t i = 0; i < SlmPFlip.size(); i++){
    SlmFlip[i] *= pow(-1., l + m);
    SlmPFlip[i] *= -pow(-1., l + m);
  }
}

TeukolskyMode::TeukolskyMode(int L, int m, int k, int n, GeodesicSource& geo): TeukolskyMode(-2, L, m, k, n, geo) {}
TeukolskyMode::TeukolskyMode(int s, int L, int m, int k, int n, GeodesicSource& geo): TeukolskyMode(s, L, m, k, n, geo.getBlackHoleSpin(), geo.getTimeFrequency(m, k, n), geo.getPolarPosition(), geo.getRadialPosition()) {}
TeukolskyMode::TeukolskyMode(int s, int L, int m, int k, int n, double a, double omega, Vector theta, Vector r):  _s(s),  _L(L), _m(m), _k(k), _n(n), _sampleSize(r.size()), _a(a), _omega(omega), _lambda(0.), _theta(theta), _Slm(theta.size(), 0.), _SlmP(theta.size(), 0.), _r(r), _Rin(r.size(), 0.), _RinP(r.size(), 0.), _Rup(r.size(), 0.), _RupP(r.size(), 0.), _ZlmIn(0.), _ZlmUp(0.) {}
TeukolskyMode::~TeukolskyMode() {}

int TeukolskyMode::generateSolutions(GeodesicSource& geo, SolutionMethod method, int samplesize){
	return generateSolutions(geo.getTimeFrequency(_m, _k, _n), geo.getTrajectoryRef(), geo.getConstantsRef(), geo.getRadialPosition(), geo.getPolarPosition(), method, samplesize);
}

int TeukolskyMode::generateSolutions(double omega, GeodesicTrajectory &traj, GeodesicConstants &geoConst, Vector r, Vector theta, SolutionMethod method, int samplesize){
	_omega = omega;
	if(std::abs(_omega) < ZERO_FREQ_LIMIT){
		// std::cout << "(TEUKOLSKY) (l, m, k, n) = ("<<_L<<", "<<_m<<", "<<_k<<", "<<_n<<") is a zero frequency mode\n";
		_omega = 0;
	}
	if(_r[0] != r[0]){
		_r = r;
	}
	if(_theta[0] != theta[0]){
		_theta = theta;
	}

	SpinWeightedHarmonic swsh(_s, _L, _m, _a*_omega, theta);
	swsh.generateSolutions();

	RadialTeukolsky teuk(_a, _s, _L, _m, _omega, swsh.getEigenvalue(), r);
	teuk.generateSolutions(method);

	return generateSolutions(swsh, teuk, traj, geoConst);
}

int TeukolskyMode::generateSolutions(SpinWeightedHarmonic& swsh, RadialTeukolsky& teuk, GeodesicTrajectory& traj, GeodesicConstants &geoConst){
	if(std::abs(_omega) < ZERO_FREQ_LIMIT){
		// std::cout << "(TEUKOLSKY) (l, m, k, n) = ("<<_L<<", "<<_m<<", "<<_k<<", "<<_n<<") is a zero frequency mode\n";
		_omega = 0;
	}
	// std::cout << "omega = " << _omega << "\n";
	// int downsampleRateR = _r.size()/_sampleSize;
	// int downsampleRateTh = _theta.size()/_sampleSize;

	// SpinWeightedHarmonic swsh(_s, _L, _m, _a*_omega, _theta);
	// swsh.generateSolutions();

	_lambda = swsh.getEigenvalue();
	_coupling = swsh.getCouplingCoefficient();
	_Slm = swsh.getSolution();
	// swsh.generateDerivatives();
	// _SlmP = swsh.getDerivative();
	// std::cout << "(TEUKOLSKY) lambda_{"<< _L <<", "<< _m <<", "<< _k <<", "<< _n <<"} = " << _lambda << "\n";
	// std::cout << "(TEUKOLSKY) S_{"<< _L <<", "<< _m <<", "<< _k <<", "<< _n <<"} = " << _Slm[0] << "\n";
	// std::cout << "(TEUKOLSKY) Sp_{"<< _L <<", "<< _m <<", "<< _k <<", "<< _n <<"} = " << _SlmP[0] << "\n";
	// std::cout << "(TEUKOLSKY) Spp_{"<< _L <<", "<< _m <<", "<< _k <<", "<< _n <<"} = " << swsh.getSecondDerivative(0) << "\n";

	// RadialTeukolsky teuk(_a, _s, _L, _m, _omega, _lambda, _r);
	// teuk.generateSolutions(AUTO);

	_Rin = teuk.getSolution(In);
	_RinP = teuk.getDerivative(In);
	// std::cout << "(TEUKOLSKY) Rin_{"<< _L <<", "<< _m <<", "<< _k <<", "<< _n <<"} = " << _Rin[0] << "\n";
	// std::cout << "(TEUKOLSKY) RinP_{"<< _L <<", "<< _m <<", "<< _k <<", "<< _n <<"} = " << _RinP[0] << "\n";
	// std::cout << "(TEUKOLSKY) RinPP_{"<< _L <<", "<< _m <<", "<< _k <<", "<< _n <<"} = " << teuk.getSecondDerivative(In)[0] << "\n";

	_Rup = teuk.getSolution(Up);
	_RupP = teuk.getDerivative(Up);
	// std::cout << "(TEUKOLSKY) Rup_{"<< _L <<", "<< _m <<", "<< _k <<", "<< _n <<"} = " << _Rup[0] << "\n";
	// std::cout << "(TEUKOLSKY) RupP_{"<< _L <<", "<< _m <<", "<< _k <<", "<< _n <<"} = " << _RupP[0] << "\n";
	// std::cout << "(TEUKOLSKY) RupPP_{"<< _L <<", "<< _m <<", "<< _k <<", "<< _n <<"} = " << teuk.getSecondDerivative(Up)[0] << "\n";

	TeukolskyAmplitudes Zlm;
	if(_s == 0){
		Zlm = scalar_amplitude(_L, _m, _k, _n, traj, geoConst, teuk, swsh);
		//
		// ComplexDerivativesMatrixStruct RinMat = {.solution = _Rin, .derivative = _RinP, .secondDerivative = _Rin};
		// ComplexDerivativesMatrixStruct RupMat = {.solution = _Rup, .derivative = _RupP, .secondDerivative = _Rup};
		// DerivativesMatrix SlmMat = {.solution = _Slm, .derivative = _SlmP, .secondDerivative = _Slm};
		// Zlm = field_amplitude(_s, _L, _m, _k, _n, geo, RinMat, RupMat, SlmMat);
	}else{
		swsh.generateDerivatives();
		_SlmP = swsh.getDerivative();

		if(geoConst.e == 0. && std::abs(geoConst.x) == 1.){
			if(std::abs(_k) > 0 || std::abs(_n) > 0){
				Zlm.in = 0.;
				Zlm.up = 0.;
			}else{
				Zlm = field_amplitude_circeq(_s, _L, _m, traj, geoConst, teuk, swsh);
			}
		}else if(geoConst.e == 0.){
			if(std::abs(_n) > 0){
				Zlm.in = 0.;
				Zlm.up = 0.;
			}else{
				Zlm = field_amplitude_sphinc(_s, _L, _m, _k, traj, geoConst, teuk, swsh);
			}
		}else if(std::abs(geoConst.x) == 1.){
			if(std::abs(_k) > 0){
				Zlm.in = 0.;
				Zlm.up = 0.;
			}else{
				Zlm = field_amplitude_ecceq(_s, _L, _m, _n, traj, geoConst, teuk, swsh);
			}
		}else{
			Zlm = field_amplitude(_s, _L, _m, _k, _n, traj, geoConst, teuk, swsh);
		}
	}

	//std::cout << "(TEUKOLSKY) Zlm^up = " << Zlm.up << "\n";
	_ZlmIn = Zlm.in;
	_ZlmUp = Zlm.up;
	_ZlmInPrecision = Zlm.inPrecision;
	_ZlmUpPrecision = Zlm.upPrecision;
	return 0;
}

int TeukolskyMode::flipSpinWeightAndFrequency(){
	// flip everything to the opposite signed spin-weight and frequency
	// this operation works for all inputs 
	double lambdaCH = _lambda + _s*(_s + 1.);

	ComplexVector inSolutionTemp = _Rin;
	ComplexVector upSolutionTemp = _Rup;
	ComplexVector inDerivativeTemp = _RinP;
	ComplexVector upDerivativeTemp = _RupP;
	flip_spin_of_radial_teukolsky_TS(_Rin, _RinP, In, _s, _m, _a, _omega, lambdaCH, _r, inSolutionTemp, inDerivativeTemp);
	flip_spin_of_radial_teukolsky_TS(_Rup, _RupP, Up, _s, _m, _a, _omega, lambdaCH, _r, upSolutionTemp, upDerivativeTemp);
	for(size_t i = 0; i < _Rin.size(); i++){
		_Rin[i] = std::conj(_Rin[i]);
		_Rup[i] = std::conj(_Rup[i]);
	}
	flip_spin_of_radial_teukolsky_amplitude_TS(_ZlmIn, In, _s, _L, _m, _k, _a, _omega, lambdaCH);
	flip_spin_of_radial_teukolsky_amplitude_TS(_ZlmUp, Up, _s, _L, _m, _k, _a, _omega, lambdaCH);
	_ZlmIn = pow(-1., _L + _k)*std::conj(_ZlmIn);
	_ZlmUp = pow(-1., _L + _k)*std::conj(_ZlmUp);
	for(size_t i = 0; i < _Slm.size(); i++){
		_Slm[i] *= pow(-1., _s + _m);
		_SlmP[i] *= pow(-1., _s + _m);
	}
	for(size_t i = 0; i < _coupling.size(); i++){
		_coupling[i] *= pow(-1., _s + _m);
	}

	_lambda = flip_eigenvalue(_s, _lambda);
	_s = flip_spin(_s);
	_omega *= -1.;
	_m *= -1.;

	return 0;
}

int TeukolskyMode::flipSpinWeight(){
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

	flip_spin_of_radial_teukolsky_amplitude_TS(_ZlmIn, In, _s, _L, _m, _k, _a, _omega, lambdaCH);
	flip_spin_of_radial_teukolsky_amplitude_TS(_ZlmUp, Up, _s, _L, _m, _k, _a, _omega, lambdaCH);
	
	flip_spin_of_spheroidal_harmonic(_Slm, _SlmP, _L, _m);
	flip_spin_of_coupling_coefficients(_coupling, _L, _m);

	_lambda = flip_eigenvalue(_s, _lambda);
	_s = flip_spin(_s);

	return 0;
}

int TeukolskyMode::getSpinWeight(){ return _s; }
int TeukolskyMode::getSpheroidalModeNumber(){ return _L; }
int TeukolskyMode::getAzimuthalModeNumber(){ return _m; }
int TeukolskyMode::getPolarModeNumber(){ return _k; }
int TeukolskyMode::getRadialModeNumber(){ return _n; }
int TeukolskyMode::getSampleSize(){ return _sampleSize; }
double TeukolskyMode::getBlackHoleSpin(){ return _a; }
double TeukolskyMode::getFrequency(){ return _omega; }
double TeukolskyMode::getHorizonFrequency(){ return _omega - 0.5*_m*_a/(1. + sqrt(1. - _a*_a)); }
double TeukolskyMode::getEigenvalue(){ return _lambda; }
Vector TeukolskyMode::getCouplingCoefficient(){ return _coupling; }
double TeukolskyMode::getCouplingCoefficient(int l){
	if(l <= getMaxCouplingModeNumber() && l >= getMinCouplingModeNumber()){
		return _coupling[l - getMinCouplingModeNumber()];
	}else{
		return 0.;
	}
}

int TeukolskyMode::getMinCouplingModeNumber(){
	return std::abs(_m) < std::abs(_s) ? std::abs(_s) : std::abs(_m);
}
int TeukolskyMode::getMaxCouplingModeNumber(){
	return _coupling.size() + getMinCouplingModeNumber() - 1;
}

Vector TeukolskyMode::getRadialPoints(){ return _r; }
ComplexVector TeukolskyMode::getHomogeneousRadialSolution(BoundaryCondition bc){
	if(bc == In){
		return _Rin;
	}else{
		return _Rup;
	}
}
ComplexVector TeukolskyMode::getHomogeneousRadialDerivative(BoundaryCondition bc){
	if(bc == In){
		return _RinP;
	}else{
		return _RupP;
	}
}
Complex TeukolskyMode::getTeukolskyAmplitude(BoundaryCondition bc){
	if(bc == In){
		return _ZlmIn;
	}else{
		return _ZlmUp;
	}
}
double TeukolskyMode::getTeukolskyAmplitudePrecision(BoundaryCondition bc){
	if(bc == In){
		return _ZlmInPrecision;
	}else{
		return _ZlmUpPrecision;
	}
}
ComplexVector TeukolskyMode::getRadialSolution(BoundaryCondition bc){
	ComplexVector ZR(_r.size());
	for(size_t i = 0; i < ZR.size(); i++){
		ZR[i] = getRadialSolution(bc, i);
	}
	return ZR;
}
ComplexVector TeukolskyMode::getRadialDerivative(BoundaryCondition bc){
	ComplexVector ZR(_r.size());
	for(size_t i = 0; i < ZR.size(); i++){
		ZR[i] = getRadialDerivative(bc, i);
	}
	return ZR;
}

Vector TeukolskyMode::getPolarPoints(){ return _theta; }
Vector TeukolskyMode::getPolarSolution(){ return _Slm; }
Vector TeukolskyMode::getPolarDerivative(){ return _SlmP; }

int TeukolskyMode::getRadialSampleNumber(){ return _r.size(); }
double TeukolskyMode::getRadialPoints(int pos){ return _r[pos]; }
Complex TeukolskyMode::getHomogeneousRadialSolution(BoundaryCondition bc, int pos){
	if(bc == In){ return _Rin[pos]; }else{ return _Rup[pos]; }
}
Complex TeukolskyMode::getHomogeneousRadialDerivative(BoundaryCondition bc, int pos){
	if(bc == In){ return _RinP[pos]; }else{ return _RupP[pos]; }
}
Complex TeukolskyMode::getRadialSolution(BoundaryCondition bc, int pos){
	if(bc == In){ return _ZlmIn*_Rin[pos]; }else{ return _ZlmUp*_Rup[pos]; }
}
Complex TeukolskyMode::getRadialDerivative(BoundaryCondition bc, int pos){
	if(bc == In){ return _ZlmIn*_RinP[pos]; }else{ return _ZlmUp*_RupP[pos]; }
}

ComplexVector TeukolskyMode::getHomogeneousSecondRadialDerivative(BoundaryCondition bc){
	ComplexVector Rpp(getHomogeneousRadialDerivative(bc));
	for(size_t i = 0; i < Rpp.size(); i++){
		Rpp[i] = teuk_secondDerivative(_a, _s, _m, _omega, _lambda, _r[i], getHomogeneousRadialSolution(bc, i), getHomogeneousRadialDerivative(bc, i));
	}
	return Rpp;
}

Complex TeukolskyMode::getHomogeneousSecondRadialDerivative(BoundaryCondition bc, int pos){
	return teuk_secondDerivative(_a, _s, _m, _omega, _lambda, _r[pos], getHomogeneousRadialSolution(bc, pos), getHomogeneousRadialDerivative(bc, pos));
}

int TeukolskyMode::getPolarSampleNumber(){ return _theta.size(); }
double TeukolskyMode::getPolarPoints(int pos){ return _theta[pos]; }
double TeukolskyMode::getPolarSolution(int pos){ return _Slm[pos]; }
double TeukolskyMode::getPolarDerivative(int pos){ return _SlmP[pos]; }
double TeukolskyMode::getPolarSecondDerivative(int pos){ 
	return Sslm_secondDerivative(_s, _L, _m, _a*_omega, _lambda, _theta[pos], _Slm[pos], _SlmP[pos]);
}

double test_teukolsky_solutions(int spin, Complex R0, Complex R1, Complex R2, double r, double a, double m, double omega, double lambda){
	double Kt = (r*r + a*a)*omega - m*a;
	double delta = r*r - 2.*r + a*a;
	double s = spin;

	Complex test = 2.*(r - 1.)*(s + 1.)/delta*R1 + ((Kt*Kt - 2.*I*s*(r - 1.)*Kt)/delta + 4.*I*s*omega*r - lambda)/delta*R0;

	return std::abs(1. + R2/test);
}

int minimum_radial_harmonic(int m, int k, GeodesicSource& geo){
	double omegaMK = geo.getTimeFrequency(m, k, 0);
	double omegaR = geo.getTimeFrequency(1);
	int n = 0;
	if(omegaMK > 0){
		while(std::abs(omegaMK + n*omegaR) > ZERO_FREQ_LIMIT && (omegaMK + n*omegaR) > 0.){
			n--;
		}
		if(std::abs(omegaMK + n*omegaR) > ZERO_FREQ_LIMIT){
			n++;
		}
	}else{
		while(std::abs(omegaMK + n*omegaR) > ZERO_FREQ_LIMIT && (omegaMK + n*omegaR) < 0.){
			n++;
		}
	}
	if(std::abs(omegaMK + n*omegaR) > ZERO_FREQ_LIMIT && omegaMK + n*omegaR < 0){
		n++;
	}

	return n;
}

int minimum_polar_harmonic_circ(int m, GeodesicSource& geo){
	double omegaPh = geo.getTimeFrequency(m, 0, 0);
	double omegaTh = geo.getTimeFrequency(2);
	int k = 0;
	while(k*omegaTh <= std::abs(omegaPh)){
		k++;
	}

	if(omegaPh > 0){
		k--;
		k *= -1;
	}

	return k;
}
