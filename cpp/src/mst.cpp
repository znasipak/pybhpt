// header mast.h

#include "mst.hpp"

///////////////////////////////////////////////////
// define constants, structures, and other types //
///////////////////////////////////////////////////

#define MST_N_INIT 5
#define MST_REL_ERROR 1.e-14

////////////////////////////
// Public class functions //
////////////////////////////

// MstSeriesData class
MstSeriesData::MstSeriesData(MstParameters &mstParameters): _mstParameters(mstParameters), _nMin(0), _nMax(0) {
	_a[0] = 0.0;
}
MstSeriesData::~MstSeriesData() {}

Complex MstSeriesData::getSeriesCoefficient(int n){
	if(n >= 0){
		return exp(getPositiveSeriesCoefficient(n));
	}else{
		return exp(getNegativeSeriesCoefficient(n));
	}
}
Complex MstSeriesData::getLogSeriesCoefficient(int n){
	if(n >= 0){
		return getPositiveSeriesCoefficient(n);
	}else{
		return getNegativeSeriesCoefficient(n);
	}
}
Complex MstSeriesData::getSeriesBasisFunction(BasisFunctionType type, int n, double x){
	Complex nu = _mstParameters.getRenormalizedAngularMomentum(), epsilon = _mstParameters.getMstEpsilon(),
		tau = _mstParameters.getMstTau(), kappa = _mstParameters.getMstKappa(), s = Complex(_mstParameters.getSpinWeight());
	Complex npnu = Complex(n) + nu;
	Complex z = epsilon*kappa*(1. - x);
	Result result(0., DBL_EPSILON);

	switch(type){
		case TypeI:{
			return hypergeo_2F1(npnu + 1. - I*tau, -npnu - I*tau, 1. - s - I*(epsilon + tau), x).getValue();
		}
		case DerivativeOfTypeI:{
			return dhypergeo_2F1(npnu + 1. - I*tau, -npnu - I*tau, 1. - s - I*(epsilon + tau), x).getValue();
		}
		case TypeII:{
			result = Complex(n)*log(2.*I*z) + log(phammer(nu + 1. + s - I*epsilon, n)) +
				log_hyper_u(npnu + 1. + s - I*epsilon, 2.*npnu + 2., -2.*I*z) - log(phammer(nu + 1. - s + I*epsilon, n));
			// std::cout << "Precision of MST "<<n<<" = " << result.getPrecision() << " for z = "<< z <<"\n";
			return exp(result.getValue());
		}
		case DerivativeOfTypeII:{
			result = exp(Complex(n)*log(2.*I*z) + log(phammer(nu + 1. + s - I*epsilon, n))- log(phammer(nu + 1. - s + I*epsilon, n)) +
				log(-epsilon*kappa/z) + log(Complex(n)) + log_hyper_u(npnu + 1. + s - I*epsilon, 2.*npnu + 2., -2.*I*z));
			result += exp(Complex(n)*log(2.*I*z) + log(phammer(nu + 1. + s - I*epsilon, n))- log(phammer(nu + 1. - s + I*epsilon, n)) +
				log(-epsilon*kappa/z) + log(-2.*I*z) + log_dhyper_u(npnu + 1. + s - I*epsilon, 2.*npnu + 2., -2.*I*z));
			//std::cout << "Precision of MST derivative "<<n<<" = " << result.getPrecision() << " for z = "<< z <<"\n";
			return result.getValue();
		}
		default:{
			std::cout << "Error: " << type << " is not a valid BasisFunctionType or has not yet been implemented.\n";
			return 0.;
		}
	}
}
Complex MstSeriesData::getLogSeriesBasisFunction(BasisFunctionType type, int n, double x){
	Complex nu = _mstParameters.getRenormalizedAngularMomentum(), epsilon = _mstParameters.getMstEpsilon(),
		tau = _mstParameters.getMstTau(), kappa = _mstParameters.getMstKappa(), s = Complex(_mstParameters.getSpinWeight());
	Complex npnu = Complex(n) + nu;
	Complex z = epsilon*kappa*(1. - x);

	switch(type){
		case TypeI:{
			return log(hypergeo_2F1(npnu + 1. - I*tau, -npnu - I*tau, 1. - s - I*(epsilon + tau), x).getValue());
		}
		case DerivativeOfTypeI:{
			return log(dhypergeo_2F1(npnu + 1. - I*tau, -npnu - I*tau, 1. - s - I*(epsilon + tau), x).getValue());
		}
		case TypeII:{
			Result result = Complex(n)*log(2.*I*z) + lphammer(nu + 1. + s - I*epsilon, n) +
				log_hyper_u(npnu + 1. + s - I*epsilon, 2.*npnu + 2., -2.*I*z) - lphammer(nu + 1. - s + I*epsilon, n);
			// std::cout << "Precision of MST "<<n<<" = " << result.getPrecision() << " for z = "<< z <<"\n";
			return result.getValue();
		}
		case DerivativeOfTypeII:{
			Result result1 = log(Complex(n)) + log_hyper_u(npnu + 1. + s - I*epsilon, 2.*npnu + 2., -2.*I*z);
			Result result2 = log(-2.*I*z) + log_dhyper_u(npnu + 1. + s - I*epsilon, 2.*npnu + 2., -2.*I*z);
			if(std::real(result1.getValue() - result2.getValue()) > -log(DBL_EPSILON)){
				return Complex(n)*log(2.*I*z) + log(phammer(nu + 1. + s - I*epsilon, n))- log(phammer(nu + 1. - s + I*epsilon, n)) +
					log(-epsilon*kappa/z) + result1.getValue();
			}else if(std::real(result2.getValue() - result1.getValue()) > -log(DBL_EPSILON)){
				return Complex(n)*log(2.*I*z) + log(phammer(nu + 1. + s - I*epsilon, n))- log(phammer(nu + 1. - s + I*epsilon, n)) +
					log(-epsilon*kappa/z) + result2.getValue();
			}else{
				return Complex(n)*log(2.*I*z) + log(phammer(nu + 1. + s - I*epsilon, n))- log(phammer(nu + 1. - s + I*epsilon, n)) +
					log(-epsilon*kappa/z) + (log(1. + exp(result1 - result2)) + result2).getValue();
			}
		}
		default:{
			std::cout << "Error: " << type << " is not a valid BasisFunctionType or has not yet been implemented.\n";
			return 0.;
		}
	}
}
Complex MstSeriesData::getSeriesTerm(BasisFunctionType type, int n, double x){
	return exp(getLogSeriesCoefficient(n) + getLogSeriesBasisFunction(type, n, x));
}

// MstSeriesWorkspace class

MstSeriesWorkspace::MstSeriesWorkspace(MstParameters &mstParameters):
	_mstParameters(mstParameters), _mstSeriesData(_mstParameters)
{
	if(std::abs(_mstParameters.getRenormalizedAngularMomentum()) == 0.){
		nu_solver(_mstParameters);
		if(std::abs(_mstParameters.getRenormalizedAngularMomentum()) == 0.){
			std::cout << "MST: (ERROR) Unable to solve for renormalized angular momentum for q = " << _mstParameters.getBlackHoleSpin() << ", s = "<< _mstParameters.getSpinWeight() <<", L = "
				<< _mstParameters.getSpinWeightedSpheroidalModeNumber() << " m = "<< _mstParameters.getAzimuthalModeNumber() <<", epsilon = "<< _mstParameters.getMstEpsilon() <<" \n";
		}
	}
}
MstSeriesWorkspace::MstSeriesWorkspace(double q, int s, int L, int m, double epsilon):
	_mstParameters(q, s, L, m, epsilon), _mstSeriesData(_mstParameters)
{
	if(std::abs(_mstParameters.getRenormalizedAngularMomentum()) == 0.){
		nu_solver(_mstParameters);
		if(std::abs(_mstParameters.getRenormalizedAngularMomentum()) == 0.){
			std::cout << "MST: (ERROR) Unable to solve for renormalized angular momentum for q = " << q << ", s = "<< s <<", L = "
				<< L << " m = "<< m <<", epsilon = "<< epsilon <<" \n";
		}
	}
}
MstSeriesWorkspace::MstSeriesWorkspace(double q, int s, int L, int m, double epsilon, double lambda):
	_mstParameters(q, s, L, m, epsilon, lambda), 	_mstSeriesData(_mstParameters)
{
	if(std::abs(_mstParameters.getRenormalizedAngularMomentum()) == 0.){
		nu_solver(_mstParameters);
		if(std::abs(_mstParameters.getRenormalizedAngularMomentum()) == 0.){
			std::cout << "MST: (ERROR) Unable to solve for renormalized angular momentum for q = " << q << ", s = "<< s <<", L = "
				<< L << " m = "<< m <<", epsilon = "<< epsilon <<" \n";
		}
	}

}
MstSeriesWorkspace::~MstSeriesWorkspace() {}

double MstSeriesWorkspace::getBlackHoleSpin(){ return _mstParameters.getBlackHoleSpin(); }
int MstSeriesWorkspace::getSpinWeight(){ return _mstParameters.getSpinWeight(); }
int MstSeriesWorkspace::getSpinWeightedSpheroidalModeNumber(){ return _mstParameters.getSpinWeightedSpheroidalModeNumber(); }
int MstSeriesWorkspace::getAzimuthalModeNumber(){ return _mstParameters.getAzimuthalModeNumber(); }
double MstSeriesWorkspace::getModeFrequency(){ return _mstParameters.getModeFrequency(); }
double MstSeriesWorkspace::getSpinWeightedSpheroidalEigenvalue(){ return _mstParameters.getSpinWeightedSpheroidalEigenvalue(); }

double MstSeriesWorkspace::getMstQ(){ return _mstParameters.getMstQ(); }
double MstSeriesWorkspace::getMstEpsilon(){ return _mstParameters.getMstEpsilon(); }
double MstSeriesWorkspace::getMstKappa(){ return _mstParameters.getMstKappa(); }
double MstSeriesWorkspace::getMstTau(){ return _mstParameters.getMstTau(); }
Complex MstSeriesWorkspace::getRenormalizedAngularMomentum(){ return _mstParameters.getRenormalizedAngularMomentum(); }

Result MstSeriesWorkspace::getSolution(BoundaryCondition bc, double r){
	if(bc == In){
		return inSeriesPrefactor(r)*sumMstSeries(TypeI, r);
	}else{
		return upSeriesPrefactor(r)*sumMstSeries(TypeII, r);
	}
}
Result MstSeriesWorkspace::getDerivative(BoundaryCondition bc, double r){
	if(bc == In){
		return _mstParameters.getMstDXDR()*(derivativeOfInSeriesPrefactor(r)*sumMstSeries(TypeI, r) + inSeriesPrefactor(r)*sumMstSeries(DerivativeOfTypeI, r));
	}else{
		return _mstParameters.getMstDXDR()*(derivativeOfUpSeriesPrefactor(r)*sumMstSeries(TypeII, r) + upSeriesPrefactor(r)*sumMstSeries(DerivativeOfTypeII, r));
	}
}
Complex MstSeriesWorkspace::getAmplitude(BoundaryCondition bc, Amplitude amplitude){
	if(bc == In){
		switch(amplitude){
			case Transmission:
				return inTransmissionAmplitude();

			case Incidence:
				return inIncidenceAmplitude();

			case Reflection:
				return inReflectionAmplitude();

			default:
				return inTransmissionAmplitude();
		}
	}else{
		switch(amplitude){
			case Transmission:
				return upTransmissionAmplitude();

			case Incidence:
				return upIncidenceAmplitude();

			case Reflection:
				return upReflectionAmplitude();

			default:
				return upTransmissionAmplitude();
		}
	}
}
Result MstSeriesWorkspace::getNormalizedSolution(BoundaryCondition bc, double r){
	return getSolution(bc, r)/getAmplitude(bc, Transmission);
}
Result MstSeriesWorkspace::getDerivativeOfNormalizedSolution(BoundaryCondition bc, double r){
	return getDerivative(bc, r)/getAmplitude(bc, Transmission);
}

/////////////////////////
// series coefficients //
/////////////////////////

int log_fn_coeff(int n, mst_coeffs* coeffs, const MstParameters &params){
	if(n == 0){
		coeffs->fn[0] = 0.;
		coeffs->nmax = n;
		coeffs->nmin = n;
		return 1;
	}else if(n < 0){
		return log_fn_coeff_neg(n, coeffs, params);
	}else{
		return log_fn_coeff_pos(n, coeffs, params);
	}
}

int log_fn_coeff_pos(int n, mst_coeffs* coeffs, const MstParameters &params){
	if(n > FN_SERIES_MAX/2. - 1){
		return 0;
	}
	if(coeffs->nmax < n-1){
		return 0;
	}
	Complex Rn = Rn_cf(n, params);

	coeffs->fn[n] = coeffs->fn[n - 1] + log(Rn);
	coeffs->nmax = n;

	return 1;
}

int log_fn_coeff_neg(int n, mst_coeffs* coeffs, const MstParameters &params){
	if(-n > FN_SERIES_MAX/2. - 1){
		return 0;
	}
	if(coeffs->nmin > n + 1){
		return 0;
	}

	Complex fnp1, Ln;
	Ln = Ln_cf(n, params);
	if(n == -1){
		fnp1 = coeffs->fn[0];
	}else{
		fnp1 = coeffs->fn[FN_SERIES_MAX + n + 1];
	}
	coeffs->fn[FN_SERIES_MAX + n] = fnp1 + log(Ln);
	coeffs->nmin = n;

	return 1;
}

int fn_coeff(int n, mst_coeffs* coeffs, const MstParameters &params){
	if(n == 0){
		coeffs->fn[0] = 1.;
		coeffs->nmax = n;
		coeffs->nmin = n;
		return 1;
	}else if(n < 0){
		return fn_coeff_neg(n, coeffs, params);
	}else{
		return fn_coeff_pos(n, coeffs, params);
	}
}

int fn_coeff_pos(int n, mst_coeffs* coeffs, const MstParameters &params){
	if(n > FN_SERIES_MAX/2. - 1){
		return 0;
	}
	if(coeffs->nmax < n-1){
		return 0;
	}
	Complex Rn = Rn_cf(n, params);

	coeffs->fn[n] = coeffs->fn[n - 1]*Rn;
	coeffs->nmax = n;

	return 1;
}

int fn_coeff_neg(int n, mst_coeffs* coeffs, const MstParameters &params){
	if(-n > FN_SERIES_MAX/2. - 1){
		return 0;
	}
	if(coeffs->nmin > n + 1){
		return 0;
	}

	Complex fnp1, Ln;
	Ln = Ln_cf(n, params);
	if(n == -1){
		fnp1 = coeffs->fn[0];
	}else{
		fnp1 = coeffs->fn[FN_SERIES_MAX + n + 1];
	}
	coeffs->fn[FN_SERIES_MAX + n] = fnp1*Ln;
	coeffs->nmin = n;

	return 1;
}

/////////////////////////////
// Private class functions //
/////////////////////////////

// MstSeriesData

int MstSeriesData::getNMin() const { return _nMin; }
int MstSeriesData::getNMax() const { return _nMax; }

Complex MstSeriesData::getPositiveSeriesCoefficient(int n){
	if(n > _nMax){
		generatePositiveSeriesCoefficient(n);
	}
	//std::cout << "a["<< n <<"] = " << _a[n] << "\n";
	return _a[n];
}

Complex MstSeriesData::getNegativeSeriesCoefficient(int n){
	if(-n > _nMin){
		generateNegativeSeriesCoefficient(n);
	}
	//std::cout << "a["<< n <<"] = " << _a[2*FN_SERIES_MAX + n + 1] << "\n";
	return _a[2*FN_SERIES_MAX + n + 1];
}

void MstSeriesData::generatePositiveSeriesCoefficient(int n){
	if(n - 1 > _nMax && n > 0){
		generatePositiveSeriesCoefficient(n - 1);
	}
	_a[n] = _a[n - 1] + log(Rn_cf(n, _mstParameters));
	_nMax = n;
}

void MstSeriesData::generateNegativeSeriesCoefficient(int n){
	if(n + 1 < _nMin && n < 0){
		generateNegativeSeriesCoefficient(n + 1);
	}
	if(n == -1){
		_a[2*FN_SERIES_MAX] = log(Ln_cf(n, _mstParameters));
	}else{
		_a[2*FN_SERIES_MAX + 1 + n] = _a[2*FN_SERIES_MAX + n + 2] + log(Ln_cf(n, _mstParameters));
	}
	_nMin = n;
}

// MstSeriesWorkspace

Result MstSeriesWorkspace::sumMstSeries(BasisFunctionType type, double r){
	double seriesMaxAmplitude = 0.;
	Complex pn = 0., pRef = 0., pPos = 0., pNeg = 0., pTot = 0.;
	double errorPos, errorNeg, errorTot, errorSum, error;
	double x = _mstParameters.getMstX(r);

	pn = _mstSeriesData.getSeriesTerm(type, 0, x);
	seriesMaxAmplitude = (std::abs(pn) < seriesMaxAmplitude)?seriesMaxAmplitude:std::abs(pn);
	pRef += pn;

	for(int n = 1; n < MST_N_INIT; n++){
		pn = _mstSeriesData.getSeriesTerm(type, n, x);
		seriesMaxAmplitude = (std::abs(pn) < seriesMaxAmplitude)?seriesMaxAmplitude:std::abs(pn);
		pRef += pn;

		pn = _mstSeriesData.getSeriesTerm(type, -n, x);
		seriesMaxAmplitude = (std::abs(pn) < seriesMaxAmplitude)?seriesMaxAmplitude:std::abs(pn);
		pRef += pn;
	}

	int n = MST_N_INIT;
	pn = _mstSeriesData.getSeriesTerm(type, n, x);
	seriesMaxAmplitude = (std::abs(pn) < seriesMaxAmplitude)?seriesMaxAmplitude:std::abs(pn);
	pPos += pn;
	errorPos = std::abs(pn/pRef);
	n++;
	while( n < FN_SERIES_MAX && errorPos > MST_REL_ERROR ){
		pn = _mstSeriesData.getSeriesTerm(type, n, x);
		seriesMaxAmplitude = (std::abs(pn) < seriesMaxAmplitude)?seriesMaxAmplitude:std::abs(pn);
		pPos += pn;
		errorPos = std::abs(pn/pRef);
		n++;
	}

	n = MST_N_INIT;
	pn = _mstSeriesData.getSeriesTerm(type, -n, x);
	seriesMaxAmplitude = (std::abs(pn) < seriesMaxAmplitude)?seriesMaxAmplitude:std::abs(pn);
	pNeg += pn;
	errorNeg = std::abs(pn/pRef);
	n++;
	while(n < FN_SERIES_MAX && errorNeg > MST_REL_ERROR ){
		pn = _mstSeriesData.getSeriesTerm(type, -n, x);
		seriesMaxAmplitude = (std::abs(pn) < seriesMaxAmplitude)?seriesMaxAmplitude:std::abs(pn);
		pNeg += pn;
		errorNeg = std::abs(pn/pRef);
		n++;
	}
	pTot = pRef + pPos + pNeg;
	// std::cout << "pTot = " << pTot << " \n";

	errorSum = std::abs(DBL_EPSILON*seriesMaxAmplitude/pTot);
	errorTot = (errorNeg < errorPos)?errorPos:errorNeg;
	errorTot *= std::abs(pRef/pTot);
	error = (errorSum < errorTot)?errorTot:errorSum;

	Result mst(pTot, error);

	return mst;
}

Complex MstSeriesWorkspace::inSeriesPrefactor(double r){
	Complex epsilon = _mstParameters.getMstEpsilon(), kappa = _mstParameters.getMstKappa(),
		tau = _mstParameters.getMstTau(), s = Complex(_mstParameters.getSpinWeight());
	double x = _mstParameters.getMstX(r);
	Complex value = exp(I*epsilon*kappa*x)*pow(-x, -s - I*(epsilon + tau)/2.)*pow(1. - x, I*(epsilon - tau)/2.);

	return value;
};

Complex MstSeriesWorkspace::derivativeOfInSeriesPrefactor(double r){
	Complex epsilon = _mstParameters.getMstEpsilon(), kappa = _mstParameters.getMstKappa(),
		tau = _mstParameters.getMstTau(), s = Complex(_mstParameters.getSpinWeight());
	double x = _mstParameters.getMstX(r);
	Complex value = exp(I*epsilon*kappa*x)*pow(-x, -s - I*(epsilon + tau)/2.)*pow(1. - x, I*(epsilon - tau)/2.);
	value *= I*(epsilon + tau + 2.*I*s*(x - 1.) - 2.*tau*x + 2.*epsilon*kappa*x*(x - 1.))/(2.*x*(x - 1.));

	return value;
};

Complex MstSeriesWorkspace::upSeriesPrefactor(double r){
	Complex epsilon = _mstParameters.getMstEpsilon(), kappa = _mstParameters.getMstKappa(),
		tau = _mstParameters.getMstTau(), s = Complex(_mstParameters.getSpinWeight()),
		nu = _mstParameters.getRenormalizedAngularMomentum();
	Complex z = epsilon*kappa*(1. - _mstParameters.getMstX(r));
	Complex value = pow(2., nu)*exp(-M_PI*epsilon)*exp(-I*M_PI*(nu + 1. + s));
	value *= exp(I*z)*pow(z, nu + 0.5*I*(epsilon + tau))*pow(z - epsilon*kappa, -s - 0.5*I*(epsilon + tau));

	return value;
}
Complex MstSeriesWorkspace::derivativeOfUpSeriesPrefactor(double r){
	Complex epsilon = _mstParameters.getMstEpsilon(), kappa = _mstParameters.getMstKappa(),
		tau = _mstParameters.getMstTau(), s = Complex(_mstParameters.getSpinWeight()),
		nu = _mstParameters.getRenormalizedAngularMomentum();
	Complex z = epsilon*kappa*(1. - _mstParameters.getMstX(r));
	Complex value = pow(2., nu)*exp(-M_PI*epsilon)*exp(-I*M_PI*(nu + 1. + s));
	value *= exp(I*z)*pow(z, nu + 0.5*I*(epsilon + tau))*pow(z - epsilon*kappa, -s - 0.5*I*(epsilon + tau));
	value *= epsilon*kappa*I*(epsilon*epsilon*kappa + 2.*I*z*(nu - s + I*z) + epsilon*kappa*(2.*z - 2.*I*nu + tau));
	value /= 2.*z*(z - epsilon*kappa);

	return value;
}

Complex MstSeriesWorkspace::inTransmissionAmplitude(){
	double epsilon = _mstParameters.getMstEpsilon(), kappa = _mstParameters.getMstKappa(), tau = _mstParameters.getMstTau();
	int s = _mstParameters.getSpinWeight();
	// int nmin = _mstSeriesData.getNMin(), nmax = _mstSeriesData.getNMax();
	Complex fn;

	Complex prefactors = pow(2.*kappa, 2*s)*exp(I*kappa*(epsilon + tau)
		*(1. + 2.*log(kappa)/(1. + kappa))/2.);
	Complex fSum = _mstSeriesData.getSeriesCoefficient(0);
	for(int n = 1; n < MST_N_INIT; n++){
		fSum += _mstSeriesData.getSeriesCoefficient(n);
		fSum += _mstSeriesData.getSeriesCoefficient(-n);
	}

	int n = MST_N_INIT;
	fn = _mstSeriesData.getSeriesCoefficient(n);
	fSum += fn;
	while(std::abs(fn/fSum) > MST_REL_ERROR){
		n++;
		fn = _mstSeriesData.getSeriesCoefficient(n);
		fSum += fn;
	}

	n = MST_N_INIT;
	fn = _mstSeriesData.getSeriesCoefficient(-n);
	fSum += fn;
	while(std::abs(fn/fSum) > MST_REL_ERROR){
		n++;
		fn = _mstSeriesData.getSeriesCoefficient(-n);
		fSum += fn;
	}

	return prefactors*fSum;
}
Complex MstSeriesWorkspace::inIncidenceAmplitude(){
	return 1.;
} // need to fix
Complex MstSeriesWorkspace::inReflectionAmplitude(){
	return 1.;
} // need to fix

Complex MstSeriesWorkspace::upTransmissionAmplitude(){
	Complex epsilon = _mstParameters.getMstEpsilon(), kappa = _mstParameters.getMstKappa(), nu = _mstParameters.getRenormalizedAngularMomentum();
	int s = _mstParameters.getSpinWeight();
	Complex fn;

	Complex prefactors = pow(epsilon/2., -1 - 2*s)*exp(I*epsilon*(log(epsilon) - 0.5*(1. - kappa)));
	prefactors *= pow(2., -1. - Complex(s) + I*epsilon);
	prefactors *= exp(-0.5*M_PI*I*(nu + 1. + Complex(s)))*exp(-0.5*M_PI*epsilon);
	Complex fSum = _mstSeriesData.getSeriesCoefficient(0);
	for(int n = 1; n < MST_N_INIT; n++){
		fSum += pow(-1., n)*_mstSeriesData.getSeriesCoefficient(n)*phammer(nu + 1. + Complex(s) - I*epsilon, n)/phammer(nu + 1. - Complex(s) + I*epsilon, n);
		fSum += pow(-1., -n)*_mstSeriesData.getSeriesCoefficient(-n)*phammer(nu + 1. + Complex(s) - I*epsilon, -n)/phammer(nu + 1. - Complex(s) + I*epsilon, -n);
	}

	int n = MST_N_INIT;
	fn = pow(-1., n)*_mstSeriesData.getSeriesCoefficient(n)*phammer(nu + 1. + Complex(s) - I*epsilon, n)/phammer(nu + 1. - Complex(s) + I*epsilon, n);
	fSum += fn;
	while(std::abs(fn/fSum) > MST_REL_ERROR){
		n++;
		fn = pow(-1., n)*_mstSeriesData.getSeriesCoefficient(n)*phammer(nu + 1. + Complex(s) - I*epsilon, n)/phammer(nu + 1. - Complex(s) + I*epsilon, n);
		fSum += fn;
	}

	n = MST_N_INIT;
	fn = pow(-1., -n)*_mstSeriesData.getSeriesCoefficient(-n)*phammer(nu + 1. + Complex(s) - I*epsilon, -n)/phammer(nu + 1. - Complex(s) + I*epsilon, -n);
	fSum += fn;
	while(std::abs(fn/fSum) > MST_REL_ERROR){
		n++;
		fn = pow(-1., -n)*_mstSeriesData.getSeriesCoefficient(-n)*phammer(nu + 1. + Complex(s) - I*epsilon, -n)/phammer(nu + 1. - Complex(s) + I*epsilon, -n);
		fSum += fn;
	}

	return prefactors*fSum;
} // need to fix
Complex MstSeriesWorkspace::upIncidenceAmplitude(){
	return 1.;
} // need to fix
Complex MstSeriesWorkspace::upReflectionAmplitude(){
	return 1.;
} // need to fix

///////////////////////////
// Asymptotic amplitudes //
///////////////////////////

Complex btrans(mst_coeffs coeffs, MstParameters params){
	double eps = params.getMstEpsilon(), kappa = params.getMstKappa(), tau = params.getMstTau();
	int s = params.getSpinWeight(), nmin = coeffs.nmin, nmax = coeffs.nmax;

	Complex prefactors = pow(2.*kappa,2.*s)*exp(I*kappa*(eps+tau)
		*(1. + 2.*log(kappa)/(1. + kappa))/2.);
	Complex fSumUp = 0., fSumDown = 0.;
	for(int n = 0; n <= nmax; n++){
		fSumUp += coeffs.fn[n];
	}
	for(int n = 1; n <= -nmin; n++){
		fSumDown += coeffs.fn[FN_SERIES_MAX-n];
	}
	std::cout << "btrans = " << prefactors*(fSumUp + fSumDown) << "\n";

	return prefactors*(fSumUp + fSumDown);
}
