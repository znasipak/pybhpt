// specialfunc.cpp

#include "specialfunc.hpp"

void set_max_compare(double& num, double comp){
	num = ( std::abs(comp) < num ) ? num : std::abs(comp);
}

//////////////////
// Result class //
//////////////////
//Result::Result(): _realValue(0.), _imaginaryValue(0.), _realAccuracy(DBL_EPSILON), _imaginaryAccuracy(DBL_EPSILON) { }
//Result::Result(Complex value): _realValue(std::real(value)), _imaginaryValue(std::imag(value)), _realAccuracy(DBL_EPSILON), _imaginaryAccuracy(DBL_EPSILON) { }
Result::Result(Complex value, Complex precision): _realValue(std::real(value)), _imaginaryValue(std::imag(value)), _realAccuracy(std::abs(std::real(precision))), _imaginaryAccuracy(std::abs(std::imag(precision))) {}

Result::~Result() {}

// define behavior under addition
Result& Result::operator+=(const Result& rhs){
	// add real parts
	_realValue += rhs.getRealValue();
	if( _realAccuracy > rhs.getRealAccuracy() ){
		_realAccuracy = _realAccuracy*sqrt(1. + pow(rhs.getRealAccuracy()/_realAccuracy, 2));
	}else{
		_realAccuracy = rhs.getRealAccuracy()*sqrt(1. + pow(_realAccuracy/rhs.getRealAccuracy(), 2));
	}

	// add imag parts
	_imaginaryValue += rhs.getImaginaryValue();
	if( _imaginaryAccuracy > rhs.getImaginaryAccuracy() ){
		_imaginaryAccuracy = _imaginaryAccuracy*sqrt(1. + pow(rhs.getImaginaryAccuracy()/_imaginaryAccuracy, 2));
	}else{
		_imaginaryAccuracy = rhs.getImaginaryAccuracy()*sqrt(1. + pow(_imaginaryAccuracy/rhs.getImaginaryAccuracy(), 2));
	}

	return *this;
}
// Result& Result::operator+=(const Complex& rhs){
// 	_realValue += std::real(rhs);
// 	_imaginaryValue += std::imag(rhs);
//
// 	return *this;
// }
Result& Result::operator+=(const Complex& rhs){
	*this += Result(rhs, rhs*DBL_EPSILON);

	return *this;
}

// define behavior under subtraction
Result& Result::operator-=(const Result& rhs){
	// subtract real parts
	_realValue -= rhs.getRealValue();
	_realAccuracy = sqrt(pow(_realAccuracy, 2) + pow(rhs.getRealAccuracy(), 2));

	// subtract imag parts
	_imaginaryValue -= rhs.getImaginaryValue();
	_imaginaryAccuracy = sqrt(pow(_imaginaryAccuracy, 2) + pow(rhs.getImaginaryAccuracy(), 2));

	return *this;
}
// Result& Result::operator-=(const Complex& rhs){
// 	_realValue -= std::real(rhs);
// 	_imaginaryValue -= std::imag(rhs);
//
// 	return *this;
// }
Result& Result::operator-=(const Complex& rhs){
	*this -= Result(rhs, rhs*DBL_EPSILON);

	return *this;
}

// define behavior under multiplication
Result& Result::operator*=(const Result& rhs){
	double temp = sqrt(pow(rhs.getRealValue(), 2)
		+ pow(rhs.getImaginaryValue()*_imaginaryAccuracy, 2)
		+ pow(rhs.getRealAccuracy()*_realValue, 2)
		+ pow(rhs.getImaginaryAccuracy()*_imaginaryValue, 2));
	_imaginaryAccuracy = sqrt(pow(rhs.getImaginaryValue()*_realAccuracy, 2)
		+ pow(rhs.getRealValue()*_imaginaryAccuracy, 2)
		+ pow(rhs.getImaginaryAccuracy()*_realValue, 2)
		+ pow(rhs.getRealAccuracy()*_imaginaryValue, 2));
	_realAccuracy = temp;

	temp = _realValue*rhs.getRealValue() - _imaginaryValue*rhs.getImaginaryValue();
	_imaginaryValue = _imaginaryValue*rhs.getRealValue() + _realValue*rhs.getImaginaryValue();
	_realValue = temp;

	return *this;
}
Result& Result::operator*=(const Complex& rhs){
	double temp = _realAccuracy*sqrt(pow(std::real(rhs), 2)
		+ pow(std::imag(rhs)*_imaginaryAccuracy/_realAccuracy, 2));
	_imaginaryAccuracy = _realAccuracy*sqrt(pow(std::imag(rhs), 2)
		+ pow(std::real(rhs)*_imaginaryAccuracy/_realAccuracy, 2));
	_realAccuracy = temp;

	temp = _realValue*std::real(rhs) - _imaginaryValue*std::imag(rhs);
	_imaginaryValue = _imaginaryValue*std::real(rhs) + _realValue*std::imag(rhs);
	_realValue = temp;

	return *this;
}
// Result& Result::operator*=(const Complex& rhs){
// 	*this *= Result(rhs, rhs*DBL_EPSILON);
//
// 	return *this;
// }

// define behavior under division
Result& Result::operator/=(const Result& rhs){
	double magnitude = std::abs(rhs.getValue());
	double temp = sqrt(pow(rhs.getRealValue()*_realAccuracy/magnitude, 2)
		+ pow(rhs.getImaginaryValue()*_imaginaryAccuracy/magnitude, 2)
		+ pow(rhs.getRealAccuracy()*rhs.getImaginaryValue()*(rhs.getImaginaryValue()*_realValue - rhs.getRealValue()*_imaginaryValue)/pow(magnitude, 3), 2)
		+ pow(rhs.getImaginaryAccuracy()*rhs.getRealValue()*(rhs.getImaginaryValue()*_realValue - rhs.getRealValue()*_imaginaryValue)/pow(magnitude, 3), 2));
	_imaginaryAccuracy = sqrt(pow(rhs.getImaginaryValue()*_realAccuracy/magnitude, 2)
		+ pow(rhs.getRealValue()*_imaginaryAccuracy/magnitude, 2)
		+ pow(rhs.getRealAccuracy()*rhs.getImaginaryValue()*(rhs.getImaginaryValue()*_imaginaryValue + rhs.getRealValue()*_realValue)/pow(magnitude, 3), 2)
		+ pow(rhs.getImaginaryAccuracy()*rhs.getRealValue()*(rhs.getImaginaryValue()*_imaginaryValue + rhs.getRealValue()*_realValue)/pow(magnitude, 3), 2));
	_realAccuracy = temp;

	temp = (_realValue*rhs.getRealValue() + _imaginaryValue*rhs.getImaginaryValue());
	_imaginaryValue = (_imaginaryValue*rhs.getRealValue() - _realValue*rhs.getImaginaryValue())/pow(magnitude, 2);
	_realValue = temp/pow(magnitude, 2);

	return *this;
}
//convert 1/(x + iy) => (x - iy)/sqrt(x^2 + y^2)
Result& Result::operator/=(const Complex& rhs){
	*this *= (std::real(rhs) - I*std::imag(rhs))/pow(std::abs(rhs), 2);
	return *this;
}
// Result& Result::operator/=(const Complex& rhs){
// 	*this *= Result(1./rhs, 1./rhs*DBL_EPSILON);
//
// 	return *this;
// }

Complex Result::getValue() const { return _realValue + I*_imaginaryValue; }
double Result::getRealValue() const { return _realValue; }
double Result::getImaginaryValue() const { return _imaginaryValue; }
//double Result::getMagnitude() const { return sqrt(_realValue*_realValue + _imaginaryValue*_imaginaryValue); }
//double Result::getPhase() const { return arg(_realValue + I*_imaginaryValue); }

Complex Result::getAccuracy() const { return _realAccuracy + I*_imaginaryAccuracy; }
double Result::getRealAccuracy() const { return _realAccuracy; }
double Result::getImaginaryAccuracy() const { return _imaginaryAccuracy; }
Complex Result::getPrecision() const {
	if(_realValue == 0 && _imaginaryValue == 0){
		return _realAccuracy + I*_imaginaryAccuracy;
	}else if(_realValue == 0){
		return _realAccuracy + I*_imaginaryAccuracy/_imaginaryValue;
	}else if(_imaginaryValue == 0){
		return _realAccuracy/_realValue + I*_imaginaryAccuracy;
	}else{
		return _realAccuracy/_realValue + I*_imaginaryAccuracy/_imaginaryValue;
	}
}

void Result::setValue(Complex value){ _realValue = std::real(value); _imaginaryValue = std::imag(value); }
void Result::setAccuracy(Complex accuracy){ _realAccuracy = std::real(accuracy); _imaginaryAccuracy = std::imag(accuracy); }

std::ostream& operator<<(std::ostream &out, const Result& result){
	out << result.getValue();
	out << " +/- " << result.getRealAccuracy() + I*result.getImaginaryAccuracy();
	return out;
}

Result operator+(Result lhs, const Result& rhs){
	lhs += rhs;
	return lhs;
}
Result operator-(Result lhs, const Result& rhs){
	lhs -= rhs;
	return lhs;
}
Result operator*(Result lhs, const Result& rhs){
	lhs *= rhs;
	return lhs;
}
Result operator/(Result lhs, const Result& rhs){
	lhs /= rhs;
	return lhs;
}
bool operator==(Result lhs, const Result& rhs){
	lhs -= rhs;
	return (std::abs(lhs.getRealValue()) <= lhs.getRealAccuracy()) && (std::abs(lhs.getImaginaryValue()) <= lhs.getImaginaryAccuracy());
}

Result operator+(Result lhs, const Complex& rhs){
	lhs += rhs;
	return lhs;
}
Result operator+(const Complex& lhs, Result rhs){
	rhs += lhs;
	return rhs;
}
Result operator-(Result lhs, const Complex& rhs){
	lhs -= rhs;
	return lhs;
}
Result operator-(const Complex& lhs, Result rhs){
	rhs -= lhs;
	rhs *= -1.;
	return rhs;
}
Result operator*(Result lhs, const Complex& rhs){
	lhs *= rhs;
	return lhs;
}
Result operator*(const Complex& lhs, Result rhs){
	rhs *= lhs;
	return rhs;
}
Result operator/(Result lhs, const Complex& rhs){
	lhs /= rhs;
	return lhs;
}
Result operator/(const Complex& lhs, Result rhs){
	rhs /= (rhs*rhs);
	rhs *= lhs;
	return rhs;
}

////////////////////////////////////
// abs, real, and imag of Results //
////////////////////////////////////

Result abs(Result result){
	return Result(std::abs(result.getValue()), std::abs(result.getAccuracy()));
}

Result real(Result result){
	return Result(result.getRealValue(), result.getRealAccuracy());
}

Result imag(Result result){
	return Result(result.getImaginaryValue(), result.getImaginaryAccuracy());
}

Result conjugate(Result result){
	Result conj = result;
	conj.setValue(result.getRealValue() - I*result.getImaginaryValue());
	return conj;
}

/////////////////////////////////////
// Log and Exp of Result functions //
/////////////////////////////////////

// this will only track the accuracy of the magnitude
Result exp(const Result& result){
	return Result(exp(result.getValue()), result.getAccuracy()*std::abs(exp(result.getValue())));
}

// this will only track the accuracy of the magnitude
Result log(const Result& result){
	return Result(log(result.getValue()), result.getAccuracy()/std::abs(result.getValue()));
}

Result pow(const Result& result, const Complex a){
	if( std::abs(a) == 0. ){
		return Result(1., 0);
	}else if( std::real(a) < 1. ){
		return pow(result, a + 1.)/result;
	}else if( std::real(a) > 1. ){
		return pow(result, a - 1.)*result;
	}
	return Result(pow(result.getValue(), a), a*result.getAccuracy()*pow(result.getValue(), a - 1.));
}

Result pow(const Result& result, const int a){
	if( a == 0 ){
		return Result(1., 0);
	}else if( a < 0 ){
		return pow(result, a + 1)/result;
	}else{
		return pow(result, a - 1)*result;
	}
	return pow(result, Complex(a));
}

Result sqrt(const Result& result){
	return pow(result, Complex(0.5));
}

/////////////////////
// Gamma functions //
/////////////////////

int is_integer(double x){
	double iter = std::abs(x);
	while(iter > 0.){
		iter -= 1;
	}
	if(iter == 0.){
		return 1;
	}else{
		return 0;
	}
}

int is_negative_integer(double x){
	return is_integer(x) && (x <= 0);
}

int is_negative_integer(Complex x){
	return is_integer(std::real(x)) && (std::real(x) <= 0) && std::abs(std::imag(x)) == 0.;
}

Complex lgamma(Complex z){

	if(is_negative_integer(z)){
		// if evaluated at pole, return residue
		return log(pow(-1, -z)/factorial(-std::real(z)));
	}

	gsl_sf_result lnr;
	gsl_sf_result arg;
	int error;

	error = gsl_sf_lngamma_complex_e(real(z), imag(z), &lnr, &arg);
	if(error == GSL_ELOSS){
		printf("Error\n");
	}

	if(isnan(std::abs(lnr.val)) && std::abs(z) < 1.e-5){
		return log(1./(z + 0.5772156649015329*z*z - 0.6558780715202539*z*z*z));
	}

	if(std::real(z) < 0 && isnan(std::abs(lnr.val))){
		return log(M_PI) - log(sin(M_PI*z)) - lgamma(1. - z);
	}

	// if(isnan(std::abs(lnr.val))){
	// 	std::cout << "Gamma("<<z<<") is nan \n";
	// }

	return I*arg.val+lnr.val;
}

Complex cgamma(Complex z){
	return exp(lgamma(z));
}

//factorial function
int factorial(int n){
	return (n < 1) ? 1 : factorial(n - 1) * n;
}
double factorial(double n){
	return (n < 1) ? 1 : factorial(n - 1) * n;
}

// chi(n) function
double chi(int n){
	return sqrt(M_PI)*exp(lgamma(0.5*n + 1.) - lgamma(0.5*n + 0.5));
}

///////////////////////
// Pochhammer symbol //
///////////////////////

Complex phammer(Complex a, int n)
{
	if( n < 0 ){
		return pow(-1., n)/phammer(1. - a, -n);
	}else if( n == 0 ){
		return 1;
	}else{
		return (a + Complex(n) - 1.)*phammer(a, n-1);
	}
}

Complex lphammer(Complex a, int n)
{
	if( n < 0 ){
		return Complex(n)*M_PI*I - lphammer(1. - a, -n);
	}else if( n == 0 ){
		return 0.;
	}else{
		return log(a + Complex(n) - 1.) + lphammer(a, n-1);
	}
}

Complex times_phammer(Complex f, Complex a, int n)
{
	Complex am1(0, 0);

	if( n < 0 ){
		am1 = 0;
	}else if( n == 0 ){
		am1 = f;
	}else{
		am1 = (a + std::complex<double>(n) - 1.)*times_phammer(f, a, n-1);
	}

	return am1;
}

Complex phammer_ratio(Complex f, Complex a, int n)
{
	Complex am1(0, 0);

	if( n < 0 ){
		am1 = 0;
	}else if( n == 0 ){
		am1 = 1;
	}else{
		am1 = (f + std::complex<double>(n) - 1.)/(a + std::complex<double>(n) - 1.)*phammer_ratio(f, a, n-1);
	}

	return am1;
}

////////////////////////////////
// Scalar spherical harmonics //
////////////////////////////////

// Coupling between scalar and derivatives of spin-weighted harmonics
double clm(const int &l, const int &m){
	if( std::abs(m) > l ) return 0.;
	// return sqrt(double(l*l - m*m)/(2.*l + 1.)/(2.*l - 1.));
	return sqrt(double(l - m)/(2.*l + 1.))*sqrt(double(l + m)/(2.*l - 1.));
}

Complex Ylm(const int &l, const int &m, const double &th, const double &ph){
	return Ylm(l, m, th)*exp(I*Complex(m)*ph);
}

double Ylm(const int &l, const int &m, const double &th){
	if( m < 0 && l >= std::abs(m) ){
		return pow(-1, m)*gsl_sf_legendre_sphPlm(l, -m, cos(th));
	}else if(l >= std::abs(m)){
		return gsl_sf_legendre_sphPlm(l, m, cos(th));
	}else{
		return 0.;
	}
}

double Ylm_derivative(const int &l, const int &m, const double &th){
	return (l*clm(l + 1, m)*Ylm(l + 1, m, th) - (l + 1)*clm(l, m)*Ylm(l - 1, m, th))/sin(th);
}

////////////////////
// Inverse cosine //
////////////////////

Complex cacos(Complex z){
	if(std::real(z) < 0){
		return M_PI - cacos(-z);
	}

	return M_PI/2. + I*log( I*z + I*sqrt(z - 1.)*sqrt(1. + z) );
}

Complex cacos(double z){ return cacos(Complex(z)); }

double cot(double th){
	return cos(th)/sin(th);
}


////////////////////////
// Elliptic integrals //
////////////////////////

double elliptic_k(double const &k){
	return boost::math::ellint_1(k);
}

double elliptic_k(double const &phi, double const &k){
	return boost::math::ellint_1(k, phi);
}

double elliptic_e(double const &k){
	return boost::math::ellint_2(k);
}

double elliptic_e(double const &phi, double const &k){
	return boost::math::ellint_2(k, phi);
}

double elliptic_pi(double const &n, double const &k){
	return boost::math::ellint_3(k, n);
}

double elliptic_pi(double const &n, double const &phi, double const &k){
	return boost::math::ellint_3(k, n, phi);
}

double jacobi_sn(double const &u, double const &k){
	return boost::math::jacobi_sn(k, u);
}

double jacobi_amp(double const &u, double const &k){
	return asin(boost::math::jacobi_sn(k, u));
}
