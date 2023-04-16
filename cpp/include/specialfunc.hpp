// specialfunc.hpp

#ifndef SF_HPP
#define SF_HPP

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include "utils.hpp"

// simple signnum function for all ordered argument types
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void set_max_compare(double& num, double comp);

//////////////////
// Result class //
//////////////////

// The Result class is an error tracking class for special functions that have been computed with some amount of numerical error that may be greater than machine epsilon
// All error is propogated via standard propogation of error formulae
class Result{
public:
	Result();
	Result(Complex value);
	Result(Complex value, Complex error);
  	~Result();
	Result& operator+=(const Result& rhs);
	Result& operator+=(const Complex& rhs);
	Result& operator-=(const Result& rhs);
	Result& operator-=(const Complex& rhs);
	Result& operator*=(const Result& rhs);
	Result& operator*=(const Complex& rhs);
	Result& operator/=(const Result& rhs);
	Result& operator/=(const Complex& rhs);

	Complex getValue() const;
	double getRealValue() const;
	double getImaginaryValue() const;
	Complex getAccuracy() const;
	Complex getPrecision() const;
	double getRealAccuracy() const;
	double getImaginaryAccuracy() const;

	void setValue(Complex value);
	void setAccuracy(Complex error);

protected:
	double _realValue;
	double _imaginaryValue;
	double _realAccuracy;
	double _imaginaryAccuracy;
};

std::ostream& operator<<(std::ostream &out, const Result& result);

Result operator+(Result lhs, const Result& rhs);
Result operator-(Result lhs, const Result& rhs);
Result operator*(Result lhs, const Result& rhs);
Result operator/(Result lhs, const Result& rhs);
bool operator==(Result lhs, const Result& rhs);

Result operator+(Result lhs, const Complex& rhs);
Result operator+(const Complex& lhs, Result rhs);
Result operator-(Result lhs, const Complex& rhs);
Result operator-(const Complex& lhs, Result rhs);
Result operator*(Result lhs, const Complex& rhs);
Result operator*(const Complex& lhs, Result rhs);
Result operator/(Result lhs, const Complex& rhs);
Result operator/(const Complex& lhs, Result rhs);

////////////////////////////////////
// abs, real, and imag of Results //
////////////////////////////////////

Result abs(Result result);
Result real(Result result);
Result imag(Result result);
Result conjugate(Result result);

/////////////////////////////////////
// Log and Exp of Result functions //
/////////////////////////////////////

Result exp(const Result& result);
Result log(const Result& result);
Result pow(const Result& result, const Complex a);
Result pow(const Result& result, const int a);
Result sqrt(const Result& result);

/////////////////////
// Gamma functions //
/////////////////////

int is_integer(double x);
int is_negative_integer(double x);
int is_negative_integer(Complex x);

Complex lgamma(Complex z);
Complex cgamma(Complex z);
int factorial(int n);
double factorial(double n);
double chi(int n);

///////////////////////
// Pochhammer symbol //
///////////////////////

Complex phammer(Complex a, int n);
Complex lphammer(Complex a, int n);
Complex times_phammer(Complex f, Complex a, int n);
Complex phammer_ratio(Complex f, Complex a, int n);

////////////////////////////////
// Scalar spherical harmonics //
////////////////////////////////

// Coupling between scalar and derivatives of spin-weighted harmonics
double clm(const int &l, const int &m);

Complex Ylm(const int &l, const int &m, const double &th, const double &ph);
double Ylm(const int &l, const int &m, const double &th);
double Ylm_derivative(const int &l, const int &m, const double &th);

////////////////////
// Inverse cosine //
////////////////////

Complex cacos(Complex f);
Complex cacos(double f);

double cot(double theta);

////////////////////////
// Elliptic integrals //
////////////////////////

double elliptic_k(double const &k);
double elliptic_k(double const &phi, double const &k);
double elliptic_e(double const &k);
double elliptic_e(double const &phi, double const &k);
double elliptic_pi(double const &n, double const &k);
double elliptic_pi(double const &n, double const &phi, double const &k);
double jacobi_sn(double const &u, double const &k);
double jacobi_amp(double const &u, double const &k);

#endif
