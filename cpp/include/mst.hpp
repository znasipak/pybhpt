// header mst.h

#ifndef MST_HPP
#define MST_HPP

#include "hypergeo_f.hpp"
#include "hypergeo_u.hpp"
#include "nusolver.hpp"

///////////////////////////////////////////////////
// define constants, structures, and other types //
///////////////////////////////////////////////////

#define FN_SERIES_MAX 200

enum BasisFunctionType {TypeI, DerivativeOfTypeI, TypeII, DerivativeOfTypeII}; // TypeI is given in (120) of LRR, TypeII in (153)
enum BoundaryCondition {In, Up};
enum Amplitude {Transmission, Incidence, Reflection};

class MstSeriesData{
public:
	MstSeriesData(MstParameters &mstParameters);
	~MstSeriesData();
	
	Complex getSeriesCoefficient(int n);
	Complex getLogSeriesCoefficient(int n);
	Complex getSeriesBasisFunction(BasisFunctionType type, int n, double x);
	Complex getLogSeriesBasisFunction(BasisFunctionType type, int n, double x);
	Complex getSeriesTerm(BasisFunctionType type, int n, double x);
	
	int getNMax() const;
	int getNMin() const;
	
protected:
	Complex getPositiveSeriesCoefficient(int n);
	Complex getNegativeSeriesCoefficient(int n);
	void generatePositiveSeriesCoefficient(int n);
	void generateNegativeSeriesCoefficient(int n);
	
	MstParameters &_mstParameters;
	Complex _a[2*FN_SERIES_MAX + 1];
	int _nMin;
	int _nMax;
};

class MstSeriesWorkspace{
public:
	MstSeriesWorkspace(double q, int s, int L, int m, double eps);
	MstSeriesWorkspace(double q, int s, int L, int m, double eps, double lambda);
	MstSeriesWorkspace(MstParameters &mstParameters);
	~MstSeriesWorkspace();
	
	double getBlackHoleSpin();
	int getSpinWeight();
	int getSpinWeightedSpheroidalModeNumber();
	int getAzimuthalModeNumber();
	double getModeFrequency();
	double getSpinWeightedSpheroidalEigenvalue();
	
	double getMstQ();
	double getMstEpsilon();
	double getMstKappa();
	double getMstTau();
	Complex getRenormalizedAngularMomentum();
	
	Result getSolution(BoundaryCondition bc, double r);
	Result getDerivative(BoundaryCondition bc, double r);
	Complex getAmplitude(BoundaryCondition bc, Amplitude amp);
	Result getNormalizedSolution(BoundaryCondition bc, double r);
	Result getDerivativeOfNormalizedSolution(BoundaryCondition bc, double r);
	
protected:
	Result sumMstSeries(BasisFunctionType type, double r);
	
	Complex inSeriesPrefactor(double r);
	Complex derivativeOfInSeriesPrefactor(double r);
	Complex inTransmissionAmplitude();
	Complex inIncidenceAmplitude();
	Complex inReflectionAmplitude();
	
	Complex upSeriesPrefactor(double r);
	Complex derivativeOfUpSeriesPrefactor(double r);
	Complex upTransmissionAmplitude();
	Complex upIncidenceAmplitude();
	Complex upReflectionAmplitude();
	
	MstParameters _mstParameters;
	MstSeriesData _mstSeriesData;
};

typedef struct mst_coeffs_struct{
	Complex fn[FN_SERIES_MAX];
	int nmin;
	int nmax;
} mst_coeffs;

/////////////////////////
// series coefficients //
/////////////////////////

int log_fn_coeff(int n, mst_coeffs* coeffs, const MstParameters &params);
int log_fn_coeff_pos(int n, mst_coeffs* coeffs, const MstParameters &params);
int log_fn_coeff_neg(int n, mst_coeffs* coeffs, const MstParameters &params);

int fn_coeff(int n, mst_coeffs* coeffs, const MstParameters &params);
int fn_coeff_pos(int n, mst_coeffs* coeffs, const MstParameters &params);
int fn_coeff_neg(int n, mst_coeffs* coeffs, const MstParameters &params);

////////////////////////////////
// MST inner series solutions //
////////////////////////////////

Complex pin(double x, MstParameters params);
Complex pinP(double x, MstParameters params);

///////////////////////////
// Asymptotic amplitudes //
///////////////////////////

Complex btrans(mst_coeffs coeffs, MstParameters params);

#endif
