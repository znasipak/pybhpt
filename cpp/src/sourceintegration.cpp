// sourceintegration.cpp

#include "sourceintegration.hpp"

#define PRECISION_THRESHOLD 1.e-2

// internal versions of the source-term coefficient generators that take the
// (mode-constant) frequency omega directly rather than recomputing it from
// the orbital frequencies on every call
static void A_coeffs_w(Complex &Ann0, Complex &Anmbar0, Complex &Ambarmbar0, Complex &Anmbar1, Complex &Ambarmbar1, Complex &Ambarmbar2, int const &m, double const &omega, GeodesicConstants &geoConstants, double const &rp, double const &thp, double const &Slm, double const &SlmP, double const &SlmPP);
static void A13_coeffs_w(Complex &All0, Complex &Alm0, Complex &Amm0, Complex &Alm1, Complex &Amm1, Complex &Amm2, int const &m, double const &omega, GeodesicConstants &geoConstants, double const &rp, double const &thp, double const &Slm, double const &SlmP, double const &SlmPP);

SummationHelper::SummationHelper(): _sum(0.), _previousSum(0.), _maxTerm(0.), _error(1.), _basePrecision(DBL_EPSILON) {}
SummationHelper::~SummationHelper() {}

void SummationHelper::add(Complex val){
	// std::cout << "val = " << val << "\n";
	if(std::abs(val) > _maxTerm) _maxTerm = std::abs(val);
	_previousSum = _sum;
	_sum += val;

	_error = _basePrecision*_maxTerm; // conservative estimate of absolute error
}

void SummationHelper::setBasePrecision(double val){
	_basePrecision = val;
}

Complex SummationHelper::getSum(){
	return _sum;
}

double SummationHelper::getMaxTerm(){
	return _maxTerm;
}

double SummationHelper::getError(){
	// provides some measure of cancellation error
	return _error;
}

double SummationHelper::getPrecision(){
	// provides some measure of cancellation error
	return getError()/std::abs(_sum);
}

int integrand_convergence(Complex old_value, Complex new_value, double eps1, double eps2){
	double eps0 = std::abs(1. - old_value/new_value);
	return ((eps0 < eps1) || (eps0 < eps2));
}

///////////////////////////////////////////////////////
// 					Field amplitudes				 //
///////////////////////////////////////////////////////

TeukolskyAmplitudes field_amplitude(int s, int L, int m, int k, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConstants, RadialTeukolsky &teuk, SpinWeightedHarmonic &swsh){
	if( std::abs(s) == 2 ){
		ComplexDerivativesMatrixStruct Rin = {.solution = teuk.getSolution(In), .derivative = teuk.getDerivative(In), .secondDerivative = teuk.getSecondDerivative(In)};
		ComplexDerivativesMatrixStruct Rup = {.solution = teuk.getSolution(Up), .derivative = teuk.getDerivative(Up), .secondDerivative = teuk.getSecondDerivative(Up)};
		DerivativesMatrix Slm = {.solution = swsh.getSolution(), .derivative = swsh.getDerivative(), .secondDerivative = swsh.getSecondDerivative()};
		return teukolsky_amplitude(s, L, m, k, n, traj, geoConstants, Rin, Rup, Slm);
	}else if( s == 0 ){
		return scalar_amplitude_generic(L, m, k, n, traj, geoConstants, teuk, swsh);
	}else{
		std::cout << "SOURCEINTEGRATION: ERROR: Source integration not yet implemented for s = " << s << " fields \n";
		TeukolskyAmplitudes Zlm = {0., 0., DBL_EPSILON, DBL_EPSILON};
		return Zlm;
	}
}

TeukolskyAmplitudes field_amplitude_circeq(int s, int L, int m, GeodesicTrajectory& traj, GeodesicConstants &geoConstants, RadialTeukolsky &teuk, SpinWeightedHarmonic &swsh){
	if( std::abs(s) == 2 ){
		ComplexDerivativesMatrixStruct Rin = {.solution = teuk.getSolution(In), .derivative = teuk.getDerivative(In), .secondDerivative = teuk.getSecondDerivative(In)};
		ComplexDerivativesMatrixStruct Rup = {.solution = teuk.getSolution(Up), .derivative = teuk.getDerivative(Up), .secondDerivative = teuk.getSecondDerivative(Up)};
		DerivativesMatrix Slm = {.solution = swsh.getSolution(), .derivative = swsh.getDerivative(), .secondDerivative = swsh.getSecondDerivative()};
		return teukolsky_amplitude_circeq(s, L, m, traj, geoConstants, Rin, Rup, Slm);
	}else if( s == 0 ){
		return scalar_amplitude_circular(L, m, 0, 0, traj, geoConstants, teuk, swsh);
	}else{
		std::cout << "SOURCEINTEGRATION: ERROR: Source integration not yet implemented for s = " << s << " fields \n";
		TeukolskyAmplitudes Zlm = {0., 0., DBL_EPSILON, DBL_EPSILON};
		return Zlm;
	}
}

TeukolskyAmplitudes field_amplitude_ecceq(int s, int L, int m, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConstants, RadialTeukolsky &teuk, SpinWeightedHarmonic &swsh){
	if( std::abs(s) == 2 ){
		ComplexDerivativesMatrixStruct Rin = {.solution = teuk.getSolution(In), .derivative = teuk.getDerivative(In), .secondDerivative = teuk.getSecondDerivative(In)};
		ComplexDerivativesMatrixStruct Rup = {.solution = teuk.getSolution(Up), .derivative = teuk.getDerivative(Up), .secondDerivative = teuk.getSecondDerivative(Up)};
		DerivativesMatrix Slm = {.solution = swsh.getSolution(), .derivative = swsh.getDerivative(), .secondDerivative = swsh.getSecondDerivative()};
		return teukolsky_amplitude_ecceq(s, L, m, n, traj, geoConstants, Rin, Rup, Slm);
	}else if( s == 0 ){
		return scalar_amplitude_equatorial(L, m, 0, n, traj, geoConstants, teuk, swsh);
	}else{
		std::cout << "SOURCEINTEGRATION: ERROR: Source integration not yet implemented for s = " << s << " fields \n";
		TeukolskyAmplitudes Zlm = {0., 0., DBL_EPSILON, DBL_EPSILON};
		return Zlm;
	}
}

TeukolskyAmplitudes field_amplitude_sphinc(int s, int L, int m, int k, GeodesicTrajectory& traj, GeodesicConstants &geoConstants, RadialTeukolsky &teuk, SpinWeightedHarmonic &swsh){
	if( std::abs(s) == 2 ){
		ComplexDerivativesMatrixStruct Rin = {.solution = teuk.getSolution(In), .derivative = teuk.getDerivative(In), .secondDerivative = teuk.getSecondDerivative(In)};
		ComplexDerivativesMatrixStruct Rup = {.solution = teuk.getSolution(Up), .derivative = teuk.getDerivative(Up), .secondDerivative = teuk.getSecondDerivative(Up)};
		DerivativesMatrix Slm = {.solution = swsh.getSolution(), .derivative = swsh.getDerivative(), .secondDerivative = swsh.getSecondDerivative()};
		return teukolsky_amplitude_sphinc(s, L, m, k, traj, geoConstants, Rin, Rup, Slm);
	}else if( s == 0 ){
		return scalar_amplitude_spherical(L, m, k, 0, traj, geoConstants, teuk, swsh);
	}else{
		std::cout << "SOURCEINTEGRATION: ERROR: Source integration not yet implemented for s = " << s << " fields \n";
		TeukolskyAmplitudes Zlm = {0., 0., DBL_EPSILON, DBL_EPSILON};
		return Zlm;
	}
}

///////////////////////////////////////
// Teukolsky s = -2 field amplitudes //
///////////////////////////////////////

TeukolskyAmplitudes teukolsky_amplitude(int s, int L, int m, int k, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConstants, const ComplexDerivativesMatrixStruct &Rin, const ComplexDerivativesMatrixStruct &Rup, const DerivativesMatrix &Slm){
	std::function<void(Complex &, Complex &, int const &, int const &, int const &, int const &, GeodesicConstants &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, Complex const &, Complex const &, Complex const &,  Complex const &, Complex const &, Complex const &, double const &, double const &, double const &)> integrand;
	std::function<void(Complex &, Complex &, int const &, int const &, int const &, int const &, GeodesicConstants &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, Complex const &, Complex const &, Complex const &,  Complex const &, Complex const &, Complex const &, double const &, double const &, double const &)> integrandRadialTurningPoint;
	std::function<void(Complex &, Complex &, int const &, int const &, int const &, int const &, GeodesicConstants &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, Complex const &, Complex const &, Complex const &,  Complex const &, Complex const &, Complex const &, double const &, double const &, double const &)> integrandPolarTurningPoint;
	std::function<void(Complex &, Complex &, int const &, int const &, int const &, int const &, GeodesicConstants &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, Complex const &, Complex const &, Complex const &,  Complex const &, Complex const &, Complex const &, double const &, double const &, double const &)> integrandRadialPolarTurningPoint;
	// std::function<Complex(int const &, int const &, int const &, int const &, GeodesicConstants &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, Complex const &, Complex const &, Complex const &, double const &, double const &, double const &)> integrandCompare;
	if(s == -2){
		integrand = teukolskyIntegrandMinus2;
		integrandRadialTurningPoint = teukolskyIntegrandMinus2RadialTurningPoint;
		integrandPolarTurningPoint = teukolskyIntegrandMinus2PolarTurningPoint;
		integrandRadialPolarTurningPoint = teukolskyIntegrandMinus2RadialPolarTurningPoint;
		// integrandCompare = teukolskyIntegrand;
	}else{
		integrand = teukolskyIntegrandPlus2;
		integrandRadialTurningPoint = teukolskyIntegrandPlus2RadialTurningPoint;
		integrandPolarTurningPoint = teukolskyIntegrandPlus2PolarTurningPoint;
		integrandRadialPolarTurningPoint = teukolskyIntegrandPlus2RadialPolarTurningPoint;
	}
	
	const ComplexVector &R0 = (Rin.solution);
	const ComplexVector &Rp0 = (Rin.derivative);
	const ComplexVector &Rpp0 = (Rin.secondDerivative);

	const ComplexVector &R1 = (Rup.solution);
	const ComplexVector &Rp1 = (Rup.derivative);
	const ComplexVector &Rpp1 = (Rup.secondDerivative);

	const Vector &S = (Slm.solution);
	const Vector &Sp = (Slm.derivative);
	const Vector &Spp = (Slm.secondDerivative);

	const Vector &rp = traj.r;
	const Vector &thp = traj.theta;
	const Vector &tR = traj.tR;
	const Vector &phiR = traj.phiR;
	const Vector &tTh = traj.tTheta;
	const Vector &phiTh = traj.phiTheta;

	Complex W = wronskian(s, geoConstants.a, rp[0], R0[0], Rp0[0], R1[0], Rp1[0]);

	int NsampleR = 8, NsampleTh = 4;
	while(NsampleR < 2*std::abs(n) + 2){
		NsampleR *= 2;
	}
	while(NsampleTh < 2*std::abs(k) + 2){
		NsampleTh *= 2;
	}
	int radialLength = rp.size(), polarLength = thp.size();
	int sampleSizeR = radialLength - 1, sampleSizeTh = polarLength - 1;
	int NsampleMaxR = 2*sampleSizeR, NsampleMaxTh = 2*sampleSizeTh;
	if(NsampleMaxR < NsampleR){
		NsampleR = NsampleMaxR;
	}

	if(NsampleMaxTh < NsampleTh){
		NsampleTh = NsampleMaxTh;
	}

	int halfSampleR = NsampleR/2, halfSampleTh = NsampleTh/2;
	int sampleDiffR = sampleSizeR/halfSampleR, sampleDiffTh = sampleSizeTh/halfSampleTh;
	double deltaQR = M_PI/double(sampleSizeR), deltaQTh = M_PI/double(sampleSizeTh);
	const double signK = (k % 2 == 0) ? 1. : -1.;
	const double signN = (n % 2 == 0) ? 1. : -1.;
	const double signKN = signK*signN;

	// first add the points at qr = 0 and qr = pi
	Complex ZlmUp, ZlmIn;
	int samplePosR = 0, samplePosTh = 0;
	double qr = 0.;
	double qth = 0.;
	Complex sumUpTerm, sumInTerm;
	SummationHelper sumUp;
	SummationHelper sumIn;
	
	// first account for the fact that we expect some numerical error in our
	// radial and polar solutions on the order of 10^-14
	sumUp.setBasePrecision(5.e-14);
	sumIn.setBasePrecision(5.e-14);

	// initial sum over fixed qr = 0. and qr = pi; qth = 0 and qth = pi
	// std::cout << "qr = " << qr/M_PI << ", qth = " << qth/M_PI << "\n";
	integrandRadialPolarTurningPoint(sumInTerm, sumUpTerm, L, m, k, n, geoConstants, 0., 0., rp[samplePosR], thp[samplePosTh], 0., 0., qr, qth,
		R0[samplePosR], Rp0[samplePosR], Rpp0[samplePosR], R1[samplePosR], Rp1[samplePosR], Rpp1[samplePosR], S[samplePosTh], Sp[samplePosTh], Spp[samplePosTh]);
	sumUp.add(0.25*sumUpTerm);
	sumIn.add(0.25*sumInTerm);

	samplePosTh = halfSampleTh*sampleDiffTh;
	qth = double(samplePosTh)*deltaQTh;
	// std::cout << "qr = " << qr/M_PI << ", qth = " << qth/M_PI << "\n";
	integrandRadialPolarTurningPoint(sumInTerm, sumUpTerm, L, m, k, n, geoConstants, 0., 0., rp[samplePosR], thp[samplePosTh], 0., 0., qr, qth,
		R0[samplePosR], Rp0[samplePosR], Rpp0[samplePosR], R1[samplePosR], Rp1[samplePosR], Rpp1[samplePosR], S[samplePosTh], Sp[samplePosTh], Spp[samplePosTh]);
	sumUp.add(signK*0.25*sumUpTerm);
	sumIn.add(signK*0.25*sumInTerm);

	samplePosR = halfSampleR*sampleDiffR;
	qr = double(samplePosR)*deltaQR;
	// std::cout << "qr = " << qr/M_PI << ", qth = " << qth/M_PI << "\n";
	integrandRadialPolarTurningPoint(sumInTerm, sumUpTerm, L, m, k, n, geoConstants, 0., 0., rp[samplePosR], thp[samplePosTh], 0., 0., qr, qth,
		R0[samplePosR], Rp0[samplePosR], Rpp0[samplePosR], R1[samplePosR], Rp1[samplePosR], Rpp1[samplePosR], S[samplePosTh], Sp[samplePosTh], Spp[samplePosTh]);
	sumUp.add(signKN*0.25*sumUpTerm);
	sumIn.add(signKN*0.25*sumInTerm);

	samplePosTh = 0;
	qth = double(samplePosTh)*deltaQTh;
	// std::cout << "qr = " << qr/M_PI << ", qth = " << qth/M_PI << "\n";
	integrandRadialPolarTurningPoint(sumInTerm, sumUpTerm, L, m, k, n, geoConstants, 0., 0., rp[samplePosR], thp[samplePosTh], 0., 0., qr, qth,
		R0[samplePosR], Rp0[samplePosR], Rpp0[samplePosR], R1[samplePosR], Rp1[samplePosR], Rpp1[samplePosR], S[samplePosTh], Sp[samplePosTh], Spp[samplePosTh]);
	sumUp.add(signN*0.25*sumUpTerm);
	sumIn.add(signN*0.25*sumInTerm);

	// initial sum over qr values with fixed qth = 0. and qth = pi
	for(int i = 1; i < halfSampleR; i++){
		samplePosR = i*sampleDiffR;
		qr = double(samplePosR)*deltaQR;

		samplePosTh = 0.;
		qth = double(samplePosTh)*deltaQTh;
		// std::cout << "qr = " << qr/M_PI << ", qth = " << qth/M_PI << "\n";
		// first sum performs integration with qth = 0
		integrandPolarTurningPoint(sumInTerm, sumUpTerm, L, m, k, n, geoConstants, tR[samplePosR], 0., rp[samplePosR], thp[samplePosTh], phiR[samplePosR], 0., qr, qth,
			R0[samplePosR], Rp0[samplePosR], Rpp0[samplePosR], R1[samplePosR], Rp1[samplePosR], Rpp1[samplePosR], S[samplePosTh], Sp[samplePosTh], Spp[samplePosTh]);
		sumUp.add(0.5*sumUpTerm);
		sumIn.add(0.5*sumInTerm);

		samplePosTh = halfSampleTh*sampleDiffTh;
		qth = double(samplePosTh)*deltaQTh;
		// std::cout << "qr = " << qr/M_PI << ", qth = " << qth/M_PI << "\n";
		// second sum performs integration with qth = pi
		integrandPolarTurningPoint(sumInTerm, sumUpTerm, L, m, k, n, geoConstants, tR[samplePosR], 0., rp[samplePosR], thp[samplePosTh], phiR[samplePosR], 0., qr, qth,
			R0[samplePosR], Rp0[samplePosR], Rpp0[samplePosR], R1[samplePosR], Rp1[samplePosR], Rpp1[samplePosR], S[samplePosTh], Sp[samplePosTh], Spp[samplePosTh]);
		sumUp.add(signK*0.5*sumUpTerm);
		sumIn.add(signK*0.5*sumInTerm);
	}

	// initial sum over qth values with fixed qr = 0. and qr = pi
	for(int j = 1; j < halfSampleTh; j++){
		samplePosTh = j*sampleDiffTh;
		qth = double(samplePosTh)*deltaQTh;

		samplePosR = 0.;
		qr = double(samplePosR)*deltaQR;
		// std::cout << "qr = " << qr/M_PI << ", qth = " << qth/M_PI << "\n";
		// std::cout << "sample position = " << samplePos << "/" << NsampleMax << "\n";
		// first sum performs integration between qr = 0 to qr = pi
		integrandRadialTurningPoint(sumInTerm, sumUpTerm, L, m, k, n, geoConstants, 0., tTh[samplePosTh], rp[samplePosR], thp[samplePosTh], 0., phiTh[samplePosTh], qr, qth,
			R0[samplePosR], Rp0[samplePosR], Rpp0[samplePosR], R1[samplePosR], Rp1[samplePosR], Rpp1[samplePosR], S[samplePosTh], Sp[samplePosTh], Spp[samplePosTh]);
		sumUp.add(0.5*sumUpTerm);
		sumIn.add(0.5*sumInTerm);

		samplePosR = halfSampleR*sampleDiffR;
		qr = double(samplePosR)*deltaQR;
		// std::cout << "qr = " << 2. - qr/M_PI << ", qth = " << qth/M_PI << "\n";
		integrandRadialTurningPoint(sumInTerm, sumUpTerm, L, m, k, n, geoConstants, 0., tTh[samplePosTh], rp[samplePosR], thp[samplePosTh], 0., phiTh[samplePosTh], qr, qth,
			R0[samplePosR], Rp0[samplePosR], Rpp0[samplePosR], R1[samplePosR], Rp1[samplePosR], Rpp1[samplePosR], S[samplePosTh], Sp[samplePosTh], Spp[samplePosTh]);
		sumUp.add(signN*0.5*sumUpTerm);
		sumIn.add(signN*0.5*sumInTerm);
	}

	// initial sum over mixed qr and qth values
	for(int i = 1; i < halfSampleR; i++){
		samplePosR = i*sampleDiffR;
		qr = double(samplePosR)*deltaQR;

		for(int j = 1; j < halfSampleTh; j++){
			samplePosTh = j*sampleDiffTh;
			qth = double(samplePosTh)*deltaQTh;
			// std::cout << "qr = " << qr/M_PI << ", qth = " << qth/M_PI << "\n";
			// std::cout << "sample position = " << samplePos << "/" << NsampleMax << "\n";
			// first sum performs integration between qr = 0 to qr = pi
			integrand(sumInTerm, sumUpTerm, L, m, k, n, geoConstants, tR[samplePosR], tTh[samplePosTh], rp[samplePosR], thp[samplePosTh], phiR[samplePosR], phiTh[samplePosTh], qr, qth,
				R0[samplePosR], Rp0[samplePosR], Rpp0[samplePosR], R1[samplePosR], Rp1[samplePosR], Rpp1[samplePosR], S[samplePosTh], Sp[samplePosTh], Spp[samplePosTh]);
			sumUp.add(sumUpTerm);
			sumIn.add(sumInTerm);
		}
	}

	ZlmUp = sumUp.getSum()/double(halfSampleR*halfSampleTh);
	ZlmIn = sumIn.getSum()/double(halfSampleR*halfSampleTh);
	Complex ZlmUpCompare = 0., ZlmInCompare = 0.;
	// std::cout << ZlmIn << ", " << ZlmUp << " for " << NsampleR << ", " << NsampleTh << "\n";

	double errorTolerance = 5.e-11;
	int convergenceTest = 0;
	while(NsampleR < NsampleMaxR && convergenceTest < 2){
		// needs to pass the convergence test twice
		if(integrand_convergence(ZlmUpCompare, ZlmUp, errorTolerance, 10.*sumUp.getPrecision()) && integrand_convergence(ZlmInCompare, ZlmIn, errorTolerance, 10.*sumIn.getPrecision())){
			convergenceTest += 1;
		}else{
			convergenceTest = 0;
		}

		Complex ZlmUpCompareTh = 0., ZlmInCompareTh = 0.;
		Complex ZlmUpTh = ZlmUp, ZlmInTh = ZlmIn;
		int convergenceTestTh = 0;
		while(NsampleTh < NsampleMaxTh && convergenceTestTh < 2){
			// needs to pass the convergence test twice
			if(integrand_convergence(ZlmUpCompareTh, ZlmUpTh, errorTolerance, 10.*sumUp.getPrecision()) && integrand_convergence(ZlmInCompareTh, ZlmInTh, errorTolerance, 10.*sumIn.getPrecision())){
				convergenceTestTh += 1;
			}else{
				convergenceTestTh = 0;
			}
			for(int j = 0; j < halfSampleTh; j++){
				samplePosTh = j*sampleDiffTh + sampleDiffTh/2;
				qth = double(samplePosTh)*deltaQTh;

				samplePosR = 0.;
				qr = double(samplePosR)*deltaQR;
				integrandRadialTurningPoint(sumInTerm, sumUpTerm, L, m, k, n, geoConstants, 0., tTh[samplePosTh], rp[samplePosR], thp[samplePosTh], 0., phiTh[samplePosTh], qr, qth,
					R0[samplePosR], Rp0[samplePosR], Rpp0[samplePosR], R1[samplePosR], Rp1[samplePosR], Rpp1[samplePosR], S[samplePosTh], Sp[samplePosTh], Spp[samplePosTh]);
				sumUp.add(0.5*sumUpTerm);
				sumIn.add(0.5*sumInTerm);

				samplePosR = halfSampleR*sampleDiffR;
				qr = double(samplePosR)*deltaQR;
				integrandRadialTurningPoint(sumInTerm, sumUpTerm, L, m, k, n, geoConstants, 0., tTh[samplePosTh], rp[samplePosR], thp[samplePosTh], 0., phiTh[samplePosTh], qr, qth,
					R0[samplePosR], Rp0[samplePosR], Rpp0[samplePosR], R1[samplePosR], Rp1[samplePosR], Rpp1[samplePosR], S[samplePosTh], Sp[samplePosTh], Spp[samplePosTh]);
				sumUp.add(signN*0.5*sumUpTerm);
				sumIn.add(signN*0.5*sumInTerm);

				for(int i = 1; i < halfSampleR; i++){
					samplePosR = i*sampleDiffR;
					qr = double(samplePosR)*deltaQR;
					integrand(sumInTerm, sumUpTerm, L, m, k, n, geoConstants, tR[samplePosR], tTh[samplePosTh], rp[samplePosR], thp[samplePosTh], phiR[samplePosR], phiTh[samplePosTh], qr, qth,
						R0[samplePosR], Rp0[samplePosR], Rpp0[samplePosR], R1[samplePosR], Rp1[samplePosR], Rpp1[samplePosR], S[samplePosTh], Sp[samplePosTh], Spp[samplePosTh]);
					sumUp.add(sumUpTerm);
					sumIn.add(sumInTerm);
				}
			}

			NsampleTh *= 2;
			halfSampleTh *= 2;
			sampleDiffTh /= 2;
			ZlmUpCompareTh = ZlmUpTh;
			ZlmInCompareTh = ZlmInTh;
			ZlmUpTh = sumUp.getSum()/double(halfSampleR*halfSampleTh);
			ZlmInTh = sumIn.getSum()/double(halfSampleR*halfSampleTh);
			// std::cout << ZlmInTh << ", " << ZlmUpTh << " for " << NsampleR << ", " << NsampleTh << "\n";
		}
		// additional qr sample points for fixed qth = 0. and qth = pi
		for(int i = 0; i < halfSampleR; i++){
			samplePosR = i*sampleDiffR + sampleDiffR/2;
			qr = double(samplePosR)*deltaQR;

			samplePosTh = 0.;
			qth = double(samplePosTh)*deltaQTh;
			// std::cout << "qr = " << qr/M_PI << ", qth = " << qth/M_PI << "\n";
			integrandPolarTurningPoint(sumInTerm, sumUpTerm, L, m, k, n, geoConstants, tR[samplePosR], 0., rp[samplePosR], thp[samplePosTh], phiR[samplePosR], 0., qr, qth,
				R0[samplePosR], Rp0[samplePosR], Rpp0[samplePosR], R1[samplePosR], Rp1[samplePosR], Rpp1[samplePosR], S[samplePosTh], Sp[samplePosTh], Spp[samplePosTh]);
			sumUp.add(0.5*sumUpTerm);
			sumIn.add(0.5*sumInTerm);

			samplePosTh = halfSampleTh*sampleDiffTh;
			qth = double(samplePosTh)*deltaQTh;
			// std::cout << "qr = " << qr/M_PI << ", qth = " << 2. - qth/M_PI << "\n";
			integrandPolarTurningPoint(sumInTerm, sumUpTerm, L, m, k, n, geoConstants, tR[samplePosR], 0., rp[samplePosR], thp[samplePosTh], phiR[samplePosR], 0., qr, qth,
				R0[samplePosR], Rp0[samplePosR], Rpp0[samplePosR], R1[samplePosR], Rp1[samplePosR], Rpp1[samplePosR], S[samplePosTh], Sp[samplePosTh], Spp[samplePosTh]);
			sumUp.add(signK*0.5*sumUpTerm);
			sumIn.add(signK*0.5*sumInTerm);

			// additional qr sample points for interior qth points between 0 and pi
			for(int j = 1; j < halfSampleTh; j++){
				samplePosTh = j*sampleDiffTh;
				qth = double(samplePosTh)*deltaQTh;
				// std::cout << "qr = " << qr/M_PI << ", qth = " << qth/M_PI << "\n";
				integrand(sumInTerm, sumUpTerm, L, m, k, n, geoConstants, tR[samplePosR], tTh[samplePosTh], rp[samplePosR], thp[samplePosTh], phiR[samplePosR], phiTh[samplePosTh], qr, qth,
					R0[samplePosR], Rp0[samplePosR], Rpp0[samplePosR], R1[samplePosR], Rp1[samplePosR], Rpp1[samplePosR], S[samplePosTh], Sp[samplePosTh], Spp[samplePosTh]);
				sumUp.add(sumUpTerm);
				sumIn.add(sumInTerm);
			}
		}

		NsampleR *= 2;
		halfSampleR *= 2;
		sampleDiffR /= 2;
		ZlmUpCompare = ZlmUp;
		ZlmInCompare = ZlmIn;
		ZlmUp = sumUp.getSum()/double(halfSampleR*halfSampleTh);
		ZlmIn = sumIn.getSum()/double(halfSampleR*halfSampleTh);
		// std::cout << ZlmIn << ", " << ZlmUp << " for " << NsampleR << ", " << NsampleTh << "\n";
		// std::cout << "Teukolsky Up amplitude = " << ZlmUp << " with "<<NsampleR*NsampleTh<<" samples \n";
		// std::cout << "Precision of Teukolsky Up amplitude = " << std::abs(1. - ZlmUpCompare/ZlmUp) << " with "<<NsampleR*NsampleTh<<" samples \n";
	}
	
	double precisionIn = std::abs(1. - ZlmInCompare/ZlmIn);
	double precisionUp = std::abs(1. - ZlmUpCompare/ZlmUp);
	if(precisionIn > PRECISION_THRESHOLD){
		precisionIn = std::max(sumIn.getPrecision(), precisionIn);
	}
	if(precisionUp > PRECISION_THRESHOLD){
		precisionUp = std::max(sumUp.getPrecision(), precisionUp);
	}

	ZlmUp *= -8.*M_PI/W/geoConstants.upsilonT;
	ZlmIn *= -8.*M_PI/W/geoConstants.upsilonT;

	TeukolskyAmplitudes Zlm = {ZlmIn, ZlmUp, precisionIn, precisionUp};

	return Zlm;
}

void weight_solutions_delta2(Complex &delta2R0, Complex &delta2Rp0, Complex &delta2Rpp0, double a, double rp, Complex R0, Complex Rp0, Complex Rpp0){
	double delta = rp*rp - 2.*rp + a*a;
	double deltaP = 2.*(rp - 1.);
	double deltaPP = 2.;
	
	delta2R0 = delta*delta*R0;
	delta2Rp0 = delta*delta*Rp0 + 2.*delta*deltaP*R0;
	delta2Rpp0 = delta*delta*Rpp0 + 4.*delta*deltaP*Rp0 + 2.*(deltaP*deltaP + delta*deltaPP)*R0;
}

void rescale_solution(ComplexVector &R, ComplexVector &Rp, ComplexVector &Rpp, double rescale){
	for(size_t i = 0; i < R.size(); i++){
		R[i] *= rescale;
		Rp[i] *= rescale;
		Rpp[i] *= rescale;
	}
}

TeukolskyAmplitudes teukolsky_amplitude_ecceq(int s, int L, int m, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConstants, const ComplexDerivativesMatrixStruct &Rin, const ComplexDerivativesMatrixStruct &Rup, const DerivativesMatrix &Slm){
	std::function<void(Complex &, Complex &, int const &, int const &, int const &, int const &, GeodesicConstants &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, Complex const &, Complex const &, Complex const &,  Complex const &, Complex const &, Complex const &, double const &, double const &, double const &)> integrand;
	std::function<void(Complex &, Complex &, int const &, int const &, int const &, int const &, GeodesicConstants &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, Complex const &, Complex const &, Complex const &,  Complex const &, Complex const &, Complex const &, double const &, double const &, double const &)> integrandTurningPoint;
	// std::function<Complex(int const &, int const &, int const &, int const &, GeodesicConstants &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, Complex const &, Complex const &, Complex const &, double const &, double const &, double const &)> integrandCompare;
	if(s == -2){
		integrand = teukolskyIntegrandMinus2PolarTurningPoint;
		integrandTurningPoint = teukolskyIntegrandMinus2RadialPolarTurningPoint;
	}else{
		integrand = teukolskyIntegrandPlus2PolarTurningPoint;
		integrandTurningPoint = teukolskyIntegrandPlus2RadialPolarTurningPoint;
	}
	const ComplexVector *R0ptr = &Rin.solution;
	const ComplexVector *Rp0ptr = &Rin.derivative;
	const ComplexVector *Rpp0ptr = &Rin.secondDerivative;

	const ComplexVector *R1ptr = &Rup.solution;
	const ComplexVector *Rp1ptr = &Rup.derivative;
	const ComplexVector *Rpp1ptr = &Rup.secondDerivative;

	double S = (Slm.solution)[0];
	double Sp = (Slm.derivative)[0];
	double Spp = (Slm.secondDerivative)[0];

	const Vector &rp = (traj.r);
	double thp = 0.5*M_PI;
	const Vector &tR = traj.tR;
	const Vector &phiR = traj.phiR;

	double rescaleIn = 1.;
	double rescaleUp = 1.;

	// local storage only used if the solutions need to be rescaled to avoid
	// an overflowing Wronskian
	ComplexVector R0loc, Rp0loc, Rpp0loc, R1loc, Rp1loc, Rpp1loc;

	Complex W = wronskian(s, geoConstants.a, rp[0], (*R0ptr)[0], (*Rp0ptr)[0], (*R1ptr)[0], (*Rp1ptr)[0]);
	if(isnan(std::abs(W)) || isinf(std::abs(W))){
		R0loc = Rin.solution;
		Rp0loc = Rin.derivative;
		Rpp0loc = Rin.secondDerivative;
		R1loc = Rup.solution;
		Rp1loc = Rup.derivative;
		Rpp1loc = Rup.secondDerivative;
		rescaleIn = 1./std::abs(R1loc[0]);
		rescaleUp = 1./std::abs(R0loc[R0loc.size() - 1]);
		rescale_solution(R0loc, Rp0loc, Rpp0loc, rescaleIn);
		rescale_solution(R1loc, Rp1loc, Rpp1loc, rescaleUp);
		R0ptr = &R0loc; Rp0ptr = &Rp0loc; Rpp0ptr = &Rpp0loc;
		R1ptr = &R1loc; Rp1ptr = &Rp1loc; Rpp1ptr = &Rpp1loc;
		W = wronskian(s, geoConstants.a, rp[0], (*R0ptr)[0], (*Rp0ptr)[0], (*R1ptr)[0], (*Rp1ptr)[0]);
	}
	const ComplexVector &R0 = *R0ptr;
	const ComplexVector &Rp0 = *Rp0ptr;
	const ComplexVector &Rpp0 = *Rpp0ptr;
	const ComplexVector &R1 = *R1ptr;
	const ComplexVector &Rp1 = *Rp1ptr;
	const ComplexVector &Rpp1 = *Rpp1ptr;

	int Nsample = 8;
	while(Nsample < 2*std::abs(n) + 2){
		Nsample *= 2;
	}
	int radialLength = rp.size();
	int sampleSize = radialLength - 1;
	int NsampleMax = 2*sampleSize;
	if(NsampleMax < Nsample){
		Nsample = NsampleMax/2;
	}
	int halfSample = Nsample/2;
	int sampleDiff = sampleSize/halfSample;
	double deltaQ = M_PI/double(sampleSize);
	const double signN = (n % 2 == 0) ? 1. : -1.;
	SummationHelper sumUp;
	SummationHelper sumIn;
	// first account for the fact that we expect some numerical error in our
	// radial and polar solutions on the order of 10^-14
	sumUp.setBasePrecision(5.e-14);
	sumIn.setBasePrecision(5.e-14);

	// first add the points at qr = 0 and qr = pi
	// Note that we need to divide by two to avoid double counting, since we are integrating from 0 to pi
	Complex ZlmUp, ZlmIn;
	int samplePos = 0.;
	double qr = samplePos*deltaQ;
	double qth = 0.;
	// double maxTerm = 0.;
	Complex sumUpTerm, sumInTerm;
	integrandTurningPoint(sumInTerm, sumUpTerm, L, m, 0, n, geoConstants, 0., 0., rp[samplePos], thp, 0., 0., qr, qth,
		R0[samplePos], Rp0[samplePos], Rpp0[samplePos], R1[samplePos], Rp1[samplePos], Rpp1[samplePos], S, Sp, Spp);
	sumUp.add(0.5*sumUpTerm);
	sumIn.add(0.5*sumInTerm);

	samplePos = halfSample*sampleDiff;
	qr = M_PI;
	integrandTurningPoint(sumInTerm, sumUpTerm, L, m, 0, n, geoConstants, 0., 0., rp[samplePos], thp, 0., 0., qr, qth,
		R0[samplePos], Rp0[samplePos], Rpp0[samplePos], R1[samplePos], Rp1[samplePos], Rpp1[samplePos], S, Sp, Spp);
	sumUp.add(signN*0.5*sumUpTerm);
	sumIn.add(signN*0.5*sumInTerm);

	for(int i = 1; i < halfSample; i++){
		samplePos = i*sampleDiff;
		qr = double(samplePos)*deltaQ;
		integrand(sumInTerm, sumUpTerm, L, m, 0, n, geoConstants, tR[samplePos], 0., rp[samplePos], thp, phiR[samplePos], 0., qr, qth,
			R0[samplePos], Rp0[samplePos], Rpp0[samplePos], R1[samplePos], Rp1[samplePos], Rpp1[samplePos], S, Sp, Spp);
		sumUp.add(sumUpTerm);
		sumIn.add(sumInTerm);
	}

	ZlmUp = sumUp.getSum()/double(halfSample);
	ZlmIn = sumIn.getSum()/double(halfSample);

	double errorTolerance = 5.e-12;
	double convergenceTest = 0;
	Complex ZlmUpCompare = 0., ZlmInCompare = 0.;
	while(convergenceTest < 2 && Nsample < NsampleMax){
		if(integrand_convergence(ZlmUpCompare, ZlmUp, errorTolerance, 10.*sumUp.getPrecision()) && integrand_convergence(ZlmInCompare, ZlmIn, errorTolerance, 10.*sumIn.getPrecision())){
			convergenceTest += 1;
		}else{
			convergenceTest = 0;
		}
		for(int i = 0; i < halfSample; i++){
			samplePos = i*sampleDiff + sampleDiff/2;
			qr = double(samplePos)*deltaQ;
			integrand(sumInTerm, sumUpTerm, L, m, 0, n, geoConstants, tR[samplePos], 0., rp[samplePos], thp, phiR[samplePos], 0., qr, qth,
				R0[samplePos], Rp0[samplePos], Rpp0[samplePos], R1[samplePos], Rp1[samplePos], Rpp1[samplePos], S, Sp, Spp);
			sumUp.add(sumUpTerm);
			sumIn.add(sumInTerm);
		}
		Nsample *= 2;
		halfSample *= 2;
		sampleDiff /= 2;
		ZlmUpCompare = ZlmUp;
		ZlmInCompare = ZlmIn;
		ZlmUp = sumUp.getSum()/double(halfSample);
		ZlmIn = sumIn.getSum()/double(halfSample);
	}
	// std::cout << Nsample << "\n";
	// std::cout << sumIn.getPrecision() << "\n";
	// std::cout << sumUp.getPrecision() << "\n";
	// std::cout << ZlmIn << "\n";
	// std::cout << ZlmUp << "\n";
	// if(sumIn.getPrecision()){

	// }
	double precisionIn = std::abs(1. - ZlmInCompare/ZlmIn);
	double precisionUp = std::abs(1. - ZlmUpCompare/ZlmUp);

	if(precisionIn > PRECISION_THRESHOLD){
		precisionIn = std::max(sumIn.getPrecision(), precisionIn);
	}
	if(precisionUp > PRECISION_THRESHOLD){
		precisionUp = std::max(sumUp.getPrecision(), precisionUp);
	}

	ZlmUp *= -rescaleUp*8.*M_PI/W/geoConstants.upsilonT;
	ZlmIn *= -rescaleIn*8.*M_PI/W/geoConstants.upsilonT;

	TeukolskyAmplitudes Zlm = {ZlmIn, ZlmUp, precisionIn, precisionUp};

	return Zlm;
}

TeukolskyAmplitudes teukolsky_amplitude_sphinc(int s, int L, int m, int k, GeodesicTrajectory& traj, GeodesicConstants &geoConstants, const ComplexDerivativesMatrixStruct &Rin, const ComplexDerivativesMatrixStruct &Rup, const DerivativesMatrix &Slm){
	std::function<void(Complex &, Complex &, int const &, int const &, int const &, int const &, GeodesicConstants &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, Complex const &, Complex const &, Complex const &,  Complex const &, Complex const &, Complex const &, double const &, double const &, double const &)> integrand;
	std::function<void(Complex &, Complex &, int const &, int const &, int const &, int const &, GeodesicConstants &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, Complex const &, Complex const &, Complex const &,  Complex const &, Complex const &, Complex const &, double const &, double const &, double const &)> integrandTurningPoint;
	if(s == -2){
		integrand = teukolskyIntegrandMinus2RadialTurningPoint;
		integrandTurningPoint = teukolskyIntegrandMinus2RadialPolarTurningPoint;
	}else{
		integrand = teukolskyIntegrandPlus2RadialTurningPoint;
		integrandTurningPoint = teukolskyIntegrandPlus2RadialPolarTurningPoint;
	}
	
	Complex R0 = (Rin.solution)[0];
	Complex Rp0 = (Rin.derivative)[0];
	Complex Rpp0 = (Rin.secondDerivative)[0];

	Complex R1 = (Rup.solution)[0];
	Complex Rp1 = (Rup.derivative)[0];
	Complex Rpp1 = (Rup.secondDerivative)[0];

	const Vector &S = (Slm.solution);
	const Vector &Sp = (Slm.derivative);
	const Vector &Spp = (Slm.secondDerivative);

	double rp = (traj.r)[0];
	const Vector &thp = (traj.theta);
	const Vector &tTh = traj.tTheta;
	const Vector &phiTh = traj.phiTheta;

	Complex W = wronskian(s, geoConstants.a, rp, R0, Rp0, R1, Rp1);

	int Nsample = 8;
	while(Nsample < 2*std::abs(k) + 2){
		Nsample *= 2;
	}
	int polarLength = thp.size();
	int sampleSize = polarLength - 1;
	int NsampleMax = 2*sampleSize;
	if(NsampleMax < Nsample){
		Nsample = NsampleMax;
	}
	int halfSample = Nsample/2;
	int sampleDiff = sampleSize/halfSample;
	double deltaQ = 2.*M_PI/double(NsampleMax);
	const double signK = (k % 2 == 0) ? 1. : -1.;
	SummationHelper sumUp;
	SummationHelper sumIn;
	// first account for the fact that we expect some numerical error in our
	// radial and polar solutions on the order of 10^-14
	sumUp.setBasePrecision(5.e-14);
	sumIn.setBasePrecision(5.e-14);

	// first add the points at qr = 0 and qr = pi
	Complex ZlmUp, ZlmIn;
	int samplePos = 0.;
	double qth = 0.;
	double qr = 0.;
	Complex sumUpTerm, sumInTerm;
	integrandTurningPoint(sumInTerm, sumUpTerm, L, m, k, 0, geoConstants, 0., 0., rp, thp[samplePos], 0., 0., qr, qth,
		R0, Rp0, Rpp0, R1, Rp1, Rpp1, S[samplePos], Sp[samplePos], Spp[samplePos]);
	sumUp.add(0.5*sumUpTerm);
	sumIn.add(0.5*sumInTerm);

	samplePos = halfSample*sampleDiff;
	qth = M_PI;
	integrandTurningPoint(sumInTerm, sumUpTerm, L, m, k, 0, geoConstants, 0., 0., rp, thp[samplePos], 0., 0., qr, qth,
		R0, Rp0, Rpp0, R1, Rp1, Rpp1, S[samplePos], Sp[samplePos], Spp[samplePos]);
	sumUp.add(signK*0.5*sumUpTerm);
	sumIn.add(signK*0.5*sumInTerm);

	for(int i = 1; i < halfSample; i++){
		samplePos = i*sampleDiff;
		qth = samplePos*deltaQ;
		// first sum performs integration between qr = 0 to qr = pi
		integrand(sumInTerm, sumUpTerm, L, m, k, 0, geoConstants, 0., tTh[samplePos], rp, thp[samplePos], 0., phiTh[samplePos], qr, qth,
			R0, Rp0, Rpp0, R1, Rp1, Rpp1, S[samplePos], Sp[samplePos], Spp[samplePos]);
		sumUp.add(sumUpTerm);
		sumIn.add(sumInTerm);
	}

	ZlmUp = sumUp.getSum()/double(halfSample);
	ZlmIn = sumIn.getSum()/double(halfSample);

	double errorTolerance = 5.e-12;
	double convergenceTest = 0;
	Complex ZlmUpCompare = 0., ZlmInCompare = 0.;
	while(convergenceTest < 2 && Nsample < NsampleMax){
		if(integrand_convergence(ZlmUpCompare, ZlmUp, errorTolerance, sumUp.getPrecision()) && integrand_convergence(ZlmInCompare, ZlmIn, errorTolerance, sumIn.getPrecision())){
			convergenceTest += 1;
		}else{
			convergenceTest = 0;
		}
		for(int i = 0; i < halfSample; i++){
			samplePos = i*sampleDiff + sampleDiff/2;
			qth = samplePos*deltaQ;
			integrand(sumInTerm, sumUpTerm, L, m, k, 0, geoConstants, 0., tTh[samplePos], rp, thp[samplePos], 0., phiTh[samplePos], qr, qth,
				R0, Rp0, Rpp0, R1, Rp1, Rpp1, S[samplePos], Sp[samplePos], Spp[samplePos]);
			sumUp.add(sumUpTerm);
			sumIn.add(sumInTerm);
		}
		Nsample *= 2;
		halfSample *= 2;
		sampleDiff /= 2;
		ZlmUpCompare = ZlmUp;
		ZlmInCompare = ZlmIn;
		ZlmUp = sumUp.getSum()/double(halfSample);
		ZlmIn = sumIn.getSum()/double(halfSample);
	}

	double precisionIn = std::abs(1. - ZlmInCompare/ZlmIn);
	double precisionUp = std::abs(1. - ZlmUpCompare/ZlmUp);

	if(precisionIn > PRECISION_THRESHOLD){
		precisionIn = std::max(sumIn.getPrecision(), precisionIn);
	}
	if(precisionUp > PRECISION_THRESHOLD){
		precisionUp = std::max(sumUp.getPrecision(), precisionUp);
	}

	ZlmUp *= -8.*M_PI/W/geoConstants.upsilonT;
	ZlmIn *= -8.*M_PI/W/geoConstants.upsilonT;

	TeukolskyAmplitudes Zlm = {ZlmIn, ZlmUp, precisionIn, precisionUp};

	return Zlm;
}

TeukolskyAmplitudes teukolsky_amplitude_circeq(int s, int L, int m, GeodesicTrajectory& traj, GeodesicConstants &geoConstants, const ComplexDerivativesMatrixStruct &Rin, const ComplexDerivativesMatrixStruct &Rup, const DerivativesMatrix &Slm){
	std::function<void(Complex &, Complex &, int const &, int const &, int const &, int const &, GeodesicConstants &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, double const &, Complex const &, Complex const &, Complex const &,  Complex const &, Complex const &, Complex const &, double const &, double const &, double const &)> integrand;
	if(s == -2){
		integrand = teukolskyIntegrandMinus2RadialPolarTurningPoint;
	}else{
		integrand = teukolskyIntegrandPlus2RadialPolarTurningPoint;
	}
	
	Complex R0 = (Rin.solution)[0];
	Complex Rp0 = (Rin.derivative)[0];
	Complex Rpp0 = (Rin.secondDerivative)[0];

	Complex R1 = (Rup.solution)[0];
	Complex Rp1 = (Rup.derivative)[0];
	Complex Rpp1 = (Rup.secondDerivative)[0];

	double S = (Slm.solution)[0];
	double Sp = (Slm.derivative)[0];
	double Spp = (Slm.secondDerivative)[0];

	double rp = (traj.r)[0];
	double thp = (traj.theta)[0];

	Complex W = wronskian(s, geoConstants.a, rp, R0, Rp0, R1, Rp1);

	Complex ZlmUp, ZlmIn;
	integrand(ZlmIn, ZlmUp, L, m, 0, 0, geoConstants, 0., 0., rp, thp, 0., 0., 0., 0.,
		R0, Rp0, Rpp0, R1, Rp1, Rpp1, S, Sp, Spp);

	ZlmUp *= -8.*M_PI/W/geoConstants.upsilonT;
	ZlmIn *= -8.*M_PI/W/geoConstants.upsilonT;

	TeukolskyAmplitudes Zlm = {ZlmIn, ZlmUp, DBL_EPSILON, DBL_EPSILON};

	return Zlm;
}

void teukolskyIntegrandMinus2(Complex &integrandIn, Complex &integrandUp, int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &tR, double const &tTh, double const &rp, double const &thp, double const &phiR, double const &phiTh, double const &qr, double const &qth, Complex const &Rin, Complex const &RinP, Complex const &RinPP,  Complex const &Rup, Complex const &RupP, Complex const &RupPP, double const &St, double const &StP, double const &StPP){	
	double freq = (m*geoConstants.upsilonPhi + k*geoConstants.upsilonTheta + n*geoConstants.upsilonR)/geoConstants.upsilonT;
	double rphase = (n*qr + freq*tR - m*phiR);
	double thphase = (k*qth + freq*tTh - m*phiTh);
	double a = geoConstants.a;

	Complex u2p, u2m, u4p, u4m;
	u_24_coeffs(u2m, u2p, u4m, u4p, geoConstants, rp, thp);
	// the four phase factors e^{i(±rphase ± thphase)} are products and
	// conjugates of just two complex exponentials
	Complex er = exp(I*rphase), eth = exp(I*thphase);
	Complex epp = er*eth, epm = er*conj(eth);
	Complex emp = conj(epm), emm = conj(epp);
	Complex Cnn, Cnmbar, Cmbarmbar;
	Cnn = u2p*u2p*(epp + epm);
	Cnn += u2m*u2m*(emp + emm);
	Cnmbar = u2p*u4p*epp;
	Cnmbar += u2p*u4m*epm;
	Cnmbar += u2m*u4p*emp;
	Cnmbar += u2m*u4m*emm;
	Cmbarmbar = u4p*u4p*(epp + emp);
	Cmbarmbar += u4m*u4m*(epm + emm);

	Complex Ann0, Anmbar0, Anmbar1, Ambarmbar0, Ambarmbar1, Ambarmbar2;
	A_coeffs_w(Ann0, Anmbar0, Ambarmbar0, Anmbar1, Ambarmbar1, Ambarmbar2, m, freq, geoConstants, rp, thp, St, StP, StPP);

	Complex prefactorR = Ann0*Cnn + Anmbar0*Cnmbar + Ambarmbar0*Cmbarmbar;
	Complex prefactorRp = Anmbar1*Cnmbar + Ambarmbar1*Cmbarmbar;
	Complex prefactorRpp = Ambarmbar2*Cmbarmbar;

	integrandUp = 0.25*( prefactorR*Rin - prefactorRp*RinP + prefactorRpp*RinPP );
	integrandIn = 0.25*( prefactorR*Rup - prefactorRp*RupP + prefactorRpp*RupPP );
}

void teukolskyIntegrandMinus2PolarTurningPoint(Complex &integrandIn, Complex &integrandUp, int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &tR, double const &tTh, double const &rp, double const &thp, double const &phiR, double const &phiTh, double const &qr, double const &qth, Complex const &Rin, Complex const &RinP, Complex const &RinPP,  Complex const &Rup, Complex const &RupP, Complex const &RupPP, double const &St, double const &StP, double const &StPP){
	double freq = (m*geoConstants.upsilonPhi + k*geoConstants.upsilonTheta + n*geoConstants.upsilonR)/geoConstants.upsilonT;
	double rphase = (n*qr + freq*tR - m*phiR);

	Complex u2p, u2m, u4;
	u_24_coeffs_PolarTurningPoint(u2m, u2p, u4, geoConstants, rp, thp);
	Complex er = exp(I*rphase), erc = conj(er);
	Complex Cnn, Cnmbar, Cmbarmbar;
	Cnn = u2p*u2p*er + u2m*u2m*erc;
	Cnmbar = u2p*u4*er + u2m*u4*erc;
	Cmbarmbar = u4*u4*(er + erc);

	Complex Ann0, Anmbar0, Anmbar1, Ambarmbar0, Ambarmbar1, Ambarmbar2;
	A_coeffs_w(Ann0, Anmbar0, Ambarmbar0, Anmbar1, Ambarmbar1, Ambarmbar2, m, freq, geoConstants, rp, thp, St, StP, StPP);

	Complex prefactorR = Ann0*Cnn + Anmbar0*Cnmbar + Ambarmbar0*Cmbarmbar;
	Complex prefactorRp = Anmbar1*Cnmbar + Ambarmbar1*Cmbarmbar;
	Complex prefactorRpp = Ambarmbar2*Cmbarmbar;

	integrandIn = 0.5*( prefactorR*Rup - prefactorRp*RupP + prefactorRpp*RupPP );
	integrandUp = 0.5*( prefactorR*Rin - prefactorRp*RinP + prefactorRpp*RinPP );
}

void teukolskyIntegrandMinus2RadialPolarTurningPoint(Complex &integrandIn, Complex &integrandUp, int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &tR, double const &tTh, double const &rp, double const &thp, double const &phiR, double const &phiTh, double const &qr, double const &qth, Complex const &Rin, Complex const &RinP, Complex const &RinPP,  Complex const &Rup, Complex const &RupP, Complex const &RupPP, double const &St, double const &StP, double const &StPP){
	double freq = (m*geoConstants.upsilonPhi + k*geoConstants.upsilonTheta + n*geoConstants.upsilonR)/geoConstants.upsilonT;
	Complex u2, u4;
	u_24_coeffs_RadialPolarTurningPoint(u2, u4, geoConstants, rp, thp);
	Complex Cnn, Cnmbar, Cmbarmbar;
	Cnn = u2*u2;
	Cnmbar = u2*u4;
	Cmbarmbar = u4*u4;

	Complex Ann0, Anmbar0, Anmbar1, Ambarmbar0, Ambarmbar1, Ambarmbar2;
	A_coeffs_w(Ann0, Anmbar0, Ambarmbar0, Anmbar1, Ambarmbar1, Ambarmbar2, m, freq, geoConstants, rp, thp, St, StP, StPP);

	Complex prefactorR = Ann0*Cnn + Anmbar0*Cnmbar + Ambarmbar0*Cmbarmbar;
	Complex prefactorRp = Anmbar1*Cnmbar + Ambarmbar1*Cmbarmbar;
	Complex prefactorRpp = Ambarmbar2*Cmbarmbar;

	integrandIn = ( prefactorR*Rup - prefactorRp*RupP + prefactorRpp*RupPP );
	integrandUp = ( prefactorR*Rin - prefactorRp*RinP + prefactorRpp*RinPP );
}

void teukolskyIntegrandMinus2RadialTurningPoint(Complex &integrandIn, Complex &integrandUp, int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &tR, double const &tTh, double const &rp, double const &thp, double const &phiR, double const &phiTh, double const &qr, double const &qth, Complex const &Rin, Complex const &RinP, Complex const &RinPP,  Complex const &Rup, Complex const &RupP, Complex const &RupPP, double const &St, double const &StP, double const &StPP){
	double freq = (m*geoConstants.upsilonPhi + k*geoConstants.upsilonTheta + n*geoConstants.upsilonR)/geoConstants.upsilonT;
	double thphase = (k*qth + freq*tTh - m*phiTh);

	Complex u2, u4p, u4m;
	u_24_coeffs_RadialTurningPoint(u2, u4m, u4p, geoConstants, rp, thp);
	Complex eth = exp(I*thphase), ethc = conj(eth);
	Complex Cnn, Cnmbar, Cmbarmbar;
	Cnn = u2*u2*(eth + ethc);
	Cnmbar = u2*u4p*eth;
	Cnmbar += u2*u4m*ethc;
	Cmbarmbar = u4p*u4p*eth;
	Cmbarmbar += u4m*u4m*ethc;

	Complex Ann0, Anmbar0, Anmbar1, Ambarmbar0, Ambarmbar1, Ambarmbar2;
	A_coeffs_w(Ann0, Anmbar0, Ambarmbar0, Anmbar1, Ambarmbar1, Ambarmbar2, m, freq, geoConstants, rp, thp, St, StP, StPP);

	Complex prefactorR = Ann0*Cnn + Anmbar0*Cnmbar + Ambarmbar0*Cmbarmbar;
	Complex prefactorRp = Anmbar1*Cnmbar + Ambarmbar1*Cmbarmbar;
	Complex prefactorRpp = Ambarmbar2*Cmbarmbar;

	integrandIn = 0.5*( prefactorR*Rup - prefactorRp*RupP + prefactorRpp*RupPP );
	integrandUp = 0.5*( prefactorR*Rin - prefactorRp*RinP + prefactorRpp*RinPP );
}

void teukolskyIntegrandPlus2(Complex &integrandIn, Complex &integrandUp, int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &tR, double const &tTh, double const &rp, double const &thp, double const &phiR, double const &phiTh, double const &qr, double const &qth, Complex const &Rin, Complex const &RinP, Complex const &RinPP,  Complex const &Rup, Complex const &RupP, Complex const &RupPP, double const &St, double const &StP, double const &StPP){	
	double freq = (m*geoConstants.upsilonPhi + k*geoConstants.upsilonTheta + n*geoConstants.upsilonR)/geoConstants.upsilonT;
	double rphase = (n*qr + freq*tR - m*phiR);
	double thphase = (k*qth + freq*tTh - m*phiTh);
	double a = geoConstants.a;
	double delta = rp*rp - 2.*rp + a*a;
	double deltaP = 2.*(rp - 1.);
	double deltaPP = 2.;

	Complex delta2Rin = delta*delta*Rin;
	Complex delta2RinP = delta*delta*RinP + 2.*delta*deltaP*Rin;
	Complex delta2RinPP = delta*delta*RinPP + 4.*delta*deltaP*RinP + 2.*(deltaP*deltaP + delta*deltaPP)*Rin;

	Complex delta2Rup = delta*delta*Rup;
	Complex delta2RupP = delta*delta*RupP + 2.*delta*deltaP*Rup;
	Complex delta2RupPP = delta*delta*RupPP + 4.*delta*deltaP*RupP + 2.*(deltaP*deltaP + delta*deltaPP)*Rup;

	Complex u1p, u1m, u3p, u3m;
	u_13_coeffs(u1m, u1p, u3m, u3p, geoConstants, rp, thp);
	// the four phase factors e^{i(±rphase ± thphase)} are products and
	// conjugates of just two complex exponentials
	Complex er = exp(I*rphase), eth = exp(I*thphase);
	Complex epp = er*eth, epm = er*conj(eth);
	Complex emp = conj(epm), emm = conj(epp);
	Complex Cll, Clm, Cmm;
	Cll = u1p*u1p*(epp + epm);
	Cll += u1m*u1m*(emp + emm);
	Clm = u1p*u3p*epp;
	Clm += u1p*u3m*epm;
	Clm += u1m*u3p*emp;
	Clm += u1m*u3m*emm;
	Cmm = u3p*u3p*(epp + emp);
	Cmm += u3m*u3m*(epm + emm);

	Complex All0, Alm0, Alm1, Amm0, Amm1, Amm2;
	A13_coeffs_w(All0, Alm0, Amm0, Alm1, Amm1, Amm2, m, freq, geoConstants, rp, thp, St, StP, StPP);

	Complex prefactorR = All0*Cll + Alm0*Clm + Amm0*Cmm;
	Complex prefactorRp = Alm1*Clm + Amm1*Cmm;
	Complex prefactorRpp = Amm2*Cmm;

	integrandUp = 0.25*( prefactorR*delta2Rin - prefactorRp*delta2RinP + prefactorRpp*delta2RinPP );
	integrandIn = 0.25*( prefactorR*delta2Rup - prefactorRp*delta2RupP + prefactorRpp*delta2RupPP );
}

void teukolskyIntegrandPlus2PolarTurningPoint(Complex &integrandIn, Complex &integrandUp, int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &tR, double const &tTh, double const &rp, double const &thp, double const &phiR, double const &phiTh, double const &qr, double const &qth, Complex const &Rin, Complex const &RinP, Complex const &RinPP,  Complex const &Rup, Complex const &RupP, Complex const &RupPP, double const &St, double const &StP, double const &StPP){	
	double freq = (m*geoConstants.upsilonPhi + k*geoConstants.upsilonTheta + n*geoConstants.upsilonR)/geoConstants.upsilonT;
	double rphase = (n*qr + freq*tR - m*phiR);

	double a = geoConstants.a;
	double delta = rp*rp - 2.*rp + a*a;
	double deltaP = 2.*(rp - 1.);
	double deltaPP = 2.;

	Complex delta2Rin = delta*delta*Rin;
	Complex delta2RinP = delta*delta*RinP + 2.*delta*deltaP*Rin;
	Complex delta2RinPP = delta*delta*RinPP + 4.*delta*deltaP*RinP + 2.*(deltaP*deltaP + delta*deltaPP)*Rin;

	Complex delta2Rup = delta*delta*Rup;
	Complex delta2RupP = delta*delta*RupP + 2.*delta*deltaP*Rup;
	Complex delta2RupPP = delta*delta*RupPP + 4.*delta*deltaP*RupP + 2.*(deltaP*deltaP + delta*deltaPP)*Rup;

	Complex u1p, u1m, u3;
	u_13_coeffs_PolarTurningPoint(u1m, u1p, u3, geoConstants, rp, thp);
	Complex er = exp(I*rphase), erc = conj(er);
	Complex Cll, Clm, Cmm;
	Cll = u1p*u1p*er;
	Cll += u1m*u1m*erc;
	Clm = u1p*u3*er;
	Clm += u1m*u3*erc;
	Cmm = u3*u3*(er + erc);

	Complex All0, Alm0, Alm1, Amm0, Amm1, Amm2;
	A13_coeffs_w(All0, Alm0, Amm0, Alm1, Amm1, Amm2, m, freq, geoConstants, rp, thp, St, StP, StPP);

	Complex prefactorR = All0*Cll + Alm0*Clm + Amm0*Cmm;
	Complex prefactorRp = Alm1*Clm + Amm1*Cmm;
	Complex prefactorRpp = Amm2*Cmm;

	integrandUp = 0.5*( prefactorR*delta2Rin - prefactorRp*delta2RinP + prefactorRpp*delta2RinPP );
	integrandIn = 0.5*( prefactorR*delta2Rup - prefactorRp*delta2RupP + prefactorRpp*delta2RupPP );
}

void teukolskyIntegrandPlus2RadialTurningPoint(Complex &integrandIn, Complex &integrandUp, int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &tR, double const &tTh, double const &rp, double const &thp, double const &phiR, double const &phiTh, double const &qr, double const &qth, Complex const &Rin, Complex const &RinP, Complex const &RinPP,  Complex const &Rup, Complex const &RupP, Complex const &RupPP, double const &St, double const &StP, double const &StPP){	
	double freq = (m*geoConstants.upsilonPhi + k*geoConstants.upsilonTheta + n*geoConstants.upsilonR)/geoConstants.upsilonT;
	double thphase = (k*qth + freq*tTh - m*phiTh);
	double a = geoConstants.a;
	double delta = rp*rp - 2.*rp + a*a;
	double deltaP = 2.*(rp - 1.);
	double deltaPP = 2.;

	Complex delta2Rin = delta*delta*Rin;
	Complex delta2RinP = delta*delta*RinP + 2.*delta*deltaP*Rin;
	Complex delta2RinPP = delta*delta*RinPP + 4.*delta*deltaP*RinP + 2.*(deltaP*deltaP + delta*deltaPP)*Rin;

	Complex delta2Rup = delta*delta*Rup;
	Complex delta2RupP = delta*delta*RupP + 2.*delta*deltaP*Rup;
	Complex delta2RupPP = delta*delta*RupPP + 4.*delta*deltaP*RupP + 2.*(deltaP*deltaP + delta*deltaPP)*Rup;

	Complex u1, u3p, u3m;
	u_13_coeffs_RadialTurningPoint(u1, u3m, u3p, geoConstants, rp, thp);
	Complex eth = exp(I*thphase), ethc = conj(eth);
	Complex Cll, Clm, Cmm;
	Cll = u1*u1*(eth + ethc);
	Clm = u1*u3p*eth;
	Clm += u1*u3m*ethc;
	Cmm = u3p*u3p*eth;
	Cmm += u3m*u3m*ethc;

	Complex All0, Alm0, Alm1, Amm0, Amm1, Amm2;
	A13_coeffs_w(All0, Alm0, Amm0, Alm1, Amm1, Amm2, m, freq, geoConstants, rp, thp, St, StP, StPP);

	Complex prefactorR = All0*Cll + Alm0*Clm + Amm0*Cmm;
	Complex prefactorRp = Alm1*Clm + Amm1*Cmm;
	Complex prefactorRpp = Amm2*Cmm;

	integrandUp = 0.5*( prefactorR*delta2Rin - prefactorRp*delta2RinP + prefactorRpp*delta2RinPP );
	integrandIn = 0.5*( prefactorR*delta2Rup - prefactorRp*delta2RupP + prefactorRpp*delta2RupPP );
}

void teukolskyIntegrandPlus2RadialPolarTurningPoint(Complex &integrandIn, Complex &integrandUp, int const &L, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &tR, double const &tTh, double const &rp, double const &thp, double const &phiR, double const &phiTh, double const &qr, double const &qth, Complex const &Rin, Complex const &RinP, Complex const &RinPP,  Complex const &Rup, Complex const &RupP, Complex const &RupPP, double const &St, double const &StP, double const &StPP){
	double freq = (m*geoConstants.upsilonPhi + k*geoConstants.upsilonTheta + n*geoConstants.upsilonR)/geoConstants.upsilonT;
	double a = geoConstants.a;
	double delta = rp*rp - 2.*rp + a*a;
	double deltaP = 2.*(rp - 1.);
	double deltaPP = 2.;

	Complex delta2Rin = delta*delta*Rin;
	Complex delta2RinP = delta*delta*RinP + 2.*delta*deltaP*Rin;
	Complex delta2RinPP = delta*delta*RinPP + 4.*delta*deltaP*RinP + 2.*(deltaP*deltaP + delta*deltaPP)*Rin;

	Complex delta2Rup = delta*delta*Rup;
	Complex delta2RupP = delta*delta*RupP + 2.*delta*deltaP*Rup;
	Complex delta2RupPP = delta*delta*RupPP + 4.*delta*deltaP*RupP + 2.*(deltaP*deltaP + delta*deltaPP)*Rup;

	Complex u1, u3;
	u_13_coeffs_RadialPolarTurningPoint(u1, u3, geoConstants, rp, thp);
	Complex Cll, Clm, Cmm;
	Cll = u1*u1;
	Clm = u1*u3;
	Cmm = u3*u3;

	Complex All0, Alm0, Alm1, Amm0, Amm1, Amm2;
	A13_coeffs_w(All0, Alm0, Amm0, Alm1, Amm1, Amm2, m, freq, geoConstants, rp, thp, St, StP, StPP);

	Complex prefactorR = All0*Cll + Alm0*Clm + Amm0*Cmm;
	Complex prefactorRp = Alm1*Clm + Amm1*Cmm;
	Complex prefactorRpp = Amm2*Cmm;

	integrandUp = ( prefactorR*delta2Rin - prefactorRp*delta2RinP + prefactorRpp*delta2RinPP );
	integrandIn = ( prefactorR*delta2Rup - prefactorRp*delta2RupP + prefactorRpp*delta2RupPP );
}

Complex A_nn_0(int const &, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &rp, double const &thp, double const &Slm, double const &SlmP, double const &SlmPP){
	double a = geoConstants.a;
	double omega = (m*geoConstants.upsilonPhi + k*geoConstants.upsilonTheta + n*geoConstants.upsilonR)/geoConstants.upsilonT;
	double varpi = rp*rp + a*a;
	double delta = varpi - 2.*rp;

	Complex rho = -1./(rp - I*a*cos(thp));
	Complex rhoBar = -1./(rp + I*a*cos(thp));

	Complex prefactor = -0.5*pow(rho, -2)*pow(rhoBar, -1)*pow(delta, -2);
	Complex L1 = -m/sin(thp) + a*omega*sin(thp) + cos(thp)/sin(thp);
	Complex L2 = -m/sin(thp) + a*omega*sin(thp) + 2.*cos(thp)/sin(thp);
	Complex L2p = m*cos(thp)/sin(thp)/sin(thp) + a*omega*cos(thp) - 2./sin(thp)/sin(thp);
	Complex L1Sp = SlmPP + L1*SlmP;
	Complex L1L2S = L1Sp + L2p*Slm + L2*SlmP + L1*L2*Slm;

	return prefactor*( pow(rho, -1)*L1L2S + 3.*I*a*sin(thp)*L1*Slm + 3.*I*a*cos(thp)*Slm + 2.*I*a*sin(thp)*SlmP - I*a*sin(thp)*L2*Slm );
}

Complex A_nmbar_0(int const &, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &rp, double const &thp, double const &Slm, double const &SlmP, double const &){
	double a = geoConstants.a;
	double omega = (m*geoConstants.upsilonPhi + k*geoConstants.upsilonTheta + n*geoConstants.upsilonR)/geoConstants.upsilonT;
	double varpi = rp*rp + a*a;
	double delta = varpi - 2.*rp;
	double K = varpi*omega - m*a;

	Complex rho = -1./(rp - I*a*cos(thp));
	Complex rhoBar = -1./(rp + I*a*cos(thp));

	Complex prefactor = pow(rho, -3)*pow(sqrt(2.)*delta, -1);
	Complex L2 = -m/sin(thp) + a*omega*sin(thp) + 2.*cos(thp)/sin(thp);
	Complex L2S = SlmP + L2*Slm;

	return prefactor*( (rho + rhoBar - I*K/delta)*L2S + (rho - rhoBar)*a*sin(thp)*K/delta*Slm );
}
Complex A_nmbar_1(int const &, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &rp, double const &thp, double const &Slm, double const &SlmP, double const &){
	double a = geoConstants.a;
	double omega = (m*geoConstants.upsilonPhi + k*geoConstants.upsilonTheta + n*geoConstants.upsilonR)/geoConstants.upsilonT;
	double varpi = rp*rp + a*a;
	double delta = varpi - 2.*rp;

	Complex rho = -1./(rp - I*a*cos(thp));
	Complex rhoBar = -1./(rp + I*a*cos(thp));

	Complex prefactor = -pow(rho, -3)*pow(sqrt(2.)*delta, -1);
	Complex L2 = -m/sin(thp) + a*omega*sin(thp) + 2.*cos(thp)/sin(thp);
	Complex L2S = SlmP + L2*Slm;

	return prefactor*( L2S + I*(rho - rhoBar)*a*sin(thp)*Slm );
}
Complex A_mbarmbar_0(int const &, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &rp, double const &thp, double const &Slm, double const &, double const &){
	double a = geoConstants.a;
	double omega = (m*geoConstants.upsilonPhi + k*geoConstants.upsilonTheta + n*geoConstants.upsilonR)/geoConstants.upsilonT;
	double varpi = rp*rp + a*a;
	double delta = varpi - 2.*rp;
	double K = varpi*omega - m*a;

	Complex rho = -1./(rp - I*a*cos(thp));
	Complex rhoBar = -1./(rp + I*a*cos(thp));

	Complex prefactor = pow(rho, -3)*rhoBar*Slm/4.;
	double deltaP = 2.*(rp - 1.);
	double Kp = 2.*omega*rp;

	return prefactor*( I*(Kp/delta - deltaP/delta*K/delta) + K/delta*K/delta + 2.*I*rho*K/delta );
}
Complex A_mbarmbar_1(int const &, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &rp, double const &thp, double const &Slm, double const &, double const &){
	double a = geoConstants.a;
	double omega = (m*geoConstants.upsilonPhi + k*geoConstants.upsilonTheta + n*geoConstants.upsilonR)/geoConstants.upsilonT;
	double varpi = rp*rp + a*a;
	double delta = varpi - 2.*rp;
	double K = varpi*omega - m*a;

	Complex rho = -1./(rp - I*a*cos(thp));
	Complex rhoBar = -1./(rp + I*a*cos(thp));

	Complex prefactor = -pow(rho, -3)*rhoBar*Slm/2.;

	return prefactor*( I*K/delta - rho );
}
Complex A_mbarmbar_2(int const &, int const &, int const &, int const &, GeodesicConstants &geoConstants, double const &rp, double const &thp, double const &Slm, double const &, double const &){
	double a = geoConstants.a;

	Complex rho = -1./(rp - I*a*cos(thp));
	Complex rhoBar = -1./(rp + I*a*cos(thp));

	Complex prefactor = -pow(rho, -3)*rhoBar*Slm/4.;

	return prefactor;
}

void A_coeffs(Complex &Ann0, Complex &Anmbar0, Complex &Ambarmbar0, Complex &Anmbar1, Complex &Ambarmbar1, Complex &Ambarmbar2, int const &, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &rp, double const &thp, double const &Slm, double const &SlmP, double const &SlmPP){
	double omega = (m*geoConstants.upsilonPhi + k*geoConstants.upsilonTheta + n*geoConstants.upsilonR)/geoConstants.upsilonT;
	A_coeffs_w(Ann0, Anmbar0, Ambarmbar0, Anmbar1, Ambarmbar1, Ambarmbar2, m, omega, geoConstants, rp, thp, Slm, SlmP, SlmPP);
}

static void A_coeffs_w(Complex &Ann0, Complex &Anmbar0, Complex &Ambarmbar0, Complex &Anmbar1, Complex &Ambarmbar1, Complex &Ambarmbar2, int const &m, double const &omega, GeodesicConstants &geoConstants, double const &rp, double const &thp, double const &Slm, double const &SlmP, double const &SlmPP){
	double a = geoConstants.a;
	double varpi = rp*rp + a*a;
	double delta = varpi - 2.*rp;
	double K = varpi*omega - m*a;
	double deltaP = 2.*(rp - 1.);
	double Kp = 2.*omega*rp;
	double cthp = cos(thp);
	double sthp = sin(thp);
	double Kdelta = K/delta;

	// q = -1/rho, so inverse powers of rho (and of rhoBar = conj(rho)) reduce
	// to products of q rather than calls to the much slower complex pow
	Complex q = rp - I*a*cthp;
	Complex q2 = q*q;
	Complex q3 = q2*q;
	Complex rho = -conj(q)/std::norm(q);
	Complex rhoBar = conj(rho);

	Complex L1 = -m/sthp + a*omega*sthp + cthp/sthp;
	Complex L2 = -m/sthp + a*omega*sthp + 2.*cthp/sthp;
	Complex L2S = SlmP + L2*Slm;
	Complex L2p = m*cthp/sthp/sthp + a*omega*cthp - 2./sthp/sthp;
	Complex L1Sp = SlmPP + L1*SlmP;
	Complex L1L2S = L1Sp + L2p*Slm + L2*SlmP + L1*L2*Slm;

	Complex q3sq2delta = q3/(M_SQRT2*delta);
	Complex q3rhoBarSlm = q3*rhoBar*Slm;

	Ann0 = q2*conj(q)*(0.5/(delta*delta))*( -q*L1L2S + 3.*I*a*sthp*L1*Slm + 3.*I*a*cthp*Slm + 2.*I*a*sthp*SlmP - I*a*sthp*L2*Slm );
	Anmbar0 = -q3sq2delta*( (rho + rhoBar - I*Kdelta)*L2S + (rho - rhoBar)*a*sthp*Kdelta*Slm );
	Anmbar1 = q3sq2delta*( L2S + I*(rho - rhoBar)*a*sthp*Slm );
	Ambarmbar0 = -0.25*q3rhoBarSlm*( I*(Kp/delta - deltaP/delta*Kdelta) + Kdelta*Kdelta + 2.*I*rho*Kdelta );
	Ambarmbar1 = 0.5*q3rhoBarSlm*( I*Kdelta - rho );
	Ambarmbar2 = 0.25*q3rhoBarSlm;
}

void A13_coeffs(Complex &All0, Complex &Alm0, Complex &Amm0, Complex &Alm1, Complex &Amm1, Complex &Amm2, int const &, int const &m, int const &k, int const &n, GeodesicConstants &geoConstants, double const &rp, double const &thp, double const &Slm, double const &SlmP, double const &SlmPP){
	double omega = (m*geoConstants.upsilonPhi + k*geoConstants.upsilonTheta + n*geoConstants.upsilonR)/geoConstants.upsilonT;
	A13_coeffs_w(All0, Alm0, Amm0, Alm1, Amm1, Amm2, m, omega, geoConstants, rp, thp, Slm, SlmP, SlmPP);
}

static void A13_coeffs_w(Complex &All0, Complex &Alm0, Complex &Amm0, Complex &Alm1, Complex &Amm1, Complex &Amm2, int const &m, double const &omega, GeodesicConstants &geoConstants, double const &rp, double const &thp, double const &Slm, double const &SlmP, double const &SlmPP){
	double a = geoConstants.a;
	double varpi = rp*rp + a*a;
	double delta = varpi - 2.*rp;
	double K = varpi*omega - m*a;
	double deltaP = 2.*(rp - 1.);
	double Kp = 2.*omega*rp;
	double cthp = cos(thp);
	double sthp = sin(thp);
	double Kdelta = K/delta;

	// q = -1/rho; inverse powers of rho and rhoBar = conj(rho) become
	// products of q, and |rho|^-2 = norm(q) is purely real
	Complex q = rp - I*a*cthp;
	double qq = std::norm(q);
	Complex rho = -conj(q)/qq;
	Complex rhoBar = conj(rho);

	Complex dLd1 = m/sthp - a*omega*sthp + cthp/sthp;
	Complex dLd2 = m/sthp - a*omega*sthp + 2.*cthp/sthp;
	Complex dLd2p = -m*cthp/sthp/sthp - a*omega*cthp - 2./sthp/sthp;
	Complex rhoOverRhoP = I*a*rho*sthp;
	Complex rhoOverRhoPP = I*a*rho*(cthp + 2.*sthp*rhoOverRhoP);

	All0 = 0.5*q*rhoBar*(SlmPP + (dLd1 + dLd2 + 2.*rhoOverRhoP)*SlmP + (dLd2p + dLd1*dLd2 - 6.*rhoOverRhoP*rhoOverRhoP + 3.*rhoOverRhoPP + (3.*dLd1 - dLd2)*rhoOverRhoP)*Slm);
	Alm0 = -M_SQRT2*q*( -1.*(rho + rhoBar + I*Kdelta)*(SlmP + dLd2*Slm) + (rho - rhoBar)*a*sthp*Kdelta*Slm);
	Alm1 = -M_SQRT2*q*( SlmP + dLd2*Slm + I*(rho - rhoBar)*a*sthp*Slm );
	Amm0 = -qq*Slm*( I*(Kp/delta - deltaP/delta*Kdelta) - Kdelta*Kdelta + 2.*I*rho*Kdelta );
	Amm1 = 2.*qq*Slm*( I*Kdelta + rho );
	Amm2 = -qq*Slm;
}

void u_13_coeffs(Complex &u1m, Complex &u1p, Complex &u3m, Complex &u3p, GeodesicConstants &geoConstants, double const &rp, double const &thp){
	double a = geoConstants.a;
	double uR = sqrt(std::abs(kerr_geo_Vr(a, geoConstants.En, geoConstants.Lz, geoConstants.Q, rp)));
	double uTheta = sqrt(std::abs(kerr_geo_Vtheta(a, geoConstants.En, geoConstants.Lz, geoConstants.Q, thp)));

	double z = cos(thp);
	Complex qb = rp + I*a*z;
	Complex rhobar = -conj(qb)/std::norm(qb);

	double Evarpi = geoConstants.En*(rp*rp + a*a) - a*geoConstants.Lz;
	double delta = rp*rp - 2.*rp + a*a;
	Complex thbracket = I*sqrt(1. - z*z)*(a*geoConstants.En - geoConstants.Lz/(1 - z*z));
	Complex rhobarSqrt2 = rhobar/M_SQRT2;
	u1p = -(Evarpi - uR)/delta;
	u1m = -(Evarpi + uR)/delta;
	u3p = rhobarSqrt2*(thbracket - uTheta);
	u3m = rhobarSqrt2*(thbracket + uTheta);
}

void u_13_coeffs_RadialTurningPoint(Complex &u1,  Complex &u3m, Complex &u3p, GeodesicConstants &geoConstants, double const &rp, double const &thp){
	double a = geoConstants.a;
	double uTheta = sqrt(std::abs(kerr_geo_Vtheta(a, geoConstants.En, geoConstants.Lz, geoConstants.Q, thp)));

	double z = cos(thp);
	Complex rhobar = -1./(rp + I*a*z);

	u1 = -(geoConstants.En*(rp*rp + a*a) - a*geoConstants.Lz)/(rp*rp - 2.*rp + a*a);
	u3p = rhobar/sqrt(2.)*(I*sqrt(1. - z*z)*(a*geoConstants.En - geoConstants.Lz/(1 - z*z)) - uTheta);
	u3m = rhobar/sqrt(2.)*(I*sqrt(1. - z*z)*(a*geoConstants.En - geoConstants.Lz/(1 - z*z)) + uTheta);
}

void u_13_coeffs_PolarTurningPoint(Complex &u1m, Complex &u1p, Complex &u3, GeodesicConstants &geoConstants, double const &rp, double const &thp){
	double a = geoConstants.a;
	double uR = sqrt(std::abs(kerr_geo_Vr(a, geoConstants.En, geoConstants.Lz, geoConstants.Q, rp)));

	double z = cos(thp);
	Complex rhobar = -1./(rp + I*a*z);

	u1p = -(geoConstants.En*(rp*rp + a*a) - a*geoConstants.Lz - uR)/(rp*rp - 2.*rp + a*a);
	u1m = -(geoConstants.En*(rp*rp + a*a) - a*geoConstants.Lz + uR)/(rp*rp - 2.*rp + a*a);
	u3 = rhobar/sqrt(2.)*(I*sqrt(1. - z*z)*(a*geoConstants.En - geoConstants.Lz/(1 - z*z)));
}

void u_13_coeffs_RadialPolarTurningPoint(Complex &u1, Complex &u3, GeodesicConstants &geoConstants, double const &rp, double const &thp){
	double a = geoConstants.a;
	double z = cos(thp);
	Complex rhobar = -1./(rp + I*a*z);

	u1 = -(geoConstants.En*(rp*rp + a*a) - a*geoConstants.Lz)/(rp*rp - 2.*rp + a*a);
	u3 = rhobar/sqrt(2.)*(I*sqrt(1. - z*z)*(a*geoConstants.En - geoConstants.Lz/(1 - z*z)));
}

void u_24_coeffs(Complex &u2m, Complex &u2p, Complex &u4m, Complex &u4p, GeodesicConstants &geoConstants, double const &rp, double const &thp){
	double a = geoConstants.a;
	double uR = sqrt(std::abs(kerr_geo_Vr(a, geoConstants.En, geoConstants.Lz, geoConstants.Q, rp)));
	double uTheta = sqrt(std::abs(kerr_geo_Vtheta(a, geoConstants.En, geoConstants.Lz, geoConstants.Q, thp)));
	double z = cos(thp);
	Complex qr0 = rp - I*a*z;
	Complex rho = -conj(qr0)/std::norm(qr0);

	double Evarpi = geoConstants.En*(rp*rp + a*a) - a*geoConstants.Lz;
	double twoSigma = 2.*(rp*rp + (a*z)*(a*z));
	Complex thbracket = I*sqrt(1. - z*z)*(a*geoConstants.En - geoConstants.Lz/(1 - z*z));
	Complex rhoSqrt2 = -rho/M_SQRT2;
	u2p = -(Evarpi + uR)/twoSigma;
	u2m = -(Evarpi - uR)/twoSigma;
	u4p = rhoSqrt2*(thbracket + uTheta);
	u4m = rhoSqrt2*(thbracket - uTheta);
}

void u_24_coeffs_PolarTurningPoint(Complex &u2m, Complex &u2p, Complex &u4, GeodesicConstants &geoConstants, double const &rp, double const &thp){
	double a = geoConstants.a;
	double uR = sqrt(std::abs(kerr_geo_Vr(a, geoConstants.En, geoConstants.Lz, geoConstants.Q, rp)));
	double z = cos(thp);
	Complex rho = -1./(rp - I*a*z);

	u2p = -(geoConstants.En*(rp*rp + a*a) - a*geoConstants.Lz + uR)/(2.*(rp*rp + pow(a*z, 2)));
	u2m = -(geoConstants.En*(rp*rp + a*a) - a*geoConstants.Lz - uR)/(2.*(rp*rp + pow(a*z, 2)));
	u4 = -rho/sqrt(2.)*(I*sqrt(1. - z*z)*(a*geoConstants.En - geoConstants.Lz/(1 - z*z)));
}

void u_24_coeffs_RadialTurningPoint(Complex &u2, Complex &u4m, Complex &u4p, GeodesicConstants &geoConstants, double const &rp, double const &thp){
	double a = geoConstants.a;
	double uTheta = sqrt(std::abs(kerr_geo_Vtheta(a, geoConstants.En, geoConstants.Lz, geoConstants.Q, thp)));
	double z = cos(thp);
	Complex rho = -1./(rp - I*a*z);

	u2 = -(geoConstants.En*(rp*rp + a*a) - a*geoConstants.Lz)/(2.*(rp*rp + pow(a*z, 2)));
	u4p = rho/sqrt(2.)*(-I*sqrt(1. - z*z)*(a*geoConstants.En - geoConstants.Lz/(1 - z*z)) - uTheta);
	u4m = rho/sqrt(2.)*(-I*sqrt(1. - z*z)*(a*geoConstants.En - geoConstants.Lz/(1 - z*z)) + uTheta);
}

void u_24_coeffs_RadialPolarTurningPoint(Complex &u2, Complex &u4, GeodesicConstants &geoConstants, double const &rp, double const &thp){
	double a = geoConstants.a;
	double z = cos(thp);
	Complex rho = -1./(rp - I*a*z);

	u2 = -(geoConstants.En*(rp*rp + a*a) - a*geoConstants.Lz)/(2.*(rp*rp + pow(a*z, 2)));
	u4 = rho/sqrt(2.)*(-I*sqrt(1. - z*z)*(a*geoConstants.En - geoConstants.Lz/(1 - z*z)));
}

double u_1(GeodesicConstants &geoConstants, double const &rp, int const &sgnUr, double const &thp, int const &){
	double uR;
	double a = geoConstants.a;
	if( sgnUr == 0) uR = 0.; else uR = sgnUr*sqrt(kerr_geo_Vr(a, geoConstants.En, geoConstants.Lz, geoConstants.Q, rp));
	if( isnan(std::abs(uR)) ) uR = 0.;

	return -(geoConstants.En*(rp*rp + a*a) - a*geoConstants.Lz - uR)/(rp*rp - 2.*rp + a*a);
}

double u_2(GeodesicConstants &geoConstants, double const &rp, int const &sgnUr, double const &thp, int const &){
	double uR;
	double a = geoConstants.a;
	if( sgnUr == 0) uR = 0.; else uR = sgnUr*sqrt(kerr_geo_Vr(a, geoConstants.En, geoConstants.Lz, geoConstants.Q, rp));
	if( isnan(std::abs(uR)) ) uR = 0.;

	return -(geoConstants.En*(rp*rp + a*a) - a*geoConstants.Lz + uR)/(2.*(rp*rp + pow(a*cos(thp), 2)));
}

Complex u_3(GeodesicConstants &geoConstants, double const &rp, int const &, double const &thp, int const &sgnUth){
	double uTheta;
	double a = geoConstants.a;
	double z = cos(thp);
	if( sgnUth == 0) uTheta = 0.; else uTheta = sgnUth*sqrt(kerr_geo_Vtheta(geoConstants.a, geoConstants.En, geoConstants.Lz, geoConstants.Q, thp));
	if( isnan(uTheta) ) uTheta = 0.;
	Complex rhobar = -1./(rp + I*a*z);

	return rhobar/sqrt(2.)*(I*sqrt(1. - z*z)*(a*geoConstants.En - geoConstants.Lz/(1 - z*z)) - uTheta);
}

Complex u_4(GeodesicConstants &geoConstants, double const &rp, int const &sgnUr, double const &thp, int const &sgnUth){
	return std::conj(u_3(geoConstants, rp, sgnUr, thp, sgnUth));
}

double C_ll(GeodesicConstants &geoConstants, double const &rp, int const &sgnUr, double const &thp, int const &sgnUth){
	double u1 = u_1(geoConstants, rp, sgnUr, thp, sgnUth);
	return u1*u1;
}

Complex C_lm(GeodesicConstants &geoConstants, double const &rp, int const &sgnUr, double const &thp, int const &sgnUth){
	double u1 = u_1(geoConstants, rp, sgnUr, thp, sgnUth);
	Complex u3 = u_3(geoConstants, rp, sgnUr, thp, sgnUth);
	return u1*u3;
}

Complex C_mm(GeodesicConstants &geoConstants, double const &rp, int const &sgnUr, double const &thp, int const &sgnUth){
	Complex u3 = u_3(geoConstants, rp, sgnUr, thp, sgnUth);
	return u3*u3;
}

double C_nn(GeodesicConstants &geoConstants, double const &rp, int const &sgnUr, double const &thp, int const &){
	double uR;
	if( sgnUr == 0) uR = 0.; else uR = sgnUr*sqrt(kerr_geo_Vr(geoConstants.a, geoConstants.En, geoConstants.Lz, geoConstants.Q, rp));
	if( isnan(std::abs(uR)) ) uR = 0.;

	return pow((geoConstants.En*(rp*rp + geoConstants.a*geoConstants.a) - geoConstants.a*geoConstants.Lz + uR)/(2.*(rp*rp + geoConstants.a*geoConstants.a*cos(thp)*cos(thp))), 2);
}

Complex C_nmbar(GeodesicConstants &geoConstants, double const &rp, int const &sgnUr, double const &thp, int const &sgnUth){
	double uR, uTheta;
	if( sgnUr == 0) uR = 0.; else uR = sgnUr*sqrt(kerr_geo_Vr(geoConstants.a, geoConstants.En, geoConstants.Lz, geoConstants.Q, rp));
	if( isnan(std::abs(uR)) ) uR = 0.;
	if( sgnUth == 0) uTheta = 0.; else uTheta = sgnUth*sqrt(kerr_geo_Vtheta(geoConstants.a, geoConstants.En, geoConstants.Lz, geoConstants.Q, thp));
	if( isnan(uTheta) ) uTheta = 0.;

	return (geoConstants.En*(rp*rp + geoConstants.a*geoConstants.a) - geoConstants.a*geoConstants.Lz + uR)/(2.*(rp*rp + geoConstants.a*geoConstants.a*cos(thp)*cos(thp)))*(-1./(rp - I*geoConstants.a*cos(thp)))*(I*sin(thp)*(geoConstants.a*geoConstants.En - geoConstants.Lz*pow(sin(thp), -2)) + uTheta)/sqrt(2.);
}

Complex C_mbarmbar(GeodesicConstants &geoConstants, double const &rp, int const &, double const &thp, int const &sgnUth){
	double uTheta;
	if( sgnUth == 0) uTheta = 0.; else uTheta = sgnUth*sqrt(kerr_geo_Vtheta(geoConstants.a, geoConstants.En, geoConstants.Lz, geoConstants.Q, thp));
	if( isnan(uTheta) ) uTheta = 0.;

	return pow((-1./(rp - I*geoConstants.a*cos(thp)))*(I*sin(thp)*(geoConstants.a*geoConstants.En - geoConstants.Lz*pow(sin(thp), -2)) + uTheta)/sqrt(2.), 2);
}

void C_coeffs(Complex &Cnn, Complex &Cnmbar, Complex &Cmbarmbar, GeodesicConstants &geoConstants, double const &rp, int const &sgnUr, double const &thp, int const &sgnUth){
	double uR, uTheta;
	if( sgnUr == 0) uR = 0.; else uR = sgnUr*sqrt(kerr_geo_Vr(geoConstants.a, geoConstants.En, geoConstants.Lz, geoConstants.Q, rp));
	if( isnan(std::abs(uR)) ) uR = 0.;
	if( sgnUth == 0) uTheta = 0.; else uTheta = sgnUth*sqrt(kerr_geo_Vtheta(geoConstants.a, geoConstants.En, geoConstants.Lz, geoConstants.Q, thp));
	if( isnan(uTheta) ) uTheta = 0.;

	Cnn = pow((geoConstants.En*(rp*rp + geoConstants.a*geoConstants.a) - geoConstants.a*geoConstants.Lz + uR)/(2.*(rp*rp + geoConstants.a*geoConstants.a*cos(thp)*cos(thp))), 2);
	Cnmbar = (geoConstants.En*(rp*rp + geoConstants.a*geoConstants.a) - geoConstants.a*geoConstants.Lz + uR)/(2.*(rp*rp + geoConstants.a*geoConstants.a*cos(thp)*cos(thp)))*(-1./(rp - I*geoConstants.a*cos(thp)))*(I*sin(thp)*(geoConstants.a*geoConstants.En - geoConstants.Lz*pow(sin(thp), -2)) + uTheta)/sqrt(2.);
	Cmbarmbar = pow((-1./(rp - I*geoConstants.a*cos(thp)))*(I*sin(thp)*(geoConstants.a*geoConstants.En - geoConstants.Lz*pow(sin(thp), -2)) + uTheta)/sqrt(2.), 2);
}

void C13_coeffs(Complex &Cll, Complex &Clm, Complex &Cmm, GeodesicConstants &geoConstants, double const &rp, int const &sgnUr, double const &thp, int const &sgnUth){
	double uR, uTheta;
	double a = geoConstants.a;
	double z = cos(thp);
	Complex rhobar = -1./(rp + I*a*z);

	if( sgnUr == 0) uR = 0.; else uR = sgnUr*sqrt(kerr_geo_Vr(a, geoConstants.En, geoConstants.Lz, geoConstants.Q, rp));
	if( isnan(std::abs(uR)) ) uR = 0.;
	if( sgnUth == 0) uTheta = 0.; else uTheta = sgnUth*sqrt(kerr_geo_Vtheta(a, geoConstants.En, geoConstants.Lz, geoConstants.Q, thp));
	if( isnan(uTheta) ) uTheta = 0.;
	Complex u1 = -(geoConstants.En*(rp*rp + a*a) - a*geoConstants.Lz - uR)/(rp*rp - 2.*rp + a*a);
	Complex u3 = rhobar/sqrt(2.)*(I*sqrt(1. - z*z)*(a*geoConstants.En - geoConstants.Lz/(1 - z*z)) - uTheta);

	Cll = u1*u1;
	Clm = u1*u3;
	Cmm = u3*u3;
}

///////////////////////////////////
// Scalar s = 0 field amplitudes //
///////////////////////////////////

TeukolskyAmplitudes scalar_amplitude_circeq(int L, int m, GeodesicTrajectory& traj, GeodesicConstants &geoConstants, ComplexDerivativesMatrixStruct Rin, ComplexDerivativesMatrixStruct Rup, DerivativesMatrix Slm){
	
	int k = 0, n = 0;

	Complex R0 = (Rin.solution)[0];
	Complex Rp0 = (Rin.derivative)[0];

	Complex R1 = (Rup.solution)[0];
	Complex Rp1 = (Rup.derivative)[0];

	double St = (Slm.solution)[0];

	double rp = (traj.r)[0];
	double thp = 0.5*M_PI;

	Complex W = scalar_wronskian(geoConstants.a, rp, R0, Rp0, R1, Rp1);

	Complex I1Up = scalar_integrand_1(L, m, k, n, geoConstants, 0., rp, 0., 0., R0);
	Complex I1In = scalar_integrand_1(L, m, k, n, geoConstants, 0., rp, 0., 0., R1);
	Complex I2 = scalar_integrand_2(L, m, k, n, geoConstants, 0., thp, 0., 0., St);

	Complex I3Up = scalar_integrand_3(L, m, k, n, geoConstants, 0., rp, 0., 0., R0);
	Complex I3In = scalar_integrand_3(L, m, k, n, geoConstants, 0., rp, 0., 0., R1);
	Complex I4 = scalar_integrand_4(L, m, k, n, geoConstants, 0., thp, 0., 0., St);

	Complex ZlmUp = -4.*M_PI/W/geoConstants.upsilonT*( I1Up*I2 + I3Up*I4 );
	Complex ZlmIn = -4.*M_PI/W/geoConstants.upsilonT*( I1In*I2 + I3In*I4 );

	TeukolskyAmplitudes Zlm = {ZlmIn, ZlmUp, DBL_EPSILON, DBL_EPSILON};

	return Zlm;
}

int scalar_integrand_I1(Complex &integrand, int m, int n, double freq, double tR, double rp, double phiR, double qr, Complex Rt){
	integrand = rp*rp*Rt*cos(n*qr + freq*tR - m*phiR);
	return 0;
}

int scalar_integrand_I2(Complex &integrand, int m, int k, double freq, double tTh, double, double phiTh, double qth, double St){
	integrand = St*cos(k*qth + freq*tTh - m*phiTh);
	return 0;
}

int scalar_integrand_I3(Complex &integrand, int m, int n, double freq, double tR, double, double phiR, double qr, Complex Rt){
	integrand = Rt*cos(n*qr + freq*tR - m*phiR);
	return 0;
}

int scalar_integrand_I4(Complex &integrand, int m, int k, double freq, double tTh, double aCosThP, double phiTh, double qth, double St){
	integrand = aCosThP*aCosThP*St*cos(k*qth + freq*tTh - m*phiTh);
	return 0;
}

int radial_integral_convergence_sum(Complex &II, int (*integrand)(Complex &, int, int, double, double, double, double, double, Complex), int m, int k, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConstants, BoundaryCondition bc, RadialTeukolsky &teuk, double errorThreshold, double errorTolerance){
	int halfSampleInit = 16;
	int halfSample = halfSampleInit;
	int argLength = traj.r.size();
	int halfSampleMax = argLength - 1;
	while(halfSample < std::abs(n) + 1){
		halfSample *= 2;
	}
	if(halfSampleMax < halfSample){
		halfSample = halfSampleMax;
	}
	int sampleDiff = halfSampleMax/halfSample;
	double deltaQ = M_PI/double(halfSampleMax);

	Complex sumTerm = 0.;
	Complex sum = sumTerm;
	double maxTerm = 0.;
	double freq = teuk.getModeFrequency();

	int samplePos = 0;
	double q = samplePos*deltaQ;
	integrand(sumTerm, m, n, freq, traj.tR[samplePos], teuk.getRadialPoints(samplePos), traj.phiR[samplePos], q, teuk.getSolution(bc, samplePos));
	sum = sumTerm;
	maxTerm = std::abs(sumTerm) < maxTerm ? maxTerm : std::abs(sumTerm);

	samplePos = halfSample*sampleDiff;
	q = samplePos*deltaQ;
	integrand(sumTerm, m, n, freq, traj.tR[samplePos], teuk.getRadialPoints(samplePos), traj.phiR[samplePos], q, teuk.getSolution(bc, samplePos));
	sum += sumTerm;
	maxTerm = std::abs(sumTerm) < maxTerm ? maxTerm : std::abs(sumTerm);

	for(int i = 1; i < halfSample; i++){
		samplePos = i*sampleDiff;
		q = samplePos*deltaQ;

		integrand(sumTerm, m, n, freq, traj.tR[samplePos], teuk.getRadialPoints(samplePos), traj.phiR[samplePos], q, teuk.getSolution(bc, samplePos));
		sum += 2.*sumTerm;
		maxTerm = std::abs(sumTerm) < maxTerm ? maxTerm : std::abs(sumTerm);
	}

	II = sum/double(2*halfSample);

	double precisionLoss = maxTerm/std::abs(II);
	double errorToleranceAdjusted = 10*DBL_EPSILON*precisionLoss;
	errorToleranceAdjusted = errorToleranceAdjusted < errorTolerance ? errorTolerance : errorToleranceAdjusted;

	Complex ICompare = 0.;
	while(halfSample < halfSampleMax && std::abs(1. - ICompare/II) > errorToleranceAdjusted){
		for(int i = 0; i < halfSample; i++){
			samplePos = i*sampleDiff + sampleDiff/2;
			q = double(samplePos)*deltaQ;

			integrand(sumTerm, m, n, freq, traj.tR[samplePos], teuk.getRadialPoints(samplePos), traj.phiR[samplePos], q, teuk.getSolution(bc, samplePos));
			sum += 2.*sumTerm;
			maxTerm = std::abs(sumTerm) < maxTerm ? maxTerm : std::abs(sumTerm);
		}
		halfSample *= 2;
		sampleDiff /= 2;
		ICompare = II;
		II = sum/double(2*halfSample);

		precisionLoss = maxTerm/std::abs(II);
		errorToleranceAdjusted = 10*DBL_EPSILON*precisionLoss;
		errorToleranceAdjusted = errorToleranceAdjusted < errorTolerance ? errorTolerance : errorToleranceAdjusted;
	}
	if(errorToleranceAdjusted > errorThreshold){
		return -1;
	}
	if(std::abs(1. - ICompare/II) > errorToleranceAdjusted){
		std::cout << "(SOURCEINT) ERROR: IR ("<<m<<","<<k<<","<<n<<") integral did not converge to expected tolerance of "<<errorToleranceAdjusted<<" within N = " << 2*halfSample << " samples. Only converged to precision of  "<<std::abs(1. - ICompare/II)<< ". \n";
		return 1;
	}

	return 0;
}

int polar_integral_convergence_sum(Complex &II, int (*integrand)(Complex &, int, int, double, double, double, double, double, double), int m, int k, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConstants, SpinWeightedHarmonic &swsh, double errorThreshold, double errorTolerance){
	int halfSampleInit = 16;
	int halfSample = halfSampleInit;
	int argLength = traj.theta.size();
	int halfSampleMax = argLength - 1;
	while(halfSample < std::abs(k) + 1){
		halfSample *= 2;
	}
	if(halfSampleMax < halfSample){
		halfSample = halfSampleMax;
	}
	int sampleDiff = halfSampleMax/halfSample;
	double deltaQ = M_PI/double(halfSampleMax);

	Complex sumTerm = 0.;
	Complex sum = sumTerm;
	double maxTerm = 0.;
	double freq = geoConstants.getTimeFrequency(m, k, n);
	double a = geoConstants.a;

	int samplePos = 0;
	double q = samplePos*deltaQ;
	integrand(sumTerm, m, k, freq, traj.tTheta[samplePos], a*cos(swsh.getArguments(samplePos)), traj.phiTheta[samplePos], q, swsh.getSolution(samplePos));
	sum = sumTerm;
	maxTerm = std::abs(sumTerm) < maxTerm ? maxTerm : std::abs(sumTerm);

	samplePos = halfSample*sampleDiff;
	q = samplePos*deltaQ;
	integrand(sumTerm, m, k, freq, traj.tTheta[samplePos], a*cos(swsh.getArguments(samplePos)), traj.phiTheta[samplePos], q, swsh.getSolution(samplePos));
	sum += sumTerm;
	maxTerm = std::abs(sumTerm) < maxTerm ? maxTerm : std::abs(sumTerm);

	for(int i = 1; i < halfSample; i++){
		samplePos = i*sampleDiff;
		q = samplePos*deltaQ;

		integrand(sumTerm, m, k, freq, traj.tTheta[samplePos], a*cos(swsh.getArguments(samplePos)), traj.phiTheta[samplePos], q, swsh.getSolution(samplePos));
		sum += 2.*sumTerm;
		maxTerm = std::abs(sumTerm) < maxTerm ? maxTerm : std::abs(sumTerm);
	}

	II = sum/double(2*halfSample);

	double precisionLoss = maxTerm/std::abs(II);
	double errorToleranceAdjusted = 10*DBL_EPSILON*precisionLoss;
	errorToleranceAdjusted = errorToleranceAdjusted < errorTolerance ? errorTolerance : errorToleranceAdjusted;

	Complex ICompare = 0.;
	while(halfSample < halfSampleMax && std::abs(1. - ICompare/II) > errorToleranceAdjusted){
		for(int i = 0; i < halfSample; i++){
			samplePos = i*sampleDiff + sampleDiff/2;
			q = double(samplePos)*deltaQ;

			integrand(sumTerm, m, k, freq, traj.tTheta[samplePos], a*cos(swsh.getArguments(samplePos)), traj.phiTheta[samplePos], q, swsh.getSolution(samplePos));
			sum += 2.*sumTerm;
			maxTerm = std::abs(sumTerm) < maxTerm ? maxTerm : std::abs(sumTerm);
		}
		halfSample *= 2;
		sampleDiff /= 2;
		ICompare = II;
		II = sum/double(2*halfSample);

		precisionLoss = maxTerm/std::abs(II);
		errorToleranceAdjusted = 10*DBL_EPSILON*precisionLoss;
		errorToleranceAdjusted = errorToleranceAdjusted < errorTolerance ? errorTolerance : errorToleranceAdjusted;
	}
	if(errorToleranceAdjusted > errorThreshold){
		return -1;
	}

	if(std::abs(1. - ICompare/II) > errorToleranceAdjusted){
		std::cout << "(SOURCEINT) ERROR: ITh ("<<m<<","<<k<<","<<n<<") integral did not converge to expected tolerance of "<<errorToleranceAdjusted<<" within N = " << 2*halfSample << " samples. Only converged to precision of  "<<std::abs(1. - ICompare/II)<< ". \n";
		return 1;
	}

	return 0;
}

TeukolskyAmplitudes scalar_amplitude(int l, int m, int k, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConstants, RadialTeukolsky &teuk, SpinWeightedHarmonic &swsh){
	if(std::abs(geoConstants.x) == 1. && geoConstants.e == 0.){
		return scalar_amplitude_circular(l, m, k, n, traj, geoConstants, teuk, swsh);
	}else if(std::abs(geoConstants.x) == 1.){
		return scalar_amplitude_equatorial(l, m, k, n, traj, geoConstants, teuk, swsh);
	}else if(geoConstants.e == 0.){
		return scalar_amplitude_spherical(l, m, k, n, traj, geoConstants, teuk, swsh);
	}
	return scalar_amplitude_generic(l, m, k, n, traj, geoConstants, teuk, swsh);
}

TeukolskyAmplitudes scalar_amplitude_generic(int, int m, int k, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConstants, RadialTeukolsky &teuk, SpinWeightedHarmonic &swsh){
	double errorThresholdR = 1.e-2;
	double errorThresholdTh = 1.e-2;
	double errorTolerance = 5.e-12;
	double upT = geoConstants.upsilonT;
	Complex W = scalar_wronskian(geoConstants.a, teuk.getRadialPoints(0), teuk.getSolution(In, 0), teuk.getDerivative(In, 0), teuk.getSolution(Up, 0), teuk.getDerivative(Up, 0));

	Complex I1Up = 0., I1In = 0., I3Up = 0., I3In = 0., I2 = 0., I4 = 0.;
	Complex ZlmUp = 0.;
	Complex ZlmIn = 0.;

	int status = polar_integral_convergence_sum(I2, scalar_integrand_I2, m, k, n, traj, geoConstants, swsh, errorThresholdTh, errorTolerance);
	if(status == -1){
		TeukolskyAmplitudes Zlm = {ZlmIn, ZlmUp, DBL_EPSILON, DBL_EPSILON};
		return Zlm;
	}
	status = polar_integral_convergence_sum(I4, scalar_integrand_I4, m, k, n, traj, geoConstants, swsh, errorThresholdTh, errorTolerance);
	if(status == -1){
		TeukolskyAmplitudes Zlm = {ZlmIn, ZlmUp, DBL_EPSILON, DBL_EPSILON};
		return Zlm;
	}

	status = radial_integral_convergence_sum(I1In, scalar_integrand_I1, m, k, n, traj, geoConstants, Up, teuk, errorThresholdR, errorTolerance);
	if(status != -1){
		if(std::abs(I4) > 0.){
			status = radial_integral_convergence_sum(I3In, scalar_integrand_I3, m, k, n, traj, geoConstants, Up, teuk, errorThresholdR, errorTolerance);
		}
		if(status != -1){
			ZlmIn = -4.*M_PI/W/upT*(I1In*I2 + I3In*I4);
		}
	}

	status = radial_integral_convergence_sum(I1Up, scalar_integrand_I1, m, k, n, traj, geoConstants, In, teuk, errorThresholdR, errorTolerance);
	if(status != -1){
		if(std::abs(I4) > 0.){
			status = radial_integral_convergence_sum(I3Up, scalar_integrand_I3, m, k, n, traj, geoConstants, In, teuk, errorThresholdR, errorTolerance);
		}
		if(status != -1){
			ZlmUp = -4.*M_PI/W/upT*(I1Up*I2 + I3Up*I4);
		}
	}

	TeukolskyAmplitudes Zlm = {ZlmIn, ZlmUp, DBL_EPSILON, DBL_EPSILON};

	return Zlm;
}

TeukolskyAmplitudes scalar_amplitude_equatorial(int, int m, int k, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConstants, RadialTeukolsky &teuk, SpinWeightedHarmonic &swsh){
	double errorThresholdR = 1.e-2;
	double errorTolerance = 5.e-12;
	double upT = geoConstants.upsilonT;
	Complex W = scalar_wronskian(geoConstants.a, teuk.getRadialPoints(0), teuk.getSolution(In, 0), teuk.getDerivative(In, 0), teuk.getSolution(Up, 0), teuk.getDerivative(Up, 0));

	Complex I1Up = 0., I1In = 0., I3Up = 0., I3In = 0., I2 = 0., I4 = 0.;
	Complex ZlmUp = 0.;
	Complex ZlmIn = 0.;

	int status = scalar_integrand_I2(I2, m, k, teuk.getModeFrequency(), traj.tTheta[0], geoConstants.a*cos(swsh.getArguments(0)), traj.getAzimuthalAccumulation(2, 0), 0, swsh.getSolution(0));
	if(status == -1){
		TeukolskyAmplitudes Zlm = {ZlmIn, ZlmUp, DBL_EPSILON, DBL_EPSILON};
		return Zlm;
	}
	status = scalar_integrand_I4(I4, m, k, teuk.getModeFrequency(), traj.tTheta[0], geoConstants.a*cos(traj.getPolarPosition(0)), traj.getAzimuthalAccumulation(2, 0), 0, swsh.getSolution(0));
	if(status == -1){
		TeukolskyAmplitudes Zlm = {ZlmIn, ZlmUp, DBL_EPSILON, DBL_EPSILON};
		return Zlm;
	}

	status = radial_integral_convergence_sum(I1In, scalar_integrand_I1, m, k, n, traj, geoConstants, Up, teuk, errorThresholdR, errorTolerance);
	if(status != -1){
		if(std::abs(I4) > 0.){
			status = radial_integral_convergence_sum(I3In, scalar_integrand_I3, m, k, n, traj, geoConstants, Up, teuk, errorThresholdR, errorTolerance);
		}
		if(status != -1){
			ZlmIn = -4.*M_PI/W/upT*(I1In*I2 + I3In*I4);
		}
	}

	status = radial_integral_convergence_sum(I1Up, scalar_integrand_I1, m, k, n, traj, geoConstants, In, teuk, errorThresholdR, errorTolerance);
	if(status != -1){
		if(std::abs(I4) > 0.){
			status = radial_integral_convergence_sum(I3Up, scalar_integrand_I3, m, k, n, traj, geoConstants, In, teuk, errorThresholdR, errorTolerance);
		}
		if(status != -1){
			ZlmUp = -4.*M_PI/W/upT*(I1Up*I2 + I3Up*I4);
		}
	}

	TeukolskyAmplitudes Zlm = {ZlmIn, ZlmUp, DBL_EPSILON, DBL_EPSILON};

	return Zlm;
}

TeukolskyAmplitudes scalar_amplitude_spherical(int, int m, int k, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConstants, RadialTeukolsky &teuk, SpinWeightedHarmonic &swsh){
	double errorThresholdTh = 1.e-2;
	double errorTolerance = 5.e-12;
	double upT = geoConstants.upsilonT;
	Complex W = scalar_wronskian(geoConstants.a, teuk.getRadialPoints(0), teuk.getSolution(In, 0), teuk.getDerivative(In, 0), teuk.getSolution(Up, 0), teuk.getDerivative(Up, 0));

	Complex I1Up = 0., I1In = 0., I3Up = 0., I3In = 0., I2 = 0., I4 = 0.;
	Complex ZlmUp = 0.;
	Complex ZlmIn = 0.;

	int status = polar_integral_convergence_sum(I2, scalar_integrand_I2, m, k, n, traj, geoConstants, swsh, errorThresholdTh, errorTolerance);
	if(status == -1){
		TeukolskyAmplitudes Zlm = {ZlmIn, ZlmUp, DBL_EPSILON, DBL_EPSILON};
		return Zlm;
	}
	status = polar_integral_convergence_sum(I4, scalar_integrand_I4, m, k, n, traj, geoConstants, swsh, errorThresholdTh, errorTolerance);
	if(status == -1){
		TeukolskyAmplitudes Zlm = {ZlmIn, ZlmUp, DBL_EPSILON, DBL_EPSILON};
		return Zlm;
	}

	status = scalar_integrand_I1(I1In, m, n, teuk.getModeFrequency(), traj.getTimeAccumulation(1, 0), traj.getRadialPosition(0), traj.getAzimuthalAccumulation(1, 0), 0., teuk.getSolution(Up, 0));
	if(status != -1){
		status = scalar_integrand_I3(I3In, m, n, teuk.getModeFrequency(), traj.getTimeAccumulation(1, 0), traj.getRadialPosition(0), traj.getAzimuthalAccumulation(1, 0), 0, teuk.getSolution(Up, 0));
		if(status != -1){
			ZlmIn = -4.*M_PI/W/upT*(I1In*I2 + I3In*I4);
		}
	}

	status = scalar_integrand_I1(I1Up, m, n, teuk.getModeFrequency(), traj.getTimeAccumulation(1, 0), traj.getRadialPosition(0), traj.getAzimuthalAccumulation(1, 0), 0., teuk.getSolution(In, 0));
	if(status != -1){
		status = scalar_integrand_I3(I3Up, m, n, teuk.getModeFrequency(), traj.getTimeAccumulation(1, 0), traj.getRadialPosition(0), traj.getAzimuthalAccumulation(1, 0), 0, teuk.getSolution(In, 0));
		if(status != -1){
			ZlmUp = -4.*M_PI/W/upT*(I1Up*I2 + I3Up*I4);
		}
	}

	TeukolskyAmplitudes Zlm = {ZlmIn, ZlmUp, DBL_EPSILON, DBL_EPSILON};

	return Zlm;
}

TeukolskyAmplitudes scalar_amplitude_circular(int, int m, int k, int n, GeodesicTrajectory& traj, GeodesicConstants &geoConstants, RadialTeukolsky &teuk, SpinWeightedHarmonic &swsh){
	double upT = geoConstants.upsilonT;
	Complex W = scalar_wronskian(geoConstants.a, teuk.getRadialPoints(0), teuk.getSolution(In, 0), teuk.getDerivative(In, 0), teuk.getSolution(Up, 0), teuk.getDerivative(Up, 0));

	Complex I1Up = 0., I1In = 0., I3Up = 0., I3In = 0., I2 = 0., I4 = 0.;
	Complex ZlmUp = 0.;
	Complex ZlmIn = 0.;

	int status = scalar_integrand_I2(I2, m, k, teuk.getModeFrequency(), traj.tTheta[0], geoConstants.a*cos(traj.getPolarPosition(0)), traj.getAzimuthalAccumulation(2, 0), 0, swsh.getSolution(0));
	if(status == -1){
		TeukolskyAmplitudes Zlm = {ZlmIn, ZlmUp, DBL_EPSILON, DBL_EPSILON};
		return Zlm;
	}
	status = scalar_integrand_I4(I4, m, k, teuk.getModeFrequency(), traj.tTheta[0], geoConstants.a*cos(traj.getPolarPosition(0)), traj.getAzimuthalAccumulation(2, 0), 0, swsh.getSolution(0));
	if(status == -1){
		TeukolskyAmplitudes Zlm = {ZlmIn, ZlmUp, DBL_EPSILON, DBL_EPSILON};
		return Zlm;
	}

	status = scalar_integrand_I1(I1In, m, n, teuk.getModeFrequency(), traj.getTimeAccumulation(1, 0), traj.getRadialPosition(0), traj.getAzimuthalAccumulation(1, 0), 0, teuk.getSolution(Up, 0));
	if(status != -1){
		status = scalar_integrand_I1(I3In, m, n, teuk.getModeFrequency(), traj.getTimeAccumulation(1, 0), traj.getRadialPosition(0), traj.getAzimuthalAccumulation(1, 0), 0, teuk.getSolution(Up, 0));
		if(status != -1){
			ZlmIn = -4.*M_PI/W/upT*(I1In*I2 + I3In*I4);
		}
	}

	status = scalar_integrand_I1(I1Up, m, n, teuk.getModeFrequency(), traj.getTimeAccumulation(1, 0), traj.getRadialPosition(0), traj.getAzimuthalAccumulation(1, 0), 0, teuk.getSolution(In, 0));
	if(status != -1){
		status = scalar_integrand_I3(I3Up, m, n, teuk.getModeFrequency(), traj.getTimeAccumulation(1, 0), traj.getRadialPosition(0), traj.getAzimuthalAccumulation(1, 0), 0, teuk.getSolution(In, 0));
		if(status != -1){
			ZlmUp = -4.*M_PI/W/upT*(I1Up*I2 + I3Up*I4);
		}
	}

	TeukolskyAmplitudes Zlm = {ZlmIn, ZlmUp, DBL_EPSILON, DBL_EPSILON};

	return Zlm;
}

Complex scalar_integrand_1(int, int m, int k, int n, GeodesicConstants &geoConstants, double tR, double rp, double phiR, double qr, Complex Rt){
	double freq = (m*geoConstants.upsilonPhi + k*geoConstants.upsilonTheta + n*geoConstants.upsilonR)/geoConstants.upsilonT;
	return pow(rp, 2)*Rt*exp(I*(n*qr + freq*tR - m*phiR));
}

Complex scalar_integrand_2(int, int m, int k, int n, GeodesicConstants &geoConstants, double tTh, double, double phiTh, double qth, double St){
	double freq = (m*geoConstants.upsilonPhi + k*geoConstants.upsilonTheta + n*geoConstants.upsilonR)/geoConstants.upsilonT;
	return St*exp(I*(k*qth + freq*tTh - m*phiTh));
}

Complex scalar_integrand_3(int, int m, int k, int n, GeodesicConstants &geoConstants, double tR, double, double phiR, double qr, Complex Rt){
	double freq = (m*geoConstants.upsilonPhi + k*geoConstants.upsilonTheta + n*geoConstants.upsilonR)/geoConstants.upsilonT;
	return Rt*exp(I*(n*qr + freq*tR - m*phiR));
}

Complex scalar_integrand_4(int, int m, int k, int n, GeodesicConstants &geoConstants, double tTh, double thp, double phiTh, double qth, double St){
	double freq = (m*geoConstants.upsilonPhi + k*geoConstants.upsilonTheta + n*geoConstants.upsilonR)/geoConstants.upsilonT;
	return pow(geoConstants.a*cos(thp), 2)*St*exp(I*(k*qth + freq*tTh - m*phiTh));
}

Complex wronskian(double a, double rp, Complex Rin, Complex RinP, Complex Rup, Complex RupP){
	return (Rin*RupP - Rup*RinP)/(rp*rp - 2.*rp + a*a);
}

Complex wronskian(int s, double a, double rp, Complex Rin, Complex RinP, Complex Rup, Complex RupP){
	return (Rin*RupP - Rup*RinP)*pow(rp*rp - 2.*rp + a*a, s + 1);
}

Complex scalar_wronskian(double a, double rp, Complex Rin, Complex RinP, Complex Rup, Complex RupP){
	return (Rin*RupP - Rup*RinP)*(rp*rp - 2.*rp + a*a);
}
