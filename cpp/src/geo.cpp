// geo.cpp

#include "geo.hpp"

#define FOURIER_EPS 1.e-13
#define FREQ_EPS 1.e-14

GeodesicTrajectory::GeodesicTrajectory(Vector tR, Vector tTheta, Vector r, Vector theta, Vector phiR, Vector phiTheta): tR(tR), tTheta(tTheta), r(r), theta(theta), phiR(phiR), phiTheta(phiTheta) {}
double GeodesicTrajectory::getTimeAccumulationRadial(int pos){ return tR[pos]; }
double GeodesicTrajectory::getTimeAccumulationPolar(int pos){ return tTheta[pos]; }
double GeodesicTrajectory::getTimeAccumulation(int j, int pos){
	switch(j){
		case 1:
			return tR[pos];

		case 2:
			return tTheta[pos];

		default:
			std::cout << "GEO: Index " << j << "out of range. Please choose 1 or 2. \n";
			return 0.;
	}
}
double GeodesicTrajectory::getRadialPosition(int pos){ return r[pos]; }
double GeodesicTrajectory::getPolarPosition(int pos){ return theta[pos]; }
double GeodesicTrajectory::getAzimuthalAccumulationRadial(int pos){ return phiR[pos]; }
double GeodesicTrajectory::getAzimuthalAccumulationPolar(int pos){ return phiTheta[pos]; }
double GeodesicTrajectory::getAzimuthalAccumulation(int j, int pos){
	switch(j){
		case 1:
			return phiR[pos];

		case 2:
			return phiTheta[pos];

		default:
			std::cout << "GEO: Index " << j << "out of range. Please choose 1 or 2. \n";
			return 0.;
	}
}

GeodesicConstants::GeodesicConstants(double a, double p, double e, double x, double En, double Lz, double Q): a(a), p(p), e(e), x(x), En(En), Lz(Lz), Q(Q) {}
GeodesicConstants::GeodesicConstants(double a, double p, double e, double x, double En, double Lz, double Q, double upsilonT, double upsilonR, double upsilonTheta, double upsilonPhi): a(a), p(p), e(e), x(x), En(En), Lz(Lz), Q(Q), upsilonT(upsilonT), upsilonR(upsilonR), upsilonTheta(upsilonTheta), upsilonPhi(upsilonPhi) {}
GeodesicConstants::GeodesicConstants(double a, double p, double e, double x, double En, double Lz, double Q, double r1, double r2, double r3, double r4, double z1, double z2, double upsilonT, double upsilonR, double upsilonTheta, double upsilonPhi, double carterR, double carterTheta, double carterPhi): a(a), p(p), e(e), x(x), En(En), Lz(Lz), Q(Q), r1(r1), r2(r2), r3(r3), r4(r4), z1(z1), z2(z2), upsilonT(upsilonT), upsilonR(upsilonR), upsilonTheta(upsilonTheta), upsilonPhi(upsilonPhi), carterR(carterR), carterTheta(carterTheta), carterPhi(carterPhi) {}


double GeodesicConstants::getTimeFrequency(int m, int k, int n){
	return (m*upsilonPhi + k*upsilonTheta + n*upsilonR)/upsilonT;
}

GeodesicSource::GeodesicSource(double a, double p, double e, double x, int Nsample){
	double En = 0., Lz = 0., Qc = 0.;
	double r1 = 0., r2 = 0., r3 = 0., r4 = 0.;
	double z1 = 0., z2 = 0.;
	double upT = 0., upR = 0., upPh = 0., upTh = 0.;
	double cR = 0., cTh = 0., cPh = 0.;
	if(e == 0. && std::abs(x) == 1.){
		En = kerr_geo_energy_circ(a*x, p);
		Lz = kerr_geo_momentum_circ(a*x, p);
		upT = kerr_geo_time_frequency_circ(a*x, p);
		upPh = x*kerr_geo_azimuthal_frequency_circ(a*x, p);
	}else{
		// auto start = std::chrono::system_clock::now();
		kerr_geo_orbital_constants(En, Lz, Qc, a, p, e, x);
		// auto end = std::chrono::system_clock::now();
		// std::chrono::duration<double> elapsed_seconds = end-start;
		// std::cout << "Orbital constants have been calculated \n";
		// std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

		// start = std::chrono::system_clock::now();
		kerr_geo_radial_roots(r1, r2, r3, r4, a, p, e, En, Lz, Qc);
		// end = std::chrono::system_clock::now();
		// elapsed_seconds = end-start;
		// std::cout << "Radial roots have been calculated \n";
		// std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

		// start = std::chrono::system_clock::now();
		kerr_geo_polar_roots(z1, z2, a, x, En, Lz, Qc);
		// end = std::chrono::system_clock::now();
		// elapsed_seconds = end-start;
		// std::cout << "Polar roots have been calculated \n";
		// std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

		// start = std::chrono::system_clock::now();
		kerr_geo_mino_frequencies(upT, upR, upTh, upPh, a, p, e, x, En, Lz, Qc, r1, r2, r3, r4, z1, z2);
		// end = std::chrono::system_clock::now();
		// elapsed_seconds = end-start;
		// std::cout << "Frequencies have been calculated \n";
		// std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

		// start = std::chrono::system_clock::now();
		kerr_geo_carter_frequencies(cR, cTh, cPh, upT, upR, upTh, upPh, a, En, Lz, Qc, z1, z2);
		// end = std::chrono::system_clock::now();
		// elapsed_seconds = end-start;
		// std::cout << "Carter frequencies have been calculated \n";
		// std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
	}
	if( std::abs(x*x + z2 - 1.) > 1.e-10 ){
		z2 = sqrt(1. - x*x);
	}

	Vector fourier_radial, fourier_psi, fourier_tr, fourier_phir;
	if(e != 0.){
		// auto start = std::chrono::system_clock::now();
		fourier_radial = mino_of_psi_fourier(a, p, e, En, r3, r4);
		// auto end = std::chrono::system_clock::now();
		// std::chrono::duration<double> elapsed_seconds = end-start;
		// std::cout << "Mino of psi has been calculated \n";
		// std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

		// start = std::chrono::system_clock::now();
		fourier_psi = kepler_phase_of_angle_fourier(fourier_radial);
		// end = std::chrono::system_clock::now();
		// elapsed_seconds = end-start;
		// std::cout << "psi of Mino has been calculated \n";
		// std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

		// start = std::chrono::system_clock::now();
		fourier_tr = tp_radial_of_angle_fourier(a, p, e, En, Lz, Qc, upR, fourier_psi);
		// end = std::chrono::system_clock::now();
		// elapsed_seconds = end-start;
		// std::cout << "tr of Mino has been calculated \n";
		// std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

		// start = std::chrono::system_clock::now();
		fourier_phir = phip_radial_of_angle_fourier(a, p, e, En, Lz, Qc, upR, fourier_psi);
		// end = std::chrono::system_clock::now();
		// elapsed_seconds = end-start;
		// std::cout << "phir of Mino has been calculated \n";
		// std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
	}

	Vector fourier_polar, fourier_chi, fourier_tz, fourier_phiz;
	if(std::abs(x) != 1.){
		if(std::abs(a) > 0){
			fourier_polar = mino_of_chi_fourier(a, En, z1, z2);
		}else{
			fourier_polar = mino_of_chi_schw_fourier(Lz, Qc);
		}
		fourier_chi = kepler_phase_of_angle_fourier(fourier_polar);
		fourier_tz = tp_polar_of_angle_fourier(a, x, En, Lz, Qc, upTh, fourier_chi);
		fourier_phiz = phip_polar_of_angle_fourier(a, x, En, Lz, Qc, upTh, fourier_chi);
	}
	if(fourier_psi.size() == 0){
		fourier_psi = Vector(3, 0.);
	}
	if(fourier_chi.size() == 0){
		fourier_chi = Vector(3, 0.);
	}
	if(fourier_tr.size() == 0){
		fourier_tr = Vector(3, 0.);
	}
	if(fourier_tz.size() == 0){
		fourier_tz = Vector(3, 0.);
	}
	if(fourier_phir.size() == 0){
		fourier_phir = Vector(3, 0.);
	}
	if(fourier_phiz.size() == 0){
		fourier_phiz = Vector(3, 0.);
	}

	int halfSample = Nsample/2 + 1;
	Vector tR(halfSample);
	Vector tTh(halfSample);
	Vector phiR(halfSample);
	Vector phiTh(halfSample);
	Vector rp(halfSample);
	Vector thetap(halfSample);

	kerr_trajectory(tR, tTh, rp, thetap, phiR, phiTh, p, e, x, fourier_tr, fourier_tz, fourier_psi, fourier_chi, fourier_phir, fourier_phiz);

	_geoConstants = GeodesicConstants(a, p, e, x, En, Lz, Qc, r1, r2, r3, r4, z1, z2, upT, upR, upTh, upPh, cR, cTh, cPh);
	_geoTrajectory = GeodesicTrajectory(tR, tTh, rp, thetap, phiR, phiTh);
	_geoCoefficients = GeodesicTrajectory(fourier_tr, fourier_tz, fourier_psi, fourier_chi, fourier_phir, fourier_phiz);
}
int GeodesicSource::getOrbitalSampleNumber(){ return 2*(_geoTrajectory.tR.size() - 1); }
double GeodesicSource::getBlackHoleSpin(){ return _geoConstants.a; }
double GeodesicSource::getSemiLatusRectum(){ return _geoConstants.p; }
double GeodesicSource::getEccentricity(){ return _geoConstants.e; }
double GeodesicSource::getInclination(){ return _geoConstants.x; }
double GeodesicSource::getOrbitalEnergy(){ return _geoConstants.En; }
double GeodesicSource::getOrbitalAngularMomentum(){ return _geoConstants.Lz; }
double GeodesicSource::getCarterConstant(){ return _geoConstants.Q; }
double GeodesicSource::getRadialRoot(int i){
	switch(i){
		case 1:
			return _geoConstants.r1;

		case 2:
			return _geoConstants.r2;

		case 3:
			return _geoConstants.r3;

		case 4:
			return _geoConstants.r4;

		default:
			std::cout << "GEO: Index " << i << " out of range. There are only four radial roots.\n";
			return 0.;
	}
}
double GeodesicSource::getPolarRoot(int i){
	switch(i){
		case 1:
			return _geoConstants.z1;

		case 2:
			return _geoConstants.z2;

		default:
			std::cout << "GEO: Index " << i << " out of range. There are only two unique polar roots.\n";
			return 0.;
	}
}
double GeodesicSource::getMinoFrequency(int mu){
	switch(mu){
		case 0:
			return _geoConstants.upsilonT;

		case 1:
			return _geoConstants.upsilonR;

		case 2:
			return _geoConstants.upsilonTheta;

		case 3:
			return _geoConstants.upsilonPhi;

		default:
			std::cout << "GEO: Frequency index " << mu << "out of range. \n";
			return 0.;
	}
}
double GeodesicSource::getTimeFrequency(int i){
	switch(i){
		case 1:
			return _geoConstants.upsilonR/_geoConstants.upsilonT;

		case 2:
			return _geoConstants.upsilonTheta/_geoConstants.upsilonT;

		case 3:
			return _geoConstants.upsilonPhi/_geoConstants.upsilonT;

		default:
			std::cout << "GEO: Frequency index " << i << "out of range. \n";
			return 0.;
	}
}
double GeodesicSource::getCarterFrequency(int i){
	switch(i){
		case 1:
			return _geoConstants.carterR;

		case 2:
			return _geoConstants.carterTheta;

		case 3:
			return _geoConstants.carterPhi;

		default:
			std::cout << "GEO: Frequency index " << i << "out of range. \n";
			return 0.;
	}
}
double GeodesicSource::getTimeFrequency(int m, int k, int n){
	return (m*_geoConstants.upsilonPhi + k*_geoConstants.upsilonTheta + n*_geoConstants.upsilonR)/_geoConstants.upsilonT;
}
double GeodesicSource::getCarterFrequency(int m, int k, int n){
	return (m*_geoConstants.carterPhi + k*_geoConstants.carterTheta + n*_geoConstants.carterR);
}

Vector GeodesicSource::getTimeAccumulation(int j){
	switch(j){
		case 1:
			return _geoTrajectory.tR;

		case 2:
			return _geoTrajectory.tTheta;

		default:
			std::cout << "GEO: Index " << j << "out of range. Please choose 1 or 2. \n";
			return Vector();
	}
}
Vector GeodesicSource::getRadialPosition(){ return _geoTrajectory.r; }
Vector GeodesicSource::getPolarPosition(){ return _geoTrajectory.theta; }
Vector GeodesicSource::getAzimuthalAccumulation(int j){
	switch(j){
		case 1:
			return _geoTrajectory.phiR;

		case 2:
			return _geoTrajectory.phiTheta;

		default:
			std::cout << "GEO: Index " << j << "out of range. Please choose 1 or 2. \n";
			return Vector();
	}
}
double GeodesicSource::getTimeAccumulation(int j, int pos){
	switch(j){
		case 1:
			return _geoTrajectory.tR[pos];

		case 2:
			return _geoTrajectory.tTheta[pos];

		default:
			std::cout << "GEO: Index " << j << "out of range. Please choose 1 or 2. \n";
			return 0.;
	}
}
double GeodesicSource::getRadialPosition(int pos){ return _geoTrajectory.r[pos]; }
double GeodesicSource::getPolarPosition(int pos){ return _geoTrajectory.theta[pos]; }
double GeodesicSource::getAzimuthalAccumulation(int j, int pos){
	switch(j){
		case 1:
			return _geoTrajectory.phiR[pos];

		case 2:
			return _geoTrajectory.phiTheta[pos];

		default:
			std::cout << "GEO: Index " << j << "out of range. Please choose 1 or 2. \n";
			return 0.;
	}
}

double GeodesicSource::getPsiRadialOfMinoTime(double lambda){
	return kepler_phase_of_angle(lambda*_geoConstants.upsilonR, _geoCoefficients.r);
}
double GeodesicSource::getPsiPolarOfMinoTime(double lambda){
	return kepler_phase_of_angle(lambda*_geoConstants.upsilonTheta, _geoCoefficients.theta);
}
double GeodesicSource::getTimePositionOfMinoTime(double lambda){
	return lambda*_geoConstants.upsilonT + phip_of_angle(lambda*_geoConstants.upsilonR, _geoCoefficients.tR) + phip_of_angle(lambda*_geoConstants.upsilonTheta, _geoCoefficients.tTheta);
}
double GeodesicSource::getRadialPositionOfMinoTime(double lambda){
	return rp_of_angle(lambda*_geoConstants.upsilonR, _geoConstants.p, _geoConstants.e, _geoCoefficients.r);
}
double GeodesicSource::getPolarPositionOfMinoTime(double lambda){
	return acos(zp_of_angle(lambda*_geoConstants.upsilonTheta, _geoConstants.x, _geoCoefficients.theta));
}
double GeodesicSource::getAzimuthalPositionOfMinoTime(double lambda){
	return lambda*_geoConstants.upsilonPhi + phip_of_angle(lambda*_geoConstants.upsilonR, _geoCoefficients.phiR) + phip_of_angle(lambda*_geoConstants.upsilonTheta, _geoCoefficients.phiTheta);
}
Vector GeodesicSource::getPositionOfMinoTime(double lambda){
	Vector xp(4);
	xp[0] = getTimePositionOfMinoTime(lambda);
	xp[1] = getRadialPositionOfMinoTime(lambda);
	xp[2] = getPolarPositionOfMinoTime(lambda);
	xp[3] = getAzimuthalPositionOfMinoTime(lambda);
	return xp;
}

Vector GeodesicSource::getPsiRadialOfMinoTime(Vector lambda){
	Vector xp(lambda.size());
	for(size_t i = 0; i < xp.size(); i++){
		xp[i] = getPsiRadialOfMinoTime(lambda[i]);
	}
	return xp;
}

Vector GeodesicSource::getPsiPolarOfMinoTime(Vector lambda){
	Vector xp(lambda.size());
	for(size_t i = 0; i < xp.size(); i++){
		xp[i] = getPsiPolarOfMinoTime(lambda[i]);
	}
	return xp;
}

Vector GeodesicSource::getTimePositionOfMinoTime(Vector lambda){
	Vector xp(lambda.size());
	for(size_t i = 0; i < xp.size(); i++){
		xp[i] = getTimePositionOfMinoTime(lambda[i]);
	}
	return xp;
}
Vector GeodesicSource::getRadialPositionOfMinoTime(Vector lambda){
	Vector xp(lambda.size());
	for(size_t i = 0; i < xp.size(); i++){
		xp[i] = getRadialPositionOfMinoTime(lambda[i]);
	}
	return xp;
}
Vector GeodesicSource::getPolarPositionOfMinoTime(Vector lambda){
	Vector xp(lambda.size());
	for(size_t i = 0; i < xp.size(); i++){
		xp[i] = getPolarPositionOfMinoTime(lambda[i]);
	}
	return xp;
}
Vector GeodesicSource::getAzimuthalPositionOfMinoTime(Vector lambda){
	Vector xp(lambda.size());
	for(size_t i = 0; i < xp.size(); i++){
		xp[i] = getAzimuthalPositionOfMinoTime(lambda[i]);
	}
	return xp;
}
struct tp_params{
	double t0;
	double upsilonT;
	double upsilonR;
	double upsilonTh;
	Vector fourierR;
	Vector fourierTh;
};

double tp_root(double x, void *params){
	struct tp_params *p = (struct tp_params *) params;
	double t0 = p->t0;
	double upsilonT = p->upsilonT;
	double upsilonR = p->upsilonR;
	double upsilonTh = p->upsilonTh;
	Vector fourierR = p->fourierR;
	Vector fourierTh = p->fourierTh;

	return upsilonT*x + tp_of_angle(upsilonR*x, fourierR) + tp_of_angle(upsilonTh*x, fourierTh) - t0;
}
double GeodesicSource::getMinoTimeOfTime(double t){
	if(t == 0.){
		return 0.;
	}
	double upsilonT = getMinoFrequency(0);
	double upsilonR = getMinoFrequency(1);
	double upsilonTh = getMinoFrequency(2);
	// double lambda_guess = t/upsilonT; recall that t \approx \Upsilon_t \times \lambda
	if(getEccentricity() == 0. && std::abs(getInclination()) == 1.){
		return t/upsilonT;
	}
	// double deltaT = (std::abs(getTimeAccumulation(1,1)) + std::abs(getTimeAccumulation(1,2)))/double(getTimeAccumulation(1).size()) + (std::abs(getTimeAccumulation(2,1)) + std::abs(getTimeAccumulation(2,2)))/double(getTimeAccumulation(2).size());
	// double lambda_lo = (t - 2.*deltaT)/upsilonT;
	// double lambda_hi = (t + 2.*deltaT)/upsilonT;
	Vector fourierR = getTimeCoefficients(1);
	Vector fourierTh = getTimeCoefficients(2);

	double fourierR_sum = 0.;
	for(size_t i = 0; i < fourierR.size(); i++){
		fourierR_sum += std::abs(fourierR[i]);
	}
	double fourierTh_sum = 0.;
	for(size_t i = 0; i < fourierTh.size(); i++){
		fourierTh_sum += std::abs(fourierTh[i]);
	}

	double lambda_lo = (t - fourierR_sum - fourierTh_sum)/upsilonT;
	double lambda_hi = (t + fourierR_sum + fourierTh_sum)/upsilonT;

	int status, out;
	int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;

	gsl_function F;
	struct tp_params params = {
		.t0 = t, // t0 is the target time
		.upsilonT = upsilonT,
		.upsilonR = upsilonR,
		.upsilonTh = upsilonTh,
		.fourierR = fourierR,
		.fourierTh = fourierTh
	};
	F.function = &tp_root;
	F.params = &params;

	double f_lo = tp_root(lambda_lo, &params);
	double f_hi = tp_root(lambda_hi, &params);

	if(f_lo*f_hi > 0.){
		std::cerr << "GEO: The function does not change sign in the interval [" << lambda_lo << ", " << lambda_hi << "]. \n";
		return 0.;
	}

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
    out = gsl_root_fsolver_set(s, &F, lambda_lo, lambda_hi);
	double lambda = 0.;

	do{
		iter++;
		status = gsl_root_fsolver_iterate (s);
		lambda = gsl_root_fsolver_root (s);
		lambda_lo = gsl_root_fsolver_x_lower (s);
		lambda_hi = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (lambda_lo, lambda_hi,
										0, 1.e-14);
	}while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free (s);

	return lambda;
}

Vector GeodesicSource::getTimeCoefficients(int j){
	switch(j){
		case 1:
			return _geoCoefficients.tR;

		case 2:
			return _geoCoefficients.tTheta;

		default:
			std::cout << "GEO: Index " << j << "out of range. Please choose 1 or 2. \n";
			return Vector();
	}
}
Vector GeodesicSource::getRadialCoefficients(){ return _geoCoefficients.r; }
Vector GeodesicSource::getPolarCoefficients(){ return _geoCoefficients.theta; }
Vector GeodesicSource::getAzimuthalCoefficients(int j){
	switch(j){
		case 1:
			return _geoCoefficients.phiR;

		case 2:
			return _geoCoefficients.phiTheta;

		default:
			std::cout << "GEO: Index " << j << "out of range. Please choose 1 or 2. \n";
			return Vector();
	}
}

GeodesicConstants GeodesicSource::getConstants(){  return _geoConstants; }
GeodesicTrajectory GeodesicSource::getTrajectory(){ return _geoTrajectory; }
GeodesicTrajectory GeodesicSource::getCoefficients(){ return _geoCoefficients; }

GeodesicConstants& GeodesicSource::getConstantsRef(){  return _geoConstants; }
GeodesicTrajectory& GeodesicSource::getTrajectoryRef(){ return _geoTrajectory; }
GeodesicTrajectory& GeodesicSource::getCoefficientsRef(){ return _geoCoefficients; }

void GeodesicSource::setConstants(double a, double p, double e, double x, double En, double Lz, double Qc, double r1, double r2, double r3, double r4, double z1, double z2, double upT, double upR, double upTh, double upPhi, double cR, double cTh, double cPhi){
	_geoConstants.a = a;
	_geoConstants.p = p;
	_geoConstants.e = e;
	_geoConstants.x = x;
	_geoConstants.En = En;
	_geoConstants.Lz = Lz;
	_geoConstants.Q = Qc;
	_geoConstants.r1 = r1;
	_geoConstants.r2 = r2;
	_geoConstants.r3 = r3;
	_geoConstants.r4 = r4;
	_geoConstants.z1 = z1;
	_geoConstants.z2 = z2;
	_geoConstants.upsilonT = upT;
	_geoConstants.upsilonR = upR;
	_geoConstants.upsilonTheta = upTh;
	_geoConstants.upsilonPhi = upPhi;
	_geoConstants.carterR = cR;
	_geoConstants.carterTheta = cTh;
	_geoConstants.carterPhi = cPhi;
}
void GeodesicSource::setTrajectory(Vector tR, Vector tTheta, Vector r, Vector theta, Vector phiR, Vector phiTheta){
	_geoTrajectory.tR = tR;
	_geoTrajectory.tTheta = tTheta;
	_geoTrajectory.r = r;
	_geoTrajectory.theta = theta;
	_geoTrajectory.phiR = phiR;
	_geoTrajectory.phiTheta = phiTheta;
}
void GeodesicSource::setCoefficients(Vector tR, Vector tTheta, Vector r, Vector theta, Vector phiR, Vector phiTheta){
	_geoCoefficients.tR = tR;
	_geoCoefficients.tTheta = tTheta;
	_geoCoefficients.r = r;
	_geoCoefficients.theta = theta;
	_geoCoefficients.phiR = phiR;
	_geoCoefficients.phiTheta = phiTheta;
}
// void GeodesicSource::save(std::string dir){
// 	int Nsample = getOrbitalSampleNumber();
// 	double a = getBlackHoleSpin();
// 	double p = getSemiLatusRectum();
// 	double e = getEccentricity();
// 	double x = getInclination();
// 	double En = getOrbitalEnergy();
// 	double Lz = getOrbitalAngularMomentum();
// 	double Qc = getCarterConstant();
// 	double r1 = getRadialRoot(1);
// 	double r2 = getRadialRoot(2);
// 	double r3 = getRadialRoot(3);
// 	double r4 = getRadialRoot(4);
// 	double z1 = getPolarRoot(1);
// 	double z2 = getPolarRoot(2);
// 	double UpT = getMinoFrequency(0);
// 	double UpR = getMinoFrequency(1);
// 	double UpTh = getMinoFrequency(2);
// 	double UpPh = getMinoFrequency(3);

// 	Vector deltaTR = getTimeAccumulation(1);
// 	Vector deltaTTh = getTimeAccumulation(2);
// 	Vector rp = getRadialPosition();
// 	Vector thp = getPolarPosition();
// 	Vector deltaPhR = getAzimuthalAccumulation(1);
// 	Vector deltaPhTh = getAzimuthalAccumulation(2);

// 	GeodesicTrajectory fourierCoefficients = getCoefficients();

// 	// if(!boost::filesystem::exists(dir)){
// 	// 	boost::filesystem::create_directory(dir);
// 	// }
// 	std::string subdir;
// 	if(dir.back() == '/'){
// 		subdir = dir + "orbit/";
// 	}else{
// 		subdir = dir + "/orbit/";
// 	}

// 	// if(!boost::filesystem::exists(subdir)){
// 	// 	boost::filesystem::create_directory(subdir);
// 	// }
// 	char buff[500];
// 	sprintf(buff, "geo_a%.5f_p%.5f_e%.5f_x%.5f.txt", a, p, e, x);

// 	std::string filepath = subdir + buff;
// 	std::cout << "Saving file to " << filepath << "\n";
// 	std::ofstream file;
// 	file.open(filepath);

// 	file << "a\tp\te\tx\tNsample\n";
// 	sprintf(buff, "%.15f\t%.15f\t%.15f\t%.15f\t%d\n\n", a, p, e, x, Nsample);
// 	file << buff;

// 	file << "En\tLz\tQc\n";
// 	sprintf(buff, "%.15f\t%.15f\t%.15f\n\n", En, Lz, Qc);
// 	file << buff;

// 	file << "r1\tr2\tr3\tr4\tz1\tz2\n";
// 	sprintf(buff, "%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n\n", r1, r2, r3, r4, z1, z2);
// 	file << buff;

// 	file << "Upsilon_t\tUpsilon_r\tUpsilon_theta\tUpsilon_phi\n";
// 	sprintf(buff, "%.15f\t%.15f\t%.15f\t%.15f\n\n", UpT, UpR, UpTh, UpPh);
// 	file << buff;

// 	// file << "t_r\tt_theta\tr_p\ttheta_p\tphi_r\tphi_theta\n";
// 	// for(size_t i = 0; i < rp.size(); i++){
// 	// 	sprintf(buff, "%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n", deltaTR[i], deltaTTh[i], rp[i], thp[i], deltaPhR[i], deltaPhTh[i]);
// 	// 	file << buff;
// 	// }
// 	// file << "\n";

// 	file << "f_tr\tf_ttheta\tf_r\tf_theta\tf_phir\tf_phitheta\n";
// 	for(size_t i = 0; i < rp.size(); i++){
// 		sprintf(buff, "%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n", fourierCoefficients.tR[i], fourierCoefficients.tTheta[i], fourierCoefficients.r[i], fourierCoefficients.theta[i], fourierCoefficients.phiR[i], fourierCoefficients.phiTheta[i]);
// 		file << buff;
// 	}

// 	file.close();
// }
////////////////////////
//   Generic orbits   //
////////////////////////

// GeodesicSource kerr_geo(double a, double p, double e, double x, int sampleN = pow(2, 10)){
// 	double En, Lz, Qc;
// 	kerr_geo_orbital_constants(En, Lz, Qc, a, p, e, x);
// 	double r1, r2, r3, r4;
// 	kerr_geo_radial_roots(r1, r2, r3, r4, a, p, e, En, Lz, Qc);
// 	//std::cout << "r1 = "<<r1<<", r2 = "<<r2<<", r3 = "<<r3<<", r4 = "<<r4<<" \n";
// 	double z1, z2;
// 	kerr_geo_polar_roots(z1, z2, a, En, Lz, Qc);
// 	//std::cout << "z1 = "<<z1<<", z2 = "<<z2<<" \n";
// 	double upT, upR, upTh, upPh;
// 	kerr_geo_mino_frequencies(upT, upR, upTh, upPh, a, p, e, x, En, Lz, Qc, r1, r2, r3, r4, z1, z2);
// 	double cR, cTh, cPh;
// 	kerr_geo_carter_frequencies(cR, cTh, cPh, upT, upR, upTh, upPh, a, p, e, x, En, Lz, Qc, r1, r2, r3, r4, z1, z2);
//
// 	Vector tpR(sampleN + 1);
// 	Vector tpTh(sampleN + 1);
//
// 	Vector phipR(sampleN + 1);
// 	Vector phipTh(sampleN + 1);
//
// 	Vector rp(sampleN + 1);
// 	rp[0] = p;
// 	Vector thetap(sampleN + 1);
// 	thetap[0] = M_PI/2.;
//
// 	GeodesicSource geo;
// 	geo.setConstants(a, p, e, x, En, Lz, Qc, upT, upR, upTh, upPh);
// 	geo.setTrajectory(tpR, tpTh, rp, thetap, phipR, phipTh);
//
// 	return geo;
// }

GeodesicSource load_geodesic(std::string file){
	int Nsample;
	double a, p, e, x;
	double En, Lz, Qc;
	double r1, r2, r3, r4, z1, z2;
	double upT, upR, upTh, upPh;
	double dtr, dtth, r, th, dpr, dpth;
	Vector ftR, ftTh, frp, fthetap, fphiR, fphiTh;

	std::istringstream lin;
	std::ifstream inFile(file);
	std::string line;
	std::getline(inFile, line);
	lin.clear();
	lin.str(line);
	while(line.front() == '#' || line.empty() || isalpha(line.front())){
		std::getline(inFile, line);
		lin.clear();
		lin.str(line);
	}
	lin >> a >> p >> e >> x >> Nsample;

	std::getline(inFile, line);
	lin.clear();
	lin.str(line);
	while(line.front() == '#' || line.empty() || isalpha(line.front())){
		std::getline(inFile, line);
		lin.clear();
		lin.str(line);
	}
	lin >> En >> Lz >> Qc;

	std::getline(inFile, line);
	lin.clear();
	lin.str(line);
	while(line.front() == '#' || line.empty() || isalpha(line.front())){
		std::getline(inFile, line);
		lin.clear();
		lin.str(line);
	}
	lin >> r1 >> r2 >> r3 >> r4 >> z1 >> z2;

	std::getline(inFile, line);
	lin.clear();
	lin.str(line);
	while(line.front() == '#' || line.empty() || isalpha(line.front())){
		std::getline(inFile, line);
		lin.clear();
		lin.str(line);
	}
	lin >> upT >> upR >> upTh >> upPh;

	std::getline(inFile, line);
	lin.clear();
	lin.str(line);
	while(line.front() == '#' || line.empty() || isalpha(line.front())){
		std::getline(inFile, line);
		lin.clear();
		lin.str(line);
	}
	lin >> dtr >> dtth >> r >> th >> dpr >> dpth;
	ftR.push_back(dtr);
	ftTh.push_back(dtth);
	frp.push_back(r);
	fthetap.push_back(th);
	fphiR.push_back(dpr);
	fphiTh.push_back(dpth);
	for (std::string line; std::getline(inFile, line); ) {
	    lin.clear();
	    lin.str(line);
	    if(lin >> dtr >> dtth >> r >> th >> dpr >> dpth){
				ftR.push_back(dtr);
				ftTh.push_back(dtth);
				frp.push_back(r);
				fthetap.push_back(th);
				fphiR.push_back(dpr);
				fphiTh.push_back(dpth);
	    }
	}

	inFile.close();

	double cR = 0., cTh = 0., cPh = 0.;
	if(e > 0. || std::abs(x) < 1.){
		kerr_geo_carter_frequencies(cR, cTh, cPh, upT, upR, upTh, upPh, a, En, Lz, Qc, z1, z2);
	}

	GeodesicSource geo;
	geo.setConstants(a, p, e, x, En, Lz, Qc, r1, r2, r3, r4, z1, z2, upT, upR, upTh, upPh, cR, cTh, cPh);
	geo.setCoefficients(ftR, ftTh, frp, fthetap, fphiR, fphiTh);

	int halfSample = Nsample/2 + 1;
	Vector tR(halfSample);
	Vector tTh(halfSample);
	Vector phiR(halfSample);
	Vector phiTh(halfSample);
	Vector rp(halfSample);
	Vector thetap(halfSample);

	kerr_trajectory(tR, tTh, rp, thetap, phiR, phiTh, p, e, x, ftR, ftTh, frp, fthetap, fphiR, fphiTh);

	geo.setTrajectory(tR, tTh, rp, thetap, phiR, phiTh);

	return geo;
}

///////////////////////
// Orbital constants //
///////////////////////

double kerr_geo_energy_circ(double a, double p, double x){
	return pow(pow(pow(p,2) - 1.*pow(a,2)*(-1. + pow(x,2)),-1)*((-3. + p)*pow(-2. + p,2)*pow(p,5) + pow(a,4)*pow(p,2)*(4. - 5.*p*(-1. + pow(x,2)) + 3.*pow(p,2)*(-1. + pow(x,2)))*(-1. + pow(x,2)) + pow(a,2)*pow(p,3)*(4. + p*(12. - 7.*pow(x,2)) - 3.*pow(p,3)*(-1. + pow(x,2)) - 4.*pow(x,2) + pow(p,2)*(-13. + 10.*pow(x,2))) - 2.*x*pow(a,5)*(-1. + pow(x,2))*pow(pow(p,3) + p*pow(a,2)*(-1. + pow(x,2)),0.5) + a*(-2.*x*pow(p,4.5)*pow(pow(p,2) + pow(a,2)*(-1. + pow(x,2)),0.5) + 4.*x*pow(p,3)*pow(pow(p,3) + p*pow(a,2)*(-1. + pow(x,2)),0.5)) + 2.*pow(a,3)*(2.*p*x*(-1. + pow(x,2))*pow(pow(p,3) + p*pow(a,2)*(-1. + pow(x,2)),0.5) - 1.*pow(x,3)*pow(pow(p,7) + pow(a,2)*pow(p,5)*(-1. + pow(x,2)),0.5)) - 1.*pow(a,6)*(pow(p,2)*(-1. + pow(x,2)) + pow(x,2) - 1.*p*(1. + 2.*pow(x,2)))*pow(-1. + pow(x,2),2))*pow(pow(-3. + p,2)*pow(p,4) - 2.*pow(a,2)*pow(p,2)*(3. + 2.*p + pow(p,2)*(-1. + pow(x,2)) - 3.*pow(x,2)) + pow(a,4)*(-1. + pow(x,2))*(-1. + pow(p,2)*(-1. + pow(x,2)) + pow(x,2) -  2.*p*(1. + pow(x,2))),-1),0.5);
}

double kerr_geo_energy_circ_2(double a, double p, double x){
	double rmax = p;
	double deltaMax = rmax*rmax - 2.*rmax + a*a;
	double ddeltaMax = 2.*rmax - 2.;

	double fmax = pow(rmax, 4) + a*a*(rmax*(rmax + 2.) + (1. - x*x)*deltaMax);
	double fmin = 4.*pow(rmax, 3) + a*a*(rmax + (rmax + 2.) + (1. - x*x)*ddeltaMax);
	double gmax = 2.*a*rmax;
	double gmin = 2.*a;
	double hmax = rmax*(rmax - 2.) + (1. - x*x)/(x*x)*deltaMax;
	double hmin = rmax + (rmax - 2.) + (1. - x*x)/(x*x)*ddeltaMax;
	double dmax = (rmax*rmax + a*a*(1. - x*x))*deltaMax;
	double dmin = 2.*rmax*deltaMax + (rmax*rmax + a*a*(1. - x*x))*ddeltaMax;

	double kap = dmax*hmin - hmax*dmin;
	double eps = dmax*gmin - gmax*dmin;
	double rho = fmax*hmin - hmax*fmin;
	double eta = fmax*gmin - gmax*fmin;
	double sig = gmax*hmin - hmax*gmin;

	return sqrt((kap*rho + 2.*eps*sig - 2.*x*sqrt(sig*(sig*eps*eps + rho*eps*kap - eta*kap*kap)/(x*x)))/(rho*rho + 4.*eta*sig));
}

double kerr_geo_energy(double a, double p, double e, double x){
	double rmax = p/(1. - e);
	double rmin = p/(1. + e);
	double deltaMax = rmax*rmax - 2.*rmax + a*a;
	double deltaMin = rmin*rmin - 2.*rmin + a*a;

	if(std::abs(e) < 1.e-14){
		return kerr_geo_energy_circ(a, p, x);
	}

	double fmax = pow(rmax, 4) + a*a*(rmax*(rmax + 2.) + (1. - x*x)*deltaMax);
	double fmin = pow(rmin, 4) + a*a*(rmin*(rmin + 2.) + (1. - x*x)*deltaMin);
	double gmax = 2.*a*rmax;
	double gmin = 2.*a*rmin;
	double hmax = rmax*(rmax - 2.) + (1. - x*x)/(x*x)*deltaMax;
	double hmin = rmin*(rmin - 2.) + (1. - x*x)/(x*x)*deltaMin;
	double dmax = (rmax*rmax + a*a*(1. - x*x))*deltaMax;
	double dmin = (rmin*rmin + a*a*(1. - x*x))*deltaMin;

	double kap = dmax*hmin - hmax*dmin;
	double eps = dmax*gmin - gmax*dmin;
	double rho = fmax*hmin - hmax*fmin;
	double eta = fmax*gmin - gmax*fmin;
	double sig = gmax*hmin - hmax*gmin;

	return sqrt((kap*rho + 2.*eps*sig - 2.*x*sqrt(sig*(sig*eps*eps + rho*eps*kap - eta*kap*kap)/(x*x)))/(rho*rho + 4.*eta*sig));
}

double kerr_geo_momentum(double a, double p, double e, double x){
	return kerr_geo_momentum(kerr_geo_energy(a, p, e, x), a, p, e, x);
}
double kerr_geo_momentum(double En, double a, double p, double e, double x){
	double rmax = p/(1. - e);
	double deltaMax = rmax*rmax - 2.*rmax + a*a;

	double fmax = pow(rmax, 4) + a*a*(rmax*(rmax + 2.) + (1. - x*x)*deltaMax);
	double gmax = 2.*a*rmax;
	double hmax = rmax*(rmax - 2.) + (1. - x*x)/(x*x)*deltaMax;
	double dmax = (rmax*rmax + a*a*(1. - x*x))*deltaMax;

	return (-En*gmax + x*sqrt((En*En*(gmax*gmax + fmax*hmax) - dmax*hmax)/(x*x)))/hmax;
}

double kerr_geo_carter(double a, double p, double e, double x){
	double En = kerr_geo_energy(a, p, e, x);
	return kerr_geo_carter(En, kerr_geo_momentum(En, a, p, e, x), a, p, e, x);
}
double kerr_geo_carter(double En, double Lz, double a, double, double, double x){
	double zmin = sqrt(1. - x*x);

	return zmin*zmin*(a*a*(1. - En*En) + Lz*Lz/(x*x));
}

void kerr_geo_orbital_constants(double &En, double &Lz, double &Qc, double a, double p, double e, double x){
	double rmax = p/(1. - e);
	double rmin = p/(1. + e);
	double deltaMax = rmax*rmax - 2.*rmax + a*a;
	double deltaMin = rmin*rmin - 2.*rmin + a*a;
	double xsq = x*x;

	if(std::abs(e) > 1.e-14){
		double fmax = pow(rmax, 4) + a*a*(rmax*(rmax + 2.) + (1. - xsq)*deltaMax);
		double fmin = pow(rmin, 4) + a*a*(rmin*(rmin + 2.) + (1. - xsq)*deltaMin);
		double gmax = 2.*a*rmax;
		double gmin = 2.*a*rmin;
		double hmax = rmax*(rmax - 2.) + (1. - xsq)/xsq*deltaMax;
		double hmin = rmin*(rmin - 2.) + (1. - xsq)/xsq*deltaMin;
		double dmax = (rmax*rmax + a*a*(1. - xsq))*deltaMax;
		double dmin = (rmin*rmin + a*a*(1. - xsq))*deltaMin;

		double kap = dmax*hmin - hmax*dmin;
		double eps = dmax*gmin - gmax*dmin;
		double rho = fmax*hmin - hmax*fmin;
		double eta = fmax*gmin - gmax*fmin;
		double sig = gmax*hmin - hmax*gmin;

		En = sqrt((kap*rho + 2.*eps*sig - 2.*x*sqrt(sig*(sig*eps*eps + rho*eps*kap - eta*kap*kap)/(x*x)))/(rho*rho + 4.*eta*sig));
		Lz = (-En*gmax + x*sqrt((En*En*(gmax*gmax + fmax*hmax) - dmax*hmax)/(x*x)))/hmax;
		Qc = (1. - x*x)*(a*a*(1. - En*En) + Lz*Lz/(x*x));
	}else{
		double fmax = 4.*pow(rmax, 3) + 2.*a*a*((2. - xsq)*rmax + xsq);
		double fmin = pow(rmin, 4) + a*a*(rmin*(rmin + 2.) + (1. - xsq)*deltaMin);
		double gmax = 2.*a;
		double gmin = 2.*a*rmin;
		double hmax = 2.*(rmax - 1.)/xsq;
		double hmin = rmin*(rmin - 2.) + (1. - xsq)/xsq*deltaMin;
		double dmax = 2.*(2.*rmax - 3.)*rmax*rmax + 2.*a*a*((2. - xsq)*rmax - (1. - xsq));
		double dmin = (rmin*rmin + a*a*(1. - xsq))*deltaMin;

		double kap = dmax*hmin - hmax*dmin;
		double eps = dmax*gmin - gmax*dmin;
		double rho = fmax*hmin - hmax*fmin;
		double eta = fmax*gmin - gmax*fmin;
		double sig = gmax*hmin - hmax*gmin;

		En = sqrt((kap*rho + 2.*eps*sig - 2.*x*sqrt(sig*(sig*eps*eps + rho*eps*kap - eta*kap*kap)/xsq))/(rho*rho + 4.*eta*sig));
		Lz = (-En*gmax + x*sqrt((En*En*(gmax*gmax + fmax*hmax) - dmax*hmax)/xsq))/hmax;
		Qc = (1. - xsq)*(a*a*(1. - En*En) + Lz*Lz/xsq);
	}
}

void kerr_geo_radial_roots(double &r1, double &r2, double &r3, double &r4, double a, double p, double e, double En, double , double Qc){
	r1 = p/(1. - e);
	r2 = p/(1. + e);

	double ApB = 2./(1. - En*En) - (r2 + r1);
	double AB = a*a*Qc/(1. - En*En)/r1/r2;
	r3 = 0.5*ApB + 0.5*sqrt(ApB*ApB - 4.*AB);
	r4 = AB/r3;
}

void kerr_geo_polar_roots(double &z1, double &z2, double a, double x, double En, double Lz, double Qc){
	z2 = sqrt(1. - x*x);
	if(std::abs(a) == 0.){
		z1 = 0.;
	}else{
		double beta = a*a*(1. - En*En);
		double alpha = Lz*Lz + Qc + beta;
		double zsqrt = sqrt(alpha*alpha - 4*Qc*beta);
		z1 = sqrt((alpha + zsqrt)/(2.*beta));
	}
}

Vector rp_psi(Vector &psi, double p, double e){
	Vector rp(psi.size());
	rp_psi(rp, psi, p, e);
	return rp;
}

void rp_psi(Vector &rp, Vector &psi, double p, double e){
	for(size_t i = 0; i < psi.size(); i++){
		rp[i] = p/(1. + e*cos(psi[i]));
	}
}

double rp_psi(double psi, double p, double e){
	return  p/(1. + e*cos(psi));
}

Vector zp_chi(Vector &chi, double x){
	Vector zp(chi.size());
	zp_chi(zp, chi, x);
	return zp;
}

void zp_chi(Vector &zp, Vector &chi, double x){
	for(size_t i = 0; i < chi.size(); i++){
		zp[i] = sqrt(1. - x*x)*cos(chi[i]);
	}
}

double zp_chi(double chi, double x){
	return sqrt(1. - x*x)*cos(chi);
}

Vector kerr_delta(Vector &r, double a){
	Vector delta(r.size());
	kerr_delta(delta, r, a);
	return delta;
}
void kerr_delta(Vector &delta, Vector &r, double a){
	for(size_t i = 0; i < r.size(); i++){
		delta[i] = kerr_delta(r[i], a);
	}
}
double kerr_delta(double r, double a){
	return r*r - 2.*r + a*a;
}

Vector kerr_varpi2(Vector &r, double a){
	Vector varpi2(r.size());
	kerr_varpi2(varpi2, r, a);
	return varpi2;
}
void kerr_varpi2(Vector &varpi2, Vector &r, double a){
	for(size_t i = 0; i < r.size(); i++){
		varpi2[i] = kerr_varpi2(r[i], a);
	}
}
double kerr_varpi2(double r, double a){
	return r*r + a*a;
}

double dmino_dpsi(double psi, double , double p, double e, double En, double r3, double const&r4){
	double p3 = r3*(1. - e), p4 = r4*(1. + e);
	return (1. - e*e)/sqrt(((p - p4) + e*(p - p4*cos(psi)))*((p - p3) - e*(p + p3*cos(psi))))/sqrt(1. - En*En);
}
double dmino_dchi(double chi, double a, double En, double z1, double z2){
	return 1./(a*sqrt(1. - En*En)*sqrt(z1*z1 - z2*z2*pow(cos(chi), 2)));
}
double dmino_dchi_schw(double Lz, double Qc){
	return 1./sqrt(Qc + pow(Lz, 2));
}
double dVtrdMino(double rp, double a, double En, double Lz, double ){
	return (En*pow(kerr_varpi2(rp, a), 2) - a*Lz*kerr_varpi2(rp, a))/kerr_delta(rp, a);
}
double dVtzdMino(double zp, double a, double En, double Lz, double ){
	return a*Lz - a*a*En*(1. - zp*zp);
}
double dVphirdMino(double rp, double a, double En, double Lz, double ){
	return a*(En*kerr_varpi2(rp, a) - a*Lz)/kerr_delta(rp, a);
}
double dVphizdMino(double zp, double a, double En, double Lz, double ){
	return Lz/(1. - zp*zp) - a*En;
}

///////////////////////////////////
// Spectral integration approach //
///////////////////////////////////

void kerr_geo_mino_frequencies(double &upT, double &upR, double &upTh, double &upPh, double a, double p, double e, double x){
	double En, Lz, Qc;
	kerr_geo_orbital_constants(En, Lz, Qc, a, p, e, x);
	//std::cout << "En = "<<En<<", Lz = "<<Lz<<", Qc = "<<Qc<<" \n";
	double r1, r2, r3, r4;
	kerr_geo_radial_roots(r1, r2, r3, r4, a, p, e, En, Lz, Qc);
	//std::cout << "r1 = "<<r1<<", r2 = "<<r2<<", r3 = "<<r3<<", r4 = "<<r4<<" \n";
	double z1, z2;
	kerr_geo_polar_roots(z1, z2, a, x, En, Lz, Qc);
	if( std::abs(x*x + z2 - 1.) > 1.e-10 ){
		z2 = sqrt(1. - x*x);
	}
	//std::cout << "z1 = "<<z1<<", z2 = "<<z2<<" \n";
	kerr_geo_mino_frequencies(upT, upR, upTh, upPh, a, p, e, x, En, Lz, Qc, r1, r2, r3, r4, z1, z2);
}

// void kerr_geo_mino_frequencies(double upT, double upR, double upTh, double upPh, double a, double p, double e, double x, double En, double Lz, double Qc){
// 	double r1, r2, r3, r4;
// 	kerr_geo_radial_roots(r1, r2, r3, r4, a, p, e, En, Lz, Qc);
// 	//std::cout << "r1 = "<<r1<<", r2 = "<<r2<<", r3 = "<<r3<<", r4 = "<<r4<<" \n";
// 	double z1, z2;
// 	kerr_geo_polar_roots(z1, z2, a, x, En, Lz, Qc);
// 	//std::cout << "z1 = "<<z1<<", z2 = "<<z2<<" \n";
// 	kerr_geo_mino_frequencies(upT, upR, upTh, upPh, a, p, e, x, En, Lz, Qc, r1, r2, r3, r4, z1, z2);
// }

void kerr_geo_mino_frequencies(double &upT, double &upR, double &upTh, double &upPh, double a, double p, double e, double , double En, double Lz, double Qc, double , double , double r3, double r4, double z1, double z2){
	kerr_geo_radial_mino_frequency(upR, a, p, e, En, r3, r4);
	if(std::abs(a) > 0){
		kerr_geo_polar_mino_frequency(upTh, a, En, z1, z2);
	}else{
		kerr_geo_polar_mino_frequency_schw(upTh, Lz, Qc);
	}
	kerr_geo_time_mino_frequency(upT, a, p, e, En, Lz, Qc, z1, z2, r3, r4, upR, upTh);
	kerr_geo_azimuthal_mino_frequency(upPh, a, p, e, En, Lz, Qc, z1, z2, r3, r4, upR, upTh);
}

void kerr_geo_carter_frequencies(double &cR, double &cTh, double &cPh, double upT, double upR, double upTh, double upPh, double a, double En, double Lz, double Qc, double z1, double z2){
	double upTTh, upPhTh;
	if(std::abs(a) > 0){
		kerr_geo_time_mino_frequency_polar(upTTh, a, En, Lz, Qc, z1, z2);
		kerr_geo_azimuthal_mino_frequency_polar(upPhTh, a, En, Lz, Qc, z1, z2);
	}else{
		kerr_geo_time_mino_frequency_polar_schw(upTTh, En, Lz, Qc, z1, z2);
		kerr_geo_azimuthal_mino_frequency_polar_schw(upPhTh, En, Lz, Qc, z1, z2);
	}
	double freqFactor = a*Lz - a*a*En - upTTh*upTh;
	double mFactor = Lz - a*En - upPhTh*upTh;
	cR = freqFactor*upR/upT;
	cTh = freqFactor*upTh/upT + upTh;
	cPh = freqFactor*upPh/upT - mFactor;
}

void kerr_geo_radial_mino_period(double &LaR, double a, double p, double e, double En, double r3, double r4){
	int Nsample = pow(2, 3);
	double sum = 0.;
	for(int i = 0; i < Nsample; i++){
		sum += dmino_dpsi(2.*M_PI*i/(Nsample), a, p, e, En, r3, r4);
	}
	LaR = 2.*M_PI*sum/Nsample;
	double periodCompare = 0;
	while(std::abs(1. - periodCompare/LaR) > FREQ_EPS && Nsample < pow(2, 15)){
		for(int i = 0; i < Nsample; i++){
			sum += dmino_dpsi(2.*M_PI*double(i)/double(Nsample) + M_PI/double(Nsample), a, p, e, En, r3, r4);
		}
		periodCompare = LaR;
		Nsample *= 2;
		LaR = 2.*M_PI*sum/Nsample;
	}
}
void kerr_geo_radial_mino_frequency(double &upR, double a, double p, double e, double En, double r3, double r4){
	kerr_geo_radial_mino_period(upR, a, p, e, En, r3, r4);
	upR = 2.*M_PI/upR;
}

void kerr_geo_polar_mino_period_schw(double &LaTh, double Lz, double Qc){
	LaTh = 2.*M_PI*dmino_dchi_schw(Lz, Qc);
}
void kerr_geo_polar_mino_frequency_schw(double &upTh, double Lz, double Qc){
	kerr_geo_polar_mino_period_schw(upTh, Lz, Qc);
	upTh = 2.*M_PI/upTh;
}

void kerr_geo_polar_mino_period(double &LaTh, double a, double En, double z1, double z2){
	int Nsample = pow(2, 3);
	double sum = 0.;
	for(int i = 0; i < Nsample; i++){
		sum += dmino_dchi(2.*M_PI*i/(Nsample), a, En, z1, z2);
	}
	LaTh = 2.*M_PI*sum/Nsample;
	double periodCompare = 0.;
	while(std::abs(1. - periodCompare/LaTh) > FREQ_EPS && Nsample < pow(2, 15)){
		for(int i = 0; i < Nsample; i++){
			sum += dmino_dchi(2.*M_PI*double(i)/double(Nsample) + M_PI/double(Nsample), a, En, z1, z2);
		}
		periodCompare = LaTh;
		Nsample *= 2;
		LaTh = 2.*M_PI*sum/Nsample;
	}
}
void kerr_geo_polar_mino_frequency(double &upTh, double a, double En, double z1, double z2){
	kerr_geo_polar_mino_period(upTh, a, En, z1, z2);
	upTh = 2.*M_PI/upTh;
}

void kerr_geo_time_mino_frequency_radial(double &upT, double a, double p, double e, double En, double Lz, double Qc, double r3, double r4){
	int Nsample = pow(2, 3);
	double sum = 0.;
	for(int i = 0; i < Nsample; i++){
		sum += dVtrdMino(rp_psi(2.*M_PI*i/(Nsample), p, e), a, En, Lz, Qc)*dmino_dpsi(2.*M_PI*i/(Nsample), a, p, e, En, r3, r4);
	}
	upT = sum/Nsample;
	double periodCompare = 0;
	while(std::abs(1. - periodCompare/upT) > FREQ_EPS && Nsample < pow(2, 15)){
		for(int i = 0; i < Nsample; i++){
			sum += dVtrdMino(rp_psi(2.*M_PI*double(i)/double(Nsample) + M_PI/double(Nsample), p, e), a, En, Lz, Qc)*dmino_dpsi(2.*M_PI*double(i)/double(Nsample) + M_PI/double(Nsample), a, p, e, En, r3, r4);
		}
		periodCompare = upT;
		Nsample *= 2;
		upT = sum/Nsample;
	}
}
void kerr_geo_time_mino_frequency_polar(double &upT, double a, double En, double Lz, double Qc, double z1, double z2){
	double x = sqrt(1. - z2*z2);
	int Nsample = pow(2, 3);
	double sum = 0.;
	for(int i = 0; i < Nsample; i++){
		sum += dVtzdMino(zp_chi(2.*M_PI*i/(Nsample), x), a, En, Lz, Qc)*dmino_dchi(2.*M_PI*i/(Nsample), a, En, z1, z2);
	}
	upT = sum/Nsample;
	double periodCompare = 0.;
	while(std::abs(1. - periodCompare/upT) > FREQ_EPS && Nsample < pow(2, 15)){
		for(int i = 0; i < Nsample; i++){
			sum += dVtzdMino(zp_chi(2.*M_PI*double(i)/double(Nsample) + M_PI/double(Nsample), x), a, En, Lz, Qc)*dmino_dchi(2.*M_PI*double(i)/double(Nsample) + M_PI/double(Nsample), a, En, z1, z2);
		}
		periodCompare = upT;
		Nsample *= 2;
		upT = sum/Nsample;
	}
}
void kerr_geo_time_mino_frequency_polar_schw(double &upT, double En, double Lz, double Qc, double , double z2){
	double x = sqrt(1. - z2*z2);
	int Nsample = pow(2, 3);
	double sum = 0.;
	for(int i = 0; i < Nsample; i++){
		sum += dVtzdMino(zp_chi(2.*M_PI*i/(Nsample), x), 0., En, Lz, Qc)*dmino_dchi_schw(Lz, Qc);
	}
	upT = sum/Nsample;
	double periodCompare = 0.;
	while(std::abs(1. - periodCompare/upT) > FREQ_EPS && Nsample < pow(2, 15)){
		for(int i = 0; i < Nsample; i++){
			sum += dVtzdMino(zp_chi(2.*M_PI*double(i)/double(Nsample) + M_PI/double(Nsample), x), 0., En, Lz, Qc)*dmino_dchi_schw(Lz, Qc);
		}
		periodCompare = upT;
		Nsample *= 2;
		upT = sum/Nsample;
	}
}
void kerr_geo_time_mino_frequency(double &upT, double a, double p, double e, double En, double Lz, double Qc, double z1, double z2, double r3, double r4, double upR, double upTh){
	double upTR, upTTh;
	kerr_geo_time_mino_frequency_radial(upTR, a, p, e, En, Lz, Qc, r3, r4);
	if(std::abs(a) > 0.){
		kerr_geo_time_mino_frequency_polar(upTTh, a, En, Lz, Qc, z1, z2);
	}else{
		kerr_geo_time_mino_frequency_polar_schw(upTTh, En, Lz, Qc, z1, z2);
	}
	upT = upTR*upR + upTTh*upTh;
}

void kerr_geo_azimuthal_mino_frequency_radial(double &upPh, double a, double p, double e, double En, double Lz, double Qc, double r3, double r4){
	int Nsample = pow(2, 3);
	double sum = 0.;
	for(int i = 0; i < Nsample; i++){
		sum += dVphirdMino(rp_psi(2.*M_PI*i/(Nsample), p, e), a,  En, Lz, Qc)*dmino_dpsi(2.*M_PI*i/(Nsample), a, p, e, En, r3, r4);
	}
	upPh = sum/Nsample;
	double periodCompare = 0.;
	while(std::abs(1. - periodCompare/upPh) > FREQ_EPS && Nsample < pow(2, 15)){
		for(int i = 0; i < Nsample; i++){
			sum += dVphirdMino(rp_psi(2.*M_PI*double(i)/double(Nsample) + M_PI/double(Nsample), p, e), a, En, Lz, Qc)*dmino_dpsi(2.*M_PI*double(i)/double(Nsample) + M_PI/double(Nsample), a, p, e, En, r3, r4);
		}
		periodCompare = upPh;
		Nsample *= 2;
		upPh = sum/Nsample;
	}
}
void kerr_geo_azimuthal_mino_frequency_polar(double &upPh, double a, double En, double Lz, double Qc, double z1, double z2){
	double x = sqrt(1. - z2*z2);
	int Nsample = pow(2, 5);
	double sum = 0.;
	for(int i = 0; i < Nsample; i++){
		sum += dVphizdMino(zp_chi(2.*M_PI*i/(Nsample), x), a, En, Lz, Qc)*dmino_dchi(2.*M_PI*i/(Nsample), a, En, z1, z2);
	}
	upPh = sum/Nsample;
	double periodCompare = 0.;
	while(std::abs(1. - periodCompare/upPh) > FREQ_EPS && Nsample < pow(2, 15)){
		for(int i = 0; i < Nsample; i++){
			sum += dVphizdMino(zp_chi(2.*M_PI*double(i)/double(Nsample) + M_PI/double(Nsample), x), a, En, Lz, Qc)*dmino_dchi(2.*M_PI*double(i)/double(Nsample) + M_PI/double(Nsample), a, En, z1, z2);
		}
		periodCompare = upPh;
		Nsample *= 2;
		upPh = sum/Nsample;
	}
}
void kerr_geo_azimuthal_mino_frequency_polar_schw(double &upPh, double En, double Lz, double Qc, double , double z2){
	double x = sqrt(1. - z2*z2);
	int Nsample = pow(2, 5);
	double sum = 0.;
	for(int i = 0; i < Nsample; i++){
		sum += dVphizdMino(zp_chi(2.*M_PI*i/(Nsample), x), 0., En, Lz, Qc)*dmino_dchi_schw(Lz, Qc);
	}
	upPh = sum/Nsample;
	double periodCompare = 0.;
	while(std::abs(1. - periodCompare/upPh) > FREQ_EPS && Nsample < pow(2, 15)){
		for(int i = 0; i < Nsample; i++){
			sum += dVphizdMino(zp_chi(2.*M_PI*double(i)/double(Nsample) + M_PI/double(Nsample), x), 0., En, Lz, Qc)*dmino_dchi_schw(Lz, Qc);
		}
		periodCompare = upPh;
		Nsample *= 2;
		upPh = sum/Nsample;
	}
}
void kerr_geo_azimuthal_mino_frequency(double &upPh, double a, double p, double e, double En, double Lz, double Qc, double z1, double z2, double r3, double r4, double upR, double upTh){
	double upPhR, upPhTh;
	kerr_geo_azimuthal_mino_frequency_radial(upPhR, a, p, e, En, Lz, Qc, r3, r4);
	if(std::abs(a) > 0.){
		kerr_geo_azimuthal_mino_frequency_polar(upPhTh, a, En, Lz, Qc, z1, z2);
	}else{
		kerr_geo_azimuthal_mino_frequency_polar_schw(upPhTh, En, Lz, Qc, z1, z2);
	}
	upPh = upPhR*upR + upPhTh*upTh;
}

Vector angle_vector(int Nsample){
	Vector angle(Nsample);
	for(int i = 0; i < Nsample; i++){
		angle[i] = 2.*M_PI*double(i)/double(Nsample);
	}
	return angle;
}
void angle_vector(Vector &angle){
	int Nsample = angle.size();
	for(int i = 0; i < Nsample; i++){
		angle[i] = 2.*M_PI*double(i)/double(Nsample);
	}
}
Vector angle_vector_half(int Nsample){
	Vector angle(Nsample);
	for(int i = 0; i < Nsample; i++){
		angle[i] = M_PI*double(i)/double(Nsample - 1);
	}
	return angle;
}
void angle_vector_half(Vector &angle){
	int Nsample = angle.size();
	for(int i = 0; i < Nsample; i++){
		angle[i] = M_PI*double(i)/double(Nsample - 1);
	}
}

Vector mino_of_psi_fourier(double a, double p, double e, double En, double r3, double r4){
	// test to see if the n = 2 harmonic converges
	int Nsample = pow(2, 3);
	int nTest = 2;
	double endpoints = dmino_dpsi(0., a, p, e, En, r3, r4) + pow(-1., nTest)*dmino_dpsi(M_PI, a, p, e, En, r3, r4);
	double sum = 0.;
	for(int i = 1; i < Nsample; i++){
		sum += 2.*dmino_dpsi(M_PI*double(i)/double(Nsample), a, p, e, En, r3, r4)*cos(nTest*M_PI*double(i)/double(Nsample));
	}
	double fourierTest = (sum + endpoints)/double(Nsample);
	double fourierCompare = 0;
	while((std::abs(fourierCompare - fourierTest) > FOURIER_EPS && std::abs(1. - fourierCompare/fourierTest) > FOURIER_EPS) && Nsample < pow(2, 15)){
		for(int i = 0; i < Nsample; i++){
			sum += 2.*dmino_dpsi((M_PI*double(i)/double(Nsample) + M_PI/double(2*Nsample)), a, p, e, En, r3, r4)*cos(nTest*(M_PI*double(i)/double(Nsample) + M_PI/double(2*Nsample)));
		}
		fourierCompare = fourierTest;
		Nsample *= 2;
		fourierTest = (sum + endpoints)/double(Nsample);
	}

	Nsample += 1;
	Vector fourier(Nsample);

	int convergeFlag = 0;
	int i = 0;
	while(i < Nsample && convergeFlag < 2){
		fourier[i] = dmino_dpsi(0., a, p, e, En, r3, r4) + pow(-1., i)*dmino_dpsi(M_PI, a, p, e, En, r3, r4);
		for(int j = 1; j < Nsample - 1; j++){
			fourier[i] += 2.*dmino_dpsi(M_PI*double(j)/double(Nsample - 1), a, p, e, En, r3, r4)*cos(M_PI*double(i*j)/double(Nsample - 1));
		}
		fourier[i] /= double(Nsample - 1);
		if(i > 0){
			// std::cout << "f_" << i << " = " << std::abs(fourier[i])/i << "\n";
			if(std::abs(fourier[i])/i < std::abs(fourier[1])*FOURIER_EPS){
				convergeFlag += 1; // if two consecutive Fourier coefficients fall below threshold, cutoff construction of higher-frequency coefficients
			}else{
				convergeFlag = 0;
			}
		}
		i++;
	}

	return fourier;
}

Vector mino_of_chi_schw_fourier(double Lz, double Qc){
	// test to see if the n = 5 harmonic converges
	int Nsample = pow(2, 3);
	int nTest = 2;
	double dminodchi = dmino_dchi_schw(Lz, Qc);
	double endpoints = dminodchi + pow(-1., nTest)*dminodchi;
	double sum = 0.;
	for(int i = 1; i < Nsample; i++){
		sum += 2.*dminodchi*cos(nTest*M_PI*double(i)/double(Nsample));
	}
	double fourierTest = (sum + endpoints)/double(Nsample);
	double fourierCompare = 0;
	while((std::abs(fourierCompare - fourierTest) > FOURIER_EPS && std::abs(1. - fourierCompare/fourierTest) > FOURIER_EPS) && Nsample < pow(2, 15)){
		for(int i = 0; i < Nsample; i++){
			sum += 2.*dminodchi*cos(nTest*(M_PI*double(i)/double(Nsample) + M_PI/double(2*Nsample)));
		}
		fourierCompare = fourierTest;
		Nsample *= 2;
		fourierTest = (sum + endpoints)/double(Nsample);
	}

	// generate fourier coefficients
	Nsample += 1;
	Vector fourier(Nsample);

	int convergeFlag = 0;
	int i = 0;
	while(i < Nsample && convergeFlag < 2){
		fourier[i] = dminodchi + pow(-1., i)*dminodchi;
		for(int j = 1; j < Nsample - 1; j++){
			fourier[i] += 2.*dminodchi*cos(M_PI*double(i*j)/double(Nsample - 1));
		}
		fourier[i] /= double(Nsample - 1);
		if(i > 0){
			if(std::abs(fourier[i])/i < std::abs(fourier[1])*FOURIER_EPS){
				convergeFlag += 1; // if two consecutive Fourier coefficients fall below threshold, cutoff construction of higher-frequency coefficients
			}else{
				convergeFlag = 0;
			}
		}
		i++;
	}

	return fourier;
}

Vector mino_of_chi_fourier(double a, double En, double z1, double z2){
	// test to see if the n = 5 harmonic converges
	int Nsample = pow(2, 3);
	int nTest = 2;
	double endpoints = dmino_dchi(0., a, En, z1, z2) + pow(-1., nTest)*dmino_dchi(M_PI, a, En, z1, z2);
	double sum = 0.;
	for(int i = 1; i < Nsample; i++){
		sum += 2.*dmino_dchi(M_PI*double(i)/double(Nsample), a, En, z1, z2)*cos(nTest*M_PI*double(i)/double(Nsample));
	}
	double fourierTest = (sum + endpoints)/double(Nsample);
	double fourierCompare = 0;
	while((std::abs(fourierCompare - fourierTest) > FOURIER_EPS && std::abs(1. - fourierCompare/fourierTest) > FOURIER_EPS) && Nsample < pow(2, 15)){
		for(int i = 0; i < Nsample; i++){
			sum += 2.*dmino_dchi((M_PI*double(i)/double(Nsample) + M_PI/double(2*Nsample)), a, En, z1, z2)*cos(nTest*(M_PI*double(i)/double(Nsample) + M_PI/double(2*Nsample)));
		}
		fourierCompare = fourierTest;
		Nsample *= 2;
		fourierTest = (sum + endpoints)/double(Nsample);
	}

	// generate fourier coefficients
	Nsample += 1;
	Vector fourier(Nsample);

	int convergeFlag = 0;
	int i = 0;
	while(i < Nsample && convergeFlag < 2){
		fourier[i] = dmino_dchi(0., a, En, z1, z2) + pow(-1., i)*dmino_dchi(M_PI, a, En, z1, z2);
		for(int j = 1; j < Nsample - 1; j++){
			fourier[i] += 2.*dmino_dchi(M_PI*double(j)/double(Nsample - 1), a, En, z1, z2)*cos(M_PI*double(i*j)/double(Nsample - 1));
		}
		fourier[i] /= double(Nsample - 1);
		if(i > 0){
			if(std::abs(fourier[i])/i < std::abs(fourier[1])*FOURIER_EPS){
				convergeFlag += 1; // if two consecutive Fourier coefficients fall below threshold, cutoff construction of higher-frequency coefficients
			}else{
				convergeFlag = 0;
			}
		}
		i++;
	}

	return fourier;
}

double integrated_dct_sum(double angle, Vector &fourier){
	double sum = 0.5*angle*fourier[0];
	int Nsample = fourier.size();
	for(int n = 1; n < Nsample - 1; n++){
		sum += fourier[n]*sin(double(n)*angle)/double(n);
	}
	sum += 0.5*fourier[Nsample - 1]*sin(double(Nsample - 1)*angle)/double(Nsample - 1);

	return sum;
}

double mino_of_kepler_phase(double angle, Vector &fourier){
	return integrated_dct_sum(angle, fourier);
}

Vector kepler_phase_of_angle_fourier(Vector &fourier){
	double freq = 2./fourier[0];
	int Nsample = pow(2, 2);
	int nTest = 2;
	double endpoints = 1. + pow(-1., nTest);
	double sum = 0.;
	for(int i = 1; i < Nsample; i++){
		sum += 2.*cos(double(nTest)*freq*mino_of_kepler_phase(M_PI*double(i)/double(Nsample), fourier));
	}
	double fourierTest = (sum + endpoints)/double(Nsample);
	double fourierCompare = 0;
	while((std::abs(fourierCompare - fourierTest) > FOURIER_EPS && std::abs(1. - fourierCompare/fourierTest) > FOURIER_EPS) && Nsample < pow(2, 15)){
		for(int i = 0; i < Nsample; i++){
			sum += 2.*cos(double(nTest)*freq*mino_of_kepler_phase(M_PI*double(i)/double(Nsample) + M_PI/double(2*Nsample), fourier));
		}
		fourierCompare = fourierTest;
		Nsample *= 2;
		fourierTest = (sum + endpoints)/double(Nsample);
	}

	// std::cout << "psi of Mino samples = " << Nsample << "\n";
	// generate fourier coefficients
	Nsample += 1;
	Vector psi_fourier(Nsample);
	Vector mino_phase(Nsample);
	for(int j = 1; j < Nsample - 1; j++){
		mino_phase[j] = mino_of_kepler_phase(M_PI*double(j)/double(Nsample - 1), fourier);
	}

	psi_fourier[0] = 2.;

	int i0 = 4;
	for(int i = 0; i < i0; i++){
		psi_fourier[i] = 1. + pow(-1., i);
		for(int j = 1; j < Nsample - 1; j++){
			psi_fourier[i] += 2.*cos(double(i)*freq*mino_phase[j]);
		}
		psi_fourier[i] /= double(Nsample - 1);
		// std::cout << "f_"<<i<<" = "<< 0.5*psi_fourier[i] << "\n";
	}
	int i = i0;
	double tol = 1.e-14;
	while(i < Nsample - 1 && (std::abs(0.5*psi_fourier[i - 1]/double(i - 1)) > tol || std::abs(0.5*psi_fourier[i - 2]/double(i - 2)) > tol)){
		psi_fourier[i] = 1. + pow(-1., i);
		for(int j = 1; j < Nsample - 1; j++){
			psi_fourier[i] += 2.*cos(double(i)*freq*mino_phase[j]);
		}
		psi_fourier[i] /= double(Nsample - 1);
		// std::cout << "f_"<<i<<" = "<< 0.5*psi_fourier[i] << "\n";
		if(std::abs(psi_fourier[i - 2]/double(i - 2)) < std::abs(psi_fourier[i]/double(i)) && std::abs(psi_fourier[i - 3]/double(i - 3)) < std::abs(psi_fourier[i - 1]/double(i - 1))){
			tol = 1.; // if the weighted (due to integration) fourier coefficients are not falling off, terminate
		}
		i++;
	}

	Vector psi_fourier_return(i);
	for(int k = 0; k < i; k++){
		psi_fourier_return[k] = psi_fourier[k];
	}

	return psi_fourier_return;
}

// double kepler_phase_of_angle(double angle, Vector &fourier){
// 	double kepler_phase = angle;
// 	for(int j = 1; j < fourier.size(); j++){
// 		kepler_phase += 2.*fourier[j]/double(j)*sin(j*angle);
// 	}
// 	return kepler_phase;
// }

double kepler_phase_of_angle(double angle, Vector &fourier){
	return integrated_dct_sum(angle, fourier);
}

double rp_of_angle(double angle, double p, double e, Vector &fourier){
	// double test1 = rp_psi(kepler_phase_of_angle(angle, fourier), p, e);
	// double test2 = rp_psi(kepler_phase_of_angle(angle, fourier), p, e);
	// if(std::abs(test2-test1) > 0.){
	// 	std::cout << "ERROR!!!\n";
	// }
	double psi = kepler_phase_of_angle(angle, fourier);
	return rp_psi(psi, p, e);
}

double zp_of_angle(double angle, double x, Vector &fourier){
	return zp_chi(kepler_phase_of_angle(angle, fourier), x);
}

Vector tp_radial_of_angle_fourier(double a, double p, double e, double En, double Lz, double Qc, double upR, Vector &fourier){
	// test to see if the n = 2 harmonic converges

	// std::cout << "Using fourier coefficients = " << fourier[0] << ", " << fourier[1] << ", " << fourier[2] << ", " << fourier[3] << ", ... \n";
	// std::cout << fourier.size() << "\n";
	int Nsample = 8;
	int nTest = 4;

	// std::cout << "Using radial points = " << rp_of_angle(M_PI*double(0)/double(Nsample), p, e, fourier) << ", " << rp_of_angle(M_PI*double(1)/double(Nsample), p, e, fourier) << ", " << rp_of_angle(M_PI*double(2)/double(Nsample), p, e, fourier) << ", " << rp_of_angle(M_PI*double(3)/double(Nsample), p, e, fourier) << ", ... \n";
	// std::cout << "Using radial points = " << rp_of_angle(M_PI*double(0)/double(Nsample), p, e, fourier) << ", " << rp_of_angle(M_PI*double(1)/double(Nsample), p, e, fourier) << ", " << rp_of_angle(M_PI*double(2)/double(Nsample), p, e, fourier) << ", " << rp_of_angle(M_PI*double(3)/double(Nsample), p, e, fourier) << ", ... \n";
	// std::cout << "Using radial points = " << rp_of_angle(M_PI*double(0)/double(Nsample), p, e, fourier) << ", " << M_PI*double(1)/double(Nsample) << ", " << p << ", " << e << ", ... \n";

	double endpoints = dVtrdMino(rp_of_angle(0., p, e, fourier), a, En, Lz, Qc) + pow(-1., nTest)*dVtrdMino(rp_of_angle(M_PI, p, e, fourier), a, En, Lz, Qc);
	double sum = 0.;
	for(int i = 1; i < Nsample; i++){
		sum += 2.*dVtrdMino(rp_of_angle(M_PI*double(i)/double(Nsample), p, e, fourier), a, En, Lz, Qc)*cos(nTest*M_PI*double(i)/double(Nsample));
	}
	double fourierTest = (sum + endpoints)/double(Nsample)/upR;
	double fourierCompare = 0;
	// std::cout << "Nsample = " << Nsample << ", test = " << fourierTest << "\n";
	// std::cout << "Nsample = " << Nsample << ", comp = " << fourierTest - fourierCompare << "\n";
	while((std::abs(fourierCompare - fourierTest) > FOURIER_EPS && std::abs(1. - fourierCompare/fourierTest) > FOURIER_EPS) && Nsample < pow(2, 15)){
		for(int i = 0; i < Nsample; i++){
			sum += 2.*dVtrdMino(rp_of_angle((M_PI*double(i)/double(Nsample) + M_PI/double(2*Nsample)), p, e, fourier), a, En, Lz, Qc)*cos(nTest*(M_PI*double(i)/double(Nsample) + M_PI/double(2*Nsample)));
		}
		fourierCompare = fourierTest;
		Nsample *= 2;
		fourierTest = (sum + endpoints)/double(Nsample)/upR;
		// std::cout << "Nsample = " << Nsample << ", test = " << fourierTest << "\n";
		// std::cout << "Nsample = " << Nsample << ", comp = " << fourierTest - fourierCompare << "\n";
	}

	// std::cout << "tr samples = " << Nsample << "\n";
	// generate fourier coefficients
	Nsample += 1;
	Vector fourierT(Nsample);
	Vector angle = angle_vector_half(Nsample);
	Vector dVtr(Nsample);
	dVtr[0] = dVtrdMino(rp_of_angle(0., p, e, fourier), a, En, Lz, Qc);
	dVtr[Nsample - 1] = dVtrdMino(rp_of_angle(M_PI, p, e, fourier), a, En, Lz, Qc);
	for(int j = 1; j < Nsample - 1; j++){
		dVtr[j] = dVtrdMino(rp_of_angle(angle[j], p, e, fourier), a, En, Lz, Qc);
	}

	for(int i = 1; i < Nsample; i++){
		fourierT[i] = dVtr[0] + pow(-1., i)*dVtr[Nsample - 1];
		for(int j = 1; j < Nsample - 1; j++){
			fourierT[i] += 2.*dVtr[j]*cos(i*angle[j]);
		}
		fourierT[i] /= upR*double(Nsample - 1);
	}
	fourierT[0] = 0.;

	return fourierT;
}

Vector tp_polar_of_angle_fourier(double a, double x, double En, double Lz, double Qc, double upTh, Vector &fourier){
	// test to see if the n = 5 harmonic converges
	int Nsample = pow(2, 3);
	int nTest = 4;
	double endpoints = dVtzdMino(zp_of_angle(0., x, fourier), a, En, Lz, Qc) + pow(-1., nTest)*dVtzdMino(zp_of_angle(M_PI, x, fourier), a, En, Lz, Qc);
	double sum = 0.;
	for(int i = 1; i < Nsample; i++){
		sum += 2.*dVtzdMino(zp_of_angle(M_PI*double(i)/double(Nsample), x, fourier), a, En, Lz, Qc)*cos(nTest*M_PI*double(i)/double(Nsample));
	}
	double fourierTest = (sum + endpoints)/double(Nsample)/upTh;
	double fourierCompare = 0;
	while((std::abs(fourierCompare - fourierTest) > FOURIER_EPS && std::abs(1. - fourierCompare/fourierTest) > FOURIER_EPS) && Nsample < pow(2, 15)){
		for(int i = 0; i < Nsample; i++){
			sum += 2.*dVtzdMino(zp_of_angle((M_PI*double(i)/double(Nsample) + M_PI/double(2*Nsample)), x, fourier), a, En, Lz, Qc)*cos(nTest*(M_PI*double(i)/double(Nsample) + M_PI/double(2*Nsample)));
		}
		fourierCompare = fourierTest;
		Nsample *= 2;
		fourierTest = (sum + endpoints)/double(Nsample)/upTh;
	}

	// generate fourier coefficients
	Nsample += 1;
	Vector fourierT(Nsample);
	Vector angle = angle_vector_half(Nsample);

	for(int i = 1; i < Nsample; i++){
		fourierT[i] = dVtzdMino(zp_of_angle(0., x, fourier), a, En, Lz, Qc) + pow(-1., i)*dVtzdMino(zp_of_angle(M_PI, x, fourier), a, En, Lz, Qc);
		for(int j = 1; j < Nsample - 1; j++){
			fourierT[i] += 2.*dVtzdMino(zp_of_angle(angle[j], x, fourier), a, En, Lz, Qc)*cos(i*angle[j]);
		}
		fourierT[i] /= upTh*double(Nsample - 1);
	}
	fourierT[0] = 0.;

	return fourierT;
}

double tp_of_angle(double angle, Vector &fourier){
	return integrated_dct_sum(angle, fourier);
}

Vector phip_radial_of_angle_fourier(double a, double p, double e, double En, double Lz, double Qc, double upR, Vector &fourier){
	// test to see if the n = 2 harmonic converges
	int Nsample = pow(2, 3);
	int nTest = 4;
	double endpoints = dVphirdMino(rp_of_angle(0., p, e, fourier), a, En, Lz, Qc) + pow(-1., nTest)*dVphirdMino(rp_of_angle(M_PI, p, e, fourier), a, En, Lz, Qc);
	double sum = 0.;
	for(int i = 1; i < Nsample; i++){
		sum += 2.*dVphirdMino(rp_of_angle(M_PI*double(i)/double(Nsample), p, e, fourier), a, En, Lz, Qc)*cos(nTest*M_PI*double(i)/double(Nsample));
	}
	double fourierTest = (sum + endpoints)/double(Nsample)/upR;
	double fourierCompare = 0;
	while((std::abs(fourierCompare - fourierTest) > FOURIER_EPS && std::abs(1. - fourierCompare/fourierTest) > FOURIER_EPS) && Nsample < pow(2, 15)){
		for(int i = 0; i < Nsample; i++){
			sum += 2.*dVphirdMino(rp_of_angle((M_PI*double(i)/double(Nsample) + M_PI/double(2*Nsample)), p, e, fourier), a, En, Lz, Qc)*cos(nTest*(M_PI*double(i)/double(Nsample) + M_PI/double(2*Nsample)));
		}
		fourierCompare = fourierTest;
		Nsample *= 2;
		fourierTest = (sum + endpoints)/double(Nsample)/upR;
	}

	// generate fourier coefficients
	Nsample += 1;
	Vector fourierPhi(Nsample);
	Vector angle = angle_vector_half(Nsample);
	Vector dVphir(Nsample);
	dVphir[0] = dVphirdMino(rp_of_angle(0., p, e, fourier), a, En, Lz, Qc);
	dVphir[Nsample - 1] = dVphirdMino(rp_of_angle(M_PI, p, e, fourier), a, En, Lz, Qc);
	for(int j = 1; j < Nsample - 1; j++){
		dVphir[j] = dVphirdMino(rp_of_angle(angle[j], p, e, fourier), a, En, Lz, Qc);
	}

	for(int i = 1; i < Nsample; i++){
		fourierPhi[i] = dVphir[0] + pow(-1., i)*dVphir[Nsample-1];
		for(int j = 1; j < Nsample - 1; j++){
			fourierPhi[i] += 2.*dVphir[j]*cos(i*angle[j]);
		}
		fourierPhi[i] /= upR*double(Nsample - 1);
	}
	fourierPhi[0] = 0.;

	return fourierPhi;
}

Vector phip_polar_of_angle_fourier(double a, double x, double En, double Lz, double Qc, double upTh, Vector &fourier){
	// test to see if the n = 5 harmonic converges
	int Nsample = pow(2, 3);
	int nTest = 4; // test needs to be on an even harmonic for polar dependence
	double endpoints = dVphizdMino(zp_of_angle(0., x, fourier), a, En, Lz, Qc) + pow(-1., nTest)*dVphizdMino(zp_of_angle(M_PI, x, fourier), a, En, Lz, Qc);
	double sum = 0.;
	for(int i = 1; i < Nsample; i++){
		sum += 2.*dVphizdMino(zp_of_angle(M_PI*double(i)/double(Nsample), x, fourier), a, En, Lz, Qc)*cos(nTest*M_PI*double(i)/double(Nsample));
	}
	double fourierTest = (sum + endpoints)/double(Nsample)/upTh;
	double fourierCompare = 0;
	while((std::abs(fourierCompare - fourierTest) > FOURIER_EPS && std::abs(1. - fourierCompare/fourierTest) > FOURIER_EPS) && Nsample < pow(2, 15)){
		for(int i = 0; i < Nsample; i++){
			sum += 2.*dVphizdMino(zp_of_angle((M_PI*double(i)/double(Nsample) + M_PI/double(2*Nsample)), x, fourier), a, En, Lz, Qc)*cos(nTest*(M_PI*double(i)/double(Nsample) + M_PI/double(2*Nsample)));
		}
		fourierCompare = fourierTest;
		Nsample *= 2;
		fourierTest = (sum + endpoints)/double(Nsample)/upTh;
	}

	// generate fourier coefficients
	Nsample += 1;
	Vector fourierPhi(Nsample);
	Vector angle = angle_vector_half(Nsample);

	for(int i = 1; i < Nsample; i++){
		fourierPhi[i] = dVphizdMino(zp_of_angle(0., x, fourier), a, En, Lz, Qc) + pow(-1., i)*dVphizdMino(zp_of_angle(M_PI, x, fourier), a, En, Lz, Qc);
		for(int j = 1; j < Nsample - 1; j++){
			fourierPhi[i] += 2.*dVphizdMino(zp_of_angle(angle[j], x, fourier), a, En, Lz, Qc)*cos(i*angle[j]);
		}
		fourierPhi[i] /= upTh*double(Nsample - 1);
	}
	fourierPhi[0] = 0.;

	return fourierPhi;
}

double phip_of_angle(double angle, Vector &fourier){
	return integrated_dct_sum(angle, fourier);
}

void kerr_trajectory(Vector& tR, Vector& tTh, Vector& rp, Vector& thetap, Vector& phiR, Vector& phiTh, double p, double e, double x, Vector fourier_tr, Vector fourier_tz, Vector fourier_psi, Vector fourier_chi, Vector fourier_phir, Vector fourier_phiz){
	int halfSample = tR.size();
	int Nsample = 2*(halfSample - 1);
	if(e == 0.){
		if(std::abs(x) == 1.){
			for(int i = 0; i < halfSample; i++){
				tR[i] = 0.;
				tTh[i] = 0.;
				phiR[i] = 0.;
				phiTh[i] = 0.;
				rp[i] = p;
				thetap[i] = 0.5*M_PI;
			}
		}else{
			for(int i = 0; i < halfSample; i++){
				tR[i] = 0.;
				tTh[i] = tp_of_angle(2.*M_PI*double(i)/double(Nsample), fourier_tz);
				phiR[i] = 0.;
				phiTh[i] = phip_of_angle(2.*M_PI*double(i)/double(Nsample), fourier_phiz);
				rp[i] = p;
				thetap[i] = acos(zp_of_angle(2.*M_PI*double(i)/double(Nsample), x, fourier_chi));
			}
		}
	}else if(std::abs(x) == 1.){
		for(int i = 0; i < halfSample; i++){
			tR[i] = tp_of_angle(2.*M_PI*double(i)/double(Nsample), fourier_tr);
			tTh[i] = 0.;
			phiR[i] = phip_of_angle(2.*M_PI*double(i)/double(Nsample), fourier_phir);
			phiTh[i] = 0.;
			rp[i] = rp_of_angle(2.*M_PI*double(i)/double(Nsample), p, e, fourier_psi);
			thetap[i] = 0.5*M_PI;
		}
	}else{
		for(int i = 0; i < halfSample; i++){
			tR[i] = tp_of_angle(2.*M_PI*double(i)/double(Nsample), fourier_tr);
			tTh[i] = tp_of_angle(2.*M_PI*double(i)/double(Nsample), fourier_tz);
			phiR[i] = phip_of_angle(2.*M_PI*double(i)/double(Nsample), fourier_phir);
			phiTh[i] = phip_of_angle(2.*M_PI*double(i)/double(Nsample), fourier_phiz);
			rp[i] = rp_of_angle(2.*M_PI*double(i)/double(Nsample), p, e, fourier_psi);
			thetap[i] = acos(zp_of_angle(2.*M_PI*double(i)/double(Nsample), x, fourier_chi));
		}
	}
}

// void mino_of_psi_test(){
// 	double a, p, e, x;
// 	a = 0.9;
// 	p = 10.;
// 	e = 0.5;
// 	x = 0.5;
// 	std::cout << std::setprecision(15);
//
// 	double En, Lz, Qc;
// 	kerr_geo_orbital_constants(En, Lz, Qc, a, p, e, x);
// 	//std::cout << "En = "<<En<<", Lz = "<<Lz<<", Qc = "<<Qc<<" \n";
// 	double r1, r2, r3, r4;
// 	kerr_geo_radial_roots(r1, r2, r3, r4, a, p, e, En, Lz, Qc);
// 	//std::cout << "r1 = "<<r1<<", r2 = "<<r2<<", r3 = "<<r3<<", r4 = "<<r4<<" \n";
// 	double z1, z2;
// 	kerr_geo_polar_roots(z1, z2, a, x, En, Lz, Qc);
// 	//std::cout << "z1 = "<<z1<<", z2 = "<<z2<<" \n";
//
// 	double upT, upR, upPh, upTh;
// 	kerr_geo_mino_frequencies(upT, upR, upTh, upPh, a, p, e, x, En, Lz, Qc, r1, r2, r3, r4, z1, z2);
//
// 	Vector fourier_radial = mino_of_psi_fourier(a, p, e, En, r3, r4);
// 	Vector fourier_polar = mino_of_chi_fourier(a, En, z1, z2);
// 	Vector fourier_psi = kepler_phase_of_angle_fourier(fourier_radial);
// 	Vector fourier_chi = kepler_phase_of_angle_fourier(fourier_polar);
// 	Vector fourier_tr = tp_radial_of_angle_fourier(a, p, e, En, Lz, Qc, upR, fourier_psi);
// 	Vector fourier_tz = tp_polar_of_angle_fourier(a, x, En, Lz, Qc, upTh, fourier_chi);
// 	Vector fourier_phir = phip_radial_of_angle_fourier(a, p, e, En, Lz, Qc, upR, fourier_psi);
// 	Vector fourier_phiz = phip_polar_of_angle_fourier(a, x, En, Lz, Qc, upTh, fourier_chi);
// 	double la = 3.;
// 	double rp = rp_of_angle(la*upR, p, e, fourier_psi);
// 	double zp = zp_of_angle(la*upTh, x, fourier_chi);
// 	double tp = upT*la + tp_of_angle(la*upR, fourier_tr) + tp_of_angle(la*upTh, fourier_tz);
// 	double phip = upPh*la + phip_of_angle(la*upR, fourier_phir) + phip_of_angle(la*upTh, fourier_phiz);
// 	std::cout << "tp(la = "<<la<<") = "<< tp <<"\n";
// 	std::cout << "rp(la = "<<la<<") = "<< rp <<"\n";
// 	std::cout << "thp(la = "<<la<<") = "<< acos(zp) <<"\n";
// 	std::cout << "phip(la = "<<la<<") = "<< phip <<"\n";
// }

GeodesicSource kerr_geo_orbit(double a, double p, double e, double x, int Nsample = pow(2, 9)){
	// if(e == 0. && std::abs(x) == 1.){
	// 	GeodesicSource geo = kerr_geo_circ(a, p, sgn(x));
	// 	GeodesicTrajectory traj = geo.getTrajectory();
	// 	double tR = geo.tR;
	// 	double tTheta = geo.tTheta;
	// 	double r = geo.r;
	// 	double theta = geo.theta;
	// 	double phiR = geo.phiR;
	// 	double phiTheta = geo.phiTheta;
	//
	// 	return geo;
	// }
	double En = 0., Lz = 0., Qc = 0.;
	double r1 = 0., r2 = 0., r3 = 0., r4 = 0.;
	double z1 = 0., z2 = 0.;
	double upT = 0., upR = 0., upPh = 0., upTh = 0.;
	double cR = 0., cTh = 0., cPh = 0.;
	if(e == 0. && std::abs(x) == 1.){
		En = kerr_geo_energy_circ(a*x, p);
		Lz = kerr_geo_momentum_circ(a*x, p);
		upT = kerr_geo_time_frequency_circ(a*x, p);
		upPh = x*kerr_geo_azimuthal_frequency_circ(a*x, p);
	}else{
		kerr_geo_orbital_constants(En, Lz, Qc, a, p, e, x);
		// std::cout << "Orbital constants have been calculated \n";
		kerr_geo_radial_roots(r1, r2, r3, r4, a, p, e, En, Lz, Qc);
		// std::cout << "Radial roots have been calculated \n";
		kerr_geo_polar_roots(z1, z2, a, x, En, Lz, Qc);
		// std::cout << "Polar roots have been calculated \n";
		kerr_geo_mino_frequencies(upT, upR, upTh, upPh, a, p, e, x, En, Lz, Qc, r1, r2, r3, r4, z1, z2);
		// std::cout << "Frequencies have been calculated \n";
		kerr_geo_carter_frequencies(cR, cTh, cPh, upT, upR, upTh, upPh, a, En, Lz, Qc, z1, z2);
		// std::cout << "Carter frequencies have been calculated \n";
	}
	if( std::abs(x*x + z2 - 1.) > 1.e-10 ){
		z2 = sqrt(1. - x*x);
	}

	Vector fourier_radial, fourier_psi, fourier_tr, fourier_phir;
	if(e != 0.){
		fourier_radial = mino_of_psi_fourier(a, p, e, En, r3, r4);
		// std::cout << "Mino of psi has been calculated \n";
		fourier_psi = kepler_phase_of_angle_fourier(fourier_radial);
		// std::cout << "psi of Mino has been calculated \n";
		fourier_tr = tp_radial_of_angle_fourier(a, p, e, En, Lz, Qc, upR, fourier_psi);
		// std::cout << "tr of Mino has been calculated \n";
		fourier_phir = phip_radial_of_angle_fourier(a, p, e, En, Lz, Qc, upR, fourier_psi);
		// std::cout << "phir of Mino has been calculated \n";
	}

	Vector fourier_polar, fourier_chi, fourier_tz, fourier_phiz;
	if(std::abs(x) != 1.){
		fourier_polar = mino_of_chi_fourier(a, En, z1, z2);
		fourier_chi = kepler_phase_of_angle_fourier(fourier_polar);
		fourier_tz = tp_polar_of_angle_fourier(a, x, En, Lz, Qc, upTh, fourier_chi);
		fourier_phiz = phip_polar_of_angle_fourier(a, x, En, Lz, Qc, upTh, fourier_chi);
	}
	if(fourier_psi.size() == 0){
		fourier_psi = Vector(3, 0.);
	}
	if(fourier_chi.size() == 0){
		fourier_chi = Vector(3, 0.);
	}
	if(fourier_tr.size() == 0){
		fourier_tr = Vector(3, 0.);
	}
	if(fourier_tz.size() == 0){
		fourier_tz = Vector(3, 0.);
	}
	if(fourier_phir.size() == 0){
		fourier_phir = Vector(3, 0.);
	}
	if(fourier_phiz.size() == 0){
		fourier_phiz = Vector(3, 0.);
	}

	int halfSample = Nsample/2 + 1;
	Vector tR(halfSample);
	Vector tTh(halfSample);
	Vector phiR(halfSample);
	Vector phiTh(halfSample);
	Vector rp(halfSample);
	Vector thetap(halfSample);

	kerr_trajectory(tR, tTh, rp, thetap, phiR, phiTh, p, e, x, fourier_tr, fourier_tz, fourier_psi, fourier_chi, fourier_phir, fourier_phiz);

	GeodesicSource geo;
	geo.setConstants(a, p, e, x, En, Lz, Qc, r1, r2, r3, r4, z1, z2, upT, upR, upTh, upPh, cR, cTh, cPh);
	geo.setTrajectory(tR, tTh, rp, thetap, phiR, phiTh);
	geo.setCoefficients(fourier_tr, fourier_tz, fourier_psi, fourier_chi, fourier_phir, fourier_phiz);

	return geo;
}

////////////////////////////////
// Special functions approach //
////////////////////////////////

double kerr_geo_radial_mino_frequency_sf(double En, double r1, double r2, double r3, double r4){
	double upR;
	kerr_geo_radial_mino_frequency_sf(upR, En, r1, r2, r3, r4);
	return upR;
}
double kerr_geo_polar_mino_frequency_sf(double a, double En, double z1, double z2){
	double upTh;
	kerr_geo_polar_mino_frequency_sf(upTh, a, En, z1, z2);
	return upTh;
}

void kerr_geo_radial_mino_frequency_sf(double &upR, double En, double r1, double r2, double r3, double r4){
	upR = 0.5*M_PI*sqrt((1. - En*En)*(r1 - r3)*(r2 - r4))/elliptic_k(sqrt((r1 - r2)*(r3 - r4)/(r1 - r3)/(r2 - r4)));
}

void kerr_geo_polar_mino_frequency_sf(double &upTh, double a, double En, double z1, double z2){
	upTh = 0.5*M_PI*sqrt(a*a*(1. - En*En))*z1/elliptic_k(z2/z1);
}

double rp_of_angle_sf(double qr, double r1, double r2, double r3, double r4){
	return (r3*(r1 - r2)*pow(jacobi_sn(elliptic_k(sqrt((r1-r2)/(r1-r3)*(r3-r4)/(r2-r4)))/M_PI*qr, sqrt((r1-r2)/(r1-r3)*(r3-r4)/(r2-r4))), 2)-r2*(r1-r3))/((r1-r2)*pow(jacobi_sn(elliptic_k(sqrt((r1-r2)/(r1-r3)*(r3-r4)/(r2-r4)))/M_PI*qr, sqrt((r1-r2)/(r1-r3)*(r3-r4)/(r2-r4))), 2)-(r1-r3));
}

double zp_of_angle_sf(double qth, double z1, double z2){
	return -z2*jacobi_sn(2.*elliptic_k(z2/z1)/M_PI*(qth + 1.5*M_PI), z2/z1);
}

/////////////////////
// Circular Orbits //
/////////////////////

GeodesicSource kerr_geo_circ(double a, double r, int sgnX){
	double En = kerr_geo_energy_circ(sgnX*a, r);
	double Lz = kerr_geo_momentum_circ(sgnX*a, r);
	double upT = kerr_geo_time_frequency_circ(sgnX*a, r);
	double upPhi = sgnX*kerr_geo_azimuthal_frequency_circ(sgnX*a, r);

	Vector tp(1);
	tp[0] = 0.;
	Vector phip(1);
	phip[0] = 0.;
	Vector rp(1);
	rp[0] = r;
	Vector thetap(1);
	thetap[0] = M_PI/2.;

	GeodesicSource geo;
	double r1 = r, r2 = r, r3 = 2./(1. - En*En) - 2.*r, r4 = 0., z1 = 0., z2 = 0.;
	geo.setConstants(a, r, 0., sgnX, En, Lz, 0., r1, r2, r3, r4, z1, z2, upT, 0., 0., upPhi, 0., 0., 0.);
	geo.setTrajectory(tp, tp, rp, thetap, phip, phip);

	return geo;
}

double kerr_geo_energy_circ(double a, double r){
	double v = 1./sqrt(r);
	return (1. - 2.*pow(v, 2) + a*pow(v, 3))/sqrt(1. - 3.*pow(v, 2) + 2.*a*pow(v, 3));
}
double kerr_geo_momentum_circ(double a, double r){
	double v = 1./sqrt(r);
	if(a < 0){
		return -v*r*(1. - 2.*a*pow(v, 3) + pow(a, 2)*pow(v, 4))/sqrt(1. - 3.*pow(v, 2) + 2.*a*pow(v, 3));
	}
	return v*r*(1. - 2.*a*pow(v, 3) + pow(a, 2)*pow(v, 4))/sqrt(1. - 3.*pow(v, 2) + 2.*a*pow(v, 3));
}
double kerr_geo_time_frequency_circ(double a, double r){
	double v = 1./sqrt(r);
	return pow(r, 2)*(1 + a*pow(v, 3))/sqrt(1. - 3.*pow(v, 2) + 2.*a*pow(v, 3));
}
double kerr_geo_azimuthal_frequency_circ(double a, double r){
	double v = 1./sqrt(r);
	// if(a < 0){
	// 	return -pow(r, 2)*pow(v, 3)/sqrt(1. - 3.*pow(v, 2) + 2.*a*pow(v, 3));
	// }
	return pow(r, 2)*pow(v, 3)/sqrt(1. - 3.*pow(v, 2) + 2.*a*pow(v, 3));
}

double kerr_geo_azimuthal_frequency_circ_time(double a, double r){
	double v = 1./sqrt(r);
	// if(a < 0){
	// 	return -pow(v, 3)/(1 + a*pow(v, 3));
	// }
	return pow(v, 3)/(1 + a*pow(v, 3));
}
double kerr_geo_azimuthal_frequency_circ_time(double a, double r, int sgnX){
	return sgnX*kerr_geo_azimuthal_frequency_circ_time(sgnX*a, r);
}
double kerr_geo_radius_circ(double a, double Omega){
	double sgnOmega = Omega/std::abs(Omega);
	return pow((1. - a*Omega)/(sgnOmega*Omega), 2./3.);
}

double kerr_geo_denergy_domega_circ(double a, double om){
	return denergy_dr(a, kerr_geo_radius_circ(a, om))*dr_domega(a, om);
}
double kerr_geo_dmomentum_domega_circ(double a, double om){
	return dmomentum_dr(a, kerr_geo_radius_circ(a, om))*dr_domega(a, om);
}

double dr_domega(double a, double om){
	return -2./(3.*pow(om, 5./3.)*pow(1. - a*om, 1./3.));
}

double denergy_dr(double a, double r){
	double v = 1./sqrt(r);
	return 0.5*(pow(v, 4) - 6.*pow(v, 6) + 8.*a*pow(v, 7) - 3.*a*a*pow(v, 8))/pow(1. + v*v*(2.*a*v - 3.), 1.5);
}

double dmomentum_dr(double a, double r){
	double v = 1./sqrt(r);
	return 0.5*v*(1 + a*pow(v, 3))*(1. - 6.*pow(v,2) + 8.*a*pow(v,3) - 3.*pow(a,2)*pow(v,4))/pow(1. - 3.*pow(v, 2) + 2.*a*pow(v, 3), 1.5);
}

double kerr_geo_VtR(double  a, double  En, double  Lz, double , double  r){
	return En*(pow(pow(r, 2) + pow(a, 2), 2))/(pow(r, 2) - 2.*r + pow(a, 2)) + a*Lz*(1. - (pow(r, 2)+ pow(a, 2))/(pow(r, 2) - 2.*r + pow(a, 2)));
}

double kerr_geo_VtTheta(double  a, double  En, double , double , double  theta){
	return -En*pow(a*sin(theta), 2);
}

double kerr_geo_Vr(double  a, double  En, double  Lz, double  Qc, double  r){
	return -(((-2. + r)*r + pow(a, 2))*(Qc + pow(-(a*En) + Lz, 2) + pow(r, 2))) + pow(a*Lz - En*(pow(a, 2) + pow(r, 2)), 2);
}

double kerr_geo_Vtheta(double  a, double  En, double  Lz, double  Qc, double  theta){
	return Qc + pow(a, 2)*(-1. + pow(En, 2))*pow(cos(theta), 2) - pow(Lz, 2)*pow(cos(theta)/sin(theta), 2);
}

double kerr_geo_VphiR(double  a, double  En, double  Lz, double , double  r){
	return -Lz*pow(a, 2)/(pow(r, 2) - 2.*r + pow(a, 2)) - a*En*(1. - (pow(r, 2)+ pow(a, 2))/(pow(r, 2) - 2.*r + pow(a, 2)));
}

double kerr_geo_VphiTheta(double , double , double  Lz, double , double  theta){
	return Lz*pow(sin(theta), -2);
}

double kerr_geo_Vz(double  a, double  En, double  Lz, double  Qc, double  z){
	return Qc - (Qc + pow(a,2)*(1. - pow(En,2)) + pow(Lz,2))*pow(z,2) + pow(a,2)*(1. - pow(En,2))*pow(z,4);
}

double kerr_geo_Vz_dz(double  a, double  En, double  Lz, double  Qc, double  z){
	return -2.*z*(Qc + pow(a,2)*(1. - pow(En,2)) + pow(Lz,2)) + 4.*pow(a,2)*(1. - pow(En,2))*pow(z,3);
}

double kerr_geo_Vz_dz2(double  a, double  En, double  Lz, double  Qc, double  z){
	return -2.*(Qc + pow(a,2)*(1. - pow(En,2)) + pow(Lz,2)) + 12.*pow(a,2)*(1. - pow(En,2))*pow(z,2);
}

double kerr_isco(double a, int sgnX){
	double z1 = 1 + pow(1 - a*a, 1./3.)*(pow(1 - a, 1./3.) + pow(1 + a, 1./3.));
	double z2 = sqrt(3*a*a + z1*z1);

	return 3 + z2 - sgnX*sqrt((3. - z1)*(3. + z1 + 2.*z2));
}

double kerr_isco_frequency(double a){
	int sgnX = int(a/std::abs(a));
	return kerr_geo_azimuthal_frequency_circ_time(std::abs(a), kerr_isco(std::abs(a), sgnX), sgnX);
}
