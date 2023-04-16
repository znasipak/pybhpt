// inspiral.cpp

#include "inspiral.hpp"

#define ALPHA_MIN 0.
#define ALPHA_MAX 1.
#define A_MAX 0.9999
#define OMEGA_MIN 2.e-3

// AdiabaticSpline::AdiabaticSpline(AdiabaticMode mode): A(mode.A.begin(), mode.A.end(), mode.xmin, mode.deltaX),
// 	Phi(mode.Phi.begin(), mode.Phi.end(), mode.xmin, mode.deltaX), EdotUp(mode.EdotUp.begin(), mode.EdotUp.end(), mode.xmin, mode.deltaX),
// 	EdotIn(mode.EdotIn.begin(), mode.EdotIn.end(), mode.xmin, mode.deltaX), _risco(mode.risco){}
// AdiabaticSpline::~AdiabaticSpline(){}
//
// Complex AdiabaticSpline::getAmplitude(double r){ return A(log_rescale(_risco, r))*exp(I*Phi(log_rescale(_risco, r))); }
// double AdiabaticSpline::getEnergyFlux(BoundaryCondition bc, double r){ if(bc == In) return EdotIn(log_rescale(_risco, r)); else return EdotUp(log_rescale(_risco, r)); }
// double AdiabaticSpline::getEnergyFlux(double r){ return EdotIn(log_rescale(_risco, r)) + EdotUp(log_rescale(_risco, r)); }

// for quasi-circular inspirals, we allow a to vary between -1 and 1 while frequencies are positive
// this is different from the geodesic module in which a is positive and frequencies can be negative

double a_of_chi_subfunc(double chi){
	return 1. - pow(chi, 3);
}
double chi_of_a_subfunc(double a){
	return pow(1. - a, 1./3.);
}

double a_of_chi(double chi){
	return 1. - pow(chi_of_a_subfunc(A_MAX) + pow(chi, 2)*(chi_of_a_subfunc(-A_MAX) - chi_of_a_subfunc(A_MAX)), 3);
}

double chi_of_a(double a){
	return pow((chi_of_a_subfunc(a) - chi_of_a_subfunc(A_MAX))/(chi_of_a_subfunc(-A_MAX) - chi_of_a_subfunc(A_MAX)), 0.5);
}
//
// double alpha_of_a_omega(double a, double omega, double oISCO){
// 	return pow(omega*(1. - a*omega)/oISCO/(1. - a*oISCO), 1./3.);
// }

// double a_of_chi(double chi){
// 	return A_MAX - 2.*A_MAX*pow(chi, 3);
// }
//
// double chi_of_a(double a){
// 	return pow((A_MAX - a)/(2.*A_MAX), 1./3.);
// }

double alpha_of_a_omega(double, double omega, double oISCO){
	if(abs(oISCO - omega) < 1.e-13){return 0.;}
	return pow(abs(pow(oISCO, 1./3.) - pow(omega, 1./3.))/(pow(oISCO, 1./3.) - pow(OMEGA_MIN, 1./3.)), 0.5);
}

double alpha_of_a_omega(double a, double omega){
	return alpha_of_a_omega(a, omega, abs(kerr_isco_frequency(a)));
}

// double omega_of_a_alpha(double a, double alpha, double oISCO){
// 	if(abs(a) < DBL_EPSILON){
// 		return pow(alpha, 3)*oISCO;
// 	}
// 	double beta = pow(alpha, 3)*oISCO*(1. - a*oISCO);
// 	return 0.5*(1. - sqrt(1. - 4.*a*beta))/a;
// }

double omega_of_a_alpha(double, double alpha, double oISCO){
	return pow(pow(oISCO, 1./3.) - pow(alpha, 2.)*(pow(oISCO, 1./3.) - pow(OMEGA_MIN, 1./3.)), 3.);
}

double omega_of_a_alpha(double a, double alpha){
	return omega_of_a_alpha(a, alpha, abs(kerr_isco_frequency(a)));
}

double omega_of_chi_alpha(double chi, double alpha, double oISCO){
	return omega_of_a_alpha(a_of_chi(chi), alpha, oISCO);
}

double omega_of_chi_alpha(double chi, double alpha){
	return omega_of_chi_alpha(chi, alpha, abs(kerr_isco_frequency(a_of_chi(chi))));
}

Vector radius_sample_a_alpha(double a, double alphaMin, double alphaMax, int sampleNum){
	double deltaAlpha = (alphaMax - alphaMin)/double(sampleNum - 1);
	double oISCO = abs(kerr_isco_frequency(a));

	Vector rPts(sampleNum);
	for(int i = 0; i < sampleNum; i++){
		rPts[i] = kerr_geo_radius_circ(a, omega_of_a_alpha(a, alphaMin + deltaAlpha*i, oISCO));
	}

	return rPts;
}

double radius_of_a_alpha(double a, double alpha, double oISCO){
	return kerr_geo_radius_circ(a, omega_of_a_alpha(a, alpha, oISCO));
}

double radius_of_a_alpha(double a, double alpha){
	return radius_of_a_alpha(a, alpha, abs(kerr_isco_frequency(a)));
}

Vector bh_spin_sample_chi(double amax, int sampleNum){
	double chiMax = chi_of_a(-amax);
	double chiMin = chi_of_a(amax);
	double deltaChi = (chiMax - chiMin)/(sampleNum - 1);

	Vector aPts(sampleNum);
	for(int i = 0; i < sampleNum; i++){
		aPts[i] = a_of_chi(chiMin + i*deltaChi);
	}

	return aPts;
}

int generate_adiabatic_circ_data(double a, int sgnX, int Lmax, const std::string& dir){
	if(!boost::filesystem::exists(dir)){
		boost::filesystem::create_directory(dir);
	}

	for(int L = 2; L <= Lmax; L++){
		for(int m = 1; m <= L; m++){
			std::cout << "INSPIRAL: Generating data for L = " << L << " m = "<< m <<" \n";
			generate_adiabatic_circ_mode_data(a, sgnX, L, m, dir);
			//std::cout << "INSPIRAL: Generating data for L = " << L << " m = "<< -m <<" \n";
			//generate_adiabatic_circ_mode_data(a, sgnX, L, -m, dir);
		}
	}

	return 0;
}

int mode_count_total(int Lmax){
	return Lmax*(Lmax + 1)/2 - 1;
}

int l_mode_from_total_mode_count(int modeNum){
	int l = 2;
	int modeCount = l;
	while(modeCount < modeNum){
		l++;
		modeCount += l;
	}
	return l;
}

int m_mode_from_total_mode_count(int modeNum){
	int l = l_mode_from_total_mode_count(modeNum);
	return modeNum - mode_count_total(l - 1);
}

int generate_adiabatic_circ_data_a_sample_parallel(int lMax, const std::string& dir){
	int modeCount = mode_count_total(lMax);
	int l, m, i;
	int sampleSize = 201;

	#pragma omp parallel shared(modeCount, dir) private(l, m, i)
	{
		#pragma omp for
			for(i = 0; i < modeCount; i++){
				l = l_mode_from_total_mode_count(i + 1);
				m = m_mode_from_total_mode_count(i + 1);
				std::cout << "Starting ("<<l<<", "<<m<<")-mode \n";
				generate_adiabatic_circ_mode_data_a_sample(l, m, sampleSize, dir);
			}
	}

	generate_adiabatic_circ_remainder_data_a_sample(lMax + 1, sampleSize, dir);

	return 0;
}

int generate_adiabatic_circ_mode_data_a_sample(int L, int m, int sampleSize, const std::string& dir){
	if(!boost::filesystem::exists(dir)){
		boost::filesystem::create_directory(dir);
	}

	std::string filename = "circ_data_alpha99_" + std::to_string(L) + "_" + std::to_string(m) + ".txt";
	std::string filepath = dir + "/" + filename;
	std::ofstream file;
	file.open(filepath);

	double alpha = 0.99;
	Vector aPts = bh_spin_sample_chi(0.999, sampleSize);
	file << "a\tr0\tomega\tAlm\tPhilm\tEdotUp\tEdotIn\n";

	for(size_t i = 0; i < aPts.size(); i++){
		double a = aPts[i];
		double r0 = radius_of_a_alpha(a, alpha);
		int sgnX = int(a/abs(a));
		GeodesicSource geo = kerr_geo_circ(abs(a), r0, sgnX);
		TeukolskyMode teukMode(-2, L, m, 0, 0, geo);
		Fluxes fluxMode = gravitational_energy_flux_mode(geo, teukMode);
		Complex amplitude = -teukMode.getTeukolskyAmplitude(Up)*pow(teukMode.getFrequency(), -2);
		double Alm = 0.;
		double Philm = 0.;
		double fluxUpOut = 0.;
		double fluxInOut = 0.;
		if( !isnan(abs(amplitude)) ){
			Alm = std::abs(amplitude);
			Philm = std::arg(amplitude);
			fluxUpOut = fluxMode.infinity;
			fluxInOut = fluxMode.horizon;
		}
		file << std::setprecision(15);
		file << a << "\t" << r0 << "\t" << teukMode.getFrequency() << "\t" << Alm << "\t" <<
			Philm  << "\t" << fluxUpOut << "\t" << fluxInOut << "\n";
	}

	file.close();

	return 0;
}

Fluxes energy_flux_remainder(int Lmin, GeodesicSource &geo){
	double EdotUp = 0.;
	double EdotIn = 0.;
	int l = Lmin;

	double EdotNorm = energy_flux_newtonian(-2, geo);
	Fluxes Edot = energy_flux_l(-2, Lmin, geo);
	double error = abs(1. - (Edot.infinity + Edot.horizon)/EdotNorm);
	EdotUp += Edot.infinity;
	EdotIn += Edot.horizon;
	l++;

	while(error > 1.e-9 && l < 60){
		Edot = energy_flux_l(-2, l, geo);
		error = abs((Edot.infinity + Edot.horizon)/EdotNorm);
		EdotUp += Edot.infinity;
		EdotIn += Edot.horizon;
		l++;
	}

	Edot.infinity = EdotUp;
	Edot.horizon = EdotIn;

	return Edot;
}

int generate_adiabatic_circ_remainder_data_a_sample(int Lmin, int sampleSize, const std::string& dir){
	if(!boost::filesystem::exists(dir)){
		boost::filesystem::create_directory(dir);
	}

	std::string filename = "circ_data_alpha99_" + std::to_string(Lmin) + ".txt";
	std::string filepath = dir + "/" + filename;
	std::ofstream file;
	file.open(filepath);

	double alpha = 0.99;
	Vector aPts = bh_spin_sample_chi(0.999, sampleSize);
	file << "a\tr0\tomega\tAlm\tPhilm\tEdotUp\tEdotIn\n";

	for(size_t i = 0; i < aPts.size(); i++){
		double a = aPts[i];
		double r0 = radius_of_a_alpha(a, alpha);
		int sgnX = int(a/abs(a));
		GeodesicSource geo = kerr_geo_circ(abs(a), r0, sgnX);
		Fluxes flux = energy_flux_remainder(Lmin, geo);

		double Alm = 0.;
		double Philm = 0.;
		double fluxUpOut = flux.infinity;
		double fluxInOut = flux.horizon;
		file << std::setprecision(15);
		file << a << "\t" << r0 << "\t" << geo.getTimeFrequency(3) << "\t" << Alm << "\t" <<
			Philm  << "\t" << fluxUpOut << "\t" << fluxInOut << "\n";
	}

	file.close();

	return 0;
}

int generate_adiabatic_circ_mode_data(double a, int sgnX, int L, int m, const std::string& dir){
	std::string filename = "circ_data_" + std::to_string(L) + "_" + std::to_string(m) + ".txt";
	std::string aDir;
	if(sgnX > 0){
		aDir = dir + "/a" + std::to_string(int(a*1000));
	}else{
		aDir = dir + "/an" + std::to_string(int(a*1000));
	}
	if(!boost::filesystem::exists(aDir)){
		boost::filesystem::create_directory(aDir);
	}


	std::string filepath = aDir + "/" + filename;
	if(!boost::filesystem::exists(filepath)){
		std::ofstream file;
		file.open(filepath);

		double risco = kerr_isco(a, sgnX);
		file << "L\tm\ta\tx\tr_isco\n";
		file << L << "\t" << m << "\t" << a << "\t" << sgnX << "\t" << risco << "\n\n";
		file << "r0\tx0\tomega\tAlm\tPhilm\tEdotUp\tEdotIn\n";

		// Vector rPts = radial_log_sample(risco);
		double rmin = risco + 0.005;
		double rmax = 30;
		Vector rPts = radial_frequency_sample(a, rmin, rmax, sgnX, 1001);
		for(size_t i = 0; i < rPts.size(); i++){
			generate_adiabatic_circ_mode_data_r0(a, sgnX, risco, rPts[i], L, m, file);
		}

		file.close();
	}

	return 0;
}

int generate_adiabatic_circ_mode_data_r0(double a, int sgnX, double risco, double r0, int L, int m, std::ofstream& file){
	double x0 = log_rescale(risco, r0);
	GeodesicSource geo = kerr_geo_circ(a, r0, sgnX);

	TeukolskyMode teukMode(L, m, 0, 0, geo);
	Fluxes fluxMode = gravitational_energy_flux_mode(geo, teukMode);
	Complex amplitude = -teukMode.getTeukolskyAmplitude(Up)*pow(teukMode.getFrequency(), -2);
	double Alm = 0.;
	double Philm = 0.;
	double fluxUpOut = 0.;
	double fluxInOut = 0.;
	if( !isnan(abs(amplitude)) ){
		Alm = std::abs(amplitude);
		Philm = std::arg(amplitude);
		fluxUpOut = fluxMode.infinity;
		fluxInOut = fluxMode.horizon;
	}
	file << std::setprecision(15);
	file << r0 << "\t" << x0 << "\t" << teukMode.getFrequency() << "\t" << Alm << "\t" <<
		Philm  << "\t" << fluxUpOut << "\t" << fluxInOut << "\n";

	return 0;
}

int generate_adiabatic_circ_data_radial_sample_parallel(double a, int lMax, int sampleSize, const std::string& dir){
	int modeCount = mode_count_total(lMax);
	int l, m, i;
	if(!boost::filesystem::exists(dir)){
		boost::filesystem::create_directory(dir);
	}

	#pragma omp parallel shared(modeCount, dir) private(l, m, i)
	{
		#pragma omp for
			for(i = 0; i < modeCount; i++){
				l = l_mode_from_total_mode_count(i + 1);
				m = m_mode_from_total_mode_count(i + 1);
				std::cout << "Starting ("<<l<<", "<<m<<")-mode \n";
				generate_adiabatic_circ_mode_data_radial_sample(a, l, m, sampleSize, dir);
			}
	}

	generate_adiabatic_circ_remainder_data_radial_sample(a, lMax + 1, sampleSize, dir);

	return 0;
}

int generate_adiabatic_circ_mode_data_radial_sample(double a, int L, int m, int sampleNum, const std::string& dir){
	std::string filename = "circ_data_" + std::to_string(L) + "_" + std::to_string(m) + ".txt";

	std::string filepath = dir + "/" + filename;
	if(!boost::filesystem::exists(filepath)){
		std::ofstream file;
		file.open(filepath);

		file << "a\tr0\tomega\tchi\talpha\tAlm\tPhilm\tEdotUp\tEdotIn\n";

		double alphaMin = ALPHA_MIN;
		double alphaMax = ALPHA_MAX;
		Vector rPts = radius_sample_a_alpha(a, alphaMin, alphaMax, sampleNum);
		for(size_t i = 0; i < rPts.size(); i++){
			generate_adiabatic_circ_mode_data_a_r0(a, rPts[i], L, m, file);
		}

		file.close();
	}

	return 0;
}

int generate_adiabatic_circ_remainder_data_radial_sample(double a, int Lmin, int sampleNum, const std::string& dir){
	std::string filename = "circ_data_" + std::to_string(Lmin) + ".txt";

	std::string filepath = dir + "/" + filename;
	std::ofstream file;
	file.open(filepath);

	file << "a\tr0\tomega\tchi\talpha\tAlm\tPhilm\tEdotUp\tEdotIn\n";

	double alphaMin = ALPHA_MIN;
	double alphaMax = ALPHA_MAX;
	Vector rPts = radius_sample_a_alpha(a, alphaMin, alphaMax, sampleNum);
	for(size_t i = 0; i < rPts.size(); i++){
		generate_adiabatic_circ_remainder_data_a_r0(a, rPts[i], Lmin, file);
	}

	file.close();

	return 0;
}

int generate_adiabatic_circ_mode_data_a_r0(double a, double r0, int L, int m, std::ofstream& file){
	int sgnX = sgn(a);
	if(sgnX == 0){
		sgnX = 1;
	}
	GeodesicSource geo = kerr_geo_circ(abs(a), r0, sgnX);

	TeukolskyMode teukMode(L, m, 0, 0, geo);
	Fluxes fluxMode = gravitational_energy_flux_mode(geo, teukMode);
	Complex amplitude = -2.*teukMode.getTeukolskyAmplitude(Up)*pow(teukMode.getFrequency(), -2);
	double Alm = 0.;
	double Philm = 0.;
	double fluxUpOut = 0.;
	double fluxInOut = 0.;
	if( !isnan(abs(amplitude)) ){
		Alm = std::abs(amplitude);
		Philm = std::arg(amplitude);
		fluxUpOut = fluxMode.infinity;
		fluxInOut = fluxMode.horizon;
	}
	double omega = abs(geo.getTimeFrequency(3));
	double alpha = alpha_of_a_omega(a, omega);
	double chi = chi_of_a(a);

	file << std::setprecision(15);
	file << a << "\t" << r0 << "\t" << omega << "\t" << chi << "\t" << alpha << "\t" << Alm << "\t" <<
		Philm  << "\t" << fluxUpOut << "\t" << fluxInOut << "\n";

	return 0;
}

int output_adiabatic_circ_spherical_mode_data_radial_sample_parallel(double a, int lmax, int sampleNum, const std::string& dir){
	int m;
	if(!boost::filesystem::exists(dir)){
		boost::filesystem::create_directory(dir);
	}

	std::cout << dir << "\n";

	#pragma omp parallel shared(a, lmax, sampleNum, dir) private(m)
	{
		#pragma omp for
			for(m = 1; m <= lmax; m++){
				std::cout << "Starting ("<<m<<")-mode \n";
				output_adiabatic_circ_spherical_mode_data_radial_sample(a, lmax, m, sampleNum, dir);
				std::cout << "("<<m<<")-mode complete \n";
			}
	}

	return 0;
}

int output_adiabatic_circ_spherical_mode_data_radial_sample(double a, int lmax, int sampleNum, const std::string& dir){
	if(!boost::filesystem::exists(dir)){
		boost::filesystem::create_directory(dir);
	}

	std::cout << dir << "\n";

	for(int m = 1; m <= lmax; m++){
		std::cout << "Starting ("<<m<<")-mode \n";
		output_adiabatic_circ_spherical_mode_data_radial_sample(a, lmax, m, sampleNum, dir);
		std::cout << "("<<m<<")-mode complete \n";
	}

	return 0;
}

int output_adiabatic_circ_spherical_mode_data_radial_sample(double a, int lmax, int m, int sampleNum, const std::string& dir){
	if(!boost::filesystem::exists(dir)){
		boost::filesystem::create_directory(dir);
	}

	int lmin = abs(m) < 2 ? 2 : abs(m);
	ComplexMatrix amplitude(sampleNum, ComplexVector(lmax - lmin + 1));
	ComplexMatrix flux(sampleNum, ComplexVector(lmax - lmin + 1));
	generate_adiabatic_circ_spherical_mode_data_radial_sample(amplitude, flux, a, m);
	double alphaMin = ALPHA_MIN;
	double alphaMax = ALPHA_MAX;
	Vector rPts = radius_sample_a_alpha(a, alphaMin, alphaMax, sampleNum);

	for(int l = lmin; l <= lmax; l++){
		std::string filename = "circ_data_" + std::to_string(l) + "_" + std::to_string(m) + ".txt";
		std::string filepath;
		if(dir.back() == '/'){
			filepath = dir + filename;
		}else{
			filepath = dir + "/" + filename;
		}
		std::ofstream file;
		file.open(filepath);

		file << "a\tr0\tomega\tchi\talpha\tAlm\tPhilm\tEdotUp\tEdotIn\n";
		for(size_t i = 0; i < amplitude.size(); i++){
			if(l - lmin >= int(amplitude[i].size())){
				std::cout << "No index for l = " << l << "\n";
			}
			double Alm = abs(amplitude[i][l - lmin]);
			double Philm = arg(amplitude[i][l - lmin]);
			double fluxUpOut = std::imag(flux[i][l - lmin]);
			double fluxInOut = std::real(flux[i][l - lmin]);
			double r0 = rPts[i];
			double omega = abs(kerr_geo_azimuthal_frequency_circ_time(a, r0));
			double alpha = alpha_of_a_omega(a, omega);
			double chi = chi_of_a(a);
			// if(i==0){std::cout << a << " " << omega << " " << alpha << "\n";}

			file << std::setprecision(14);
			file << a << "\t" << rPts[i] << "\t" << omega << "\t" << chi << "\t" << alpha << "\t" << Alm << "\t" <<
				Philm  << "\t" << fluxUpOut << "\t" << fluxInOut << "\n";
		}

		file.close();
	}

	return 0;
}

int generate_adiabatic_circ_spherical_mode_data_radial_sample(ComplexMatrix &amplitude, ComplexMatrix &flux, double a, int m){
	int sampleNum = amplitude.size();
	double alphaMin = ALPHA_MIN;
	double alphaMax = ALPHA_MAX;
	Vector rPts = radius_sample_a_alpha(a, alphaMin, alphaMax, sampleNum);
	for(size_t i = 0; i < rPts.size(); i++){
		generate_adiabatic_circ_spherical_mode_data_a_r0(amplitude[i], flux[i], a, rPts[i], m);
	}

	return 0;
}

int generate_adiabatic_circ_spherical_mode_data_a_r0(ComplexVector &amplitude, ComplexVector &flux, double a, double r0, int m){
	int sgnX = sgn(a);
	if(sgnX == 0){
		sgnX = 1;
	}
	GeodesicSource geo = kerr_geo_circ(abs(a), r0, sgnX);

	int lmin = abs(m) < 2 ? 2 : abs(m);
	int lmax = amplitude.size() - 1 + lmin;
	int L = lmin;
	TeukolskyMode teukMode(-2, L, m, 0, 0, geo);
	teukMode.generateSolutions(geo);

	int couplingMin = L;
	while(teukMode.getCouplingCoefficient(couplingMin) > 1.e-10 && couplingMin >= lmin){
		couplingMin -= 1;
	}

	// while the spheroidal L mode is below 50, but the coupling between the spheroidal L mode and the max spherical l mode is strong (> 1.e-10)
	while(L < 50 && couplingMin <= lmax){
		// first calculate complex amplitude of h(t) on a spheroidal basis
		Complex ampContribution = -2.*teukMode.getTeukolskyAmplitude(Up)*pow(teukMode.getFrequency(), -2);
		for(int l = lmin; l <= lmax; l++){
			if( !isnan(abs(ampContribution)) ){
				// divide spheroidal contributions on a spherical basis
				// std::cout << "b^" << l << "_" << L << " = " << teukMode.getCouplingCoefficient(l) << "\n";
				amplitude[l - lmin] += teukMode.getCouplingCoefficient(l)*ampContribution;
			}
		}
		// calculate the flux on a spheroidal basis
		if(L <= lmax){
			Fluxes fluxMode = gravitational_energy_flux_mode(geo, teukMode);
			flux[L - lmin] = fluxMode.horizon + I*fluxMode.infinity;
		}

		// move onto the next spheroidal mode
		L += 1;
		teukMode = TeukolskyMode(-2, L, m, 0, 0, geo);
		teukMode.generateSolutions(geo);
		couplingMin = L;
		while(abs(teukMode.getCouplingCoefficient(couplingMin)) > 1.e-10 && couplingMin >= lmin){
			couplingMin -= 1;
		}
	}

	return 0;
}

int generate_adiabatic_circ_remainder_data_a_r0(double a, double r0, int Lmin, std::ofstream& file){
	int sgnX = sgn(a);
	if(sgnX == 0){
		sgnX = 1;
	}
	GeodesicSource geo = kerr_geo_circ(abs(a), r0, sgnX);
	Fluxes flux = energy_flux_remainder(Lmin, geo);

	double Alm = 0.;
	double Philm = 0.;
	double fluxUpOut = flux.infinity;
	double fluxInOut = flux.horizon;

	double omega = abs(geo.getTimeFrequency(3));
	double alpha = alpha_of_a_omega(a, abs(omega));
	double chi = chi_of_a(a);

	file << std::setprecision(15);
	file << a << "\t" << r0 << "\t" << omega << "\t" << chi << "\t" << alpha << "\t" << Alm << "\t" <<
		Philm  << "\t" << fluxUpOut << "\t" << fluxInOut << "\n";

	return 0;
}

Vector radial_frequency_sample(double a, double rmin, double rmax, double sgnX, int sampleNum){
	double omegaMax = kerr_geo_azimuthal_frequency_circ_time(a, rmin, sgnX);
	double omegaMin = kerr_geo_azimuthal_frequency_circ_time(a, rmax, sgnX);
	double deltaOmega = (omegaMax - omegaMin)/(sampleNum - 1);

	Vector rPts(sampleNum);
	for(int i = 0; i < sampleNum; i++){
		rPts[i] = kerr_geo_radius_circ(a, omegaMin + i*deltaOmega);
	}

	return rPts;
}

double log_rescale(double r0, double r){
	return log((r + r0)/(r - r0));
}

double exp_rescale(double r0, double x){
	return r0*(1. + exp(-x))/(1. - exp(-x));
}

AdiabaticMode read_adiabatic_circ_mode(int L, int m, const std::string& dir){
	std::string filename = "/circ_data_" + std::to_string(L) + "_" + std::to_string(m) + ".txt";
	if(m == 0){
		filename = "/circ_data_" + std::to_string(L) + ".txt";
	}
	std::string filepath = dir + filename;

	double spin, chi, r0, freq, alpha, Alm, Philm, EdotUp, EdotIn;
	Vector rVec, alphaVec, freqVec, AlmVec, PhilmVec, EdotUpVec, EdotInVec;
	std::istringstream lin;
	std::ifstream inFile(filepath);
	for (std::string line; std::getline(inFile, line); ) {
	    lin.clear();
	    lin.str(line);
	    if(lin >> spin >> r0 >> freq >> chi >> alpha >> Alm >> Philm >> EdotUp >> EdotIn){
	        rVec.push_back(r0);
					alphaVec.push_back(alpha);
					freqVec.push_back(freq);
					AlmVec.push_back(Alm);
					PhilmVec.push_back(Philm);
					EdotUpVec.push_back(EdotUp);
					EdotInVec.push_back(EdotIn);
					// if(isnan(EdotUp)){
					// 	std::cout << "(inspiral) EdotUp is NaN for ("<<L<<", "<<m<<")-mode.\n";
					// }
					// if(isnan(EdotIn)){
					// 	std::cout << "(inspiral) EdotUp is NaN for ("<<L<<", "<<m<<")-mode.\n";
					// }
	    }
	}

	double freqISCO = abs(kerr_isco_frequency(spin));

	AdiabaticMode mode = {
		.spin = spin,
		.chi = chi,
		.freqISCO = freqISCO,
		.r0 = rVec,
		.alpha = alphaVec,
		.freq = freqVec,
		.A = AlmVec,
		.Phi = PhilmVec,
		.EdotUp = EdotUpVec,
		.EdotIn = EdotInVec
	};

	return mode;
}

// double dalpha_domega_of_a_omega(const double &a, const double &omega, const double &oISCO){
// 	return (1. - 2.*a*omega)*pow(oISCO*(1. - a*oISCO)*pow(omega*(1. - a*omega), 2), -1./3.)/3.;
// }

// double dalpha_domega_of_a_omega(const double &a, const double &omega, const double &oISCO){
// 	double alpha = alpha_of_a_omega(a, omega, oISCO);
// 	return -2.*alpha/(3.*(oISCO - omega));
// }

// double domega_dalpha_of_a_omega(const double &a, const double &omega, const double &oISCO){
// 	return -1.5*(oISCO - OMEGA_MIN)*sqrt(alpha_of_a_omega(a, omega, oISCO));
// }

double domega_dalpha_of_a_omega(const double &, const double &omega, const double &oISCO){
	if(abs(oISCO - omega) < 1.e-13){return 0.;}
	return -6.*pow((pow(oISCO, 1./3.) - pow(OMEGA_MIN, 1./3.))*(pow(oISCO, 1./3.) - pow(omega, 1./3.)), 0.5)*pow(omega, 2./3.);
}

double energy_flux_of_omega_spline(double omega, SplineData spline){
	return gsl_spline_eval(spline.spline, alpha_of_a_omega(spline.a, omega, spline.rescale), spline.acc)*energy_flux_newtonian(-2, omega);
}

double energy_flux_of_alpha_spline(double alpha, SplineData spline){
	return gsl_spline_eval(spline.spline, alpha, spline.acc)*energy_flux_newtonian(-2, omega_of_a_alpha(spline.a, alpha, spline.rescale));
}

int dt_dAlpha(double alpha, const double*, double dtdalpha[], void* parameters){
  SplineData spline = *(SplineData *)parameters;
	double o = omega_of_a_alpha(spline.a, alpha, spline.rescale);

	double Edot = energy_flux_of_alpha_spline(alpha, spline);
	// double Ldot = Edot/o;
	double dOmega_dE = 1./kerr_geo_denergy_domega_circ(spline.a, o);
	// double dOmega_dL = 1./kerr_geo_dmomentum_domega_circ(spline.a, o);
	// double dOmega_dalpha = 1./dalpha_domega_of_a_omega(spline.a, o, spline.rescale);
	double dOmega_dalpha = domega_dalpha_of_a_omega(spline.a, o, spline.rescale);

	// dtdalpha[0] = -dOmega_dalpha/(dOmega_dE*Edot + dOmega_dL*Ldot);
	dtdalpha[0] = -dOmega_dalpha/(dOmega_dE*Edot);

  return GSL_SUCCESS;
}

int dPhi_dAlpha(double alpha, const double*, double dPhidalpha[], void* parameters){
  SplineData spline = *(SplineData *)parameters;
	double o = omega_of_a_alpha(spline.a, alpha, spline.rescale);

	double Edot = energy_flux_of_alpha_spline(alpha, spline);
	// double Ldot = Edot/o;
	double dOmega_dE = 1./kerr_geo_denergy_domega_circ(spline.a, o);
	// double dOmega_dL = 1./kerr_geo_dmomentum_domega_circ(spline.a, o);
	// double dOmega_dalpha = 1./dalpha_domega_of_a_omega(spline.a, o, spline.rescale);
	double dOmega_dalpha = domega_dalpha_of_a_omega(spline.a, o, spline.rescale);

	// dPhidalpha[0] = -o*dOmega_dalpha/(dOmega_dE*Edot + dOmega_dL*Ldot);
	dPhidalpha[0] = -o*dOmega_dalpha/(dOmega_dE*Edot);

  return GSL_SUCCESS;
}

int dt_dOmega(double omega, const double*, double dtdomega[], void* parameters){
  SplineData spline = *(SplineData *)parameters;

	double Edot = energy_flux_of_omega_spline(omega, spline);
	double dOmega_dE = 1./kerr_geo_denergy_domega_circ(spline.a, omega);

	dtdomega[0] = -1./(dOmega_dE*Edot);

  return GSL_SUCCESS;
}

int dPhi_dOmega(double omega, const double*, double dPhidomega[], void* parameters){
  SplineData spline = *(SplineData *)parameters;

	double Edot = energy_flux_of_omega_spline(omega, spline);
	double dOmega_dE = 1./kerr_geo_denergy_domega_circ(spline.a, omega);

	dPhidomega[0] = -omega/(dOmega_dE*Edot);

  return GSL_SUCCESS;
}

int dr_dOmega(double omega, const double*, double drdomega[], void* parameters){
  SplineData spline = *(SplineData *)parameters;
	double r = kerr_geo_radius_circ(spline.a, omega);

	double Edot = energy_flux_of_omega_spline(omega, spline);
	double dE_dr = denergy_dr(spline.a, r);
	double dr_dOmega = dr_domega(spline.a, omega);

	drdomega[0] = -dr_dOmega*dE_dr*Edot;

  return GSL_SUCCESS;
}

int dr_dt(double, const double r[], double drdt[], void* parameters){
  SplineData spline = *(SplineData *)parameters;
	double omega = kerr_geo_azimuthal_frequency_circ_time(spline.a, r[0]);
	if(alpha_of_a_omega(spline.a, omega, spline.rescale) > ALPHA_MAX){
		std::cout << "shift down\n";
		omega = spline.rescale;
	}
	if(alpha_of_a_omega(spline.a, omega, spline.rescale) < ALPHA_MIN){
		std::cout << "shift up\n";
		omega = omega_of_a_alpha(spline.a, ALPHA_MIN, spline.rescale);
	}

	double Edot = energy_flux_of_omega_spline(omega, spline);
	double dr_dE = 1./denergy_dr(spline.a, r[0]);

	drdt[0] = -dr_dE*Edot;

  return GSL_SUCCESS;
}

Vector alpha_sample(double alphaMin, double alphaMax, int sampleNum){
	double deltaAlpha = (alphaMax - alphaMin)/(sampleNum - 1);

	Vector alphaPts(sampleNum);
	for(int i = 0; i < sampleNum; i++){
		alphaPts[i] = alphaMin + deltaAlpha*i;
	}

	return alphaPts;
}

int generate_adiabatic_circ_trajectory(const std::string& dir, int outputSampleSize){
	AdiabaticMode insp = read_adiabatic_circ_mode(2, 1, dir);
	double a = insp.spin;
	int sampleNum = insp.r0.size();
	double omegaISCO = insp.freqISCO;
	if(abs(a) < 1.e-6){
		a = 0.;
		omegaISCO = kerr_isco_frequency(a);
	}

	double alpha[sampleNum];
	double EdotTot[sampleNum];
	double omega;

	if(insp.alpha[0] < insp.alpha[1]){
		for(int i = 0; i < sampleNum; i++){
			alpha[i] = insp.alpha[i];
			EdotTot[i] = 0.;
		}
	}else{
		for(int i = 0; i < sampleNum; i++){
			alpha[i] = insp.alpha[sampleNum - 1 - i];
			EdotTot[i] = 0.;
		}
	}
	alpha[0] = ALPHA_MIN;
	alpha[sampleNum - 1] = ALPHA_MAX;

	std::cout << "Start summing fluxes \n";
	int lmax = 15;

	for(int l = 2; l <= lmax; l++){
		for(int m = 1; m <= l; m++){
			insp = read_adiabatic_circ_mode(l, m, dir);
			if(insp.alpha[0] < insp.alpha[1]){
				for(int i = 0; i < sampleNum; i++){
					omega = omega_of_a_alpha(a, insp.alpha[i], omegaISCO);
					EdotTot[i] += 2.*(insp.EdotIn[i] + insp.EdotUp[i])/energy_flux_newtonian(-2, omega);
					if(isnan(EdotTot[i]) || isinf(EdotTot[i])){
						std::cout << l << "\t" << m << "\n";
					}
				}
			}else{
				for(int i = 0; i < sampleNum; i++){
					omega = omega_of_a_alpha(a, insp.alpha[sampleNum - 1 - i], omegaISCO);
					EdotTot[i] += 2.*(insp.EdotIn[sampleNum - 1 - i] + insp.EdotUp[sampleNum - 1 - i])/energy_flux_newtonian(-2, omega);
					if(isnan(EdotTot[i]) || isinf(EdotTot[i])){
						std::cout << l << "\t" << m << "\n";
					}
				}
			}
		}
	}

	std::cout << "Add flux remainders \n";

	insp = read_adiabatic_circ_mode(lmax + 1, 0, dir);
	if(insp.alpha[0] < insp.alpha[1]){
		for(int i = 0; i < sampleNum; i++){
			omega = omega_of_a_alpha(a, insp.alpha[i], omegaISCO);
			EdotTot[i] += 2.*(insp.EdotIn[i] + insp.EdotUp[i])/energy_flux_newtonian(-2, omega);
		}
	}else{
		for(int i = 0; i < sampleNum; i++){
			omega = omega_of_a_alpha(a, insp.alpha[sampleNum - 1 - i], omegaISCO);
			EdotTot[i] += 2.*(insp.EdotIn[sampleNum - 1 - i] + insp.EdotUp[sampleNum - 1 - i])/energy_flux_newtonian(-2, omega);
		}
	}

	std::cout << "Finished summing fluxes \n";

	gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, sampleNum);

  gsl_spline_init(spline, alpha, EdotTot, sampleNum);

	SplineData splineParams = {
		.spline = spline,
		.acc = acc,
		.a = a,
		.rescale = omegaISCO
	};

	int sampleNumHigh = outputSampleSize;
	Vector tOfAlpha(sampleNumHigh);
	Vector PhiOfAlpha(sampleNumHigh);
	double alpha_min = alpha[0];
	double alpha_max = alpha[sampleNum - 1];
	Vector alphaSample = alpha_sample(alpha_min, alpha_max, sampleNumHigh);
	if(omega_of_a_alpha(a, alpha_min) < omega_of_a_alpha(a, alpha_max)){
		alphaSample = alpha_sample(alpha_max, alpha_min, sampleNumHigh);
	}

	double alpha0 = alphaSample[0];
	double t0 = 0.;
	double Phi0 = 0.;

	std::cout << "Integrate ODEs \n";

	// double o = omega_of_a_alpha(splineParams.a, 0.5, splineParams.rescale);
	//
	// double Edot = energy_flux_of_alpha_spline(0.5, splineParams);
	// // double Ldot = Edot/o;
	// double dOmega_dE = 1./kerr_geo_denergy_domega_circ(splineParams.a, o);
	// // double dOmega_dL = 1./kerr_geo_dmomentum_domega_circ(spline.a, o);
	// double dOmega_dalpha = 1./dalpha_domega_of_a_omega(splineParams.a, o, splineParams.rescale);

	integrate_first_order_real_ode(tOfAlpha, alphaSample, &dt_dAlpha, t0, alpha0, &splineParams);
	integrate_first_order_real_ode(PhiOfAlpha, alphaSample, &dPhi_dAlpha, Phi0, alpha0, &splineParams);

	double nanFlag = 0;
	for(int i = 1; i < sampleNumHigh; i++){
		if(tOfAlpha[i] >= tOfAlpha[i-1]){
			std::cout << "(INSPIRAL) Error: Time values are not strictly increasing at step "<< i <<" \n";
			std::cout << "(INSPIRAL) Error: t_"<<i<<" = "<<-tOfAlpha[i]<<", t_"<<i-1<<" = "<<-tOfAlpha[i-1]<<" \n";
		}
		if(isnan(tOfAlpha[i])){
			nanFlag = 1;
		}
	}
	if(nanFlag){
		std::cout << "(INSPIRAL) Error: Integration returned NaNs \n";
		std::cout << "(INSPIRAL) Error: Initial frequency = "<<omega_of_a_alpha(a, alphaSample[0], omegaISCO)<<" \n";
		std::cout << "(INSPIRAL) Error: Derivative at initial integration step = "<<domega_dalpha_of_a_omega(a, omega_of_a_alpha(a, alphaSample[0], omegaISCO), omegaISCO)<<" \n";
		std::cout << "(INSPIRAL) Error: Derivative at next alpha value = "<<domega_dalpha_of_a_omega(a, omega_of_a_alpha(a, alphaSample[1], omegaISCO), omegaISCO)<<" \n";
	}

	double beta[sampleNumHigh];
	double alphaArray[sampleNumHigh];
	for(int i = 0; i < sampleNumHigh; i++){
		beta[i] = sqrt(log(1. - tOfAlpha[i]));
		alphaArray[i] = alphaSample[i];
	}
	std::cout << "Resample frequency as a function of time \n";

	gsl_interp_accel *acc_2 = gsl_interp_accel_alloc();
  gsl_spline *spline_2 = gsl_spline_alloc(gsl_interp_cspline, sampleNumHigh);
  gsl_spline_init(spline_2, beta, alphaArray, sampleNumHigh);

	Vector betaEvenSpacing(sampleNumHigh);
	double time_before_merge = tOfAlpha.back();
	double betaMax = sqrt(log(1. - time_before_merge));
	betaEvenSpacing[0] = beta[0];
	for(int i = 1; i < sampleNumHigh - 1; i++){
		betaEvenSpacing[i] = betaMax*double(i)/double(sampleNumHigh - 1);
	}
	betaEvenSpacing[sampleNumHigh - 1] = beta[sampleNumHigh - 1];

	Vector omegaOfBeta(sampleNumHigh);
	for(int i = 0; i < sampleNumHigh; i++){
		omegaOfBeta[i] = omega_of_a_alpha(a, gsl_spline_eval(spline_2, betaEvenSpacing[i], acc_2), omegaISCO);
	}

	std::cout << "Generate dense sampling of trajectories \n";

	std::string interp_dir = "interpolation_data";
	std::string interp_dir_full = dir + "/" + interp_dir;
	std::string t_file = "time.txt";
	std::string t_path = interp_dir_full + "/" + t_file;
	std::string phase_file = "phase.txt";
	std::string phase_path = interp_dir_full + "/" + phase_file;
	std::string flux_file = "EdotNorm.txt";
	std::string flux_path = interp_dir_full + "/" + flux_file;
	std::string freq_file = "freq.txt";
	std::string freq_path = interp_dir_full + "/" + freq_file;

	if(!boost::filesystem::exists(dir + "/" + interp_dir)){
		boost::filesystem::create_directory(dir + "/" + interp_dir);
	}

	std::ofstream file;
	file.open(t_path);

	std::cout << "Save time of omega \n";
	char buff[500];
	for(int i = 0; i < sampleNumHigh; i++){
		sprintf(buff, "%.14e\t%.14e\t%.14e\n", chi_of_a(a), alphaSample[i],  abs(tOfAlpha[i]));
		file << buff;
	}

	file.close();

	file.open(phase_path);

	std::cout << "Save Phi of omega \n";

	for(int i = 0; i < sampleNumHigh; i++){
		sprintf(buff, "%.14e\t%.14e\t%.14e\n", chi_of_a(a), alphaSample[i],  abs(PhiOfAlpha[i]));
		file << buff;
	}

	file.close();

	file.open(flux_path);

	std::cout << "Save flux of omega \n";

	for(int i = 0; i < sampleNumHigh; i++){
		sprintf(buff, "%.14e\t%.14e\t%.14e\n", chi_of_a(a), alphaSample[i],  gsl_spline_eval(spline, alphaSample[i], acc));
		file << buff;
	}

	file.close();

	file.open(freq_path);

	std::cout << "Save omega of time \n";

	for(int i = 0; i < sampleNumHigh; i++){
		sprintf(buff, "%.14e\t%.14e\t%.14e\n", chi_of_a(a), betaEvenSpacing[i],  omegaOfBeta[i]);
		file << buff;
	}

	file.close();

	gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
	gsl_spline_free(spline_2);
  gsl_interp_accel_free(acc_2);

	return 0;
}

TrajectoryData read_trajectory_data(const std::string& filename){
	double chi, alpha, data;
	Vector chiVec, alphaVec, dataVec;
	std::istringstream lin;
	std::ifstream inFile(filename);
	for (std::string line; std::getline(inFile, line); ) {
	    lin.clear();
	    lin.str(line);
	    if(lin >> chi >> alpha >> data){
					chiVec.push_back(chi);
					alphaVec.push_back(alpha);
					dataVec.push_back(data);
	    }
	}

	TrajectoryData trajData = {
		.x = chiVec,
		.y = alphaVec,
		.z = dataVec
	};

	return trajData;
}

void compile_and_output_trajectory_data(std::vector<std::string> inputDirs, const std::string& outputDir){
	if(!boost::filesystem::exists(outputDir)){
		boost::filesystem::create_directory(outputDir);
	}

	std::string outFile;
	if(outputDir.back() == '/'){
		outFile = outputDir + "trajectory.txt";
	}else{
		outFile = outputDir + "/trajectory.txt";
	}

	std::ofstream file;
	file.open(outFile);

	int fileNum = inputDirs.size();
	for(int i = 0; i < fileNum; i++){
		std::string inputDir = inputDirs[i];
		std::string fileEdot, fileTime, filePhi, fileOmega;
		if(inputDir.back() == '/'){
			fileEdot = inputDir + "interpolation_data/EdotNorm.txt";
			fileTime = inputDir + "interpolation_data/time.txt";
			filePhi = inputDir + "interpolation_data/phase.txt";
			fileOmega = inputDir + "interpolation_data/freq.txt";
		}else{
			fileEdot = inputDir + "/interpolation_data/EdotNorm.txt";
			fileTime = inputDir + "/interpolation_data/time.txt";
			filePhi = inputDir + "/interpolation_data/phase.txt";
			fileOmega = inputDir + "/interpolation_data/freq.txt";
		}
		TrajectoryData Edot = read_trajectory_data(fileEdot);
		TrajectoryData time = read_trajectory_data(fileTime);
		TrajectoryData phase = read_trajectory_data(filePhi);
		TrajectoryData freq = read_trajectory_data(fileOmega);
		int sampleNum = Edot.x.size();
		if(i == 0){
			file << std::setprecision(14);
			file << "chiN\talphaN\n";
			file << fileNum << "\t" << sampleNum << "\n";
			file << "chi\talpha\tflux\ttime\tphase\tbeta\tomega\n";
		}

		char buff[500];
		for(int i = 0; i < sampleNum; i++){
			sprintf(buff, "%.14e\t%.14e\t%.14e\t%.14e\t%.14e\t%.14e\t%.14e\n", Edot.x[i], Edot.y[i], Edot.z[i], time.z[i], phase.z[i], freq.y[i], freq.z[i]);
			file << buff;
		}
	}

	file.close();
}

void compile_and_output_mode_amplitude_data(std::vector<std::string> inputDirs, const std::string& outputDir, int lmax, int radialSampleN){
	if(!boost::filesystem::exists(outputDir)){
		boost::filesystem::create_directory(outputDir);
	}

	for(int l = 2; l <= lmax; l++){
		for(int m = 1; m <= l; m++){
			std::string outFile;
			char buff[60];
			sprintf(buff, "circ_data_%d_%d.txt", l, m);
			std::string filename = buff;
			if(outputDir.back() == '/'){
				outFile = outputDir + filename;
			}else{
				outFile = outputDir + "/" + filename;
			}

			std::ofstream file;
			file.open(outFile);

			int fileNum = inputDirs.size();
			for(int i = 0; i < fileNum; i++){
				std::string inputDir = inputDirs[i];
				// std::string fileLMmode;
				// if(inputDir.back() == '/'){
				// 	fileLMmode = inputDir + filename;
				// }else{
				// 	fileLMmode = inputDir + "/" + filename;
				// }
				AdiabaticMode mode = read_adiabatic_circ_mode(l, m, inputDir);
				if(i == 0){
					file << std::setprecision(14);
					file << "chiN\talphaN\n";
					file << fileNum << "\t" << radialSampleN << "\n";
					file << "chi\talpha\tamp\tphase\n";
				}

				Vector phase = unwrap_phase(mode.Phi);
				gsl_interp_accel *acc = gsl_interp_accel_alloc();
				gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, mode.Phi.size());
				gsl_interp_accel *acc_2 = gsl_interp_accel_alloc();
				gsl_spline *spline_2 = gsl_spline_alloc(gsl_interp_cspline, mode.Phi.size());

				Vector alpha = mode.alpha;
				Vector amp = mode.A;
				if(alpha[0] > alpha[1]){
					std::reverse(alpha.begin(), alpha.end());
					std::reverse(amp.begin(), amp.end());
					std::reverse(phase.begin(), phase.end());
				}
				alpha[0] = ALPHA_MIN;
				alpha[alpha.size() - 1] = ALPHA_MAX;

				gsl_spline_init(spline, &alpha[0], &amp[0], amp.size());
				gsl_spline_init(spline_2, &alpha[0], &phase[0], phase.size());

				char buff[500];
				double delta_alpha = (ALPHA_MAX - ALPHA_MIN)/(radialSampleN - 1);
				for(int i = 0; i < radialSampleN; i++){
					double alphaPoint = ALPHA_MIN + i*delta_alpha;
					// I added this additional phase adjustment because I noticed a jump in the phase when going from a > 0 to a < 0. I think for a < 0, I was actually storing negative m-modes due to the switch in the frequency parametrization. So I think adding a negative sign should give smoother variations in the unwrapped phase
					sprintf(buff, "%.14e\t%.14e\t%.14e\t%.14e\n", mode.chi, alphaPoint, gsl_spline_eval(spline, alphaPoint, acc), sgn(mode.spin)*gsl_spline_eval(spline_2, alphaPoint, acc_2));
					file << buff;
				}

				gsl_spline_free(spline);
				gsl_interp_accel_free(acc);
				gsl_spline_free(spline_2);
				gsl_interp_accel_free(acc_2);
			}

			file.close();
		}
	}
}

Vector bh_spin_sample(double amax, int sampleN){
	double delta_a = 2.*amax/(sampleN - 1);
	Vector aVec(sampleN);
	for(int i = 0; i < sampleN; i++){
		aVec[i] = amax - i*delta_a;
	}
	return aVec;
}

void generate_adiabatic_inspiral_data_2d(std::string mainDir, int sampleSize_spin, int sampleSize_radius){
	int lmax = 15;
	int sampleSize_radius_base = 257; // this seems to be good enough sampling to get interpolated fluxes with fractional errors <1e-7

	Vector spinList = bh_spin_sample_chi(A_MAX, sampleSize_spin);
	std::vector<std::string> spinDirs(spinList.size());
	if(!(mainDir.back() == '/')){
		mainDir += "/";
	}
	if(!boost::filesystem::exists(mainDir)){
		boost::filesystem::create_directory(mainDir);
	}
	for(size_t i = 0; i < spinList.size(); i++){
	// for(size_t i = 0; i < 1; i++){
		std::string dir;
		if(spinList[i] > 0){
			dir = mainDir + "a" + std::to_string(int(abs(spinList[i])*1000000));
		}else{
			dir = mainDir + "an" + std::to_string(int(abs(spinList[i])*1000000));
		}
		spinDirs[i] = dir;
		std::cout << "a = " << spinList[i] << ", chi = " << chi_of_a(spinList[i]) << "\n";
		if(!boost::filesystem::exists(dir)){
			std::cout << "(INSPIRAL) Output spherical mode data to lmax = "<<lmax<<"\n";
			output_adiabatic_circ_spherical_mode_data_radial_sample_parallel(spinList[i], lmax, sampleSize_radius_base, dir);
			std::cout << "(INSPIRAL) Calculate flux remainder to l-modes > lmax \n";
			generate_adiabatic_circ_remainder_data_radial_sample(spinList[i], lmax + 1, sampleSize_radius_base, dir);
		}
		std::cout << "(INSPIRAL) Generating inspiral data from fluxes \n";
		generate_adiabatic_circ_trajectory(dir, sampleSize_radius);
	}

	compile_and_output_trajectory_data(spinDirs, mainDir + "trajectory/");
	compile_and_output_mode_amplitude_data(spinDirs, mainDir + "mode_amplitudes/", lmax, sampleSize_radius);
}

void generate_adiabatic_inspiral_data_1d(std::string mainDir, double spin, int sampleSize_radius){
	int lmax = 15;
	int sampleSize_radius_base = pow(2, 9) + 1; // this seems to be good enough sampling to get interpolated fluxes with fractional errors <1e-7

	if(!(mainDir.back() == '/')){
		mainDir += "/";
	}
	if(!boost::filesystem::exists(mainDir)){
		boost::filesystem::create_directory(mainDir);
	}

	std::string dir;
	if(spin > 0){
		dir = mainDir + "a" + std::to_string(int(abs(spin)*1000000));
	}else{
		dir = mainDir + "an" + std::to_string(int(abs(spin)*1000000));
	}
	if(!boost::filesystem::exists(dir)){
		std::cout << "(INSPIRAL) Output spherical mode data to lmax = "<<lmax<<"\n";
		output_adiabatic_circ_spherical_mode_data_radial_sample_parallel(spin, lmax, sampleSize_radius_base, dir);
		std::cout << "(INSPIRAL) Calculate flux remainder to l-modes > lmax \n";
		generate_adiabatic_circ_remainder_data_radial_sample(spin, lmax + 1, sampleSize_radius_base, dir);
	}
	std::cout << "(INSPIRAL) Generating inspiral data from fluxes \n";
	generate_adiabatic_circ_trajectory(dir, sampleSize_radius);
}

int integrate_first_order_real_ode(Vector &psi, const Vector &x, int (*sys)(double, const double*, double*, void*), double psi0, const double x0, void *params){
	size_t dim = 1;
	double h = abs(1.e-6*x0);
	if(abs(h) == 0.){
		h = 1.e-6;
	}
	if(x[1] - x[0] < 0){
		h *= -1.;
	}

	gsl_odeiv2_system gsl_sys = {sys, NULL, dim, params};

	double abs_error = 1.e-12;
	const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
	gsl_odeiv2_step* s = gsl_odeiv2_step_alloc(T, dim);
	gsl_odeiv2_control* c = gsl_odeiv2_control_y_new(abs_error, abs_error);
	gsl_odeiv2_evolve* e = gsl_odeiv2_evolve_alloc(dim);

	int stepNum = x.size(), status;
	int stepCount = 0;
	double xi = x0;
	double y[1];
	y[0] = psi0;
	for(int i = 0; i < stepNum; i++){
		if( h > 0 ){
			while(xi < x[i]){
				status = gsl_odeiv2_evolve_apply(e, c, s, &gsl_sys, &xi, x[i], &h, y);
				stepCount += 1;
			}
		}else{
			while(xi > x[i]){
				status = gsl_odeiv2_evolve_apply(e, c, s, &gsl_sys, &xi, x[i], &h, y);
				stepCount += 1;
			}
		}
		psi[i] = y[0];
	}
	// std::cout << "Step count = " << stepCount << "\n";

	gsl_odeiv2_evolve_free(e);
	gsl_odeiv2_control_free(c);
	gsl_odeiv2_step_free(s);

	return status;
}
