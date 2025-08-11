#include "resflux.hpp"

#define MODEMAX 100
#define FLUX_EPSILON 5.e-14
#define AMP_EPSILON 5.e-30
#define FLUX_LMAX 100
#define COUNT_MAX 1000
#define ZERO_FREQ_MAX 1.e-8

int minimize_resonant_mode_numbers(int& k, int& n, int nth, int nr, int Nres){
	k = round((nth*Nres)/(nr*nr + nth*nth));
	n = Nres - nth*k;
	while(std::abs(n) % nr != 0){
		k--;
		n = Nres - nth*k;
	}
	n = n/nr;
	return 0;
}

int minimize_resonant_polar_mode_number(int& k, int& n, int l, int m, int nth, int nr, int Nres){
	k = l - m;
	n = Nres - nth*k;
	while(std::abs(n) % nr != 0){
		k--;
		n = Nres - nth*k;
	}
	n = n/nr;
	return 0;
}

double horizonAmplitude(double a, int m, double omega, double lambda){
	return (pow(lambda + 2., 2) + 4*m*a*omega - pow(2*a*omega, 2))*(pow(lambda, 2) + 36.*m*a*omega - pow(6.*a*omega, 2))
		+ (2.*lambda + 3)*(96.*pow(a*omega, 2) - 48.*m*a*omega) + pow(12.*omega, 2)*(1. - a*a);
}

double horizonConstant(double a, int m, double omega, double lambda){
	double rplus = 1. + sqrt(1. - a*a);
	double k = omega - m*a/(2.*rplus);
	double epsilon = (rplus - 1.)/(4.*rplus);
	double Alm = horizonAmplitude(a, m, omega, lambda);

	return 256.*pow(2.*rplus, 5)*k*(k*k + 4.*epsilon*epsilon)*(k*k + 16.*epsilon*epsilon)*pow(omega, 3)/Alm;
}

Fluxes gravitational_flux_mode(double a, int m, double omega, double lambda, double horizonAmplitude, double infinityAmplitude){
	double prefactor = pow(2.*omega, -2)/M_PI;
	Fluxes eflux(2.*prefactor*infinityAmplitude, 2.*horizonConstant(a, m, omega, lambda)*prefactor*horizonAmplitude);
	return eflux;
}

Fluxes scalar_flux_mode(double a, int m, double omega, double horizonAmplitude, double infinityAmplitude){
	double prefactor = 0.25/M_PI;
	// (r_+^2 + a^2)/M^2 = 2 r_+/M factor due to normalization at the horizon
	double rplus = 1. + sqrt(1. - a*a);
	double horizonPrefactor = 2.*(1. + sqrt(1. - pow(a, 2)));
	double horizonFrequency = omega - m*a/(2.*rplus);
	Fluxes eflux(2.*prefactor*omega*infinityAmplitude, 2.*horizonPrefactor*prefactor*horizonFrequency*horizonAmplitude);
	return eflux;
}

Fluxes flux_mode(int s, double a, int m, double omega, double lambda, double horizonAmplitude, double infinityAmplitude){
	if(s == 0){
		return scalar_flux_mode(a, m, omega, horizonAmplitude, infinityAmplitude);
	}else{
		return gravitational_flux_mode(a, m, omega, lambda, horizonAmplitude, infinityAmplitude);
	}
}

int minimum_resonant_harmonic(int m, int nth, int nr, GeodesicSource& geo){
	double omegaM = geo.getTimeFrequency(m, 0, 0);
	double omegaRes = geo.getTimeFrequency(1)/nr;
	double omegaResCheck = geo.getTimeFrequency(2)/nth;
	if(std::abs(omegaRes - omegaResCheck) > ZERO_FREQ_MAX ){std::cout << "(RESFLUX): Frequencies do not form a "<<nr<<":"<<nth<<" \n"; }
	int Nres = 0;
	if(omegaM > 0){
		while(std::abs(omegaM + Nres*omegaRes) > ZERO_FREQ_MAX && (omegaM + Nres*omegaRes) > 0.){
			Nres--;
		}
		if(std::abs(omegaM + Nres*omegaRes) > ZERO_FREQ_MAX){
			Nres++;
		}
	}else{
		while(std::abs(omegaM + Nres*omegaRes) > ZERO_FREQ_MAX && (omegaM + Nres*omegaRes) < 0.){
			Nres++;
		}
	}
	if(std::abs(omegaM + Nres*omegaRes) > ZERO_FREQ_MAX && (omegaM + Nres*omegaRes) < 0){
		Nres++;
	}

	return Nres;
}

ResonantFluxList res_flux_l(int s, int L, int nth, int nr, GeodesicSource& geo){
	ResonantFluxList fluxMode;
	int Kmax = fluxMode.Edot.size();

	Vector efluxInf(Kmax, 0.);
	Vector efluxHor(Kmax, 0.);
	double efluxPreviousInf;
	double efluxPreviousHor;
	double efluxRef = 0.5*energy_flux_newtonian(s, geo);

	Vector lfluxInf(Kmax, 0.);
	Vector lfluxHor(Kmax, 0.);
	double lfluxPreviousInf;
	double lfluxPreviousHor;
	double lfluxRef = 0.5*angular_momentum_flux_newtonian(s, geo);

	Vector qfluxInf(Kmax, 0.);
	Vector qfluxHor(Kmax, 0.);
	double qfluxPreviousInf;
	double qfluxPreviousHor;
	double qfluxRef = 0.5*carter_flux_newtonian(s, geo);
	if(qfluxRef == 0){
		qfluxRef = 1.;
	}
	int convergenceCriteria = 2;
	int PRINT_INF_MODES = 1;
	int PRINT_HOR_MODES = 0;

	int m = L;
	fluxMode = res_flux_lm(s, L, m, nth, nr, geo);
	efluxPreviousInf = fluxMode.Edot[0].infinity;
	efluxPreviousHor = fluxMode.Edot[0].horizon;
	lfluxPreviousInf = fluxMode.Ldot[0].infinity;
	lfluxPreviousHor = fluxMode.Ldot[0].horizon;
	qfluxPreviousInf = fluxMode.Qdot[0].infinity;
	qfluxPreviousHor = fluxMode.Qdot[0].horizon;
	for(int K = 0; K < Kmax; K++){
		efluxInf[K] += fluxMode.Edot[K].infinity;
		efluxHor[K] += fluxMode.Edot[K].horizon;
		lfluxInf[K] += fluxMode.Ldot[K].infinity;
		lfluxHor[K] += fluxMode.Ldot[K].horizon;
		qfluxInf[K] += fluxMode.Qdot[K].infinity;
		qfluxHor[K] += fluxMode.Qdot[K].horizon;
	}
	if(PRINT_INF_MODES){
		std::cout << "Edot_infinity ("<<L<<", "<<m<<")-mode = "<< efluxPreviousInf << " \n";
		// std::cout << "Ldot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< lfluxPreviousInf << " \n";
		// std::cout << "Qdot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< qfluxPreviousInf << " \n";
	}
	if(PRINT_HOR_MODES){
		std::cout << "Edot_horizon ("<<L<<", "<<m<<")-mode = "<< efluxPreviousHor << " \n";
	}
	m--;
	int convergenceCheckE = 0;
	int convergenceCheckL = 0;
	int convergenceCheckQ = 0;
	while(m > L - 2 && m > 0){
		fluxMode = res_flux_lm(s, L, m, nth, nr, geo);
		if(std::abs(fluxMode.Edot[0].infinity/efluxRef) < FLUX_EPSILON/10 && (fluxMode.Edot[0].infinity < efluxPreviousInf || fluxMode.Edot[0].infinity < 1.e-18) && std::abs(fluxMode.Edot[0].horizon/efluxRef) < FLUX_EPSILON/10 && (fluxMode.Edot[0].horizon < efluxPreviousHor || fluxMode.Edot[0].horizon < 1.e-18)){
			convergenceCheckE += 1;
		}else{
			if(convergenceCheckE < COUNT_MAX){
				convergenceCheckE = 0;
			}
		}
		if(std::abs(fluxMode.Ldot[0].infinity/lfluxRef) < FLUX_EPSILON/10 && (fluxMode.Ldot[0].infinity < lfluxPreviousInf || fluxMode.Ldot[0].infinity < 1.e-18) && std::abs(fluxMode.Ldot[0].horizon/lfluxRef) < FLUX_EPSILON/10 && (fluxMode.Ldot[0].horizon < lfluxPreviousHor || fluxMode.Ldot[0].horizon < 1.e-18)){
			convergenceCheckL += 1;
		}else{
			if(convergenceCheckL < COUNT_MAX){
				convergenceCheckL = 0;
			}
		}
		if(std::abs(fluxMode.Qdot[0].infinity/qfluxRef) < FLUX_EPSILON/10 && (fluxMode.Qdot[0].infinity < qfluxPreviousInf || fluxMode.Qdot[0].infinity < 1.e-18) && std::abs(fluxMode.Qdot[0].horizon/qfluxRef) < FLUX_EPSILON/10 && (fluxMode.Qdot[0].horizon < qfluxPreviousHor || fluxMode.Qdot[0].horizon < 1.e-18)){
			convergenceCheckQ += 1;
		}else{
			if(convergenceCheckQ < COUNT_MAX){
				convergenceCheckQ = 0;
			}
		}
		efluxPreviousInf = fluxMode.Edot[0].infinity;
		efluxPreviousHor = fluxMode.Edot[0].horizon;
		lfluxPreviousInf = fluxMode.Ldot[0].infinity;
		lfluxPreviousHor = fluxMode.Ldot[0].horizon;
		qfluxPreviousInf = fluxMode.Qdot[0].infinity;
		qfluxPreviousHor = fluxMode.Qdot[0].horizon;

		for(int K = 0; K < Kmax; K++){
			if(convergenceCheckE < convergenceCriteria){
				efluxInf[K] += fluxMode.Edot[K].infinity;
				efluxHor[K] += fluxMode.Edot[K].horizon;
			}else{
				convergenceCheckE = COUNT_MAX;
			}
			if(convergenceCheckL < convergenceCriteria){
				lfluxInf[K] += fluxMode.Ldot[K].infinity;
				lfluxHor[K] += fluxMode.Ldot[K].horizon;
			}else{
				convergenceCheckL = COUNT_MAX;
			}
			if(convergenceCheckQ < convergenceCriteria){
				qfluxInf[K] += fluxMode.Qdot[K].infinity;
				qfluxHor[K] += fluxMode.Qdot[K].horizon;
			}else{
				convergenceCheckQ = COUNT_MAX;
			}
		}
		if(PRINT_INF_MODES){
			std::cout << "Edot_infinity ("<<L<<", "<<m<<")-mode = "<< efluxPreviousInf << " \n";
			// std::cout << "Ldot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< lfluxPreviousInf << " \n";
			// std::cout << "Qdot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< qfluxPreviousInf << " \n";
		}
		if(PRINT_HOR_MODES){
			std::cout << "Edot_horizon ("<<L<<", "<<m<<")-mode = "<< efluxPreviousHor << " \n";
		}
		m--;
	}

	int minM = -L;
	if(geo.getEccentricity() < DBL_EPSILON && 1. - std::abs(geo.getInclination()) < DBL_EPSILON){
		minM = 1;
	}

	while(m >= minM && (convergenceCheckE < convergenceCriteria || convergenceCheckL < convergenceCriteria || convergenceCheckQ < convergenceCriteria)){
		fluxMode = res_flux_lm(s, L, m, nth, nr, geo);
		if(std::abs(fluxMode.Edot[0].infinity/efluxRef) < FLUX_EPSILON/10 && (fluxMode.Edot[0].infinity < efluxPreviousInf || fluxMode.Edot[0].infinity < 1.e-18) && std::abs(fluxMode.Edot[0].horizon/efluxRef) < FLUX_EPSILON/10 && (fluxMode.Edot[0].horizon < efluxPreviousHor || fluxMode.Edot[0].horizon < 1.e-18)){
			convergenceCheckE += 1;
		}else{
			if(convergenceCheckE < COUNT_MAX){
				convergenceCheckE = 0;
			}
		}
		if(std::abs(fluxMode.Ldot[0].infinity/lfluxRef) < FLUX_EPSILON/10 && (fluxMode.Ldot[0].infinity < lfluxPreviousInf || fluxMode.Ldot[0].infinity < 1.e-18) && std::abs(fluxMode.Ldot[0].horizon/lfluxRef) < FLUX_EPSILON/10 && (fluxMode.Ldot[0].horizon < lfluxPreviousHor || fluxMode.Ldot[0].horizon < 1.e-18)){
			convergenceCheckL += 1;
		}else{
			if(convergenceCheckL < COUNT_MAX){
				convergenceCheckL = 0;
			}
		}
		if(std::abs(fluxMode.Qdot[0].infinity/qfluxRef) < FLUX_EPSILON/10 && (fluxMode.Qdot[0].infinity < qfluxPreviousInf || fluxMode.Qdot[0].infinity < 1.e-18) && std::abs(fluxMode.Qdot[0].horizon/qfluxRef) < FLUX_EPSILON/10 && (fluxMode.Qdot[0].horizon < qfluxPreviousHor || fluxMode.Qdot[0].horizon < 1.e-18)){
			convergenceCheckQ += 1;
		}else{
			if(convergenceCheckQ < COUNT_MAX){
				convergenceCheckQ = 0;
			}
		}
		efluxPreviousInf = fluxMode.Edot[0].infinity;
		efluxPreviousHor = fluxMode.Edot[0].horizon;
		lfluxPreviousInf = fluxMode.Ldot[0].infinity;
		lfluxPreviousHor = fluxMode.Ldot[0].horizon;
		qfluxPreviousInf = fluxMode.Qdot[0].infinity;
		qfluxPreviousHor = fluxMode.Qdot[0].horizon;

		for(int K = 0; K < Kmax; K++){
			if(convergenceCheckE < convergenceCriteria){
				efluxInf[K] += fluxMode.Edot[K].infinity;
				efluxHor[K] += fluxMode.Edot[K].horizon;
			}else{
				convergenceCheckE = COUNT_MAX;
			}
			if(convergenceCheckL < convergenceCriteria){
				lfluxInf[K] += fluxMode.Ldot[K].infinity;
				lfluxHor[K] += fluxMode.Ldot[K].horizon;
			}else{
				convergenceCheckL = COUNT_MAX;
			}
			if(convergenceCheckQ < convergenceCriteria){
				qfluxInf[K] += fluxMode.Qdot[K].infinity;
				qfluxHor[K] += fluxMode.Qdot[K].horizon;
			}else{
				convergenceCheckQ = COUNT_MAX;
			}
		}
		if(PRINT_INF_MODES){
			std::cout << "Edot_infinity ("<<L<<", "<<m<<")-mode = "<< efluxPreviousInf << " \n";
			// std::cout << "Ldot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< lfluxPreviousInf << " \n";
			// std::cout << "Qdot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< qfluxPreviousInf << " \n";
		}
		if(PRINT_HOR_MODES){
			std::cout << "Edot_horizon ("<<L<<", "<<m<<")-mode = "<< efluxPreviousHor << " \n";
		}
		m--;
	}

	ResonantFluxList fluxes;
	for(int K = 0; K < Kmax; K++){
		fluxes.Edot[K].infinity = efluxInf[K];
		fluxes.Edot[K].horizon = efluxHor[K];
		fluxes.Ldot[K].infinity = lfluxInf[K];
		fluxes.Ldot[K].horizon = lfluxHor[K];
		fluxes.Qdot[K].infinity = qfluxInf[K];
		fluxes.Qdot[K].horizon = qfluxHor[K];
	}

	return fluxes;
}

ResonantFluxList res_flux_lm(int s, int L, int m, int nth, int nr, GeodesicSource& geo){
	ResonantFluxList fluxMode;
	int Kmax = fluxMode.Edot.size();

	Vector efluxInf(Kmax, 0.);
	Vector efluxHor(Kmax, 0.);
	double efluxPreviousInf;
	double efluxPreviousHor;
	double efluxRef = 0.5*energy_flux_newtonian(s, geo);

	Vector lfluxInf(Kmax, 0.);
	Vector lfluxHor(Kmax, 0.);
	double lfluxPreviousInf;
	double lfluxPreviousHor;
	double lfluxRef = 0.5*angular_momentum_flux_newtonian(s, geo);

	Vector qfluxInf(Kmax, 0.);
	Vector qfluxHor(Kmax, 0.);
	double qfluxPreviousInf;
	double qfluxPreviousHor;
	double qfluxRef = 0.5*carter_flux_newtonian(s, geo);
	if(qfluxRef == 0){
		qfluxRef = 1.;
	}
	// factor of 0.5 is a fudge factor because Newtonian fluxes usually overestimate
	int convergenceCriteria = 10;

	int NPeak = nth*(L - m);
	int NMax = NPeak + 100;
	// calculate the minimum radial harmonic that still gives us a postive frequency
	int NMin = minimum_resonant_harmonic(m, nth, nr, geo);
	if(NMin >= NPeak){
		NPeak = NMin + 1;
	}

	int PRINT_INF_MODES = 0;
	int PRINT_HOR_MODES = 0;

	int Nres = NPeak;
	fluxMode = res_flux_lmN(s, L, m, Nres, nth, nr, geo);
	efluxPreviousInf = fluxMode.Edot[0].infinity;
	efluxPreviousHor = fluxMode.Edot[0].horizon;
	lfluxPreviousInf = fluxMode.Ldot[0].infinity;
	lfluxPreviousHor = fluxMode.Ldot[0].horizon;
	qfluxPreviousInf = fluxMode.Qdot[0].infinity;
	qfluxPreviousHor = fluxMode.Qdot[0].horizon;
	for(int K = 0; K < Kmax; K++){
		efluxInf[K] += fluxMode.Edot[K].infinity;
		efluxHor[K] += fluxMode.Edot[K].horizon;
		lfluxInf[K] += fluxMode.Ldot[K].infinity;
		lfluxHor[K] += fluxMode.Ldot[K].horizon;
		qfluxInf[K] += fluxMode.Qdot[K].infinity;
		qfluxHor[K] += fluxMode.Qdot[K].horizon;
	}
	if(PRINT_INF_MODES){
		std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< efluxPreviousInf << " \n";
		// std::cout << "Ldot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< lfluxPreviousInf << " \n";
		// std::cout << "Qdot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< qfluxPreviousInf << " \n";
	}
	if(PRINT_HOR_MODES){
		std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< efluxPreviousHor << " \n";
	}
	Nres++;

	int convergenceCheckE = 0;
	int convergenceCheckL = 0;
	int convergenceCheckQ = 0;
	while(Nres <= NMax && (convergenceCheckE < convergenceCriteria || convergenceCheckL < convergenceCriteria || convergenceCheckQ < convergenceCriteria)){
		fluxMode = res_flux_lmN(s, L, m, Nres, nth, nr, geo);
		if(std::abs(fluxMode.Edot[0].infinity/efluxRef) < FLUX_EPSILON/10 && (fluxMode.Edot[0].infinity < efluxPreviousInf || fluxMode.Edot[0].infinity < 1.e-18) && std::abs(fluxMode.Edot[0].horizon/efluxRef) < FLUX_EPSILON/10 && (fluxMode.Edot[0].horizon < efluxPreviousHor || fluxMode.Edot[0].horizon < 1.e-18)){
			convergenceCheckE += 1;
		}else{
			if(convergenceCheckE < COUNT_MAX){
				convergenceCheckE = 0;
			}
		}
		if(std::abs(fluxMode.Ldot[0].infinity/lfluxRef) < FLUX_EPSILON/10 && (fluxMode.Ldot[0].infinity < lfluxPreviousInf || fluxMode.Ldot[0].infinity < 1.e-18) && std::abs(fluxMode.Ldot[0].horizon/lfluxRef) < FLUX_EPSILON/10 && (fluxMode.Ldot[0].horizon < lfluxPreviousHor || fluxMode.Ldot[0].horizon < 1.e-18)){
			convergenceCheckL += 1;
		}else{
			if(convergenceCheckL < COUNT_MAX){
				convergenceCheckL = 0;
			}
		}
		if(std::abs(fluxMode.Qdot[0].infinity/qfluxRef) < FLUX_EPSILON/10 && (fluxMode.Qdot[0].infinity < qfluxPreviousInf || fluxMode.Qdot[0].infinity < 1.e-18) && std::abs(fluxMode.Qdot[0].horizon/qfluxRef) < FLUX_EPSILON/10 && (fluxMode.Qdot[0].horizon < qfluxPreviousHor || fluxMode.Qdot[0].horizon < 1.e-18)){
			convergenceCheckQ += 1;
		}else{
			if(convergenceCheckQ < COUNT_MAX){
				convergenceCheckQ = 0;
			}
		}
		efluxPreviousInf = fluxMode.Edot[0].infinity;
		efluxPreviousHor = fluxMode.Edot[0].horizon;
		lfluxPreviousInf = fluxMode.Ldot[0].infinity;
		lfluxPreviousHor = fluxMode.Ldot[0].horizon;
		qfluxPreviousInf = fluxMode.Qdot[0].infinity;
		qfluxPreviousHor = fluxMode.Qdot[0].horizon;

		for(int K = 0; K < Kmax; K++){
			if(convergenceCheckE < convergenceCriteria){
				efluxInf[K] += fluxMode.Edot[K].infinity;
				efluxHor[K] += fluxMode.Edot[K].horizon;
			}else{
				convergenceCheckE = COUNT_MAX;
			}
			if(convergenceCheckL < convergenceCriteria){
				lfluxInf[K] += fluxMode.Ldot[K].infinity;
				lfluxHor[K] += fluxMode.Ldot[K].horizon;
			}else{
				convergenceCheckL = COUNT_MAX;
			}
			if(convergenceCheckQ < convergenceCriteria){
				qfluxInf[K] += fluxMode.Qdot[K].infinity;
				qfluxHor[K] += fluxMode.Qdot[K].horizon;
			}else{
				convergenceCheckQ = COUNT_MAX;
			}
		}
		if(PRINT_INF_MODES){
			std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< efluxPreviousInf << " \n";
			// std::cout << "Ldot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< lfluxPreviousInf << " \n";
			// std::cout << "Qdot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< qfluxPreviousInf << " \n";
		}
		if(PRINT_HOR_MODES){
			std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< efluxPreviousHor << " \n";
		}
		Nres++;
	}

	Nres = NPeak - 1;
	fluxMode = res_flux_lmN(s, L, m, Nres, nth, nr, geo);
	efluxPreviousInf = fluxMode.Edot[0].infinity;
	efluxPreviousHor = fluxMode.Edot[0].horizon;
	lfluxPreviousInf = fluxMode.Ldot[0].infinity;
	lfluxPreviousHor = fluxMode.Ldot[0].horizon;
	qfluxPreviousInf = fluxMode.Qdot[0].infinity;
	qfluxPreviousHor = fluxMode.Qdot[0].horizon;
	for(int K = 0; K < Kmax; K++){
		efluxInf[K] += fluxMode.Edot[K].infinity;
		efluxHor[K] += fluxMode.Edot[K].horizon;
		lfluxInf[K] += fluxMode.Ldot[K].infinity;
		lfluxHor[K] += fluxMode.Ldot[K].horizon;
		qfluxInf[K] += fluxMode.Qdot[K].infinity;
		qfluxHor[K] += fluxMode.Qdot[K].horizon;
	}
	if(PRINT_INF_MODES){
		std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< efluxPreviousInf << " \n";
		// std::cout << "Ldot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< lfluxPreviousInf << " \n";
		// std::cout << "Qdot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< qfluxPreviousInf << " \n";
	}
	if(PRINT_HOR_MODES){
		std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< efluxPreviousHor << " \n";
	}
	Nres--;
	convergenceCheckE = 0;
	convergenceCheckL = 0;
	convergenceCheckQ = 0;
	while(Nres >= NMin && (convergenceCheckE < convergenceCriteria || convergenceCheckL < convergenceCriteria || convergenceCheckQ < convergenceCriteria)){
		fluxMode = res_flux_lmN(s, L, m, Nres, nth, nr, geo);
		if(std::abs(fluxMode.Edot[0].infinity/efluxRef) < FLUX_EPSILON/10 && (fluxMode.Edot[0].infinity < efluxPreviousInf || fluxMode.Edot[0].infinity < 1.e-18) && std::abs(fluxMode.Edot[0].horizon/efluxRef) < FLUX_EPSILON/10 && (fluxMode.Edot[0].horizon < efluxPreviousHor || fluxMode.Edot[0].horizon < 1.e-18)){
			convergenceCheckE += 1;
		}else{
			if(convergenceCheckE < COUNT_MAX){
				convergenceCheckE = 0;
			}
		}
		if(std::abs(fluxMode.Ldot[0].infinity/lfluxRef) < FLUX_EPSILON/10 && (fluxMode.Ldot[0].infinity < lfluxPreviousInf || fluxMode.Ldot[0].infinity < 1.e-18) && std::abs(fluxMode.Ldot[0].horizon/lfluxRef) < FLUX_EPSILON/10 && (fluxMode.Ldot[0].horizon < lfluxPreviousHor || fluxMode.Ldot[0].horizon < 1.e-18)){
			convergenceCheckL += 1;
		}else{
			if(convergenceCheckL < COUNT_MAX){
				convergenceCheckL = 0;
			}
		}
		if(std::abs(fluxMode.Qdot[0].infinity/qfluxRef) < FLUX_EPSILON/10 && (fluxMode.Qdot[0].infinity < qfluxPreviousInf || fluxMode.Qdot[0].infinity < 1.e-18) && std::abs(fluxMode.Qdot[0].horizon/qfluxRef) < FLUX_EPSILON/10 && (fluxMode.Qdot[0].horizon < qfluxPreviousHor || fluxMode.Qdot[0].horizon < 1.e-18)){
			convergenceCheckQ += 1;
		}else{
			if(convergenceCheckQ < COUNT_MAX){
				convergenceCheckQ = 0;
			}
		}
		efluxPreviousInf = fluxMode.Edot[0].infinity;
		efluxPreviousHor = fluxMode.Edot[0].horizon;
		lfluxPreviousInf = fluxMode.Ldot[0].infinity;
		lfluxPreviousHor = fluxMode.Ldot[0].horizon;
		qfluxPreviousInf = fluxMode.Qdot[0].infinity;
		qfluxPreviousHor = fluxMode.Qdot[0].horizon;

		for(int K = 0; K < Kmax; K++){
			if(convergenceCheckE < convergenceCriteria){
				efluxInf[K] += fluxMode.Edot[K].infinity;
				efluxHor[K] += fluxMode.Edot[K].horizon;
			}else{
				convergenceCheckE = COUNT_MAX;
			}
			if(convergenceCheckL < convergenceCriteria){
				lfluxInf[K] += fluxMode.Ldot[K].infinity;
				lfluxHor[K] += fluxMode.Ldot[K].horizon;
			}else{
				convergenceCheckL = COUNT_MAX;
			}
			if(convergenceCheckQ < convergenceCriteria){
				qfluxInf[K] += fluxMode.Qdot[K].infinity;
				qfluxHor[K] += fluxMode.Qdot[K].horizon;
			}else{
				convergenceCheckQ = COUNT_MAX;
			}
		}
		if(PRINT_INF_MODES){
			std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< efluxPreviousInf << " \n";
			// std::cout << "Ldot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< lfluxPreviousInf << " \n";
			// std::cout << "Qdot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< qfluxPreviousInf << " \n";
		}
		if(PRINT_HOR_MODES){
			std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< efluxPreviousHor << " \n";
		}
		Nres--;
	}
	ResonantFluxList fluxes;
	for(int K = 0; K < Kmax; K++){
		fluxes.Edot[K].infinity = efluxInf[K];
		fluxes.Edot[K].horizon = efluxHor[K];
		fluxes.Ldot[K].infinity = lfluxInf[K];
		fluxes.Ldot[K].horizon = lfluxHor[K];
		fluxes.Qdot[K].infinity = qfluxInf[K];
		fluxes.Qdot[K].horizon = qfluxHor[K];
	}

	return fluxes;
}

ResonantFluxList res_flux_lmN(int s, int L, int m, int Nres, int nth, int nr, GeodesicSource& geo){
	ResonantFluxList fluxes;
	int k0 = 0, n0 = 0;
	minimize_resonant_polar_mode_number(k0, n0, L, m, nth, nr, Nres);
	// k0 = Nres/nth;

	double a = geo.getBlackHoleSpin();
	double omega = geo.getTimeFrequency(m, k0, n0);
	if(std::abs(omega) < ZERO_FREQ_MAX){
		return fluxes;
	}

	RadialTeukolsky Rt(a, s, L, m, omega, geo.getRadialPosition());
	Rt.generateSolutions(AUTO);
	SpinWeightedHarmonic Slm(s, L, m, a*omega, geo.getPolarPosition());
	Slm.generateSolutions();
	double lambda = Slm.getEigenvalue();

	FieldAmplitudeVector amps = resonant_flux_mode_coefficients(s, L, m, k0, n0, nth, nr, Slm, Rt, geo);

	int modeMax = amps.size()/2;
	int Kmax = fluxes.Edot.size();
	for(int K = 0; K < Kmax; K++){
		int Kres = nth*K;
		double baseAmplitudeInf = 0.;
		double EnAmplitudeInf = 0.;
		double LzAmplitudeInf = 0.;
		double QcAmplitudeInf = 0.;

		double baseAmplitudeHor = 0.;
		double EnAmplitudeHor = 0.;
		double LzAmplitudeHor = 0.;
		double QcAmplitudeHor = 0.;

		double baseAmplitudeMix = 0.;
		double MmknMix = 0.;
		double QcAmplitudeMix = 0.;
		double rp = Rt.getRadialPoints(0);
		Complex Bup = Rt.getSolution(In, 0)*std::conj(Rt.getDerivative(Up, 0)) - std::conj(Rt.getSolution(Up, 0))*Rt.getDerivative(In, 0);
		Bup *= (rp*rp - 2.*rp + pow(Rt.getBlackHoleSpin(), 2))/(-2.*I*omega);

		double rpComp = Rt.getRadialPoints(10);
		Complex BupComp = Rt.getSolution(In, 10)*std::conj(Rt.getDerivative(Up, 10)) - std::conj(Rt.getSolution(Up, 10))*Rt.getDerivative(In, 10);
		BupComp *= (rpComp*rpComp - 2.*rpComp + pow(Rt.getBlackHoleSpin(), 2))/(-2.*I*omega);
		if(std::abs(1. - Bup/BupComp) > 1.e-8){
			std::cout << "(RESFLUX) Error: Bup calculation failed for ("<<m<<", "<<Nres<<")-mode\n";
		}


		int PRINT_INF_MODES = 0;
		int PRINT_HOR_MODES = 0;

		double Mmkn = 0.;
		if( std::abs(nth*k0 + nr*n0 - Nres) > ZERO_FREQ_MAX ){
			std::cout << "(RESFLUX) ERROR: Initial mode numbers ("<<k0<<", "<<n0<<") do not match frequency. \n";
		}
		for(int j = 0; j < modeMax; j++){
			int kres = nth*k0 + j;
			int nres = nr*n0 - j;
			//int nres = npres - Kres;
			//int kres = Nres - npres + Kres;
			if((kres % nth) == 0 && (nres % nr) == 0){
				// if(K == 2){std::cout << "("<< Nres <<", "<<kres/nth<<", "<<nres/nr<<") \n";}
				Mmkn = geo.getCarterFrequency(m, (kres + Kres)/nth, (nres - Kres)/nr) + geo.getCarterFrequency(m, kres/nth, nres/nr);
				MmknMix = geo.getCarterFrequency(m, (kres + Kres)/nth, (nres - Kres)/nr) - geo.getCarterFrequency(m, kres/nth, nres/nr);
			}else{
				Mmkn = 0.;
				MmknMix = 0.;
			}

			int jp = j + Kres;
			if(jp < 0){
				jp = -jp;
				// if(K == 2 && (kres % nth) == 0 && (nres % nr) == 0){
				// 	std::cout << "j = "<<j<<", jp = "<<2*modeMax - jp<<" and C_("<<L<<", "<<m<<", "<<kres/nth<<", "<<nres/nr<<") x C_("<<L<<", "<<m<<", "<<(nth*k0 - jp)/nth<<", "<<(nr*n0 + jp)/nr<<")\n";
				// }
				baseAmplitudeInf = std::real(amps[j].up*std::conj(amps[2*modeMax - jp].up));
				baseAmplitudeHor = std::real(amps[j].in*std::conj(amps[2*modeMax - jp].in));
				baseAmplitudeMix = std::imag(Bup*amps[j].in*std::conj(amps[2*modeMax - jp].up)) - std::imag(std::conj(Bup)*amps[j].up*std::conj(amps[2*modeMax - jp].in));
			}else{
				// if(K == 2 && (kres % nth) == 0 && (nres % nr) == 0){
				// 	std::cout << "j = "<<j<<", jp = "<<jp<<" and C_("<<L<<", "<<m<<", "<<kres/nth<<", "<<nres/nr<<") x C_("<<L<<", "<<m<<", "<<(nth*k0 + jp)/nth<<", "<<(nr*n0 - jp)/nr<<")\n";
				// }
				baseAmplitudeInf = std::real(amps[j].up*std::conj(amps[jp].up));
				baseAmplitudeHor = std::real(amps[j].in*std::conj(amps[jp].in));
				baseAmplitudeMix = std::imag(Bup*amps[j].in*std::conj(amps[jp].up)) - std::imag(std::conj(Bup)*amps[j].up*std::conj(amps[jp].in));
				// if(K == 2 && (kres % nth) == 0 && (nres % nr) == 0 && std::abs(Nres) < 30 && std::abs(baseAmplitudeMix) > 0.){
				// 	std::cout << Bup*amps[j].in*std::conj(amps[jp].up) << " vs " << std::conj(Bup)*amps[j].up*std::conj(amps[jp].in) << "\n";
				// }
			}
			EnAmplitudeInf += omega*baseAmplitudeInf;
			LzAmplitudeInf += m*baseAmplitudeInf;
			QcAmplitudeInf += Mmkn*baseAmplitudeInf;

			EnAmplitudeHor += omega*baseAmplitudeHor;
			LzAmplitudeHor += m*baseAmplitudeHor;
			QcAmplitudeHor += Mmkn*baseAmplitudeHor;

			QcAmplitudeMix += MmknMix*baseAmplitudeMix;

			if(std::abs(baseAmplitudeInf) > 0. && std::abs(Mmkn) == 0.){
				std::cout << "(RESFLUX) Error: Zero frequency Carter frequency for ("<<m<<", "<<kres/nth<<", "<<nres/nr<<") \n";
			}

			if(PRINT_INF_MODES && K == 0){
				std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<kres/double(nth)<<", "<<nres/double(nr)<<")-mode amplitude = "<< baseAmplitudeInf << " \n";
				// std::cout << "Ldot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< lfluxPreviousInf << " \n";
				// std::cout << "Qdot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< qfluxPreviousInf << " \n";
			}
			if(PRINT_HOR_MODES && K == 0){
				std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<kres/double(nth)<<", "<<nres/double(nr)<<")-mode amplitude = "<< baseAmplitudeHor << " \n";
			}

			if(j > 0){
				kres = nth*k0 - j;
				nres = nr*n0 + j;
				// nres = npres - Kres;
				// kres = Nres - npres + Kres;
				if((kres % nth) == 0 && (nres % nr) == 0){
					// if(K == 2){std::cout << "("<< Nres <<", "<<kres/nth<<", "<<nres/nr<<") \n";}
					Mmkn = geo.getCarterFrequency(m, (kres + Kres)/nth, (nres - Kres)/nr) + geo.getCarterFrequency(m, kres/nth, nres/nr);
					MmknMix = geo.getCarterFrequency(m, (kres + Kres)/nth, (nres - Kres)/nr) - geo.getCarterFrequency(m, kres/nth, nres/nr);
				}else{
					Mmkn = 0.;
					MmknMix = 0.;
				}

				jp = -j + Kres;
				if(jp < 0){
					jp = -jp;
					// if(K == 2 && (kres % nth) == 0 && (nres % nr) == 0){
					// 	std::cout << "j = "<<2*modeMax - j<<", jp = "<<2*modeMax - jp<<" and C_("<<L<<", "<<m<<", "<<kres/nth<<", "<<nres/nr<<") x C_("<<L<<", "<<m<<", "<<(nth*k0 - jp)/nth<<", "<<(nr*n0 + jp)/nr<<")\n";
					// }
					baseAmplitudeInf = std::real(amps[2*modeMax - j].up*std::conj(amps[2*modeMax - jp].up));
					baseAmplitudeHor = std::real(amps[2*modeMax - j].in*std::conj(amps[2*modeMax - jp].in));
					baseAmplitudeMix = std::imag(Bup*amps[2*modeMax - j].in*std::conj(amps[2*modeMax - jp].up)) - std::imag(std::conj(Bup)*amps[2*modeMax - j].up*std::conj(amps[2*modeMax - jp].in));
				}else{
					// if(K == 2 && (kres % nth) == 0 && (nres % nr) == 0){
					// 	std::cout << "j = "<<2*modeMax - j<<", jp = "<<jp<<" and C_("<<L<<", "<<m<<", "<<kres/nth<<", "<<nres/nr<<") x C_("<<L<<", "<<m<<", "<<(nth*k0 + jp)/nth<<", "<<(nr*n0 - jp)/nr<<")\n";
					// }
					baseAmplitudeInf = std::real(amps[2*modeMax - j].up*std::conj(amps[jp].up));
					baseAmplitudeHor = std::real(amps[2*modeMax - j].in*std::conj(amps[jp].in));
					baseAmplitudeMix = std::imag(Bup*amps[2*modeMax - j].in*std::conj(amps[jp].up)) - std::imag(std::conj(Bup)*amps[2*modeMax - j].up*std::conj(amps[jp].in));
				}
				EnAmplitudeInf += omega*baseAmplitudeInf;
				LzAmplitudeInf += m*baseAmplitudeInf;
				QcAmplitudeInf += Mmkn*baseAmplitudeInf;

				EnAmplitudeHor += omega*baseAmplitudeHor;
				LzAmplitudeHor += m*baseAmplitudeHor;
				QcAmplitudeHor += Mmkn*baseAmplitudeHor;

				QcAmplitudeMix += MmknMix*baseAmplitudeMix;

				if(std::abs(baseAmplitudeInf) > 0. && std::abs(Mmkn) == 0.){
					std::cout << "(RESFLUX) Error: Zero frequency Carter frequency for ("<<m<<", "<<kres/nth<<", "<<nres/nr<<") \n";
				}

				if(PRINT_INF_MODES && K == 0){
					std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<kres/double(nth)<<", "<<nres/double(nr)<<")-mode amplitude = "<< baseAmplitudeInf << " \n";
					// std::cout << "Ldot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< lfluxPreviousInf << " \n";
					// std::cout << "Qdot_infinity ("<<L<<", "<<m<<", "<<Nres<<")-mode = "<< qfluxPreviousInf << " \n";
				}
				if(PRINT_HOR_MODES && K == 0){
					std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<kres/double(nth)<<", "<<nres/double(nr)<<")-mode amplitude = "<< baseAmplitudeHor << " \n";
				}
			}
		}
		fluxes.Edot[K] = flux_mode(s, a, m, omega, lambda, EnAmplitudeHor, EnAmplitudeInf);
		fluxes.Ldot[K] = flux_mode(s, a, m, omega, lambda, LzAmplitudeHor, LzAmplitudeInf);
		fluxes.Qdot[K] = flux_mode(s, a, m, omega, lambda, QcAmplitudeHor, QcAmplitudeInf);
		// fluxes.Qdot[K] = flux_mode(s, a, m, omega, lambda, 0., QcAmplitudeMix);
	}

	return fluxes;
}

FieldAmplitudeVector resonant_flux_mode_coefficients(int s, int L, int m, int k0, int n0, int nth, int nr, SpinWeightedHarmonic& Slm, RadialTeukolsky& Rt, GeodesicSource& geo){
	FieldAmplitudeVector resAmplitudes(2*MODEMAX);

	int kres = nth*k0;
	int nres = nr*n0;
	int Nres = nres + kres;

	int j = 0;
	int jmaxCounter = 0;
	double ampMag = 1.;
	double ampMagPrev = ampMag;
	for(int j = 0; j < MODEMAX && (ampMag + ampMagPrev) > AMP_EPSILON; j++){
		kres = nth*k0 + j;
		nres = nr*n0 - j;
		if((nres % nr) == 0 && (kres % nth) == 0){
			int k = kres/nth;
			int n = nres/nr;
			// std::cout << "j = "<<j<<" and (k, n) = ("<< k <<", "<< n <<")\n";
			TeukolskyMode teuk(s, L, m, k, n, geo);
			teuk.generateSolutions(Slm, Rt, geo.getTrajectory(), geo.getConstants());
			// teuk.generateSolutions(geo);
			resAmplitudes[j].up = teuk.getTeukolskyAmplitude(Up);
			resAmplitudes[j].in = teuk.getTeukolskyAmplitude(In);
			ampMagPrev = ampMag;
			ampMag = std::abs(resAmplitudes[j].up);
		}
	}
	if(ampMag > AMP_EPSILON){
		std::cout << "(RESFLUX) Error: Amplitudes did not converge for ("<<L<<", "<<m<<", "<<Nres<<")-mode \n";
	}

	j = 0;
	jmaxCounter = 0;
	ampMag = 1.;
	ampMagPrev = ampMag;
	for(int j = 1; j < MODEMAX && (ampMag + ampMagPrev) > AMP_EPSILON; j++){
		kres = nth*k0 - j;
		nres = nr*n0 + j;
		// std::cout << "(nr * n, nth * k) = ("<< nrn <<", "<< nthk <<")\n";
		if((nres % nr) == 0 && (kres % nth) == 0){
			int k = kres/nth;
			int n = nres/nr;
			// std::cout << "j = "<< 2*MODEMAX - j <<" and (k, n) = ("<< k <<", "<< n <<")\n";
			TeukolskyMode teuk(s, L, m, k, n, geo);
			teuk.generateSolutions(Slm, Rt, geo.getTrajectory(), geo.getConstants());
			// teuk.generateSolutions(geo);
			resAmplitudes[2*MODEMAX - j].up = teuk.getTeukolskyAmplitude(Up);
			resAmplitudes[2*MODEMAX - j].in = teuk.getTeukolskyAmplitude(In);
			ampMagPrev = ampMag;
			ampMag = std::abs(resAmplitudes[2*MODEMAX - j].up);
		}
	}
	if(ampMag > AMP_EPSILON){
		std::cout << "(RESFLUX) Error: Amplitudes did not converge for ("<<L<<", "<<m<<", "<<Nres<<")-mode \n";
	}

	return resAmplitudes;
}