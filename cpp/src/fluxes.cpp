// fluxes.cpp

#include "fluxes.hpp"

#define FLUX_EPSILON 5.e-10
#define FLUX_LMAX 100
#define COUNT_MAX 1000
#define ZERO_FREQ_MAX 1.e-8

double energy_flux_newtonian(GeodesicSource& geo){
	return energy_flux_newtonian(-2, geo);
}

double energy_flux_newtonian(int s, double omegaPhi){
	double v = pow(std::abs(omegaPhi), 1./3.);
	if(s == -2){
		return 32./5.*pow(v, 10);
	}else{ // s = 0 case
		return 1./3.*pow(v, 8); // found scalar formula in Bini et al. arXiv:1610.02235
	}
}

double energy_flux_newtonian(int s, double omegaPhi, double e){
	return energy_flux_newtonian(s, omegaPhi)*pow(1. - e*e, 1.5);
}

double angular_momentum_flux_newtonian(int s, double omegaPhi){
	double v = pow(std::abs(omegaPhi), 1./3.);
	if(s == -2){
		return 32./5.*pow(v, 7);
	}else{ // s = 0 case
		return 1./3.*pow(v, 5); // found scalar formula in Bini et al. arXiv:1610.02235
	}
}

double angular_momentum_flux_newtonian(int s, double omegaPhi, double e){
	return angular_momentum_flux_newtonian(s, omegaPhi)*pow(1. - e*e, 1.5);
}

double carter_flux_newtonian(int s, double omegaPhi){
	double v = pow(std::abs(omegaPhi), 1./3.);
	if(s == -2){
		return 64./5.*pow(v, 6);
	}else{ // s = 0 case
		return 1./3.*pow(v, 4); // inferred from gravitational formulae
	}
}

double carter_flux_newtonian(int s, double omegaPhi, double e, double x){
	return carter_flux_newtonian(s, omegaPhi)*pow(1. - e*e, 1.5)*(1. - x*x);
}

double energy_flux_newtonian(int s, GeodesicSource& geo){
	return energy_flux_newtonian(s, geo.getTimeFrequency(3), geo.getEccentricity());
}

double angular_momentum_flux_newtonian(int s, GeodesicSource& geo){
	return angular_momentum_flux_newtonian(s, geo.getTimeFrequency(3), geo.getEccentricity());
}

double carter_flux_newtonian(int s, GeodesicSource& geo){
	return carter_flux_newtonian(s, geo.getTimeFrequency(3), geo.getEccentricity(), geo.getInclination());
}

FluxList fluxes(GeodesicSource &geo){
	return fluxes(-2, geo);
}

FluxList fluxes(int s, GeodesicSource& geo){
	double efluxInf = 0.;
	double efluxHor = 0.;
	double efluxRef = energy_flux_newtonian(s, geo);

	double lfluxInf = 0.;
	double lfluxHor = 0.;
	double lfluxRef = angular_momentum_flux_newtonian(s, geo);

	double qfluxInf = 0.;
	double qfluxHor = 0.;
	double qfluxRef = carter_flux_newtonian(s, geo);

	int L = std::abs(s);
	FluxList fluxMode = flux_l(s, L, geo);
	efluxInf += fluxMode.Edot.infinity;
	efluxHor += fluxMode.Edot.horizon;
	lfluxInf += fluxMode.Ldot.infinity;
	lfluxHor += fluxMode.Ldot.horizon;
	qfluxInf += fluxMode.Qdot.infinity;
	qfluxHor += fluxMode.Qdot.horizon;

	double errorE = std::abs(fluxMode.Edot.infinity/efluxRef);
	double errorL = std::abs(fluxMode.Ldot.infinity/lfluxRef);
	double errorQ = std::abs(fluxMode.Qdot.infinity/qfluxRef);

	std::cout << "Edot_infinity ("<<L<<")-mode = "<< fluxMode.Edot.infinity << " \n";
	std::cout << "Edot_horizon ("<<L<<")-mode = "<< fluxMode.Edot.horizon << " \n";
	std::cout << "Edot_infinity sum = "<< efluxInf << " \n";
	std::cout << "Edot_horizon sum = "<< efluxHor << " \n";

	L++;
	while(L <= FLUX_LMAX && (errorE > FLUX_EPSILON || errorL > FLUX_EPSILON || errorQ > FLUX_EPSILON)){
		fluxMode = flux_l(s, L, geo);
		efluxInf += fluxMode.Edot.infinity;
		efluxHor += fluxMode.Edot.horizon;
		lfluxInf += fluxMode.Ldot.infinity;
		lfluxHor += fluxMode.Ldot.horizon;
		qfluxInf += fluxMode.Qdot.infinity;
		qfluxHor += fluxMode.Qdot.horizon;

		errorE = std::abs(fluxMode.Edot.infinity/efluxRef);
		errorL = std::abs(fluxMode.Ldot.infinity/lfluxRef);
		errorQ = std::abs(fluxMode.Qdot.infinity/qfluxRef);

		std::cout << "Edot_infinity ("<<L<<")-mode = "<< fluxMode.Edot.infinity << " \n";
		std::cout << "Edot_horizon ("<<L<<")-mode = "<< fluxMode.Edot.horizon << " \n";
		std::cout << "Edot_infinity sum = "<< efluxInf << " \n";
		std::cout << "Edot_horizon sum = "<< efluxHor << " \n";
		L++;
	}
	FluxList fluxes;
	fluxes.Edot.infinity = efluxInf;
	fluxes.Edot.horizon = efluxHor;
	fluxes.Ldot.infinity = lfluxInf;
	fluxes.Ldot.horizon = lfluxHor;
	fluxes.Qdot.infinity = qfluxInf;
	fluxes.Qdot.horizon = qfluxHor;

	return fluxes;
}

FluxList flux_l(int s, int L, GeodesicSource& geo){
	double efluxInf = 0.;
	double efluxHor = 0.;
	double efluxPreviousInf, efluxPreviousHor;
	double efluxRef = 0.5*energy_flux_newtonian(s, geo);

	double lfluxInf = 0.;
	double lfluxHor = 0.;
	double lfluxPreviousInf, lfluxPreviousHor;
	double lfluxRef = 0.5*angular_momentum_flux_newtonian(s, geo);

	double qfluxInf = 0.;
	double qfluxHor = 0.;
	double qfluxPreviousInf, qfluxPreviousHor;
	double qfluxRef = 0.5*carter_flux_newtonian(s, geo);
	if(qfluxRef == 0){
		qfluxRef = 1.;
	}

	FluxList fluxMode;
	int convergenceCriteria = 2;

	int m = L;
	fluxMode = flux_lm(s, L, m, geo);
	efluxPreviousInf = fluxMode.Edot.infinity;
	efluxPreviousHor = fluxMode.Edot.horizon;
	lfluxPreviousInf = fluxMode.Ldot.infinity;
	lfluxPreviousHor = fluxMode.Ldot.horizon;
	qfluxPreviousInf = fluxMode.Qdot.infinity;
	qfluxPreviousHor = fluxMode.Qdot.horizon;
	efluxInf += efluxPreviousInf;
	efluxHor += efluxPreviousHor;
	lfluxInf += lfluxPreviousInf;
	lfluxHor += lfluxPreviousHor;
	qfluxInf += qfluxPreviousInf;
	qfluxHor += qfluxPreviousHor;
	std::cout << "Edot_infinity ("<<L<<", "<<m<<")-mode = "<< fluxMode.Edot.infinity << " \n";
	// std::cout << "Edot_horizon ("<<L<<", "<<m<<")-mode = "<< fluxMode.Edot.horizon << " \n";
	std::cout << "Ldot_infinity ("<<L<<", "<<m<<")-mode = "<< fluxMode.Ldot.infinity << " \n";
	std::cout << "Qdot_infinity ("<<L<<", "<<m<<")-mode = "<< fluxMode.Qdot.infinity << " \n";
	m--;
	int convergenceCheckE = 0;
	int convergenceCheckL = 0;
	int convergenceCheckQ = 0;
	while(m > L - 2 && m > 0){
		fluxMode = flux_lm(s, L, m, geo);
		if(std::abs(fluxMode.Edot.infinity/efluxRef) < FLUX_EPSILON && (fluxMode.Edot.infinity < efluxPreviousInf || fluxMode.Edot.infinity < 1.e-18) && std::abs(fluxMode.Edot.horizon/efluxRef) < FLUX_EPSILON && (fluxMode.Edot.horizon < efluxPreviousHor || fluxMode.Edot.horizon < 1.e-18)){
			convergenceCheckE += 1;
		}else{
			convergenceCheckE = 0;
		}
		if(std::abs(fluxMode.Ldot.infinity/lfluxRef) < FLUX_EPSILON && (fluxMode.Ldot.infinity < lfluxPreviousInf || fluxMode.Ldot.infinity < 1.e-18) && std::abs(fluxMode.Ldot.horizon/lfluxRef) < FLUX_EPSILON && (fluxMode.Ldot.horizon < lfluxPreviousHor || fluxMode.Ldot.horizon < 1.e-18)){
			convergenceCheckL += 1;
		}else{
			convergenceCheckL = 0;
		}
		if(std::abs(fluxMode.Qdot.infinity/qfluxRef) < FLUX_EPSILON && (fluxMode.Qdot.infinity < qfluxPreviousInf || fluxMode.Qdot.infinity < 1.e-18) && std::abs(fluxMode.Qdot.horizon/qfluxRef) < FLUX_EPSILON && (fluxMode.Qdot.horizon < qfluxPreviousHor || fluxMode.Qdot.horizon < 1.e-18)){
			convergenceCheckQ += 1;
		}else{
			convergenceCheckQ = 0;
		}

		efluxPreviousInf = fluxMode.Edot.infinity;
		efluxPreviousHor = fluxMode.Edot.horizon;
		lfluxPreviousInf = fluxMode.Ldot.infinity;
		lfluxPreviousHor = fluxMode.Ldot.horizon;
		qfluxPreviousInf = fluxMode.Qdot.infinity;
		qfluxPreviousHor = fluxMode.Qdot.horizon;
		if(convergenceCheckE < convergenceCriteria){
			efluxInf += efluxPreviousInf;
			efluxHor += efluxPreviousHor;
		}
		if(convergenceCheckL < convergenceCriteria){
			lfluxInf += lfluxPreviousInf;
			lfluxHor += lfluxPreviousHor;
		}
		if(convergenceCheckQ < convergenceCriteria){
			qfluxInf += qfluxPreviousInf;
			qfluxHor += qfluxPreviousHor;
		}
		std::cout << "Edot_infinity ("<<L<<", "<<m<<")-mode = "<< fluxMode.Edot.infinity << " \n";
		// std::cout << "Edot_horizon ("<<L<<", "<<m<<")-mode = "<< fluxMode.Edot.horizon << " \n";
		std::cout << "Ldot_infinity ("<<L<<", "<<m<<")-mode = "<< fluxMode.Ldot.infinity << " \n";
		std::cout << "Qdot_infinity ("<<L<<", "<<m<<")-mode = "<< fluxMode.Qdot.infinity << " \n";
		m--;
	}

	int minM = -L;
	if(geo.getEccentricity() < DBL_EPSILON && 1. - std::abs(geo.getInclination()) < DBL_EPSILON){
		minM = 1;
	}

	while(m >= minM && (convergenceCheckE < convergenceCriteria || convergenceCheckL < convergenceCriteria || convergenceCheckQ < convergenceCriteria)){
		fluxMode = flux_lm(s, L, m, geo);
		if(std::abs(fluxMode.Edot.infinity/efluxRef) < FLUX_EPSILON && (fluxMode.Edot.infinity < efluxPreviousInf || fluxMode.Edot.infinity < 1.e-18) && std::abs(fluxMode.Edot.horizon/efluxRef) < FLUX_EPSILON && (fluxMode.Edot.horizon < efluxPreviousHor || fluxMode.Edot.horizon < 1.e-18)){
			convergenceCheckE += 1;
		}else{
			convergenceCheckE = 0;
		}
		if(std::abs(fluxMode.Ldot.infinity/lfluxRef) < FLUX_EPSILON && (fluxMode.Ldot.infinity < lfluxPreviousInf || fluxMode.Ldot.infinity < 1.e-18) && std::abs(fluxMode.Ldot.horizon/lfluxRef) < FLUX_EPSILON && (fluxMode.Ldot.horizon < lfluxPreviousHor || fluxMode.Ldot.horizon < 1.e-18)){
			convergenceCheckL += 1;
		}else{
			convergenceCheckL = 0;
		}
		if(std::abs(fluxMode.Qdot.infinity/qfluxRef) < FLUX_EPSILON && (fluxMode.Qdot.infinity < qfluxPreviousInf || fluxMode.Qdot.infinity < 1.e-18) && std::abs(fluxMode.Qdot.horizon/qfluxRef) < FLUX_EPSILON && (fluxMode.Qdot.horizon < qfluxPreviousHor || fluxMode.Qdot.horizon < 1.e-18)){
			convergenceCheckQ += 1;
		}else{
			convergenceCheckQ = 0;
		}

		efluxPreviousInf = fluxMode.Edot.infinity;
		efluxPreviousHor = fluxMode.Edot.horizon;
		lfluxPreviousInf = fluxMode.Ldot.infinity;
		lfluxPreviousHor = fluxMode.Ldot.horizon;
		qfluxPreviousInf = fluxMode.Qdot.infinity;
		qfluxPreviousHor = fluxMode.Qdot.horizon;
		if(convergenceCheckE < convergenceCriteria){
			efluxInf += efluxPreviousInf;
			efluxHor += efluxPreviousHor;
		}
		if(convergenceCheckL < convergenceCriteria){
			lfluxInf += lfluxPreviousInf;
			lfluxHor += lfluxPreviousHor;
		}
		if(convergenceCheckQ < convergenceCriteria){
			qfluxInf += qfluxPreviousInf;
			qfluxHor += qfluxPreviousHor;
		}
		std::cout << "Edot_infinity ("<<L<<", "<<m<<")-mode = "<< fluxMode.Edot.infinity << " \n";
		// std::cout << "Edot_horizon ("<<L<<", "<<m<<")-mode = "<< fluxMode.Edot.horizon << " \n";
		std::cout << "Ldot_infinity ("<<L<<", "<<m<<")-mode = "<< fluxMode.Ldot.infinity << " \n";
		std::cout << "Qdot_infinity ("<<L<<", "<<m<<")-mode = "<< fluxMode.Qdot.infinity << " \n";
		m--;
	}

	FluxList fluxes;
	fluxes.Edot.infinity = efluxInf;
	fluxes.Edot.horizon = efluxHor;
	fluxes.Ldot.infinity = lfluxInf;
	fluxes.Ldot.horizon = lfluxHor;
	fluxes.Qdot.infinity = qfluxInf;
	fluxes.Qdot.horizon = qfluxHor;

	return fluxes;
}

FluxList flux_lm(int s, int L, int m, GeodesicSource& geo){
	if(geo.getEccentricity() == 0. && std::abs(geo.getInclination()) == 1.){
		if(m == 0){
			FluxList fluxes;
			return fluxes;
		}
		TeukolskyMode teukMode(s, L, m, 0, 0, geo);
		return flux_mode(s, geo, teukMode);
	}else{
		return flux_lm_sum(s, L, m, geo);
	}
}

FluxList flux_lm_sum(int s, int L, int m, GeodesicSource& geo){
	double efluxInf = 0.;
	double efluxHor = 0.;
	double efluxPreviousInf, efluxPreviousHor;
	double efluxRef = 0.5*energy_flux_newtonian(s, geo);

	double lfluxInf = 0.;
	double lfluxHor = 0.;
	double lfluxPreviousInf, lfluxPreviousHor;
	double lfluxRef = 0.5*angular_momentum_flux_newtonian(s, geo);

	double qfluxInf = 0.;
	double qfluxHor = 0.;
	double qfluxPreviousInf, qfluxPreviousHor;
	double qfluxRef = 0.5*carter_flux_newtonian(s, geo);
	if(qfluxRef == 0){
		qfluxRef = 1.;
	}
	// factor of 0.5 is a fudge factor because Newtonian fluxes usually overestimate

	FluxList fluxMode;
	int convergenceCriteria = 10;

	if(std::abs(geo.getInclination()) == 1.){
		return flux_lmk(s, L, m, 0, geo);
	}

	int kPeak = L - m;
	int kMin = kPeak - 100;
	if(geo.getEccentricity() == 0.){
		// for circular orbits don't sum over k values that give you a negative frequency
		kMin = minimum_polar_harmonic_circ(m, geo);
		if(kPeak < kMin){
			kPeak = kMin + 1;
		}
	}
	int kMax = kPeak + 100;
	int PRINT_INF_MODES = 1;
	int PRINT_HOR_MODES = 0;

	int k = kPeak;
	fluxMode = flux_lmk(s, L, m, k, geo);
	efluxPreviousInf = fluxMode.Edot.infinity;
	efluxPreviousHor = fluxMode.Edot.horizon;
	lfluxPreviousInf = fluxMode.Ldot.infinity;
	lfluxPreviousHor = fluxMode.Ldot.horizon;
	qfluxPreviousInf = fluxMode.Qdot.infinity;
	qfluxPreviousHor = fluxMode.Qdot.horizon;
	efluxInf += efluxPreviousInf;
	efluxHor += efluxPreviousHor;
	lfluxInf += lfluxPreviousInf;
	lfluxHor += lfluxPreviousHor;
	qfluxInf += qfluxPreviousInf;
	qfluxHor += qfluxPreviousHor;
	if(PRINT_INF_MODES){
		std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<k<<")-mode = "<< efluxPreviousInf << " \n";
	}
	if(PRINT_HOR_MODES){
		std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<k<<")-mode = "<< efluxPreviousHor << " \n";
	}
	k++;

	int convergenceCheckE = 0;
	int convergenceCheckL = 0;
	int convergenceCheckQ = 0;
	while(k <= kMax && (convergenceCheckE < convergenceCriteria || convergenceCheckL < convergenceCriteria || convergenceCheckQ < convergenceCriteria)){
		fluxMode = flux_lmk(s, L, m, k, geo);
		if(std::abs(fluxMode.Edot.infinity/efluxRef) < FLUX_EPSILON/10 && (fluxMode.Edot.infinity < efluxPreviousInf || fluxMode.Edot.infinity < 1.e-18) && std::abs(fluxMode.Edot.horizon/efluxRef) < FLUX_EPSILON/10 && (fluxMode.Edot.horizon < efluxPreviousHor || fluxMode.Edot.horizon < 1.e-18)){
			convergenceCheckE += 1;
		}else{
			if(convergenceCheckE < COUNT_MAX){
				convergenceCheckE = 0;
			}
		}
		if(std::abs(fluxMode.Ldot.infinity/lfluxRef) < FLUX_EPSILON/10 && (fluxMode.Ldot.infinity < lfluxPreviousInf || fluxMode.Ldot.infinity < 1.e-18) && std::abs(fluxMode.Ldot.horizon/lfluxRef) < FLUX_EPSILON/10 && (fluxMode.Ldot.horizon < lfluxPreviousHor || fluxMode.Ldot.horizon < 1.e-18)){
			convergenceCheckL += 1;
		}else{
			if(convergenceCheckL < COUNT_MAX){
				convergenceCheckL = 0;
			}
		}
		if(std::abs(fluxMode.Qdot.infinity/qfluxRef) < FLUX_EPSILON/10 && (fluxMode.Qdot.infinity < qfluxPreviousInf || fluxMode.Qdot.infinity < 1.e-18) && std::abs(fluxMode.Qdot.horizon/qfluxRef) < FLUX_EPSILON/10 && (fluxMode.Qdot.horizon < qfluxPreviousHor || fluxMode.Qdot.horizon < 1.e-18)){
			convergenceCheckQ += 1;
		}else{
			if(convergenceCheckQ < COUNT_MAX){
				convergenceCheckQ = 0;
			}
		}
		efluxPreviousInf = fluxMode.Edot.infinity;
		efluxPreviousHor = fluxMode.Edot.horizon;
		lfluxPreviousInf = fluxMode.Ldot.infinity;
		lfluxPreviousHor = fluxMode.Ldot.horizon;
		qfluxPreviousInf = fluxMode.Qdot.infinity;
		qfluxPreviousHor = fluxMode.Qdot.horizon;
		if(convergenceCheckE < convergenceCriteria){
			efluxInf += efluxPreviousInf;
			efluxHor += efluxPreviousHor;
		}else{
			convergenceCheckE = COUNT_MAX;
		}
		if(convergenceCheckL < convergenceCriteria){
			lfluxInf += lfluxPreviousInf;
			lfluxHor += lfluxPreviousHor;
		}else{
			convergenceCheckL = COUNT_MAX;
		}
		if(convergenceCheckQ < convergenceCriteria){
			qfluxInf += qfluxPreviousInf;
			qfluxHor += qfluxPreviousHor;
		}else{
			convergenceCheckQ = COUNT_MAX;
		}
		if(PRINT_INF_MODES){
			std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<k<<")-mode = "<< efluxPreviousInf << " \n";
		}
		if(PRINT_HOR_MODES){
			std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<k<<")-mode = "<< efluxPreviousHor << " \n";
		}
		k++;
	}

	k = kPeak - 1;
	fluxMode = flux_lmk(s, L, m, k, geo);
	efluxPreviousInf = fluxMode.Edot.infinity;
	efluxPreviousHor = fluxMode.Edot.horizon;
	lfluxPreviousInf = fluxMode.Ldot.infinity;
	lfluxPreviousHor = fluxMode.Ldot.horizon;
	qfluxPreviousInf = fluxMode.Qdot.infinity;
	qfluxPreviousHor = fluxMode.Qdot.horizon;
	efluxInf += efluxPreviousInf;
	efluxHor += efluxPreviousHor;
	lfluxInf += lfluxPreviousInf;
	lfluxHor += lfluxPreviousHor;
	qfluxInf += qfluxPreviousInf;
	qfluxHor += qfluxPreviousHor;
	if(PRINT_INF_MODES){
		std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<k<<")-mode = "<< efluxPreviousInf << " \n";
	}
	if(PRINT_HOR_MODES){
		std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<k<<")-mode = "<< efluxPreviousHor << " \n";
	}
	k--;
	convergenceCheckE = 0;
	convergenceCheckL = 0;
	convergenceCheckQ = 0;
	while(k >= kMin && (convergenceCheckE < convergenceCriteria || convergenceCheckL < convergenceCriteria || convergenceCheckQ < convergenceCriteria)){
		fluxMode = flux_lmk(s, L, m, k, geo);
		if(std::abs(fluxMode.Edot.infinity/efluxRef) < FLUX_EPSILON/10 && (fluxMode.Edot.infinity < efluxPreviousInf || fluxMode.Edot.infinity < 1.e-18) && std::abs(fluxMode.Edot.horizon/efluxRef) < FLUX_EPSILON/10 && (fluxMode.Edot.horizon < efluxPreviousHor || fluxMode.Edot.horizon < 1.e-18)){
			convergenceCheckE += 1;
		}else{
			if(convergenceCheckE < COUNT_MAX){
				convergenceCheckE = 0;
			}
		}
		if(std::abs(fluxMode.Ldot.infinity/lfluxRef) < FLUX_EPSILON/10 && (fluxMode.Ldot.infinity < lfluxPreviousInf || fluxMode.Ldot.infinity < 1.e-18) && std::abs(fluxMode.Ldot.horizon/lfluxRef) < FLUX_EPSILON/10 && (fluxMode.Ldot.horizon < lfluxPreviousHor || fluxMode.Ldot.horizon < 1.e-18)){
			convergenceCheckL += 1;
		}else{
			if(convergenceCheckL < COUNT_MAX){
				convergenceCheckL = 0;
			}
		}
		if(std::abs(fluxMode.Qdot.infinity/qfluxRef) < FLUX_EPSILON/10 && (fluxMode.Qdot.infinity < qfluxPreviousInf || fluxMode.Qdot.infinity < 1.e-18) && std::abs(fluxMode.Qdot.horizon/qfluxRef) < FLUX_EPSILON/10 && (fluxMode.Qdot.horizon < qfluxPreviousHor || fluxMode.Qdot.horizon < 1.e-18)){
			convergenceCheckQ += 1;
		}else{
			if(convergenceCheckQ < COUNT_MAX){
				convergenceCheckQ = 0;
			}
		}
		efluxPreviousInf = fluxMode.Edot.infinity;
		efluxPreviousHor = fluxMode.Edot.horizon;
		lfluxPreviousInf = fluxMode.Ldot.infinity;
		lfluxPreviousHor = fluxMode.Ldot.horizon;
		qfluxPreviousInf = fluxMode.Qdot.infinity;
		qfluxPreviousHor = fluxMode.Qdot.horizon;
		if(convergenceCheckE < convergenceCriteria){
			efluxInf += efluxPreviousInf;
			efluxHor += efluxPreviousHor;
		}else{
			convergenceCheckE = COUNT_MAX;
		}
		if(convergenceCheckL < convergenceCriteria){
			lfluxInf += lfluxPreviousInf;
			lfluxHor += lfluxPreviousHor;
		}else{
			convergenceCheckL = COUNT_MAX;
		}
		if(convergenceCheckQ < convergenceCriteria){
			qfluxInf += qfluxPreviousInf;
			qfluxHor += qfluxPreviousHor;
		}else{
			convergenceCheckQ = COUNT_MAX;
		}
		if(PRINT_INF_MODES){
			std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<k<<")-mode = "<< efluxPreviousInf << " \n";
		}
		if(PRINT_HOR_MODES){
			std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<k<<")-mode = "<< efluxPreviousHor << " \n";
		}
		k--;
	}
	FluxList fluxes;
	fluxes.Edot.infinity = efluxInf;
	fluxes.Edot.horizon = efluxHor;
	fluxes.Ldot.infinity = lfluxInf;
	fluxes.Ldot.horizon = lfluxHor;
	fluxes.Qdot.infinity = qfluxInf;
	fluxes.Qdot.horizon = qfluxHor;

	return fluxes;
}

FluxList flux_lmk(int s, int L, int m, int k, GeodesicSource& geo){
	if(std::abs(geo.getInclination()) == 1. && std::abs(k) > 0){
		FluxList fluxes;
		return fluxes;
	}else if(geo.getEccentricity() == 0.){
		//std::cout << "(FLUXES) Restricting to spherical inclined orbits. \n";
		TeukolskyMode teukMode(s, L, m, k, 0, geo);
		return flux_mode(s, geo, teukMode);
	}else{
		return flux_lmk_sum(s, L, m, k, geo);
	}
}

FluxList flux_lmk_sum(int s, int L, int m, int k, GeodesicSource& geo){
	double efluxInf = 0.;
	double efluxHor = 0.;
	double efluxPreviousInf, efluxPreviousHor;
	double efluxRef = 0.5*energy_flux_newtonian(s, geo);

	double lfluxInf = 0.;
	double lfluxHor = 0.;
	double lfluxPreviousInf, lfluxPreviousHor;
	double lfluxRef = 0.5*angular_momentum_flux_newtonian(s, geo);

	double qfluxInf = 0.;
	double qfluxHor = 0.;
	double qfluxPreviousInf, qfluxPreviousHor;
	double qfluxRef = 0.5*carter_flux_newtonian(s, geo);
	if(qfluxRef == 0){
		qfluxRef = 1.;
	}

	FluxList fluxMode;
	int convergenceCriteria = 10;

	// guess where the energy peaks in the radial n-modes
	int n0 = radial_n_mode_max(L, geo.getEccentricity());
	if(geo.getInclination() < 0){
		n0 *= -1;
	}
	// give some width to the region we will initially evaluate
	// in case we are slightly off from the actual peak value
	int nWidth = 2;
	nWidth += std::abs(round(0.2*n0));
	int nInit = n0 - nWidth;

	// calculate the minimum radial harmonic that still gives us a postive frequency
	int nMin = minimum_radial_harmonic(m, k, geo);

	if(nMin >= nInit){
		nInit = nMin + 1;
	}
	int nMax = nInit + nWidth + 200;
	int PRINT_INF_MODES = 0;
	int PRINT_HOR_MODES = 0;

	int n = nInit;
	TeukolskyMode teukMode(s, L, m, k, n, geo);
	fluxMode = flux_mode(s, geo, teukMode);
	efluxPreviousInf = fluxMode.Edot.infinity;
	efluxPreviousHor = fluxMode.Edot.horizon;
	lfluxPreviousInf = fluxMode.Ldot.infinity;
	lfluxPreviousHor = fluxMode.Ldot.horizon;
	qfluxPreviousInf = fluxMode.Qdot.infinity;
	qfluxPreviousHor = fluxMode.Qdot.horizon;
	efluxInf += efluxPreviousInf;
	efluxHor += efluxPreviousHor;
	lfluxInf += lfluxPreviousInf;
	lfluxHor += lfluxPreviousHor;
	qfluxInf += qfluxPreviousInf;
	qfluxHor += qfluxPreviousHor;
	if(PRINT_INF_MODES){
		std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<k<<", "<<n<<")-mode = "<< fluxMode.Edot.infinity << " \n";
	}
	if(PRINT_HOR_MODES){
		std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<k<<", "<<n<<")-mode = "<< fluxMode.Edot.horizon << " \n";
	}
	n++;

	// search original region without testing for convergence to make sure we have calculate the peak value
	while(n <= nInit + nWidth){
		teukMode = TeukolskyMode(s, L, m, k, n, geo);
		fluxMode = flux_mode(s, geo, teukMode);
		efluxPreviousInf = fluxMode.Edot.infinity;
		efluxPreviousHor = fluxMode.Edot.horizon;
		lfluxPreviousInf = fluxMode.Ldot.infinity;
		lfluxPreviousHor = fluxMode.Ldot.horizon;
		qfluxPreviousInf = fluxMode.Qdot.infinity;
		qfluxPreviousHor = fluxMode.Qdot.horizon;
		efluxInf += efluxPreviousInf;
		efluxHor += efluxPreviousHor;
		lfluxInf += lfluxPreviousInf;
		lfluxHor += lfluxPreviousHor;
		qfluxInf += qfluxPreviousInf;
		qfluxHor += qfluxPreviousHor;
		if(PRINT_INF_MODES){
			std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<k<<", "<<n<<")-mode = "<< fluxMode.Edot.infinity << " \n";
		}
		if(PRINT_HOR_MODES){
			std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<k<<", "<<n<<")-mode = "<< fluxMode.Edot.horizon << " \n";
		}
		n++;
	}

	// start testing for convergence of n-mode sum up to nMax
	int convergenceCheckE = 0;
	int convergenceCheckL = 0;
	int convergenceCheckQ = 0;
	while(n <= nMax && (convergenceCheckE < convergenceCriteria || convergenceCheckL < convergenceCriteria || convergenceCheckQ < convergenceCriteria)){
		teukMode = TeukolskyMode(s, L, m, k, n, geo);
		fluxMode = flux_mode(s, geo, teukMode);
		if(std::abs(fluxMode.Edot.infinity/efluxRef) < FLUX_EPSILON/10 && (fluxMode.Edot.infinity < efluxPreviousInf || fluxMode.Edot.infinity < 1.e-18) && std::abs(fluxMode.Edot.horizon/efluxRef) < FLUX_EPSILON/10 && (fluxMode.Edot.horizon < efluxPreviousHor || fluxMode.Edot.horizon < 1.e-18)){
			convergenceCheckE += 1;
		}else{
			if(convergenceCheckE < COUNT_MAX){
				convergenceCheckE = 0;
			}
		}
		if(std::abs(fluxMode.Ldot.infinity/lfluxRef) < FLUX_EPSILON/10 && (fluxMode.Ldot.infinity < lfluxPreviousInf || fluxMode.Ldot.infinity < 1.e-18) && std::abs(fluxMode.Ldot.horizon/lfluxRef) < FLUX_EPSILON/10 && (fluxMode.Ldot.horizon < lfluxPreviousHor || fluxMode.Ldot.horizon < 1.e-18)){
			convergenceCheckL += 1;
		}else{
			if(convergenceCheckL < COUNT_MAX){
				convergenceCheckL = 0;
			}
		}
		if(std::abs(fluxMode.Qdot.infinity/qfluxRef) < FLUX_EPSILON/10 && (fluxMode.Qdot.infinity < qfluxPreviousInf || fluxMode.Qdot.infinity < 1.e-18) && std::abs(fluxMode.Qdot.horizon/qfluxRef) < FLUX_EPSILON/10 && (fluxMode.Qdot.horizon < qfluxPreviousHor || fluxMode.Qdot.horizon < 1.e-18)){
			convergenceCheckQ += 1;
		}else{
			if(convergenceCheckQ < COUNT_MAX){
				convergenceCheckQ = 0;
			}
		}
		efluxPreviousInf = fluxMode.Edot.infinity;
		efluxPreviousHor = fluxMode.Edot.horizon;
		lfluxPreviousInf = fluxMode.Ldot.infinity;
		lfluxPreviousHor = fluxMode.Ldot.horizon;
		qfluxPreviousInf = fluxMode.Qdot.infinity;
		qfluxPreviousHor = fluxMode.Qdot.horizon;
		if(convergenceCheckE < convergenceCriteria){
			efluxInf += efluxPreviousInf;
			efluxHor += efluxPreviousHor;
		}else{
			convergenceCheckE = COUNT_MAX;
		}
		if(convergenceCheckL < convergenceCriteria){
			lfluxInf += lfluxPreviousInf;
			lfluxHor += lfluxPreviousHor;
		}else{
			convergenceCheckL = COUNT_MAX;
		}
		if(convergenceCheckQ < convergenceCriteria){
			qfluxInf += qfluxPreviousInf;
			qfluxHor += qfluxPreviousHor;
		}else{
			convergenceCheckQ = COUNT_MAX;
		}
		if(PRINT_INF_MODES){
			std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<k<<", "<<n<<")-mode = "<< fluxMode.Edot.infinity << " \n";
		}
		if(PRINT_HOR_MODES){
			std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<k<<", "<<n<<")-mode = "<< fluxMode.Edot.horizon << " \n";
		}
		n++;
	}

	// start testing for convergence of n-mode sum down to nMin
	convergenceCheckE = 0;
	convergenceCheckL = 0;
	convergenceCheckQ = 0;
	n = nInit - 1;

	teukMode = TeukolskyMode(s, L, m, k, n, geo);
	fluxMode = flux_mode(s, geo, teukMode);
	efluxPreviousInf = fluxMode.Edot.infinity;
	efluxPreviousHor = fluxMode.Edot.horizon;
	lfluxPreviousInf = fluxMode.Ldot.infinity;
	lfluxPreviousHor = fluxMode.Ldot.horizon;
	qfluxPreviousInf = fluxMode.Qdot.infinity;
	qfluxPreviousHor = fluxMode.Qdot.horizon;
	efluxInf += efluxPreviousInf;
	efluxHor += efluxPreviousHor;
	lfluxInf += lfluxPreviousInf;
	lfluxHor += lfluxPreviousHor;
	qfluxInf += qfluxPreviousInf;
	qfluxHor += qfluxPreviousHor;
	if(PRINT_INF_MODES){
		std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<k<<", "<<n<<")-mode = "<< fluxMode.Edot.infinity << " \n";
	}
	if(PRINT_HOR_MODES){
		std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<k<<", "<<n<<")-mode = "<< fluxMode.Edot.horizon << " \n";
	}
	n--;

	while(n >= nMin && (convergenceCheckE < convergenceCriteria || convergenceCheckL < convergenceCriteria || convergenceCheckQ < convergenceCriteria)){
		teukMode = TeukolskyMode(s, L, m, k, n, geo);
		fluxMode = flux_mode(s, geo, teukMode);
		if(std::abs(fluxMode.Edot.infinity/efluxRef) < FLUX_EPSILON/10 && (fluxMode.Edot.infinity < efluxPreviousInf || fluxMode.Edot.infinity < 1.e-18) && std::abs(fluxMode.Edot.horizon/efluxRef) < FLUX_EPSILON/10 && (fluxMode.Edot.horizon < efluxPreviousHor || fluxMode.Edot.horizon < 1.e-18)){
			convergenceCheckE += 1;
		}else{
			if(convergenceCheckE < COUNT_MAX){
				convergenceCheckE = 0;
			}
		}
		if(std::abs(fluxMode.Ldot.infinity/lfluxRef) < FLUX_EPSILON/10 && (fluxMode.Ldot.infinity < lfluxPreviousInf || fluxMode.Ldot.infinity < 1.e-18) && std::abs(fluxMode.Ldot.horizon/lfluxRef) < FLUX_EPSILON/10 && (fluxMode.Ldot.horizon < lfluxPreviousHor || fluxMode.Ldot.horizon < 1.e-18)){
			convergenceCheckL += 1;
		}else{
			if(convergenceCheckL < COUNT_MAX){
				convergenceCheckL = 0;
			}
		}
		if(std::abs(fluxMode.Qdot.infinity/qfluxRef) < FLUX_EPSILON/10 && (fluxMode.Qdot.infinity < qfluxPreviousInf || fluxMode.Qdot.infinity < 1.e-18) && std::abs(fluxMode.Qdot.horizon/qfluxRef) < FLUX_EPSILON/10 && (fluxMode.Qdot.horizon < qfluxPreviousHor || fluxMode.Qdot.horizon < 1.e-18)){
			convergenceCheckQ += 1;
		}else{
			if(convergenceCheckQ < COUNT_MAX){
				convergenceCheckQ = 0;
			}
		}
		efluxPreviousInf = fluxMode.Edot.infinity;
		efluxPreviousHor = fluxMode.Edot.horizon;
		lfluxPreviousInf = fluxMode.Ldot.infinity;
		lfluxPreviousHor = fluxMode.Ldot.horizon;
		qfluxPreviousInf = fluxMode.Qdot.infinity;
		qfluxPreviousHor = fluxMode.Qdot.horizon;
		if(convergenceCheckE < convergenceCriteria){
			efluxInf += efluxPreviousInf;
			efluxHor += efluxPreviousHor;
		}else{
			convergenceCheckE = COUNT_MAX;
		}
		if(convergenceCheckL < convergenceCriteria){
			lfluxInf += lfluxPreviousInf;
			lfluxHor += lfluxPreviousHor;
		}else{
			convergenceCheckL = COUNT_MAX;
		}
		if(convergenceCheckQ < convergenceCriteria){
			qfluxInf += qfluxPreviousInf;
			qfluxHor += qfluxPreviousHor;
		}else{
			convergenceCheckQ = COUNT_MAX;
		}
		if(PRINT_INF_MODES){
			std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<k<<", "<<n<<")-mode = "<< fluxMode.Edot.infinity << " \n";
		}
		if(PRINT_HOR_MODES){
			std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<k<<", "<<n<<")-mode = "<< fluxMode.Edot.horizon << " \n";
		}
		n--;
	}
	FluxList fluxes;
	fluxes.Edot.infinity = efluxInf;
	fluxes.Edot.horizon = efluxHor;
	fluxes.Ldot.infinity = lfluxInf;
	fluxes.Ldot.horizon = lfluxHor;
	fluxes.Qdot.infinity = qfluxInf;
	fluxes.Qdot.horizon = qfluxHor;

	return fluxes;
}

double angular_momentum_flux_factor(int m, double omega){
	return m/omega;
}

double carter_flux_factor(GeodesicSource& geo, int m, int k, int n, double omega){
	double Mmkn = geo.getCarterFrequency(m, k, n);
	return 2.*Mmkn/omega;
}

double angular_momentum_flux_factor(GeodesicSource&, TeukolskyMode& teukMode){
	return angular_momentum_flux_factor(teukMode.getAzimuthalModeNumber(), teukMode.getFrequency());
}

double carter_flux_factor(GeodesicSource& geo, TeukolskyMode& teukMode){
	return carter_flux_factor(geo, teukMode.getAzimuthalModeNumber(), teukMode.getPolarModeNumber(), teukMode.getRadialModeNumber(), teukMode.getFrequency());
}

FluxList flux_mode(int s, GeodesicSource& geo, TeukolskyMode& teukMode, int include_minus_m){
	Fluxes Edot;
	double mFactor = 1.;
	if(include_minus_m == 1){
		mFactor = 2.;
	}
	if(s == 0){
		Edot = scalar_energy_flux_mode(geo, teukMode);
	}else{
		Edot = gravitational_energy_flux_mode(geo, teukMode);
	}
	double LdotPrefactor = 0.;
	double QdotPrefactor = 0.;
	if(std::abs(teukMode.getFrequency()) > 1.e-8){
		LdotPrefactor = angular_momentum_flux_factor(geo, teukMode);
		QdotPrefactor = carter_flux_factor(geo, teukMode);
	}
	FluxList fluxes;
	fluxes.Edot.infinity = mFactor*Edot.infinity;
	fluxes.Edot.horizon = mFactor*Edot.horizon;
	fluxes.Ldot.infinity = mFactor*LdotPrefactor*Edot.infinity;
	fluxes.Ldot.horizon = mFactor*LdotPrefactor*Edot.horizon;
	fluxes.Qdot.infinity = mFactor*QdotPrefactor*Edot.infinity;
	fluxes.Qdot.horizon = mFactor*QdotPrefactor*Edot.horizon;

	return fluxes;
}

Fluxes energy_flux_mode(int s, GeodesicSource& geo, TeukolskyMode& teukMode){
	if(s == 0){
		return scalar_energy_flux_mode(geo, teukMode);
	}else{
		return gravitational_energy_flux_mode(geo, teukMode);
	}
}

Fluxes scalar_energy_flux_mode(double a, double gamma, double omega, double horizonAmplitude, double infinityAmplitude){
	double prefactor = 0.25*omega/M_PI;
	// (r_+^2 + a^2)/M^2 = 2 r_+/M factor due to normalization at the horizon
	double horizonPrefactor = 2.*(1. + sqrt(1. - pow(a, 2)));
	Fluxes eflux = {prefactor*omega*infinityAmplitude, horizonPrefactor*prefactor*gamma*horizonAmplitude};
	return eflux;
}

Fluxes scalar_energy_flux_mode(GeodesicSource& geo, TeukolskyMode& teukMode){
	if(std::abs(teukMode.getFrequency()) > ZERO_FREQ_MAX){
		teukMode.generateSolutions(geo);
		return scalar_energy_flux_mode(geo.getBlackHoleSpin(), teukMode.getHorizonFrequency(), teukMode.getFrequency(), pow(std::abs(teukMode.getTeukolskyAmplitude(In)), 2), pow(std::abs(teukMode.getTeukolskyAmplitude(Up)), 2));
	}else{
		Fluxes eflux = {0., 0.};
		return eflux;
	}
}

Fluxes gravitational_energy_flux_mode(GeodesicSource& geo, TeukolskyMode& teukMode){
	if(std::abs(teukMode.getFrequency()) > ZERO_FREQ_MAX){
		teukMode.generateSolutions(geo);
		double prefactor = pow(2.*teukMode.getFrequency(), -2)/M_PI;
		Fluxes eflux = {prefactor*pow(std::abs(teukMode.getTeukolskyAmplitude(Up)), 2), horizonConstant(teukMode)*prefactor*pow(std::abs(teukMode.getTeukolskyAmplitude(In)), 2)};
		return eflux;
	}else{
		// std::cout << "Skipping (l, m, k, n) = (" << teukMode.getSpheroidalModeNumber() << ", " << teukMode.getAzimuthalModeNumber() << ", " << teukMode.getPolarModeNumber() << ", " << teukMode.getRadialModeNumber() << ")\n";
		Fluxes eflux = {0., 0.};
		return eflux;
	}
}

double horizonConstant(TeukolskyMode teuk){
	return horizonConstant(teuk.getBlackHoleSpin(), teuk.getSpheroidalModeNumber(), teuk.getAzimuthalModeNumber(), teuk.getFrequency(), teuk.getEigenvalue());
}

double horizonConstant(double a, int L, int m, double omega, double lambda){
	double rplus = 1. + sqrt(1. - a*a);
	double k = omega - m*a/(2.*rplus);
	double epsilon = (rplus - 1.)/(4.*rplus);
	double Alm = horizonAmplitude(a, L, m, omega, lambda);

	return 256.*pow(2.*rplus, 5)*k*(k*k + 4.*epsilon*epsilon)*(k*k + 16.*epsilon*epsilon)*pow(omega, 3)/Alm;
}

double horizonAmplitude(double a, int, int m, double omega, double lambda){
	return (pow(lambda + 2., 2) + 4*m*a*omega - pow(2*a*omega, 2))*(pow(lambda, 2) + 36.*m*a*omega - pow(6.*a*omega, 2))
		+ (2.*lambda + 3)*(96.*pow(a*omega, 2) - 48.*m*a*omega) + pow(12.*omega, 2)*(1. - a*a);
}

////////////////////////
/// ENERGY FLUX CASE ///
////////////////////////

Fluxes energy_flux(GeodesicSource& geo){
	return energy_flux(-2, geo);
}

Fluxes energy_flux(int s, GeodesicSource& geo){
	double efluxInf = 0.;
	double efluxHor = 0.;
	double efluxRef = energy_flux_newtonian(s, geo);

	int L = std::abs(s);
	Fluxes fluxMode = energy_flux_l(s, L, geo);
	efluxInf += fluxMode.infinity;
	efluxHor += fluxMode.horizon;
	double error = std::abs(fluxMode.infinity/efluxRef);
	// std::cout << "Edot_infinity ("<<L<<")-mode = "<< fluxMode.infinity << " \n";
	// std::cout << "Edot_horizon ("<<L<<")-mode = "<< fluxMode.horizon << " \n";
	// std::cout << "Edot_infinity sum = "<< efluxInf << " \n";
	// std::cout << "Edot_horizon sum = "<< efluxHor << " \n";
	L++;
	while(L <= FLUX_LMAX && error > FLUX_EPSILON){
		fluxMode = energy_flux_l(s, L, geo);
		efluxInf += fluxMode.infinity;
		efluxHor += fluxMode.horizon;
		error = std::abs(fluxMode.infinity/efluxRef);
		// std::cout << "Edot_infinity ("<<L<<")-mode = "<< fluxMode.infinity << " \n";
		// std::cout << "Edot_horizon ("<<L<<")-mode = "<< fluxMode.horizon << " \n";
		// std::cout << "Edot_infinity sum = "<< efluxInf << " \n";
		// std::cout << "Edot_horizon sum = "<< efluxHor << " \n";
		L++;
	}
	Fluxes eflux = {efluxInf, efluxHor};

	return eflux;
}

Fluxes energy_flux_alt(int s, GeodesicSource& geo){
	double efluxInf = 0.;
	double efluxHor = 0.;
	double efluxRef = energy_flux_newtonian(s, geo);

	int m = 0;
	Fluxes fluxMode = energy_flux_m(s, m, geo);
	efluxInf += fluxMode.infinity;
	efluxHor += fluxMode.horizon;
	double error = std::abs(fluxMode.infinity/efluxRef);
	// std::cout << "Edot_infinity ("<<m<<")-mode = "<< fluxMode.infinity << " \n";
	// std::cout << "Edot_horizon ("<<m<<")-mode = "<< fluxMode.horizon << " \n";
	// std::cout << "Edot_infinity sum = "<< efluxInf << " \n";
	// std::cout << "Edot_horizon sum = "<< efluxHor << " \n";
	m++;
	while(m <= FLUX_LMAX && error > FLUX_EPSILON){
		fluxMode = energy_flux_m(s, m, geo);
		efluxInf += fluxMode.infinity;
		efluxHor += fluxMode.horizon;
		error = std::abs(fluxMode.infinity/efluxRef);
		// std::cout << "Edot_infinity ("<<m<<")-mode = "<< fluxMode.infinity << " \n";
		// std::cout << "Edot_horizon ("<<m<<")-mode = "<< fluxMode.horizon << " \n";
		// std::cout << "Edot_infinity sum = "<< efluxInf << " \n";
		// std::cout << "Edot_horizon sum = "<< efluxHor << " \n";
		m++;
	}
	Fluxes eflux = {efluxInf, efluxHor};

	return eflux;
}

Fluxes energy_flux_l(int s, int L, GeodesicSource& geo){
	double efluxInf = 0.;
	double efluxHor = 0.;
	double efluxPreviousInf, efluxPreviousHor;
	double efluxRef = energy_flux_newtonian(s, geo);
	Fluxes fluxMode;
	int convergenceCriteria = 2;

	int m = L;
	fluxMode = energy_flux_lm(s, L, m, geo);
	// if(std::abs(m) > 0){
	// 	efluxInf += 2.*fluxMode.infinity;
	// 	efluxHor += 2.*fluxMode.horizon;
	// }else{
	// 	efluxInf += fluxMode.infinity;
	// 	efluxHor += fluxMode.horizon;
	// }
	efluxInf += 2.*fluxMode.infinity;
	efluxHor += 2.*fluxMode.horizon;
	// std::cout << "Edot_infinity ("<<L<<", "<<m<<")-mode = "<< fluxMode.infinity << " \n";
	// std::cout << "Edot_horizon ("<<L<<", "<<m<<")-mode = "<< fluxMode.horizon << " \n";
	m--;
	int convergenceCheck = 0;
	while(m > L - 2 && m > 0){
		efluxPreviousInf = 2.*fluxMode.infinity;
		efluxPreviousHor = 2.*fluxMode.horizon;
		fluxMode = energy_flux_lm(s, L, m, geo);
		if(std::abs(fluxMode.infinity/efluxRef) < FLUX_EPSILON && (fluxMode.infinity < efluxPreviousInf || fluxMode.infinity < 1.e-18) && std::abs(fluxMode.horizon/efluxRef) < FLUX_EPSILON && (fluxMode.horizon < efluxPreviousHor || fluxMode.horizon < 1.e-18)){
			convergenceCheck += 1;
		}else{
			convergenceCheck = 0;
		}
		efluxInf += 2.*fluxMode.infinity;
		efluxHor += 2.*fluxMode.horizon;
		// std::cout << "Edot_infinity ("<<L<<", "<<m<<")-mode = "<< fluxMode.infinity << " \n";
		// std::cout << "Edot_horizon ("<<L<<", "<<m<<")-mode = "<< fluxMode.horizon << " \n";
		m--;
	}

	int minM = -L;
	if(geo.getEccentricity() < DBL_EPSILON && 1. - std::abs(geo.getInclination()) < DBL_EPSILON){
		minM = 1;
	}

	while(m >= minM && convergenceCheck < convergenceCriteria){
		efluxPreviousInf = 2.*fluxMode.infinity;
		efluxPreviousHor = 2.*fluxMode.horizon;
		fluxMode = energy_flux_lm(s, L, m, geo);
		if(std::abs(fluxMode.infinity/efluxRef) < FLUX_EPSILON && (fluxMode.infinity < efluxPreviousInf || fluxMode.infinity < 1.e-18) && std::abs(fluxMode.horizon/efluxRef) < FLUX_EPSILON && (fluxMode.horizon < efluxPreviousHor || fluxMode.horizon < 1.e-18)){
			convergenceCheck += 1;
		}else{
			convergenceCheck = 0;
		}
		efluxInf += 2.*fluxMode.infinity;
		efluxHor += 2.*fluxMode.horizon;
		// std::cout << "Edot_infinity ("<<L<<", "<<m<<")-mode = "<< fluxMode.infinity << " \n";
		// std::cout << "Edot_horizon ("<<L<<", "<<m<<")-mode = "<< fluxMode.horizon << " \n";
		m--;
	}
	// if(m == 0 && convergenceCheck < convergenceCriteria){
	// 	fluxMode = energy_flux_lm(s, L, m, geo);
	// 	efluxInf += 2.*fluxMode.infinity;
	// 	efluxHor += 2.*fluxMode.horizon;
	// 	std::cout << "Edot_infinity ("<<L<<", "<<m<<")-mode = "<< fluxMode.infinity << " \n";
	// 	std::cout << "Edot_horizon ("<<L<<", "<<m<<")-mode = "<< fluxMode.horizon << " \n";
	// }

	Fluxes eflux = {efluxInf, efluxHor};

	return eflux;
}

Fluxes energy_flux_m(int s, int m, GeodesicSource& geo){
	double efluxInf = 0.;
	double efluxHor = 0.;
	double efluxRef = energy_flux_newtonian(s, geo);
	Fluxes fluxMode;
	int convergenceCriteria = 2;

	int L = std::abs(m) < std::abs(s) ? std::abs(s) : std::abs(m);

	fluxMode = energy_flux_lm(s, L, m, geo);
	efluxInf += fluxMode.infinity;
	efluxHor += fluxMode.horizon;
	// std::cout << "Edot_infinity ("<<L<<", "<<m<<")-mode = "<< fluxMode.infinity << " \n";
	// std::cout << "Edot_horizon ("<<L<<", "<<m<<")-mode = "<< fluxMode.horizon << " \n";
	L++;

	double efluxPreviousInf = fluxMode.infinity;
	double efluxPreviousHor = fluxMode.horizon;
	int convergenceCheck = 0;
	while(L < FLUX_LMAX && convergenceCheck < convergenceCriteria){
		fluxMode = energy_flux_lm(s, L, m, geo);
		if(std::abs(fluxMode.infinity/efluxRef) < FLUX_EPSILON && (fluxMode.infinity < efluxPreviousInf || fluxMode.infinity < 1.e-18) && std::abs(fluxMode.horizon/efluxRef) < FLUX_EPSILON && (fluxMode.horizon < efluxPreviousHor || fluxMode.horizon < 1.e-18)){
			convergenceCheck += 1;
		}else{
			convergenceCheck = 0;
		}
		efluxPreviousInf = fluxMode.infinity;
		efluxInf += efluxPreviousInf;
		efluxPreviousHor = fluxMode.horizon;
		efluxHor += efluxPreviousHor;
		// std::cout << "Edot_infinity ("<<L<<", "<<m<<")-mode = "<< fluxMode.infinity << " \n";
		// std::cout << "Edot_horizon ("<<L<<", "<<m<<")-mode = "<< fluxMode.horizon << " \n";
		L++;
	}
	Fluxes eflux = {efluxInf, efluxHor};

	return eflux;
}

Fluxes energy_flux_lm(int s, int L, int m, GeodesicSource& geo){
	if(geo.getEccentricity() == 0. && std::abs(geo.getInclination()) == 1.){
		if(m == 0){
			Fluxes energy_flux = {0., 0.};
			return energy_flux;
		}
		TeukolskyMode teukMode(s, L, m, 0, 0, geo);
		return energy_flux_mode(s, geo, teukMode);
	}else{
		return energy_flux_lm_sum(s, L, m, geo);
	}
}

Fluxes energy_flux_lm_sum(int s, int L, int m, GeodesicSource& geo){
	double efluxInf = 0.;
	double efluxHor = 0.;
	double efluxRef = energy_flux_newtonian(s, geo);
	double efluxPreviousInf = 0., efluxPreviousHor = 0.;
	Fluxes fluxMode;
	int convergenceCriteria = 3;
	if(std::abs(geo.getInclination()) == 1.){
		return energy_flux_lmk(s, L, m, 0, geo);
	}

	int kPeak = L - m;
	int kMin = kPeak - 100;
	if(geo.getEccentricity() == 0.){
		// for circular orbits don't sum over k values that give you a negative frequency
		kMin = minimum_polar_harmonic_circ(m, geo);
		if(kPeak < kMin){
			kPeak = kMin + 1;
		}
	}
	int kMax = kPeak + 100;

	int k = kPeak;
	fluxMode = energy_flux_lmk(s, L, m, k, geo);
	efluxPreviousInf = fluxMode.infinity;
	efluxInf += efluxPreviousInf;
	efluxPreviousHor = fluxMode.horizon;
	efluxHor += efluxPreviousHor;
	// std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<k<<")-mode = "<< efluxPreviousInf << " \n";
	// std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<k<<")-mode = "<< efluxPreviousHor << " \n";
	k++;

	int convergenceCheck = 0;
	while(k <= kMax && convergenceCheck < convergenceCriteria){
		fluxMode = energy_flux_lmk(s, L, m, k, geo);
		if(std::abs(fluxMode.infinity/efluxRef) < FLUX_EPSILON/10 && (fluxMode.infinity < efluxPreviousInf || fluxMode.infinity < 1.e-18) && std::abs(fluxMode.horizon/efluxRef) < FLUX_EPSILON/10 && (std::abs(fluxMode.horizon) < std::abs(efluxPreviousHor) || std::abs(fluxMode.horizon) < 1.e-18)){
			convergenceCheck += 1;
		}else{
			convergenceCheck = 0;
		}
		efluxPreviousInf = fluxMode.infinity;
		efluxInf += efluxPreviousInf;
		efluxPreviousHor = fluxMode.horizon;
		efluxHor += efluxPreviousHor;
		// std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<k<<")-mode = "<< efluxPreviousInf << " \n";
		// std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<k<<")-mode = "<< efluxPreviousHor << " \n";
		k++;
	}

	k = kPeak - 1;
	fluxMode = energy_flux_lmk(s, L, m, k, geo);
	efluxPreviousInf = fluxMode.infinity;
	efluxInf += efluxPreviousInf;
	efluxPreviousHor = fluxMode.horizon;
	efluxHor += efluxPreviousHor;
	// std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<k<<")-mode = "<< efluxPreviousInf << " \n";
	// std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<k<<")-mode = "<< efluxPreviousHor << " \n";
	k--;
	convergenceCheck = 0;
	while(k >= kMin && convergenceCheck < convergenceCriteria){
		fluxMode = energy_flux_lmk(s, L, m, k, geo);
		if(std::abs(fluxMode.infinity/efluxRef) < FLUX_EPSILON/10 && (fluxMode.infinity < efluxPreviousInf || fluxMode.infinity < 1.e-18) && std::abs(fluxMode.horizon/efluxRef) < FLUX_EPSILON/10 && (std::abs(fluxMode.horizon) < std::abs(efluxPreviousHor) || std::abs(fluxMode.horizon) < 1.e-18)){
			convergenceCheck += 1;
		}else{
			convergenceCheck = 0;
		}
		efluxPreviousInf = fluxMode.infinity;
		efluxInf += efluxPreviousInf;
		efluxPreviousHor = fluxMode.horizon;
		efluxHor += efluxPreviousHor;
		// std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<k<<")-mode = "<< efluxPreviousInf << " \n";
		// std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<k<<")-mode = "<< efluxPreviousHor << " \n";
		k--;
	}
	Fluxes eflux = {efluxInf, efluxHor};

	return eflux;
}

Fluxes energy_flux_lmk(int s, int L, int m, int k, GeodesicSource& geo){
	if(geo.getEccentricity() == 0. && std::abs(geo.getInclination()) == 1.){
		Fluxes flux = {0., 0.};
		return flux;
	}else if(geo.getEccentricity() == 0.){
		//std::cout << "(FLUXES) Restricting to spherical inclined orbits. \n";
		TeukolskyMode teukMode(s, L, m, k, 0, geo);
		return energy_flux_mode(s, geo, teukMode);
	}else{
		return energy_flux_lmk_sum(s, L, m, k, geo);
	}
}

Fluxes energy_flux_lmk_sum(int s, int L, int m, int k, GeodesicSource& geo){
	double efluxInf = 0.;
	double efluxHor = 0.;
	double efluxRef = energy_flux_newtonian(s, geo);
	double efluxPreviousInf = 0., efluxPreviousHor = 0.;
	Fluxes fluxMode;
	int convergenceCriteria = 5;

	// guess where the energy peaks in the radial n-modes
	int n0 = radial_n_mode_max(L, geo.getEccentricity());
	if(geo.getInclination() < 0){
		n0 *= -1;
	}
	// give some width to the region we will initially evaluate
	// in case we are slightly off from the actual peak value
	int nWidth = 2;
	nWidth += std::abs(round(0.2*n0));
	int nInit = n0 - nWidth;

	// calculate the minimum radial harmonic that still gives us a postive frequency
	int nMin = minimum_radial_harmonic(m, k, geo);
	if(nMin > nInit){
		nInit = nMin + 1;
	}
	int nMax = nInit + 2*nWidth + 200;

	int n = nInit;
	TeukolskyMode teukMode(s, L, m, k, n, geo);
	fluxMode = energy_flux_mode(s, geo, teukMode);
	efluxPreviousInf = fluxMode.infinity;
	efluxInf += efluxPreviousInf;
	efluxPreviousHor = fluxMode.horizon;
	efluxHor += efluxPreviousHor;
	// std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<k<<", "<<n<<")-mode = "<< fluxMode.infinity << " \n";
	// std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<k<<", "<<n<<")-mode = "<< fluxMode.horizon << " \n";
	n++;

	// search original region without testing for convergence to make sure we have calculate the peak value
	while(n <= nInit + 2*nWidth){
		teukMode = TeukolskyMode(s, L, m, k, n, geo);
		fluxMode = energy_flux_mode(s, geo, teukMode);
		efluxPreviousInf = fluxMode.infinity;
		efluxInf += efluxPreviousInf;
		efluxPreviousHor = fluxMode.horizon;
		efluxHor += efluxPreviousHor;
		// std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<k<<", "<<n<<")-mode = "<< fluxMode.infinity << " \n";
		// std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<k<<", "<<n<<")-mode = "<< fluxMode.horizon << " \n";
		n++;
	}

	// start testing for convergence of n-mode sum up to nMax
	int convergenceCheck = 0;
	while(n <= nMax && convergenceCheck < convergenceCriteria){
		teukMode = TeukolskyMode(s, L, m, k, n, geo);
		fluxMode = energy_flux_mode(s, geo, teukMode);
		if(std::abs(fluxMode.infinity/efluxRef) < FLUX_EPSILON/10 && (fluxMode.infinity < efluxPreviousInf || fluxMode.infinity < 1.e-18) && std::abs(fluxMode.horizon/efluxRef) < FLUX_EPSILON/10 && (std::abs(fluxMode.horizon) < std::abs(efluxPreviousHor) || std::abs(fluxMode.horizon) < 1.e-18)){
			convergenceCheck += 1;
		}else{
			convergenceCheck = 0;
		}
		efluxPreviousInf = fluxMode.infinity;
		efluxInf += efluxPreviousInf;
		efluxPreviousHor = fluxMode.horizon;
		efluxHor += efluxPreviousHor;
		// std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<k<<", "<<n<<")-mode = "<< fluxMode.infinity << " \n";
		// std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<k<<", "<<n<<")-mode = "<< fluxMode.horizon << " \n";
		n++;
	}

	// start testing for convergence of n-mode sum down to nMin
	convergenceCheck = 0;
	n = nInit - 1;
	while(n >= nMin && convergenceCheck < convergenceCriteria){
		teukMode = TeukolskyMode(s, L, m, k, n, geo);
		fluxMode = energy_flux_mode(s, geo, teukMode);
		if(std::abs(fluxMode.infinity/efluxRef) < FLUX_EPSILON/10 && (fluxMode.infinity < efluxPreviousInf || fluxMode.infinity < 1.e-18) && std::abs(fluxMode.horizon/efluxRef) < FLUX_EPSILON/10 && (std::abs(fluxMode.horizon) < std::abs(efluxPreviousHor) || std::abs(fluxMode.horizon) < 1.e-18)){
			convergenceCheck += 1;
		}else{
			convergenceCheck = 0;
		}
		efluxPreviousInf = fluxMode.infinity;
		efluxInf += efluxPreviousInf;
		efluxPreviousHor = fluxMode.horizon;
		efluxHor += efluxPreviousHor;
		// std::cout << "Edot_infinity ("<<L<<", "<<m<<", "<<k<<", "<<n<<")-mode = "<< fluxMode.infinity << " \n";
		// std::cout << "Edot_horizon ("<<L<<", "<<m<<", "<<k<<", "<<n<<")-mode = "<< fluxMode.horizon << " \n";
		n--;
	}
	Fluxes eflux = {efluxInf, efluxHor};

	return eflux;
}

//predictor for where the n-mode sum peaks in n
int radial_n_mode_max(int l, double e){
	if(l == 2){
		return round(exp(0.5 - 1.5*log(1. - e)));
	}else if(l == 3){
		return round(exp(1.0 - 1.5*log(1. - e)));
	}else if(l < 2){
  	return 0;
	}else{
		return round(exp(1.0 - 1.5*log(1. - e)))*(l - 2) - round(exp(0.5 - 1.5*log(1. - e)))*(l - 3);
	}
}
