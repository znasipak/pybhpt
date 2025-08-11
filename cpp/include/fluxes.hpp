// fluxes.hpp

#ifndef FLUXES_HPP
#define FLUXES_HPP

#include "teukolsky.hpp"

class Fluxes{
public:
	Fluxes(): infinity(0.), horizon(0.) {}
	Fluxes(double inf, double hor): infinity(inf), horizon(hor) {}
	double infinity;
	double horizon;
};

class FluxList{
public:
	FluxList() {}
	Fluxes Edot;
	Fluxes Ldot;
	Fluxes Qdot;
};

double energy_flux_newtonian(GeodesicSource& geo);
double energy_flux_newtonian(int s, double omegaPhiCirc);

double energy_flux_newtonian(int s, GeodesicSource& geo);
double angular_momentum_flux_newtonian(int s, GeodesicSource& geo);
double carter_flux_newtonian(int s, GeodesicSource& geo);

FluxList fluxes(GeodesicSource& geo);
FluxList fluxes(int s, GeodesicSource& geo);
FluxList flux_l(int s, int L, GeodesicSource& geo);
FluxList flux_lm(int s, int L, int m, GeodesicSource& geo);
FluxList flux_lm_sum(int s, int L, int m, GeodesicSource& geo);
FluxList flux_lmk(int s, int L, int m, int k, GeodesicSource& geo);
FluxList flux_lmk_sum(int s, int L, int m, int k, GeodesicSource& geo);
double angular_momentum_flux_factor(GeodesicSource&, TeukolskyMode& teukMode);
double carter_flux_factor(GeodesicSource& geo, TeukolskyMode& teukMode);
FluxList flux_mode(int s, GeodesicSource& geo, TeukolskyMode& teukMode, int include_minus_m = 1);

Fluxes energy_flux(GeodesicSource& geo);
Fluxes energy_flux(int s, GeodesicSource& geo);
Fluxes energy_flux_l(int s, int L, GeodesicSource& geo);
Fluxes energy_flux_m(int s, int m, GeodesicSource& geo);
Fluxes energy_flux_lm(int s, int L, int m, GeodesicSource& geo);
Fluxes energy_flux_lm_sum(int s, int L, int m, GeodesicSource& geo);
Fluxes energy_flux_lmk(int s, int L, int m, int k, GeodesicSource& geo);
Fluxes energy_flux_lmk_sum(int s, int L, int m, int k, GeodesicSource& geo);
Fluxes energy_flux_mode(int s, GeodesicSource& geo, TeukolskyMode& teukMode);

Fluxes scalar_energy_flux_mode(GeodesicSource& geo, TeukolskyMode& teukMode);
Fluxes gravitational_energy_flux_mode(GeodesicSource& geo, TeukolskyMode& teukMode);

double horizonConstant(TeukolskyMode teuk);
double horizonConstant(double a, int L, int m, double omega, double lambda);
double horizonAmplitude(double a, int L, int m, double omega, double lambda);
int radial_n_mode_max(int l, double e);

#endif
