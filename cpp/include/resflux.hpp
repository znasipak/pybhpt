#ifndef RESFLUXES_HPP
#define RESFLUXES_HPP

#include "fluxes.hpp"

struct ResonantFluxList{
  ResonantFluxList(): Edot(9), Ldot(9), Qdot(9) {}
	std::vector<Fluxes> Edot;
	std::vector<Fluxes> Ldot;
	std::vector<Fluxes> Qdot;
};

struct FieldAmplitudes{
  FieldAmplitudes() {}
  Complex up;
  Complex in;
};

typedef std::vector<FieldAmplitudes> FieldAmplitudeVector;

Fluxes flux_mode(int s, double a, int m, double omega, double lambda, double horizonAmplitude, double infinityAmplitude);
ResonantFluxList res_flux_l(int s, int L, int nth, int nr, GeodesicSource& geo);
ResonantFluxList res_flux_lm(int s, int L, int m, int nth, int nr, GeodesicSource& geo);
ResonantFluxList res_flux_lmN(int s, int L, int m, int Nres, int nth, int nr, GeodesicSource& geo);
FieldAmplitudeVector resonant_flux_mode_coefficients(int s, int L, int m, int k0, int n0, int nth, int nr, SpinWeightedHarmonic& Slm, RadialTeukolsky& Rt, GeodesicSource& geo);

#endif // RESFLUXES_HPP
