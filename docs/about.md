# About

`pybhpt` is a collection of numerical tools for analyzing perturbations of Kerr spacetime, particularly the self-forces and metric-perturbations experienced by small bodies moving in a Kerr background.

## Subpackages
- [`pybhpt.geo`](pybhpt.geo): generates bound periodic timelike geodesics in Kerr spacetime
- [`pybhpt.radial`](pybhpt.radial): calculates homogeneous solutions of the radial Teukolsky equation
- [`pybhpt.swsh`](pybhpt.swsh): constructs the spin-weighted spheroidal harmonics
- [`pybhpt.teuk`](pybhpt.teuk): evaluates inhomogeneous solutions (Teukolsky amplitudes) of the radial Teukolsky equation due to a point-particle on a bound timelike Kerr geodesic
- [`pybhpt.flux`](pybhpt.flux): produces gravitational wave fluxes sourced by a point-particle on a generic bound timelike Kerr geodesic
- [`pybhpt.hertz`](pybhpt.hertz): solves for the Hertz potentials for the CCK and AAB metric reconstruction procedures
- [`pybhpt.metric`](pybhpt.metric): produces coefficients needed to reconstruct the metric from the Hertz potentials
- [`pybhpt.redshift`](pybhpt.redshift): computes the generalized Detweiler redshift invariant in a variety of gauges

One can find out more information about each module by exploring the User Guides or clicking on the subpackages, which are linked to the API.

## References

Theoretical background for the code and explanations of the numerical methods used within are summarized in the references below: 

- Z. Nasipak, *Metric reconstruction and the Hamiltonian for eccentric, precessing binaries in the small-mass-ratio limit* (2025) [arXiv:2507.07746](https://arxiv.org/abs/2507.07746)
- Z. Nasipak, *An adiabatic gravitational waveform model for compact objects undergoing quasi-circular inspirals into rotating massive black holes*, Phys. Rev. D 109, 044020 (2024) [arXiv:2310.19706](https://arxiv.org/abs/2310.19706)
- Z. Nasipak, *Adiabatic evolution due to the conservative scalar self-force during orbital resonances*, Phys. Rev. D 106, 064042 (2022) [arXiv:2207.02224](https://arxiv.org/abs/2207.02224)

## Authors
Zachary Nasipak