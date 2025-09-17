# Welcome to pybhpt's documentation

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

One can find out more information about each module by exploring the [User Guides](#contents) or clicking on the subpackages, which are linked to the API. References and author information are provided in the [About](about) section.

## Quick Installation

Tagged releases of `pybhpt` are available as wheel packages for macOS and 64-bit Linux on [PyPI](https://pypi.org/project/pybhpt). Install using `pip`:
```
python3 -m pip install pybhpt
```
Developers can compile from source using the instructions in [Installation](installation). User guides and documentation are provided below.

## Contents

```{toctree}
:maxdepth: 1
:caption: User Guides

about
installation
notebooks/geodesics.ipynb
notebooks/radial.ipynb
notebooks/swsh.ipynb
notebooks/teuk.ipynb
notebooks/fluxes.ipynb
notebooks/waveform.ipynb
notebooks/tutorial.ipynb
```

```{toctree}
:maxdepth: 1
:caption: Background

background/geo
background/teuk
background/radial
background/swsh
```

```{toctree}
:maxdepth: 1
:caption: API Reference

pybhpt.geo
pybhpt.radial
pybhpt.swsh
pybhpt.teuk
pybhpt.hertz
pybhpt.flux
pybhpt.metric
pybhpt.redshift
```

## Authors
Zachary Nasipak