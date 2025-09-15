# pybhpt

A python package for solving problems in black hole perturbation theory. 

`pybhpt` is a collection of numerical tools for analyzing perturbations of Kerr spacetime, particularly the self-forces and metric-perturbations experienced by small bodies moving in a Kerr background. Subpackages include: 

- `pybhpt.geodesic`: a module that generates bound periodic timelike geodesics in Kerr spacetime
- `pybhpt.radial`: a module that calculates homogeneous solutions of the radial Teukolsky equation
- `pybhpt.swsh`: a module that constructs the spin-weighted spheroidal harmonics
- `pybhpt.teuk`: a module that evaluates the inhomogeneous solutions (Teukolsky amplitudes) of the radial Teukolsky equation due to a point-particle on a bound timelike Kerr geodesic
- `pybhpt.flux`: a module that produces the gravitational wave fluxes sourced by a point-particle on a generic bound timelike Kerr geodesic
- `pybhpt.hertz`: a module that solves for the Hertz potentials for the CCK and AAB metric reconstruction procedures
- `pybhpt.metric`: a module that produces the coefficients needed to reconstruct the metric from the Hertz potentials
- `pybhpt.redshift`: a module that computes the generalized Detweiler redshift invariant in a variety of gauges

See the [Documentation](https://pybhpt.readthedocs.io/en/latest/) pages for more information about the package, including User Guides and API. References and author information can be found at the bottom of the README.

## Quick Installation

Tagged releases of `pybhpt` are available as wheel packages for macOS and 64-bit Linux on [PyPI](https://pypi.org/project/pybhpt). Install using `pip`:
```
python3 -m pip install pybhpt
```

## Developer Installation

`pybhpt` relies on several dependencies to install and run, namely a C/C++ compiler (e.g., `g++`), `gsl`, `boost`, `Cython`, `numpy`, `scipy`, and `python >= 3.8`. To reduce package conflicts and ensure that the proper dependencies are installed, we recommend using [conda](https://docs.conda.io/en/latest/) (paricularly through [Miniforge](https://github.com/conda-forge/miniforge)) and its virtual environments.

To create a conda environment `pybhpt-env` with just the necessary dependencies use `environment.yml`:
```
conda env create -f environment.yml
conda activate pybhpt-env
```
For an environment with the extended recommended software dependencies, one can replace `environment.yml` with `environment-extended.yml` or follow the extended `pip` install instructions below. 

Next clone the `pybhpt` repository from GitHub:
```
git clone https://github.com/znasipak/pybhpt.git
cd pybhpt
```
The `boost` repository is included as a submodule under `extern/boost`. To activate the submodule:
```
git submodule update --init --recursive extern/boost
```
To include the Boost headers in the compilation path:
```
cd extern/boost
chmod +x bootstrap.sh
./bootstrap.sh
./b2 headers
```
Finally, we recommend installing the package via `pip`:
```
pip install .
```
or to get the basic development dependencies:
```
pip install '.[dev, testing]'
```
To get all recommended software:
```
pip install '.[dev, testing, extended]'
```
To ensure the package installed properly, run the unit tests via `pytest`:
```
pytest .
```

## Conda Environments with Jupyter

To run the code in a Jupyter notebook, we recommend `pip` installing the following dependencies:
```
pip install ipykernel matplotlib
```
These dependencies are already included in `environment-extended.yml` or in the optional `pip install` dependencies. To make the environment accessible within Jupyter:
```
python -m ipykernel install --user --name=pybhpt-env
```

## Uninstalling

If the package is installed using `pip`, then one can easily uninstall the package by executing
```
pip uninstall pybhpt
```
Developers may also need to remove the directories `build/` and `pybhpt.egg.info/` from the main repository.

## Troubleshooting compiling from source

If there are problems compiling `pybhpt` from source, it may be because `cmake` cannot identify the correct compiler. To fix this issue, one can try to explicitly download the necessary compiler into their conda environment.

To include the necessary compiler on macOS Intel:
```
conda install clang_osx-64 clangxx_osx-64
```
To include the necessary compiler on macOS arm/silicon:
```
conda install clang_osx-arm64 clangxx_osx-arm64
```
To include the necessary compiler on Linux:
```
conda install gcc_linux-64 gxx_linux-64
```

## References

Theoretical background for the code and explanations of the numerical methods used within are summarized in the references below: 

- Z. Nasipak, *Metric reconstruction and the Hamiltonian for eccentric, precessing binaries in the small-mass-ratio limit* (2025) [arXiv:2507.07746](https://arxiv.org/abs/2507.07746)
- Z. Nasipak, *An adiabatic gravitational waveform model for compact objects undergoing quasi-circular inspirals into rotating massive black holes*, Phys. Rev. D 109, 044020 (2024) [arXiv:2310.19706](https://arxiv.org/abs/2310.19706)
- Z. Nasipak, *Adiabatic evolution due to the conservative scalar self-force during orbital resonances*, Phys. Rev. D 106, 064042 (2022) [arXiv:2207.02224](https://arxiv.org/abs/2207.02224)

## Authors

Zachary Nasipak