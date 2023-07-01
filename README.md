# pybhpt

A python package for solving problems in black hole perturbation theory

`pybhpt` is a collection of numerical tools for analyzing perturbations of Kerr spacetime, particularly the self-forces experienced by small bodies moving in a Kerr background. Subpackages include: 

- `pybhpt.geodesic`: a module that generates bound timelike geodesics in Kerr spacetime
- `pybhpt.radial`: a module that calculates homogeneous solutions of the radial Teukolsky equation
- `pybhpt.swsh`: a module that constructs the spin-weighted spheroidal harmonics
- `pybhpt.teuk`: a module that evaluates the inhomogeneous solutions (Teukolsky amplitudes) of the radial Teukolsky equation due to a point-particle on a bound timelike Kerr geodesic
- `pybhpt.flux`: a module that produces the gravitational wave fluxes sourced by a point-particle on a generic bound timelike Kerr geodesic
- `pybhpt.hertz`: a module that solves for the Hertz potentials for the CCK and AAB metric reconstruction procedures
- `pybhpt.metric`: a module that produces the coefficients needed to reconstruct the metric from the Hertz potentials
- `pybhpt.redshift`: a module that computes the generalized Detweiler redshift invariant in a variety of gauges

> **Warning**
> This project is still under heavy active development. Please use at your own risk. If you plan on using this code for research purposes, please reach out to the authors first.

## Installation

PyBHPT relies on several dependencies to install and run, namely a C/C++ compiler (e.g., `g++`), `gsl`, `BLAS`, `boost`, `Cython`, `numpy`, and `python >= 3.7`.
To reduce package conflicts and ensure that the proper dependencies are installed,
we recommend using Anaconda and its virtual environments.

Create a conda environment `pybhpt-env` (or whatever name you would like)
with the necessary dependencies to install `pybhpt`:
```
$ conda create -n pybhpt-env -c conda-forge gsl boost-cpp Cython numpy scipy python=3.7
$ conda activate pybhpt-env
```
To include the necessary compiler on MACOSX run:
```
$ conda install clang_osx-64 clangxx_osx-64
```
To include the necessary compiler on linux run:
```
$ conda install gcc_linux-64 gxx_linux-64
```
Next clone the :code:`pybhpt` repository from GitHub:
```
$ git clone https://github.com/znasipak/pybhpt.git
$ cd pybhpt
```
Finally, we recommend installing the package via `pip`:
```
$ pip install .
```

## Conda Environments with Jupyter

To run the code in a Jupyter notebook, we recommend `pip` installing the following dependencies:
```
$ pip install ipykernel matplotlib
```
One can then make the environment accessible within Jupyter by running
```
$ python -m ipykernel install --user --name=pybhpt-env
```

## Uninstalling

If the package is installed using `pip`, then one can easily uninstall the package by executing
```
$ pip uninstall pybhpt
```
To clean the repository, one will also need to remove the directories `build/` and `pybhpt.egg.info/` from the main repository 

## Authors

Zachary Nasipak