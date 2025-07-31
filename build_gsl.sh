#!/bin/bash
set -euo pipefail

GSL_VERSION=2.7.1
INSTALL_PREFIX="$(pwd)/gsl-install"

echo "Downloading GSL ${GSL_VERSION}..."
curl -fsSL -O https://ftp.gnu.org/gnu/gsl/gsl-${GSL_VERSION}.tar.gz
tar -xzf gsl-${GSL_VERSION}.tar.gz
cd gsl-${GSL_VERSION}

echo "Configuring GSL..."
./configure --prefix="${INSTALL_PREFIX}" --disable-shared --enable-static
make -j$(nproc || sysctl -n hw.ncpu)
make install

# Export so CMake and compiler can find it
export CMAKE_PREFIX_PATH="${INSTALL_PREFIX}:${CMAKE_PREFIX_PATH:-}"
export CMAKE_INCLUDE_PATH="${INSTALL_PREFIX}/include:${CMAKE_INCLUDE_PATH:-}"
export CMAKE_LIBRARY_PATH="${INSTALL_PREFIX}/lib:${CMAKE_LIBRARY_PATH:-}"
export CPATH="${INSTALL_PREFIX}/include:${CPATH:-}"
export LIBRARY_PATH="${INSTALL_PREFIX}/lib:${LIBRARY_PATH:-}"

echo "GSL installed to ${INSTALL_PREFIX}"

#!/usr/bin/env bash
set -e

# Define installation prefix
INSTALL_DIR="$(pwd)/gsl-install"

# Build and install GSL
curl -LO https://ftp.gnu.org/gnu/gsl/gsl-2.7.tar.gz
tar -xzf gsl-2.7.tar.gz
cd gsl-2.7
./configure --prefix="${INSTALL_DIR}"
make -j$(nproc)
make install

# Output CMAKE_PREFIX_PATH for later use
echo "GSL installed to ${INSTALL_DIR}"
echo "CMAKE_PREFIX_PATH=${INSTALL_DIR}" >> $GITHUB_ENV
