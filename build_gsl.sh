#!/bin/bash
set -e  # stop on first error

# Version and install directory
GSL_VERSION=2.8
INSTALL_DIR="${INSTALL_DIR:-/opt/gsl-install}"

# Detect number of cores (Linux vs macOS)
if command -v nproc >/dev/null; then
    CORES=$(nproc)
else
    CORES=$(sysctl -n hw.ncpu)
fi

# Download source
mkdir -p /tmp/gsl-src
cd /tmp/gsl-src
curl -LO https://ftp.gnu.org/gnu/gsl/gsl-${GSL_VERSION}.tar.gz
tar -xzf gsl-${GSL_VERSION}.tar.gz
cd gsl-${GSL_VERSION}

# Build & install
echo "Configuring GSL ${GSL_VERSION}..."
./configure --prefix="${INSTALL_DIR}" > /dev/null

echo "Building with ${CORES} cores..."
make -j"${CORES}" > /dev/null

echo "Installing..."
make install > /dev/null

echo "GSL ${GSL_VERSION} installed to ${INSTALL_DIR}"
