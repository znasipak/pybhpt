#!/bin/bash
set -euo pipefail

GSL_VERSION=2.7.1
INSTALL_DIR="$(pwd)/gsl-install"

echo "Downloading GSL ${GSL_VERSION}..."
curl -fsSL -O https://ftp.gnu.org/gnu/gsl/gsl-${GSL_VERSION}.tar.gz
tar -xzf gsl-${GSL_VERSION}.tar.gz
cd gsl-${GSL_VERSION}

echo "Configuring GSL..."
./configure --prefix="${INSTALL_DIR}" --disable-shared --enable-static
make -j$(nproc || sysctl -n hw.ncpu)
make install

# Output CMAKE_PREFIX_PATH for later use
echo "GSL installed to ${INSTALL_DIR}"
echo "CMAKE_PREFIX_PATH=${INSTALL_DIR}" >> $GITHUB_ENV

