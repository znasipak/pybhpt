#!/bin/bash

GSL_VERSION=2.8
INSTALL_DIR="/opt/gsl-install"

mkdir -p /tmp/gsl-src
cd /tmp/gsl-src
curl -LO https://ftp.gnu.org/gnu/gsl/gsl-${GSL_VERSION}.tar.gz
tar -xzf gsl-${GSL_VERSION}.tar.gz
cd gsl-${GSL_VERSION}

echo "Configuring..."
./configure --prefix="${INSTALL_DIR}" > /dev/null

echo "Building..."
make -j$(nproc) > /dev/null

echo "Installing..."
make install > /dev/null

echo "GSL ${GSL_VERSION} installed to ${INSTALL_DIR}"
