#!/bin/bash

GSL_VERSION=2.8
INSTALL_DIR="/opt/gsl-${GSL_VERSION}"

mkdir -p /tmp/gsl-src
cd /tmp/gsl-src
curl -LO https://ftp.gnu.org/gnu/gsl/gsl-${GSL_VERSION}.tar.gz
tar -xzf gsl-${GSL_VERSION}.tar.gz
cd gsl-${GSL_VERSION}
./configure --prefix="${INSTALL_DIR}"
make -j$(nproc)
make install

echo "CMAKE_PREFIX_PATH=${INSTALL_DIR}" >> $GITHUB_ENV


