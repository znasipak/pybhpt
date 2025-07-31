#!/bin/bash

wget https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz
tar -xzf gsl-2.7.1.tar.gz
cd gsl-2.7.1
./configure --prefix=/usr/local
make -j$(nproc)
sudo make install

