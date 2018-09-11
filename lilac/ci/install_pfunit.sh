#!/usr/bin/env bash

set -e
set -x

cd ./external/pfunit

# set environemnt variables
export F90=gfortran
export F90_VENDOR=GNU

mkdir -p build
cd build
cmake ..
make install INSTALL_DIR=/usr/pfunit

cd -
