#!/usr/bin/env bash

set -e
set -x

cd ${HOME}

git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git ${HOME}/deps/pfunit
cd deps/pfunit

# set environemnt variables
export F90=gfortran
export F90_VENDOR=GNU

mkdir build
cd build
cmake ..
make install INSTALL_DIR=/usr

cd ${TRAVIS_BUILD_DIR}
