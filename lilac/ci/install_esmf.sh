#!/usr/bin/env bash
set -e
set -x

cd external/esmf

export FC="gfortran"
export ESMF_COMPILER="gfortran"

export ESMF_DIR=$PWD
export ESMF_INSTALL_PREFIX=/usr/
make -j4 lib
make install
export ESMFMKFILE=${ESMF_INSTALL_PREFIX}/esmf.mk

cd -