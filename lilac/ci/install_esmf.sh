#!/usr/bin/env bash
set -e
set -x

cd ./external/esmf

export FC="mpif90"

export ESMF_DIR=$PWD
export ESMF_COMM="mpich3"
export ESMF_COMPILER="gfortran"
export ESMF_INSTALL_PREFIX="/usr/local"
export ESMF_INSTALL_LIBDIR="/usr/local/lib"
export ESMF_INSTALL_MODDIR="/usr/local/mod"
export ESMF_INSTALL_BINDIR="/usr/local/bin"
export ESMF_INSTALL_DOCDIR="/usr/local/doc"
export ESMFMKFILE="${ESMF_INSTALL_LIBDIR}/esmf.mk"

make -j4 lib
make install
# make install check

cd -