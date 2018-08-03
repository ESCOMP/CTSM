#!/usr/bin/env bash
set -e
set -x

cd ${HOME}

export FC="/usr/bin/gfortran"
export ESMF_COMPILER="/usr/bin/gfortran"

git clone https://git.code.sf.net/p/esmf/esmf deps/esmf
cd deps/esmf
export ESMF_DIR=$PWD
export ESMF_INSTALL_PREFIX=/usr/
make -j4 lib
make install
export ESMFMKFILE=${ESMF_INSTALL_PREFIX}/esmf.mk

cd ${TRAVIS_BUILD_DIR}
