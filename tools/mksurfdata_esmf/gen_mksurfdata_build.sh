#! /bin/bash -f

export MACH="cheyenne"
export COMPILER="intel"
export MPILIB="mpt"

cwd=`pwd`
rm -rf bld
mkdir bld
cd bld
../../../cime/tools/configure --macros-format CMake --machine $MACH --compiler $COMPILER --mpilib $MPILIB
ls -l
source ./.env_mach_specific.sh
CC=mpicc FC=mpif90 cmake -DCMAKE_BUILD_TYPE=debug ../src
make VERBOSE=1

