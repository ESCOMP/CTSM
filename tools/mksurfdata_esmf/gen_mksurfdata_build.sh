#! /bin/bash -f

# Define what machine to use that's been ported to cime
export MACH="cheyenne"
export COMPILER="intel"
export MPILIB="mpt"

# Create a build directory
cwd=`pwd`
rm -rf bld
mkdir bld
cd bld
# Run the cime configure tool to figure out what modules need to be loaded
../../../cime/tools/configure --macros-format CMake --machine $MACH --compiler $COMPILER --mpilib $MPILIB
source ./.env_mach_specific.sh
# Build the cmake files
CC=mpicc FC=mpif90 cmake -DCMAKE_BUILD_TYPE=debug ../src
# Build the actual executable
make VERBOSE=1

