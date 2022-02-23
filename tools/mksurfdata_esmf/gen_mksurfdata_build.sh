#! /bin/bash -f

module load cmake

export MACH="cheyenne"
export COMPILER="intel"
export OS="LINUX"
export ESMFMKFILE="/glade/work/oehmke/ESMF/scalable_mesh_from_file/lib/libg/Linux.intel.64.mpt.default/esmf.mk"

cwd=`pwd`
rm -rf bld
mkdir bld
cd bld
../../../cime/tools/configure --macros-format CMake --machine cheyenne --compiler intel --mpilib mpt
ls -l
source ./.env_mach_specific.sh
CC=mpicc FC=mpif90 cmake -DCMAKE_BUILD_TYPE=debug ../src
make VERBOSE=1

