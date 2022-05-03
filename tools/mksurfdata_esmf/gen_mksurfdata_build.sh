#! /bin/bash -f

hostname=`hostname --short`

# Define what machine to use that's been ported to cime
case $hostname in
  ##cheyenne
  cheyenne* | r* )
      export MACH="cheyenne"
      ;;
  ## Other machines
  *)
      export MACH="$hostname"
      ;;
esac

# Create a build directory
echo "cime Machine is: $MACH..."
if [ -d "bld" ]; then
  echo "An existing bld directory exists already remove it, if you want to do a clean build..."
  exit 1
fi
cwd=`pwd`
rm -rf bld
mkdir bld
cd bld
# Run the cime configure tool to figure out what modules need to be loaded
echo "Run cime configure for machine $MACH..."
# You can specify the non-default compiler and mpi-library by adding --compiler and --mpilib settings
if [ -z "$COMPILER" ] && [ -z "$MPILIB" ]; then
  echo "configure for the default MPI-library and compiler..."
  ../../../cime/tools/configure --macros-format CMake --machine $MACH
else
  echo "configure for the specific MPILIB=$MPILIB and COMPILER=$COMPILER..."
  ../../../cime/tools/configure --macros-format CMake --machine $MACH --compiler $COMPILER --mpilib $MPILIB
fi

if [ $? != 0 ]; then
  echo "Error doing configure for machine name: $MACH"
  exit 1
fi
. .env_mach_specific.sh
echo "COMPILER = $COMPILER, MPILIB = $MPILIB, DEBUG = $DEBUG, OS = $OS"
if [ -z "$PIO" ]; then
  echo "The PIO directory for the PIO build is required and was not set in the configure"
  echo "Make sure a PIO build is provided for $MACH_$COMPILER with $MPILIB in config_machines"
  exit 1
fi
# Build the cmake files
echo "Do the cmake build..."
CC=mpicc FC=mpif90 cmake -DCMAKE_BUILD_TYPE=debug ../src
if [ $? != 0 ]; then
  echo "Error doing cmake for $MACH $MPILIB $COMPILER"
  exit 1
fi
# Build the actual executable
echo "Build the mksurfdata_esmf build..."
make VERBOSE=1
if [ $? != 0 ]; then
  echo "Error doing make for $MACH $MPILIB $COMPILER"
  exit 1
fi
echo "\n\n\nSuccessfully created mksurfdata_esmf executable for: ${MACH}_${COMPILER} for $MPILIB library"
