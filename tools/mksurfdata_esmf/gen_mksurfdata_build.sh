#! /bin/bash -f

#----------------------------------------------------------------------
# Usage subroutine
usage() {
  echo ""
  echo "***********************************************************************"
  echo "usage:"
  echo "./gen_mksurfdata_build.sh"
  echo ""
  echo "valid arguments: "
  echo "[-h|--help]  "
  echo "     Displays this help message"
  echo "[-v|--verbose]  "
  echo "     Run in verbose mode"
  echo "[-b|--blddir <blddir>]  "
  echo "     Overwrites default, which is /tool_bld in the same directory as ./gen_mksurfdata_build.sh"
  echo "[-m|--machine <machine>]  "
  echo "     Overwrites default MACH"
  echo "***********************************************************************"
}


# Current working directory: the location of ./gen_mksurfdata_build.sh
cwd=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Default settings
verbose="No"
blddir=$cwd/tool_bld  # may overwrite this default with command-line option (below)

# Define what machine to use that's been ported to cime
# May overwrite this default with command-line option --machine
hostname=`hostname --short`
case $hostname in
  ##cheyenne
  cheyenne* | r* )
      export MACH="cheyenne"
      pio_iotype=1
      ;;
  ##casper
  casper* )
      export MACH="casper"
      pio_iotype=1
      ;;
  ##izumi
  izumi*)
      export MACH="izumi"
      pio_iotype=2
      ;;
  ##Hobart
  hobart*)
      export MACH="hobart"
      pio_iotype=2
      ;;
  ## Other machines
  ## Assumption: pnetcdf is off; therefore, pio_iotype = 2
  *)
      export MACH="$hostname"
      pio_iotype=2
      ;;
esac

# Parse command-line options
while [ $# -gt 0 ]; do
   case $1 in
       -h|--help )
           usage
           exit 0
           ;;
       -v|--verbose )
           verbose="Yes"
           ;;
       -b|--blddir )
           blddir=$2
           shift
           ;;
       -m|--machine )
           MACH=$2
           shift
           ;;
       * )
           echo "ERROR:: invalid argument sent in: $2"
           usage
           exit 1
           ;;
   esac
   shift
done

# Create /tool_bld directory
echo "cime Machine is: $MACH..."
if [ -d "$blddir" ]; then
  echo "A /tool_bld directory exists; remove it to do a clean build..."
  echo " or if you want to use the existing build do (assuming using bash):"
  echo "cd $blddir"
  echo ". .env_mach_specific.sh"
  echo "make"
  exit 1
fi
mkdir $blddir
cd $blddir

# Write pio_iotype to file with name pio_iotype.txt
pio_iotype_filepath=../pio_iotype.txt  # one up from /tool_bld
echo 'VALUE OF pio_iotype WRITTEN BY gen_mksurfdata_build.sh AND USED BY mksurfdata (i.e. THE FORTRAN EXECUTABLE):' > $pio_iotype_filepath
echo $pio_iotype >> $pio_iotype_filepath

# Run the cime configure tool to figure out what modules need to be loaded
echo "Run cime configure for machine $MACH..."
# You can specify the non-default compiler and mpi-library by adding --compiler and --mpilib settings
if [ -z "$COMPILER" ] || [ -z "$MPILIB" ]; then
  echo "configure for the default MPI-library and compiler..."
  $cwd/../../cime/CIME/scripts/configure --macros-format CMake --machine $MACH
else
  echo "configure for the specific MPILIB=$MPILIB and COMPILER=$COMPILER..."
  $cwd/../../cime/CIME/scripts/configure --macros-format CMake --machine $MACH --compiler $COMPILER --mpilib $MPILIB
fi

if [ $? != 0 ]; then
  echo "Error doing configure for machine name: $MACH"
  exit 1
fi
. ./.env_mach_specific.sh
echo "COMPILER = $COMPILER, MPILIB = $MPILIB, DEBUG = $DEBUG, OS = $OS"
if [ -z "$PIO" ]; then
  echo "The PIO directory for the PIO build is required and was not set in the configure"
  echo "Make sure a PIO build is provided for $MACH_$COMPILER with $MPILIB in config_machines"
  exit 1
fi

# Build the cmake files
echo "Do the cmake build..."
CC=mpicc FC=mpif90 cmake -DCMAKE_BUILD_TYPE=debug $cwd/src
if [ $? != 0 ]; then
  echo "Error doing cmake for $MACH $MPILIB $COMPILER"
  exit 1
fi

# Build the executable
echo "Build the mksurfdata_esmf build..."
make VERBOSE=1
if [ $? != 0 ]; then
  echo "Error doing make for $MACH $MPILIB $COMPILER"
  exit 1
fi
echo ""
echo ""
echo ""
echo "Successfully created mksurfdata_esmf executable for: ${MACH}_${COMPILER} for $MPILIB library"
