# This file is for user convenience only and is not used by the model
# Changes to this file will be ignored and overwritten
# Changes to the environment should be made in env_mach_specific.xml
# Run ./case.setup --reset to regenerate this file
. /glade/u/apps/ch/opt/lmod/7.5.3/lmod/lmod/init/sh
module purge 
module load ncarenv/1.3 python/3.7.9 cmake/3.22.0 intel/19.1.1 esmf_libs mkl
module use /glade/campaign/cesm/cesmdata/cseg/PROGS/modulefiles/esmfpkgs/intel/19.1.1/
module load esmf-8.4.1b02-ncdfio-mpt-O mpt/2.25 netcdf-mpi/4.9.0 pnetcdf/1.12.3 ncarcompilers/0.5.0 pio/2.5.10
export OMP_STACKSIZE=1024M
export TMPDIR=/glade/derecho/scratch/erik
export MPI_TYPE_DEPTH=16
export MPI_USE_ARRAY=None
export COMPILER=intel
export MPILIB=mpt
export DEBUG=FALSE
export OS=LINUX
