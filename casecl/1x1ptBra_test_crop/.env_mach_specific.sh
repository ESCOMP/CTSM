module purge 
module load ncarenv/1.2
module load intel/17.0.1
module load esmf_libs
module load mkl
module load esmf-7.0.0-defio-mpi-O
module load mpt/2.15f
module load netcdf-mpi/4.4.1.1
module load pnetcdf/1.8.0
module load ncarcompilers/0.4.1
export OMP_STACKSIZE=256M
export TMPDIR=/glade/scratch/lfjiang
export MPI_TYPE_DEPTH=16