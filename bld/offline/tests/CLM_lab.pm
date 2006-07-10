#
#	CLM_lab.pm			Erik Kluzek
#
#	Default settings needed to run the scripts at different labs.
#
#	$Id$
#
use strict;
use lab_default;
use Env qw(LOGNAME HOME MODEL_CFGDIR SPMD_NODES SPMD_CPUS_ON_NODES);

#
# List of labs setup with default values
#
@CLM_run::LAB_list = ("ncar", "ornl", "dao", "nersc", "llnl", "default");
#
# Location for case directory
#
my %CASE_DIR_ncar;
my %CASE_DIR_ornl;
my %CASE_DIR_dao;
my %CASE_DIR_nersc;
my %CASE_DIR_llnl;
my %CASE_DIR_default;
$CASE_DIR_ncar{'default'} = "/ptmp/$LOGNAME/clmtest";
$CASE_DIR_ncar{'solaris'} = "/net/flagstaf/export/scratch/$LOGNAME";
$CASE_DIR_ncar{'linux'} = "/scratch/cluster/$LOGNAME";
$CASE_DIR_default{'default'} = "\$MODEL_CFGDIR";
$CASE_DIR_default{'dec_osf'} = "/tmp/$LOGNAME";
$CASE_DIR_llnl{'aix'} = "/p/gb1/$LOGNAME";
$CASE_DIR_llnl{'default'} = "/nfs/tmp0/$LOGNAME";
$CASE_DIR_ornl{'aix'} = "/tmp/work/$LOGNAME";
$CASE_DIR_ornl{'linux'} = "/u1/$LOGNAME";
$CASE_DIR_ornl{'default'} = "/scratch/scr101/$LOGNAME";
$CASE_DIR_dao{'default'} = "/scratch/$LOGNAME";
$CASE_DIR_nersc{'default'} = "/scratch/scratchdirs/$LOGNAME";
$CLM_run::CASE_DIR_default = lab_default->new( \@CLM_run::LAB_list, (
                                                'ncar', \%CASE_DIR_ncar, 
                                                'ornl', \%CASE_DIR_ornl,
                                                'dao', \%CASE_DIR_dao,
                                                'llnl', \%CASE_DIR_llnl,
                                                'nersc', \%CASE_DIR_nersc,
                                                'default', \%CASE_DIR_default
                                              ) );
#
# Location to build
#
my %BUILD_DIR_ncar;
my %BUILD_DIR_ornl;
my %BUILD_DIR_llnl;
my %BUILD_DIR_nersc;
my %BUILD_DIR_default;
$BUILD_DIR_default{'aix'} = "\${SCRIPT_DIR}";
$BUILD_DIR_default{'default'} = "\${CASE_DIR}";
$BUILD_DIR_ncar{'default'} = "\${CASE_DIR}";
$BUILD_DIR_ornl{'default'} = "\${CASE_DIR}";
$BUILD_DIR_llnl{'default'} = "\${CASE_DIR}";
$BUILD_DIR_nersc{'aix'} = "\${CASE_DIR}";
$BUILD_DIR_nersc{'default'} = "\${SCRIPT_DIR}";
$CLM_run::BUILD_DIR_default = lab_default->new( \@CLM_run::LAB_list, (
                                                'ncar', \%BUILD_DIR_ncar,
                                                'ornl', \%BUILD_DIR_ornl,
                                                'llnl', \%BUILD_DIR_llnl,
                                                'nersc', \%BUILD_DIR_nersc,
                                                'default', \%BUILD_DIR_default
                                              ) );
#
# Location to put log-files
#
my %LOG_DIR_nersc;
my %LOG_DIR_default;
$LOG_DIR_default{'default'} = "\${MODEL_EXEDIR}";
$LOG_DIR_nersc{'default'} = "\${MODEL_EXEDIR}";
$CLM_run::LOG_DIR_default = lab_default->new( \@CLM_run::LAB_list, (
                                                'nersc', \%LOG_DIR_nersc,
                                                'default', \%LOG_DIR_default
                                              ) );
#
# Location to build model
#
my %MODEL_BLDDIR_default;
my %MODEL_BLDDIR_nersc;
my %MODEL_BLDDIR_ncar;
$MODEL_BLDDIR_default{'aix'} = "\${BUILD_DIR}/\${CASE}obj";
$MODEL_BLDDIR_default{'default'} = "\${BUILD_DIR}/\${CASE}/obj";
$MODEL_BLDDIR_ncar{'default'} = "\${BUILD_DIR}/\${CASE}/obj";
$MODEL_BLDDIR_nersc{'default'} = "\${BUILD_DIR}/\${CASE}/obj";
$CLM_run::MODEL_BLDDIR_default = lab_default->new( \@CLM_run::LAB_list, (
                                                'ncar', \%MODEL_BLDDIR_ncar,
                                                'nersc', \%MODEL_BLDDIR_nersc,
                                                'default', \%MODEL_BLDDIR_default
                                              ) );
#
# Whether to run in SPMD mode or not
#
my %SPMD_default;
my %SPMD_ncar;
my %SPMD_ornl;
$SPMD_default{'default'} = "TRUE";
$SPMD_default{'linux'} = "FALSE";
$SPMD_default{'darwin'} = "FALSE";
$SPMD_default{'irix'} = "FALSE";
$SPMD_ncar{'default'} = "TRUE";
$SPMD_ncar{'dec_osf'} = "FALSE";
$SPMD_ncar{'linux'} = "FALSE";
$SPMD_ncar{'irix'} = "FALSE";
$SPMD_ornl{'default'} = "TRUE";
$CLM_run::SPMD_default = lab_default->new( \@CLM_run::LAB_list, (
                                                'ornl', \%SPMD_ornl,
                                                'ncar', \%SPMD_ncar,
                                                'default', \%SPMD_default
                                              ) );
#
# Name of the GNU-make command and the number of parallel processes to use
#
my %GNUMAKE_dao;
my %GNUMAKE_ncar;
my %GNUMAKE_ornl;
my %GNUMAKE_nersc;
my %GNUMAKE_default;
$GNUMAKE_ncar{'default'} = "gmake -j 2";
$GNUMAKE_ncar{'solaris'} = "gmake ";
$GNUMAKE_ncar{'irix'} = "gmake -j 8";
$GNUMAKE_ncar{'dec_osf'} = "gmake -j 4";
$GNUMAKE_ncar{'aix'} = "gmake -j 4";
$GNUMAKE_ncar{'darwin'} = "make";
$GNUMAKE_ornl{'linux'} = "gmake -j 2";
$GNUMAKE_ornl{'default'} = "gmake -j 8";
$GNUMAKE_nersc{'default'} = "gmake -j 16";
$GNUMAKE_dao{'default'} = "gmake";
$GNUMAKE_dao{'aix'} = "gmake -j 8";
$GNUMAKE_dao{'irix'} = "gmake -j 8";
$GNUMAKE_dao{'linux'} = "gmake -j 2";
$GNUMAKE_default{'default'} = "gmake";
$GNUMAKE_default{'darwin'} = "make";
$CLM_run::GNUMAKE_default = lab_default->new( \@CLM_run::LAB_list, (
                                                 'dao', \%GNUMAKE_dao, 
                                                 'ncar', \%GNUMAKE_ncar, 
                                                 'ornl', \%GNUMAKE_ornl, 
                                                 'nersc', \%GNUMAKE_nersc, 
                                                 'default', \%GNUMAKE_default
                                            ) );
#
# Number of shared memory CPU's to run with
# (If set to 0, will unset these to use default values)
#
my %SHMEM_CPUS_ncar;
my %SHMEM_CPUS_dao;
my %SHMEM_CPUS_ornl;
my %SHMEM_CPUS_llnl;
my %SHMEM_CPUS_default;
$SHMEM_CPUS_ncar{'default'} = 2;
$SHMEM_CPUS_ncar{'irix'}    = 4;
$SHMEM_CPUS_ncar{'linux'}   = 2;
$SHMEM_CPUS_ncar{'aix'}     = 0;
$SHMEM_CPUS_ncar{'dec_osf'} = 4;
$SHMEM_CPUS_dao{'default'}  = 4;
$SHMEM_CPUS_dao{'irix'}     = 4;
$SHMEM_CPUS_dao{'linux'}    = 2;
$SHMEM_CPUS_dao{'aix'}      = 0;
$SHMEM_CPUS_llnl{'default'} = 2;
$SHMEM_CPUS_llnl{'aix'}     = 4;
$SHMEM_CPUS_ornl{'default'} = 2;
$SHMEM_CPUS_ornl{'unicosmp'} = 4;
$SHMEM_CPUS_default{'default'} = 2;
$SHMEM_CPUS_default{'aix'}     = 0;
$SHMEM_CPUS_default{'darwin'}  = 1;
$CLM_run::SHMEM_CPUS_default = lab_default->new( \@CLM_run::LAB_list, (
                                                   'dao', \%SHMEM_CPUS_dao, 
                                                   'ncar', \%SHMEM_CPUS_ncar, 
                                                   'ornl', \%SHMEM_CPUS_ornl, 
                                                   'llnl', \%SHMEM_CPUS_llnl, 
                                                   'default', \%SHMEM_CPUS_default
                                               ) );
#
# Number of nodes when running with MPI
#
my %SPMD_NODES_dao;
my %SPMD_NODES_ncar;
my %SPMD_NODES_ornl;
my %SPMD_NODES_llnl;
my %SPMD_NODES_default;
$SPMD_NODES_dao{'default'}  = 2;
$SPMD_NODES_dao{'aix'}      = 4;
$SPMD_NODES_ncar{'default'} = 2;
$SPMD_NODES_ncar{'aix'}     = 4;
$SPMD_NODES_ncar{'dec_osf'} = 4;
$SPMD_NODES_ornl{'default'} = 2;
$SPMD_NODES_ornl{'aix'} = 4;
$SPMD_NODES_ornl{'unicosmp'} = 4;
$SPMD_NODES_llnl{'default'} = 2;
$SPMD_NODES_llnl{'aix'}     = 4;
$SPMD_NODES_default{'default'} = 2;
$CLM_run::SPMD_NODES_default = lab_default->new( \@CLM_run::LAB_list, (
                                                  'dao', \%SPMD_NODES_dao,
                                                  'ncar', \%SPMD_NODES_ncar,
                                                  'ornl', \%SPMD_NODES_ornl,
                                                  'llnl', \%SPMD_NODES_llnl,
                                                  'default', \%SPMD_NODES_default
                                               ) );
#
# SPMD Run command
#
my %SPMD_RUNCMND_default;
my %SPMD_RUNCMND_llnl;
$SPMD_RUNCMND_default{'default'} = "mpirun -np \$SPMD_NODES";
$SPMD_RUNCMND_default{'solaris'} = "mpirun -machinefile machine -np \$SPMD_NODES";
$SPMD_RUNCMND_default{'aix'}     = "poe";
$SPMD_RUNCMND_default{'dec_osf'} = "prun -N \$SPMD_NODES -c \$SPMD_CPUS_ON_NODE sh -c";
$SPMD_RUNCMND_default{'unicosmp'} = "aprun -A -n \$SPMD_NODES -d \$SHMEM_CPUS";
$SPMD_RUNCMND_llnl{'default'}    = "mpirun -np \$SPMD_NODES ";
$SPMD_RUNCMND_llnl{'aix'}        = "poe";
# $SPMD_RUNCMND_llnl{'dec_osf'}    = "dmpirun -np \$SPMD_NODES ";
$SPMD_RUNCMND_llnl{'dec_osf'} = "prun -p pdebug -N \$SPMD_NODES -c \$SPMD_CPUS_ON_NODE sh -c";
$CLM_run::SPMD_RUNCMND_default = lab_default->new( \@CLM_run::LAB_list, (
                                                  'llnl', \%SPMD_RUNCMND_llnl, 
                                                  'default', \%SPMD_RUNCMND_default
                                               ) );
#
# NetCDF library location
#
my %LIB_NETCDF_default;
my %LIB_NETCDF_ncar;
my %LIB_NETCDF_ornl;
my %LIB_NETCDF_nersc;
my %LIB_NETCDF_dao;
$LIB_NETCDF_ncar{'solaris'} = "/contrib/lib";
$LIB_NETCDF_ncar{'linux'}   = "/usr/local/netcdf/lib";
$LIB_NETCDF_ncar{'irix'}    = "/usr/local/lib64/r4i4";
$LIB_NETCDF_ncar{'default'} = "/usr/local/lib64/r4i4";
$LIB_NETCDF_nersc{'default'}= "$ENV{'NETCDF_DIR'}" . "/lib";
$LIB_NETCDF_ornl{'unicosmp'} = "/apps/netcdf/3.5.1/x1_pe53_r4/lib";
$LIB_NETCDF_ornl{'linux'} = "/usr/local/netcdf-3.5.1_pgi/lib";
$LIB_NETCDF_ornl{'default'} = "/apps/netcdf/prod/rs_aix5_64_r4/lib";
$LIB_NETCDF_dao{'default'} = "/usr/local/lib";
$LIB_NETCDF_dao{'irix'}    = "/ford1/local/IRIX64/netcdf/lib";
$LIB_NETCDF_dao{'linux'}   = "/usr/local/netcdf-3.4_Lahey/lib";
$LIB_NETCDF_default{'default'} = "/usr/local/lib";
$CLM_run::LIB_NETCDF_default = lab_default->new( \@CLM_run::LAB_list, (
                                                  'ncar', \%LIB_NETCDF_ncar, 
                                                  'ornl', \%LIB_NETCDF_ornl, 
                                                  'dao', \%LIB_NETCDF_dao, 
                                                  'nersc', \%LIB_NETCDF_nersc, 
                                                  'default', \%LIB_NETCDF_default
                                               ) );
#
# NetCDF includes location
#
my %INC_NETCDF_default;
my %INC_NETCDF_dao;
my %INC_NETCDF_nersc;
my %INC_NETCDF_ornl;
my %INC_NETCDF_ncar;
$INC_NETCDF_ncar{'solaris'} = "/contrib/include";
$INC_NETCDF_ncar{'linux'}   = "/usr/local/netcdf/include";
$INC_NETCDF_ncar{'default'} = "/usr/local/include";
$INC_NETCDF_nersc{'default'}= "$ENV{'NETCDF_DIR'}" . "/include";
$INC_NETCDF_ornl{'unicosmp'} = "/apps/netcdf/3.5.1/x1_pe53_r4/include";
$INC_NETCDF_ornl{'linux'} = "/usr/local/netcdf-3.5.1_pgi/include";
$INC_NETCDF_ornl{'default'} = "/apps/netcdf/prod/rs_aix5_64_r4/include";
$INC_NETCDF_dao{'default'}  = "/usr/local/include";
$INC_NETCDF_dao{'irix'}     = "/ford1/local/IRIX64/netcdf/include";
$INC_NETCDF_dao{'linux'}    = "/usr/local/netcdf-3.4_Lahey/include";
$INC_NETCDF_default{'default'} = "/usr/local/include";
$CLM_run::INC_NETCDF_default = lab_default->new( \@CLM_run::LAB_list, (
                                                  'ncar', \%INC_NETCDF_ncar, 
                                                  'dao', \%INC_NETCDF_dao, 
                                                  'ornl', \%INC_NETCDF_ornl, 
                                                  'nersc', \%INC_NETCDF_nersc, 
                                                  'default', \%INC_NETCDF_default
                                               ) );
#
# MPI includes location
#
my %INC_MPI_default;
my %INC_MPI_dao;
my %INC_MPI_ncar;
my %INC_MPI_ornl;
$INC_MPI_ncar{'solaris'} = "/contrib/include";
$INC_MPI_ncar{'linux'}   = "/usr/local/mpich/include";
$INC_MPI_ncar{'irix'}   = "$ENV{'MPT_SGI'}" . "/usr/include";
$INC_MPI_ncar{'dec_osf'} = "/usr/include";
$INC_MPI_ncar{'default'} = "/usr/local/include";
$INC_MPI_dao{'default'}  = "/usr/local/include";
$INC_MPI_dao{'irix'}     = "$ENV{'MPT_SGI'}" . "/usr/include";
$INC_MPI_dao{'linux'}    = "/usr/local/mpich-1.2.1-ffc/include";
$INC_MPI_ornl{'linux'} = "/usr/local/mpich_pgi/include";
$INC_MPI_ornl{'default'} = "/usr/include";
$INC_MPI_default{'default'} = "/usr/local/include";
$CLM_run::INC_MPI_default = lab_default->new( \@CLM_run::LAB_list, (
                                                  'ncar', \%INC_MPI_ncar,
                                                  'dao', \%INC_MPI_dao,
                                                  'ornl', \%INC_MPI_ornl,
                                                  'default', \%INC_MPI_default
                                               ) );
#
# MPI library location
#
my %LIB_MPI_default;
my %LIB_MPI_ncar;
my %LIB_MPI_ornl;
my %LIB_MPI_dao;
$LIB_MPI_ncar{'solaris'} = "/contrib/lib";
$LIB_MPI_ncar{'linux'}   = "/usr/local/mpich/lib";
$LIB_MPI_ncar{'dec_osf'} = "/usr/lib";
$LIB_MPI_ncar{'irix'}    = "$ENV{'MPT_SGI'}" . "/usr/lib64";
$LIB_MPI_ncar{'default'} = "/usr/local/lib32/r4i4";
$LIB_MPI_ornl{'linux'} = "/usr/local/mpich_pgi/lib";
$LIB_MPI_ornl{'default'} = "/usr/local/lib64/r4i4";
$LIB_MPI_dao{'default'}  = "/usr/local/lib";
$LIB_MPI_dao{'irix'}     = "$ENV{'MPT_SGI'}" . "/usr/lib64";
$LIB_MPI_dao{'linux'}    = "/usr/local/mpich-1.2.1-ffc/lib";
$LIB_MPI_default{'default'} = "/usr/local/lib";
$CLM_run::LIB_MPI_default = lab_default->new( \@CLM_run::LAB_list, (
                                                  'ncar', \%LIB_MPI_ncar,
                                                  'ornl', \%LIB_MPI_ornl,
                                                  'dao', \%LIB_MPI_dao,
                                                  'default', \%LIB_MPI_default
                                               ) );

#
# Command to use to archive log files
# (Don't set if doesn't exist)
#
my %ARCHIVE_ncar;
my %ARCHIVE_default;
$ARCHIVE_ncar{'default'} = "msrcp";
$ARCHIVE_default{'default'} = " ";
$CLM_run::ARCHIVE_default = lab_default->new( \@CLM_run::LAB_list, (
                                                'ncar', \%ARCHIVE_ncar, 
                                                'default', \%ARCHIVE_default
                                            ) );

#
# Model data directory
my %MODEL_DATDIR_default;
my %MODEL_DATDIR_ncar;
my %MODEL_DATDIR_ornl;
$MODEL_DATDIR_ncar{'default'} = "/fs/cgd/csm/inputdata";
$MODEL_DATDIR_ncar{'darwin'} = "$HOME/inputdata";
$MODEL_DATDIR_ornl{'default'} = "/spin/proj/ccsm/inputdata";
$MODEL_DATDIR_default{'default'} = "/fs/cgd/csm/inputdata";
$MODEL_DATDIR_default{'darwin'} = "$HOME/inputdata";
$CLM_run::MODEL_DATDIR_default = lab_default->new( \@CLM_run::LAB_list, (
						'ncar', \%MODEL_DATDIR_ncar,
						'ornl', \%MODEL_DATDIR_ornl,
						'default', \%MODEL_DATDIR_default
					) );


1   # To make use or require happy
