#!/bin/csh -fv
#=======================================================================
#
#  test_batch.csh
#
#  This is a batch submission script for test-model.pl for NCAR and ORNL
#  platforms and batch queues. At other sites this file can be copied
#  and customized for the batch ques available there.  Also the options
#  sent to "test-model.pl" below can be customized for how you would
#  like to run it.
#
# Note: Make sure the number of nodes and shared-memory CPU's used given
# 	in the batch submission agrees with the settings below.
# Note: To use any of the directives below, replace the leading '##'
#       with a single '#' as appropriate.
# Note: Make sure to set the LAB environment variable as appropriate.
#
#-----------------------------------------------------------------------
# NCAR Linux Cluster: anchorage
# Usage: qsub test_batch.csh
#-----------------------------------------------------------------------
##PBS -N test-model
##PBS -j oe
##PBS -q default
##PBS -l nodes=2
##PBS -l walltime=06:00:00
##PBS -m ae
##PBS -V
#
#-----------------------------------------------------------------------
# NCAR IBM SP: bluesky
# Usage: env SCRIPT_DIR=`pwd` llsubmit test_batch.csh
#-----------------------------------------------------------------------
# @ class            = share
# @ node             = 1
# @ tasks_per_node   = 2
# @ network.MPI      = csss,shared,us
# @ output           = test-model.aix.log
# @ error            = test-model.aix.err
# @ node_usage       = shared
# @ job_type         = parallel
# @ account_no       = 93300006
# @ environment      = COPY_ALL
# @ queue
#
#-----------------------------------------------------------------------
# NCAR SGI: tempest
# Usage: qsub test_batch.csh
#-----------------------------------------------------------------------
##QSUB -q ded_16
##QSUB -l mpp_p=4
##QSUB -mb -me -eo
#
#-----------------------------------------------------------------------
# ORNL Cray X1E: phoenix
# Usage: qsub -A acctname test_batch.csh
#-----------------------------------------------------------------------
##PBS -N test-model
##PBS -j oe
##PBS -q batch
##PBS -l walltime=12:00:00,mppe=16,mem=128gb
##PBS -m ae
##PBS -V
#
#-----------------------------------------------------------------------
# ORNL IBM SP: cheetah
# Usage: env SCRIPT_DIR=`pwd` llsubmit test_batch.csh
#-----------------------------------------------------------------------
## @ class       = batch
## @ node        = 1
## @ network.MPI = csss,shared,us
## @ wall_clock_limit = 03:00:00
## @ output      = test-model.aix.log
## @ error       = test-model.aix.err
## @ node_usage  = shared
## @ job_type    = parallel
## @ tasks_per_node = 2
## @ environment = COPY_ALL
## @ queue
#
#-----------------------------------------------------------------------
# ORNL Linux Cluster: penguin
# Usage: qsub test_batch.csh
#-----------------------------------------------------------------------
##PBS -N test-model
##PBS -l nodes=2
##PBS -l walltime=03:00:00
##PBS -j oe
##PBS -m ae
#
#=======================================================================
set OS = `uname -s`
switch ( $OS )
  case UNICOS/mp:
  case Linux:
     setenv SCRIPT_DIR $PBS_O_WORKDIR
     setenv ROOTDIR $SCRIPT_DIR
     breaksw
  case IRIX64:
     setenv SCRIPT_DIR $QSUB_WORKDIR
     setenv ROOTDIR $SCRIPT_DIR
     breaksw
  default:
     if ( ! $?SCRIPT_DIR )then
       echo "ERROR:: The SCRIPT_DIR env variable is not set\!"
       echo "   Set SCRIPT_DIR to the location of test-model.pl"
       echo "   On IBM SP use either:"
       echo "      env SCRIPT_DIR=`pwd` llsubmit $0"
       echo "      setenv SCRIPT_DIR `pwd`; llsubmit $0"
       exit
     endif
     breaksw
endsw
echo "Changing directory to $SCRIPT_DIR"
cd $SCRIPT_DIR
switch ( $OS )
  case AIX:
     setenv SPMD_NODES 4
     setenv SHMEM_CPUS 4
     echo "Set SPMD_NODES to $SPMD_NODES"
     echo "Set SHMEM_CPUS to $SHMEM_CPUS"
     breaksw
  case IRIX64:
     setenv SPMD_NODES 2
     echo "Set SPMD_NODES to $SPMD_NODES"
     setenv SHMEM_CPUS 2
     echo "Set SHMEM_CPUS to $SHMEM_CPUS"
     breaksw
  case UNICOS/mp:
     # Use the Programming Environment which will be used for production runs
     source /opt/modules/modules/init/csh
     module purge
     module load open
     module load PrgEnv.5407
     module unload mpt
     module load mpt.2.4.0.6
     module load pbs
     module list
     #
     setenv OMP_THREAD_STACK_SIZE 60000000
     #
     setenv SPMD_NODES 4
     echo "Set SPMD_NODES to $SPMD_NODES"
     setenv SHMEM_CPUS 4
     echo "Set SHMEM_CPUS to $SHMEM_CPUS"
     breaksw
  case Linux:
     setenv SPMD_NODES 2
     echo "Set SPMD_NODES to $SPMD_NODES"
     cp -p $PBS_NODEFILE machines.$OS
     setenv SPMD_RUNCMND "mpirun -machinefile $SCRIPT_DIR/machines.$OS -np \$SPMD_NODES"
     breaksw
  default:
    echo "Use default values for number of nodes and shared memory CPUs"
    breaksw
endsw

# CHANGE: Set to the appropriate LAB value
setenv LAB "ncar"
#
echo "Starting test-model.pl"

 ./test-model.pl -res T31cnall -c /fis/cgd/ccr/tcraig/fm/clm3_expa_53
 ./test-model.pl -res T31      -c /fis/cgd/ccr/tcraig/fm/clm3_expa_53
 ./test-model.pl -res T31cn    -c /fis/cgd/ccr/tcraig/fm/clm3_expa_53
 ./test-model.pl -res T31casa  -c /fis/cgd/ccr/tcraig/fm/clm3_expa_53
#./test-model.pl -res T31dgvm  


