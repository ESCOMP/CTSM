#! /bin/csh -f

##-------------------------------------------------------------------------------
# The following script creates a CLM offline model executable
# at T42 model resolution with RTM activated and DGVM and VOCs
# deactivated, and runs the model for one day
##-------------------------------------------------------------------------------

##-------------------------------------------------------------------------------
## USER MODIFICATION:
## Define batch queue settings, necessary script environment variables,
## and  build namelist. Some of the environment variables will be used by 
## the Makefile
##-------------------------------------------------------------------------------

##-------------------------------------------------------------------------------
## 1. Set batch system options if necessary
##    (these are provided below only for NCAR IBM SP, NCAR SGI,
##     NCAR-CGD Linux Cluster, ORNL Cray X1, and NEC SX6)
##-------------------------------------------------------------------------------

##-------------------------------------------------------------------------------
## ------------ NCAR IBM SP: blackforest ------------
##-------------------------------------------------------------------------------
## Setting LoadLeveler options for batch queue submission.
## @ class is the queue in which to run the job.
##     To see a list of available queues, type "llclass" interactively.
## @ node is the number of nodes. @tasks_per_node should be set to 1 because
##     of the hybrid OpenMP/MPI configuration of CLM.  
## @ output and @ error are the names of files written to the directory from
##     which the script is submitted containing STDOUT and STDERR respectively
## @ network.MPI: Has to do with network connection between nodes.  Best to leave alone.
## @ job_type = parallel declares that multiple nodes will be used.
## @ node_usage = not_shared acquires dedicated access to nodes for the job.
## @ queue tells load leveler to submit the job

#@ class          = com_reg
#@ node           = 2
#@ tasks_per_node = 1
#@ output         = clm3.$(jobid)
#@ error          = clm3.$(jobid)
#@ job_type       = parallel
#@ network.MPI    = csss,not_shared,us
#@ node_usage     = not_shared
## @ account_no   = 00000000
## @ wall_clock_limit = 3800
#@ queue
##-------------------------------------------------------------------------------

##-------------------------------------------------------------------------------
## ------------ NCAR SGI O2K: chinoook ------------
##-------------------------------------------------------------------------------
## Setting NQS options for batch queue submission.
## -q  is the queue in which to run the job
## -l  on NCAR's SGI this number should be 4 greater than the number of CPUs 
##     specified on the -q line
## -lT is the time limit in seconds
## -eo says to include both stdout and stderr in the output file
## NOTE: if using a ded_32, make sure that NTASKS or NTHREADS is set to 32 below
## or set the number of threads or tasks to the appropriate queue setting

#QSUB -q ded_32          # batch queue
#QSUB -l mpp_p=36        # max number of CPUs accessible
#QSUB -mb -me -eo        # combine stderr & stout
#QSUB -lT 20000          # CPU time limits
#QSUB -eo
#QSUB
##-------------------------------------------------------------------------------

##-------------------------------------------------------------------------------
## ------------ NCAR CGD Linux Cluster: anchorage -----------
##-------------------------------------------------------------------------------
### Job name 
### #PBS -N clm3
###
### Declare job non-rerunable 
### #PBS -r n 
###
### Declare STDERR and STDOUT files 
### #PBS -e clm3_08_nsT42cam-pbs.err 
### #PBS -o clm3_08_nsT42cam-pbs.log 
###
### Mail to user 
### #PBS -m ae 
###
### Queue name (small, medium, long, verylong) 
### #PBS -q small
### 
### Request 4 nodes: 
### #PBS -l nodes=4 
###
### Request 4 nodes with 2 CPUs on each node (total of 8 CPUs): 
### #PBS -l nodes=1:ppn=2 
###
### Request 2 nodes with 4 CPUs on each node that have "hyper" turned on: 
### #PBS -l nodes=2:hyper:ppn=4 

#PBS -N clm3
#PBS -r n 
#PBS -e clm3.err 
#PBS -o clm3.log 
#PBS -m ae 
#PBS -q small
#PBS -l nodes=4 
##-------------------------------------------------------------------------------
 
##-------------------------------------------------------------------------------
## ------------ ORNL Cray X1: phoenix ------------
##-------------------------------------------------------------------------------
##PBS -N clm3
##PBS -j oe
##PBS -q batch
##PBS -l walltime=2:00:00,mppe=16,mem=8gb
##PBS -m ae
##PBS -V
##-------------------------------------------------------------------------------

##-------------------------------------------------------------------------------
## ------------ NEC SX6 --------------------------
##-------------------------------------------------------------------------------
##@$-q hnsxlg
##@$-s /bin/csh
##@$-eo
##@$-ro
##-------------------------------------------------------------------------------

##-------------------------------------------------------------------------------
# 2. Set model resolution 
##-------------------------------------------------------------------------------

setenv LSMLON       128
setenv LSMLAT       64
setenv MAXPATCH_PFT numpft+1

##-------------------------------------------------------------------------------
# 3. Set root directory for CLM source distribution
##-------------------------------------------------------------------------------

setenv ROOTDIR /fs/cgd/data0/mvertens/src/clm/clm_exp/clm3_expa_15

##-------------------------------------------------------------------------------
# 4. Set directory to build the executable in
##-------------------------------------------------------------------------------

setenv MODEL_EXEDIR /scratch/cluster/$LOGNAME/pftdyn

##-------------------------------------------------------------------------------
# 5. Set root directory for input data
##-------------------------------------------------------------------------------

setenv CSMDATA /fs/cgd/csm/inputdata/lnd/clm2

##-------------------------------------------------------------------------------
# 6. Determine if MPI and/or OpenMP 
##-------------------------------------------------------------------------------

# SPMD TRUE implies that MPI will be enabled
# --If SPMD is TRUE , NTASKS is the number of MPI tasks
# --If SPMD is FALSE, NTASKS is ignored
# SMP  TRUE imples  that OpenMP will be enabled
# -- If SMP is TRUE , NTHREADS is the number of OpenMP threads
# -- If SMP is FALSE, NTHREADS is ignored

# Note: it is highly recommended to set SMP to FALSE for Linux systems
# since it is not supported under pgf90 and it does not work robustly
# under pgf90

# NOTE:
# In SPMD mode    : NTHREADS must be set equal to (CPUs-per-node / tasks_per_node)
# In non-SPMD mode: NTHREADS must be set equal to (CPUs-per-node)

setenv SPMD     FALSE
setenv NTASKS   1 

setenv SMP      FALSE
setenv NTHREADS 1 


# force OpenMP threading to be disabled, regardless of the setting above
setenv OS `uname -s`
setenv HOST `uname -n`
if ($OS == 'UNICOS/mp') then
   # Need gmake alias for Cray X1 (but not when cross-compiling on barracuda1 or robin)
   alias gmake /opt/open/open/bin/make
endif
if ($HOST == barracuda1.ccs.ornl.gov) setenv OS "UNICOS/mp"
if ($HOST == robin.ccs.ornl.gov) setenv OS "UNICOS/mp"
if ($OS == Linux) setenv SMP FALSE

##-------------------------------------------------------------------------------
# 7. Set required library and include paths
##-------------------------------------------------------------------------------

## ---SGI---
if ($OS == 'IRIX64') then

  # SGI-specific stuff. MIPSpro is required for f90 (NCAR specific)
  # NOTE: $MPT_SGI is obtained from invoking module mpt
  source /opt/modules/modules/init/csh        
  module purge
  module load MIPSpro mpt nqe modules
  setenv INC_NETCDF  /usr/local/include
  setenv LIB_NETCDF  /usr/local/lib64/r4i4
  if ($SPMD == TRUE) then
     setenv INC_MPI $MPT_SGI/usr/include
     setenv LIB_MPI $MPT_SGI/usr/lib64
  endif

## --- IBM SP --
else if ($OS == 'AIX') then

  # LIB_MPI and INC_MPI do not need to be set
  setenv INC_NETCDF /usr/local/include
  setenv LIB_NETCDF /usr/local/lib64/r4i4

## --- Linux ---
else if ($OS == Linux) then

  # set linux_fort_compiler, the linux fortran compiler [pgf90, lf95] 
  # set the linux machine name (anchorage)

  set linux_fort_compiler = lf95
  set linux_machine = anchorage

  if ($linux_machine == 'anchorage') then 

     if ($linux_fort_compiler == 'pgf90') then 
        setenv INC_NETCDF /usr/local/netcdf-3.5.0/include
        setenv LIB_NETCDF /usr/local/netcdf-3.5.0/lib
        if ($SPMD == TRUE) then
           setenv INC_MPI /usr/local/mpich/include
           setenv LIB_MPI /usr/local/mpich/lib
           set linux_mpirun_cmnd = /usr/local/mpich/bin/mpirun
        endif
     endif
     if ($linux_fort_compiler == 'lf95') then
        setenv INC_NETCDF /usr/local/netcdf-3.5.0/include
        setenv LIB_NETCDF /usr/local/netcdf-3.5.0/lib
        if ($SPMD == TRUE) then
           setenv INC_MPI /usr/local/mpich-lf95/include
           setenv LIB_MPI /usr/local/mpich-lf95/lib
           set linux_mpirun_cmnd = /usr/local/mpich-lf95/bin/mpirun
        endif
     endif
    
  endif

## --- Cray X1 ---
else if ($OS == 'UNICOS/mp') then

  source /opt/modules/modules/init/csh
  module purge
  module load open
  #
  # PE5105
  #module load PrgEnv.5105
  #module unload mpt.2.3.0.2
  #module load mpt.2.3.0.3
  # PE5302
  module load PrgEnv.5302
  module unload mpt
  module load mpt.2.4.0.2
  #
  module load pbs
  module list

  # LIB_MPI and INC_MPI do not need to be set
  # PE5105
  #setenv INC_NETCDF /apps/netcdf/3.5.1/x1_pe51_r4/include
  #setenv LIB_NETCDF /apps/netcdf/3.5.1/x1_pe51_r4/lib
  # PE5302
  setenv INC_NETCDF /apps/netcdf/3.5.1/x1_pe53_r4/include
  setenv LIB_NETCDF /apps/netcdf/3.5.1/x1_pe53_r4/lib

## --- NEC SX6 ---
else if ($OS == 'SUPER-UX') then

  # LIB_MPI and INC_MPI do not need to be set
  setenv INC_NETCDF ${HOME}/SX/include
  setenv LIB_NETCDF ${HOME}/SX/lib

endif


##-------------------------------------------------------------------------------
# 8. Create an input parameter namelist file      
##-------------------------------------------------------------------------------

mkdir -p $MODEL_EXEDIR; cd $MODEL_EXEDIR
cat >! lnd.stdin << EOF
 &clmexp
 caseid         = 'pftdyn'
 ctitle         = 'pftdyn'
 finidat        = ' '
 fsurdat        = '/fs/cgd/data0/mvertens/src/clm/clm_exp/clm3_expa_15/bld/offline/surface-data.128x064.nc' 
 fpftdyn        = '/fs/cgd/data0/mvertens/src/clm/clm_exp/clm3_expa_15/bld/offline/surface-data.dynpft.128x064.nc'
 fpftcon        = '$CSMDATA/pftdata/pft-physiology-cn16.c040719'
 frivinp_rtm    = '$CSMDATA/rtmdata/rdirc.05'
 offline_atmdir = '$CSMDATA/NCEPDATA'
 nsrest         =  0
 nelapse        =  48
 dtime          =  1800
 start_ymd      =  19980101
 start_tod      =  0
 irad           = -1
 wrtdia         = .true.
 mss_irt        =  0
 hist_dov2xy    = .true.
 hist_nhtfrq    =  -24
 hist_mfilt     =  1
 hist_crtinic   = 'MONTHLY'
 /
EOF

##-------------------------------------------------------------------------------
# 9. Determine which model cpp options will be enabled
##-------------------------------------------------------------------------------

##------------------------------------------------------
## build preproc.h in ./obj directory
## define OFFLINE if offline mode
## define LSMLON  to number of longitudes on land grid
## define LSMLAT  to number of latitudes  on land grid
## define RTM     if using RTM river routing
## define VOC     if want to use VOCMod.F90
## define DGVM    if want to use DGVM mode
##------------------------------------------------------
if !( -d $MODEL_EXEDIR/obj ) mkdir -p $MODEL_EXEDIR/obj; cd $MODEL_EXEDIR/obj

\cat >! .tmp << EOF; cmp -s .tmp preproc.h || mv -f .tmp preproc.h
#define OFFLINE
#define LSMLON       $LSMLON
#define LSMLAT       $LSMLAT
#define MAXPATCH_PFT $MAXPATCH_PFT
#undef  DUST
#undef  RTM
#undef  VOC
#undef  DGVM
#define CN
#define SUPLN
#define SUNSHA
#define STOMATA2
EOF

##=======================================================================
## In general do not need to edit below this point
##=======================================================================

##------------------------------------------------------
## build misc.h in ./obj directory
## define SPMD if running in message passing mode
##------------------------------------------------------
if ($SPMD == TRUE) then
   set spmd = "#define SPMD"   
else
   set spmd = "#undef  SPMD"
endif
\cat >! .tmp << EOF; cmp -s .tmp misc.h || mv -f .tmp misc.h
$spmd
EOF

##------------------------------------------------------
## build Filepath 
##------------------------------------------------------
\cat >! .tmp << EOF; cmp -s .tmp Filepath || mv -f .tmp Filepath
$ROOTDIR/src/main
$ROOTDIR/src/biogeophys
$ROOTDIR/src/biogeochem
$ROOTDIR/src/riverroute
$ROOTDIR/src/csm_share/shr
$ROOTDIR/src/utils/timing
$ROOTDIR/src/utils/esmf_wrf_timemgr
$ROOTDIR/src/mksrfdata
EOF

touch $MODEL_EXEDIR/compile_log.clm
echo '------------------------------------------------------------------' >>& $MODEL_EXEDIR/compile_log.clm
date                                                                      >>& $MODEL_EXEDIR/compile_log.clm
echo '------------------------------------------------------------------' >>& $MODEL_EXEDIR/compile_log.clm

##------------------------------------------------------
## build executable
##------------------------------------------------------

## Set debugging flag to TRUE or FALSE
setenv DEBUG TRUE

## Following settings are linux-specific
if ($OS == Linux) then

  ## Following is set to transfer data from all open units in big  endian  format
  if ($linux_fort_compiler == lf95) then
      setenv USER_FC lf95
      setenv FORT90L -Wl,-T
  endif
  if ($linux_fort_compiler == 'pgf90') then 
      setenv USER_FC pgf90
  endif
endif

echo ' '
echo ' compiling clm model'
echo ' The following Makefile environment is used'
               echo ' ROOTDIR      is ' $ROOTDIR
               echo ' MODEL_EXEDIR is ' $MODEL_EXEDIR 
               echo ' SMP          is ' $SMP 
               echo ' SPMD         is ' $SPMD
               echo ' DEBUG        is ' $DEBUG
               echo ' INC_NETCDF   is ' $INC_NETCDF
               echo ' LIB_NETCDF   is ' $LIB_NETCDF
if ($?USER_FC) then
               echo ' USER_FC      is ' $USER_FC
endif 
if ($?INC_MPI) then
               echo ' INC_MPI      is ' $INC_MPI
endif
if ($?LIB_MPI) then
               echo ' LIB_MPI      is ' $LIB_MPI
endif
if ($?LD_LIBRARY_PATH) then
               echo ' LD_LBRARY_PATH      is ' $LD_LIBRARY_PATH
endif
if ($?linux_machine) then
               echo ' Linux machine name  is ' $linux_machine
endif
if ($?linux_mpirun_cmnd) then
               echo ' Linux mpirun cmnd   is ' $linux_mpirun_cmnd
endif

# On NEC SX6 compiling and linking shold be done in a separate script 
# since cross compilers are used to do the builds. The user will have to do
# this independently since this is not officially supported in this script.

if ($OS == 'UNICOS/mp') then
   alias gmake /opt/open/open/bin/make
endif
gmake -j8 -f $ROOTDIR/bld/offline/Makefile >>& $MODEL_EXEDIR/compile_log.clm 
if ($status != 0) then
  echo GMAKE ERROR failed: see $MODEL_EXEDIR/compile_log.clm
  exit  5
else
  echo ' Compiled CLM model successfully'
endif

## go into executable directory
#------------------------------------------------------
cd $MODEL_EXEDIR; setenv LID  "`date +%y%m%d-%H%M%S`"    

##------------------------------------------------------
## run executable
##------------------------------------------------------
echo ' executing clm model'
echo ' using input model data from ' $CSMDATA

## Do our best to get sufficient stack memory
limit stacksize unlimited

## SGI
#------------------------------------------------------
if ($OS == 'IRIX64') then
  ## MP_SLAVE_STACKSIZE sets the size of the thread stack
  setenv    TRAP_FPE "UNDERFL=FLUSH_ZERO; OVERFL=ABORT,TRACE; DIVZERO=ABORT,TRACE; INVALID=ABORT,TRACE"
  setenv    OMP_DYNAMIC   FALSE
  setenv   _DSM_PLACEMENT ROUND_ROBIN
  setenv   _DSM_WAIT      SPIN
  setenv    MPC_GANG      OFF
  setenv MP_SLAVE_STACKSIZE 40000000
  setenv MP_SET_NUMTHREADS $NTHREADS  
  if ($SPMD == TRUE) then
     mpirun -np $NTASKS clm < lnd.stdin >&! clm.log.$LID || echo "CLM run failed" && exit 1 
  else
                        clm < lnd.stdin >&! clm.log.$LID || echo "CLM run failed" && exit 1
  endif
endif

## IBM
#------------------------------------------------------
if ($OS == 'AIX') then
  ## POE Environment.  Set these for interactive jobs.  
  ## They're ignored by LoadLeveler for batch jobs
  ## MP_NODES is the number of nodes.  XLSMPOPTS sets the size of the thread stack
  ## Note - the number of threads here is not set by NTHREADS by rather
  ## is set by the number of tasks per node divided by the number 
  ## of processors per node
  setenv MP_EUILIB us
  setenv MP_NODES $NTASKS
  setenv MP_TASKS_PER_NODE 1
  setenv MP_RMPOOL 1
  setenv XLSMPOPTS "stack=86000000"
  setenv OMP_DYNAMIC false
  setenv OMP_NUM_THREADS $NTHREADS 
  if ($SPMD == TRUE) then
      poe clm < lnd.stdin >&! clm.log.$LID || echo "CLM run failed" && exit 1
  else
          clm < lnd.stdin >&! clm.log.$LID || echo "CLM run failed" && exit 1
  endif 
endif

## Linux
#------------------------------------------------------
if ($OS == 'Linux') then

  if ($SPMD == TRUE) then

     if ($?PBS_NODEFILE) then 
       #
       # This job's working directory 
       echo `date` 
       echo Running offline CLM3 on host `hostname` 
       echo Time is `date` 
       echo Directory is `pwd` 
       echo This job runs on the following processors: 
       #
       # Define number of processors for the MPI command 
       set NPROCS=`wc -l < $PBS_NODEFILE` 
       cat "$PBS_NODEFILE" 
       echo This job has allocated $NPROCS nodes 
       #
       # Run the parallel MPI executable "clm" 
       echo "`date` MPIRUN - Start" 
       $linux_mpirun_cmnd -v -machinefile $PBS_NODEFILE -np $NPROCS clm < lnd.stdin >&! clm.log.$LID || echo "CLM run failed" && exit 1
       echo "`date` MPIRUN - END"                       
     else
       echo " running in interactive mode"
       $linux_mpirun_cmnd -v -np $NTASKS clm < lnd.stdin >&! clm.log.$LID || echo "CLM run failed" && exit 1
       echo "`date` MPIRUN - END"                       
     endif   
     echo "`date` MPIRUN - END"                       
  else     
     clm < lnd.stdin >&! clm.log.$LID || echo "CLM run failed" && exit 1
  endif

endif

## Cray X1
#------------------------------------------------------
if ($OS == 'UNICOS/mp') then
   setenv TRACEBK 100
   env OMP_NUM_THREADS=$NTHREADS aprun -A -d $NTHREADS -n $NTASKS clm < lnd.stdin >&! clm.log.$LID || echo "CLM run failed" && exit 1
endif

## NEC SX6
#------------------------------------------------------
if ($OS == 'SUPER-UX') then
  cd ${RUNDIR}
  setenv MPISUSPEND ON
  setenv F_SETBUF6 0
  if ($NTASKS > 1) then
     mpirun -np $NTASKS ${MODEL_EXEDIR}/clm < ${MODEL_EXEDIR}/lnd.stdin >&! clm.log.$LID || echo "CLM run failed" && exit 1
  else
     ${MODEL_EXEDIR}/clm < ${MODEL_EXEDIR}/lnd.stdin >&! clm.log.$LID || echo "CLM run failed" && exit 1
  endif
endif

wait
echo "`date` -- CLM EXECUTION HAS FINISHED" 
