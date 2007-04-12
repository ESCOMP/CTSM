#! /usr/bin/csh -f

#-----------------------------------------------------------------------
## IBM
##------------
##
## This is an example script to build and run the default CLM configuration
## on an IBM SP. 
##
## Setting LoadLeveler options for batch queue submission.
## @ class is the queue in which to run the job.
##     To see a list of available queues, type "llclass" interactively.
## @ node is the number of nodes. @tasks_per_node should be set to 1 because
##     of the hybrid OpenMP/MPI configuration of CLM.
## @ output and @error are the names of file written to the directory from
##     which the script is submitted containing STDOUT and STDERR respectively
## @ job_type = parallel declares that multiple nodes will be used.
## @ network.MPI: Has to do with network connection between nodes.  Best to leave alone.
## @ node_usage = not_shared acquires dedicated access to nodes for the job.
## @ queue tells load leveler to submit the job
## @ environment tells load leveler to export all Environment variables

#@ class          = csl_rg8
#@ node           = 2
#@ tasks_per_node = 2
#@ output         = out.$(jobid)
#@ error          = out.$(jobid)
#@ job_type       = parallel
#@ network.MPI    = csss,not_shared,us
#@ node_usage     = not_shared
#@ environment = COPY_ALL
#@ queue

## Setting LSF options for batch queue submission.
#BSUB -a poe                    # use LSF openmp elim
#BSUB -x                        # exclusive use of node (not_shared)
#BSUB -n 8                      # total tasks and threads (processors) needed
#BSUB -R "span[ptile=8]"        # max number of tasks (MPI) per node
#BSUB -o out.%J                 # output filename
#BSUB -e out.%J                 # error filename
#BSUB -q regular                # queue
#BSUB -W 0:10                   # wall clock limit

## POE Environment.  Set these for interactive jobs.  They're ignored by LoadLeveler
## MP_NODES is the number of nodes.  
setenv MP_NODES 1
setenv MP_TASKS_PER_NODE 8
setenv MP_EUILIB us
setenv MP_RMPOOL 1

unsetenv MP_PROCS

setenv MP_STDINMODE 0

# should be set equal to (CPUs-per-node / tasks_per_node)
setenv OMP_NUM_THREADS 4

## suggestion from Jim Edwards to reintroduce XLSMPOPTS on 11/13/03
setenv XLSMPOPTS "stack=256000000"
setenv AIXTHREAD_SCOPE S
setenv MALLOCMULTIHEAP true
setenv OMP_DYNAMIC false
## Do our best to get sufficient stack memory
limit stacksize unlimited

## netCDF stuff
setenv INC_NETCDF /usr/local/include
setenv LIB_NETCDF /usr/local/lib64/r4i4

## ROOT OF CLM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CLM distribution.
## (the root directory contains the subdirectory "models")
set clmroot   = /fis/cgd/.......

## ROOT OF CLM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CLM distribution.
setenv CSMDATA /fs/cgd/csm/inputdata/lnd/clm2

## Default configuration settings:
## By default both spmd and smp are off
set spmd     = on       # settings are [on   | off       ] (default is off)
set smp      = off      # settings are [on   | off       ] (default is off)
set maxpft   = 4        # settings are 4->17               (default is 4)
set bgc      = none     # settings are [none | cn | casa ] (default is none)
set supln    = off      # settings are [on   | off       ] (default is off)
set rtm      = off      # settings are [on   | off       ] (default is off)
set dust     = off      # settings are [on   | off       ] (default is off)   
set voc      = off      # settings are [on   | off       ] (default is off)   

## $wrkdir is a working directory where the model will be built and run.
## $blddir is the directory where model will be compiled.
## $rundir is the directory where the model will be run.
## $cfgdir is the directory containing the CLM configuration scripts.
set case    = clmrun
set wrkdir  = /ptmp/$LOGNAME
set blddir  = $wrkdir/$case/bld
set rundir  = $wrkdir/$case
set cfgdir  = $clmroot/bld
set usr_src = $clmroot/bld/empty

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## Build (or re-build) executable
set flags = "-maxpft $maxpft -bgc $bgc -supln $supln -rtm $rtm -voc $voc -dust $dust -usr_src $usr_src"
if ($spmd == on ) set flags = "$flags -spmd"
if ($spmd == off) set flags = "$flags -nospmd"
if ($smp  == on ) set flags = "$flags -smp"
if ($smp  == off) set flags = "$flags -nosmp"
echo "cd $blddir"
cd $blddir                  || echo "cd $blddir failed" && exit 1
if ( ! -f $blddir/config_cache.xml ) then
    echo "flags to configure are $flags"
    $cfgdir/configure $flags    || echo "configure failed" && exit 1
    echo "Building CLM in $blddir ..."
    gmake -j8 >&! MAKE.out      || echo "CLM build failed: see $blddir/MAKE.out" && exit 1
else
    echo "Re-building CLM in $blddir ..."
    rm -f Depends
    gmake -j8 >&! REMAKE.out      || echo "CLM build failed: see $blddir/REMAKE.out" && exit 1
endif

## Create the namelist
cd $rundir                      || echo "cd $blddir failed" && exit 1

cat >! lnd.stdin << EOF
 &clm_inparm
 caseid         = '$case'
 ctitle         = '$case'
 finidat        = ' '
 fsurdat        = "$CSMDATA/surfdata/surfdata_048x096_061108.nc"
 fatmgrid       = "$CSMDATA/griddata/griddata_48x96_060829.nc"
 fatmlndfrc     = "$CSMDATA/griddata/fracdata_48x96_gx3v5_060829.nc"
 fpftcon        = '$CSMDATA/pftdata/pft-physiology.c070207'
 fndepdat       = "$CSMDATA/ndepdata/1890/regrid_ndep_clm.nc"
 frivinp_rtm    = "$CSMDATA/rtmdata/rdirc.05.061026"
 offline_atmdir = "$CSMDATA/NCEPDATA"
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
 &prof_inparm
 /
EOF

## Run CLM 

cd $rundir                    || echo "cd $rundir failed" && exit 1
echo "running CLM in $rundir"

setenv LID "`date +%y%m%d-%H%M%S`"

if ($spmd == on) then

  #uncomment the line below to run on bluesky with MPI
  #poe $blddir/clm >&! clm.log.$LID              || echo "CLM run failed" && exit 1

  #uncomment the line below to run on bluevista with MPI
  #mpirun.lsf $blddir/clm >&! clm.log.$LID       || echo "CLM run failed" && exit 1
else 
  $blddir/clm  >&! clm.log.$LID                  || echo "CLM run failed" && exit 1
endif

exit 0
