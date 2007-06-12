#! /usr/bin/csh -f

#-----------------------------------------------------------------------
## IBM
##------------
##
## This is an example script to build and run the default CLM configuration
## on an IBM SP bluevista
## Setting LSF options for batch queue submission.
#BSUB -a poe                    # use LSF openmp elim
#BSUB -x                        # exclusive use of node (not_shared)
#BSUB -n 32                     # total tasks and threads (processors) needed
#BSUB -R "span[ptile=16]"       # max number of tasks (MPI) per node
#BSUB -o out.%J                 # output filename
#BSUB -e out.%J                 # error filename
#BSUB -q premium                # queue
#BSUB -W 0:30                   # wall clock limit
#BSUB -P 93300006

## POE Environment.  Set these for interactive jobs.  They're ignored by LSF
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
setenv  CLMDATA  /fs/cgd/csm/inputdata/lnd/clm2
setenv DATMDATA  /fs/cgd/csm/inputdata/atm/datm7

## Default configuration settings:
## By default both spmd and smp are off
set spmd     = on       # settings are [on   | off       ] (default is off)
set smp      = off      # settings are [on   | off       ] (default is off)
set maxpft   = 4        # settings are 4->17               (default is 4)
set bgc      = none     # settings are [none | cn | casa ] (default is none)
set supln    = off      # settings are [on   | off       ] (default is off)
set rtm      = on       # settings are [on   | off       ] (default is off)
set dust     = off      # settings are [on   | off       ] (default is off)   
set voc      = off      # settings are [on   | off       ] (default is off)   

## $wrkdir is a working directory where the model will be built and run.
## $blddir is the directory where model will be compiled.
## $rundir is the directory where the model will be run.
## $cfgdir is the directory containing the CLM configuration scripts.
set case    = I_1.9x2.5_gx1v4_seq
set wrkdir  = /ptmp/$LOGNAME
set blddir  = $wrkdir/$case/bld
set rundir  = $wrkdir/$case/run
set cfgdir  = $clmroot/bld

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## Create the namelist
cd $rundir                      || echo "cd $blddir failed" && exit 1

cat >! drv_in << EOF
&ccsm_inparm
 case_name	= '$case'
 mss_irt	= 0
 restart_pfile	= './rpointer.drv'
 start_type	= "startup"
/
&timemgr_inparm
 atm_cpl_dt	= 3600
 orb_iyear_ad	= 1990
 start_ymd	= 00010101
 stop_n		= 1
 stop_option	= 'nyears'
 restart_option	= 'yearly'
/
EOF

cat >! lnd_in << EOF
 &clm_inparm
 finidat        = ' '
 fsurdat        = "$CLMDATA/surfdata/surfdata_1.9x2.5_061130.nc"
 fatmgrid       = "$CLMDATA/griddata/griddata_1.9x2.5_060404.nc"
 fatmlndfrc     = "$CLMDATA/griddata/fracdata_1.9x2.5_gx1v4_060808.nc"
 fpftcon        = '$CLMDATA/pftdata/pft-physiology.c061129'
 frivinp_rtm    = "$CLMDATA/rtmdata/rdirc.05.061026"
 dtime          =  1800
 irad           = -1
 wrtdia         = .false.
 mss_irt        =  0
 hist_dov2xy    = .true.
 hist_nhtfrq    =  0
 hist_mfilt     =  1
 hist_crtinic   = 'YEARLY'
 /
EOF

cat >! datm_dshr_in << EOF
 &dshr_nml
 restPfile      = 'rpointer.atm '
 dataMode       = 'CLMNCEP'
 domainFile     = 'domain.lnd.1.9x2.5_gx1v4.060810.nc'
 streams        = 'clmncep.T62.stream.060215.txt 1 2003 2004'
 infoDbug       = 2
 /
EOF

cat >! datm_in << EOF
 &datm_nml
  tn460_factorFn = "unused  "
 /
EOF

## If an executable doesn't exist, build one.
set flags = "-mode ccsm_seq -maxpft $maxpft -bgc $bgc -supln $supln -rtm $rtm -voc $voc -dust $dust"
if ($spmd == on ) set flags = "$flags -spmd"
if ($spmd == off) set flags = "$flags -nospmd"
if ($smp  == on ) set flags = "$flags -smp"
if ($smp  == off) set flags = "$flags -nosmp"
if ( ! -x $blddir/clm ) then
    echo "cd $blddir"
    cd $blddir                  || echo "cd $blddir failed" && exit 1
    echo "blddir is $blddir"
    echo "cfgdir is $cfgdir"
    echo "flags to configure are $flags"
    $cfgdir/configure $flags    || echo "configure failed" && exit 1
    echo "building CLM in $blddir ..."
    rm -f Depends
    gmake -j8 >&! MAKE.out      || echo "CLM build failed: see $blddir/MAKE.out" && exit 1
endif

## Prestage the relevant data
cd $rundir
if !(-f clmncep.T62.stream.060215.txt     ) cp $DATMDATA/CLMNCEP/clmncep.T62.stream.060215.txt .
if !(-f domain.lnd.1.9x2.5_gx1v4.060810.nc) cp $DATMDATA/domain.lnd.1.9x2.5_gx1v4.060810.nc    .  
if !(-f domain.T62.050609.nc              ) cp $DATMDATA/domain.T62.050609.nc                  .

## Run CLM utilizing sequential driver

cd $rundir                    || echo "cd $rundir failed" && exit 1
echo "running CLM in $rundir"

setenv LID "`date +%y%m%d-%H%M%S`"

cd $rundir
if ($spmd == on) then
  mpirun.lsf $blddir/clm >&! clm.log.$LID       || echo "CLM run failed" && exit 1
else 
  $blddir/clm  >&! clm.log.$LID                  || echo "CLM run failed" && exit 1
endif

exit 0
