#! /usr/bin/csh -f

#-----------------------------------------------------------------------
## SGI
##------------
##
## This is an example script to build and run the default CLM configuration on
## an SGI.  
##
## Set NQS options 
## -q  is the queue in which to run the job
## -l  on NCAR's SGI this number should be 4 greater than the number of CPUs 
##       specified on the -q line
## -lT is the time limit in seconds
## -eo says to include both stdout and stderr in the output file

#QSUB -q ded_32
#QSUB -l mpp_p=36
#QSUB -lT 20000
#QSUB -eo

## Number of OMP threads.  Only need to set if requesting
## fewer than the maximum number for the queue.
#setenv OMP_NUM_THREADS 32

## Set runtime variables for optimal execution on a DSM
setenv OMP_DYNAMIC FALSE
setenv _DSM_PLACEMENT ROUND_ROBIN
setenv _DSM_WAIT SPIN
setenv MPC_GANG    OFF

## MP_SLAVE_STACKSIZE sets the size of the thread stack
setenv MP_SLAVE_STACKSIZE 40000000
## Do our best to get sufficient stack memory
limit stacksize unlimited

## SGI-specific stuff. MIPSpro is required for f90
source /opt/modules/modules/init/csh
module purge
module load MIPSpro modules nqe mpt

## netCDF stuff
setenv INC_NETCDF /usr/local/include
setenv LIB_NETCDF /usr/local/lib64/r4i4

## mpi stuff
setenv INC_MPI /opt/mpt/mpt/usr/include
setenv LIB_MPI /opt/mpt/mpt/usr/lib64

## ROOT OF CLM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CLM distribution.
## (the root directory contains the subdirectory "models")
set clmroot   = /fis/cgd/...

## ROOT OF CLM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CLM distribution.
setenv CSMDATA /fs/cgd/csm/inputdata/lnd/clm2

## Default configuration settings:
## By default smp is on
set smp      = on       # settings are [on   | off       ] (default is off)
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
set case     = clmrun
set wrkdir   = /ptmp/$LOGNAME
set blddir   = $wrkdir/$case/bld
set rundir   = $wrkdir/$case
set cfgdir   = $clmroot/bld
set usr_src  = $clmroot/bld/empty

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## If an executable doesn't exist, build one.
set flags = "-maxpft $maxpft -bgc $bgc -supln $supln -rtm $rtm -voc $voc -dust $dust"
if ($smp == on) set flags = "$flags -smp"
if ( ! -x $blddir/clm ) then
    echo "cd $blddir"
    cd $blddir                  || echo "cd $blddir failed" && exit 1
    echo "flags to configure are $flags"
    $cfgdir/configure $flags    || echo "configure failed" && exit 1
    echo "building CLM in $blddir ..."
    rm -f Depends
    gmake -j8 >&! MAKE.out      || echo "CLM build failed: see $blddir/MAKE.out" && exit 1
endif

## Create the namelist
cd $rundir                      || echo "cd $blddir failed" && exit 1

cat >! lnd.stdin << EOF
 &clm_inparm
 caseid         = $case
 ctitle         = $case
 finidat        = ' '
 fsurdat        = "$CSMDATA/surfdata/surfdata_048x096_061108.nc"
 fatmgrid       = "$CSMDATA/griddata/griddata_48x96_060829.nc"
 fatmlndfrc     = "$CSMDATA/griddata/fracdata_48x96_gx3v5_060829.nc"
 fpftcon        = '$CSMDATA/pftdata/pft-physiology.c061129'
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
EOF

## Run ClM
cd $rundir                      || echo "cd $rundir failed" && exit 1
echo "running CLM in $rundir"
$blddir/clm                     || echo "CLM run failed" && exit 1

exit 0
