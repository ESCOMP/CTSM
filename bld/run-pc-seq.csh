#! /bin/tcsh -f
#
#=======================================================================
#
#  run-pc.csh
#
#  Generic batch submission script for PC-linux using PBS.  
#
#-----------------------------------------------------------------------
# Batch options for machine with PBS batch system. (anchorage)
# Usage for Lahey compiler (default): 
#   qsub run-pc.csh
# Usage for pgf90 compiler with pgcc: 
#   env OPT_BLD=pgf90-pgcc qsub run-pc.csh
# Usage for pgf90 compiler with gcc: 
#   env OPT_BLD=pgf90-gcc qsub run-pc.csh
#-----------------------------------------------------------------------
#
# Name of the queue (CHANGE THIS if needed)
#PBS -q long
# Number of nodes & procs/node - must have ecc memory (CHANGE THIS if needed)
#PBS -l nodes=2:ppn=2
# output file base name
#PBS -N run-pc
# Put standard error and standard out in same file
#PBS -j oe
# Export all Environment variables
#PBS -V
# End of options
#=======================================================================

set OS = `uname -s`;
switch ( $OS )
  case Linux:
     if ( ! $?PBS_JOBID ) then
       echo "${0}: ERROR::  This batch script must be submitted via PBS";
       echo "${0}:          on a Linux machine\!";
       exit;
     else
       echo "${0}: Running CLM on Linux using PBS";
     endif

     setenv PGI /usr/local/pgi
     setenv LAHEY /usr/local/lf95
     setenv LD_LIBRARY_PATH ${PGI}/linux86/6.1/lib:${LAHEY}/lib:/cluster/torque/lib:${LD_LIBRARY_PATH}

     if ( ! $?OPT_BLD ) then
       setenv USER_FC lf95
       setenv LAHEY /usr/local/lf9562
       setenv INC_NETCDF /usr/local/netcdf-3.6.1-gcc-3.4.6-lf9562/include
       setenv LIB_NETCDF /usr/local/netcdf-3.6.1-gcc-3.4.6-lf9562/lib
       set mpich = /usr/local/mpich-1.2.7p1-gcc-g++-3.4.6-3-lf9562
       setenv INC_MPI ${mpich}/include
       setenv LIB_MPI ${mpich}/lib
       setenv PATH ${LAHEY}/bin:${mpich}/bin:${PATH}
     else
       if ( $OPT_BLD == "pgf90-gcc" ) then
         setenv USER_CC gcc
       endif
       setenv INC_MPI /usr/local/mpich-pgi-pgcc-pghf/include
       setenv LIB_MPI /usr/local/mpich-pgi-pgcc-pghf/lib
       setenv INC_NETCDF /usr/local/netcdf-pgi-hpf-cc/include
       setenv LIB_NETCDF /usr/local/netcdf-pgi-hpf-cc/lib
       setenv PGI /usr/local/pgi-pgcc-pghf
       setenv PATH ${PGI}/linux86/6.1/bin:${PATH}
     endif
     setenv MOD_NETCDF $INC_NETCDF

     breaksw;
  default:
    echo "${0}: This script meant for running CLM on Linux machines";    exit;
endsw

## set this equal to #nodes X #ppn
set procs = 4

## Do our best to get sufficient stack memory
limit stacksize unlimited

## ROOT OF CLM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CLM distribution.
## (the root directory contains the subdirectory "models")
set clmroot   = /fis/cgd/.....

## ROOT OF CLM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CLM distribution.
setenv CLMDATA   /fs/cgd/csm/inputdata/lnd/clm2
setenv DATMDATA  /fs/cgd/csm/inputdata/atm/datm7

## Default configuration settings:
## By default spmd is off
set spmd     = on       # settings are [on   | off       ] (default is off)
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
set case     = I_1.9x2.5_gx1v4_seq
set wrkdir   = /ptmp/$LOGNAME
set blddir   = $wrkdir/$case/bld
set rundir   = $wrkdir/$case/run
set cfgdir   = $clmroot/bld

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## If an executable doesn't exist, build one.
set flags = "-mode ccsm_seq -maxpft $maxpft -bgc $bgc -supln $supln -rtm $rtm -voc $voc -dust $dust "
if ($spmd == on ) set flags = "$flags -spmd"
if ($spmd == off) set flags = "$flags -nospmd"
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

## Prestage the relevant data
cd $rundir
if !(-f clmncep.T62.stream.060215.txt     ) cp $DATMDATA/CLMNCEP/clmncep.T62.stream.060215.txt .
if !(-f domain.lnd.1.9x2.5_gx1v4.060810.nc) cp $DATMDATA/domain.lnd.1.9x2.5_gx1v4.060810.nc    .  
if !(-f domain.T62.050609.nc              ) cp $DATMDATA/domain.T62.050609.nc                  .

## Run CLM
cd $rundir                      || echo "cd $rundir failed" && exit 1
echo "running CLM in $rundir"

setenv LID "`date +%y%m%d-%H%M%S`"

if ($spmd == on) then
  mpiexec -n $procs $blddir/clm >&! clm.log.$LID || echo "CLM run failed" && exit 1
else
  $blddir/clm >&! clm.log.$LID                   || echo "CLM run failed" && exit 1
endif

exit 0
