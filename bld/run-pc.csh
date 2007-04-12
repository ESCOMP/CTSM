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
     if ( ! $?OPT_BLD ) then
       setenv USER_FC lf95
       setenv LAHEY /usr/local/lf9562
       set netcdf = /usr/local/netcdf-3.6.1beta3-gcc-4.0.2-g77-lf9562
       setenv INC_NETCDF ${netcdf}/include
       setenv LIB_NETCDF ${netcdf}/lib
       set mpich = /usr/local/mpich-1.2.7p1-gcc-g++-4.0.2-8-lf9562
       setenv INC_MPI ${mpich}/include
       setenv LIB_MPI ${mpich}/lib
       setenv PATH ${LAHEY}/bin:${mpich}/bin:${PATH}
       setenv MOD_NETCDF $INC_NETCDF
       setenv LD_LIBRARY_PATH "${LAHEY}/lib:/lib/i686:/usr/local/lib"
     else
       if ( $OPT_BLD == "pgf90-gcc" ) then
         setenv USER_CC gcc
       endif
       setenv INC_MPI /usr/local/mpich-pgi-pgcc-pghf/include
       setenv LIB_MPI /usr/local/mpich-pgi-pgcc-pghf/lib
       setenv INC_NETCDF /usr/local/netcdf-pgi-hpf-cc/include
       setenv MOD_NETCDF $INC_NETCDF
       setenv LIB_NETCDF /usr/local/netcdf-pgi-hpf-cc/lib
       setenv PGI /usr/local/pgi-pgcc-pghf
       setenv PATH ${PGI}/linux86/6.1/bin:${PATH}
       setenv LD_LIBRARY_PATH "${PGI}/linux86/lib:${PGI}/linux86/liblf:/lib/i686:/usr/local/lib"

     endif
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
set clmroot   = /fis/cgd/...

## ROOT OF CLM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CLM distribution.
setenv CSMDATA /fs/cgd/csm/inputdata/lnd/clm2

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
set case     = clmrun
set wrkdir   = /ptmp/$LOGNAME
set blddir   = $wrkdir/$case/bld
set rundir   = $wrkdir/$case
set cfgdir   = $clmroot/bld
set usr_src  = $clmroot/bld/empty

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## Build (or re-build) executable
set flags = "-maxpft $maxpft -bgc $bgc -supln $supln -rtm $rtm -voc $voc -dust $dust -usr_src $usr_src"
if ($spmd == on ) set flags = "$flags -spmd"
if ($spmd == off) set flags = "$flags -nospmd"
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
 caseid         = $case
 ctitle         = $case
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
cd $rundir                      || echo "cd $rundir failed" && exit 1
echo "running CLM in $rundir"

if ($spmd == on) then
  mpiexec -n $procs $blddir/clm || echo "CLM run failed" && exit 1
else
  $blddir/clm                   || echo "CLM run failed" && exit 1
endif

exit 0
