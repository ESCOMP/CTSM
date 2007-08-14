#! /bin/tcsh -f
#
#=======================================================================
#
#  run-pc.csh
#
#  Generic batch submission script for PC-linux using PBS.  
#
#-----------------------------------------------------------------------
# Batch options for machine with PBS batch system. (bangkok/calgary)
# Usage for Lahey compiler (default): 
#   qsub run-pc.csh
# Usage for pgf90 compiler with pgcc: 
#   env OPT_BLD=PGI qsub run-pc.csh
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
       echo "${0}: WARNING: Running CLM interactively -- resubmit job to batch que"
       sleep 5
     else
       echo "${0}: Running CLM on Linux using PBS";
     endif
     if ( ! $?OPT_BLD ) then
       setenv USER_FC lf95
       setenv LAHEY /usr/local/lf9562
       set netcdf = /usr/local/netcdf-gcc-lf95
       setenv INC_NETCDF ${netcdf}/include
       setenv LIB_NETCDF ${netcdf}/lib
       set mpich = /usr/local/mpich-gcc-g++-lf95
       setenv INC_MPI ${mpich}/include
       setenv LIB_MPI ${mpich}/lib
       setenv PATH ${LAHEY}/bin:${mpich}/bin:${PATH}
       setenv MOD_NETCDF $INC_NETCDF
       setenv LD_LIBRARY_PATH "${LAHEY}/lib:/lib/i686:/usr/local/lib"
     else
       set mpich = /usr/local/mpich-pgi-pgcc-pghf
       setenv INC_MPI ${mpich}/include
       setenv LIB_MPI ${mpich}/lib
       set netcdf = /usr/local/netcdf-pgi-hpf-cc
       setenv INC_NETCDF ${netcdf}/include
       setenv LIB_NETCDF ${netcdf}/lib
       setenv MOD_NETCDF $INC_NETCDF
       setenv PGI /usr/local/pgi-pgcc-pghf
       setenv PATH ${PGI}/linux86/6.1/bin:${mpich}/bin:${PATH}
       setenv LD_LIBRARY_PATH "${PGI}/linux86/6.1/lib:${LD_LIBRARY_PATH}"

     endif
     breaksw;
  default:
    echo "${0}: This script meant for running CLM on Linux machines";    exit 2;
endsw

## set this equal to #nodes X #ppn
set procs = 4

## Do our best to get sufficient stack memory
limit stacksize unlimited

## ROOT OF CLM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CLM distribution.
## (the root directory contains the subdirectory "src")
set clmroot   = /fs/cgd/...                # (MAKE SURE YOU CHANGE THIS!!!)

## ROOT OF CLM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CLM distribution.
setenv CSMDATA /fs/cgd/csm/inputdata                # (MAKE SURE YOU CHANGE THIS!!!)

## Configuration settings:

set mode     = offline  # settings are [offline | ccsm_seq ] (default is offline)
set spmd     = on       # settings are [on   | off         ] (default is off)
set maxpft   = 4        # settings are 4->17                 (default is 4)
set bgc      = none     # settings are [none | cn | casa   ] (default is none)
set supln    = off      # settings are [on   | off         ] (default is off)
set dust     = off      # settings are [on   | off         ] (default is off)   
set voc      = off      # settings are [on   | off         ] (default is off)   
set rtm      = off      # settings are [on   | off         ] (default is off)   
## IF YOU CHANGE ANY OF THE CONFIGURATION SETTINGS -- DELETE THE $blddir/config_cache.xml AND RESUBMIT
## (See below)
#--------------------------------------------------------------------------------------------

## Run time settings:
set res        = 48x96    # settings are [48x96   | 64x128  | 4x5  | 10x15 | 1.9x2.5 etc. ]
set mask       = default  # settings are [default | USGS    | navy | gx3v5 | gx1v5   etc. ]
set sim_year   = default  # settings are [default | 1890    | 2000 | 2100                 ]
set start_type = arb_ic   # settings are [arb_ic  | startup | continue | branch           ] (default is arb_ic)
set runlen     = 2d       # settings are [ integer<sdy> where s=step, d=days, y=years     ] (default is 2d)

## $wrkdir is a working directory where the model will be built and run.
## $blddir is the directory where model will be compiled.
## $rundir is the directory where the model will be run.
## $cfgdir is the directory containing the CLM configuration scripts.
set case     = clmrun
set wrkdir   = /ptmp/$LOGNAME
set blddir   = $wrkdir/$case/bld
set rundir   = $wrkdir/$case
set cfgdir   = $clmroot/bld
set usr_src  = $clmroot/bld/usr.src

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 3
mkdir -p $blddir                || echo "cannot create $blddir" && exit 3

## Build (or re-build) executable
set flags = "-maxpft $maxpft -bgc $bgc -supln $supln -voc $voc -rtm $rtm -dust $dust -usr_src $usr_src"
if ($spmd == on ) set flags = "$flags -spmd"
if ($spmd == off) set flags = "$flags -nospmd"

echo "cd $blddir"
cd $blddir                  || echo "cd $blddir failed" && exit 3

set config = "$blddir/config_cache.xml"
## Check if config_cache.xml file exists -- if so just run make -- if NOT then run configure.
## IF YOU CHANGE ANY OF THE CONFIGURATION SETTINGS -- DELETE THE $blddir/config_cache.xml AND RESUBMIT
#--------------------------------------------------------------------------------------------
if ( ! -f $config ) then
    echo "flags to configure are $flags"
    $cfgdir/configure $flags    || echo "configure failed" && exit 4
    echo "Building CLM in $blddir ..."
    gmake -j8 >&! MAKE.out      || echo "CLM build failed: see $blddir/MAKE.out" && exit 5
else
    echo "Re-building CLM in $blddir ..."
    rm -f Depends
    gmake -j8 >&! REMAKE.out      || echo "CLM build failed: see $blddir/REMAKE.out" && exit 5
endif

## Create the namelist
cd $rundir                      || echo "cd $blddir failed" && exit 6

#
# If you want to include a specific input file to clm simply add it to the clm_inparm namelist below
# e.g.
#    finidat = '$CSMDATA/lnd/clm2/inidata_3.1/offline/clmi_0000-01-01_064x128_c070403.nc'
#
cat >! lndinput << EOF
 &drv_in
 start_ymd      =  19980101
 start_tod      =  0
 /
 &clm_inparm
 dtime          =  1800
 irad           = -1
 wrtdia         = .true.
 mss_irt        =  0
 hist_dov2xy    = .true.
 hist_nhtfrq    =  -24
 hist_mfilt     =  1
 hist_crtinic   = 'MONTHLY'
 /
EOF
set bnflags="-case $case -start_type $start_type -config $config -mask $mask -sim_year $sim_year -infile lndinput -runlength $runlen" 
$cfgdir/build-namelist $bnflags    || echo "build-namelist failed" && exit 1

## Run CLM
cd $rundir                      || echo "cd $rundir failed" && exit 6
setenv LID "`date +%y%m%d-%H%M%S`"

echo "running CLM in $rundir log output to $rundir/clm.log.$LID"


if ($spmd == on) then
  mpirun -np $procs $blddir/clm >&! clm.log.$LID || echo "CLM run failed" && exit 7
else
  $blddir/clm >&! clm.log.$LID                   || echo "CLM run failed" && exit 7
endif

exit 0
