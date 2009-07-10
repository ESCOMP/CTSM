#! /usr/bin/csh -f

#-----------------------------------------------------------------------
## IBM
##------------
##
## This is an example script to build and run the default CLM configuration
## on an IBM SP.  This is setup to run on NCAR's machine bluefire.
##
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
## NOTE: FOR PRODUCTION RUNS -- USE THE create_newcase SCRIPT in ../../../../scripts!!!
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
##
## To submit this on bluefire:
##
## bsub < run-ibm.csh
##
## To run interatively: 
##
##    first set spmd=off and smp=on below then
##    ./run-ibm.csh
##

## Setting LSF options for batch queue submission.
#BSUB -a poe                    # use poe for multiprocessing
#BSUB -x                        # exclusive use of node (not_shared)
## Number of tasks and tasks per node (CHANGE THIS IF YOU TURN smp on)
#BSUB -n 32                     # total number of MPI-tasks (processors) needed
#BSUB -R "span[ptile=64]"       # max number of tasks (MPI) per node
#BSUB -o out.%J                 # output filename
#BSUB -e out.%J                 # error filename
#BSUB -q regular                # queue
#BSUB -W 0:10                   # wall clock limit
#BSUB -P xxxxxxxx               # Project number to charge to (MAKE SURE YOU CHANGE THIS!!!)

#===========================================================================
#=================== THINGS MOST COMONLY CHANGED ===========================

## Configuration settings:
set spmd     = off      # settings are [on   | off         ] (default is on) (off for interactive)
set smp      = on       # settings are [on   | off         ] (default is off)
set maxpft   = numpft+1 # settings are 4->17                 (default is 17)
set bgc      = none     # settings are [none | cn | casa   ] (default is none)
set supln    = off      # settings are [on   | off         ] (default is off)
set dust     = on       # settings are [on   | off         ] (default is on)   
set voc      = on       # settings are [on   | off         ] (default is on)   
set rtm      = on       # settings are [on   | off         ] (default is on)   

## IF YOU CHANGE ANY OF THE CONFIGURATION SETTINGS -- DELETE THE $blddir/config_cache.xml AND RESUBMIT
## (see below)
#--------------------------------------------------------------------------------------------
## Run time settings:
## May also make changes to namelist in build-namelist section below:
set res        = 4x5        # settings are [48x96   | 64x128  | 4x5  | 10x15 | 1.9x2.5 etc.      ]
set mask       = default    # settings are [default | USGS    | navy | gx3v5 | gx1v6   etc.      ]
set sim_year   = default    # settings are [default | 1850    | 2000                             ]
set start_type = arb_ic     # settings are [cold    | arb_ic  | startup | continue | branch      ] (default is arb_ic)
                            # for branch type you need to enter needed files by hand, and make sure they are available
set runlen     = 2d         # settings are [ integer<sdy> where s=cpling-step, d=days, y=years   ] (default is 2d)
set start_ymd  = 19980101   # Start date [yyyymmdd]
set cycle_init = 1998       # Initial year to use atm data from
set cycle_beg  = 1948       # Begining year to cycle through for input atm data
set cycle_end  = 2004       # Ending   year to cycle through for input atm data
#--------------------------------------------------------------------------------------------
## Locations of important directories:
##

## netCDF stuff (MAKE SURE YOU CHANGE THIS -- IF NOT RUNNING AT NCAR)
setenv INC_NETCDF /usr/local/include
setenv LIB_NETCDF /usr/local/lib64/r4i4

## CCSM machine name
set ccsm_mach = bluefire

## ROOT OF CLM DISTRIBUTION - probably needs to be customized (changed by create_newcase)
## Contains the source code for the CLM distribution.
## (the root directory for CLM contains the subdirectory "src")
set clmroot   = /fis/cgd/.......       # (MAKE SURE YOU CHANGE THIS!!!)

## ROOT OF CLM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CLM distribution.
setenv CSMDATA /fs/cgd/csm/inputdata                # (MAKE SURE YOU CHANGE THIS!!!)

## Location of datm data - needs to be customized unless running at NCAR.
## Contains the location for the datm7 input data
setenv datm_data_dir /cgd/tss/atm_forcing.datm7.Qian.T62.c080727  # (MAKE SURE YOU CHANGE THIS!!!)

## $wrkdir  is a working directory where the model will be built and run.
## $blddir  is the directory where model will be compiled.
## $rundir  is the directory where the model will be run.
## $cfgdir  is the directory containing the CLM configuration scripts.
## $casdir  is the directory this script is found in.
## $macdir  is the directory the ccsm machines files are found in
set case    = clmrun                  # changed by create_newcase
set wrkdir  = /ptmp/$LOGNAME          # changed by create_newcase
set blddir  = $wrkdir/$case/bld
set rundir  = $wrkdir/$case
set cfgdir  = $clmroot/bld
set casdir  = $clmroot/bld            # changed by create_newcase
set macdir  = $clmroot/../../../scripts/ccsm_utils/Machines

# Number of threads to use:
#
# should be set equal to (CPUs-per-node / tasks_per_node)
# Only activated if smp=on above
setenv OMP_NUM_THREADS 32

#=================== END OF THINGS MOST COMONLY CHANGED ====================
#===========================================================================

## Invoke ccsm machine vars
source $macdir/env_machopts.$ccsm_mach

## Do our best to get sufficient stack memory
limit stacksize unlimited

echo "This script is deprecated."
echo "NOTE: FOR PRODUCTION RUNS -- USE THE create_newcase SCRIPT in ../../../../scripts\!\!\!"


## Ensure that run and build directories exist
mkdir -p $rundir                    || echo "cannot create $rundir"        && exit 1
mkdir -p $rundir/timing/checkpoints || echo "cannot create $rundir/timing" && exit 1
mkdir -p $blddir                    || echo "cannot create $blddir"        && exit 1

## Build (or re-build) executable
set flags = "-maxpft $maxpft -bgc $bgc -supln $supln -voc $voc -rtm $rtm -dust $dust "
set flags = "$flags -mach $ccsm_mach"
if ($spmd == on ) set flags = "$flags -spmd"
if ($spmd == off) set flags = "$flags -nospmd"
if ($smp  == on ) set flags = "$flags -smp"
if ($smp  == off) set flags = "$flags -nosmp"

echo "cd $blddir"
cd $blddir                  || echo "cd $blddir failed" && exit 1

set config="$blddir/config_cache.xml"
## Check if config_cache.xml file exists -- if so just run make -- if NOT then run configure.
## IF YOU CHANGE ANY OF THE CONFIGURATION SETTINGS -- DELETE THE $blddir/config_cache.xml AND RESUBMIT
#--------------------------------------------------------------------------------------------
if ( ! -f $config ) then
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

#
# If you want to include a specific input file to clm simply add it to the clm_inparm namelist below
# e.g.
#    finidat = '$CSMDATA/lnd/clm2/inidata_3.1/offline/clmi_0000-01-01_064x128_c070403.nc'
#
set dtime = 1800
cat >! lndinput << EOF
 &drv_in
 start_ymd      =  $start_ymd
 start_tod      =  0
 atm_cpl_dt     =  $dtime
 /
 &clm_inparm
 dtime          =  $dtime
 wrtdia         = .true.
 hist_dov2xy    = .true.
 hist_nhtfrq    =  -24
 hist_mfilt     =  1
 hist_crtinic   = 'NONE'
 /
EOF
set bnflags="-case $case -start_type $start_type -config $config -mask $mask -sim_year $sim_year -infile lndinput -runlength $runlen"
set bnflags="$bnflags -datm_data_dir $datm_data_dir -res $res -csmdata $CSMDATA"
set bnflags="$bnflags -cycle_init $cycle_init -cycle_beg_year $cycle_beg -cycle_end_year $cycle_end"
$cfgdir/build-namelist $bnflags    || echo "build-namelist failed" && exit 1

## Run CLM 

cd $rundir                    || echo "cd $rundir failed" && exit 1
setenv LID "`date +%y%m%d-%H%M%S`"

echo "running CLM in $rundir log file out to $rundir/clm.log.$LID"

if ($spmd == on && $smp == on) then
  setenv TARGET_CPU_LIST "-1"
  set launch = "/usr/local/bin/hybrid_launch"
else if ($spmd == on) then
  setenv TARGET_CPU_RANGE "-1"
  set launch = "/usr/local/bin/launch"
else
  set launch = ""
endif
if ( ! -f "$launch" ) set launch = ""

if ($spmd == on) then

  mpirun.lsf $launch $blddir/clm >&! clm.log.$LID       || echo "CLM run failed" && exit 1
  set runstatus = $status
else 
  $blddir/clm  >&! clm.log.$LID                 || echo "CLM run failed" && exit 1
  set runstatus = $status
endif

exit 0
