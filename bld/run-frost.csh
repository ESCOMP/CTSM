#! /bin/tcsh -f
#
#=======================================================================
#
#  run-frost.csh
#
#  Generic batch submission script for frost using LSF.  
#
#-----------------------------------------------------------------------
# Run this job interactively on frost,
#   ./run-frost.csh
#-----------------------------------------------------------------------
#
#XXXX -n 32
#XXXX -t 00:40:00
#
#=======================================================================

setenv INC_NETCDF /contrib/bgl/netcdf-3.6.1/include
setenv LIB_NETCDF /contrib/bgl/netcdf-3.6.1/lib

setenv INC_MPI /bgl/BlueLight/ppcfloor/bglsys/include
setenv LIB_MPI /bgl/BlueLight/ppcfloor/bglsys/lib

## ROOT OF CLM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CLM distribution.
## (the root directory contains the subdirectory "src")
set clmroot   = ~/fmf_clm3_expa_88

## ROOT OF CLM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CLM distribution.
setenv CSMDATA /ptmp/tcraig/inputdata/lnd/clm2

## Default configuration settings:
set spmd     = on       # settings are [on   | off       ] (default is off)
set maxpft   = 17       # settings are 4->17               (default is 4)
set bgc      = none     # settings are [none | cn | casa ] (default is none)
set supln    = off      # settings are [on   | off       ] (default is off)
set rtm      = off      # settings are [on   | off       ] (default is off)
set dust     = off      # settings are [on   | off       ] (default is off)   
set voc      = off      # settings are [on   | off       ] (default is off)   

## $wrkdir is a working directory where the model will be built and run.
## $blddir is the directory where model will be compiled.
## $rundir is the directory where the model will be run.
## $cfgdir is the directory containing the CLM configuration scripts.
set case     = clmrun05
set wrkdir   = /ptmp/$LOGNAME
set blddir   = $wrkdir/$case/bld
set rundir   = $wrkdir/$case
set cfgdir   = $clmroot/bld
set usr_src  = $clmroot/bld/usr.src


## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## Build (or re-build) executable
set flags = "-maxpft $maxpft -bgc $bgc -supln $supln -rtm $rtm -voc $voc -dust $dust -usr_src $usr_src -target_os bgl"
if ($spmd == on) set flags = "$flags -spmd"
echo "cd $blddir"
cd $blddir                  || echo "cd $blddir failed" && exit 1
set config="$blddir/config_cache.xml"
if ( ! -f "$config" )then
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

set fatmgrid       = `$cfgdir/queryDefaultNamelist.pl -res $res -silent -config $config -csmdata $CSMDATA -var fatmgrid`
set fsurdat        = `$cfgdir/queryDefaultNamelist.pl -res $res -silent -config $config -csmdata $CSMDATA -var fsurdat`
set fatmlndfrc     = `$cfgdir/queryDefaultNamelist.pl -res $res -silent -config $config -csmdata $CSMDATA -var fatmlndfrc`
set fpftcon        = `$cfgdir/queryDefaultNamelist.pl -silent -config $config -csmdata $CSMDATA -var fpftcon`
set frivinp_rtm    = `$cfgdir/queryDefaultNamelist.pl -silent -config $config -csmdata $CSMDATA -var frivinp_rtm`
set offline_atmdir = `$cfgdir/queryDefaultNamelist.pl -silent -config $config -csmdata $CSMDATA -var offline_atmdir`


cat >! lnd.stdin << EOF
 &clm_inparm
 caseid         = '$case'
 ctitle         = '$case'
 finidat        = ' '
 $fsurdat
 $fatmgrid
 $fatmlndfrc
 $fpftcon
 $frivinp_rtm
 $offline_atmdir
 nsrest         =  0
 nelapse        =  48
 dtime          =  1800
 start_ymd      =  19980101
 start_tod      =  0
 irad           = -1
 wrtdia         = .true.
 mss_irt        =  0
 hist_dov2xy    = .true.
 hist_nhtfrq    =  1000
 hist_mfilt     =  1
 hist_crtinic   = 'MONTHLY'
 rest_flag      = .false.
 /
 &prof_inparm
 /
EOF

## Run CLM 
cd $rundir                              || echo "cd $rundir failed" && exit 1
echo "running CLM in $rundir"

cqsub -e LOGNAME=$LOGNAME -n 32 -t 00:20:00 $blddir/clm
#$blddir/clm                   || echo "CLM run failed" && exit 1


exit 0
