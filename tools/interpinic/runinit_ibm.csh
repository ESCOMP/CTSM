#! /usr/bin/csh -f

#-----------------------------------------------------------------------
## IBM
##------------
##
## This is a script to do short 1-step initial runs with CLM in order
## to create initial condition files that can then be interpolated into
## to create files to start from.
##
## This is setup to run on NCAR's machine bluefire.
##
## To run this on bluefire:
##
##    ./runinit_ibm.csh
##
#BSUB -P 03010999           # project number
#BSUB -x                    # exclusive use of node (only enable for non-share que)
#BSUB -n 1                  # use 1 to 64 processors on one node
#BSUB -R "span[ptile=1]"    # only 1 task per node
#BSUB -o %J.out             # ouput filename
#BSUB -e %J.err             # input filename
#BSUB -J runinit_ibm        # job name
#BSUB -W 6:00               # Wall clock time limit
#BSUB -q regular            # queue


#===========================================================================
#=================== CONFIGURATION THAT ARE SET  ===========================

## Configuration settings:
set spmd     = off      # settings are [on   | off         ] (default is on) (off for interactive)
set smp      = on       # settings are [on   | off         ] (default is off)
set supln    = off      # settings are [on   | off         ] (default is off)
set dust     = on       # settings are [on   | off         ] (default is off)
set seaslt   = on       # settings are [on   | off         ] (default is off)
set voc      = off      # settings are [on   | off         ] (default is off)   
set rtm      = on       # settings are [on   | off         ] (default is off)   

#--------------------------------------------------------------------------------------------
## Run time settings that are set for all cases:
## May also make changes to namelist in build-namelist section below:
set start_type = cold       # settings are [cold    | arb_ic  | startup | continue | branch      ] (default is arb_ic)
set runlen     = 1s         # settings are [ integer<sdy> where s=cpling-step, d=days, y=years   ] (default is 2d)
#--------------------------------------------------------------------------------------------
## Locations of important directories:
##

## netCDF stuff
setenv INC_NETCDF /usr/local/include
setenv LIB_NETCDF /usr/local/lib64/r4i4

## ROOT OF CLM DISTRIBUTION
## Contains the source code for the CLM distribution.
## (the root directory for CLM contains the subdirectory "src")
## UTILROOT is the root of the CCSM tools directory
set curdir    = `pwd`
set clmroot   = $curdir/../..
setenv UTILROOT $clmroot/../../../scripts/ccsm_utils
set ccsm_mach = "bluefire"

## ROOT OF CLM DATA DISTRIBUTION
## Contains the initial and boundary data for the CLM distribution.
setenv CSMDATA /fs/cgd/csm/inputdata

## Location of datm data
## Contains the location for the datm7 input data
setenv datm_data_dir /cgd/tss/atm_forcing.datm7.Qian.T62.c080727

## $wrkdir  is a working directory where the model will be built and run.
## $cfgdir  is the directory containing the CLM configuration scripts.
set wrkdir  = /ptmp/$LOGNAME          # changed by create_newcase
set cfgdir  = $clmroot/bld

# Number of threads to use:
#
# should be set equal to (CPUs-per-node / tasks_per_node)
# Only activated if smp=on above
setenv OMP_NUM_THREADS 64

## env variables
set CCSMUSER=$USER
source $UTILROOT/Machines/env_machopts.$ccsm_mach

## Do our best to get sufficient stack memory
limit stacksize unlimited

set sdate = "c"`date +%y%m%d`

## Create a script file that copies files to proper place and also keeps the list of filenames
set datalist = clm.input_data_list
if ( -f "$datalist" ) mv -f $datalist ${datalist}.previous
touch $datalist
echo '#! /bin/csh -f'              >> $datalist
echo "set DIN_LOC_ROOT = $CSMDATA" >> $datalist
chmod +x $datalist
set svnrepo = "https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata"
set svnmesg = 'Update finidat files with interpinic $HeadURL$ $Id$'
set initdir = "lnd/clm2/initdata"

# make interpinic in local directory

cd $curdir
gmake clean
gmake -j 65 OPT=TRUE SMP=TRUE interpinic

date

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#------------------ Loop over different and configuration types -------------------------
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#foreach bgc ( "none" "cn" )
foreach bgc ( "cn" )
   set maxpft   = "numpft+1"
   if (      "$bgc" == "none" )then
      set compset = "IQ"
   else if ( "$bgc" == "casa" )then
      set compset = "IQCASA"
   else if ( "$bgc" == "cn"   )then
      set compset = "IQCN"
   endif
   set bldcase = "clmbld_$bgc"
   set blddir  = $wrkdir/$bldcase/bld
   mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

   ## Build (or re-build) executable
   set flags = "-maxpft $maxpft -bgc $bgc -supln $supln -voc $voc -rtm $rtm -dust $dust "
   set flags = "$flags -prog_seasalt $seaslt -mach $ccsm_mach -ccsm_bld on -mode ccsm_seq"
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

   #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   #---------------- Loop over different resolutions and run-time configurations ---------
   #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   foreach res ( "1.9x2.5" "10x15" "4x5" "0.9x1.25" "2.5x3.33" "0.47x0.63" "48x96" )
      set fres    = `echo $res | tr '.' '_'`
      set case    = "clmrun_${bgc}_${fres}"
      set rundir  = $wrkdir/$case
      mkdir -p $rundir/timing         || echo "cannot create $rundir" && exit 1
      if ( "$res" == "48x96" || "$res" == "4x5" )then
         set masks = ( "gx3v7" )
      else if ( "$res" == "1.9x2.5" ) then
         set masks = ( "gx1v6" )
      else if ( "$res" == "1.9x2.5" || "$res" == "0.47x0.63" || "$res" == "0.23x0.31" \
             || "$res" == "0.9x1.25" || "$res" == "64x128" || "$res" == "48x96" )then
         set masks = ( "gx1v6" )
      else
         set masks = ( "USGS" )
      endif

      set sim_years = ("1850" "2000" )

      if ( "$res" == "0.47x0.63" && "$bgc" == "cn" ) then
         set outnc = "outnc_large_files = .true."
      else
         set outnc = " "
      endif

      #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      #---------------- Loop over different masks -------------------------------------------
      #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      foreach mask ( $masks )


         #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         #---------------- Loop over different sim_years   -------------------------------------
         #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

         foreach sim_year ( $sim_years )

            if (      "$sim_year" == "1850" )then
               set start_times = ( 19481231 )
            else if ( "$sim_year" == "2000" )then
               set start_times = ( 19991231 )
            endif


            #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            #---------------- Loop over different start times -------------------------------------
            #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   
            foreach start_ymd ( $start_times )
   
               @ cycle_init   = $start_ymd / 10000
               set cycle_beg  = $cycle_init
               @ cycle_end    = $cycle_init + 1
   
               ## Create the namelist
               cd $rundir                      || echo "cd $blddir failed" && exit 1
      
               cat >! lndinput << EOF
 &drv_in
 start_ymd      = $start_ymd
 start_tod      = 84600
 atm_cpl_dt     = 1800
 $outnc
 /
 &clm_inparm
 dtime          =  1800
 hist_crtinic   = 'NONE'
 /
EOF
               set bnflags="-drv_case $case -clm_start_type $start_type -config $config -mask $mask -sim_year $sim_year -infile lndinput"
               set bnflags="$bnflags -datm_data_dir $datm_data_dir -res $res -csmdata $CSMDATA -drv_runlength $runlen"
               set bnflags="$bnflags -datm_cycle_init $cycle_init -datm_cycle_beg_year $cycle_beg -datm_cycle_end_year $cycle_end"
               set bnflags="$bnflags -test"
               echo "Build-namelist flags: $bnflags"
               $cfgdir/build-namelist $bnflags    || echo "build-namelist failed" && exit 1

               ## Run CLM 
   
               setenv LID "`date +%y%m%d-%H%M%S`"

               echo "running CLM in $rundir log file out to $rundir/clm.log.$LID"

               $blddir/clm  >&! clm.log.$LID                 || echo "CLM run failed" && exit 1
               tail -20 clm.log.$LID | grep SUCCESS
               if ( $status ) then
                  tail clm.log.$LID
                  echo "CLM run failed"
                  exit 1
               endif
        
               echo "Copy $case.clm?.r.*.nc to $curdir"
               set enddate = `echo $case.clm?.r.*.nc | awk -F. '{print $4}' | awk -F- '{print $1"-"$2"-"$3}'`
               set outfile = "clmi.${compset}.${enddate}_${res}_${mask}_simyr${sim_year}_${sdate}.nc"
               echo "mv $case.clm?.r.*.nc $curdir/$outfile"
               mv $case.clm?.r.*.nc $curdir/$outfile
               if ( $status ) then
                  echo "could NOT copy file to interpinic directory"
                  exit 1
               endif

               # Run interpinic on resulting file
               cd $curdir
               set finidat = `../../bld/queryDefaultNamelist.pl -res 0.9x1.25 -options mask=gx1v6,sim_year=$sim_year -onlyfiles -justvalue -var finidat -s -config $config`

               echo "Run interpinic to interpolate from $finidat to $outfile"

               ./interpinic -i $finidat -o $outfile > interpinic.log.$LID || echo "interpinic failed" && exit 1
               tail interpinic.log.$LID | grep Success
               if ( $status ) then
                  tail interpinic.log.$LID
                  echo "Interpinic failed"
                  exit 1
               endif
               cat >> $datalist << EOF
if ( "\$1" != "import" ) \cp -p $outfile \$DIN_LOC_ROOT/$initdir/$outfile
#file = \$DIN_LOC_ROOT/$initdir/$outfile
if ( "\$1" == "import" ) svn import -m '$svnmesg' \$DIN_LOC_ROOT/$initdir/$outfile $svnrepo/$initdir/$outfile
EOF
            end
         end

      end

   end
end
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#------------------ End Loop over different resolutions and configuration types ----------
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

date

echo "Successfully ran interpinic on needed resolutions and configuration types"

exit 0
