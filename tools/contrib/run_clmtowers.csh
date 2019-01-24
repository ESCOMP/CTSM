#!/bin/csh -f
#
# run_clmtowers.v1.csh
# Purpose: This script will run any number of flux tower sites. You will need to
#          run the script for each spinup and post spinup set of simulations (i.e., 
#          for BGC on, run it separately for AD, PostAD, and post spinup simulations;
#          for BGC off, run it separately for spinup and post spinup simulations)
#          You will need to do two things:
#          1. Copy this script into $Clm_Tag_Dir/components/clm/tools/shared/PTCLM
#             where $Clm_Tag_Dir is the location of your clm tag
#          2. Set up a directory structure where you can put any sourcemods you might want.
#             These sourcemods will be copied into the appropriate case directory.
#             The structure is up to you but here is an example:
#             cd $Clm_Tag_Dir/components/clm/tools/shared/PTCLM
#             mkdir SourceMods
#             mkdir SourceMods/clm4_5 
#             mkdir SourceMods/clm5_0
#             mkdir SourceMods/clm4_5/BASE  ; This might contain any sourcemods that you want
#                                           ;  to use in your clm4_5 control experiment
#             mkdir SourceMods/clm5_0/BASE  ; This might contain any sourcemods that you want
#                                           ;  to use in your clm5_0 control experiment
#             mkdir SourceMods/clm4_5/EXP1  ; This might contain any sourcemods that you want
#                                           ;  to use in your first clm4_5 experiment
#             mkdir SourceMods/clm5_0/EXP1  ; This might contain any sourcemods that you want
#                                           ;  to use in your first clm5_0 experiment
# Author: Keith Oleson
# Last Revised: Aug 25 2015
# Last CLM Tag that this worked on: clm4_5_1_r111
# Warning: This script is complicated and does not have good (any) error checking currently.
#          You might want to come see me for a quick tutorial before using this.
#
# Assumes that PTCLMmkdata has already been run for the tower sites chosen below
# (Surface datasets and shell commands have already been created)
#

set pwd=`pwd`

# =================================Start User Mods================================
# Pick a compset (these are the only two compsets supported)
#set compset = I1PTCLM45
set compset = I1PTCLM50
if ($compset == I1PTCLM45) then
  set model = clm4_5
else
  set model = clm5_0
endif

# Set location of your run directories
set rundata = /glade/scratch/oleson
# Set the location of your CLM tag
set Clm_Tag_Dir   = /glade/p/work/oleson/clm4_5_1_r111
# Set the location of your surface datasets and shell commands. This will not necessarily be in the 
# same location as the CLM tag that you are running above
set User_Mods_Dir = /glade/p/work/oleson/clm4_5_1_r111

# What sites to run?
# These are the sites that can be evaluated with some combination of level 2 data and synthesis (gap-filled) data
#set sites     = ( US-Var US-Bo1 US-UMB US-Brw US-ARM US-Ho1 US-Me4 US-Me2 US-Dk3 US-NR1 DE-Tha ES-ES1 FL-Hyy CA-Man BR-Sa3 BR-Sa1 IT-Cpz US-Dk2 US-MOz US-WCr US-MMS US-Ha1 BE-Vie IT-Col CA-Let US-FPe FL-Kaa US-IB1 US-Ne3 CA-Qfo BR-Sa1LBA BR-Sa3LBA BR-Ma1LBA BR-RJA BR-Ji1 )
#set startyear = ( 2000   1996   1999   1998   2000   1996   1996   2002   1998   1998   1998   1999   1997   1994   2001   2002   2001   2003   2004   1998   1999   1991   1997   1996   1998   2000   2000   2005   2001   2004   2002      2001      2000      2000   1999 )
#set endyear   = ( 2007   2008   2006   2006   2007   2004   2000   2010   2005   2007   2003   2005   2005   2003   2003   2004   2005   2005   2007   2006   2007   2006   2005   2001   2007   2007   2005   2007   2006   2006   2004      2003      2005      2002   2001 )
# Or you could just do one site
set sites     = ( US-Var )
set startyear = ( 2000 )
set endyear   = ( 2007 )

# USER MODS FOR BGC ON
# For BGC on, the sequence of simulations is AD spinup (1000 years),
#                                            PostAD spinup (500 years),
#                                            post spinup (the number of tower years with atmospheric forcing)
# For BGC on, AD spinup is            SPINUP_P1=TRUE,  SPINUP_P2=FALSE
#             PostAD spinup is        SPINUP_P1=TRUE,  SPINUP_P2=TRUE
#             post spinup is          SPINUP_P1=FALSE, SPINUP_P2=FALSE
#set BGC   = "ON"
#setenv SPINUP_P1 "FALSE"
#setenv SPINUP_P2 "FALSE"

# For BGC on, you could use these for either type of spinup (AD or PostAD)
#set newcase   = spinADbgc
#set clonecase = spinPostADbgc
# For BGC on, you could use these for post spinup (the number of tower years with atmospheric forcing)
# You should change the "t111" for the clonecase to whatever tag you are using (e.g., t111 is used here to 
#   denote clm4_5_1_r111) and/or add some designation for your particular experiment with that tag (e.g.,
#   t111exp1wspinbgc)
#set newcase   = spinPostADbgc
#set clonecase = t111wspinbgc

# USER MODS FOR BGC OFF
# For BGC off, the sequence of simulations is normal spinup (32 years)
#                                             post spinup (the number of tower years with atmospheric forcing)
# For BGC off, normal spinup is SPINUP_P1=true,  SPINUP_P2=false
#              post spinup is   SPINUP_P1=false, SPINUP_P2=false 
set BGC   = "OFF"
setenv SPINUP_P1 "FALSE"
setenv SPINUP_P2 "FALSE"

# For BGC off, use these for either normal spinup or post spinup
# You should change the "t111" for the clonecase to whatever tag you are using (e.g., t111 is used here to 
#   denote clm4_5_1_r111) and/or add some designation for your particular experiment with that tag (e.g., 
#   t111exp1wspinsp)
set newcase   = spinsp
set clonecase = t111wspinsp

# These sourcemods will be copied into every case directory (you will need to setup a 
# directory structure for your sourcemods, see instructions at top of script)
set sourcemods_dir = {$Clm_Tag_Dir}/components/clm/tools/shared/PTCLM/SourceMods/
set sourcemods = $sourcemods_dir{$model}/BASE/*.F90 

# =================================End User Mods================================

@ cnt = 1
foreach mysite ( $sites )
  @ numyears = $endyear[$cnt] - $startyear[$cnt] + 1
  if ($SPINUP_P1 == TRUE) then
    @ numfour = $numyears / 4
    # If have three years or less (numfour = 0) just repeat first year  
    # unless first year is leap year then use next year.
    # Since just using one year that is not a leap year we can use 
    # an alignyear of 1 and endyear is the startyear
    if ( $numfour == 0 ) then
      if ( $startyear[$cnt] % 4 == 0 ) then
        @ startyears = $startyear[$cnt] + 1
        @ endyears  = $startyears
      else
        @ endyears = $startyear[$cnt]
        @ startyears = $endyears
      endif
      @ alignyear = 1
    endif
    if ( $numfour != 0 ) then
      @ startyears = $startyear[$cnt]
      @ endyears  = $startyear[$cnt] + $numfour * 4 - 1
      @ alignyear = $startyear[$cnt]
    endif
    echo $endyear[$cnt]
    echo $endyears
    echo $startyears
    echo $alignyear
  endif
  cd {$Clm_Tag_Dir}/cime/scripts
  if ($SPINUP_P1 == FALSE) then
    set casename = ${clonecase}_${mysite}_$compset
    ./create_clone -case $casename -clone ${newcase}_${mysite}_${compset}
  else
    if ($BGC == ON && $SPINUP_P2 == TRUE) then
      set casename = ${clonecase}_${mysite}_$compset
      echo $casename
      ./create_clone -case $casename -clone ${newcase}_${mysite}_${compset}
    else
      set casename = ${newcase}_${mysite}_$compset
      ./create_newcase -user_mods_dir {$User_Mods_Dir}/components/clm/tools/shared/PTCLM/mydatafiles/1x1pt_${mysite} -case $casename -mach yellowstone_intel -compset $compset -res CLM_USRDAT
    endif
  endif
  cd {$Clm_Tag_Dir}/cime/scripts/${casename}
  if ($SPINUP_P1 == FALSE) then
    rm -f cesm.stderr*
    rm -f cesm.stdout*
    rm -f STATUS.out
    ./xmlchange -file env_run.xml -id STOP_OPTION -val nyears
    ./xmlchange -file env_run.xml -id STOP_N -val $numyears
    ./xmlchange -file env_run.xml -id RUN_STARTDATE -val $startyear[$cnt]-01-01
    ./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_ALIGN -val $startyear[$cnt]
    ./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_START -val $startyear[$cnt]
    ./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_END -val $endyear[$cnt]
    ./xmlchange -file env_build.xml -id CALENDAR -val GREGORIAN
    if ($BGC == ON) then
      ./xmlchange -file env_run.xml -id CLM_BLDNML_OPTS -val "-mask navy -bgc bgc -bgc_spinup off"
    endif
  else
    if ($BGC == ON && $SPINUP_P2 == TRUE) then
      rm -f cesm.stderr*
      rm -f cesm.stdout*
      rm -f STATUS.out
      ./xmlchange -file env_run.xml -id STOP_OPTION -val nyears
      ./xmlchange -file env_run.xml -id STOP_N -val 500
      ./xmlchange -file env_run.xml -id CLM_BLDNML_OPTS -val "-mask navy -bgc bgc -bgc_spinup off"
    else
      ./xmlchange -file env_run.xml -id STOP_OPTION -val nyears
      if ($BGC == ON) then
        ./xmlchange -file env_run.xml -id STOP_N -val 1000
      else
        ./xmlchange -file env_run.xml -id STOP_N -val 32
      endif
      if ($alignyear == 1) then
        ./xmlchange -file env_run.xml -id RUN_STARTDATE -val 000{$alignyear}-01-01
      else
        ./xmlchange -file env_run.xml -id RUN_STARTDATE -val $startyear[$cnt]-01-01
      endif
      ./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_ALIGN -val $alignyear
      ./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_START -val $startyears
      ./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_END -val $endyears
      if ($alignyear == 1) then
        ./xmlchange -file env_build.xml -id CALENDAR -val NO_LEAP
      endif
      if ($BGC == ON) then
        ./xmlchange -file env_run.xml -id CLM_BLDNML_OPTS -val "-mask navy -bgc bgc -bgc_spinup on"
      endif
    endif
  endif
  if ($mysite == BR-Sa1LBA || $mysite == BR-Sa3LBA || $mysite == BR-Ma1LBA || $mysite == BR-RJA || $mysite == BR-Ji1) then
    if ($SPINUP_P1 == FALSE) then
      rm -f user_datm.streams.txt.CLM1PT.CLM_USRDAT
    endif
  endif
  ./cesm_setup
  if ( $status != 0 )then
     echo "CESM_SETUP FAIL $status" >> ./STATUS.out
     exit -1
  else
     echo "CESM_SETUP PASS" >> ./STATUS.out
  endif
  if ($SPINUP_P1 == TRUE && $BGC == ON) then
    sed "/BSUB -R/d;s/BSUB -P P93300641/BSUB -P P93300041/g;s/BSUB -W 2:00/BSUB -W 23:59/g" ./${casename}.run > tmp.run
  else
    sed "/BSUB -R/d;s/BSUB -P P93300641/BSUB -P P93300041/g" ./${casename}.run > tmp.run
  endif
  mv tmp.run ./${casename}.run
  chmod u+x ./${casename}.run
  cp $sourcemods SourceMods/src.clm
  echo $mysite
  if ($mysite == BR-Sa1LBA || $mysite == BR-Sa3LBA || $mysite == BR-Ma1LBA || $mysite == BR-RJA || $mysite == BR-Ji1) then
    cp CaseDocs/datm.streams.txt.CLM1PT.CLM_USRDAT ./user_datm.streams.txt.CLM1PT.CLM_USRDAT
    chmod u+wx ./user_datm.streams.txt.CLM1PT.CLM_USRDAT
    sed "s/RH   /QBOT /g" user_datm.streams.txt.CLM1PT.CLM_USRDAT > tmp.user_datm.streams.txt.CLM1PT.CLM_USRDAT
    sed "s/  rh/  shum/g" tmp.user_datm.streams.txt.CLM1PT.CLM_USRDAT > tmp2.user_datm.streams.txt.CLM1PT.CLM_USRDAT
    rm -f tmp.user_datm.streams.txt.CLM1PT.CLM_USRDAT
    mv tmp2.user_datm.streams.txt.CLM1PT.CLM_USRDAT ./user_datm.streams.txt.CLM1PT.CLM_USRDAT
  endif
  if ($SPINUP_P1 == FALSE) then
    sed "/taxmode = 'cycle','cycle'/d" user_nl_datm > tmp.user_nl_datm
    mv tmp.user_nl_datm ./user_nl_datm
    sed "s/hist_nhtfrq = 0/hist_nhtfrq = 0,1/g" ./user_nl_clm > tmp.user_nl_clm
    sed "s/hist_mfilt  = 1200/hist_mfilt  = 1,350400/g" ./tmp.user_nl_clm > tmp2.user_nl_clm
    rm -f tmp.user_nl_clm
    sed "/finidat/d" ./tmp2.user_nl_clm > tmp3.user_nl_clm
    rm -f tmp2.user_nl_clm
    if ($BGC == ON) then
      echo " hist_fincl2  = 'FSDS','FLDS','FSR','FSA','FIRE','FIRA','FSH','FCTR','FCEV','FGEV','FGR','FGR12','FSM','TSOI','COSZEN','RAIN','SNOW','H2OSOI','WA','ZWT','GPP','NEE','ELAI','BTRAN','TV','RSSUN','RSSHA','FSH_G','RHAF','RH_LEAF','RH','T10','TG','SABG','SABV'" >> tmp3.user_nl_clm
    else
      echo " hist_fincl2  = 'FSDS','FLDS','FSR','FSA','FIRE','FIRA','FSH','FCTR','FCEV','FGEV','FGR','FGR12','FSM','TSOI','COSZEN','RAIN','SNOW','H2OSOI','WA','ZWT','ELAI','BTRAN','FPSN','TV','RSSUN','RSSHA','FSH_G','RHAF','RH_LEAF','RH','T10','TG','SABG','SABV'" >> tmp3.user_nl_clm
    endif
    set finidat=`ls -1 $rundata/${newcase}_${mysite}_${compset}/run/${newcase}_${mysite}_${compset}.clm?.r.*.nc | tail -1`
    echo $finidat
    echo " finidat = '$finidat'" >> tmp3.user_nl_clm
    mv tmp3.user_nl_clm ./user_nl_clm
  else
    if ($BGC == ON && $SPINUP_P2 == TRUE) then 
      set finidat=`ls -1 $rundata/${newcase}_${mysite}_${compset}/run/${newcase}_${mysite}_${compset}.clm?.r.*.nc | tail -1`
      echo $finidat
      echo " finidat = '$finidat'" >> user_nl_clm
    else
      echo "taxmode = 'cycle','cycle'" >> user_nl_datm
      if ($BGC == ON) then
        sed "s/hist_mfilt  = 1200/hist_mfilt  = 12000/g" ./user_nl_clm > tmp.user_nl_clm
        mv tmp.user_nl_clm ./user_nl_clm
      endif
    endif
  endif
  ./{$casename}.build
  if ( $status != 0 )then
     echo "BUILD FAIL $status" >> ./STATUS.out
     exit -1
  else
     echo "BUILD PASS" >> ./STATUS.out
  endif
  ./${casename}.submit
  if ( $status != 0 )then
     echo "SUBMIT FAIL $status" >> ./STATUS.out
     exit -1
  else
     echo "SUBMIT PASS" >> ./STATUS.out
  endif
  cd $pwd
  @ cnt++
end
