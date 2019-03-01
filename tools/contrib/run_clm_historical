#!/bin/csh

#########################################################################################
#
#  - Execute this script to do a CLM historical simulation from 1850 - 2014.  This 
#    script will complete all the changes required at year 1901 to deal with the
#    fact that met forcing data does not go back to 1850.
#
#  - Unmodified script will do the following. 
#      Part 1: Simulation 1: 1850 - 1870 (21 years) using repeated 1901-1920 forcing
#      Part 2: Simulation 2: 1871 - 1900 (30 years) using repeated 1901-1920 forcing
#      Part 3: Simulation 3+4+5+6: 1901-1988 (four 22 year) simulations using 1901-1988 forcing
#      Part 4: Simulation 7: 1989-2004 (one 16 year) branch simulation w/daily output using 1989-2004 forcing
#      Part 5: Simulation 8: 2005-2014 (one 10 year) branch simulation w/daily & subdaily output using 2005-2014 forcing
#
#  - Script assumes that simulation can run at least 30 years within a 12 hour block on 
#    Cheyenne.  To find the timing in an equivalent sample run, look in the timing 
#    directory and grep as follows  > grep 'simulated_years/day' cesm_timing*
#
#  - In the env_batch.xml file for the case.run group ensure the following: <entry id="JOB_WALLCLOCK_TIME" value="12:00:00">
#    that way you can use up to 12:00 hours of wall-clock computer time per run block
#
#  - This script assumes that env_mach_pes.xml has been setup and case.setup has already been run
#
#  - This script makes use of user_nl_datm1901-1920 and user_nl_datm1901-2014
#
#  - Before submitting script, make a copy of your modified or unmodified user_nl_clm file
#    into "original_user_nl_clm".  This should only contain namelist items that will not change throughout
#    the run.
#    Create a file called user_nl_clm_histdaily that contains the desired history output namelist items
#    for the 1989-2004 simulation
#    Create a file called user_nl_clm_histsubdaily that contains the desired history output namelist items
#    for the 2005-2014 simulation
#
#  - The atm data files start in 1901, so with :
#    ALIGN year of 1901, (this is in units of RUN or simulation years)
#    START year of 1901, (this is in units of FORCE'ing data years)
#
#    RUN Year   : 1850 ... 1860 1861 ... 1870 ... 1880 1881 ... 1890 ... 1900 1901 ... 2014
#    FORCE Year : 1910 ... 1920 1901 ... 1910 ... 1920 1901 ... 1910 ... 1920 1901 ... 2014
#
#  - The script could be broken up into several parts if you want to check the initial set of 
#    simulations.
#
# - Written by Andrew Slater - Late July, 2015; aslater@kryos.colorado.edu
# - Modified by Dave Lawrence, August, 2015
# - Updated with better check that run has also been archived to short term - Dave Lawrence
#      October, 2015
# - Updated to adjust for the fact that the model is now slower - Keith Oleson December, 2016
#      ./run_clm_historical.v5.csh ! > & run_historical.out &
# - Extend to 2014 and obtain daily & sub-daily output for end of run - Keith Oleson January, 2018
#      ./run_clm_historical.v6.csh ! > & run_historical.out &
# - Modify history output for CMIP6 - Keith Oleson January, 2019
#      ./run_clm_historical.v7.csh ! > & run_historical.out &
#########################################################################################

#########################################################################################
# PART 1
#########################################################################################
#
#    This portion does the initial setup and the initial 21 year run (1850-1870)
#
#########################################################################################

# --- CASENAME is your case name
set CASENAME = 'clm50_release-clm5.0.15_2deg_GSWP3V1_hist'

# --- Set the user namelist file.
cp original_user_nl_clm user_nl_clm

# --- Ensure that the env_run.xml file has the correct content
# Since this particular run won't go for 51 years in 12 hours (as in v4 script), set it up for 21 years first
./xmlchange RUN_TYPE=startup
./xmlchange RUN_STARTDATE=1850-01-01
./xmlchange STOP_OPTION=nyears
./xmlchange STOP_N=21
./xmlchange CONTINUE_RUN=FALSE
./xmlchange RESUBMIT=0
./xmlchange DATM_CLMNCEP_YR_ALIGN=1901
./xmlchange DATM_CLMNCEP_YR_START=1901
./xmlchange DATM_CLMNCEP_YR_END=1920

# need to use user_nl_datm files to get years right
cp user_nl_datm1901-1920 user_nl_datm

# --- Check that you end up using the correct env_run.xml file
set nenvr = `ls -1 env_run*.xml | wc -l`
if ($nenvr > 1) then
  echo 'There is more than one file of the type env_run*.xml'
  echo 'There should only be one such file'
  exit
endif

# --- If you have not already built the code, then do so now
#./case.clean_build
qcmd -- ./case.build

# --- Now submit the job and let it run
./case.submit

#########################################################################################
# PART 2
#########################################################################################
#
#    This portion checks to see if the 1850-1870 portion of the run is done (or it waits
#    10 minutes before checking again). 
#
#    This will then start the run from 1871 (hence the CONTINUE_RUN=TRUE) and do one 30 year 
#    simulation
#
#    The new values for env_run.xml are put in place
#    Then submit the job
#
#########################################################################################


set WDIR = '/glade/scratch/'$USER'/'$CASENAME'/run/'
set DONE_RUNA = 0
set DONE_ARCHIVE = 0
set RESTART_FILE = $WDIR$CASENAME'.clm2.r.1871-01-01-00000.nc'

# --- Check if the first set of simulations have completed and the data archived (every ten minutes)
while ($DONE_RUNA == 0)
  if (-e $RESTART_FILE) then
    set DONE_RUNA = 1
    echo '1850-1870 run is complete'
    while ($DONE_ARCHIVE == 0) 
       set nh0 = `ls -l $WDIR/*clm2.h0.* | egrep -c '^-'`
       echo $nh0
       if ($nh0 == 1) then
          set DONE_ARCHIVE = 1 
          echo 'Files have been archived'
       else 
          sleep 600
          date
       endif
  else
    sleep 600
    date
  endif
end

# These are the proper settings to let this script continue the run through 1900
./xmlchange STOP_N=30
./xmlchange CONTINUE_RUN=TRUE

# --- Now submit the job and let it run
./case.submit

#########################################################################################
# PART 3
#########################################################################################
#
#    This portion checks to see if the 1871-1900 portion of the run is done (or it waits
#    10 minutes before checking again). It then removes (or rather moves and renames) the 
#    datm files so that the model will use the full array of data from 1901-2014.
#    This part runs with forcing data files that actually exist for 1901-2014
#
#    This will start the run from 1901 (hence the CONTINUE_RUN=TRUE) and do four 22 year 
#    simulations: 1901 + 4*22 - 1 = 1988  (minus 1 because we do 1901)
#
#    The new values for env_run.xml are put in place
#    Then submit the job
#
#########################################################################################


set WDIR = '/glade/scratch/'$USER'/'$CASENAME'/run/'
set DDIR = $WDIR'restart_dump/'
set DONE_RUNA = 0
set DONE_ARCHIVE = 0
set RESTART_FILE = $WDIR$CASENAME'.clm2.r.1901-01-01-00000.nc'

# --- Check if the first set of simulations have completed and the data archived (every ten minutes)
while ($DONE_RUNA == 0)
  if (-e $RESTART_FILE) then
    set DONE_RUNA = 1
    echo '1850-1900 run is complete'
    while ($DONE_ARCHIVE == 0) 
       set nh0 = `ls -l $WDIR/*clm2.h0.* | egrep -c '^-'`
       echo $nh0
       if ($nh0 == 1) then
          set DONE_ARCHIVE = 1 
          echo 'Files have been archived'
       else 
          sleep 600
          date
       endif
  else
    sleep 600
    date
  endif
end

# --- If the first two sets of simulations are done, move the datm files and compress them
if (! -d $DDIR) then
  mkdir $DDIR
endif
mv -i $WDIR$CASENAME.datm.rs1*.bin $DDIR
gzip $DDIR$CASENAME*.bin

# Since this particular run won't go for 55 years in 12 hours, do this in four 22 year chunks, thus
# we have to resubmit the job 3 times.
./xmlchange STOP_OPTION=nyears
./xmlchange STOP_N=22
./xmlchange DATM_CLMNCEP_YR_ALIGN=1901
./xmlchange DATM_CLMNCEP_YR_START=1901
./xmlchange DATM_CLMNCEP_YR_END=2014
./xmlchange CONTINUE_RUN=TRUE
./xmlchange RESUBMIT=3

# need to use user_nl_datm files to get years right
cp user_nl_datm1901-2014 user_nl_datm

# --- Check that you end up using the correct env_run.xml file
set nenvr = `ls -1 env_run*.xml | wc -l`
if ($nenvr > 1) then
  echo 'There is more than one file of the type env_run*.xml'
  echo 'There should only be one such file'
  exit
endif

# --- Now submit the job and let it run
./case.submit

#########################################################################################
# PART 4
#########################################################################################
#
#    This portion checks to see if the 1901-1988 part of the run is complete
#    and then runs the model for 1989-2004 as a branch run to get daily output
#
#########################################################################################

set DONE_RUNA = 0
set DONE_ARCHIVE = 0
set RESTART_FILE = $WDIR$CASENAME'.clm2.r.1989-01-01-00000.nc'

# --- Check if the second set of simulations have completed and the data archived (every ten minutes)
while ($DONE_RUNA == 0)
  if (-e $RESTART_FILE) then
    set DONE_RUNA = 1
    echo '1901-1989 run is complete'
    while ($DONE_ARCHIVE == 0)
       set nh0 = `ls -l $WDIR/*clm2.h0.* | egrep -c '^-'`
       echo $nh0
       if ($nh0 == 1) then
          set DONE_ARCHIVE = 1
          echo 'Files have been archived'
       else
          sleep 600
          date
       endif
  else
    sleep 600
    date
  endif
end

# --- Ensure that the env_run.xml file has the correct content
./xmlchange RUN_TYPE=branch
./xmlchange RUN_REFCASE={$CASENAME}
./xmlchange RUN_REFDATE=1989-01-01
./xmlchange STOP_OPTION=nyears
./xmlchange STOP_N=16
./xmlchange CONTINUE_RUN=FALSE
./xmlchange RESUBMIT=0

# --- Add in the daily output streams
# --- Reset the user namelist file.
cp original_user_nl_clm user_nl_clm

# --- Add in the daily history output items
cat user_nl_clm_histdaily >> user_nl_clm

# --- Now submit the job and let it run
./case.submit

#########################################################################################
# PART 5
#########################################################################################
#
#    This portion checks to see if the 1989-2004 part of the run is complete
#    and then runs the model for 2005-2014 as a branch run to get daily and subdaily output
#
#########################################################################################

set DONE_RUNA = 0
set DONE_ARCHIVE = 0
set RESTART_FILE = $WDIR$CASENAME'.clm2.r.2005-01-01-00000.nc'

# --- Check if the second set of simulations have completed and the data archived (every ten minutes)
while ($DONE_RUNA == 0)
  if (-e $RESTART_FILE) then
    set DONE_RUNA = 1
    echo '1989-2004 run is complete'
    while ($DONE_ARCHIVE == 0)
       set nh0 = `ls -l $WDIR/*clm2.h0.* | egrep -c '^-'`
       echo $nh0
       if ($nh0 == 1) then
          set DONE_ARCHIVE = 1
          echo 'Files have been archived'
       else
          sleep 600
          date
       endif
  else
    sleep 600
    date
  endif
end

# --- Ensure that the env_run.xml file has the correct content
./xmlchange RUN_TYPE=branch
./xmlchange RUN_REFCASE={$CASENAME}
./xmlchange RUN_REFDATE=2005-01-01
./xmlchange STOP_OPTION=nyears
./xmlchange STOP_N=10
./xmlchange CONTINUE_RUN=FALSE
./xmlchange RESUBMIT=0

# --- Add in the subdaily output streams
# --- Reset the user namelist file.
cp original_user_nl_clm user_nl_clm

# --- Add in the daily history output items
cat user_nl_clm_histdaily >> user_nl_clm

# --- Add in the subdaily history output items
cat user_nl_clm_histsubdaily >> user_nl_clm

# --- Now submit the job and let it run
./case.submit

#########################################################################################
#
#    This portion checks to see if the 2005-2014 part of the run is complete
#    and ends the script
#
#########################################################################################

set DONE_RUNA = 0
set DONE_ARCHIVE = 0
set RESTART_FILE = $WDIR$CASENAME'.clm2.r.2015-01-01-00000.nc'

# --- Check if the second set of simulations have completed and the data archived (every ten minutes)
while ($DONE_RUNA == 0)
  if (-e $RESTART_FILE) then
    set DONE_RUNA = 1
    echo '2005-2014 run is complete'
    while ($DONE_ARCHIVE == 0)
       set nh0 = `ls -l $WDIR/*clm2.h0.* | egrep -c '^-'`
       echo $nh0
       if ($nh0 == 1) then
          set DONE_ARCHIVE = 1
          echo 'Files have been archived'
       else
          sleep 600
          date
       endif
  else
    sleep 600
    date
  endif
end
