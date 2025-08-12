#!/bin/csh

#########################################################################################
#
#  - Execute this script to do a CLM historical simulation from 1850 - 2023.  This 
#    script will complete all the changes required at year 1901 to deal with the
#    fact that met forcing data does not go back to 1850.
#
#  - The CASENAME should be the only thing the user needs to change.
#
#  - Run this script in the background on a Derecho login node like this:
#    ./run_clm_historical.v11.csh >&! run_clm_historical.out &
#    CAUTION: Note the hostname (echo $HOST) and the PID for the script and the "sleep" job (identify them using ps -u $USER)
#    E.g., 
#      PID TTY          TIME CMD
#      581 ?        00:00:00 sleep
#    30684 ?        00:00:00 run_clm_histori
#    If the simulation crashes at any point, you'll need to kill the script and the "sleep" job (kill -9 $PID) 
#    before restarting the simulation. You will have to run the remainder of the simulation manually.
#
#  - Unmodified script will do the following. 
#      Part 1: Simulation 1: 1850 - 1870 (21 years) using repeated 1901-1920 forcing
#      Part 2: Simulation 2: 1871 - 1900 (30 years) using repeated 1901-1920 forcing
#      Part 3: Simulation 3+4+5: 1901-1978 (three 26 year) simulations using 1901-1978 forcing
#      Part 4: Simulation 7: 1979-1999 (one 21 year) branch simulation w/daily output using 1979-1999 forcing
#      Part 5: Simulation 8: 2000-2023 (one 24 year) branch simulation w/daily & subdaily output using 2000-2023 forcing
#      Note that this approach generates year 1979 and 2000 restart files needed for other specific compsets.
#
#  - Script assumes that simulation can run at least 30 years within a 12 hour block on 
#    Derecho.  To find the timing in an equivalent sample run, look in the timing 
#    directory and grep as follows  > grep 'simulated_years/day' cesm_timing*
#
#  - In the env_batch.xml file for the case.run group ensure the following: <entry id="JOB_WALLCLOCK_TIME" value="12:00:00">
#    that way you can use up to 12:00 hours of wall-clock computer time per run block
#
#  - This script assumes that env_mach_pes.xml has been setup and case.setup has already been run
#
#  - Before submitting script, make a copy of your modified or unmodified user_nl_clm file
#    into "original_user_nl_clm".  This should only contain namelist items that will not change throughout
#    the run. You can start with this example:
#    /glade/u/home/oleson/run_hist_1850_files/BGC/original_user_nl_clm_ctsm5.3.0
#    Create a file called user_nl_clm_histdaily that contains the desired history output namelist items
#    for the 1979-2023 simulation. You can start with this example:
#    /glade/u/home/oleson/run_hist_1850_files/BGC/user_nl_clm_histdaily_ctsm
#    Create a file called user_nl_clm_histsubdaily that contains the desired history output namelist items
#    for the 2001-2023 simulation. You can start with this example:
#    /glade/u/home/oleson/run_hist_1850_files/BGC/user_nl_clm_histsubdaily_ctsm
#
#  - The atm data files start in 1901, so with :
#    ALIGN year of 1901, (this is in units of RUN or simulation years)
#    START year of 1901, (this is in units of FORCE'ing data years)
#
#    RUN Year   : 1850 ... 1860 1861 ... 1870 ... 1880 1881 ... 1890 ... 1900 1901 ... 2023
#    FORCE Year : 1910 ... 1920 1901 ... 1910 ... 1920 1901 ... 1910 ... 1920 1901 ... 2023
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
# - Modify to run with nuopc - Keith Oleson September, 2022
#      ./run_clm_historical.v9.csh ! > & run_historical.out &
# - Modify to produce restart files for 1979 and 2000 and run with CRUv7 extension to 2023 - Keith Oleson September, 2024
#      ./run_clm_historical.v10.csh ! > & run_clm_historical.out &
# - Modify to look for either h0 or h0a files - Keith Oleson July, 2025
#      ./run_clm_historical.v11.csh ! > & run_clm_historical.out &
#########################################################################################

#########################################################################################
# PART 1
#########################################################################################
#
#    This portion does the initial setup and the initial 21 year run (1850-1870)
#
#########################################################################################

# --- CASENAME is your case name
set CASENAME = 'ctsm530_f19_PPE_TESTv11_hist'

# --- Set to either h0 (ctsm5.3.061 or earlier) or h0a (ctsm5.3.062 or later)
set HIST_EXT = 'h0'

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
./xmlchange DATM_YR_ALIGN=1901
./xmlchange DATM_YR_START=1901
./xmlchange DATM_YR_END=1920
./xmlchange DATM_SKIP_RESTART_READ=FALSE

# --- Check that you end up using the correct env_run.xml file
set nenvr = `ls -1 env_run*.xml | wc -l`
if ($nenvr > 1) then
  echo 'There is more than one file of the type env_run*.xml'
  echo 'There should only be one such file'
  exit
endif

# --- If you have not already built the code, then do so now
#./case.build --clean-all
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


set WDIR = '/glade/derecho/scratch/'$USER'/'$CASENAME'/run/'
set DONE_RUNA = 0
set DONE_ARCHIVE = 0
set RESTART_FILE = $WDIR$CASENAME'.clm2.r.1871-01-01-00000.nc'

# --- Check if the first set of simulations have completed and the data archived (every ten minutes)
while ($DONE_RUNA == 0)
  if (-e $RESTART_FILE) then
    set DONE_RUNA = 1
    echo '1850-1870 run is complete'
    while ($DONE_ARCHIVE == 0) 
       set nh0 = `ls -l $WDIR/*clm2.{$HIST_EXT}.* | egrep -c '^-'`
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
#    10 minutes before checking again). It then sets DATM_SKIP_RESTART_READ to TRUE
#    so that the model will use the full array of forcing data from 1901-2023.
#
#    This will start the run from 1901 (hence the CONTINUE_RUN=TRUE) and do three 26 year 
#    simulations: 1901 + 3*26 - 1 = 1978  (minus 1 because we do 1901)
#    This will give us a restart file at the beginning of 1979.
#
#    The new values for env_run.xml are put in place
#    Then submit the job
#
#########################################################################################


set WDIR = '/glade/derecho/scratch/'$USER'/'$CASENAME'/run/'
set DONE_RUNA = 0
set DONE_ARCHIVE = 0
set RESTART_FILE = $WDIR$CASENAME'.clm2.r.1901-01-01-00000.nc'

# --- Check if the first set of simulations have completed and the data archived (every ten minutes)
while ($DONE_RUNA == 0)
  if (-e $RESTART_FILE) then
    set DONE_RUNA = 1
    echo '1850-1900 run is complete'
    while ($DONE_ARCHIVE == 0) 
       set nh0 = `ls -l $WDIR/*clm2.{$HIST_EXT}.* | egrep -c '^-'`
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

# Since this particular run won't go for 78 years in 12 hours, do this in three 26 year chunks, thus
# we have to resubmit the job 2 times.
./xmlchange DATM_SKIP_RESTART_READ=TRUE
./xmlchange STOP_OPTION=nyears
./xmlchange STOP_N=26
./xmlchange DATM_YR_ALIGN=1901
./xmlchange DATM_YR_START=1901
./xmlchange DATM_YR_END=2023
./xmlchange CONTINUE_RUN=TRUE
./xmlchange RESUBMIT=2

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
#    This portion checks to see if the 1901-1978 part of the run is complete
#    and then runs the model for 1979-1999 as a branch run to get daily output
#    and to get a restart file at the beginning of 2000
#
#########################################################################################

set DONE_RUNA = 0
set DONE_ARCHIVE = 0
set RESTART_FILE = $WDIR$CASENAME'.clm2.r.1979-01-01-00000.nc'

# --- Check if the second set of simulations have completed and the data archived (every ten minutes)
while ($DONE_RUNA == 0)
  if (-e $RESTART_FILE) then
    set DONE_RUNA = 1
    echo '1901-1978 run is complete'
    while ($DONE_ARCHIVE == 0)
       set nh0 = `ls -l $WDIR/*clm2.{$HIST_EXT}.* | egrep -c '^-'`
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
./xmlchange RUN_REFDATE=1979-01-01
./xmlchange STOP_OPTION=nyears
./xmlchange STOP_N=21
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
#    This portion checks to see if the 1979-1999 part of the run is complete
#    and then runs the model for 2000-2023 as a branch run to get daily and subdaily output
#
#########################################################################################

set DONE_RUNA = 0
set DONE_ARCHIVE = 0
set RESTART_FILE = $WDIR$CASENAME'.clm2.r.2000-01-01-00000.nc'

# --- Check if the second set of simulations have completed and the data archived (every ten minutes)
while ($DONE_RUNA == 0)
  if (-e $RESTART_FILE) then
    set DONE_RUNA = 1
    echo '1979-1999 run is complete'
    while ($DONE_ARCHIVE == 0)
       set nh0 = `ls -l $WDIR/*clm2.{$HIST_EXT}.* | egrep -c '^-'`
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
./xmlchange RUN_REFDATE=2000-01-01
./xmlchange STOP_OPTION=nyears
./xmlchange STOP_N=24
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
#    This portion checks to see if the 2005-2023 part of the run is complete
#    and ends the script
#
#########################################################################################

set DONE_RUNA = 0
set DONE_ARCHIVE = 0
set RESTART_FILE = $WDIR$CASENAME'.clm2.r.2024-01-01-00000.nc'

# --- Check if the second set of simulations have completed and the data archived (every ten minutes)
while ($DONE_RUNA == 0)
  if (-e $RESTART_FILE) then
    set DONE_RUNA = 1
    echo '2005-2023 run is complete'
    while ($DONE_ARCHIVE == 0)
       set nh0 = `ls -l $WDIR/*clm2.{$HIST_EXT}.* | egrep -c '^-'`
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
