#!/bin/csh

#########################################################################################
#
#  - Execute this script to do a CLM historical simulation from 1850 - 2010.  This 
#    script will complete all the changes required at year 1901 to deal with the
#    fact that met forcing data does not go back to 1850.
#
#  - Unmodified script will do the following. 
#      Simulation 1: 1850 - 1901 (51 years) using repeated 1901-1920 forcing
#      Simulation 2+3: 1901-2010 (two 55 year) simulations
#
#  - Script assumes that can simulation > 60 years within a 12 hour block on 
#    Yellowstone.  To find the timing in an equivalent sample run, look in the timing 
#    directory and grep as follows  > grep 'simulated_years/day' cesm_timing*
#
#  - In the $CASENAME.run file ensure the following:  #BSUB -W 12:00
#    that way you can use up to 12:00 hours of wall-clock computer time per run block
#
#  - Before submitting script, make a copy of your modified or unmodified user_nl_clm file i
#    into "original_user_nl_clm".  If you wish to make additional changes to user_nl_clm
#    via this script, you can do that by editing thee Set the user namelist file section below.
#
#  - The data files start in 1901 for CRUNCEP, so with :
#    ALIGN year of 1901, (this is in units of RUN or simulation years)
#    START year of 1901, (this is in units of FORCE'ing data years)
#
#    RUN Year   : 1850 ... 1860 1861 ... 1870 ... 1880 1881 ... 1890 ... 1899 1900 1901
#    FORCE Year : 1910 ... 1920 1901 ... 1910 ... 1920 1901 ... 1910 ... 1919 1920 1901
#
#  - The script could be broken up into two parts if you want to check the initial set of 
#    simulations.
#
# - Written by Andrew Slater - Late July, 2015; aslater@kryos.colorado.edu
# - Modified by Dave Lawrence, August, 2015
# - Updated with better check that run has also been archived to short term - Dave Lawrence
#      October, 2015
# - ./run_clm_historical.v3.csh ! > & run_historical.out &
#########################################################################################

# --- CASENAME is your case name
set CASENAME = 'clm5_respn05r162_2degGSWP3_hist_allN_v01'

# --- Set the user namelist file. i.e. add the snow alterations for my run. Others may have 
#      diferent needs
#set UNLCLM = 'script_user_nl_clm'
#cat original_user_nl_clm > $UNLCLM

# Examples of how to modify user_nl_clm using this script.  Uncomment and edit as desired.
#echo "hist_empty_htapes = .true." >> $UNLCLM
#echo "hist_fincl1 = 'TSA'" >> $UNLCLM
#cat $UNLCLM > user_nl_clm


# --- Check that you end up using the correct env_run.xml file
set nenvr = `ls -1 env_run*.xml | wc -l`
if ($nenvr > 1) then
  echo 'There is more than one file of the type env_run*.xml'
  echo 'There should only be one such file'
  exit
endif


# --- Ensure that the env_run.xml file has the correct content
./xmlchange RUN_TYPE=startup
./xmlchange RUN_STARTDATE=1850-01-01
./xmlchange STOP_OPTION=nyears
./xmlchange STOP_N=51
./xmlchange DATM_CLMNCEP_YR_ALIGN=1901
./xmlchange DATM_CLMNCEP_YR_START=1901
./xmlchange DATM_CLMNCEP_YR_END=1920
./xmlchange CONTINUE_RUN=FALSE
./xmlchange RESUBMIT=0

# --- If you have not already built the code, then do so now
#./case.clean_build
#./case.build

# --- Now submit the job and let it run
./case.submit


#########################################################################################
# PART 2
#########################################################################################
#
#    This portion checks to see if the 1850-1901 portion of the run is done (or it waits
#    10 minutes before checking again). It then removes (or rather moves and renames) the 
#    datm files so that the model will use the full array of data from 1901-2010.
#    This part runs with forcing data files that actually exist
#
#    This will start the run from 1901 (hence the CONTINUE=TRUE) and do two 55 year 
#    simulations: 1901 + 2*55 - 1 = 2010  (minus 1 becase we do 1901)
#
#    The new values for env_run.xml are put in place
#    The values for user_nl_clm are put in place
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
          sleep 300
          date
       endif
  else
    sleep 300
    date
  endif
end

# --- If the first set of simulations are done,  move the datm files and compress them
if (! -d $DDIR) then
  mkdir $DDIR
endif
mv -i $WDIR$CASENAME.datm.rs1*.bin $DDIR
gzip $DDIR$CASENAME*.bin


# --- Note that the name list file (user_nl_clm) has not changed from before


# --- Check that you end up using the correct env_run.xml file
set nenvr = `ls -1 env_run*.xml | wc -l`
if ($nenvr > 1) then
  echo 'There is more than one file of the type env_run*.xml'
  echo 'There should only be one such file'
  exit
endif


# --- Ensure that the env_run.xml file has the correct content
#     We want to do a 110 year simulation, but in two 55 year chunks, thus
#     we have to resubmit the job once.
./xmlchange STOP_OPTION=nyears
./xmlchange STOP_N=55
./xmlchange DATM_CLMNCEP_YR_ALIGN=1901
./xmlchange DATM_CLMNCEP_YR_START=1901
./xmlchange DATM_CLMNCEP_YR_END=2010
./xmlchange CONTINUE_RUN=TRUE
./xmlchange RESUBMIT=1

./case.submit

