#!/bin/bash

# function to abort tether when handling errors
kill_tether () {
    cd $WDIR
    rm commands.txt
    exit 1
}


# import a bunch of variables from config file
#    e.g. $WDIR $PI_COMPSET 
config=$1
source <(grep = $1)


# create case
caseroot=$WDIR"/"$CASE2
cd $CDIR
./create_newcase --case $caseroot --compset $PI_COMPSET --res $GRID --project $PROJECT --mach derecho --run-unsupported


# setup case
cd $caseroot
./case.setup


# finding the latest restart from previous case
# this code is likely very brittle
ARCHIVE="/glade/derecho/scratch/"$USER"/archive"
REFCASE=$CASE1
REFREST=$ARCHIVE"/"$REFCASE"/rest"
last_date=$(ls $REFREST | tail -n1)
REFDIR=$REFREST"/"$last_date
REFDATE=${last_date%-*} 


# some xml changes
./xmlchange RUN_TYPE=hybrid
./xmlchange STOP_N=1
./xmlchange STOP_OPTION=nmonths
./xmlchange JOB_PRIORITY=premium
./xmlchange JOB_WALLCLOCK_TIME="01:00:00"
./xmlchange RUN_REFCASE=$REFCASE
./xmlchange GET_REFCASE="True"
./xmlchange RUN_REFDIR=$REFDIR
./xmlchange RUN_REFDATE=$REFDATE


#build the case (DO NOT SUBMIT!!!)
./case.build


# tether commands
cd $WDIR
rm commands.txt
echo $CASE2 > case.txt

