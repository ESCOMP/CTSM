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
caseroot=$WDIR"/"$CASE1
cd $CDIR
./create_newcase --case $caseroot --compset $PI_COMPSET --res $GRID --project $PROJECT --mach derecho --run-unsupported


# setup case
cd $caseroot
./case.setup


# some xml changes
./xmlchange STOP_N=1
./xmlchange STOP_OPTION=nmonths
./xmlchange JOB_PRIORITY=premium
./xmlchange JOB_WALLCLOCK_TIME="01:00:00"

#build the case (DO NOT SUBMIT!!!)
./case.build


# tether commands
cd $WDIR
echo "./part2.sh example1.config" > commands.txt
echo $CASE1 > case.txt

