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
caseroot=$WDIR"/"$CASE_AD
cd $CDIR
./create_newcase --case $caseroot --compset $PI_COMPSET --res $GRID --project $PROJECT --mach derecho --run-unsupported


# setup case
cd $caseroot
./case.setup


# some xml changes
./xmlchange JOB_PRIORITY=premium
./xmlchange STOP_OPTION=nyears
./xmlchange STOP_N=80
./xmlchange CLM_ACCELERATED_SPINUP=on
./xmlchange MOSART_MODE=NULL

# cp in namelist
ad_namelists=$NLDIR"/AD"
if [ -d $ad_namelists ]; then
    cd $ad_namelists
    if [ -f user_nl_clm ]; then
	cp user_nl_* $caseroot
    else
	echo "ERROR: no user_nl_clm in "$ad_namelists
	kill_tether
    fi
else
    echo "ERROR: "$ad_namelists" is not a directory"
    kill_tether
fi


# build case
#   do not submit! tether.sh will submit
cd $caseroot
./case.build


# tether commands
cd $WDIR
echo "./stabAD.sh spinup.config" > commands.txt
echo $CASE_AD > case.txt
