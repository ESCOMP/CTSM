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
source <(grep = $config)


# create case
caseroot=$WDIR"/"$CASE_SASU
cd $CDIR
./create_newcase --case $caseroot --compset $PI_COMPSET --res $GRID --project $PROJECT --mach derecho --run-unsupported


# setup case
cd $caseroot
./case.setup


# finding the latest restart from AD
# this code is likely very brittle
ARCHIVE="/glade/derecho/scratch/"$USER"/archive"
REFCASE=$CASE_AD
REFREST=$ARCHIVE"/"$REFCASE"/rest"
last_date=$(ls $REFREST | tail -n1)
REFDIR=$REFREST"/"$last_date
REFDATE=${last_date%-*} 


# xml changes
./xmlchange JOB_PRIORITY=premium
./xmlchange STOP_OPTION=nyears
./xmlchange STOP_N=80
./xmlchange MOSART_MODE=NULL
./xmlchange RUN_TYPE=hybrid
./xmlchange CLM_ACCELERATED_SPINUP=sasu
./xmlchange RUN_REFCASE=$REFCASE
./xmlchange GET_REFCASE="True"
./xmlchange RUN_REFDIR=$REFDIR
./xmlchange RUN_REFDATE=$REFDATE


# cp in namelist
namelists=$NLDIR"/SASU"
if [ -d $namelists ]; then
    cd $namelists
    if [ -f user_nl_clm ]; then
	cp user_nl_* $caseroot
    else
	echo "ERROR: no user_nl_clm in "$namelists
	kill_tether
    fi
else
    echo "ERROR: "$namelists" is not a directory"
    kill_tether
fi


# build case
#   do not submit! tether.sh will submit
cd $caseroot
./case.build


# tether commands
cd $WDIR
echo "./stabSASU.sh spinup.config" > commands.txt
echo $CASE_SASU > case.txt
