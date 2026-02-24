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


# run the spinup stability script
python $SDIR"/spinup_stability.py" SASU.yml
status=$?
echo "status: "$status


# proceed accordingly
if [[ "$status" == "11" ]]; then
    echo "needs more spinup"
    cd $CASE_SASU
    ./xmlchange CONTINUE_RUN=True
    ./xmlchange STOP_N=20
    ./xmlchange JOB_WALLCLOCK_TIME=4:00:00 --subgroup case.run
    cd $WDIR
    echo $CASE_SASU>case.txt
    echo "./stabSASU.sh spinup.config">commands.txt
elif [[ "$status" == "0" ]]; then
    echo "spinup appears sufficient"
    ./setupND.sh spinup.config
else
    echo "something looks wrong, halting tether"
    echo $TDIR"spinup_stability.py SASU.yml experienced an error"
    kill_tether
fi
