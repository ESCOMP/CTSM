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
python $TDIR"/spinup_stability.py" ND.yml
status=$?
echo "status: "$status


# proceed accordingly
if [[ "$status" == "11" ]]; then
    echo "needs more spinup"
    cd $CASE_ND
    ./xmlchange CONTINUE_RUN=True
    ./xmlchange STOP_N=20
    ./xmlchange JOB_WALLCLOCK_TIME=4:00:00 --subgroup case.run
    cd $WDIR
    echo $CASE_ND>case.txt
    echo "./stabND.sh spinup.config">commands.txt
elif [[ "$status" == "0" ]]; then
    echo "spinup appears sufficient"
    rm commands.txt  # sequence complete, exit tether
else
    echo "something looks wrong, halting tether"
    echo $TDIR"spinup_stability.py ND.yml experienced an error"
    kill_tether
fi
