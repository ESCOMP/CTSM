

SDIR="/glade/work/djk2120/ctsm_tether/tools/tether/old_spinup/"
WDIR="/glade/u/home/djk2120/vp/sims/"
case=$(<case.txt)

python $SDIR"spinup_stability.py" PAD.yml
status=$?
echo "status: "$status

if [[ "$status" == "11" ]]; then
    echo "PAD needs more spinup"
    cd $case
    ./xmlchange CONTINUE_RUN=True
    ./xmlchange STOP_N=20
    ./xmlchange JOB_WALLCLOCK_TIME="4:00:00"
    cd $WDIR
    echo $case>case.txt
    echo "stabPAD.sh">commands.txt
elif [[ "$status" == "0" ]]; then
    echo "PAD spinup appears sufficient"
    ./setupHIST.sh
    rm commands.txt
else
    echo "something looks wrong, halting tether"
    rm commands.txt
fi
