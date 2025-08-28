
case="I1850.f45_g37.txd.part1"
wdir=$(pwd)

#configuing the new case
CIME="/glade/work/djk2120/cesm2.1.5/cime/scripts" #my checkout feel free to replace with yours
COMPSET=I1850Clm50Sp
GRID=f45_g37
PROJECT=P93300041

#creating the case
$CIME/create_newcase --case $case --compset $COMPSET --res $GRID --project $PROJECT --mach derecho --run-unsupported

#working on the case a little bit
cd $case
./case.setup
./xmlchange STOP_N=1
./xmlchange STOP_OPTION=nmonths
./xmlchange JOB_PRIORITY=premium
./xmlchange JOB_WALLCLOCK_TIME="01:00:00"

#build the case (DO NOT SUBMIT!!!)
./case.build

#tether will submit the case in case.txt
cd $wdir
echo $case > case.txt

#tell tether what to do after this case
echo "bash part2.sh" > commands.txt
