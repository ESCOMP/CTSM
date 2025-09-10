#!/bin/bash

WDIR="/glade/u/home/djk2120/vp/sims/ICLM50Bgc.CPLHIST.default/"
NAMELISTS="/glade/u/home/djk2120/vp/scripts/namelists/"
CLONENAME="I1850Clm50Bgc.CPLHIST.default.AD"
CASENAME="I1850Clm50Bgc.CPLHIST.default.pAD"
PROJECT=P93300041
ARCHIVE="/glade/derecho/scratch/djk2120/archive/"
CESMROOT="/glade/work/djk2120/cesm2.1.5/"

curdir=$(pwd)
CASEROOTBASE=$WDIR
cloneroot=$CASEROOTBASE$CLONENAME
caseroot=$CASEROOTBASE$CASENAME
echo $caseroot


cd $CESMROOT/cime/scripts
./create_clone --case $caseroot --clone $cloneroot
cd $caseroot
./case.setup


./xmlchange RUN_TYPE=hybrid
./xmlchange CONTINUE_RUN=False
./xmlchange PROJECT=$PROJECT
./xmlchange JOB_PRIORITY="regular"
./xmlchange JOB_WALLCLOCK_TIME="12:00:00"
./xmlchange RUN_STARTDATE="0001-01-01"
./xmlchange STOP_N=100,STOP_OPTION=nyears
./xmlchange CLM_ACCELERATED_SPINUP=off

# finding the latest restart from AD
# this code is likely very brittle
./xmlchange RUN_REFCASE=$CLONENAME
./xmlchange GET_REFCASE="True"
REFREST=$ARCHIVE$CLONENAME"/rest"
last_date=$(ls $REFREST | tail -n1)
REFDIR=$REFREST"/"$last_date
REFDATE=${last_date%-*} 
./xmlchange RUN_REFDIR=$REFDIR
./xmlchange RUN_REFDATE=$REFDATE


cp $NAMELISTS"pAD/"* ./
./case.build

cd $WDIR
echo $CASENAME>case.txt
echo "./stabPAD.sh">commands.txt





