#!/bin/bash

WDIR="/glade/u/home/djk2120/vp/sims/ICLM50Bgc.CPLHIST.default/"
NAMELISTS="/glade/u/home/djk2120/vp/scripts/namelists/"
CASENAME="I1850Clm50Bgc.CPLHIST.default.AD"
PROJECT=P93300041

COMPSET=1850_DATM%CPLHIST_CLM50%BGC-CROP_SICE_SOCN_MOSART_CISM2%NOEVOLVE_SWAV
GRID=f09_f09_mg17
CESMROOT="/glade/work/djk2120/cesm2.1.5/"

curdir=$(pwd)
CASEROOTBASE=$WDIR
caseroot=$CASEROOTBASE$CASENAME
echo $caseroot

cd $CESMROOT/cime/scripts
./create_newcase --case $caseroot --compset $COMPSET --res $GRID --project $PROJECT --mach derecho --run-unsupported
cd $caseroot
./case.setup

./xmlchange RUN_TYPE=hybrid
./xmlchange PROJECT=$PROJECT
./xmlchange JOB_PRIORITY="premium"
./xmlchange RUN_STARTDATE="0001-01-01"
./xmlchange RUN_REFCASE="b.e21.BHISTcmip6.f09_g17.LE2-1001.001"
./xmlchange RUN_REFDIR="/glade/campaign/cesm/collections/CESM2-LE/restarts/b.e21.BHISTcmip6.f09_g17.LE2-1001.001/rest/1880-01-01-00000"
./xmlchange RUN_REFDATE="1880-01-01"
./xmlchange GET_REFCASE="True"
./xmlchange DATM_MODE="CPLHIST"
./xmlchange DATM_PRESAERO="cplhist"
./xmlchange DATM_TOPO="cplhist"
./xmlchange DATM_CPLHIST_CASE="f.e21.FHIST_BGC.f09_f09.ersstv5.cplhist"
./xmlchange DATM_CPLHIST_DIR="/glade/derecho/scratch/djk2120/archive/f.e21.FHIST_BGC.f09_f09.ersstv5.cplhist/cpl/proc/"
./xmlchange DATM_CPLHIST_YR_ALIGN="1880"
./xmlchange DATM_CPLHIST_YR_START="1880"
./xmlchange DATM_CPLHIST_YR_END="1899"
./xmlchange STOP_N=100,STOP_OPTION=nyears
./xmlchange CLM_ACCELERATED_SPINUP=on


cp $NAMELISTS"AD/"* ./
./case.build

cd $WDIR
echo $CASENAME>case.txt
echo "./stabAD.sh">commands.txt
