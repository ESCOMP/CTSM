#! /bin/bash

# Script to compile CLM with case-specific settings
# For standalone CLM or coupling with COSMO (for coupling set COMPILER=MY_COMPILER-oasis)
# Domain can be global or EURO-CORDEX (set DOMAIN=eur, requires domain and mapping files) 

set -e # failing commands will cause the shell script to exit


#==========================================
# Case settings
#==========================================

echo "*** Setting up case ***"

date=`date +'%Y%m%d-%H%M'` # get current date and time
startdate=`date +'%Y-%m-%d %H:%M:%S'`
COMPSET=I2000Clm50SpGs # for CCLM2 (use stub glacier component for regional domain!)
RES=hcru_hcru # hcru_hcru for CCLM2-0.44, f09_g17 to test glob (inputdata downloaded)
DOMAIN=eur # eur for CCLM2 (EURO-CORDEX), sa for South-America, glob otherwise 

CODE=clm5.0 # clm5.0 for official release, clm5.0_features for Ronny's version, CTSMdev for latest 
COMPILER=nvhpc # gnu for gcc, or nvhpc; setting to gnu-oasis or nvhpc-oasis will: (1) use different compiler config, (2) copy oasis source code to CASEDIR
DRIVER=mct # default is mct, using nuopc requires ESMF installation
EXP=cclm2_${date} # custom case name
CASENAME=$CODE.$COMPILER.$COMPSET.$RES.$DOMAIN.$EXP
MACH=pizdaint
QUEUE=normal # USER_REQUESTED_QUEUE, overrides default JOB_QUEUE
WALLTIME="01:00:00" # USER_REQUESTED_WALLTIME, overrides default JOB_WALLCLOCK_TIME
PROJ=$(basename "$(dirname "${PROJECT}")") # extract project name (sm61/sm62)
NTASKS=24
NSUBMIT=0 # partition into smaller chunks, excludes the first submission
let "NCORES = $NTASKS * 12"
STARTDATE="2004-01-01"
NYEARS=1

# Set directories
export CLMROOT=$PWD # CLM code base directory on $PROJECT where this script is located
export CCLM2ROOT=$CLMROOT/.. # CCLM2 code base directory on $PROJECT where CLM, OASIS and COSMO are located
export CASEDIR=$SCRATCH/CCLM2_cases/$CASENAME # case directory on scratch
export CESMDATAROOT=$SCRATCH/CCLM2_inputdata # inputdata directory on scratch (to reuse, includes downloads and preprocessed EURO-CORDEX files)
export CESMOUTPUTROOT=$SCRATCH/CCLM2_output/$CASENAME # output directory on scratch

# Log output (use "tee" to send output to both screen and $outfile)
logfile=$SCRATCH/CCLM2_logs/${CASENAME}_mylogfile.log
mkdir -p "$(dirname "$logfile")" && touch "$logfile" # create parent/child directories and logfile
cp $CLMROOT/$BASH_SOURCE $SCRATCH/CCLM2_logs/${CASENAME}_myjobscipt.sh # copy this script to logs
print_log() {
    output="$1"
    echo "${output}" | tee -a $logfile
}

print_log "*** Case at: ${CASEDIR} ***"
print_log "*** Case settings: compset ${COMPSET}, resolution ${RES}, domain ${DOMAIN}, compiler ${COMPILER} ***"
print_log "*** Logfile at: ${logfile} ***"

# Sync inputdata on scratch because scratch will be cleaned every month (change inputfiles on $PROJECT!)
print_log "*** Syncing inputdata on scratch  ***"

#rsync -rv --ignore-existing /project/$PROJ/shared/CCLM2_inputdata/ $CESMDATAROOT/ | tee -a $logfile
sbatch --account=$PROJ --export=ALL,PROJ=$PROJ transfer_clm_inputdata.sh # xfer job to prevent overflowing the loginnode


#==========================================
# Load modules and find spack_oasis
#==========================================

# Load modules: now done through $USER/.cime/config_machines.xml
# print_log "*** Loading modules ***"
# daint-gpu (although CLM will run on cpus)
# PrgEnv-xxx (also switch compiler version if needed)
# cray-mpich
# cray-netcdf-hdf5parallel
# cray-hdf5-parallel
# cray-parallel-netcdf

#module list | tee -a $logfile

# Find spack_oasis installation (used in .cime/config_compilers.xml)
if [[ $COMPILER =~ "oasis" ]]; then
    print_log "*** Finding spack_oasis ***"
    export OASIS_PATH=$(spack location -i oasis%gcc) # e.g. /project/sm61/psieber/spack-install/oasis/master/gcc/24obfvejulxnpfxiwatzmtcddx62pikc
    print_log "*** OASIS at: ${OASIS_PATH} ***"
fi

print_log "*** LD_LIBRARY_PATH: ${LD_LIBRARY_PATH} ***"


#==========================================
# Create case
#==========================================

print_log "*** Creating CASE: ${CASENAME} ***"

cd $CLMROOT/cime/scripts
./create_newcase --case $CASEDIR --compset $COMPSET --res $RES --mach $MACH --compiler $COMPILER --driver $DRIVER --project $PROJ --run-unsupported | tee -a $logfile


#==========================================
# Modify namelists
#==========================================

print_log "*** Modifying env_*.xml  ***"
cd $CASEDIR

# Set directory structure
./xmlchange RUNDIR="$CASEDIR/run" # by defaut, RUNDIR is $SCRATCH/$CASENAME/run
./xmlchange EXEROOT="$CASEDIR/bld"

# Change job settings (env_batch.xml or env_workflow.xml). Do this here to change for both case.run and case.st_archive
./xmlchange JOB_QUEUE=$QUEUE --force
./xmlchange JOB_WALLCLOCK_TIME=$WALLTIME

# Set run start/stop options and DATM forcing (env_run.xml)
./xmlchange RUN_TYPE=startup
./xmlchange RESUBMIT=$NSUBMIT
./xmlchange RUN_STARTDATE=$STARTDATE
./xmlchange STOP_OPTION=nyears,STOP_N=$NYEARS
./xmlchange NCPL_BASE_PERIOD="day",ATM_NCPL=48 # coupling freq default 30min = day,48
if [ $DRIVER = mct ] ; then
    ./xmlchange DATM_CLMNCEP_YR_START=2004,DATM_CLMNCEP_YR_END=2004,DATM_CLMNCEP_YR_ALIGN=2004 # in clm5.0 with mct
fi
if [ $DRIVER = nuopc ] ; then
    ./xmlchange DATM_YR_START=2004,DATM_YR_END=2004,DATM_YR_ALIGN=2004 # new variable names in CTSMdev; probably only with nuopc driver
fi

# Set the number of cores and nodes (env_mach_pes.xml)
./xmlchange COST_PES=$NCORES
./xmlchange NTASKS_CPL=-$NTASKS
./xmlchange NTASKS_ATM=-$NTASKS
./xmlchange NTASKS_OCN=-$NTASKS
./xmlchange NTASKS_WAV=-$NTASKS
./xmlchange NTASKS_GLC=-$NTASKS
./xmlchange NTASKS_ICE=-$NTASKS
./xmlchange NTASKS_ROF=-$NTASKS
./xmlchange NTASKS_LND=-$NTASKS 

# If parallel netcdf is used, PIO_VERSION="2" (have not gotten this to work!)
#./xmlchange PIO_VERSION="1" # 1 is default in clm5.0, 2 is default in CTSMdev

# Activate debug mode (env_build.xml)
#./xmlchange DEBUG=TRUE
#./xmlchange INFO_DBUG=2 # Change amount of output

# Additional options
#./xmlchange CLM_BLDNML_OPTS="-irrig .true." -append # switch on irrigation
#./xmlchange RTM_MODE="NULL" # switch off river routing
#./xmlchange CCSM_BGC=CO2A,CLM_CO2_TYPE=diagnostic,DATM_CO2_TSERIES=20tr # set transient CO2

#./xmlchange CLM_NAMELIST_OPTS="use_init_interp=.false. # Ronny sets interp to false, not sure about this

# Domain and mapping files for limited spatial extent
if [ $DOMAIN == eur ]; then
    ./xmlchange LND_DOMAIN_PATH="$CESMDATAROOT/CCLM2_EUR_inputdata/domain"
    ./xmlchange LND_DOMAIN_FILE="domain_EU-CORDEX_0.5_correctedlons.nc"
    ./xmlchange LND2ROF_FMAPNAME="$CESMDATAROOT/CCLM2_EUR_inputdata/mapping/map_360x720_nomask_to_0.5x0.5_nomask_aave_da_c130103.nc"
    ./xmlchange ROF2LND_FMAPNAME="$CESMDATAROOT/CCLM2_EUR_inputdata/mapping/map_0.5x0.5_nomask_to_360x720_nomask_aave_da_c120830.nc"
    ./xmlchange LND2GLC_FMAPNAME="$CESMDATAROOT/CCLM2_EUR_inputdata/mapping/map_360x720_TO_gland4km_aave.170429.nc"
    ./xmlchange LND2GLC_SMAPNAME="$CESMDATAROOT/CCLM2_EUR_inputdata/mapping/map_360x720_TO_gland4km_aave.170429.nc"
    ./xmlchange GLC2LND_FMAPNAME="$CESMDATAROOT/CCLM2_EUR_inputdata/mapping/map_gland4km_TO_360x720_aave.170429.nc"
    ./xmlchange GLC2LND_SMAPNAME="$CESMDATAROOT/CCLM2_EUR_inputdata/mapping/map_gland4km_TO_360x720_aave.170429.nc"
    ./xmlchange MOSART_MODE=NULL # turn off MOSART for the moment because it runs globally
fi
if [ $DOMAIN == sa ]; then
    ./xmlchange LND_DOMAIN_PATH="$CESMDATAROOT/CCLM2_SA_inputdata/domain"
    ./xmlchange LND_DOMAIN_FILE="domain.lnd.360x720_SA-CORDEX_cruncep.100429.nc"
    ./xmlchange LND2ROF_FMAPNAME="$CESMDATAROOT/CCLM2_EUR_inputdata/mapping/map_360x720_nomask_to_0.5x0.5_nomask_aave_da_c130103.nc"
    ./xmlchange ROF2LND_FMAPNAME="$CESMDATAROOT/CCLM2_EUR_inputdata/mapping/map_0.5x0.5_nomask_to_360x720_nomask_aave_da_c120830.nc"
    ./xmlchange LND2GLC_FMAPNAME="$CESMDATAROOT/CCLM2_EUR_inputdata/mapping/map_360x720_TO_gland4km_aave.170429.nc"
    ./xmlchange LND2GLC_SMAPNAME="$CESMDATAROOT/CCLM2_EUR_inputdata/mapping/map_360x720_TO_gland4km_aave.170429.nc"
    ./xmlchange GLC2LND_FMAPNAME="$CESMDATAROOT/CCLM2_EUR_inputdata/mapping/map_gland4km_TO_360x720_aave.170429.nc"
    ./xmlchange GLC2LND_SMAPNAME="$CESMDATAROOT/CCLM2_EUR_inputdata/mapping/map_gland4km_TO_360x720_aave.170429.nc"
    ./xmlchange MOSART_MODE=NULL # turn off MOSART for the moment because it runs globally
fi

# Still global
./xmlchange LND2ROF_FMAPNAME="$CESMDATAROOT/CCLM2_EUR_inputdata/mapping/map_360x720_nomask_to_0.5x0.5_nomask_aave_da_c130103.nc"
./xmlchange ROF2LND_FMAPNAME="$CESMDATAROOT/CCLM2_EUR_inputdata/mapping/map_0.5x0.5_nomask_to_360x720_nomask_aave_da_c120830.nc"

# Not needed for stub components (?)
#./xmlchange LND2GLC_FMAPNAME="$CESMDATAROOT/CCLM2_EUR_inputdata/mapping/map_360x720_TO_gland4km_aave.170429.nc"
#./xmlchange LND2GLC_SMAPNAME="$CESMDATAROOT/CCLM2_EUR_inputdata/mapping/map_360x720_TO_gland4km_aave.170429.nc"
#./xmlchange GLC2LND_FMAPNAME="$CESMDATAROOT/CCLM2_EUR_inputdata/mapping/map_gland4km_TO_360x720_aave.170429.nc"
#./xmlchange GLC2LND_SMAPNAME="$CESMDATAROOT/CCLM2_EUR_inputdata/mapping/map_gland4km_TO_360x720_aave.170429.nc"

# ESMF interface and time manager (env_build.xml)
#./xmlchange -file env_build.xml -id COMP_INTERFACE -val "mct" # mct is default in clm5.0, nuopc is default in CTSMdev (requires ESMF installation); adding --driver mct to create_newcase creates the case with everything needed
#./xmlchange -file env_build.xml -id USE_ESMF_LIB -val "FALSE" # FALSE is default in clm5.0; since cesm1_2 ESMF is no longer necessary to run with calendar=gregorian
#./xmlchange -file env_build.xml -id ESMF_LIBDIR -val ".../lib/libO/Linux.pgi.64.mpiuni.default" # path to ESMF library; can be set in config_machines.xml


#==========================================
# Modify user namelists
#==========================================
print_log "*** Modifying unser_nl_*.xml  ***"

# Output frequency and averaging (example)
# hist_empty_htapes = .true. # turn off all default output on h0
# hist_fincl1 or hist_fexcl1 # include or exclude selected variables
# hist_nhtfrq # output frequency
# hist_mfilt # number of values per file
# hist_avgflag_pertape # averaging over the output interval
# hist_dov2xy # true for 2D, false for 1D vector
# hist_type1d_pertape # '' for 2D and no averaging (i.e. PFT output), 'COL' for columns, 'LAND' for land-units, 'GRID' for grid-cells

# Commented out during testing to avoid lots of output
: '
cat > user_nl_clm << EOF
hist_fincl1 = 'TG', 'QAF', 'TAF', 'UAF'
hist_fincl2 = 'QAF', 'TAF', 'UAF', 'VPD_CAN', 'TLAI', 'FCEV', 'FCTR', 'TG', 'TSOI', 'TSOI_10CM', 'TSA', 'Q2M', 'VPD'
hist_fincl3 = 'QAF', 'TAF', 'UAF', 'VPD_CAN', 'TLAI', 'FCEV', 'FCTR', 'TG', 'TSOI', 'TSOI_10CM', 'TSA', 'Q2M', 'VPD'
hist_fincl4 = 'QAF', 'TAF', 'UAF', 'VPD_CAN', 'TLAI', 'FCEV', 'FCTR', 'TG', 'TSOI', 'TSOI_10CM', 'TSA', 'Q2M', 'VPD'

hist_nhtfrq = 0, -24, -24, -6 
hist_mfilt  = 12, 365, 365, 4 
hist_avgflag_pertape = 'A','M','X','I' 
hist_dov2xy = .true.,.true.,.true.,.false. 
hist_type1d_pertape = '','','',''
EOF

print_log "*** Output frequency and averaging  ***"
print_log "h0: default + selected variables, monthly values (0), yearly file (12 vals per file), average over the output interval (A)"
print_log "h1: selected variables, daily values (-24), yearly file (365 vals per file), min over the output interval (M)"
print_log "h2: selected variables, daily values (-24), yearly file (365 vals per file), max over the output interval (X)"
print_log "h3: selected variables, 6-hourly values (-6), daily file (4 vals per file), instantaneous at the output interval (I) by PFT"
'

# Params file: can be exchanged for newer versions
cat > user_nl_clm << EOF
paramfile = "$CESMDATAROOT/CCLM2_EUR_inputdata/CLM5params/clm5_params.cpbiomass.c190103.nc"
EOF

# Surface data (domain-specific), specify to re-use downloaded files
if [ $DOMAIN == eur ]; then
cat > user_nl_clm << EOF
fsurdat = "$CESMDATAROOT/CCLM2_EUR_inputdata/surfdata/surfdata_0.5x0.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190418.nc"
EOF
cat > user_nl_datm << EOF
domainfile = "$CESMDATAROOT/CCLM2_EUR_inputdata/domain/domain_EU-CORDEX_0.5_correctedlons.nc"
EOF
fi

if [ $DOMAIN == sa ]; then
cat > user_nl_clm << EOF
fsurdat = "$CESMDATAROOT/CCLM2_SA_inputdata/surfdata/surfdata_360x720cru_SA-CORDEX_16pfts_Irrig_CMIP6_simyr2000_c170824.nc"
EOF
cat > user_nl_datm << EOF
domainfile = "$CESMDATAROOT/CCLM2_SA_inputdata/domain/domain.lnd.360x720_SA-CORDEX_cruncep.100429.nc"
EOF
fi

if [ $DOMAIN == glob ]; then
cat > user_nl_clm << EOF
fsurdat = "$CESMDATAROOT/CCLM2_EUR_inputdata/surfdata/surfdata_360x720cru_16pfts_simyr2000_c170428.nc"
EOF
fi

# Namelist options available in Ronny's code
if [ $CODE == clm5.0_features ]; then
cat > user_nl_clm << EOF
use_biomass_heat_storage = .true.
use_individual_pft_soil_column = .true.
zetamaxstable = 100.0d00
EOF
fi

# Namelist options available in CTSMdev (?)
if [ $CODE == CTSMdev ]; then
cat > user_nl_clm << EOF
use_biomass_heat_storage = .true.
z0param_method = 'Meier2022'
zetamaxstable = 100.0d00
use_z0mg_2d = .true.
use_z0m_snowmelt = .true.
flanduse_timeseries=''
EOF
fi


#==========================================
# Set up the case
#==========================================

print_log "*** Running case.setup ***"
./case.setup -r | tee -a $logfile

#print_log "*** Downloading missing inputdata (if needed) ***"
#print_log "*** Consider transferring new data to PROJECT, e.g. rsync -av ${SCRATCH}/CCLM2_inputdata /project/${PROJ}/shared/CCLM2_inputdata ***"
#./check_input_data --download

#==========================================
# For OASIS coupling: before building, add the additional routines for OASIS interface in your CASEDIR on scratch
#==========================================

if [[ $COMPILER =~ "oasis" ]]; then
    print_log "*** Adding OASIS routines ***"
    ln -sf $CCLM2ROOT/cesm2_oas/src/oas/* SourceMods/src.drv/
    rm SourceMods/src.drv/oas_clm_vardef.F90
    ln -sf $CCLM2ROOT/cesm2_oas/src/drv/* SourceMods/src.drv/
    ln -sf $CCLM2ROOT/cesm2_oas/src/oas/oas_clm_vardef.F90 SourceMods/src.share/
    ln -sf $CCLM2ROOT/cesm2_oas/src/datm/* SourceMods/src.datm/
fi


#==========================================
# Build
#==========================================

print_log "*** Building case ***"
./case.build --clean-all | tee -a $logfile

if [ $CODE == clm5.0_features ]; then
    ./case.build --skip-provenance-check | tee -a $logfile # needed with Ronny's old code base
else
    ./case.build | tee -a $logfile
fi

print_log "*** Finished building new case in ${CASEDIR} ***"


#==========================================
# Preview and submit job
#==========================================

print_log "*** Preview the run ***"
./preview_run | tee -a $logfile

print_log "*** Submitting job ***"
./case.submit -a "-C gpu" | tee -a $logfile

# fails for clm5.0_features because tasks-per-node evaluates to float (12.0) with python 3. Cannot find where the calculation is made. Can also not override it like this:
# ./case.submit -a "-C gpu -p normal --ntasks-per-node 12" 
# or by setting in config_batch.xml

squeue --user=$USER | tee -a $logfile
#less CaseStatus

enddate=`date +'%Y-%m-%d %H:%M:%S'`
duration=$SECONDS
print_log "Started at: $startdate"
print_log "Finished at: $enddate"
print_log "Duration to create, setup, build, submit: $(($duration / 60)) min $(($duration % 60)) sec"

print_log "*** Check the job: squeue --user=${USER} ***"
print_log "*** Check the case: in ${CASEDIR}, run less CaseStatus ***"
print_log "*** Output at: ${CESMOUTPUTROOT}, run less CaseStatus ***"


#==========================================
# Copy final CaseStatus to logs
#==========================================

# Notes:
#env_case = model version, components, resolution, machine, compiler [do not modify]
#env_mach_pes = NTASKS, number of MPI tasks (or nodes if neg. values) [modify before setup]
#env_mach_specific = controls machine specific environment [modify before setup]
#env_build = component settings [modify before build]
#env_batch = batch job settings [modify any time]
#env_run = run settings incl runtype, coupling, pyhsics/sp/bgc and output [modify any time]
#env_workflow = walltime, queue, project
