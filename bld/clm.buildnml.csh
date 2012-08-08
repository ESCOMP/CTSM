#! /bin/csh -f

if !(-d $CASEBUILD/clmconf) mkdir -p $CASEBUILD/clmconf
cd $CASEBUILD/clmconf  

#--------------------------------------------------------------------
# Invoke clm configure - output will go in CASEBUILD/clmconf
#--------------------------------------------------------------------

# Set configure-time settings specific for regional grids
# This includes such things as turning RTM off
set config_opts=" "
if ($LND_GRID == "reg" && $GRID != "CLM_USRDAT" ) then
   set config_opts=" -sitespf_pt $GRID"
endif
if ("$CCSM_COMPSET" =~ P* || "$CCSM_COMPSET" =~ R* ) then
   set config_opts=" -sitespf_pt $LND_GRID"
endif

if ($COMP_INTERFACE == 'MCT' ) setenv COMP MCT
if ($COMP_INTERFACE == 'ESMF') setenv COMP ESMF

$CODEROOT/lnd/clm/bld/configure  $config_opts -comp_intf $COMP \
    $CLM_CONFIG_OPTS -usr_src $CASEROOT/SourceMods/src.clm || exit -1 

#--------------------------------------------------------------------
# Create clm.buildnml.csh
#--------------------------------------------------------------------

if ($RUN_TYPE == startup ) then
   if ($CLM_FORCE_COLDSTART == on) then
     set START_TYPE = "cold"
   else
     set START_TYPE = "default"
   endif
else
   if ($RUN_TYPE == hybrid ) then
     set START_TYPE = "startup"
   else
     set START_TYPE = $RUN_TYPE
   endif
endif

set RESOLUTION = $LND_GRID
if ($LND_GRID == reg ) then
   if ( $GRID == CLM_USRDAT ) then
      set RESOLUTION = $CLM_USRDAT_NAME
   else
      set RESOLUTION = $GRID
   endif
endif

set mask = " "
if ($OCN_GRID == "reg") then
   # Don't set mask in this case
else if ($ATM_GRID != $OCN_GRID) then
   set mask = "-mask $OCN_GRID"
endif

set clm_usr_name = ""
if ( "$CLM_USRDAT_NAME" != "UNSET" ) then
    set clm_usr_name = "-clm_usr_name $CLM_USRDAT_NAME"
endif

set glc_opts = ""
if ("$COMP_GLC" != "sglc" )then
   set glc_opts = "-glc_grid $GLC_GRID -glc_smb .$GLC_SMB. "
endif

set ignore = "-ignore_ic_date"
if ($RUN_STARTDATE =~ *-01-01* || $RUN_STARTDATE =~ *-09-01*) then
  set ignore = "-ignore_ic_year"
endif

set usecase = " "
if ($CLM_NML_USE_CASE != "UNSET") set usecase = "-use_case $CLM_NML_USE_CASE"

set clm_startfile = " "
if ( $RUN_TYPE == "hybrid" || $RUN_TYPE == "branch" ) then
   set clm_startfile = "-clm_startfile ${RUN_REFCASE}.clm2%inst_string.r.${RUN_REFDATE}-${RUN_REFTOD}.nc"
endif

set rtm = " "
if ( "$PTS_MODE" == TRUE   ) set rtm = "-rtm off"

set default_lnd_in_filename = "lnd_in"

set inst_counter = 1
while ($inst_counter <= $NINST_LND)

set instval=""
if ($NINST_LND > 1) then
   set inst_string = $inst_counter
   if ($inst_counter <= 999) set inst_string = "0$inst_string"
   if ($inst_counter <=  99) set inst_string = "0$inst_string"
   if ($inst_counter <=   9) set inst_string = "0$inst_string"
   set inst_string = "_${inst_string}"    
   set instval="-inst_string $inst_string"
else
   set inst_string = ""       
endif
set lnd_in_filename = ${default_lnd_in_filename}${inst_string}

if (-e $CASEROOT/user_nl_clm${inst_string}) then
  $UTILROOT/Tools/user_nlcreate -user_nl_file $CASEROOT/user_nl_clm${inst_string} -namelist_name clm_inparm >! $CASEBUILD/clmconf/cesm_namelist  || exit -2
endif

if (-e $CASEBUILD/clm.input_data_list) rm $CASEBUILD/clm.input_data_list

setenv INST_STRING $inst_string
$CODEROOT/lnd/clm/bld/build-namelist -infile $CASEBUILD/clmconf/cesm_namelist \
    -csmdata $DIN_LOC_ROOT  \
    -inputdata $CASEBUILD/clm.input_data_list \
    -namelist "&clm_inparm $CLM_NAMELIST_OPTS /" $instval \
    $ignore $usecase $clm_usr_name $glc_opts -co2_type $CLM_CO2_TYPE \
    -res $RESOLUTION $mask -clm_start_type $START_TYPE $clm_startfile \
    -l_ncpl $LND_NCPL -lnd_frac "${LND_DOMAIN_PATH}/${LND_DOMAIN_FILE}" \
    -rtm_tstep 10800 $rtm -glc_nec $GLC_NEC -co2_ppmv $CCSM_CO2_PPMV \
    -config $CASEBUILD/clmconf/config_cache.xml $CLM_BLDNML_OPTS
    
if (-d ${RUNDIR}) then
  cp $CASEBUILD/clmconf/lnd_in ${RUNDIR}/$lnd_in_filename || exit -2
endif

@ inst_counter = $inst_counter + 1

end


