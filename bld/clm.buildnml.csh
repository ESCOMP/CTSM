#! /bin/csh -f

if !(-d $CASEBUILD/clmconf) mkdir -p $CASEBUILD/clmconf

#--------------------------------------------------------------------
# Invoke clm configure - output will go in CASEBUILD/clmconf
#--------------------------------------------------------------------

set config_opts=" "
if ($LND_GRID == "reg" && $GRID != "CLM_USRDAT" ) then
   set config_opts=" -sitespf_pt $GRID"
endif
if ("$CCSM_COMPSET" =~ P* || "$CCSM_COMPSET" =~ R* ) then
   set config_opts=" -sitespf_pt $LND_GRID"
endif

cd $CASEBUILD/clmconf  
$CODEROOT/lnd/clm/bld/configure  $config_opts -comp_intf $COMP_INTERFACE \
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
set clmusr     = ""
if ($LND_GRID == reg ) then
   if ( $GRID == CLM_USRDAT ) then
      set RESOLUTION = $CLM_USRDAT_NAME
      set clmusr     = " -clm_usr_name $CLM_USRDAT_NAME"
   else
      set RESOLUTION = $GRID
   endif
endif

set default_lnd_in_filename = "lnd_in"

set inst_counter = 1
while ($inst_counter <= $NINST_LND)

if ($NINST_LND > 1) then
   set inst_string = $inst_counter
   if ($inst_counter <= 999) set inst_string = "0$inst_string"
   if ($inst_counter <=  99) set inst_string = "0$inst_string"
   if ($inst_counter <=   9) set inst_string = "0$inst_string"
   set inst_string = "_${inst_string}"    
else
   set inst_string = ""       
endif
set lnd_in_filename = ${default_lnd_in_filename}${inst_string}

setenv INST_STRING $inst_string

cd $CASEBUILD/clmconf  

if (-e $CASEBUILD/clm.input_data_list) rm $CASEBUILD/clm.input_data_list

if (-e $CASEROOT/user_nl_clm${inst_string}) then
  $UTILROOT/Tools/user_nlcreate -user_nl_file $CASEROOT/user_nl_clm${inst_string} \
    -namelist_name clm_inparm >! $CASEBUILD/clmconf/cesm_namelist  || exit -2
endif

set glc_opts = ""
if ("$COMP_GLC" != "sglc" )then
   set glc_opts = "-glc_grid $GLC_GRID -glc_smb .$GLC_SMB. "
endif

set usecase = " "
if ($CLM_NML_USE_CASE != "UNSET") set usecase = "-use_case $CLM_NML_USE_CASE"

set clm_startfile = " "
if ( $RUN_TYPE == "hybrid" || $RUN_TYPE == "branch" ) then
   set clm_startfile = "-clm_startfile ${RUN_REFCASE}.clm2.r.${RUN_REFDATE}-${RUN_REFTOD}.nc"
endif

$CODEROOT/lnd/clm/bld/build-namelist -infile $CASEBUILD/clmconf/cesm_namelist \
    -csmdata $DIN_LOC_ROOT  \
    -inputdata $CASEBUILD/clm.input_data_list \
    -namelist "&clm_inparm $CLM_NAMELIST_OPTS /" $usecase $glc_opts \
    -res $RESOLUTION $clmusr -clm_start_type $START_TYPE $clm_startfile \
    -l_ncpl $LND_NCPL -lnd_frac "${LND_DOMAIN_PATH}/${LND_DOMAIN_FILE}" \
    -glc_nec $GLC_NEC -co2_ppmv $CCSM_CO2_PPMV \
    -config $CASEBUILD/clmconf/config_cache.xml $CLM_BLDNML_OPTS || exit -3
    
if (-d ${RUNDIR}) then
  cp $CASEBUILD/clmconf/lnd_in ${RUNDIR}/$lnd_in_filename || exit -2
  # Only copy drv_flds_in namelist file if one doesn't already exist
  if ( ! -f "${RUNDIR}/drv_flds_in" ) cp $CASEBUILD/clmconf/drv_flds_in ${RUNDIR}/. >& /dev/null
endif

@ inst_counter = $inst_counter + 1

end


