module lnd_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for CTSM
  !----------------------------------------------------------------------------

  use ESMF
  use NUOPC                  , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC                  , only : NUOPC_CompFilterPhaseMap, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
  use NUOPC_Model            , only : model_routine_SS           => SetServices
  use NUOPC_Model            , only : SetVM
  use NUOPC_Model            , only : model_label_Advance        => label_Advance
  use NUOPC_Model            , only : model_label_DataInitialize => label_DataInitialize
  use NUOPC_Model            , only : model_label_SetRunClock    => label_SetRunClock
  use NUOPC_Model            , only : model_label_Finalize       => label_Finalize
  use NUOPC_Model            , only : NUOPC_ModelGet
  use shr_kind_mod           , only : r8 => shr_kind_r8, cl=>shr_kind_cl
  use shr_sys_mod            , only : shr_sys_abort
  use shr_file_mod           , only : shr_file_getlogunit, shr_file_setlogunit
  use shr_orb_mod            , only : shr_orb_decl, shr_orb_params, SHR_ORB_UNDEF_REAL, SHR_ORB_UNDEF_INT
  use shr_cal_mod            , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_ymd2date
  use spmdMod                , only : masterproc, mpicom, spmd_init
  use decompMod              , only : bounds_type, ldecomp, get_proc_bounds
  use domainMod              , only : ldomain
  use controlMod             , only : control_setNL
  use clm_varorb             , only : eccen, obliqr, lambm0, mvelpp
  use clm_varctl             , only : inst_index, inst_suffix, inst_name
  use clm_varctl             , only : single_column, clm_varctl_set, iulog
  use clm_varctl             , only : nsrStartup, nsrContinue, nsrBranch
  use clm_varcon             , only : re
  use clm_time_manager       , only : set_timemgr_init, advance_timestep
  use clm_time_manager       , only : set_nextsw_cday, update_rad_dtime
  use clm_time_manager       , only : get_nstep, get_step_size
  use clm_time_manager       , only : get_curr_date, get_curr_calday
  use clm_initializeMod      , only : initialize1, initialize2
  use clm_driver             , only : clm_drv
  use nuopc_shr_methods      , only : chkerr, state_setscalar, state_getscalar, state_diagnose, alarmInit
  use nuopc_shr_methods      , only : set_component_logging, get_component_instance, log_clock_advance
  use perf_mod               , only : t_startf, t_stopf, t_barrierf
  use netcdf                 , only : nf90_open, nf90_nowrite, nf90_noerr, nf90_close, nf90_strerror
  use netcdf                 , only : nf90_inq_dimid, nf90_inq_varid, nf90_get_var
  use netcdf                 , only : nf90_inquire_dimension, nf90_inquire_variable
  use lnd_import_export      , only : advertise_fields, realize_fields, import_fields, export_fields
  use lnd_comp_shr           , only : mesh, model_meshfile, model_clock

  implicit none
  private ! except

  ! Module routines
  public  :: SetServices
  public  :: SetVM
  private :: InitializeP0
  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelSetRunClock
  private :: ModelAdvance
  private :: ModelFinalize
  private :: clm_orbital_init
  private :: clm_orbital_update

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  character(len=CL)      :: flds_scalar_name = ''
  integer                :: flds_scalar_num = 0
  integer                :: flds_scalar_index_nx = 0
  integer                :: flds_scalar_index_ny = 0
  integer                :: flds_scalar_index_nextsw_cday = 0

  logical                :: glc_present
  logical                :: rof_prognostic
  integer, parameter     :: dbug = 0
  character(*),parameter :: modName =  "(lnd_comp_nuopc)"
  character(*),parameter :: u_FILE_u = &
       __FILE__

  character(len=CL)      :: orb_mode        ! attribute - orbital mode
  integer                :: orb_iyear       ! attribute - orbital year
  integer                :: orb_iyear_align ! attribute - associated with model year
  real(R8)               :: orb_obliq       ! attribute - obliquity in degrees
  real(R8)               :: orb_mvelp       ! attribute - moving vernal equinox longitude
  real(R8)               :: orb_eccen       ! attribute and update-  orbital eccentricity

  character(len=*) , parameter :: orb_fixed_year       = 'fixed_year'
  character(len=*) , parameter :: orb_variable_year    = 'variable_year'
  character(len=*) , parameter :: orb_fixed_parameters = 'fixed_parameters'

!===============================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
         specRoutine=ModelSetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine SetServices

  !===============================================================================

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries

    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeP0

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM)      :: vm
    integer            :: lmpicom
    integer            :: ierr
    integer            :: n
    integer            :: localpet
    integer            :: compid      ! component id
    integer            :: shrlogunit  ! original log unit
    character(len=CL)  :: cvalue
    character(len=CL)  :: logmsg
    logical            :: isPresent, isSet
    logical            :: cism_evolve
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    character(len=*), parameter :: format = "('("//trim(subname)//") :',A)"
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------------------------------------------------------------
    ! generate local mpi comm
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, localpet=localpet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! initialize CTSM MPI info
    !----------------------------------------------------------------------------

    call mpi_comm_dup(lmpicom, mpicom, ierr)

    !----------------------------------------------------------------------------
    ! reset shr logging to my log file
    !----------------------------------------------------------------------------

    call set_component_logging(gcomp, localPet==0, iulog, shrlogunit, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Note still need compid for those parts of the code that use the data model
    ! functionality through subroutine calls
    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) compid  ! convert from string to integer

    call spmd_init(mpicom, compid)

    !----------------------------------------------------------------------------
    ! determine instance information
    !----------------------------------------------------------------------------

    call get_component_instance(gcomp, inst_suffix, inst_index, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    inst_name = 'LND'

    !----------------------------------------------------------------------------
    ! advertise fields
    !----------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       flds_scalar_name = trim(cvalue)
       call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldName')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue, *) flds_scalar_num
       write(logmsg,*) flds_scalar_num
       call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldCount')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nx
       write(logmsg,*) flds_scalar_index_nx
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nx = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxGridNX')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_ny
       write(logmsg,*) flds_scalar_index_ny
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ny = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxGridNY')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxNextSwCday", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nextsw_cday
       write(logmsg,*) flds_scalar_index_nextsw_cday
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nextsw_cday = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxNextSwCday')
    endif

    call NUOPC_CompAttributeGet(gcomp, name='ROF_model', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (trim(cvalue) == 'srof' .or. trim(cvalue) == 'drof') then
       rof_prognostic = .false.
    else
       rof_prognostic = .true.
    end if

    call NUOPC_CompAttributeGet(gcomp, name='GLC_model', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (trim(cvalue) == 'sglc') then
       glc_present = .false.
    else
       glc_present = .true.
       cism_evolve = .true.
       call NUOPC_CompAttributeGet(gcomp, name="cism_evolve", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (isPresent .and. isSet) then
          call ESMF_LogWrite(trim(subname)//' cism_evolve = '//trim(cvalue), ESMF_LOGMSG_INFO)
          if (trim(cvalue) .eq. '.true.') then
             cism_evolve = .true.
          else if (trim(cvalue) .eq. '.false.') then
             cism_evolve = .false.
          else
             call shr_sys_abort(subname//'Could not determine cism_evolve value '//trim(cvalue))
          endif
       else
          call shr_sys_abort(subname//'Need to set cism_evolve if glc is present')
       endif
    end if

    if (masterproc) then
       write(iulog,*)' rof_prognostic = ',rof_prognostic
       write(iulog,*)' glc_present    = ',glc_present
       if (glc_present) write(iulog,*)' cism_evolve    = ',cism_evolve
    end if

    call advertise_fields(gcomp, flds_scalar_name, glc_present, cism_evolve, rof_prognostic, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)
    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)
  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    use clm_instMod, only : lnd2atm_inst, lnd2glc_inst, water_inst
!$  use omp_lib, only : omp_set_num_threads
    use ESMF, only : ESMF_VM, ESMF_VMGet
    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Mesh)         :: gridmesh              ! temporary esmf mesh
    type(ESMF_DistGrid)     :: DistGrid              ! esmf global index space descriptor
    type(ESMF_VM)           :: vm
    type(ESMF_Time)         :: currTime              ! Current time
    type(ESMF_Time)         :: startTime             ! Start time
    type(ESMF_Time)         :: refTime               ! Ref time
    type(ESMF_TimeInterval) :: timeStep              ! Model timestep
    type(ESMF_Calendar)     :: esmf_calendar         ! esmf calendar
    type(ESMF_CalKind_Flag) :: esmf_caltype          ! esmf calendar type
    integer                 :: ref_ymd               ! reference date (YYYYMMDD)
    integer                 :: ref_tod               ! reference time of day (sec)
    integer                 :: yy,mm,dd              ! Temporaries for time query
    integer                 :: start_ymd             ! start date (YYYYMMDD)
    integer                 :: start_tod             ! start time of day (sec)
    integer                 :: curr_ymd              ! Start date (YYYYMMDD)
    integer                 :: curr_tod              ! Start time of day (sec)
    integer                 :: dtime_sync            ! coupling time-step from the input synchronization clock
    integer                 :: localPet
    integer                 :: localpecount
    integer, pointer        :: gindex(:)             ! global index space for land and ocean points
    integer, pointer        :: gindex_lnd(:)         ! global index space for just land points
    integer, pointer        :: gindex_ocn(:)         ! global index space for just ocean points
    integer, pointer        :: mask(:)               ! local land/ocean mask
    character(ESMF_MAXSTR)  :: cvalue                ! config data
    integer                 :: nlnd, nocn            ! local size ofarrays
    integer                 :: g,n                   ! indices
    real(r8)                :: scmlat                ! single-column latitude
    real(r8)                :: scmlon                ! single-column longitude
    real(r8)                :: nextsw_cday           ! calday from clock of next radiation computation
    character(len=CL)       :: caseid                ! case identifier name
    character(len=CL)       :: ctitle                ! case description title
    character(len=CL)       :: starttype             ! start-type (startup, continue, branch, hybrid)
    character(len=CL)       :: calendar              ! calendar type name
    character(len=CL)       :: hostname              ! hostname of machine running on
    character(len=CL)       :: model_version         ! Model version
    character(len=CL)       :: username              ! user running the model
    integer                 :: nsrest                ! ctsm restart type
    logical                 :: brnch_retain_casename ! flag if should retain the case name on a branch start type
    integer                 :: lbnum                 ! input to memory diagnostic
    type(bounds_type)       :: bounds                ! bounds
    integer                 :: shrlogunit            ! original log unit
    real(r8)                :: mesh_lon, mesh_lat, mesh_area
    real(r8)                :: tolerance_latlon = 1.e-5
    real(r8)                :: tolerance_area   = 1.e-3
    integer                 :: spatialDim
    integer                 :: numOwnedElements
    real(R8), pointer       :: ownedElemCoords(:)
    real(r8), pointer       :: areaPtr(:)
    type(ESMF_Field)        :: areaField
    integer                 :: dimid_ni, dimid_nj, dimid_nv
    integer                 :: ncid, ierr
    integer                 :: ni, nj, nv
    integer                 :: varid_xv, varid_yv
    real(r8), allocatable   :: xv(:,:,:), yv(:,:,:)
    integer                 :: maxIndex(2)
    real(r8)                :: mincornerCoord(2)
    real(r8)                :: maxcornerCoord(2)
    type(ESMF_Grid)         :: lgrid
    character(len=*),parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (iulog)
#if (defined _MEMTRACE)
    if (masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_comp_nuopc_InitializeRealize:start::',lbnum)
    endif
#endif

    call ESMF_GridCompGet(gcomp, vm=vm, localPet=localPet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, pet=localPet, peCount=localPeCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

!$  call omp_set_num_threads(localPeCount)

    !----------------------
    ! Obtain attribute values
    !----------------------

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) caseid
    ctitle= trim(caseid)

    call NUOPC_CompAttributeGet(gcomp, name='scmlon', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlon

    call NUOPC_CompAttributeGet(gcomp, name='scmlat', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlat

    call NUOPC_CompAttributeGet(gcomp, name='single_column', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) single_column

    call NUOPC_CompAttributeGet(gcomp, name='brnch_retain_casename', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) brnch_retain_casename

    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) starttype

    call NUOPC_CompAttributeGet(gcomp, name='model_version', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) model_version

    call NUOPC_CompAttributeGet(gcomp, name='hostname', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) hostname

    call NUOPC_CompAttributeGet(gcomp, name='username', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) username

    !TODO: the following strings must not be hard-wired - must have module variables
    if (     trim(starttype) == trim('startup')) then
       nsrest = nsrStartup
    else if (trim(starttype) == trim('continue') ) then
       nsrest = nsrContinue
    else if (trim(starttype) == trim('branch')) then
       nsrest = nsrBranch
    else
       call shr_sys_abort( subname//' ERROR: unknown starttype' )
    end if

    !----------------------
    ! Consistency check on namelist filename
    !----------------------

    call control_setNL("lnd_in"//trim(inst_suffix))

    !----------------------
    ! Get properties from clock
    !----------------------

    call ESMF_ClockGet( clock, &
         currTime=currTime, startTime=startTime, refTime=RefTime, &
         timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=yy, mm=mm, dd=dd, s=curr_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,curr_ymd)

    call ESMF_TimeGet( startTime, yy=yy, mm=mm, dd=dd, s=start_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,start_ymd)

    call ESMF_TimeGet( refTime, yy=yy, mm=mm, dd=dd, s=ref_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,ref_ymd)

    call ESMF_TimeGet( currTime, calkindflag=esmf_caltype, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar = shr_cal_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar = shr_cal_gregorian
    else
       call shr_sys_abort( subname//'ERROR:: bad calendar for ESMF' )
    end if

    call ESMF_TimeIntervalGet( timeStep, s=dtime_sync, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (masterproc) then
       write(iulog,*)'dtime = ', dtime_sync
    end if

    !----------------------
    ! Initialize module orbital values and update orbital
    !----------------------

    call clm_orbital_init(gcomp, iulog, masterproc, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call clm_orbital_update(clock, iulog, masterproc, eccen, obliqr, lambm0, mvelpp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Initialize CTSM time manager
    !----------------------

    ! Note that we assume that CTSM's internal dtime matches the coupling time step.
    ! i.e., we currently do NOT allow sub-cycling within a coupling time step.
    call set_timemgr_init(       &
         calendar_in=calendar,   &
         start_ymd_in=start_ymd, &
         start_tod_in=start_tod, &
         ref_ymd_in=ref_ymd,     &
         ref_tod_in=ref_tod,     &
         dtime_in=dtime_sync)

    !----------------------------------------------------------------------------
    ! Set model clock in lnd_comp_shr
    !----------------------------------------------------------------------------

    model_clock = clock

    !----------------------
    ! Read namelist, grid and surface data
    !----------------------

    ! set default values for run control variables
    call clm_varctl_set(&
         caseid_in=caseid, ctitle_in=ctitle,                     &
         brnch_retain_casename_in=brnch_retain_casename,         &
         single_column_in=single_column, scmlat_in=scmlat, scmlon_in=scmlon, &
         nsrest_in=nsrest, &
         version_in=model_version, &
         hostname_in=hostname, &
         username_in=username)

    ! note that the memory for gindex_ocn will be allocated in the following call
    call initialize1(dtime=dtime_sync, gindex_ocn=gindex_ocn)

    ! If no land then abort for now
    ! TODO: need to handle the case of noland with CMEPS
    ! if ( noland ) then
    !    call shr_sys_abort(trim(subname)//"ERROR: Currently cannot handle case of single column with non-land")
    ! end if

    ! obtain global index array for just land points which includes mask=0 or ocean points
    call get_proc_bounds( bounds )
    nlnd = bounds%endg - bounds%begg + 1
    allocate(gindex_lnd(nlnd))
    do g = bounds%begg,bounds%endg
       n = 1 + (g - bounds%begg)
       gindex_lnd(n) = ldecomp%gdc2glo(g)
    end do

    ! create a global index that includes both land and ocean points
    nocn = size(gindex_ocn)
    allocate(gindex(nlnd + nocn))
    allocate(mask(nlnd + nocn))
    do n = 1,nlnd+nocn
       if (n <= nlnd) then
          gindex(n) = gindex_lnd(n)
          mask(n) = 1
       else
          gindex(n) = gindex_ocn(n-nlnd)
          mask(n) = 0
       end if
    end do

    ! create distGrid from global index array
    DistGrid = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    deallocate(gindex)

    !--------------------------------
    ! generate the mesh and realize fields
    !--------------------------------

    ! determine if the mesh will be created or read in
    call NUOPC_CompAttributeGet(gcomp, name='mesh_lnd', value=model_meshfile, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (single_column) model_meshfile = 'create_mesh'

    if (trim(model_meshfile) == 'create_mesh') then
       ! get the datm grid from the domain file
       call NUOPC_CompAttributeGet(gcomp, name='domain_lnd', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! open file
       ierr = nf90_open(cvalue, NF90_NOWRITE, ncid)
       call nc_check_err(ierr, 'nf90_open', trim(cvalue))
       ! get dimension ids
       ierr = nf90_inq_dimid(ncid, 'ni', dimid_ni)
       call nc_check_err(ierr, 'nf90_inq_dimid for ni', trim(cvalue))
       ierr = nf90_inq_dimid(ncid, 'nj', dimid_nj)
       call nc_check_err(ierr, 'nf90_inq_dimid for nj', trim(cvalue))
       ierr = nf90_inq_dimid(ncid, 'nv', dimid_nv)
       call nc_check_err(ierr, 'nf90_inq_dimid for nv', trim(cvalue))
       ! get dimension values
       ierr = nf90_inquire_dimension(ncid, dimid_ni, len=ni)
       call nc_check_err(ierr, 'nf90_inq_dimension for ni', trim(cvalue))
       ierr = nf90_inquire_dimension(ncid, dimid_nj, len=nj)
       call nc_check_err(ierr, 'nf90_inq_dimension for nj', trim(cvalue))
       ierr = nf90_inquire_dimension(ncid, dimid_nv, len=nv)
       call nc_check_err(ierr, 'nf90_inq_dimension for nv', trim(cvalue))
       ! get variable ids
       ierr = nf90_inq_varid(ncid, 'xv', varid_xv)
       call nc_check_err(ierr, 'nf90_inq_varid for xv', trim(cvalue))
       ierr = nf90_inq_varid(ncid, 'yv', varid_yv)
       call nc_check_err(ierr, 'nf90_inq_varid for yv', trim(cvalue))
       ! allocate memory for variables and get variable values
       allocate(xv(nv,ni,nj), yv(nv,ni,nj))
       ierr = nf90_get_var(ncid, varid_xv, xv)
       call nc_check_err(ierr, 'nf90_get_var for xv', trim(cvalue))
       ierr = nf90_get_var(ncid, varid_yv, yv)
       call nc_check_err(ierr, 'nf90_get_var for yv', trim(cvalue))
       ! close file
       ierr = nf90_close(ncid)
       call nc_check_err(ierr, 'nf90_close', trim(cvalue))
       ! create the grid
       maxIndex(1)       = ni          ! number of lons
       maxIndex(2)       = nj          ! number of lats
       mincornerCoord(1) = xv(1,1,1)   ! min lon
       mincornerCoord(2) = yv(1,1,1)   ! min lat
       maxcornerCoord(1) = xv(3,ni,nj) ! max lon
       maxcornerCoord(2) = yv(3,ni,nj) ! max lat
       deallocate(xv,yv)
       lgrid = ESMF_GridCreateNoPeriDimUfrm (maxindex=maxindex, &
            mincornercoord=mincornercoord, maxcornercoord= maxcornercoord, &
            staggerloclist=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! create the mesh from the grid
       mesh =  ESMF_MeshCreate(lgrid, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! TODO: is the mask by default set to 1 if created from a grid?
       ! reset the global mask (which is 1) to the land/ocean mask
       !
       ! Currently, this call requires that the information has
       ! already been added to the mesh during creation. For example,
       ! you can only change the element mask information, if the mesh
       ! was initially created with element masking.
       !!! call ESMF_MeshSet(mesh, elementMask=mask, rc=rc)
       !!! if (ChkErr(rc,__LINE__,u_FILE_u)) return

       deallocate(mask)

    else

       ! read in the mesh from the file
       mesh = ESMF_MeshCreate(filename=trim(model_meshfile), fileformat=ESMF_FILEFORMAT_ESMFMESH, &
            elementDistgrid=Distgrid, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (masterproc) then
          write(iulog,*)'mesh file for domain is ',trim(model_meshfile)
       end if

       ! Determine the areas on the mesh
       areaField = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name='mesh_areas', meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldRegridGetArea(areaField, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(areaField, farrayPtr=areaPtr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    ! realize the actively coupled fields
    call realize_fields(gcomp,  mesh, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Finish initializing ctsm
    !--------------------------------

    call initialize2()

    !--------------------------------
    ! Check that lats, lons and areas on mesh are the same as those internal to ctsm
    ! obtain mesh lats and lons
    !--------------------------------

    if (trim(model_meshfile) /= 'create_mesh') then
       call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(ownedElemCoords(spatialDim*numOwnedElements))
       call ESMF_MeshGet(mesh, ownedElemCoords=ownedElemCoords)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do g = bounds%begg,bounds%endg
          n = 1 + (g - bounds%begg)
          mesh_lon = ownedElemCoords(2*n-1)
          mesh_lat = ownedElemCoords(2*n)
          mesh_area = areaPtr(n)
          if (abs(mesh_lon - ldomain%lonc(g)) > tolerance_latlon) then
             write(6,100)'ERROR: clm_lon, mesh_lon, diff_lon = ',&
                  ldomain%lonc(g), mesh_lon, abs(mesh_lon - ldomain%lonc(g))
             !call shr_sys_abort()
          end if
          if (abs(mesh_lat - ldomain%latc(g)) > tolerance_latlon) then
             write(6,100)'ERROR: clm_lat, mesh_lat, diff_lat = ',&
                  ldomain%latc(g), mesh_lat, abs(mesh_lat - ldomain%latc(g))
             !call shr_sys_abort()
          end if
          if (abs(mesh_area - ldomain%area(g)/(re*re)) > tolerance_area) then
             write(6,100)'ERROR: clm_area, mesh_area, diff_area = ',&
                  ldomain%area(g)/(re*re), mesh_area, abs(mesh_area - ldomain%area(g)/(re*re))
             !call shr_sys_abort()
          end if
       end do
100    format(a,3(d13.5,2x))
    end if

    !--------------------------------
    ! Create land export state
    !--------------------------------

    call export_fields(gcomp, bounds, glc_present, rof_prognostic, &
         water_inst%waterlnd2atmbulk_inst, lnd2atm_inst, lnd2glc_inst, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get calendar day of nextsw calculation
    if (nsrest == nsrStartup) then
       call ESMF_ClockGet( clock, currTime=currTime, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_TimeGet( currTime, dayOfYear_r8=nextsw_cday, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call set_nextsw_cday(nextsw_cday)

    ! Set scalars in export state
    call State_SetScalar(dble(ldomain%ni), flds_scalar_index_nx, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call State_SetScalar(dble(ldomain%nj), flds_scalar_index_ny, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       call State_diagnose(exportState, subname//':ExportState',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call shr_file_setLogUnit (shrlogunit)

#if (defined _MEMTRACE)
    if(masterproc) then
       write(iulog,*) TRIM(subname) // ':end::'
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_comp_nuopc_InitializeRealize:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine InitializeRealize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)

    !------------------------
    ! Run CTSM
    !------------------------

    use clm_instMod, only : water_inst, atm2lnd_inst, glc2lnd_inst, lnd2atm_inst, lnd2glc_inst
    !$  use omp_lib, only : omp_set_num_threads
    use ESMF, only : ESMF_VM, ESMF_VMGet

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables:
    type(ESMF_Clock)       :: clock
    type(ESMF_Alarm)       :: alarm
    type(ESMF_Time)        :: currTime
    type(ESMF_Time)        :: nextTime
    type(ESMF_State)       :: importState, exportState
    type(ESMF_VM)          :: vm
    character(ESMF_MAXSTR) :: cvalue
    character(ESMF_MAXSTR) :: case_name      ! case name
    integer                :: ymd            ! CTSM current date (YYYYMMDD)
    integer                :: yr             ! CTSM current year
    integer                :: mon            ! CTSM current month
    integer                :: day            ! CTSM current day
    integer                :: tod            ! CTSM current time of day (sec)
    integer                :: ymd_sync       ! Sync date (YYYYMMDD)
    integer                :: yr_sync        ! Sync current year
    integer                :: mon_sync       ! Sync current month
    integer                :: day_sync       ! Sync current day
    integer                :: tod_sync       ! Sync current time of day (sec)
    integer                :: dtime          ! time step increment (sec)
    integer                :: nstep          ! time step index
    integer                :: localPet
    integer                :: localpecount
    logical                :: rstwr          ! .true. ==> write restart file before returning
    logical                :: nlend          ! .true. ==> last time-step
    logical                :: dosend         ! true => send data back to driver
    logical                :: doalb          ! .true. ==> do albedo calculation on this time step
    real(r8)               :: nextsw_cday    ! calday from clock of next radiation computation
    real(r8)               :: caldayp1       ! ctsm calday plus dtime offset
    integer                :: lbnum          ! input to memory diagnostic
    integer                :: g,i            ! counters
    real(r8)               :: calday         ! calendar day for nstep
    real(r8)               :: declin         ! solar declination angle in radians for nstep
    real(r8)               :: declinp1       ! solar declination angle in radians for nstep+1
    real(r8)               :: eccf           ! earth orbit eccentricity factor
    type(bounds_type)      :: bounds         ! bounds
    character(len=32)      :: rdate          ! date char string for restart file names
    integer                :: shrlogunit ! original log unit
    character(len=*),parameter  :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call ESMF_GridCompGet(gcomp, vm=vm, localPet=localPet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, pet=localPet, peCount=localPeCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !$  call omp_set_num_threads(localPeCount)

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (iulog)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_comp_nuopc_ModelAdvance:start::',lbnum)
    endif
#endif

    !--------------------------------
    ! Query the Component for its clock, importState and exportState
    !--------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Determine time of next atmospheric shortwave calculation
    !--------------------------------

    call State_GetScalar(importState, &
         flds_scalar_index_nextsw_cday, nextsw_cday, &
         flds_scalar_name, flds_scalar_num, rc)
    call set_nextsw_cday( nextsw_cday )

    !----------------------
    ! Get orbital values
    !----------------------


    !--------------------------------
    ! Unpack import state
    !--------------------------------

    call t_startf ('lc_lnd_import')

    call get_proc_bounds(bounds)
    call import_fields( gcomp, bounds, glc_present, rof_prognostic, &
         atm2lnd_inst, glc2lnd_inst, water_inst%wateratm2lndbulk_inst, rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call t_stopf ('lc_lnd_import')

    !--------------------------------
    ! Run model
    !--------------------------------

    dtime = get_step_size()
    dosend = .false.
    do while(.not. dosend)

       ! TODO: This is currently hard-wired - is there a better way for nuopc?
       ! Note that the model clock is updated at the end of the time step not at the beginning
       nstep = get_nstep()
       if (nstep > 0) then
          dosend = .true.
       end if

       !--------------------------------
       ! Determine doalb based on nextsw_cday sent from atm model
       !--------------------------------

       caldayp1 = get_curr_calday(offset=dtime)

       if (nstep == 0) then
          doalb = .false.
       else if (nstep == 1) then
          doalb = (abs(nextsw_cday- caldayp1) < 1.e-10_r8)
       else
          doalb = (nextsw_cday >= -0.5_r8)
       end if
       call update_rad_dtime(doalb)

       !--------------------------------
       ! Determine if time to write restart
       !--------------------------------

       call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          rstwr = .true.
          call ESMF_AlarmRingerOff( alarm, rc=rc )
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          rstwr = .false.
       endif

       !--------------------------------
       ! Determine if time to stop
       !--------------------------------

       call ESMF_ClockGetAlarm(clock, alarmname='alarm_stop', alarm=alarm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          nlend = .true.
          call ESMF_AlarmRingerOff( alarm, rc=rc )
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          nlend = .false.
       endif

       !--------------------------------
       ! Run CTSM
       !--------------------------------

       call t_barrierf('sync_ctsm_run1', mpicom)

       call t_startf ('shr_orb_decl')

       ! Note - the orbital inquiries set the values in clm_varorb via the module use statements
       call  clm_orbital_update(clock, iulog, masterproc, eccen, obliqr, lambm0, mvelpp, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       calday = get_curr_calday()
       call shr_orb_decl( calday     , eccen, mvelpp, lambm0, obliqr, declin  , eccf )
       call shr_orb_decl( nextsw_cday, eccen, mvelpp, lambm0, obliqr, declinp1, eccf )
       call t_stopf ('shr_orb_decl')

       call t_startf ('ctsm_run')

       ! Restart File - use nexttimestr rather than currtimestr here since that is the time at the end of
       ! the timestep and is preferred for restart file names

       call ESMF_ClockGetNextTime(clock, nextTime=nextTime, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(nexttime, yy=yr_sync, mm=mon_sync, dd=day_sync, s=tod_sync, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr_sync, mon_sync, day_sync, tod_sync

       call clm_drv(doalb, nextsw_cday, declinp1, declin, rstwr, nlend, rdate, rof_prognostic)

       call t_stopf ('ctsm_run')

       !--------------------------------
       ! Pack export state
       !--------------------------------

       call t_startf ('lc_lnd_export')

       call export_fields(gcomp, bounds, glc_present, rof_prognostic, &
            water_inst%waterlnd2atmbulk_inst, lnd2atm_inst, lnd2glc_inst, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call t_stopf ('lc_lnd_export')

       !--------------------------------
       ! Advance ctsm time step
       !--------------------------------

       call t_startf ('lc_ctsm2_adv_timestep')
       call advance_timestep()
       call t_stopf ('lc_ctsm2_adv_timestep')

    end do

    ! Check that internal clock is in sync with master clock
    ! Note that the driver clock has not been updated yet - so at this point
    ! CTSM is actually 1 coupling intervals ahead of the driver clock

    call get_curr_date( yr, mon, day, tod, offset=-2*dtime )
    ymd = yr*10000 + mon*100 + day
    tod = tod

    call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=yr_sync, mm=mon_sync, dd=day_sync, s=tod_sync, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yr_sync, mon_sync, day_sync, ymd_sync)

    if ( (ymd /= ymd_sync) .or. (tod /= tod_sync) ) then
       write(iulog,*)'ctsm ymd=',ymd     ,' ctsm tod= ',tod
       write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       call ESMF_LogWrite(subname//" CTSM clock not in sync with Master Sync clock",ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       call State_diagnose(exportState,subname//':ES',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (masterproc) then
          call log_clock_advance(clock, 'CTSM', iulog, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    !--------------------------------
    ! Reset shr logging to my original values
    !--------------------------------

    call shr_file_setLogUnit (shrlogunit)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_comp_nuopc_ModelAdvance:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)         :: mclock, dclock
    type(ESMF_Time)          :: mcurrtime, dcurrtime
    type(ESMF_Time)          :: mstoptime
    type(ESMF_TimeInterval)  :: mtimestep, dtimestep
    character(len=256)       :: cvalue
    character(len=256)       :: restart_option ! Restart option units
    integer                  :: restart_n      ! Number until restart interval
    integer                  :: restart_ymd    ! Restart date (YYYYMMDD)
    type(ESMF_ALARM)         :: restart_alarm
    character(len=256)       :: stop_option    ! Stop option units
    integer                  :: stop_n         ! Number until stop interval
    integer                  :: stop_ymd       ! Stop date (YYYYMMDD)
    type(ESMF_ALARM)         :: stop_alarm
    character(len=128)       :: name
    integer                  :: alarmcount
    character(len=*),parameter :: subname=trim(modName)//':(ModelSetRunClock) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart and stop alarms
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then

       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for' // trim(name), ESMF_LOGMSG_INFO)

       !----------------
       ! Restart alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------
       ! Stop alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="stop_option", value=stop_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="stop_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_n

       call NUOPC_CompAttributeGet(gcomp, name="stop_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_ymd

       call alarmInit(mclock, stop_alarm, stop_option, &
            opt_n   = stop_n,           &
            opt_ymd = stop_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_stop', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(stop_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelSetRunClock

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(*), parameter :: F00   = "('(lnd_comp_nuopc) ',8a)"
    character(*), parameter :: F91   = "('(lnd_comp_nuopc) ',73('-'))"
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    !-------------------------------------------------------------------------------

    !--------------------------------
    ! Finalize routine
    !--------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    if (masterproc) then
       write(iulog,F91)
       write(iulog,F00) 'CTSM: end of main integration loop'
       write(iulog,F91)
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelFinalize

  !===============================================================================

  subroutine clm_orbital_init(gcomp, logunit, mastertask, rc)

    !----------------------------------------------------------
    ! Initialize orbital related values
    !----------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp) , intent(in)    :: gcomp
    integer             , intent(in)    :: logunit
    logical             , intent(in)    :: mastertask
    integer             , intent(out)   :: rc              ! output error

    ! local variables
    character(len=CL) :: msgstr          ! temporary
    character(len=CL) :: cvalue          ! temporary
    character(len=*) , parameter :: subname = "(clm_orbital_init)"
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine orbital attributes from input
    call NUOPC_CompAttributeGet(gcomp, name="orb_mode", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_mode

    call NUOPC_CompAttributeGet(gcomp, name="orb_iyear", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_iyear

    call NUOPC_CompAttributeGet(gcomp, name="orb_iyear_align", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_iyear_align

    call NUOPC_CompAttributeGet(gcomp, name="orb_obliq", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_obliq

    call NUOPC_CompAttributeGet(gcomp, name="orb_eccen", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_eccen

    call NUOPC_CompAttributeGet(gcomp, name="orb_mvelp", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_mvelp

    ! Error checks
    if (trim(orb_mode) == trim(orb_fixed_year)) then
       orb_obliq = SHR_ORB_UNDEF_REAL
       orb_eccen = SHR_ORB_UNDEF_REAL
       orb_mvelp = SHR_ORB_UNDEF_REAL
       if (orb_iyear == SHR_ORB_UNDEF_INT) then
          if (mastertask) then
             write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
             write(logunit,*) trim(subname),' ERROR: fixed_year settings = ',orb_iyear
             write (msgstr, *) ' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          end if
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       endif
    elseif (trim(orb_mode) == trim(orb_variable_year)) then
       orb_obliq = SHR_ORB_UNDEF_REAL
       orb_eccen = SHR_ORB_UNDEF_REAL
       orb_mvelp = SHR_ORB_UNDEF_REAL
       if (orb_iyear == SHR_ORB_UNDEF_INT .or. orb_iyear_align == SHR_ORB_UNDEF_INT) then
          if (mastertask) then
             write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
             write(logunit,*) trim(subname),' ERROR: variable_year settings = ',orb_iyear, orb_iyear_align
             write (msgstr, *) subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          end if
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       endif
    elseif (trim(orb_mode) == trim(orb_fixed_parameters)) then
       !-- force orb_iyear to undef to make sure shr_orb_params works properly
       orb_iyear = SHR_ORB_UNDEF_INT
       orb_iyear_align = SHR_ORB_UNDEF_INT
       if (orb_eccen == SHR_ORB_UNDEF_REAL .or. &
           orb_obliq == SHR_ORB_UNDEF_REAL .or. &
           orb_mvelp == SHR_ORB_UNDEF_REAL) then
          if (mastertask) then
             write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
             write(logunit,*) trim(subname),' ERROR: orb_eccen = ',orb_eccen
             write(logunit,*) trim(subname),' ERROR: orb_obliq = ',orb_obliq
             write(logunit,*) trim(subname),' ERROR: orb_mvelp = ',orb_mvelp
             write (msgstr, *) subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          end if
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       endif
    else
       write (msgstr, *) subname//' ERROR: invalid orb_mode '//trim(orb_mode)
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       rc = ESMF_FAILURE
       return  ! bail out
    endif

  end subroutine clm_orbital_init

  !===============================================================================

  subroutine clm_orbital_update(clock, logunit,  mastertask, eccen, obliqr, lambm0, mvelpp, rc)

    !----------------------------------------------------------
    ! Update orbital settings
    !----------------------------------------------------------

    ! input/output variables
    type(ESMF_Clock) , intent(in)    :: clock
    integer          , intent(in)    :: logunit
    logical          , intent(in)    :: mastertask
    real(R8)         , intent(inout) :: eccen  ! orbital eccentricity
    real(R8)         , intent(inout) :: obliqr ! Earths obliquity in rad
    real(R8)         , intent(inout) :: lambm0 ! Mean long of perihelion at vernal equinox (radians)
    real(R8)         , intent(inout) :: mvelpp ! moving vernal equinox longitude of perihelion plus pi (radians)
    integer          , intent(out)   :: rc     ! output error

    ! local variables
    type(ESMF_Time)   :: CurrTime ! current time
    integer           :: year     ! model year at current time
    integer           :: orb_year ! orbital year for current orbital computation
    character(len=CL) :: msgstr   ! temporary
    logical           :: lprint
    logical           :: first_time = .true.
    character(len=*) , parameter :: subname = "(clm_orbital_update)"
    !-------------------------------------------

    if (trim(orb_mode) == trim(orb_variable_year)) then
       call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(CurrTime, yy=year, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       orb_year = orb_iyear + (year - orb_iyear_align)
       lprint = mastertask
    else
       orb_year = orb_iyear
       if (first_time) then
          lprint = mastertask
          first_time = .false.
       else
          lprint = .false.
       end if
    end if

    eccen = orb_eccen
    call shr_orb_params(orb_year, eccen, orb_obliq, orb_mvelp, obliqr, lambm0, mvelpp, lprint)

    if ( eccen  == SHR_ORB_UNDEF_REAL .or. obliqr == SHR_ORB_UNDEF_REAL .or. &
         mvelpp == SHR_ORB_UNDEF_REAL .or. lambm0 == SHR_ORB_UNDEF_REAL) then
       write (msgstr, *) subname//' ERROR: orb params incorrect'
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       return  ! bail out
    endif

  end subroutine clm_orbital_update

  !===============================================================================

  subroutine nc_check_err(ierror, description, filename)
    integer     , intent(in) :: ierror
    character(*), intent(in) :: description
    character(*), intent(in) :: filename

    if (ierror /= nf90_noerr) then
       write (*,'(6a)') 'ERROR ', trim(description),'. NetCDF file : "', trim(filename),&
            '". Error message:', trim(nf90_strerror(ierror))
       call shr_sys_abort()
    endif
  end subroutine nc_check_err

end module lnd_comp_nuopc
