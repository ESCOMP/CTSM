module lnd_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for CTSM
  !----------------------------------------------------------------------------

  use ESMF                   , only : ESMF_SUCCESS, ESMF_GridComp, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_GridCompSetEntryPoint
  use ESMF                   , only : ESMF_METHOD_INITIALIZE, ESMF_FAILURE, ESMF_Time, ESMF_LogSetError, ESMF_RC_NOT_VALID
  use ESMF                   , only : ESMF_Clock, ESMF_State, ESMF_Field, ESMF_MAXSTR, ESMF_LOGMSG_WARNING,ESMF_RC_ARG_BAD
  use ESMF                   , only : ESMF_LOGFOUNDERROR, ESMF_LOGERR_PASSTHRU, ESMF_CALKIND_FLAG, ESMF_TIMEINTERVAL
  use ESMF                   , only : ESMF_LogSetError, ESMF_FieldGet, ESMF_ClockGet, ESMF_GridCompGet, ESMF_ClockGetNextTime
  use ESMF                   , only : ESMF_AlarmRingerOff, ESMF_TimeIntervalGet,  ESMF_TimeGet, ESMF_StateGet
  use ESMF                   , only : ESMF_MethodRemove, ESMF_VM, ESMF_VMGet, ESMF_CALKIND_NOLEAP, ESMF_CALKIND_GREGORIAN
  use ESMF                   , only : ESMF_ALARMLIST_ALL, ESMF_ALARM, ESMF_ALARMISRINGING,  ESMF_ClockGetAlarm, ESMF_ClockGetAlarmList
  use ESMF                   , only : ESMF_AlarmSet, ESMF_ClockAdvance
  use ESMF                   , only : operator(==), operator(+)
  use ESMF                   , only : ESMF_AlarmIsCreated, ESMF_LOGMSG_ERROR, ESMF_ClockSet
  use NUOPC                  , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC                  , only : NUOPC_CompFilterPhaseMap, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
  use NUOPC                  , only : NUOPC_CompGet, Nuopc_IsAtTime, Nuopc_GetAttribute
  use NUOPC_Model            , only : model_routine_SS           => SetServices
  use NUOPC_Model            , only : SetVM
  use NUOPC_Model            , only : model_label_Advance        => label_Advance
  use NUOPC_Model            , only : model_label_DataInitialize => label_DataInitialize
  use NUOPC_Model            , only : model_label_SetRunClock    => label_SetRunClock
  use NUOPC_Model            , only : model_label_Finalize       => label_Finalize
  use NUOPC_Model            , only : label_CheckImport

  use NUOPC_Model            , only : NUOPC_ModelGet
  use shr_kind_mod           , only : r8 => shr_kind_r8, cl=>shr_kind_cl
  use shr_sys_mod            , only : shr_sys_abort
  use shr_log_mod            , only : shr_log_setLogUnit, shr_log_getLogUnit
  use shr_orb_mod            , only : shr_orb_decl, shr_orb_params, SHR_ORB_UNDEF_REAL, SHR_ORB_UNDEF_INT
  use shr_cal_mod            , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_ymd2date
  use spmdMod                , only : masterproc, mpicom, spmd_init
  use controlMod             , only : control_setNL, control_init, control_print, NLFilename
  use clm_varorb             , only : eccen, obliqr, lambm0, mvelpp
  use clm_varctl             , only : inst_index, inst_suffix, inst_name
  use clm_varctl             , only : single_column, clm_varctl_set, iulog
  use clm_varctl             , only : nsrStartup, nsrContinue, nsrBranch
  use clm_varctl             , only : FL => fname_len
  use clm_time_manager       , only : set_timemgr_init, advance_timestep
  use clm_time_manager       , only : update_rad_dtime
  use clm_time_manager       , only : get_nstep, get_step_size
  use clm_time_manager       , only : get_curr_date, get_curr_calday
  use clm_initializeMod      , only : initialize1, initialize2
  use nuopc_shr_methods      , only : chkerr, state_setscalar, state_getscalar, state_diagnose, alarmInit
  use nuopc_shr_methods      , only : set_component_logging, get_component_instance, log_clock_advance
  use lnd_import_export      , only : advertise_fields, realize_fields, import_fields, export_fields
  use lnd_comp_shr           , only : mesh, model_meshfile, model_clock
  use perf_mod               , only : t_startf, t_stopf, t_barrierf

  implicit none
  private ! except

  ! Module public routines
  public  :: SetServices         ! Setup the pointers to the function calls for the different models phases (initialize, run, finalize)
  public  :: SetVM               ! Set the virtual machine description of the paralell model (both MPI and OpenMP)

  ! Module private routines
  private :: InitializeP0        ! Phase zero of initialization
  private :: InitializeAdvertise ! Advertise the fields that can be passed
  private :: InitializeRealize   ! Realize the list of fields that will be exchanged
  private :: ModelSetRunClock    ! Set the run clock
  private :: ModelAdvance        ! Advance the model
  private :: ModelFinalize       ! Finalize the model
  private :: clm_orbital_init    ! Initialize the orbital information
  private :: clm_orbital_update  ! Update the orbital information
  private :: CheckImport

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
  logical                :: atm_prognostic
  integer, parameter     :: dbug = 0
  character(*),parameter :: modName =  "(lnd_comp_nuopc)"

  character(len=CL)      :: orb_mode        ! attribute - orbital mode
  integer                :: orb_iyear       ! attribute - orbital year
  integer                :: orb_iyear_align ! attribute - associated with model year
  real(R8)               :: orb_obliq       ! attribute - obliquity in degrees
  real(R8)               :: orb_mvelp       ! attribute - moving vernal equinox longitude
  real(R8)               :: orb_eccen       ! attribute and update-  orbital eccentricity

  logical                :: scol_valid      ! if single_column, does point have a mask of zero

  integer                :: nthrds          ! Number of threads per task in this component

  character(len=*) , parameter :: orb_fixed_year       = 'fixed_year'
  character(len=*) , parameter :: orb_variable_year    = 'variable_year'
  character(len=*) , parameter :: orb_fixed_parameters = 'fixed_parameters'

  character(len=*) , parameter :: startup_run  = 'startup'
  character(len=*) , parameter :: continue_run = 'continue'
  character(len=*) , parameter :: branch_run   = 'branch'

  logical :: write_restart_at_endofrun = .false.

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)
    ! Setup the pointers to the function calls for the different models phases (initialize, run, finalize)
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

    call NUOPC_CompSpecialize(gcomp, specLabel=label_CheckImport, &
         specRoutine=CheckImport, rc=rc)
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

    ! Phase zero initialization
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

    ! Advertise the fields that can be exchanged
    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM)     :: vm
    integer           :: lmpicom
    integer           :: ierr
    integer           :: n
    integer           :: localPet    ! local PET (Persistent Execution Threads) (both MPI tasks and OpenMP threads)
    integer           :: compid      ! component id
    integer           :: shrlogunit  ! original log unit
    character(len=CL) :: cvalue
    character(len=CL) :: logmsg
    logical           :: cism_evolve
    character(len=CL) :: atm_model
    character(len=CL) :: rof_model
    character(len=CL) :: glc_model
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
    ! functionality through subroutine calls (MCTID just means the Model ComonenT IDentification number)
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

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    flds_scalar_name = trim(cvalue)
    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue, *) flds_scalar_num
    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_scalar_index_nx
    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_scalar_index_ny
    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxNextSwCday", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_scalar_index_nextsw_cday
    call NUOPC_CompAttributeGet(gcomp, name='ROF_model', value=rof_model, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (trim(rof_model) == 'srof' .or. trim(rof_model) == 'drof') then
       rof_prognostic = .false.
    else
       rof_prognostic = .true.
    end if
    call NUOPC_CompAttributeGet(gcomp, name='ATM_model', value=atm_model, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (trim(atm_model) == 'satm' .or. trim(atm_model) == 'datm') then
       atm_prognostic = .false.
    else
       atm_prognostic = .true.
    end if
    call NUOPC_CompAttributeGet(gcomp, name='GLC_model', value=glc_model, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (trim(glc_model) == 'sglc') then
       glc_present = .false.
    else
       glc_present = .true.
    end if
    if (.not. glc_present) then
       cism_evolve = .false.
    else
       call NUOPC_CompAttributeGet(gcomp, name="cism_evolve", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (trim(cvalue) == '.true.') then
          cism_evolve = .true.
       else if (trim(cvalue) == '.false.') then
          cism_evolve = .false.
       else
          call shr_sys_abort(subname//'Could not determine cism_evolve value '//trim(cvalue))
       endif
    end if

    if (masterproc) then
       write(iulog,'(a   )')' atm component                 = '//trim(atm_model)
       write(iulog,'(a   )')' rof component                 = '//trim(rof_model)
       write(iulog,'(a   )')' glc component                 = '//trim(glc_model)
       write(iulog,'(a,L2)')' atm_prognostic                = ',atm_prognostic
       write(iulog,'(a,L2)')' rof_prognostic                = ',rof_prognostic
       write(iulog,'(a,L2)')' glc_present                   = ',glc_present
       if (glc_present) then
          write(iulog,'(a,L2)')' cism_evolve                  = ',cism_evolve
       end if
       write(iulog,'(a   )')' flds_scalar_name              = '//trim(flds_scalar_name)
       write(iulog,'(a,i8)')' flds_scalar_num               = ',flds_scalar_num
       write(iulog,'(a,i8)')' flds_scalar_index_nx          = ',flds_scalar_index_nx
       write(iulog,'(a,i8)')' flds_scalar_index_ny          = ',flds_scalar_index_ny
       write(iulog,'(a,i8)')' flds_scalar_index_nextsw_cday = ',flds_scalar_index_nextsw_cday
    end if

    !----------------------
    ! Set the namelist filename
    !----------------------
    call control_setNL("lnd_in"//trim(inst_suffix))


    call advertise_fields(gcomp, flds_scalar_name, glc_present, cism_evolve, rof_prognostic, atm_prognostic, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! reset shr logging to original values
    !----------------------------------------------------------------------------
    call shr_log_setLogUnit(shrlogunit)
    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine InitializeAdvertise

  !===============================================================================
  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    ! Realize the list of fields that will be exchanged
    !$  use omp_lib, only : omp_set_num_threads
    use ESMF                      , only : ESMF_VM, ESMF_VMGet
    use clm_instMod               , only : lnd2atm_inst, lnd2glc_inst, water_inst
    use domainMod                 , only : ldomain
    use decompMod                 , only : bounds_type, get_proc_bounds
    use lnd_set_decomp_and_domain , only : lnd_set_decomp_and_domain_from_readmesh
    use lnd_set_decomp_and_domain , only : lnd_set_mesh_for_single_column
    use lnd_set_decomp_and_domain , only : lnd_set_decomp_and_domain_for_single_column

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM)           :: vm                    ! Virtual machine, description of parallel procesors being used (both MPI and OpenMP)
    type(ESMF_Time)         :: currTime              ! Current time
    type(ESMF_Time)         :: startTime             ! Start time
    type(ESMF_Time)         :: refTime               ! Ref time
    type(ESMF_TimeInterval) :: timeStep              ! Model timestep
    type(ESMF_CalKind_Flag) :: esmf_caltype          ! esmf calendar type
    integer                 :: ref_ymd               ! reference date (YYYYMMDD)
    integer                 :: ref_tod               ! reference time of day (sec)
    integer                 :: yy,mm,dd              ! Temporaries for time query
    integer                 :: start_ymd             ! start date (YYYYMMDD)
    integer                 :: start_tod             ! start time of day (sec)
    integer                 :: curr_ymd              ! Start date (YYYYMMDD)
    integer                 :: curr_tod              ! Start time of day (sec)
    integer                 :: dtime_sync            ! coupling time-step from the input synchronization clock
    integer                 :: localPet              ! local PET (Persistent Execution Threads) (both MPI tasks and OpenMP threads)
    integer                 :: localPeCount          ! Number of local Processors
    character(len=CL)       :: starttype             ! start-type (startup, continue, branch, hybrid)
    character(len=CL)       :: calendar              ! calendar type name
    logical                 :: brnch_retain_casename ! flag if should retain the case name on a branch start type
    integer                 :: nsrest                ! ctsm restart type
    integer                 :: lbnum                 ! input to memory diagnostic
    integer                 :: shrlogunit            ! original log unit
    integer                 :: n, ni, nj             ! Indices
    character(len=CL)       :: cvalue                ! config data
    character(len=FL)       :: meshfile_mask         ! filename of mesh file with land mask
    character(len=CL)       :: ctitle                ! case description title
    character(len=CL)       :: caseid                ! case identifier name
    real(r8)                :: scol_lat              ! single-column latitude
    real(r8)                :: scol_lon              ! single-column longitude
    real(r8)                :: scol_area             ! single-column area
    real(r8)                :: scol_frac             ! single-column frac
    integer                 :: scol_mask             ! single-column mask
    real(r8)                :: scol_spval            ! single-column special value to indicate it isn't set
    character(len=FL)       :: single_column_lnd_domainfile   ! domain filename to use for single-column mode (i.e. SCAM)
    type(bounds_type)      :: bounds                          ! bounds
    type(ESMF_Field)        :: lfield                         ! Land field read in
    character(CL) ,pointer  :: lfieldnamelist(:) => null()    ! Land field namelist item sent with land field
    integer                 :: fieldCount                     ! Number of fields on export state
    integer                 :: rank                           ! Rank of field (1D or 2D)
    real(r8), pointer       :: fldptr1d(:)                    ! 1D field pointer
    real(r8), pointer       :: fldptr2d(:,:)                  ! 2D field pointer
    logical                 :: isPresent                      ! If attribute is present
    logical                 :: isSet                          ! If attribute is present and also set
    character(len=CL)       :: model_version                  ! Model version
    character(len=CL)       :: hostname                       ! hostname of machine running on
    character(len=CL)       :: username                       ! user running the model
    character(len=*),parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------------------------------------------------------------
    ! Single column logic - if mask is zero for nearest neighbor search then
    ! set all export state fields to zero and return
    !----------------------------------------------------------------------------

    ! If single_column is true - used single_column_domainfile to
    ! obtain nearest neighbor values for scol_lon and scol_lat
    ! If single_column is false and scol_lon and scol_lat are not equal to scol_spval then
    ! use scol_lon and scol_lat directly

    call NUOPC_CompAttributeGet(gcomp, name='scol_lon', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scol_lon
    call NUOPC_CompAttributeGet(gcomp, name='scol_lat', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scol_lat
    call NUOPC_CompAttributeGet(gcomp, name='single_column_lnd_domainfile', value=single_column_lnd_domainfile, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! TODO: there is a problem retrieving scol_spval from the driver - for now
    ! hard-wire scol_spval - this needs to be fixed
    scol_spval = -999._r8
    ! call NUOPC_CompAttributeGet(gcomp, name='scol_spval', value=cvalue, rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! read(cvalue,*) scol_spval
    scol_valid = .true.
    if (scol_lon > scol_spval .and. scol_lat > scol_spval) then
       single_column = (trim(single_column_lnd_domainfile) /= 'UNSET')

       call NUOPC_CompAttributeGet(gcomp, name='scol_lndmask', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) scol_mask

       call NUOPC_CompAttributeGet(gcomp, name='scol_lndfrac', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) scol_frac

       call lnd_set_mesh_for_single_column(scol_lon, scol_lat, mesh, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       scol_valid = (scol_mask == 1)
       if (.not. scol_valid) then
          write(iulog,'(a)')' single column mode point does not contain any land - will set all export data to 0'
          ! if single column is not valid - set all export state fields to zero and return
          call realize_fields(importState, exportState, mesh, flds_scalar_name, flds_scalar_num, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          allocate(lfieldnamelist(fieldCount))
          call ESMF_StateGet(exportState, itemNameList=lfieldnamelist, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          do n = 1, fieldCount
             if (trim(lfieldnamelist(n)) /= flds_scalar_name) then
                call ESMF_StateGet(exportState, itemName=trim(lfieldnamelist(n)), field=lfield, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                call ESMF_FieldGet(lfield, rank=rank, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                if (rank == 2) then
                   call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   fldptr2d(:,:) = 0._r8
                else
                   call ESMF_FieldGet(lfield, farrayPtr=fldptr1d, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   fldptr1d(:) = 0._r8
                end if
             end if
          enddo
          deallocate(lfieldnamelist)
          ! *******************
          ! *** RETURN HERE ***
          ! *******************
          RETURN
       else
          write(iulog,'(a,3(f10.5,2x))')' single column mode scol_lon/scol_lat/scol_frac is ',&
               scol_lon,scol_lat,scol_frac
       end if
    else
       single_column = .false.
    end if

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_log_getLogUnit (shrlogunit)
    call shr_log_setLogUnit (iulog)
#if (defined _MEMTRACE)
    if (masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_comp_nuopc_InitializeRealize:start::',lbnum)
    endif
#endif
    !----------------------------------------------------------------------------
    ! Initialize component threading
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, localPet=localPet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, pet=localPet, peCount=localPeCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if(localPeCount == 1) then
       call NUOPC_CompAttributeGet(gcomp, "nthreads", value=cvalue, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       read(cvalue,*) nthrds
    else
       nthrds = localPeCount
    endif

    !$  call omp_set_num_threads(nthrds)

    !----------------------
    ! Get properties from clock
    !----------------------
    call ESMF_ClockGet( clock, currTime=currTime, startTime=startTime, refTime=RefTime, &
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
    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) caseid
    ctitle= trim(caseid)
    call NUOPC_CompAttributeGet(gcomp, name='model_version', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) model_version

    ! Note that we assume that CTSM's internal dtime matches the coupling time step.
    ! i.e., we currently do NOT allow sub-cycling within a coupling time step.
    call set_timemgr_init(       &
         calendar_in=calendar,   &
         start_ymd_in=start_ymd, &
         start_tod_in=start_tod, &
         ref_ymd_in=ref_ymd,     &
         ref_tod_in=ref_tod,     &
         dtime_in=dtime_sync)

    ! Set model clock in lnd_comp_shr
    model_clock = clock

    call NUOPC_CompAttributeGet(gcomp, name="write_restart_at_endofrun", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) .eq. '.true.') write_restart_at_endofrun = .true.
    end if
    ! ---------------------
    ! Initialize first phase of ctsm
    ! ---------------------
    call NUOPC_CompAttributeGet(gcomp, name='hostname', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) hostname
    call NUOPC_CompAttributeGet(gcomp, name='username', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) username
    call NUOPC_CompAttributeGet(gcomp, name='brnch_retain_casename', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) brnch_retain_casename
    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) starttype

    if (     trim(starttype) == trim(startup_run)) then
       nsrest = nsrStartup
    else if (trim(starttype) == trim(continue_run)) then
       nsrest = nsrContinue
    else if (trim(starttype) == trim(branch_run)) then
       nsrest = nsrBranch
    else
       call shr_sys_abort( subname//' ERROR: unknown starttype' )
    end if

    ! set default values for run control variables
    call clm_varctl_set(&
         caseid_in=caseid, ctitle_in=ctitle,                     &
         brnch_retain_casename_in=brnch_retain_casename,         &
         single_column_in=single_column, scmlat_in=scol_lat, scmlon_in=scol_lon, &
         nsrest_in=nsrest, &
         version_in=model_version, &
         hostname_in=hostname, &
         username_in=username)

    call initialize1(dtime=dtime_sync)

    ! ---------------------
    ! Create ctsm decomp and domain info
    ! ---------------------
    if (scol_lon > scol_spval .and. scol_lat > scol_spval) then
       call lnd_set_decomp_and_domain_for_single_column(scol_lon, scol_lat, scol_mask, scol_frac)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call NUOPC_CompAttributeGet(gcomp, name='mesh_lnd', value=model_meshfile, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeGet(gcomp, name='mesh_mask', value=meshfile_mask, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call lnd_set_decomp_and_domain_from_readmesh(driver='cmeps', vm=vm, &
            meshfile_lnd=model_meshfile, meshfile_mask=meshfile_mask, mesh_ctsm=mesh, ni=ni, nj=nj, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! ---------------------
    ! Realize the actively coupled fields
    ! ---------------------
    call realize_fields(importState, exportState, mesh, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Finish initializing ctsm
    ! ---------------------
    call ESMF_ClockGet(clock, currTime=currtime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call initialize2(ni, nj, currtime)

    !--------------------------------
    ! Create land export state
    !--------------------------------
    call get_proc_bounds(bounds)
    call export_fields(gcomp, bounds, glc_present, rof_prognostic, &
         water_inst%waterlnd2atmbulk_inst, lnd2atm_inst, lnd2glc_inst, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

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

    call shr_log_setLogUnit (shrlogunit)

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

    !$  use omp_lib, only : omp_set_num_threads
    use ESMF        , only : ESMF_VM, ESMF_VMGet
    use clm_instMod , only : water_inst, atm2lnd_inst, glc2lnd_inst, lnd2atm_inst, lnd2glc_inst
    use decompMod   , only : bounds_type, get_proc_bounds
    use clm_driver  , only : clm_drv

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
    integer                :: localPet       ! local PET (Persistent Execution Threads) (both MPI tasks and OpenMP threads)
    integer                :: localPeCount   ! Number of local Processors
    logical                :: rstwr          ! .true. ==> write restart file before returning
    logical                :: nlend          ! .true. ==> last time-step
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

    !--------------------------------
    ! Single column logic if nearest neighbor point has a mask of zero
    !--------------------------------

    if (single_column .and. .not. scol_valid) then
       RETURN
    end if

    !$  call omp_set_num_threads(nthrds)

    !--------------------------------
    ! Reset share log units
    !--------------------------------

    call shr_log_getLogUnit (shrlogunit)
    call shr_log_setLogUnit (iulog)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_comp_nuopc_ModelAdvance:start::',lbnum)
    endif
#endif

    !--------------------------------
    ! Query the Component for its clock, importState and exportState and vm
    !--------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Determine time of next atmospheric shortwave calculation
    !--------------------------------

    call State_GetScalar(importState, &
         flds_scalar_index_nextsw_cday, nextsw_cday, &
         flds_scalar_name, flds_scalar_num, rc)

    ! Get proc bounds
    call get_proc_bounds(bounds)

    !--------------------------------
    ! Unpack import state
    !--------------------------------

    call t_startf ('lc_lnd_import')
    call import_fields( gcomp, bounds, glc_present, rof_prognostic, &
         atm2lnd_inst, glc2lnd_inst, water_inst%wateratm2lndbulk_inst, rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf ('lc_lnd_import')

    !--------------------------------
    ! Run model
    !--------------------------------

    dtime = get_step_size()

    ! TODO: This is currently hard-wired - is there a better way for nuopc?
    ! Note that the model clock is updated at the end of the time step not at the beginning
    nstep = get_nstep()

    !--------------------------------
    ! Determine doalb based on nextsw_cday sent from atm model
    !--------------------------------

    caldayp1 = get_curr_calday(offset=dtime, reuse_day_365_for_day_366=.true.)

    if (nstep == 1) then
       doalb = (abs(nextsw_cday- caldayp1) < 1.e-10_r8)
    else
       doalb = (nextsw_cday >= -0.5_r8)
    end if

    call update_rad_dtime(doalb)

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
    ! Determine if time to write restart
    !--------------------------------
    rstwr = .false.
    if (nlend .and. write_restart_at_endofrun) then
       rstwr = .true.
    else 
       call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (ESMF_AlarmIsCreated(alarm, rc=rc)) then
          if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             rstwr = .true.
             call ESMF_AlarmRingerOff( alarm, rc=rc )
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          endif
       endif
    end if

    !--------------------------------
    ! Run CTSM
    !--------------------------------

    ! call ESMF_VMBarrier(vm, rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call t_startf ('shr_orb_decl')
    ! Note - the orbital inquiries set the values in clm_varorb via the module use statements
    call  clm_orbital_update(clock, iulog, masterproc, eccen, obliqr, lambm0, mvelpp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    calday = get_curr_calday(reuse_day_365_for_day_366=.true.)
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

    call shr_log_setLogUnit (shrlogunit)

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
    if (.not. scol_valid) return

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
       call ESMF_LogWrite(subname//'setting alarms for ' // trim(name), ESMF_LOGMSG_INFO)

       !----------------------------------------------------------------------------------
       ! Stop alarm
       ! MUST be set before the restart alarm in case restarts happen at the stop alarm
       !----------------------------------------------------------------------------------
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

       !----------------------------------------------------------------------------------
       ! Restart alarm
       ! MUST be set after the stop alarm in case restarts happen at the stop alarm
       !----------------------------------------------------------------------------------
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

  subroutine CheckImport(gcomp, rc)
    type(ESMF_GridComp) :: gcomp
    integer, intent(out) :: rc
    character(len=*) , parameter :: subname = "("//__FILE__//":CheckImport)"

    ! This is the routine that enforces the explicit time dependence on the
    ! import fields. This simply means that the timestamps on the Fields in the
    ! importState are checked against the currentTime on the Component's
    ! internalClock. Consequenty, this model starts out with forcing fields
    ! at the current time as it does its forward step from currentTime to
    ! currentTime + timeStep.

    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_Time)               :: time
    type(ESMF_State)              :: importState
    logical                       :: allCurrent
    type(ESMF_Field), allocatable :: fieldList(:)
    integer                       :: i
    character(ESMF_MAXSTR)        :: fieldName
    character(ESMF_MAXSTR)        :: name
    character(ESMF_MAXSTR)        :: valueString

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    if (single_column .and. .not. scol_valid) then
       RETURN
    end if
    ! The remander of this should be equivalent to the NUOPC internal routine
    ! from NUOPC_ModeBase.F90

    ! query the component for info
    call NUOPC_CompGet(gcomp, name=name, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! query the Component for its clock and importState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! get the current time out of the clock
    call ESMF_ClockGet(clock, currTime=time, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! check that Fields in the importState show correct timestamp
    allCurrent = NUOPC_IsAtTime(importState, time, fieldList=fieldList, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (.not.allCurrent) then
      !TODO: introduce and use INCOMPATIBILITY return codes!!!!
      do i=1, size(fieldList)
        call ESMF_FieldGet(fieldList(i), name=fieldName, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return

        call NUOPC_GetAttribute(fieldList(i), name="StandardName", &
             value=valueString, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return

        call ESMF_LogWrite(trim(name)//": Field '"//trim(fieldName)//&
          "' in the importState is not at the expected time. StandardName: "&
          //trim(valueString), ESMF_LOGMSG_WARNING)
      enddo
      deallocate(fieldList)
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
        msg="NUOPC INCOMPATIBILITY DETECTED: Import Fields not at current time", &
        line=__LINE__, file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    endif

  end subroutine CheckImport
end module lnd_comp_nuopc
