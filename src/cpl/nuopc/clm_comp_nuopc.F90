module clm_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for CTSM
  !----------------------------------------------------------------------------

  use shr_kind_mod          , only : R8=>SHR_KIND_R8, CXX=>SHR_KIND_CXX, CL=>SHR_KIND_CL
  use shr_sys_mod           , only : shr_sys_abort
  use shr_file_mod          , only : shr_file_getlogunit, shr_file_setlogunit
  use shr_file_mod          , only : shr_file_getloglevel, shr_file_setloglevel
  use shr_file_mod          , only : shr_file_setIO, shr_file_getUnit
  use shr_string_mod        , only : shr_string_listGetNum
  use shr_orb_mod           , only : shr_orb_decl
  use shr_cal_mod           , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_ymd2date
  use esmFlds               , only : fldListFr, fldListTo, complnd, compname
  use esmFlds               , only : flds_scalar_name, flds_scalar_num
  use esmFlds               , only : flds_scalar_index_nx, flds_scalar_index_ny
  use esmFlds               , only : flds_scalar_index_nextsw_cday
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Realize
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Concat
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Getnumflds
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Getfldinfo
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_chkerr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_Clock_TimePrint
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_SetScalar
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_GetScalar
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_Diagnose
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_ArrayToState
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_StateToArray

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS        => SetServices,       &
    model_label_Advance     => label_Advance,     &
    model_label_SetRunClock => label_SetRunClock, &
    model_label_Finalize    => label_Finalize

  use lnd_import_export , only : lnd_import, lnd_export
  use clm_cpl_indices   , only : clm_cpl_indices_set
  use perf_mod          , only : t_startf, t_stopf, t_barrierf
  use spmdMod           , only : masterproc, mpicom, spmd_init
  use decompMod         , only : bounds_type, ldecomp, get_proc_bounds
  use domainMod         , only : ldomain
  use controlMod        , only : control_setNL
  use clm_varorb        , only : eccen, obliqr, lambm0, mvelpp
  use clm_varctl        , only : inst_index, inst_suffix, inst_name
  use clm_varctl        , only : single_column, clm_varctl_set, iulog
  use clm_varctl        , only : nsrStartup, nsrContinue, nsrBranch
  use clm_time_manager  , only : set_timemgr_init, advance_timestep
  use clm_time_manager  , only : set_nextsw_cday, update_rad_dtime
  use clm_time_manager  , only : get_nstep, get_step_size
  use clm_time_manager  , only : get_curr_date, get_curr_calday
  use clm_initializeMod , only : initialize1, initialize2
  use clm_initializeMod , only : atm2lnd_inst
  use clm_initializeMod , only : glc2lnd_inst
  use clm_initializeMod , only : lnd2atm_inst
  use clm_initializeMod , only : lnd2glc_inst
  use clm_driver        , only : clm_drv

  implicit none
  private ! except

  ! Module routines
  public  :: SetServices
  private :: InitializeP0
  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelSetRunClock
  private :: ModelAdvance
  private :: ModelFinalize

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------
 
  character(CXX)             :: flds_l2x = ''
  character(CXX)             :: flds_x2l = ''
  real(r8), allocatable      :: x2l(:,:)
  real(r8), allocatable      :: l2x(:,:)
  integer                    :: nflds_l2x
  integer                    :: nflds_x2l
  character(len=*),parameter :: grid_option = "mesh" ! grid_de, grid_arb, grid_reg, mesh
  integer, parameter         :: dbug = 10
  integer                    :: dbrc
  integer                    :: compid               ! component id

  !----- formats -----
  character(*),parameter :: modName =  "(clm_comp_nuopc)"
  character(*),parameter :: u_FILE_u = __FILE__

  !===============================================================================
  contains
  !===============================================================================

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
         specRoutine=ModelSetRunClock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine SetServices

  !===============================================================================

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries

    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv01p"/), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeP0

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM)          :: vm
    integer                :: lmpicom
    character(ESMF_MAXSTR) :: cvalue
    logical                :: exists
    integer                :: lsize      ! local array size
    integer                :: ierr       ! error code
    integer                :: shrlogunit ! original log unit
    integer                :: shrloglev  ! original log level
    integer                :: n,nflds       
    logical                :: isPresent
    character(len=512)     :: diro
    character(len=512)     :: logfile
    logical                :: activefld
    character(ESMF_MAXSTR) :: stdname, shortname
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    character(len=*), parameter :: format = "('("//trim(subname)//") :',A)"
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------------------------
    ! generate local mpi comm
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! initialize CLM MPI stuff
    !----------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) compid  ! convert from string to integer

    call mpi_comm_dup(lmpicom, mpicom, ierr)

    call spmd_init( mpicom, compid )

    !----------------------------------------------------------------------------
    ! determine instance information
    !----------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name="inst_name", value=inst_name, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name="inst_index", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) inst_index 

    call ESMF_AttributeGet(gcomp, name="inst_suffix", isPresent=isPresent, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name="inst_suffix", value=inst_suffix, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       inst_suffix = ''
    end if

    !----------------------------------------------------------------------------
    ! reset shr logging to my log file
    !----------------------------------------------------------------------------

    if (masterproc) then
       call NUOPC_CompAttributeGet(gcomp, name="diro", value=diro, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeGet(gcomp, name="logfile", value=logfile, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       if (len_trim(logfile) > 0) then
          iulog = shr_file_getUnit()
          open(iulog,file=trim(diro)//"/"//trim(logfile))
       else
          iulog = shrlogunit
       endif
    endif

    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

    !--------------------------------
    ! create import and export field list 
    !--------------------------------

    call shr_nuopc_fldList_Concat(fldListFr(complnd), fldListTo(complnd), flds_l2x, flds_x2l, flds_scalar_name)

    !--------------------------------
    ! advertise import and export fields
    !--------------------------------

    nflds = shr_nuopc_fldList_Getnumflds(fldListFr(complnd))
    do n = 1,nflds
       call shr_nuopc_fldList_Getfldinfo(fldListFr(complnd), n, activefld, stdname, shortname)
       if (activefld) then
          call NUOPC_Advertise(exportState, standardName=stdname, shortname=shortname, name=shortname, &
               TransferOfferGeomObject='will provide', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       call ESMF_LogWrite(subname//':Fr_'//trim(compname(complnd))//': '//trim(shortname), ESMF_LOGMSG_INFO)
    end do

    nflds = shr_nuopc_fldList_Getnumflds(fldListTo(complnd))
    do n = 1,nflds
       call shr_nuopc_fldList_Getfldinfo(fldListTo(complnd), n, activefld, stdname, shortname)
       if (activefld) then
          call NUOPC_Advertise(importState, standardName=stdname, shortname=shortname, name=shortname, &
               TransferOfferGeomObject='will provide', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       call ESMF_LogWrite(subname//':To_'//trim(compname(complnd))//': '//trim(shortname), ESMF_LOGMSG_INFO)
    end do

    !----------------------------------------------------------------------------
    ! reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogLevel(shrloglev)
    call shr_file_setLogUnit (shrlogunit)

  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Time)         :: currTime              ! Current time
    type(ESMF_Time)         :: startTime             ! Start time
    type(ESMF_Time)         :: stopTime              ! Stop time
    type(ESMF_Time)         :: refTime               ! Ref time
    type(ESMF_TimeInterval) :: timeStep              ! Model timestep
    type(ESMF_Calendar)     :: esmf_calendar         ! esmf calendar     
    type(ESMF_CalKind_Flag) :: esmf_caltype          ! esmf calendar type
    integer                 :: ref_ymd               ! reference date (YYYYMMDD)
    integer                 :: ref_tod               ! reference time of day (sec)
    integer                 :: yy,mm,dd              ! Temporaries for time query 
    integer                 :: start_ymd             ! start date (YYYYMMDD)
    integer                 :: start_tod             ! start time of day (sec)
    integer                 :: stop_ymd              ! stop date (YYYYMMDD)
    integer                 :: stop_tod              ! stop time of day (sec)
    integer                 :: curr_ymd              ! Start date (YYYYMMDD)
    integer                 :: curr_tod              ! Start time of day (sec)
    integer                 :: dtime_sync            ! coupling time-step from the input synchronization clock
    integer                 :: dtime_clm             ! clm time-step
    integer , allocatable   :: gindex(:)
    real(r8), pointer       :: lat(:)
    real(r8), pointer       :: lon(:)
    real(r8), pointer       :: elemCoords(:,:)
    real(r8), pointer       :: elemCornerCoords(:,:,:)
    real(r8)                :: dx,dy
    character(ESMF_MAXSTR)  :: cvalue
    character(ESMF_MAXSTR)  :: convCIM, purpComp
    type(ESMF_Grid)         :: Egrid
    type(ESMF_Mesh)         :: Emesh
    integer                 :: shrlogunit            ! original log unit
    integer                 :: shrloglev             ! original log level
    type(ESMF_VM)           :: vm
    integer                 :: n, m
    logical                 :: connected             ! is field connected?
    integer                 :: lsize                 ! local size ofarrays
    integer                 :: g,i,j                 ! indices
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
    integer                 :: nsrest                ! clm restart type
    logical                 :: brnch_retain_casename ! flag if should retain the case name on a branch start type
    integer                 :: lbnum                 ! input to memory diagnostic
    type(bounds_type)       :: bounds                ! bounds
    integer                 :: numg
    character(len=*),parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

    !----------------------
    ! Determine field indices of import/export arrays
    !----------------------

    call clm_cpl_indices_set(flds_x2l, flds_l2x)
    nflds_l2x = shr_string_listGetNum(flds_l2x)
    nflds_x2l = shr_string_listGetNum(flds_x2l)

#if (defined _MEMTRACE)
    if (masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','clm_comp_nuopc_InitializeRealize:start::',lbnum)
    endif
#endif

    !----------------------
    ! Obtain orbital values
    !----------------------

    ! Note - the orbital inquiries set the values in clm_varorb via the module use statements

    call NUOPC_CompAttributeGet(gcomp, name='orb_eccen', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) eccen

    call NUOPC_CompAttributeGet(gcomp, name='orb_obliqr', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) obliqr

    call NUOPC_CompAttributeGet(gcomp, name='orb_lambm0', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) lambm0

    call NUOPC_CompAttributeGet(gcomp, name='orb_mvelpp', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) mvelpp

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) caseid

    !TODO: case_desc does not appear in the esm_AddAttributes in esm.F90
    ! just hard-wire from now - is this even needed?
    ! call NUOPC_CompAttributeGet(gcomp, name='case_desc', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    ! read(cvalue,*) ctitle
    ctitle='UNSET'

    call NUOPC_CompAttributeGet(gcomp, name='scmlon', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlon

    call NUOPC_CompAttributeGet(gcomp, name='scmlat', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlat

    call NUOPC_CompAttributeGet(gcomp, name='single_column', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) single_column

    call NUOPC_CompAttributeGet(gcomp, name='brnch_retain_casename', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) brnch_retain_casename

    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) starttype

    call NUOPC_CompAttributeGet(gcomp, name='model_version', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) model_version

    call NUOPC_CompAttributeGet(gcomp, name='hostname', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) hostname

    call NUOPC_CompAttributeGet(gcomp, name='username', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
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
         currTime=currTime, startTime=startTime, stopTime=stopTime, refTime=RefTime, &
         timeStep=timeStep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=yy, mm=mm, dd=dd, s=curr_tod, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,curr_ymd)

    call ESMF_TimeGet( startTime, yy=yy, mm=mm, dd=dd, s=start_tod, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,start_ymd)

    call ESMF_TimeGet( stopTime, yy=yy, mm=mm, dd=dd, s=stop_tod, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,stop_ymd)

    call ESMF_TimeGet( refTime, yy=yy, mm=mm, dd=dd, s=ref_tod, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,ref_ymd)

    call ESMF_TimeGet( currTime, calkindflag=esmf_caltype, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar = shr_cal_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar = shr_cal_gregorian
    else
       call shr_sys_abort( subname//'ERROR:: bad calendar for ESMF' )
    end if

    call set_timemgr_init(       &
         calendar_in=calendar,   &
         start_ymd_in=start_ymd, &
         start_tod_in=start_tod, &
         ref_ymd_in=ref_ymd,     &
         ref_tod_in=ref_tod,     &
         stop_ymd_in=stop_ymd,   &
         stop_tod_in=stop_tod)

    call clm_varctl_set(&
         caseid_in=caseid, ctitle_in=ctitle,                     &
         brnch_retain_casename_in=brnch_retain_casename,         &
         single_column_in=single_column, scmlat_in=scmlat, scmlon_in=scmlon, &
         nsrest_in=nsrest, &
         version_in=model_version, &
         hostname_in=hostname, &
         username_in=username)

    !----------------------
    ! Read namelist, grid and surface data
    !----------------------

    call initialize1( )

    !----------------------
    ! Determine field arays
    !----------------------

    call get_proc_bounds( bounds )
    lsize = bounds%endg - bounds%begg + 1

    allocate(l2x(nflds_l2x,lsize))
    allocate(x2l(nflds_x2l,lsize))

    l2x(:,:)  = 0._r8
    x2l(:,:)  = 0._r8

    !--------------------------------
    ! generate the mesh
    ! grid_option specifies grid or mesh
    !--------------------------------

    allocate(gindex(lsize))
    allocate(elemCoords(2,lsize))          ! (lon+lat) * n_gridcells
    allocate(elemCornerCoords(2,4,lsize))  ! (lon+lat) * n_corners * n_gridcells

    do g = bounds%begg,bounds%endg
       n = 1 + (g - bounds%begg)
       gindex(n) = ldecomp%gdc2glo(g)
       elemCoords(1,n) = ldomain%lonc(g)
       elemCoords(2,n) = ldomain%latc(g)
       write(6,*)'CLM: g,lon,lat= ',g,ldomain%lonc(g),ldomain%latc(g)
       ! TBD
       ! tcraig, clm does not define corner values and there is no info about grid sizes (ie. dx, dy)
       ! anywhere so make something up for now.  this has to be fixed if weights are generated on the fly!
       ! someone from clm has to define the corner lon and lat for use here.
       ! corners are defined counterclockwise
       do m = 1,4
          if (m == 1 .or. m == 4) dx = -0.05
          if (m == 2 .or. m == 3) dx =  0.05
          if (m == 1 .or. m == 2) dy = -0.05
          if (m == 3 .or. m == 4) dy =  0.05
          elemCornerCoords(1,m,n) = ldomain%lonc(g) + dx
          elemCornerCoords(2,m,n) = ldomain%latc(g) + dy
       enddo
    end do

    Emesh = ESMF_MeshCreate(parametricDim=2, &
       coordSys=ESMF_COORDSYS_SPH_DEG, &
       elementIds=gindex, &
       elementType=ESMF_MESHELEMTYPE_QUAD, &
       elementCoords=elemCoords, &
       elementCornerCoords=elemCornerCoords, &
       rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    deallocate(gindex)
    deallocate(elemCoords)
    deallocate(elemCornerCoords)

    !--------------------------------
    ! realize the actively coupled fields
    !--------------------------------

    call shr_nuopc_fldList_Realize(importState, fldListTo(complnd), flds_scalar_name, flds_scalar_num, &
         mesh=Emesh, tag=subname//':ciceImport', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_fldList_Realize(exportState, fldListFr(complnd), flds_scalar_name, flds_scalar_num, &
         mesh=Emesh, tag=subname//':ciceExport', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Finish initializing clm
    !--------------------------------

    call initialize2()

    !--------------------------------
    ! Check that clm internal dtime aligns with clm coupling interval
    !--------------------------------

    call ESMF_ClockGet( clock, timeStep=timeStep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeIntervalGet( timeStep, s=dtime_sync, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    dtime_clm = get_step_size()

    if (masterproc) then
       write(iulog,*)'dtime_sync= ',dtime_sync,' dtime_clm= ',dtime_clm,' mod = ',mod(dtime_sync,dtime_clm)
    end if
    if (mod(dtime_sync,dtime_clm) /= 0) then
       write(iulog,*)'clm dtime ',dtime_clm,' and clock dtime ',dtime_sync,' never align'
       rc = ESMF_FAILURE
       return
    end if

    !--------------------------------
    ! Create land export state
    !--------------------------------

    call lnd_export(bounds, lnd2atm_inst, lnd2glc_inst, l2x, flds_l2x)

    !--------------------------------
    ! Get calendar day of nextsw calculation
    !--------------------------------

    ! TODO: for nuopc at startup will set nextsw_cday to the currtime - and this will simplify the whole
    ! initialization dependence with the atmosphere - does this make sense?

    if (nsrest == nsrStartup) then
       call ESMF_ClockGet( clock, currTime=currTime, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_TimeGet( currTime, dayOfYear_r8=nextsw_cday, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call set_nextsw_cday(nextsw_cday)

    !--------------------------------
    ! Pack export state
    ! Copy from l2x to exportState
    ! Set the coupling scalars
    !--------------------------------

    call shr_nuopc_grid_ArrayToState(l2x, flds_l2x, exportState, grid_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(dble(ldomain%ni), flds_scalar_index_nx, exportState, mpicom, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(dble(ldomain%nj), flds_scalar_index_ny, exportState, mpicom, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       call shr_nuopc_methods_State_diagnose(exportState,subname//':ES',rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

#ifdef USE_ESMF_METADATA
    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"
    call ESMF_AttributeAdd(comp, convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ShortName", "CLM", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", "Community Land Model", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Description", "Community Land Model", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2017", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "Terrestrial", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Name", "TBD", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "EmailAddress", TBD, convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", convention=convCIM, purpose=purpComp, rc=rc)
#endif

    call shr_file_setLogLevel(shrloglev)
    call shr_file_setLogUnit (shrlogunit)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

#if (defined _MEMTRACE)
    if(masterproc) then
       write(iulog,*) TRIM(Sub) // ':end::'
       lbnum=1
       call memmon_dump_fort('memmon.out','clm_comp_nuopc_InitializeRealize:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

  end subroutine InitializeRealize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Run CTSM

    ! local variables:
    type(ESMF_Clock)       :: clock
    type(ESMF_Alarm)       :: alarm
    type(ESMF_Time)        :: currTime
    type(ESMF_State)       :: importState, exportState
    character(ESMF_MAXSTR) :: cvalue
    integer                :: shrlogunit     ! original log unit
    integer                :: shrloglev      ! original log level
    character(ESMF_MAXSTR) :: case_name      ! case name
    integer                :: ymd            ! CLM current date (YYYYMMDD)
    integer                :: yr             ! CLM current year
    integer                :: mon            ! CLM current month
    integer                :: day            ! CLM current day
    integer                :: tod            ! CLM current time of day (sec)
    integer                :: ymd_sync       ! Sync date (YYYYMMDD)
    integer                :: yr_sync        ! Sync current year
    integer                :: mon_sync       ! Sync current month
    integer                :: day_sync       ! Sync current day
    integer                :: tod_sync       ! Sync current time of day (sec)
    integer                :: dtime          ! time step increment (sec)
    integer                :: nstep          ! time step index
    logical                :: rstwr          ! .true. ==> write restart file before returning
    logical                :: nlend          ! .true. ==> last time-step
    logical                :: dosend         ! true => send data back to driver
    logical                :: doalb          ! .true. ==> do albedo calculation on this time step
    logical                :: rof_prognostic ! .true. => running with a prognostic ROF model
    logical                :: glc_present    ! .true. => running with a non-stub GLC model
    real(r8)               :: nextsw_cday    ! calday from clock of next radiation computation
    real(r8)               :: caldayp1       ! clm calday plus dtime offset
    integer                :: lbnum          ! input to memory diagnostic
    integer                :: g,i            ! counters
    real(r8)               :: calday         ! calendar day for nstep
    real(r8)               :: declin         ! solar declination angle in radians for nstep
    real(r8)               :: declinp1       ! solar declination angle in radians for nstep+1
    real(r8)               :: eccf           ! earth orbit eccentricity factor
    type(bounds_type)      :: bounds         ! bounds
    character(len=32)      :: rdate          ! date char string for restart file names
    character(len=*),parameter  :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogLevel(max(shrloglev,1))
    call shr_file_setLogUnit (iulog)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','clm_comp_nuopc_ModelAdvance:start::',lbnum)
    endif
#endif

    !--------------------------------
    ! Query the Component for its clock, importState and exportState
    !--------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 1) then
      call shr_nuopc_methods_Clock_TimePrint(clock,subname//'clock',rc=rc)
    endif

    !--------------------------------
    ! Determine time of next atmospheric shortwave calculation
    !--------------------------------

    call shr_nuopc_methods_State_GetScalar(importState, &
         flds_scalar_index_nextsw_cday, nextsw_cday, mpicom, &
         flds_scalar_name, flds_scalar_num, rc)
    call set_nextsw_cday( nextsw_cday )

    !----------------------
    ! Obtain orbital values
    !----------------------

    call NUOPC_CompAttributeGet(gcomp, name='orb_eccen', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) eccen

    call NUOPC_CompAttributeGet(gcomp, name='orb_obliqr', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) obliqr

    call NUOPC_CompAttributeGet(gcomp, name='orb_lambm0', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) lambm0

    call NUOPC_CompAttributeGet(gcomp, name='orb_mvelpp', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) mvelpp

    ! TODO: Determine if we're running with prognostic ROF and GLC models. 
    ! This won't change throughout the run, but we can't count on this being set in 
    ! initialization, so need to get it in the run method.
    ! This seems to imply an order dependence - for now simply set this to false for debugging - since
    ! are running without a river model for the first test

    rof_prognostic = .false. !THIS MUST BE FIXED before a prognostic river model is added with a nuopc cap
    glc_present    = .false. !THIS MUST BE FIXED before a prognostic ice-sheet model is added with a nuopc cap

    !--------------------------------
    ! Unpack export state
    !--------------------------------

    call get_proc_bounds(bounds)

    call t_startf ('lc_lnd_import')

    call shr_nuopc_grid_StateToArray(importState, x2l, flds_x2l, grid_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call lnd_import( bounds, &
         x2l=x2l, &
         flds_x2l=flds_x2l, &
         glc_present=glc_present, &
         atm2lnd_inst=atm2lnd_inst, &
         glc2lnd_inst=glc2lnd_inst )

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

       call ESMF_ClockGetAlarm(clock, alarmname='seq_timemgr_alarm_restart', alarm=alarm, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          rstwr = .true.
          call ESMF_AlarmRingerOff( alarm, rc=rc )
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          rstwr = .false.
       endif

       !--------------------------------
       ! Determine if time to stop
       !--------------------------------

       call ESMF_ClockGetAlarm(clock, alarmname='seq_timemgr_alarm_stop', alarm=alarm, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          nlend = .true.
          call ESMF_AlarmRingerOff( alarm, rc=rc )
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          nlend = .false.
       endif

       !--------------------------------
       ! Run clm
       !--------------------------------

       call t_barrierf('sync_clm_run1', mpicom)

       call t_startf ('shr_orb_decl')
       calday = get_curr_calday()
       call shr_orb_decl( calday     , eccen, mvelpp, lambm0, obliqr, declin  , eccf )
       call shr_orb_decl( nextsw_cday, eccen, mvelpp, lambm0, obliqr, declinp1, eccf )
       call t_stopf ('shr_orb_decl')

       call t_startf ('clm_run')

       call ESMF_ClockGet( clock, currTime=currTime )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet( currTime, yy=yr_sync, mm=mon_sync, dd=day_sync, s=tod_sync, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr_sync, mon_sync, day_sync, tod_sync

       call clm_drv(doalb, nextsw_cday, declinp1, declin, rstwr, nlend, rdate, rof_prognostic)

       call t_stopf ('clm_run')

       !--------------------------------
       ! Pack export state
       !--------------------------------

       call t_startf ('lc_lnd_export')

       ! create array l2x
       call lnd_export(bounds, lnd2atm_inst, lnd2glc_inst, l2x, flds_l2x)

       ! map l2x to export state
       call shr_nuopc_grid_ArrayToState(l2x, flds_l2x, exportState, grid_option, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call t_stopf ('lc_lnd_export')

       !--------------------------------
       ! Advance clm time step
       !--------------------------------

       call t_startf ('lc_clm2_adv_timestep')
       call advance_timestep()
       call t_stopf ('lc_clm2_adv_timestep')

    end do

    ! Check that internal clock is in sync with master clock
    ! Note that the driver clock has not been updated yet - so at this point
    ! CLM is actually 1 coupling intervals ahead of the driver clock

    call get_curr_date( yr, mon, day, tod, offset=-2*dtime )
    ymd = yr*10000 + mon*100 + day
    tod = tod

    call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=yr_sync, mm=mon_sync, dd=day_sync, s=tod_sync, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yr_sync, mon_sync, day_sync, ymd_sync)

    if ( (ymd /= ymd_sync) .and. (tod /= tod_sync) ) then
       write(iulog,*)' clm ymd=',ymd     ,'  clm tod= ',tod
       write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       rc = ESMF_FAILURE
       call ESMF_LogWrite(subname//" CLM clock not in sync with Master Sync clock",ESMF_LOGMSG_ERROR, rc=dbrc)
    end if

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       call shr_nuopc_methods_State_diagnose(exportState,subname//':ES',rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_ClockPrint(clock, options="currTime", preString="------>Advancing LND from: ", rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_ClockPrint(clock, options="stopTime", preString="--------------------------------> to: ", rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    !--------------------------------
    ! Reset shr logging to my original values
    !--------------------------------

    call shr_file_setLogLevel(shrloglev)
    call shr_file_setLogUnit (shrlogunit)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','clm_comp_nuopc_ModelAdvance:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)
    implicit none
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)         :: mclock, dclock
    type(ESMF_Time)          :: mcurrtime, dcurrtime
    type(ESMF_Time)          :: mstoptime
    type(ESMF_TimeInterval)  :: mtimestep, dtimestep
    character(len=128)       :: mtimestring, dtimestring
    type(ESMF_Alarm),pointer :: alarmList(:)
    type(ESMF_Alarm)         :: dalarm
    integer                  :: alarmcount, n
    character(len=*),parameter :: subname=trim(modName)//':(ModelSetRunClock) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! check that the current time in the model and driver are the same
    !--------------------------------

    if (mcurrtime /= dcurrtime) then
      call ESMF_TimeGet(dcurrtime, timeString=dtimestring, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_TimeGet(mcurrtime, timeString=mtimestring, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_LogWrite(subname//" ERROR in time consistency; "//trim(dtimestring)//" ne "//trim(mtimestring),  &
           ESMF_LOGMSG_ERROR, rc=dbrc)
      rc=ESMF_Failure
      return
    endif

    !--------------------------------                                                                                 
    ! force model clock currtime and timestep to match driver and set stoptime                                        
    !--------------------------------                                                                                 

    mstoptime = mcurrtime + dtimestep

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! copy alarms from driver to model clock if model clock has no alarms (do this only once!)
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then
      call ESMF_ClockGetAlarmList(dclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      allocate(alarmList(alarmCount))
      call ESMF_ClockGetAlarmList(dclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmList=alarmList, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      do n = 1, alarmCount
         ! call ESMF_AlarmPrint(alarmList(n), rc=rc)
         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
         dalarm = ESMF_AlarmCreate(alarmList(n), rc=rc)
         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
         call ESMF_AlarmSet(dalarm, clock=mclock, rc=rc)
         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      enddo

      deallocate(alarmList)
    endif

    !--------------------------------                                                                                 
    ! Advance model clock to trigger alarms then reset model clock back to currtime                                   
    !--------------------------------                                                                                 

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ModelSetRunClock

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(*), parameter :: F00   = "('(clm_comp_nuopc) ',8a)"
    character(*), parameter :: F91   = "('(clm_comp_nuopc) ',73('-'))"
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    !-------------------------------------------------------------------------------

    !--------------------------------
    ! Finalize routine
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    if (masterproc) then
       write(iulog,F91)
       write(iulog,F00) 'CLM: end of main integration loop'
       write(iulog,F91)
    end if

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ModelFinalize

  !===============================================================================

end module clm_comp_nuopc
