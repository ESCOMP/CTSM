module lilac_mod

  !-----------------------------------------------------------------------
  ! This is the driver for running CTSM, the ESMF lilac atm cap, and
  ! optionally the MOSART river model that is put in place to ensure
  ! that the host atmosphere does not need to know about ESMF
  !-----------------------------------------------------------------------

  ! External libraries
  use ESMF

  ! shr code routines
  use shr_sys_mod   , only : shr_sys_abort
  use shr_kind_mod  , only : r8 => shr_kind_r8

  ! lilac routines and data
  use lilac_io      , only : lilac_io_init
  use lilac_time    , only : lilac_time_clockinit, lilac_time_alarminit
  use lilac_time    , only : lilac_time_restart_write, lilac_time_restart_read
  use lilac_atmaero , only : lilac_atmaero_init, lilac_atmaero_interp
  use lilac_atmcap  , only : lilac_atmcap_init_vars
  use lilac_history , only : lilac_history_init
  use lilac_history , only : lilac_history_write
  use lilac_methods , only : chkerr
  use lilac_constants, only : logunit
  use lilac_driver_pio_mod, only : lilac_driver_pio_init1, lilac_driver_pio_init2
  use ctsm_LilacCouplingFields, only : create_a2l_field_list, create_l2a_field_list
  use ctsm_LilacCouplingFields, only : complete_a2l_field_list, complete_l2a_field_list
  use ctsm_LilacCouplingFields, only : a2l_fields

  ! lilac register phaes
  use lilac_atmcap  , only : lilac_atmcap_register
  use lilac_cpl     , only : cpl_atm2lnd_register, cpl_lnd2atm_register
  use lilac_cpl     , only : cpl_lnd2rof_register, cpl_rof2lnd_register

  ! ctsm register
  use lnd_comp_esmf , only : lnd_register ! ctsm routine

  ! mosart register
  use rof_comp_esmf , only : rof_register ! mosart routine

  implicit none

  public :: lilac_init1
  public :: lilac_init2
  public :: lilac_run
  public :: lilac_final

  ! Gridded components and states in gridded components
  type(ESMF_GridComp)        :: atm_gcomp
  type(ESMF_GridComp)        :: lnd_gcomp
  type(ESMF_GridComp)        :: rof_gcomp

  ! Coupler components
  type(ESMF_CplComp)         :: cpl_atm2lnd_comp
  type(ESMF_CplComp)         :: cpl_lnd2atm_comp
  type(ESMF_CplComp)         :: cpl_lnd2rof_comp
  type(ESMF_CplComp)         :: cpl_rof2lnd_comp

  ! States
  type(ESMF_State)           :: atm2cpl_state, cpl2atm_state ! on atm mesh (1 field bundle)
  type(ESMF_State)           :: lnd2cpl_state, cpl2lnd_state ! on lnd mesh (2 field bundles)
  type(ESMF_State)           :: rof2cpl_state, cpl2rof_state ! on rof mesh (1 field bundle)

  ! Clock, TimeInterval, and Times
  type(ESMF_Clock)           :: lilac_clock
  type(ESMF_Calendar),target :: lilac_calendar
  type(ESMF_Alarm)           :: lilac_restart_alarm
  type(ESMF_Alarm)           :: lilac_stop_alarm

  ! Coupling to mosart is now set to .false. by default
  logical                    :: couple_to_river = .false.

  integer                    :: mytask
  character(ESMF_MAXSTR)     :: starttype

  character(*) , parameter   :: modname     = "lilac_mod"
  character(*), parameter    :: u_FILE_u = &
       __FILE__

!========================================================================
contains
!========================================================================

  subroutine lilac_init1()

    ! --------------------------------------------------------------------------------
    ! This is called by the host atmosphere. This is phase 1 of the lilac initialization.
    !
    ! Indices defined in lilac_coupling_fields (lilac_a2l_* and lilac_l2a_*) are not
    ! valid until this is called.
    ! --------------------------------------------------------------------------------

    call create_a2l_field_list()
    call create_l2a_field_list()

  end subroutine lilac_init1


  subroutine lilac_init2(mpicom, atm_global_index, atm_lons, atm_lats, &
       atm_global_nx, atm_global_ny, atm_calendar, atm_timestep, &
       atm_start_year, atm_start_mon, atm_start_day, atm_start_secs, &
       starttype_in, fields_needed_from_data)

    ! --------------------------------------------------------------------------------
    ! This is called by the host atmosphere. This is phase 2 of the lilac initialization.
    ! --------------------------------------------------------------------------------

    ! input/output variables
    integer          , intent(inout) :: mpicom  ! input commiunicator from atm
    integer          , intent(in)    :: atm_global_index(:)
    real(r8)         , intent(in)    :: atm_lons(:)
    real(r8)         , intent(in)    :: atm_lats(:)
    integer          , intent(in)    :: atm_global_nx
    integer          , intent(in)    :: atm_global_ny
    character(len=*) , intent(in)    :: atm_calendar
    integer          , intent(in)    :: atm_timestep
    integer          , intent(in)    :: atm_start_year !(yyyy)
    integer          , intent(in)    :: atm_start_mon  !(mm)
    integer          , intent(in)    :: atm_start_day
    integer          , intent(in)    :: atm_start_secs
    character(len=*) , intent(in)    :: starttype_in

    ! List of field indices that need to be read from data, because the host atmosphere
    ! isn't going to provide them. These should be indices given in
    ! ctsm_LilacCouplingFields (lilac_a2l_Faxa_bcphidry). This can be an empty list if no
    ! fields need to be read from data.
    integer          , intent(in)    :: fields_needed_from_data(:)

    ! local variables
    character(ESMF_MAXSTR)      :: caseid
    logical                     :: create_esmf_pet_files
    type(ESMF_LogKind_Flag)     :: logkindflag
    type(ESMF_TimeInterval)     :: timeStep
    type(ESMF_Time)             :: startTime
    integer                     :: yy,mm,dd,sec
    integer                     :: lsize
    type(ESMF_State)            :: importState, exportState
    type(ESMF_VM)               :: vm
    integer                     :: user_rc, rc
    character(len=ESMF_MAXSTR)  :: cname   !components or cpl names
    integer                     :: ierr
    integer                     :: n, i
    integer                     :: fileunit
    integer, parameter          :: debug = 1   !-- internal debug level
    character(len=*), parameter :: subname=trim(modname)//': [lilac_init] '

    ! initialization of pio
    integer           :: compids(1) = (/1/)         ! for pio_init2 - array with component ids
    character(len=32) :: compLabels(1) = (/'LND'/)  ! for pio_init2
    character(len=64) :: comp_name(1) = (/'LND'/)   ! for pio_init2
    logical           :: comp_iamin(1) = (/.true./) ! for pio init2
    !------------------------------------------------------------------------

    namelist /lilac_run_input/ caseid, create_esmf_pet_files

    ! Initialize return code
    rc = ESMF_SUCCESS

    !-------------------------------------------------------------------------
    ! Set module variable starttype
    !-------------------------------------------------------------------------
    starttype = starttype_in

    ! ------------------------------------------------------------------------
    ! Read main namelist
    ! ------------------------------------------------------------------------

    ! Initialize variables in case not set in namelist (but we expect them to be set)
    caseid = 'UNSET'
    create_esmf_pet_files = .false.

    open(newunit=fileunit, status="old", file="lilac_in")
    read(fileunit, lilac_run_input, iostat=ierr)
    if (ierr > 0) then
       call shr_sys_abort(trim(subname) // 'error reading in lilac_run_input')
    end if
    close(fileunit)

    if (create_esmf_pet_files) then
       logkindflag = ESMF_LOGKIND_MULTI
    else
       logkindflag = ESMF_LOGKIND_MULTI_ON_ERROR
    end if

    ! ------------------------------------------------------------------------
    ! Complete setup of field lists started in lilac_init1, now that we know the number
    ! of atm points.
    ! ------------------------------------------------------------------------
    call complete_a2l_field_list(size(atm_global_index), fields_needed_from_data)
    call complete_l2a_field_list(size(atm_global_index))

    !-------------------------------------------------------------------------
    ! Initialize pio with first initialization
    ! AFTER call to MPI_init (which is in the host atm driver) and
    ! BEFORE call to ESMF_Initialize
    !-------------------------------------------------------------------------
    call lilac_driver_pio_init1(ncomps=1, nlfilename="lilac_in", Global_Comm=mpicom)

    !-------------------------------------------------------------------------
    ! Initialize ESMF, set the default calendar and log type.
    !-------------------------------------------------------------------------

    ! NOTE: the default calendar is set to GREGORIAN and is reset below in the initialization of
    ! the lilac clock
    call ESMF_Initialize(mpiCommunicator=mpicom, defaultCalKind=ESMF_CALKIND_GREGORIAN, &
         logkindflag=logkindflag, logappendflag=.false., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogSet(flush=.true.)
    call ESMF_LogWrite(subname//".........................", ESMF_LOGMSG_INFO)
    call ESMF_LogWrite(subname//"Initializing ESMF        ", ESMF_LOGMSG_INFO)

    call ESMF_VMGetGlobal(vm=vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=mytask,  rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------------------------------------------------
    ! Initialize PIO with second initialization
    !-------------------------------------------------------------------------
    call lilac_driver_pio_init2(compids, compLabels, comp_iamin, (/mpicom/), (/mytask/))
    call ESMF_LogWrite(subname//"initialized lilac_driver_pio_init2 ...", ESMF_LOGMSG_INFO)

    !-------------------------------------------------------------------------
    ! Initial lilac atmosphere cap module variables
    !-------------------------------------------------------------------------
    ! This must be done BEFORE the atmcap initialization
    call lilac_atmcap_init_vars(atm_global_index, atm_lons, atm_lats, atm_global_nx, atm_global_ny)

    !-------------------------------------------------------------------------
    ! Create Gridded and Coupler Components
    !-------------------------------------------------------------------------
    cname = " LILAC atm cap "
    atm_gcomp = ESMF_GridCompCreate(name=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error lilac atmcap initialization')
    call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       write(logunit,*) trim(subname) // "lilac atm cap gridded component created"
    end if

    cname = " CTSM "
    lnd_gcomp = ESMF_GridCompCreate(name=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error lilac ctsm initialization')
    call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       write(logunit,*) trim(subname) // " ctsm gridded component created"
    end if

    if (couple_to_river) then
       cname = " MOSART "
       rof_gcomp = ESMF_GridCompCreate(name=cname, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error lilac mosart initialization')
       call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
       if (mytask == 0) then
          write(logunit,*) trim(subname) // " mosart gridded component created"
       end if
    end if

    cname = "Coupler from atmosphere to land"
    cpl_atm2lnd_comp = ESMF_CplCompCreate(name=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error lilac cpl_a2l initialization')
    call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       write(logunit,*) trim(subname) // " coupler component (atmosphere to land) created"
    end if

    cname = "Coupler from land to atmosphere"
    cpl_lnd2atm_comp = ESMF_CplCompCreate(name=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error lilac cpl_l2a initialization')
    call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       write(logunit,*) trim(subname) // " coupler component (land to atmosphere) created"
    end if

    if (couple_to_river) then
       cname = "Coupler from river to land"
       cpl_rof2lnd_comp = ESMF_CplCompCreate(name=cname, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error lilac cpl_r2l initialization')
       call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
       if (mytask == 0) then
          write(logunit,*) trim(subname) // " coupler component (river to land) created"
       end if

       cname = "Coupler from land to river"
       cpl_lnd2rof_comp = ESMF_CplCompCreate(name=cname, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error lilac cpl_l2r initialization')
       call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
       if (mytask == 0) then
          write(logunit,*) trim(subname) // " coupler component (land to river) created"
       end if
    end if

    !-------------------------------------------------------------------------
    ! Register gridded and coupler components
    !-------------------------------------------------------------------------

    ! Register section -- set services -- atmcap
    call ESMF_GridCompSetServices(atm_gcomp, userRoutine=lilac_atmcap_register, userRc=user_rc, rc=rc)
    if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort('atm_gcomp register failure')
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('atm_gcomp register failure')
    call ESMF_LogWrite(subname//"  atmos SetServices finished!", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       write(logunit,*) trim(subname) // " lilac atm cap setservices finished"
    end if

    ! Register section -- set services -- ctsm
    call ESMF_GridCompSetServices(lnd_gcomp, userRoutine=lnd_register, userRc=user_rc, rc=rc)
    if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort('lnd_gcomp register failure')
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('lnd_gcomp register failure')
    call ESMF_LogWrite(subname//"CSTM SetServices finished!", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       write(logunit,*) trim(subname) // " CTSM setservices finished"
    end if

    if (couple_to_river) then
       ! Register section -- set services -- mosart
       call ESMF_GridCompSetServices(rof_gcomp, userRoutine=rof_register, userRc=user_rc, rc=rc)
       if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort('rof_gcomp register failure')
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('rof_gcomp register failure')
       call ESMF_LogWrite(subname//"MOSART SetServices finished!", ESMF_LOGMSG_INFO)
       if (mytask == 0) then
          write(logunit,*) trim(subname) // " CTSM setservices finished"
       end if
    end if

    ! Register section -- set services -- coupler atmosphere to land
    call ESMF_CplCompSetServices(cpl_atm2lnd_comp, userRoutine=cpl_atm2lnd_register, userRc=user_rc, rc=rc)
    if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort('cpl_atm2lnd_comp register failure')
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('cpl_atm2lnd_comp register failure')
    call ESMF_LogWrite(subname//"Coupler from atmosphere to land  SetServices finished!", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       write(logunit,*) trim(subname) // " coupler from atmosphere to land setservices finished"
    end if

    if (couple_to_river) then
       ! Register section -- set services -- river to land
       call ESMF_CplCompSetServices(cpl_rof2lnd_comp, userRoutine=cpl_rof2lnd_register, userRc=user_rc, rc=rc)
       if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort('cpl_rof2lnd_comp register failure')
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('cpl_rof2lnd_comp register failure')
       call ESMF_LogWrite(subname//"Coupler from river to land  SetServices finished!", ESMF_LOGMSG_INFO)
       if (mytask == 0) then
          write(logunit,*) trim(subname) // " coupler from river to land setservices finished"
       end if
    end if

    ! Register section -- set services -- coupler land to atmosphere
    call ESMF_CplCompSetServices(cpl_lnd2atm_comp, userRoutine=cpl_lnd2atm_register, userRc=user_rc, rc=rc)
    if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort('cpl_lnd2atm_comp register failure')
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('cpl_lnd2atm_comp register failure')
    call ESMF_LogWrite(subname//"Coupler from land to atmosphere SetServices finished!", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       write(logunit,*) trim(subname) // " coupler from land to atmosphere setservices finished"
    end if

    if (couple_to_river) then
       ! Register section -- set services -- coupler land to river
       call ESMF_CplCompSetServices(cpl_lnd2rof_comp, userRoutine=cpl_lnd2rof_register, userRc=user_rc, rc=rc)
       if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort('cpl_lnd2rof_comp register failure')
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('cpl_lnd2rof_comp register failure')
       call ESMF_LogWrite(subname//"Coupler from land to river SetServices finished!", ESMF_LOGMSG_INFO)
       if (mytask == 0) then
          write(logunit,*) trim(subname) // " coupler from land to river setservices finished"
       end if
    end if

    !-------------------------------------------------------------------------
    !  Create and initialize the lilac_clock, alarms and calendar
    !-------------------------------------------------------------------------

    call lilac_time_clockInit(caseid, starttype, atm_calendar, atm_timestep, &
       atm_start_year, atm_start_mon, atm_start_day, atm_start_secs, &
       lilac_clock, rc)

    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing clock")
    call ESMF_LogWrite(subname//"lilac_clock initialized", ESMF_LOGMSG_INFO)

    ! -------------------------------------------------------------------------
    ! Initialize LILAC gridded components
    ! First Create the empty import and export states used to pass data
    ! between components. (these are module variables)
    ! -------------------------------------------------------------------------

    ! Create import and export states for atm_gcomp
    atm2cpl_state = ESMF_StateCreate(name='state_from_atm', stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    cpl2atm_state = ESMF_StateCreate(name='state_to_atm', stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Initialze lilac_atm gridded component
    call ESMF_GridCompInitialize(atm_gcomp, importState=cpl2atm_state, exportState=atm2cpl_state, &
         clock=lilac_clock, userRc=user_rc, rc=rc)
    if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing atmcap")
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing atmcap")
    call ESMF_LogWrite(subname//"lilac_atm gridded component initialized", ESMF_LOGMSG_INFO)

    ! Create import and export states for lnd_gcomp (i.e. CTSM)
    cpl2lnd_state = ESMF_StateCreate(name='state_to_land', stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lnd2cpl_state = ESMF_StateCreate(name='state_fr_land', stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Add caseid and starttype as attributes of cpl2lnd_state
    call ESMF_AttributeSet(cpl2lnd_state, name="caseid", value=trim(caseid), rc=rc) 
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_AttributeSet(cpl2lnd_state, name="starttype", value=trim(starttype), rc=rc) 
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Initialze CTSM Gridded Component
    call ESMF_GridCompInitialize(lnd_gcomp, importState=cpl2lnd_state, exportState=lnd2cpl_state, &
         clock=lilac_clock, userRc=user_rc, rc=rc)
    if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing ctsm")
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing ctsm")
    call ESMF_LogWrite(subname//"CTSM gridded component initialized", ESMF_LOGMSG_INFO)

    if (couple_to_river) then
       ! Create import and export states for rof_gcomp (i.e. MOSART)
       cpl2rof_state = ESMF_StateCreate(name='state_to_river', stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       rof2cpl_state = ESMF_StateCreate(name='state_fr_river', stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Initialize MOSART Gridded Component
       call ESMF_GridCompInitialize(rof_gcomp, importState=cpl2rof_state, exportState=rof2cpl_state, &
            clock=lilac_clock, userRc=user_rc, rc=rc)
       if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing mosart")
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing mosart")
       call ESMF_LogWrite(subname//"MOSART gridded component initialized", ESMF_LOGMSG_INFO)
    end if

    ! -------------------------------------------------------------------------
    ! Initialize LILAC coupler components
    ! -------------------------------------------------------------------------

    ! Note that the lnd2cpl_state and cpl2lnd_state are each made up of 2 field bundles,
    ! one for the river and one for the atm -
    ! The following fills in the atm field bundle in cpl2lnd_state
    call ESMF_CplCompInitialize(cpl_atm2lnd_comp, importState=atm2cpl_state, exportState=cpl2lnd_state, &
         clock=lilac_clock, userRc=user_rc, rc=rc)
    if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing cpl_atm2lnd component")
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing cpl_atm2lnd component")
    call ESMF_LogWrite(subname//"coupler :: cpl_atm2lnd_comp initialized", ESMF_LOGMSG_INFO)

    ! The following maps the atm field bundle in lnd2cpl_state to the atm mesh
    call ESMF_CplCompInitialize(cpl_lnd2atm_comp, importState=lnd2cpl_state, exportState=cpl2atm_state, &
         clock=lilac_clock, userRc=user_rc, rc=rc)
    if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing cpl_lnd2atm component")
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing cpl_lnd2atm component")
    call ESMF_LogWrite(subname//"coupler :: cpl_lnd2atm_comp initialized", ESMF_LOGMSG_INFO)

    if (couple_to_river) then
       ! The following maps the rof field bundle in lnd2cpl_state to the rof mesh
       call ESMF_CplCompInitialize(cpl_lnd2rof_comp, importState=lnd2cpl_state, exportState=cpl2rof_state, &
            clock=lilac_clock, userRc=user_rc, rc=rc)
       if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing cpl_lnd2rof component")
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing cpl_lnd2rof component")
       call ESMF_LogWrite(subname//"coupler :: cpl_atm2lnd_comp initialized", ESMF_LOGMSG_INFO)

       ! The following fills in the rof field bundle in cpl2lnd_state
       call ESMF_CplCompInitialize(cpl_rof2lnd_comp, importState=rof2cpl_state, exportState=cpl2lnd_state, &
            clock=lilac_clock, userRc=user_rc, rc=rc)
       if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing cpl_lnd2atm component")
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing cpl_lnd2atm component")
       call ESMF_LogWrite(subname//"coupler :: cpl_lnd2atm_comp initialized", ESMF_LOGMSG_INFO)
    end if

    if (mytask == 0) then
       write(logunit,*) trim(subname) // "finished lilac initialization"
    end if

    !-------------------------------------------------------------------------
    ! Initialize atmaero stream data (using share strearm capability from CIME)
    !-------------------------------------------------------------------------

    call lilac_atmaero_init(atm2cpl_state, lilac_clock, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing lilac_atmaero_init")

    !-------------------------------------------------------------------------
    ! Initialize lilac_io_mod module data
    !-------------------------------------------------------------------------

    call lilac_io_init()
    call ESMF_LogWrite(subname//"initialized lilac io ...", ESMF_LOGMSG_INFO)

    !-------------------------------------------------------------------------
    ! Initialize lilac history output
    !-------------------------------------------------------------------------

    call lilac_history_init(lilac_clock, caseid, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing lilac_history_init")
    call ESMF_LogWrite(subname//"initialized lilac history output ...", ESMF_LOGMSG_INFO)

  end subroutine lilac_init2

  !========================================================================

  subroutine lilac_run(write_restarts_now, stop_now)

    ! input/output variables
    logical, intent(in) :: write_restarts_now ! if true, CTSM will write restarts at end of time step
    logical, intent(in) :: stop_now           ! if true, CTSM will do some finalization at end of time step

    ! local variables
    type(ESMF_Alarm)            :: lilac_history_alarm
    type(ESMF_Alarm)            :: lilac_restart_alarm
    type(ESMF_State)            :: importState, exportState
    integer                     :: user_rc, rc
    character(len=*), parameter :: subname=trim(modname)//': [lilac_run] '
    !------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (mytask == 0) then
       write(logunit,*)  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       write(logunit,*)  "               Lilac Run               "
       write(logunit,*)  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    end if

    ! Note that the lilac caps for ctsm and possible mosart will
    ! listen to the restart and stop alarms on the lilac clock

    ! Set the clock restart alarm if restart_alarm_ringing is true
    if (write_restarts_now) then
       ! Turn on lilac restart alarm (this will be needed by ctsm)
       call ESMF_ClockGetAlarm(lilac_clock, 'lilac_restart_alarm', lilac_restart_alarm, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in obtaining lilac_restart_alarm")
       call ESMF_AlarmRingerOn(lilac_restart_alarm, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running lilac atm_cap")
       call ESMF_LogWrite(subname//"lilac restart alarm is ringing", ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("Error in querying lilac restart alarm ring")

       ! Write out lilac restart output if lilac_restart_alarm is ringing
       call lilac_time_restart_write(lilac_clock, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in restart write")
    end if

    ! Set the clock stop alarm if stop_alarm_ringing is true
    if (stop_now) then
       call ESMF_ClockGetAlarm(lilac_clock, 'lilac_stop_alarm', lilac_stop_alarm, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in obtaining lilac_stop_alarm")
       call ESMF_AlarmRingerOn(lilac_stop_alarm, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running lilac atm_cap")
       call ESMF_LogWrite(subname//"lilac stop alarm is ringing", ESMF_LOGMSG_INFO)
    end if

    ! Run lilac atmcap - update the cpl2atm_state
    call ESMF_LogWrite(subname//"running lilac atmos_cap", ESMF_LOGMSG_INFO)
    if (mytask == 0) write(logunit,*) "Running atmos_cap gridded component , rc =", rc
    call ESMF_GridCompRun(atm_gcomp, importState=cpl2atm_state, exportState=atm2cpl_state, &
         clock=lilac_clock, userRc=user_rc, rc=rc)
    if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running lilac atm_cap")
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running lilac atm_cap")

    ! Update prescribed aerosols  atm2cpl_a_state
    call lilac_atmaero_interp(lilac_clock, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running lilac_atmaero_interp")

    ! Make sure all atm2lnd fields have been set
    call a2l_fields%check_all_set()

    ! Run cpl_atm2lnd
    call ESMF_LogWrite(subname//"running cpl_atm2lnd_comp ", ESMF_LOGMSG_INFO)
    if (mytask == 0) write(logunit,*) "Running coupler component..... cpl_atm2lnd_comp"
    call ESMF_CplCompRun(cpl_atm2lnd_comp, importState=atm2cpl_state, exportState=cpl2lnd_state, &
         clock=lilac_clock, userRc=user_rc, rc=rc)
    if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running cpl_atm2lnd")
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running cpl_atm2lnd")

    ! Run ctsm
    ! Write ctsm restart file if lilac_restart_alarm is ringing
    ! Finalize ctsm if lilac_stop_alarm is ringing
    call ESMF_LogWrite(subname//"running ctsm", ESMF_LOGMSG_INFO)
    if (mytask == 0) write(logunit,*) "Running ctsm"
    call ESMF_GridCompRun(lnd_gcomp,  importState=cpl2lnd_state, exportState=lnd2cpl_state, &
         clock=lilac_clock, userRc=user_rc, rc=rc)
    if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running ctsm")
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running ctsm")

    ! Run cpl_lnd2atm
    call ESMF_LogWrite(subname//"running cpl_lnd2atm_comp ", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       write(logunit,*) "Running coupler component..... cpl_lnd2atm_comp , rc =", rc
    end if
    call ESMF_CplCompRun(cpl_lnd2atm_comp, importState=lnd2cpl_state, exportState=cpl2atm_state, &
         clock=lilac_clock, userRc=user_rc, rc=rc)
    if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in cpl_lnd2atm")
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in cpl_lnd2atm")

    if (couple_to_river) then
       ! Run cpl_lnd2rof
       call ESMF_LogWrite(subname//"running cpl_lnd2rof_comp ", ESMF_LOGMSG_INFO)
       if (mytask == 0) write(logunit,*) "Running coupler component..... cpl_lnd2rof_comp"
       call ESMF_CplCompRun(cpl_lnd2rof_comp, importState=lnd2cpl_state, exportState=cpl2rof_state, &
            clock=lilac_clock, userRc=user_rc, rc=rc)
       if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running cpl_lnd2rof")
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running cpl_lnd2rof")

       ! Run mosart
       ! Write mosart restart file if lilac_restart_alarm is ringing
       ! Finalize mosart if lilac_stop_alarm is ringing
       call ESMF_LogWrite(subname//"running mosart", ESMF_LOGMSG_INFO)
       if (mytask == 0) write(logunit,*) "Running mosart"
       call ESMF_GridCompRun(rof_gcomp, importState=cpl2rof_state, exportState=rof2cpl_state, &
            clock=lilac_clock, userRc=user_rc, rc=rc)
       if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running rof")
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running rof")

       ! Run cpl_rof2lnd
       ! TODO: uncommenting this needs to be tested
       ! call ESMF_LogWrite(subname//"running cpl_rof2lnd_comp ", ESMF_LOGMSG_INFO)
       ! if (mytask == 0) write(logunit,*) "Running coupler component..... cpl_rof2lnd_comp"
       ! call ESMF_CplCompRun(cpl_rof2lnd_comp, importState=rof2cpl_state, exportState=cpl2lnd_state, &
       !      clock=lilac_clock, userRc=user_rc, rc=rc)
       ! if (chkerr(user_rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running cpl_rof2lnd")
       ! if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running cpl_rof2lnd")
    end if

    ! Write out lilac history output if lilac_history_alarm is ringing
    call ESMF_ClockGetAlarm(lilac_clock, 'lilac_history_alarm', lilac_history_alarm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in obtaining lilac_history_alarm")
    if (ESMF_AlarmIsRinging(lilac_history_alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("Error in querying lilac history alarm ring")
       call ESMF_AlarmRingerOff( lilac_history_alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("Error in turning ringer off in lilac history alarm")
       if (couple_to_river) then
          call lilac_history_write(atm2cpl_state, cpl2atm_state, lnd2cpl_state, cpl2lnd_state, &
               rof2cpl_state=rof2cpl_state, cpl2rof_state=cpl2rof_state, clock=lilac_clock, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in history write")
       else
          call lilac_history_write(atm2cpl_state, cpl2atm_state, lnd2cpl_state, cpl2lnd_state, &
               clock=lilac_clock, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in history write")
       end if
    end if
       
    if (write_restarts_now) then
       call ESMF_ClockGetAlarm(lilac_clock, 'lilac_restart_alarm', lilac_restart_alarm, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in obtaining lilac_restart_alarm")
       call ESMF_AlarmRingerOff( lilac_restart_alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Reset atm2lnd provided flags for next time step
    call a2l_fields%reset_provided()

    ! Advance the lilac clock at the end of the time step
    call ESMF_ClockAdvance(lilac_clock, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in advancing time step")
    call ESMF_LogWrite(subname//"time is icremented now... (ClockAdvance)", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       write(logunit,*) "time is icremented now... (ClockAdvance) , rc =", rc
    end if

  end subroutine lilac_run

!========================================================================

  subroutine lilac_final( )

    ! local variables
    type(ESMF_State)            :: importState, exportState
    integer                     :: rc, user_rc
    character(len=*), parameter :: subname=trim(modname)//': [lilac_final] '
    !------------------------------------------------------------------------

    ! Initialize return code
    rc = ESMF_SUCCESS

    if (mytask == 0) then
       write(logunit,*)  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       write(logunit,*)  "        Lilac Finalizing               "
       write(logunit,*)  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    end if

    ! Gridded Component Finalizing!     ---       atmosphere
    call ESMF_GridCompFinalize(atm_gcomp, importState=cpl2atm_state, exportState=atm2cpl_state, clock=lilac_clock, &
         userRc=user_rc, rc=rc)
    if (chkerr(user_rc,__LINE__,u_FILE_u)) return
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//"atmos_cap or atm_gcomp is running", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       write(logunit,*) "Finalizing atmos_cap gridded component , rc =", rc
    end if

    ! Coupler component Finalizing     ---    coupler atmos to land
    call ESMF_CplCompFinalize(cpl_atm2lnd_comp, importState=atm2cpl_state, exportState=cpl2lnd_state, clock=lilac_clock, &
         userRc=user_rc, rc=rc)
    if (chkerr(user_rc,__LINE__,u_FILE_u)) return
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//"running cpl_atm2lnd_comp ", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       write(logunit,*) "Finalizing coupler component..... cpl_atm2lnd_comp , rc =", rc
    end if

    ! Gridded Component Finalizing!     ---       land
    call ESMF_GridCompFinalize(lnd_gcomp,  importState=cpl2lnd_state, exportState=lnd2cpl_state, clock=lilac_clock, &
         userRc=user_rc, rc=rc)
    if (chkerr(user_rc,__LINE__,u_FILE_u)) return
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//"lnd_cap or lnd_gcomp is running", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       write(logunit,*) "Finalizing lnd_cap gridded component , rc =", rc
    end if

    ! Coupler component Finalizing     ---    coupler land to atmos
    call ESMF_CplCompFinalize(cpl_lnd2atm_comp, importState=cpl2lnd_state, exportState=cpl2atm_state, clock=lilac_clock, &
         userRc=user_rc, rc=rc)
    if (chkerr(user_rc,__LINE__,u_FILE_u)) return
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//"running cpl_lnd2atm_comp ", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       write(logunit,*) "Finalizing coupler component..... cpl_lnd2atm_comp , rc =", rc
    end if

    ! Then clean them up
    call ESMF_LogWrite(subname//".........................", ESMF_LOGMSG_INFO)
    call ESMF_LogWrite(subname//"destroying all states    ", ESMF_LOGMSG_INFO)

    if (mytask == 0) then
       write(logunit,*) "ready to destroy all states"
    end if
    call ESMF_StateDestroy(atm2cpl_state , rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_StateDestroy(cpl2atm_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_StateDestroy(lnd2cpl_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_StateDestroy(cpl2lnd_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    if (couple_to_river) then
       call ESMF_StateDestroy(rof2cpl_state, rc=rc)
       if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
       call ESMF_StateDestroy(cpl2rof_state, rc=rc)
       if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    end if

    call ESMF_LogWrite(subname//"destroying all components    ", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       write(logunit,*) "ready to destroy all components"
    end if

    call ESMF_GridCompDestroy(atm_gcomp, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_GridCompDestroy(lnd_gcomp, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    if (couple_to_river) then
       call ESMF_GridCompDestroy(rof_gcomp, rc=rc)
       if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    end if

    call ESMF_CplCompDestroy(cpl_atm2lnd_comp, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_CplCompDestroy(cpl_lnd2atm_comp, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    call ESMF_LogWrite(subname//".........................", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       write(logunit,*) "end of Lilac Finalization routine"
    end if

    ! Finalize ESMF; keep mpi alive so that atmosphere can do any finalization needed
    ! before it calls MPI_Finalize
    call ESMF_Finalize  (endflag=ESMF_END_KEEPMPI)

  end subroutine lilac_final

end module lilac_mod
