module lilac_mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This is the driver for running CTSM and the ESMF lilac atm cap that
  ! is put in place to ensure that the host atmosphere does not need to
  ! know about ESMF
  !-----------------------------------------------------------------------

  ! External libraries
  use ESMF
  use mct_mod        , only : mct_world_init

  ! shr code routines
  use shr_pio_mod   , only : shr_pio_init1, shr_pio_init2
  use shr_sys_mod   , only : shr_sys_abort

  ! lilac routines
  use lilac_io      , only : lilac_io_init
  use lilac_utils   , only : lilac_init_lnd2atm, lilac_init_atm2lnd
  use lilac_utils   , only : gindex_atm, atm_mesh_filename
  use lilac_cpl     , only : cpl_atm2lnd_register, cpl_lnd2atm_register
  use lilac_cpl     , only : cpl_lnd2rof_register, cpl_rof2lnd_register
  use lilac_atmcap  , only : lilac_atmos_register
  use lilac_atmaero , only : lilac_atmaero_init
  use lilac_atmaero , only : lilac_atmaero_interp
  use lilac_history , only : lilac_history_write
  use lilac_methods , only : chkerr

  ! ctsm routines
  use lnd_comp_esmf , only : lnd_register ! ctsm routine

  ! mosart routines
  use rof_comp_esmf , only : rof_register ! mosart routine

  implicit none

  public :: lilac_init
  public :: lilac_run

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

  character(*) , parameter   :: modname     = "lilac_mod"
  character(*), parameter    :: u_FILE_u = &
       __FILE__

!========================================================================
contains
!========================================================================

  subroutine lilac_init(mpicom, atm_mesh_file, atm_global_index, atm_lons, atm_lats, &
       atm_calendar, atm_timestep, &
       atm_start_year, atm_start_mon, atm_start_day, atm_start_secs, &
       atm_stop_year, atm_stop_mon, atm_stop_day, atm_stop_secs)

    ! --------------------------------------------------------------------------------
    ! This is called by the host atmosphere
    ! --------------------------------------------------------------------------------

    ! input/output variables
    integer          , intent(inout) :: mpicom  ! input commiunicator from atm
    character(len=*) , intent(in)    :: atm_mesh_file
    integer          , intent(in)    :: atm_global_index(:)
    real             , intent(in)    :: atm_lons(:)
    real             , intent(in)    :: atm_lats(:)
    character(len=*) , intent(in)    :: atm_calendar
    integer          , intent(in)    :: atm_timestep
    integer          , intent(in)    :: atm_start_year !(yyyy)
    integer          , intent(in)    :: atm_start_mon  !(mm)
    integer          , intent(in)    :: atm_start_day
    integer          , intent(in)    :: atm_start_secs
    integer          , intent(in)    :: atm_stop_year  !(yyyy)
    integer          , intent(in)    :: atm_stop_mon   !(mm)
    integer          , intent(in)    :: atm_stop_day
    integer          , intent(in)    :: atm_stop_secs

    ! local variables
    type(ESMF_TimeInterval)     :: timeStep
    type(ESMF_Time)             :: startTime
    type(ESMF_Time)             :: stopTime
    type(ESMF_Alarm)            :: EAlarm_stop, EAlarm_rest
    integer                     :: yy,mm,dd,sec
    integer                     :: lsize
    type(ESMF_State)            :: importState, exportState
    type(ESMF_VM)               :: vm
    integer                     :: rc
    character(len=ESMF_MAXSTR)  :: cname   !components or cpl names
    integer                     :: ierr
    integer                     :: n, i
    integer                     :: fileunit
    integer, parameter          :: debug = 1   !-- internal debug level
    character(len=*), parameter :: subname=trim(modname)//': [lilac_init] '

    ! initialization of mct and pio
    integer           :: ncomps = 1                 ! for mct
    integer, pointer  :: mycomms(:)                 ! for mct
    integer, pointer  :: myids(:)                   ! for mct
    integer           :: compids(1) = (/1/)         ! for pio_init2 - array with component ids
    integer           :: comms(1)                   ! for both mct and pio_init2 - array with mpicoms
    character(len=32) :: compLabels(1) = (/'LND'/)  ! for pio_init2
    character(len=64) :: comp_name(1) = (/'LND'/)   ! for pio_init2
    logical           :: comp_iamin(1) = (/.true./) ! for pio init2
    !------------------------------------------------------------------------

    ! Initialize return code
    rc = ESMF_SUCCESS

    !-------------------------------------------------------------------------
    ! Initialize pio with first initialization
    ! AFTER call to MPI_init (which is in the host atm driver) and 
    ! BEFORE call to ESMF_Initialize
    !-------------------------------------------------------------------------
    call shr_pio_init1(ncomps=1, nlfilename="lilac_in", Global_Comm=mpicom)

    !-------------------------------------------------------------------------
    ! Initialize ESMF, set the default calendar and log type.
    !-------------------------------------------------------------------------

    ! TODO: cannot assume that the calendar is always gregorian unless CTSM assumes this as well
    ! Need to coordinate the calendar info between lilac and the host component
    call ESMF_Initialize(mpiCommunicator=mpicom, defaultCalKind=ESMF_CALKIND_GREGORIAN, &
         logappendflag=.false., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogSet(flush=.true.)
    call ESMF_LogWrite(subname//".........................", ESMF_LOGMSG_INFO)
    call ESMF_LogWrite(subname//"Initializing ESMF        ", ESMF_LOGMSG_INFO)

    call ESMF_VMGetGlobal(vm=vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=mytask,  rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------------------------------------------------
    ! Initialize MCT (this is needed for data model functionality)
    !-------------------------------------------
    allocate(mycomms(1), myids(1))
    mycomms = (/mpicom/) ; myids = (/1/)
    call mct_world_init(ncomps, mpicom, mycomms, myids)
    call ESMF_LogWrite(subname//"initialized mct ...  ", ESMF_LOGMSG_INFO)

    !-------------------------------------------------------------------------
    ! Initialize PIO with second initialization
    !-------------------------------------------------------------------------
    call shr_pio_init2(compids, compLabels, comp_iamin, (/mpicom/), (/mytask/))
    call ESMF_LogWrite(subname//"initialized shr_pio_init2 ...", ESMF_LOGMSG_INFO)

    !-------------------------------------------------------------------------
    ! Initial lilac_utils module variables
    !-------------------------------------------------------------------------

    ! Initialize gindex_atm
    lsize = size(atm_global_index)
    allocate(gindex_atm(lsize))
    gindex_atm(:) = atm_global_index(:)

    ! Initialize atm_mesh_filename
    atm_mesh_filename = atm_mesh_file

    ! Initialize datatypes atm2lnd and lnd2atm
    ! This must be done BEFORE the atmcap initialization - since the dataptr attributes
    ! are only needed to initialize the atmcap field bundles
    call lilac_init_atm2lnd(lsize)
    call lilac_init_lnd2atm(lsize)

    !-------------------------------------------------------------------------
    ! Create Gridded Component --  lilac atmos_cap
    !-------------------------------------------------------------------------
    cname = " LILAC atm cap "
    atm_gcomp = ESMF_GridCompCreate(name=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error lilac atmcap initialization')
    call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // "lilac atm cap gridded component created"
    end if

    !-------------------------------------------------------------------------
    ! Create Gridded Component -- CTSM land
    !-------------------------------------------------------------------------
    cname = " CTSM "
    lnd_gcomp = ESMF_GridCompCreate(name=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error lilac ctsm initialization')
    call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // " ctsm gridded component created"
    end if

    !-------------------------------------------------------------------------
    ! Create Gridded Component -- MOSART river
    !-------------------------------------------------------------------------
    cname = " MOSART "
    rof_gcomp = ESMF_GridCompCreate(name=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error lilac mosart initialization')
    call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // " mosart gridded component created"
    end if

    !-------------------------------------------------------------------------
    ! Create Coupling Component! --- Coupler  from atmos  to land
    !-------------------------------------------------------------------------
    cname = "Coupler from atmosphere to land"
    cpl_atm2lnd_comp = ESMF_CplCompCreate(name=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error lilac cpl_a2l initialization')
    call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // " coupler component (atmosphere to land) created"
    end if

    !-------------------------------------------------------------------------
    ! Create Coupling  Component!  -- Coupler from land to atmos
    !-------------------------------------------------------------------------
    cname = "Coupler from land to atmosphere"
    cpl_lnd2atm_comp = ESMF_CplCompCreate(name=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error lilac cpl_l2a initialization')
    call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // " coupler component (land to atmosphere) created"
    end if

    !-------------------------------------------------------------------------
    ! Create Coupling Component! --- Coupler  from rof  to land
    !-------------------------------------------------------------------------
    cname = "Coupler from river to land"
    cpl_rof2lnd_comp = ESMF_CplCompCreate(name=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error lilac cpl_r2l initialization')
    call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // " coupler component (atmosphere to land) created"
    end if

    !-------------------------------------------------------------------------
    ! Create Coupling  Component!  -- Coupler from land to atmos
    !-------------------------------------------------------------------------
    cname = "Coupler from land to river"
    cpl_lnd2rof_comp = ESMF_CplCompCreate(name=cname, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error lilac cpl_l2r initialization')
    call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // " coupler component (land to atmosphere) created"
    end if

    !-------------------------------------------------------------------------
    ! Register section -- set services -- atmcap
    !-------------------------------------------------------------------------
    call ESMF_GridCompSetServices(atm_gcomp, userRoutine=lilac_atmos_register, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('atm_gcomp register failure')
    call ESMF_LogWrite(subname//"  atmos SetServices finished!", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // " lilac atm cap setservices finished"
    end if

    !-------------------------------------------------------------------------
    ! Register section -- set services -- ctsm
    !-------------------------------------------------------------------------
    call ESMF_GridCompSetServices(lnd_gcomp, userRoutine=lnd_register, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('lnd_gcomp register failure')
    call ESMF_LogWrite(subname//"CSTM SetServices finished!", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // " CTSM setservices finished"
    end if

    !-------------------------------------------------------------------------
    ! Register section -- set services -- mosart
    !-------------------------------------------------------------------------
    call ESMF_GridCompSetServices(rof_gcomp, userRoutine=rof_register, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('rof_gcomp register failure')
    call ESMF_LogWrite(subname//"MOSART SetServices finished!", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // " CTSM setservices finished"
    end if

    !-------------------------------------------------------------------------
    ! Register section -- set services -- coupler atmosphere to land
    !-------------------------------------------------------------------------
    call ESMF_CplCompSetServices(cpl_atm2lnd_comp, userRoutine=cpl_atm2lnd_register, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('cpl_atm2lnd_comp register failure')
    call ESMF_LogWrite(subname//"Coupler from atmosphere to land  SetServices finished!", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // " coupler from atmosphere to land setservices finished"
    end if

    !-------------------------------------------------------------------------
    ! Register section -- set services -- river to land
    !-------------------------------------------------------------------------
    call ESMF_CplCompSetServices(cpl_rof2lnd_comp, userRoutine=cpl_rof2lnd_register, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('cpl_rof2lnd_comp register failure')
    call ESMF_LogWrite(subname//"Coupler from river to land  SetServices finished!", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // " coupler from river to land setservices finished"
    end if

    !-------------------------------------------------------------------------
    ! Register section -- set services -- coupler land to atmosphere
    !-------------------------------------------------------------------------
    call ESMF_CplCompSetServices(cpl_lnd2atm_comp, userRoutine=cpl_lnd2atm_register, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('cpl_lnd2atm_comp register failure')
    call ESMF_LogWrite(subname//"Coupler from land to atmosphere SetServices finished!", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // " coupler from land to atmosphere setservices finished"
    end if

    !-------------------------------------------------------------------------
    ! Register section -- set services -- coupler land to river
    !-------------------------------------------------------------------------
    call ESMF_CplCompSetServices(cpl_lnd2rof_comp, userRoutine=cpl_lnd2rof_register, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('cpl_lnd2rof_comp register failure')
    call ESMF_LogWrite(subname//"Coupler from land to river SetServices finished!", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // " coupler from land to river setservices finished"
    end if

    !-------------------------------------------------------------------------
    !  Create and initialize the lilac_clock and calendar
    !-------------------------------------------------------------------------

    if (trim(atm_calendar) == 'NOLEAP') then
       lilac_calendar = ESMF_CalendarCreate(name='NOLEAP', calkindflag=ESMF_CALKIND_NOLEAP, rc=rc )
    else if (trim(atm_calendar) == 'GREGORIAN') then
       lilac_calendar = ESMF_CalendarCreate(name='NOLEAP', calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc )
    else
       ! TODO: add supported calendars here
    end if

    call ESMF_TimeIntervalSet(TimeStep, s=atm_timestep, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeSet(StartTime, yy=atm_start_year, mm=atm_start_mon, dd=atm_start_day , s=atm_start_secs, &
          calendar=lilac_calendar, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeSet(StopTime , yy=atm_stop_year , mm=atm_stop_mon , dd=atm_stop_day  , s=atm_stop_secs , &
          calendar=lilac_calendar, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    lilac_clock = ESMF_ClockCreate(name='lilac_clock', TimeStep=TimeStep, startTime=StartTime, &
         RefTime=StartTime, stopTime=stopTime, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (mytask == 0) then
       print *, trim(subname) // "---------------------------------------"
       call ESMF_ClockPrint (lilac_clock, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_CalendarPrint (lilac_calendar , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       print *, trim(subname) // "---------------------------------------"
    end if

    ! Add a restart alarm to the clock
    lilac_restart_alarm = ESMF_AlarmCreate(lilac_clock, ringTime=StopTime, name='lilac_restart_alarm', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in initializing restart alarm')

    ! Add a stop alarm to the clock
    lilac_stop_alarm = ESMF_AlarmCreate(lilac_clock, ringTime=StopTime, name='lilac_stop_alarm', rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in initializing stop alarm')

    ! -------------------------------------------------------------------------
    ! Initialze lilac_atm gridded component
    ! First Create the empty import and export states used to pass data
    ! between components. (these are module variables)
    ! -------------------------------------------------------------------------

    atm2cpl_state = ESMF_StateCreate(name='state_from_atm', stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    cpl2atm_state = ESMF_StateCreate(name='state_to_atm', stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_GridCompInitialize(atm_gcomp, importState=cpl2atm_state, exportState=atm2cpl_state, &
         clock=lilac_clock, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing atmcap")
    call ESMF_LogWrite(subname//"lilac_atm gridded component initialized", ESMF_LOGMSG_INFO)

    ! -------------------------------------------------------------------------
    ! Initialze CTSM Gridded Component
    ! First Create the empty import and export states used to pass data
    ! between components. (these are module variables)
    ! -------------------------------------------------------------------------

    cpl2lnd_state = ESMF_StateCreate(name='state_to_land', stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lnd2cpl_state = ESMF_StateCreate(name='state_fr_land', stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_GridCompInitialize(lnd_gcomp, importState=cpl2lnd_state, exportState=lnd2cpl_state, &
         clock=lilac_clock, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing ctsm")
    call ESMF_LogWrite(subname//"CTSM gridded component initialized", ESMF_LOGMSG_INFO)

    ! -------------------------------------------------------------------------
    ! Initialize MOSART Gridded Component
    ! First Create the empty import and export states used to pass data
    ! between components. (these are module variables)
    ! -------------------------------------------------------------------------

    if (couple_to_river) then
       cpl2rof_state = ESMF_StateCreate(name='state_to_river', stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       rof2cpl_state = ESMF_StateCreate(name='state_fr_river', stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_GridCompInitialize(rof_gcomp, importState=cpl2rof_state, exportState=rof2cpl_state, &
            clock=lilac_clock, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing mosart")
       call ESMF_LogWrite(subname//"MOSART gridded component initialized", ESMF_LOGMSG_INFO)
    end if

    ! -------------------------------------------------------------------------
    ! Initialze LILAC coupler components
    ! -------------------------------------------------------------------------

    ! Note that the lnd2cpl_state and cpl2lnd_state are each made up of 2 field bundles, 
    ! one for the river and one for the atm - 

    ! The following fills in the atm field bundle in cpl2lnd_state
    call ESMF_CplCompInitialize(cpl_atm2lnd_comp, importState=atm2cpl_state, exportState=cpl2lnd_state, &
         clock=lilac_clock, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing cpl_atm2lnd component")
    call ESMF_LogWrite(subname//"coupler :: cpl_atm2lnd_comp initialized", ESMF_LOGMSG_INFO)

    ! The following maps the atm field bundle in lnd2cpl_state to the atm mesh
    call ESMF_CplCompInitialize(cpl_lnd2atm_comp, importState=lnd2cpl_state, exportState=cpl2atm_state, &
         clock=lilac_clock, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing cpl_lnd2atm component")
    call ESMF_LogWrite(subname//"coupler :: cpl_lnd2atm_comp initialized", ESMF_LOGMSG_INFO)

    if (couple_to_river) then
       ! The following maps the rof field bundle in lnd2cpl_state to the rof mesh
       call ESMF_CplCompInitialize(cpl_lnd2rof_comp, importState=lnd2cpl_state, exportState=cpl2rof_state, &
            clock=lilac_clock, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing cpl_lnd2rof component")
       call ESMF_LogWrite(subname//"coupler :: cpl_atm2lnd_comp initialized", ESMF_LOGMSG_INFO)

       ! The following fills in the rof field bundle in cpl2lnd_state
       call ESMF_CplCompInitialize(cpl_rof2lnd_comp, importState=rof2cpl_state, exportState=cpl2lnd_state, &
            clock=lilac_clock, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing cpl_lnd2atm component")
       call ESMF_LogWrite(subname//"coupler :: cpl_lnd2atm_comp initialized", ESMF_LOGMSG_INFO)
    end if

    if (mytask == 0) then
       print *, trim(subname) // "finished lilac initialization"
    end if

    !-------------------------------------------------------------------------
    ! Initialize lilac_io_mod module data
    !-------------------------------------------------------------------------

    call lilac_io_init()
    call ESMF_LogWrite(subname//"initialized lilac_io ...", ESMF_LOGMSG_INFO)

    !-------------------------------------------------------------------------
    ! Initialize atmaero stream data (using share strearm capability from CIME)
    !-------------------------------------------------------------------------

    call lilac_atmaero_init(atm2cpl_state, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in initializing lilac_atmaero_init")

  end subroutine lilac_init

  !========================================================================

  subroutine lilac_run(restart_alarm_is_ringing, stop_alarm_is_ringing)

    ! input/output variables
    logical, intent(in) :: restart_alarm_is_ringing
    logical, intent(in) :: stop_alarm_is_ringing

    ! local variables
    type(ESMF_State)            :: importState, exportState
    integer                     :: rc
    character(len=*), parameter :: subname=trim(modname)//': [lilac_run] '
    !------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (mytask == 0) then
       print *,  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       print *,  "               Lilac Run               "
       print *,  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    end if

    ! Set the clock restart alarm if restart_alarm_ringing is true
    if (restart_alarm_is_ringing) then
       call ESMF_AlarmRingerOn(lilac_restart_alarm, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running lilac atm_cap")
    end if

    ! Set the clock stop alarm if stop_alarm_ringing is true
    if (stop_alarm_is_ringing) then
       call ESMF_AlarmRingerOn(lilac_stop_alarm, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running lilac atm_cap")
    end if

    ! Run lilac atmcap - update the cpl2atm_state
    call ESMF_LogWrite(subname//"running lilac atmos_cap", ESMF_LOGMSG_INFO)
    if (mytask == 0) print *, "Running atmos_cap gridded component , rc =", rc
    call ESMF_GridCompRun(atm_gcomp, importState=cpl2atm_state, exportState=atm2cpl_state, &
         clock=lilac_clock, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running lilac atm_cap")

    ! Update prescribed aerosols  atm2cpl_a_state
    call lilac_atmaero_interp(atm2cpl_state, lilac_clock, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running lilac_atmaero_interp")

    ! Run cpl_atm2lnd
    call ESMF_LogWrite(subname//"running cpl_atm2lnd_comp ", ESMF_LOGMSG_INFO)
    if (mytask == 0) print *, "Running coupler component..... cpl_atm2lnd_comp"
    call ESMF_CplCompRun(cpl_atm2lnd_comp, importState=atm2cpl_state, exportState=cpl2lnd_state, &
         clock=lilac_clock, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running cpl_atm2lnd")

    ! Run ctsm
    call ESMF_LogWrite(subname//"running ctsm", ESMF_LOGMSG_INFO)
    if (mytask == 0) print *, "Running ctsm"
    call ESMF_GridCompRun(lnd_gcomp,  importState=cpl2lnd_state, exportState=lnd2cpl_state, &
         clock=lilac_clock, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running ctsm")

    ! Run cpl_lnd2atm
    call ESMF_LogWrite(subname//"running cpl_lnd2atm_comp ", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "Running coupler component..... cpl_lnd2atm_comp , rc =", rc
    end if
    call ESMF_CplCompRun(cpl_lnd2atm_comp, importState=lnd2cpl_state, exportState=cpl2atm_state, &
         clock=lilac_clock, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in cpl_lnd2atm")

    if (couple_to_river) then
       ! Run cpl_lnd2rof 
       call ESMF_LogWrite(subname//"running cpl_lnd2rof_comp ", ESMF_LOGMSG_INFO)
       if (mytask == 0) print *, "Running coupler component..... cpl_lnd2rof_comp"
       call ESMF_CplCompRun(cpl_lnd2rof_comp, importState=lnd2cpl_state, exportState=cpl2rof_state, &
            clock=lilac_clock, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running cpl_lnd2rof")

       ! Run mosart
       call ESMF_LogWrite(subname//"running mosart", ESMF_LOGMSG_INFO)
       if (mytask == 0) print *, "Running mosart"
       call ESMF_GridCompRun(rof_gcomp, importState=cpl2rof_state, exportState=rof2cpl_state, &
            clock=lilac_clock, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running ctsm")

       ! Run cpl_rof2lnd
       ! TODO: uncommenting this needs to be tested
       ! call ESMF_LogWrite(subname//"running cpl_rof2lnd_comp ", ESMF_LOGMSG_INFO)
       ! if (mytask == 0) print *, "Running coupler component..... cpl_rof2lnd_comp"
       ! call ESMF_CplCompRun(cpl_rof2lnd_comp, importState=rof2cpl_state, exportState=cpl2lnd_state, &
       !      clock=lilac_clock, rc=rc)
       ! if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in running cpl_rof2lnd")
    end if

    ! Write out history output
    if (couple_to_river) then
       call lilac_history_write(atm2cpl_state, cpl2atm_state, lnd2cpl_state, cpl2lnd_state, &
            rof2cpl_state=rof2cpl_state, cpl2rof_state=cpl2rof_state, clock=lilac_clock, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in history write")
    else
       call lilac_history_write(atm2cpl_state, cpl2atm_state, lnd2cpl_state, cpl2lnd_state, &
            clock=lilac_clock, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in history write")
    end if

    ! Advance the time at the end of the time step
    call ESMF_ClockAdvance(lilac_clock, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort("lilac error in advancing time step")
    call ESMF_LogWrite(subname//"time is icremented now... (ClockAdvance)", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "time is icremented now... (ClockAdvance) , rc =", rc
    end if

  end subroutine lilac_run

!========================================================================

  subroutine lilac_final( )

    ! local variables
    type(ESMF_State)            :: importState, exportState
    integer                     :: rc, userRC
    character(len=*), parameter :: subname=trim(modname)//': [lilac_final] '
    !------------------------------------------------------------------------

    ! Initialize return code
    rc = ESMF_SUCCESS

    if (mytask == 0) then
       print *,  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       print *,  "        Lilac Finalizing               "
       print *,  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    end if

    !-------------------------------------------------------------------------
    ! Gridded Component Finalizing!     ---       atmosphere
    !-------------------------------------------------------------------------
    call ESMF_GridCompFinalize(atm_gcomp, importState=cpl2atm_state, exportState=atm2cpl_state, clock=lilac_clock, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//"atmos_cap or atm_gcomp is running", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "Finalizing atmos_cap gridded component , rc =", rc
    end if

    !-------------------------------------------------------------------------
    ! Coupler component Finalizing     ---    coupler atmos to land
    !-------------------------------------------------------------------------
    call ESMF_CplCompFinalize(cpl_atm2lnd_comp, importState=atm2cpl_state, exportState=cpl2lnd_state, clock=lilac_clock, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//"running cpl_atm2lnd_comp ", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "Finalizing coupler component..... cpl_atm2lnd_comp , rc =", rc
    end if

    !-------------------------------------------------------------------------
    ! Gridded Component Finalizing!     ---       land
    !-------------------------------------------------------------------------
    call ESMF_GridCompFinalize(lnd_gcomp,  importState=cpl2lnd_state, exportState=lnd2cpl_state, clock=lilac_clock, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//"lnd_cap or lnd_gcomp is running", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "Finalizing lnd_cap gridded component , rc =", rc
    end if

    !-------------------------------------------------------------------------
    ! Coupler component Finalizing     ---    coupler land to atmos
    !-------------------------------------------------------------------------
    call ESMF_CplCompFinalize(cpl_lnd2atm_comp, importState=cpl2lnd_state, exportState=cpl2atm_state, clock=lilac_clock, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//"running cpl_lnd2atm_comp ", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "Finalizing coupler component..... cpl_lnd2atm_comp , rc =", rc
    end if

    ! Then clean them up
    call ESMF_LogWrite(subname//".........................", ESMF_LOGMSG_INFO)
    call ESMF_LogWrite(subname//"destroying all states    ", ESMF_LOGMSG_INFO)

    if (mytask == 0) then
       print *, "ready to destroy all states"
    end if
    call ESMF_StateDestroy(atm2cpl_state , rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_StateDestroy(cpl2atm_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_StateDestroy(lnd2cpl_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_StateDestroy(cpl2lnd_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_StateDestroy(rof2cpl_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_StateDestroy(cpl2rof_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    call ESMF_LogWrite(subname//"destroying all components    ", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "ready to destroy all components"
    end if

    call ESMF_GridCompDestroy(atm_gcomp, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_GridCompDestroy(lnd_gcomp, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_GridCompDestroy(rof_gcomp, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    call ESMF_CplCompDestroy(cpl_atm2lnd_comp, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_CplCompDestroy(cpl_lnd2atm_comp, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    call ESMF_LogWrite(subname//".........................", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "end of Lilac Finalization routine"
    end if

    ! Finalize ESMF
    call ESMF_Finalize  ( )

  end subroutine lilac_final

end module lilac_mod
