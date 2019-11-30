module lilac_mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This is the driver for running CTSM and the ESMF lilac atm cap that
  ! is put in place to ensure that the host atmosphere does not need to
  ! know about ESMF
  !-----------------------------------------------------------------------

  use ESMF
  implicit none

  public :: lilac_init
  public :: lilac_run

  ! Gridded components and states in gridded components
  type(ESMF_GridComp) :: atm_gcomp
  type(ESMF_GridComp) :: lnd_gcomp

  ! Coupler components
  type(ESMF_CplComp)  :: cpl_atm2lnd_comp
  type(ESMF_CplComp)  :: cpl_lnd2atm_comp

  ! States
  type(ESMF_State)    :: atm2lnd_l_state, atm2lnd_a_state
  type(ESMF_State)    :: lnd2atm_a_state, lnd2atm_l_state

  ! Clock, TimeInterval, and Times
  type(ESMF_Clock)           :: lilac_clock
  type(ESMF_Calendar),target :: lilac_calendar
  type(ESMF_Alarm)           :: lilac_restart_alarm
  type(ESMF_Alarm)           :: lilac_stop_alarm

  character(*) , parameter   :: modname     = "lilac_mod"

  integer :: mytask

!========================================================================
contains
!========================================================================

  subroutine lilac_init(atm_mesh_file, atm_global_index, atm_lons, atm_lats, &
       atm_calendar, atm_timestep, &
       atm_start_year, atm_start_mon, atm_start_day, atm_start_secs, &
       atm_stop_year, atm_stop_mon, atm_stop_day, atm_stop_secs)

    ! --------------------------------------------------------------------------------
    ! This is called by the host atmosphere
    ! --------------------------------------------------------------------------------

    use lilac_io      , only : lilac_io_init
    use lilac_utils   , only : lilac_init_lnd2atm, lilac_init_atm2lnd
    use lilac_utils   , only : gindex_atm, atm_mesh_filename
    use lilac_cpl     , only : cpl_atm2lnd_register, cpl_lnd2atm_register
    use lilac_atmcap  , only : lilac_atmos_register
    use lnd_comp_esmf , only : lnd_register !ctsm routine
    use shr_pio_mod   , only : shr_pio_init1
    use shr_sys_mod   , only : shr_sys_abort

    ! input/output variables
    character(len=*) , intent(in) :: atm_mesh_file
    integer          , intent(in) :: atm_global_index(:)
    real             , intent(in) :: atm_lons(:)
    real             , intent(in) :: atm_lats(:)
    character(len=*) , intent(in) :: atm_calendar
    integer          , intent(in) :: atm_timestep 
    integer          , intent(in) :: atm_start_year !(yyyy)
    integer          , intent(in) :: atm_start_mon  !(mm)
    integer          , intent(in) :: atm_start_day
    integer          , intent(in) :: atm_start_secs
    integer          , intent(in) :: atm_stop_year  !(yyyy)
    integer          , intent(in) :: atm_stop_mon   !(mm)
    integer          , intent(in) :: atm_stop_day
    integer          , intent(in) :: atm_stop_secs

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
    integer                     :: mpic   ! mpi communicator 
    integer                     :: n, i
    integer                     :: fileunit
    integer, parameter          :: debug = 1   !-- internal debug level
    character(len=*), parameter :: subname=trim(modname)//': [lilac_init] '
    !------------------------------------------------------------------------

    ! Initialize return code
    rc = ESMF_SUCCESS

    !-------------------------------------------------------------------------
    ! Initialize ESMF, set the default calendar and log type.
    !-------------------------------------------------------------------------

    ! TODO: cannot assume that the calendar is always gregorian unless CTSM assumes this as well
    ! Need to coordinate the calendar info between lilac and the host component
    call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, logappendflag=.false., rc=rc)  
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_LogSet(flush=.true.)

    call ESMF_LogWrite(subname//".........................", ESMF_LOGMSG_INFO)
    call ESMF_LogWrite(subname//"Initializing ESMF        ", ESMF_LOGMSG_INFO)

    !-------------------------------------------------------------------------
    ! Initialize pio with first initialization
    !-------------------------------------------------------------------------

    ! Initialize pio (needed by CTSM) - TODO: this should be done within CTSM not here

    call ESMF_VMGetGlobal(vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_VMGet(vm, localPet=mytask, mpiCommunicator=mpic, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call shr_pio_init1(ncomps=1, nlfilename="lilac_in", Global_Comm=mpic)

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
    call lilac_init_atm2lnd(lsize)
    call lilac_init_lnd2atm(lsize)

    !-------------------------------------------------------------------------
    ! Create Gridded Component --  lilac atmos_cap
    !-------------------------------------------------------------------------
    cname = " LILAC atm cap "
    atm_gcomp = ESMF_GridCompCreate(name=cname, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // "lilac atm cap gridded component created"
    end if

    !-------------------------------------------------------------------------
    ! Create Gridded Component -- CTSM land
    !-------------------------------------------------------------------------
    cname = " CTSM "
    lnd_gcomp = ESMF_GridCompCreate(name=cname, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // " ctsm gridded component created"
    end if

    !-------------------------------------------------------------------------
    ! Create Coupling Component! --- Coupler  from atmos  to land
    !-------------------------------------------------------------------------
    cname = "Coupler from atmosphere to land"
    cpl_atm2lnd_comp = ESMF_CplCompCreate(name=cname, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // " coupler component (atmosphere to land) created"
    end if

    !-------------------------------------------------------------------------
    ! Create Coupling  Component!  -- Coupler from land to atmos
    !-------------------------------------------------------------------------
    cname = "Coupler from land to atmosphere"
    cpl_lnd2atm_comp = ESMF_CplCompCreate(name=cname, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // " coupler component (land to atmosphere) created"
    end if

    !-------------------------------------------------------------------------
    ! Register section -- set services -- atmos_cap
    !-------------------------------------------------------------------------
    call ESMF_GridCompSetServices(atm_gcomp, userRoutine=lilac_atmos_register, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"  atmos SetServices finished!", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // " lilac atm cap setservices finished"
    end if

    !-------------------------------------------------------------------------
    ! Register section -- set services -- land cap
    !-------------------------------------------------------------------------
    call ESMF_GridCompSetServices(lnd_gcomp, userRoutine=lnd_register, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"land SetServices finished!", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // " CTSM setservices finished"
    end if

    !-------------------------------------------------------------------------
    ! Register section -- set services -- coupler atmosphere to land
    !-------------------------------------------------------------------------
    call ESMF_CplCompSetServices(cpl_atm2lnd_comp, userRoutine=cpl_atm2lnd_register, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"Coupler from atmosphere to land  SetServices finished!", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // " coupler from atmosphere to land setservices finished"
    end if

    !-------------------------------------------------------------------------
    ! Register section -- set services -- coupler land to atmosphere
    !-------------------------------------------------------------------------
    call ESMF_CplCompSetServices(cpl_lnd2atm_comp, userRoutine=cpl_lnd2atm_register, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"Coupler from land to atmosphere SetServices finished!", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // " coupler from land to atmosphere setservices finished"
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
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_TimeSet(StartTime, yy=atm_start_year, mm=atm_start_mon, dd=atm_start_day , s=atm_start_secs, &
          calendar=lilac_calendar, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_TimeSet(StopTime , yy=atm_stop_year , mm=atm_stop_mon , dd=atm_stop_day  , s=atm_stop_secs , &
          calendar=lilac_calendar, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    lilac_clock = ESMF_ClockCreate(name='lilac_clock', TimeStep=TimeStep, startTime=StartTime, &
         RefTime=StartTime, stopTime=stopTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (mytask == 0) then
       print *, trim(subname) // "---------------------------------------"
       call ESMF_ClockPrint (lilac_clock, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       call ESMF_CalendarPrint (lilac_calendar , rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       print *, trim(subname) // "---------------------------------------"
    end if

    ! Add a restart alarm to the clock
    lilac_restart_alarm = ESMF_AlarmCreate(lilac_clock, ringTime=StopTime, name='lilac_restart_alarm', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call shr_sys_abort('error in initializing restart alarm')
    end if

    ! Add a stop alarm to the clock
    lilac_stop_alarm = ESMF_AlarmCreate(lilac_clock, ringTime=StopTime, name='lilac_stop_alarm', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call shr_sys_abort('error in initializing stop alarm')
    end if

    ! -------------------------------------------------------------------------
    ! Initialze lilac_atm gridded component
    ! First Create the empty import and export states used to pass data
    ! between components. (these are module variables)
    ! -------------------------------------------------------------------------

    atm2lnd_a_state = ESMF_StateCreate(name='atm_state_on_atm_mesh', stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    lnd2atm_a_state = ESMF_StateCreate(name='lnd_state_on_lnd_mesh', stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridCompInitialize(atm_gcomp, importState=lnd2atm_a_state, exportState=atm2lnd_a_state, clock=lilac_clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call shr_sys_abort("lilac error in initializing atmcap")
    end if
    call ESMF_LogWrite(subname//"lilac_atm gridded component initialized", ESMF_LOGMSG_INFO)

    ! -------------------------------------------------------------------------
    ! Initialze CTSM Gridded Component
    ! First Create the empty import and export states used to pass data
    ! between components. (these are module variables)
    ! -------------------------------------------------------------------------

    atm2lnd_l_state = ESMF_StateCreate(name='atm_state_on_lnd_mesh', stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    lnd2atm_l_state = ESMF_StateCreate(name='lnd_state_on_atm_mesh', stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridCompInitialize(lnd_gcomp, importState=atm2lnd_l_state, exportState=lnd2atm_l_state, clock=lilac_clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call shr_sys_abort("lilac error in initializing ctsm")
    end if
    call ESMF_LogWrite(subname//"CTSM gridded component initialized", ESMF_LOGMSG_INFO)

    ! -------------------------------------------------------------------------
    ! Initialze LILAC coupler components
    ! -------------------------------------------------------------------------

    call ESMF_CplCompInitialize(cpl_atm2lnd_comp, importState=atm2lnd_a_state, exportState=atm2lnd_l_state, clock=lilac_clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call shr_sys_abort("lilac error in initializing cpl_lnd2atm coupler component")
    end if
    call ESMF_LogWrite(subname//"coupler :: cpl_atm2lnd_comp initialized", ESMF_LOGMSG_INFO)

    call ESMF_CplCompInitialize(cpl_lnd2atm_comp, importState=lnd2atm_l_state, exportState=lnd2atm_a_state, clock=lilac_clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call shr_sys_abort("lilac error in initializing cpl_atm2lnd coupler component")
    end if
    call ESMF_LogWrite(subname//"coupler :: cpl_lnd2atm_comp initialized", ESMF_LOGMSG_INFO)

    if (mytask == 0) then
       print *, trim(subname) // "finished lilac initialization"
    end if

    !-------------------------------------------------------------------------
    ! Initialize lilac_io_mod module data
    !-------------------------------------------------------------------------

    call lilac_io_init()
    call ESMF_LogWrite(subname//"initialized lilac_io ...", ESMF_LOGMSG_INFO)

  end subroutine lilac_init

  !========================================================================

  subroutine lilac_run(restart_alarm_is_ringing, stop_alarm_is_ringing)

    use shr_sys_mod  , only : shr_sys_abort
    use lilac_history, only : lilac_history_write

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
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call shr_sys_abort("lilac error in running lilac atm_cap")
       end if
    end if

    ! Set the clock stop alarm if stop_alarm_ringing is true
    if (stop_alarm_is_ringing) then
       call ESMF_AlarmRingerOn(lilac_stop_alarm, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call shr_sys_abort("lilac error in running lilac atm_cap")
       end if
    end if

    ! Run lilac atmcap
    call ESMF_LogWrite(subname//"running lilac atmos_cap", ESMF_LOGMSG_INFO)
    if (mytask == 0) print *, "Running atmos_cap gridded component , rc =", rc
    call ESMF_GridCompRun(atm_gcomp, importState=lnd2atm_a_state, exportState=atm2lnd_a_state, &
         clock=lilac_clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call shr_sys_abort("lilac error in running lilac atm_cap")
    end if

    ! Run cpl_atm2lnd
    call ESMF_LogWrite(subname//"running cpl_atm2lnd_comp ", ESMF_LOGMSG_INFO)
    if (mytask == 0) print *, "Running coupler component..... cpl_atm2lnd_comp"
    call ESMF_CplCompRun(cpl_atm2lnd_comp, importState=atm2lnd_a_state, exportState=atm2lnd_l_state, &
         clock=lilac_clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call shr_sys_abort("lilac error in running cpl_atm2lnd")
    end if

    ! Run ctsm
    call ESMF_LogWrite(subname//"running ctsm", ESMF_LOGMSG_INFO)
    if (mytask == 0) print *, "Running ctsm"
    call ESMF_GridCompRun(lnd_gcomp,  importState=atm2lnd_l_state, exportState=lnd2atm_l_state, &
         clock=lilac_clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call shr_sys_abort("lilac error in running ctsm")
    end if

    ! Run cpl_lnd2atm
    call ESMF_LogWrite(subname//"running cpl_lnd2atm_comp ", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "Running coupler component..... cpl_lnd2atm_comp , rc =", rc
    end if
    call ESMF_CplCompRun(cpl_lnd2atm_comp, importState=lnd2atm_l_state, exportState=lnd2atm_a_state, &
         clock=lilac_clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call shr_sys_abort("lilac error in cpl_lnd2atm")
    end if

    ! Write out history output
    call lilac_history_write(atm2lnd_a_state, atm2lnd_l_state, lnd2atm_l_state, lnd2atm_a_state, &
         lilac_clock, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call shr_sys_abort("lilac error in history write")
    end if

    ! Advance the time at the end of the time step
    call ESMF_ClockAdvance(lilac_clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call shr_sys_abort("lilac error in advancing time step")
    end if
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
    call ESMF_GridCompFinalize(atm_gcomp, importState=lnd2atm_a_state, exportState=atm2lnd_a_state, clock=lilac_clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"atmos_cap or atm_gcomp is running", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "Finalizing atmos_cap gridded component , rc =", rc
    end if

    !-------------------------------------------------------------------------
    ! Coupler component Finalizing     ---    coupler atmos to land
    !-------------------------------------------------------------------------
    call ESMF_CplCompFinalize(cpl_atm2lnd_comp, importState=atm2lnd_a_state, exportState=atm2lnd_l_state, clock=lilac_clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"running cpl_atm2lnd_comp ", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "Finalizing coupler component..... cpl_atm2lnd_comp , rc =", rc
    end if

    !-------------------------------------------------------------------------
    ! Gridded Component Finalizing!     ---       land
    !-------------------------------------------------------------------------
    call ESMF_GridCompFinalize(lnd_gcomp,  importState=atm2lnd_l_state, exportState=lnd2atm_l_state, clock=lilac_clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"lnd_cap or lnd_gcomp is running", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "Finalizing lnd_cap gridded component , rc =", rc
    end if

    !-------------------------------------------------------------------------
    ! Coupler component Finalizing     ---    coupler land to atmos
    !-------------------------------------------------------------------------
    call ESMF_CplCompFinalize(cpl_lnd2atm_comp, importState=lnd2atm_l_state, exportState=lnd2atm_a_state, clock=lilac_clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
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
    call ESMF_StateDestroy(atm2lnd_a_state , rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_StateDestroy(atm2lnd_l_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_StateDestroy(lnd2atm_a_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_StateDestroy(lnd2atm_l_state, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    call ESMF_LogWrite(subname//"destroying all components    ", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "ready to destroy all components"
    end if

    call ESMF_GridCompDestroy(atm_gcomp, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_GridCompDestroy(lnd_gcomp, rc=rc)
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
