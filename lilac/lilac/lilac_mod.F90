module lilac_mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  !-----------------------------------------------------------------------

  use ESMF
  implicit none

  public :: lilac_init
  public :: lilac_run

  ! Clock, TimeInterval, and Times
  type(ESMF_Clock)           :: clock
  type(ESMF_TimeInterval)    :: timeStep
  type(ESMF_Time)            :: startTime
  type(ESMF_Time)            :: stopTime
  type(ESMF_Alarm)           :: EAlarm_stop, EAlarm_rest
  type(ESMF_Calendar),target :: calendar
  integer                    :: yy,mm,dd,sec

  ! Gridded Components and Coupling Components
  type(ESMF_GridComp)        :: atm_gcomp
  type(ESMF_GridComp)        :: lnd_gcomp
  type(ESMF_CplComp)         :: cpl_atm2lnd_comp
  type(ESMF_CplComp)         :: cpl_lnd2atm_comp
  type(ESMF_State)           :: atm2lnd_l_state, atm2lnd_a_state
  type(ESMF_State)           :: lnd2atm_a_state, lnd2atm_l_state
  character(*) , parameter   :: modname     = "lilac_mod"

  integer :: mytask

!========================================================================
contains
!========================================================================

  subroutine lilac_init(lsize)

    ! --------------------------------------------------------------------------------
    ! This is called by the host atmosphere
    ! --------------------------------------------------------------------------------

    use lilac_utils   ,  only : lilac_init_lnd2atm, lilac_init_atm2lnd
    use lilac_cpl     ,  only : cpl_atm2lnd_register, cpl_lnd2atm_register
    use lilac_atmcap  ,  only : atmos_register
    use lnd_comp_esmf ,  only : lnd_register !ctsm routine
    use shr_pio_mod   ,  only : shr_pio_init1

    ! input/output variables
    integer, intent(in) :: lsize

    ! local variables
    type(ESMF_State)            :: importState, exportState
    type(ESMF_VM)               :: vm
    integer                     :: rc
    character(len=ESMF_MAXSTR)  :: cname   !components or cpl names
    integer                     :: COMP_COMM
    integer                     :: ierr
    integer                     :: mpic   ! mpi communicator 
    integer                     :: n, i
    integer                     :: fileunit
    integer, parameter          :: debug = 1   !-- internal debug level
    character(len=*), parameter :: subname=trim(modname)//': [lilac_init] '

    ! Namelist and related variables
    integer ::  s_month, s_day, s_hour, s_min
    integer ::  e_month, e_day, e_hour, e_min
    namelist /input/ s_month, s_day, s_hour, s_min, e_month, e_day, e_hour, e_min
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

    ! Initialize pio (needed by CTSM) - TODO: this should be done within CTSM not here

    call ESMF_VMGetGlobal(vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_VMGet(vm, localPet=mytask, mpiCommunicator=mpic, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call shr_pio_init1(ncomps=1, nlfilename="drv_in", Global_Comm=mpic)  ! TODO: make the filename lilac_in

    ! Initialize atm2lnd and lnd2atm data types
    call lilac_init_atm2lnd(lsize)
    call lilac_init_lnd2atm(lsize)

    !-------------------------------------------------------------------------
    ! Read in configuration data -- namelist.input from host atmosphere(wrf)
    !-------------------------------------------------------------------------

    ! Read in namelist file ...

    if (mytask == 0) then
       print *,  "---------------------------------------"
    end if

    ! TODO: put checks for error below
    ! TODO: only the master lilac proc should read the namelist file and do a broadcast to the
    ! other processors
    open(newunit=fileunit, status="old", file="namelist_lilac",  action="read", iostat=rc)
    read(fileunit, input)
    close(fileunit)

    !-------------------------------------------------------------------------
    ! Create Gridded Component!  --  atmosphere ( atmos_cap)
    !-------------------------------------------------------------------------
    cname = " Atmosphere or Atmosphere Cap"
    atm_gcomp = ESMF_GridCompCreate(name=cname, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "Atmosphere Gridded Component Created!"
    end if

    !-------------------------------------------------------------------------
    ! Create Gridded Component!  --- CTSM land ( land_capX )
    !-------------------------------------------------------------------------
    cname = " Land ctsm "
    lnd_gcomp = ESMF_GridCompCreate(name=cname, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "  Land (ctsm) Gridded Component Created!"
    end if

    !-------------------------------------------------------------------------
    ! Create Coupling Component! --- Coupler  from atmos  to land
    !-------------------------------------------------------------------------
    cname = "Coupler from atmosphere to land"
    cpl_atm2lnd_comp = ESMF_CplCompCreate(name=cname, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "1st Coupler Component (atmosphere to land ) Created!"
    end if

    !-------------------------------------------------------------------------
    ! Create Coupling  Component!  -- Coupler from land to atmos
    !-------------------------------------------------------------------------
    cname = "Coupler from land to atmosphere"
    cpl_lnd2atm_comp = ESMF_CplCompCreate(name=cname, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"Created "//trim(cname)//" component", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "2nd Coupler Component (land to atmosphere) Created!"
    end if

    !-------------------------------------------------------------------------
    ! Register section -- set services -- atmos_cap
    !-------------------------------------------------------------------------
    call ESMF_GridCompSetServices(atm_gcomp, userRoutine=atmos_register, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"  atmos SetServices finished!", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "  Atmosphere Gridded Component SetServices finished!"
    end if

    !-------------------------------------------------------------------------
    ! Register section -- set services -- land cap
    !-------------------------------------------------------------------------
    call ESMF_GridCompSetServices(lnd_gcomp, userRoutine=lnd_register, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"land SetServices finished!", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "Land  Gridded Component SetServices finished!"
    end if

    !-------------------------------------------------------------------------
    ! Register section -- set services -- coupler atmosphere to land
    !-------------------------------------------------------------------------
    call ESMF_CplCompSetServices(cpl_atm2lnd_comp, userRoutine=cpl_atm2lnd_register, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"Coupler from atmosphere to land  SetServices finished!", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "Coupler from atmosphere to land SetServices finished!"
    end if

    !-------------------------------------------------------------------------
    ! Register section -- set services -- coupler land to atmosphere
    !-------------------------------------------------------------------------
    call ESMF_CplCompSetServices(cpl_lnd2atm_comp, userRoutine=cpl_lnd2atm_register, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"Coupler from land to atmosphere SetServices finished!", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "Coupler from land to atmosphere SetServices finished!"
    end if
    !-------------------------------------------------------------------------
    !  Create and initialize a clock!
    !  Clock is initialized here from namelist.input from WRF..... still we
    !  are looping over time from host atmosphere
    !-------------------------------------------------------------------------
    calendar = ESMF_CalendarCreate(name='lilac_drv_NOLEAP', calkindflag=ESMF_CALKIND_NOLEAP, rc=rc )
    call ESMF_TimeIntervalSet(TimeStep, s=2, rc=rc) ! time step every 2second
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    !call ESMF_TimeSet(startTime, yy=2003, mm=s_month, dd=s_day, h=s_hour, m=s_min, s=0, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    !call ESMF_TimeSet(stopTime,  yy=2003, mm=e_month, dd=e_day, h=e_hour, m=e_min, s=0, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    !clock = ESMF_ClockCreate(timeStep=timeStep, startTime=startTime, stopTime=stopTime, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_TimeSet(StartTime, yy=2000, mm=1, dd=1 , s=0, calendar=Calendar, rc=rc)
    call ESMF_TimeSet(StopTime , yy=2000, mm=03, dd=01, s=0, calendar=Calendar, rc=rc)
    !call ESMF_TimeIntervalSet(TimeStep, s=3600, rc=rc)
    call ESMF_TimeIntervalSet(TimeStep, s=1800, rc=rc)
    clock = ESMF_ClockCreate(name='lilac_drv_EClock', TimeStep=TimeStep, startTime=StartTime, &
         RefTime=StartTime, stopTime=stopTime, rc=rc)

    if (mytask == 0) then
       print *,  "---------------------------------------"
    end if
    !call ESMF_ClockPrint (clock, rc=rc)
    if (mytask == 0) then
       print *,  "======================================="
    end if
    !call ESMF_CalendarPrint ( calendar , rc=rc)
    if (mytask == 0) then
       print *,  "---------------------------------------"
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

    call ESMF_GridCompInitialize(atm_gcomp, importState=lnd2atm_a_state, exportState=atm2lnd_a_state, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"atmos_cap or atm_gcomp initialized", ESMF_LOGMSG_INFO)

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

    call ESMF_GridCompInitialize(lnd_gcomp, importState=atm2lnd_l_state, exportState=lnd2atm_l_state, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"lnd_cap or lnd_gcomp initialized", ESMF_LOGMSG_INFO)

    call ESMF_LogWrite(subname//"CTSM gridded component initialized", ESMF_LOGMSG_INFO)

    ! -------------------------------------------------------------------------
    ! Initialze LILAC coupler components
    ! -------------------------------------------------------------------------

    call ESMF_CplCompInitialize(cpl_atm2lnd_comp, importState=atm2lnd_a_state, exportState=atm2lnd_l_state, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"coupler :: cpl_atm2lnd_comp initialized", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "coupler :: cpl_atm2lnd_comp initialize finished" !, rc =", rc
    end if

    call ESMF_CplCompInitialize(cpl_lnd2atm_comp, importState=lnd2atm_l_state, exportState=lnd2atm_a_state, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"coupler :: cpl_lnd2atm_comp initialized", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "coupler :: cpl_lnd2atm_comp initialize finished" !, rc =", rc
    end if

  end subroutine lilac_init

  !========================================================================

  subroutine lilac_run( )

    ! local variables
    type(ESMF_State)            :: importState, exportState
    integer                     :: rc, userRC
    type (ESMF_Clock)           :: local_clock
    character(len=*), parameter :: subname=trim(modname)//': [lilac_run] '
    !------------------------------------------------------------------------

    ! Initialize return code
    rc = ESMF_SUCCESS

    if (mytask == 0) then
       print *,  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       print *,  "               Lilac Run               "
       print *,  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    end if

    !-------------------------------------------------------------------------
    ! Create a local clock from the general clock!
    !-------------------------------------------------------------------------

    local_clock = ESMF_ClockCreate(clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (mytask == 0) then
       print *, "Run Loop Start time"
    end if
    !call ESMF_ClockPrint(local_clock, options="currtime string", rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    !-------------------------------------------------------------------------
    ! We are running components in this order:
    ! 1- atmos_cap   2- cpl_atm2lnd! 3- lnd_cap   4- cpl_lnd2atm
    !-------------------------------------------------------------------------

    ! if we want to loop through clock in atmos cap.
    !do while (.NOT. ESMF_ClockIsStopTime(local_clock, rc=rc))
    call ESMF_GridCompRun(atm_gcomp, importState=lnd2atm_a_state, exportState=atm2lnd_a_state, &
         clock=local_clock, rc=rc, userRC=userRC)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=userRC, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite(subname//"atmos_cap or atm_gcomp is running", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "Running atmos_cap gridded component , rc =", rc
    end if

    call ESMF_CplCompRun(cpl_atm2lnd_comp, importState=atm2lnd_a_state, exportState=atm2lnd_l_state, &
         clock=local_clock, rc=rc , userRC=userRC)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=userRC, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"running cpl_atm2lnd_comp ", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "Running coupler component..... cpl_atm2lnd_comp , rc =", rc
    end if

    call ESMF_GridCompRun(lnd_gcomp,  importState=atm2lnd_l_state, exportState=lnd2atm_l_state, &
         clock=local_clock, rc=rc, userRC=userRC)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=userRC, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"lnd_cap or lnd_gcomp is running", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "Running lnd_cap gridded component  , rc =", rc
    end if

    call ESMF_CplCompRun(cpl_lnd2atm_comp, importState=lnd2atm_l_state, exportState=lnd2atm_a_state, &
         clock=local_clock, rc=rc, userRC=userRC)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=userRC, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"running cpl_lnd2atm_comp ", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "Running coupler component..... cpl_lnd2atm_comp , rc =", rc
    end if

    ! Advance the time
    call ESMF_ClockAdvance(local_clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
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
    call ESMF_GridCompFinalize(atm_gcomp, importState=lnd2atm_a_state, exportState=atm2lnd_a_state, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"atmos_cap or atm_gcomp is running", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "Finalizing atmos_cap gridded component , rc =", rc
    end if

    !-------------------------------------------------------------------------
    ! Coupler component Finalizing     ---    coupler atmos to land
    !-------------------------------------------------------------------------
    call ESMF_CplCompFinalize(cpl_atm2lnd_comp, importState=atm2lnd_a_state, exportState=atm2lnd_l_state, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"running cpl_atm2lnd_comp ", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "Finalizing coupler component..... cpl_atm2lnd_comp , rc =", rc
    end if

    !-------------------------------------------------------------------------
    ! Gridded Component Finalizing!     ---       land
    !-------------------------------------------------------------------------
    call ESMF_GridCompFinalize(lnd_gcomp,  importState=atm2lnd_l_state, exportState=lnd2atm_l_state, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//"lnd_cap or lnd_gcomp is running", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "Finalizing lnd_cap gridded component , rc =", rc
    end if

    !-------------------------------------------------------------------------
    ! Coupler component Finalizing     ---    coupler land to atmos
    !-------------------------------------------------------------------------
    call ESMF_CplCompFinalize(cpl_lnd2atm_comp, importState=lnd2atm_l_state, exportState=lnd2atm_a_state, clock=clock, rc=rc)
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
