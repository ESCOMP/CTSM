!Khoda
module LilacMod
use ESMF
use lilac_utils
!use DummyAtmos
use DummyAtmos, only : x2a_fields
use DummyAtmos, only : a2x_fields
use DummyAtmos, only : atmos_register

use cpl_mod



implicit none

   ! Clock, TimeInterval, and Times
   type(ESMF_Clock)           :: clock
   type(ESMF_TimeInterval)    :: timeStep
   type(ESMF_Time)            :: startTime
   type(ESMF_Time)            :: stopTime
   type(ESMF_Alarm)           :: EAlarm_stop, EAlarm_rest
   type(ESMF_Calendar),target :: calendar
   integer                    :: yy,mm,dd,sec

   character(*), parameter :: modname =  "(LilacMod)"

  !===============================================================================

  public :: lilac_init
  public :: lilac_run
  contains

  subroutine lilac_init( dum_var1, dum_var2)
    ! modules
    implicit none

    real, dimension(10) :: dum_var1
    real, dimension(10) :: dum_var2

    ! Component, and State
    type(ESMF_GridComp)      :: dummy_atmos_comp
    type(ESMF_GridComp)      :: dummy_land_comp
    type(ESMF_CplComp)      :: cpl_atm2lnd_comp
    type(ESMF_CplComp)      :: cpl_lnd2atm_comp


    type(ESMF_State)         :: coupledFlowState ! the coupled flow State
    type(ESMF_Mesh)          :: Emesh
    character(len=*), parameter :: subname=trim(modname)//':(lilac_init) '
    type(ESMF_State)         :: importState, exportState
    !character(len=*)        :: atm_mesh_filepath
    
    ! local variables
    integer                               :: rc, urc
    character(len=ESMF_MAXSTR)            :: cname1, cname2, cname3, cname4
    !integer, parameter     :: fldsMax = 100

    !-------------------------------------------------------------------------
    ! Initialize ESMF, set the default calendar and log type.
    !-------------------------------------------------------------------------
    call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    print *, "----------------------------"
    print *, "lilac Demo Application Start"

    !-------------------------------------------------------------------------
    ! Create Gridded Component!     --- dummy atmosphere 
    !-------------------------------------------------------------------------
    cname1 = "Dummy Atmosphere"

    dummy_atmos_comp = ESMF_GridCompCreate(name=cname1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Created "//trim(cname1)//" component", ESMF_LOGMSG_INFO)
    print *, "Dummy Atmosphere Gridded Component Created!"

    !-------------------------------------------------------------------------
    ! Create Gridded Component!   --- Coupler  from atmosphere  to land
    !-------------------------------------------------------------------------
    cname2 = "Coupler from atmosphere to land"
    cpl_atm2lnd_comp = ESMF_CplCompCreate(name=cname2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Created "//trim(cname2)//" component", ESMF_LOGMSG_INFO)
    print *, "1st Coupler Gridded Component (atmosphere to land ) Created!"

    !-------------------------------------------------------------------------
    ! Create Gridded Component!   --- dummy land (land cap)
    !-------------------------------------------------------------------------
    cname3 = "Dummy Land"

    dummy_land_comp = ESMF_GridCompCreate(name=cname3, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Created "//trim(cname3)//" component", ESMF_LOGMSG_INFO)
    print *, "Dummy Land  Gridded Component Created!"

    !-------------------------------------------------------------------------
    ! Create Gridded Component!  -- Coupler from land to atmos
    !-------------------------------------------------------------------------
    cname4 = "Coupler from land to atmosphere"
    cpl_lnd2atm_comp = ESMF_CplCompCreate(name=cname4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Created "//trim(cname4)//" component", ESMF_LOGMSG_INFO)
    print *, "2nd Coupler Gridded Component (land to atmosphere) Created!"


    ! ========================================================================
    !-------------------------------------------------------------------------
    ! Register section -- set services -- dummy atmosphere
    !-------------------------------------------------------------------------
    call ESMF_GridCompSetServices(dummy_atmos_comp, userRoutine=atmos_register, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"dummy atmos SetServices finished!", ESMF_LOGMSG_INFO)
    print *, "Dummy Atmosphere Gridded Component SetServices finished!"
    !-------------------------------------------------------------------------
    ! Register section -- set services  -- coupler atmosphere to land 
    !-------------------------------------------------------------------------
    call ESMF_CplCompSetServices(cpl_atm2lnd_comp, userRoutine=atmos_register, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Coupler from atmosphere to land  SetServices finished!", ESMF_LOGMSG_INFO)
    print *, "Coupler from atmosphere to land SetServices finished!"
    !-------------------------------------------------------------------------
    ! Register section -- set services -- dummy land
    !-------------------------------------------------------------------------
    call ESMF_GridCompSetServices(dummy_land_comp, userRoutine=atmos_register, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"dummy land SetServices finished!", ESMF_LOGMSG_INFO)
    print *, "Dummy Land  Gridded Component SetServices finished!"
    !-------------------------------------------------------------------------
    ! Register section -- set services -- coupler land to atmosphere
    !-------------------------------------------------------------------------
    call ESMF_CplCompSetServices(cpl_lnd2atm_comp, userRoutine=atmos_register, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Coupler from land to atmosphere SetServices finished!", ESMF_LOGMSG_INFO)
    print *, "Coupler from land to atmosphere SetServices finished!"

    ! ========================================================================

    !-------------------------------------------------------------------------
    !  Create and initialize a clock!
    ! ????? Should I create a clock here or in driver?
    !-------------------------------------------------------------------------
    calendar = ESMF_CalendarCreate(name='lilac_drv_NOLEAP', calkindflag=ESMF_CALKIND_NOLEAP, rc=rc )
    call ESMF_TimeSet(StartTime, yy=2000, mm=1, dd=1, s=0, calendar=Calendar, rc=rc)
    call ESMF_TimeSet(StopTime , yy=2000, mm=1, dd=10, s=0, calendar=Calendar, rc=rc)
    call ESMF_TimeIntervalSet(TimeStep, s=3600, rc=rc)
    clock = ESMF_ClockCreate(name='lilac_drv_EClock', TimeStep=TimeStep, startTime=StartTime, RefTime=StartTime, stopTime=stopTime, rc=rc)

    !print *, 
    !clock = ESMF_ClockCreate(timeStep=timeStep, startTime=startTime, stopTime=stopTime, rc=rc)
    !EClock = ESMF_ClockCreate(name='lilac_drv_EClock', TimeStep=TimeStep, startTime=StartTime, RefTime=StartTime, stopTime=stopTime, rc=rc)

    !-------------------------------------------------------------------------
    !  Atmosphere Initialization....
    !-------------------------------------------------------------------------
    call ESMF_GridCompInitialize(dummy_atmos_comp, &
        importState=importState, exportState=exportState, &
        clock=clock, rc=rc)
    !, dum_var1= dum_var1, dum_var2= dum_var2)
    !call ESMF_GridCompInitialize(self%land_comp, importState=self%land_import, exportState=self%land_export, rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out



  end subroutine lilac_init

  subroutine lilac_run(dum_var1, dum_var2)

    use DummyAtmos, only : x2a_fields
    use DummyAtmos, only : a2x_fields

    real, dimension(:) :: dum_var1  ! from host atm
    real, dimension(:) :: dum_var2 ! to host atm

    integer                :: n, num

    !integer, parameter     :: fldsMax = 100
    integer                :: fldsToLnd_num = 0
    integer                :: fldsFrLnd_num = 0

    type (fld_list_type)   :: fldsToLnd(fldsMax)
    type (fld_list_type)   :: fldsFrLnd(fldsMax)
    !-----------------------------------------
    !-----------------------------------------
    type(ESMF_State)     :: importState, exportState

    !search through fldlist array to find the right fldist object to do the copy - say its index N

    !x2a_fields(n)%datafld1d(:) = dum_var_input(:)

    !call ESMF_CplCompRun(cpl_atm2lnd, rc=rc)

    !call ESMF_GridCompRun(lndcomp, rc=rc)

    !call ESMF_CplCompRun(cpl_lnd2atm, rc=rc)

    !dum_var_output(:) = a2x_fields(N)%datafld1d(:)

  end subroutine lilac_run

end module LilacMod

