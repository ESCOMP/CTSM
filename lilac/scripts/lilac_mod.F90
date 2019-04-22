module lilac_mod
use ESMF
use lilac_utils
!use DummyAtmos
use atmos_cap, only : atmos_register
use lnd_cap,   only : lnd_register
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

    type atm2lnd_data1d_type
       real(r8), pointer :: uwind(:)
       real(r8), pointer :: vwind(:)
       real(r8), pointer :: tbot(:)
    end type atm2lnd_data1d_type

    type atm2lnd_data2d_type
       real(r8), pointer :: uwind(:,:)
       real(r8), pointer :: vwind(:,:)
       real(r8), pointer :: tbot(:,:)
    end type atm2lnd_data2d_type

    type lnd2atm_data1d_type
       real(r8), pointer :: lwup(:)
       real(r8), pointer :: taux(:)
       real(r8), pointer :: tauy(:)
    end type lnd2atm_data1d_type

    type lnd2atm_data2d_type
       real(r8), pointer :: lwup(:,:)
       real(r8), pointer :: taux(:,:)
       real(r8), pointer :: tauy(:,:)
    end type lnd2atm_data2d_type

  !===============================================================================

  public :: lilac_init
  public :: lilac_run
  contains

  subroutine lilac_init( atm2lnd1d, atm2lnd2d, lnd2atm1d, lnd2atm2d)

    use atmos_cap, only : a2l_fields !****
    use atmos_cap, only : l2a_fields !****

    ! modules
    implicit none

    type(atm2lnd_data1d_type), intent(in), optional :: atm2lnd1d
    type(atm2lnd_data2d_type), intent(in), optional :: atm2lnd2d
    type(lnd2atm_data1d_type), intent(in), optional :: lnd2atm1d
    type(lnd2atm_data2d_type), intent(in), optional :: lnd2atm2d

    ! Component, and State
    type(ESMF_GridComp)      :: dummy_atmos_comp
    type(ESMF_GridComp)      :: dummy_land_comp
    type(ESMF_CplComp)       :: cpl_atm2lnd_comp
    type(ESMF_CplComp)       :: cpl_lnd2atm_comp


    type(ESMF_State)         :: coupledFlowState ! the coupled flow State
    type(ESMF_Mesh)          :: Emesh
    character(len=*), parameter :: subname=trim(modname)//':(lilac_init) '
    type(ESMF_State)         :: importState, exportState
    type(ESMF_State)         :: atm2lnd_l_state , atm2lnd_a_state
    type(ESMF_State)         :: lnd2atm_a_state, lnd2atm_l_state

    !character(len=*)        :: atm_mesh_filepath   !!! For now this is hard
    !coded in the atmos init

    ! local variables
    integer                               :: rc, urc
    character(len=ESMF_MAXSTR)            :: cname1, cname2, cname3, cname4
    !integer, parameter                   :: fldsMax = 100

    !-------------------------------------------------------------------------
    ! Initialize ESMF, set the default calendar and log type.
    !-------------------------------------------------------------------------
    call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    print *, "----------------------------"
    print *, "lilac Demo Application Start"

    !-------------------------------------------------------------------------
    ! Create Field lists -- Basically create a list of fields and add a default
    ! value to them.
    !-------------------------------------------------------------------------
    fldsFrCpl_num = 1
    fldsToCpl_num = 1

    if (present(var1_a2l_1d)) then
       a2l_fields%name = 'uwind'
       a2l_fields%arrptr1d => atm2lnd_data1d%uwind
       call create_fldlists(flds_a2l, fldsfldsToCpl, fldsToCpl_num, fldsFrCpl_num, 'mesh')
    else
       a2l_fields%name = 'name'
       a2l_fields%arrptr2d => var1_a2l_2d
       call create_fldlists(fldsFrCpl, fldsToCpl, fldsToCpl_num, fldsFrCpl_num, 'grid')
    end if

    !-------------------------------------------------------------------------
    ! Create Gridded Component!     --- dummy atmosphere 
    !-------------------------------------------------------------------------
    cname1 = "Dummy Atmosphere or Atmosphere Cap"

    dummy_atmos_comp = ESMF_GridCompCreate(name=cname1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Created "//trim(cname1)//" component", ESMF_LOGMSG_INFO)
    print *, "Dummy Atmosphere Gridded Component Created!"

    !-------------------------------------------------------------------------
    ! Create Gridded Component!   --- dummy land (land cap)
    !-------------------------------------------------------------------------
    cname2 = "Dummy Land or Land Cap"

    dummy_land_comp = ESMF_GridCompCreate(name=cname2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Created "//trim(cname2)//" component", ESMF_LOGMSG_INFO)
    print *, "Dummy Land  Gridded Component Created!"

    !-------------------------------------------------------------------------
    ! Create Gridded Component!   --- Coupler  from atmos  to land
    !-------------------------------------------------------------------------
    cname3 = "Coupler from atmosphere to land"
    cpl_atm2lnd_comp = ESMF_CplCompCreate(name=cname3, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Created "//trim(cname3)//" component", ESMF_LOGMSG_INFO)
    print *, "1st Coupler Gridded Component (atmosphere to land ) Created!"

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
    ! Register section -- set services -- land cap
    !-------------------------------------------------------------------------
    call ESMF_GridCompSetServices(dummy_land_comp, userRoutine=lnd_register, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"land SetServices finished!", ESMF_LOGMSG_INFO)
    print *, "Land  Gridded Component SetServices finished!"
    !-------------------------------------------------------------------------
    ! Register section -- set services  -- coupler atmosphere to land 
    !-------------------------------------------------------------------------
    call ESMF_CplCompSetServices(cpl_atm2lnd_comp, userRoutine=cpl_atm2lnd_register, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Coupler from atmosphere to land  SetServices finished!", ESMF_LOGMSG_INFO)
    print *, "Coupler from atmosphere to land SetServices finished!"
    !-------------------------------------------------------------------------
    ! Register section -- set services -- coupler land to atmosphere
    !-------------------------------------------------------------------------
    call ESMF_CplCompSetServices(cpl_lnd2atm_comp, userRoutine=cpl_lnd2atm_register, rc=rc)
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
    ! Create the necessary import and export states used to pass data
    !  between components.
    !-------------------------------------------------------------------------

    ! following 4 states are lilac module variables

    atm2lnd_a_state = ESMF_StateCreate(name=cname1, stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    atm2lnd_l_state = ESMF_StateCreate(name=cname1, stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    lnd2atm_a_state = ESMF_StateCreate(name=cname2, stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    lnd2atm_l_state = ESMF_StateCreate(name=cname2, stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out



    call ESMF_LogWrite(subname//"Empty import and export states are created!!", ESMF_LOGMSG_INFO)
    print *, "Empty import and export states are created!!"

    !, dum_var1= dum_var1, dum_var2= dum_var2)
    ! returns a valid state_to_lnd_atm and an empty state_from_land_atmgrid
    ! ========================================================================
    !-------------------------------------------------------------------------
    ! Grid Componenet Initialization -- 1- atmos cap 2- lnd cap 3- cpl_atm2lnd
    ! 4- cpl_lnd2atm
    !-------------------------------------------------------------------------

    call ESMF_GridCompInitialize(dummy_atmos_comp, importState=lnd2atm_a_state, exportState=atm2lnd_a_state, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"atmos_cap or dummy_atmos_comp initialized", ESMF_LOGMSG_INFO)
    print *, "atmos_cap initialize finished, rc =", rc


    call ESMF_GridCompInitialize(dummy_land_comp       , importState=atm2lnd_l_state, exportState=lnd2atm_l_state, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"lnd_cap or dummy_land_comp initialized", ESMF_LOGMSG_INFO)
    print *, "lnd_cap initialize finished, rc =", rc


    ! All 4 states that are module variables are no longer empty - have been initialized

    call ESMF_CplCompInitialize(cpl_atm2lnd_comp, importState=atm2lnd_a_state, exportState=atm2lnd_l_state, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"coupler :: cpl_atm2lnd_comp initialized", ESMF_LOGMSG_INFO)
    print *, "coupler :: cpl_atm2lnd_comp initialize finished, rc =", rc

    call ESMF_CplCompInitialize(cpl_lnd2atm_comp, importState=lnd2atm_l_state, exportState=lnd2atm_a_state, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"coupler :: cpl_lnd2atm_comp initialized", ESMF_LOGMSG_INFO)
    print *, "coupler :: cpl_lnd2atm_comp initialize finished, rc =", rc


  end subroutine lilac_init

  subroutine lilac_run(dum_var1, dum_var2)

    use atmos_cap, only : x2a_fields
    use atmos_cap, only : a2x_fields

    real, dimension(:,:) :: dum_var1  ! from host atm
    real, dimension(:,:) :: dum_var2 ! to host atm

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



end module lilac_mod

