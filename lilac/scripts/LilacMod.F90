!Khoda
module LilacMod
use ESMF
use lilac_utils
!use DummyAtmos
use DummyAtmos, only : x2a_fields
use DummyAtmos, only : a2x_fields
use DummyAtmos, only : atmos_register
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

  !public :: lilac_init
  contains

  subroutine lilac_init( dum_var1, dum_var2)
    ! modules
    implicit none

    real, dimension(10) :: dum_var1
    real, dimension(10) :: dum_var2

    ! Component, and State
    type(ESMF_GridComp) :: dummy_atmos_comp  ! the coupled flow Component
    type(ESMF_State)    :: coupledFlowState ! the coupled flow State
    type(ESMF_Mesh)      :: Emesh
    character(len=*), parameter :: subname=trim(modname)//':(lilac_init) '
    type(ESMF_State)     :: importState, exportState
    !character(len=*)    :: atm_mesh_filepath
    
    ! local variables
    integer :: rc, urc
    character(len=ESMF_MAXSTR)            :: cname1, cname2, cplname

    !-------------------------------------------------------------------------
    ! Initialize ESMF, set the default calendar and log type.
    !-------------------------------------------------------------------------
    call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    print *, "----------------------------"
    print *, "lilac Demo Application Start"

    !-------------------------------------------------------------------------
    ! Create Gridded Component!
    !-------------------------------------------------------------------------
    cname1 = "Dummy Atmosphere"
    
    ! Create dummy atmosphere gridded component.
    dummy_atmos_comp = ESMF_GridCompCreate(name=cname1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Created "//trim(cname1)//" component", ESMF_LOGMSG_INFO)
    print *, "Dummy Atmosphere Gridded Component Created!"

    !-------------------------------------------------------------------------
    ! Register section -- set services
    !-------------------------------------------------------------------------
    call ESMF_GridCompSetServices(dummy_atmos_comp, userRoutine=atmos_register, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"dummy atmos SetServices finished!", ESMF_LOGMSG_INFO)
    print *, "Dummy Atmosphere Gridded Component SetServices finished!"

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

  subroutine lilac_run(dum_var_input, dum_var_output)

    use DummyAtmos, only : x2a_fields
    use DummyAtmos, only : a2x_fields

    real, dimension(:) :: dum_var_input ! from host atm
    real, dimension(:) :: dum_var_output ! to host atm

    integer                :: n, num

    integer, parameter     :: fldsMax = 100
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

    dum_var_output(:) = a2x_fields(N)%datafld1d(:)

  end subroutine lilac_run





  !subroutine fldlist_add(num, fldlist, stdname, default_value, units)
  !  integer,                     intent(inout) :: num
  !  type(fld_list_type),         intent(inout) :: fldlist(:)
  !  character(len=*),            intent(in)    :: stdname
  !  real, optional,              intent(in)    :: default_value
  !  character(len=*),  optional, intent(in)    :: units

    ! local variables
  !  integer :: rc
  !  character(len=*), parameter :: subname='(fldlist_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information
 !   num = num + 1
 !   if (num > fldsMax) then
 !      call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=rc) return
 !   endif
 !   fldlist(num)%stdname = trim(stdname)
 !   if(present(default_value)) then
 !      fldlist(num)%default_value = default_value
 !   else
 !      fldlist(num)%default_value = 0.
 !   end if
 !   if(present(units)) then
 !      fldlist(num)%units = trim(units)
 !   else
 !      fldlist(num)%units = ""
 !   end if

 ! end subroutine fldlist_add


subroutine fldlist_add_dumb(num, fldlist, stdname)
    integer,                    intent(inout) :: num
    type(fld_list_type),        intent(inout) :: fldlist(:)
    character(len=*),           intent(in)    :: stdname

    ! local variables
    integer :: rc
    integer :: dbrc
    character(len=*), parameter :: subname='(lnd_import_export:fldlist_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
       call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=dbrc)
       return
    endif
    fldlist(num)%stdname = trim(stdname)

  end subroutine fldlist_add_dumb
 subroutine create_fldlists_dumb(fldsFrCpl, fldsToCpl, fldsToCpl_num, fldsFrCpl_num)
    type(fld_list_type),        intent(inout) :: fldsFrCpl(:)
    type(fld_list_type),        intent(inout) :: fldsToCpl(:)
    !integer, intent(out)                     :: fldsToCpl_num = 0
    !integer, intent(out)                     :: fldsFrCpl_num = 0

    ! import fields
    ! call fldlist_add(fldsFrCpl_num, fldsFrCpl, trim(flds_scalar_name))

    integer :: fldsFrCpl_num, fldsToCpl_num

    ! land states
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_lfrin'      )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_t'          )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_tref'       )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_qref'       )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_avsdr'      )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_anidr'      )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_avsdf'      )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_anidf'      )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_snowh'      )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_u10'        )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_fv'         )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_ram1'       )

    ! fluxes to atm
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_taux'     )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_tauy'     )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_lat'      )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_sen'      )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_lwup'     )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_evap'     )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_swnet'    )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_flxdst1'  )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_flxdst2'  )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_flxdst3'  )
    call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_flxdst4'  )

    ! call fldlist_add(fldsToCpl_num, fldsToCpl, trim(flds_scalar_name))

    ! from atm
    call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_z', default_value=30.0, units='m')
    ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_topo')
    call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_u', default_value=0.0, units='m/s')
    call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_v', default_value=0.0, units='m/s')
    call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_ptem', default_value=280.0, units='degK')
    call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_pbot', default_value=100100.0, units='Pa')
    call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_tbot', default_value=280.0, units='degK')
    call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_shum', default_value=0.0004, units='kg/kg')
    !call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_methane'   )

    call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_lwdn', default_value=200.0, units='W/m2')
    call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_rainc', default_value=4.0e-8, units='kg/m2s')
    call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_rainl', default_value=3.0e-8, units='kg/m2s')
    call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_snowc', default_value=1.0e-8, units='kg/m2s')
    call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_snowl', default_value=2.0e-8, units='kg/m2s')
    call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_swndr', default_value=100.0, units='W/m2')
    call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_swvdr', default_value=90.0, units='W/m2')
    call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_swndf', default_value=20.0, units='W/m2')
    call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_swvdf', default_value=40.0, units='W/m2')
    ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_bcphidry')
    ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_bcphodry')
    ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_bcphiwet')
    ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_ocphidry')
    ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_ocphodry')
    ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_ocphiwet')
    ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstdry1' )
    ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstdry2' )
    ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstdry3' )
    ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstdry4' )
    ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstwet1' )
    ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstwet2' )
    ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstwet3' )
    ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstwet4' )

    ! more: https://github.com/mvertens/ctsm/blob/ae02ffe25dbc4a85c769c9137b5b3d50f2843e89/src/cpl/nuopc/lnd_import_export.F90#L131
  end subroutine create_fldlists_dumb

end module LilacMod

