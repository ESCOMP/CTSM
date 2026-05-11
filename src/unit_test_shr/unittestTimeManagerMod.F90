module unittestTimeManagerMod

  ! This module provides wrappers to the clm_time_manager, which facilitate configuring
  ! the time manager as desired for each unit test.
  !
  ! In the setup for a test, the following should be done:
  !
  ! (1) call unittest_timemgr_setup
  !
  ! (2) optionally (if the unit test needs a specific date/time): call
  !     unittest_timemgr_set_curr_date
  !
  ! (3) optionally (if the unit test needs a specific time step number): call
  !     unittest_timemgr_set_nstep
  !
  ! In the teardown for any test that called unittest_timemgr_init, the following should
  ! be done:
  !
  ! (1) call unittest_timemgr_teardown
  !
  !
  ! Note that there are still some test-specific routines in clm_time_manager. Those
  ! include (a) routines that have info that is closely tied to info already in
  ! clm_time_manager (e.g., timemgr_reset, which needs to reset all module data
  ! defined in clm_time_manager), and/or (b) routines that modify data that are private to
  ! clm_time_manager. The routines in this unittest-specific file, in contrast, tend to be
  ! higher-level wrappers.

  use shr_sys_mod, only : shr_sys_abort

  implicit none
  private
  save

  ! Public routines
  public :: unittest_timemgr_setup         ! do the initial setup of the time manager
  public :: unittest_timemgr_set_curr_date ! set the current date
  public :: unittest_timemgr_teardown      ! tear down the time manager at the end of a test
  public :: unittest_timemgr_set_curr_year ! set the current year, keeping other date components unchanged
  public :: unittest_timemgr_set_nstep     ! set the time step number

contains

  !-----------------------------------------------------------------------
  subroutine unittest_timemgr_setup(dtime, use_gregorian_calendar)
    !
    ! !DESCRIPTION:
    ! Set up the time manager for each unit test.
    !
    ! Should be called once for every test that uses the time manager.
    !
    ! !USES:
    use ESMF, only : ESMF_Initialize, ESMF_IsInitialized, ESMF_SUCCESS
    use clm_time_manager, only : set_timemgr_init, timemgr_init, NO_LEAP_C, GREGORIAN_C
    !
    ! !ARGUMENTS:
    integer, intent(in), optional :: dtime  ! time step (seconds)
    logical, intent(in), optional :: use_gregorian_calendar
    !
    ! !LOCAL VARIABLES:
    integer :: l_dtime  ! local version of dtime
    logical :: l_use_gregorian_calendar  ! local version of use_gregorian_calendar
    character(len=:), allocatable :: calendar
    logical :: esmf_is_initialized
    integer :: rc ! return code

    integer, parameter :: dtime_default = 1800  ! time step (seconds)
    
    ! Set ymd values to be year N, month 1, day 1
    integer, parameter :: start_ymd = 10101
    integer, parameter :: ref_ymd = start_ymd
    integer, parameter :: perpetual_ymd = start_ymd

    character(len=*), parameter :: subname = 'unittest_timemgr_setup'
    !-----------------------------------------------------------------------
    
    if (present(dtime)) then
       l_dtime = dtime
    else
       l_dtime = dtime_default
    end if

    if (present(use_gregorian_calendar)) then
       l_use_gregorian_calendar = use_gregorian_calendar
    else
       l_use_gregorian_calendar = .false.
    end if

    esmf_is_initialized = ESMF_IsInitialized(rc=rc)
    if (rc /= ESMF_SUCCESS) then
       call shr_sys_abort(subname//': Error in ESMF_IsInitialized')
    end if
    if (.not. esmf_is_initialized) then
       call ESMF_Initialize(rc=rc)
       if (rc /= ESMF_SUCCESS) then
          call shr_sys_abort(subname//': Error in ESMF_Initialize')
       end if
    end if

    if (l_use_gregorian_calendar) then
       calendar = GREGORIAN_C
    else
       calendar = NO_LEAP_C
    end if

    call set_timemgr_init( &
         calendar_in = calendar, &
         start_ymd_in = start_ymd, &
         start_tod_in = 0, &
         ref_ymd_in = ref_ymd, &
         ref_tod_in = 0, &
         perpetual_run_in = .false., &
         perpetual_ymd_in = perpetual_ymd, &
         dtime_in = l_dtime)

    call timemgr_init()

  end subroutine unittest_timemgr_setup

  !-----------------------------------------------------------------------
  subroutine unittest_timemgr_set_curr_date(yr, mon, day, tod)
    !
    ! !DESCRIPTION:
    ! Set the current model date in the time manager. This is the time at the END of the
    ! time step.
    !
    ! Note that a side effect of this subroutine is that the time step count is
    ! incremented by 1 (because of the method used by for_test_set_curr_date).
    !
    ! !USES:
    use clm_time_manager, only : for_test_set_curr_date
    !
    ! !ARGUMENTS:
    integer, intent(in) :: yr  ! year
    integer, intent(in) :: mon ! month
    integer, intent(in) :: day ! day of month
    integer, intent(in) :: tod ! time of day (seconds past 0Z)
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'unittest_timemgr_set_curr_date'
    !-----------------------------------------------------------------------
    
    call for_test_set_curr_date(yr, mon, day, tod)

  end subroutine unittest_timemgr_set_curr_date

  !-----------------------------------------------------------------------
  subroutine unittest_timemgr_set_curr_year(yr)
    !
    ! !DESCRIPTION:
    ! Set the current model year, keeping other date components unchanged
    !
    ! !USES:
    use clm_time_manager, only : get_curr_date
    !
    ! !ARGUMENTS:
    integer, intent(in) :: yr ! new year
    !
    ! !LOCAL VARIABLES:
    integer :: curr_yr
    integer :: curr_mon
    integer :: curr_day
    integer :: curr_tod

    character(len=*), parameter :: subname = 'unittest_timemgr_set_curr_year'
    !-----------------------------------------------------------------------
    
    call get_curr_date(curr_yr, curr_mon, curr_day, curr_tod)
    call unittest_timemgr_set_curr_date(yr, curr_mon, curr_day, curr_tod)

  end subroutine unittest_timemgr_set_curr_year

  !-----------------------------------------------------------------------
  subroutine unittest_timemgr_set_nstep(nstep)
    !
    ! !DESCRIPTION:
    ! Set the time step number
    !
    ! !USES:
    use clm_time_manager, only : advance_timestep
    !
    ! !ARGUMENTS:
    integer, intent(in) :: nstep
    !
    ! !LOCAL VARIABLES:
    integer :: n

    character(len=*), parameter :: subname = 'unittest_timemgr_set_nstep'
    !-----------------------------------------------------------------------

    do n = 2, nstep
       call advance_timestep()
    end do

  end subroutine unittest_timemgr_set_nstep



  !-----------------------------------------------------------------------
  subroutine unittest_timemgr_teardown
    !
    ! !DESCRIPTION:
    ! Tear down the time manager from each unit test.
    !
    ! Should be called once at the end of every test that set up the time manager.
    !
    ! !USES:
    use ESMF, only : ESMF_Finalize, ESMF_SUCCESS
    use clm_time_manager, only : timemgr_reset
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: rc ! return code
    
    character(len=*), parameter :: subname = 'unittest_timemgr_teardown'
    !-----------------------------------------------------------------------
    
    call timemgr_reset()

    ! If this is the end of the executable, we should call
    ! ESMF_Finalize. But the timemgr setup and teardown routines can
    ! be called multiple times within a single unit test executable,
    ! and it's an error to re-call ESMF_Initialize after calling
    ! ESMF_Finalize. So for now we just won't attempt to do an
    ! ESMF_Finalize.

  end subroutine unittest_timemgr_teardown


end module unittestTimeManagerMod
