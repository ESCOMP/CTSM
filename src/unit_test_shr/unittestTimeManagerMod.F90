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

  implicit none
  private
  save

  ! Public routines
  public :: unittest_timemgr_setup         ! do the initial setup of the time manager
  public :: unittest_timemgr_set_curr_date ! set the current date
  public :: unittest_timemgr_teardown      ! tear down the time manager at the end of a test
  public :: unittest_timemgr_set_curr_year ! set the current year, keeping other date components unchanged

contains

  !-----------------------------------------------------------------------
  subroutine unittest_timemgr_setup(dtime)
    !
    ! !DESCRIPTION:
    ! Set up the time manager for each unit test.
    !
    ! Should be called once for every test that uses the time manager.
    !
    ! !USES:
    use ESMF, only : ESMF_Initialize, ESMF_SUCCESS
    use clm_time_manager, only : set_timemgr_init, timemgr_init, NO_LEAP_C
    !
    ! !ARGUMENTS:
    integer, intent(in), optional :: dtime  ! time step (seconds)
    !
    ! !LOCAL VARIABLES:
    integer :: l_dtime  ! local version of dtime
    integer :: rc ! return code

    integer, parameter :: dtime_default = 1800  ! time step (seconds)
    
    ! Set ymd values to be year N, month 1, day 1
    integer, parameter :: start_ymd = 10101
    integer, parameter :: ref_ymd = start_ymd
    integer, parameter :: stop_ymd = 20101
    integer, parameter :: perpetual_ymd = start_ymd

    ! Set current time to be at the start of year 1
    integer, parameter :: curr_yr = 1
    integer, parameter :: curr_mon = 1
    integer, parameter :: curr_day = 1
    integer, parameter :: curr_tod = 0

    character(len=*), parameter :: subname = 'unittest_timemgr_setup'
    !-----------------------------------------------------------------------
    
    if (present(dtime)) then
       l_dtime = dtime
    else
       l_dtime = dtime_default
    end if

    call ESMF_Initialize(rc=rc)
    if (rc /= ESMF_SUCCESS) then
       stop 'Error in ESMF_Initialize'
    end if

    call set_timemgr_init( &
         calendar_in = NO_LEAP_C, &
         start_ymd_in = start_ymd, &
         start_tod_in = 0, &
         ref_ymd_in = ref_ymd, &
         ref_tod_in = 0, &
         stop_ymd_in = stop_ymd, &
         stop_tod_in = 0, &
         perpetual_run_in = .false., &
         perpetual_ymd_in = perpetual_ymd, &
         nelapse_in = 1, &
         dtime_in = l_dtime)

    call timemgr_init()

    call unittest_timemgr_set_curr_date( &
         yr = curr_yr, &
         mon = curr_mon, &
         day = curr_day, &
         tod = curr_tod)

  end subroutine unittest_timemgr_setup

  !-----------------------------------------------------------------------
  subroutine unittest_timemgr_set_curr_date(yr, mon, day, tod)
    !
    ! !DESCRIPTION:
    ! Set the current model date in the time manager.
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

    call ESMF_Finalize(rc=rc)
    if (rc /= ESMF_SUCCESS) then
       stop 'Error in ESMF_Finalize'
    end if

  end subroutine unittest_timemgr_teardown


end module unittestTimeManagerMod
