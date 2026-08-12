module clm_time_manager

#include "shr_assert.h"
   use shr_kind_mod, only: r8 => shr_kind_r8
   use shr_sys_mod , only: shr_sys_abort
   use spmdMod     , only: masterproc
   use clm_varctl  , only: iulog
   use clm_varcon  , only: isecspday
   use ESMF        , only: ESMF_Clock, ESMF_Calendar, ESMF_MAXSTR, ESMF_Time, ESMF_TimeInterval
   use ESMF        , only: ESMF_TimeIntervalSet, ESMF_TimeSet, ESMF_TimeGet, ESMF_ClockGet
   use ESMF        , only: operator(==), operator(+), operator(<=), operator(>=)
   use ESMF        , only: operator(>), operator(<), operator(-)
   use ESMF        , only: ESMF_KIND_I8, ESMF_TimeIntervalGet

   implicit none
   private

   ! Public methods

   public ::&
        set_timemgr_init,         &! setup startup values
        timemgr_init,             &! time manager initialization, called always
        timemgr_restart_io,       &! read/write time manager restart info and restart time manager
        timemgr_restart,          &! check that time manager is setup coorectly upcon restart
        timemgr_datediff,         &! calculate difference between two time instants
        advance_timestep,         &! increment timestep number
        get_curr_ESMF_Time,       &! get current time in terms of the ESMF_Time
        get_step_size,            &! return step size in seconds
        get_step_size_real,       &! return step size in seconds, real-valued
        get_rad_step_size,        &! return radiation step size in seconds
        get_nstep,                &! return timestep number
        get_nstep_since_startup_or_lastDA_restart_or_pause, &! return number of timesteps since restart was modified
        get_curr_date,            &! return date components at end of current timestep
        get_prev_date,            &! return date components at beginning of current timestep
        get_start_date,           &! return components of the start date
        get_driver_start_ymd,     &! return year/month/day (as integer in YYYYMMDD format) of driver start date
        get_ref_date,             &! return components of the reference date
        get_perp_date,            &! return components of the perpetual date, and current time of day
        get_curr_time,            &! return components of elapsed time since reference date at end of current timestep
        get_prev_time,            &! return components of elapsed time since reference date at beg of current timestep
        get_curr_calday,          &! return calendar day at end of current timestep
        get_prev_calday,          &! return calendar day at beginning of current timestep
        get_calday,               &! return calendar day from input date
        get_calendar,             &! return calendar
        get_doy_tomorrow,         &! return next day of year
        get_average_days_per_year,&! return the average number of days per year for the given calendar
        get_curr_days_per_year,   &! return the days per year for year as of the end of the current time step
        get_prev_days_per_year,   &! return the days per year for year as of the beginning of the current time step
        get_curr_yearfrac,        &! return the fractional position in the current year, as of the end of the current timestep
        get_prev_yearfrac,        &! return the fractional position in the current year, as of the beginning of the current timestep
        get_rest_date,            &! return the date from the restart file
        get_local_timestep_time,  &! return the local time for the input longitude to the nearest time-step
        get_local_time,           &! return the local time for the input longitude
        is_first_step,            &! return true on first step of initial run
        is_first_restart_step,    &! return true on first step of restart or branch run
        is_first_step_of_this_run_segment, &! return true on first step of any run segment (initial, restart or branch run)
        is_beg_curr_day,          &! return true on first timestep in current day
        is_end_curr_day,          &! return true on last timestep in current day
        is_end_curr_month,        &! return true on last timestep in current month
        is_beg_curr_year,         &! return true on first timestep in current year
        is_end_curr_year,         &! return true on last timestep in current year
        is_perpetual,             &! return true if perpetual calendar is in use
        is_doy_in_interval,       &! return true if day of year is in the provided interval
        is_today_in_doy_interval, &! return true if today's day of year is in the provided interval
        is_near_local_noon,       &! return true if near local noon
        is_restart,               &! return true if this is a restart run
        update_rad_dtime,         &! track radiation interval via nstep
        update_DA_nstep,          &! update the Data Assimulation time step
        timemgr_reset              ! reset values to their defaults, and free memory

   ! Public methods, but just to support unit testing:
   public :: for_test_set_curr_date  ! set the current date and time

   ! Public parameter data
   character(len=*), public, parameter :: NO_LEAP_C   = 'NO_LEAP'
   character(len=*), public, parameter :: GREGORIAN_C = 'GREGORIAN'

   ! Private module data

   ! Private data for input

   character(len=ESMF_MAXSTR), save ::&
        calendar   = NO_LEAP_C        ! Calendar to use in date calculations.
   integer,  parameter :: uninit_int = -999999999
   real(r8), parameter :: uninit_r8  = -999999999.0

   ! We'll use this really big year to effectively mean infinitely into the future.
   integer,  parameter :: really_big_year = 999999999

   ! Input
   integer, save ::&
        dtime          = uninit_int,  &! timestep in seconds
        dtime_rad      = uninit_int,  &! radiation interval in seconds
        nstep_rad_prev = uninit_int    ! radiation interval in seconds

   ! Input from CESM driver
   integer, save ::&
        start_ymd     = uninit_int,  &! starting date for run in yearmmdd format
        start_tod     = 0,           &! starting time of day for run in seconds
        ref_ymd       = uninit_int,  &! reference date for time coordinate in yearmmdd format
        ref_tod       = 0             ! reference time of day for time coordinate in seconds
   type(ESMF_Calendar), target, save   :: tm_cal       ! calendar
   type(ESMF_Clock),    save   :: tm_clock     ! model clock
   integer,             save   :: tm_clock_step_size_sec ! Cache of clock timestep.
   type(ESMF_Time),     save   :: tm_perp_date ! perpetual date

   ! Data required to restart time manager (only set if timemgr_restart_io is called):
   integer, save :: rst_step_sec          = uninit_int ! timestep size seconds
   integer, save :: rst_start_ymd         = uninit_int ! start date
   integer, save :: rst_start_tod         = uninit_int ! start time of day
   integer, save :: rst_ref_ymd           = uninit_int ! reference date
   integer, save :: rst_ref_tod           = uninit_int ! reference time of day
   integer, save :: rst_curr_ymd          = uninit_int ! current date
   integer, save :: rst_curr_tod          = uninit_int ! current time of day

   integer, save :: rst_nstep_rad_prev                 ! nstep of previous radiation call
   integer, save :: perpetual_ymd         = uninit_int ! Perpetual calendar date (YYYYMMDD)
   logical, save :: tm_first_restart_step = .false.    ! true for first step of a restart or branch run
   logical, save :: tm_perp_calendar      = .false.    ! true when using perpetual calendar
   logical, save :: timemgr_set           = .false.    ! true when timemgr initialized

   !
   ! The time-step number of startup or last Data Assimulation (DA) restart or pause
   integer, save :: DA_nstep = 0 ! Last step number that state was modified externally (by DA)

   ! Private module methods

   private :: timemgr_spmdbcast
   private :: init_calendar
   private :: init_clock
   private :: timemgr_print
   private :: TimeGetymd
   private :: check_timemgr_initialized

   character(len=*), parameter, private :: sourcefile = &
        __FILE__

   !=========================================================================================
contains
  !=========================================================================================

  subroutine set_timemgr_init( calendar_in,      start_ymd_in,     start_tod_in, ref_ymd_in,        &
       ref_tod_in, perpetual_run_in, perpetual_ymd_in, dtime_in )

    !---------------------------------------------------------------------------------
    ! set time manager startup values
    ! 
    ! Arguments
    character(len=*), optional, intent(IN) :: calendar_in       ! Calendar type
    integer         , optional, intent(IN) :: start_ymd_in      ! Start date       (YYYYMMDD)
    integer         , optional, intent(IN) :: start_tod_in      ! Start time of day (sec)
    integer         , optional, intent(IN) :: ref_ymd_in        ! Reference date   (YYYYMMDD)
    integer         , optional, intent(IN) :: ref_tod_in        ! Reference time of day (sec)
    logical         , optional, intent(IN) :: perpetual_run_in  ! If in perpetual mode or not
    integer         , optional, intent(IN) :: perpetual_ymd_in  ! Perpetual date   (YYYYMMDD)
    integer         , optional, intent(IN) :: dtime_in          ! Time-step (sec)
    !
    character(len=*), parameter :: sub = 'clm::set_timemgr_init'

    if ( timemgr_set ) call shr_sys_abort( sub//":: timemgr_init already called" )
    if (present(calendar_in)      ) calendar         = trim(calendar_in)
    if (present(start_ymd_in)     ) start_ymd        = start_ymd_in
    if (present(start_tod_in)     ) start_tod        = start_tod_in
    if (present(ref_ymd_in)       ) ref_ymd          = ref_ymd_in
    if (present(ref_tod_in)       ) ref_tod          = ref_tod_in
    if (present(perpetual_run_in) )then
       tm_perp_calendar = perpetual_run_in
       if ( tm_perp_calendar ) then
          if ( .not. present(perpetual_ymd_in) .or. perpetual_ymd == uninit_int) &
               call shr_sys_abort( sub//":: perpetual_run set but NOT perpetual_ymd" )
          perpetual_ymd    = perpetual_ymd_in
       end if
    end if
    if (present(dtime_in)         ) dtime            = dtime_in

  end subroutine set_timemgr_init

  !=========================================================================================

  subroutine timemgr_init(curr_date_in )

    use clm_varctl, only : nsrest, nsrContinue, nsrBranch

    type(ESMF_Time), intent(in), optional :: curr_date_in 

    !---------------------------------------------------------------------------------
    ! Initialize the ESMF time manager from the sync clock
    ! 
    ! Arguments
    !
    character(len=*), parameter :: sub = 'clm::timemgr_init'
    integer :: rc                            ! return code
    type(ESMF_Time) :: start_date            ! start date for run
    type(ESMF_Time) :: ref_date              ! reference date for time coordinate
    type(ESMF_Time) :: curr_date             ! temporary date used in logic
    type(ESMF_TimeInterval) :: day_step_size ! day step size
    type(ESMF_TimeInterval) :: step_size     ! timestep size
    !---------------------------------------------------------------------------------
    call timemgr_spmdbcast( )

    ! Initalize calendar 

    call init_calendar()

    ! Initalize start date.

    if ( start_ymd == uninit_int ) then
       write(iulog,*)sub,': start_ymd must be specified '
       call shr_sys_abort
    end if
    if ( start_tod == uninit_int ) then
       write(iulog,*)sub,': start_tod must be specified '
       call shr_sys_abort
    end if
    start_date = TimeSetymd( start_ymd, start_tod, "start_date" )

    ! Initialize current date
    if(present(curr_date_in)) then
       curr_date = curr_date_in
    else
       curr_date = start_date
    endif

    call ESMF_TimeIntervalSet( step_size, s=dtime, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet: setting step_size')

    call ESMF_TimeIntervalSet( day_step_size, d=1, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet: setting day_step_size')

    ! Initalize reference date for time coordinate.

    if ( ref_ymd /= uninit_int ) then
       ref_date = TimeSetymd( ref_ymd, ref_tod, "ref_date" )
    else
       ref_date = start_date
    end if

    ! Initialize clock

    call init_clock( start_date, ref_date, curr_date)

    ! Initialize date used for perpetual calendar day calculation.

    if (tm_perp_calendar) then
       tm_perp_date = TimeSetymd( perpetual_ymd, 0, "tm_perp_date" )
    end if

    ! Advance time step to start at nstep=1
    if (nsrest /= nsrContinue .and. nsrest /= nsrBranch) then
       call advance_timestep()
    end if

    ! Print configuration summary to log file (stdout).

    if (masterproc) call timemgr_print()

    timemgr_set = .true.

  end subroutine timemgr_init

  !=========================================================================================

  subroutine init_clock( start_date, ref_date, curr_date )

    !---------------------------------------------------------------------------------
    ! Purpose: Initialize the clock based on the start_date, ref_date and curr_date
    !
    use ESMF       , only : ESMF_ClockCreate, ESMF_ClockAdvance, esmf_clockiscreated

    type(ESMF_Time), intent(in) :: start_date  ! start date for run
    type(ESMF_Time), intent(in) :: ref_date    ! reference date for time coordinate
    type(ESMF_Time), intent(in) :: curr_date   ! current date (equal to start_date)
    !
    character(len=*), parameter :: sub = 'clm::init_clock'
    type(ESMF_Time)             :: stop_date         ! stop date for run
    type(ESMF_TimeInterval)     :: step_size         ! timestep size
    type(ESMF_Time)             :: current           ! current date (from clock)
    integer                     :: yr, mon, day, tod ! Year, month, day, and second as integers
    integer                     :: rc                ! return code
    !---------------------------------------------------------------------------------

    call ESMF_TimeIntervalSet( step_size, s=dtime, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet: setting step_size')

    ! We don't use a stop time in the CTSM clock. Instead, we set the clock to
    ! effectively have a stop time infinitely far into the future, and rely on other
    ! mechanisms to tell CTSM when to stop. If we were always using the real ESMF
    ! library, we could avoid setting the stopTime on the clock. But the ESMF time
    ! manager included in cime appears to require stopTime.
    call ESMF_TimeSet(stop_date, yy=really_big_year, mm=12, dd=31, s=0, &
         calendar=tm_cal, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet: setting step_size')

    ! Error check 

    if ( stop_date <= start_date ) then
       write(iulog,*)sub, ': Assumed stop date is earlier than start date: '
       call ESMF_TimeGet( start_date, yy=yr, mm=mon, dd=day, s=tod )
       write(iulog,*) ' Start date (yr, mon, day, tod): ', yr, mon, day, tod
       call ESMF_TimeGet( stop_date, yy=yr, mm=mon, dd=day, s=tod )
       write(iulog,*) ' Stop date  (yr, mon, day, tod): ', yr, mon, day, tod
       call shr_sys_abort
    end if
    if ( stop_date <= curr_date ) then
       write(iulog,*)sub, ': Assumed stop date is earlier than current date: '
       call ESMF_TimeGet( curr_date, yy=yr, mm=mon, dd=day, s=tod )
       write(iulog,*) ' Current date (yr, mon, day, tod): ', yr, mon, day, tod
       call ESMF_TimeGet( stop_date, yy=yr, mm=mon, dd=day, s=tod )
       write(iulog,*) ' Stop date    (yr, mon, day, tod): ', yr, mon, day, tod
       call shr_sys_abort
    end if

    ! Initialize the clock

    
    tm_clock = ESMF_ClockCreate(name="CLM Time-manager clock", timeStep=step_size, startTime=start_date, &
         stopTime=stop_date, refTime=ref_date, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_ClockCreate')

    ! Advance clock to the current time (in case of a restart)

    call ESMF_ClockGet(tm_clock, currTime=current, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_ClockGet')
    do while( curr_date > current )
       call ESMF_ClockAdvance( tm_clock, rc=rc )
       call chkrc(rc, sub//': error return from ESMF_ClockAdvance')
       call ESMF_ClockGet(tm_clock, currTime=current )
       call chkrc(rc, sub//': error return from ESMF_ClockGet')
    end do


    ! Cache step size, we query it a lot.
    call ESMF_TimeIntervalGet(step_size, s=tm_clock_step_size_sec, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_ClockTimeIntervalGet')

  end subroutine init_clock

  !=========================================================================================

  function TimeSetymd( ymd, tod, desc )
    !---------------------------------------------------------------------------------
    !
    ! Set the time by an integer as YYYYMMDD and integer seconds in the day
    !
    integer, intent(in) :: ymd            ! Year, month, day YYYYMMDD
    integer, intent(in) :: tod            ! Time of day in seconds
    character(len=*), intent(in) :: desc  ! Description of time to set

    type(ESMF_Time) :: TimeSetymd    ! Return value

    character(len=*), parameter :: sub = 'clm::TimeSetymd'
    integer :: yr, mon, day          ! Year, month, day as integers
    integer :: rc                    ! return code
    !---------------------------------------------------------------------------------

    if ( (ymd < 0) .or. (tod < 0) .or. (tod > isecspday) )then
       write(iulog,*) sub//': error yymmdd is a negative number or time-of-day out of bounds', &
            ymd, tod
       call shr_sys_abort
    end if
    yr  = ymd / 10000
    mon = (ymd - yr*10000) / 100
    day =  ymd - yr*10000 - mon*100
    call ESMF_TimeSet( TimeSetymd, yy=yr, mm=mon, dd=day, s=tod, &
         calendar=tm_cal, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_TimeSet: setting '//trim(desc))
  end function TimeSetymd

  !=========================================================================================

  integer function TimeGetymd( date, tod )
    !
    ! Get the date and time of day in ymd from ESMF Time.
    !
    type(ESMF_Time), intent(inout) :: date ! Input date to convert to ymd
    integer, intent(out), optional :: tod  ! Time of day in seconds

    character(len=*), parameter :: sub = 'clm::TimeGetymd'
    integer :: yr, mon, day
    integer :: rc                          ! return code

    call ESMF_TimeGet( date, yy=yr, mm=mon, dd=day, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_TimeGet')
    TimeGetymd = yr*10000 + mon*100 + day
    if ( present( tod ) )then
       call ESMF_TimeGet( date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
       call chkrc(rc, sub//': error return from ESMF_TimeGet')
    end if
    if ( yr < 0 )then
       write(iulog,*) sub//': error year is less than zero', yr
       call shr_sys_abort
    end if
  end function TimeGetymd

  !=========================================================================================

  subroutine timemgr_restart_io( ncid, flag )

    !---------------------------------------------------------------------------------
    ! Read/Write information needed on restart to a netcdf file. 
    use ncdio_pio, only: ncd_int, file_desc_t
    use restUtilMod
    !
    ! Arguments
    type(file_desc_t), intent(inout) :: ncid  ! netcdf id
    character(len=*), intent(in) :: flag  ! 'read' or 'write'
    !
    ! Local variables
    character(len=*), parameter :: sub = 'clm::timemgr_restart'
    integer :: rc                  ! return code
    logical :: readvar             ! determine if variable is on initial file
    type(ESMF_Time) :: start_date  ! start date for run
    type(ESMF_Time) :: ref_date    ! reference date for run
    type(ESMF_Time) :: curr_date   ! date of data in restart file
    integer :: rst_caltype         ! calendar type
    integer, parameter :: noleap = 1
    integer, parameter :: gregorian = 2
    character(len=len(calendar)) :: cal
    !---------------------------------------------------------------------------------

    if (flag == 'write') then
       rst_nstep_rad_prev  = nstep_rad_prev
    end if
    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_nstep_rad_prev', xtype=ncd_int,  &
         long_name='previous_radiation_nstep', units='unitless positive integer', &
         ifill_value=uninit_int, &
         interpinic_flag='skip', readvar=readvar, data=rst_nstep_rad_prev)
    if (flag == 'read') then
       nstep_rad_prev = rst_nstep_rad_prev
    end if

    if (flag == 'write') then
       cal = to_upper(calendar)
       if ( trim(cal) == NO_LEAP_C ) then
          rst_caltype = noleap
       else if ( trim(cal) == GREGORIAN_C ) then
          rst_caltype = gregorian
       else
          call shr_sys_abort(sub//'ERROR: unrecognized calendar specified= '//trim(calendar))
       end if
    end if
    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_type', xtype=ncd_int,  &
         long_name='calendar type', units='unitless', flag_meanings=(/ "NO_LEAP_C", "GREGORIAN" /), &
         flag_values=(/ noleap, gregorian /), ifill_value=uninit_int, &
         interpinic_flag='skip', readvar=readvar, data=rst_caltype)
    if (flag == 'read') then
       if ( rst_caltype == noleap ) then
          calendar = NO_LEAP_C
       else if ( rst_caltype == gregorian ) then
          calendar = GREGORIAN_C
       else
          write(iulog,*)sub,': unrecognized calendar type in restart file: ',rst_caltype
          call shr_sys_abort( sub//'ERROR: bad calendar type in restart file')
       end if
    end if

    if (flag == 'write') then
       call ESMF_ClockGet( tm_clock, startTime=start_date, currTime=curr_date, refTime=ref_date, rc=rc )
       call chkrc(rc, sub//': error return from ESMF_ClockGet')
       rst_step_sec  = dtime
       rst_start_ymd = TimeGetymd( start_date, tod=rst_start_tod )
       rst_ref_ymd   = TimeGetymd( ref_date,   tod=rst_ref_tod   )
       rst_curr_ymd  = TimeGetymd( curr_date,  tod=rst_curr_tod  )
    end if
    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_step_sec', xtype=ncd_int, &
         long_name='seconds component of timestep size', units='sec',         &
         nvalid_range=(/0,isecspday/), ifill_value=uninit_int,                &
         interpinic_flag='skip', readvar=readvar, data=rst_step_sec)
    if ((flag == 'read') .and. ( rst_step_sec < 0 .or. rst_step_sec > isecspday )) then
       call shr_sys_abort( sub//'ERROR: timemgr_rst_step_sec out of range')
    end if

    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_start_ymd', xtype=ncd_int, &
         long_name='start date', units='YYYYMMDD', ifill_value=uninit_int,     &
         interpinic_flag='skip', readvar=readvar, data=rst_start_ymd)

    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_start_tod', xtype=ncd_int, &
         long_name='start time of day', units='sec',                           &
         nvalid_range=(/0,isecspday/), ifill_value=uninit_int,                 &
         interpinic_flag='skip', readvar=readvar, data=rst_start_tod)
    if ((flag == 'read') .and. ( rst_start_tod < 0 .or. rst_start_tod > isecspday )) then
       call shr_sys_abort( sub//'ERROR: timemgr_rst_strart_tod out of range')
    end if

    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_ref_ymd', xtype=ncd_int,   &
         long_name='reference date', units='YYYYMMDD', ifill_value=uninit_int, &
         interpinic_flag='skip', readvar=readvar, data=rst_ref_ymd)

    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_ref_tod', xtype=ncd_int,   &
         long_name='reference time of day', units='sec',                       &
         nvalid_range=(/0,isecspday/), ifill_value=uninit_int,                 &
         interpinic_flag='skip', readvar=readvar, data=rst_ref_tod)
    if ((flag == 'read') .and. ( rst_start_tod < 0 .or. rst_start_tod > isecspday )) then
       call shr_sys_abort( sub//'ERROR: timemgr_rst_ref_tod out of range')
    end if

    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_curr_ymd', xtype=ncd_int,  &
         long_name='current date', units='YYYYMMDD', ifill_value=uninit_int,   &
         interpinic_flag='skip', readvar=readvar, data=rst_curr_ymd)

    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_curr_tod', xtype=ncd_int,  &
         long_name='current time of day', units='sec',                         &
         nvalid_range=(/0,isecspday/), ifill_value=uninit_int,                 &
         interpinic_flag='skip', readvar=readvar, data=rst_curr_tod)
    if ((flag == 'read') .and. ( rst_curr_tod < 0 .or. rst_curr_tod > isecspday )) then
       call shr_sys_abort( sub//'ERROR: timemgr_rst_ref_ymd out of range')
    end if

  end subroutine timemgr_restart_io

  !=========================================================================================

  subroutine timemgr_restart()

    !---------------------------------------------------------------------------------
    ! On restart do some checkcing to make sure time is synchronized with the clock from CESM.
    ! Set a couple of variables, and advance the clock, so time is aligned properly.
   !
    ! timemgr_init MIST be called before this
    !

    character(len=*), parameter :: sub = 'clm::timemgr_restart'
    integer :: rc                            ! return code
    integer :: yr, mon, day, tod             ! Year, month, day, and second as integers
    type(ESMF_Time) :: start_date            ! start date for run
    type(ESMF_Time) :: ref_date              ! reference date for run
    type(ESMF_Time) :: curr_date             ! date of data in restart file
    type(ESMF_TimeInterval) :: day_step_size ! day step size
    type(ESMF_TimeInterval) :: step_size     ! timestep size
    !---------------------------------------------------------------------------------
    ! Check that timemgr_init was already called
    if ( .not. check_timemgr_initialized(sub) ) return

    ! Initialize the timestep

    dtime = rst_step_sec

    ! Check start date from restart info

    if (rst_start_ymd .ne. start_ymd .or. rst_start_tod .ne. start_tod) then
       call shr_sys_abort(sub//'ERROR: mismatch in start date with restart file')
    endif

    if (rst_ref_ymd .ne. ref_ymd .or. rst_ref_tod .ne. ref_tod) then
       call shr_sys_abort(sub//'ERROR: mismatch in reference date with restart file')
    endif

    call ESMF_TimeIntervalSet( day_step_size, d=1, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet: setting day_step_size')

    ! Initialize nstep_rad_prev from restart info

    nstep_rad_prev = rst_nstep_rad_prev

    ! Initialize ref date from restart info

    ! Advance the timestep.  
    ! Data from the restart file corresponds to the last timestep of the previous run.

    call advance_timestep()

    ! Set flag that this is the first timestep of the restart run.

    tm_first_restart_step = .true.

  end subroutine timemgr_restart

  !=========================================================================================

  subroutine init_calendar( )

    !---------------------------------------------------------------------------------
    ! Initialize calendar
    use ESMF        , only : ESMF_CalKind_Flag, ESMF_CALKIND_NOLEAP
    use ESMF        , only : ESMF_CALKIND_GREGORIAN, ESMF_CalendarCreate
    !
    ! Local variables
    !
    character(len=*), parameter :: sub = 'clm::init_calendar'
    type(ESMF_CalKind_Flag) :: cal_type        ! calendar type
    character(len=len(calendar)) :: caltmp
    integer :: rc                              ! return code
    !---------------------------------------------------------------------------------

    caltmp = to_upper(calendar)
    if ( trim(caltmp) == NO_LEAP_C ) then
       cal_type = ESMF_CALKIND_NOLEAP
    else if ( trim(caltmp) == GREGORIAN_C ) then
       cal_type = ESMF_CALKIND_GREGORIAN
    else
       write(iulog,*)sub,': unrecognized calendar specified: ',calendar
       call shr_sys_abort
    end if
    tm_cal = ESMF_CalendarCreate( name=caltmp, calkindflag=cal_type, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_CalendarSet')
  end subroutine init_calendar

  !=========================================================================================

  subroutine timemgr_print()

    !---------------------------------------------------------------------------------
    character(len=*), parameter :: sub = 'clm::timemgr_print'
    integer :: rc
    integer :: yr, mon, day
    integer :: &                   ! Data required to restart time manager:
         nstep     = uninit_int,  &! current step number
         step_sec  = uninit_int,  &! timestep size seconds
         start_yr  = uninit_int,  &! start year
         start_mon = uninit_int,  &! start month
         start_day = uninit_int,  &! start day of month
         start_tod = uninit_int,  &! start time of day
         ref_yr    = uninit_int,  &! reference year
         ref_mon   = uninit_int,  &! reference month
         ref_day   = uninit_int,  &! reference day of month
         ref_tod   = uninit_int,  &! reference time of day
         curr_yr   = uninit_int,  &! current year
         curr_mon  = uninit_int,  &! current month
         curr_day  = uninit_int,  &! current day of month
         curr_tod  = uninit_int    ! current time of day
    integer(ESMF_KIND_I8) :: step_no
    type(ESMF_Time) :: start_date! start date for run
    type(ESMF_Time) :: curr_date ! date of data in restart file
    type(ESMF_Time) :: ref_date  ! reference date
    type(ESMF_TimeInterval) :: step ! Time-step
    !---------------------------------------------------------------------------------

    call ESMF_ClockGet( tm_clock, startTime=start_date, currTime=curr_date, &
         refTime=ref_date, timeStep=step, &
         advanceCount=step_no, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_ClockGet')
    nstep = step_no

    write(iulog,*)' ******** CLM Time Manager Configuration ********'

    call ESMF_TimeIntervalGet( step, s=step_sec, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_TimeIntervalGet')

    call ESMF_TimeGet( start_date, yy=start_yr, mm=start_mon, dd=start_day, &
         s=start_tod, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_TimeGet')
    call ESMF_TimeGet( ref_date, yy=ref_yr, mm=ref_mon, dd=ref_day, s=ref_tod, &
         rc=rc )
    call chkrc(rc, sub//': error return from ESMF_TimeGet')
    call ESMF_TimeGet( curr_date, yy=curr_yr, mm=curr_mon, dd=curr_day, &
         s=curr_tod, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_TimeGet')

    write(iulog,*)'  Calendar type:            ',trim(calendar)
    write(iulog,*)'  Timestep size (seconds):  ', step_sec
    write(iulog,*)'  Start date (yr mon day tod):     ', start_yr, start_mon, &
         start_day, start_tod
    write(iulog,*)'  Reference date (yr mon day tod): ', ref_yr, ref_mon, &
         ref_day, ref_tod
    write(iulog,*)'  Current step number:      ', nstep
    write(iulog,*)'  Current date (yr mon day tod):   ', curr_yr, curr_mon, &
         curr_day, curr_tod

    if ( tm_perp_calendar ) then
       call ESMF_TimeGet( tm_perp_date, yy=yr, mm=mon, dd=day, rc=rc )
       call chkrc(rc, sub//': error return from ESMF_TimeGet')
       write(iulog,*)'  Use perpetual diurnal cycle date (yr mon day): ', &
            yr, mon, day
    end if

    write(iulog,*)' ************************************************'

  end subroutine timemgr_print

  !=========================================================================================

  subroutine advance_timestep()

    ! Increment the timestep number.
    use ESMF       , only : ESMF_ClockAdvance

    character(len=*), parameter :: sub = 'clm::advance_timestep'
    integer :: rc

    call ESMF_ClockAdvance( tm_clock, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_ClockAdvance')

    tm_first_restart_step = .false.

  end subroutine advance_timestep

  !=========================================================================================

  function get_curr_ESMF_Time( )

    ! Return the current time as ESMF_Time

    type(ESMF_Time) :: get_curr_ESMF_Time
    character(len=*), parameter :: sub = 'clm::get_curr_ESMF_Time'
    integer :: rc

    if ( .not. check_timemgr_initialized(sub) ) return

    call ESMF_ClockGet( tm_clock, currTime=get_curr_ESMF_Time, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_ClockGet')

  end function get_curr_ESMF_Time

  !=========================================================================================

  integer function get_step_size()

    ! Return the step size in seconds.

    character(len=*), parameter :: sub = 'clm::get_step_size'

    if ( .not. check_timemgr_initialized(sub) ) return

    get_step_size = tm_clock_step_size_sec

  end function get_step_size

  !=========================================================================================

  real(r8) function get_step_size_real()

    ! Return the step size in seconds, as a real value

    get_step_size_real = real(get_step_size(), r8)

  end function get_step_size_real

  !=========================================================================================

  subroutine update_DA_nstep()
     ! Update the Data Assimulation time-step to the current time step, since DA has been done
     DA_nstep = get_nstep()
  end subroutine update_DA_nstep

  !=========================================================================================

  subroutine update_rad_dtime(doalb)
    !---------------------------------------------------------------------------------
    ! called only on doalb timesteps to save off radiation nsteps
    ! 
    ! Local Arguments
    logical,intent(in) ::  doalb
    integer :: dtime,nstep

    if (doalb) then 

       dtime=get_step_size()
       nstep = get_nstep()

       if (nstep_rad_prev == uninit_int ) then
          dtime_rad = dtime
          nstep_rad_prev = nstep
       else
          dtime_rad = (nstep - nstep_rad_prev) * dtime
          nstep_rad_prev = nstep
       endif
    end if
  end subroutine update_rad_dtime

  !=========================================================================================

  integer function get_rad_step_size()

    character(len=*), parameter :: sub = 'clm::get_rad_step_size'

    if ( .not. check_timemgr_initialized(sub) ) return

    if (nstep_rad_prev == uninit_int ) then
       get_rad_step_size=get_step_size()
    else
       get_rad_step_size=dtime_rad
    end if

  end function get_rad_step_size

  !=========================================================================================

  integer function get_nstep()

    ! Return the timestep number.

    character(len=*), parameter :: sub = 'clm::get_nstep'
    integer :: rc
    integer(ESMF_KIND_I8) :: step_no

    if ( .not. check_timemgr_initialized(sub) ) return

    call ESMF_ClockGet(tm_clock, advanceCount=step_no, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_ClockGet')

    get_nstep = step_no

  end function get_nstep

  !=========================================================================================

  integer function get_nstep_since_startup_or_lastDA_restart_or_pause()

    ! Return the number of time-steps since the restart file was modified

    character(len=*), parameter :: sub = 'clm::get_nstep_since_rest_mod'

    get_nstep_since_startup_or_lastDA_restart_or_pause = get_nstep() - DA_nstep

  end function get_nstep_since_startup_or_lastDA_restart_or_pause

  !=========================================================================================

  subroutine get_curr_date(yr, mon, day, tod, offset)

    !-----------------------------------------------------------------------------------------
    ! Return date components valid at end of current timestep with an optional
    ! offset (positive or negative) in seconds.

    integer, intent(out) ::&
         yr,    &! year
         mon,   &! month
         day,   &! day of month
         tod     ! time of day (seconds past 0Z)

    integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
    ! Positive for future times, negative 
    ! for previous times.

    character(len=*), parameter :: sub = 'clm::get_curr_date'
    integer :: rc
    type(ESMF_Time) :: date
    type(ESMF_TimeInterval) :: off
    !-----------------------------------------------------------------------------------------

    if ( .not. check_timemgr_initialized(sub) ) return

    call ESMF_ClockGet( tm_clock, currTime=date, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_ClockGet')

    if (present(offset)) then
       if (offset > 0) then
          call ESMF_TimeIntervalSet( off, s=offset, rc=rc )
          call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
          date = date + off
       else if (offset < 0) then
          call ESMF_TimeIntervalSet( off, s=-offset, rc=rc )
          call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
          date = date - off
       end if
    end if

    call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_TimeGet')

  end subroutine get_curr_date

  !=========================================================================================

  subroutine get_perp_date(yr, mon, day, tod, offset)

    !-----------------------------------------------------------------------------------------
    ! Return time of day valid at end of current timestep and the components
    ! of the perpetual date (with an optional offset (positive or negative) in seconds.

    integer, intent(out) ::&
         yr,    &! year
         mon,   &! month
         day,   &! day of month
         tod     ! time of day (seconds past 0Z)

    integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
    ! Positive for future times, negative 
    ! for previous times.

    character(len=*), parameter :: sub = 'clm::get_perp_date'
    integer :: rc
    type(ESMF_Time) :: date
    type(ESMF_TimeInterval) :: DelTime
    !-----------------------------------------------------------------------------------------

    if ( .not. check_timemgr_initialized(sub) ) return

    call ESMF_ClockGet( tm_clock, currTime=date, rc=rc )
    ! Get time of day add it to perpetual date
    ! Get year, month, day so that seconds are time-of-day rather than since start time
    call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_TimeGet')
    call ESMF_TimeIntervalSet(DelTime, s=tod, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
    date = tm_perp_date + DelTime
    if ( present(offset) )then
       call ESMF_TimeIntervalSet(DelTime, s=offset, rc=rc)
       call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
       date = date + DelTime
    end if
    ! Get time of day from the result
    ! Get year, month, day so that seconds are time-of-day rather than since start time
    call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)

    ! Get the date from the fixed perpetual date (in case it overflows to next day)
    call ESMF_TimeGet(tm_perp_date, yy=yr, mm=mon, dd=day, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_TimeGet')

  end subroutine get_perp_date

  !=========================================================================================

  subroutine get_prev_date(yr, mon, day, tod)

    ! Return date components valid at beginning of current timestep.

    ! Arguments
    integer, intent(out) ::&
         yr,    &! year
         mon,   &! month
         day,   &! day of month
         tod     ! time of day (seconds past 0Z)

    ! Local variables
    character(len=*), parameter :: sub = 'clm::get_prev_date'
    integer :: rc
    type(ESMF_Time) :: date
    !-----------------------------------------------------------------------------------------

    if ( .not. check_timemgr_initialized(sub) ) return

    call ESMF_ClockGet(tm_clock, prevTime=date, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_ClockGet')

    call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_TimeGet')

  end subroutine get_prev_date

  !=========================================================================================

  subroutine get_start_date(yr, mon, day, tod)

    ! Return date components valid at beginning of initial run.

    ! Arguments
    integer, intent(out) ::&
         yr,    &! year
         mon,   &! month
         day,   &! day of month
         tod     ! time of day (seconds past 0Z)

    ! Local variables
    character(len=*), parameter :: sub = 'clm::get_start_date'
    integer :: rc
    type(ESMF_Time) :: date
    !-----------------------------------------------------------------------------------------

    if ( .not. check_timemgr_initialized(sub) ) return

    call ESMF_ClockGet(tm_clock, startTime=date, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_ClockGet')

    call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_TimeGet')

  end subroutine get_start_date

  !=========================================================================================

  integer function get_driver_start_ymd( tod )

    ! Return date of start of simulation from driver (i.e. NOT from restart file)
    ! Note: get_start_date gets you the date from the beginning of the simulation
    !       on the restart file.

    ! Arguments
    integer, optional, intent(out) ::&
         tod     ! time of day (seconds past 0Z)

    ! Local variables
    character(len=*), parameter :: sub = 'clm::get_driver_start_ymd'
    !-----------------------------------------------------------------------------------------

    if ( .not. check_timemgr_initialized(sub) ) return

    if ( start_ymd == uninit_int )then
       call shr_sys_abort( sub//': error driver start date is NOT set yet' )
    end if
    if ( start_ymd < 101 .or. start_ymd > 99991231 )then
       call shr_sys_abort( sub//': error driver start date is invalid' )
    end if
    if ( present(tod) )then
       tod = start_tod
       if ( (tod < 0) .or. (tod > isecspday) )then
          call shr_sys_abort( sub//': error driver start tod is invalid' )
       end if
    end if
    get_driver_start_ymd = start_ymd

  end function get_driver_start_ymd

  !=========================================================================================

  subroutine get_ref_date(yr, mon, day, tod)

    ! Return date components of the reference date.

    ! Arguments
    integer, intent(out) ::&
         yr,    &! year
         mon,   &! month
         day,   &! day of month
         tod     ! time of day (seconds past 0Z)

    ! Local variables
    character(len=*), parameter :: sub = 'clm::get_ref_date'
    integer :: rc
    type(ESMF_Time) :: date
    !-----------------------------------------------------------------------------------------

    if ( .not. check_timemgr_initialized(sub) ) return

    call ESMF_ClockGet(tm_clock, refTime=date, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_ClockGet')

    call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_TimeGet')

  end subroutine get_ref_date

  !=========================================================================================

  subroutine get_curr_time(days, seconds)

    ! Return time components valid at end of current timestep.
    ! Current time is the time interval between the current date and the reference date.

    ! Arguments
    integer, intent(out) ::&
         days,   &! number of whole days in time interval
         seconds  ! remaining seconds in time interval

    ! Local variables
    character(len=*), parameter :: sub = 'clm::get_curr_time'
    integer :: rc
    type(ESMF_Time) :: cdate, rdate
    type(ESMF_TimeInterval) :: diff
    !-----------------------------------------------------------------------------------------

    if ( .not. check_timemgr_initialized(sub) ) return

    call ESMF_ClockGet( tm_clock, currTime=cdate, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_ClockGet')

    call ESMF_ClockGet( tm_clock, refTime=rdate, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_ClockGet')

    diff = cdate - rdate

    call ESMF_TimeIntervalGet(diff, d=days, s=seconds, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_TimeIntervalGet')

  end subroutine get_curr_time

  !=========================================================================================

  subroutine get_prev_time(days, seconds)

    ! Return time components valid at beg of current timestep.
    ! prev time is the time interval between the prev date and the reference date.

    ! Arguments
    integer, intent(out) ::&
         days,   &! number of whole days in time interval
         seconds  ! remaining seconds in time interval

    ! Local variables
    character(len=*), parameter :: sub = 'clm::get_prev_time'
    integer :: rc
    type(ESMF_Time) :: date, ref_date
    type(ESMF_TimeInterval) :: diff
    !-----------------------------------------------------------------------------------------

    if ( .not. check_timemgr_initialized(sub) ) return

    call ESMF_ClockGet(tm_clock, prevTime=date, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_ClockGet for prevTime')
    call ESMF_ClockGet(tm_clock, refTime=ref_date, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_ClockGet for refTime')
    diff = date - ref_date
    call ESMF_TimeIntervalGet( diff, d=days, s=seconds, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_TimeintervalGet')

  end subroutine get_prev_time

  !=========================================================================================

  function get_curr_calday(offset, reuse_day_365_for_day_366)

    ! Return calendar day at end of current timestep with optional offset.
    ! Calendar day 1.0 = 0Z on Jan 1.

    ! Arguments
    integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
    ! Positive for future times, negative
    ! for previous times.

    ! If present and true, then day 366 (i.e., the last day of the year on leap years when
    ! using a Gregorian calendar) reuses day 365. Note that this leads to non-monotonic
    ! values throughout the year. This is needed in situations where the calday is used
    ! in code that assumes a 365 day year and won't work right for day 366, such as
    ! shr_orb_decl.
    logical, optional, intent(in) :: reuse_day_365_for_day_366

    ! Return value
    real(r8) :: get_curr_calday

    ! Local variables
    character(len=*), parameter :: sub = 'clm::get_curr_calday'
    integer :: rc
    logical :: l_reuse_day_365_for_day_366  ! local version of reuse_day_365_for_day_366
    type(ESMF_Time) :: date
    type(ESMF_TimeInterval) :: off, diurnal
    integer :: year, month, day, tod
    !-----------------------------------------------------------------------------------------

    if ( .not. check_timemgr_initialized(sub) ) return

    if (present(reuse_day_365_for_day_366)) then
       l_reuse_day_365_for_day_366 = reuse_day_365_for_day_366
    else
       l_reuse_day_365_for_day_366 = .false.
    end if

    call ESMF_ClockGet( tm_clock, currTime=date, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_ClockGet')

    if (present(offset)) then
       if (offset > 0) then
          call ESMF_TimeIntervalSet( off, s=offset, rc=rc )
          call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
          date = date + off
       else if (offset < 0) then
          call ESMF_TimeIntervalSet( off, s=-offset, rc=rc )
          call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
          date = date - off
       end if
    end if

    if ( tm_perp_calendar ) then
       call ESMF_TimeGet(date, yy=year, mm=month, dd=day, s=tod, rc=rc)
       call chkrc(rc, sub//': error return from ESMF_TimeGet')
       call ESMF_TimeIntervalSet( diurnal, s=tod, rc=rc )
       call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
       date = tm_perp_date + diurnal
    end if

    call ESMF_TimeGet( date, dayOfYear_r8=get_curr_calday, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_TimeGet')
    !----------------------------------------------------------------------------------------!
    !!!!!!!!!!!!!! WARNING HACK TO ENABLE Gregorian CALENDAR WITH SHR_ORB !!!!!!!!!!!!!!!!!!!!
    !!!! The following hack fakes day 366 by reusing day 365. This is just because the  !!!!!!
    !!!! current shr_orb_decl calculation can't handle days > 366.                      !!!!!!
    !!!!       Dani Bundy-Coleman and Erik Kluzek Aug/2008                              !!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (l_reuse_day_365_for_day_366) then
       if ( (get_curr_calday > 366.0) .and. (get_curr_calday <= 367.0) .and. &
            (trim(calendar) == GREGORIAN_C) )then
          get_curr_calday = get_curr_calday - 1.0_r8
       end if
    end if
    !!!!!!!!!!!!!! END HACK TO ENABLE Gregorian CALENDAR WITH SHR_ORB !!!!!!!!!!!!!!!!!!!!!!!!
    !----------------------------------------------------------------------------------------!
    if ( (get_curr_calday < 1.0) .or. (get_curr_calday > 367.0) )then
       write(iulog,*) sub, ' = ', get_curr_calday
       if ( present(offset) ) write(iulog,*) 'offset = ', offset
       call shr_sys_abort( sub//': error get_curr_calday out of bounds' )
    end if

  end function get_curr_calday

  !=========================================================================================

  function get_prev_calday(reuse_day_365_for_day_366)

    ! Return calendar day at beginning of current timestep.
    ! Calendar day 1.0 = 0Z on Jan 1.

    ! If present and true, then day 366 (i.e., the last day of the year on leap years when
    ! using a Gregorian calendar) reuses day 365. Note that this leads to non-monotonic
    ! values throughout the year. This is needed in situations where the calday is used
    ! in code that assumes a 365 day year and won't work right for day 366, such as
    ! shr_orb_decl.
    logical, optional, intent(in) :: reuse_day_365_for_day_366

    ! Return value
    real(r8) :: get_prev_calday

    character(len=*), parameter :: sub = 'clm::get_prev_calday'
    !-----------------------------------------------------------------------------------------

    if ( .not. check_timemgr_initialized(sub) ) return

    get_prev_calday = get_curr_calday( &
         offset = -dtime, &
         reuse_day_365_for_day_366 = reuse_day_365_for_day_366)

  end function get_prev_calday

  !=========================================================================================

  function get_calday(ymd, tod, reuse_day_365_for_day_366)

    ! Return calendar day corresponding to specified time instant.
    ! Calendar day 1.0 = 0Z on Jan 1.

    ! If the current run is using a Gregorian calendar, then the year is important, in
    ! that it determines whether or not we're in a leap year for the sake of determining
    ! the calendar day.

    ! Arguments
    integer, intent(in) :: &
         ymd,   &! date in yearmmdd format
         tod     ! time of day (seconds past 0Z)

    ! If present and true, then day 366 (i.e., the last day of the year on leap years when
    ! using a Gregorian calendar) reuses day 365. Note that this leads to non-monotonic
    ! values throughout the year. This is needed in situations where the calday is used
    ! in code that assumes a 365 day year and won't work right for day 366, such as
    ! shr_orb_decl.
    logical, optional, intent(in) :: reuse_day_365_for_day_366

    ! Return value
    real(r8) :: get_calday

    ! Local variables
    character(len=*), parameter :: sub = 'clm::get_calday'
    integer :: rc                 ! return code
    logical :: l_reuse_day_365_for_day_366  ! local version of reuse_day_365_for_day_366
    type(ESMF_Time) :: date
    !-----------------------------------------------------------------------------------------

    if (present(reuse_day_365_for_day_366)) then
       l_reuse_day_365_for_day_366 = reuse_day_365_for_day_366
    else
       l_reuse_day_365_for_day_366 = .false.
    end if

    date = TimeSetymd( ymd, tod, "get_calday" )
    call ESMF_TimeGet( date, dayOfYear_r8=get_calday, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_TimeGet')
    !----------------------------------------------------------------------------------------!
!!!!!!!!!!!!!! WARNING HACK TO ENABLE Gregorian CALENDAR WITH SHR_ORB !!!!!!!!!!!!!!!!!!!!
!!!! The following hack fakes day 366 by reusing day 365. This is just because the  !!!!!!
!!!! current shr_orb_decl calculation can't handle days > 366.                      !!!!!!
!!!!       Dani Bundy-Coleman and Erik Kluzek Aug/2008                              !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (l_reuse_day_365_for_day_366) then
       if ( (get_calday > 366.0) .and. (get_calday <= 367.0) .and. &
            (trim(calendar) == GREGORIAN_C) )then
          get_calday = get_calday - 1.0_r8
       end if
    end if
!!!!!!!!!!!!!! END HACK TO ENABLE Gregorian CALENDAR WITH SHR_ORB !!!!!!!!!!!!!!!!!!!!!!!!
    !----------------------------------------------------------------------------------------!
    if ( (get_calday < 1.0) .or. (get_calday > 367.0) )then
       write(iulog,*) sub, ' = ', get_calday
       call shr_sys_abort( sub//': error calday out of range' )
    end if

  end function get_calday

  !=========================================================================================

  function get_calendar()

    ! Return calendar

    ! Return value
    character(len=ESMF_MAXSTR) :: get_calendar

    get_calendar = calendar

  end function get_calendar

  !=========================================================================================

  function get_doy_tomorrow(doy_today) result(doy_tomorrow)

    !---------------------------------------------------------------------------------
    ! Given a day of the year (doy_today), return the next day of the year

    integer, intent(in) :: doy_today
    integer             :: doy_tomorrow
    integer             :: days_in_year
    character(len=*), parameter :: sub = 'clm::get_doy_tomorrow'

    ! Use get_prev_days_per_year() instead of get_curr_days_per_year() because the latter, in the last timestep of a year, actually returns the number of days in the NEXT year.
    days_in_year = get_prev_days_per_year()

    if ( doy_today < 1 .or. doy_today > days_in_year )then
       write(iulog,*) 'doy_today    = ', doy_today
       write(iulog,*) 'days_in_year = ', days_in_year
       call shr_sys_abort( sub//': error doy_today out of range' )
    end if

    if (doy_today == days_in_year) then
        doy_tomorrow = 1
    else
        doy_tomorrow = doy_today + 1
    end if
  end function get_doy_tomorrow

  !=========================================================================================

  real(r8) function get_average_days_per_year()

    !---------------------------------------------------------------------------------
    ! Get the average number of days per year for the given calendar.
    !
    ! This should be used, for example, when converting a parameter from units of
    ! per-year to units of per-second (so that the parameter will have a fixed, constant
    ! value rather than a slightly different value on leap years vs. non-leap years).

    real(r8) :: avg_days_per_year
    real(r8) :: curr_days_per_year

    real(r8), parameter :: days_per_year_noleap = 365._r8

    ! From the definition of ESMF_CALKIND_GREGORIAN in
    ! https://earthsystemmodeling.org/docs/release/latest/ESMF_refdoc/node6.html: "In the
    ! Gregorian calendar every fourth year is a leap year in which February has 29 and not
    ! 28 days; however, years divisible by 100 are not leap years unless they are also
    ! divisible by 400." This results in an average number of days per year of 365.2425.
    real(r8), parameter :: days_per_year_gregorian = 365.2425_r8

    character(len=*), parameter :: subname = 'get_average_days_per_year'
    !---------------------------------------------------------------------------------

    ! BUG(wjs, 2022-02-01, ESCOMP/CTSM#1624) Ideally we would use ESMF_CalendarGet here,
    ! but that currently isn't possible (see notes in issue 1624 for details)
    if (to_upper(calendar) == NO_LEAP_C) then
       avg_days_per_year = days_per_year_noleap
    else if (to_upper(calendar) == GREGORIAN_C) then
       avg_days_per_year = days_per_year_gregorian
    else
       call shr_sys_abort(subname//' ERROR: unrecognized calendar specified= '//trim(calendar))
    end if

    ! Paranoia: Since we're using a hard-coded value, let's make sure that the user hasn't
    ! done some customizations to the calendar that change the days per year from what we
    ! expect: Compare the hard-coded value with the number of days per year in the
    ! current year, which comes from the actual ESMF calendar; the two should be close.
    ! (This check can be removed once we address issue 1624, making the results of this
    ! function depend on the actual ESMF calendar instead of a hard-coded value.)
    curr_days_per_year = get_curr_days_per_year()
    if (abs(avg_days_per_year - curr_days_per_year) > 1._r8) then
       write(iulog,*) 'ERROR: hard-coded average days per year differs by more than expected'
       write(iulog,*) 'from current days per year. Are you using a non-standard calendar?'
       write(iulog,*) 'avg_days_per_year (hard-coded)          = ', avg_days_per_year
       write(iulog,*) 'curr_days_per_year (from ESMF calendar) = ', curr_days_per_year
       write(iulog,*) 'You can fix this by changing the hard-coded parameters in '//subname
       write(iulog,*) 'in file: '//sourcefile
       call shr_sys_abort(subname//' ERROR: hard-coded average days per year differs by more than expected')
    end if

    get_average_days_per_year = avg_days_per_year

  end function get_average_days_per_year

  !=========================================================================================

  integer function get_curr_days_per_year( offset )

    !---------------------------------------------------------------------------------
    ! Get the number of days per year for the year as of the end of the current time step
    ! (or offset from the end of the current time step, if offset is provided).
    !
    ! For the last time step of the year, note that the this will give the number of days
    ! per year in the about-to-start year, not in the just-finishing year.

    !
    ! Arguments
    integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
    ! Positive for future times, negative 
    ! for previous times.

    character(len=*), parameter :: sub = 'clm::get_curr_days_per_year'
    integer         :: yr, mon, day, tod ! current date year, month, day and time-of-day
    type(ESMF_Time) :: eDate             ! ESMF date
    integer         :: rc                ! ESMF return code
    !---------------------------------------------------------------------------------

    if ( .not. check_timemgr_initialized(sub) ) return

    if ( present(offset) )then
       call get_curr_date(yr, mon, day, tod, offset )
    else
       call get_curr_date(yr, mon, day, tod )
    end if
    eDate = TimeSetymd( ymd=yr*10000+1231, tod=0, desc="end of year" )
    call ESMF_TimeGet( eDate, dayOfYear=get_curr_days_per_year, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_TimeGet')

  end function get_curr_days_per_year

  !=========================================================================================

  integer function get_prev_days_per_year()

    !---------------------------------------------------------------------------------
    ! Get the number of days per year for the year as of the beginning of the current time step
    !
    ! For the last time step of the year, note that the this will give the number of days
    ! per year in the just-finishing year.

    character(len=*), parameter :: sub = 'clm::get_prev_days_per_year'

    if ( .not. check_timemgr_initialized(sub) ) return

    get_prev_days_per_year = get_curr_days_per_year(offset = -dtime)

  end function get_prev_days_per_year

  !=========================================================================================

  function get_curr_yearfrac( offset )

    !---------------------------------------------------------------------------------
    ! Get the fractional position in the current year, as of the end of the current
    ! timestep. This is 0 at midnight on Jan 1, and 1 at the end of Dec 31.

    !
    ! Arguments
    real(r8) :: get_curr_yearfrac  ! function result
    
    integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
    ! Positive for future times, negative 
    ! for previous times.

    character(len=*), parameter :: sub = 'clm::get_curr_yearfrac'
    real(r8) :: cday               ! current calendar day (1.0 = 0Z on Jan 1)
    real(r8) :: days_per_year      ! days per year

    if ( .not. check_timemgr_initialized(sub) ) return

    cday          = get_curr_calday(offset=offset)
    days_per_year = get_curr_days_per_year(offset=offset)

    get_curr_yearfrac = (cday - 1._r8)/days_per_year

  end function get_curr_yearfrac

  !=========================================================================================

  function get_prev_yearfrac()

    !---------------------------------------------------------------------------------
    ! Get the fractional position in the current year, as of the beginning of the current
    ! timestep. This is 0 at midnight on Jan 1, and 1 at the end of Dec 31.

    !
    ! Arguments
    real(r8) :: get_prev_yearfrac  ! function result
    
    character(len=*), parameter :: sub = 'clm::get_prev_yearfrac'
    
    if ( .not. check_timemgr_initialized(sub) ) return
    
    get_prev_yearfrac = get_curr_yearfrac(offset = -dtime)

  end function get_prev_yearfrac


  !=========================================================================================

  subroutine get_rest_date(ncid, yr)

    !---------------------------------------------------------------------------------
    ! Get the date from the restart file.
    !
    ! Currently just returns the year (because the month & day are harder to extract, and
    ! currently aren't needed).
    use ncdio_pio, only: ncd_io, file_desc_t
    !
    ! Arguments
    type(file_desc_t) , intent(inout) :: ncid ! netcdf id for the restart file
    integer           , intent(out)   :: yr   ! year from restart file

    integer :: ymd     ! yyyymmdd from the restart file
    logical :: readvar ! whether the variable was read from the file
    
    integer, parameter :: year_mask = 10000  ! divide by this to get year from ymd

    character(len=*), parameter :: subname = 'get_rest_date'
    !-----------------------------------------------------------------------
    
    ! Get the date (yyyymmdd) from restart file.
    ! Note that we cannot simply use the rst_curr_ymd module variable, because that isn't
    ! set under some circumstances
    call ncd_io(varname='timemgr_rst_curr_ymd', data=ymd, &
         ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) then
       call shr_sys_abort(subname//' ERROR: timemgr_rst_curr_ymd not found on restart file')
    end if
    
    ! Extract the year
    yr = ymd / year_mask
  end subroutine get_rest_date

  !=========================================================================================

  integer function get_local_timestep_time( londeg, offset )

    !---------------------------------------------------------------------------------
    ! Get the local time for this longitude that is evenly divisible by the time-step
    !
    ! uses
    use clm_varcon, only: degpsec, isecspday
    ! Arguments
    real(r8)         , intent(in) :: londeg  ! Longitude in degrees
    integer, optional, intent(in) :: offset  ! Offset from current time in seconds (either sign)

    ! Local variables
    integer  :: yr, mon, day    ! year, month, day, unused
    integer  :: secs            ! seconds into the day
    real(r8) :: lon             ! positive longitude
    integer  :: offset_sec      ! offset seconds (either 0 for current time or -dtime for previous time)
    !---------------------------------------------------------------------------------
    if ( present(offset) ) then
       offset_sec = offset
    else
       offset_sec = 0
    end if
    SHR_ASSERT( londeg >= -180.0_r8, "londeg must be greater than -180" )
    SHR_ASSERT( londeg <= 360.0_r8,  "londeg must be less than 360" )
    call  get_curr_date(yr, mon, day, secs, offset=offset_sec )
    lon = londeg
    if ( lon < 0.0_r8 ) lon = lon + 360.0_r8
    get_local_timestep_time  = secs + nint((lon/degpsec)/real(dtime,r8))*dtime
    get_local_timestep_time  = mod(get_local_timestep_time,isecspday)
  end function get_local_timestep_time

  !=========================================================================================

  integer function get_local_time( londeg, starttime, offset )

    !---------------------------------------------------------------------------------
    ! Get the local time for this longitude
    !
    ! uses
    use clm_varcon, only: degpsec, isecspday
    ! Arguments
    real(r8)         , intent(in) :: londeg       ! Longitude in degrees
    integer, optional, intent(in) :: starttime    ! Start time (sec)
    integer, optional, intent(in) :: offset       ! Offset from current time in seconds (either sign)

    ! Local variables
    integer  :: yr, mon, day    ! year, month, day, unused
    integer  :: secs            ! seconds into the day
    integer  :: start           ! start seconds
    integer  :: offset_sec      ! offset seconds (either 0 for current time or -dtime for previous time)
    real(r8) :: lon             ! positive longitude
    !---------------------------------------------------------------------------------
    if ( present(starttime) ) then
       start = starttime
    else
       start = 0
    end if
    if ( present(offset) ) then
       offset_sec = offset
    else
       offset_sec = 0
    end if
    SHR_ASSERT( start >= 0,            "starttime must be greater than or equal to zero" )
    SHR_ASSERT( start <= isecspday,    "starttime must be less than or equal to number of seconds in a day" )
    SHR_ASSERT( londeg >= -180.0_r8,   "londeg must be greater than -180" )
    SHR_ASSERT( londeg <= 360.0_r8,    "londeg must be less than 360" )
    SHR_ASSERT( (offset_sec == 0) .or. (offset_sec == -dtime), "offset must be zero or negative time-step" )
    call  get_curr_date(yr, mon, day, secs, offset=offset_sec )
    lon = londeg
    if ( lon < 0.0_r8 ) lon = lon + 360.0_r8
    get_local_time  = modulo(secs + nint(londeg/degpsec), isecspday)
    get_local_time  = modulo(get_local_time - start,isecspday)
  end function get_local_time

  !=========================================================================================

  logical function is_near_local_noon( londeg, deltasec )

    !---------------------------------------------------------------------------------
    ! Is this longitude near it's local noon?
    !
    ! uses
    use clm_varcon, only: degpsec, isecspday
    ! Arguments
    real(r8), intent(in) :: londeg   ! Longitude in degrees
    integer , intent(in) :: deltasec ! Number of seconds before or after local noon
    
    ! Local variables
    integer :: local_secs                         ! Local time in seconds
    integer, parameter :: noonsec = isecspday / 2 ! seconds at local noon
    !---------------------------------------------------------------------------------
    SHR_ASSERT( deltasec < noonsec, "deltasec must be less than 12 hours" )
    local_secs = get_local_timestep_time( londeg )

    if ( local_secs >= (noonsec - deltasec) .and. local_secs <= (noonsec + deltasec)) then
       is_near_local_noon = .true.
    else
       is_near_local_noon = .false.
    end if

    !---------------------------------------------------------------------------------
  end function is_near_local_noon

  !=========================================================================================
 
  function is_beg_curr_day()
 
     ! Return true if current timestep is first timestep in current day.
     
     ! Return value
     logical :: is_beg_curr_day
  
     ! Local variables
     integer ::&
        yr,    &! year
        mon,   &! month
        day,   &! day of month
        tod     ! time of day (seconds past 0Z)
 
    character(len=*), parameter :: sub = 'clm::is_beg_curr_day'

     if ( .not. check_timemgr_initialized(sub) ) return

     call get_curr_date(yr, mon, day, tod)
     is_beg_curr_day = ( tod == dtime )
 
  end function is_beg_curr_day

  !=========================================================================================

  function is_end_curr_day()

    !---------------------------------------------------------------------------------
    ! Return true if current timestep is last timestep in current day.

    ! Return value
    logical :: is_end_curr_day

    ! Local variables
    integer ::&
         yr,    &! year
         mon,   &! month
         day,   &! day of month
         tod     ! time of day (seconds past 0Z)

    character(len=*), parameter :: sub = 'clm::is_end_curr_day'
    !---------------------------------------------------------------------------------

    if ( .not. check_timemgr_initialized(sub) ) return

    call get_curr_date(yr, mon, day, tod)
    is_end_curr_day = (tod == 0)

  end function is_end_curr_day

  !=========================================================================================

  logical function is_end_curr_month()

    !---------------------------------------------------------------------------------
    ! Return true if current timestep is last timestep in current month.

    ! Local variables
    integer ::&
         yr,    &! year
         mon,   &! month
         day,   &! day of month
         tod     ! time of day (seconds past 0Z)

    character(len=*), parameter :: sub = 'clm::is_end_curr_month'
    !---------------------------------------------------------------------------------

    if ( .not. check_timemgr_initialized(sub) ) return

    call get_curr_date(yr, mon, day, tod)
    is_end_curr_month = (day == 1  .and.  tod == 0)

  end function is_end_curr_month

  !-----------------------------------------------------------------------
  logical function is_beg_curr_year()
    !
    ! !DESCRIPTION:
    ! Return true if current timestep is first timestep in current year.
    !
    ! !LOCAL VARIABLES:
    integer ::&
         yr,    &! year
         mon,   &! month
         day,   &! day of month
         tod     ! time of day (seconds past 0Z)

    character(len=*), parameter :: subname = 'is_beg_curr_year'
    !-----------------------------------------------------------------------

    if ( .not. check_timemgr_initialized(subname) ) return

    call get_curr_date(yr, mon, day, tod)
    is_beg_curr_year = (mon == 1 .and. day == 1 .and. tod == dtime)
    
  end function is_beg_curr_year

  !-----------------------------------------------------------------------
  logical function is_end_curr_year()
    !
    ! !DESCRIPTION:
    ! Return true if current timestep is last timestep in current year.
    !
    ! !LOCAL VARIABLES:
    integer ::&
         yr,    &! year
         mon,   &! month
         day,   &! day of month
         tod     ! time of day (seconds past 0Z)

    character(len=*), parameter :: subname = 'is_end_curr_year'
    !-----------------------------------------------------------------------

    if ( .not. check_timemgr_initialized(subname) ) return

    call get_curr_date(yr, mon, day, tod)
    is_end_curr_year = (mon == 1 .and. day == 1 .and. tod == 0)
    
  end function is_end_curr_year


  !=========================================================================================

  logical function is_first_step()

    !---------------------------------------------------------------------------------
    ! Return true on first step of startup and hybrid runs.

    ! Local variables
    character(len=*), parameter :: sub = 'clm::is_first_step'
    integer :: rc
    integer :: nstep
    integer(ESMF_KIND_I8) :: step_no
    !---------------------------------------------------------------------------------

    if ( .not. check_timemgr_initialized(sub) ) return

    call ESMF_ClockGet( tm_clock, advanceCount=step_no, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_ClockGet')
    nstep = step_no
    is_first_step = (nstep == 1)

  end function is_first_step
  !=========================================================================================

  logical function is_first_restart_step()

    ! Return true on first step of restart or branch run only.
    character(len=*), parameter :: sub = 'clm::is_first_restart_step'

    if ( .not. check_timemgr_initialized(sub) ) return

    is_first_restart_step = tm_first_restart_step

  end function is_first_restart_step

  !=========================================================================================
  
  logical function is_first_step_of_this_run_segment()

    ! Return true if this is the first step of this run segment. This will be true for
    ! the first step of a startup, restart or branch run.
    character(len=*), parameter :: sub = 'clm::is_first_step_of_this_run_segment'

    if ( .not. check_timemgr_initialized(sub) ) return

    is_first_step_of_this_run_segment = (is_first_step() .or. is_first_restart_step())

  end function is_first_step_of_this_run_segment

  !=========================================================================================

  logical function is_perpetual()

    ! Return true on last timestep.
    character(len=*), parameter :: sub = 'clm::is_perpetual'

    if ( .not. check_timemgr_initialized(sub) ) return

    is_perpetual = tm_perp_calendar

  end function is_perpetual

  !=========================================================================================

  logical function is_doy_in_interval(start, end, doy)

    ! Return true if day of year is in the provided interval.
    ! Does not treat leap years differently from normal years.
    ! Arguments
    integer, intent(in) :: start ! start of interval (day of year)
    integer, intent(in) :: end ! end of interval (day of year)
    integer, intent(in) :: doy ! day of year to query
    
    ! Local variables
    logical :: window_crosses_newyear

    character(len=*), parameter :: sub = 'clm::is_doy_in_interval'

    window_crosses_newyear = end < start

    if (window_crosses_newyear .and. &
        (doy >= start .or. doy <= end)) then
       is_doy_in_interval = .true.
    else if (.not. window_crosses_newyear .and. &
        (doy >= start .and. doy <= end)) then
       is_doy_in_interval = .true.
    else
       is_doy_in_interval = .false.
    end if
    
  end function is_doy_in_interval

  !=========================================================================================

  logical function is_today_in_doy_interval(start, end)

    ! Return true if today's day of year is in the provided interval.
    ! Does not treat leap years differently from normal years.
    ! Arguments
    integer, intent(in) :: start ! start of interval (day of year)
    integer, intent(in) :: end ! end of interval (day of year)

    ! Local variable(s)
    integer :: doy_today

    character(len=*), parameter :: sub = 'clm::is_today_in_doy_interval'

    ! Get doy of beginning of current timestep
    doy_today = get_prev_calday()

    is_today_in_doy_interval = is_doy_in_interval(start, end, doy_today)

  end function is_today_in_doy_interval

  !=========================================================================================

  subroutine timemgr_datediff(ymd1, tod1, ymd2, tod2, days)

    ! Calculate the difference (ymd2,tod2) - (ymd1,tod1) and return the result in days.
    ! Arguments
    integer, intent(in) ::&
         ymd1,    &! date1 in yyyymmdd format
         tod1,    &! time of day relative to date1 (seconds past 0Z)
         ymd2,    &! date2 in yyyymmdd format
         tod2      ! time of day relative to date2 (seconds past 0Z)

    real(r8) :: days ! (ymd2,tod2)-(ymd1,tod1) in days

    ! Local variables
    character(len=*), parameter :: sub = 'clm::timemgr_datediff'
    integer :: rc   ! return code

    type(ESMF_Time) :: date1
    type(ESMF_Time) :: date2
    type(ESMF_TimeInterval) :: diff
    !-----------------------------------------------------------------------------------------

    date1 = TimeSetymd( ymd1, tod1, "date1" )
    date2 = TimeSetymd( ymd2, tod2, "date2" )
    diff = date2 - date1
    call ESMF_TimeIntervalGet( diff, d_r8=days, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_TimeIntervalGet')
    days = days + 1.0_r8

  end subroutine timemgr_datediff

  !=========================================================================================

  subroutine chkrc(rc, mes)
    use ESMF        , only : ESMF_SUCCESS
    integer, intent(in)          :: rc   ! return code from time management library
    character(len=*), intent(in) :: mes  ! error message
    if ( rc == ESMF_SUCCESS ) return
    write(iulog,*) mes
    call shr_sys_abort ('CHKRC')
  end subroutine chkrc

  !=========================================================================================

  function to_upper(str)

    !---------------------------------------------------------------------------------
    ! Convert character string to upper case. Use achar and iachar intrinsics
    ! to ensure use of ascii collating sequence.
    !
    ! !INPUT PARAMETERS:
    character(len=*), intent(in) :: str ! String to convert to upper case
    ! !RETURN VALUE:
    character(len=len(str))      :: to_upper
    ! !LOCAL VARIABLES:
    integer :: i                ! Index
    integer :: aseq             ! ascii collating sequence
    character(len=1) :: ctmp    ! Character temporary
    !---------------------------------------------------------------------------------

    do i = 1, len(str)
       ctmp = str(i:i)
       aseq = iachar(ctmp)
       if ( aseq >= 97  .and.  aseq <= 122 ) ctmp = achar(aseq - 32)
       to_upper(i:i) = ctmp
    end do

  end function to_upper

  !=========================================================================================

  logical function is_restart( )
    ! Determine if restart run
    use clm_varctl, only : nsrest, nsrContinue
    if (nsrest == nsrContinue) then
       is_restart = .true.
    else
       is_restart = .false.
    end if
  end function is_restart

  !=========================================================================================

  subroutine timemgr_spmdbcast( )

    use spmdMod    , only : mpicom, MPI_INTEGER
    use shr_mpi_mod, only : shr_mpi_bcast

    integer :: ier

    call shr_mpi_bcast (dtime, mpicom)

  end subroutine timemgr_spmdbcast

  !=========================================================================================
  
  logical function check_timemgr_initialized(caller)
    !
    ! !DESCRIPTION:
    ! Checks if the time manager has been initialized. If not, aborts with an error
    ! message.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: caller ! name of calling routine
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'check_timemgr_initialized'
    !-----------------------------------------------------------------------
    
    if (.not. timemgr_set) then
       call shr_sys_abort(trim(caller)//":: Time manager has not been initialized")
       check_timemgr_initialized = .false.
    else
       check_timemgr_initialized = .true.
    end if
    
  end function check_timemgr_initialized

  !-----------------------------------------------------------------------
  subroutine timemgr_reset()
    !
    ! !DESCRIPTION:
    ! Reset time manager module data to default values.
    !
    ! All unit tests that modify the time manager should call this routine in their
    ! teardown section.
    !
    ! It is safe to call this subroutine even if the time manager hasn't been initialized.
    ! (In this case, this reset routine won't do anything.)
    !
    ! Note: we could probably get away with doing much less resetting than is currently
    ! done here. For example, we could simply set timemgr_set = .false., and deallocate
    ! anything that needs deallocation. That would provide the benefit of less
    ! maintenance, at the cost of slightly less robustness (in case some variable isn't
    ! set in the initialization of a unit test, either because the unit test forgets to
    ! call the time manager initialization method, or because the initialization method
    ! does not explicitly initialize all variables).
    !
    ! !USES:
    use ESMF      , only : ESMF_ClockDestroy
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: rc ! return code

    character(len=*), parameter :: sub = 'timemgr_reset'
    !-----------------------------------------------------------------------
    
    ! ------------------------------------------------------------------------
    ! The values in the following section should match the initialization values given in
    ! the variable declarations at the top of the module. 
    !
    ! Note: it would be easier to ensure this match if we introduced a time manager
    ! derived type, which had default initialization of its components. Then this routine
    ! could simply set to time manager instance to a new instance of the derived type.
    ! ------------------------------------------------------------------------

    if (.not. timemgr_set) then
       ! If the time manager hasn't been initialized, then we don't need to do anything.
       ! This logic makes it safe to call this reset routine even in cases where the time
       ! manager hasn't been initialized.
       return
    end if

    calendar = NO_LEAP_C

    dtime          = uninit_int
    dtime_rad      = uninit_int
    nstep_rad_prev = uninit_int
    
    start_ymd = uninit_int
    start_tod = 0
    ref_ymd   = uninit_int
    ref_tod   = 0

    rst_step_sec  = uninit_int
    rst_start_ymd = uninit_int
    rst_start_tod = uninit_int
    rst_ref_ymd   = uninit_int
    rst_ref_tod   = uninit_int
    rst_curr_ymd  = uninit_int
    rst_curr_tod  = uninit_int

    ! note that rst_nstep_rad_prev is NOT initialized in its declaration
    rst_nstep_rad_prev    = uninit_int
    perpetual_ymd         = uninit_int
    tm_first_restart_step = .false.
    tm_perp_calendar      = .false.
    timemgr_set           = .false.

    ! ------------------------------------------------------------------------
    ! Reset other module-level variables to some reasonable default, to ensure that they
    ! don't carry over any state from one unit test to the next.
    ! ------------------------------------------------------------------------

    ! Reset tm_cal
    call init_calendar()

    ! Reset portions of the clock. Note that this does not fully reset the clock, and so
    ! there is still the potential for information in the clock to carry over to the next
    ! unit test if the next test does not properly initialize things.
    call ESMF_ClockDestroy(tm_clock, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_ClockDestroy')

    ! Note that we do NOT currently reset tm_perp_date, because it's unclear what that
    ! should be reset to. Thus, there is potential for its information to carry over to
    ! the next unit test if the next test does not properly initialize things.

  end subroutine timemgr_reset

  ! ========================================================================
  ! The following routines are meant to be used just in unit tests
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine for_test_set_curr_date(yr, mon, day, tod)
    !
    ! !DESCRIPTION:
    ! Sets the current date - i.e., the date at the end of the time step
    !
    ! This is done in a way that mimics what would happen if the clock advanced from the
    ! previous time step to the specified time (so, in particular, sets the previous time
    ! as if that's what happened).
    !
    ! Note that, because of the method used to do this setting, it is an error to try to
    ! call this with the earliest possible time (yr,mon,day,tod = 1,1,1,0): instead, it
    ! needs to be called with a time at least one time step later.
    !
    ! Also note that an unavoidable side-effect of this method is that the time step count
    ! is incremented by 1.
    !
    ! *** Should only be used in unit tests!!! ***
    !
    ! !USES:
    use ESMF    , only : ESMF_ClockSet, ESMF_ClockAdvance
    !
    ! !ARGUMENTS:
    integer, intent(in) :: yr  ! year
    integer, intent(in) :: mon ! month
    integer, intent(in) :: day ! day of month
    integer, intent(in) :: tod ! time of day (seconds past 0Z)
    !
    ! !LOCAL VARIABLES:
    type(ESMF_Time) :: input_time ! ESMF Time corresponding to the inputs
    type(ESMF_TimeInterval) :: interval_dtime
    type(ESMF_Time) :: input_time_minus_dtime
    integer :: rc ! return code

    character(len=*), parameter :: sub = 'for_test_set_curr_date'
    !-----------------------------------------------------------------------
    
    ! Rather than simply setting the clock to the specified date, we instead set it to one
    ! time step before the specified date, then advance the clock by a time step. This is
    ! needed so that the clock's previous time is set as if we reached the specified time
    ! through a typical one-time-step advance of the clock, rather than via an arbitrary
    ! jump.

    ! Because of this method of setting the clock, though, it is an error to call this
    ! with yr,mon,day,tod = 1,1,1,0; catch that common error here and give a meaningful
    ! error message.
    if (yr == 1 .and. mon == 1 .and. day == 1 .and. tod == 0) then
       call shr_sys_abort(sub//': need to use a time later than yr,mon,day,tod = 1,1,1,0')
    end if

    call ESMF_TimeSet(input_time, yy=yr, mm=mon, dd=day, s=tod, &
         calendar=tm_cal, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_TimeSet')

    call ESMF_TimeIntervalSet(interval_dtime, s=dtime, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')

    input_time_minus_dtime = input_time - interval_dtime
    
    call ESMF_ClockSet(tm_clock, CurrTime=input_time_minus_dtime, rc=rc)
    call chkrc(rc, sub//': error return from ESMF_ClockSet')

    ! Note that this ClockAdvance call increments the time step count; this seems like an
    ! unavoidable side effect of this technique of starting with an earlier time and
    ! advancing the clock (which we do in order to get the clock's previous time set
    ! correctly).
    call ESMF_ClockAdvance( tm_clock, rc=rc )
    call chkrc(rc, sub//': error return from ESMF_ClockAdvance')

  end subroutine for_test_set_curr_date


end module clm_time_manager
