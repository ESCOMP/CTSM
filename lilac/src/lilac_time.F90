module lilac_time

  use ESMF
  use shr_kind_mod  , only : cx=>shr_kind_cx, cs=>shr_kind_cs, cl=>shr_kind_cl, r8=>shr_kind_r8
  use shr_sys_mod   , only : shr_sys_abort
  use shr_cal_mod   , only : shr_cal_ymd2date
  use lilac_io      , only : lilac_io_write, lilac_io_wopen, lilac_io_enddef
  use lilac_io      , only : lilac_io_close, lilac_io_date2yyyymmdd, lilac_io_sec2hms
  use lilac_methods , only : chkerr
  use lilac_constants, only : logunit
  use netcdf        , only : nf90_open, nf90_nowrite, nf90_noerr
  use netcdf        , only : nf90_inq_varid, nf90_get_var, nf90_close

  implicit none
  private    ! default private

  public :: lilac_time_clockinit     ! initialize the lilac clock
  public :: lilac_time_alarminit     ! initialize an alarm
  public :: lilac_time_restart_write ! only writes the time info
  public :: lilac_time_restart_read  ! only reads the time info

  ! Clock and alarm options
  character(len=*), private, parameter :: &
       optNone           = "none"      , &
       optNever          = "never"     , &
       optNSteps         = "nsteps"    , &
       optNSeconds       = "nseconds"  , &
       optNMinutes       = "nminutes"  , &
       optNHours         = "nhours"    , &
       optNDays          = "ndays"     , &
       optNMonths        = "nmonths"   , &
       optNYears         = "nyears"    

  ! Module data
  character(len=ESMF_MAXSTR)  :: caseid
  type(ESMF_Calendar)         :: lilac_calendar
  integer                     :: mytask
  integer, parameter          :: SecPerDay = 86400 ! Seconds per day

  ! We'll use this really big year to effectively mean infinitely into the future.
  integer,  parameter :: really_big_year = 999999999

  character(len=*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine lilac_time_clockInit(caseid_in, starttype, atm_calendar, atm_timestep, &
       atm_start_year, atm_start_mon, atm_start_day, atm_start_secs, &
       lilac_clock, rc)

    ! -------------------------------------------------
    ! Initialize the lilac clock
    ! -------------------------------------------------

    ! input/output variables
    character(len=*)    , intent(in)    :: caseid_in
    character(len=*)    , intent(in)    :: starttype
    character(len=*)    , intent(in)    :: atm_calendar
    integer             , intent(in)    :: atm_timestep
    integer             , intent(in)    :: atm_start_year !(yyyy)
    integer             , intent(in)    :: atm_start_mon  !(mm)
    integer             , intent(in)    :: atm_start_day
    integer             , intent(in)    :: atm_start_secs
    type(ESMF_Clock)    , intent(inout) :: lilac_clock
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_Alarm)        :: lilac_restart_alarm
    type(ESMF_Alarm)        :: lilac_stop_alarm
    type(ESMF_Clock)        :: clock
    type(ESMF_VM)           :: vm
    type(ESMF_Time)         :: StartTime           ! Start time
    type(ESMF_Time)         :: CurrTime            ! Current time
    type(ESMF_Time)         :: Clocktime           ! Loop time
    type(ESMF_TimeInterval) :: TimeStep            ! Clock time-step
    type(ESMF_TimeInterval) :: TimeStep_advance
    integer                 :: start_ymd           ! Start date (YYYYMMDD)
    integer                 :: start_tod           ! Start time of day (seconds)
    integer                 :: curr_ymd            ! Current ymd (YYYYMMDD)
    integer                 :: curr_tod            ! Current tod (seconds)
    character(len=CL)       :: restart_file
    character(len=CL)       :: restart_pfile
    integer                 :: yr, mon, day, secs  ! Year, month, day, seconds as integers
    integer                 :: unitn               ! unit number
    integer                 :: ierr                ! Return code
    integer                 :: tmp(2)              ! Array for Broadcast
    character(len=*), parameter :: subname = '(lilactime_clockInit): '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    caseid = trim(caseid_in)

    ! ------------------------------
    ! get my task
    ! ------------------------------

    call ESMF_VMGetCurrent(vm=vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------------
    ! create lilac_calendar
    ! ------------------------------

    if (trim(atm_calendar) == 'NOLEAP') then
       lilac_calendar = ESMF_CalendarCreate(name='NOLEAP', calkindflag=ESMF_CALKIND_NOLEAP, rc=rc )
    else if (trim(atm_calendar) == 'GREGORIAN') then
       lilac_calendar = ESMF_CalendarCreate(name='GREGORIAN', calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc )
    else
       call shr_sys_abort(trim(subname)//'ERROR: only NOLEAP and GREGORIAN calendars currently supported')
    end if
    call ESMF_CalendarSetDefault(lilac_calendar, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort(trim(subname)//'ERROR: default calendar set error')

    ! ------------------------------
    ! create and initialize lilac_clock
    ! ------------------------------

    call ESMF_TimeIntervalSet(TimeStep, s=atm_timestep, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeSet(StartTime, yy=atm_start_year, mm=atm_start_mon, dd=atm_start_day , s=atm_start_secs, &
         calendar=lilac_calendar, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create the lilac clock
    ! NOTE: the reference time is set to the start time
    ! NOTE: no stop time is given. Stopping will be determined via an argument passed to lilac_run.
    lilac_clock = ESMF_ClockCreate(name='lilac_clock', TimeStep=TimeStep, startTime=StartTime, RefTime=StartTime, &
         rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort(trim(subname)//'error initializing lilac clock')

    ! ------------------------------
    ! For a continue run - obtain current time from the lilac restart file and 
    ! advance the clock to the current time for a continue run
    ! ------------------------------

    if (starttype == 'continue') then

       ! Read the pointer file to obtain the restart file, read the restart file for curr_ymd and curr_tod
       ! and then convert this to an esmf current time (currtime)
       if (mytask == 0) then
          restart_pfile = 'rpointer.lilac'
          call ESMF_LogWrite(trim(subname)//" reading rpointer file = "//trim(restart_pfile), ESMF_LOGMSG_INFO)
          open(newunit=unitn, file=restart_pfile, form='FORMATTED', status='old',iostat=ierr)
          if (ierr < 0) call shr_sys_abort(trim(subname)//' ERROR rpointer file open returns error')
          read(unitn,'(a)', iostat=ierr) restart_file
          if (ierr < 0) call shr_sys_abort(trim(subname)//' ERROR rpointer file read returns error') 
          close(unitn)
          call ESMF_LogWrite(trim(subname)//" read driver restart from "//trim(restart_file), ESMF_LOGMSG_INFO)

          ! Read the restart file on mastertask and then broadcast the data
          call lilac_time_restart_read(restart_file, start_ymd, start_tod, curr_ymd, curr_tod, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          tmp(1) = curr_ymd  ; tmp(2) = curr_tod
       endif
       call ESMF_VMBroadcast(vm, tmp, 4, 0, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       curr_ymd  = tmp(1) ; curr_tod  = tmp(2)

       ! Determine current time
       yr  = int(curr_ymd/10000)
       mon = int( mod(curr_ymd,10000)/ 100)
       day = mod(curr_ymd, 100)
       call ESMF_TimeSet( currtime, yy=yr, mm=mon, dd=day, s=curr_tod, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Determine the current time from the lilac clock
       call ESMF_ClockGet(lilac_clock, currtime=clocktime, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Compute the time step difference from the current time from the restart file to the current lilac clock time
       ! (which is really just the start time)

       TimeStep_advance = currtime - clocktime
       call ESMF_TimeIntervalGet(timestep_advance, s=secs, rc=rc)
       if (mytask == 0) write(logunit,*)'DEBUG: time step advance is ',secs

       ! Advance the clock to the current time (in case of a restart)
       call ESMF_ClockAdvance (lilac_clock, timestep=timestep_advance, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in initializing restart alarm')

    end if

    ! Write out diagnostic info
    if (mytask == 0) then
       print *, trim(subname) // "---------------------------------------"
       call ESMF_CalendarPrint (lilac_calendar , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_ClockPrint (lilac_clock, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       print *, trim(subname) // "---------------------------------------"
    end if

    ! Add a restart alarm and stop alarm to the clock.
    !
    ! These alarms are initially set up to never go off, but they are turned on by
    ! arguments passed to lilac_run.
    !
    ! NOTE: The history alarm will be added in lilac_history_init and can go off multiple times during the run

    call lilac_time_alarmInit(lilac_clock, lilac_restart_alarm, 'lilac_restart_alarm', &
         option = optNever, opt_n = -1, rc = rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call lilac_time_alarmInit(lilac_clock, lilac_stop_alarm, 'lilac_stop_alarm', &
         option = optNever, opt_n = -1, rc = rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine lilac_time_clockInit

!===============================================================================

  subroutine lilac_time_alarmInit( clock, alarm, alarmname, option, opt_n, rc)

    ! Setup an alarm in a clock
    ! The ringtime sent to AlarmCreate MUST be the next alarm
    ! time.  If you send an arbitrary but proper ringtime from the
    ! past and the ring interval, the alarm will always go off on the
    ! next clock advance and this will cause serious problems.  Even
    ! if it makes sense to initialize an alarm with some reference
    ! time and the alarm interval, that reference time has to be
    ! advance forward to be >= the current time.  In the logic below
    ! we set an appropriate "NextAlarm" and then we make sure to
    ! advance it properly based on the ring interval.

    ! input/output variables
    type(ESMF_Clock) , intent(inout) :: clock     ! clock
    type(ESMF_Alarm) , intent(inout) :: alarm     ! alarm
    character(len=*) , intent(in)    :: alarmname ! alarm name
    character(len=*) , intent(in)    :: option    ! alarm option
    integer          , intent(in)    :: opt_n     ! alarm freq (ignored for option of optNone or optNever)
    integer          , intent(inout) :: rc        ! Return code

    ! local variables
    type(ESMF_Calendar)     :: cal              ! calendar
    integer                 :: lymd             ! local ymd
    integer                 :: ltod             ! local tod
    integer                 :: cyy,cmm,cdd,csec ! time info
    logical                 :: update_nextalarm ! update next alarm
    type(ESMF_Time)         :: CurrTime         ! Current Time
    type(ESMF_Time)         :: NextAlarm        ! Next restart alarm time
    type(ESMF_TimeInterval) :: AlarmInterval    ! Alarm interval
    integer                 :: sec
    character(len=*), parameter :: subname = '(lilac_time_alarmInit): '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(CurrTime, yy=cyy, mm=cmm, dd=cdd, s=csec, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! initial guess of next alarm, this will be updated below
    NextAlarm = CurrTime

    ! Get calendar from clock
    call ESMF_ClockGet(clock, calendar=cal)

    ! Determine inputs for call to create alarm
    if (trim(option) /= optNone .and. trim(option) /= optNever) then
       if (opt_n <= 0) call shr_sys_abort(subname//trim(option)//' invalid opt_n - must be > 0')
    end if
    select case (trim(option))

    case (optNone, optNever)
       call ESMF_TimeIntervalSet(AlarmInterval, yy=really_big_year, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=really_big_year, mm=12, dd=1, s=0, calendar=cal, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .false.

    case (optNSteps)
       call ESMF_ClockGet(clock, TimeStep=AlarmInterval, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNSeconds)
       call ESMF_TimeIntervalSet(AlarmInterval, s=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNMinutes)
       call ESMF_TimeIntervalSet(AlarmInterval, s=60, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNHours)
       call ESMF_TimeIntervalSet(AlarmInterval, s=3600, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNDays)
       call ESMF_TimeIntervalSet(AlarmInterval, d=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNMonths)
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNYears)
       call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case default
       call ESMF_LogWrite(subname//'unknown option '//trim(option), ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return

    end select

    ! -------------------------------------------------
    ! AlarmInterval and NextAlarm should be set 
    ! -------------------------------------------------

    ! advance Next Alarm so it won't ring on first timestep for
    ! most options above. go back one alarminterval just to be careful

    if (update_nextalarm) then
       NextAlarm = NextAlarm - AlarmInterval
       do while (NextAlarm <= CurrTime)
          NextAlarm = NextAlarm + AlarmInterval
       enddo
    endif

    alarm = ESMF_AlarmCreate( name=alarmname, clock=clock, ringTime=NextAlarm, &
         ringInterval=AlarmInterval, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//'created lilac alarm '//trim(alarmname), ESMF_LOGMSG_INFO)

  end subroutine lilac_time_alarmInit

!===============================================================================

  subroutine lilac_time_restart_write(clock, rc)

    ! -------------------------------------------------
    ! Write lilac restart time info
    ! -------------------------------------------------

    ! Input/output variables
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM)              :: vm
    type(ESMF_Time)            :: currtime, starttime, nexttime
    type(ESMF_TimeInterval)    :: timediff       ! Used to calculate curr_time
    type(ESMF_Calendar)        :: calendar
    character(len=64)          :: currtimestr, nexttimestr
    integer                    :: i,j,m,n,n1,ncnt
    integer                    :: curr_ymd       ! Current date YYYYMMDD
    integer                    :: curr_tod       ! Current time-of-day (s)
    integer                    :: start_ymd      ! Starting date YYYYMMDD
    integer                    :: start_tod      ! Starting time-of-day (s)
    integer                    :: next_ymd       ! Starting date YYYYMMDD
    integer                    :: next_tod       ! Starting time-of-day (s)
    integer                    :: nx,ny          ! global grid size
    integer                    :: yr,mon,day,sec ! time units
    real(R8)                   :: dayssince      ! Time interval since reference time
    integer                    :: unitn          ! unit number
    character(ESMF_MAXSTR)     :: time_units     ! units of time variable
    character(ESMF_MAXSTR)     :: restart_file   ! Local path to restart filename
    character(ESMF_MAXSTR)     :: restart_pfile  ! Local path to restart pointer filename
    character(ESMF_MAXSTR)     :: freq_option    ! freq_option setting (ndays, nsteps, etc)
    integer                    :: freq_n         ! freq_n setting relative to freq_option
    logical                    :: write_header   ! true => write netcdf header
    logical                    :: write_data     ! true => write netcdf data
    character(len=*), parameter :: subname='(lilac_time_phases_restart_write)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm=vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(clock, currtime=currtime, starttime=starttime, calendar=calendar, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGetNextTime(clock, nextTime=nexttime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(currtime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(currtimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
    call ESMF_LogWrite(trim(subname)//": currtime = "//trim(currtimestr), ESMF_LOGMSG_INFO, rc=rc)

    call ESMF_TimeGet(nexttime,yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(nexttimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
    call ESMF_LogWrite(trim(subname)//": nexttime = "//trim(nexttimestr), ESMF_LOGMSG_INFO, rc=rc)

    timediff = nexttime - starttime
    call ESMF_TimeIntervalGet(timediff, d=day, s=sec, rc=rc)
    dayssince = day + sec/real(SecPerDay,R8)

    call ESMF_TimeGet(starttime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    call shr_cal_ymd2date(yr,mon,day,start_ymd)
    start_tod = sec
    time_units = 'days since '//trim(lilac_io_date2yyyymmdd(start_ymd))//' '//lilac_io_sec2hms(start_tod, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(nexttime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    call shr_cal_ymd2date(yr,mon,day,next_ymd)
    next_tod = sec

    call ESMF_TimeGet(currtime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    call shr_cal_ymd2date(yr,mon,day,curr_ymd)
    curr_tod = sec

    !---------------------------------------
    ! Write restart file
    ! Use nexttimestr rather than currtimestr here since that is the time at the end of
    ! the timestep and is preferred for restart file names
    !---------------------------------------

    write(restart_file,"(5a)") trim(caseid),'.lilac','.r.', trim(nexttimestr),'.nc'
    call ESMF_LogWrite(subname//"lilac restart file is "//trim(restart_file), ESMF_LOGMSG_INFO)

    if (mytask == 0) then
       open(newunit=unitn, file="rpointer.lilac", form='FORMATTED')
       write(unitn,'(a)') trim(restart_file)
       close(unitn)
       call ESMF_LogWrite(trim(subname)//" wrote lilac restart pointer file rpointer.lilac", ESMF_LOGMSG_INFO)
    endif

    call ESMF_LogWrite(trim(subname)//": writing "//trim(restart_file), ESMF_LOGMSG_INFO)
    call lilac_io_wopen(restart_file, vm, mytask, clobber=.true.)

    do m = 1,2
       if (m == 1) then
          write_header = .true.  ;  write_data = .false.
       else 
          write_header = .false. ;  write_data = .true.
       endif

       if (write_data) then
          call lilac_io_enddef(restart_file)
       end if

       call ESMF_LogWrite(trim(subname)//": time "//trim(time_units), ESMF_LOGMSG_INFO)
       call lilac_io_write(restart_file, iam=mytask, &
            time_units=time_units, calendar=calendar, time_val=dayssince, &
            whead=write_header, wdata=write_data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Write out next ymd/tod in place of curr ymd/tod because
       ! the currently the restart represents the time at end of
       ! the current timestep and that is where we want to start the next run.

       call lilac_io_write(restart_file, mytask, next_ymd , 'curr_ymd' , whead=write_header, wdata=write_data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call lilac_io_write(restart_file, mytask, next_tod , 'curr_tod' , whead=write_header, wdata=write_data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call lilac_io_write(restart_file, mytask, start_ymd , 'start_ymd' , whead=write_header, wdata=write_data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call lilac_io_write(restart_file, mytask, start_tod , 'start_tod' , whead=write_header, wdata=write_data, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    call lilac_io_close(restart_file, mytask, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(trim(subname)//": closing "//trim(restart_file), ESMF_LOGMSG_INFO)

  end subroutine lilac_time_restart_write

  !===============================================================================

  subroutine lilac_time_restart_read(restart_file, start_ymd, start_tod, curr_ymd, curr_tod, rc)

    ! -------------------------------------------------
    ! Read the restart time info needed to initialize the clock
    ! -------------------------------------------------

    ! input/output variables
    character(len=*), intent(in) :: restart_file
    integer, intent(out)         :: start_ymd           ! Current ymd (YYYYMMDD)
    integer, intent(out)         :: start_tod           ! Current tod (seconds)
    integer, intent(out)         :: curr_ymd            ! Current ymd (YYYYMMDD)
    integer, intent(out)         :: curr_tod            ! Current tod (seconds)
    integer, intent(out)         :: rc

    ! local variables
    integer                 :: status, ncid, varid ! netcdf stuff
    character(CL)           :: tmpstr              ! temporary
    character(len=*), parameter :: subname = "(lilac_time_restart_read)"
    !----------------------------------------------------------------

    ! use netcdf here since it's serial
    status = nf90_open(restart_file, NF90_NOWRITE, ncid)
    if (status /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_open')

    status = nf90_inq_varid(ncid, 'curr_ymd', varid)
    if (status /= nf90_NoErr) call shr_sys_abort('ERROR: nf90_inq_varid curr_ymd')
    status = nf90_get_var(ncid, varid, curr_ymd)
    if (status /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_get_var curr_ymd')
    status = nf90_inq_varid(ncid, 'curr_tod', varid)
    if (status /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_inq_varid curr_tod')
    status = nf90_get_var(ncid, varid, curr_tod)
    if (status /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_get_var curr_tod')

    status = nf90_inq_varid(ncid, 'start_ymd', varid)
    if (status /= nf90_NoErr) call shr_sys_abort('ERROR: nf90_inq_varid start_ymd')
    status = nf90_get_var(ncid, varid, start_ymd)
    if (status /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_get_var start_ymd')
    status = nf90_inq_varid(ncid, 'start_tod', varid)
    if (status /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_inq_varid start_tod')
    status = nf90_get_var(ncid, varid, start_tod)
    if (status /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_get_var start_tod')


    status = nf90_close(ncid)
    if (status /= nf90_NoErr) call shr_sys_abort(' ERROR: nf90_close')

    write(tmpstr,*) trim(subname)//" read start_ymd  = ",start_ymd
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
    write(tmpstr,*) trim(subname)//" read start_tod  = ",start_tod
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

    write(tmpstr,*) trim(subname)//" read curr_ymd  = ",curr_ymd
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
    write(tmpstr,*) trim(subname)//" read curr_tod  = ",curr_tod
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

  end subroutine lilac_time_restart_read

end module lilac_time
