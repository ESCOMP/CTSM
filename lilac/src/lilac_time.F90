module lilac_time

  use ESMF
  use shr_kind_mod    , only : cx=>shr_kind_cx, cs=>shr_kind_cs, cl=>shr_kind_cl, r8=>shr_kind_r8
  use lilac_constants , only : dbug_flag => lilac_constants_dbug_flag
  use lilac_methods   , only : chkerr 

  implicit none
  private    ! default private

  public  :: lilac_time_alarmInit  ! initialize an alarm

  ! Clock and alarm options
  character(len=*), private, parameter :: &
       optNONE           = "none"      , &
       optNever          = "never"     , &
       optNSteps         = "nsteps"    , &
       optNStep          = "nstep"     , &
       optNSeconds       = "nseconds"  , &
       optNSecond        = "nsecond"   , &
       optNMinutes       = "nminutes"  , &
       optNMinute        = "nminute"   , &
       optNHours         = "nhours"    , &
       optNHour          = "nhour"     , &
       optNDays          = "ndays"     , &
       optNDay           = "nday"      , &
       optNMonths        = "nmonths"   , &
       optNMonth         = "nmonth"    , &
       optNYears         = "nyears"    , &
       optNYear          = "nyear"     , &
       optMonthly        = "monthly"   , &
       optYearly         = "yearly"    , &
       optIfdays0        = "ifdays0"   , &
       optGLCCouplingPeriod = "glc_coupling_period"

  ! Module data
  integer, parameter          :: SecPerDay = 86400 ! Seconds per day
  character(len=*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine lilac_time_alarmInit( clock, alarm, option, &
       opt_n, opt_ymd, opt_tod, RefTime, alarmname, rc)

    ! !DESCRIPTION: Setup an alarm in a clock
    ! Notes: The ringtime sent to AlarmCreate MUST be the next alarm
    ! time.  If you send an arbitrary but proper ringtime from the
    ! past and the ring interval, the alarm will always go off on the
    ! next clock advance and this will cause serious problems.  Even
    ! if it makes sense to initialize an alarm with some reference
    ! time and the alarm interval, that reference time has to be
    ! advance forward to be >= the current time.  In the logic below
    ! we set an appropriate "NextAlarm" and then we make sure to
    ! advance it properly based on the ring interval.

    ! input/output variables
    type(ESMF_Clock)            , intent(inout) :: clock     ! clock
    type(ESMF_Alarm)            , intent(inout) :: alarm     ! alarm
    character(len=*)            , intent(in)    :: option    ! alarm option
    integer          , optional , intent(in)    :: opt_n     ! alarm freq
    integer          , optional , intent(in)    :: opt_ymd   ! alarm ymd
    integer          , optional , intent(in)    :: opt_tod   ! alarm tod (sec)
    type(ESMF_Time)  , optional , intent(in)    :: RefTime   ! ref time
    character(len=*) , optional , intent(in)    :: alarmname ! alarm name
    integer                     , intent(inout) :: rc        ! Return code

    ! local variables
    type(ESMF_Calendar)     :: cal              ! calendar
    integer                 :: lymd             ! local ymd
    integer                 :: ltod             ! local tod
    integer                 :: cyy,cmm,cdd,csec ! time info
    character(len=64)       :: lalarmname       ! local alarm name
    logical                 :: update_nextalarm ! update next alarm
    type(ESMF_Time)         :: CurrTime         ! Current Time
    type(ESMF_Time)         :: NextAlarm        ! Next restart alarm time
    type(ESMF_TimeInterval) :: AlarmInterval    ! Alarm interval
    integer                 :: sec
    character(len=*), parameter :: subname = '(lilac_time_alarmInit): '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lalarmname = 'alarm_unknown'
    if (present(alarmname)) lalarmname = trim(alarmname)
    ltod = 0
    if (present(opt_tod)) ltod = opt_tod
    lymd = -1
    if (present(opt_ymd)) lymd = opt_ymd

    call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(CurrTime, yy=cyy, mm=cmm, dd=cdd, s=csec, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! initial guess of next alarm, this will be updated below
    if (present(RefTime)) then
       NextAlarm = RefTime
    else
       NextAlarm = CurrTime
    endif

    ! Get calendar from clock
    call ESMF_ClockGet(clock, calendar=cal)

    ! Determine inputs for call to create alarm
    selectcase (trim(option))

    case (optNONE)
       call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=9999, mm=12, dd=1, s=0, calendar=cal, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .false.

    case (optNever)
       call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=9999, mm=12, dd=1, s=0, calendar=cal, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .false.

    case (optIfdays0)
       if (.not. present(opt_ymd)) then
          call ESMF_LogWrite(subname//trim(option)//' requires opt_ymd', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(subname//trim(option)//' requires opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0)  then
          call ESMF_LogWrite(subname//trim(option)//' invalid opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=opt_n, s=0, calendar=cal, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .true.

    case (optNSteps)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(subname//trim(option)//' requires opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(subname//trim(option)//' invalid opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_ClockGet(clock, TimeStep=AlarmInterval, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNStep)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(subname//trim(option)//' requires opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0)  then
          call ESMF_LogWrite(subname//trim(option)//' invalid opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_ClockGet(clock, TimeStep=AlarmInterval, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNSeconds)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(subname//trim(option)//' requires opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(subname//trim(option)//' invalid opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNSecond)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(subname//trim(option)//' requires opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(subname//trim(option)//' invalid opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNMinutes)
       call ESMF_TimeIntervalSet(AlarmInterval, s=60, rc=rc)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(subname//trim(option)//' requires opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(subname//trim(option)//' invalid opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNMinute)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(subname//trim(option)//' requires opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(subname//trim(option)//' invalid opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=60, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNHours)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(subname//trim(option)//' requires opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(subname//trim(option)//' invalid opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=3600, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNHour)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(subname//trim(option)//' requires opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(subname//trim(option)//' invalid opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=3600, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNDays)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(subname//trim(option)//' requires opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(subname//trim(option)//' invalid opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, d=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNDay)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(subname//trim(option)//' requires opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(subname//trim(option)//' invalid opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, d=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNMonths)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(subname//trim(option)//' requires opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(subname//trim(option)//' invalid opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNMonth)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(subname//trim(option)//' requires opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(subname//trim(option)//' invalid opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optMonthly)
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=1, s=0, calendar=cal, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .true.

    case (optNYears)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(subname//trim(option)//' requires opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(subname//trim(option)//' invalid opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNYear)
       if (.not.present(opt_n)) then
          call ESMF_LogWrite(subname//trim(option)//' requires opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       if (opt_n <= 0) then
          call ESMF_LogWrite(subname//trim(option)//' invalid opt_n', ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optYearly)
       call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=1, dd=1, s=0, calendar=cal, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .true.

    case default
       call ESMF_LogWrite(subname//'unknown option '//trim(option), ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return

    end select

    ! --------------------------------------------------------------------------------
    ! --- AlarmInterval and NextAlarm should be set ---
    ! --------------------------------------------------------------------------------

    ! --- advance Next Alarm so it won't ring on first timestep for
    ! --- most options above. go back one alarminterval just to be careful

    if (update_nextalarm) then
       NextAlarm = NextAlarm - AlarmInterval
       do while (NextAlarm <= CurrTime)
          NextAlarm = NextAlarm + AlarmInterval
       enddo
    endif

    alarm = ESMF_AlarmCreate( name=lalarmname, clock=clock, ringTime=NextAlarm, &
         ringInterval=AlarmInterval, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine lilac_time_alarmInit

  !===============================================================================

  subroutine lilac_time_read_restart(restart_file, &
       start_ymd, start_tod, ref_ymd, ref_tod, curr_ymd, curr_tod, rc)

    use netcdf , only : nf90_open, nf90_nowrite, nf90_noerr
    use netcdf , only : nf90_inq_varid, nf90_get_var, nf90_close
    use ESMF   , only : ESMF_LogWrite, ESMF_LOGMSG_INFO

    ! input/output variables
    character(len=*), intent(in) :: restart_file
    integer, intent(out)         :: ref_ymd             ! Reference date (YYYYMMDD)
    integer, intent(out)         :: ref_tod             ! Reference time of day (seconds)
    integer, intent(out)         :: start_ymd           ! Start date (YYYYMMDD)
    integer, intent(out)         :: start_tod           ! Start time of day (seconds)
    integer, intent(out)         :: curr_ymd            ! Current ymd (YYYYMMDD)
    integer, intent(out)         :: curr_tod            ! Current tod (seconds)
    integer, intent(out)         :: rc

    ! local variables
    integer                 :: status, ncid, varid ! netcdf stuff
    character(CL)           :: tmpstr              ! temporary
    character(len=*), parameter :: subname = "(lilac_time_read_restart)"
    !----------------------------------------------------------------

    ! use netcdf here since it's serial
    status = nf90_open(restart_file, NF90_NOWRITE, ncid)
    if (status /= nf90_NoErr) then
       print *,__FILE__,__LINE__,trim(restart_file)
       call ESMF_LogWrite(trim(subname)//' ERROR: nf90_open', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    endif
    status = nf90_inq_varid(ncid, 'start_ymd', varid)
    if (status /= nf90_NoErr) then
       call ESMF_LogWrite(trim(subname)//' ERROR: nf90_inq_varid start_ymd', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if
    status = nf90_get_var(ncid, varid, start_ymd)
    if (status /= nf90_NoErr) then
       call ESMF_LogWrite(trim(subname)//' ERROR: nf90_get_var start_ymd', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if
    status = nf90_inq_varid(ncid, 'start_tod', varid)
    if (status /= nf90_NoErr) then
       call ESMF_LogWrite(trim(subname)//' ERROR: nf90_inq_varid start_tod', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if
    status = nf90_get_var(ncid, varid, start_tod)
    if (status /= nf90_NoErr) then
       call ESMF_LogWrite(trim(subname)//' ERROR: nf90_get_var start_tod', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if
    status = nf90_inq_varid(ncid, 'ref_ymd', varid)
    if (status /= nf90_NoErr) then
       call ESMF_LogWrite(trim(subname)//' ERROR: nf90_inq_varid ref_ymd', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if
    status = nf90_get_var(ncid, varid, ref_ymd)
    if (status /= nf90_NoErr) then
       call ESMF_LogWrite(trim(subname)//' ERROR: nf90_get_var ref_ymd', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if
    status = nf90_inq_varid(ncid, 'ref_tod', varid)
    if (status /= nf90_NoErr) then
       call ESMF_LogWrite(trim(subname)//' ERROR: nf90_inq_varid ref_tod', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if
    status = nf90_get_var(ncid, varid, ref_tod)
    if (status /= nf90_NoErr) then
       call ESMF_LogWrite(trim(subname)//' ERROR: nf90_get_var ref_tod', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if
    status = nf90_inq_varid(ncid, 'curr_ymd', varid)
    if (status /= nf90_NoErr) then
       call ESMF_LogWrite(trim(subname)//' ERROR: nf90_inq_varid curr_ymd', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if
    status = nf90_get_var(ncid, varid, curr_ymd)
    if (status /= nf90_NoErr) then
       call ESMF_LogWrite(trim(subname)//' ERROR: nf90_get_var curr_ymd', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if
    status = nf90_inq_varid(ncid, 'curr_tod', varid)
    if (status /= nf90_NoErr) then
       call ESMF_LogWrite(trim(subname)//' ERROR: nf90_inq_varid curr_tod', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if
    status = nf90_get_var(ncid, varid, curr_tod)
    if (status /= nf90_NoErr) then
       call ESMF_LogWrite(trim(subname)//' ERROR: nf90_get_var curr_tod', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if
    status = nf90_close(ncid)
    if (status /= nf90_NoErr) then
       call ESMF_LogWrite(trim(subname)//' ERROR: nf90_close', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if

    write(tmpstr,*) trim(subname)//" read start_ymd = ",start_ymd
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
    write(tmpstr,*) trim(subname)//" read start_tod = ",start_tod
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
    write(tmpstr,*) trim(subname)//" read ref_ymd   = ",ref_ymd
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
    write(tmpstr,*) trim(subname)//" read ref_tod   = ",ref_tod
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
    write(tmpstr,*) trim(subname)//" read curr_ymd  = ",curr_ymd
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
    write(tmpstr,*) trim(subname)//" read curr_tod  = ",curr_tod
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

  end subroutine lilac_time_read_restart

end module lilac_time


