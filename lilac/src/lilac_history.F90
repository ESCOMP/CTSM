module lilac_history

  !-----------------------------------------------------------------------------
  ! Mediator Phases
  !-----------------------------------------------------------------------------

  use ESMF
  use shr_kind_mod    , only : cx=>shr_kind_cx, cs=>shr_kind_cs, cl=>shr_kind_cl, r8=>shr_kind_r8
  use lilac_constants , only : dbug_flag       => lilac_constants_dbug_flag
  use lilac_constants , only : SecPerDay       => lilac_constants_SecPerDay
  use lilac_methods   , only : FB_reset        => lilac_methods_FB_reset
  use lilac_methods   , only : FB_diagnose     => lilac_methods_FB_diagnose
  use lilac_methods   , only : FB_GetFldPtr    => lilac_methods_FB_GetFldPtr
  use lilac_methods   , only : FB_accum        => lilac_methods_FB_accum
  use lilac_methods   , only : chkerr
  use lilac_time      , only : alarmInit       => lilac_time_alarmInit
  use lilac_io        , only : lilac_io_write, lilac_io_wopen, lilac_io_enddef
  use lilac_io        , only : lilac_io_close, lilac_io_date2yyyymmdd, lilac_io_sec2hms
  use lilac_io        , only : lilac_io_ymd2date

  implicit none
  private

  public :: lilac_history_alarm_init
  public :: lilac_history_write

  type(ESMF_Alarm)        :: AlarmHist
  type(ESMF_Alarm)        :: AlarmHistAvg
  character(*), parameter :: u_FILE_u  = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine lilac_history_alarm_init(clock, hist_freq_n, hist_freq_option, rc)

    ! Initialize lilac history file alarm

    ! input/output variables
    type(ESMF_Clock)               :: clock  ! lilac clock
    integer          , intent(in)  :: hist_freq_n
    character(len=*) , intent(in)  :: hist_freq_option
    integer          , intent(out) :: rc

    ! local variables
    type(ESMF_Time)         :: currtime
    type(ESMF_Time)         :: reftime
    type(ESMF_Time)         :: starttime
    type(ESMF_Time)         :: nexttime
    type(ESMF_Calendar)     :: calendar       ! calendar type
    character(len=64)       :: currtimestr
    character(CS)           :: histavg_option ! Histavg option units
    integer                 :: yr,mon,day,sec ! time units
    character(CL)           :: freq_option    ! freq_option setting (ndays, nsteps, etc)
    integer                 :: freq_n         ! freq_n setting relative to freq_option
    integer                 :: iam
    character(len=*), parameter :: subname='(lilac_history_alarm_init)'
    !---------------------------------------

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)
    endif
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the clock info
    !---------------------------------------

    call ESMF_ClockGet(clock, currtime=currtime, reftime=reftime, starttime=starttime, calendar=calendar, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(currtime,yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(currtimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": currtime = "//trim(currtimestr), ESMF_LOGMSG_INFO, rc=rc)
    endif

    !---------------------------------------
    ! Initialize thie history alarm
    !---------------------------------------

    call alarmInit(clock, AlarmHist, option=freq_option, opt_n=freq_n, RefTime=RefTime, alarmname='history', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine lilac_history_alarm_init

  !===============================================================================

  subroutine lilac_history_write(atm2lnd_a_state, atm2lnd_l_state, lnd2atm_l_state, lnd2atm_a_state, clock, rc)

    ! Write lilac history file

    ! input/output variables
    type(ESMF_State) :: atm2lnd_a_state
    type(ESMF_State) :: atm2lnd_l_state
    type(ESMF_State) :: lnd2atm_l_state
    type(ESMF_State) :: lnd2atm_a_state
    type(ESMF_Clock) :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_FieldBundle)  :: c2a_fb , a2c_fb, c2l_fb, l2c_fb
    type(ESMF_VM)           :: vm
    type(ESMF_Time)         :: currtime
    type(ESMF_Time)         :: reftime
    type(ESMF_Time)         :: starttime
    type(ESMF_Time)         :: nexttime
    type(ESMF_TimeInterval) :: timediff       ! Used to calculate curr_time
    type(ESMF_Calendar)     :: calendar       ! calendar type
    character(len=64)       :: currtimestr
    character(len=64)       :: nexttimestr
    character(CS)           :: histavg_option ! Histavg option units
    integer                 :: i,j,m,n,n1,ncnt
    integer                 :: start_ymd      ! Starting date YYYYMMDD
    integer                 :: start_tod      ! Starting time-of-day (s)
    integer                 :: nx,ny          ! global grid size
    integer                 :: yr,mon,day,sec ! time units
    real(r8)                :: rval           ! real tmp value
    real(r8)                :: dayssince      ! Time interval since reference time
    integer                 :: fk             ! index
    character(CL)           :: time_units     ! units of time variable
    character(CL)           :: case_name      ! case name
    character(CL)           :: hist_file      ! Local path to history filename
    character(CS)           :: cpl_inst_tag   ! instance tag
    character(CL)           :: freq_option    ! freq_option setting (ndays, nsteps, etc)
    integer                 :: freq_n         ! freq_n setting relative to freq_option
    logical                 :: alarmIsOn      ! generic alarm flag
    real(r8)                :: tbnds(2)       ! CF1.0 time bounds
    logical                 :: whead,wdata    ! for writing restart/history cdf files
    integer                 :: dbrc
    integer                 :: iam
    logical,save            :: first_call = .true.
    character(len=*), parameter :: subname='(lilac_history_write)'
    !---------------------------------------

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)
    endif
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! --- Get the communicator and localpet
    !---------------------------------------

    call ESMF_VMGetGlobal(vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_VMGet(vm, localPet=iam, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! --- Get the clock info
    !---------------------------------------

    call ESMF_ClockGet(clock, currtime=currtime, reftime=reftime, starttime=starttime, calendar=calendar, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGetNextTime(clock, nextTime=nexttime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(currtime,yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(currtimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": currtime = "//trim(currtimestr), ESMF_LOGMSG_INFO, rc=rc)
    endif

    call ESMF_TimeGet(nexttime,yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(nexttimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": nexttime = "//trim(nexttimestr), ESMF_LOGMSG_INFO, rc=rc)
    endif
    timediff = nexttime - reftime
    call ESMF_TimeIntervalGet(timediff, d=day, s=sec, rc=rc)
    dayssince = day + sec/real(SecPerDay,R8)

    call ESMF_TimeGet(reftime, yy=yr, mm=mon, dd=day, s=sec, rc=rc)
    call lilac_io_ymd2date(yr,mon,day,start_ymd)
    start_tod = sec
    time_units = 'days since ' // trim(lilac_io_date2yyyymmdd(start_ymd)) // ' ' // lilac_io_sec2hms(start_tod, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! --- History Alarms
    !---------------------------------------

    ! if (ESMF_AlarmIsRinging(AlarmHist, rc=rc)) then
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !    alarmIsOn = .true.
    !    call ESMF_AlarmRingerOff( AlarmHist, rc=rc )
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! else
    !    alarmisOn = .false.
    ! endif
    ! hard-wire for now
    alarmisOn = .true.
    case_name = 'test_lilac'

    !---------------------------------------
    ! --- History File
    ! Use nexttimestr rather than currtimestr here since that is the time at the end of
    ! the timestep and is preferred for history file names
    !---------------------------------------

    if (alarmIsOn) then
       write(hist_file,"(6a)") trim(case_name), '.cpl.hi.',trim(nexttimestr),'.nc'
       call ESMF_LogWrite(trim(subname)//": write "//trim(hist_file), ESMF_LOGMSG_INFO, rc=rc)

       call lilac_io_wopen(hist_file, vm, iam, clobber=.true.)

       do m = 1,2
          whead=.false.
          wdata=.false.
          if (m == 1) then
             whead=.true.
          elseif (m == 2) then
             wdata=.true.
             call lilac_io_enddef(hist_file)
          endif

          tbnds = dayssince

          call ESMF_LogWrite(trim(subname)//": time "//trim(time_units), ESMF_LOGMSG_INFO, rc=rc)
          if (tbnds(1) >= tbnds(2)) then
             call lilac_io_write(hist_file, iam, &
                  time_units=time_units, calendar=calendar, time_val=dayssince, &
                  whead=whead, wdata=wdata, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             call lilac_io_write(hist_file, iam, &
                  time_units=time_units, calendar=calendar, time_val=dayssince, &
                  whead=whead, wdata=wdata, tbnds=tbnds, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          endif

          nx = 72 ! hard-wire for now
          ny = 46 ! hard-wire for now

          call ESMF_StateGet(atm2lnd_a_state, 'a2c_fb', a2c_fb) ! from atm
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call lilac_io_write(hist_file, iam, a2c_fb, &
               nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre='a2c_from_atm', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          call ESMF_StateGet(atm2lnd_l_state, 'c2l_fb', c2l_fb) ! to land
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call lilac_io_write(hist_file, iam, c2l_fb, &
               nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre='c2l_to_land', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          call ESMF_StateGet(lnd2atm_l_state, 'l2c_fb', l2c_fb)  ! from land 
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call lilac_io_write(hist_file, iam, c2l_fb, &
               nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre='l2c_from_land', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          call ESMF_StateGet(lnd2atm_a_state, 'c2a_fb', c2a_fb) ! to atm
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call lilac_io_write(hist_file, iam, c2l_fb, &
               nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre='c2a_to_atm', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       enddo

       call lilac_io_close(hist_file, iam, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    endif

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    first_call = .false.

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=rc)
    endif

  end subroutine lilac_history_write

  !===============================================================================

end module lilac_history
