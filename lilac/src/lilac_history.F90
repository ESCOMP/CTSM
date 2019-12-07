module lilac_history

  !-----------------------------------------------------------------------------
  ! LILAC history output
  !-----------------------------------------------------------------------------

  use ESMF
  use shr_kind_mod    , only : cx=>shr_kind_cx, cs=>shr_kind_cs, cl=>shr_kind_cl, r8=>shr_kind_r8
  use shr_sys_mod     , only : shr_sys_abort
  use lilac_atmcap    , only : atm_nx       => atm_global_nx
  use lilac_atmcap    , only : atm_ny       => atm_global_ny
  use lilac_constants , only : dbug_flag    => lilac_constants_dbug_flag
  use lilac_constants , only : SecPerDay    => lilac_constants_SecPerDay
  use lilac_time      , only : alarmInit    => lilac_time_alarmInit
  use lilac_methods   , only : FB_reset     => lilac_methods_FB_reset
  use lilac_methods   , only : FB_diagnose  => lilac_methods_FB_diagnose
  use lilac_methods   , only : FB_GetFldPtr => lilac_methods_FB_GetFldPtr
  use lilac_methods   , only : FB_accum     => lilac_methods_FB_accum
  use lilac_methods   , only : chkerr
  use lilac_io        , only : lilac_io_write, lilac_io_wopen, lilac_io_enddef
  use lilac_io        , only : lilac_io_close, lilac_io_date2yyyymmdd, lilac_io_sec2hms
  use lilac_io        , only : lilac_io_ymd2date

  implicit none
  private

  public :: lilac_history_init
  public :: lilac_history_write

  character(CL)           :: histfile_prefix
  type(ESMF_Alarm)        :: AlarmHist
  type(ESMF_Alarm)        :: AlarmHistAvg
  character(*), parameter :: u_FILE_u  = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine lilac_history_init(clock, rc)

    ! ------------------------------------------
    ! Initialize lilac history file alarm
    ! ------------------------------------------

    ! input/output variables
    type(ESMF_Clock)               :: clock  ! lilac clock
    integer          , intent(out) :: rc

    ! local variables
    type(ESMF_Time)         :: currtime
    type(ESMF_Time)         :: reftime
    type(ESMF_Time)         :: starttime
    type(ESMF_Time)         :: nexttime
    type(ESMF_Calendar)     :: calendar       ! calendar type
    character(len=64)       :: currtimestr
    integer                 :: yr,mon,day,sec ! time units
    character(CL)           :: freq_option    ! freq_option setting (ndays, nsteps, etc)
    integer                 :: freq_n         ! freq_n setting relative to freq_option
    integer                 :: iam
    integer                 :: fileunit
    integer                 :: ierr
    character(CL)           :: caseid
    character(CS)           :: starttype
    character(CS)           :: lilac_histfreq_option
    integer                 :: lilac_histfreq_n
    character(len=*), parameter :: subname='(lilac_history_init)'
    !---------------------------------------

    namelist /lilac_run_input/ caseid, starttype
    namelist /lilac_io_input/ lilac_histfreq_n, lilac_histfreq_option

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)
    endif
    rc = ESMF_SUCCESS

    ! read in caseid
    open(newunit=fileunit, status="old", file="lilac_in")
    read(fileunit, lilac_run_input, iostat=ierr)
    if (ierr > 0) then
       call shr_sys_abort(trim(subname) // 'error reading in lilac_run_input')
    end if
    close(fileunit)
    write(histfile_prefix,"(2a)") trim(caseid),'.lilac.hi.'

    ! read in history file output frequencies
    open(newunit=fileunit, status="old", file="lilac_in")
    read(fileunit, lilac_io_input, iostat=ierr)
    if (ierr > 0) then
       call shr_sys_abort(trim(subname) // 'error reading in lilac_io_input')
    end if
    close(fileunit)

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

    call alarmInit(clock, AlarmHist, option=lilac_histfreq_option, opt_n=lilac_histfreq_n, &
         RefTime=RefTime, alarmname='lilac_history_alarm', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine lilac_history_init

  !===============================================================================

  subroutine lilac_history_write(atm2cpl_state, cpl2atm_state, &
       lnd2cpl_state, cpl2lnd_state, rof2cpl_state, cpl2rof_state, clock, rc)

    ! ------------------------------
    ! Write lilac history file
    ! ------------------------------

    ! input/output variables
    type(ESMF_State)           :: atm2cpl_state
    type(ESMF_State)           :: cpl2atm_state
    type(ESMF_State)           :: lnd2cpl_state
    type(ESMF_State)           :: cpl2lnd_state
    type(ESMF_State), optional :: rof2cpl_state
    type(ESMF_State), optional :: cpl2rof_state
    type(ESMF_Clock)           :: clock
    integer, intent(out)       :: rc

    ! local variables
    type(ESMF_FieldBundle)      :: c2a_fb, a2c_fb
    type(ESMF_FieldBundle)      :: c2l_fb_atm, c2l_fb_rof, l2c_fb_atm, l2c_fb_rof
    type(ESMF_FieldBundle)      :: c2r_fb, r2c_fb
    type(ESMF_VM)               :: vm
    type(ESMF_Time)             :: currtime
    type(ESMF_Time)             :: reftime
    type(ESMF_Time)             :: starttime
    type(ESMF_Time)             :: nexttime
    type(ESMF_TimeInterval)     :: timediff       ! Used to calculate curr_time
    type(ESMF_Calendar)         :: calendar       ! calendar type
    character(len=CS)           :: currtimestr
    character(len=CS)           :: nexttimestr
    integer                     :: i,j,m,n
    integer                     :: start_ymd      ! Starting date YYYYMMDD
    integer                     :: start_tod      ! Starting time-of-day (s)
    integer                     :: nx,ny          ! global grid size
    integer                     :: yr,mon,day,sec ! time units
    real(r8)                    :: rval           ! real tmp value
    real(r8)                    :: dayssince      ! Time interval since reference time
    character(CL)               :: time_units     ! units of time variable
    character(CL)               :: hist_file      ! Local path to history filename
    character(CL)               :: freq_option    ! freq_option setting (ndays, nsteps, etc)
    integer                     :: freq_n         ! freq_n setting relative to freq_option
    logical                     :: alarmIsOn      ! generic alarm flag
    real(r8)                    :: tbnds(2)       ! CF1.0 time bounds
    logical                     :: whead,wdata    ! for writing restart/history cdf files
    character(CL)               :: cvalue
    integer                     :: lnd_nx, lnd_ny
    integer                     :: rof_nx, rof_ny
    integer                     :: iam
    logical,save                :: first_call = .true.
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

    if (ESMF_AlarmIsRinging(AlarmHist, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       alarmIsOn = .true.
       call ESMF_AlarmRingerOff( AlarmHist, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       alarmIsOn = .false.
    endif

    ! TODO: remove this hard-wiring
    alarmisOn = .true.

    !---------------------------------------
    ! --- History File
    ! Use nexttimestr rather than currtimestr here since that is the time at the end of
    ! the timestep and is preferred for history file names
    !---------------------------------------

    if (alarmIsOn) then
       write(hist_file,"(3a)") trim(histfile_prefix),trim(nexttimestr),'.nc'
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

          ! obtain global sizes for lnd_nx and lnd_ny
          call ESMF_AttributeGet(lnd2cpl_state, name="lnd_nx", value=cvalue, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
          read(cvalue,*) lnd_nx
          call ESMF_LogWrite(subname//"got attribute lnd_nx "//trim(cvalue), ESMF_LOGMSG_INFO)
          call ESMF_AttributeGet(lnd2cpl_state, name="lnd_ny", value=cvalue, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
          read(cvalue,*) lnd_ny
          call ESMF_LogWrite(subname//"got attribute lnd_nx "//trim(cvalue), ESMF_LOGMSG_INFO)

          call ESMF_StateGet(atm2cpl_state, 'a2c_fb', a2c_fb)         ! from atm
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call lilac_io_write(hist_file, iam, a2c_fb, &
               nx=atm_nx, ny=atm_ny, nt=1, whead=whead, wdata=wdata, pre='atm_to_cpl', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          call ESMF_StateGet(cpl2atm_state, 'c2a_fb', c2a_fb)         ! to atm
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call lilac_io_write(hist_file, iam, c2a_fb, &
               nx=atm_nx, ny=atm_ny, nt=1, whead=whead, wdata=wdata, pre='cpl_to_atm', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          call ESMF_StateGet(lnd2cpl_state, 'l2c_fb_atm', l2c_fb_atm) ! from lnd
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call lilac_io_write(hist_file, iam, l2c_fb_atm, &
               nx=lnd_nx, ny=lnd_ny, nt=1, whead=whead, wdata=wdata, pre='lnd_to_cpl_atm', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          call ESMF_StateGet(lnd2cpl_state, 'l2c_fb_rof', l2c_fb_rof) ! from lnd
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call lilac_io_write(hist_file, iam, l2c_fb_rof, &
               nx=lnd_nx, ny=lnd_ny, nt=1, whead=whead, wdata=wdata, pre='lnd_to_cpl_rof', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          call ESMF_StateGet(cpl2lnd_state, 'c2l_fb_atm', c2l_fb_atm) ! to lnd
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call lilac_io_write(hist_file, iam, c2l_fb_atm, &
               nx=lnd_nx, ny=lnd_ny, nt=1, whead=whead, wdata=wdata, pre='cpl_to_lnd_atm', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          if (present(rof2cpl_state) .and. present(cpl2rof_state)) then
             ! obtain global sizes for rof_nx and rof_ny
             call ESMF_AttributeGet(rof2cpl_state, name="rof_nx", value=cvalue, rc=rc)
             if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
             read(cvalue,*) rof_nx
             call ESMF_LogWrite(subname//"got attribute rof_nx "//trim(cvalue), ESMF_LOGMSG_INFO)
             call ESMF_AttributeGet(rof2cpl_state, name="rof_ny", value=cvalue, rc=rc)
             if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
             read(cvalue,*) rof_ny
             call ESMF_LogWrite(subname//"got attribute rof_nx "//trim(cvalue), ESMF_LOGMSG_INFO)

             call ESMF_StateGet(rof2cpl_state, 'r2c_fb', r2c_fb)      ! from rof
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call lilac_io_write(hist_file, iam, r2c_fb, &
                  nx=rof_nx, ny=rof_ny, nt=1, whead=whead, wdata=wdata, pre='rof_to_cpl', rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             call ESMF_StateGet(cpl2rof_state, 'c2r_fb', c2r_fb)      ! to rof 
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call lilac_io_write(hist_file, iam, c2r_fb, &
                  nx=rof_nx, ny=rof_ny, nt=1, whead=whead, wdata=wdata, pre='cpl_to_rof', rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             call ESMF_StateGet(cpl2lnd_state, 'c2l_fb_rof', c2l_fb_rof) ! to rof 
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call lilac_io_write(hist_file, iam, c2l_fb_rof, &
                  nx=lnd_nx, ny=lnd_ny, nt=1, whead=.true., wdata=wdata, pre='cpl_to_lnd_rof', rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

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

end module lilac_history
