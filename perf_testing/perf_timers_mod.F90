module perf_timers_mod

  !-----------------------------------------------------------------------
  ! Lightweight per-label wall-clock timer for the standalone perf-testing
  ! drivers in perf_testing/. Use:
  !
  !   call perf_timer_start('init_sminn_tot')
  !   ! ... loop body ...
  !   call perf_timer_stop('init_sminn_tot')
  !
  !   call perf_timer_print(6)            ! table to stdout
  !   call perf_timer_dump_csv(unit_csv)  ! one row per label
  !
  ! Gated by the cpp macro INNER_TIMING. When undefined, all public
  ! routines are empty no-ops and there is zero per-call overhead in
  ! the science code. Backed by the Fortran intrinsic system_clock
  ! (no external library dependencies).
  !-----------------------------------------------------------------------

  implicit none
  private

  public :: perf_timer_start
  public :: perf_timer_stop
  public :: perf_timer_print
  public :: perf_timer_dump_csv
  public :: perf_timer_reset

#ifdef INNER_TIMING

  integer, parameter :: r8 = selected_real_kind(12)
  integer, parameter :: max_timers = 64
  integer, parameter :: max_label  = 32

  type :: timer_t
     character(len=max_label) :: label   = ''
     integer(kind=8)          :: t_start = 0_8
     integer(kind=8)          :: t_total = 0_8
     integer                  :: ncalls  = 0
     logical                  :: in_use  = .false.
  end type timer_t

  type(timer_t)  , save :: timers(max_timers)
  integer        , save :: ntimers = 0
  integer(kind=8), save :: t_rate  = 0_8
  logical        , save :: rate_initialized = .false.

#endif

contains

  !-----------------------------------------------------------------------
  subroutine perf_timer_start(label)
    character(len=*), intent(in) :: label
#ifdef INNER_TIMING
    integer :: idx
    if (.not. rate_initialized) then
       call system_clock(count_rate=t_rate)
       rate_initialized = .true.
    end if
    idx = find_or_add_timer(label)
    call system_clock(timers(idx)%t_start)
#endif
  end subroutine perf_timer_start

  !-----------------------------------------------------------------------
  subroutine perf_timer_stop(label)
    character(len=*), intent(in) :: label
#ifdef INNER_TIMING
    integer         :: idx
    integer(kind=8) :: t_now
    call system_clock(t_now)
    idx = find_or_add_timer(label)
    timers(idx)%t_total = timers(idx)%t_total + (t_now - timers(idx)%t_start)
    timers(idx)%ncalls  = timers(idx)%ncalls + 1
#endif
  end subroutine perf_timer_stop

  !-----------------------------------------------------------------------
  subroutine perf_timer_print(unit)
    integer, intent(in) :: unit
#ifdef INNER_TIMING
    integer  :: i
    real(r8) :: total_s, per_call_s
    write(unit,'(a)') '--- per-loop wall-clock timers ---'
    write(unit,'(a32,2x,a14,2x,a10,2x,a14)') &
         'label', 'total (s)', 'ncalls', 'per call (s)'
    do i = 1, ntimers
       if (.not. timers(i)%in_use) cycle
       total_s = real(timers(i)%t_total, r8) / real(t_rate, r8)
       if (timers(i)%ncalls > 0) then
          per_call_s = total_s / real(timers(i)%ncalls, r8)
       else
          per_call_s = 0.0_r8
       end if
       write(unit,'(a32,2x,es14.6,2x,i10,2x,es14.6)') &
            trim(timers(i)%label), total_s, timers(i)%ncalls, per_call_s
    end do
#endif
  end subroutine perf_timer_print

  !-----------------------------------------------------------------------
  subroutine perf_timer_dump_csv(unit)
    integer, intent(in) :: unit
#ifdef INNER_TIMING
    integer  :: i
    real(r8) :: total_s, per_call_s
    write(unit,'(a)') 'label,total_s,ncalls,per_call_s'
    do i = 1, ntimers
       if (.not. timers(i)%in_use) cycle
       total_s = real(timers(i)%t_total, r8) / real(t_rate, r8)
       if (timers(i)%ncalls > 0) then
          per_call_s = total_s / real(timers(i)%ncalls, r8)
       else
          per_call_s = 0.0_r8
       end if
       write(unit,'(a,",",es24.16,",",i0,",",es24.16)') &
            trim(timers(i)%label), total_s, timers(i)%ncalls, per_call_s
    end do
#endif
  end subroutine perf_timer_dump_csv

  !-----------------------------------------------------------------------
  subroutine perf_timer_reset()
#ifdef INNER_TIMING
    integer :: i
    do i = 1, max_timers
       timers(i)%label   = ''
       timers(i)%t_start = 0_8
       timers(i)%t_total = 0_8
       timers(i)%ncalls  = 0
       timers(i)%in_use  = .false.
    end do
    ntimers = 0
#endif
  end subroutine perf_timer_reset

#ifdef INNER_TIMING

  !-----------------------------------------------------------------------
  function find_timer(label) result(idx)
    character(len=*), intent(in) :: label
    integer :: idx
    integer :: i
    idx = 0
    do i = 1, ntimers
       if (timers(i)%in_use .and. trim(timers(i)%label) == trim(label)) then
          idx = i
          return
       end if
    end do
  end function find_timer

  !-----------------------------------------------------------------------
  function find_or_add_timer(label) result(idx)
    character(len=*), intent(in) :: label
    integer :: idx
    idx = find_timer(label)
    if (idx == 0) then
       if (ntimers >= max_timers) then
          write(*,'(a,a,a,i0,a)') &
               'perf_timers_mod: too many timers (adding "', trim(label), &
               '"); raise max_timers (currently ', max_timers, ')'
          stop 1
       end if
       ntimers = ntimers + 1
       timers(ntimers)%label  = label
       timers(ntimers)%in_use = .true.
       idx = ntimers
    end if
  end function find_or_add_timer

#endif

end module perf_timers_mod
