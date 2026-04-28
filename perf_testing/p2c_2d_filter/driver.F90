program p2c_2d_filter_driver

  !-----------------------------------------------------------------------
  ! Standalone timing harness for p2c_2d_filter.
  !
  ! Allocates synthetic patch- and column-level arrays, calls the routine
  ! niters times under a system_clock wrapper, prints config + timings +
  ! checksum, and writes the result to last_run.txt. If baseline_checksum.txt
  ! is present and its parameter set matches the current run, compares the
  ! checksum and prints MATCH or MISMATCH.
  !
  ! Usage:
  !   ./driver [ncol [npatch_per_col [nlev [numfc [niters]]]]]
  !
  ! Defaults (when args omitted): 8000 16 10 8000 1
  !-----------------------------------------------------------------------

  use p2c_2d_filter_mod, only : r8, p2c_2d_filter

  implicit none

  ! Run parameters
  integer :: ncol           = 8000
  integer :: npatch_per_col = 16
  integer :: nlev           = 10
  integer :: numfc          = -1     ! -1 sentinel -> default to ncol
  integer :: niters         = 1

  ! Arrays (pointer where the routine expects pointer)
  integer , pointer :: patchi(:)        => null()
  integer , pointer :: patchf(:)        => null()
  logical , pointer :: active(:)        => null()
  real(r8), pointer :: wtcol(:)         => null()
  real(r8), pointer :: patcharr(:,:)    => null()
  real(r8), pointer :: colarr(:,:)      => null()
  integer , allocatable :: filterc(:)

  integer  :: npatch_total
  integer  :: c, p, j, iter
  integer(kind=8) :: t_start, t_end, t_rate
  real(r8) :: elapsed_s, per_call_s, checksum

  call parse_args()
  if (numfc < 0) numfc = ncol

  npatch_total = ncol * npatch_per_col

  call allocate_arrays()
  call fill_inputs()

  ! Time niters calls
  call system_clock(count_rate=t_rate)
  call system_clock(t_start)
  do iter = 1, niters
     call p2c_2d_filter(nlev, numfc, filterc, &
                        patchi, patchf,       &
                        active, wtcol,        &
                        patcharr, colarr)
  end do
  call system_clock(t_end)

  elapsed_s  = real(t_end - t_start, r8) / real(t_rate, r8)
  per_call_s = elapsed_s / real(niters, r8)
  checksum   = sum(colarr)

  call report(elapsed_s, per_call_s, checksum)
  call write_last_run(checksum)
  call compare_to_baseline(checksum)

contains

  !---------------------------------------------------------------------
  subroutine parse_args()
    integer :: nargs
    character(len=64) :: arg

    nargs = command_argument_count()
    if (nargs >= 1) then; call get_command_argument(1, arg); read(arg,*) ncol;           end if
    if (nargs >= 2) then; call get_command_argument(2, arg); read(arg,*) npatch_per_col; end if
    if (nargs >= 3) then; call get_command_argument(3, arg); read(arg,*) nlev;           end if
    if (nargs >= 4) then; call get_command_argument(4, arg); read(arg,*) numfc;          end if
    if (nargs >= 5) then; call get_command_argument(5, arg); read(arg,*) niters;         end if
  end subroutine parse_args

  !---------------------------------------------------------------------
  subroutine allocate_arrays()
    allocate(patchi(ncol))
    allocate(patchf(ncol))
    allocate(active(npatch_total))
    allocate(wtcol(npatch_total))
    allocate(patcharr(npatch_total, nlev))
    allocate(colarr(ncol, nlev))
    allocate(filterc(numfc))
  end subroutine allocate_arrays

  !---------------------------------------------------------------------
  subroutine fill_inputs()
    integer :: c, p, j, fc

    ! Contiguous patch ranges per column: column c owns patches
    !   (c-1)*npatch_per_col + 1 .. c*npatch_per_col
    do c = 1, ncol
       patchi(c) = (c - 1) * npatch_per_col + 1
       patchf(c) =  c      * npatch_per_col
    end do

    ! All patches active, equal weights
    active = .true.
    wtcol  = 1.0_r8 / real(npatch_per_col, r8)

    ! Reproducible non-trivial input data
    do j = 1, nlev
       do p = 1, npatch_total
          patcharr(p, j) = real(p + j, r8)
       end do
    end do

    ! Filter selects the first numfc columns
    do fc = 1, numfc
       filterc(fc) = fc
    end do

    colarr = 0.0_r8
  end subroutine fill_inputs

  !---------------------------------------------------------------------
  subroutine report(elapsed_s, per_call_s, checksum)
    real(r8), intent(in) :: elapsed_s, per_call_s, checksum

    write(*,'(a)') '=== p2c_2d_filter standalone driver ==='
    write(*,'(a,i0)')      '  ncol            = ', ncol
    write(*,'(a,i0)')      '  npatch_per_col  = ', npatch_per_col
    write(*,'(a,i0)')      '  nlev            = ', nlev
    write(*,'(a,i0)')      '  numfc           = ', numfc
    write(*,'(a,i0)')      '  niters          = ', niters
    write(*,'(a,i0)')      '  npatch_total    = ', npatch_total
    write(*,'(a,es14.6)')  '  elapsed (s)     = ', elapsed_s
    write(*,'(a,es14.6)')  '  per call (s)    = ', per_call_s
    write(*,'(a,es24.16)') '  checksum        = ', checksum
  end subroutine report

  !---------------------------------------------------------------------
  subroutine write_last_run(checksum)
    real(r8), intent(in) :: checksum
    integer :: u
    open(newunit=u, file='last_run.txt', status='replace', action='write')
    write(u,'(a,i0)')      'ncol            ', ncol
    write(u,'(a,i0)')      'npatch_per_col  ', npatch_per_col
    write(u,'(a,i0)')      'nlev            ', nlev
    write(u,'(a,i0)')      'numfc           ', numfc
    write(u,'(a,i0)')      'niters          ', niters
    write(u,'(a,es24.16)') 'checksum        ', checksum
    close(u)
  end subroutine write_last_run

  !---------------------------------------------------------------------
  subroutine compare_to_baseline(checksum)
    real(r8), intent(in) :: checksum
    integer  :: u, ios
    logical  :: exists
    character(len=64) :: key
    integer  :: b_ncol, b_npatch_per_col, b_nlev, b_numfc, b_niters
    real(r8) :: b_checksum, tol, diff

    inquire(file='baseline_checksum.txt', exist=exists)
    if (.not. exists) then
       write(*,'(a)') '  baseline        = (no baseline_checksum.txt found; skipping compare)'
       return
    end if

    open(newunit=u, file='baseline_checksum.txt', status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(*,'(a)') '  baseline        = (could not open baseline_checksum.txt)'
       return
    end if
    read(u,*,iostat=ios) key, b_ncol;           if (ios /= 0) goto 99
    read(u,*,iostat=ios) key, b_npatch_per_col; if (ios /= 0) goto 99
    read(u,*,iostat=ios) key, b_nlev;           if (ios /= 0) goto 99
    read(u,*,iostat=ios) key, b_numfc;          if (ios /= 0) goto 99
    read(u,*,iostat=ios) key, b_niters;         if (ios /= 0) goto 99
    read(u,*,iostat=ios) key, b_checksum;       if (ios /= 0) goto 99
    close(u)

    if (b_ncol /= ncol .or. b_npatch_per_col /= npatch_per_col .or. &
        b_nlev /= nlev .or. b_numfc /= numfc) then
       write(*,'(a)') '  baseline        = (param set differs; skipping compare)'
       return
    end if

    tol  = 1.0e-10_r8 * max(abs(b_checksum), 1.0_r8)
    diff = abs(checksum - b_checksum)
    if (diff <= tol) then
       write(*,'(a,es14.6,a)') '  baseline        = MATCH (|diff| = ', diff, ')'
    else
       write(*,'(a,es24.16)')  '  baseline value  = ', b_checksum
       write(*,'(a,es14.6,a,es14.6,a)') &
           '  baseline        = MISMATCH (|diff| = ', diff, ', tol = ', tol, ')'
    end if
    return

 99 continue
    close(u)
    write(*,'(a)') '  baseline        = (parse error in baseline_checksum.txt; skipping compare)'
  end subroutine compare_to_baseline

end program p2c_2d_filter_driver
