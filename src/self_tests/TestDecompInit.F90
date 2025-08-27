module TestDecompInit

  ! ------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module contains tests of decomp_init

#include "shr_assert.h"
  use shr_kind_mod, only : r8 => shr_kind_r8, CX => shr_kind_cx
  use Assertions, only : assert_equal
  use clm_varctl, only : iulog
  use abortutils, only : endrun, endrun_init, get_last_endrun_msg
  use spmdMod, only : masterproc, npes
  use decompInitMod, only : decompInit_lnd, clump_pproc
  use decompMod

  implicit none
  private
  save

  ! Public routines

  public :: test_decomp_init

  ! Module data used in various tests

  integer, parameter :: ni = 10, nj = 5
  integer :: amask(ni*nj)

  integer :: default_npes
  integer :: default_clump_pproc

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine test_decomp_init()
    !
    ! !DESCRIPTION:
    ! Drive tests of decomp_init
    !
    ! NOTE(wjs, 2020-10-15) Currently, endrun is called when any test assertion fails. I
    ! thought about changing this so that, instead, a counter is incremented for each
    ! failure, then at the end of the testing (in the higher-level self-test driver),
    ! endrun is called if this counter is greater than 0. The benefit of this is that we'd
    ! see all test failures, not just the first failure. To do that, we'd need to change
    ! the assertions here to increment a counter rather than aborting. However, I'm not
    ! spending the time to make this change for now because (1) I'm not sure how much
    ! value we'd get from it; (2) even if we made that change, it's still very possible
    ! for test code to abort for reasons other than assertions, if something goes wrong
    ! inside decomp_init or pio; and (3) some tests here are dependent on earlier tests (for
    ! example, the reads depend on the writes having worked), so a failure in an early
    ! phase could really muck things up for later testing phases. Migrating to a
    ! pFUnit-based unit test would solve this problem, since each pFUnit test is
    ! independent, though would prevent us from being able to have dependent tests the
    ! way we do here (where reads depend on earlier writes), for better or for worse.
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

    default_npes = npes
    default_clump_pproc = clump_pproc
    call write_to_log('start_test_decomp_init')

    call write_to_log('test_check_nclumps')
    call test_check_nclumps()
    call write_to_log('test_decompInit_lnd_abort_on_bad_clump_pproc')
    call test_decompInit_lnd_abort_on_bad_clump_pproc()
    call write_to_log('test_decompInit_lnd_abort_on_too_big_clump_pproc')
    call test_decompInit_lnd_abort_on_too_big_clump_pproc()
    call write_to_log('test_decompInit_lnd_abort_when_npes_too_large')
    call test_decompInit_lnd_abort_when_npes_too_large()
    call write_to_log('test_decompInit_lnd_abort_on_too_small_nsegspc')
    call test_decompInit_lnd_abort_on_too_small_nsegspc()
    call write_to_log('test_decompInit_lnd_check_sizes')
    call test_decompInit_lnd_check_sizes()

  end subroutine test_decomp_init

  !-----------------------------------------------------------------------
  subroutine setup()
     use clm_varctl, only : nsegspc

     clump_pproc = default_clump_pproc
     nsegspc = 20
     npes = default_npes
     amask(:) = 1 ! Set all to land

  end subroutine setup

  !-----------------------------------------------------------------------
  subroutine test_decompInit_lnd_abort_on_bad_clump_pproc()
     character(len=CX) :: expected_msg, actual_msg

     call setup()
     call endrun_init( .true. )  ! Do not abort on endrun for self-tests
     clump_pproc = 0
     call write_to_log('decompInit_lnd with clump_pproc=0 should abort')
     call decompInit_lnd( ni, nj, amask )
     call write_to_log('check expected abort message')
     expected_msg = 'clump_pproc must be greater than 0'
     actual_msg = get_last_endrun_msg()
     call endrun_init( .false. )   ! Turn back on to abort on the assert
     call write_to_log('call assert_equal to check the abort message')
     call assert_equal( &
           expected=expected_msg, actual=actual_msg, &
           msg='decompInit_lnd did not abort with clump_pproc=0' )
     call clean()
  end subroutine test_decompInit_lnd_abort_on_bad_clump_pproc

  !-----------------------------------------------------------------------
  subroutine test_decompInit_lnd_abort_on_too_big_clump_pproc()
     character(len=CX) :: expected_msg, actual_msg

     call setup()
     call endrun_init( .true. )  ! Do not abort on endrun for self-tests
     amask(:) = 1 ! Set all to land
     clump_pproc = (ni * nj + 1) / npes
     call write_to_log('decompInit_lnd with clump_pproc too large should abort')
     call decompInit_lnd( ni, nj, amask )
     call write_to_log('check expected abort message')
     expected_msg = 'Number of clumps exceeds number of land grid cells'
     actual_msg = get_last_endrun_msg()
     call endrun_init( .false. )   ! Turn back on to abort on the assert
     call write_to_log('call assert_equal to check the abort message')
     call assert_equal( &
           expected=expected_msg, actual=actual_msg, &
           msg='decompInit_lnd did not abort with clump_pproc too large' )
     call clean()
  end subroutine test_decompInit_lnd_abort_on_too_big_clump_pproc

  !-----------------------------------------------------------------------
  subroutine test_decompInit_lnd_check_sizes()
     use decompMod, only : get_proc_bounds
     type(bounds_type) :: bounds

     integer :: expected_endg, expected_numg

     call setup()
     expected_numg = ni*nj
     if ( expected_numg < npes )then
        call endrun( msg="npes is too large for this test", file=sourcefile, line=__LINE__ )
     end if
     if ( modulo( expected_numg, npes ) /= 0 )then
        call endrun( msg="npes does not evenly divide into numg so this test will not work", file=sourcefile, line=__LINE__ )
     end if
     expected_endg = ni*nj / npes
     amask(:) = 1 ! Set all to land
     call decompInit_lnd( ni, nj, amask )
     call get_proc_bounds(bounds, allow_errors=.true.)
     call assert_equal( bounds%begg, 1, msg='begg is not as expected' )
     call assert_equal( bounds%endg, expected_endg, msg='endg is not as expected' )
     call clean()
  end subroutine test_decompInit_lnd_check_sizes

  !-----------------------------------------------------------------------
  subroutine test_decompInit_lnd_abort_when_npes_too_large()
     character(len=CX) :: expected_msg, actual_msg

     call setup()
     ! NOTE: This is arbitrarily modifying the NPES value -- so it MUST be reset set the END!
     npes = ni*nj + 1

     call endrun_init( .true. )  ! Do not abort on endrun for self-tests
     amask(:) = 1 ! Set all to land
     call write_to_log('decompInit_lnd with npes too large should abort')
     call decompInit_lnd( ni, nj, amask )
     call write_to_log('check expected abort message')
     expected_msg = 'Number of processes exceeds number of land grid cells'
     actual_msg = get_last_endrun_msg()
     call endrun_init( .false. )   ! Turn back on to abort on the assert
     call write_to_log('call assert_equal to check the abort message')
     call assert_equal( &
           expected=expected_msg, actual=actual_msg, &
           msg='decompInit_lnd did not abort with npes too large' )

     ! NOTE: Return npes to its original value
     npes = default_npes
     call clean()
  end subroutine test_decompInit_lnd_abort_when_npes_too_large

  !-----------------------------------------------------------------------
  subroutine test_decompInit_lnd_abort_on_too_small_nsegspc()
     use clm_varctl, only : nsegspc
     character(len=CX) :: expected_msg, actual_msg

     call setup()
     call endrun_init( .true. )  ! Do not abort on endrun for self-tests
     amask(:) = 1 ! Set all to land
     nsegspc = 0
     call write_to_log('decompInit_lnd with nsegspc too small should abort')
     call decompInit_lnd( ni, nj, amask )
     call write_to_log('check expected abort message')
     expected_msg = 'Number of segments per clump (nsegspc) is less than 1 and can NOT be'
     actual_msg = get_last_endrun_msg()
     call endrun_init( .false. )   ! Turn back on to abort on the assert
     call write_to_log('call assert_equal to check the abort message')
     call assert_equal( &
           expected=expected_msg, actual=actual_msg, &
           msg='decompInit_lnd did not abort with too nsegspc too small' )
     call clean()
  end subroutine test_decompInit_lnd_abort_on_too_small_nsegspc

  !-----------------------------------------------------------------------
  subroutine test_check_nclumps()
    integer :: expected_nclumps

    call setup()
    call endrun_init( .true. )  ! Do not abort on endrun for self-tests
    expected_nclumps = npes / clump_pproc
    call assert_equal(expected=expected_nclumps, actual=nclumps, &
         msg='nclumps are not as expected')
    call endrun_init( .false. )
    call clean()
  end subroutine test_check_nclumps

  !-----------------------------------------------------------------------
  subroutine write_to_log(msg)
    !
    ! !DESCRIPTION:
    ! Write a message to the log file, just from the masterproc
    !
    use shr_sys_mod, only : shr_sys_flush
    ! !ARGUMENTS:
    character(len=*), intent(in) :: msg
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'write_to_log'
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,'(a)') msg
       call shr_sys_flush(iulog)  ! Flush the I/O buffers always
    end if

  end subroutine write_to_log

  !-----------------------------------------------------------------------
  subroutine clean
    !
    ! !DESCRIPTION:
    ! Do end-of-testing cleanup after each test
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------
    call decompmod_clean()

  end subroutine clean


end module TestDecompInit
