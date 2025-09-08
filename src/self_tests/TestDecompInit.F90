module TestDecompInit

  ! ------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module contains tests of decomp_init

#include "shr_assert.h"
  use shr_kind_mod, only : r8 => shr_kind_r8, CX => shr_kind_cx
  use Assertions, only : assert_equal
  use clm_varctl, only : iulog
  use abortutils, only : endrun, endrun_init, get_last_endrun_msg
  use spmdMod, only : masterproc, npes, iam
  use decompInitMod, only : decompInit_lnd, clump_pproc, decompInit_clumps
  use clm_InstMod, only : glc_behavior
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
    ! !USERS:
    use decompInitMod, only : decompInit_clumps, decompInit_glcp
    use domainMod, only : ldomain
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer, allocatable :: model_amask(:)
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
    call write_to_log('test_decompInit_clump_gcell_info_correct')
    call test_decompInit_clump_gcell_info_correct()
    ! Comment out for now -- needs some work
    !call write_to_log('test_decompMod_get_clump_bounds_correct')
    !call test_decompMod_get_clump_bounds_correct()

    !
    ! Call the decompInit initialization series a last time so that decompMod data can still be used
    !
    !allocate( model_amask(ldomain%ni*ldomain%nj) )
    !model_amask(:) = 1
    !call decompInit_lnd( ldomain%ni, ldomain%nj, model_amask )
    !call decompInit_clumps(ldomain%ni, ldomain%nj, glc_behavior)
    !call decompInit_glcp(ldomain%ni, ldomain%nj, glc_behavior)
    !deallocate( model_amask )

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

     if ( npes > 1 ) return   ! error checking testing only works seriallly
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

     if ( npes > 1 ) return   ! error checking testing only works seriallly
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

     if ( npes > 1 ) return   ! error checking testing only works seriallly
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

     if ( npes > 1 ) return   ! error checking testing only works seriallly
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
  subroutine test_decompMod_get_clump_bounds_correct()
    ! Some testing for get_clump_bounds
    use decompMod, only : get_clump_bounds, bounds_type
    use unittestSimpleSubgridSetupsMod, only :  setup_ncells_single_veg_patch
    use unittestSubgridMod, only : unittest_subgrid_teardown
    use pftconMod, only : noveg
    type(bounds_type) :: bounds
    integer :: expected_begg, expected_endg, expected_numg, gcell_per_task
    integer :: iclump

    call setup()
     ! Now setup a singple grid that's just the full test with every point a single baresoil patch
     call setup_ncells_single_veg_patch( ncells=ni*nj, pft_type=noveg )
    clump_pproc = 1  ! Ensure we are just doing this for one clump per proc for now
    expected_numg = ni*nj
    if ( expected_numg < npes )then
       call endrun( msg="npes is too large for this test", file=sourcefile, line=__LINE__ )
    end if
    if ( modulo( expected_numg, npes ) /= 0 )then
       call endrun( msg="npes does not evenly divide into numg so this test will not work", file=sourcefile, line=__LINE__ )
    end if
    gcell_per_task = expected_numg / npes
    expected_begg = gcell_per_task * iam + 1
    expected_endg = expected_begg + gcell_per_task
    amask(:) = 1 ! Set all to land
    call decompInit_lnd( ni, nj, amask )
    call decompInit_clumps( ni, nj, glc_behavior )
    iclump = 1 ! Clump is just 1 since there's only one clump per task
    call get_clump_bounds(iclump, bounds)
    call assert_equal( bounds%begg, expected_begg, msg='begg is not as expected' )
    call assert_equal( bounds%endg, expected_endg, msg='endg is not as expected' )
    ! Other subgrtid level information will be the same -- since there's only one landunit, column, and patch per gridcell
    call assert_equal( bounds%begl, expected_begg, msg='begl is not as expected' )
    call assert_equal( bounds%endl, expected_endg, msg='endl is not as expected' )
    call assert_equal( bounds%begc, expected_begg, msg='begc is not as expected' )
    call assert_equal( bounds%endc, expected_endg, msg='endc is not as expected' )
    call assert_equal( bounds%begp, expected_begg, msg='begp is not as expected' )
    call assert_equal( bounds%endp, expected_endg, msg='endp is not as expected' )
    call unittest_subgrid_teardown( )
    call clean()
  end subroutine test_decompMod_get_clump_bounds_correct

  !-----------------------------------------------------------------------
  subroutine test_decompInit_clump_gcell_info_correct()
    ! Some testing for get_clump_bounds
    use decompMod, only : clumps
    use decompMod, only : get_proc_bounds
    type(bounds_type) :: bounds
    integer :: expected_gcells, iclump, g, beg_global_index, gcell_per_task
    integer :: expected_begg, expected_endg

    call setup()
    expected_gcells = ni*nj
    if ( expected_gcells < npes )then
       call endrun( msg="npes is too large for this test", file=sourcefile, line=__LINE__ )
    end if
    if ( modulo( expected_gcells, npes ) /= 0 )then
       call endrun( msg="npes does not evenly divide into gcell so this test will not work", file=sourcefile, line=__LINE__ )
    end if
    gcell_per_task = expected_gcells / npes
    expected_begg = gcell_per_task * iam + 1
    expected_endg = expected_begg + gcell_per_task
    amask(:) = 1 ! Set all to land
    call decompInit_lnd( ni, nj, amask )
    ! When clump_pproc is one clumps will be the same as PE
    if ( clump_pproc == 1 ) then
       call assert_equal( nclumps, npes, msg='nclumps should match number of processors when clump_pproc is 1' )
    else
       call assert_equal( nclumps/clump_pproc, npes, msg='nclumps divided by clump_pproc should match number of processors when clump_pproc > 1' )
    end if
    do iclump = 1, nclumps
       call assert_equal( clumps(iclump)%owner, iclump-1, msg='clumps owner is not correct' )
       call assert_equal( clumps(iclump)%ncells, gcell_per_task, msg='clumps ncells is not correct' )
    end do
    call clean()
  end subroutine test_decompInit_clump_gcell_info_correct

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
