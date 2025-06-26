module SelfTestDriver

  ! ------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module contains a top-level driver to the self-test code.
  !
  ! See the README file in this directory for a high-level overview of these self-tests.

  implicit none
  private
  save

  ! Public routines

  public :: self_test_driver

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine self_test_driver(bounds)
    !
    ! !DESCRIPTION:
    ! Top-level driver to the self-test code
    !
    ! This subroutine should be called all the time, but each set of self tests is only
    ! run if the appropriate flag is set.
    !
    ! !USES:
    use clm_varctl, only : for_testing_run_ncdiopio_tests, for_testing_run_decomp_init_tests
    use clm_varctl, only : for_testing_exit_after_self_tests, iulog
    use decompMod, only : bounds_type
    use TestNcdioPio, only : test_ncdio_pio
    use ESMF, only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_Finalize
    use shr_sys_mod, only : shr_sys_flush
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'self_test_driver'
    integer :: ntests = 0
    !-----------------------------------------------------------------------

    if (for_testing_run_ncdiopio_tests) then
       ntests = ntests + 1
       call test_ncdio_pio(bounds)
    end if
    if (for_testing_run_decomp_init_tests) then
       ntests = ntests + 1
    end if
    if (for_testing_exit_after_self_tests) then
      if ( ntests == 0 )then
          write(iulog,*) 'WARNING: You are exiting after self tests were run -- but no self tests were run.'
        else
          write(iulog,*) 'Exiting after running ', ntests, ' self tests.'
        end if
        call shr_sys_flush(iulog)
        call ESMF_LogWrite(' exiting after running self tests', ESMF_LOGMSG_INFO)
        call ESMF_Finalize()
    end if

  end subroutine self_test_driver

end module SelfTestDriver
