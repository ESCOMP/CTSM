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
  subroutine self_test_driver
    !
    ! !DESCRIPTION:
    ! Top-level driver to the self-test code
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'self_test_driver'
    !-----------------------------------------------------------------------

  end subroutine self_test_driver

end module SelfTestDriver
