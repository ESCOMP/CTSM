module SelfTestDriver

  ! ------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module contains a top-level driver to the self-test code.
  !
  ! See the README file in this directory for a high-level overview of these self-tests.

  use decompMod, only : bounds_type
  use TestNcdioPio, only : test_ncdio_pio

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
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'self_test_driver'
    !-----------------------------------------------------------------------

    call test_ncdio_pio(bounds)

  end subroutine self_test_driver

end module SelfTestDriver
