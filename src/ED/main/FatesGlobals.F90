module FatesGlobals
  ! NOTE(bja, 201608) This is a temporary hack module to store global
  ! data used inside fates. It's main use it to explicitly call out
  ! global data that needs to be dealt with, but doesn't have an
  ! immediately obvious home.

  implicit none

  integer, private :: fates_log_
  logical, private :: fates_global_verbose_

  public :: FatesGlobalsInit
  public :: fates_log
  public :: fates_global_verbose

contains

  subroutine FatesGlobalsInit(log_unit, global_verbose)

    implicit none

    integer, intent(in) :: log_unit
    logical, intent(in) :: global_verbose

    fates_log_ = log_unit
    fates_global_verbose_ = global_verbose

  end subroutine FatesGlobalsInit

  integer function fates_log()
    fates_log = fates_log_
  end function fates_log

  logical function fates_global_verbose()
    fates_global_verbose = fates_global_verbose_
  end function fates_global_verbose

end module FatesGlobals
