module unittestArrayMod

  ! Provides utility functions for working with array inputs to (or outputs from)
  ! subroutines.
  !
  ! The routines here assume that the subgrid structure has been set up via
  ! unittestSubgridMod.

  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
  use unittestSubgridMod, only : bounds

  implicit none
  private
  save

  public :: grc_array_logical  ! create a gridcell-level array of logical values
  public :: col_array          ! create a column-level array of r8 values

contains

  !-----------------------------------------------------------------------
  pure function grc_array_logical(val) result(grc_array)
    !
    ! !DESCRIPTION:
    ! Creates a gridcell-level array of logical values. All elements are set to val.
    !
    ! If val is provided, all elements of the array are set to val; if not, all elements
    ! are set to false.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    logical, allocatable :: grc_array(:) ! function result
    logical, intent(in)  :: val ! all elements in grc_array are set to val
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'grc_array_logical'
    !-----------------------------------------------------------------------

    allocate(grc_array(bounds%begg:bounds%endg))
    grc_array(:) = val

  end function grc_array_logical


  !-----------------------------------------------------------------------
  pure function col_array(val)
    !
    ! !DESCRIPTION:
    ! Creates a column-level array of r8 values.
    !
    ! If val is provided, all elements of the array are set to val; if not, all elements
    ! are set to NaN.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), allocatable :: col_array(:) ! function result
    real(r8), intent(in), optional :: val ! if provided, all elements in col_array are set to val
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'col_array'
    !-----------------------------------------------------------------------

    allocate(col_array(bounds%begc:bounds%endc))
    if (present(val)) then
       col_array(:) = val
    else
       col_array(:) = nan
    end if

  end function col_array


end module unittestArrayMod
