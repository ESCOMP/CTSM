module initInterpMultilevelCopy

  ! ------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module defines a class for handling multi-level fields by simply copying the
  ! source to the destination, assuming the same number of levels in each.
  !
  ! !USES:
#include "shr_assert.h" 

  use shr_kind_mod             , only : r8 => shr_kind_r8
  use initInterpMultilevelBase , only : interp_multilevel_type

  implicit none
  private
  save

  ! Public types

  public :: interp_multilevel_copy_type

  type, extends(interp_multilevel_type) :: interp_multilevel_copy_type
     ! COMPILER_BUG(wjs, 2015-10-20, intel15.0.1) intel has problems creating a pointer to
     ! a class without any data components. Thus, including this unused dummy_var to make
     ! intel happy.
     integer :: dummy_var
   contains
     procedure :: check_npts
     procedure :: interp_multilevel
     procedure :: get_description
  end type interp_multilevel_copy_type

  interface interp_multilevel_copy_type
     module procedure constructor
  end interface interp_multilevel_copy_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  ! ========================================================================
  ! Constructors
  ! ========================================================================

  !-----------------------------------------------------------------------
  type(interp_multilevel_copy_type) function constructor()
    !
    ! !DESCRIPTION:
    ! Creates a new interp_multilevel_copy_type object
    !-----------------------------------------------------------------------

    ! Nothing needs to be done

  end function constructor

  ! ========================================================================
  ! Public methods
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine check_npts(this, npts, varname)
    !
    ! !DESCRIPTION:
    ! Checks the number of destination points, to ensure that this interpolator is
    ! appropriate for this variable. This should be called once for each variable.
    !
    ! This version accepts any number of points, because it has no point-based metadata.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(interp_multilevel_copy_type), intent(in) :: this
    integer, intent(in) :: npts             ! number of dest points (on this processor)
    character(len=*), intent(in) :: varname ! variable name (for diagnostic output)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'check_npts'
    !-----------------------------------------------------------------------

    return

  end subroutine check_npts


  !-----------------------------------------------------------------------
  subroutine interp_multilevel(this, data_dest, data_source, index_dest, &
                               scale_by_thickness)
    !
    ! !DESCRIPTION:
    ! Interpolates a multi-level field from source to dest, for a single point.
    !
    ! This version requires that data_dest and data_source be the same size, and it simply
    ! copies source to dest.
    !
    ! !ARGUMENTS:
    class(interp_multilevel_copy_type), intent(in) :: this
    real(r8) , intent(inout) :: data_dest(:)
    real(r8) , intent(in)    :: data_source(:)
    integer  , intent(in)    :: index_dest
    logical  , intent(in)    :: scale_by_thickness
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'interp_multilevel'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((size(data_source) == size(data_dest)), sourcefile, __LINE__)
    SHR_ASSERT_FL((.not. scale_by_thickness), sourcefile, __LINE__)

    ! Note that it's safe to do whole-array assignment here because we never decompose
    ! along the level dimension (in contrast to the spatial dimension, where you need to
    ! specify explicit bounds).
    data_dest(:) = data_source(:)

  end subroutine interp_multilevel

  !-----------------------------------------------------------------------
  pure function get_description(this) result(description)
    !
    ! !DESCRIPTION:
    ! Returns a text description of this interpolator
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    character(len=:), allocatable :: description  ! function result
    class(interp_multilevel_copy_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_description'
    !-----------------------------------------------------------------------
    
    description = 'Copy levels'

  end function get_description

end module initInterpMultilevelCopy

