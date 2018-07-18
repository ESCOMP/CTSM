module initInterpMultilevelBase

  ! ------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module defines a base class for handling multi-level fields
  !
  ! Usage: For a given variable:
  !
  !   First call check_npts with the number of destination points, to ensure that this
  !   interpolator is appropriate for this variable.
  !
  !   Then call interp_multilevel once for each destination point
  !
  ! Note that these interpolators only care about the destination point, not the source
  ! point. Any information on the source grid is expected to have already been
  ! interpolated to the destination grid. (This is needed for the sake of memory
  ! scalability: Information on the source grid is not decomposed across processors, so
  ! results in large amounts of memory usage.)
  !
  ! !USES:

  implicit none
  private
  save

  ! Public types

  public :: interp_multilevel_type

  type, abstract :: interp_multilevel_type
   contains
     procedure(check_npts_interface), deferred        :: check_npts
     procedure(interp_multilevel_interface), deferred :: interp_multilevel
     procedure(get_description_interface), deferred   :: get_description  ! get text description of interpolator
  end type interp_multilevel_type

  abstract interface

     subroutine check_npts_interface(this, npts, varname)
       ! !DESCRIPTION:
       ! Checks the number of destination points, to ensure that this interpolator is
       ! appropriate for this variable. This should be called once for each variable.
       !
       ! Aborts if there is a mismatch.
       !
       ! !USES:
       import:: interp_multilevel_type
       !
       ! !ARGUMENTS:
       class(interp_multilevel_type), intent(in) :: this
       integer, intent(in) :: npts             ! number of dest points (on this processor)
       character(len=*), intent(in) :: varname ! variable name (for diagnostic output)
     end subroutine check_npts_interface

     subroutine interp_multilevel_interface(this, &
          data_dest, data_source, index_dest)
       ! !DESCRIPTION:
       ! Interpolates a multi-level field from source to dest, for a single point.
       !
       ! data_dest and data_source give values for all levels, for one destination or
       ! source point. index_dest gives the spatial index (e.g., column index) of the dest
       ! point; this is needed for some types of interpolation (to find the appropriate
       ! metadata), but is ignored by others.  This index should be 1-based (i.e., if the
       ! lower bounds in the caller are not 1, they should be adjusted so that the first
       ! index is 1).
       !
       ! !USES:
       use shr_kind_mod, only : r8 => shr_kind_r8
       import :: interp_multilevel_type
       !
       ! !ARGUMENTS:
       class(interp_multilevel_type), intent(in) :: this
       real(r8) , intent(inout) :: data_dest(:)
       real(r8) , intent(in)    :: data_source(:)
       integer  , intent(in)    :: index_dest
     end subroutine interp_multilevel_interface

     pure function get_description_interface(this) result(description)
       ! !DESCRIPTION
       ! Returns a text description of this interpolator
       !
       ! !USES:
       import :: interp_multilevel_type
       !
       ! !ARGUMENTS:
       character(len=:), allocatable :: description  ! function result
       class(interp_multilevel_type), intent(in) :: this
     end function get_description_interface
  end interface

end module initInterpMultilevelBase
