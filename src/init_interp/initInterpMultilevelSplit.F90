module initInterpMultilevelSplit

  ! ------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module defines a class for handling multi-level fields by doing two different
  ! interpolations: One for some first set of levels, and a different one for some second
  ! set of levels.
  !
  ! !USES:
#include "shr_assert.h"

  use shr_kind_mod             , only : r8 => shr_kind_r8
  use shr_log_mod              , only : errMsg => shr_log_errMsg
  use abortutils               , only : endrun
  use clm_varctl               , only : iulog
  use initInterpMultilevelBase , only : interp_multilevel_type

  implicit none
  private
  save

  ! Public types

  public :: interp_multilevel_split_type

  type, extends(interp_multilevel_type) :: interp_multilevel_split_type
     private
     class(interp_multilevel_type), pointer :: interpolator_first_levels  => null()
     class(interp_multilevel_type), pointer :: interpolator_second_levels => null()
     integer :: num_second_levels_source
     integer :: num_second_levels_dest
   contains
     ! Public methods
     procedure :: check_npts
     procedure :: interp_multilevel
     procedure :: get_description
  end type interp_multilevel_split_type

  ! Constructor
  ! NOTE(wjs, 2015-10-23) This is given a different name from
  ! interp_multilevel_split_type because some compilers (in particular intel 15.0.1) had
  ! trouble distinguishing between the user-defined constructor and the default structure
  ! constructor. I'm not sure if this is a compiler bug or not.
  public :: create_interp_multilevel_split_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  ! ========================================================================
  ! Constructors
  ! ========================================================================

  !-----------------------------------------------------------------------
  function create_interp_multilevel_split_type( &
       interpolator_first_levels, interpolator_second_levels, &
       num_second_levels_source, num_second_levels_dest) &
       result(this)
    !
    ! !DESCRIPTION:
    ! Construct an interp_multilevel_split_type object.
    !
    ! interpolator_first_levels gives the interpolator for the first set of levels (1:nsrc
    ! and 1:ndst); interpolator_second_levels gives the interpolator for the second set of
    ! levels (nsrc+1:msrc and ndst+1:mdst).
    !
    ! You must specify the number of levels in the *second* set of levels for the source
    ! and dest (i.e., the number of levels that are used by interpolator_second_levels).
    ! (In principle, the number of levels in the *first* set of levels can then vary from
    ! one call to interp_multilevel to another.)
    !
    ! NOTE(wjs, 2015-10-23) This is given a different name from
    ! interp_multilevel_split_type because some compilers (in particular intel 15.0.1) had
    ! trouble distinguishing between the user-defined constructor and the default structure
    ! constructor. I'm not sure if this is a compiler bug or not.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(interp_multilevel_split_type) :: this  ! function result
    class(interp_multilevel_type), target, intent(in) :: interpolator_first_levels
    class(interp_multilevel_type), target, intent(in) :: interpolator_second_levels
    integer, intent(in) :: num_second_levels_source
    integer, intent(in) :: num_second_levels_dest
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'create_interp_multilevel_split_type'
    !-----------------------------------------------------------------------

    if (num_second_levels_source <= 0) then
       write(iulog,*) "For interp_multilevel_split_type, num_second_levels_source must be > 0"
       write(iulog,*) "num_second_levels_source = ", num_second_levels_source
       call endrun(msg="num_second_levels_source must be > 0 "//errMsg(sourcefile, __LINE__))
    end if
    if (num_second_levels_dest <= 0) then
       write(iulog,*) "For interp_multilevel_split_type, num_second_levels_dest must be > 0"
       write(iulog,*) "num_second_levels_dest = ", num_second_levels_dest
       call endrun(msg="num_second_levels_dest must be > 0 "//errMsg(sourcefile, __LINE__))
    end if

    this%interpolator_first_levels => interpolator_first_levels
    this%interpolator_second_levels => interpolator_second_levels
    this%num_second_levels_source = num_second_levels_source
    this%num_second_levels_dest = num_second_levels_dest

  end function create_interp_multilevel_split_type

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
    ! !USES:
    !
    ! !ARGUMENTS:
    class(interp_multilevel_split_type), intent(in) :: this
    integer, intent(in) :: npts             ! number of dest points (on this processor)
    character(len=*), intent(in) :: varname ! variable name (for diagnostic output)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'check_npts'
    !-----------------------------------------------------------------------

    call this%interpolator_first_levels%check_npts(npts, varname)
    call this%interpolator_second_levels%check_npts(npts, varname)
  end subroutine check_npts

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
    class(interp_multilevel_split_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_description'
    !-----------------------------------------------------------------------

    description = 'Split levels: ' // &
         trim(this%interpolator_first_levels%get_description()) // &
         ' + ' // &
         trim(this%interpolator_second_levels%get_description())

  end function get_description

  !-----------------------------------------------------------------------
  subroutine interp_multilevel(this, data_dest, data_source, index_dest, scale_by_thickness)
    !
    ! !DESCRIPTION:
    ! Interpolates a multi-level field from source to dest, for a single point.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(interp_multilevel_split_type), intent(in) :: this
    real(r8) , intent(inout) :: data_dest(:)
    real(r8) , intent(in)    :: data_source(:)
    integer  , intent(in)    :: index_dest
    logical  , intent(in)    :: scale_by_thickness
    !
    ! !LOCAL VARIABLES:
    integer :: num_first_levels_dest
    integer :: num_first_levels_source

    character(len=*), parameter :: subname = 'interp_multilevel'
    !-----------------------------------------------------------------------

    num_first_levels_dest = size(data_dest) - this%num_second_levels_dest
    num_first_levels_source = size(data_source) - this%num_second_levels_source

    if (num_first_levels_source <= 0) then
       write(iulog,*) "For interp_multilevel_split_type, num_first_levels_source must be > 0"
       write(iulog,*) "num_first_levels_source = ", num_first_levels_source
       call endrun(msg="num_first_levels_source must be > 0 "//errMsg(sourcefile, __LINE__))
    end if
    if (num_first_levels_dest <= 0) then
       write(iulog,*) "For interp_multilevel_split_type, num_first_levels_dest must be > 0"
       write(iulog,*) "num_first_levels_dest = ", num_first_levels_dest
       call endrun(msg="num_first_levels_dest must be > 0 "//errMsg(sourcefile, __LINE__))
    end if

    call this%interpolator_first_levels%interp_multilevel( &
         data_dest = data_dest(1:num_first_levels_dest), &
         data_source = data_source(1:num_first_levels_source), &
         index_dest = index_dest, &
         scale_by_thickness = scale_by_thickness)

    call this%interpolator_second_levels%interp_multilevel( &
         data_dest = data_dest((num_first_levels_dest+1):size(data_dest)), &
         data_source = data_source((num_first_levels_source+1):size(data_source)), &
         index_dest = index_dest, &
         scale_by_thickness = scale_by_thickness)

  end subroutine interp_multilevel


end module initInterpMultilevelSplit
