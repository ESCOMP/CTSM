module initInterpMultilevelSnow

  ! ------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module defines a class for handling multi-level snow fields.
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

  public :: interp_multilevel_snow_type

  type, extends(interp_multilevel_type) :: interp_multilevel_snow_type
     private
     character(len=:), allocatable :: num_layers_name
     integer, allocatable :: num_snow_layers_source(:)  ! number of active snow layers for each point on the source grid
     integer :: npts_source  ! number of spatial points on source grid
   contains
     procedure :: check_npts
     procedure :: interp_multilevel
     procedure :: get_description
  end type interp_multilevel_snow_type

  interface interp_multilevel_snow_type
     module procedure constructor
  end interface interp_multilevel_snow_type

contains

  ! ========================================================================
  ! Constructors
  ! ========================================================================

  !-----------------------------------------------------------------------
  type(interp_multilevel_snow_type) function constructor(num_snow_layers_source, &
       num_layers_name)
    !
    ! !DESCRIPTION:
    ! Creates a new interp_multilevel_snow_type object
    !
    ! !USES:
    !
    ! !ARGUMENTS:

    ! number of active snow layers for each point on the source grid
    integer, intent(in) :: num_snow_layers_source(:)

    ! name of variable giving number of snow layers (just used for identification purposes)
    character(len=*), intent(in) :: num_layers_name
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'constructor'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL(num_snow_layers_source >= 0, errMsg(__FILE__, __LINE__))

    constructor%npts_source = size(num_snow_layers_source)

    allocate(constructor%num_snow_layers_source(constructor%npts_source))
    constructor%num_snow_layers_source(:) = num_snow_layers_source(:)

    constructor%num_layers_name = trim(num_layers_name)

  end function constructor

  ! ========================================================================
  ! Public methods
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine check_npts(this, npts_source, npts_dest, varname)
    !
    ! !DESCRIPTION:
    ! Checks the number of source and destination points, to ensure that this
    ! interpolator is appropriate for this variable. This should be called once for
    ! each variable.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(interp_multilevel_snow_type), intent(in) :: this
    integer, intent(in) :: npts_source      ! number of source points
    integer, intent(in) :: npts_dest        ! number of dest points (on this processor)
    character(len=*), intent(in) :: varname ! variable name (for diagnostic output)
    !
    ! !LOCAL VARIABLES:
    logical :: mismatch

    character(len=*), parameter :: subname = 'check_npts'
    !-----------------------------------------------------------------------

    mismatch = .false.

    if (npts_source /= this%npts_source) then
       write(iulog,*) subname//' ERROR: mismatch in number of source points for ', &
            trim(varname)
       write(iulog,*) 'Number of source points: ', npts_source
       write(iulog,*) 'Expected number of source points: ', this%npts_source
       mismatch = .true.
    end if

    ! Note: this class doesn't care about the number of destination points

    if (mismatch) then
       call endrun(msg=subname//' ERROR: mismatch in number of points for '//&
            trim(varname) // ' ' // errMsg(__FILE__, __LINE__))
    end if

  end subroutine check_npts

  !-----------------------------------------------------------------------
  subroutine interp_multilevel(this, data_dest, data_source, index_dest, index_source)
    !
    ! !DESCRIPTION:
    ! Interpolates a multi-level field from source to dest, for a single point.
    !
    ! This version does an "interpolation" (really a copy with offsets) appropriate for
    ! snow variables. This is based on the number of EXISTING snow layers in each source
    ! point. If the destination has at least as many levels as the number of existing snow
    ! layers in the source, then all existing snow layers are copied to the destination,
    ! with non-existing levels set to 0. If the destination has fewer levels than the
    ! number of existing snow layers, then the top N levels of the source are copied to
    ! the destination (where N is the number of levels in the destination).
    !
    ! index_source is used in this version, in order to match each point with its number
    ! of existing snow layers; index_source should be 1-based.
    !
    ! !ARGUMENTS:
    class(interp_multilevel_snow_type), intent(in) :: this
    real(r8) , intent(inout) :: data_dest(:)
    real(r8) , intent(in)    :: data_source(:)
    integer  , intent(in)    :: index_dest
    integer  , intent(in)    :: index_source
    !
    ! !LOCAL VARIABLES:
    integer :: num_snow_layers_source  ! number of existing snow layers at this source point
    integer :: top_snow_layer_source
    integer :: top_snow_layer_dest
    integer :: dest_layer
    integer :: source_layer

    character(len=*), parameter :: subname = 'interp_multilevel'
    !-----------------------------------------------------------------------

    num_snow_layers_source = this%num_snow_layers_source(index_source)
    SHR_ASSERT(num_snow_layers_source <= size(data_source), errMsg(__FILE__, __LINE__))

    ! Determine the index in the source for the top layer that has snow. If there is snow
    ! in every layer, then top_snow_layer_source will be 1. If there is no snow,
    ! top_snow_layer_source will be size(data_source) + 1.
    top_snow_layer_source = size(data_source) - num_snow_layers_source + 1

    ! If there are N snow layers in the source, there will be N snow layers in the dest...
    top_snow_layer_dest = top_snow_layer_source + (size(data_dest) - size(data_source))

    ! ... but there cannot be more snow layers than the size allows; if we cannot fit all
    ! snow layers from the source in the dest, we will copy the *top* snow layers from the
    ! source into the dest
    top_snow_layer_dest = max(top_snow_layer_dest, 1)

    do dest_layer = 1, (top_snow_layer_dest - 1)
       data_dest(dest_layer) = 0._r8
    end do

    source_layer = top_snow_layer_source
    do dest_layer = top_snow_layer_dest, size(data_dest)
       data_dest(dest_layer) = data_source(source_layer)
       source_layer = source_layer + 1
    end do

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
    class(interp_multilevel_snow_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_description'
    !-----------------------------------------------------------------------
    
    description = 'Copy snow-covered levels using '//this%num_layers_name

  end function get_description

end module initInterpMultilevelSnow
