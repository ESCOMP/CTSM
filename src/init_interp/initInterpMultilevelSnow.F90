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

     ! Number of active snow layers on the source grid, regridded to the destination grid
     !
     ! Thus, num_snow_layers_source(i) gives the number of active snow layers on the source
     ! grid for the source grid point that maps to destination point i.
     integer, allocatable :: num_snow_layers_source(:)
   contains
     procedure :: check_npts
     procedure :: interp_multilevel
     procedure :: get_description
  end type interp_multilevel_snow_type

  interface interp_multilevel_snow_type
     module procedure constructor
  end interface interp_multilevel_snow_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

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

    ! Number of active snow layers on the source grid, regridded to the destination grid
    !
    ! Thus, num_snow_layers_source(i) gives the number of active snow layers on the source
    ! grid for the source grid point that maps to destination point i.
    integer, intent(in) :: num_snow_layers_source(:)

    ! name of variable giving number of snow layers (just used for identification purposes)
    character(len=*), intent(in) :: num_layers_name
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'constructor'
    !-----------------------------------------------------------------------

    ! We do not check validity of num_snow_layers_source here (i.e., confirming that it
    ! is >= 0 everywhere) because it's okay for it to be invalid for points that are
    ! never invoked. This is the case for destination points with no corresponding source
    ! point (e.g., inactive destination points).

    allocate(constructor%num_snow_layers_source(size(num_snow_layers_source)))
    constructor%num_snow_layers_source(:) = num_snow_layers_source(:)

    constructor%num_layers_name = trim(num_layers_name)

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
    ! !USES:
    !
    ! !ARGUMENTS:
    class(interp_multilevel_snow_type), intent(in) :: this
    integer, intent(in) :: npts             ! number of dest points (on this processor)
    character(len=*), intent(in) :: varname ! variable name (for diagnostic output)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'check_npts'
    !-----------------------------------------------------------------------

    if (npts /= size(this%num_snow_layers_source)) then
       write(iulog,*) subname//' ERROR: mismatch in number of dest points for ', &
            trim(varname)
       write(iulog,*) 'Number of dest points: ', npts
       write(iulog,*) 'Expected number of dest points: ', &
            size(this%num_snow_layers_source)
       call endrun(msg=subname//' ERROR: mismatch in number of points for '//&
            trim(varname) // ' ' // errMsg(sourcefile, __LINE__))
    end if

  end subroutine check_npts

  !-----------------------------------------------------------------------
  subroutine interp_multilevel(this, data_dest, data_source, index_dest)
    !
    ! !DESCRIPTION:
    ! Interpolates a multi-level field from source to dest, for a single point.
    !
    ! This version does an "interpolation" (really a copy with offsets) appropriate for
    ! snow variables.
    !
    ! For most snow variables, we only need to copy data from EXISTING snow layers.
    ! However, there are a few variables (specifically, absorbed radiation fluxes) where
    ! we need to copy data even from non-existing snow layers in order to get bit-for-bit
    ! behavior upon interpolation. Thus, to be safe, we copy as many levels as possible.
    ! The algorithm is as follows:
    !
    ! - If (number of destination levels) >= (number of source levels), then copy all
    !   source levels to the destination. If there are more destination levels than
    !   source levels, then the lowest-index destination levels will be set to 0.
    !
    ! - If (number of destination levels) < (number of source levels), but (number of
    !   destination levels) >= (number of EXISTING snow levels in source), then copy all
    !   existing snow levels to the destination, plus as many non-existing levels as will
    !   fit.
    !
    ! - If (number of destination levels) < (number of EXISTING snow levels in source),
    !   then copy the top N existing snow levels to the destination, where N is the
    !   number of destination levels.
    !
    ! index_dest is used in this version, in order to match each point with its number
    ! of existing snow layers; index_dest should be 1-based.
    !
    ! !ARGUMENTS:
    class(interp_multilevel_snow_type), intent(in) :: this
    real(r8) , intent(inout) :: data_dest(:)
    real(r8) , intent(in)    :: data_source(:)
    integer  , intent(in)    :: index_dest
    !
    ! !LOCAL VARIABLES:
    integer :: num_source         ! total number of source layers
    integer :: num_dest           ! total number of dest layers
    integer :: num_snow_layers_source  ! number of EXISTING snow layers at this source point
    integer :: top_layer_source   ! top source layer copied to dest
    integer :: top_layer_dest     ! top dest layer receiving data from source
    integer :: source_layer       ! current source layer in copy
    integer :: dest_layer         ! current dest layer in copy

    character(len=*), parameter :: subname = 'interp_multilevel'
    !-----------------------------------------------------------------------

    num_source = size(data_source)
    num_dest = size(data_dest)
    num_snow_layers_source = this%num_snow_layers_source(index_dest)
    SHR_ASSERT_FL(num_snow_layers_source >= 0, sourcefile, __LINE__)
    SHR_ASSERT_FL(num_snow_layers_source <= num_source, sourcefile, __LINE__)

    if (num_dest >= num_source) then
       ! Copy all source layers to dest (even non-existent snow layers)
       top_layer_source = 1
       top_layer_dest = 1 + (num_dest - num_source)
    else if (num_dest >= num_snow_layers_source) then
       ! Copy all existing source layers to dest, plus as many non-existing layers as will fit
       top_layer_dest = 1
       top_layer_source = 1 + (num_source - num_dest)
    else
       ! num_dest < num_snow_layers_source. Copy as many existing layers from source as
       ! will fit in dest, starting at the top of the snow pack.
       top_layer_dest = 1
       top_layer_source = num_source - num_snow_layers_source + 1
    end if

    do dest_layer = 1, (top_layer_dest - 1)
       data_dest(dest_layer) = 0._r8
    end do

    source_layer = top_layer_source
    do dest_layer = top_layer_dest, num_dest
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
