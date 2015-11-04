module initInterpBounds

  ! ------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module defines a class for storing and querying bounds information for initInterp
  !
  ! !USES:
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use abortutils      , only : endrun

  implicit none
  private
  save

  ! Public types

  public :: interp_bounds_type

  type :: interp_bounds_type
     private
     integer :: begp  ! beginning patch-level index
     integer :: endp  ! ending patch-level index
     integer :: begc  ! beginning col-level index
     integer :: endc  ! ending col-level index
     integer :: begl  ! beginning landunit-level index
     integer :: endl  ! ending landunit-level index
   contains
     procedure :: get_begp
     procedure :: get_endp
     procedure :: get_begc
     procedure :: get_endc
     procedure :: get_begl
     procedure :: get_endl
     procedure :: get_beg  ! get beginning index for a given subgrid level
     procedure :: get_end  ! get ending index for a given subgrid level
  end type interp_bounds_type

  interface interp_bounds_type
     module procedure constructor
  end interface interp_bounds_type

contains

  ! ========================================================================
  ! Constructors
  ! ========================================================================

  !-----------------------------------------------------------------------
  function constructor(begp, endp, begc, endc, begl, endl) result(this)
    !
    ! !DESCRIPTION:
    ! Create an interp_bounds_type instance
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(interp_bounds_type) :: this  ! function result
    integer, intent(in) :: begp, endp
    integer, intent(in) :: begc, endc
    integer, intent(in) :: begl, endl
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'constructor'
    !-----------------------------------------------------------------------

    this%begp = begp
    this%endp = endp
    this%begc = begc
    this%endc = endc
    this%begl = begl
    this%endl = endl

  end function constructor

  ! ========================================================================
  ! Public methods
  ! ========================================================================

  integer function get_begp(this)
    class(interp_bounds_type), intent(in) :: this
    get_begp = this%begp
  end function get_begp

  integer function get_endp(this)
    class(interp_bounds_type), intent(in) :: this
    get_endp = this%endp
  end function get_endp

  integer function get_begc(this)
    class(interp_bounds_type), intent(in) :: this
    get_begc = this%begc
  end function get_begc

  integer function get_endc(this)
    class(interp_bounds_type), intent(in) :: this
    get_endc = this%endc
  end function get_endc

  integer function get_begl(this)
    class(interp_bounds_type), intent(in) :: this
    get_begl = this%begl
  end function get_begl

  integer function get_endl(this)
    class(interp_bounds_type), intent(in) :: this
    get_endl = this%endl
  end function get_endl

  !-----------------------------------------------------------------------
  function get_beg(this, subgrid_level) result(beg_index)
    !
    ! !DESCRIPTION:
    ! Get beginning index for a given subgrid level
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer :: beg_index  ! function result
    class(interp_bounds_type), intent(in) :: this
    character(len=*), intent(in) :: subgrid_level  ! 'pft', 'column' or 'landunit'
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_beg'
    !-----------------------------------------------------------------------

    select case (subgrid_level)
    case('pft')
       beg_index = this%begp
    case('column')
       beg_index = this%begc
    case('landunit')
       beg_index = this%begl
    case default
       call endrun(msg=subname//' ERROR: Unknown subgrid level: '//trim(subgrid_level)// &
            errMsg(__FILE__, __LINE__))
    end select

  end function get_beg

  !-----------------------------------------------------------------------
  function get_end(this, subgrid_level) result(end_index)
    !
    ! !DESCRIPTION:
    ! Get ending index for a given subgrid level
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer :: end_index  ! function result
    class(interp_bounds_type), intent(in) :: this
    character(len=*), intent(in) :: subgrid_level  ! 'pft', 'column' or 'landunit'
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_end'
    !-----------------------------------------------------------------------

    select case (subgrid_level)
    case('pft')
       end_index = this%endp
    case('column')
       end_index = this%endc
    case('landunit')
       end_index = this%endl
    case default
       call endrun(msg=subname//' ERROR: Unknown subgrid level: '//trim(subgrid_level)// &
            errMsg(__FILE__, __LINE__))
    end select

  end function get_end

end module initInterpBounds
