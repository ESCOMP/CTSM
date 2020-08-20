module ctsm_NcdIoDim
  ! This module is specific to using the ncdio_pio_fake version of ncdio_pio. This
  ! provides a derived type for holding a single dimension from a fake netcdf file, and
  ! associated methods for working with this derived type.

  implicit none
  private
  save

  public :: ctsm_NcdIoDim_type

  integer, parameter, public :: max_name = 256   ! max length for a dimension name

  ! This type store a single dimension in a fake file
  type :: ctsm_NcdIoDim_type
     private
     character(len=max_name) :: dimname  ! dimension name
     integer :: dimlen  ! dimension length

   contains
     procedure :: get_dimname  ! get the dimension name
     procedure :: get_dimlen   ! get the dimension length
  end type ctsm_NcdIoDim_type

  interface ctsm_NcdIoDim_type
     module procedure constructor
  end interface ctsm_NcdIoDim_type

contains

  !-----------------------------------------------------------------------
  type(ctsm_NcdIoDim_type) function constructor(dimname, dimlen)
    ! Create a new object of type ctsm_NcdIoDim_type

    character(len=*) , intent(in) :: dimname  ! dimension name
    integer          , intent(in) :: dimlen   ! dimension length

    constructor%dimname = dimname
    constructor%dimlen = dimlen
  end function constructor

  !-----------------------------------------------------------------------
  character(len=max_name) function get_dimname(this)
    ! Get the name associated with this dimension
    class(ctsm_NcdIoDim_type), intent(in) :: this

    get_dimname = this%dimname
  end function get_dimname

  !-----------------------------------------------------------------------
  integer function get_dimlen(this)
    ! Get the length associated with this dimension
    class(ctsm_NcdIoDim_type), intent(in) :: this

    get_dimlen = this%dimlen
  end function get_dimlen

end module ctsm_NcdIoDim
