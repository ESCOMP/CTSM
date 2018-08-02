module WaterInfoBulkType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a class for working with information describing a given water instance, such
  ! as building history and restart field names.
  !
  ! This version is used for the bulk variables.
  !
  ! !USES:
  !
  use WaterInfoBaseType, only : water_info_base_type
  !
  implicit none
  private

  type, extends(water_info_base_type), public :: water_info_bulk_type
     private
   contains
     procedure, public :: fname  ! Get a history/restart field name
     procedure, public :: lname  ! Get a history/restart long name
  end type water_info_bulk_type

  interface water_info_bulk_type
     module procedure constructor
  end interface water_info_bulk_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  function constructor() result (this)
    ! Create a water_info_bulk_type object
    type(water_info_bulk_type) :: this  ! function result

    ! nothing to do
  end function constructor

  !-----------------------------------------------------------------------
  pure function fname(this, basename)
    !
    ! !DESCRIPTION:
    ! Get a history/restart field name
    !
    ! basename gives the base name of the history/restart field
    !
    ! !ARGUMENTS:
    character(len=:), allocatable :: fname  ! function result
    class(water_info_bulk_type), intent(in) :: this
    character(len=*), intent(in) :: basename
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'fname'
    !-----------------------------------------------------------------------

    fname = trim(basename)

  end function fname

  !-----------------------------------------------------------------------
  pure function lname(this, basename)
    !
    ! !DESCRIPTION:
    ! Get a history/restart long name
    !
    ! basename gives the base name of the history/restart long name
    !
    ! !ARGUMENTS:
    character(len=:), allocatable :: lname  ! function result
    class(water_info_bulk_type), intent(in) :: this
    character(len=*), intent(in) :: basename
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'lname'
    !-----------------------------------------------------------------------

    lname = trim(basename)

  end function lname

end module WaterInfoBulkType
