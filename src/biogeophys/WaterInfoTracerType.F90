module WaterInfoTracerType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a class for working with information describing a given water instance, such
  ! as building history and restart field names.
  !
  ! This version is used for water tracers.
  !
  ! !USES:
  !
  use WaterInfoBaseType, only : water_info_base_type
  !
  implicit none
  private

  type, extends(water_info_base_type), public :: water_info_tracer_type
     private
     character(len=:), allocatable :: tracer_name
   contains
     procedure, public :: fname  ! Get a history/restart field name for this tracer
     procedure, public :: lname  ! Get a history/restart long name for this tracer
  end type water_info_tracer_type

  interface water_info_tracer_type
     module procedure constructor
  end interface water_info_tracer_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  function constructor(tracer_name) result(this)
    ! Create a water_info_tracer_type object
    type(water_info_tracer_type) :: this  ! function result
    character(len=*), intent(in) :: tracer_name

    this%tracer_name = trim(tracer_name)
  end function constructor

  !-----------------------------------------------------------------------
  pure function fname(this, basename)
    !
    ! !DESCRIPTION:
    ! Get a history/restart field name for this tracer
    !
    ! basename gives the base name of the history/restart field
    !
    ! !ARGUMENTS:
    character(len=:), allocatable :: fname  ! function result
    class(water_info_tracer_type), intent(in) :: this
    character(len=*), intent(in) :: basename
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'fname'
    !-----------------------------------------------------------------------

    fname = trim(basename) // '_' // this%tracer_name

  end function fname

  !-----------------------------------------------------------------------
  pure function lname(this, basename)
    !
    ! !DESCRIPTION:
    ! Get a history/restart long name for this tracer
    !
    ! basename gives the base name of the history/restart long name
    !
    ! !ARGUMENTS:
    character(len=:), allocatable :: lname  ! function result
    class(water_info_tracer_type), intent(in) :: this
    character(len=*), intent(in) :: basename
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'lname'
    !-----------------------------------------------------------------------

    lname = this%tracer_name // ' ' // trim(basename)

  end function lname

end module WaterInfoTracerType
