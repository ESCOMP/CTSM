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
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use WaterInfoBaseType, only : water_info_base_type
  !
  implicit none
  private

  type, extends(water_info_base_type), public :: water_info_tracer_type
     private
     character(len=:), allocatable :: tracer_name

     ! If true, this tracer is received from and sent to the coupler. If false, this
     ! tracer is just used internally in CTSM, and is set to some fixed ratio times the
     ! bulk water.
     logical :: communicated_with_coupler
     logical :: included_in_consistency_check
   contains
     procedure, public :: fname  ! Get a history/restart field name for this tracer
     procedure, public :: lname  ! Get a history/restart long name for this tracer
     procedure, public :: is_communicated_with_coupler
     procedure, public :: is_included_in_consistency_check
  end type water_info_tracer_type

  interface water_info_tracer_type
     module procedure constructor
  end interface water_info_tracer_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  function constructor(tracer_name, ratio, included_in_consistency_check, &
                       communicated_with_coupler) result(this)
    ! Create a water_info_tracer_type object
    type(water_info_tracer_type) :: this  ! function result
    character(len=*), intent(in) :: tracer_name
    real(r8), intent(in)         :: ratio
    logical, intent(in)          :: included_in_consistency_check

    ! If true, this tracer is received from and sent to the coupler. If false, this tracer
    ! is just used internally in CTSM, and is set to some fixed ratio times the bulk
    ! water.
    logical, intent(in)          :: communicated_with_coupler

    this%tracer_name = trim(tracer_name)
    this%communicated_with_coupler = communicated_with_coupler
    this%included_in_consistency_check = included_in_consistency_check
    call this%set_metadata(ratio = ratio)
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

  !-----------------------------------------------------------------------
  pure function is_communicated_with_coupler(this) result(coupled)
    !
    ! !DESCRIPTION:
    ! Returns true if this tracer is received from and sent to the coupler. Returns false
    ! if this tracer is just used internally in CTSM, and is set to some fixed ratio
    ! times the bulk water.
    !
    ! !ARGUMENTS:
    logical :: coupled  ! function result
    class(water_info_tracer_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'is_communicated_with_coupler'
    !-----------------------------------------------------------------------

    coupled = this%communicated_with_coupler

  end function is_communicated_with_coupler

  pure function is_included_in_consistency_check(this) result(included)

    logical :: included  ! function result
    class(water_info_tracer_type), intent(in) :: this

    included = this%included_in_consistency_check

  end function is_included_in_consistency_check


end module WaterInfoTracerType
