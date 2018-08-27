module WaterInfoIsotopeType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a class for working with information describing a given water instance, such
  ! as building history and restart field names.
  !
  ! This version is used for water isotopes. These are a subclass of general water
  ! tracers, but have some other information that is particular to isotopes.
  !
  ! !USES:
  !
  use WaterInfoTracerType, only : water_info_tracer_type
  !
  implicit none
  private

  type, extends(water_info_tracer_type), public :: water_info_isotope_type
     private

     ! Eventually, I envision this including a number of isotope-specific fields, like
     ! molecular weight, diffusivity ratio, etc. I envision these being private, with
     ! public getters. (That way, it's hidden to the user whether a given quantity is
     ! stored itself or calculated from other quantities. Also, that way we don't need to
     ! worry about someone accidentally re-setting one of these fields.) I'd like this
     ! class to be immutable: once it is created, its fields should never change - that
     ! way we don't need to worry about it possibly being changed when we pass it around
     ! via pointers.
  end type water_info_isotope_type

  interface water_info_isotope_type
     module procedure constructor
  end interface water_info_isotope_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  function constructor(tracer_name) result(this)
    ! Create a water_info_isotope_type object
    !
    ! Eventually, this will either (a) accept various arguments specifying information
    ! about this isotope (molecular weight, diffusivity ratio, etc.), or (b) look up this
    ! information from a lookup table defined here or elsewhere, based on the tracer_name.
    type(water_info_isotope_type) :: this  ! function result
    character(len=*), intent(in) :: tracer_name

    this%water_info_tracer_type = water_info_tracer_type(tracer_name)
  end function constructor

end module WaterInfoIsotopeType
