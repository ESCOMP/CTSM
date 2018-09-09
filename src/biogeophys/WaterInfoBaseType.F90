module WaterInfoBaseType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a base class for working with information describing a given water instance
  ! (bulk or tracer), such as building history and restart field names.
  !
  ! !USES:
  !
  use shr_kind_mod            , only : r8 => shr_kind_r8

  implicit none
  private

  ! !PUBLIC TYPES:

  type, abstract, public :: water_info_base_type

     private
     real(r8) :: ratio

   contains
     ! Get a history/restart field name for this tracer (or bulk)
     procedure(fname_interface), public, deferred :: fname

     ! Get a history/restart long name for this tracer (or bulk)
     procedure(lname_interface), public, deferred :: lname
     procedure :: get_ratio
     procedure :: set_metadata

  end type water_info_base_type

  abstract interface
     pure function fname_interface(this, basename) result(fname)
       ! Get a history/restart field name for this tracer (or bulk)
       !
       ! basename gives the base name of the history/restart field
       import :: water_info_base_type

       character(len=:)            , allocatable :: fname
       class(water_info_base_type) , intent(in)  :: this
       character(len=*)            , intent(in)  :: basename
     end function fname_interface

     pure function lname_interface(this, basename) result(lname)
       ! Get a history/restart long name for this tracer (or bulk)
       !
       ! basename gives the base name of the history/restart long name
       import :: water_info_base_type

       character(len=:)            , allocatable :: lname
       class(water_info_base_type) , intent(in)  :: this
       character(len=*)            , intent(in)  :: basename
     end function lname_interface
  end interface

contains

  subroutine set_metadata(this,ratio) 
     ! !ARGUMENTS:
     class(water_info_base_type) , intent(inout)  :: this
     real(r8), intent(in), optional :: ratio

     if (present(ratio)) then
        this%ratio = ratio
     end if
  end subroutine

  function get_ratio(this) result(ratio)
     ! !ARGUMENTS:
     class(water_info_base_type) , intent(in)  :: this
     real(r8) :: ratio  ! function result

     ratio = this%ratio
  end function

end module WaterInfoBaseType
