module WaterInfoBaseType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a base class for working with information describing a given water instance
  ! (bulk or tracer), such as building history and restart field names.
  !
  ! !USES:
  !
  implicit none
  private

  ! !PUBLIC TYPES:

  type, abstract, public :: water_info_base_type
   contains
     ! Get a history/restart field name for this tracer (or bulk)
     procedure(fname_interface), public, deferred :: fname

     ! Get a history/restart long name for this tracer (or bulk)
     procedure(lname_interface), public, deferred :: lname
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
end module WaterInfoBaseType
