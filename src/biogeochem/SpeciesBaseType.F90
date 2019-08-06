module SpeciesBaseType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a base class for working with chemical species, such as building history and
  ! restart field names.
  !
  ! !USES:
  !
  implicit none
  private

  ! !PUBLIC TYPES:

  type, abstract, public :: species_base_type
   contains
     ! Get a history field name for this species
     procedure(hist_fname_interface), public, deferred :: hist_fname

     ! Get a restart field name for this species
     procedure(rest_fname_interface), public, deferred :: rest_fname

     ! Get the full species name
     procedure(get_species_interface), public, deferred :: get_species

     ! Return true if this is an isotope, false if not
     procedure(is_isotope_interface), public, deferred :: is_isotope
  end type species_base_type

  abstract interface
     pure function hist_fname_interface(this, basename, suffix) result(fname)
       ! Get a history field name for this species
       !
       ! basename gives the base name of the history field
       !
       ! suffix, if provided, gives a suffix that appears after all species information
       ! in the field name
       import :: species_base_type

       character(len=:)         , allocatable :: fname  ! function result
       class(species_base_type) , intent(in)  :: this
       character(len=*)         , intent(in)  :: basename
       character(len=*)         , optional, intent(in) :: suffix
     end function hist_fname_interface

     function rest_fname_interface(this, basename, suffix) result(fname)
       ! Get a restart field name for this species
       !
       ! basename gives the base name of the restart field
       !
       ! suffix, if provided, gives a suffix that appears after all species information
       ! in the field name
       import :: species_base_type

       character(len=:)         , allocatable :: fname  ! function result
       class(species_base_type) , intent(in)  :: this
       character(len=*)         , intent(in)  :: basename
       character(len=*)         , optional, intent(in) :: suffix
     end function rest_fname_interface

     pure function get_species_interface(this) result(species_name)
       ! Get the full species name
       import :: species_base_type

       character(len=:), allocatable :: species_name
       class(species_base_type) , intent(in)  :: this
     end function get_species_interface

     pure function is_isotope_interface(this) result(is_isotope)
       ! Return true if this is an isotope, false if not
       import :: species_base_type

       logical :: is_isotope ! function result
       class(species_base_type), intent(in) :: this
     end function is_isotope_interface
  end interface

end module SpeciesBaseType
