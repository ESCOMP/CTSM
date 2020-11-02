module ctsm_LilacLnd2AtmFieldListType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a class and related methods for a list of lilac fields sent from lnd -> atm.
  !
  ! (Note: this is very similar to LilacAtm2LndFieldListType. However, between the fact
  ! that (1) they have different sets of supported methods, and (2) the use of the
  ! dynamic vector, it seemed to make more sense to have totally separate classes rather
  ! than trying to share code between the two.)
  !
  ! To set up this list (lilac_lnd2atm_field_list_type):
  !
  ! - Initialize it by calling the 'init' method
  !
  ! - Add variables with add_var
  !
  ! - When done adding variables, call complete_setup
  !   - Note that you cannot access or perform any operations on any of the fields until
  !     this is done!
  !
  ! To use this list (after complete_setup has been called):
  !
  ! - Query number of fields with num_fields
  !
  ! - Extract data from a field with get_field
  !
  ! !USES:

  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : OOBMsg => shr_log_OOBMsg
  use shr_sys_mod    , only : shr_sys_abort
  use lilac_constants, only : field_index_unset, logunit

  implicit none
  private

  ! !PRIVATE TYPES:

  type, private :: lilac_lnd2atm_field_type
     private

     ! Metadata set initially in initialization
     character(len=:), allocatable :: fieldname
     character(len=:), allocatable :: units

     ! Metadata set later in initialization
     logical :: required_by_atm ! whether this field is actually required by the atmosphere

     ! Data set each time step
     real(r8), pointer :: dataptr(:)
  end type lilac_lnd2atm_field_type

  ! Define a dynamic vector for lilac_lnd2atm_field_type
#define VECTOR_NAME lilac_lnd2atm_field_vector
#define TYPE_NAME type(lilac_lnd2atm_field_type)
#define THROW(string) call shr_sys_abort(string)
#include "dynamic_vector_typedef.inc"

  !
  ! !PUBLIC TYPES:
  type, public :: lilac_lnd2atm_field_list_type
     private
     type(lilac_lnd2atm_field_vector) :: field_vec
     type(lilac_lnd2atm_field_type), allocatable :: fields(:)
   contains
     ! Methods for setting up the list:
     procedure, public :: init
     procedure, public :: add_var
     procedure, public :: complete_setup

     ! Methods to query or set data:
     procedure, public :: num_fields ! return the number of fields
     procedure, public :: get_field ! get data for one field
     procedure, public :: get_fieldname  ! get the field name for a given field
     procedure, public :: get_units      ! get the units for a given field
     procedure, public :: get_dataptr    ! get a pointer to the data for a given field

     ! Private methods:
     procedure, private :: check_field_index ! check whether a field index is valid
  end type lilac_lnd2atm_field_list_type

  interface lilac_lnd2atm_field_type
     module procedure new_lilac_lnd2atm_field_type
  end interface lilac_lnd2atm_field_type

contains

  ! Complete the dynamic vector definition.
#include "dynamic_vector_procdef.inc"

  !-----------------------------------------------------------------------
  function new_lilac_lnd2atm_field_type(fieldname, units) result(this)
    !
    ! !DESCRIPTION:
    ! Initialize a new lilac_lnd2atm_field_type object
    !
    ! !ARGUMENTS:
    type(lilac_lnd2atm_field_type) :: this  ! function result
    character(len=*), intent(in) :: fieldname
    character(len=*), intent(in) :: units
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'new_lilac_lnd2atm_field_type'
    !-----------------------------------------------------------------------

    this%fieldname = fieldname
    this%units = units

    ! Assume true until told otherwise
    this%required_by_atm = .true.

    nullify(this%dataptr)

  end function new_lilac_lnd2atm_field_type

  !-----------------------------------------------------------------------
  subroutine init(this)
    !
    ! !DESCRIPTION:
    ! Initialize a new lilac_lnd2atm_field_list_type object
    !
    ! !ARGUMENTS:
    class(lilac_lnd2atm_field_list_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'init'
    !-----------------------------------------------------------------------

    this%field_vec = lilac_lnd2atm_field_vector()

  end subroutine init

  !-----------------------------------------------------------------------
  subroutine add_var(this, fieldname, units, field_index)
    !
    ! !DESCRIPTION:
    ! Add the given field to this list
    !
    ! Also set field_index to be the index of this field in list. For the sake of error
    ! checking, field_index should be initialized to field_index_unset before the call to
    ! this subroutine.
    !
    ! !ARGUMENTS:
    class(lilac_lnd2atm_field_list_type), intent(inout) :: this
    character(len=*), intent(in) :: fieldname
    character(len=*), intent(in) :: units
    integer, intent(inout) :: field_index
    !
    ! !LOCAL VARIABLES:
    type(lilac_lnd2atm_field_type) :: one_field

    character(len=*), parameter :: subname = 'add_var'
    !-----------------------------------------------------------------------

    if (allocated(this%fields)) then
       write(logunit,*) subname//' ERROR: this%fields is already allocated.'
       write(logunit,*) 'fieldname = ', trim(fieldname)
       write(logunit,*) 'This is likely a sign that you are trying to add a variable'
       write(logunit,*) 'after complete_setup has already been called.'
       call shr_sys_abort('Attempt to call '//subname//' after complete_setup was called')
    end if

    if (field_index /= field_index_unset) then
       write(logunit,*) subname//' ERROR: attempt to add var with a field index that has already been set.'
       write(logunit,*) 'fieldname, field_index = ', trim(fieldname), field_index
       call shr_sys_abort('Attempt to add var with a field index that has already been set')
    end if

    one_field = lilac_lnd2atm_field_type( &
         fieldname = trim(fieldname), &
         units = trim(units))

    call this%field_vec%push_back(one_field)

    field_index = this%field_vec%vsize()
  end subroutine add_var

  !-----------------------------------------------------------------------
  subroutine complete_setup(this, data_size)
    !
    ! !DESCRIPTION:
    ! Finalize the creation of this field list; this includes allocating the data arrays for each field
    !
    ! !ARGUMENTS:
    class(lilac_lnd2atm_field_list_type), intent(inout) :: this
    integer, intent(in) :: data_size  ! number of points in each field (assumed to be the same for all fields)
    !
    ! !LOCAL VARIABLES:
    integer :: i

    character(len=*), parameter :: subname = 'complete_setup'
    !-----------------------------------------------------------------------

    call this%field_vec%move_out(this%fields)

    do i = 1, this%num_fields()
       allocate(this%fields(i)%dataptr(data_size))
    end do

  end subroutine complete_setup

  !-----------------------------------------------------------------------
  function num_fields(this)
    !
    ! !DESCRIPTION:
    ! Return the number of fields
    !
    ! !ARGUMENTS:
    integer :: num_fields  ! function result
    class(lilac_lnd2atm_field_list_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'num_fields'
    !-----------------------------------------------------------------------

    if (.not. allocated(this%fields)) then
       write(logunit,*) subname//' ERROR: this%fields has not yet been allocated'
       write(logunit,*) 'This is likely a sign that you are trying to call num_fields'
       write(logunit,*) 'before complete_setup has been called.'
       call shr_sys_abort('Attempt to get number of fields before complete_setup was called')
    end if

    num_fields = size(this%fields)

  end function num_fields

  !-----------------------------------------------------------------------
  subroutine get_field(this, field_index, data)
    !
    ! !DESCRIPTION:
    ! Get data for the given field
    !
    ! !ARGUMENTS:
    class(lilac_lnd2atm_field_list_type), intent(in) :: this
    integer, intent(in) :: field_index
    real(r8), intent(out) :: data(:)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_field'
    !-----------------------------------------------------------------------

    call this%check_field_index(field_index, subname)

    if (size(data) /= size(this%fields(field_index)%dataptr)) then
       call shr_sys_abort(subname//' field size mismatch for '//trim(this%fields(field_index)%fieldname))
    end if

    data(:) = this%fields(field_index)%dataptr(:)

  end subroutine get_field

  !-----------------------------------------------------------------------
  function get_fieldname(this, field_index) result(fieldname)
    !
    ! !DESCRIPTION:
    ! Get the field name for a given field
    !
    ! (This will already be trimmed - no further trimming is needed.)
    !
    ! !ARGUMENTS:
    character(len=:), allocatable :: fieldname  ! function result
    class(lilac_lnd2atm_field_list_type), intent(in) :: this
    integer, intent(in) :: field_index
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_fieldname'
    !-----------------------------------------------------------------------

    call this%check_field_index(field_index, subname)

    fieldname = this%fields(field_index)%fieldname

  end function get_fieldname

  !-----------------------------------------------------------------------
  function get_units(this, field_index) result(units)
    !
    ! !DESCRIPTION:
    ! Get the units for a given field
    !
    ! (This will already be trimmed - no further trimming is needed.)
    !
    ! !ARGUMENTS:
    character(len=:), allocatable :: units  ! function result
    class(lilac_lnd2atm_field_list_type), intent(in) :: this
    integer, intent(in) :: field_index
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_units'
    !-----------------------------------------------------------------------

    call this%check_field_index(field_index, subname)

    units = this%fields(field_index)%units

  end function get_units

  !-----------------------------------------------------------------------
  function get_dataptr(this, field_index) result(dataptr)
    !
    ! !DESCRIPTION:
    ! Get a pointer to the data for a given field
    !
    ! !ARGUMENTS:
    real(r8), pointer :: dataptr(:)  ! function result
    class(lilac_lnd2atm_field_list_type), intent(in) :: this
    integer, intent(in) :: field_index
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_dataptr'
    !-----------------------------------------------------------------------

    call this%check_field_index(field_index, subname)

    dataptr => this%fields(field_index)%dataptr

  end function get_dataptr

  !-----------------------------------------------------------------------
  subroutine check_field_index(this, field_index, caller)
    !
    ! !DESCRIPTION:
    ! Check the provided field_index for validity. If not valid, aborts.
    !
    ! !ARGUMENTS:
    class(lilac_lnd2atm_field_list_type), intent(in) :: this
    integer, intent(in) :: field_index
    character(len=*), intent(in) :: caller  ! name of caller, for error messages
    !
    ! !LOCAL VARIABLES:
    integer :: nfields

    character(len=*), parameter :: subname = 'check_field_index'
    !-----------------------------------------------------------------------

    if (field_index == field_index_unset) then
       call shr_sys_abort(caller//':'//subname//' attempt to set field for unset field index')
    end if

    if (field_index < 1 .or. field_index > this%num_fields()) then
       write(logunit,*) caller//':'//subname//' ERROR: field_index out of bounds'
       nfields = this%num_fields()
       write(logunit,*) 'field_index, num_fields = ', field_index, nfields
       call shr_sys_abort(caller//':'//subname//' field_index out of bounds')
    end if

  end subroutine check_field_index

end module ctsm_LilacLnd2AtmFieldListType
