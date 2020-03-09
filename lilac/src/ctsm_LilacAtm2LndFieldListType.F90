module ctsm_LilacAtm2LndFieldListType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a class and related methods for a list of lilac fields sent from atm -> lnd.
  !
  ! (Note: this is very similar to LilacLnd2AtmFieldListType. However, between the fact
  ! that (1) they have different sets of supported methods, and (2) the use of the
  ! dynamic vector, it seemed to make more sense to have totally separate classes rather
  ! than trying to share code between the two.)
  !
  ! To set up this list (lilac_atm2lnd_field_list_type):
  !
  ! - Initialize it by calling the 'init' method
  !
  ! - Add variables with add_var
  !
  ! - When done adding variables, call complete_setup
  !   - Note that you cannot access or perform any operations on any of the fields until
  !     this is done!
  !
  ! - Set which fields are needed from data with set_needed_from_data
  !
  ! To use this list (after complete_setup has been called), here is the workflow:
  !
  ! - Query number of fields with num_fields
  !
  ! - Set fields each time step with set_field
  !
  ! - After fields are set, check that all required fields have been set by calling check_all_set
  !
  ! - At the end of each time step, call reset_provided
  !
  ! !USES:

  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : OOBMsg => shr_log_OOBMsg
  use shr_sys_mod    , only : shr_sys_abort
  use lilac_constants, only : field_index_unset, logunit

  implicit none
  private

  ! !PRIVATE TYPES:

  type, private :: lilac_atm2lnd_field_type
     private

     ! Metadata set initially in initialization
     character(len=:), allocatable :: fieldname
     character(len=:), allocatable :: units
     logical :: available_from_data  ! whether this field can be obtained from data if not provided by the sending component
     logical :: can_be_time_const    ! if true, it's okay for this field to be set just once, in the first time step, keeping its same value for the entire run (e.g., a landfrac field that doesn't vary in time)

     ! Metadata set later in initialization
     logical :: needed_from_data ! whether the host atmosphere wants LILAC to read this field from data
     logical :: required_by_lnd  ! whether this field is actually required by the land

     ! Data set each time step
     real(r8), pointer :: dataptr(:)
     logical :: provided_this_time   ! whether this variable has been set this time step
     logical :: provided_ever        ! whether this variable has ever been set
  end type lilac_atm2lnd_field_type

  ! Define a dynamic vector for lilac_atm2lnd_field_type
#define VECTOR_NAME lilac_atm2lnd_field_vector
#define TYPE_NAME type(lilac_atm2lnd_field_type)
#define THROW(string) call shr_sys_abort(string)
#include "dynamic_vector_typedef.inc"

  !
  ! !PUBLIC TYPES:
  type, public :: lilac_atm2lnd_field_list_type
     private
     type(lilac_atm2lnd_field_vector) :: field_vec
     type(lilac_atm2lnd_field_type), allocatable :: fields(:)
   contains
     ! Methods for setting up the list:
     procedure, public :: init
     procedure, public :: add_var
     procedure, public :: complete_setup

     ! Methods to query or set data:
     procedure, public :: num_fields ! return the number of fields
     procedure, public :: set_needed_from_data ! dictate that the given fields need to be read from data
     procedure, public :: is_needed_from_data ! query whether the given field is needed from data
     procedure, public :: set_field ! set data for one field
     procedure, public :: check_all_set  ! check to ensure that all required fields have been set this time
     procedure, public :: reset_provided ! reset the provided_this_time variable for all fields
     procedure, public :: get_fieldname  ! get the field name for a given field
     procedure, public :: get_units      ! get the units for a given field
     procedure, public :: get_dataptr    ! get a pointer to the data for a given field (this should be treated as read-only!)

     ! Private methods:
     procedure, private :: check_field_index ! check whether a field index is valid
  end type lilac_atm2lnd_field_list_type

  interface lilac_atm2lnd_field_type
     module procedure new_lilac_atm2lnd_field_type
  end interface lilac_atm2lnd_field_type

contains

  ! Complete the dynamic vector definition.
#include "dynamic_vector_procdef.inc"

  !-----------------------------------------------------------------------
  function new_lilac_atm2lnd_field_type(fieldname, units, available_from_data, can_be_time_const) result(this)
    !
    ! !DESCRIPTION:
    ! Initialize a new lilac_atm2lnd_field_type object
    !
    ! !ARGUMENTS:
    type(lilac_atm2lnd_field_type) :: this  ! function result
    character(len=*), intent(in) :: fieldname
    character(len=*), intent(in) :: units
    logical, intent(in) :: available_from_data ! whether this field can be obtained from data if not provided by the sending component
    logical, intent(in) :: can_be_time_const   ! if true, it's okay for this field to be set just once, in the first time step, keeping its same value for the entire run (e.g., a landfrac field that doesn't vary in time)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'new_lilac_atm2lnd_field_type'
    !-----------------------------------------------------------------------

    this%fieldname = fieldname
    this%units = units
    this%available_from_data = available_from_data
    this%can_be_time_const = can_be_time_const

    ! Assume false until told otherwise
    this%needed_from_data = .false.

    ! Assume true until told otherwise
    this%required_by_lnd = .true.

    nullify(this%dataptr)
    this%provided_this_time = .false.
    this%provided_ever      = .false.

  end function new_lilac_atm2lnd_field_type

  !-----------------------------------------------------------------------
  subroutine init(this)
    !
    ! !DESCRIPTION:
    ! Initialize a new lilac_atm2lnd_field_list_type object
    !
    ! !ARGUMENTS:
    class(lilac_atm2lnd_field_list_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'init'
    !-----------------------------------------------------------------------

    this%field_vec = lilac_atm2lnd_field_vector()

  end subroutine init

  !-----------------------------------------------------------------------
  subroutine add_var(this, fieldname, units, available_from_data, can_be_time_const, field_index)
    !
    ! !DESCRIPTION:
    ! Add the given field to this list
    !
    ! Also set field_index to be the index of this field in list. For the sake of error
    ! checking, field_index should be initialized to field_index_unset before the call to
    ! this subroutine.
    !
    ! !ARGUMENTS:
    class(lilac_atm2lnd_field_list_type), intent(inout) :: this
    character(len=*), intent(in) :: fieldname
    character(len=*), intent(in) :: units
    logical, intent(in) :: available_from_data ! whether this field can be obtained from data if not provided by the sending component
    logical, intent(in) :: can_be_time_const   ! if true, it's okay for this field to be set just once, in the first time step, keeping its same value for the entire run (e.g., a landfrac field that doesn't vary in time)
    integer, intent(inout) :: field_index
    !
    ! !LOCAL VARIABLES:
    type(lilac_atm2lnd_field_type) :: one_field

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

    one_field = lilac_atm2lnd_field_type( &
         fieldname = trim(fieldname), &
         units = trim(units), &
         available_from_data = available_from_data, &
         can_be_time_const = can_be_time_const)

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
    class(lilac_atm2lnd_field_list_type), intent(inout) :: this
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
    class(lilac_atm2lnd_field_list_type), intent(in) :: this
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
  subroutine set_needed_from_data(this, fields_needed_from_data)
    !
    ! !DESCRIPTION:
    ! Dictate that the given fields need to be read from data
    !
    ! !ARGUMENTS:
    class(lilac_atm2lnd_field_list_type), intent(inout) :: this
    integer, intent(in) :: fields_needed_from_data(:)  ! vector of field indices that need to be read from data
    !
    ! !LOCAL VARIABLES:
    integer :: i
    integer :: field_index

    character(len=*), parameter :: subname = 'set_needed_from_data'
    !-----------------------------------------------------------------------

    do i = 1, size(fields_needed_from_data)
       field_index = fields_needed_from_data(i)
       call this%check_field_index(field_index, subname)

       if (this%fields(field_index)%needed_from_data) then
          call shr_sys_abort(subname//' attempt to set needed_from_data on field for which it has already been set: '//&
               this%fields(field_index)%fieldname)
       end if

       if (.not. this%fields(field_index)%available_from_data) then
          call shr_sys_abort(subname//' attempt to set needed_from_data on field not available from data: '//&
               this%fields(field_index)%fieldname)
       end if

       this%fields(field_index)%needed_from_data = .true.
    end do

  end subroutine set_needed_from_data

  !-----------------------------------------------------------------------
  function is_needed_from_data(this, field_index)
    !
    ! !DESCRIPTION:
    ! Query whether the given field is needed from data
    !
    ! !ARGUMENTS:
    logical :: is_needed_from_data  ! function result
    class(lilac_atm2lnd_field_list_type), intent(in) :: this
    integer, intent(in) :: field_index
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'is_needed_from_data'
    !-----------------------------------------------------------------------

    call this%check_field_index(field_index, subname)

    is_needed_from_data = this%fields(field_index)%needed_from_data

  end function is_needed_from_data

  !-----------------------------------------------------------------------
  subroutine set_field(this, field_index, data)
    !
    ! !DESCRIPTION:
    ! Set data for the given field
    !
    ! It is an error to try to set a field that has already been set this time
    !
    ! !ARGUMENTS:
    class(lilac_atm2lnd_field_list_type), intent(inout) :: this
    integer, intent(in) :: field_index
    real(r8), intent(in) :: data(:)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'set_field'
    !-----------------------------------------------------------------------

    call this%check_field_index(field_index, subname)

    if (size(data) /= size(this%fields(field_index)%dataptr)) then
       call shr_sys_abort(subname//' field size mismatch for '//trim(this%fields(field_index)%fieldname))
    end if

    if (this%fields(field_index)%provided_this_time) then
       ! This can typically happen in one of two ways:
       ! - A component tries to re-set a field that has already been set
       ! - reset_provided hasn't been called in between times
       write(logunit,*) subname//' ERROR: attempt to set an already-set field: ', this%fields(field_index)%fieldname
       if (this%fields(field_index)%needed_from_data) then
          write(logunit,*) "This field was marked as being needed from data."
          write(logunit,*) "A possible cause of this error is that it is being set by both"
          write(logunit,*) "the host atmosphere and LILAC's internal data atmosphere."
       end if
       call shr_sys_abort(subname//' attempt to set an already-set field: '//this%fields(field_index)%fieldname)
    end if

    this%fields(field_index)%dataptr(:) = data(:)
    this%fields(field_index)%provided_this_time = .true.
    this%fields(field_index)%provided_ever = .true.

  end subroutine set_field

  !-----------------------------------------------------------------------
  subroutine check_all_set(this)
    !
    ! !DESCRIPTION:
    ! Check to ensure that all required fields have been set this time
    !
    ! !ARGUMENTS:
    class(lilac_atm2lnd_field_list_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:
    integer :: i

    character(len=*), parameter :: subname = 'check_all_set'
    !-----------------------------------------------------------------------

    do i = 1, this%num_fields()
       if (this%fields(i)%required_by_lnd) then
          if (.not. this%fields(i)%provided_ever) then
             call shr_sys_abort(trim(this%fields(i)%fieldname)//' required but never provided')
          end if
          if (.not. this%fields(i)%can_be_time_const) then
             if (.not. this%fields(i)%provided_this_time) then
                call shr_sys_abort(trim(this%fields(i)%fieldname)//' required but not provided this time')
             end if
          end if
       end if
    end do

  end subroutine check_all_set

  !-----------------------------------------------------------------------
  subroutine reset_provided(this)
    !
    ! !DESCRIPTION:
    ! Reset the provided_this_time variable for all fields
    !
    ! !ARGUMENTS:
    class(lilac_atm2lnd_field_list_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:
    integer :: i

    character(len=*), parameter :: subname = 'reset_provided'
    !-----------------------------------------------------------------------

    do i = 1, this%num_fields()
       this%fields(i)%provided_this_time = .false.
    end do

  end subroutine reset_provided

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
    class(lilac_atm2lnd_field_list_type), intent(in) :: this
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
    class(lilac_atm2lnd_field_list_type), intent(in) :: this
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
    ! This should be treated as read-only! Setting data should be done via the provided
    ! methods in this class.
    !
    ! !ARGUMENTS:
    real(r8), pointer :: dataptr(:)  ! function result
    class(lilac_atm2lnd_field_list_type), intent(in) :: this
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
    class(lilac_atm2lnd_field_list_type), intent(in) :: this
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

end module ctsm_LilacAtm2LndFieldListType
