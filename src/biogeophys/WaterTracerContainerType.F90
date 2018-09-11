module WaterTracerContainerType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a class and related methods for a container of pointers to all water tracer
  ! variables. This allows looping through all variables to do some operation on each
  ! variable in turn.
  !
  ! To use the container (water_tracer_container_type):
  !
  ! - Initialize it by calling the 'init' method
  !
  ! - Add variables with add_var
  !
  ! - When done adding variables, call complete_setup()
  !
  ! - Then use get_num_vars to get the number of variables, and use any of the other
  !   routines defined in water_tracer_container_type to extract data
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : OOBMsg => shr_log_OOBMsg
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use abortutils     , only : endrun
  use decompMod      , only : bounds_type, get_beg, get_end
  use clm_varctl     , only : iulog
  !
  implicit none
  private

  ! !PRIVATE TYPES:

  ! A single water tracer
  !
  ! Note that multi-level tracers are represented via multiple instances of this type
  type, private :: water_tracer_type
     private
     real(r8), pointer :: data(:)
     character(len=:), allocatable :: description  ! description of the variable (typically, the variable name, and optionally some other information like the level index)
     integer :: subgrid_level  ! one of the level codes defined in decompMod
  end type water_tracer_type

  ! Define a dynamic vector for water_tracer_type
#define VECTOR_NAME water_tracer_vector
#define TYPE_NAME type(water_tracer_type)
#define THROW(string) call endrun(msg=string)
#include "dynamic_vector_typedef.inc"

  !
  ! !PUBLIC TYPES:
  type, public :: water_tracer_container_type
     private
     type(water_tracer_vector) :: tracer_vec
     type(water_tracer_type), allocatable :: tracers(:)
   contains
     procedure, public :: init
     procedure, public :: add_var
     procedure, public :: complete_setup
     procedure, public :: get_num_vars
     procedure, public :: get_description
     procedure, public :: get_bounds
     procedure, public :: get_data
  end type water_tracer_container_type

  interface water_tracer_type
     module procedure new_water_tracer_type
  end interface water_tracer_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  ! Complete the dynamic vector definition.
#include "dynamic_vector_procdef.inc"

  !-----------------------------------------------------------------------
  function new_water_tracer_type(data, begi, description, subgrid_level) result(this)
    !
    ! !DESCRIPTION:
    ! Create a water_tracer_type object
    !
    ! !ARGUMENTS:
    type(water_tracer_type) :: this  ! function result
    integer, intent(in)          :: begi           ! beginning index of data array
    real(r8), target, intent(in) :: data(begi:)
    character(len=*), intent(in) :: description
    integer         , intent(in) :: subgrid_level  ! one of the levels defined in decompMod  
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'new_water_tracer_type'
    !-----------------------------------------------------------------------

    this%data => data
    this%description = description
    this%subgrid_level = subgrid_level

  end function new_water_tracer_type

  !-----------------------------------------------------------------------
  subroutine init(this)
    !
    ! !DESCRIPTION:
    ! Initialize a new water_tracer_container_type object
    !
    ! !ARGUMENTS:
    class(water_tracer_container_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'init'
    !-----------------------------------------------------------------------

    this%tracer_vec = water_tracer_vector()

  end subroutine init

  !-----------------------------------------------------------------------
  subroutine add_var(this, var, begi, description, subgrid_level)
    !
    ! !DESCRIPTION:
    ! Add the given variable to this container
    !
    ! !ARGUMENTS:
    class(water_tracer_container_type), intent(inout) :: this
    integer, intent(in)          :: begi           ! beginning index of var array
    real(r8), target, intent(in) :: var(begi:)
    character(len=*), intent(in) :: description
    integer         , intent(in) :: subgrid_level  ! one of the levels defined in decompMod
    !
    ! !LOCAL VARIABLES:
    type(water_tracer_type) :: one_tracer

    character(len=*), parameter :: subname = 'add_var'
    !-----------------------------------------------------------------------

    if (allocated(this%tracers)) then
       write(iulog,*) subname//' ERROR: this%tracers is already allocated.'
       write(iulog,*) 'This is likely a sign that you are trying to add a variable'
       write(iulog,*) 'after complete_setup has already been called.'
       call endrun(msg='Attempt to call '//subname//' after complete_setup was called', &
            additional_msg=errMsg(sourcefile, __LINE__))
    end if

    one_tracer = water_tracer_type( &
         data = var, &
         begi = begi, &
         description = description, &
         subgrid_level = subgrid_level)

    call this%tracer_vec%push_back(one_tracer)

  end subroutine add_var

  !-----------------------------------------------------------------------
  subroutine complete_setup(this)
    !
    ! !DESCRIPTION:
    ! Finalize the creation of this container
    !
    ! This must be called when all add_var calls are complete
    !
    ! !ARGUMENTS:
    class(water_tracer_container_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'complete_setup'
    !-----------------------------------------------------------------------

    call this%tracer_vec%move_out(this%tracers)

  end subroutine complete_setup

  !-----------------------------------------------------------------------
  function get_num_vars(this) result(num_vars)
    !
    ! !DESCRIPTION:
    ! Returns the number of variables held here
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer :: num_vars  ! function result
    class(water_tracer_container_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_num_vars'
    !-----------------------------------------------------------------------

    ! We could put this check in all of the getters. But typically, this get_num_vars
    ! routine will be called before any of the others, so for efficiency, we just have
    ! the check here.
    if (.not. allocated(this%tracers)) then
       write(iulog,*) subname//' ERROR: this%tracers is not yet allocated.'
       write(iulog,*) 'This is likely a sign that complete_setup was not called.'
       call endrun(msg='Attempt to call '//subname//' without calling complete_setup', &
            additional_msg=errMsg(sourcefile, __LINE__))
    end if

    num_vars = size(this%tracers)

  end function get_num_vars

  !-----------------------------------------------------------------------
  function get_description(this, var_num) result(description)
    !
    ! !DESCRIPTION:
    ! Returns the description of the variable with the given index
    !
    ! !ARGUMENTS:
    character(len=:), allocatable :: description  ! function result
    class(water_tracer_container_type), intent(in) :: this
    integer, intent(in) :: var_num
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_description'
    !-----------------------------------------------------------------------

    description = this%tracers(var_num)%description

  end function get_description

  !-----------------------------------------------------------------------
  subroutine get_bounds(this, var_num, bounds, begi, endi)
    !
    ! !DESCRIPTION:
    ! Gets the bounds of the variable with the given index, based on its subgrid level
    !
    ! !ARGUMENTS:
    class(water_tracer_container_type), intent(in) :: this
    integer, intent(in) :: var_num
    type(bounds_type), intent(in) :: bounds
    integer, intent(out) :: begi
    integer, intent(out) :: endi
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_bounds'
    !-----------------------------------------------------------------------

    begi = get_beg(bounds, this%tracers(var_num)%subgrid_level)
    endi = get_end(bounds, this%tracers(var_num)%subgrid_level)

  end subroutine get_bounds

  !-----------------------------------------------------------------------
  subroutine get_data(this, var_num, data)
    !
    ! !DESCRIPTION:
    ! Set a pointer to the data associated with the given variable
    !
    ! Assumes that the 'data' pointer is not currently allocated. (Otherwise this will
    ! result in a memory leak. It is okay for the data pointer to be previously
    ! associated with something else, though, as long as it doesn't require deallocation
    ! before being associated with something new.)
    !
    ! !ARGUMENTS:
    class(water_tracer_container_type), intent(in) :: this
    integer, intent(in) :: var_num
    real(r8), pointer, intent(out) :: data(:)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_data'
    !-----------------------------------------------------------------------

    data => this%tracers(var_num)%data

  end subroutine get_data

end module WaterTracerContainerType
