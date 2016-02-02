module dynColumnStateUpdaterMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Class for adjusting column-level state variables due to transient column areas.
  !
  ! In each time step, the object should be set up with:
  !
  !    - call column_state_updater%set_old_weights (before dyn subgrid weight updates)
  !
  !    - call column_state_updater%set_new_weights (after dyn subgrid weight updates)
  !
  ! Then it can be used to update each state variable with a call to one of the variants of:
  !
  !    - call column_state_updater%update_column_state_*
  !
  ! The following methods are available for state updates:
  !
  !    - update_column_state_no_special_handling
  !        No special handling is done for any columns. This method is appropriate for
  !        state variables that have valid values for all landunits.
  !
  !    - update_column_state_fill_special_using_natveg
  !        Columns on any "special" landunit contribute mass-per-unit-area equal to the
  !        current state on the first natural vegetation column on its grid cell. This
  !        method is appropriate for state variables that are not tracked on special
  !        landunits, when you want their value to remain unchanged when a vegetated
  !        landunit takes over area from a special landunit.
  !
  !    - update_column_state_fill_using_fixed_values
  !        Each specially-treated landunit is given its own constant value for the state
  !        variable, which is the mass-per-unit-area contributed when a column on that
  !        landunit shrinks. This method is appropriate when you want to treat certain
  !        landunits as contributing some fixed amount of mass-per-unit-area when they
  !        shrink - this fixed amount can be 0 or some non-0 value.
  !
  ! For methods other than update_column_state_no_special_handling, an additional inout
  ! argument (non_conserved_mass_grc) accumulates the non-conserved mass due to shrinking
  ! or growing specially-treated columns.
  !
  ! !USES:
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use abortutils           , only : endrun
  use clm_varctl           , only : iulog  
  use clm_varcon           , only : namec, spval
  use decompMod            , only : bounds_type, BOUNDS_LEVEL_PROC
  use ColumnType           , only : col
  use LandunitType         , only : lun
  use landunit_varcon      , only : max_lunit
  use dynColumnTemplateMod , only : template_col_from_natveg_array, TEMPLATE_NONE_FOUND
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  public :: column_state_updater_type

  type column_state_updater_type
     private
     real(r8), allocatable :: cwtgcell_old(:)  ! old column weights on the grid cell
     real(r8), allocatable :: cwtgcell_new(:)  ! new column weights on the grid cell

     ! (cwtgcell_new - cwtgcell_old) from last call to set_new_weights
     real(r8), allocatable :: area_gained_col(:)

     ! For each column, a 'template' column determined as: the first active column on the
     ! natural veg landunit in the same grid cell as the target column. 'active' is
     ! determined at the time of the call to set_old_weights, so that we consider whether
     ! a column was active in the previous time step, rather than newly-active.
     integer , allocatable :: natveg_template_col(:)
   contains
     ! Public routines
     procedure, public :: set_old_weights     ! set weights before dyn subgrid updates
     procedure, public :: set_new_weights     ! set weights after dyn subgrid updates

     ! Various ways to update a column-level state variable based on changing column areas:
     procedure, public :: update_column_state_no_special_handling
     procedure, public :: update_column_state_fill_special_using_natveg
     procedure, public :: update_column_state_fill_using_fixed_values

     ! Private routines
     procedure, private :: update_column_state  ! do the work of updating a column state
  end type column_state_updater_type

  interface column_state_updater_type
     module procedure constructor  ! initialize a column_state_updater_type object
  end interface column_state_updater_type

  ! !PUBLIC VARIABLES:
  ! For update_column_state_fill_using_fixed_values, any landunit with
  ! FILLVAL_USE_EXISTING_VALUE will use the existing value in the state variable
  real(r8), parameter, public :: FILLVAL_USE_EXISTING_VALUE = spval

contains

  ! ========================================================================
  ! Constructors
  ! ========================================================================

  !-----------------------------------------------------------------------
  function constructor(bounds)
    !
    ! !DESCRIPTION:
    ! Initialize a column_state_updater_type object
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    type(column_state_updater_type) :: constructor  ! function result
    type(bounds_type), intent(in)   :: bounds       ! processor bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'constructor'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, errMsg(__FILE__, __LINE__))

    allocate(constructor%cwtgcell_old(bounds%begc:bounds%endc))
    constructor%cwtgcell_old(:) = nan
    allocate(constructor%cwtgcell_new(bounds%begc:bounds%endc))
    constructor%cwtgcell_new(:) = nan
    allocate(constructor%area_gained_col(bounds%begc:bounds%endc))
    constructor%area_gained_col(:) = nan
    allocate(constructor%natveg_template_col(bounds%begc:bounds%endc))
    constructor%natveg_template_col(:) = TEMPLATE_NONE_FOUND

  end function constructor

  ! ========================================================================
  ! Public methods
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine set_old_weights(this, bounds)
    !
    ! !DESCRIPTION:
    ! Set subgrid weights before dyn subgrid updates
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(column_state_updater_type) , intent(inout) :: this
    type(bounds_type)                , intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: c

    character(len=*), parameter :: subname = 'set_old_weights'
    !-----------------------------------------------------------------------

    do c = bounds%begc, bounds%endc
       this%cwtgcell_old(c) = col%wtgcell(c)
    end do

    call template_col_from_natveg_array(bounds, col%active(bounds%begc:bounds%endc), &
         this%natveg_template_col(bounds%begc:bounds%endc))

  end subroutine set_old_weights

  !-----------------------------------------------------------------------
  subroutine set_new_weights(this, bounds)
    !
    ! !DESCRIPTION:
    ! Set subgrid weights after dyn subgrid updates
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(column_state_updater_type) , intent(inout) :: this
    type(bounds_type)                , intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: c, g

    character(len=*), parameter :: subname = 'set_new_weights'
    !-----------------------------------------------------------------------

    do c = bounds%begc, bounds%endc
       g = col%gridcell(c)
       this%cwtgcell_new(c)     = col%wtgcell(c)
       this%area_gained_col(c)  = this%cwtgcell_new(c) - this%cwtgcell_old(c)
    end do

  end subroutine set_new_weights

  !-----------------------------------------------------------------------
  subroutine update_column_state_no_special_handling(this, bounds, var)
    !
    ! !DESCRIPTION:
    ! Adjust the values of a column-level state variable due to changes in subgrid
    ! weights.
    !
    ! This method does no special handling of any columns.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(column_state_updater_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(inout) :: var( bounds%begc: ) ! column-level variable
    !
    ! !LOCAL VARIABLES:
    real(r8) :: vals_input(bounds%begc:bounds%endc)
    logical  :: vals_input_valid(bounds%begc:bounds%endc)
    logical  :: has_prognostic_state(bounds%begc:bounds%endc)
    real(r8) :: non_conserved_mass(bounds%begg:bounds%endg)
    character(len=:), allocatable :: err_msg
    real(r8), parameter :: conservation_tolerance = 1.e-12_r8

    character(len=*), parameter :: subname = 'update_column_state_no_special_handling'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(var) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))

    vals_input(bounds%begc:bounds%endc) = var(bounds%begc:bounds%endc)
    vals_input_valid(bounds%begc:bounds%endc) = .true.
    has_prognostic_state(bounds%begc:bounds%endc) = .true.
    non_conserved_mass(bounds%begg:bounds%endg) = 0._r8

    call this%update_column_state(&
         bounds = bounds, &
         vals_input = vals_input(bounds%begc:bounds%endc), &
         vals_input_valid = vals_input_valid(bounds%begc:bounds%endc), &
         has_prognostic_state = has_prognostic_state(bounds%begc:bounds%endc), &
         var = var(bounds%begc:bounds%endc), &
         non_conserved_mass = non_conserved_mass(bounds%begg:bounds%endg))

    ! Since there is no special handling in this routine, the non_conserved_mass variable
    ! should not have any accumulation. We allow for roundoff-level accumulation in case
    ! non-conserved mass is determined in a way that is prone to roundoff-level errors.
    err_msg = subname//': ERROR: failure to conserve mass when using no special handling'
    SHR_ASSERT_ALL(abs(non_conserved_mass(bounds%begg:bounds%endg)) < conservation_tolerance, err_msg)

  end subroutine update_column_state_no_special_handling

  !-----------------------------------------------------------------------
  subroutine update_column_state_fill_special_using_natveg(this, bounds, var, &
       non_conserved_mass_grc)
    !
    ! !DESCRIPTION:
    ! Adjust the values of a column-level state variable due to changes in subgrid
    ! weights.
    !
    ! In this method, any shrinking column on a special landunit contributes state equal
    ! to the first natural vegetation column on its grid cell.
    !
    ! For the sake of determining non-conserved mass, special landunits are treated as
    ! having a value of 0 - i.e., any non-zero quantity in var(c) (for c in a special
    ! landunit) is ignored. Thus, even though a shrinking special column can contribute a
    ! non-zero value, this non-zero state is treated as having been created out of thin
    ! air for conservation-tracking purposes.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(column_state_updater_type), intent(in) :: this
    type(bounds_type) , intent(in)    :: bounds
    real(r8)          , intent(inout) :: var( bounds%begc: ) ! column-level variable

    ! Mass lost (per unit of grid cell area) from each grid cell due to changing area of
    ! special landunits. When a special landunit grows, it throws away all mass that would
    ! be added to it; when a special landunit shrinks, it pulls some virtual mass out of
    ! thin air (as given by the quantity in the vegetated column in that grid cell); this
    ! variable accounts for both of these cases. Positive denotes mass lost from the grid
    ! cell, negative denotes mass gained by the grid cell.
    real(r8)          , intent(inout) :: non_conserved_mass_grc( bounds%begg: )
    !
    ! !LOCAL VARIABLES:
    integer  :: c, l
    integer  :: template_col
    real(r8) :: vals_input(bounds%begc:bounds%endc)
    logical  :: vals_input_valid(bounds%begc:bounds%endc)
    logical  :: has_prognostic_state(bounds%begc:bounds%endc)

    character(len=*), parameter :: subname = 'update_column_state_fill_special_using_natveg'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(var) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(non_conserved_mass_grc) == (/bounds%begg/)), errMsg(__FILE__, __LINE__))

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          has_prognostic_state(c) = .false.

          template_col = this%natveg_template_col(c)
          if (template_col == TEMPLATE_NONE_FOUND) then
             vals_input(c) = spval
             vals_input_valid(c) = .false.
          else
             vals_input(c) = var(template_col)
             vals_input_valid(c) = .true.
          end if
       else
          has_prognostic_state(c) = .true.
          vals_input(c) = var(c)
          vals_input_valid(c) = .true.
       end if
    end do

    call this%update_column_state(&
         bounds = bounds, &
         vals_input = vals_input(bounds%begc:bounds%endc), &
         vals_input_valid = vals_input_valid(bounds%begc:bounds%endc), &
         has_prognostic_state = has_prognostic_state(bounds%begc:bounds%endc), &
         var = var(bounds%begc:bounds%endc), &
         non_conserved_mass = non_conserved_mass_grc(bounds%begg:bounds%endg))

  end subroutine update_column_state_fill_special_using_natveg

  !-----------------------------------------------------------------------
  subroutine update_column_state_fill_using_fixed_values(this, bounds, var, &
       landunit_values, non_conserved_mass_grc)
    !
    ! !DESCRIPTION:
    ! Adjust the values of a column-level state variable due to changes in subgrid
    ! weights.
    !
    ! In this method, any column in landunit type i is assumed to have state equal to
    ! landunit_values(i). This is used when such a column is shrinking. However, if
    ! landunit_values(i) is FILLVAL_USE_EXISTING_VALUE, then columns in landunit type i
    ! keep their existing value in the state variable (var).
    !
    ! Landunits with landunit_values(i) = FILLVAL_USE_EXISTING_VALUE can accept mass,
    ! others cannot.
    !
    ! For example: If landunit_values(istwet) = 30, then any shrinking wetland column
    ! provides a mass-per-unit-area of 30 for this state variable, and growing wetland
    ! columns cannot except mass. On the other hand, if landunit_values(istwet) =
    ! FILLVAL_USE_EXISTING_VALUE, then shrinking wetland columns provide a
    ! mass-per-unit-area equal to var for the given column, and growing wetland columns
    ! can accept mass for this state variable.
    !
    ! For the sake of determining non-conserved mass, specially-treated landunits (those
    ! with landunit_values(i) /= FILLVAL_USE_EXISTING_VALUE) are treated as having a value
    ! of 0 - i.e., any non-zero quantity in var(c) (for c in a specially-treated landunit)
    ! is ignored. Thus, even though a shrinking special column can contribute a
    ! non-zero value, this non-zero state is treated as having been created out of thin
    ! air for conservation-tracking purposes.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(column_state_updater_type), intent(in) :: this
    type(bounds_type) , intent(in)    :: bounds
    real(r8)          , intent(inout) :: var( bounds%begc: ) ! column-level variable
    real(r8)          , intent(in)    :: landunit_values(:)  ! value to use as input for each landunit type

    ! Mass lost (per unit of grid cell area) from each grid cell due to changing area of
    ! specially-treated landunits - any landunit with landunit_values(i) /=
    ! FILLVAL_USE_EXISTING_VALUE. When a specially-treated landunit grows, it throws away
    ! all mass that would be added to it; when a specially-treated landunit shrinks, it
    ! pulls some virtual mass out of thin air (as given by landunit_values); this variable
    ! accounts for both of these cases. Positive denotes mass lost from the grid cell,
    ! negative denotes mass gained by the grid cell.
    real(r8)          , intent(inout) :: non_conserved_mass_grc( bounds%begg: )
    !
    ! !LOCAL VARIABLES:
    character(len=:), allocatable :: err_msg
    integer  :: c, l
    integer  :: ltype
    real(r8) :: my_fillval
    real(r8) :: vals_input(bounds%begc:bounds%endc)
    logical  :: vals_input_valid(bounds%begc:bounds%endc)
    logical  :: has_prognostic_state(bounds%begc:bounds%endc)

    character(len=*), parameter :: subname = 'update_column_state_fill_using_fixed_values'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(var) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))
    err_msg = subname//': must provide values for all landunits'
    SHR_ASSERT((size(landunit_values) == max_lunit), err_msg)
    SHR_ASSERT_ALL((ubound(non_conserved_mass_grc) == (/bounds%begg/)), errMsg(__FILE__, __LINE__))

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       ltype = lun%itype(l)
       my_fillval = landunit_values(ltype)

       if (my_fillval == FILLVAL_USE_EXISTING_VALUE) then
          vals_input(c) = var(c)
          vals_input_valid(c) = .true.
          has_prognostic_state(c) = .true.
       else
          vals_input(c) = my_fillval
          vals_input_valid(c) = .true.
          has_prognostic_state(c) = .false.
       end if
    end do

    call this%update_column_state( &
         bounds = bounds, &
         vals_input = vals_input(bounds%begc:bounds%endc), &
         vals_input_valid = vals_input_valid(bounds%begc:bounds%endc), &
         has_prognostic_state = has_prognostic_state(bounds%begc:bounds%endc), &
         var = var(bounds%begc:bounds%endc), &
         non_conserved_mass = non_conserved_mass_grc(bounds%begg:bounds%endg))

  end subroutine update_column_state_fill_using_fixed_values



  ! ========================================================================
  ! Private methods
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine update_column_state(this, bounds, &
       vals_input, vals_input_valid, has_prognostic_state, &
       var, non_conserved_mass)
    !
    ! !DESCRIPTION:
    ! Do the work of updating a column-level state variable due to changes in subgrid
    ! weights.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(column_state_updater_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds

    ! value used as input for each column, if that column is shrinking (can differ from
    ! var if we're doing special handling of that column)
    real(r8), intent(in) :: vals_input( bounds%begc: )

    ! whether each item in vals_input_valid is valid. An entry can be invalid if there is
    ! was no way to derive an input for that column.
    logical, intent(in) :: vals_input_valid( bounds%begc: )

    ! whether each column simulates the given variable (which, among other things,
    ! determines whether it can accept mass of this variable)
    logical, intent(in) :: has_prognostic_state( bounds%begc: )

    ! column-level variable of interest, updated in-place
    real(r8), intent(inout) :: var( bounds%begc: )

    ! mass lost (per unit of grid cell area) from each grid cell, if doing special
    ! handling that leads mass to not be conserved; this can happen due to growing columns
    ! where has_prognostic_state is false, or due to shrinking columns where
    ! has_prognostic_state is false but vals_input /= 0. Positive denotes mass lost from
    ! the grid cell, negative denotes mass gained by the grid cell.
    real(r8), intent(inout) :: non_conserved_mass( bounds%begg: )
    !
    ! !LOCAL VARIABLES:
    integer  :: c, g

    ! whether vals_input /= var in the columns where it should be equal
    logical :: bad_vals_input(bounds%begc:bounds%endc)

    ! fractional area lost from this column
    real(r8) :: area_lost

    ! area-weighted amount lost from a given column
    real(r8) :: area_weighted_loss

    ! area-weighted amount lost from decreasing weights, in each grid cell
    ! ((mass-per-unit-area) * (fractional area lost))
    real(r8) :: total_loss_grc(bounds%begg:bounds%endg)

    ! total area lost by columns decreasing in area (fractional area of grid cell)
    ! Note that this should exactly equal the total area gained by columns increasing in area
    real(r8) :: total_area_lost_grc(bounds%begg:bounds%endg)

    ! amount of state gain needed per unit area of area gained
    real(r8) :: gain_per_unit_area_grc(bounds%begg:bounds%endg)

    ! mass gained by a given column ((mass-gained-per-unit-area) * (fractional area gained))
    real(r8) :: mass_gained

    character(len=:), allocatable :: message

    character(len=*), parameter :: subname = 'update_column_state'
    !-----------------------------------------------------------------------

    ! ------------------------------------------------------------------------
    ! Error-checking on inputs
    ! ------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(var) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(vals_input) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(has_prognostic_state) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(non_conserved_mass) == (/bounds%endg/)), errMsg(__FILE__, __LINE__))

    ! For the sake of conservation - including the calculation of non_conserved_mass - we
    ! assume that vals_input == var wherever has_prognostic_state is .true. We ensure that
    ! is the case here.
    where(has_prognostic_state .and. vals_input_valid)
       bad_vals_input = (vals_input /= var)
    elsewhere
       ! where has_prognostic_state is false, vals_input can be anything
       bad_vals_input = .false.
    end where
    message = subname//': ERROR: where has_prognostic_state is true, vals_input must equal var'
    SHR_ASSERT_ALL(.not. bad_vals_input, message)

    ! ------------------------------------------------------------------------
    ! Begin main work
    ! ------------------------------------------------------------------------

    ! Determine the total mass loss for each grid cell, along with the gross area loss
    ! (which should match the gross area gain)
    total_loss_grc(bounds%begg:bounds%endg) = 0._r8
    total_area_lost_grc(bounds%begg:bounds%endg) = 0._r8
    do c = bounds%begc, bounds%endc
       g = col%gridcell(c)
       if (this%area_gained_col(c) < 0._r8) then
          if (.not. vals_input_valid(c)) then
             write(iulog,*) subname//' ERROR: shrinking column without valid input value'
             call endrun(decomp_index=c, clmlevel=namec, msg=errMsg(__FILE__, __LINE__))
          end if
          area_lost = -1._r8 * this%area_gained_col(c)
          total_area_lost_grc(g) = total_area_lost_grc(g) + area_lost
          area_weighted_loss = area_lost * vals_input(c)
          total_loss_grc(g) = total_loss_grc(g) + area_weighted_loss

          if (.not. has_prognostic_state(c)) then
             ! If a column doesn't model this state variable, then its vals_input value is
             ! really some fictitious quantity. So we track how much of this fictitious
             ! quantity we added to the system.
             non_conserved_mass(g) = non_conserved_mass(g) - area_weighted_loss
          end if
       end if
    end do

    ! Determine the mass loss per unit area for each grid cell. We essentially lump all of
    ! the loss together in a "loss" pool in each grid cell, so that we can then
    ! distribute that loss amongst the growing columns.
    do g = bounds%begg, bounds%endg
       if (total_area_lost_grc(g) > 0._r8) then
          gain_per_unit_area_grc(g) = total_loss_grc(g) / total_area_lost_grc(g)
       else
          gain_per_unit_area_grc(g) = 0._r8
       end if
    end do

    ! Distribute gain to growing columns
    do c = bounds%begc, bounds%endc
       g = col%gridcell(c)
       if (this%area_gained_col(c) > 0._r8) then
          mass_gained = this%area_gained_col(c) * gain_per_unit_area_grc(g)
          if (has_prognostic_state(c)) then
             var(c) = (this%cwtgcell_old(c) * var(c) + mass_gained) / this%cwtgcell_new(c)
          else
             non_conserved_mass(g) = non_conserved_mass(g) + mass_gained
          end if
       end if
    end do

  end subroutine update_column_state

end module dynColumnStateUpdaterMod
