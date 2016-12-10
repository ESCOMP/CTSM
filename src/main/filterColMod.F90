module filterColMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a type to hold column-level filters, along with factory methods to help create
  ! a column-level filter
  !
  ! To loop over the filter, use code like this:
  ! do fc = 1, myfilter%num
  !    c = myfilter%indices(fc)
  !    ...
  ! end do
  !
  ! !USES:
#include "shr_assert.h"
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use decompMod    , only : bounds_type
  use GridcellType , only : grc
  use LandunitType , only : lun
  use ColumnType   , only : col
  use clm_varcon   , only : ispval
  use clm_varctl   , only : iulog

  ! !PUBLIC TYPES:
  implicit none
  private
  save

  type, public :: filter_col_type
     integer :: num  ! number of points in the filter
     integer, allocatable :: indices(:)  ! column indices included in the filter
   contains
     procedure :: equals_filter
     generic :: operator(==) => equals_filter
  end type filter_col_type

  ! !PUBLIC ROUTINES:

  ! Create an empty filter
  public :: col_filter_empty

  ! Create a filter from an array of indices. This is mainly useful for unit testing.
  public :: col_filter_from_index_array

  ! Create a filter from a column-level logical array
  public :: col_filter_from_logical_array

  ! Create a filter from a column-level logical array, but including only active points
  public :: col_filter_from_logical_array_active_only

  ! Create a filter that contains one or more landunit type(s) of interest
  public :: col_filter_from_ltypes

  ! Create a filter from a landunit-level logical array
  public :: col_filter_from_lunflags

  ! Create a filter from a gridcell-level logical array and an array of landunit type(s)
  ! of interest
  public :: col_filter_from_grcflags_ltypes

  ! !PRIVATE ROUTINES:

  ! Whether a given column should be included in the filter based on the active flag
  private :: include_based_on_active

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  ! TODO(wjs, 2016-04-07) If repeated reallocation of the indices arrays (every time a
  ! filter is recreated - each time through the run loop) is a performance issue, then we
  ! could rewrite the creation functions to instead be subroutines that act on an existing
  ! filter object: I think this would involve replacing calls to col_filter_empty with
  ! something like filter%reset_filter; this would only allocate the indices array if it
  ! is not already allocated.

  !-----------------------------------------------------------------------
  function col_filter_empty(bounds) result(filter)
    !
    ! !DESCRIPTION:
    ! Initialize a filter object
    !
    ! !ARGUMENTS:
    type(filter_col_type) :: filter  ! function result
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'col_filter_empty'
    !-----------------------------------------------------------------------

    filter%num = 0
    allocate(filter%indices(bounds%endc - bounds%begc + 1))

  end function col_filter_empty

  !-----------------------------------------------------------------------
  function col_filter_from_index_array(bounds, indices_col) result(filter)
    !
    ! !DESCRIPTION:
    ! Create a filter from an array of indices.
    !
    ! This is mainly useful for unit testing.
    !
    ! !ARGUMENTS:
    type(filter_col_type) :: filter  ! function result
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: indices_col(:)  ! column-level array of indices to include in filter
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'col_filter_from_index_array'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL(indices_col >= bounds%begc, errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL(indices_col <= bounds%endc, errMsg(sourcefile, __LINE__))

    filter = col_filter_empty(bounds)

    filter%num = size(indices_col)
    filter%indices(1:filter%num) = indices_col

  end function col_filter_from_index_array


  !-----------------------------------------------------------------------
  function col_filter_from_logical_array(bounds, logical_col) result(filter)
    !
    ! !DESCRIPTION:
    ! Create a column-level filter from a column-level logical array.
    !
    ! This version does not consider whether a column is active: it simply includes any
    ! column 'c' for which logical_col(c) is .true.
    !
    ! !ARGUMENTS:
    type(filter_col_type) :: filter  ! function result
    type(bounds_type), intent(in) :: bounds
    logical, intent(in) :: logical_col(bounds%begc:) ! column-level logical array
    !
    ! !LOCAL VARIABLES:
    integer :: c

    character(len=*), parameter :: subname = 'col_filter_from_logical_array'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(logical_col) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    filter = col_filter_empty(bounds)

    do c = bounds%begc, bounds%endc
       if (logical_col(c)) then
          filter%num = filter%num + 1
          filter%indices(filter%num) = c
       end if
    end do

  end function col_filter_from_logical_array

  !-----------------------------------------------------------------------
  function col_filter_from_logical_array_active_only(bounds, logical_col) result(filter)
    !
    ! !DESCRIPTION:
    ! Create a column-level filter from a column-level logical array. Only include active
    ! points in the filter: even if the logical array is true for a given column, that
    ! column is excluded if it is inactive.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(filter_col_type) :: filter  ! function result
    type(bounds_type), intent(in) :: bounds
    logical, intent(in) :: logical_col(bounds%begc:) ! column-level logical array
    !
    ! !LOCAL VARIABLES:
    integer :: c

    character(len=*), parameter :: subname = 'col_filter_from_logical_array_active_only'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(logical_col) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    filter = col_filter_empty(bounds)

    do c = bounds%begc, bounds%endc
       if (col%active(c)) then
          if (logical_col(c)) then
             filter%num = filter%num + 1
             filter%indices(filter%num) = c
          end if
       end if
    end do

  end function col_filter_from_logical_array_active_only


  !-----------------------------------------------------------------------
  function col_filter_from_ltypes(bounds, ltypes, include_inactive) &
       result(filter)
    !
    ! !DESCRIPTION:
    ! Create a column-level filter that includes one or more landunit type(s) of interest
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(filter_col_type) :: filter  ! function result
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: ltypes(:)  ! landunit type(s) of interest
    logical, intent(in) :: include_inactive ! whether inactive points should be included in the filter
    !
    ! !LOCAL VARIABLES:
    integer :: c
    integer :: l

    character(len=*), parameter :: subname = 'col_filter_from_ltypes'
    !-----------------------------------------------------------------------

    filter = col_filter_empty(bounds)

    do c = bounds%begc, bounds%endc
       if (include_based_on_active(c, include_inactive)) then
          l = col%landunit(c)
          if (any(ltypes(:) == lun%itype(l))) then
             filter%num = filter%num + 1
             filter%indices(filter%num) = c
          end if
       end if
    end do

  end function col_filter_from_ltypes

  !-----------------------------------------------------------------------
  function col_filter_from_lunflags(bounds, lunflags, include_inactive) &
       result(filter)
    !
    ! !DESCRIPTION:
    ! Create a column-level filter from a landunit-level logical array.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(filter_col_type) :: filter  ! function result
    type(bounds_type), intent(in) :: bounds
    logical, intent(in) :: lunflags(bounds%begl:)  ! landunit-level logical array
    logical, intent(in) :: include_inactive ! whether inactive points should be included in the filter
    !
    ! !LOCAL VARIABLES:
    integer :: c
    integer :: l

    character(len=*), parameter :: subname = 'col_filter_from_lunflags'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(lunflags) == (/bounds%endl/)), errMsg(sourcefile, __LINE__))

    filter = col_filter_empty(bounds)

    do c = bounds%begc, bounds%endc
       if (include_based_on_active(c, include_inactive)) then
          l = col%landunit(c)
          if (lunflags(l)) then
             filter%num = filter%num + 1
             filter%indices(filter%num) = c
          end if
       end if
    end do

  end function col_filter_from_lunflags


  !-----------------------------------------------------------------------
  function col_filter_from_grcflags_ltypes(bounds, grcflags, ltypes, include_inactive) &
       result(filter)
    !
    ! !DESCRIPTION:
    ! Create a column-level filter from a gridcell-level logical array and an array of
    ! landunit type(s) of interest. The filter will contain all columns for which (a)
    ! grcflags is true for the gridcell containing this column, and (b) the landunit type
    ! for the landunit containing this column is one of the types in ltypes.
    !
    ! !ARGUMENTS:
    type(filter_col_type) :: filter  ! function result
    type(bounds_type), intent(in) :: bounds
    logical, intent(in) :: grcflags(bounds%begg:)  ! gridcell-level logical array
    integer, intent(in) :: ltypes(:)  ! landunit type(s) of interest
    logical, intent(in) :: include_inactive ! whether inactive points should be included in the filter
    !
    ! !LOCAL VARIABLES:
    integer :: g  ! gridcell index
    integer :: l  ! landunit index
    integer :: c  ! column index
    integer :: i  ! array index
    integer :: ltype  ! landunit type

    character(len=*), parameter :: subname = 'col_filter_from_grcflags_ltypes'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(grcflags) == (/bounds%endg/)), errMsg(sourcefile, __LINE__))

    filter = col_filter_empty(bounds)
    
    ! This loops over g then l then c rather than just looping over all columns, because
    ! this is likely more efficient for sparse filters (e.g., sparse grcflags or uncommon
    ! ltypes).
    do g = bounds%begg, bounds%endg
       if (grcflags(g)) then
          do i = 1, size(ltypes)
             ltype = ltypes(i)
             l = grc%landunit_indices(ltype, g)
             if (l == ispval) then
                cycle
             end if

             do c = lun%coli(l), lun%colf(l)
                if (include_based_on_active(c, include_inactive)) then
                   filter%num = filter%num + 1
                   filter%indices(filter%num) = c
                end if
             end do  ! c
          end do  ! i = 1, size(ltypes)
       end if  ! grcflags(g)
    end do  ! g

  end function col_filter_from_grcflags_ltypes

  !-----------------------------------------------------------------------
  pure function include_based_on_active(c, include_inactive) result(include_point)
    !
    ! !DESCRIPTION:
    ! Returns true if the given column should be included in a filter based on its active
    ! flag
    !
    ! !ARGUMENTS:
    logical :: include_point  ! function result
    integer, intent(in) :: c  ! column index
    logical, intent(in) :: include_inactive ! whether inactive points are included in this filter
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'include_based_on_active'
    !-----------------------------------------------------------------------

    ! This code is written to avoid the check of col%active if include_inactive is true.
    ! This is needed in the case of filters that are created in initialization, before
    ! the active flags are set.
    if (include_inactive) then
       include_point = .true.
    else if (col%active(c)) then
       include_point = .true.
    else
       include_point = .false.
    end if

  end function include_based_on_active


  !-----------------------------------------------------------------------
  function equals_filter(this, other) result(equal)
    !
    ! !DESCRIPTION:
    ! Returns true if the two filters are equal.
    !
    ! If they differ, prints some information about how they differ.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    logical :: equal  ! function result
    class(filter_col_type), intent(in) :: this
    class(filter_col_type), intent(in) :: other
    !
    ! !LOCAL VARIABLES:
    integer :: i

    character(len=*), parameter :: subname = 'equals_filter'
    !-----------------------------------------------------------------------

    equal = .true.

    if (this%num /= other%num) then
       equal = .false.
       write(iulog,*) ' '
       write(iulog,'(a, i0, a, i0)') 'equals_filter false: Sizes differ: ', &
            this%num, ' /= ', other%num
    else
       do i = 1, this%num
          if (this%indices(i) /= other%indices(i)) then
             equal = .false.
             write(iulog,*) ' '
             write(iulog,'(a, i0, a, i0, a, i0)') &
                  'equals_filter false: Values differ; first difference at ', &
                  i, ': ', this%indices(i), ' /= ', other%indices(i)
             exit
          end if
       end do
    end if

  end function equals_filter


end module filterColMod
