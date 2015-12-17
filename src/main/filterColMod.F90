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
  use decompMod    , only : bounds_type
  use GridcellType , only : grc
  use LandunitType , only : lun
  use clm_varcon   , only : ispval

  ! !PUBLIC TYPES:
  implicit none
  private
  save

  type, public :: filter_col_type
     integer :: num  ! number of points in the filter
     integer, allocatable :: indices(:)  ! column indices included in the filter
  end type filter_col_type

  ! !PUBLIC ROUTINES:

  ! create a filter from a gridcell-level logical array and an array of landunit type(s)
  ! of interest
  public :: col_filter_from_grcflags_ltypes

contains

  !-----------------------------------------------------------------------
  function col_filter_from_grcflags_ltypes(bounds, grcflags, ltypes) result(filter)
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
    !
    ! !LOCAL VARIABLES:
    integer :: g  ! gridcell index
    integer :: l  ! landunit index
    integer :: c  ! column index
    integer :: i  ! array index
    integer :: ltype  ! landunit type

    character(len=*), parameter :: subname = 'col_filter_from_grcflags_ltypes'
    !-----------------------------------------------------------------------

    filter%num = 0
    allocate(filter%indices(bounds%endc - bounds%begc + 1))

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
                filter%num = filter%num + 1
                filter%indices(filter%num) = c
             end do  ! c
          end do  ! i = 1, size(ltypes)
       end if  ! grcflags(g)
    end do  ! g

  end function col_filter_from_grcflags_ltypes

end module filterColMod
