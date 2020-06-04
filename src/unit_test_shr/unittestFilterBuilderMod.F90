module unittestFilterBuilderMod

  ! This module builds simple filters that can be used as inputs to routines that require
  ! a filter.

  implicit none
  private
  save

  public :: filter_empty  ! build a filter that doesn't contain any points

  public :: filter_from_range ! build a filter that includes all points between a start and end index

contains

  !-----------------------------------------------------------------------
  subroutine filter_empty(numf, filter)
    !
    ! !DESCRIPTION:
    ! Build a filter that doesn't contain any points
    !
    ! Allocates the 'filter' argument
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(out) :: numf ! number of points in the filter
    integer, allocatable, intent(out) :: filter(:)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'filter_empty'
    !-----------------------------------------------------------------------

    numf = 0
    allocate(filter(0))

  end subroutine filter_empty


  !-----------------------------------------------------------------------
  subroutine filter_from_range(start, end, numf, filter)
    !
    ! !DESCRIPTION:
    ! Build a filter that includes all points between a start and end index
    !
    ! Allocates the 'filter' argument
    !
    ! !ARGUMENTS:
    integer, intent(in) :: start ! start index to include
    integer, intent(in) :: end   ! end index to include
    integer, intent(out) :: numf ! number of points in the filter
    integer, allocatable, intent(out) :: filter(:)
    !
    ! !LOCAL VARIABLES:
    integer :: i

    character(len=*), parameter :: subname = 'filter_from_range'
    !-----------------------------------------------------------------------
    
    numf = end - start + 1
    numf = max(numf, 0)

    allocate(filter(numf))
    do i = 1, numf
       filter(i) = start + i - 1
    end do
    
  end subroutine filter_from_range

end module unittestFilterBuilderMod
