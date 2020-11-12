! This module contains logical functions to find point(s) of interest in CLM.
! Note the abbreviation: p.o.i = point of interest

! Its typical use will be something like:
!    if (poi_c(c)) then
!       write(iulog,*) ...

! Look for comments about "customize" to see what to customize.

module point_of_interest

  use GridcellType, only : grc
  use LandunitType, only : lun
  use ColumnType, only : col
  use shr_kind_mod  , only : r8 => shr_kind_r8
  
  implicit none
  save
  private

  ! Customize: can define other levels like poi_p, etc.
  public :: poi_c

contains
  
  logical function poi_c(c)
    integer, intent(in) :: c

    integer :: g, l

    g = col%gridcell(c)
    l = col%landunit(c)

    poi_c = .false.

    ! Customize this conditional
    if (at_poi(g) .and. lun%itype(l) == 1) then
       poi_c = .true.
    end if

  end function poi_c



  logical function at_poi(g)
    integer, intent(in) :: g

    ! Customize these parameters (adding more blocks if necessary)
    real(r8), parameter :: poi_lon = 237.5_r8
    real(r8), parameter :: poi_lat = -72.94737_r8

    real(r8), parameter :: poi_tol = 0.01_r8      ! tolerance on check of lat/lon

    
    if ( abs(grc%londeg(g) - poi_lon) < poi_tol .and. &
         abs(grc%latdeg(g) - poi_lat) < poi_tol) then
       at_poi = .true.
    else
       at_poi = .false.
    end if

  end function at_poi

end module point_of_interest
