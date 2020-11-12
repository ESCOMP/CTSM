module point_of_interest

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module contains logical functions to find point(s) of interest in CTSM. Note the
  ! abbreviation: p.o.i = point of interest.
  !
  ! The idea is: It is common, when debugging, to want to print the values of various
  ! variables for all patches or columns of certain landunit types within a certain grid
  ! cell of interest. This module helps create the logical functions needed to do this.
  !
  ! This module is compiled into every CTSM build, but is not invoked by default. To use
  ! it:
  !
  ! (1) Customize the code here by changing the functions below and/or adding new
  !     functions. Look for comments about "customize" to see what to customize.
  !
  ! (2) Add calls in the code
  !
  !     Its typical use will be something like:
  !
  !       do fc = 1, num_nolakec
  !          c = filter_nolakec(fc)
  !
  !          ! Various code here, maybe setting foo and bar variables
  !
  !          if (poi_c(c)) then
  !             write(iulog,*) 'DEBUG: foo, bar = ', foo(c), bar(c)
  !          end if
  !       end do
  !
  ! !USES:
  use GridcellType, only : grc
  use LandunitType, only : lun
  use ColumnType, only : col
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use landunit_varcon, only : istsoil

  implicit none
  save
  private

  ! Customize: can define other levels like poi_p, etc.
  public :: poi_c

contains
  
  !-----------------------------------------------------------------------
  logical function poi_c(c)
    ! This function can be used in a column-level loop to find columns with a given
    ! landunit type(s) within the grid cell(s) of interest.
    integer, intent(in) :: c

    integer :: g, l
    !-----------------------------------------------------------------------

    g = col%gridcell(c)
    l = col%landunit(c)

    poi_c = .false.

    ! Customize this conditional; it is currently set up to flag columns in the natural
    ! vegetated landunit (istsoil) of the target grid cell(s).
    if (at_poi(g) .and. lun%itype(l) == istsoil) then
       poi_c = .true.
    end if

  end function poi_c

  !-----------------------------------------------------------------------
  logical function at_poi(g)
    integer, intent(in) :: g

    ! Customize these parameters (adding more blocks if necessary). These give the
    ! longitude and latitude of the grid cell of interest.
    real(r8), parameter :: poi_lon = 237.5_r8
    real(r8), parameter :: poi_lat = -72.94737_r8

    real(r8), parameter :: poi_tol = 0.01_r8      ! tolerance on check of lat/lon
    !-----------------------------------------------------------------------
    
    if ( abs(grc%londeg(g) - poi_lon) < poi_tol .and. &
         abs(grc%latdeg(g) - poi_lat) < poi_tol) then
       at_poi = .true.
    else
       at_poi = .false.
    end if

  end function at_poi

end module point_of_interest
