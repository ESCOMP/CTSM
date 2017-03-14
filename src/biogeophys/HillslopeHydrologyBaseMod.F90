module HillslopeHydrologyBaseMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate geomorphological quantities for hillslope columns.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use spmdMod        , only : masterproc
  use abortutils     , only : endrun
  use clm_varctl     , only : iulog
  use decompMod      , only : bounds_type

  implicit none
  private   
  save

  !-----------------------------------------------------------------------

  ! !PUBLIC TYPES:
  type, abstract, public :: hillslope_geomorphology_type
     private
     ! variable declarations

   contains
     ! procedure declarations
     procedure (Init_interface)          , deferred :: Init
     procedure (hcol_width_interface)    , deferred :: hcol_width
     procedure (hcol_elevation_interface), deferred :: hcol_elevation
     procedure (hcol_slope_interface)    , deferred :: hcol_slope
     procedure :: hcol_distance
     procedure :: hcol_area

  end type hillslope_geomorphology_type

  !-----------------------------------------------------------------------

  abstract interface

     subroutine Init_interface(this,bounds,fsurdat)
       !
       ! !DESCRIPTION:
       ! Initialize hillslope geomorphology
       !
       ! !USES:
       
       use decompMod,   only : bounds_type
       import hillslope_geomorphology_type
              
       !
       ! !ARGUMENTS:
       class(hillslope_geomorphology_type), intent(inout) :: this
       type(bounds_type), intent(in)    :: bounds
       character(len=*) , intent(in)    :: fsurdat   ! surface data file name
     end subroutine Init_interface
     
     function hcol_width_interface(this,x,alpha,beta,hill_length,hill_width,hill_height) result(width)
       !
       ! !DESCRIPTION:
       ! Returns width of hillslope column.
       !
       ! !USES:
       use shr_kind_mod   , only : r8 => shr_kind_r8
       import hillslope_geomorphology_type
       !
       ! !ARGUMENTS:
       class(hillslope_geomorphology_type) , intent(in) :: this
       real(r8) :: width                   ! function result
       real(r8), intent(in) :: x           ! distance along hillslope
       real(r8), intent(in) :: alpha       ! profile curvature parameter
       real(r8), intent(in) :: beta        ! plan curvature parameter
       real(r8), intent(in) :: hill_length ! total hillslope length
       real(r8), intent(in) :: hill_width  ! total hillslope width
       real(r8), intent(in) :: hill_height ! total hillslope height

     end function hcol_width_interface

     !-----------------------------------------------------------------------
     function hcol_slope_interface(this,xtop,xbottom,alpha, hill_length, hill_height) result(slope)
       !
       ! !DESCRIPTION:
       ! Returns mean along-hillslope slope of hillslope column 
       !
       ! !USES:
       use shr_kind_mod   , only : r8 => shr_kind_r8
       import hillslope_geomorphology_type
       !
       ! !ARGUMENTS:
       class(hillslope_geomorphology_type) , intent(in) :: this
       real(r8) :: slope                   ! function result
       real(r8), intent(in) :: xtop        ! distance to upper edge of column
       real(r8), intent(in) :: xbottom     ! distance to lower edge of column
       real(r8), intent(in) :: alpha       ! hillslope profile curvature parameter
       real(r8), intent(in) :: hill_length ! total hillslope length
       real(r8), intent(in) :: hill_height ! total hillslope height
     end function hcol_slope_interface
     
     !-----------------------------------------------------------------------

     function hcol_elevation_interface(this,xtop,xbottom,alpha,beta,hill_length,hill_width,hill_height) result(elev)
       !
       ! !DESCRIPTION:
       ! Returns mean elevation of column (relative to hillslope bottom).
       ! Area-weighted mean elevation is calculated by 
       ! numerically integrating using hcol_width function.    !
       ! !USES:
       use shr_kind_mod   , only : r8 => shr_kind_r8
       import hillslope_geomorphology_type
       !
       ! !ARGUMENTS:
       class(hillslope_geomorphology_type) , intent(in) :: this
       real(r8) :: elev                    ! function result
       real(r8), intent(in) :: xtop        ! upper integration limit
       real(r8), intent(in) :: xbottom     ! lower integration limit
       real(r8), intent(in) :: alpha       ! profile curvature parameter
       real(r8), intent(in) :: beta        ! plan curvature parameter
       real(r8), intent(in) :: hill_length ! total hillslope length
       real(r8), intent(in) :: hill_width  ! total hillslope width
       real(r8), intent(in) :: hill_height ! total hillslope height

     end function hcol_elevation_interface

  end interface

! PRIVATE 
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  function hcol_distance(this, c, ctop, cbottom, hill_length) result(x)
    !
    ! !DESCRIPTION:
    ! Returns distance from the bottom of the hillslope of the 
    ! column's node.  Assumes hilltop to hillbottom column 
    ! ordering based on lun%coli, lun%colf (see initHillslopeMod).
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(hillslope_geomorphology_type) , intent(in) :: this
    real(r8) :: x                       ! function result
    integer, intent(in) :: c            ! current column
    integer, intent(in) :: ctop         ! hillslope top column
    integer, intent(in) :: cbottom      ! hillslope bottom column
    real(r8), intent(in) :: hill_length ! total hillslope length
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'hcol_distance'
    !-----------------------------------------------------------------------

    x = hill_length * (real(cbottom - c,r8) +0.5_r8) &
         / real(cbottom - ctop + 1,r8) 

  end function hcol_distance

  !-----------------------------------------------------------------------
  function hcol_area(this,xtop,xbottom,alpha,beta,hill_length,hill_width,hill_height) result(area)
    !
    ! !DESCRIPTION:
    ! Returns area of a hillslope column.  Area is calculated by 
    ! numerically integrating using hcol_width function.
    ! 
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(hillslope_geomorphology_type) , intent(in) :: this
    real(r8) :: area                    ! function result
    real(r8), intent(in) :: xtop        ! upper integration limit
    real(r8), intent(in) :: xbottom     ! lower integration limit
    real(r8), intent(in) :: alpha       ! profile curvature parameter
    real(r8), intent(in) :: beta        ! plan curvature parameter
    real(r8), intent(in) :: hill_length ! total hillslope length
    real(r8), intent(in) :: hill_width  ! total hillslope width
    real(r8), intent(in) :: hill_height ! total hillslope height
    !
    ! !LOCAL VARIABLES:
    integer  :: n
    integer, parameter  :: ndiv = 100
    real(r8) :: x
    real(r8) :: dx
    character(len=*), parameter :: subname = 'hcol_area'
    !-----------------------------------------------------------------------
    ! surface area of column
    dx = (xtop - xbottom)/real(ndiv)

    area = 0._r8
    do n = 0, ndiv-1
       x = xbottom + (n+0.5)*dx
       area = area + dx * this%hcol_width(x,alpha,beta,hill_length,hill_width,hill_height)
    enddo

  end function hcol_area
             
end module HillslopeHydrologyBaseMod
