module HillslopeHydrologyMod

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

  implicit none
  private   
  save

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: hcol_distance
  public :: hcol_width
  public :: hcol_area
  public :: hcol_elevation
  public :: hcol_slope


  ! !PRIVATE MEMBER FUNCTIONS:
!  private :: 

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  function hcol_distance(c, ctop, cbottom, hill_length) result(x)
    !
    ! !DESCRIPTION:
    ! Returns distance from the bottom of the hillslope of the 
    ! column's node.  Assumes hilltop to hillbottom column 
    ! ordering based on lun%coli, lun%colf (see initHillslopeMod).
    !
    ! !USES:
    !
    ! !ARGUMENTS:
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

  function hcol_width(x,alpha,beta,hill_length,hill_width,hill_height) result(width)
    !
    ! !DESCRIPTION:
    ! Returns width of hillslope column.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8) :: width                   ! function result
    real(r8), intent(in) :: x           ! distance along hillslope
    real(r8), intent(in) :: alpha       ! profile curvature parameter
    real(r8), intent(in) :: beta        ! plan curvature parameter
    real(r8), intent(in) :: hill_length ! total hillslope length
    real(r8), intent(in) :: hill_width  ! total hillslope width
    real(r8), intent(in) :: hill_height ! total hillslope height
    !
    ! !LOCAL VARIABLES:
    real(r8) :: eps = 1.e-6_r8
    real(r8) :: x0,y0,yl                ! integration limits
    character(len=*), parameter :: subname = 'hcol_width'
    !-----------------------------------------------------------------------

    ! width function has special case for n = 2
    ! in this implementation, integration limits depend on sign of beta
    if (abs(alpha - 2._r8) < eps) then

       ! function blows up for x0=0; integration limits set by trial and error
       if(beta < 0._r8) then
          y0 = hill_width/2._r8
          yl = 0.1_r8
          x0=hill_length *(yl/y0)**(-hill_height/(beta*hill_length**2))
       else
          x0 = 0.2_r8
          y0 = (hill_width/2._r8)&
               *(x0/hill_length)**(beta*hill_length**2/hill_height)
       endif

       ! compiler does not like log(zero)
       if (x == 0._r8) then
          if (beta <  0._r8) then
             width = hill_width/2._r8
          else
             width = eps
          endif
       else
          width = y0*(x/x0)**(beta*hill_length**2/hill_height)
       endif
    else 
       ! alpha /= 2 case, x0 equals zero
       y0 = hill_width/2._r8
       if(beta > 0._r8) then
          y0 = y0 * exp(-(2._r8*beta*hill_length**2) &
               / (hill_height*(2._r8 - alpha)*alpha))
       endif
       ! compiler does not like zero to a negative power.
       if (x == 0._r8) then
          width = y0
       else
          width = y0*exp((((2._r8*beta*hill_length**2) &
               / (hill_height*(2._r8 - alpha)*alpha)) &
               * (x/hill_length)**(2._r8-alpha)))
       endif
    endif
    ! hillslope width is twice integral [0:x]
    width = 2._r8 * width 

  end function hcol_width

  function hcol_area(xtop,xbottom,alpha,beta,hill_length,hill_width,hill_height) result(area)
    !
    ! !DESCRIPTION:
    ! Returns area of a hillslope column.  Area is calculated by 
    ! numerically integrating using hcol_width function.
    ! 
    !
    ! !USES:
    !
    ! !ARGUMENTS:
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
       area = area + dx * hcol_width(x,alpha,beta,hill_length,hill_width,hill_height)
    enddo

  end function hcol_area
             
  function hcol_elevation(xtop,xbottom,alpha,beta,hill_length,hill_width,hill_height) result(elev)
    !
    ! !DESCRIPTION:
    ! Returns mean elevation of column (relative to hillslope bottom).
    ! Area-weighted mean elevation is calculated by 
    ! numerically integrating using hcol_width function.    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8) :: elev                    ! function result
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
    real(r8) :: x, y
    real(r8) :: dx
    real(r8) :: dA
    real(r8) :: area
    character(len=*), parameter :: subname = 'hcol_elevation'
    !-----------------------------------------------------------------------
    ! mean elevation of column relative to hillslope bottom
    ! elevation is first integrated analytically in across-slope direction
    ! then summed in along-slope direction

    dx = (xtop - xbottom)/real(ndiv)

    elev = 0._r8
    area = 0._r8
    do n = 0, ndiv-1
       x = xbottom + (n+0.5)*dx
       y = hcol_width(x,alpha,beta,hill_length,hill_width,hill_height)
       dA = dx * y
       area = area + dA
       elev = elev + dx * (hill_height*y*(x/hill_length)**alpha &
            + beta*y**3/12._r8)
    enddo
    elev = elev / area

  end function hcol_elevation

   function hcol_slope(xtop,xbottom,alpha, hill_length, hill_height) result(slope)
    !
    ! !DESCRIPTION:
    ! Returns mean along-hillslope slope of hillslope column 
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8) :: slope                   ! function result
    real(r8), intent(in) :: xtop        ! distance to upper edge of column
    real(r8), intent(in) :: xbottom     ! distance to lower edge of column
    real(r8), intent(in) :: alpha       ! hillslope profile curvature parameter
    real(r8), intent(in) :: hill_length ! total hillslope length
    real(r8), intent(in) :: hill_height ! total hillslope height
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'hcol_slope'
    !-----------------------------------------------------------------------

    ! mean along-hill slope of column
    slope = hill_height &
         * ((xtop/hill_length)**alpha &
         - (xbottom/hill_length)**alpha) &
         / (xtop - xbottom)

  end function hcol_slope

end module HillslopeHydrologyMod
