#include <misc.h>
#include <preproc.h>

module clm_varsur

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varsur
!
! !DESCRIPTION:
! Module containing 2-d surface boundary data information
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varpar  , only : lsmlon, lsmlat
!
! !PUBLIC TYPES:
  implicit none
  save
!
! land model grid
!
  integer  :: numlon(lsmlat)                 ! longitude points for each latitude strip
  integer  :: landmask(lsmlon,lsmlat)        ! land mask: 1 = land. 0 = ocean
  real(r8) :: landfrac(lsmlon,lsmlat)        ! fractional land
  real(r8) :: latixy(lsmlon,lsmlat)          ! latitude of grid cell (degrees)
  real(r8) :: longxy(lsmlon,lsmlat)          ! longitude of grid cell (degrees)
  real(r8) :: area(lsmlon,lsmlat)            ! grid cell area (km**2)
  real(r8) :: lats(lsmlat+1)                 ! grid cell latitude, southern edge (degrees)
  real(r8) :: lonw(lsmlon+1,lsmlat)          ! grid cell longitude, western edge (degrees)
  real(r8) :: lsmedge(4)                     ! North,East,South,West edges of grid (deg)
  logical  :: fullgrid  = .true.             ! true => no grid reduction towards poles
!
! surface boundary data
!
  real(r8) :: pctspec(lsmlon,lsmlat)         ! percent of special landunits wrt gridcell
  logical  :: all_pfts_on_srfdat = .false.   ! true=>old format dataset 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

end module clm_varsur
