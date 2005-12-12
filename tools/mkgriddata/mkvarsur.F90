module mkvarsur

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkvarsur
!
! !DESCRIPTION:
! Module containing 2-d surface boundary data information
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! land model grid
!
  real(r8) :: spval = 1.e36                     ! special value
  integer , allocatable :: numlon(:)            ! longitude points for each latitude strip
  real(r8), allocatable :: latixy(:,:)          ! latitude of grid cell (degrees)
  real(r8), allocatable :: longxy(:,:)          ! longitude of grid cell (degrees)
  real(r8), allocatable :: area(:,:)            ! grid cell area (km**2)
  real(r8), allocatable :: lats(:)              ! grid cell latitude, southern edge (degrees)
  real(r8), allocatable :: lonw(:,:)            ! grid cell longitude, western edge (degrees)
  integer , allocatable :: landmask(:,:)        ! land mask: 1 = land. 0 = ocean
  real(r8), allocatable :: landfrac(:,:)        ! fractional land
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

end module mkvarsur


