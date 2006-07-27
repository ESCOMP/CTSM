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
  use domainMod   , only : domain_type
!
! !PUBLIC TYPES:
  implicit none
  save
!
! land model grid
!
  real(r8) :: spval = 1.e36                     ! special value

  type(domain_type) :: ldomain

! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

end module mkvarsur


