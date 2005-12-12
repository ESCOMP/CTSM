module mkvarpar

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varpar
!
! !DESCRIPTION:
! Module containing CLM parameters
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
  integer, parameter :: nlevsoi = 10  ! number of soil layers
  integer, parameter :: numpft  = 16  ! number of plant types
  integer, parameter :: noveg   = 0   ! value for non-vegetated pft
!
!EOP
!-----------------------------------------------------------------------

end module mkvarpar
