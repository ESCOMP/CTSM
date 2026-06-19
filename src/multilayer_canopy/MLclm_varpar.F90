module MLclm_varpar

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing multilayer canopy model parameters
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  save

  ! Parameters for multilayer canopy

  integer, parameter :: nlevmlcan = 100     ! Number of layers in multilayer canopy model
  integer, parameter :: nleaf = 2           ! Number of leaf types (sunlit and shaded)
  integer, parameter :: isun = 1            ! Sunlit leaf index
  integer, parameter :: isha = 2            ! Shaded leaf index

end module MLclm_varpar
