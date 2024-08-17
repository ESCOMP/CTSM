module mkpftConstantsMod

  !-----------------------------------------------------------------------
  ! Constants used by mkpft and related code
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  ! public data members: 
  integer, parameter, public :: maxpft = 78   ! maximum # of PFT

  integer, public    :: num_natpft = -1       ! number of PFTs on the natural vegetation
                                              ! landunit, NOT including bare ground
                                              ! (includes generic crops for runs with
                                              ! create_crop_landunit=false)

  integer, public    :: num_cft               ! number of CFTs on the crop landunit
  integer, public    :: natpft_lb             ! lower bound for natural pft arrays
  integer, public    :: natpft_ub             ! upper bound for natural pft arrays
  integer, public    :: cft_lb                ! lower bound for cft arrays
  integer, public    :: cft_ub                ! upper bound for cft arrays
  
  integer, parameter, public :: baregroundindex = 0  ! index of bare ground in a natural pft array
  
  ! The following is NOT set as a parameter so that it can be overridden in unit tests
  integer, public :: c3cropindex = 15
  integer, public :: c3irrcropindex = 16

end module mkpftConstantsMod
