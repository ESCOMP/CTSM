module reweightMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Top level driver for things that happen when subgrid weights are changed. This is in
  ! a separate module from subgridWeightsMod in order to keep subgridWeightsMod lower-
  ! level - and particularly to break its dependency on filterMod.
  !
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod   , only : r8 => shr_kind_r8
  !
  ! PUBLIC TYPES:
  implicit none
  save

  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: reweight_wrapup               ! do modifications and error-checks after modifying subgrid weights

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine reweight_wrapup(bounds, glc_behavior)
    !
    ! !DESCRIPTION:
    ! Do additional modifications and error-checks that should be done after modifying subgrid
    ! weights
    !
    ! This should be called whenever any weights change (e.g., patch weights on the column,
    ! landunit weights on the grid cell, etc.).
    !
    ! !USES:
    use filterMod         , only : setFilters
    use subgridWeightsMod , only : set_active, check_weights
    use decompMod         , only : bounds_type, BOUNDS_LEVEL_CLUMP
    use glcBehaviorMod    , only : glc_behavior_type
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds                      ! clump bounds
    type(glc_behavior_type), intent(in) :: glc_behavior
    !------------------------------------------------------------------------

    SHR_ASSERT_FL(bounds%level == BOUNDS_LEVEL_CLUMP, sourcefile, __LINE__)

    call set_active(bounds, glc_behavior)
    call check_weights(bounds, active_only=.false.)
    call check_weights(bounds, active_only=.true.)
    call setFilters(bounds, glc_behavior)

  end subroutine reweight_wrapup

end module reweightMod
