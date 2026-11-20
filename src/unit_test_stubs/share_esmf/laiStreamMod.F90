module laiStreamMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Read LAI from stream
  !
  ! !USES:
  use decompMod        , only : bounds_type
  use abortutils       , only : endrun
  use clm_varctl       , only : iulog
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: lai_init    ! position datasets for LAI
  public :: lai_advance ! Advance the LAI streams (outside of a Open-MP threading loop)
  public :: lai_interp  ! interpolates between two years of LAI data (when LAI streams

  character(len=*), parameter :: sourcefile = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine lai_init(bounds)
    !
    ! Initialize data stream information for LAI.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds          ! bounds
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

  end subroutine lai_init

  !================================================================
  subroutine lai_advance( bounds )
    !
    ! Advance LAI streams
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

  end subroutine lai_advance

  !================================================================
  subroutine lai_interp(bounds, canopystate_inst)
    !
    ! Interpolate data stream information for Lai.
    !
    ! !USES:
    use CanopyStateType  , only : canopystate_type
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    type(canopystate_type) , intent(inout) :: canopystate_inst
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

  end subroutine lai_interp

end module LaiStreamMod
