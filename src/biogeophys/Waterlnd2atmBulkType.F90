module Waterlnd2atmBulkType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type containing water lnd2atm variables that just apply to bulk
  ! water. Note that this type extends the base waterlnd2atm_type, so the full
  ! waterlnd2atmbulk_type contains the union of the fields defined here and the fields
  ! defined in waterlnd2atm_type.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use clm_varpar     , only : nlevgrnd
  use WaterLnd2atmType , only : waterlnd2atm_type
  use WaterInfoBaseType, only : water_info_base_type
  use WaterTracerContainerType, only : water_tracer_container_type
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, extends(waterlnd2atm_type), public :: waterlnd2atmbulk_type

     real(r8), pointer :: h2osoi_vol_grc     (:,:) ! volumetric soil water (0~watsat, m3/m3, nlevgrnd) (for dust model)
   contains
     procedure, public  :: InitBulk
     procedure, private :: InitBulkAllocate
     procedure, private :: InitBulkHistory
     procedure, private :: InitBulkCold

  end type waterlnd2atmbulk_type


  character(len=*), parameter, private :: sourcefile = &
       __FILE__
 !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine InitBulk(this, bounds, info, vars)

    class(waterlnd2atmbulk_type), intent(inout) :: this
    type(bounds_type) , intent(in) :: bounds
    class(water_info_base_type), intent(in), target :: info
    type(water_tracer_container_type), intent(inout) :: vars


    call this%Init(bounds, info, vars)

    call this%InitBulkAllocate(bounds)

    call this%InitBulkHistory(bounds)

    call this%InitBulkCold(bounds)

  end subroutine InitBulk

  !------------------------------------------------------------------------
  subroutine InitBulkAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(waterlnd2atmbulk_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ival  = 0.0_r8  ! initial value
    integer :: begg, endg
    !------------------------------------------------------------------------

    begg = bounds%begg; endg= bounds%endg

    allocate(this%h2osoi_vol_grc     (begg:endg,1:nlevgrnd)) ; this%h2osoi_vol_grc     (:,:) = ival

  end subroutine InitBulkAllocate

  !------------------------------------------------------------------------
  subroutine InitBulkHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize history vars
    !
    ! !ARGUMENTS:
    class(waterlnd2atmbulk_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    !------------------------------------------------------------------------

    ! Nothing to do for now

  end subroutine InitBulkHistory

  !-----------------------------------------------------------------------
  subroutine InitBulkCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize cold start conditions
    !
    ! !ARGUMENTS:
    class(waterlnd2atmbulk_type), intent(inout) :: this
    type(bounds_type)     , intent(in)          :: bounds
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

    ! Nothing to do for now

  end subroutine InitBulkCold

end module Waterlnd2atmBulkType
