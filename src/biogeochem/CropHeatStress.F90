module CropHeatStress

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !MODULE: CNPhenologyMod
  !
  ! !DESCRIPTION:
  ! Module for implementing the effect of heat stress on crop production
  ! Added by SdR. 
  !
  ! !USES:
  !------------------------------------------------------------------------
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use shr_sys_mod       , only : shr_sys_flush
  use shr_infnan_mod    , only : isnan => shr_infnan_isnan, isinf => shr_infnan_isinf, nan => shr_infnan_nan, assignment(=)
  use clm_varcon, only : spval
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: crop_heatstress_ndays       ! SdR: checks for number of days above Tcrit for crop heat stress
  !
  ! !PUBLIC FOR UNIT TESTING
  real(r8), public, parameter :: tcrit = 306._r8
  real(r8), public, parameter :: HS_ndays_min = 3._r8

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  subroutine crop_heatstress_reset(HS_ndays, heatwave_crop)

    ! !DESCRIPTION
    ! Unsets variables related to crop heat stress

    ! !ARGUMENTS:
    real(r8), intent(inout)    :: HS_ndays ! number of crop heat stress days (ndays) should be integer at final implementation
    real(r8), intent(inout)    :: heatwave_crop ! keep track if heatwave is activated

    HS_ndays = 0._r8
    heatwave_crop = 0._r8

  end subroutine crop_heatstress_reset

  !------------------------------------------------------------------------
  subroutine crop_heatstress_ndays(HS_ndays, heatwave_crop, t_veg_day,croplive)

    ! !DESCRIPTION:
    ! added by SdR for heat stress implementation
    ! function to keep track of critical temperature for crops that is exceeded at the end of each day
    ! needs a minimum of 3 consecutive days before heat stress is activated

    ! !ARGUMENTS:
    real(r8),        intent(inout)    :: HS_ndays              ! number of crop heat stress days (ndays) should be integer at final implementation
    real(r8),        intent(inout)    :: heatwave_crop         ! keep track if heatwave is activated
    real(r8),        intent(in)       :: t_veg_day
    logical,        intent(in)        :: croplive              !check if crop is alive

    !----------------------------------------------------------------------

    ! No heat stress if crop isn't alive
    if (.not. croplive) then
       call crop_heatstress_reset(HS_ndays, heatwave_crop)
       return
    end if

    ! Don't do anything if it's not a real temperature
    if (t_veg_day > spval / 1000._r8) then
       return
    end if

    ! check if tcrit is exceeded and count days
    if (t_veg_day >= tcrit) then
         HS_ndays = HS_ndays + 1.0_r8
    else
         HS_ndays = 0.0_r8
    end if

    ! check if a heatwave is occurring
    if (isinf(HS_ndays)) then
        heatwave_crop = 4.0_r8
    else if (HS_ndays > 1000000000000) then
        heatwave_crop = nan
    else if (HS_ndays >= (HS_ndays_min - 0.2_r8)) then
         heatwave_crop = 1.0_r8 
    else
         heatwave_crop = 0.0_r8
    end if

  end subroutine crop_heatstress_ndays


end module CropHeatStress