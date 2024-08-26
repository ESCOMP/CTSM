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
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: crop_heatstress_ndays       ! SdR: checks for number of days above Tcrit for crop heat stress


  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

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


    ! !LOCAL VARIABLES:
    real(r8) :: tcrit				! placeholder for for inputparameter
    real(r8) :: HS_ndays_min            ! minimum required number of heat stress days

    !----------------------------------------------------------------------

    tcrit = 306.0_r8
    HS_ndays_min = 3.0_r8

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
    else if (HS_ndays >= (HS_ndays_min + 0.2_r8) .and. &
                                                   croplive) then
         heatwave_crop = 1.0_r8 
    else
         heatwave_crop = 0.0_r8
    end if

  end subroutine crop_heatstress_ndays


end module CropHeatStress