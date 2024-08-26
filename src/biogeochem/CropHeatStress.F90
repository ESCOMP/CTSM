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
  use shr_infnan_mod    , only : isnan => shr_infnan_isnan, isinf => shr_infnan_isinf
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
  subroutine crop_heatstress_ndays(HS_ndays, heatwave_crop, t_veg_day)

    ! !DESCRIPTION:
    ! added by SdR for heat stress implementation
    ! function to keep track of critical temperature for crops that is exceeded at the end of each day
    ! needs a minimum of 3 consecutive days before heat stress is activated

    ! !ARGUMENTS:
    real(r8),        intent(inout)    :: HS_ndays         ! number of crop heat stress days (ndays) should be integer at final implementation
    real(r8),        intent(inout)    :: heatwave_crop         ! keep track if heatwave is activated
    real(r8),        intent(in)       :: t_veg_day


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
    if (isnan(HS_ndays)) then
         heatwave_crop = 4.0_r8
    else if (isinf(HS_ndays)) then
         heatwave_crop = 3.0_r8
    else if (HS_ndays >= 100) then
             heatwave_crop = 2.5_r8
    else if (HS_ndays >= 10) then
             heatwave_crop = 2.0_r8
    else if (HS_ndays >= (HS_ndays_min + 0.2_r8) .and. &
                                                   .not. isinf(HS_ndays)) then
         heatwave_crop = 1.0_r8 
    else
         heatwave_crop = 0.0_r8
    end if

  end subroutine crop_heatstress_ndays

  subroutine calc_HS_factor(HS_ndays, HS_factor, t_veg_day)

    ! !DESCRIPTION:
    ! function to calculate heat stress instensity by applying a factor to bglfr (increasing LAI decline)

    ! !ARGUMENTS:
    real(r8),        intent(inout)     :: HS_factor         ! keep track if heatwave is activated
    real(r8),        intent(in)        :: t_veg_day         ! daily vegetation temperature (Kelvin)
    real(r8),        intent(in)        :: HS_ndays          ! number of crop heat stress days (ndays) should be integer at final implementation


    ! !LOCAL VARIABLES:
    real(r8) :: tcrit, tmax            ! placeholders for input parameters (will be seperate development)
    integer  :: day_min, day_max

    !-----------------------------------------------------------------------

    tcrit = 300
    tmax  = 310
    day_min = 3
    day_max = 14

    !check  if stress occurs
    if (HS_ndays < day_min) then
       HS_factor = 1._r8
    else if (HS_ndays == day_min) then
       ! onset heatwave
         HS_factor = 1.05_r8
    else if (HS_ndays > day_min) then
       if (t_veg_day < tmax .and. t_veg_day > tcrit) then
            HS_factor = 1._r8 + 0.5_r8 * ((t_veg_day - tcrit) / (tmax - tcrit))
         else if (t_veg_day > tmax .and. HS_ndays < day_max) then
               HS_factor = 1.5_r8
         else if (t_veg_day > tmax .and. HS_ndays > day_max) then
               HS_factor = 1.7_r8
         end if
    end if

end subroutine calc_HS_factor 


end module CropHeatStress