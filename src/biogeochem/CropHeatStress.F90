module CropHeatStress

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !MODULE: CropHeatStress
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
  public  :: calc_HS_factor              ! SdR: calculates heat stress magnitude affecting leaf area decline (grainfill phase)
  public  :: crop_heatstress_reset 

  !
  ! !PUBLIC FOR UNIT TESTING
  real(r8), public, parameter :: tcrit = 302.15_r8
  real(r8), public, parameter :: tmax = 318.15_r8
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
  subroutine crop_heatstress_ndays(HS_ndays, heatwave_crop, t_veg_day, croplive)

    ! !DESCRIPTION:
    ! added by SdR for heat stress implementation
    ! function to keep track of critical temperature for crops that is exceeded at the end of each day
    ! needs a minimum of 3 consecutive days before heat stress is activated

    ! !ARGUMENTS:
    real(r8),        intent(inout)    :: HS_ndays              ! number of crop heat stress days (ndays) should be integer at final implementation
    real(r8),        intent(inout)    :: heatwave_crop         ! keep track if heatwave is activated
    real(r8),        intent(in)       :: t_veg_day
    logical,         intent(in)       :: croplive              ! crop between sowing and harvest

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
    if (HS_ndays >= (HS_ndays_min - 0.2_r8)) then
         heatwave_crop = 1.0_r8 
    else
         heatwave_crop = 0.0_r8
    end if

  end subroutine crop_heatstress_ndays


  subroutine calc_HS_factor(HS_factor, HS_ndays, t_veg_day, croplive)

    ! !DESCRIPTION:
    ! function to calculate heat stress instensity by applying a factor to bglfr (increasing LAI decline). function based on Apsim-Nwheat model Asseng et al. 2011:https://doi.org/10.1111/j.1365-2486.2010.02262.x
    ! WIP: different tcrit, tmax values for different crops or climatologies

    ! !ARGUMENTS:
    real(r8),        intent(inout)     :: HS_factor         ! keep track if heatwave is activated
    real(r8),        intent(in)        :: HS_ndays          ! number of crop heat stress days (ndays) should be integer at final implementation
    real(r8),        intent(in)        :: t_veg_day         ! daily vegetation temperature (Kelvin)
    logical,         intent(in)        :: croplive          ! crop between sowing and harvest

    ! !LOCAL VARIABLES:
    integer  :: day_min

    !-----------------------------------------------------------------------

    day_min = 3

    !check  if stress occurs
    if (HS_ndays == day_min .and. croplive .and. t_veg_day >= tcrit) then
       ! onset heatwave
       HS_factor = 1.5_r8
    else if (HS_ndays > day_min .and. croplive .and. t_veg_day >= tcrit) then
      if (t_veg_day <= tmax ) then
          HS_factor = 4 - (1 - (t_veg_day - tcrit)/2)
      else
          HS_factor = 11._r8
      end if
    else
       HS_factor = 1._r8
    end if


  end subroutine calc_HS_factor 


end module CropHeatStress
