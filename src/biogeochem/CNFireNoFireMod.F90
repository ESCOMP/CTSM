module ctsm_CNFireNoFireMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! module for fire dynamics with fire explicitly turned off
  !
  ! !USES:
  use shr_kind_mod                       , only : r8 => shr_kind_r8, CL => shr_kind_CL
  use ctsm_Decomp                          , only : bounds_type
  use ctsm_Atm2LndType                        , only : atm2lnd_type
  use ctsm_CNVegStateType                     , only : cnveg_state_type
  use ctsm_CNVegCarbonStateType               , only : cnveg_carbonstate_type
  use ctsm_CNVegCarbonFluxType                , only : cnveg_carbonflux_type
  use ctsm_CNVegNitrogenStateType             , only : cnveg_nitrogenstate_type
  use ctsm_CNVegNitrogenFluxType              , only : cnveg_nitrogenflux_type
  use ctsm_EnergyFluxType                     , only : energyflux_type
  use ctsm_SaturatedExcessRunoff           , only : saturated_excess_runoff_type
  use ctsm_WaterDiagnosticBulkType                     , only : waterdiagnosticbulk_type
  use ctsm_WaterAtm2LndBulkType                     , only : wateratm2lndbulk_type
  use ctsm_FireMethodType                     , only : fire_method_type
  use ctsm_CNFireBaseMod                      , only : cnfire_base_type
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  public :: cnfire_nofire_type
  !
  type, extends(cnfire_base_type) :: cnfire_nofire_type
    private
  contains
    !
    ! !PUBLIC MEMBER FUNCTIONS:
    procedure, public :: need_lightning_and_popdens
    procedure, public :: CNFireArea    ! Calculate fire area
  end type cnfire_nofire_type

contains

  !-----------------------------------------------------------------------
  function need_lightning_and_popdens(this)
    ! !ARGUMENTS:
    class(cnfire_nofire_type), intent(in) :: this
    logical :: need_lightning_and_popdens  ! function result
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'need_lightning_and_popdens'
    !-----------------------------------------------------------------------

    need_lightning_and_popdens = .false.
  end function need_lightning_and_popdens

  !-----------------------------------------------------------------------
  subroutine CNFireArea (this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       atm2lnd_inst, energyflux_inst, saturated_excess_runoff_inst, &
       waterdiagnosticbulk_inst, wateratm2lndbulk_inst, &
       cnveg_state_inst, cnveg_carbonstate_inst, totlitc_col, decomp_cpools_vr_col, t_soi17cm_col)
    !
    ! !DESCRIPTION:
    ! Computes column-level burned area 
    !
    ! !USES:
    use ctsm_SubgridAve                      , only : p2c
    !
    ! !ARGUMENTS:
    class(cnfire_nofire_type)                             :: this
    type(bounds_type)                     , intent(in)    :: bounds 
    integer                               , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                               , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                               , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                               , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(atm2lnd_type)                    , intent(in)    :: atm2lnd_inst
    type(energyflux_type)                 , intent(in)    :: energyflux_inst
    type(saturated_excess_runoff_type)    , intent(in)    :: saturated_excess_runoff_inst
    type(waterdiagnosticbulk_type)                 , intent(in)    :: waterdiagnosticbulk_inst
    type(wateratm2lndbulk_type)                 , intent(in)    :: wateratm2lndbulk_inst
    type(cnveg_state_type)                , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)          , intent(inout) :: cnveg_carbonstate_inst
    real(r8)                              , intent(in)    :: totlitc_col(bounds%begc:)
    real(r8)                              , intent(in)    :: decomp_cpools_vr_col(bounds%begc:,1:,1:)
    real(r8)                              , intent(in)    :: t_soi17cm_col(bounds%begc:)
    !
    ! !LOCAL VARIABLES:
    integer  :: c,fc   ! index variables
    !-----------------------------------------------------------------------

    associate(                                                                      & 
         cropf_col          => cnveg_state_inst%cropf_col                      , & ! Input:  [real(r8) (:)     ]  cropland fraction in veg column                   
         baf_crop           => cnveg_state_inst%baf_crop_col                   , & ! Output: [real(r8) (:)     ]  burned area fraction for cropland (/sec)  
         baf_peatf          => cnveg_state_inst%baf_peatf_col                  , & ! Output: [real(r8) (:)     ]  burned area fraction for peatland (/sec)  
         fbac               => cnveg_state_inst%fbac_col                       , & ! Output: [real(r8) (:)     ]  total burned area out of conversion (/sec)
         fbac1              => cnveg_state_inst%fbac1_col                      , & ! Output: [real(r8) (:)     ]  burned area out of conversion region due to land use fire
         lfc                => cnveg_state_inst%lfc_col                        , & ! Output: [real(r8) (:)     ]  conversion area frac. of BET+BDT that haven't burned before
         leafc              => cnveg_carbonstate_inst%leafc_patch              , & ! Input:  [real(r8) (:)     ]  (gC/m2) leaf C                                    
         leafc_col          => cnveg_carbonstate_inst%leafc_col                , & ! Output: [real(r8) (:)     ]  leaf carbon at column level 
         farea_burned       => cnveg_state_inst%farea_burned_col                 & ! Output: [real(r8) (:)     ]  total fractional area burned (/sec)
         )
 
      !pft to column average 
      call p2c(bounds, num_soilc, filter_soilc, &
           leafc(bounds%begp:bounds%endp), &
           leafc_col(bounds%begc:bounds%endc))
     !
     ! begin column loop to calculate fractional area affected by fire
     !
     do fc = 1, num_soilc
        c = filter_soilc(fc)

        ! zero out the fire area

        farea_burned(c) = 0._r8
        baf_crop(c)     = 0._r8
        baf_peatf(c)    = 0._r8
        fbac(c)         = 0._r8
        fbac1(c)        = 0._r8
        cropf_col(c)    = 0._r8 
        lfc(c)          = 0._r8
        ! with NOFIRE, tree carbon is still removed in landuse change regions by the
        ! landuse code
     end do  ! end of column loop

   end associate

 end subroutine CNFireArea

end module ctsm_CNFireNoFireMod
