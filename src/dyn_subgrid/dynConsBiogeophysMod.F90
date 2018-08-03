module dynConsBiogeophysMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Handle conservation of biogeophysical quantities (water & energy) with dynamic land
  ! cover.
  !
  ! !USES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use UrbanParamsType   , only : urbanparams_type
  use EnergyFluxType    , only : energyflux_type
  use SoilHydrologyType , only : soilhydrology_type  
  use SoilStateType     , only : soilstate_type
  use TemperatureType   , only : temperature_type
  use WaterFluxBulkType     , only : waterfluxbulk_type
  use WaterStateBulkType    , only : waterstatebulk_type
  use WaterDiagnosticBulkType    , only : waterdiagnosticbulk_type
  use WaterBalanceType    , only : waterbalance_type
  use TotalWaterAndHeatMod, only : ComputeLiqIceMassNonLake, ComputeLiqIceMassLake
  use TotalWaterAndHeatMod, only : ComputeHeatNonLake, ComputeHeatLake
  use TotalWaterAndHeatMod, only : AdjustDeltaHeatForDeltaLiq
  use TotalWaterAndHeatMod, only : heat_base_temp
  use subgridAveMod     , only : p2c, c2g
  use LandunitType      , only : lun
  use ColumnType        , only : col                
  use PatchType         , only : patch                
  use clm_varcon        , only : tfrz, cpliq, hfus
  use dynSubgridControlMod, only : get_for_testing_zero_dynbal_fluxes
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  !
  public :: dyn_hwcontent_init            ! compute grid-level heat and water content, before land cover change
  public :: dyn_hwcontent_final           ! compute grid-level heat and water content, after land cover change; also compute dynbal fluxes
  !
  ! !PRIVATE MEMBER FUNCTIONS
  private :: dyn_water_content            ! compute gridcell total liquid and ice water contents
  private :: dyn_heat_content             ! compute gridcell total heat contents

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !---------------------------------------------------------------------------

contains

  !---------------------------------------------------------------------------
  subroutine dyn_hwcontent_init(bounds, &
       num_nolakec, filter_nolakec, &
       num_lakec, filter_lakec, &
       urbanparams_inst, soilstate_inst, soilhydrology_inst, &
       waterstatebulk_inst, waterdiagnosticbulk_inst, waterbalancebulk_inst, &
       waterfluxbulk_inst, temperature_inst, energyflux_inst)
    !
    ! !DESCRIPTION:
    ! Initialize variables used for dyn_hwcontent, and compute grid cell-level heat
    ! and water content before land cover change
    !
    ! Should be called BEFORE any subgrid weight updates this time step
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds  
    integer                  , intent(in)    :: num_nolakec
    integer                  , intent(in)    :: filter_nolakec(:)
    integer                  , intent(in)    :: num_lakec
    integer                  , intent(in)    :: filter_lakec(:)
    type(urbanparams_type)   , intent(in)    :: urbanparams_inst
    type(soilstate_type)     , intent(in)    :: soilstate_inst
    type(soilhydrology_type) , intent(in)    :: soilhydrology_inst
    type(waterstatebulk_type)    , intent(inout) :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)    , intent(inout) :: waterdiagnosticbulk_inst
    type(waterbalance_type)    , intent(inout) :: waterbalancebulk_inst
    type(waterfluxbulk_type)     , intent(inout) :: waterfluxbulk_inst
    type(temperature_type)   , intent(inout) :: temperature_inst
    type(energyflux_type)    , intent(inout) :: energyflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: g   ! grid cell index

    !-------------------------------------------------------------------------------
    
    call dyn_water_content(bounds, &
         num_nolakec, filter_nolakec, &
         num_lakec, filter_lakec, &
         soilhydrology_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, &
         liquid_mass = waterbalancebulk_inst%liq1_grc(bounds%begg:bounds%endg), &
         ice_mass    = waterbalancebulk_inst%ice1_grc(bounds%begg:bounds%endg))

    call dyn_heat_content( bounds, &
         num_nolakec, filter_nolakec, &
         num_lakec, filter_lakec, &
         urbanparams_inst, soilstate_inst, soilhydrology_inst, &
         temperature_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, &
         heat_grc = temperature_inst%heat1_grc(bounds%begg:bounds%endg), &
         liquid_water_temp_grc = temperature_inst%liquid_water_temp1_grc(bounds%begg:bounds%endg))

  end subroutine dyn_hwcontent_init

  !---------------------------------------------------------------------------
  subroutine dyn_hwcontent_final(bounds, &
       num_nolakec, filter_nolakec, &
       num_lakec, filter_lakec, &
       urbanparams_inst, soilstate_inst, soilhydrology_inst, &
       waterstatebulk_inst, waterdiagnosticbulk_inst, waterbalancebulk_inst, &
       waterfluxbulk_inst, temperature_inst, energyflux_inst)
    !
    ! !DESCRIPTION:
    ! Compute grid cell-level heat and water content after land cover change, and compute
    ! the dynbal fluxes
    !
    ! Should be called AFTER all subgrid weight updates this time step
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds  
    integer                  , intent(in)    :: num_nolakec
    integer                  , intent(in)    :: filter_nolakec(:)
    integer                  , intent(in)    :: num_lakec
    integer                  , intent(in)    :: filter_lakec(:)
    type(urbanparams_type)   , intent(in)    :: urbanparams_inst
    type(soilstate_type)     , intent(in)    :: soilstate_inst
    type(soilhydrology_type) , intent(in)    :: soilhydrology_inst
    type(waterstatebulk_type)    , intent(inout) :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)    , intent(inout) :: waterdiagnosticbulk_inst
    type(waterbalance_type)    , intent(inout) :: waterbalancebulk_inst
    type(waterfluxbulk_type)     , intent(inout) :: waterfluxbulk_inst
    type(temperature_type)   , intent(inout) :: temperature_inst
    type(energyflux_type)    , intent(inout) :: energyflux_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: begg, endg
    integer  :: g     ! grid cell index
    real(r8) :: delta_liq(bounds%begg:bounds%endg)  ! change in gridcell h2o liq content
    real(r8) :: delta_ice(bounds%begg:bounds%endg)  ! change in gridcell h2o ice content
    real(r8) :: delta_heat(bounds%begg:bounds%endg) ! change in gridcell heat content
    !---------------------------------------------------------------------------

    begg = bounds%begg
    endg = bounds%endg

    call dyn_water_content(bounds, &
         num_nolakec, filter_nolakec, &
         num_lakec, filter_lakec, &
         soilhydrology_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, &
         liquid_mass = waterbalancebulk_inst%liq2_grc(bounds%begg:bounds%endg), &
         ice_mass    = waterbalancebulk_inst%ice2_grc(bounds%begg:bounds%endg))

    call dyn_heat_content( bounds,                                &
         num_nolakec, filter_nolakec, &
         num_lakec, filter_lakec, &
         urbanparams_inst, soilstate_inst, soilhydrology_inst, &
         temperature_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, &
         heat_grc = temperature_inst%heat2_grc(bounds%begg:bounds%endg), &
         liquid_water_temp_grc = temperature_inst%liquid_water_temp2_grc(bounds%begg:bounds%endg))

    if (get_for_testing_zero_dynbal_fluxes()) then
       do g = begg, endg
          delta_liq(g) = 0._r8
          delta_ice(g) = 0._r8
          delta_heat(g) = 0._r8
       end do
    else
       do g = begg, endg
          delta_liq(g)  = waterbalancebulk_inst%liq2_grc(g) - waterbalancebulk_inst%liq1_grc(g)
          delta_ice(g)  = waterbalancebulk_inst%ice2_grc(g) - waterbalancebulk_inst%ice1_grc(g)
          delta_heat(g) = temperature_inst%heat2_grc(g) - temperature_inst%heat1_grc(g)
       end do
    end if

    call AdjustDeltaHeatForDeltaLiq( &
         bounds, &
         delta_liq = delta_liq(bounds%begg:bounds%endg), &
         liquid_water_temp1 = temperature_inst%liquid_water_temp1_grc(bounds%begg:bounds%endg), &
         liquid_water_temp2 = temperature_inst%liquid_water_temp2_grc(bounds%begg:bounds%endg), &
         delta_heat = delta_heat(bounds%begg:bounds%endg))

    call waterfluxbulk_inst%qflx_liq_dynbal_dribbler%set_curr_delta(bounds, &
         delta_liq(begg:endg))
    call waterfluxbulk_inst%qflx_liq_dynbal_dribbler%get_curr_flux(bounds, &
         waterfluxbulk_inst%qflx_liq_dynbal_grc(begg:endg))

    call waterfluxbulk_inst%qflx_ice_dynbal_dribbler%set_curr_delta(bounds, &
         delta_ice(begg:endg))
    call waterfluxbulk_inst%qflx_ice_dynbal_dribbler%get_curr_flux(bounds, &
         waterfluxbulk_inst%qflx_ice_dynbal_grc(begg:endg))

    call energyflux_inst%eflx_dynbal_dribbler%set_curr_delta(bounds, &
         delta_heat(begg:endg))
    call energyflux_inst%eflx_dynbal_dribbler%get_curr_flux(bounds, &
         energyflux_inst%eflx_dynbal_grc(begg:endg))

  end subroutine dyn_hwcontent_final

  !-----------------------------------------------------------------------
  subroutine dyn_water_content(bounds, &
       num_nolakec, filter_nolakec, &
       num_lakec, filter_lakec, &
       soilhydrology_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, &
       liquid_mass, ice_mass)
    !
    ! !DESCRIPTION:
    ! Compute gridcell total liquid and ice water contents
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds  
    integer                  , intent(in)    :: num_nolakec
    integer                  , intent(in)    :: filter_nolakec(:)
    integer                  , intent(in)    :: num_lakec
    integer                  , intent(in)    :: filter_lakec(:)
    type(soilhydrology_type) , intent(in)    :: soilhydrology_inst
    type(waterstatebulk_type)    , intent(in)    :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)    , intent(in)    :: waterdiagnosticbulk_inst
    real(r8)                 , intent(out)   :: liquid_mass( bounds%begg: ) ! kg m-2
    real(r8)                 , intent(out)   :: ice_mass( bounds%begg: )    ! kg m-2
    !
    ! !LOCAL VARIABLES:
    real(r8) :: liquid_mass_col(bounds%begc:bounds%endc) ! kg m-2
    real(r8) :: ice_mass_col(bounds%begc:bounds%endc)    ! kg m-2

    character(len=*), parameter :: subname = 'dyn_water_content'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(liquid_mass) == (/bounds%endg/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(ice_mass) == (/bounds%endg/)), errMsg(sourcefile, __LINE__))

    call ComputeLiqIceMassNonLake(bounds, num_nolakec, filter_nolakec, &
         soilhydrology_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, &
         liquid_mass_col(bounds%begc:bounds%endc), &
         ice_mass_col(bounds%begc:bounds%endc))

    call ComputeLiqIceMassLake(bounds, num_lakec, filter_lakec, &
         waterstatebulk_inst, &
         liquid_mass_col(bounds%begc:bounds%endc), &
         ice_mass_col(bounds%begc:bounds%endc))

    call c2g(bounds, &
         carr = liquid_mass_col(bounds%begc:bounds%endc), &
         garr = liquid_mass(bounds%begg:bounds%endg), &
         c2l_scale_type = 'urbanf', &
         l2g_scale_type = 'unity')

    call c2g(bounds, &
         carr = ice_mass_col(bounds%begc:bounds%endc), &
         garr = ice_mass(bounds%begg:bounds%endg), &
         c2l_scale_type = 'urbanf', &
         l2g_scale_type = 'unity')

  end subroutine dyn_water_content


  !---------------------------------------------------------------------------
  subroutine dyn_heat_content(bounds, &
       num_nolakec, filter_nolakec, &
       num_lakec, filter_lakec, &
       urbanparams_inst, soilstate_inst, soilhydrology_inst, &
       temperature_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, &
       heat_grc, liquid_water_temp_grc)

    ! !DESCRIPTION:
    ! Compute grid-level heat and water content to track conservation with respect to
    ! dynamic land cover.
    !
    ! Heat content is computed relative to a baseline of 0 C. So temperatures above 0 C
    ! lead to a positive heat content, temperatures below 0 C lead to a negative heat
    ! content. For water, the baseline is considered to be ice at 0 C, so for liquid water
    ! we include the latent heat of fusion.
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)  :: bounds  
    integer                  , intent(in)  :: num_nolakec
    integer                  , intent(in)  :: filter_nolakec(:)
    integer                  , intent(in)  :: num_lakec
    integer                  , intent(in)  :: filter_lakec(:)
    type(urbanparams_type)   , intent(in)  :: urbanparams_inst
    type(soilstate_type)     , intent(in)  :: soilstate_inst
    type(soilhydrology_type) , intent(in)  :: soilhydrology_inst
    type(temperature_type)   , intent(in)  :: temperature_inst
    type(waterstatebulk_type)    , intent(in)  :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)    , intent(in)  :: waterdiagnosticbulk_inst

    real(r8)                 , intent(out) :: heat_grc( bounds%begg: ) ! total heat content for each grid cell [J/m^2]
    real(r8)                 , intent(out) :: liquid_water_temp_grc( bounds%begg: ) ! weighted average liquid water temperature for each grid cell (K)

    !
    ! !LOCAL VARIABLES:
    integer  :: g

    real(r8) :: heat_col(bounds%begc:bounds%endc)  ! sum of heat content for all columns [J/m^2]
    real(r8) :: heat_liquid_col(bounds%begc:bounds%endc) ! sum of heat content for all columns: liquid water, excluding latent heat [J/m^2]
    real(r8) :: cv_liquid_col(bounds%begc:bounds%endc) ! sum of liquid heat capacity for all columns [J/(m^2 K)]

    real(r8) :: heat_liquid_grc(bounds%begg:bounds%endg) ! heat_liquid_col averaged to grid cell [J/m^2]
    real(r8) :: cv_liquid_grc(bounds%begg:bounds%endg) ! cv_liquid_col averaged to grid cell [J/(m^2 K)]
    !-------------------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(heat_grc) == (/bounds%endg/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(liquid_water_temp_grc) == (/bounds%endg/)), errMsg(sourcefile, __LINE__))

    call ComputeHeatNonLake(bounds, num_nolakec, filter_nolakec, &
         urbanparams_inst, soilstate_inst, &
         temperature_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, soilhydrology_inst, &
         heat = heat_col(bounds%begc:bounds%endc), &
         heat_liquid = heat_liquid_col(bounds%begc:bounds%endc), &
         cv_liquid = cv_liquid_col(bounds%begc:bounds%endc))

    call ComputeHeatLake(bounds, num_lakec, filter_lakec, &
         soilstate_inst, temperature_inst, waterstatebulk_inst, &
         heat = heat_col(bounds%begc:bounds%endc), &
         heat_liquid = heat_liquid_col(bounds%begc:bounds%endc), &
         cv_liquid = cv_liquid_col(bounds%begc:bounds%endc))

    call c2g(bounds, &
         carr = heat_col(bounds%begc:bounds%endc), &
         garr = heat_grc(bounds%begg:bounds%endg), &
         c2l_scale_type = 'urbanf', &
         l2g_scale_type = 'unity')

    call c2g(bounds, &
         carr = heat_liquid_col(bounds%begc:bounds%endc), &
         garr = heat_liquid_grc(bounds%begg:bounds%endg), &
         c2l_scale_type = 'urbanf', &
         l2g_scale_type = 'unity')

    call c2g(bounds, &
         carr = cv_liquid_col(bounds%begc:bounds%endc), &
         garr = cv_liquid_grc(bounds%begg:bounds%endg), &
         c2l_scale_type = 'urbanf', &
         l2g_scale_type = 'unity')

    do g = bounds%begg, bounds%endg
       if (cv_liquid_grc(g) > 0._r8) then
          liquid_water_temp_grc(g) = &
               (heat_liquid_grc(g) / cv_liquid_grc(g)) + heat_base_temp
       else
          ! 0 or negative water mass in this grid cell: set an arbitrary temperature
          liquid_water_temp_grc(g) = tfrz
       end if
    end do

  end subroutine dyn_heat_content

end module dynConsBiogeophysMod
