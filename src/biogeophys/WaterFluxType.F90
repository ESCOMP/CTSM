module WaterFluxType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type containing water fluxes that apply to both bulk water and
  ! water tracers.
  !
  ! !USES:
  use shr_kind_mod   , only: r8 => shr_kind_r8
  use clm_varpar     , only : nlevsno, nlevsoi
  use clm_varcon     , only : spval
  use decompMod      , only : bounds_type
  use decompMod      , only : BOUNDS_SUBGRID_PATCH, BOUNDS_SUBGRID_COLUMN, BOUNDS_SUBGRID_GRIDCELL
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use AnnualFluxDribbler, only : annual_flux_dribbler_type, annual_flux_dribbler_gridcell
  use WaterInfoBaseType, only : water_info_base_type
  use WaterTracerContainerType, only : water_tracer_container_type
  use WaterTracerUtils, only : AllocateVar1d, AllocateVar2d
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  type, public :: waterflux_type

     class(water_info_base_type), pointer :: info

     ! water fluxes are in units or mm/s

     real(r8), pointer :: qflx_through_snow_patch  (:)   ! patch canopy throughfall of snow (mm H2O/s)
     real(r8), pointer :: qflx_through_liq_patch  (:)    ! patch canopy throughfal of liquid (rain+irrigation) (mm H2O/s)
     real(r8), pointer :: qflx_intercepted_snow_patch(:) ! patch canopy interception of snow (mm H2O/s)
     real(r8), pointer :: qflx_intercepted_liq_patch(:)  ! patch canopy interception of liquid (rain+irrigation) (mm H2O/s)
     real(r8), pointer :: qflx_snocanfall_patch(:)       ! patch rate of excess canopy snow falling off canopy (mm H2O/s)
     real(r8), pointer :: qflx_liqcanfall_patch(:)       ! patch rate of excess canopy liquid falling off canopy (mm H2O/s)
     real(r8), pointer :: qflx_snow_unload_patch(:)      ! patch rate of canopy snow unloading (mm H2O/s)
     real(r8), pointer :: qflx_liq_grnd_col        (:)   ! col liquid (rain+irrigation) on ground after interception (mm H2O/s) [+]
     real(r8), pointer :: qflx_snow_grnd_col       (:)   ! col snow on ground after interception (mm H2O/s) [+]
     real(r8), pointer :: qflx_sub_snow_patch      (:)   ! patch sublimation rate from snow pack (mm H2O /s) [+]
     real(r8), pointer :: qflx_sub_snow_col        (:)   ! col sublimation rate from snow pack (mm H2O /s) [+]
     real(r8), pointer :: qflx_evap_soi_patch      (:)   ! patch soil evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_soi_col        (:)   ! col soil evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_veg_patch      (:)   ! patch vegetation evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_veg_col        (:)   ! col vegetation evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_can_patch      (:)   ! patch evaporation from leaves and stems (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_can_col        (:)   ! col evaporation from leaves and stems (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_tot_patch      (:)   ! patch pft_qflx_evap_soi + pft_qflx_evap_veg + qflx_tran_veg
     real(r8), pointer :: qflx_evap_tot_col        (:)   ! col col_qflx_evap_soi + col_qflx_evap_veg + qflx_tran_veg
     real(r8), pointer :: qflx_evap_grnd_patch     (:)   ! patch ground surface evaporation rate (mm H2O/s) [+]
     real(r8), pointer :: qflx_evap_grnd_col       (:)   ! col ground surface evaporation rate (mm H2O/s) [+]

     ! In the snow capping parametrization excess mass above h2osno_max is removed.  A breakdown of mass into liquid 
     ! and solid fluxes is done, these are represented by qflx_snwcp_liq_col and qflx_snwcp_ice_col. 
     real(r8), pointer :: qflx_snwcp_liq_col       (:)   ! col excess liquid h2o due to snow capping (outgoing) (mm H2O /s)
     real(r8), pointer :: qflx_snwcp_ice_col       (:)   ! col excess solid h2o due to snow capping (outgoing) (mm H2O /s)
     real(r8), pointer :: qflx_snwcp_discarded_liq_col(:) ! col excess liquid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s)
     real(r8), pointer :: qflx_snwcp_discarded_ice_col(:) ! col excess solid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s)
     real(r8), pointer :: qflx_glcice_col(:)              ! col net flux of new glacial ice (growth - melt) (mm H2O/s), passed to GLC; only valid inside the do_smb_c filter
     real(r8), pointer :: qflx_glcice_frz_col (:)         ! col ice growth (positive definite) (mm H2O/s); only valid inside the do_smb_c filter
     real(r8), pointer :: qflx_glcice_melt_col(:)         ! col ice melt (positive definite) (mm H2O/s); only valid inside the do_smb_c filter
     real(r8), pointer :: qflx_glcice_dyn_water_flux_col(:) ! col water flux needed for balance check due to glc_dyn_runoff_routing (mm H2O/s) (positive means addition of water to the system); valid for all columns

     real(r8), pointer :: qflx_tran_veg_patch      (:)   ! patch vegetation transpiration (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_tran_veg_col        (:)   ! col vegetation transpiration (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_dew_snow_patch      (:)   ! patch surface dew added to snow pack (mm H2O /s) [+]
     real(r8), pointer :: qflx_dew_snow_col        (:)   ! col surface dew added to snow pack (mm H2O /s) [+]
     real(r8), pointer :: qflx_dew_grnd_patch      (:)   ! patch ground surface dew formation (mm H2O /s) [+]
     real(r8), pointer :: qflx_dew_grnd_col        (:)   ! col ground surface dew formation (mm H2O /s) [+] (+ = to atm); usually eflx_bot >= 0)

     real(r8), pointer :: qflx_infl_col            (:)   ! col infiltration (mm H2O /s)
     real(r8), pointer :: qflx_surf_col            (:)   ! col total surface runoff (mm H2O /s)
     real(r8), pointer :: qflx_drain_col           (:)   ! col sub-surface runoff (mm H2O /s)
     real(r8), pointer :: qflx_drain_perched_col   (:)   ! col sub-surface runoff from perched wt (mm H2O /s)                                                                                                      
     real(r8), pointer :: qflx_top_soil_col        (:)   ! col net water input into soil from top (mm/s)
     real(r8), pointer :: qflx_floodc_col          (:)   ! col flood water flux at column level
     real(r8), pointer :: qflx_sl_top_soil_col     (:)   ! col liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)
     real(r8), pointer :: qflx_snomelt_col         (:)   ! col snow melt (mm H2O /s)
     real(r8), pointer :: qflx_qrgwl_col           (:)   ! col qflx_surf at glaciers, wetlands, lakes
     real(r8), pointer :: qflx_runoff_col          (:)   ! col total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
     real(r8), pointer :: qflx_runoff_r_col        (:)   ! col Rural total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
     real(r8), pointer :: qflx_runoff_u_col        (:)   ! col urban total runoff (qflx_drain+qflx_surf) (mm H2O /s) 
     real(r8), pointer :: qflx_rsub_sat_col        (:)   ! col soil saturation excess [mm/s]
     real(r8), pointer :: qflx_snofrz_lyr_col      (:,:) ! col snow freezing rate (positive definite) (col,lyr) [kg m-2 s-1]
     real(r8), pointer :: qflx_snofrz_col          (:)   ! col column-integrated snow freezing rate (positive definite) (col) [kg m-2 s-1]
     real(r8), pointer :: qflx_snow_drain_col      (:)   ! col drainage from snow pack
     real(r8), pointer :: qflx_ice_runoff_snwcp_col(:)   ! col solid runoff from snow capping (mm H2O /s)
     real(r8), pointer :: qflx_ice_runoff_xs_col   (:)   ! col solid runoff from excess ice in soil (mm H2O /s)

     real(r8), pointer :: qflx_h2osfc_to_ice_col   (:)   ! col conversion of h2osfc to ice
     real(r8), pointer :: qflx_snow_h2osfc_col     (:)   ! col snow falling on surface water

     ! Dynamic land cover change
     real(r8), pointer :: qflx_liq_dynbal_grc      (:)   ! grc liq dynamic land cover change conversion runoff flux
     real(r8), pointer :: qflx_ice_dynbal_grc      (:)   ! grc ice dynamic land cover change conversion runoff flux

     real(r8), pointer :: qflx_sfc_irrig_col        (:)   ! col surface irrigation flux (mm H2O/s) [+]             
     real(r8), pointer :: qflx_gw_uncon_irrig_col   (:)   ! col unconfined groundwater irrigation flux (mm H2O/s)
     real(r8), pointer :: qflx_gw_uncon_irrig_lyr_col(:,:) ! col unconfined groundwater irrigation flux, separated by layer (mm H2O/s)
     real(r8), pointer :: qflx_gw_con_irrig_col     (:)   ! col confined groundwater irrigation flux (mm H2O/s)
     real(r8), pointer :: qflx_irrig_drip_patch     (:)   ! patch drip irrigation
     real(r8), pointer :: qflx_irrig_sprinkler_patch(:)   ! patch sprinkler irrigation

     ! Objects that help convert once-per-year dynamic land cover changes into fluxes
     ! that are dribbled throughout the year
     type(annual_flux_dribbler_type) :: qflx_liq_dynbal_dribbler
     type(annual_flux_dribbler_type) :: qflx_ice_dynbal_dribbler

   contains
     
     procedure, public  :: Init
     procedure, public  :: Restart      
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     

  end type waterflux_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, info, tracer_vars)

    class(waterflux_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds  
    class(water_info_base_type), intent(in), target :: info
    type(water_tracer_container_type), intent(inout) :: tracer_vars

    this%info => info

    call this%InitAllocate(bounds, tracer_vars)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds, tracer_vars)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(waterflux_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds  
    type(water_tracer_container_type), intent(inout) :: tracer_vars
    !
    ! !LOCAL VARIABLES:
    !------------------------------------------------------------------------

    call AllocateVar1d(var = this%qflx_through_snow_patch, name = 'qflx_through_snow_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_PATCH)
    call AllocateVar1d(var = this%qflx_through_liq_patch, name = 'qflx_through_liq_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_PATCH)
    call AllocateVar1d(var = this%qflx_intercepted_snow_patch, name = 'qflx_intercepted_snow_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_PATCH)
    call AllocateVar1d(var = this%qflx_intercepted_liq_patch, name = 'qflx_intercepted_liq_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_PATCH)
    call AllocateVar1d(var = this%qflx_snocanfall_patch, name = 'qflx_snocanfall_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_PATCH)
    call AllocateVar1d(var = this%qflx_liqcanfall_patch, name = 'qflx_liqcanfall_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_PATCH)
    call AllocateVar1d(var = this%qflx_snow_unload_patch, name = 'qflx_snow_unload_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_PATCH)
    call AllocateVar1d(var = this%qflx_sub_snow_patch, name = 'qflx_sub_snow_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_PATCH, &
         ival = 0.0_r8)
    call AllocateVar1d(var = this%qflx_tran_veg_patch, name = 'qflx_tran_veg_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_PATCH)
    call AllocateVar1d(var = this%qflx_dew_grnd_patch, name = 'qflx_dew_grnd_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_PATCH)
    call AllocateVar1d(var = this%qflx_dew_snow_patch, name = 'qflx_dew_snow_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_PATCH)

    call AllocateVar1d(var = this%qflx_liq_grnd_col, name = 'qflx_liq_grnd_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_snow_grnd_col, name = 'qflx_snow_grnd_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_sub_snow_col, name = 'qflx_sub_snow_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN, &
         ival = 0.0_r8)
    call AllocateVar1d(var = this%qflx_snwcp_liq_col, name = 'qflx_snwcp_liq_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_snwcp_ice_col, name = 'qflx_snwcp_ice_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_snwcp_discarded_liq_col, name = 'qflx_snwcp_discarded_liq_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_snwcp_discarded_ice_col, name = 'qflx_snwcp_discarded_ice_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_glcice_col, name = 'qflx_glcice_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_glcice_frz_col, name = 'qflx_glcice_frz_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_glcice_melt_col, name = 'qflx_glcice_melt_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_glcice_dyn_water_flux_col, name = 'qflx_glcice_dyn_water_flux_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_tran_veg_col, name = 'qflx_tran_veg_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_evap_veg_col, name = 'qflx_evap_veg_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_evap_can_col, name = 'qflx_evap_can_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_evap_soi_col, name = 'qflx_evap_soi_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_evap_tot_col, name = 'qflx_evap_tot_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_evap_grnd_col, name = 'qflx_evap_grnd_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_dew_grnd_col, name = 'qflx_dew_grnd_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_dew_snow_col, name = 'qflx_dew_snow_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_evap_veg_patch, name = 'qflx_evap_veg_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_PATCH)
    call AllocateVar1d(var = this%qflx_evap_can_patch, name = 'qflx_evap_can_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_PATCH)
    call AllocateVar1d(var = this%qflx_evap_soi_patch, name = 'qflx_evap_soi_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_PATCH)
    call AllocateVar1d(var = this%qflx_evap_tot_patch, name = 'qflx_evap_tot_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_PATCH)
    call AllocateVar1d(var = this%qflx_evap_grnd_patch, name = 'qflx_evap_grnd_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_PATCH)

    call AllocateVar1d(var = this%qflx_infl_col, name = 'qflx_infl_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_surf_col, name = 'qflx_surf_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_drain_col, name = 'qflx_drain_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_drain_perched_col, name = 'qflx_drain_perched_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_top_soil_col, name = 'qflx_top_soil_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_snomelt_col, name = 'qflx_snomelt_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_snofrz_col, name = 'qflx_snofrz_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar2d(var = this%qflx_snofrz_lyr_col, name = 'qflx_snofrz_lyr_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN, &
         dim2beg = -nlevsno+1, dim2end = 0)
    call AllocateVar1d(var = this%qflx_snow_drain_col, name = 'qflx_snow_drain_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_ice_runoff_snwcp_col, name = 'qflx_ice_runoff_snwcp_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_ice_runoff_xs_col, name = 'qflx_ice_runoff_xs_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_qrgwl_col, name = 'qflx_qrgwl_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_floodc_col, name = 'qflx_floodc_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_sl_top_soil_col, name = 'qflx_sl_top_soil_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_runoff_col, name = 'qflx_runoff_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_runoff_r_col, name = 'qflx_runoff_r_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_runoff_u_col, name = 'qflx_runoff_u_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_rsub_sat_col, name = 'qflx_rsub_sat_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)

    call AllocateVar1d(var = this%qflx_h2osfc_to_ice_col, name = 'qflx_h2osfc_to_ice_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qflx_snow_h2osfc_col, name = 'qflx_snow_h2osfc_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)

    call AllocateVar1d(var = this%qflx_liq_dynbal_grc, name = 'qflx_liq_dynbal_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL)
    call AllocateVar1d(var = this%qflx_ice_dynbal_grc, name = 'qflx_ice_dynbal_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL)

    call AllocateVar1d(var = this%qflx_sfc_irrig_col, name = 'qflx_sfc_irrig_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)

    call AllocateVar1d(var = this%qflx_gw_uncon_irrig_col, name = 'qflx_gw_uncon_irrig_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)

    call AllocateVar2d(var = this%qflx_gw_uncon_irrig_lyr_col, name = 'qflx_gw_uncon_irrig_lyr_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN, &
         dim2beg = 1, dim2end = nlevsoi)

    call AllocateVar1d(var = this%qflx_gw_con_irrig_col, name = 'qflx_gw_con_irrig_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)

    call AllocateVar1d(var = this%qflx_irrig_drip_patch, name = 'qflx_irrig_drip_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_PATCH)

    call AllocateVar1d(var = this%qflx_irrig_sprinkler_patch, name = 'qflx_irrig_sprinkler_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_PATCH)
    
    this%qflx_liq_dynbal_dribbler = annual_flux_dribbler_gridcell( &
         bounds = bounds, &
         name = this%info%fname('qflx_liq_dynbal'), &
         units = 'mm H2O')

    this%qflx_ice_dynbal_dribbler = annual_flux_dribbler_gridcell( &
         bounds = bounds, &
         name = this%info%fname('qflx_ice_dynbal'), &
         units = 'mm H2O')

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod , only : hist_addfld1d, hist_addfld2d, no_snow_normal
    !
    ! !ARGUMENTS:
    class(waterflux_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    integer           :: begg, endg
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    this%qflx_through_liq_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QDIRECT_THROUGHFALL'), &
         units='mm/s', &
         avgflag='A', &
         long_name=this%info%lname('direct throughfall of liquid (rain + above-canopy irrigation)'), &
         ptr_patch=this%qflx_through_liq_patch, c2l_scale_type='urbanf', default='inactive')

    this%qflx_through_snow_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QDIRECT_THROUGHFALL_SNOW'), &
         units='mm/s', &
         avgflag='A', &
         long_name=this%info%lname('direct throughfall of snow'), &
         ptr_patch=this%qflx_through_snow_patch, c2l_scale_type='urbanf', default='inactive')

    this%qflx_liqcanfall_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QDRIP'), &
         units='mm/s', &
         avgflag='A', &
         long_name=this%info%lname('rate of excess canopy liquid falling off canopy'), &
         ptr_patch=this%qflx_liqcanfall_patch, c2l_scale_type='urbanf', default='inactive')

    this%qflx_snocanfall_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QDRIP_SNOW'), &
         units='mm/s', &
         avgflag='A', &
         long_name=this%info%lname('rate of excess canopy snow falling off canopy'), &
         ptr_patch=this%qflx_snocanfall_patch, c2l_scale_type='urbanf', default='inactive')

    this%qflx_snow_unload_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QSNOUNLOAD'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('canopy snow unloading'), &
         ptr_patch=this%qflx_snow_unload_patch, c2l_scale_type='urbanf')

    this%qflx_top_soil_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QTOPSOIL'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('water input to surface'), &
         ptr_col=this%qflx_top_soil_col, c2l_scale_type='urbanf', default='inactive')

    this%qflx_infl_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QINFL'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('infiltration'), &
         ptr_col=this%qflx_infl_col, c2l_scale_type='urbanf')

    this%qflx_surf_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QOVER'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('total surface runoff (includes QH2OSFC)'), &
         ptr_col=this%qflx_surf_col, c2l_scale_type='urbanf')

    this%qflx_qrgwl_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QRGWL'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname( &
         'surface runoff at glaciers (liquid only), wetlands, lakes; also includes melted ice runoff from QSNWCPICE'), &
         ptr_col=this%qflx_qrgwl_col, c2l_scale_type='urbanf')

    this%qflx_drain_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QDRAI'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('sub-surface drainage'), &
         ptr_col=this%qflx_drain_col, c2l_scale_type='urbanf')

    this%qflx_drain_perched_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QDRAI_PERCH'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('perched wt drainage'), &
         ptr_col=this%qflx_drain_perched_col, c2l_scale_type='urbanf')

    this%qflx_liq_dynbal_grc(begg:endg) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_LIQ_DYNBAL'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('liq dynamic land cover change conversion runoff flux'), &
         ptr_lnd=this%qflx_liq_dynbal_grc)     

    this%qflx_ice_dynbal_grc(begg:endg) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_ICE_DYNBAL'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('ice dynamic land cover change conversion runoff flux'), &
         ptr_lnd=this%qflx_ice_dynbal_grc)

    this%qflx_runoff_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QRUNOFF'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('total liquid runoff not including correction for land use change'), &
         ptr_col=this%qflx_runoff_col, c2l_scale_type='urbanf')

    call hist_addfld1d ( &
         fname=this%info%fname('QRUNOFF_ICE'), &
         units='mm/s', avgflag='A', &
         long_name=this%info%lname('total liquid runoff not incl corret for LULCC (ice landunits only)'), &
         ptr_col=this%qflx_runoff_col, c2l_scale_type='urbanf', l2g_scale_type='ice')

    this%qflx_runoff_u_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QRUNOFF_U'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('Urban total runoff'), &
         ptr_col=this%qflx_runoff_u_col, set_nourb=spval, c2l_scale_type='urbanf', default='inactive')

    this%qflx_runoff_r_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QRUNOFF_R'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('Rural total runoff'), &
         ptr_col=this%qflx_runoff_r_col, set_spec=spval, default='inactive')

    this%qflx_snomelt_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QSNOMELT'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('snow melt rate'), &
         ptr_col=this%qflx_snomelt_col, c2l_scale_type='urbanf')

    call hist_addfld1d ( &
         fname=this%info%fname('QSNOMELT_ICE'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('snow melt (ice landunits only)'), &
         ptr_col=this%qflx_snomelt_col, c2l_scale_type='urbanf', l2g_scale_type='ice')


    this%qflx_snofrz_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QSNOFRZ'), &
         units='kg/m2/s', &
         avgflag='A', &
         long_name=this%info%lname('column-integrated snow freezing rate'), &
         ptr_col=this%qflx_snofrz_col, set_lake=spval, c2l_scale_type='urbanf')

    call hist_addfld1d ( &
         fname=this%info%fname('QSNOFRZ_ICE'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('column-integrated snow freezing rate (ice landunits only)'), &
         ptr_col=this%qflx_snofrz_col, c2l_scale_type='urbanf', l2g_scale_type='ice')

    this%qflx_snofrz_lyr_col(begc:endc,-nlevsno+1:0) = spval
    data2dptr => this%qflx_snofrz_lyr_col(begc:endc,-nlevsno+1:0)
    call hist_addfld2d ( &
         fname=this%info%fname('SNO_FRZ'),  &
         units='kg/m2/s', type2d='levsno', &
         avgflag='A', &
         long_name=this%info%lname('snow freezing rate in each snow layer'), &
         ptr_col=data2dptr, c2l_scale_type='urbanf',no_snow_behavior=no_snow_normal, default='inactive')

    call hist_addfld2d ( &
         fname=this%info%fname('SNO_FRZ_ICE'),  &
         units='mm/s', type2d='levsno', &
         avgflag='A', &
         long_name=this%info%lname('snow freezing rate in each snow layer (ice landunits only)'), &
         ptr_col=data2dptr, c2l_scale_type='urbanf',no_snow_behavior=no_snow_normal, &
         l2g_scale_type='ice', default='inactive')

    this%qflx_snow_drain_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_SNOW_DRAIN'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('drainage from snow pack'), &
         ptr_col=this%qflx_snow_drain_col, c2l_scale_type='urbanf')

    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_SNOW_DRAIN_ICE'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('drainage from snow pack melt (ice landunits only)'), &
         ptr_col=this%qflx_snow_drain_col, c2l_scale_type='urbanf', l2g_scale_type='ice')

    this%qflx_evap_soi_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QSOIL'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname( 'Ground evaporation (soil/snow evaporation + soil/snow sublimation - dew)'), &
         ptr_patch=this%qflx_evap_soi_patch, c2l_scale_type='urbanf')

    call hist_addfld1d ( &
         fname=this%info%fname('QSOIL_ICE'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('Ground evaporation (ice landunits only)'), &
         ptr_patch=this%qflx_evap_soi_patch, c2l_scale_type='urbanf', l2g_scale_type='ice')

    this%qflx_evap_can_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QVEGE'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('canopy evaporation'), &
         ptr_patch=this%qflx_evap_can_patch, set_lake=0._r8, c2l_scale_type='urbanf')

    this%qflx_tran_veg_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QVEGT'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('canopy transpiration'), &
         ptr_patch=this%qflx_tran_veg_patch, set_lake=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d ( &
         fname=this%info%fname('QSNOEVAP'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('evaporation from snow'), &
         ptr_patch=this%qflx_tran_veg_patch, set_lake=0._r8, c2l_scale_type='urbanf')

    this%qflx_snwcp_liq_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QSNOCPLIQ'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('excess liquid h2o due to snow capping not including correction for land use change'), &
         ptr_col=this%qflx_snwcp_liq_col, c2l_scale_type='urbanf')

    this%qflx_snwcp_ice_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QSNWCPICE'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('excess solid h2o due to snow capping not including correction for land use change'), &
         ptr_col=this%qflx_snwcp_ice_col, c2l_scale_type='urbanf')

    this%qflx_glcice_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QICE'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('ice growth/melt'), &
         ptr_col=this%qflx_glcice_col, l2g_scale_type='ice')

    this%qflx_glcice_frz_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QICE_FRZ'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('ice growth'), &
         ptr_col=this%qflx_glcice_frz_col, l2g_scale_type='ice')

    this%qflx_glcice_melt_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QICE_MELT'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('ice melt'), &
         ptr_col=this%qflx_glcice_melt_col, l2g_scale_type='ice')

    this%qflx_liq_grnd_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_LIQ_GRND'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('liquid (rain+irrigation) on ground after interception'), &
         ptr_col=this%qflx_liq_grnd_col, default='inactive', c2l_scale_type='urbanf')

    this%qflx_snow_grnd_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_SNOW_GRND'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('snow on ground after interception'), &
         ptr_col=this%qflx_snow_grnd_col, default='inactive', c2l_scale_type='urbanf')

    this%qflx_evap_grnd_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_EVAP_GRND'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('ground surface evaporation'), &
         ptr_patch=this%qflx_evap_grnd_patch, default='inactive', c2l_scale_type='urbanf')

    this%qflx_evap_veg_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_EVAP_VEG'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('vegetation evaporation'), &
         ptr_patch=this%qflx_evap_veg_patch, default='inactive', c2l_scale_type='urbanf')

    this%qflx_evap_tot_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_EVAP_TOT'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('qflx_evap_soi + qflx_evap_can + qflx_tran_veg'), &
         ptr_patch=this%qflx_evap_tot_patch, c2l_scale_type='urbanf')

    this%qflx_dew_grnd_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_DEW_GRND'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('ground surface dew formation'), &
         ptr_patch=this%qflx_dew_grnd_patch, c2l_scale_type='urbanf')

    this%qflx_sub_snow_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_SUB_SNOW'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('sublimation rate from snow pack (also includes bare ice sublimation from glacier columns)'), &
         ptr_patch=this%qflx_sub_snow_patch, c2l_scale_type='urbanf')

    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_SUB_SNOW_ICE'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('sublimation rate from snow pack (also includes bare ice sublimation from glacier columns) '// &
         '(ice landunits only)'), &
         ptr_patch=this%qflx_sub_snow_patch, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')

    this%qflx_dew_snow_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_DEW_SNOW'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('surface dew added to snow pacK'), &
         ptr_patch=this%qflx_dew_snow_patch, c2l_scale_type='urbanf')

    this%qflx_rsub_sat_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QDRAI_XS'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('saturation excess drainage'), &
         ptr_col=this%qflx_rsub_sat_col, c2l_scale_type='urbanf')

    this%qflx_h2osfc_to_ice_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QH2OSFC_TO_ICE'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('surface water converted to ice'), &
         ptr_col=this%qflx_h2osfc_to_ice_col, default='inactive')

    this%qflx_sfc_irrig_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QIRRIG_FROM_SURFACE'), &
         units='mm/s', &
         avgflag='A', &
         long_name=this%info%lname('water added through surface water irrigation'), &
         ptr_col=this%qflx_sfc_irrig_col)

    this%qflx_gw_uncon_irrig_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QIRRIG_FROM_GW_UNCONFINED'), &
         units='mm/s', &
         avgflag='A', &
         long_name=this%info%lname('water added through unconfined groundwater irrigation'), &
         ptr_col=this%qflx_gw_uncon_irrig_col)

    this%qflx_gw_con_irrig_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QIRRIG_FROM_GW_CONFINED'), &
         units='mm/s', &
         avgflag='A', &
         long_name=this%info%lname('water added through confined groundwater irrigation'), &
         ptr_col=this%qflx_gw_con_irrig_col)

    this%qflx_irrig_drip_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QIRRIG_DRIP'), &
         units='mm/s', &
         avgflag='A', &
         long_name=this%info%lname('water added via drip irrigation'), &
         ptr_patch=this%qflx_irrig_drip_patch, default='inactive')

    this%qflx_irrig_sprinkler_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QIRRIG_SPRINKLER'), &
         units='mm/s', &
         avgflag='A', &
         long_name=this%info%lname('water added via sprinkler irrigation'), &
         ptr_patch=this%qflx_irrig_sprinkler_patch, default='inactive')

  end subroutine InitHistory
  
  
  

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !USES:
    use landunit_varcon, only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    class(waterflux_type), intent(in) :: this
    type(bounds_type) , intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: c,l
    !-----------------------------------------------------------------------

    this%qflx_snocanfall_patch(bounds%begp:bounds%endp)       = 0.0_r8
    this%qflx_liqcanfall_patch(bounds%begp:bounds%endp)       = 0.0_r8
    this%qflx_snow_unload_patch(bounds%begp:bounds%endp)      = 0.0_r8

    this%qflx_evap_grnd_patch(bounds%begp:bounds%endp)        = 0.0_r8
    this%qflx_dew_grnd_patch (bounds%begp:bounds%endp)        = 0.0_r8
    this%qflx_dew_snow_patch (bounds%begp:bounds%endp)        = 0.0_r8

    this%qflx_sfc_irrig_col (bounds%begc:bounds%endc)         = 0.0_r8
    this%qflx_gw_uncon_irrig_col (bounds%begc:bounds%endc)    = 0.0_r8
    this%qflx_gw_uncon_irrig_lyr_col(bounds%begc:bounds%endc,:) = 0.0_r8
    this%qflx_gw_con_irrig_col (bounds%begc:bounds%endc)      = 0.0_r8
    this%qflx_irrig_drip_patch (bounds%begp:bounds%endp)      = 0.0_r8
    this%qflx_irrig_sprinkler_patch (bounds%begp:bounds%endp) = 0.0_r8
    
    this%qflx_evap_grnd_col(bounds%begc:bounds%endc) = 0.0_r8
    this%qflx_dew_grnd_col (bounds%begc:bounds%endc) = 0.0_r8
    this%qflx_dew_snow_col (bounds%begc:bounds%endc) = 0.0_r8
    this%qflx_snow_drain_col(bounds%begc:bounds%endc)  = 0._r8

    ! This variable only gets set in the hydrology filter; need to initialize it to 0 for
    ! the sake of columns outside this filter
    this%qflx_ice_runoff_xs_col(bounds%begc:bounds%endc) = 0._r8

    ! Initialize qflx_glcice_dyn_water_flux_col to 0 for all columns because we want this
    ! flux to remain 0 for columns where is is never set, including non-glacier columns.
    !
    ! Other qflx_glcice fluxes intentionally remain unset (spval) outside the do_smb
    ! filter, so that they are flagged as missing value outside that filter.
    this%qflx_glcice_dyn_water_flux_col(bounds%begc:bounds%endc) = 0._r8

    ! needed for CNNLeaching 
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%qflx_drain_col(c) = 0._r8
          this%qflx_surf_col(c)  = 0._r8
       end if
    end do

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !USES:
    use ncdio_pio, only : file_desc_t, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(waterflux_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    ! needed for SNICAR
    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('qflx_snow_drain')//':'//this%info%fname('qflx_snow_melt'), &
         xtype=ncd_double,  &
         dim1name='column', &
         long_name=this%info%lname('drainage from snow column'), &
         units='mm/s', &
         interpinic_flag='interp', readvar=readvar, data=this%qflx_snow_drain_col)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize qflx_snow_drain to zero
       this%qflx_snow_drain_col(bounds%begc:bounds%endc) = 0._r8
    endif

    call this%qflx_liq_dynbal_dribbler%Restart(bounds, ncid, flag)
    call this%qflx_ice_dynbal_dribbler%Restart(bounds, ncid, flag)

  end subroutine Restart

end module WaterFluxType
