module WaterfluxType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  !
  ! !USES:
  use shr_kind_mod   , only: r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : nlevsno, nlevsoi
  use clm_varcon     , only : spval
  use decompMod      , only : bounds_type
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use PatchType      , only : patch   
  use CNSharedParamsMod           , only : use_fun             
  use AnnualFluxDribbler, only : annual_flux_dribbler_type, annual_flux_dribbler_gridcell
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  type, public :: waterflux_type

     ! water fluxes are in units or mm/s

     real(r8), pointer :: qflx_prec_grnd_patch     (:)   ! patch water onto ground including canopy runoff [kg/(m2 s)]
     real(r8), pointer :: qflx_prec_grnd_col       (:)   ! col water onto ground including canopy runoff [kg/(m2 s)]
     real(r8), pointer :: qflx_rain_grnd_patch     (:)   ! patch rain on ground after interception (mm H2O/s) [+]
     real(r8), pointer :: qflx_rain_grnd_col       (:)   ! col rain on ground after interception (mm H2O/s) [+]
     real(r8), pointer :: qflx_snow_grnd_patch     (:)   ! patch snow on ground after interception (mm H2O/s) [+]
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
     real(r8), pointer :: qflx_phs_neg_col         (:)   ! col sum of negative hydraulic redistribution fluxes (mm H2O/s) [+]

     ! In the snow capping parametrization excess mass above h2osno_max is removed.  A breakdown of mass into liquid 
     ! and solid fluxes is done, these are represented by qflx_snwcp_liq_col and qflx_snwcp_ice_col. 
     real(r8), pointer :: qflx_snwcp_liq_col       (:)   ! col excess liquid h2o due to snow capping (outgoing) (mm H2O /s)
     real(r8), pointer :: qflx_snwcp_ice_col       (:)   ! col excess solid h2o due to snow capping (outgoing) (mm H2O /s)
     real(r8), pointer :: qflx_snwcp_discarded_liq_col(:) ! col excess liquid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s)
     real(r8), pointer :: qflx_snwcp_discarded_ice_col(:) ! col excess solid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s)

     real(r8), pointer :: qflx_tran_veg_patch      (:)   ! patch vegetation transpiration (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_tran_veg_col        (:)   ! col vegetation transpiration (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_dew_snow_patch      (:)   ! patch surface dew added to snow pack (mm H2O /s) [+]
     real(r8), pointer :: qflx_dew_snow_col        (:)   ! col surface dew added to snow pack (mm H2O /s) [+]
     real(r8), pointer :: qflx_dew_grnd_patch      (:)   ! patch ground surface dew formation (mm H2O /s) [+]
     real(r8), pointer :: qflx_dew_grnd_col        (:)   ! col ground surface dew formation (mm H2O /s) [+] (+ = to atm); usually eflx_bot >= 0)
     real(r8), pointer :: qflx_prec_intr_patch     (:)   ! patch interception of precipitation [mm/s]
     real(r8), pointer :: qflx_prec_intr_col       (:)   ! col interception of precipitation [mm/s]
     real(r8), pointer :: qflx_snowindunload_patch (:)   ! patch canopy snow wind unloading (mm H2O /s)
     real(r8), pointer :: qflx_snowindunload_col   (:)   ! col canopy snow wind unloading (mm H2O /s)
     real(r8), pointer :: qflx_snotempunload_patch (:)   ! patch canopy snow temp unloading (mm H2O /s) 
     real(r8), pointer :: qflx_snotempunload_col   (:)   ! col canopy snow temp unloading (mm H2O /s) 

     real(r8), pointer :: qflx_ev_snow_patch       (:)   ! patch evaporation heat flux from snow       (mm H2O/s) [+ to atm]
     real(r8), pointer :: qflx_ev_snow_col         (:)   ! col evaporation heat flux from snow         (mm H2O/s) [+ to atm]
     real(r8), pointer :: qflx_ev_soil_patch       (:)   ! patch evaporation heat flux from soil       (mm H2O/s) [+ to atm]
     real(r8), pointer :: qflx_ev_soil_col         (:)   ! col evaporation heat flux from soil         (mm H2O/s) [+ to atm]
     real(r8), pointer :: qflx_ev_h2osfc_patch     (:)   ! patch evaporation heat flux from soil       (mm H2O/s) [+ to atm]
     real(r8), pointer :: qflx_ev_h2osfc_col       (:)   ! col evaporation heat flux from soil         (mm H2O/s) [+ to atm]

     real(r8), pointer :: qflx_adv_col             (:,:) ! col advective flux across different soil layer interfaces [mm H2O/s] [+ downward]
     real(r8), pointer :: qflx_rootsoi_col         (:,:) ! col root and soil water exchange [mm H2O/s] [+ into root]
     real(r8), pointer :: qflx_infl_col            (:)   ! col infiltration (mm H2O /s)
     real(r8), pointer :: qflx_surf_col            (:)   ! col surface runoff (mm H2O /s)
     real(r8), pointer :: qflx_drain_col           (:)   ! col sub-surface runoff (mm H2O /s)
     real(r8), pointer :: qflx_top_soil_col        (:)   ! col net water input into soil from top (mm/s)
     real(r8), pointer :: qflx_h2osfc_to_ice_col   (:)   ! col conversion of h2osfc to ice
     real(r8), pointer :: qflx_h2osfc_surf_col     (:)   ! col surface water runoff
     real(r8), pointer :: qflx_snow_h2osfc_col     (:)   ! col snow falling on surface water
     real(r8), pointer :: qflx_drain_perched_col   (:)   ! col sub-surface runoff from perched wt (mm H2O /s)
     real(r8), pointer :: qflx_deficit_col         (:)   ! col water deficit to keep non-negative liquid water content (mm H2O)   
     real(r8), pointer :: qflx_floodc_col          (:)   ! col flood water flux at column level
     real(r8), pointer :: qflx_sl_top_soil_col     (:)   ! col liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)
     real(r8), pointer :: qflx_snomelt_col         (:)   ! col snow melt (mm H2O /s)
     real(r8), pointer :: qflx_snomelt_lyr_col     (:,:) ! col snow melt in each layer (mm H2O /s)
     real(r8), pointer :: qflx_snow_drain_col      (:)   ! col drainage from snow pack
     real(r8), pointer :: qflx_qrgwl_col           (:)   ! col qflx_surf at glaciers, wetlands, lakes
     real(r8), pointer :: qflx_runoff_rain_to_snow_conversion_col(:) ! col runoff flux from rain-to-snow conversion, when this conversion leads to immediate runoff rather than snow (mm H2O /s)
     real(r8), pointer :: qflx_runoff_col          (:)   ! col total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
     real(r8), pointer :: qflx_runoff_r_col        (:)   ! col Rural total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
     real(r8), pointer :: qflx_runoff_u_col        (:)   ! col urban total runoff (qflx_drain+qflx_surf) (mm H2O /s) 
     real(r8), pointer :: qflx_ice_runoff_snwcp_col(:)   ! col solid runoff from snow capping (mm H2O /s)
     real(r8), pointer :: qflx_ice_runoff_xs_col   (:)   ! col solid runoff from excess ice in soil (mm H2O /s)
     real(r8), pointer :: qflx_rsub_sat_col        (:)   ! col soil saturation excess [mm/s]
     real(r8), pointer :: qflx_snofrz_lyr_col      (:,:) ! col snow freezing rate (positive definite) (col,lyr) [kg m-2 s-1]
     real(r8), pointer :: qflx_snofrz_col          (:)   ! col column-integrated snow freezing rate (positive definite) (col) [kg m-2 s-1]
     real(r8), pointer :: qflx_drain_vr_col        (:,:) ! col liquid water losted as drainage (m /time step)
     real(r8), pointer :: snow_sources_col         (:)   ! col snow sources (mm H2O/s)
     real(r8), pointer :: snow_sinks_col           (:)   ! col snow sinks (mm H2O/s)

     ! Dynamic land cover change
     real(r8), pointer :: qflx_liq_dynbal_grc      (:)   ! grc liq dynamic land cover change conversion runoff flux
     real(r8), pointer :: qflx_ice_dynbal_grc      (:)   ! grc ice dynamic land cover change conversion runoff flux

     ! Objects that help convert once-per-year dynamic land cover changes into fluxes
     ! that are dribbled throughout the year
     type(annual_flux_dribbler_type) :: qflx_liq_dynbal_dribbler
     type(annual_flux_dribbler_type) :: qflx_ice_dynbal_dribbler

     ! ET accumulation
     real(r8), pointer :: AnnEt                    (:)   ! Annual average ET flux mmH20/s


   contains
 
     
     
     procedure, public  :: Init
     procedure, public  :: Restart      
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     
     procedure, public  :: InitAccBuffer
     procedure, public  :: InitAccVars
     procedure, public  :: UpdateAccVars

  end type waterflux_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(waterflux_type) :: this
    type(bounds_type), intent(in)    :: bounds  

    call this%InitAllocate(bounds) ! same as "call initAllocate_type(hydro, bounds)"
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(waterflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    allocate(this%qflx_prec_intr_patch     (begp:endp))              ; this%qflx_prec_intr_patch     (:)   = nan
    allocate(this%qflx_prec_grnd_patch     (begp:endp))              ; this%qflx_prec_grnd_patch     (:)   = nan
    allocate(this%qflx_rain_grnd_patch     (begp:endp))              ; this%qflx_rain_grnd_patch     (:)   = nan
    allocate(this%qflx_snow_grnd_patch     (begp:endp))              ; this%qflx_snow_grnd_patch     (:)   = nan
    allocate(this%qflx_sub_snow_patch      (begp:endp))              ; this%qflx_sub_snow_patch      (:)   = 0.0_r8
    allocate(this%qflx_tran_veg_patch      (begp:endp))              ; this%qflx_tran_veg_patch      (:)   = nan

    allocate(this%qflx_snowindunload_patch (begp:endp))              ; this%qflx_snowindunload_patch (:)   = nan
    allocate(this%qflx_snowindunload_col   (begp:endp))              ; this%qflx_snowindunload_col   (:)   = nan
    allocate(this%qflx_snotempunload_patch (begp:endp))              ; this%qflx_snotempunload_patch (:)   = nan
    allocate(this%qflx_snotempunload_col   (begp:endp))              ; this%qflx_snotempunload_col   (:)   = nan
	
    allocate(this%qflx_dew_grnd_patch      (begp:endp))              ; this%qflx_dew_grnd_patch      (:)   = nan
    allocate(this%qflx_dew_snow_patch      (begp:endp))              ; this%qflx_dew_snow_patch      (:)   = nan

    allocate(this%qflx_prec_intr_col       (begc:endc))              ; this%qflx_prec_intr_col       (:)   = nan
    allocate(this%qflx_prec_grnd_col       (begc:endc))              ; this%qflx_prec_grnd_col       (:)   = nan
    allocate(this%qflx_rain_grnd_col       (begc:endc))              ; this%qflx_rain_grnd_col       (:)   = nan
    allocate(this%qflx_snow_grnd_col       (begc:endc))              ; this%qflx_snow_grnd_col       (:)   = nan
    allocate(this%qflx_sub_snow_col        (begc:endc))              ; this%qflx_sub_snow_col        (:)   = 0.0_r8
    allocate(this%qflx_snwcp_liq_col       (begc:endc))              ; this%qflx_snwcp_liq_col       (:)   = nan
    allocate(this%qflx_snwcp_ice_col       (begc:endc))              ; this%qflx_snwcp_ice_col       (:)   = nan
    allocate(this%qflx_snwcp_discarded_liq_col(begc:endc))           ; this%qflx_snwcp_discarded_liq_col(:) = nan
    allocate(this%qflx_snwcp_discarded_ice_col(begc:endc))           ; this%qflx_snwcp_discarded_ice_col(:) = nan
    allocate(this%qflx_tran_veg_col        (begc:endc))              ; this%qflx_tran_veg_col        (:)   = nan
    allocate(this%qflx_evap_veg_col        (begc:endc))              ; this%qflx_evap_veg_col        (:)   = nan
    allocate(this%qflx_evap_can_col        (begc:endc))              ; this%qflx_evap_can_col        (:)   = nan
    allocate(this%qflx_evap_soi_col        (begc:endc))              ; this%qflx_evap_soi_col        (:)   = nan
    allocate(this%qflx_evap_tot_col        (begc:endc))              ; this%qflx_evap_tot_col        (:)   = nan
    allocate(this%qflx_evap_grnd_col       (begc:endc))              ; this%qflx_evap_grnd_col       (:)   = nan
    allocate(this%qflx_dew_grnd_col        (begc:endc))              ; this%qflx_dew_grnd_col        (:)   = nan
    allocate(this%qflx_dew_snow_col        (begc:endc))              ; this%qflx_dew_snow_col        (:)   = nan
    allocate(this%qflx_evap_veg_patch      (begp:endp))              ; this%qflx_evap_veg_patch      (:)   = nan
    allocate(this%qflx_evap_can_patch      (begp:endp))              ; this%qflx_evap_can_patch      (:)   = nan
    allocate(this%qflx_evap_soi_patch      (begp:endp))              ; this%qflx_evap_soi_patch      (:)   = nan
    allocate(this%qflx_evap_tot_patch      (begp:endp))              ; this%qflx_evap_tot_patch      (:)   = nan
    allocate(this%qflx_evap_grnd_patch     (begp:endp))              ; this%qflx_evap_grnd_patch     (:)   = nan
    allocate(this%qflx_phs_neg_col         (begc:endc))              ; this%qflx_phs_neg_col       (:)   = nan

    allocate( this%qflx_ev_snow_patch      (begp:endp))              ; this%qflx_ev_snow_patch       (:)   = nan
    allocate( this%qflx_ev_snow_col        (begc:endc))              ; this%qflx_ev_snow_col         (:)   = nan
    allocate( this%qflx_ev_soil_patch      (begp:endp))              ; this%qflx_ev_soil_patch       (:)   = nan
    allocate( this%qflx_ev_soil_col        (begc:endc))              ; this%qflx_ev_soil_col         (:)   = nan
    allocate( this%qflx_ev_h2osfc_patch    (begp:endp))              ; this%qflx_ev_h2osfc_patch     (:)   = nan
    allocate( this%qflx_ev_h2osfc_col      (begc:endc))              ; this%qflx_ev_h2osfc_col       (:)   = nan

    allocate(this%qflx_drain_vr_col      (begc:endc,1:nlevsoi))      ; this%qflx_drain_vr_col        (:,:) = nan
    allocate(this%qflx_adv_col             (begc:endc,0:nlevsoi))    ; this%qflx_adv_col             (:,:) = nan
    allocate(this%qflx_rootsoi_col         (begc:endc,1:nlevsoi))    ; this%qflx_rootsoi_col         (:,:) = nan
    allocate(this%qflx_infl_col            (begc:endc))              ; this%qflx_infl_col            (:)   = nan
    allocate(this%qflx_surf_col            (begc:endc))              ; this%qflx_surf_col            (:)   = nan
    allocate(this%qflx_drain_col           (begc:endc))              ; this%qflx_drain_col           (:)   = nan
    allocate(this%qflx_top_soil_col        (begc:endc))              ; this%qflx_top_soil_col        (:)   = nan
    allocate(this%qflx_h2osfc_to_ice_col   (begc:endc))              ; this%qflx_h2osfc_to_ice_col   (:)   = nan
    allocate(this%qflx_h2osfc_surf_col     (begc:endc))              ; this%qflx_h2osfc_surf_col     (:)   = nan
    allocate(this%qflx_snow_h2osfc_col     (begc:endc))              ; this%qflx_snow_h2osfc_col     (:)   = nan
    allocate(this%qflx_snomelt_col         (begc:endc))              ; this%qflx_snomelt_col         (:)   = nan
    allocate(this%qflx_snomelt_lyr_col     (begc:endc,-nlevsno+1:0)) ; this%qflx_snomelt_lyr_col     (:,:) = nan
    allocate(this%qflx_snow_drain_col      (begc:endc))              ; this%qflx_snow_drain_col      (:)   = nan
    allocate(this%qflx_snofrz_col          (begc:endc))              ; this%qflx_snofrz_col          (:)   = nan
    allocate(this%qflx_snofrz_lyr_col      (begc:endc,-nlevsno+1:0)) ; this%qflx_snofrz_lyr_col      (:,:) = nan
    allocate(this%qflx_qrgwl_col           (begc:endc))              ; this%qflx_qrgwl_col           (:)   = nan
    allocate(this%qflx_runoff_rain_to_snow_conversion_col(begc:endc)); this%qflx_runoff_rain_to_snow_conversion_col(:) = nan
    allocate(this%qflx_drain_perched_col   (begc:endc))              ; this%qflx_drain_perched_col   (:)   = nan
    allocate(this%qflx_deficit_col         (begc:endc))              ; this%qflx_deficit_col         (:)   = nan
    allocate(this%qflx_floodc_col          (begc:endc))              ; this%qflx_floodc_col          (:)   = nan
    allocate(this%qflx_sl_top_soil_col     (begc:endc))              ; this%qflx_sl_top_soil_col     (:)   = nan
    allocate(this%qflx_runoff_col          (begc:endc))              ; this%qflx_runoff_col          (:)   = nan
    allocate(this%qflx_runoff_r_col        (begc:endc))              ; this%qflx_runoff_r_col        (:)   = nan
    allocate(this%qflx_runoff_u_col        (begc:endc))              ; this%qflx_runoff_u_col        (:)   = nan
    allocate(this%qflx_ice_runoff_snwcp_col(begc:endc))              ; this%qflx_ice_runoff_snwcp_col(:)   = nan
    allocate(this%qflx_ice_runoff_xs_col   (begc:endc))              ; this%qflx_ice_runoff_xs_col   (:)   = nan
    allocate(this%qflx_rsub_sat_col        (begc:endc))              ; this%qflx_rsub_sat_col        (:)   = nan
    allocate(this%snow_sources_col         (begc:endc))              ; this%snow_sources_col         (:)   = nan   
    allocate(this%snow_sinks_col           (begc:endc))              ; this%snow_sinks_col           (:)   = nan   

    allocate(this%qflx_liq_dynbal_grc      (begg:endg))              ; this%qflx_liq_dynbal_grc      (:)   = nan
    allocate(this%qflx_ice_dynbal_grc      (begg:endg))              ; this%qflx_ice_dynbal_grc      (:)   = nan
    allocate(this%AnnET                    (begc:endc))              ; this%AnnET                    (:)   = nan

    this%qflx_liq_dynbal_dribbler = annual_flux_dribbler_gridcell( &
         bounds = bounds, &
         name = 'qflx_liq_dynbal', &
         units = 'mm H2O')

    this%qflx_ice_dynbal_dribbler = annual_flux_dribbler_gridcell( &
         bounds = bounds, &
         name = 'qflx_ice_dynbal', &
         units = 'mm H2O')

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use clm_varctl  , only : use_cn
    use histFileMod , only : hist_addfld1d, hist_addfld2d, no_snow_normal
    !
    ! !ARGUMENTS:
    class(waterflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    integer           :: begg, endg
    character(10)     :: active
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    this%qflx_top_soil_col(begc:endc) = spval
    call hist_addfld1d (fname='QTOPSOIL',  units='mm/s',  &
         avgflag='A', long_name='water input to surface', &
         ptr_col=this%qflx_top_soil_col, c2l_scale_type='urbanf', default='inactive')

    this%qflx_infl_col(begc:endc) = spval
    call hist_addfld1d (fname='QINFL',  units='mm/s',  &
         avgflag='A', long_name='infiltration', &
         ptr_col=this%qflx_infl_col, c2l_scale_type='urbanf')

    this%qflx_surf_col(begc:endc) = spval
    call hist_addfld1d (fname='QOVER',  units='mm/s',  &
         avgflag='A', long_name='surface runoff', &
         ptr_col=this%qflx_surf_col, c2l_scale_type='urbanf')

    this%qflx_qrgwl_col(begc:endc) = spval
    call hist_addfld1d (fname='QRGWL',  units='mm/s',  &
         avgflag='A', &
         long_name='surface runoff at glaciers (liquid only), wetlands, lakes; also includes melted ice runoff from QSNWCPICE', &
         ptr_col=this%qflx_qrgwl_col, c2l_scale_type='urbanf')

    this%qflx_runoff_rain_to_snow_conversion_col(begc:endc) = spval
    call hist_addfld1d (fname='QRUNOFF_RAIN_TO_SNOW_CONVERSION', units='mm/s', &
         avgflag='A', &
         long_name='liquid runoff from rain-to-snow conversion when this conversion leads to immediate runoff', &
         ptr_col=this%qflx_runoff_rain_to_snow_conversion_col, c2l_scale_type='urbanf')

    this%qflx_drain_col(begc:endc) = spval
    call hist_addfld1d (fname='QDRAI',  units='mm/s',  &
         avgflag='A', long_name='sub-surface drainage', &
         ptr_col=this%qflx_drain_col, c2l_scale_type='urbanf')

    this%qflx_liq_dynbal_grc(begg:endg) = spval
    call hist_addfld1d (fname='QFLX_LIQ_DYNBAL',  units='mm/s',  &  
         avgflag='A', long_name='liq dynamic land cover change conversion runoff flux', &              
         ptr_lnd=this%qflx_liq_dynbal_grc)     

    this%qflx_ice_dynbal_grc(begg:endg) = spval
    call hist_addfld1d (fname='QFLX_ICE_DYNBAL',  units='mm/s',  &
         avgflag='A', long_name='ice dynamic land cover change conversion runoff flux', &                                   
         ptr_lnd=this%qflx_ice_dynbal_grc)

    this%qflx_runoff_col(begc:endc) = spval
    call hist_addfld1d (fname='QRUNOFF',  units='mm/s',  &
         avgflag='A', &
         long_name='total liquid runoff not including correction for land use change', &
         ptr_col=this%qflx_runoff_col, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QRUNOFF_ICE', units='mm/s', avgflag='A', &
         long_name='total liquid runoff not incl corret for LULCC (ice landunits only)', &
         ptr_col=this%qflx_runoff_col, c2l_scale_type='urbanf', l2g_scale_type='ice')

    this%qflx_runoff_u_col(begc:endc) = spval
    call hist_addfld1d (fname='QRUNOFF_U', units='mm/s',  &
         avgflag='A', long_name='Urban total runoff', &
         ptr_col=this%qflx_runoff_u_col, set_nourb=spval, c2l_scale_type='urbanf', default='inactive')

    this%qflx_runoff_r_col(begc:endc) = spval
    call hist_addfld1d (fname='QRUNOFF_R', units='mm/s',  &
         avgflag='A', long_name='Rural total runoff', &
         ptr_col=this%qflx_runoff_r_col, set_spec=spval, default='inactive')

    this%qflx_snow_drain_col(begc:endc) = spval
    call hist_addfld1d (fname='QFLX_SNOW_DRAIN',  units='mm/s',  &
         avgflag='A', long_name='drainage from snow pack', &
         ptr_col=this%qflx_snow_drain_col, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QFLX_SNOW_DRAIN_ICE', units='mm/s',  &
         avgflag='A', long_name='drainage from snow pack melt (ice landunits only)', &
         ptr_col=this%qflx_snow_drain_col, c2l_scale_type='urbanf', l2g_scale_type='ice')

    this%qflx_snomelt_col(begc:endc) = spval
    call hist_addfld1d (fname='QSNOMELT',  units='mm/s',  &
         avgflag='A', long_name='snow melt rate', &
         ptr_col=this%qflx_snomelt_col, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QSNOMELT_ICE', units='mm/s',  &
         avgflag='A', long_name='snow melt (ice landunits only)', &
         ptr_col=this%qflx_snomelt_col, c2l_scale_type='urbanf', l2g_scale_type='ice')

    this%qflx_snomelt_lyr_col(begc:endc,-nlevsno+1:0) = spval
    data2dptr => this%qflx_snomelt_lyr_col(begc:endc,-nlevsno+1:0)
    call hist_addfld2d (fname='SNO_MELT',  units='mm/s', type2d='levsno', &
         avgflag='A', long_name='snow melt rate in each snow layer', &
         ptr_col=data2dptr, c2l_scale_type='urbanf',no_snow_behavior=no_snow_normal, default='inactive')

    call hist_addfld2d (fname='SNO_MELT_ICE',  units='mm/s', type2d='levsno', &
         avgflag='A', long_name='snow melt rate in each snow layer (ice landunits only)', &
         ptr_col=data2dptr, c2l_scale_type='urbanf',no_snow_behavior=no_snow_normal, &
         l2g_scale_type='ice', default='inactive')

    this%qflx_snofrz_col(begc:endc) = spval
    call hist_addfld1d (fname='QSNOFRZ', units='kg/m2/s', &
         avgflag='A', long_name='column-integrated snow freezing rate', &
         ptr_col=this%qflx_snofrz_col, set_lake=spval, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QSNOFRZ_ICE', units='mm/s',  &
         avgflag='A', long_name='column-integrated snow freezing rate (ice landunits only)', &
         ptr_col=this%qflx_snofrz_col, c2l_scale_type='urbanf', l2g_scale_type='ice')

    this%qflx_snofrz_lyr_col(begc:endc,-nlevsno+1:0) = spval
    data2dptr => this%qflx_snofrz_lyr_col(begc:endc,-nlevsno+1:0)
    call hist_addfld2d (fname='SNO_FRZ',  units='kg/m2/s', type2d='levsno', &
         avgflag='A', long_name='snow freezing rate in each snow layer', &
         ptr_col=data2dptr, c2l_scale_type='urbanf',no_snow_behavior=no_snow_normal, default='inactive')

    call hist_addfld2d (fname='SNO_FRZ_ICE',  units='mm/s', type2d='levsno', &
         avgflag='A', long_name='snow freezing rate in each snow layer (ice landunits only)', &
         ptr_col=data2dptr, c2l_scale_type='urbanf',no_snow_behavior=no_snow_normal, &
         l2g_scale_type='ice', default='inactive')

    this%qflx_h2osfc_to_ice_col(begc:endc) = spval
    call hist_addfld1d (fname='QH2OSFC_TO_ICE',  units='mm/s',  &
         avgflag='A', long_name='surface water converted to ice', &
         ptr_col=this%qflx_h2osfc_to_ice_col, default='inactive')

    this%qflx_prec_intr_patch(begp:endp) = spval
    call hist_addfld1d (fname='QINTR', units='mm/s',  &
         avgflag='A', long_name='interception', &
         ptr_patch=this%qflx_prec_intr_patch, set_lake=0._r8)

    this%qflx_prec_grnd_patch(begp:endp) = spval
    call hist_addfld1d (fname='QDRIP', units='mm/s',  &
         avgflag='A', long_name='throughfall', &
         ptr_patch=this%qflx_prec_grnd_patch, c2l_scale_type='urbanf')

    this%qflx_evap_soi_patch(begp:endp) = spval
    call hist_addfld1d (fname='QSOIL', units='mm/s',  &
         avgflag='A', long_name= 'Ground evaporation (soil/snow evaporation + soil/snow sublimation - dew)', &
         ptr_patch=this%qflx_evap_soi_patch, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QSOIL_ICE', units='mm/s',  &
         avgflag='A', long_name='Ground evaporation (ice landunits only)', &
         ptr_patch=this%qflx_evap_soi_patch, c2l_scale_type='urbanf', l2g_scale_type='ice')

    call hist_addfld2d (fname='QROOTSINK',  units='mm/s', type2d='levsoi', &
         avgflag='A', long_name='water flux from soil to root in each soil-layer', &
         ptr_col=this%qflx_rootsoi_col, set_spec=spval, l2g_scale_type='veg', default='inactive')

    this%qflx_evap_can_patch(begp:endp) = spval
    call hist_addfld1d (fname='QVEGE', units='mm/s',  &
         avgflag='A', long_name='canopy evaporation', &
         ptr_patch=this%qflx_evap_can_patch, set_lake=0._r8, c2l_scale_type='urbanf')

    this%qflx_tran_veg_patch(begp:endp) = spval
    call hist_addfld1d (fname='QVEGT', units='mm/s',  &
         avgflag='A', long_name='canopy transpiration', &
         ptr_patch=this%qflx_tran_veg_patch, set_lake=0._r8, c2l_scale_type='urbanf')

    this%qflx_ev_snow_patch(begp:endp) = spval
    call hist_addfld1d (fname='QSNOEVAP', units='mm/s',  &
         avgflag='A', long_name='evaporation from snow', &
         ptr_patch=this%qflx_ev_snow_patch, set_lake=0._r8, c2l_scale_type='urbanf')

    this%qflx_snowindunload_patch(begp:endp) = spval
    call hist_addfld1d (fname='QSNO_WINDUNLOAD', units='mm/s',  &
         avgflag='A', long_name='canopy snow wind unloading', &
         ptr_patch=this%qflx_snowindunload_patch, set_lake=0._r8, c2l_scale_type='urbanf')

    this%qflx_snotempunload_patch(begp:endp) = spval
    call hist_addfld1d (fname='QSNO_TEMPUNLOAD', units='mm/s',  &
         avgflag='A', long_name='canopy snow temp unloading', &
         ptr_patch=this%qflx_snotempunload_patch, set_lake=0._r8, c2l_scale_type='urbanf')

    this%qflx_snwcp_liq_col(begc:endc) = spval
    call hist_addfld1d (fname='QSNOCPLIQ', units='mm H2O/s', &
         avgflag='A', long_name='excess liquid h2o due to snow capping not including correction for land use change', &
         ptr_col=this%qflx_snwcp_liq_col, c2l_scale_type='urbanf')

    this%qflx_snwcp_ice_col(begc:endc) = spval
    call hist_addfld1d (fname='QSNWCPICE', units='mm H2O/s', &
         avgflag='A', long_name='excess solid h2o due to snow capping not including correction for land use change', &
         ptr_col=this%qflx_snwcp_ice_col, c2l_scale_type='urbanf')

    this%qflx_rain_grnd_patch(begp:endp) = spval
    call hist_addfld1d (fname='QFLX_RAIN_GRND', units='mm H2O/s', &
         avgflag='A', long_name='rain on ground after interception', &
         ptr_patch=this%qflx_rain_grnd_patch, default='inactive', c2l_scale_type='urbanf')

    this%qflx_snow_grnd_patch(begp:endp) = spval
    call hist_addfld1d (fname='QFLX_SNOW_GRND', units='mm H2O/s', &
         avgflag='A', long_name='snow on ground after interception', &
         ptr_patch=this%qflx_snow_grnd_patch, default='inactive', c2l_scale_type='urbanf')

    this%qflx_evap_grnd_patch(begp:endp) = spval
    call hist_addfld1d (fname='QFLX_EVAP_GRND', units='mm H2O/s', &
         avgflag='A', long_name='ground surface evaporation', &
         ptr_patch=this%qflx_evap_grnd_patch, default='inactive', c2l_scale_type='urbanf')

    this%qflx_evap_veg_patch(begp:endp) = spval
    call hist_addfld1d (fname='QFLX_EVAP_VEG', units='mm H2O/s', &
         avgflag='A', long_name='vegetation evaporation', &
         ptr_patch=this%qflx_evap_veg_patch, default='inactive', c2l_scale_type='urbanf')

    this%qflx_evap_tot_patch(begp:endp) = spval
    call hist_addfld1d (fname='QFLX_EVAP_TOT', units='mm H2O/s', &
         avgflag='A', long_name='qflx_evap_soi + qflx_evap_can + qflx_tran_veg', &
         ptr_patch=this%qflx_evap_tot_patch, c2l_scale_type='urbanf')

    this%qflx_dew_grnd_patch(begp:endp) = spval
    call hist_addfld1d (fname='QFLX_DEW_GRND', units='mm H2O/s', &
         avgflag='A', long_name='ground surface dew formation', &
         ptr_patch=this%qflx_dew_grnd_patch, c2l_scale_type='urbanf')

    this%qflx_sub_snow_patch(begp:endp) = spval
    call hist_addfld1d (fname='QFLX_SUB_SNOW', units='mm H2O/s', &
         avgflag='A', &
         long_name='sublimation rate from snow pack (also includes bare ice sublimation from glacier columns)', &
         ptr_patch=this%qflx_sub_snow_patch, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='QFLX_SUB_SNOW_ICE', units='mm H2O/s', &
         avgflag='A', &
         long_name='sublimation rate from snow pack (also includes bare ice sublimation from glacier columns) '// &
         '(ice landunits only)', &
         ptr_patch=this%qflx_sub_snow_patch, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')

    this%qflx_dew_snow_patch(begp:endp) = spval
    call hist_addfld1d (fname='QFLX_DEW_SNOW', units='mm H2O/s', &
         avgflag='A', long_name='surface dew added to snow pacK', &
         ptr_patch=this%qflx_dew_snow_patch, c2l_scale_type='urbanf')

    this%qflx_h2osfc_surf_col(begc:endc) = spval
    call hist_addfld1d (fname='QH2OSFC',  units='mm/s',  &
         avgflag='A', long_name='surface water runoff', &
         ptr_col=this%qflx_h2osfc_surf_col)

    this%qflx_drain_perched_col(begc:endc) = spval
    call hist_addfld1d (fname='QDRAI_PERCH',  units='mm/s',  &
         avgflag='A', long_name='perched wt drainage', &
         ptr_col=this%qflx_drain_perched_col, c2l_scale_type='urbanf')

    this%qflx_rsub_sat_col(begc:endc) = spval
    call hist_addfld1d (fname='QDRAI_XS',  units='mm/s',  &
         avgflag='A', long_name='saturation excess drainage', &
         ptr_col=this%qflx_rsub_sat_col, c2l_scale_type='urbanf')

    this%qflx_phs_neg_col(begc:endc) = spval
    call hist_addfld1d (fname='QPHSNEG',  units='mm/s',  &
         avgflag='A', long_name='net negative hydraulic redistribution flux', &
         ptr_col=this%qflx_phs_neg_col, default='inactive')

    ! As defined here, snow_sources - snow_sinks will equal the change in h2osno at any
    ! given time step but only if there is at least one snow layer (for all landunits 
    ! except lakes).  Also note that monthly average files of snow_sources and snow_sinks
    ! sinks must be weighted by number of days in the month to diagnose, for example, an 
    ! annual value of the change in h2osno. 

    this%snow_sources_col(begc:endc) = spval
    call hist_addfld1d (fname='SNOW_SOURCES',  units='mm/s',  &
         avgflag='A', long_name='snow sources (liquid water)', &
         ptr_col=this%snow_sources_col, c2l_scale_type='urbanf')

    this%snow_sinks_col(begc:endc) = spval
    call hist_addfld1d (fname='SNOW_SINKS',  units='mm/s',  &
         avgflag='A', long_name='snow sinks (liquid water)', &
         ptr_col=this%snow_sinks_col, c2l_scale_type='urbanf')
         
    this%AnnET(begc:endc) = spval
    call hist_addfld1d (fname='AnnET',  units='mm/s',  &
         avgflag='A', long_name='Annual ET', &
         ptr_col=this%AnnET, c2l_scale_type='urbanf', default='inactive')

  end subroutine InitHistory
  
  
  
   !-----------------------------------------------------------------------
    subroutine InitAccBuffer (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for all required module accumulated fields
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    !
    ! !USES 
    use clm_varcon  , only : spval
    use accumulMod  , only : init_accum_field
    !
    ! !ARGUMENTS:
    class(waterflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !---------------------------------------------------------------------

    if (use_fun) then
   
       call init_accum_field (name='AnnET', units='MM H2O/S', &
            desc='365-day running mean of total ET', accum_type='runmean', accum_period=-365, &
            subgrid_type='column', numlev=1, init_value=0._r8)

    end if

  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
    !
     subroutine InitAccVars (this, bounds)
    ! !DESCRIPTION:
    ! Initialize module variables that are associated with
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file 
    ! is read in and the accumulation buffer is obtained)
    !
    ! !USES 
    use accumulMod       , only : extract_accum_field
    use clm_time_manager , only : get_nstep
    !
    ! !ARGUMENTS:
    class(waterflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: begc, endc
    integer  :: nstep
    integer  :: ier
    real(r8), pointer :: rbufslp(:)  ! temporary
    !---------------------------------------------------------------------
    begc = bounds%begc; endc = bounds%endc

    ! Allocate needed dynamic memory for single level patch field
    allocate(rbufslp(begc:endc), stat=ier)

    ! Determine time step
    nstep = get_nstep()

    if (use_fun) then
       call extract_accum_field ('AnnET', rbufslp, nstep)
       this%AnnET(begc:endc) = rbufslp(begc:endc)
    end if

    deallocate(rbufslp)

  end subroutine InitAccVars
  
  
  !-----------------------------------------------------------------------
  subroutine UpdateAccVars (this, bounds)
    !
    ! USES
    use clm_time_manager, only : get_nstep
    use accumulMod      , only : update_accum_field, extract_accum_field
    !
    ! !ARGUMENTS:
    class(waterflux_type)                 :: this
    type(bounds_type)      , intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: g,c,p                     ! indices
    integer :: dtime                     ! timestep size [seconds]
    integer :: nstep                     ! timestep number
    integer :: ier                       ! error status
    integer :: begc, endc
    real(r8), pointer :: rbufslp(:)      ! temporary single level - patch level
    !---------------------------------------------------------------------

    begc = bounds%begc; endc = bounds%endc

    nstep = get_nstep()

    ! Allocate needed dynamic memory for single level patch field

    allocate(rbufslp(begc:endc), stat=ier)
    
    do c = begc,endc
       rbufslp(c) = this%qflx_evap_tot_col(c)
    end do
    if (use_fun) then
       ! Accumulate and extract AnnET (accumulates total ET as 365-day running mean)
       call update_accum_field  ('AnnET', rbufslp, nstep)
       call extract_accum_field ('AnnET', this%AnnET, nstep)
    
    end if

    deallocate(rbufslp)
    
  end subroutine UpdateAccVars


  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !USES:
    use landunit_varcon, only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    class(waterflux_type) :: this
    type(bounds_type) , intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l
    !-----------------------------------------------------------------------

    this%qflx_evap_grnd_patch(bounds%begp:bounds%endp) = 0.0_r8
    this%qflx_dew_grnd_patch (bounds%begp:bounds%endp) = 0.0_r8
    this%qflx_dew_snow_patch (bounds%begp:bounds%endp) = 0.0_r8

    this%qflx_evap_grnd_col(bounds%begc:bounds%endc) = 0.0_r8
    this%qflx_dew_grnd_col (bounds%begc:bounds%endc) = 0.0_r8
    this%qflx_dew_snow_col (bounds%begc:bounds%endc) = 0.0_r8

    this%qflx_phs_neg_col(bounds%begc:bounds%endc)   = 0.0_r8

    this%qflx_h2osfc_surf_col(bounds%begc:bounds%endc) = 0._r8
    this%qflx_snow_drain_col(bounds%begc:bounds%endc)  = 0._r8

    ! This variable only gets set in the hydrology filter; need to initialize it to 0 for
    ! the sake of columns outside this filter
    this%qflx_ice_runoff_xs_col(bounds%begc:bounds%endc) = 0._r8

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
    class(waterflux_type)            :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    ! needed for SNICAR
    call restartvar(ncid=ncid, flag=flag, varname='qflx_snofrz_lyr', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer ice freezing rate', units='kg m-2 s-1', &
         interpinic_flag='interp', readvar=readvar, data=this%qflx_snofrz_lyr_col)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize qflx_snofrz_lyr to zero
       this%qflx_snofrz_lyr_col(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    endif

    call restartvar(ncid=ncid, flag=flag, varname='qflx_snow_drain:qflx_snow_melt', xtype=ncd_double,  &
         dim1name='column', &
         long_name='drainage from snow column', units='mm/s', &
         interpinic_flag='interp', readvar=readvar, data=this%qflx_snow_drain_col)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize qflx_snow_drain to zero
       this%qflx_snow_drain_col(bounds%begc:bounds%endc) = 0._r8
    endif
    
    call this%qflx_liq_dynbal_dribbler%Restart(bounds, ncid, flag)
    call this%qflx_ice_dynbal_dribbler%Restart(bounds, ncid, flag)

  end subroutine Restart

end module WaterfluxType
