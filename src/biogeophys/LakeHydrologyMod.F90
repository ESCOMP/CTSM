module LakeHydrologyMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculation of Lake Hydrology. Full hydrology, aerosol deposition, etc. of snow layers is
  ! done. However, there is no infiltration, and the water budget is balanced with 
  ! qflx_qrgwl. Lake water mass is kept constant. The soil is simply maintained at
  ! volumetric saturation if ice melting frees up pore space. Likewise, if the water
  ! portion alone at some point exceeds pore capacity, it is reduced. This is consistent
  ! with the possibility of initializing the soil layer with excess ice.
  ! 
  ! If snow layers are present over an unfrozen lake, and the top layer of the lake
  ! is capable of absorbing the latent heat without going below freezing, 
  ! the snow-water is runoff and the latent heat is subtracted from the lake.
  !
  ! Minimum snow layer thickness for lakes has been increased to avoid instabilities with 30 min timestep.
  ! Also frost / dew is prevented from being added to top snow layers that have already melted during the phase change step.
  !
  ! ! USES
#include "shr_assert.h"
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use decompMod            , only : bounds_type
  use clm_varpar           , only : nlevsno, nlevgrnd, nlevsoi
  use ColumnType           , only : col                
  use PatchType            , only : patch                
  use atm2lndType          , only : atm2lnd_type
  use AerosolMod           , only : aerosol_type
  use EnergyFluxType       , only : energyflux_type
  use FrictionVelocityMod  , only : frictionvel_type
  use LakeStateType        , only : lakestate_type
  use SoilStateType        , only : soilstate_type
  use TemperatureType      , only : temperature_type
  use WaterType            , only : water_type
  use TotalWaterAndHeatMod , only : ComputeWaterMassLake
  use perf_mod             , only : t_startf, t_stopf
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: LakeHydrology              ! Calculates soil/snow hydrology
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: SumFlux_FluxesOntoGround    ! Compute summed fluxes onto ground, for bulk or one tracer

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine LakeHydrology(bounds, &
       num_lakec, filter_lakec, num_lakep, filter_lakep, &
       num_shlakesnowc, filter_shlakesnowc, num_shlakenosnowc, filter_shlakenosnowc, &
       scf_method, water_inst, &
       atm2lnd_inst, temperature_inst, soilstate_inst, &
       energyflux_inst, aerosol_inst, lakestate_inst, topo_inst)
    !
    ! !DESCRIPTION:
    ! WARNING: This subroutine assumes lake columns have one and only one pft.
    !
    ! Sequence is:
    !  LakeHydrology:
    !    Do needed tasks from CanopyHydrology, Biogeophysics2, & top of SoilHydrology.
    !    -> SnowWater:             change of snow mass and snow water onto soil
    !    -> SnowCompaction:        compaction of snow layers
    !    -> CombineSnowLayers:     combine snow layers that are thinner than minimum
    !    -> DivideSnowLayers:      subdivide snow layers that are thicker than maximum
    !
    !    Add water to soil if melting has left it with open pore space.
    !    If snow layers are found above a lake with unfrozen top layer, whose top
    !    layer has enough heat to melt all the snow ice without freezing, do so
    !    and eliminate the snow layers.
    !    Cleanup and do water balance.
    !
    ! !USES:
    use clm_varcon      , only : denh2o, denice, spval, hfus, tfrz, cpliq, cpice
    use clm_varctl      , only : iulog
    use clm_time_manager, only : get_step_size_real
    use SnowHydrologyMod, only : UpdateQuantitiesForNewSnow, InitializeExplicitSnowPack
    use SnowHydrologyMod, only : SnowCompaction, CombineSnowLayers, SnowWater
    use SnowHydrologyMod, only : ZeroEmptySnowLayers, BuildSnowFilter, SnowCapping
    use SnowHydrologyMod, only : DivideSnowLayers
    use TopoMod         , only : topo_type
    use SnowCoverFractionBaseMod, only : snow_cover_fraction_base_type
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds  
    integer                , intent(in)    :: num_lakec               ! number of column lake points in column filter
    integer                , intent(in)    :: filter_lakec(:)         ! column filter for lake points
    integer                , intent(in)    :: num_lakep               ! number of pft lake points in column filter
    integer                , intent(in)    :: filter_lakep(:)         ! patch filter for lake points
    integer                , intent(out)   :: num_shlakesnowc         ! number of column snow points
    integer                , intent(out)   :: filter_shlakesnowc(:)   ! column filter for snow points
    integer                , intent(out)   :: num_shlakenosnowc       ! number of column non-snow points
    integer                , intent(out)   :: filter_shlakenosnowc(:) ! column filter for non-snow points
    class(snow_cover_fraction_base_type), intent(in) :: scf_method
    type(water_type)       , intent(inout) :: water_inst
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(temperature_type) , intent(inout) :: temperature_inst
    type(soilstate_type)   , intent(in)    :: soilstate_inst
    type(energyflux_type)  , intent(inout) :: energyflux_inst
    type(aerosol_type)     , intent(inout) :: aerosol_inst
    type(lakestate_type)   , intent(inout) :: lakestate_inst
    class(topo_type)   , intent(in)    :: topo_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: p,fp,g,l,c,j,fc,jtop                            ! indices
    integer  :: i                                               ! index of water tracer or bulk
    real(r8) :: dtime                                           ! land model time step (sec)
    real(r8) :: dz_snowf                                        ! layer thickness rate change due to precipitation [mm/s]
    real(r8) :: bifall(bounds%begc:bounds%endc)                 ! bulk density of newly fallen dry snow [kg/m3]
    real(r8) :: fracsnow(bounds%begp:bounds%endp)               ! frac of precipitation that is snow
    real(r8) :: fracrain(bounds%begp:bounds%endp)               ! frac of precipitation that is rain
    real(r8) :: qflx_evap_soi_lim                               ! temporary evap_soi limited by top snow layer content [mm/s]
    real(r8) :: h2osno_temp                                     ! temporary h2osno [kg/m^2]
    real(r8) :: h2osno_total(bounds%begc:bounds%endc)           ! total snow water (mm H2O)
    real(r8) :: sumsnowice(bounds%begc:bounds%endc)             ! sum of snow ice if snow layers found above unfrozen lake [kg/m&2]
    logical  :: unfrozen(bounds%begc:bounds%endc)               ! true if top lake layer is unfrozen with snow layers above
    real(r8) :: heatrem                                         ! used in case above [J/m^2]
    real(r8) :: heatsum(bounds%begc:bounds%endc)                ! used in case above [J/m^2]
    real(r8) :: qflx_dew_minus_sub_snow                         ! qflx_soliddew_to_top_layer - qflx_solidevap_from_top_layer [mm/s]
    real(r8), parameter :: frac_sno_small = 1.e-6_r8            ! small value of frac_sno used when initiating a snow pack due to frost
    real(r8), parameter :: snow_bd = 250._r8                    ! assumed snow bulk density (for lakes w/out resolved snow layers) [kg/m^3]
    ! Should only be used for frost below.
    !-----------------------------------------------------------------------

    associate( &
         b_waterstate_inst  => water_inst%waterstatebulk_inst, &
         b_waterdiagnostic_inst => water_inst%waterdiagnosticbulk_inst, &
         b_waterbalance_inst => water_inst%waterbalancebulk_inst, &
         b_waterflux_inst => water_inst%waterfluxbulk_inst, &
         b_wateratm2lnd_inst => water_inst%wateratm2lndbulk_inst, &
         
         begc => bounds%begc, &
         endc => bounds%endc  &
         )

    ! TODO(wjs, 2019-08-07) Once this subroutine has been modularized so that it is just
    ! wrapper that calls the actual science routines, we should be able to get rid of
    ! this associate statement.
    associate(                                                            & 
         pcolumn              =>  patch%column                            , & ! Input:  [integer  (:)   ]  pft's column index                       
         pgridcell            =>  patch%gridcell                          , & ! Input:  [integer  (:)   ]  pft's gridcell index                     
         cgridcell            =>  col%gridcell                          , & ! Input:  [integer  (:)   ]  column's gridcell                        
         clandunit            =>  col%landunit                          , & ! Input:  [integer  (:)   ]  column's landunit                        
         dz_lake              =>  col%dz_lake                           , & ! Input:  [real(r8) (:,:) ]  layer thickness for lake (m)          
         z                    =>  col%z                                 , & ! Input:  [real(r8) (:,:) ]  layer depth  (m)                      
         dz                   =>  col%dz                                , & ! Input:  [real(r8) (:,:) ]  layer thickness depth (m)             
         zi                   =>  col%zi                                , & ! Input:  [real(r8) (:,:) ]  interface depth (m)                   
         snl                  =>  col%snl                               , & ! Input:  [integer  (:)   ]  number of snow layers                    
         
         forc_rain            =>  b_wateratm2lnd_inst%forc_rain_downscaled_col , & ! Input:  [real(r8) (:)   ]  rain rate [mm/s]                        
         forc_snow            =>  b_wateratm2lnd_inst%forc_snow_downscaled_col , & ! Input:  [real(r8) (:)   ]  snow rate [mm/s]                        
         qflx_floodg          =>  b_wateratm2lnd_inst%forc_flood_grc           , & ! Input:  [real(r8) (:)   ]  gridcell flux of flood water from RTM   
         
         watsat               =>  soilstate_inst%watsat_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         
         t_lake               =>  temperature_inst%t_lake_col           , & ! Input:  [real(r8) (:,:) ]  lake temperature (Kelvin)             
         t_grnd               =>  temperature_inst%t_grnd_col           , & ! Input:  [real(r8) (:)   ]  ground temperature (Kelvin)             
         t_soisno             =>  temperature_inst%t_soisno_col         , & ! Output: [real(r8) (:,:) ]  snow temperature (Kelvin)             
         dTdz_top             =>  temperature_inst%dTdz_top_col         , & ! Output: [real(r8) (:)   ]  temperature gradient in top layer K m-1] !TOD 
         snot_top             =>  temperature_inst%snot_top_col         , & ! Output: [real(r8) (:)   ]  snow temperature in top layer [K]  !TODO
         t_sno_mul_mss        =>  temperature_inst%t_sno_mul_mss_col     , & ! Output: [real(r8) (:)   ]  col snow temperature multiplied by layer mass, layer sum (K * kg/m2) 
         
         begwb                =>  b_waterbalance_inst%begwb_col             , & ! Input:  [real(r8) (:)   ]  water mass begining of the time step    
         endwb                =>  b_waterbalance_inst%endwb_col             , & ! Output: [real(r8) (:)   ]  water mass end of the time step         
         snw_rds              =>  b_waterdiagnostic_inst%snw_rds_col           , & ! Output: [real(r8) (:,:) ]  effective snow grain radius (col,lyr) [microns, m^-6] 
         snw_rds_top          =>  b_waterdiagnostic_inst%snw_rds_top_col       , & ! Output: [real(r8) (:)   ]  effective snow grain size, top layer [microns] 
         h2osno_top           =>  b_waterdiagnostic_inst%h2osno_top_col        , & ! Output: [real(r8) (:)   ]  mass of snow in top layer [kg]    
         sno_liq_top          =>  b_waterdiagnostic_inst%sno_liq_top_col       , & ! Output: [real(r8) (:)   ]  liquid water fraction in top snow layer [frc] 
         frac_sno             =>  b_waterdiagnostic_inst%frac_sno_col          , & ! Output: [real(r8) (:)   ]
         frac_sno_eff         =>  b_waterdiagnostic_inst%frac_sno_eff_col      , & ! Output: [real(r8) (:)   ]  needed for snicar code                  
         frac_iceold          =>  b_waterdiagnostic_inst%frac_iceold_col       , & ! Output: [real(r8) (:,:) ]  fraction of ice relative to the tot water
         snow_depth           =>  b_waterdiagnostic_inst%snow_depth_col        , & ! Output: [real(r8) (:)   ]  snow height (m)                         
         h2osno_no_layers     => b_waterstate_inst%h2osno_no_layers_col        , & ! Output: [real(r8) (:)   ]  snow that is not resolved into layers (kg/m2)
         snowice              =>  b_waterdiagnostic_inst%snowice_col           , & ! Output: [real(r8) (:)   ]  average snow ice lens                   
         snowliq              =>  b_waterdiagnostic_inst%snowliq_col           , & ! Output: [real(r8) (:)   ]  average snow liquid water               
         h2osoi_ice           =>  b_waterstate_inst%h2osoi_ice_col        , & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)                      
         h2osoi_liq           =>  b_waterstate_inst%h2osoi_liq_col        , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)                  
         h2osoi_vol           =>  b_waterstate_inst%h2osoi_vol_col        , & ! Output: [real(r8) (:,:) ]  volumetric soil water [m3/m3]         
         
         qflx_floodc          =>  b_waterflux_inst%qflx_floodc_col        , & ! Output: [real(r8) (:)   ]  column flux of flood water from RTM     
         qflx_liq_grnd        =>  b_waterflux_inst%qflx_liq_grnd_col      , & ! Output: [real(r8) (:)   ]  liquid on ground after interception (mm H2O/s) [+]
         qflx_evap_tot        =>  b_waterflux_inst%qflx_evap_tot_patch    , & ! Output: [real(r8) (:)   ]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg
         qflx_evap_soi        =>  b_waterflux_inst%qflx_evap_soi_patch    , & ! Output: [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)
         qflx_solidevap_from_top_layer => b_waterflux_inst%qflx_solidevap_from_top_layer_patch, & ! Output: [real(r8) (:)   ]  rate of ice evaporated from top soil or snow layer (sublimation) (mm H2O /s) [+]
         qflx_liqevap_from_top_layer   => b_waterflux_inst%qflx_liqevap_from_top_layer_patch  , & ! Output: [real(r8) (:)   ]  rate of liquid water evaporated from top soil or snow layer (mm H2O/s) [+]
         qflx_soliddew_to_top_layer    =>  b_waterflux_inst%qflx_soliddew_to_top_layer_patch  , & ! Output: [real(r8) (:)   ]  rate of solid water deposited on top soil or snow layer (frost) (mm H2O /s) [+]
         qflx_liqdew_to_top_layer      => b_waterflux_inst%qflx_liqdew_to_top_layer_patch     , & ! Output: [real(r8) (:)   ]  rate of liquid water deposited on top soil or snow layer (dew) (mm H2O /s) [+]
         qflx_snomelt         =>  b_waterflux_inst%qflx_snomelt_col       , & ! Output: [real(r8) (:)   ]  snow melt (mm H2O /s)
         qflx_snomelt_lyr     =>  b_waterflux_inst%qflx_snomelt_lyr_col   , & ! Output: [real(r8) (:)   ]  snow melt in each layer (mm H2O /s)
         qflx_liqevap_from_top_layer_col   => b_waterflux_inst%qflx_liqevap_from_top_layer_col  , & ! Output: [real(r8) (:)   ]  rate of liquid water evaporated from top soil or snow layer (mm H2O/s) [+]
         qflx_liqdew_to_top_layer_col      => b_waterflux_inst%qflx_liqdew_to_top_layer_col     , & ! Output: [real(r8) (:)   ]  rate of liquid water deposited on top soil or snow layer (dew) (mm H2O /s) [+]
         qflx_soliddew_to_top_layer_col    =>  b_waterflux_inst%qflx_soliddew_to_top_layer_col  , & ! Output: [real(r8) (:)   ]  rate of solid water deposited on top soil or snow layer (frost) (mm H2O /s) [+]
         qflx_solidevap_from_top_layer_col => b_waterflux_inst%qflx_solidevap_from_top_layer_col, & ! Output: [real(r8) (:)   ]  rate of ice evaporated from top soil or snow layer (sublimation) (mm H2O /s) [+]
         qflx_ev_snow         =>  b_waterflux_inst%qflx_ev_snow_patch     , & ! Output: [real(r8) (:)   ]  evaporation flux from snow (mm H2O/s) [+ to atm]
         qflx_ev_snow_col     =>  b_waterflux_inst%qflx_ev_snow_col       , & ! Output: [real(r8) (:)   ]  evaporation flux from snow (mm H2O/s) [+ to atm]
         qflx_evap_tot_col    =>  b_waterflux_inst%qflx_evap_tot_col      , & ! Output: [real(r8) (:)   ]  pft quantity averaged to the column (assuming one pft)
         qflx_snwcp_ice       =>  b_waterflux_inst%qflx_snwcp_ice_col     , & ! Output: [real(r8) (:)   ]  excess solid h2o due to snow capping (outgoing) (mm H2O /s) [+]
         qflx_snwcp_discarded_ice => b_waterflux_inst%qflx_snwcp_discarded_ice_col, & ! Input: [real(r8) (:)   ]  excess solid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s) [+]
         qflx_snwcp_discarded_liq => b_waterflux_inst%qflx_snwcp_discarded_liq_col, & ! Input: [real(r8) (:)   ]  excess liquid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s) [+]
         qflx_drain_perched   =>  b_waterflux_inst%qflx_drain_perched_col , & ! Output: [real(r8) (:)   ]  perched wt sub-surface runoff (mm H2O /s) !TODO - move this to somewhere else
         qflx_snow_drain      =>  b_waterflux_inst%qflx_snow_drain_col    , & ! Output: [real(r8) (:)   ]  drainage from snow pack                          
         qflx_rsub_sat        =>  b_waterflux_inst%qflx_rsub_sat_col      , & ! Output: [real(r8) (:)   ]  soil saturation excess [mm h2o/s]        
         qflx_surf            =>  b_waterflux_inst%qflx_surf_col          , & ! Output: [real(r8) (:)   ]  surface runoff (mm H2O /s)              
         qflx_drain           =>  b_waterflux_inst%qflx_drain_col         , & ! Output: [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)          
         qflx_infl            =>  b_waterflux_inst%qflx_infl_col          , & ! Output: [real(r8) (:)   ]  infiltration (mm H2O /s)                
         qflx_qrgwl           =>  b_waterflux_inst%qflx_qrgwl_col         , & ! Output: [real(r8) (:)   ]  qflx_surf at glaciers, wetlands, lakes  
         qflx_runoff          =>  b_waterflux_inst%qflx_runoff_col        , & ! Output: [real(r8) (:)   ]  total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
         qflx_ice_runoff_snwcp => b_waterflux_inst%qflx_ice_runoff_snwcp_col, & ! Output: [real(r8) (:)] solid runoff from snow capping (mm H2O /s)
         qflx_rain_plus_snomelt => b_waterflux_inst%qflx_rain_plus_snomelt_col , & ! Output: [real(r8) (:)   ] rain plus snow melt falling on the soil (mm/s)
         qflx_top_soil        =>  b_waterflux_inst%qflx_top_soil_col      , & ! Output: [real(r8) (:)   ]  net water input into soil from top (mm/s)
         
         eflx_snomelt         =>  energyflux_inst%eflx_snomelt_col      , & ! Output: [real(r8) (:)   ]  snow melt heat flux (W/m**2)
         eflx_sh_tot          =>  energyflux_inst%eflx_sh_tot_patch     , & ! Output: [real(r8) (:)   ]  total sensible heat flux (W/m**2) [+ to atm]
         eflx_sh_grnd         =>  energyflux_inst%eflx_sh_grnd_patch    , & ! Output: [real(r8) (:)   ]  sensible heat flux from ground (W/m**2) [+ to atm]
         eflx_soil_grnd       =>  energyflux_inst%eflx_soil_grnd_patch  , & ! Output: [real(r8) (:)   ]  heat flux into snow / lake (W/m**2) [+ = into soil]
         eflx_gnet            =>  energyflux_inst%eflx_gnet_patch       , & ! Output: [reay(r8) (:)   ]  net heat flux into ground (W/m**2)      
         eflx_grnd_lake       =>  energyflux_inst%eflx_grnd_lake_patch  , & ! Output: [real(r8) (:)   ]  net heat flux into lake / snow surface, excluding light transmission (W/m**2)
         
         lake_icefrac         =>  lakestate_inst%lake_icefrac_col         & ! Output: [real(r8) (:,:) ]  mass fraction of lake layer that is frozen
         )

    ! BUG(wjs, 2019-07-12, ESCOMP/ctsm#762) This is needed so that temporary tracer
    ! consistency checks later in this routine pass. Remove this block once code before
    ! this point is fully tracerized.
    if (water_inst%DoConsistencyCheck()) then
       call water_inst%ResetCheckedTracers(bounds)
       call water_inst%TracerConsistencyCheck(bounds, 'start of LakeHydrology')
    end if

    ! Determine step size
    dtime = get_step_size_real()

    ! Compute "summed" (really just copies here) fluxes onto "ground" (really the lake
    ! surface here), for bulk water and each tracer. (Subroutine name mimics the one in
    ! CanopyHydrologyMod.)
    do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
       associate(w => water_inst%bulk_and_tracers(i))
       call SumFlux_FluxesOntoGround(bounds, &
            num_lakec, filter_lakec, &
            ! Inputs
            forc_snow      = w%wateratm2lnd_inst%forc_snow_downscaled_col(begc:endc), &
            forc_rain      = w%wateratm2lnd_inst%forc_rain_downscaled_col(begc:endc), &
            ! Outputs
            qflx_snow_grnd = w%waterflux_inst%qflx_snow_grnd_col(begc:endc), &
            qflx_liq_grnd  = w%waterflux_inst%qflx_liq_grnd_col(begc:endc))
       end associate
    end do

    ! Determine snow height and snow water

    call UpdateQuantitiesForNewSnow(bounds, num_lakec, filter_lakec, &
         scf_method, atm2lnd_inst, water_inst)

    call InitializeExplicitSnowPack(bounds, num_lakec, filter_lakec, &
         atm2lnd_inst, temperature_inst, aerosol_inst, water_inst)

    ! TODO(wjs, 2019-08-01) Eventually move this down, merging this with later tracer
    ! consistency checks. If/when we remove calls to TracerConsistencyCheck from this
    ! module, remember to also remove 'use perf_mod' at the top.
    if (water_inst%DoConsistencyCheck()) then
       call t_startf("tracer_consistency_check")
       call water_inst%TracerConsistencyCheck(bounds, 'after initial snow stuff in LakeHydrology')
       call t_stopf("tracer_consistency_check")
    end if

    ! Calculate sublimation and dew, adapted from HydrologyLake and Biogeophysics2.

    do fp = 1,num_lakep
       p = filter_lakep(fp)
       c = pcolumn(p)
       jtop = snl(c)+1

       qflx_liqevap_from_top_layer(p)   = 0._r8
       qflx_solidevap_from_top_layer(p) = 0._r8
       qflx_soliddew_to_top_layer(p)    = 0._r8
       qflx_liqdew_to_top_layer(p)      = 0._r8
       qflx_ev_snow(p)                  = qflx_evap_soi(p)

       if (jtop <= 0) then ! snow layers
          j = jtop
          ! Assign ground evaporation to sublimation from soil ice or to dew
          ! on snow or ground

          if (qflx_evap_soi(p) >= 0._r8) then
             ! for evaporation partitioning between liquid evap and ice sublimation, 
             ! use the ratio of liquid to (liquid+ice) in the top layer to determine split
             ! Since we're not limiting evap over lakes, but still can't remove more from top
             ! snow layer than there is there, create temp. limited evap_soi.
             qflx_evap_soi_lim = min(qflx_evap_soi(p), (h2osoi_liq(c,j)+h2osoi_ice(c,j))/dtime)
             qflx_ev_snow(p) = qflx_evap_soi_lim
             if ((h2osoi_liq(c,j)+h2osoi_ice(c,j)) > 0._r8) then
                qflx_liqevap_from_top_layer(p) = max(qflx_evap_soi_lim*(h2osoi_liq(c,j)/ &
                     (h2osoi_liq(c,j)+h2osoi_ice(c,j))), 0._r8)
             else
                qflx_liqevap_from_top_layer(p) = 0._r8
             end if
             qflx_solidevap_from_top_layer(p) = qflx_evap_soi_lim - qflx_liqevap_from_top_layer(p)     
          else
             ! if (t_grnd(c) < tfrz) then
             ! Causes rare blowup when thin snow layer should completely melt and has a high temp after thermal physics,
             ! but then is not eliminated in SnowHydrology because of this added frost. Also see below removal of
             ! completely melted single snow layer.
             if (t_grnd(c) < tfrz .and. t_soisno(c,j) < tfrz) then
                qflx_soliddew_to_top_layer(p) = abs(qflx_evap_soi(p))
                ! If top layer is only snow layer, SnowHydrology won't eliminate it if dew is added.
             else if (j < 0 .or. (t_grnd(c) == tfrz .and. t_soisno(c,j) == tfrz)) then
                qflx_liqdew_to_top_layer(p) = abs(qflx_evap_soi(p))
             end if
          end if

       else ! No snow layers
          if (qflx_evap_soi(p) >= 0._r8) then
             ! Sublimation: do not allow for more sublimation than there is snow
             ! after melt.  Remaining surface evaporation used for infiltration.
             qflx_solidevap_from_top_layer(p) = min(qflx_evap_soi(p), h2osno_no_layers(c)/dtime)
             qflx_liqevap_from_top_layer(p) = qflx_evap_soi(p) - qflx_solidevap_from_top_layer(p)
          else
             if (t_grnd(c) < tfrz-0.1_r8) then
                qflx_soliddew_to_top_layer(p) = abs(qflx_evap_soi(p))
             else
                qflx_liqdew_to_top_layer(p) = abs(qflx_evap_soi(p))
             end if
          end if

          ! Update snow pack for dew & sub.

          h2osno_temp = h2osno_no_layers(c)
          qflx_dew_minus_sub_snow = -qflx_solidevap_from_top_layer(p)+qflx_soliddew_to_top_layer(p)
          h2osno_no_layers(c) = h2osno_no_layers(c) + qflx_dew_minus_sub_snow*dtime
          h2osno_no_layers(c) = max(h2osno_no_layers(c), 0._r8)
          if (qflx_dew_minus_sub_snow > 0._r8) then
             ! If we're accumulating snow from dew, then ensure that we have at least a
             ! small, non-zero frac_sno. (It complicates the code too much to call
             ! UpdateSnowDepthAndFrac for this purpose - see
             ! <https://github.com/ESCOMP/CTSM/issues/827#issuecomment-546163067>.)
             if (frac_sno(c) <= 0._r8) then
                frac_sno(c) = frac_sno_small
             end if
          else if (qflx_dew_minus_sub_snow < 0._r8) then
             ! If we're losing snow from sublimation, and this has caused the snow pack
             ! to completely vanish, then ensure that frac_sno is reset to 0.
             if (h2osno_no_layers(c) == 0._r8) then
                frac_sno(c) = 0._r8
             end if
          end if
          if (h2osno_temp > 0._r8) then
             ! Assume that snow bulk density remains the same as before
             ! NOTE (SSR, 2023-11-08): Small h2osno_temp can cause unrealistically high snow depths: see https://github.com/ESCOMP/CTSM/issues/2227. Suggested fix there is to replace this line with
             !     snow_depth(c) = h2osno_no_layers(c) * min (snow_depth(c)/h2osno_temp, 1._r8/50._r8)
             ! where 50 kg/m3 is suggested as a lower limit for snow density.
             ! As this bug seemingly has never been encountered in CTSM, we are not yet implementing the fix.
             snow_depth(c) = snow_depth(c) * h2osno_no_layers(c) / h2osno_temp
          else
             ! Assume a constant snow bulk density = 250.
             snow_depth(c) = h2osno_no_layers(c)/snow_bd
          end if

       end if
    end do

    ! Since frac_sno may have been updated above, recalculate frac_sno_eff accordingly
    call scf_method%CalcFracSnoEff(bounds, num_lakec, filter_lakec, &
         ! Inputs
         lun_itype_col = col%lun_itype(begc:endc), &
         urbpoi        = col%urbpoi(begc:endc), &
         frac_sno      = frac_sno(begc:endc), &
         ! Outputs
         frac_sno_eff  = frac_sno_eff(begc:endc))

    ! patch averages must be done here -- BEFORE SNOW CALCULATIONS AS THEY USE IT.
    ! for output to history tape and other uses
    ! (note that pft2col is called before LakeHydrology, so we can't use that routine
    ! to do these column -> pft averages)
    do fp = 1,num_lakep
       p = filter_lakep(fp)
       c = pcolumn(p)

       qflx_evap_tot_col(c)  = qflx_evap_tot(p)
       qflx_liqevap_from_top_layer_col(c)   = qflx_liqevap_from_top_layer(p)
       qflx_liqdew_to_top_layer_col(c)      = qflx_liqdew_to_top_layer(p)
       qflx_soliddew_to_top_layer_col(c)    = qflx_soliddew_to_top_layer(p)
       qflx_solidevap_from_top_layer_col(c) = qflx_solidevap_from_top_layer(p)
       qflx_ev_snow_col(c)                  = qflx_ev_snow(p)
    enddo

    ! BUG(wjs, 2019-07-12, ESCOMP/ctsm#762) This is needed so that we can test the
    ! tracerization of the following snow stuff without having tracerized all of the above
    ! code. Remove this block once code before this point is fully tracerized.
    if (water_inst%DoConsistencyCheck()) then
       call water_inst%ResetCheckedTracers(bounds)
       call water_inst%TracerConsistencyCheck(bounds, 'before main snow code in LakeHydrology')
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Determine initial snow/no-snow filters (will be modified possibly by
    ! routines CombineSnowLayers and DivideSnowLayers below)

    call BuildSnowFilter(bounds, num_lakec, filter_lakec, &
         num_shlakesnowc, filter_shlakesnowc, num_shlakenosnowc, filter_shlakenosnowc)

    ! Determine the change of snow mass and the snow water onto soil

    call SnowWater(bounds, &
         num_shlakesnowc, filter_shlakesnowc, num_shlakenosnowc, filter_shlakenosnowc, &
         atm2lnd_inst, aerosol_inst, water_inst)

    call SnowCapping(bounds, num_lakec, filter_lakec, num_shlakesnowc, filter_shlakesnowc, &
         topo_inst, aerosol_inst, water_inst)

    ! Natural compaction and metamorphosis.

    call SnowCompaction(bounds, num_shlakesnowc, filter_shlakesnowc, &
         scf_method, &
         temperature_inst, b_waterstate_inst, b_waterdiagnostic_inst, atm2lnd_inst)

    ! Combine thin snow elements

    call CombineSnowLayers(bounds, num_shlakesnowc, filter_shlakesnowc, &
         aerosol_inst, temperature_inst, water_inst)

    ! Divide thick snow elements

    call DivideSnowLayers(bounds, num_shlakesnowc, filter_shlakesnowc, &
         aerosol_inst, temperature_inst, water_inst, is_lake=.true.)

    ! Set empty snow layers to zero
    call ZeroEmptySnowLayers(bounds, num_shlakesnowc, filter_shlakesnowc, &
         col, water_inst, temperature_inst)

    ! TODO(wjs, 2019-09-16) Eventually move this down, merging this with later tracer
    ! consistency checks. If/when we remove calls to TracerConsistencyCheck from this
    ! module, remember to also remove 'use perf_mod' at the top.
    if (water_inst%DoConsistencyCheck()) then
       call t_startf("tracer_consistency_check")
       call water_inst%TracerConsistencyCheck(bounds, 'LakeHydrology: after main snow code')
       call t_stopf("tracer_consistency_check")
    end if

    ! Recompute h2osno_total for possible updates in the above snow routines
    call b_waterstate_inst%CalculateTotalH2osno(bounds, num_lakec, filter_lakec, &
         caller = 'LakeHydrology-2', &
         h2osno_total = h2osno_total(bounds%begc:bounds%endc))

    ! Determine soil hydrology
    ! Here this consists only of making sure that soil is saturated even as it melts and
    ! pore space opens up. Conversely, if excess ice is melting and the liquid water exceeds the
    ! saturation value, then remove water.

    do j = 1,nlevsoi  !nlevgrnd
       ! changed to nlevsoi on 8/11/10 to make consistent with non-lake bedrock
       do fc = 1, num_lakec
          c = filter_lakec(fc)

          h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)
          ! Could have changed during phase change! (Added 8/11/10)

          if (h2osoi_vol(c,j) < watsat(c,j)) then
             h2osoi_liq(c,j) = (watsat(c,j)*dz(c,j) - h2osoi_ice(c,j)/denice)*denh2o
             ! h2osoi_vol will be updated below, and this water addition will come from qflx_qrgwl
          else if (h2osoi_liq(c,j) > watsat(c,j)*denh2o*dz(c,j)) then
             h2osoi_liq(c,j) = watsat(c,j)*denh2o*dz(c,j)
             ! Another way to do this would be: if h2osoi_vol > watsat then remove min(h2osoi_liq,
             !(h2osoi_vol-watsat)*dz*denh2o) from h2osoi_liq.  The question is whether the excess ice
             ! melts first or last (or simultaneously) to the pore ice.  Because excess ice is often in chunks,
             ! requiring greater convergence of heat to melt, assume it melts last.
             ! This will also improve the initialization behavior or in an occasionally warm year, the excess ice
             ! won't start going away if a layer is briefly at freezing.

             ! Allow up to 10% excess ice over watsat in refreezing soil,
             ! e.g. heaving soil.  (As with > 10% excess ice modeling, and for the lake water,
             ! the thermal conductivity will be adjusted down to compensate for the fact that the nominal dz is smaller
             ! than the real soil volume.)  The current solution is consistent but perhaps unrealistic in real soils,
             ! where slow drainage may occur during freezing; drainage is only assumed to occur here when >10% excess
             ! ice melts. The latter is more likely to be permanent rather than seasonal anyway. Attempting to remove the
             ! ice volume after some has already frozen during the timestep would not conserve energy unless this were
             ! incorporated into the ice stream.

          end if

       end do
    end do
!!!!!!!!!!

    ! Check for single completely unfrozen snow layer over lake.  Modeling this ponding is unnecessary and
    ! can cause instability after the timestep when melt is completed, as the temperature after melt can be
    ! excessive because the fluxes were calculated with a fixed ground temperature of freezing, but the 
    ! phase change was unable to restore the temperature to freezing.

    do fp = 1, num_lakep
       p = filter_lakep(fp)
       c = pcolumn(p)

       j = 0

       if (snl(c) == -1) then 
          if (h2osoi_ice(c,j) > 0._r8 .and. t_soisno(c,j) > tfrz) then

             ! Take extra heat of layer and release to sensible heat in order 
             ! to maintain energy conservation.
             heatrem           = (cpliq*h2osoi_liq(c,j))*(t_soisno(c,j) - tfrz)
             t_soisno(c,j)     = tfrz
             eflx_sh_tot(p)    = eflx_sh_tot(p) + heatrem/dtime
             eflx_sh_grnd(p)   = eflx_sh_grnd(p) + heatrem/dtime
             eflx_soil_grnd(p) = eflx_soil_grnd(p) - heatrem/dtime
             eflx_gnet(p)      = eflx_gnet(p) - heatrem/dtime
             eflx_grnd_lake(p) = eflx_grnd_lake(p) - heatrem/dtime
          else if (h2osoi_ice(c,j) == 0._r8) then
             ! Remove layer
             ! Take extra heat of layer and release to sensible heat in order 
             ! to maintain energy conservation.
             heatrem             = cpliq*h2osoi_liq(c,j)*(t_soisno(c,j) - tfrz)
             eflx_sh_tot(p)      = eflx_sh_tot(p) + heatrem/dtime
             eflx_sh_grnd(p)     = eflx_sh_grnd(p) + heatrem/dtime
             eflx_soil_grnd(p)   = eflx_soil_grnd(p) - heatrem/dtime
             eflx_gnet(p)        = eflx_gnet(p) - heatrem/dtime
             eflx_grnd_lake(p)   = eflx_grnd_lake(p) - heatrem/dtime
             qflx_snow_drain(c)  = qflx_snow_drain(c) + h2osno_total(c)/dtime
             snl(c)              = 0
             h2osno_no_layers(c) = 0._r8
             h2osno_total(c)     = 0._r8
             snow_depth(c)       = 0._r8
             ! Rest of snow layer book-keeping will be done below.
          else
             eflx_grnd_lake(p) = eflx_gnet(p)
          end if
       else
          eflx_grnd_lake(p) = eflx_gnet(p)
       end if
    end do

    ! Check for snow layers above lake with unfrozen top layer.  Mechanically,
    ! the snow will fall into the lake and melt or turn to ice.  If the top layer has
    ! sufficient heat to melt the snow without freezing, then that will be done.
    ! Otherwise, the top layer will undergo freezing, but only if the top layer will
    ! not freeze completely.  Otherwise, let the snow layers persist and melt by diffusion.

    do fc = 1, num_lakec
       c = filter_lakec(fc)

       if (t_lake(c,1) > tfrz .and. lake_icefrac(c,1) == 0._r8 .and. snl(c) < 0) then
          unfrozen(c) = .true.
       else
          unfrozen(c) = .false.
       end if
    end do

    do j = -nlevsno+1,0
       do fc = 1, num_lakec
          c = filter_lakec(fc)

          if (unfrozen(c)) then
             if (j == -nlevsno+1) then
                sumsnowice(c) = 0._r8
                heatsum(c) = 0._r8
             end if
             if (j >= snl(c)+1) then
                sumsnowice(c) = sumsnowice(c) + h2osoi_ice(c,j)
                heatsum(c) = heatsum(c) + h2osoi_ice(c,j)*cpice*(tfrz - t_soisno(c,j)) &
                     + h2osoi_liq(c,j)*cpliq*(tfrz - t_soisno(c,j))
             end if
          end if
       end do
    end do

    do fc = 1, num_lakec
       c = filter_lakec(fc)

       if (unfrozen(c)) then
          heatsum(c) = heatsum(c) + sumsnowice(c)*hfus
          heatrem = (t_lake(c,1) - tfrz)*cpliq*denh2o*dz_lake(c,1) - heatsum(c)

          if (heatrem + denh2o*dz_lake(c,1)*hfus > 0._r8) then            
             ! Remove snow and subtract the latent heat from the top layer.
             qflx_snomelt(c) = qflx_snomelt(c) + sumsnowice(c)/dtime
             eflx_snomelt(c) = eflx_snomelt(c) + sumsnowice(c)*hfus/dtime 

             ! Update melt per layer. Note that sumsnowice = sum(h2osoi_ice), where the
             ! sum is taken over layers snl(c)+1 to 0. Thus, this code partitions the
             ! above addition to qflx_snomelt (which is based on sumsnowice).
             do j = snl(c)+1,0
                qflx_snomelt_lyr(c,j) = qflx_snomelt_lyr(c,j) + h2osoi_ice(c,j) / dtime
             end do

             ! update incidental drainage from snow pack for this case
             qflx_snow_drain(c) = qflx_snow_drain(c) + h2osno_total(c)/dtime

             h2osno_no_layers(c) = 0._r8
             snow_depth(c) = 0._r8
             snl(c) = 0
             ! The rest of the bookkeeping for the removed snow will be done below.
             if (heatrem > 0._r8) then ! simply subtract the heat from the layer
                t_lake(c,1) = t_lake(c,1) - heatrem/(cpliq*denh2o*dz_lake(c,1))
             else !freeze part of the layer
                t_lake(c,1) = tfrz
                lake_icefrac(c,1) = -heatrem/(denh2o*dz_lake(c,1)*hfus)
             end if
          end if
       end if
    end do

    ! Set empty snow layers to zero
    call ZeroEmptySnowLayers(bounds, num_shlakesnowc, filter_shlakesnowc, &
         col, water_inst, temperature_inst)

    ! Build new snow filter

    call BuildSnowFilter(bounds, num_lakec, filter_lakec, &
         num_shlakesnowc, filter_shlakesnowc, num_shlakenosnowc, filter_shlakenosnowc)

    ! Vertically average t_soisno and sum of h2osoi_liq and h2osoi_ice
    ! over all snow layers for history output

    do fc = 1, num_lakec
       c = filter_lakec(fc)
       snowice(c) = 0._r8
       snowliq(c) = 0._r8
    end do

    do j = -nlevsno+1, 0
       do fc = 1, num_shlakesnowc
          c = filter_shlakesnowc(fc)
          if (j >= snl(c)+1) then
             snowice(c) = snowice(c) + h2osoi_ice(c,j)
             snowliq(c) = snowliq(c) + h2osoi_liq(c,j)
          end if
       end do
    end do

    ! Snow internal temperature
    ! See description in HydrologyNoDrainageMod

    do fc = 1, num_lakec
       c = filter_lakec(fc)
       t_sno_mul_mss(c) = 0._r8
    end do

    do j = -nlevsno+1, 0
       do fc = 1, num_shlakesnowc
          c = filter_shlakesnowc(fc)
          if (j >= snl(c)+1) then
             t_sno_mul_mss(c) = t_sno_mul_mss(c) + h2osoi_ice(c,j) * t_soisno(c,j)
             t_sno_mul_mss(c) = t_sno_mul_mss(c) + h2osoi_liq(c,j) * tfrz
          end if
       end do
    end do

    ! Determine ending water balance and volumetric soil water

    call ComputeWaterMassLake(bounds, num_lakec, filter_lakec, &
         b_waterstate_inst, lakestate_inst, &
         add_lake_water_and_subtract_dynbal_baselines = .false., &
         water_mass = endwb(bounds%begc:bounds%endc))

    do j = 1, nlevgrnd
       do fc = 1, num_lakec
          c = filter_lakec(fc)
          h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)
       end do
    end do

    do fp = 1,num_lakep
       p = filter_lakep(fp)
       c = pcolumn(p)
       g = pgridcell(p)

       qflx_drain_perched(c) = 0._r8
       qflx_rsub_sat(c)      = 0._r8
       qflx_infl(c)          = 0._r8
       qflx_surf(c)          = 0._r8
       qflx_drain(c)         = 0._r8

       ! Insure water balance using qflx_qrgwl
       ! qflx_snwcp_ice(c) has been computed in routine SnowCapping
       qflx_qrgwl(c)     = forc_rain(c) + forc_snow(c) - qflx_evap_tot(p) - qflx_snwcp_ice(c) - &
            qflx_snwcp_discarded_ice(c) - qflx_snwcp_discarded_liq(c) - &
            (endwb(c)-begwb(c))/dtime + qflx_floodg(g)
       qflx_floodc(c)    = qflx_floodg(g)
       qflx_runoff(c)    = qflx_drain(c) + qflx_qrgwl(c)
       qflx_rain_plus_snomelt(c) = qflx_liq_grnd(c) + qflx_snow_drain(c)
       qflx_top_soil(c)  = qflx_rain_plus_snomelt(c)
       qflx_ice_runoff_snwcp(c) = qflx_snwcp_ice(c)

    enddo

    ! top-layer diagnostics
    do fc = 1, num_shlakesnowc
       c = filter_shlakesnowc(fc)
       h2osno_top(c)  = h2osoi_ice(c,snl(c)+1) + h2osoi_liq(c,snl(c)+1)
    end do

    ! Zero variables in columns without snow
    do fc = 1, num_shlakenosnowc
       c = filter_shlakenosnowc(fc)

       h2osno_top(c)      = 0._r8
       snw_rds(c,:)       = 0._r8

       ! top-layer diagnostics (spval is not averaged when computing history fields)
       snot_top(c)        = spval
       dTdz_top(c)        = spval
       snw_rds_top(c)     = spval
       sno_liq_top(c)     = spval
    end do

    end associate

    end associate

  end subroutine LakeHydrology

  !-----------------------------------------------------------------------
  subroutine SumFlux_FluxesOntoGround(bounds, &
       num_lakec, filter_lakec, &
       forc_snow, forc_rain, &
       qflx_snow_grnd, qflx_liq_grnd)
    !
    ! !DESCRIPTION:
    ! Compute "summed" (really just copies here) fluxes onto "ground" (really the lake
    ! surface here), for bulk or one tracer.
    !
    ! (Subroutine name mimics the one in CanopyHydrologyMod.)
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_lakec
    integer, intent(in) :: filter_lakec(:)

    real(r8) , intent(in)    :: forc_rain( bounds%begc: )         ! atm rain rate [mm/s]
    real(r8) , intent(in)    :: forc_snow( bounds%begc: )         ! atm snow rate [mm/s]

    real(r8) , intent(inout) :: qflx_snow_grnd( bounds%begc: )    ! snow on ground after interception (mm H2O/s)
    real(r8) , intent(inout) :: qflx_liq_grnd( bounds%begc: )     ! liquid on ground after interception (mm H2O/s)
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c

    character(len=*), parameter :: subname = 'SumFlux_FluxesOntoGround'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(forc_rain, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(forc_snow, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_snow_grnd, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_liq_grnd, 1) == bounds%endc), sourcefile, __LINE__)

    do fc = 1, num_lakec
       c = filter_lakec(fc)

       qflx_snow_grnd(c) = forc_snow(c)
       qflx_liq_grnd(c)  = forc_rain(c)
    end do

  end subroutine SumFlux_FluxesOntoGround

end module LakeHydrologyMod
