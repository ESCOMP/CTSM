module BalanceCheckMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Water and energy balance check.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use decompMod          , only : bounds_type
  use abortutils         , only : endrun
  use clm_varctl         , only : iulog
  use clm_varcon         , only : namep, namec
  use clm_varpar         , only : nlevsoi
  use GetGlobalValuesMod , only : GetGlobalIndex
  use atm2lndType        , only : atm2lnd_type
  use EnergyFluxType     , only : energyflux_type
  use SolarAbsorbedType  , only : solarabs_type
  use SoilHydrologyType  , only : soilhydrology_type  
  use WaterStateBulkType     , only : waterstatebulk_type
  use WaterDiagnosticBulkType     , only : waterdiagnosticbulk_type
  use WaterBalanceType     , only : waterbalance_type
  use WaterFluxBulkType      , only : waterfluxbulk_type
  use TotalWaterAndHeatMod, only : ComputeWaterMassNonLake, ComputeWaterMassLake
  use GridcellType       , only : grc                
  use LandunitType       , only : lun                
  use ColumnType         , only : col                
  use PatchType          , only : patch                
  use landunit_varcon    , only : istdlak, istsoil,istcrop,istwet,istice_mec
  use column_varcon      , only : icol_roof, icol_sunwall, icol_shadewall
  use column_varcon      , only : icol_road_perv, icol_road_imperv
  use clm_varcon         , only : aquifer_water_baseline
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:

  public :: BeginWaterBalance        ! Initialize water balance check
  public :: BalanceCheck             ! Water and energy balance check

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine BeginWaterBalance(bounds, &
       num_nolakec, filter_nolakec, num_lakec, filter_lakec, &
       soilhydrology_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, waterbalancebulk_inst)
    !
    ! !DESCRIPTION:
    ! Initialize column-level water balance at beginning of time step
    !
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds     
    integer                   , intent(in)    :: num_nolakec          ! number of column non-lake points in column filter
    integer                   , intent(in)    :: filter_nolakec(:)    ! column filter for non-lake points
    integer                   , intent(in)    :: num_lakec            ! number of column lake points in column filter
    integer                   , intent(in)    :: filter_lakec(:)      ! column filter for lake points
    type(soilhydrology_type)  , intent(inout) :: soilhydrology_inst
    type(waterstatebulk_type)     , intent(inout) :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)     , intent(inout) :: waterdiagnosticbulk_inst
    type(waterbalance_type)     , intent(inout) :: waterbalancebulk_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c, j, fc                  ! indices
    !-----------------------------------------------------------------------

    associate(                                               & 
         zi           =>    col%zi                         , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m) 
         zwt          =>    soilhydrology_inst%zwt_col     , & ! Input:  [real(r8) (:)   ]  water table depth (m)                   
         wa           =>    soilhydrology_inst%wa_col      , & ! Output: [real(r8) (:)   ]  water in the unconfined aquifer (mm)    
         begwb        =>    waterbalancebulk_inst%begwb_col        & ! Output: [real(r8) (:)   ]  water mass begining of the time step    
         )

   do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       if (col%hydrologically_active(c)) then
          if(zwt(c) <= zi(c,nlevsoi)) then
             wa(c) = aquifer_water_baseline
          end if
       end if
    end do

    call ComputeWaterMassNonLake(bounds, num_nolakec, filter_nolakec, &
         soilhydrology_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, begwb(bounds%begc:bounds%endc))

    call ComputeWaterMassLake(bounds, num_lakec, filter_lakec, &
         waterstatebulk_inst, begwb(bounds%begc:bounds%endc))

    end associate 

  end subroutine BeginWaterBalance

   !-----------------------------------------------------------------------
   subroutine BalanceCheck( bounds, &
        atm2lnd_inst, solarabs_inst, waterfluxbulk_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, waterbalancebulk_inst, &
        energyflux_inst, canopystate_inst)
     !
     ! !DESCRIPTION:
     ! This subroutine accumulates the numerical truncation errors of the water
     ! and energy balance calculation. It is helpful to see the performance of
     ! the process of integration.
     !
     ! The error for energy balance:
     !
     ! error = abs(Net radiation - change of internal energy - Sensible heat
     !             - Latent heat)
     !
     ! The error for water balance:
     !
     ! error = abs(precipitation - change of water storage - evaporation - runoff)
     !
     ! !USES:
     use clm_varcon        , only : spval
     use clm_time_manager  , only : get_step_size, get_nstep
     use clm_time_manager  , only : get_nstep_since_startup_or_lastDA_restart_or_pause
     use clm_instMod       , only : surfalb_inst
     use CanopyStateType   , only : canopystate_type
     use subgridAveMod
     !
     ! !ARGUMENTS:
     type(bounds_type)     , intent(in)    :: bounds  
     type(atm2lnd_type)    , intent(in)    :: atm2lnd_inst
     type(solarabs_type)   , intent(in)    :: solarabs_inst
     type(waterfluxbulk_type)  , intent(inout) :: waterfluxbulk_inst
     type(waterstatebulk_type) , intent(inout) :: waterstatebulk_inst
     type(waterdiagnosticbulk_type) , intent(inout) :: waterdiagnosticbulk_inst
     type(waterbalance_type) , intent(inout) :: waterbalancebulk_inst
     type(energyflux_type) , intent(inout) :: energyflux_inst
     type(canopystate_type), intent(inout) :: canopystate_inst
     !
     ! !LOCAL VARIABLES:
     integer  :: p,c,l,g,fc                             ! indices
     real(r8) :: dtime                                  ! land model time step (sec)
     integer  :: nstep                                  ! time step number
     integer  :: DAnstep                                ! time step number since last Data Assimilation (DA)
     logical  :: found                                  ! flag in search loop
     integer  :: indexp,indexc,indexl,indexg            ! index of first found in search loop
     real(r8) :: forc_rain_col(bounds%begc:bounds%endc) ! column level rain rate [mm/s]
     real(r8) :: forc_snow_col(bounds%begc:bounds%endc) ! column level snow rate [mm/s]
     !-----------------------------------------------------------------------

     associate(                                                                   & 
          volr                    =>    atm2lnd_inst%volr_grc                   , & ! Input:  [real(r8) (:)   ]  river water storage (m3)                 
          forc_solad              =>    atm2lnd_inst%forc_solad_grc             , & ! Input:  [real(r8) (:,:) ]  direct beam radiation (vis=forc_sols , nir=forc_soll )
          forc_solai              =>    atm2lnd_inst%forc_solai_grc             , & ! Input:  [real(r8) (:,:) ]  diffuse radiation     (vis=forc_solsd, nir=forc_solld)
          forc_rain               =>    atm2lnd_inst%forc_rain_downscaled_col   , & ! Input:  [real(r8) (:)   ]  rain rate [mm/s]
          forc_snow               =>    atm2lnd_inst%forc_snow_downscaled_col   , & ! Input:  [real(r8) (:)   ]  snow rate [mm/s]
          forc_lwrad              =>    atm2lnd_inst%forc_lwrad_downscaled_col  , & ! Input:  [real(r8) (:)   ]  downward infrared (longwave) radiation (W/m**2)

          h2osno                  =>    waterstatebulk_inst%h2osno_col              , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)                     
          h2osno_old              =>    waterbalancebulk_inst%h2osno_old_col          , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O) at previous time step
          frac_sno_eff            =>    waterdiagnosticbulk_inst%frac_sno_eff_col        , & ! Input:  [real(r8) (:)   ]  effective snow fraction                 
          frac_sno                =>    waterdiagnosticbulk_inst%frac_sno_col            , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
          snow_depth              =>    waterdiagnosticbulk_inst%snow_depth_col          , & ! Input:  [real(r8) (:)   ]  snow height (m)                         
          begwb                   =>    waterbalancebulk_inst%begwb_col               , & ! Input:  [real(r8) (:)   ]  water mass begining of the time step    
          errh2o                  =>    waterbalancebulk_inst%errh2o_col              , & ! Output: [real(r8) (:)   ]  water conservation error (mm H2O)       
          errh2osno               =>    waterbalancebulk_inst%errh2osno_col           , & ! Output: [real(r8) (:)   ]  error in h2osno (kg m-2)                
          endwb                   =>    waterbalancebulk_inst%endwb_col               , & ! Output: [real(r8) (:)   ]  water mass end of the time step         
          total_plant_stored_h2o_col => waterdiagnosticbulk_inst%total_plant_stored_h2o_col, & ! Input: [real(r8) (:)   ]  water mass in plant tissues (kg m-2)
          qflx_rootsoi_col        =>    waterfluxbulk_inst%qflx_rootsoi_col         , & ! Input   [real(r8) (:)   ]  water loss in soil layers to root uptake (mm H2O/s)
                                                                                    !                            (ie transpiration demand, often = transpiration)
          qflx_rain_grnd_col      =>    waterfluxbulk_inst%qflx_rain_grnd_col       , & ! Input:  [real(r8) (:)   ]  rain on ground after interception (mm H2O/s) [+]
          qflx_snow_grnd_col      =>    waterfluxbulk_inst%qflx_snow_grnd_col       , & ! Input:  [real(r8) (:)   ]  snow on ground after interception (mm H2O/s) [+]
          qflx_evap_soi           =>    waterfluxbulk_inst%qflx_evap_soi_col        , & ! Input:  [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)
          qflx_snwcp_liq          =>    waterfluxbulk_inst%qflx_snwcp_liq_col       , & ! Input:  [real(r8) (:)   ]  excess liquid h2o due to snow capping (outgoing) (mm H2O /s) [+]`
          qflx_snwcp_ice          =>    waterfluxbulk_inst%qflx_snwcp_ice_col       , & ! Input:  [real(r8) (:)   ]  excess solid h2o due to snow capping (outgoing) (mm H2O /s) [+]`
          qflx_snwcp_discarded_liq =>   waterfluxbulk_inst%qflx_snwcp_discarded_liq_col, & ! Input:  [real(r8) (:)   ]  excess liquid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s) [+]`
          qflx_snwcp_discarded_ice =>   waterfluxbulk_inst%qflx_snwcp_discarded_ice_col, & ! Input:  [real(r8) (:)   ]  excess solid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s) [+]`
          qflx_evap_tot           =>    waterfluxbulk_inst%qflx_evap_tot_col        , & ! Input:  [real(r8) (:)   ]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg
          qflx_dew_snow           =>    waterfluxbulk_inst%qflx_dew_snow_col        , & ! Input:  [real(r8) (:)   ]  surface dew added to snow pack (mm H2O /s) [+]
          qflx_sub_snow           =>    waterfluxbulk_inst%qflx_sub_snow_col        , & ! Input:  [real(r8) (:)   ]  sublimation rate from snow pack (mm H2O /s) [+]
          qflx_evap_grnd          =>    waterfluxbulk_inst%qflx_evap_grnd_col       , & ! Input:  [real(r8) (:)   ]  ground surface evaporation rate (mm H2O/s) [+]
          qflx_dew_grnd           =>    waterfluxbulk_inst%qflx_dew_grnd_col        , & ! Input:  [real(r8) (:)   ]  ground surface dew formation (mm H2O /s) [+]
          qflx_prec_grnd          =>    waterfluxbulk_inst%qflx_prec_grnd_col       , & ! Input:  [real(r8) (:)   ]  water onto ground including canopy runoff [kg/(m2 s)]
          qflx_snow_h2osfc        =>    waterfluxbulk_inst%qflx_snow_h2osfc_col     , & ! Input:  [real(r8) (:)   ]  snow falling on surface water (mm/s)    
          qflx_h2osfc_to_ice      =>    waterfluxbulk_inst%qflx_h2osfc_to_ice_col   , & ! Input:  [real(r8) (:)   ]  conversion of h2osfc to ice             
          qflx_drain_perched      =>    waterfluxbulk_inst%qflx_drain_perched_col   , & ! Input:  [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)          
          qflx_floodc             =>    waterfluxbulk_inst%qflx_floodc_col          , & ! Input:  [real(r8) (:)   ]  total runoff due to flooding            
          qflx_snow_drain         =>    waterfluxbulk_inst%qflx_snow_drain_col      , & ! Input:  [real(r8) (:)   ]  drainage from snow pack                         
          qflx_surf               =>    waterfluxbulk_inst%qflx_surf_col            , & ! Input:  [real(r8) (:)   ]  surface runoff (mm H2O /s)              
          qflx_qrgwl              =>    waterfluxbulk_inst%qflx_qrgwl_col           , & ! Input:  [real(r8) (:)   ]  qflx_surf at glaciers, wetlands, lakes  
          qflx_drain              =>    waterfluxbulk_inst%qflx_drain_col           , & ! Input:  [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)          
          qflx_runoff             =>    waterfluxbulk_inst%qflx_runoff_col          , & ! Input:  [real(r8) (:)   ]  total runoff (mm H2O /s)                
          qflx_ice_runoff_snwcp   =>    waterfluxbulk_inst%qflx_ice_runoff_snwcp_col, & ! Input:  [real(r8) (:)   ] solid runoff from snow capping (mm H2O /s)
          qflx_ice_runoff_xs      =>    waterfluxbulk_inst%qflx_ice_runoff_xs_col   , & ! Input:  [real(r8) (:)   ] solid runoff from excess ice in soil (mm H2O /s)
          qflx_top_soil           =>    waterfluxbulk_inst%qflx_top_soil_col        , & ! Input:  [real(r8) (:)   ]  net water input into soil from top (mm/s)
          qflx_sl_top_soil        =>    waterfluxbulk_inst%qflx_sl_top_soil_col     , & ! Input:  [real(r8) (:)   ]  liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)
          qflx_liq_dynbal         =>    waterfluxbulk_inst%qflx_liq_dynbal_grc      , & ! Input:  [real(r8) (:)   ]  liq runoff due to dynamic land cover change (mm H2O /s)
          qflx_ice_dynbal         =>    waterfluxbulk_inst%qflx_ice_dynbal_grc      , & ! Input:  [real(r8) (:)   ]  ice runoff due to dynamic land cover change (mm H2O /s)
          snow_sources            =>    waterfluxbulk_inst%snow_sources_col         , & ! Output: [real(r8) (:)   ]  snow sources (mm H2O /s)  
          snow_sinks              =>    waterfluxbulk_inst%snow_sinks_col           , & ! Output: [real(r8) (:)   ]  snow sinks (mm H2O /s)    

          qflx_irrig              =>    waterfluxbulk_inst%qflx_irrig_col          , & ! Input:  [real(r8) (:)   ]  irrigation flux (mm H2O /s)             

          qflx_glcice_dyn_water_flux => waterfluxbulk_inst%qflx_glcice_dyn_water_flux_col, & ! Input: [real(r8) (:)]  water flux needed for balance check due to glc_dyn_runoff_routing (mm H2O/s) (positive means addition of water to the system)

          eflx_lwrad_out          =>    energyflux_inst%eflx_lwrad_out_patch    , & ! Input:  [real(r8) (:)   ]  emitted infrared (longwave) radiation (W/m**2)
          eflx_lwrad_net          =>    energyflux_inst%eflx_lwrad_net_patch    , & ! Input:  [real(r8) (:)   ]  net infrared (longwave) rad (W/m**2) [+ = to atm]
          eflx_sh_tot             =>    energyflux_inst%eflx_sh_tot_patch       , & ! Input:  [real(r8) (:)   ]  total sensible heat flux (W/m**2) [+ to atm]
          eflx_lh_tot             =>    energyflux_inst%eflx_lh_tot_patch       , & ! Input:  [real(r8) (:)   ]  total latent heat flux (W/m**2)  [+ to atm]
          eflx_soil_grnd          =>    energyflux_inst%eflx_soil_grnd_patch    , & ! Input:  [real(r8) (:)   ]  soil heat flux (W/m**2) [+ = into soil] 
          eflx_wasteheat_patch    =>    energyflux_inst%eflx_wasteheat_patch    , & ! Input:  [real(r8) (:)   ]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
          eflx_heat_from_ac_patch =>    energyflux_inst%eflx_heat_from_ac_patch , & ! Input:  [real(r8) (:)   ]  sensible heat flux put back into canyon due to removal by AC (W/m**2)
          eflx_traffic_patch      =>    energyflux_inst%eflx_traffic_patch      , & ! Input:  [real(r8) (:)   ]  traffic sensible heat flux (W/m**2)     
          eflx_dynbal             =>    energyflux_inst%eflx_dynbal_grc         , & ! Input:  [real(r8) (:)   ]  energy conversion flux due to dynamic land cover change(W/m**2) [+ to atm]
          errsoi_col              =>    energyflux_inst%errsoi_col              , & ! Output: [real(r8) (:)   ]  column-level soil/lake energy conservation error (W/m**2)
          errsol                  =>    energyflux_inst%errsol_patch            , & ! Output: [real(r8) (:)   ]  solar radiation conservation error (W/m**2)
          errseb                  =>    energyflux_inst%errseb_patch            , & ! Output: [real(r8) (:)   ]  surface energy conservation error (W/m**2)
          errlon                  =>    energyflux_inst%errlon_patch            , & ! Output: [real(r8) (:)   ]  longwave radiation conservation error (W/m**2)

          sabg_soil               =>    solarabs_inst%sabg_soil_patch           , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by soil (W/m**2)
          sabg_snow               =>    solarabs_inst%sabg_snow_patch           , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by snow (W/m**2)
          sabg_chk                =>    solarabs_inst%sabg_chk_patch            , & ! Input:  [real(r8) (:)   ]  sum of soil/snow using current fsno, for balance check
          fsa                     =>    solarabs_inst%fsa_patch                 , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed (total) (W/m**2)
          fsr                     =>    solarabs_inst%fsr_patch                 , & ! Input:  [real(r8) (:)   ]  solar radiation reflected (W/m**2)      
          sabv                    =>    solarabs_inst%sabv_patch                , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by vegetation (W/m**2)
          sabg                    =>    solarabs_inst%sabg_patch                , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by ground (W/m**2)
          
          elai                    =>    canopystate_inst%elai_patch             , & ! Input:  [real(r8) (:,:)]  
          esai                    =>    canopystate_inst%esai_patch             , & ! Input:  [real(r8) (:,:)]  

          fabd                    =>    surfalb_inst%fabd_patch                 , & ! Input:  [real(r8) (:,:)]  flux absorbed by canopy per unit direct flux
          fabi                    =>    surfalb_inst%fabi_patch                 , & ! Input:  [real(r8) (:,:)]  flux absorbed by canopy per unit indirect flux
          albd                    =>    surfalb_inst%albd_patch                 , & ! Output: [real(r8) (:,:)]  surface albedo (direct)
          albi                    =>    surfalb_inst%albi_patch                 , & ! Output: [real(r8) (:,:)]  surface albedo (diffuse)
          ftdd                    =>    surfalb_inst%ftdd_patch                 , & ! Input:  [real(r8) (:,:)]  down direct flux below canopy per unit direct flux
          ftid                    =>    surfalb_inst%ftid_patch                 , & ! Input:  [real(r8) (:,:)]  down diffuse flux below canopy per unit direct flux
          ftii                    =>    surfalb_inst%ftii_patch                 , & ! Input:  [real(r8) (:,:)]  down diffuse flux below canopy per unit diffuse flux

          netrad                  =>    energyflux_inst%netrad_patch              & ! Output: [real(r8) (:)   ]  net radiation (positive downward) (W/m**2)
          )

       ! Get step size and time step

       nstep = get_nstep()
       DAnstep = get_nstep_since_startup_or_lastDA_restart_or_pause()
       dtime = get_step_size()

       ! Determine column level incoming snow and rain
       ! Assume no incident precipitation on urban wall columns (as in CanopyHydrologyMod.F90).

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          l = col%landunit(c)       

          if (col%itype(c) == icol_sunwall .or.  col%itype(c) == icol_shadewall) then
             forc_rain_col(c) = 0.
             forc_snow_col(c) = 0.
          else
             forc_rain_col(c) = forc_rain(c)
             forc_snow_col(c) = forc_snow(c)
          end if
       end do

       ! Water balance check

       do c = bounds%begc, bounds%endc

          ! add qflx_drain_perched and qflx_flood
          if (col%active(c)) then

             errh2o(c) = endwb(c) - begwb(c) &
                  - (forc_rain_col(c)        &
                  + forc_snow_col(c)         &
                  + qflx_floodc(c)           &
                  + qflx_irrig(c)            &
                  + qflx_glcice_dyn_water_flux(c) &
                  - qflx_evap_tot(c)         &
                  - qflx_surf(c)             &
                  - qflx_qrgwl(c)            &
                  - qflx_drain(c)            &
                  - qflx_drain_perched(c)    &
                  - qflx_ice_runoff_snwcp(c) &
                  - qflx_ice_runoff_xs(c)    &
                  - qflx_snwcp_discarded_liq(c) &
                  - qflx_snwcp_discarded_ice(c)) * dtime

          else

             errh2o(c) = 0.0_r8

          end if

       end do

       found = .false.
       do c = bounds%begc, bounds%endc
          if (abs(errh2o(c)) > 1.e-9_r8) then
             found = .true.
             indexc = c
          end if
       end do

       if ( found ) then

          write(iulog,*)'WARNING:  water balance error ',&
               ' nstep= ',nstep, &
               ' local indexc= ',indexc,&
               ! ' global indexc= ',GetGlobalIndex(decomp_index=indexc, clmlevel=namec), &
               ' errh2o= ',errh2o(indexc)

          if ((col%itype(indexc) == icol_roof .or. &
               col%itype(indexc) == icol_road_imperv .or. &
               col%itype(indexc) == icol_road_perv) .and. &
               abs(errh2o(indexc)) > 1.e-5_r8 .and. (DAnstep > 2) ) then

             write(iulog,*)'clm urban model is stopping - error is greater than 1e-5 (mm)'
             write(iulog,*)'nstep                 = ',nstep
             write(iulog,*)'errh2o                = ',errh2o(indexc)
             write(iulog,*)'forc_rain             = ',forc_rain_col(indexc)*dtime
             write(iulog,*)'forc_snow             = ',forc_snow_col(indexc)*dtime
             write(iulog,*)'endwb                 = ',endwb(indexc)
             write(iulog,*)'begwb                 = ',begwb(indexc)
             write(iulog,*)'qflx_evap_tot         = ',qflx_evap_tot(indexc)*dtime
             write(iulog,*)'qflx_irrig            = ',qflx_irrig(indexc)*dtime
             write(iulog,*)'qflx_surf             = ',qflx_surf(indexc)*dtime
             write(iulog,*)'qflx_qrgwl            = ',qflx_qrgwl(indexc)*dtime
             write(iulog,*)'qflx_drain            = ',qflx_drain(indexc)*dtime
             write(iulog,*)'qflx_ice_runoff_snwcp = ',qflx_ice_runoff_snwcp(indexc)*dtime
             write(iulog,*)'qflx_ice_runoff_xs    = ',qflx_ice_runoff_xs(indexc)*dtime
             write(iulog,*)'qflx_snwcp_discarded_ice = ',qflx_snwcp_discarded_ice(indexc)*dtime
             write(iulog,*)'qflx_snwcp_discarded_liq = ',qflx_snwcp_discarded_liq(indexc)*dtime
             write(iulog,*)'qflx_rootsoi_col(1:nlevsoil)  = ',qflx_rootsoi_col(indexc,:)*dtime
             write(iulog,*)'total_plant_stored_h2o_col = ',total_plant_stored_h2o_col(indexc)
             write(iulog,*)'deltawb          = ',endwb(indexc)-begwb(indexc)
             write(iulog,*)'deltawb/dtime    = ',(endwb(indexc)-begwb(indexc))/dtime
             write(iulog,*)'deltaflux        = ',forc_rain_col(indexc)+forc_snow_col(indexc) - (qflx_evap_tot(indexc) + &
                  qflx_surf(indexc)+qflx_drain(indexc))

             write(iulog,*)'clm model is stopping'
             call endrun(decomp_index=indexc, clmlevel=namec, msg=errmsg(sourcefile, __LINE__))

          else if (abs(errh2o(indexc)) > 1.e-5_r8 .and. (DAnstep > 2) ) then

             write(iulog,*)'clm model is stopping - error is greater than 1e-5 (mm)'
             write(iulog,*)'nstep                 = ',nstep
             write(iulog,*)'errh2o                = ',errh2o(indexc)
             write(iulog,*)'forc_rain             = ',forc_rain_col(indexc)*dtime
             write(iulog,*)'forc_snow             = ',forc_snow_col(indexc)*dtime
             write(iulog,*)'total_plant_stored_h2o_col = ',total_plant_stored_h2o_col(indexc)
             write(iulog,*)'endwb                 = ',endwb(indexc)
             write(iulog,*)'begwb                 = ',begwb(indexc)
             
             write(iulog,*)'qflx_evap_tot         = ',qflx_evap_tot(indexc)*dtime
             write(iulog,*)'qflx_irrig            = ',qflx_irrig(indexc)*dtime
             write(iulog,*)'qflx_surf             = ',qflx_surf(indexc)*dtime
             write(iulog,*)'qflx_qrgwl            = ',qflx_qrgwl(indexc)*dtime
             write(iulog,*)'qflx_drain            = ',qflx_drain(indexc)*dtime
             write(iulog,*)'qflx_drain_perched    = ',qflx_drain_perched(indexc)*dtime
             write(iulog,*)'qflx_flood            = ',qflx_floodc(indexc)*dtime
             write(iulog,*)'qflx_ice_runoff_snwcp = ',qflx_ice_runoff_snwcp(indexc)*dtime
             write(iulog,*)'qflx_ice_runoff_xs    = ',qflx_ice_runoff_xs(indexc)*dtime
             write(iulog,*)'qflx_glcice_dyn_water_flux = ', qflx_glcice_dyn_water_flux(indexc)*dtime
             write(iulog,*)'qflx_snwcp_discarded_ice = ',qflx_snwcp_discarded_ice(indexc)*dtime
             write(iulog,*)'qflx_snwcp_discarded_liq = ',qflx_snwcp_discarded_liq(indexc)*dtime
             write(iulog,*)'qflx_rootsoi_col(1:nlevsoil)  = ',qflx_rootsoi_col(indexc,:)*dtime
             write(iulog,*)'clm model is stopping'
             call endrun(decomp_index=indexc, clmlevel=namec, msg=errmsg(sourcefile, __LINE__))
          end if
       end if

       ! Snow balance check

       do c = bounds%begc,bounds%endc
          if (col%active(c)) then
             g = col%gridcell(c)
             l = col%landunit(c)

             ! As defined here, snow_sources - snow_sinks will equal the change in h2osno at 
             ! any given time step but only if there is at least one snow layer.  h2osno 
             ! also includes snow that is part of the soil column (an initial snow layer is 
             ! only created if h2osno > 10mm).

             if (col%snl(c) < 0) then
                snow_sources(c) = qflx_prec_grnd(c) + qflx_dew_snow(c) + qflx_dew_grnd(c)
                snow_sinks(c)  = qflx_sub_snow(c) + qflx_evap_grnd(c) + qflx_snow_drain(c) &
                     + qflx_snwcp_ice(c) + qflx_snwcp_liq(c) &
                     + qflx_snwcp_discarded_ice(c) + qflx_snwcp_discarded_liq(c) &
                     + qflx_sl_top_soil(c)

                if (lun%itype(l) == istdlak) then 
                   snow_sources(c) = qflx_snow_grnd_col(c) &
                        + frac_sno_eff(c) * (qflx_rain_grnd_col(c) &
                        +  qflx_dew_snow(c) + qflx_dew_grnd(c) ) 
                   snow_sinks(c)   = frac_sno_eff(c) * (qflx_sub_snow(c) + qflx_evap_grnd(c) ) &
                        + qflx_snwcp_ice(c) + qflx_snwcp_liq(c)  &
                        + qflx_snwcp_discarded_ice(c) + qflx_snwcp_discarded_liq(c)  &
                        + qflx_snow_drain(c)  + qflx_sl_top_soil(c)
                endif

                 if (col%itype(c) == icol_road_perv .or. lun%itype(l) == istsoil .or. &
                      lun%itype(l) == istcrop .or. lun%itype(l) == istwet .or. &
                      lun%itype(l) == istice_mec) then
                   snow_sources(c) = (qflx_snow_grnd_col(c) - qflx_snow_h2osfc(c) ) &
                          + frac_sno_eff(c) * (qflx_rain_grnd_col(c) &
                          +  qflx_dew_snow(c) + qflx_dew_grnd(c) ) + qflx_h2osfc_to_ice(c)
                   snow_sinks(c) = frac_sno_eff(c) * (qflx_sub_snow(c) + qflx_evap_grnd(c)) &
                          + qflx_snwcp_ice(c) + qflx_snwcp_liq(c) &
                          + qflx_snwcp_discarded_ice(c) + qflx_snwcp_discarded_liq(c) &
                          + qflx_snow_drain(c) + qflx_sl_top_soil(c)
                endif

                errh2osno(c) = (h2osno(c) - h2osno_old(c)) - (snow_sources(c) - snow_sinks(c)) * dtime
             else
                snow_sources(c) = 0._r8
                snow_sinks(c) = 0._r8
                errh2osno(c) = 0._r8
             end if

          end if
       end do

       found = .false.
       do c = bounds%begc,bounds%endc
          if (col%active(c)) then
             if (abs(errh2osno(c)) > 1.0e-9_r8) then
                found = .true.
                indexc = c
             end if
          end if
       end do
       if ( found ) then
          write(iulog,*)'WARNING:  snow balance error '
          write(iulog,*)'nstep= ',nstep, &
               ' local indexc= ',indexc, &
               ! ' global indexc= ',GetGlobalIndex(decomp_index=indexc, clmlevel=namec), &
               ' col%itype= ',col%itype(indexc), &
               ' lun%itype= ',lun%itype(col%landunit(indexc)), &
               ' errh2osno= ',errh2osno(indexc)

          if (abs(errh2osno(indexc)) > 1.e-5_r8 .and. (DAnstep > 2) ) then
             write(iulog,*)'clm model is stopping - error is greater than 1e-5 (mm)'
             write(iulog,*)'nstep              = ',nstep
             write(iulog,*)'errh2osno          = ',errh2osno(indexc)
             write(iulog,*)'snl                = ',col%snl(indexc)
             write(iulog,*)'snow_depth         = ',snow_depth(indexc)
             write(iulog,*)'frac_sno_eff       = ',frac_sno_eff(indexc)
             write(iulog,*)'h2osno             = ',h2osno(indexc)
             write(iulog,*)'h2osno_old         = ',h2osno_old(indexc)
             write(iulog,*)'snow_sources       = ',snow_sources(indexc)*dtime
             write(iulog,*)'snow_sinks         = ',snow_sinks(indexc)*dtime
             write(iulog,*)'qflx_prec_grnd     = ',qflx_prec_grnd(indexc)*dtime
             write(iulog,*)'qflx_snow_grnd_col = ',qflx_snow_grnd_col(indexc)*dtime
             write(iulog,*)'qflx_rain_grnd_col = ',qflx_rain_grnd_col(indexc)*dtime
             write(iulog,*)'qflx_sub_snow      = ',qflx_sub_snow(indexc)*dtime
             write(iulog,*)'qflx_snow_drain    = ',qflx_snow_drain(indexc)*dtime
             write(iulog,*)'qflx_evap_grnd     = ',qflx_evap_grnd(indexc)*dtime
             write(iulog,*)'qflx_top_soil      = ',qflx_top_soil(indexc)*dtime
             write(iulog,*)'qflx_dew_snow      = ',qflx_dew_snow(indexc)*dtime
             write(iulog,*)'qflx_dew_grnd      = ',qflx_dew_grnd(indexc)*dtime
             write(iulog,*)'qflx_snwcp_ice     = ',qflx_snwcp_ice(indexc)*dtime
             write(iulog,*)'qflx_snwcp_liq     = ',qflx_snwcp_liq(indexc)*dtime
             write(iulog,*)'qflx_snwcp_discarded_ice = ',qflx_snwcp_discarded_ice(indexc)*dtime
             write(iulog,*)'qflx_snwcp_discarded_liq = ',qflx_snwcp_discarded_liq(indexc)*dtime
             write(iulog,*)'qflx_sl_top_soil   = ',qflx_sl_top_soil(indexc)*dtime
             write(iulog,*)'clm model is stopping'

             call endrun(decomp_index=indexc, clmlevel=namec, msg=errmsg(sourcefile, __LINE__))
          end if
       end if

       ! Energy balance checks

       do p = bounds%begp, bounds%endp
          if (patch%active(p)) then
             c = patch%column(p)
             l = patch%landunit(p)
             g = patch%gridcell(p)

             ! Solar radiation energy balance
             ! Do not do this check for an urban patch since it will not balance on a per-column
             ! level because of interactions between columns and since a separate check is done
             ! in the urban radiation module
             if (.not. lun%urbpoi(l)) then
                errsol(p) = fsa(p) + fsr(p) &
                     - (forc_solad(g,1) + forc_solad(g,2) + forc_solai(g,1) + forc_solai(g,2))
             else
                errsol(p) = spval
             end if

             ! Longwave radiation energy balance
             ! Do not do this check for an urban patch since it will not balance on a per-column
             ! level because of interactions between columns and since a separate check is done
             ! in the urban radiation module
             if (.not. lun%urbpoi(l)) then
                errlon(p) = eflx_lwrad_out(p) - eflx_lwrad_net(p) - forc_lwrad(c)
             else
                errlon(p) = spval
             end if

             ! Surface energy balance
             ! Changed to using (eflx_lwrad_net) here instead of (forc_lwrad - eflx_lwrad_out) because
             ! there are longwave interactions between urban columns (and therefore patches). 
             ! For surfaces other than urban, (eflx_lwrad_net) equals (forc_lwrad - eflx_lwrad_out),
             ! and a separate check is done above for these terms.

             if (.not. lun%urbpoi(l)) then
                errseb(p) = sabv(p) + sabg_chk(p) + forc_lwrad(c) - eflx_lwrad_out(p) &
                     - eflx_sh_tot(p) - eflx_lh_tot(p) - eflx_soil_grnd(p)
             else
                errseb(p) = sabv(p) + sabg(p) &
                     - eflx_lwrad_net(p) &
                     - eflx_sh_tot(p) - eflx_lh_tot(p) - eflx_soil_grnd(p) &
                     + eflx_wasteheat_patch(p) + eflx_heat_from_ac_patch(p) + eflx_traffic_patch(p)
             end if
             !TODO MV - move this calculation to a better place - does not belong in BalanceCheck 
             netrad(p) = fsa(p) - eflx_lwrad_net(p) 
          end if
       end do

       ! Solar radiation energy balance check

       found = .false.
       do p = bounds%begp, bounds%endp
          if (patch%active(p)) then
             if ( (errsol(p) /= spval) .and. (abs(errsol(p)) > 1.e-7_r8) ) then
                found = .true.
                indexp = p
                indexg = patch%gridcell(indexp)
             end if
          end if
       end do
       if ( found  .and. (DAnstep > 2) ) then
          write(iulog,*)'WARNING:: BalanceCheck, solar radiation balance error (W/m2)'
          write(iulog,*)'nstep         = ',nstep
          write(iulog,*)'errsol        = ',errsol(indexp)
          if (abs(errsol(indexp)) > 1.e-5_r8 ) then
             write(iulog,*)'clm model is stopping - error is greater than 1e-5 (W/m2)'
             write(iulog,*)'fsa           = ',fsa(indexp)
             write(iulog,*)'fsr           = ',fsr(indexp)
             write(iulog,*)'forc_solad(1) = ',forc_solad(indexg,1)
             write(iulog,*)'forc_solad(2) = ',forc_solad(indexg,2)
             write(iulog,*)'forc_solai(1) = ',forc_solai(indexg,1)
             write(iulog,*)'forc_solai(2) = ',forc_solai(indexg,2)
             write(iulog,*)'forc_tot      = ',forc_solad(indexg,1)+forc_solad(indexg,2) &
               +forc_solai(indexg,1)+forc_solai(indexg,2)
             write(iulog,*)'clm model is stopping'
             call endrun(decomp_index=indexp, clmlevel=namep, msg=errmsg(sourcefile, __LINE__))
          end if
       end if

       ! Longwave radiation energy balance check

       found = .false.
       do p = bounds%begp, bounds%endp
          if (patch%active(p)) then
             if ( (errlon(p) /= spval) .and. (abs(errlon(p)) > 1.e-7_r8) ) then
                found = .true.
                indexp = p
             end if
          end if
       end do
       if ( found  .and. (DAnstep > 2) ) then
          write(iulog,*)'WARNING: BalanceCheck: longwave energy balance error (W/m2)' 
          write(iulog,*)'nstep        = ',nstep 
          write(iulog,*)'errlon       = ',errlon(indexp)
          if (abs(errlon(indexp)) > 1.e-5_r8 ) then
             write(iulog,*)'clm model is stopping - error is greater than 1e-5 (W/m2)'
             call endrun(decomp_index=indexp, clmlevel=namep, msg=errmsg(sourcefile, __LINE__))
          end if
       end if

       ! Surface energy balance check

       found = .false.
       do p = bounds%begp, bounds%endp
          if (patch%active(p)) then
             if (abs(errseb(p)) > 1.e-7_r8 ) then
                found = .true.
                indexp = p
                indexc = patch%column(indexp)
                indexg = patch%gridcell(indexp)
             end if
          end if
       end do
       if ( found  .and. (DAnstep > 2) ) then
          write(iulog,*)'WARNING: BalanceCheck: surface flux energy balance error (W/m2)'
          write(iulog,*)'nstep          = ' ,nstep
          write(iulog,*)'errseb         = ' ,errseb(indexp)
          if (abs(errseb(indexp)) > 1.e-5_r8 ) then
             write(iulog,*)'clm model is stopping - error is greater than 1e-5 (W/m2)'
             write(iulog,*)'sabv           = ' ,sabv(indexp)

             write(iulog,*)'sabg           = ' ,sabg(indexp), ((1._r8- frac_sno(indexc))*sabg_soil(indexp) + &
                  frac_sno(indexc)*sabg_snow(indexp)),sabg_chk(indexp)

             write(iulog,*)'forc_tot      = '  ,forc_solad(indexg,1) + forc_solad(indexg,2) + &
                  forc_solai(indexg,1) + forc_solai(indexg,2)

             write(iulog,*)'eflx_lwrad_net = ' ,eflx_lwrad_net(indexp)
             write(iulog,*)'eflx_sh_tot    = ' ,eflx_sh_tot(indexp)
             write(iulog,*)'eflx_lh_tot    = ' ,eflx_lh_tot(indexp)
             write(iulog,*)'eflx_soil_grnd = ' ,eflx_soil_grnd(indexp)
             write(iulog,*)'fsa fsr = '        ,fsa(indexp),    fsr(indexp)
             write(iulog,*)'fabd fabi = '      ,fabd(indexp,:), fabi(indexp,:)
             write(iulog,*)'albd albi = '      ,albd(indexp,:), albi(indexp,:)
             write(iulog,*)'ftii ftdd ftid = ' ,ftii(indexp,:), ftdd(indexp,:),ftid(indexp,:)
             write(iulog,*)'elai esai = '      ,elai(indexp),   esai(indexp)      
             write(iulog,*)'clm model is stopping'
             call endrun(decomp_index=indexp, clmlevel=namep, msg=errmsg(sourcefile, __LINE__))
          end if
       end if

       ! Soil energy balance check

       found = .false.
       do c = bounds%begc,bounds%endc
          if (col%active(c)) then
             if (abs(errsoi_col(c)) > 1.0e-5_r8 ) then
                found = .true.
                indexc = c
             end if
          end if
       end do
       if ( found ) then
          write(iulog,*)'WARNING: BalanceCheck: soil balance error (W/m2)'
          write(iulog,*)'nstep         = ',nstep
          write(iulog,*)'errsoi_col    = ',errsoi_col(indexc)
          if (abs(errsoi_col(indexc)) > 1.e-4_r8 .and. (DAnstep > 2) ) then
             write(iulog,*)'clm model is stopping'
             call endrun(decomp_index=indexc, clmlevel=namec, msg=errmsg(sourcefile, __LINE__))
          end if
       end if

     end associate

   end subroutine BalanceCheck

end module BalanceCheckMod
