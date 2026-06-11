module SoilFluxesMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Updates surface fluxes based on the new ground temperature.
  !
  ! !USES:
  use shr_kind_mod	, only : r8 => shr_kind_r8
  use shr_log_mod	, only : errMsg => shr_log_errMsg
  use decompMod		, only : bounds_type
  use abortutils	, only : endrun
  use perf_mod		, only : t_startf, t_stopf
  use clm_varctl	, only : iulog, use_nvp  ! [PORTED by Hui Tang: use_nvp for NVP snow-burial guard]
  use clm_varpar	, only : nlevsno, nlevgrnd, nlevurb
  use atm2lndType	, only : atm2lnd_type
  use CanopyStateType   , only : canopystate_type
  use EnergyFluxType    , only : energyflux_type
  use SolarAbsorbedType , only : solarabs_type
  use TemperatureType   , only : temperature_type
  use WaterStateBulkType    , only : waterstatebulk_type
  use WaterDiagnosticBulkType    , only : waterdiagnosticbulk_type
  use WaterFluxBulkType     , only : waterfluxbulk_type
  use LandunitType	, only : lun                
  use ColumnType	, only : col                
  use PatchType		, only : patch                
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilFluxes   ! Calculate soil/snow and ground temperatures
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SoilFluxes (bounds, num_urbanl, filter_urbanl, &
       num_urbanp, filter_urbanp, &
       num_nolakec, filter_nolakec, num_nolakep, filter_nolakep, &
       atm2lnd_inst, solarabs_inst, temperature_inst, canopystate_inst, &
       waterstatebulk_inst, waterdiagnosticbulk_inst, energyflux_inst, waterfluxbulk_inst)            
    !
    ! !DESCRIPTION:
    ! Update surface fluxes based on the new ground temperature
    !
    ! !USES:
    use clm_time_manager , only : get_step_size_real
    use clm_varcon       , only : hvap, cpair, grav, vkc, tfrz, sb 
    use landunit_varcon  , only : istsoil, istcrop
    use column_varcon    , only : icol_roof, icol_sunwall, icol_shadewall, icol_road_perv
    use subgridAveMod    , only : p2c
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds    
    integer                , intent(in)    :: num_nolakec                      ! number of column non-lake points in column filter
    integer                , intent(in)    :: filter_nolakec(:)                ! column filter for non-lake points
    integer                , intent(in)    :: num_urbanl                       ! number of urban landunits in clump
    integer                , intent(in)    :: filter_urbanl(:)                 ! urban landunit filter
    integer                , intent(in)    :: num_urbanp                       ! number of urban pfts in clump
    integer                , intent(in)    :: filter_urbanp(:)                 ! urban pft filter
    integer                , intent(in)    :: num_nolakep                      ! number of column non-lake points in pft filter
    integer                , intent(in)    :: filter_nolakep(:)                ! patch filter for non-lake points
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(solarabs_type)    , intent(in)    :: solarabs_inst
    type(temperature_type) , intent(in)    :: temperature_inst
    type(canopystate_type) , intent(in)    :: canopystate_inst
    type(waterstatebulk_type)  , intent(in)    :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)  , intent(in)    :: waterdiagnosticbulk_inst
    type(waterfluxbulk_type)   , intent(inout) :: waterfluxbulk_inst
    type(energyflux_type)  , intent(inout) :: energyflux_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,g,j,pi,l                                       ! indices
    integer  :: fc,fp                                              ! lake filtered column and pft indices
    real(r8) :: dtime                                              ! land model time step (sec)
    real(r8) :: tinc(bounds%begc:bounds%endc)                      ! temperature difference of two time step
    real(r8) :: eflx_lwrad_del(bounds%begp:bounds%endp)            ! update due to eflx_lwrad
    real(r8) :: t_grnd0(bounds%begc:bounds%endc)                   ! t_grnd of previous time step
    real(r8) :: lw_grnd
    real(r8) :: frac_nvp_eff                                       ! [PORTED by Hui Tang: effective NVP fraction for LW weighting]
    real(r8) :: frac_soil                                          ! [PORTED by Hui Tang: exposed bare-soil fraction (1-fsno-fh2osfc-fnvp)]
    real(r8) :: qflx_evap_grnd_eff                                 ! [PORTED by Hui Tang: per-surface ground evap total for energy/diagnostic consistency]
    real(r8) :: evaporation_limit                                  ! top layer moisture available for evaporation
    real(r8) :: evaporation_demand                                   ! evaporative demand
    real(r8) :: heat_store_diag                                      ! [PORTED by Hui Tang: errsoi diagnostic - heat storage sum]
    real(r8) :: wgt                                                  ! [PORTED by Hui Tang: errsoi per-layer diagnostic - applied frac weight]
    real(r8) :: eflx_soil_grnd_nvp                                   ! [PORTED by Hui Tang: VERIFY-ONLY candidate NVP-consistent errsoi input flux (W/m2)]
    real(r8) :: errsoi_test                                          ! [PORTED by Hui Tang: VERIFY-ONLY candidate errsoi residual using NVP-consistent input (W/m2)]
    !-----------------------------------------------------------------------

    associate(                                                                & 
         eflx_sh_stem            => energyflux_inst%eflx_sh_stem_patch      , & ! Output: [real(r8) (:)   ]  sensible heat flux from stems (W/m**2) [+ to atm]
         eflx_h2osfc_to_snow_col => energyflux_inst%eflx_h2osfc_to_snow_col , & ! Input:  [real(r8) (:)   ]  col snow melt to h2osfc heat flux (W/m**2)

         forc_lwrad              => atm2lnd_inst%forc_lwrad_downscaled_col  , & ! Input:  [real(r8) (:)   ]  downward infrared (longwave) radiation (W/m**2)

         frac_veg_nosno          => canopystate_inst%frac_veg_nosno_patch   , & ! Input:  [integer (:)    ]  fraction of veg not covered by snow (0/1 now) [-]

         frac_sno_eff            => waterdiagnosticbulk_inst%frac_sno_eff_col        , & ! Input:  [real(r8) (:)   ]  eff. fraction of ground covered by snow (0 to 1)
         frac_h2osfc             => waterdiagnosticbulk_inst%frac_h2osfc_col         , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)
         h2osoi_ice              => waterstatebulk_inst%h2osoi_ice_col          , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2) (new)                
         h2osoi_liq              => waterstatebulk_inst%h2osoi_liq_col          , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2) (new)            
         sabg_soil               => solarabs_inst%sabg_soil_patch           , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by soil (W/m**2)
         sabg_snow               => solarabs_inst%sabg_snow_patch           , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by snow (W/m**2)
         sabg                    => solarabs_inst%sabg_patch                , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by ground (W/m**2)
         ! [PORTED by Hui Tang: NVP errsoi fix - solar by layer and NVP sensible heat]
         sabg_lyr                => solarabs_inst%sabg_lyr_patch            , & ! Input:  [real(r8) (:,:) ]  solar radiation absorbed per snow/soil layer (W/m**2)
         eflx_sh_nvp             => energyflux_inst%eflx_sh_nvp_patch      , & ! Input:  [real(r8) (:)   ]  sensible heat flux from NVP (W/m**2) [+ to atm]

         emg                     => temperature_inst%emg_col                , & ! Input:  [real(r8) (:)   ]  ground emissivity                       
!         emv                     => temperature_inst%emv_patch              , & ! Input:  [real(r8) (:)   ]  vegetation emissivity
!         t_veg                   => temperature_inst%t_veg_patch            , & ! Output: [real(r8) (:)   ]  vegetation temperature (Kelvin) 
         t_skin_patch            => temperature_inst%t_skin_patch           , & ! Output: [real(r8) (:)   ]  patch skin temperature (K)
         t_h2osfc                => temperature_inst%t_h2osfc_col           , & ! Input:  [real(r8) (:)   ]  surface water temperature               
         tssbef                  => temperature_inst%t_ssbef_col            , & ! Input:  [real(r8) (:,:) ]  soil/snow temperature before update   
         t_h2osfc_bef            => temperature_inst%t_h2osfc_bef_col       , & ! Input:  [real(r8) (:)   ]  saved surface water temperature         
         t_grnd                  => temperature_inst%t_grnd_col             , & ! Input:  [real(r8) (:)   ]  ground temperature (Kelvin)             
         t_soisno                => temperature_inst%t_soisno_col           , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)             
         xmf                     => temperature_inst%xmf_col                , & ! Input:  [real(r8) (:)   ]  
         xmf_h2osfc              => temperature_inst%xmf_h2osfc_col         , & ! Input:  [real(r8) (:)   ]  
         fact                    => temperature_inst%fact_col               , & ! Input:  [real(r8) (:)   ]  
         c_h2osfc                => temperature_inst%c_h2osfc_col           , & ! Input:  [real(r8) (:)   ]  

         htvp                    => energyflux_inst%htvp_col                , & ! Input:  [real(r8) (:)   ]  latent heat of vapor of water (or sublimation) [j/kg]
         eflx_building_heat_errsoi=> energyflux_inst%eflx_building_heat_errsoi_col  , & ! Input: [real(r8) (:)] heat flux to interior surface of walls and roof for errsoi check (W m-2)
         eflx_wasteheat_patch    => energyflux_inst%eflx_wasteheat_patch    , & ! Input:  [real(r8) (:)   ]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
         eflx_ventilation_patch  => energyflux_inst%eflx_ventilation_patch  , & ! Input:  [real(r8) (:)   ]  sensible heat flux from building ventilation (W/m**2)
         eflx_heat_from_ac_patch => energyflux_inst%eflx_heat_from_ac_patch , & ! Input:  [real(r8) (:)   ]  sensible heat flux put back into canyon due to removal by AC (W/m**2)
         eflx_traffic_patch      => energyflux_inst%eflx_traffic_patch      , & ! Input:  [real(r8) (:)   ]  traffic sensible heat flux (W/m**2)     
         dlrad                   => energyflux_inst%dlrad_patch             , & ! Input:  [real(r8) (:)   ]  downward longwave radiation below the canopy [W/m2]
         ulrad                   => energyflux_inst%ulrad_patch             , & ! Input:  [real(r8) (:)   ]  upward longwave radiation above the canopy [W/m2]
         cgrnds                  => energyflux_inst%cgrnds_patch            , & ! Input:  [real(r8) (:)   ]  deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
         cgrndl                  => energyflux_inst%cgrndl_patch            , & ! Input:  [real(r8) (:)   ]  deriv of soil latent heat flux wrt soil temp [w/m**2/k]
         
         qflx_evap_can           => waterfluxbulk_inst%qflx_evap_can_patch      , & ! Output: [real(r8) (:)   ]  evaporation from leaves and stems (mm H2O/s) (+ = to atm)
         qflx_evap_soi           => waterfluxbulk_inst%qflx_evap_soi_patch      , & ! Output: [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)
         qflx_evap_veg           => waterfluxbulk_inst%qflx_evap_veg_patch      , & ! Output: [real(r8) (:)   ]  vegetation evaporation (mm H2O/s) (+ = to atm)
         qflx_tran_veg           => waterfluxbulk_inst%qflx_tran_veg_patch      , & ! Input:  [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)
         qflx_evap_tot           => waterfluxbulk_inst%qflx_evap_tot_patch      , & ! Output: [real(r8) (:)   ]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg
         qflx_liqevap_from_top_layer   => waterfluxbulk_inst%qflx_liqevap_from_top_layer_patch  , & ! Output: [real(r8) (:)   ]  rate of liquid water evaporated from top soil or snow layer (mm H2O/s) [+]
         qflx_solidevap_from_top_layer => waterfluxbulk_inst%qflx_solidevap_from_top_layer_patch, & ! Output: [real(r8) (:)   ]  rate of ice evaporated from top soil or snow layer (sublimation) (mm H2O /s) [+]
         qflx_liqdew_to_top_layer      => waterfluxbulk_inst%qflx_liqdew_to_top_layer_patch     , & ! Output: [real(r8) (:)   ]  rate of liquid water deposited on top soil or snow layer (dew) (mm H2O /s) [+]
         qflx_soliddew_to_top_layer    => waterfluxbulk_inst%qflx_soliddew_to_top_layer_patch   , & ! Output: [real(r8) (:)   ]  rate of solid water deposited on top soil or snow layer (frost) (mm H2O /s) [+]
         qflx_ev_snow            => waterfluxbulk_inst%qflx_ev_snow_patch       , & ! In/Out: [real(r8) (:)   ]  evaporation flux from snow (mm H2O/s) [+ to atm]
         qflx_ev_soil            => waterfluxbulk_inst%qflx_ev_soil_patch       , & ! In/Out: [real(r8) (:)   ]  evaporation flux from soil (mm H2O/s) [+ to atm]
         qflx_ev_h2osfc          => waterfluxbulk_inst%qflx_ev_h2osfc_patch     , & ! In/Out: [real(r8) (:)   ]  evaporation flux from h2osfc (mm H2O/s) [+ to atm]
         ! [PORTED by Hui Tang: NVP evaporation flux linearization correction]
         qflx_ev_nvp             => waterfluxbulk_inst%qflx_ev_nvp_patch        , & ! In/Out: [real(r8) (:)   ]  evaporation flux from NVP (mm H2O/s) [+ to atm]
         
         eflx_sh_grnd            => energyflux_inst%eflx_sh_grnd_patch      , & ! Output: [real(r8) (:)   ]  sensible heat flux from ground (W/m**2) [+ to atm]
         eflx_sh_veg             => energyflux_inst%eflx_sh_veg_patch       , & ! Output: [real(r8) (:)   ]  sensible heat flux from leaves (W/m**2) [+ to atm]
         eflx_soil_grnd          => energyflux_inst%eflx_soil_grnd_patch    , & ! Output: [real(r8) (:)   ]  soil heat flux (W/m**2) [+ = into soil] 
         eflx_soil_grnd_u        => energyflux_inst%eflx_soil_grnd_u_patch  , & ! Output: [real(r8) (:)   ]  urban soil heat flux (W/m**2) [+ = into soil]
         eflx_soil_grnd_r        => energyflux_inst%eflx_soil_grnd_r_patch  , & ! Output: [real(r8) (:)   ]  rural soil heat flux (W/m**2) [+ = into soil]
         eflx_sh_tot             => energyflux_inst%eflx_sh_tot_patch       , & ! Output: [real(r8) (:)   ]  total sensible heat flux (W/m**2) [+ to atm]
         eflx_sh_tot_u           => energyflux_inst%eflx_sh_tot_u_patch     , & ! Output: [real(r8) (:)   ]  urban total sensible heat flux (W/m**2) [+ to atm]
         eflx_sh_tot_r           => energyflux_inst%eflx_sh_tot_r_patch     , & ! Output: [real(r8) (:)   ]  rural total sensible heat flux (W/m**2) [+ to atm]
         eflx_lh_tot             => energyflux_inst%eflx_lh_tot_patch       , & ! Output: [real(r8) (:)   ]  total latent heat flux (W/m**2)  [+ to atm]
         eflx_lh_tot_u           => energyflux_inst%eflx_lh_tot_u_patch     , & ! Output: [real(r8) (:)   ]  urban total latent heat flux (W/m**2)  [+ to atm]
         eflx_lh_tot_r           => energyflux_inst%eflx_lh_tot_r_patch     , & ! Output: [real(r8) (:)   ]  rural total latent heat flux (W/m**2)  [+ to atm]
         eflx_lwrad_out          => energyflux_inst%eflx_lwrad_out_patch    , & ! Output: [real(r8) (:)   ]  emitted infrared (longwave) radiation (W/m**2)
         eflx_lwrad_net          => energyflux_inst%eflx_lwrad_net_patch    , & ! Output: [real(r8) (:)   ]  net infrared (longwave) rad (W/m**2) [+ = to atm]
         eflx_lwrad_net_r        => energyflux_inst%eflx_lwrad_net_r_patch  , & ! Output: [real(r8) (:)   ]  rural net infrared (longwave) rad (W/m**2) [+ = to atm]
         eflx_lwrad_out_r        => energyflux_inst%eflx_lwrad_out_r_patch  , & ! Output: [real(r8) (:)   ]  rural emitted infrared (longwave) rad (W/m**2)
         eflx_lwrad_net_u        => energyflux_inst%eflx_lwrad_net_u_patch  , & ! Output: [real(r8) (:)   ]  urban net infrared (longwave) rad (W/m**2) [+ = to atm]
         eflx_lwrad_out_u        => energyflux_inst%eflx_lwrad_out_u_patch  , & ! Output: [real(r8) (:)   ]  urban emitted infrared (longwave) rad (W/m**2)
         eflx_lh_vege            => energyflux_inst%eflx_lh_vege_patch      , & ! Output: [real(r8) (:)   ]  veg evaporation heat flux (W/m**2) [+ to atm]
         eflx_lh_vegt            => energyflux_inst%eflx_lh_vegt_patch      , & ! Output: [real(r8) (:)   ]  veg transpiration heat flux (W/m**2) [+ to atm]
         eflx_lh_grnd            => energyflux_inst%eflx_lh_grnd_patch      , & ! Output: [real(r8) (:)   ]  ground evaporation heat flux (W/m**2) [+ to atm]
         errsoi_col              => energyflux_inst%errsoi_col              , & ! Output: [real(r8) (:)   ]  column-level soil/lake energy conservation error (W/m**2)
         errsoi_patch            => energyflux_inst%errsoi_patch              & ! Output: [real(r8) (:)   ]  patch-level soil/lake energy conservation error (W/m**2)
         )

      ! Get step size

      dtime = get_step_size_real()

      call t_startf('bgp2_loop_1')
      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         j = col%snl(c)+1

         ! Calculate difference in soil temperature from last time step, for
         ! flux corrections

         if (col%snl(c) < 0) then
            ! [PORTED by Hui Tang: Phase 1c RESTORE (2026-06-11) — NVP-weighted t_grnd0 for snl<0,
            !  mirroring the BiogeophysPreFluxCalcsMod snl<0 restore. tinc = t_grnd - t_grnd0 must use
            !  the NVP-weighted blend on both sides so the LW linearization emg*sb*t_grnd0^3*4*tinc in
            !  eflx_soil_grnd matches the solve (which now applies the NVP surface flux at j=0).]
            if (use_nvp .and. col%nvp_layer_active(c)) then
               frac_nvp_eff = min(1._r8 - frac_h2osfc(c) - frac_sno_eff(c), &
                                  max(0._r8, col%frac_nvp(c) - frac_sno_eff(c)))
               frac_soil    = max(0._r8, 1._r8 - frac_sno_eff(c) - frac_h2osfc(c) - frac_nvp_eff)
               t_grnd0(c) = frac_sno_eff(c) * tssbef(c,col%snl(c)+1) &
                           + frac_nvp_eff    * tssbef(c,0) &
                           + frac_soil       * tssbef(c,1) &
                           + frac_h2osfc(c)  * t_h2osfc_bef(c)
            else
               t_grnd0(c) = frac_sno_eff(c) * tssbef(c,col%snl(c)+1) &
                    + (1 - frac_sno_eff(c) - frac_h2osfc(c)) * tssbef(c,1) &
                    + frac_h2osfc(c) * t_h2osfc_bef(c)
            end if
         else
            ! [PORTED by Hui Tang: include NVP layer temperature in t_grnd0 for snl==0 NVP columns.
            !  Mirrors BiogeophysPreFluxCalcsMod snl==0 branch so tinc = t_grnd - t_grnd0 is consistent.]
            if (use_nvp .and. col%nvp_layer_active(c)) then
               frac_nvp_eff = min(1._r8 - frac_h2osfc(c) - frac_sno_eff(c), &
                                  max(0._r8, col%frac_nvp(c) - frac_sno_eff(c)))
               frac_soil    = max(0._r8, 1._r8 - frac_sno_eff(c) - frac_h2osfc(c) - frac_nvp_eff)
               t_grnd0(c) = frac_nvp_eff   * tssbef(c,0) &
                           + frac_soil      * tssbef(c,1) &
                           + frac_h2osfc(c) * t_h2osfc_bef(c)
            else
               t_grnd0(c) = (1 - frac_h2osfc(c)) * tssbef(c,1) + frac_h2osfc(c) * t_h2osfc_bef(c)
            end if
         endif

         tinc(c) = t_grnd(c) - t_grnd0(c)

      end do

      ! Correct fluxes to present soil temperature

      do fp = 1,num_nolakep
         p = filter_nolakep(fp)
         c = patch%column(p)
         eflx_sh_grnd(p) = eflx_sh_grnd(p) + tinc(c)*cgrnds(p)
         qflx_evap_soi(p) = qflx_evap_soi(p) + tinc(c)*cgrndl(p)

         ! Set ev_soil, ev_h2osfc, ev_snow for urban landunits here
         l = patch%landunit(p)
         if (lun%urbpoi(l)) then
            qflx_ev_soil(p) = 0._r8
            qflx_ev_h2osfc(p) = 0._r8
            qflx_ev_nvp(p) = 0._r8
            qflx_ev_snow(p) = qflx_evap_soi(p)
         else
            ! [PORTED by Hui Tang: NaN diagnostic — identify which term makes qflx_ev_snow NaN]
            if ((qflx_ev_snow(p) /= qflx_ev_snow(p)) .or. &
                (tinc(c)*cgrndl(p) /= tinc(c)*cgrndl(p))) then
               write(iulog,*) "NaN DIAGNOSTIC SoilFluxesMod: before qflx_ev_snow linearization"
               write(iulog,*) "  p, c                  = ", p, c
               write(iulog,*) "  qflx_ev_snow(p)       = ", qflx_ev_snow(p)
               write(iulog,*) "  tinc(c)                = ", tinc(c)
               write(iulog,*) "  cgrndl(p)              = ", cgrndl(p)
               write(iulog,*) "  tinc*cgrndl            = ", tinc(c)*cgrndl(p)
               write(iulog,*) "  t_grnd(c)              = ", t_grnd(c)
               write(iulog,*) "  t_grnd0(c)             = ", t_grnd0(c)
               write(iulog,*) "  tssbef(c,snl+1)        = ", tssbef(c,col%snl(c)+1)
               write(iulog,*) "  frac_sno_eff(c)        = ", frac_sno_eff(c)
               !call endrun(subgrid_index=p, subgrid_level=subgrid_level_patch, &
               !     msg="NaN in qflx_ev_snow or tinc*cgrndl in SoilFluxesMod")
            end if
            qflx_ev_snow(p) = qflx_ev_snow(p) + tinc(c)*cgrndl(p)
            qflx_ev_soil(p) = qflx_ev_soil(p) + tinc(c)*cgrndl(p)
            qflx_ev_h2osfc(p) = qflx_ev_h2osfc(p) + tinc(c)*cgrndl(p)
            ! [PORTED by Hui Tang: apply linearization correction to NVP evaporation diagnostic]
            ! Skip when NVP is fully buried (frac_nvp_eff <= 0): qflx_ev_nvp was zeroed in
            ! BareGroundFluxesMod/CanopyFluxesMod and must remain zero to avoid a water
            ! balance error (non-zero qflx_ev_nvp_col with no corresponding water removal).
            ! [PORTED by Hui Tang: gate on exposed NVP fraction instead of the binary snl<-1, so
            !  partial snow cover keeps the correction wherever NVP is still exposed. frac_nvp_eff
            !  (frac_sno_eff based, matching this module) is computed locally since it is not yet set.]
            frac_nvp_eff = min(1._r8 - frac_h2osfc(c) - frac_sno_eff(c), max(0._r8, &
                               col%frac_nvp(c) - frac_sno_eff(c)))
            if (use_nvp .and. frac_nvp_eff <= 0._r8) then
               qflx_ev_nvp(p) = 0._r8
            else
               ! [PORTED by Hui Tang: linearize qflx_ev_nvp with the NVP layer's OWN
               !  temperature increment instead of the bulk tinc. tinc = t_grnd(post, NVP-weighted)
               !  - t_grnd0(pre, soil-only) conflates the temporal increment with a soil-vs-NVP basis
               !  difference (~frac_nvp_eff*(t_nvp-t_soil)); for the thin, thermally-decoupled moss
               !  that spurious term is several K and can flip qflx_ev_nvp negative (unphysical dew in
               !  summer). When NVP is active (jbot_sno=-1) layer 0 IS the NVP layer for ALL snow
               !  states (snow bottoms at j=-1, not 0), so t_soisno(c,0)/tssbef(c,0) are the NVP
               !  post/pre-solve temperatures (t_nvp_col = t_soisno(c,0), SoilTemperatureMod:597).
               !  cgrndl (bulk raiw*dqgdT) is retained — only the increment is corrected. Non-NVP
               !  columns keep the standard bulk tinc correction.]

               ! The correct correction uses the NVP layer's own temperature increment and derivative:
               ! qflx_ev_nvp += (t_soisno(c,0) − tssbef(c,0)) · cgrndl_nvp, where cgrndl_nvp = raiw_nvp·hr_nvp·qsatgdT_nvp.

               qflx_ev_nvp(p) = qflx_ev_nvp(p) + (t_soisno(c,0) - tssbef(c,0))*cgrndl(p)

            end if
         endif
      end do

      ! Partition evaporation into liquid and solid
      do fp = 1, num_nolakep
         p = filter_nolakep(fp)
         c = patch%column(p)
         l = patch%landunit(p)
         j = col%snl(c)+1

         qflx_liqevap_from_top_layer(p)   = 0._r8
         qflx_solidevap_from_top_layer(p) = 0._r8
         qflx_soliddew_to_top_layer(p)    = 0._r8
         qflx_liqdew_to_top_layer(p)      = 0._r8

         ! Partition the evaporation from snow/soil surface into liquid evaporation,
         ! solid evaporation (sublimation), liquid dew, or solid dew.  Note that the variables
         ! affected here are all related to the snow subgrid patch only because of the use of qflx_ev_snow.
         ! In the situations where there are snow layers or there is snow without an explicit snow layer,
         ! the partitioned variables will represent the components of snow evaporation
         ! (qflx_ev_snow = qflx_liqevap_from_top_layer + qflx_solidevap_from_top_layer
         ! - qflx_liqdew_to_top_layer - qflx_soliddew_to_top_layer).
         ! In the case of no snow, qflx_ev_snow has already been set equal to qflx_ev_soil (the evaporation
         ! from the subgrid soil patch) and the partitioned variables will then represent evaporation from the
         ! subgrid soil patch.
         ! In the case of urban columns (and lake columns - see LakeHydrologyMod), there are no subgrid
         ! patches and qflx_evap_soi is used. qflx_evap_soi = qflx_liqevap_from_top_layer + qflx_solidevap_from_top_layer
         ! - qflx_liqdew_to_top_layer - qflx_soliddew_to_top_layer.
         if (.not. lun%urbpoi(l)) then
            if (qflx_ev_snow(p) >= 0._r8) then
               ! for evaporation partitioning between liquid evap and ice sublimation,
               ! use the ratio of liquid to (liquid+ice) in the top layer to determine split
               if ((h2osoi_liq(c,j)+h2osoi_ice(c,j)) > 0._r8) then
                  qflx_liqevap_from_top_layer(p) = max(qflx_ev_snow(p)*(h2osoi_liq(c,j)/ &
                       (h2osoi_liq(c,j)+h2osoi_ice(c,j))), 0._r8)
               else
                  qflx_liqevap_from_top_layer(p) = 0._r8
               end if
               qflx_solidevap_from_top_layer(p) = qflx_ev_snow(p) - qflx_liqevap_from_top_layer(p)
            else
               if (t_grnd(c) < tfrz) then
                  qflx_soliddew_to_top_layer(p) = abs(qflx_ev_snow(p))
               else
                  qflx_liqdew_to_top_layer(p) = abs(qflx_ev_snow(p))
               end if
            end if

         else ! Urban columns

            if (qflx_evap_soi(p) >= 0._r8) then
               ! for evaporation partitioning between liquid evap and ice sublimation,
               ! use the ratio of liquid to (liquid+ice) in the top layer to determine split
               if ((h2osoi_liq(c,j)+h2osoi_ice(c,j)) > 0._r8) then
                  qflx_liqevap_from_top_layer(p) = max(qflx_evap_soi(p)*(h2osoi_liq(c,j)/ &
                       (h2osoi_liq(c,j)+h2osoi_ice(c,j))), 0._r8)
               else
                  qflx_liqevap_from_top_layer(p) = 0._r8
               end if
               qflx_solidevap_from_top_layer(p) = qflx_evap_soi(p) - qflx_liqevap_from_top_layer(p)
            else
               if (t_grnd(c) < tfrz) then
                  qflx_soliddew_to_top_layer(p) = abs(qflx_evap_soi(p))
               else
                  qflx_liqdew_to_top_layer(p) = abs(qflx_evap_soi(p))
               end if
            end if

         end if

      end do

      ! Constrain evaporation from snow to be <= available moisture
      do fp = 1,num_nolakep
         p = filter_nolakep(fp)
         c = patch%column(p)
         j = col%snl(c)+1
         ! snow layers; assumes for j < 1 that frac_sno_eff > 0
         if (j < 1) then
            ! Defining the limitation uniformly for all patches is more 
            ! strict than absolutely necessary.  This definition assumes 
            ! each patch is spatially distinct and may remove all the snow
            ! on its patch, but may not remove snow from adjacent patches. 
            evaporation_limit = (h2osoi_ice(c,j)+h2osoi_liq(c,j))/(frac_sno_eff(c)*dtime)
            if (qflx_ev_snow(p) > evaporation_limit) then
               evaporation_demand = qflx_ev_snow(p)
               qflx_ev_snow(p)    = evaporation_limit
               qflx_evap_soi(p)   = qflx_evap_soi(p) - frac_sno_eff(c)*(evaporation_demand - evaporation_limit)
               qflx_liqevap_from_top_layer(p)   = max(h2osoi_liq(c,j)/(frac_sno_eff(c)*dtime), 0._r8)
               qflx_solidevap_from_top_layer(p) = max(h2osoi_ice(c,j)/(frac_sno_eff(c)*dtime), 0._r8)
               ! conserve total energy flux
               eflx_sh_grnd(p) = eflx_sh_grnd(p) + frac_sno_eff(c)*(evaporation_demand - evaporation_limit)*htvp(c)
            endif
         endif
         
         ! top soil layer for urban columns (excluding pervious road, which 
         ! shouldn't be limited here b/c it uses the uses the soilwater
         ! equations, while the other urban columns do not)
         if (lun%urbpoi(patch%landunit(p)) .and. (col%itype(c)/=icol_road_perv) .and. (j == 1)) then
            evaporation_limit = (h2osoi_ice(c,j)+h2osoi_liq(c,j))/dtime
            if (qflx_evap_soi(p) > evaporation_limit) then
               evaporation_demand = qflx_evap_soi(p)
               qflx_evap_soi(p)   = evaporation_limit
               qflx_ev_snow(p)    = qflx_evap_soi(p)
               qflx_liqevap_from_top_layer(p)   = max(h2osoi_liq(c,j)/dtime, 0._r8)
               qflx_solidevap_from_top_layer(p) = max(h2osoi_ice(c,j)/dtime, 0._r8)
               ! conserve total energy flux
               eflx_sh_grnd(p) = eflx_sh_grnd(p) +(evaporation_demand -evaporation_limit)*htvp(c)
            endif
         endif
         
         ! limit only solid evaporation (sublimation) from top soil layer
         ! (liquid evaporation from soil should not be limited)
         if (j==1 .and. frac_h2osfc(c) < 1._r8) then
            evaporation_limit = h2osoi_ice(c,j)/(dtime*(1._r8 - frac_h2osfc(c)))
            if (qflx_solidevap_from_top_layer(p) >= evaporation_limit) then
               evaporation_demand = qflx_solidevap_from_top_layer(p)
               qflx_solidevap_from_top_layer(p) &
                    = evaporation_limit
               qflx_liqevap_from_top_layer(p)  &
                    = qflx_liqevap_from_top_layer(p)  &
                    + (evaporation_demand - evaporation_limit)
            endif
         endif

      enddo
      
      call t_stopf('bgp2_loop_1')
      call t_startf('bgp2_loop_2')

      do fp = 1,num_nolakep
         p = filter_nolakep(fp)
         c = patch%column(p)
         l = patch%landunit(p)
         g = patch%gridcell(p)
         j = col%snl(c)+1

         ! Ground heat flux
         
         if (.not. lun%urbpoi(l)) then
            ! [PORTED by Hui Tang: fix lw_grnd for NVP — area-weighted LW emission including NVP]
            ! Standard formula uses frac_sno_eff/bare-soil/h2osfc fractions summing to 1.
            ! When NVP is active (snl=0 or snl<0), NVP occupies frac_nvp_eff of the non-snow,
            ! non-water ground and emits LW from tssbef(c,0).  Four fractions sum to 1:
            !   snow (frac_sno_eff) + NVP (frac_nvp_eff) + bare soil + surface water (frac_h2osfc)
            ! When snl=0: frac_sno_eff=0 so the snow term vanishes and the formula reduces to
            !   three terms (NVP + bare soil + water).
            if (use_nvp .and. col%nvp_layer_active(c)) then
               frac_nvp_eff = min(col%frac_nvp(c), max(0._r8, 1._r8 - frac_sno_eff(c) - frac_h2osfc(c)))
               lw_grnd = frac_sno_eff(c)                                    * tssbef(c,col%snl(c)+1)**4 &
                    + frac_nvp_eff                                           * tssbef(c,0)**4 &
                    + (1._r8 - frac_sno_eff(c) - frac_nvp_eff - frac_h2osfc(c)) * tssbef(c,1)**4 &
                    + frac_h2osfc(c)                                         * t_h2osfc_bef(c)**4
            else
               lw_grnd=(frac_sno_eff(c)*tssbef(c,col%snl(c)+1)**4 &
                    +(1._r8-frac_sno_eff(c)-frac_h2osfc(c))*tssbef(c,1)**4 &
                    +frac_h2osfc(c)*t_h2osfc_bef(c)**4)
            end if

            eflx_soil_grnd(p) = ((1._r8- frac_sno_eff(c))*sabg_soil(p) + frac_sno_eff(c)*sabg_snow(p)) + dlrad(p) &
                 + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) &
                 - emg(c)*sb*lw_grnd - emg(c)*sb*t_grnd0(c)**3*(4._r8*tinc(c)) &
                 - (eflx_sh_grnd(p)+qflx_evap_soi(p)*htvp(c))
            ! [PORTED by Hui Tang: NVP solar absorption in eflx_soil_grnd]
            ! SurfaceRadiationMod removes sabg_lyr(p,0) from sabg_soil (snl=0: line 781;
            ! snl<0: sabg_soil=sabg_lyr(p,1) only); add it back so NVP absorbed solar is
            ! counted as column energy input in both the exposed (snl=0) and buried (snl<0) cases.
            ! NOTE: sabg_lyr(p,0) must NOT be multiplied by frac_nvp_eff.  NVPBeerLawAbsorptance
            ! (FatesRadiationDriveMod.F90) already includes nvp_frac in the absorptance:
            !   fabd_nvp = nvp_frac * (1 - exp(-k * lai_nvp))
            ! so sabg_lyr(p,0) = fabd_nvp * trd + fabi_nvp * tri is already per unit ground area.
            ! Applying frac_nvp_eff again would double-count the NVP coverage fraction.
            if (use_nvp .and. col%nvp_layer_active(c)) then
               eflx_soil_grnd(p) = eflx_soil_grnd(p) + sabg_lyr(p,0)
            end if

            if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
               eflx_soil_grnd_r(p) = eflx_soil_grnd(p)
            end if
         else
            ! For all urban columns we use the net longwave radiation (eflx_lwrad_net) since
            ! the term (emg*sb*tssbef(col%snl+1)**4) is not the upward longwave flux because of 
            ! interactions between urban columns.

            eflx_lwrad_del(p) = 4._r8*emg(c)*sb*t_grnd0(c)**3*tinc(c)

            ! Include transpiration term because needed for pervious road
            ! and wasteheat and traffic flux
            eflx_soil_grnd(p) = sabg(p) + dlrad(p) &
                 - eflx_lwrad_net(p) - eflx_lwrad_del(p) &
                 - (eflx_sh_grnd(p) + qflx_evap_soi(p)*htvp(c) + qflx_tran_veg(p)*hvap) &
                 + eflx_wasteheat_patch(p) + eflx_heat_from_ac_patch(p) + eflx_traffic_patch(p) &
                 + eflx_ventilation_patch(p)
            eflx_soil_grnd_u(p) = eflx_soil_grnd(p)
         end if

         ! Total fluxes (vegetation + ground)

         eflx_sh_tot(p) = eflx_sh_veg(p) + eflx_sh_grnd(p)
         if (.not. lun%urbpoi(l)) eflx_sh_tot(p) = eflx_sh_tot(p) + eflx_sh_stem(p)
         
         qflx_evap_tot(p) = qflx_evap_veg(p) + qflx_evap_soi(p)

         eflx_lh_tot(p)= hvap*qflx_evap_veg(p) + htvp(c)*qflx_evap_soi(p)
         if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
            eflx_lh_tot_r(p)= eflx_lh_tot(p)
            eflx_sh_tot_r(p)= eflx_sh_tot(p)
         else if (lun%urbpoi(l)) then
            eflx_lh_tot_u(p)= eflx_lh_tot(p)
            eflx_sh_tot_u(p)= eflx_sh_tot(p)
         end if

         ! Variables needed by history tape

         qflx_evap_can(p)  = qflx_evap_veg(p) - qflx_tran_veg(p)
         eflx_lh_vege(p)   = (qflx_evap_veg(p) - qflx_tran_veg(p)) * hvap
         eflx_lh_vegt(p)   = qflx_tran_veg(p) * hvap
         eflx_lh_grnd(p)   = qflx_evap_soi(p) * htvp(c)

      end do
      call t_stopf('bgp2_loop_2')
      call t_startf('bgp2_loop_3')

      ! Soil Energy balance check

      do fp = 1,num_nolakep
         p = filter_nolakep(fp)
         c = patch%column(p)
         errsoi_patch(p) = eflx_soil_grnd(p) - xmf(c) - xmf_h2osfc(c) &
              - frac_h2osfc(c)*(t_h2osfc(c)-t_h2osfc_bef(c)) &
              *(c_h2osfc(c)/dtime)
         errsoi_patch(p) =  errsoi_patch(p)+eflx_h2osfc_to_snow_col(c) 
         ! For urban sunwall, shadewall, and roof columns, the "soil" energy balance check
         ! must include the heat flux from the interior of the building.
         if (col%itype(c)==icol_sunwall .or. col%itype(c)==icol_shadewall .or. col%itype(c)==icol_roof) then
            errsoi_patch(p) = errsoi_patch(p) + eflx_building_heat_errsoi(c) 
         end if
      end do
      do j = -nlevsno+1,nlevgrnd
         do fp = 1,num_nolakep
            p = filter_nolakep(fp)
            c = patch%column(p)

! Do this for perv and imperv road for -nlevsno+1,nlevgrnd
            if (col%itype(c) /= icol_sunwall .and. col%itype(c) /= icol_shadewall &
                 .and. col%itype(c) /= icol_roof) then
               ! area weight heat absorbed by snow layers
               if (j >= col%snl(c)+1 .and. j < 1) errsoi_patch(p) = errsoi_patch(p) &
                    - frac_sno_eff(c)*(t_soisno(c,j)-tssbef(c,j))/fact(c,j)
               if (j >= 1) errsoi_patch(p) = errsoi_patch(p) &
                    - (t_soisno(c,j)-tssbef(c,j))/fact(c,j)
            end if
         end do
      end do

! Do this for sunwall, shadewall, roof but for -nlevsno+1,nlevurb
      do j = -nlevsno+1,nlevurb
         do fp = 1,num_urbanp
            p = filter_urbanp(fp)
            c = patch%column(p)

            if (col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall &
                 .or. col%itype(c) == icol_roof) then
            ! area weight heat absorbed by snow layers
               if (j >= col%snl(c)+1 .and. j < 1) errsoi_patch(p) = errsoi_patch(p) &
                   - frac_sno_eff(c)*(t_soisno(c,j)-tssbef(c,j))/fact(c,j)
               if (j >= 1) errsoi_patch(p) = errsoi_patch(p) &
                   - (t_soisno(c,j)-tssbef(c,j))/fact(c,j)
            end if
         end do
      end do

      ! [PORTED by Hui Tang: NVP errsoi fix]
      ! snl=0: j=0 is the atmospheric boundary (top BC = hs_nvp) and is skipped entirely
      ! by the loops above (j>=snl+1=1 fails for j=0) → subtract full heat storage.
      ! snl<0: j=0 is included in the snow loop with frac_sno_eff weight, but cv(c,0) is
      ! per unit column area (unlike snow layers where cv is per unit snow area) → subtract
      ! the missing (1-frac_sno_eff) fraction to account for the full column-area heat storage.
      if (use_nvp) then
         do fp = 1, num_nolakep
            p = filter_nolakep(fp)
            c = patch%column(p)
            if (col%nvp_layer_active(c) .and. col%snl(c) == 0) then
               ! snl=0: j=0 skipped entirely; subtract full heat storage (weight = 1.0).
               errsoi_patch(p) = errsoi_patch(p) &
                    - (t_soisno(c,0) - tssbef(c,0)) / fact(c,0)
            else if (col%nvp_layer_active(c) .and. col%snl(c) < 0) then
               ! snl<0: j=0 included with frac_sno_eff; cv(c,0) is per unit column area,
               ! so subtract the missing (1-frac_sno_eff) fraction.
               ! [PORTED by Hui Tang: NVP under-snow errsoi correction]
               errsoi_patch(p) = errsoi_patch(p) &
                    - (1.0_r8 - frac_sno_eff(c)) * (t_soisno(c,0) - tssbef(c,0)) / fact(c,0)
            end if
         end do
      end if

      ! [PORTED by Hui Tang: errsoi diagnostic - decompose terms when error is large]
      do fp = 1, num_nolakep
         p = filter_nolakep(fp)
         c = patch%column(p)
         if (abs(errsoi_patch(p)) > 0.5_r8) then
            heat_store_diag = 0._r8
            do j = -nlevsno+1, nlevgrnd
               if (col%itype(c) /= icol_sunwall .and. col%itype(c) /= icol_shadewall &
                    .and. col%itype(c) /= icol_roof) then
                  if (j >= col%snl(c)+1 .and. j < 1) heat_store_diag = heat_store_diag + &
                       frac_sno_eff(c)*(t_soisno(c,j)-tssbef(c,j))/fact(c,j)
                  if (j >= 1) heat_store_diag = heat_store_diag + &
                       (t_soisno(c,j)-tssbef(c,j))/fact(c,j)
               end if
            end do
            ! [PORTED by Hui Tang: NVP j=0 correction to heat_store_diag — mirrors errsoi fix]
            ! Weight is 1.0 not frac_nvp_eff: cv(c,0) is per unit column area so
            ! (t_0-tbef_0)/fact(0) is already W/m2 column area.
            ! snl=0: j=0 skipped by loop (j>=1 and j<1 impossible) → add full term.
            ! snl<0: j=0 included with frac_sno_eff → add missing (1-frac_sno_eff) fraction.
            if (use_nvp .and. col%nvp_layer_active(c)) then
               if (col%snl(c) == 0) then
                  heat_store_diag = heat_store_diag + &
                       (t_soisno(c,0) - tssbef(c,0)) / fact(c,0)
               else if (col%snl(c) < 0) then
                  heat_store_diag = heat_store_diag + &
                       (1.0_r8 - frac_sno_eff(c)) * (t_soisno(c,0) - tssbef(c,0)) / fact(c,0)
               end if
            end if

            write(iulog,*) '[ERRSOI DBG] p=',p,' c=',c,' snl=',col%snl(c)
            write(iulog,*) '  errsoi_patch     =', errsoi_patch(p)
            write(iulog,*) '  eflx_soil_grnd   =', eflx_soil_grnd(p)
            write(iulog,*) '  xmf (phase chg)  =', xmf(c)
            write(iulog,*) '  xmf_h2osfc       =', xmf_h2osfc(c)
            write(iulog,*) '  heat_store_sum   =', heat_store_diag
            write(iulog,*) '  eflx_h2osfc_snow =', eflx_h2osfc_to_snow_col(c)
            write(iulog,*) '  frac_h2osfc_term =', &
                 frac_h2osfc(c)*(t_h2osfc(c)-t_h2osfc_bef(c))*(c_h2osfc(c)/dtime)
            write(iulog,*) '  expected_errsoi  =', &
                 eflx_soil_grnd(p) - xmf(c) - xmf_h2osfc(c) - heat_store_diag &
                 + eflx_h2osfc_to_snow_col(c) &
                 - frac_h2osfc(c)*(t_h2osfc(c)-t_h2osfc_bef(c))*(c_h2osfc(c)/dtime)
         end if
      end do

      call t_stopf('bgp2_loop_3')
      call t_startf('bgp2_loop_4')

      ! Outgoing long-wave radiation from vegetation + ground
      ! For conservation we put the increase of ground longwave to outgoing
      ! For urban patches, ulrad=0 and (1-fracveg_nosno)=1, and eflx_lwrad_out and eflx_lwrad_net 
      ! are calculated in UrbanRadiation. The increase of ground longwave is added directly 
      ! to the outgoing longwave and the net longwave.

      do fp = 1,num_nolakep
         p = filter_nolakep(fp)
         c = patch%column(p)
         l = patch%landunit(p)
         g = patch%gridcell(p)
         j = col%snl(c)+1

         if (.not. lun%urbpoi(l)) then
            ! [PORTED by Hui Tang: NVP-aware lw_grnd for eflx_lwrad_out — mirror the eflx_soil_grnd
            !  lw_grnd (lines ~387-390). Without this, eflx_lwrad_out emits the NVP fraction at soil
            !  temperature tssbef(c,1) while eflx_soil_grnd emits it at NVP temperature tssbef(c,0),
            !  leaving errseb = -emg*sb*frac_nvp_eff*(tssbef(c,1)**4 - tssbef(c,0)**4) (~-3 W/m2 at
            !  partial snow cover). At full snow frac_nvp_eff=0 and this reduces to the standard form.]
            if (use_nvp .and. col%nvp_layer_active(c)) then
               frac_nvp_eff = min(1._r8 - frac_h2osfc(c) - frac_sno_eff(c), max(0._r8, col%frac_nvp(c) - frac_sno_eff(c)))
               lw_grnd = frac_sno_eff(c)                                    * tssbef(c,col%snl(c)+1)**4 &
                    + frac_nvp_eff                                           * tssbef(c,0)**4 &
                    + (1._r8 - frac_sno_eff(c) - frac_nvp_eff - frac_h2osfc(c)) * tssbef(c,1)**4 &
                    + frac_h2osfc(c)                                         * t_h2osfc_bef(c)**4
            else
               lw_grnd=(frac_sno_eff(c)*tssbef(c,col%snl(c)+1)**4 &
                    +(1._r8-frac_sno_eff(c)-frac_h2osfc(c))*tssbef(c,1)**4 &
                    +frac_h2osfc(c)*t_h2osfc_bef(c)**4)
            end if

            eflx_lwrad_out(p) = ulrad(p) &
                 + (1-frac_veg_nosno(p))*(1.-emg(c))*forc_lwrad(c) &
                 + (1-frac_veg_nosno(p))*emg(c)*sb*lw_grnd &
                 + 4._r8*emg(c)*sb*t_grnd0(c)**3*tinc(c)


            ! Calculate the skin temperature as a weighted sum of all the surface contributions (surface water table, snow, etc...)
            ! Note: This is the bare ground calculation of skin temperature
            !       The Urban and Vegetation are done in other place.  Urban=Later in this function Veg=CanopyFluxMod
!            t_skin_patch(p) = ((1._r8 - emv(p))*(1-frac_veg_nosno(p)) * sqrt(sqrt(lw_grnd)))  +  emv(p)*t_veg(p)
!            if( frac_veg_nosno(p).eq.0 ) then
!               t_skin_patch(p) = ((1._r8 - emv(p))*(1-frac_veg_nosno(p)) * sqrt(sqrt(lw_grnd)))  +  &
!                                           emv(p) *   frac_veg_nosno(p)  * t_veg(p)
!            end if
             if(frac_veg_nosno(p).eq.0)  t_skin_patch(p) = sqrt(sqrt(lw_grnd))

            eflx_lwrad_net(p) = eflx_lwrad_out(p) - forc_lwrad(c)
            if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
               eflx_lwrad_net_r(p) = eflx_lwrad_out(p) - forc_lwrad(c)
               eflx_lwrad_out_r(p) = eflx_lwrad_out(p)
            end if
         else
            eflx_lwrad_out(p) = eflx_lwrad_out(p) + eflx_lwrad_del(p)
            eflx_lwrad_net(p) = eflx_lwrad_net(p) + eflx_lwrad_del(p)
            eflx_lwrad_net_u(p) = eflx_lwrad_net_u(p) + eflx_lwrad_del(p)
            eflx_lwrad_out_u(p) = eflx_lwrad_out(p)
         end if
      end do

      ! lake balance for errsoi is not over pft
      ! therefore obtain column-level radiative temperature

      call p2c(bounds, num_nolakec, filter_nolakec, &
           errsoi_patch(bounds%begp:bounds%endp), &
           errsoi_col(bounds%begc:bounds%endc))

      ! Assign column-level t_soisno(snl+1) to t_skin for each urban pft
      do fp = 1, num_urbanp
         p = filter_urbanp(fp)         
         c = patch%column(p)
         
         t_skin_patch(p) = t_soisno(c,col%snl(c)+1)
  
      end do

      call t_stopf('bgp2_loop_4')

    end associate 

  end subroutine SoilFluxes

end module SoilFluxesMod

