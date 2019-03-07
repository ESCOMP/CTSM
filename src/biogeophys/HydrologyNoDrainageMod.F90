Module HydrologyNoDrainageMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate snow and soil temperatures including phase change
  !
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use clm_varctl        , only : iulog, use_vichydro, use_fates
  use clm_varcon        , only : e_ice, denh2o, denice, rpi, spval
  use CLMFatesInterfaceMod, only : hlm_fates_interface_type
  use atm2lndType       , only : atm2lnd_type
  use AerosolMod        , only : aerosol_type
  use EnergyFluxType    , only : energyflux_type
  use TemperatureType   , only : temperature_type
  use SoilHydrologyType , only : soilhydrology_type
  use SoilStateType     , only : soilstate_type
  use SaturatedExcessRunoffMod, only : saturated_excess_runoff_type
  use InfiltrationExcessRunoffMod, only : infiltration_excess_runoff_type
  use IrrigationMod, only : irrigation_type
  use WaterType , only : water_type
  use WaterFluxBulkType     , only : waterfluxbulk_type
  use WaterStateBulkType    , only : waterstatebulk_type
  use WaterDiagnosticBulkType    , only : waterdiagnosticbulk_type
  use CanopyStateType   , only : canopystate_type
  use LandunitType      , only : lun                
  use ColumnType        , only : col                
  use TopoMod, only : topo_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: CalcAndWithdrawIrrigationFluxes  ! Calculates irrigation withdrawal fluxes and withdraws from groundwater
  public  :: HydrologyNoDrainage    ! Calculates soil/snow hydrology without drainage
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CalcAndWithdrawIrrigationFluxes(bounds, &
       num_soilc, filter_soilc, &
       num_soilp, filter_soilp, &
       soilhydrology_inst, soilstate_inst, &
       irrigation_inst, &
       water_inst)
    !
    ! !DESCRIPTION:
    ! Calculates irrigation withdrawal fluxes and withdraws from groundwater
    !
    ! !USES:
    use SoilHydrologyMod       , only : WithdrawGroundwaterIrrigation
    !
    ! !ARGUMENTS:
    type(bounds_type)              , intent(in)    :: bounds
    integer                        , intent(in)    :: num_soilc            ! number of points in filter_soilc
    integer                        , intent(in)    :: filter_soilc(:)      ! column filter for soil points
    integer                        , intent(in)    :: num_soilp            ! number of points in filter_soilp
    integer                        , intent(in)    :: filter_soilp(:)      ! patch filter for soil points
    type(soilhydrology_type)       , intent(in)    :: soilhydrology_inst
    type(soilstate_type)           , intent(in)    :: soilstate_inst
    type(irrigation_type)          , intent(inout) :: irrigation_inst
    type(water_type)               , intent(inout) :: water_inst
    !
    ! !LOCAL VARIABLES:
    integer :: i  ! tracer index

    character(len=*), parameter :: subname = 'CalcAndWithdrawIrrigationFluxes'
    !-----------------------------------------------------------------------

    ! Calculate irrigation flux
    call irrigation_inst%CalcIrrigationFluxes(bounds, num_soilc, &
         filter_soilc, num_soilp, filter_soilp, &
         soilhydrology_inst, soilstate_inst, &
         water_inst)

    ! Remove groundwater irrigation
    if (irrigation_inst%UseGroundwaterIrrigation()) then
       do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
          call WithdrawGroundwaterIrrigation(bounds, num_soilc, filter_soilc, &
               water_inst%bulk_and_tracers(i)%waterflux_inst, &
               water_inst%bulk_and_tracers(i)%waterstate_inst)
       end do
    end if

  end subroutine CalcAndWithdrawIrrigationFluxes

  !-----------------------------------------------------------------------
  subroutine HydrologyNoDrainage(bounds, &
       num_nolakec, filter_nolakec, &
       num_hydrologyc, filter_hydrologyc, &
       num_urbanc, filter_urbanc, &
       num_snowc, filter_snowc, &
       num_nosnowc, filter_nosnowc, &
       clm_fates, &
       atm2lnd_inst, soilstate_inst, energyflux_inst, temperature_inst, &
       waterfluxbulk_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, &
       soilhydrology_inst, saturated_excess_runoff_inst, infiltration_excess_runoff_inst, &
       aerosol_inst, canopystate_inst, soil_water_retention_curve, topo_inst)
    !
    ! !DESCRIPTION:
    ! This is the main subroutine to execute the calculation of soil/snow
    ! hydrology
    !
    ! !USES:
    use clm_varcon           , only : denh2o, denice, hfus, grav, tfrz
    use landunit_varcon      , only : istwet, istsoil, istcrop, istdlak 
    use column_varcon        , only : icol_roof, icol_road_imperv, icol_road_perv, icol_sunwall
    use column_varcon        , only : icol_shadewall
    use clm_varctl           , only : use_cn
    use clm_varpar           , only : nlevgrnd, nlevsno, nlevsoi, nlevurb
    use clm_time_manager     , only : get_step_size, get_nstep
    use SnowHydrologyMod     , only : SnowCompaction, CombineSnowLayers, DivideSnowLayers, SnowCapping
    use SnowHydrologyMod     , only : SnowWater, BuildSnowFilter 
    use SoilHydrologyMod     , only : CLMVICMap, SetSoilWaterFractions
    use SoilHydrologyMod     , only : SetQflxInputs, RouteInfiltrationExcess, UpdateH2osfc
    use SoilHydrologyMod     , only : Infiltration, TotalSurfaceRunoff
    use SoilHydrologyMod     , only : UpdateUrbanPonding
    use SoilHydrologyMod     , only : WaterTable, PerchedWaterTable
    use SoilHydrologyMod     , only : ThetaBasedWaterTable, RenewCondensation
    use SoilWaterMovementMod , only : SoilWater 
    use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
    use SoilWaterMovementMod , only : use_aquifer_layer
    use SoilWaterPlantSinkMod , only : Compute_EffecRootFrac_And_VertTranSink

    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds               
    integer                  , intent(in)    :: num_nolakec          ! number of column non-lake points in column filter
    integer                  , intent(in)    :: filter_nolakec(:)    ! column filter for non-lake points
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
    integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
    integer                  , intent(inout) :: num_snowc            ! number of column snow points
    integer                  , intent(inout) :: filter_snowc(:)      ! column filter for snow points
    integer                  , intent(inout) :: num_nosnowc          ! number of column non-snow points
    integer                  , intent(inout) :: filter_nosnowc(:)    ! column filter for non-snow points
    type(hlm_fates_interface_type), intent(inout) :: clm_fates
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_inst
    type(soilstate_type)     , intent(inout) :: soilstate_inst
    type(energyflux_type)    , intent(in)    :: energyflux_inst
    type(temperature_type)   , intent(inout) :: temperature_inst
    type(waterfluxbulk_type)     , intent(inout) :: waterfluxbulk_inst
    type(waterstatebulk_type)    , intent(inout) :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)    , intent(inout) :: waterdiagnosticbulk_inst
    type(aerosol_type)       , intent(inout) :: aerosol_inst
    type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
    type(saturated_excess_runoff_type), intent(inout) :: saturated_excess_runoff_inst
    type(infiltration_excess_runoff_type), intent(inout) :: infiltration_excess_runoff_inst
    type(canopystate_type)   , intent(inout) :: canopystate_inst
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    class(topo_type)   , intent(in)    :: topo_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: g,l,c,j,fc                    ! indices
    real(r8) :: dtime                         ! land model time step (sec)
    real(r8) :: psi,vwc,fsattmp,psifrz        ! temporary variables for soilpsi calculation
    real(r8) :: watdry                        ! temporary
    real(r8) :: rwat(bounds%begc:bounds%endc) ! soil water wgted by depth to maximum depth of 0.5 m
    real(r8) :: swat(bounds%begc:bounds%endc) ! same as rwat but at saturation
    real(r8) :: rz(bounds%begc:bounds%endc)   ! thickness of soil layers contributing to rwat (m)
    real(r8) :: tsw                           ! volumetric soil water to 0.5 m
    real(r8) :: stsw                          ! volumetric soil water to 0.5 m at saturation
    real(r8) :: fracl                         ! fraction of soil layer contributing to 10cm total soil water
    real(r8) :: s_node                        ! soil wetness (-)
    real(r8) :: icefrac(bounds%begc:bounds%endc,1:nlevsoi)
    !-----------------------------------------------------------------------
    
    associate(                                                          & 
         z                  => col%z                                  , & ! Input:  [real(r8) (:,:) ]  layer depth  (m)                      
         dz                 => col%dz                                 , & ! Input:  [real(r8) (:,:) ]  layer thickness depth (m)             
         zi                 => col%zi                                 , & ! Input:  [real(r8) (:,:) ]  interface depth (m)                   
         snl                => col%snl                                , & ! Input:  [integer  (:)   ]  number of snow layers                    
         ctype              => col%itype                              , & ! Input:  [integer  (:)   ]  column type                              

         t_h2osfc           => temperature_inst%t_h2osfc_col          , & ! Input:  [real(r8) (:)   ]  surface water temperature               
         dTdz_top           => temperature_inst%dTdz_top_col          , & ! Output: [real(r8) (:)   ]  temperature gradient in top layer (col) [K m-1] !
         snot_top           => temperature_inst%snot_top_col          , & ! Output: [real(r8) (:)   ]  snow temperature in top layer (col) [K]  
         t_soisno           => temperature_inst%t_soisno_col          , & ! Output: [real(r8) (:,:) ]  soil temperature (Kelvin)             
         t_grnd             => temperature_inst%t_grnd_col            , & ! Output: [real(r8) (:)   ]  ground temperature (Kelvin)             
         t_grnd_u           => temperature_inst%t_grnd_u_col          , & ! Output: [real(r8) (:)   ]  Urban ground temperature (Kelvin)       
         t_grnd_r           => temperature_inst%t_grnd_r_col          , & ! Output: [real(r8) (:)   ]  Rural ground temperature (Kelvin)       
         tsl                => temperature_inst%tsl_col               , & ! Output: [real(r8) (:)   ]  temperature of near-surface soil layer (Kelvin)
         t_soi_10cm         => temperature_inst%t_soi10cm_col         , & ! Output: [real(r8) (:)   ]  soil temperature in top 10cm of soil (Kelvin)
         tsoi17             => temperature_inst%t_soi17cm_col         , & ! Output: [real(r8) (:)   ]  soil temperature in top 17cm of soil (Kelvin) 
         t_sno_mul_mss      => temperature_inst%t_sno_mul_mss_col     , & ! Output: [real(r8) (:)   ]  col snow temperature multiplied by layer mass, layer sum (K * kg/m2) 

         snow_depth         => waterdiagnosticbulk_inst%snow_depth_col         , & ! Input:  [real(r8) (:)   ]  snow height of snow covered area (m)     
         snowdp             => waterdiagnosticbulk_inst%snowdp_col             , & ! Input:  [real(r8) (:)   ]  area-averaged snow height (m)       
         frac_sno_eff       => waterdiagnosticbulk_inst%frac_sno_eff_col       , & ! Input:  [real(r8) (:)   ]  eff.  snow cover fraction (col) [frc]    
         frac_h2osfc        => waterdiagnosticbulk_inst%frac_h2osfc_col        , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)
         snw_rds            => waterdiagnosticbulk_inst%snw_rds_col            , & ! Output: [real(r8) (:,:) ]  effective snow grain radius (col,lyr) [microns, m^-6] 
         snw_rds_top        => waterdiagnosticbulk_inst%snw_rds_top_col        , & ! Output: [real(r8) (:)   ]  effective snow grain size, top layer(col) [microns] 
         sno_liq_top        => waterdiagnosticbulk_inst%sno_liq_top_col        , & ! Output: [real(r8) (:)   ]  liquid water fraction in top snow layer (col) [frc] 
         snowice            => waterdiagnosticbulk_inst%snowice_col            , & ! Output: [real(r8) (:)   ]  average snow ice lens                   
         snowliq            => waterdiagnosticbulk_inst%snowliq_col            , & ! Output: [real(r8) (:)   ]  average snow liquid water               
         snow_persistence   => waterstatebulk_inst%snow_persistence_col   , & ! Output: [real(r8) (:)   ]  counter for length of time snow-covered
         h2osoi_liqice_10cm => waterdiagnosticbulk_inst%h2osoi_liqice_10cm_col , & ! Output: [real(r8) (:)   ]  liquid water + ice lens in top 10cm of soil (kg/m2)
         h2osoi_ice         => waterstatebulk_inst%h2osoi_ice_col         , & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)                      
         h2osoi_liq         => waterstatebulk_inst%h2osoi_liq_col         , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)                  
         h2osoi_ice_tot     => waterdiagnosticbulk_inst%h2osoi_ice_tot_col     , & ! Output: [real(r8) (:)   ]  vertically summed ice lens (kg/m2)
         h2osoi_liq_tot     => waterdiagnosticbulk_inst%h2osoi_liq_tot_col     , & ! Output: [real(r8) (:)   ]  vertically summed liquid water (kg/m2)   
         h2osoi_vol         => waterstatebulk_inst%h2osoi_vol_col         , & ! Output: [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
         h2osno_top         => waterdiagnosticbulk_inst%h2osno_top_col         , & ! Output: [real(r8) (:)   ]  mass of snow in top layer (col) [kg]    
         wf                 => waterdiagnosticbulk_inst%wf_col                 , & ! Output: [real(r8) (:)   ]  soil water as frac. of whc for top 0.05 m 
         wf2                => waterdiagnosticbulk_inst%wf2_col                , & ! Output: [real(r8) (:)   ]  soil water as frac. of whc for top 0.17 m 

         watsat             => soilstate_inst%watsat_col              , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         sucsat             => soilstate_inst%sucsat_col              , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)             
         bsw                => soilstate_inst%bsw_col                 , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"              
         smp_l              => soilstate_inst%smp_l_col               , & ! Input:  [real(r8) (:,:) ]  soil matrix potential [mm]            
         smpmin             => soilstate_inst%smpmin_col              , & ! Input:  [real(r8) (:)   ]  restriction for min of soil potential (mm)
         soilpsi            => soilstate_inst%soilpsi_col               & ! Output: [real(r8) (:,:) ]  soil water potential in each soil layer (MPa)
         )

      ! Determine step size

      dtime = get_step_size()

      ! Determine initial snow/no-snow filters (will be modified possibly by
      ! routines CombineSnowLayers and DivideSnowLayers below

      call BuildSnowFilter(bounds, num_nolakec, filter_nolakec, &
           num_snowc, filter_snowc, num_nosnowc, filter_nosnowc)

      ! Determine the change of snow mass and the snow water onto soil

      call SnowWater(bounds, num_snowc, filter_snowc, num_nosnowc, filter_nosnowc, &
           atm2lnd_inst, waterfluxbulk_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, aerosol_inst)

      ! mapping soilmoist from CLM to VIC layers for runoff calculations
      if (use_vichydro) then
         call CLMVICMap(bounds, num_hydrologyc, filter_hydrologyc, &
              soilhydrology_inst, waterstatebulk_inst)
      end if

      call SetSoilWaterFractions(bounds, num_hydrologyc, filter_hydrologyc, &
           soilhydrology_inst, soilstate_inst, waterstatebulk_inst)

      call saturated_excess_runoff_inst%SaturatedExcessRunoff(&
           bounds, num_hydrologyc, filter_hydrologyc, lun, col, &
           soilhydrology_inst, soilstate_inst, waterfluxbulk_inst)

      call SetQflxInputs(bounds, num_hydrologyc, filter_hydrologyc, &
           waterfluxbulk_inst, waterdiagnosticbulk_inst)

      call infiltration_excess_runoff_inst%InfiltrationExcessRunoff( &
           bounds, num_hydrologyc, filter_hydrologyc, &
           soilhydrology_inst, soilstate_inst, saturated_excess_runoff_inst, waterfluxbulk_inst, &
           waterdiagnosticbulk_inst)

      call RouteInfiltrationExcess(bounds, num_hydrologyc, filter_hydrologyc, &
           waterfluxbulk_inst, soilhydrology_inst)

      call UpdateH2osfc(bounds, num_hydrologyc, filter_hydrologyc, &
           infiltration_excess_runoff_inst, &
           energyflux_inst, soilhydrology_inst, &
           waterfluxbulk_inst, waterstatebulk_inst, waterdiagnosticbulk_inst)

      call Infiltration(bounds, num_hydrologyc, filter_hydrologyc, &
           waterfluxbulk_inst)

      call TotalSurfaceRunoff(bounds, num_hydrologyc, filter_hydrologyc, &
           num_urbanc, filter_urbanc, &
           waterfluxbulk_inst, soilhydrology_inst, &
           waterstatebulk_inst)

      call UpdateUrbanPonding(bounds, num_urbanc, filter_urbanc, &
           waterstatebulk_inst, soilhydrology_inst, waterfluxbulk_inst)

      call Compute_EffecRootFrac_And_VertTranSink(bounds, num_hydrologyc, &
           filter_hydrologyc, soilstate_inst, canopystate_inst, waterfluxbulk_inst, energyflux_inst)
      
      if ( use_fates ) then
         call clm_fates%ComputeRootSoilFlux(bounds, num_hydrologyc, filter_hydrologyc, soilstate_inst, waterfluxbulk_inst)
      end if
      
      call SoilWater(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc, &
           soilhydrology_inst, soilstate_inst, waterfluxbulk_inst, waterstatebulk_inst, temperature_inst, &
           canopystate_inst, energyflux_inst, soil_water_retention_curve)

      if (use_vichydro) then
         ! mapping soilmoist from CLM to VIC layers for runoff calculations
         call CLMVICMap(bounds, num_hydrologyc, filter_hydrologyc, &
              soilhydrology_inst, waterstatebulk_inst)
      end if

      if (use_aquifer_layer()) then
         call WaterTable(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc, &
              soilhydrology_inst, soilstate_inst, temperature_inst, waterstatebulk_inst, &
              waterdiagnosticbulk_inst, waterfluxbulk_inst)
      else

         call PerchedWaterTable(bounds, num_hydrologyc, filter_hydrologyc, &
              num_urbanc, filter_urbanc, soilhydrology_inst, soilstate_inst, &
              temperature_inst, waterstatebulk_inst, waterfluxbulk_inst) 

         call ThetaBasedWaterTable(bounds, num_hydrologyc, filter_hydrologyc, &
              num_urbanc, filter_urbanc, soilhydrology_inst, soilstate_inst, &
              waterstatebulk_inst, waterfluxbulk_inst) 

         call RenewCondensation(bounds, num_hydrologyc, filter_hydrologyc, &
              num_urbanc, filter_urbanc,&
              soilhydrology_inst, soilstate_inst, &
              waterstatebulk_inst, waterdiagnosticbulk_inst, waterfluxbulk_inst)
         
      endif

      ! Snow capping
      call SnowCapping(bounds, num_nolakec, filter_nolakec, num_snowc, filter_snowc, &
           aerosol_inst, waterfluxbulk_inst, waterstatebulk_inst, topo_inst)

      ! Natural compaction and metamorphosis.
      call SnowCompaction(bounds, num_snowc, filter_snowc, &
           temperature_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, atm2lnd_inst)

      ! Combine thin snow elements
      call CombineSnowLayers(bounds, num_snowc, filter_snowc, &
           aerosol_inst, temperature_inst, waterfluxbulk_inst, waterstatebulk_inst, waterdiagnosticbulk_inst)

      ! Divide thick snow elements
      call DivideSnowLayers(bounds, num_snowc, filter_snowc, &
           aerosol_inst, temperature_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, is_lake=.false.)

      ! Set empty snow layers to zero
      do j = -nlevsno+1,0
         do fc = 1, num_snowc
            c = filter_snowc(fc)
            if (j <= snl(c) .and. snl(c) > -nlevsno) then
               h2osoi_ice(c,j) = 0._r8
               h2osoi_liq(c,j) = 0._r8
               t_soisno(c,j)  = 0._r8
               dz(c,j)    = 0._r8
               z(c,j)     = 0._r8
               zi(c,j-1)  = 0._r8
            end if
         end do
      end do
       
      ! Build new snow filter

      call BuildSnowFilter(bounds, num_nolakec, filter_nolakec, &
           num_snowc, filter_snowc, num_nosnowc, filter_nosnowc)

      ! For columns where snow exists, accumulate 'time-covered-by-snow' counters.
      ! Otherwise, re-zero counter, since it is bareland

      do fc = 1, num_snowc
         c = filter_snowc(fc)
         snow_persistence(c) = snow_persistence(c) + dtime
      end do
      do fc = 1, num_nosnowc
         c = filter_nosnowc(fc)
         snow_persistence(c) = 0._r8
      enddo

      ! Vertically average t_soisno and sum of h2osoi_liq and h2osoi_ice
      ! over all snow layers for history output

      do fc = 1, num_nolakec
         c = filter_nolakec(fc)
         snowice(c) = 0._r8
         snowliq(c) = 0._r8
      end do

      do j = -nlevsno+1, 0
         do fc = 1, num_snowc
            c = filter_snowc(fc)
            if (j >= snl(c)+1) then
               snowice(c) = snowice(c) + h2osoi_ice(c,j)
               snowliq(c) = snowliq(c) + h2osoi_liq(c,j)
            end if
         end do
      end do

      ! Calculate column average snow depth
      do c = bounds%begc,bounds%endc
         snowdp(c) = snow_depth(c) * frac_sno_eff(c)
      end do

      ! Calculate snow internal temperature
      !
      ! Snow internal (or: integrated) temperature is the average temperature of the entire 
      ! snowpack, weighted by the layer mass. In a formula this reads:
      !
      ! SIT = [ Sum_t Sum_i w(t,i) * T(t,i) ] / 
      !       [ Sum_t Sum_i w(t,i) ]
      !
      ! where
      !
      ! t = time
      ! i = layer index
      ! w(t,i) = layer mass or weight (kg /m2) 
      ! T(t,i) = layer temperature (K)
      ! 
      ! SIT can be calculated offline from two components: the nominator and denominator, which are output
      ! separately.
      ! 
      ! Both the nominator and denominator are scaled with a factor 1/Nstep, the number of time samples, 
      ! to make them independent of the number of time steps that were used in the averaging. 
      ! This is simply implemented as avgflag='A' in any calls to the history output routines. 
      !
      ! Snow packs without layers are not taken into account, as they have no temperature.
      !
      ! The denominator Sum_t Sum_i w(t,i) is already output as SNOWICE and SNOWLIQ (mass of 
      ! snow in layered snowpacks). Note that these != H2OSNO which does account for layerless snowpacks.
      !
      ! The nominator Sum_t Sum_i w(t,i) * T(t,i) is computed and stored as t_sno_mul_mss

      do fc = 1, num_nolakec
         c = filter_nolakec(fc)
         t_sno_mul_mss(c) = 0._r8
      end do

      do j = -nlevsno+1, 0
         do fc = 1, num_snowc
            c = filter_snowc(fc)
            if (j >= snl(c)+1) then
               t_sno_mul_mss(c) = t_sno_mul_mss(c) + h2osoi_ice(c,j) * t_soisno(c,j)
               t_sno_mul_mss(c) = t_sno_mul_mss(c) + h2osoi_liq(c,j) * tfrz
            end if
         end do
      end do

      ! Determine ground temperature, ending water balance and volumetric soil water
      ! Calculate temperature of near-surface soil layer
      ! Calculate soil temperature and total water (liq+ice) in top 10cm of soil
      ! Calculate soil temperature and total water (liq+ice) in top 17cm of soil
      do fc = 1, num_nolakec
         c = filter_nolakec(fc)
         l = col%landunit(c)
         if (.not. lun%urbpoi(l)) then
            t_soi_10cm(c) = 0._r8
            tsoi17(c) = 0._r8
            h2osoi_liqice_10cm(c) = 0._r8
            h2osoi_liq_tot(c) = 0._r8
            h2osoi_ice_tot(c) = 0._r8
         end if
      end do
      do j = 1, nlevsoi
         do fc = 1, num_nolakec
            c = filter_nolakec(fc)
            l = col%landunit(c)
            if (.not. lun%urbpoi(l)) then
               if (j == 1) then
                  tsl(c) = t_soisno(c,j)
               end if
               ! soil T at top 17 cm added by F. Li and S. Levis
               if (zi(c,j) <= 0.17_r8) then
                  fracl = 1._r8
                  tsoi17(c) = tsoi17(c) + t_soisno(c,j)*dz(c,j)*fracl
               else
                  if (zi(c,j) > 0.17_r8 .and. zi(c,j-1) < 0.17_r8) then 
                     fracl = (0.17_r8 - zi(c,j-1))/dz(c,j)
                     tsoi17(c) = tsoi17(c) + t_soisno(c,j)*dz(c,j)*fracl
                  end if
               end if

               if (zi(c,j) <= 0.1_r8) then
                  fracl = 1._r8
                  t_soi_10cm(c) = t_soi_10cm(c) + t_soisno(c,j)*dz(c,j)*fracl
                  h2osoi_liqice_10cm(c) = h2osoi_liqice_10cm(c) + &
                       (h2osoi_liq(c,j)+h2osoi_ice(c,j))* &
                       fracl
               else
                  if (zi(c,j) > 0.1_r8 .and. zi(c,j-1) < 0.1_r8) then
                     fracl = (0.1_r8 - zi(c,j-1))/dz(c,j)
                     t_soi_10cm(c) = t_soi_10cm(c) + t_soisno(c,j)*dz(c,j)*fracl
                     h2osoi_liqice_10cm(c) = h2osoi_liqice_10cm(c) + &
                          (h2osoi_liq(c,j)+h2osoi_ice(c,j))* &
                          fracl
                  end if
               end if

               h2osoi_liq_tot(c) = h2osoi_liq_tot(c) + h2osoi_liq(c,j)
               h2osoi_ice_tot(c) = h2osoi_ice_tot(c) + h2osoi_ice(c,j)

            end if
         end do
      end do

      ! TODO - if this block of code is moved out of here - the SoilHydrology 
      ! will NOT effect t_grnd, t_grnd_u or t_grnd_r

      do fc = 1, num_nolakec

         c = filter_nolakec(fc)
         l = col%landunit(c)

         ! t_grnd is weighted average of exposed soil and snow
         if (snl(c) < 0) then
            t_grnd(c) = frac_sno_eff(c) * t_soisno(c,snl(c)+1) &
                 + (1 - frac_sno_eff(c)- frac_h2osfc(c)) * t_soisno(c,1) &
                 + frac_h2osfc(c) * t_h2osfc(c)
         else
            t_grnd(c) = (1 - frac_h2osfc(c)) * t_soisno(c,1) + frac_h2osfc(c) * t_h2osfc(c)
         endif

         if (lun%urbpoi(l)) then
            t_grnd_u(c) = t_soisno(c,snl(c)+1)
         else
            t_soi_10cm(c) = t_soi_10cm(c)/0.1_r8
            tsoi17(c) =  tsoi17(c)/0.17_r8         ! F. Li and S. Levis
         end if
         if (lun%itype(l)==istsoil .or. lun%itype(l)==istcrop) then
            t_grnd_r(c) = t_soisno(c,snl(c)+1)
         end if

      end do

      do j = 1, nlevgrnd
         do fc = 1, num_nolakec
            c = filter_nolakec(fc)
            if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
                 .or. ctype(c) == icol_roof) .and. j > nlevurb) then
            else
               h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)
            end if
         end do
      end do

!      if (use_cn) then
         ! Update soilpsi.
         ! ZMS: Note this could be merged with the following loop updating smp_l in the future.
         do j = 1, nlevgrnd
            do fc = 1, num_hydrologyc
               c = filter_hydrologyc(fc)

               if (h2osoi_liq(c,j) > 0._r8) then

                  vwc = h2osoi_liq(c,j)/(dz(c,j)*denh2o)

                  ! the following limit set to catch very small values of 
                  ! fractional saturation that can crash the calculation of psi

                  ! use the same contants used in the supercool so that psi for frozen soils is consistent
                  fsattmp = max(vwc/watsat(c,j), 0.001_r8)
                  psi = sucsat(c,j) * (-9.8e-6_r8) * (fsattmp)**(-bsw(c,j))  ! Mpa
                  soilpsi(c,j) = min(max(psi,-15.0_r8),0._r8)

               else 
                  soilpsi(c,j) = -15.0_r8
               end if
            end do
         end do
!      end if

      ! Update smp_l for history and for ch4Mod.
      ! ZMS: Note, this form, which seems to be the same as used in SoilWater, DOES NOT distinguish between
      ! ice and water volume, in contrast to the soilpsi calculation above. It won't be used in ch4Mod if
      ! t_soisno <= tfrz, though.
      do j = 1, nlevgrnd
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)

            s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
            s_node = min(1.0_r8, s_node)

            smp_l(c,j) = -sucsat(c,j)*s_node**(-bsw(c,j))
            smp_l(c,j) = max(smpmin(c), smp_l(c,j))
         end do
      end do

 !     if (use_cn) then
         ! Available soil water up to a depth of 0.05 m.
         ! Potentially available soil water (=whc) up to a depth of 0.05 m.
         ! Water content as fraction of whc up to a depth of 0.05 m.

         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            rwat(c) = 0._r8
            swat(c) = 0._r8
            rz(c)   = 0._r8
         end do

         do j = 1, nlevgrnd
            do fc = 1, num_hydrologyc
               c = filter_hydrologyc(fc)
               !if (z(c,j)+0.5_r8*dz(c,j) <= 0.5_r8) then
               if (z(c,j)+0.5_r8*dz(c,j) <= 0.05_r8) then
                  watdry = watsat(c,j) * (316230._r8/sucsat(c,j)) ** (-1._r8/bsw(c,j))
                  rwat(c) = rwat(c) + (h2osoi_vol(c,j)-watdry) * dz(c,j)
                  swat(c) = swat(c) + (watsat(c,j)    -watdry) * dz(c,j)
                  rz(c) = rz(c) + dz(c,j)
               end if
            end do
         end do

         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            if (rz(c) /= 0._r8) then
               tsw  = rwat(c)/rz(c)
               stsw = swat(c)/rz(c)
            else
               watdry = watsat(c,1) * (316230._r8/sucsat(c,1)) ** (-1._r8/bsw(c,1))
               tsw = h2osoi_vol(c,1) - watdry
               stsw = watsat(c,1) - watdry
            end if
            wf(c) = tsw/stsw
         end do

         do j = 1, nlevgrnd
            do fc = 1, num_hydrologyc
               c = filter_hydrologyc(fc)
               if (z(c,j)+0.5_r8*dz(c,j) <= 0.17_r8) then
                  watdry = watsat(c,j) * (316230._r8/sucsat(c,j)) ** (-1._r8/bsw(c,j))
                  rwat(c) = rwat(c) + (h2osoi_vol(c,j)-watdry) * dz(c,j)
                  swat(c) = swat(c) + (watsat(c,j)    -watdry) * dz(c,j)
                  rz(c) = rz(c) + dz(c,j)
               end if
            end do
         end do

         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            if (rz(c) /= 0._r8) then
               tsw  = rwat(c)/rz(c)
               stsw = swat(c)/rz(c)
            else
               watdry = watsat(c,1) * (316230._r8/sucsat(c,1)) ** (-1._r8/bsw(c,1))
               tsw = h2osoi_vol(c,1) - watdry
               stsw = watsat(c,1) - watdry
            end if
            wf2(c) = tsw/stsw
         end do
  !    end if

      ! top-layer diagnostics
      do fc = 1, num_snowc
         c = filter_snowc(fc)
         h2osno_top(c)  = h2osoi_ice(c,snl(c)+1) + h2osoi_liq(c,snl(c)+1)
      enddo

      ! Zero variables in columns without snow
      do fc = 1, num_nosnowc
         c = filter_nosnowc(fc)
            
         h2osno_top(c)      = 0._r8
         snw_rds(c,:)       = 0._r8

         ! top-layer diagnostics (spval is not averaged when computing history fields)
         snot_top(c)        = spval
         dTdz_top(c)        = spval
         snw_rds_top(c)     = spval
         sno_liq_top(c)     = spval
      end do

    end associate

  end subroutine HydrologyNoDrainage

end Module HydrologyNoDrainageMod
