module CNDriverMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Ecosystem dynamics: phenology, vegetation
  !
  ! !USES:
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use clm_varctl                      , only : use_c13, use_c14, use_fates, use_fates_bgc
  use dynSubgridControlMod            , only : get_do_harvest, get_do_grossunrep
  use decompMod                       , only : bounds_type
  use perf_mod                        , only : t_startf, t_stopf
  use clm_varctl                      , only : use_nitrif_denitrif, use_nguardrail
  use clm_varctl                      , only : use_crop, use_crop_agsys, use_cn
  use SoilBiogeochemDecompCascadeConType, only : mimics_decomp, century_decomp, decomp_method
  use CNSharedParamsMod               , only : use_fun
  use CNVegStateType                  , only : cnveg_state_type
  use CNVegCarbonStateType            , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType             , only : cnveg_carbonflux_type
  use CNVegNitrogenStateType          , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType           , only : cnveg_nitrogenflux_type
  use CNProductsMod                   , only : cn_products_type
  use SoilBiogeochemStateType         , only : soilbiogeochem_state_type
  use SoilBiogeochemCarbonStateType   , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemCarbonFluxType    , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
  use CNDVType                        , only : dgvs_type
  use CanopyStateType                 , only : canopystate_type
  use SoilStateType                   , only : soilstate_type
  use TemperatureType                 , only : temperature_type
  use WaterStateBulkType                  , only : waterstatebulk_type
  use WaterDiagnosticBulkType                  , only : waterdiagnosticbulk_type
  use WaterFluxBulkType                   , only : waterfluxbulk_type
  use Wateratm2lndBulkType                   , only : wateratm2lndbulk_type
  use atm2lndType                     , only : atm2lnd_type
  use SoilStateType                   , only : soilstate_type
  use TemperatureType                 , only : temperature_type 
  use PhotosynthesisMod               , only : photosyns_type
  use ch4Mod                          , only : ch4_type
  use EnergyFluxType                  , only : energyflux_type
  use SaturatedExcessRunoffMod        , only : saturated_excess_runoff_type
  use ActiveLayerMod                  , only : active_layer_type
  use SoilWaterRetentionCurveMod      , only : soil_water_retention_curve_type
  use CLMFatesInterfaceMod            , only : hlm_fates_interface_type
  use CropReprPoolsMod                , only : nrepr
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNDriverInit         ! Ecosystem dynamics: initialization
  public :: CNDriverNoLeaching   ! Ecosystem dynamics: phenology, vegetation, before doing N leaching
  public :: CNDriverLeaching     ! Ecosystem dynamics: phenology, vegetation, doing N leaching
  public :: CNDriverSummarizeStates
  public :: CNDriverSummarizeFluxes
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNDriverInit(bounds, NLFilename, cnfire_method)
    !
    ! !DESCRIPTION:
    ! Initialzation of the CN Ecosystem dynamics.
    !
    ! !USES:
    use CNSharedParamsMod           , only : use_fun
    use CNPhenologyMod              , only : CNPhenologyInit
    use FireMethodType              , only : fire_method_type
    use SoilBiogeochemCompetitionMod, only : SoilBiogeochemCompetitionInit
    !
    ! !ARGUMENTS:
    type(bounds_type)                      , intent(in)    :: bounds      
    character(len=*)                       , intent(in)    :: NLFilename     ! Namelist filename
    class(fire_method_type)                , intent(inout) :: cnfire_method 
    !-----------------------------------------------------------------------
    call SoilBiogeochemCompetitionInit(bounds)
    if(use_cn)then
       call CNPhenologyInit(bounds)
       call cnfire_method%FireInit(bounds, NLFilename)
    end if
  end subroutine CNDriverInit

  !-----------------------------------------------------------------------
  subroutine CNDriverNoLeaching(bounds,                                                    &
       num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp,                     &
       num_pcropp, filter_pcropp, num_soilnopcropp, filter_soilnopcropp,                   &
       num_actfirec, filter_actfirec, num_actfirep, filter_actfirep,                       &
       num_exposedvegp, filter_exposedvegp, num_noexposedvegp, filter_noexposedvegp,       &
       cnveg_state_inst,                                                                   &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst,                                      &
       c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst,                              &
       c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst,                              &
       cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst,                                  &
       c_products_inst, c13_products_inst, c14_products_inst, n_products_inst,             &
       soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst,                    &
       c13_soilbiogeochem_carbonflux_inst, c13_soilbiogeochem_carbonstate_inst,            &
       c14_soilbiogeochem_carbonflux_inst, c14_soilbiogeochem_carbonstate_inst,            &
       soilbiogeochem_state_inst,                                                          &
       soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst,                &
       active_layer_inst, clm_fates,                                                       &
       atm2lnd_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, waterfluxbulk_inst,    &
       wateratm2lndbulk_inst, canopystate_inst, soilstate_inst, temperature_inst,          &
       soil_water_retention_curve, crop_inst, ch4_inst,                                    &
       dgvs_inst, photosyns_inst, saturated_excess_runoff_inst, energyflux_inst,           &
       nutrient_competition_method, cnfire_method, dribble_crophrv_xsmrpool_2atm)
    !
    ! !DESCRIPTION:
    ! The core CN code is executed here. Calculates fluxes for maintenance
    ! respiration, decomposition, allocation, phenology, and growth respiration.
    ! These routines happen on the radiation time step so that canopy structure
    ! stays synchronized with albedo calculations.
    !
    ! !USES:
    use clm_varpar                        , only: nlevdecomp_full, nvegcpool, nvegnpool 
    use clm_varpar                        , only: nlevdecomp, ndecomp_cascade_transitions, ndecomp_pools
    use subgridAveMod                     , only: p2c
    use CropType                          , only: crop_type
    use CNAllocationMod                   , only: calc_gpp_mr_availc, calc_crop_allocation_fractions
    use CNAllocationMod                   , only: calc_allometry
    use CNNDynamicsMod                    , only: CNNDeposition,CNNFixation, CNNFert, CNSoyfix,CNFreeLivingFixation
    use CNMRespMod                        , only: CNMResp
    use CNFUNMod                          , only: CNFUNInit  !, CNFUN 
    use CNPhenologyMod                    , only: CNPhenology
    use CNGRespMod                        , only: CNGResp
    use FireMethodType                    , only: fire_method_type
    use CNCIsoFluxMod                     , only: CIsoFlux1, CIsoFlux2, CIsoFlux2h, CIsoFlux2g, CIsoFlux3
    use CNC14DecayMod                     , only: C14Decay
    use CNCStateUpdate1Mod                , only: CStateUpdate1,CStateUpdate0
    use CNCStateUpdate2Mod                , only: CStateUpdate2, CStateUpdate2h, CStateUpdate2g
    use CNCStateUpdate3Mod                , only: CStateUpdate3
    use CNNStateUpdate1Mod                , only: NStateUpdate1
    use CNNStateUpdate2Mod                , only: NStateUpdate2, NStateUpdate2h, NStateUpdate2g
    use CNGapMortalityMod                 , only: CNGapMortality
    use CNSharedParamsMod                 , only: use_fun
    use dynHarvestMod                     , only: CNHarvest
    use dynGrossUnrepMod                  , only: CNGrossUnrep
    use SoilBiogeochemDecompCascadeMIMICSMod, only: decomp_rates_mimics
    use SoilBiogeochemDecompCascadeBGCMod , only: decomp_rate_constants_bgc
    use SoilBiogeochemCompetitionMod      , only: SoilBiogeochemCompetition
    use SoilBiogeochemDecompMod           , only: SoilBiogeochemDecomp
    use SoilBiogeochemLittVertTranspMod   , only: SoilBiogeochemLittVertTransp
    use SoilBiogeochemPotentialMod        , only: SoilBiogeochemPotential 
    use SoilBiogeochemVerticalProfileMod  , only: SoilBiogeochemVerticalProfile
    use SoilBiogeochemNitrifDenitrifMod   , only: SoilBiogeochemNitrifDenitrif
    use SoilBiogeochemNStateUpdate1Mod    , only: SoilBiogeochemNStateUpdate1
    use NutrientCompetitionMethodMod      , only: nutrient_competition_method_type
    use CNPrecisionControlMod             , only: CNPrecisionControl
    !
    ! !ARGUMENTS:
    type(bounds_type)                       , intent(in)    :: bounds  
    integer                                 , intent(in)    :: num_bgc_soilc        ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_bgc_soilc(:)  ! filter for soil columns
    integer                                 , intent(in)    :: num_bgc_vegp         ! number of veg patches in filter
    integer                                 , intent(in)    :: filter_bgc_vegp(:)   ! filter for veg patches
    integer                                 , intent(out)   :: num_actfirep         ! number of soil patches on fire in filter
    integer                                 , intent(out)   :: filter_actfirep(:)   ! filter for soil patches on fire
    integer                                 , intent(out)   :: num_actfirec         ! number of soil columns on fire in filter
    integer                                 , intent(out)   :: filter_actfirec(:)   ! filter for soil columns on fire
    integer                                 , intent(in)    :: num_pcropp        ! number of prog. crop patches in filter
    integer                                 , intent(in)    :: filter_pcropp(:)  ! filter for prognostic crop patches
    integer                                 , intent(in)    :: num_soilnopcropp       ! number of non-prog. crop soil patches in filter
    integer                                 , intent(in)    :: filter_soilnopcropp(:) ! filter for non-prog. crop soil patches
    integer                                 , intent(in)    :: num_exposedvegp        ! number of points in filter_exposedvegp
    integer                                 , intent(in)    :: filter_exposedvegp(:)  ! patch filter for non-snow-covered veg
    integer                                 , intent(in)    :: num_noexposedvegp       ! number of points in filter_noexposedvegp
    integer                                 , intent(in)    :: filter_noexposedvegp(:) ! patch filter where frac_veg_nosno is 0 
    type(cnveg_state_type)                  , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonflux_type)             , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)            , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)             , intent(inout) :: c13_cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)            , intent(inout) :: c13_cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)             , intent(inout) :: c14_cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)            , intent(inout) :: c14_cnveg_carbonstate_inst
    type(cnveg_nitrogenflux_type)           , intent(inout) :: cnveg_nitrogenflux_inst
    type(cnveg_nitrogenstate_type)          , intent(inout) :: cnveg_nitrogenstate_inst
    type(cn_products_type)                  , intent(inout) :: c_products_inst
    type(cn_products_type)                  , intent(inout) :: c13_products_inst
    type(cn_products_type)                  , intent(inout) :: c14_products_inst
    type(cn_products_type)                  , intent(inout) :: n_products_inst
    type(soilbiogeochem_state_type)         , intent(inout) :: soilbiogeochem_state_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: c13_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c13_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: c14_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c14_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(active_layer_type)                 , intent(in)    :: active_layer_inst
    type(atm2lnd_type)                      , intent(in)    :: atm2lnd_inst
    type(waterstatebulk_type)                   , intent(in)    :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)                   , intent(in)    :: waterdiagnosticbulk_inst
    type(waterfluxbulk_type)                    , intent(inout)    :: waterfluxbulk_inst
    type(wateratm2lndbulk_type)                    , intent(inout)    :: wateratm2lndbulk_inst
    type(canopystate_type)                  , intent(inout)    :: canopystate_inst
    type(soilstate_type)                    , intent(inout) :: soilstate_inst
    type(temperature_type)                  , intent(inout) :: temperature_inst
    class(soil_water_retention_curve_type)  , intent(in)    :: soil_water_retention_curve
    type(crop_type)                         , intent(inout) :: crop_inst
    type(ch4_type)                          , intent(in)    :: ch4_inst
    type(dgvs_type)                         , intent(inout) :: dgvs_inst
    type(photosyns_type)                    , intent(in)    :: photosyns_inst
    type(saturated_excess_runoff_type)      , intent(in)    :: saturated_excess_runoff_inst
    type(energyflux_type)                   , intent(in)    :: energyflux_inst
    class(nutrient_competition_method_type) , intent(inout) :: nutrient_competition_method
    class(fire_method_type)                 , intent(inout) :: cnfire_method
    logical                                 , intent(in)    :: dribble_crophrv_xsmrpool_2atm
    type(hlm_fates_interface_type)          , intent(inout) :: clm_fates
    !
    ! !LOCAL VARIABLES:
    real(r8):: cn_decomp_pools(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
    real(r8):: p_decomp_cpool_loss(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_cascade_transitions) !potential C loss from one pool to another
    real(r8):: pmnf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_cascade_transitions) !potential mineral N flux, from one pool to another
    real(r8):: p_decomp_npool_to_din(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_cascade_transitions)  ! potential flux to dissolved inorganic N
    real(r8):: p_decomp_cn_gain(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)  ! C:N ratio of the flux gained by the receiver pool
    integer :: begp,endp
    integer :: begc,endc

    integer :: dummy_to_make_pgi_happy
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    associate(                                                                    &
         laisun                    => canopystate_inst%laisun_patch             , & ! Input:  [real(r8) (:)   ]  sunlit projected leaf area index        
         laisha                    => canopystate_inst%laisha_patch             , & ! Input:  [real(r8) (:)   ]  shaded projected leaf area index        
         frac_veg_nosno            => canopystate_inst%frac_veg_nosno_patch     , & ! Input:  [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]
         frac_veg_nosno_alb        => canopystate_inst%frac_veg_nosno_alb_patch , & ! Output: [integer  (:) ] frac of vegetation not covered by snow [-]         
         tlai                      => canopystate_inst%tlai_patch               , & ! Input:  [real(r8) (:) ]  one-sided leaf area index, no burying by snow     
         tsai                      => canopystate_inst%tsai_patch               , & ! Input:  [real(r8) (:)   ]  one-sided stem area index, no burying by snow     
         elai                      => canopystate_inst%elai_patch               , & ! Output: [real(r8) (:) ] one-sided leaf area index with burying by snow    
         esai                      => canopystate_inst%esai_patch               , & ! Output: [real(r8) (:) ] one-sided stem area index with burying by snow    
         htop                      => canopystate_inst%htop_patch               , & ! Output: [real(r8) (:) ] canopy top (m)                                     
         hbot                      => canopystate_inst%hbot_patch                 & ! Output: [real(r8) (:) ] canopy bottom (m)                                  
      )

    ! --------------------------------------------------
    ! zero the column-level C and N fluxes
    ! --------------------------------------------------

    call t_startf('CNZero')
    ! COMPILER_BUG(wjs, 2014-11-29, pgi 14.7) Without this, the filter is full of garbage
    ! in some situations
    call t_startf('CNZero-soilbgc-cflux')
    dummy_to_make_pgi_happy = ubound(filter_bgc_soilc, 1)
    call soilbiogeochem_carbonflux_inst%SetValues( &
         num_bgc_soilc, filter_bgc_soilc, 0._r8)
    if ( use_c13 ) then
       call c13_soilbiogeochem_carbonflux_inst%SetValues( &
            num_bgc_soilc, filter_bgc_soilc, 0._r8)
    end if
    if ( use_c14 ) then
       call c14_soilbiogeochem_carbonflux_inst%SetValues( &
            num_bgc_soilc, filter_bgc_soilc, 0._r8)
    end if
    call t_stopf('CNZero-soilbgc-cflux')

    if(num_bgc_vegp>0)then
       call t_startf('CNZero-vegbgc-cflux')
       call cnveg_carbonflux_inst%SetValues( &
            nvegcpool,&
            num_bgc_vegp, filter_bgc_vegp, 0._r8, &
            num_bgc_soilc, filter_bgc_soilc, 0._r8)
       if ( use_c13 ) then
          call c13_cnveg_carbonflux_inst%SetValues( &
               nvegcpool,&
               num_bgc_vegp, filter_bgc_vegp, 0._r8, &
               num_bgc_soilc, filter_bgc_soilc, 0._r8)
       end if
       if ( use_c14 ) then
          call c14_cnveg_carbonflux_inst%SetValues( &
               nvegcpool,&
               num_bgc_vegp, filter_bgc_vegp, 0._r8, &
               num_bgc_soilc, filter_bgc_soilc, 0._r8)
       end if
       call t_stopf('CNZero-vegbgc-cflux')
       
       call t_startf('CNZero-vegbgc-nflux')
       call cnveg_nitrogenflux_inst%SetValues( &
            nvegnpool, &
            num_bgc_vegp, filter_bgc_vegp, 0._r8, &
            num_bgc_soilc, filter_bgc_soilc, 0._r8)
       call t_stopf('CNZero-vegbgc-nflux')
    end if
       
    call t_startf('CNZero-soilbgc-nflux')
    call soilbiogeochem_nitrogenflux_inst%SetValues( &
         num_bgc_soilc, filter_bgc_soilc, 0._r8)

    call t_stopf('CNZero-soilbgc-nflux')
    call t_stopf('CNZero')
    
    ! --------------------------------------------------
    ! Nitrogen Deposition, Fixation and Respiration
    ! --------------------------------------------------

    call t_startf('CNDeposition')
    call CNNDeposition(bounds, &
         atm2lnd_inst, soilbiogeochem_nitrogenflux_inst)
    call t_stopf('CNDeposition')
    
    if(use_fun)then
        call t_startf('CNFLivFixation')
        call CNFreeLivingFixation( num_bgc_soilc, filter_bgc_soilc, &
             waterfluxbulk_inst, soilbiogeochem_nitrogenflux_inst)
        call t_stopf('CNFLivFixation')
    else
       call t_startf('CNFixation')
       call CNNFixation( num_bgc_soilc, filter_bgc_soilc, &
            cnveg_carbonflux_inst, soilbiogeochem_nitrogenflux_inst, &
            clm_fates, bounds%clump_index)
       call t_stopf('CNFixation')
    end if
  

    if (use_crop) then
       call CNNFert(bounds, num_bgc_soilc,filter_bgc_soilc, &
            cnveg_nitrogenflux_inst, soilbiogeochem_nitrogenflux_inst)

       if (.not. use_fun) then  ! if FUN is active, then soy fixation handled by FUN
          call  CNSoyfix (bounds, num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
               waterdiagnosticbulk_inst, crop_inst, cnveg_state_inst, cnveg_nitrogenflux_inst , &
               soilbiogeochem_state_inst, soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)
       end if
    end if

    call t_startf('CNMResp')
    call CNMResp(bounds, num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
         canopystate_inst, soilstate_inst, temperature_inst, photosyns_inst, &
         cnveg_carbonflux_inst, cnveg_nitrogenstate_inst)
    call t_stopf('CNMResp')

    !--------------------------------------------
    ! Soil Biogeochemistry
    !--------------------------------------------

    call t_startf('SoilBiogeochem')
    call t_startf('DecompRate')
    if (decomp_method == century_decomp) then
       call decomp_rate_constants_bgc(bounds, num_bgc_soilc, filter_bgc_soilc, &
            soilstate_inst, temperature_inst, ch4_inst, soilbiogeochem_carbonflux_inst, &
            cnveg_state_inst%idop_patch)
    else if (decomp_method == mimics_decomp) then
       call decomp_rates_mimics(bounds, num_bgc_soilc, filter_bgc_soilc, &
            num_bgc_vegp, filter_bgc_vegp, clm_fates, &
            soilstate_inst, temperature_inst, cnveg_carbonflux_inst, ch4_inst, &
            soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst, &
            cnveg_state_inst%idop_patch)
    end if
    call t_stopf('DecompRate')

    call t_startf('SoilBiogeochemPotential')
    ! calculate potential decomp rates and total immobilization demand (previously inlined in CNDecompAlloc)
    call SoilBiogeochemPotential (bounds, num_bgc_soilc, filter_bgc_soilc,                                                    &
         soilbiogeochem_state_inst, soilbiogeochem_carbonstate_inst, soilbiogeochem_carbonflux_inst,                  &
         soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst,                                         &
         cn_decomp_pools=cn_decomp_pools(begc:endc,1:nlevdecomp,1:ndecomp_pools), & 
         p_decomp_cpool_loss=p_decomp_cpool_loss(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions), &
         p_decomp_cn_gain=p_decomp_cn_gain(begc:endc,1:nlevdecomp,1:ndecomp_pools), &
         pmnf_decomp_cascade=pmnf_decomp_cascade(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions), &
         p_decomp_npool_to_din=p_decomp_npool_to_din(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions))
    call t_stopf('SoilBiogeochemPotential')

    ! calculate vertical profiles for distributing soil and litter C and N
    ! (previously subroutine decomp_vertprofiles called from CNDecompAlloc)
    call SoilBiogeochemVerticalProfile(bounds, num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
         active_layer_inst, soilstate_inst,soilbiogeochem_state_inst)

    ! calculate nitrification and denitrification rates (previously subroutine nitrif_denitrif called from CNDecompAlloc)
    if (use_nitrif_denitrif) then 
       call SoilBiogeochemNitrifDenitrif(bounds, num_bgc_soilc, filter_bgc_soilc, &
            soilstate_inst, waterstatebulk_inst, temperature_inst, ch4_inst, &
            soilbiogeochem_carbonflux_inst, soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)
    end if
    call t_stopf('SoilBiogeochem')

    !--------------------------------------------
    ! Resolve the competition between plants and soil heterotrophs 
    ! for available soil mineral N resource 
    !--------------------------------------------

    call t_startf('CNDecompAlloc')

    ! Jinyun Tang: at this stage, the plant_nutrient_demand only calculates the plant ntirgeon demand.
    ! Assume phosphorus dynamics will be included in the future. Also, I consider plant_nutrient_demand
    ! as a generic interface to call actual nutrient calculation from different aboveground plantbgc. 
    ! Right now it is assumed the plant nutrient demand is summarized into columnwise demand, and the 
    ! nutrient redistribution after uptake is done by the plant bgc accordingly. 
    ! When nutrient competition is required to be done at cohort level both plant_nutrient_demand and 
    ! do_nutrient_competition should be modified, but that modification should not significantly change 
    ! the current interface.

    !RF: moved ths call to before nutrient_demand, so that croplive didn't change half way through crop N cycle. 
    if(num_bgc_vegp>0)then
       if ( use_fun) then
          call t_startf('CNPhenology_phase1')
          call CNPhenology (bounds, num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, &
               filter_bgc_vegp, num_pcropp, filter_pcropp, &
               waterdiagnosticbulk_inst, wateratm2lndbulk_inst, temperature_inst, atm2lnd_inst, &
               crop_inst, canopystate_inst, soilstate_inst, dgvs_inst, &
               cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst, &
               cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
               c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst, &
               leaf_prof_patch=soilbiogeochem_state_inst%leaf_prof_patch(begp:endp,1:nlevdecomp_full), &
               froot_prof_patch=soilbiogeochem_state_inst%froot_prof_patch(begp:endp,1:nlevdecomp_full), &
               phase=1)
          call t_stopf('CNPhenology_phase1')
          
          call t_startf('CNFUNInit')
          call CNFUNInit(bounds,cnveg_state_inst,cnveg_carbonstate_inst,cnveg_nitrogenstate_inst)
          call t_stopf('CNFUNInit')
          
       end if
       
       call t_startf('cnalloc')
       call calc_gpp_mr_availc( &
            bounds, num_bgc_vegp, filter_bgc_vegp, &
            crop_inst, photosyns_inst, canopystate_inst, &
            cnveg_carbonstate_inst, cnveg_carbonflux_inst, &
            c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst)
       
       if (.not. use_crop_agsys) then
          call calc_crop_allocation_fractions(bounds, num_pcropp, filter_pcropp, &
               crop_inst, cnveg_state_inst)
       end if
       
       call calc_allometry(num_bgc_vegp, filter_bgc_vegp, &
            cnveg_carbonflux_inst, cnveg_state_inst)
       call t_stopf('cnalloc')
    end if
    
     call t_startf('calc_plant_nutrient_demand')
     ! We always call calc_plant_nutrient_demand for natural veg patches, but only call
     ! it for crop patches if NOT running with AgSys (since AgSys calculates the relevant
     ! output variables in its own way).
     call nutrient_competition_method%calc_plant_nutrient_demand ( &
          bounds,                                                          &
          num_soilnopcropp, filter_soilnopcropp, .false.,                  &
          crop_inst, canopystate_inst,                                     &
          cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst, &
          cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst,               &
          soilbiogeochem_carbonflux_inst, soilbiogeochem_nitrogenstate_inst, &
          energyflux_inst)
     if (.not. use_crop_agsys) then
        call nutrient_competition_method%calc_plant_nutrient_demand ( &
             bounds,                                                          &
             num_pcropp, filter_pcropp, .true.,                               &
             crop_inst, canopystate_inst,                                     &
             cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst, &
             cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst,               &
             soilbiogeochem_carbonflux_inst, soilbiogeochem_nitrogenstate_inst, &
             energyflux_inst)
     end if

     ! get the column-averaged plant_ndemand (needed for following call to SoilBiogeochemCompetition)

     if(num_bgc_vegp>0)then
        call p2c(bounds, num_bgc_soilc, filter_bgc_soilc,                    &
             cnveg_nitrogenflux_inst%plant_ndemand_patch(begp:endp), &
             soilbiogeochem_state_inst%plant_ndemand_col(begc:endc))
     else
        ! With FATES N coupling, we will have a call to fill
        ! this in on the filter_bgc_soilc
        soilbiogeochem_state_inst%plant_ndemand_col(begc:endc) = 0._r8
     end if

     call t_stopf('calc_plant_nutrient_demand')

     ! resolve plant/heterotroph competition for mineral N 
 
   
     call t_startf('soilbiogeochemcompetition')
     call SoilBiogeochemCompetition (bounds, num_bgc_soilc, filter_bgc_soilc,num_bgc_vegp, filter_bgc_vegp, &
                                     p_decomp_cn_gain, pmnf_decomp_cascade, waterstatebulk_inst, &
                                     waterfluxbulk_inst,temperature_inst,soilstate_inst,cnveg_state_inst,          &
                                     cnveg_carbonstate_inst               ,&
                                     cnveg_carbonflux_inst,cnveg_nitrogenstate_inst,cnveg_nitrogenflux_inst,   &
                                     soilbiogeochem_carbonflux_inst,&
                                     soilbiogeochem_state_inst,soilbiogeochem_nitrogenstate_inst,              &
                                     soilbiogeochem_nitrogenflux_inst,canopystate_inst)
     call t_stopf('soilbiogeochemcompetition')

    ! distribute the available N between the competing patches  on the basis of 
    ! relative demand, and allocate C and N to new growth and storage

    call t_startf('calc_plant_nutrient_competition')
    call nutrient_competition_method%calc_plant_nutrient_competition ( &
         bounds, num_bgc_vegp, filter_bgc_vegp, &
         cnveg_state_inst, crop_inst, canopystate_inst, &
         cnveg_carbonstate_inst, cnveg_carbonflux_inst, &
         c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst, &
         cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
         soilbiogeochem_nitrogenstate_inst, &
         fpg_col=soilbiogeochem_state_inst%fpg_col(begc:endc))
    call t_stopf('calc_plant_nutrient_competition')

    call t_stopf('CNDecompAlloc')

    !--------------------------------------------
    ! Calculate litter and soil decomposition rate
    !--------------------------------------------

    ! Calculation of actual immobilization and decomp rates, following
    ! resolution of plant/heterotroph  competition for mineral N (previously inlined in CNDecompAllocation in CNDecompMod)

    call t_startf('SoilBiogeochemDecomp')

    call SoilBiogeochemDecomp (bounds, num_bgc_soilc, filter_bgc_soilc,                                                       &
         soilbiogeochem_state_inst, soilbiogeochem_carbonstate_inst, soilbiogeochem_carbonflux_inst,                  &
         soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst,                                         &
         cn_decomp_pools=cn_decomp_pools(begc:endc,1:nlevdecomp,1:ndecomp_pools),                       & 
         p_decomp_cpool_loss=p_decomp_cpool_loss(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions), &
         pmnf_decomp_cascade=pmnf_decomp_cascade(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions), &
         p_decomp_npool_to_din=p_decomp_npool_to_din(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions))

    call t_stopf('SoilBiogeochemDecomp')

    !--------------------------------------------
    ! Phenology
    !--------------------------------------------

    ! CNphenology needs to be called after above calls, since it depends on current
    ! time-step fluxes to new growth on the lastlitterfall timestep in deciduous systems
    if(num_bgc_vegp>0)then
       call t_startf('CNPhenology')
       if ( .not. use_fun ) then
          call CNPhenology (bounds, num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, &
               filter_bgc_vegp, num_pcropp, filter_pcropp, &
               waterdiagnosticbulk_inst, wateratm2lndbulk_inst, temperature_inst, atm2lnd_inst, &
               crop_inst, canopystate_inst, soilstate_inst, dgvs_inst, &
               cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst, &
               cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
               c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst, &
               leaf_prof_patch=soilbiogeochem_state_inst%leaf_prof_patch(begp:endp,1:nlevdecomp_full), &
               froot_prof_patch=soilbiogeochem_state_inst%froot_prof_patch(begp:endp,1:nlevdecomp_full), &
               phase=1)
       end if
       call CNPhenology (bounds, num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, &
            filter_bgc_vegp, num_pcropp, filter_pcropp, &
            waterdiagnosticbulk_inst, wateratm2lndbulk_inst, temperature_inst, atm2lnd_inst, &
            crop_inst, canopystate_inst, soilstate_inst, dgvs_inst, &
            cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst, &
            cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
            c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst, &
            leaf_prof_patch=soilbiogeochem_state_inst%leaf_prof_patch(begp:endp,1:nlevdecomp_full), &
            froot_prof_patch=soilbiogeochem_state_inst%froot_prof_patch(begp:endp,1:nlevdecomp_full), &
            phase=2)
       
       call t_stopf('CNPhenology')
    end if
    !--------------------------------------------
    ! Growth respiration
    !--------------------------------------------

    call t_startf('CNGResp')

    call CNGResp(num_bgc_vegp, filter_bgc_vegp,&
         cnveg_carbonflux_inst, canopystate_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst)  
         
    call t_stopf('CNGResp')

    !--------------------------------------------------------------------------
    ! CNUpdate0
    ! The state updates are still called for the matrix solution (use_matrixn
    ! and use_soil_matrixcn) but most of the state updates are done after
    ! the matrix multiply in VegMatrix and SoilMatrix.
    !--------------------------------------------------------------------------

    call t_startf('CNUpdate0')

    call CStateUpdate0(num_bgc_vegp, filter_bgc_vegp, &
         cnveg_carbonflux_inst, cnveg_carbonstate_inst)

    if ( use_c13 ) then
       call CStateUpdate0(num_bgc_vegp, filter_bgc_vegp, &
            c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst)
    end if

    if ( use_c14 ) then
       call CStateUpdate0(num_bgc_vegp, filter_bgc_vegp, &
            c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst)
    end if

    call t_stopf('CNUpdate0')

    if ( use_nguardrail .and. num_bgc_vegp>0 ) then
       call t_startf('CNPrecisionControl')
       call CNPrecisionControl(bounds, num_bgc_vegp, filter_bgc_vegp, &
            cnveg_carbonstate_inst, c13_cnveg_carbonstate_inst, &
            c14_cnveg_carbonstate_inst, cnveg_nitrogenstate_inst)
       call t_stopf('CNPrecisionControl')
    end if
    !--------------------------------------------
    ! Update1
    !--------------------------------------------

    call t_startf('CNUpdate1')

    ! Set the carbon isotopic flux variables (except for gap-phase mortality and fire fluxes)
    if ( use_c13 ) then

       call CIsoFlux1(num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp,              &
            soilbiogeochem_state_inst,                                               &
            soilbiogeochem_carbonflux_inst,  soilbiogeochem_carbonstate_inst,        &
            cnveg_carbonflux_inst, cnveg_carbonstate_inst,                           &
            c13_soilbiogeochem_carbonflux_inst, c13_soilbiogeochem_carbonstate_inst, &
            c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst,                   &
            isotope='c13')
    end if
    if ( use_c14 ) then
       call CIsoFlux1(num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp,              &
            soilbiogeochem_state_inst,                                               &
            soilbiogeochem_carbonflux_inst,  soilbiogeochem_carbonstate_inst,        &
            cnveg_carbonflux_inst, cnveg_carbonstate_inst,                           &
            c14_soilbiogeochem_carbonflux_inst, c14_soilbiogeochem_carbonstate_inst, &
            c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst,                   &
            isotope='c14')
    end if

    ! Update all prognostic carbon state variables (except for gap-phase mortality and fire fluxes)
    call CStateUpdate1( num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
         crop_inst, cnveg_carbonflux_inst, cnveg_carbonstate_inst, &
         soilbiogeochem_carbonflux_inst, dribble_crophrv_xsmrpool_2atm, &
         clm_fates, bounds%clump_index)
    if ( use_c13 ) then
       call CStateUpdate1(num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
            crop_inst, c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst, &
            c13_soilbiogeochem_carbonflux_inst, dribble_crophrv_xsmrpool_2atm, &
            clm_fates, bounds%clump_index)
    end if
    if ( use_c14 ) then
       call CStateUpdate1(num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
            crop_inst, c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst, &
            c14_soilbiogeochem_carbonflux_inst, dribble_crophrv_xsmrpool_2atm, &
            clm_fates, bounds%clump_index)
    end if

    ! Update all prognostic nitrogen state variables (except for gap-phase mortality and fire fluxes)
    call NStateUpdate1(num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
         cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst, &
         clm_fates, bounds%clump_index)

    call t_stopf('CNUpdate1')

    if ( use_nguardrail .and. num_bgc_vegp>0 ) then
       call t_startf('CNPrecisionControl')
       call CNPrecisionControl(bounds, num_bgc_vegp, filter_bgc_vegp, &
            cnveg_carbonstate_inst, c13_cnveg_carbonstate_inst, &
            c14_cnveg_carbonstate_inst, cnveg_nitrogenstate_inst)
       call t_stopf('CNPrecisionControl')
    end if

    call t_startf('SoilBiogeochemStateUpdate1')
    call SoilBiogeochemNStateUpdate1(num_bgc_soilc, filter_bgc_soilc,  &
         soilbiogeochem_state_inst, soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst)
    call t_stopf('SoilBiogeochemStateUpdate1')


    !--------------------------------------------
    ! Calculate vertical mixing of soil and litter pools
    !--------------------------------------------

    call t_startf('SoilBiogeochemLittVertTransp')

    call SoilBiogeochemLittVertTransp(bounds, num_bgc_soilc, filter_bgc_soilc,            &
         active_layer_inst, soilbiogeochem_state_inst,                             &
         soilbiogeochem_carbonstate_inst, soilbiogeochem_carbonflux_inst,         &
         c13_soilbiogeochem_carbonstate_inst, c13_soilbiogeochem_carbonflux_inst, &
         c14_soilbiogeochem_carbonstate_inst, c14_soilbiogeochem_carbonflux_inst, &
         soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)

    call t_stopf('SoilBiogeochemLittVertTransp')

    !--------------------------------------------
    ! Calculate the gap mortality carbon and nitrogen fluxes
    !--------------------------------------------

    if_bgc_vegp1: if(num_bgc_vegp>0)then

       call t_startf('CNGapMortality')

       call CNGapMortality (bounds, num_bgc_vegp, filter_bgc_vegp,                                                   &
            dgvs_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst,  soilbiogeochem_nitrogenflux_inst,          &
            cnveg_carbonflux_inst, cnveg_nitrogenflux_inst,  canopystate_inst,                                       &   
            leaf_prof_patch=soilbiogeochem_state_inst%leaf_prof_patch(begp:endp, 1:nlevdecomp_full),   &
            froot_prof_patch=soilbiogeochem_state_inst%froot_prof_patch(begp:endp, 1:nlevdecomp_full), & 
            croot_prof_patch=soilbiogeochem_state_inst%croot_prof_patch(begp:endp, 1:nlevdecomp_full), &
            stem_prof_patch=soilbiogeochem_state_inst%stem_prof_patch(begp:endp, 1:nlevdecomp_full))   

       call t_stopf('CNGapMortality')

       !--------------------------------------------------------------------------
       ! Update2 (gap mortality)
       ! The state updates are still called for the matrix solution (use_matrixn
       ! and use_soil_matrixcn) but most of the state updates are done after
       ! the matrix multiply in VegMatrix and SoilMatrix.
       !--------------------------------------------------------------------------

       call t_startf('CNUpdate2')

       ! Set the carbon isotopic fluxes for gap mortality
       if ( use_c13 ) then
          call CIsoFlux2(num_bgc_vegp, filter_bgc_vegp,                                  &
               soilbiogeochem_state_inst, cnveg_carbonflux_inst, cnveg_carbonstate_inst, &
               iso_cnveg_carbonflux_inst=c13_cnveg_carbonflux_inst,                      &
               iso_cnveg_carbonstate_inst=c13_cnveg_carbonstate_inst,                    &
               isotope='c13')
       end if
       if ( use_c14 ) then
          call CIsoFlux2(num_bgc_vegp, filter_bgc_vegp,                                  &
               soilbiogeochem_state_inst, cnveg_carbonflux_inst, cnveg_carbonstate_inst, &
               iso_cnveg_carbonflux_inst=c14_cnveg_carbonflux_inst,                      &
               iso_cnveg_carbonstate_inst=c14_cnveg_carbonstate_inst,                    &
               isotope='c14')
       end if

       ! Update all the prognostic carbon state variables affected by gap-phase mortality fluxes
       call CStateUpdate2(num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
            cnveg_carbonflux_inst, cnveg_carbonstate_inst, soilbiogeochem_carbonstate_inst, &
            soilbiogeochem_carbonflux_inst)
       if ( use_c13 ) then
          call CStateUpdate2(num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
               c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst, c13_soilbiogeochem_carbonstate_inst, &
               c13_soilbiogeochem_carbonflux_inst)
       end if
       if ( use_c14 ) then
          call CStateUpdate2(num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
               c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst, c14_soilbiogeochem_carbonstate_inst, &
               c14_soilbiogeochem_carbonflux_inst)
       end if

       ! Update all the prognostic nitrogen state variables affected by gap-phase mortality fluxes
       call NStateUpdate2(num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
            cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst,soilbiogeochem_nitrogenstate_inst, &
            soilbiogeochem_nitrogenflux_inst)

       !--------------------------------------------------------------------------
       ! Update2h (harvest)
       ! The state updates are still called for the matrix solution (use_matrixn
       ! and use_soil_matrixcn) but most of the state updates are done after
       ! the matrix multiply in VegMatrix and SoilMatrix.
       !--------------------------------------------------------------------------

       ! Set harvest mortality routine
       if (get_do_harvest()) then
          call CNHarvest(num_bgc_vegp, filter_bgc_vegp, &
               soilbiogeochem_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, &
               cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
       end if

       if ( use_c13 ) then
          call CIsoFlux2h(num_bgc_vegp, filter_bgc_vegp,                      &
               soilbiogeochem_state_inst,                                     &
               cnveg_carbonflux_inst, cnveg_carbonstate_inst,                 &
               c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst,         &                         
               isotope='c13')
       end if
       if ( use_c14 ) then
          call CIsoFlux2h(num_bgc_vegp, filter_bgc_vegp,                      &
               soilbiogeochem_state_inst,                                     &
               cnveg_carbonflux_inst, cnveg_carbonstate_inst,                 &
               c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst,         &                         
               isotope='c14')
       end if

       call CStateUpdate2h( num_bgc_soilc, filter_bgc_soilc,  num_bgc_vegp, filter_bgc_vegp, &
            cnveg_carbonflux_inst, cnveg_carbonstate_inst, soilbiogeochem_carbonstate_inst, &
            soilbiogeochem_carbonflux_inst)
       if ( use_c13 ) then
          call CStateUpdate2h(num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
            c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst, c13_soilbiogeochem_carbonstate_inst, &
            c13_soilbiogeochem_carbonflux_inst)
       end if
       if ( use_c14 ) then
          call CStateUpdate2h(num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
               c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst, c14_soilbiogeochem_carbonstate_inst, &
               c14_soilbiogeochem_carbonflux_inst)
       end if

       call NStateUpdate2h(num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
            cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, soilbiogeochem_nitrogenstate_inst, &
            soilbiogeochem_nitrogenflux_inst)

       !--------------------------------------------
       ! Update2g (gross unrepresented landcover change)
       !--------------------------------------------

       ! Set gross unrepresented landcover change mortality routine 
       if (get_do_grossunrep()) then
          call CNGrossUnrep(num_bgc_vegp, filter_bgc_vegp, &
               soilbiogeochem_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, &
               cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
       end if

       if ( use_c13 ) then
          call CIsoFlux2g(num_bgc_vegp, filter_bgc_vegp,                      &
               soilbiogeochem_state_inst,                                     &
               cnveg_carbonflux_inst, cnveg_carbonstate_inst,                 &
               c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst,         &                         
               isotope='c13')
       end if
       if ( use_c14 ) then
          call CIsoFlux2g(num_bgc_vegp, filter_bgc_vegp,                      &
               soilbiogeochem_state_inst,                                     &
               cnveg_carbonflux_inst, cnveg_carbonstate_inst,                 &
               c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst,         &                         
               isotope='c14')
       end if

       call CStateUpdate2g( num_bgc_soilc, filter_bgc_soilc,  num_bgc_vegp, filter_bgc_vegp, &
            cnveg_carbonflux_inst, cnveg_carbonstate_inst, &
            soilbiogeochem_carbonstate_inst, soilbiogeochem_carbonflux_inst)
       if ( use_c13 ) then
          call CStateUpdate2g(num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
               c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst, &
               c13_soilbiogeochem_carbonstate_inst, c13_soilbiogeochem_carbonflux_inst)
       end if
       if ( use_c14 ) then
          call CStateUpdate2g(num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
               c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst, &
               c14_soilbiogeochem_carbonstate_inst, c14_soilbiogeochem_carbonflux_inst)
       end if

       call NStateUpdate2g(num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
            cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, &
            soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)

       call t_stopf('CNUpdate2')

    end if if_bgc_vegp1
       
    if ( use_nguardrail .and. num_bgc_vegp>0 ) then

       call t_startf('CNPrecisionControl')
       call CNPrecisionControl(bounds, num_bgc_vegp, filter_bgc_vegp, &
            cnveg_carbonstate_inst, c13_cnveg_carbonstate_inst, &
            c14_cnveg_carbonstate_inst, cnveg_nitrogenstate_inst)
       call t_stopf('CNPrecisionControl')
       
    end if
    
    !--------------------------------------------
    ! Calculate loss fluxes from wood products pools
    ! and update product pool state variables
    !--------------------------------------------

    call t_startf('CNWoodProducts')
    
    call c_products_inst%SetValues(bounds,0._r8)
    if (use_c13) call c13_products_inst%SetValues(bounds,0._r8)
    if (use_c14) call c14_products_inst%SetValues(bounds,0._r8)
    call n_products_inst%SetValues(bounds,0._r8)
          
    if(use_fates_bgc) then
       call clm_fates%wrap_WoodProducts(bounds, num_bgc_soilc, filter_bgc_soilc, c_products_inst, n_products_inst)
    end if

    if_bgc_vegp2: if(num_bgc_vegp>0)then
       call c_products_inst%UpdateProducts(bounds, &
            num_bgc_vegp, filter_bgc_vegp, &
            dwt_wood_product_gain_patch = cnveg_carbonflux_inst%dwt_wood_productc_gain_patch(begp:endp), &
            gru_wood_product_gain_patch = cnveg_carbonflux_inst%gru_wood_productc_gain_patch(begp:endp), &
            wood_harvest_patch = cnveg_carbonflux_inst%wood_harvestc_patch(begp:endp), &
            dwt_crop_product_gain_patch = cnveg_carbonflux_inst%dwt_crop_productc_gain_patch(begp:endp), &
            crop_harvest_to_cropprod_patch = cnveg_carbonflux_inst%crop_harvestc_to_cropprodc_patch(begp:endp))
 
       if (use_c13) then
          call c13_products_inst%UpdateProducts(bounds, &
               num_bgc_vegp, filter_bgc_vegp, &
               dwt_wood_product_gain_patch = c13_cnveg_carbonflux_inst%dwt_wood_productc_gain_patch(begp:endp), &
               gru_wood_product_gain_patch = c13_cnveg_carbonflux_inst%gru_wood_productc_gain_patch(begp:endp), &
               wood_harvest_patch = c13_cnveg_carbonflux_inst%wood_harvestc_patch(begp:endp), &
               dwt_crop_product_gain_patch = c13_cnveg_carbonflux_inst%dwt_crop_productc_gain_patch(begp:endp), &
               crop_harvest_to_cropprod_patch = c13_cnveg_carbonflux_inst%crop_harvestc_to_cropprodc_patch(begp:endp))
       end if
    
       if (use_c14) then
          call c14_products_inst%UpdateProducts(bounds, &
               num_bgc_vegp, filter_bgc_vegp, &
               dwt_wood_product_gain_patch = c14_cnveg_carbonflux_inst%dwt_wood_productc_gain_patch(begp:endp), &
               gru_wood_product_gain_patch = c14_cnveg_carbonflux_inst%gru_wood_productc_gain_patch(begp:endp), &
               wood_harvest_patch = c14_cnveg_carbonflux_inst%wood_harvestc_patch(begp:endp), &
               dwt_crop_product_gain_patch = c14_cnveg_carbonflux_inst%dwt_crop_productc_gain_patch(begp:endp), &
               crop_harvest_to_cropprod_patch = c14_cnveg_carbonflux_inst%crop_harvestc_to_cropprodc_patch(begp:endp))
       end if

       call n_products_inst%UpdateProducts(bounds, &
            num_bgc_vegp, filter_bgc_vegp, &
            dwt_wood_product_gain_patch = cnveg_nitrogenflux_inst%dwt_wood_productn_gain_patch(begp:endp), &
            gru_wood_product_gain_patch = cnveg_nitrogenflux_inst%gru_wood_productn_gain_patch(begp:endp), &
            wood_harvest_patch = cnveg_nitrogenflux_inst%wood_harvestn_patch(begp:endp), &
            dwt_crop_product_gain_patch = cnveg_nitrogenflux_inst%dwt_crop_productn_gain_patch(begp:endp), &
            crop_harvest_to_cropprod_patch = cnveg_nitrogenflux_inst%crop_harvestn_to_cropprodn_patch(begp:endp))
       
    end if if_bgc_vegp2
 
    call c_products_inst%ComputeProductSummaryVars(bounds)
    if (use_c13) call c13_products_inst%ComputeProductSummaryVars(bounds)
    if (use_c14) call c14_products_inst%ComputeProductSummaryVars(bounds)
    call n_products_inst%ComputeProductSummaryVars(bounds)


    call c_products_inst%ComputeSummaryVars(bounds)
    if (use_c13) call c13_products_inst%ComputeSummaryVars(bounds)
    if (use_c14) call c14_products_inst%ComputeSummaryVars(bounds)
    call n_products_inst%ComputeSummaryVars(bounds)
    
    call t_stopf('CNWoodProducts')
       
    !--------------------------------------------
    ! Calculate fire area and fluxes
    !--------------------------------------------

    if_bgc_vegp3: if(num_bgc_vegp>0)then

       call t_startf('CNFire')
       call cnfire_method%CNFireArea(bounds, num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
            num_exposedvegp, filter_exposedvegp, num_noexposedvegp, filter_noexposedvegp, &
            atm2lnd_inst, energyflux_inst, saturated_excess_runoff_inst, waterdiagnosticbulk_inst, wateratm2lndbulk_inst, &
            waterstatebulk_inst, soilstate_inst, soil_water_retention_curve, &
            crop_inst, cnveg_state_inst, cnveg_carbonstate_inst, &
            totlitc_col=soilbiogeochem_carbonstate_inst%totlitc_col(begc:endc), &
            decomp_cpools_vr_col=soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools), &
            t_soi17cm_col=temperature_inst%t_soi17cm_col(begc:endc))

       call cnfire_method%CNFireFluxes(bounds, num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp,                        &
            num_actfirec, filter_actfirec, num_actfirep, filter_actfirep,                                                             &
            dgvs_inst, cnveg_state_inst,                                                                                              &
            cnveg_carbonstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst,                         &
            soilbiogeochem_carbonflux_inst,                                                                                           &
            leaf_prof_patch=soilbiogeochem_state_inst%leaf_prof_patch(begp:endp, 1:nlevdecomp_full),                                  &
            froot_prof_patch=soilbiogeochem_state_inst%froot_prof_patch(begp:endp, 1:nlevdecomp_full),                                &
            croot_prof_patch=soilbiogeochem_state_inst%croot_prof_patch(begp:endp, 1:nlevdecomp_full),                                &
            stem_prof_patch=soilbiogeochem_state_inst%stem_prof_patch(begp:endp, 1:nlevdecomp_full),                                  &
            totsomc_col=soilbiogeochem_carbonstate_inst%totsomc_col(begc:endc),                                                       &
            decomp_cpools_vr_col=soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools),   &
            decomp_npools_vr_col=soilbiogeochem_nitrogenstate_inst%decomp_npools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools), &
            somc_fire_col=soilbiogeochem_carbonflux_inst%somc_fire_col(begc:endc))
       call t_stopf('CNFire')


       !--------------------------------------------------------------------------
       ! Update3
       ! The state updates are still called for the matrix solution (use_matrixn
       ! and use_soil_matrixcn) but most of the state updates are done after
       ! the matrix multiply in VegMatrix and SoilMatrix.
       !--------------------------------------------------------------------------

       call t_startf('CNUpdate3')
       if ( use_c13 ) then
          call CIsoFlux3(num_bgc_vegp, filter_bgc_vegp,                                   &
               soilbiogeochem_state_inst , soilbiogeochem_carbonstate_inst,         &
               cnveg_carbonflux_inst, cnveg_carbonstate_inst,                       &
               c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst,               &
               c13_soilbiogeochem_carbonstate_inst, &
               isotope='c13')
       end if
       if ( use_c14 ) then
          call CIsoFlux3(num_bgc_vegp, filter_bgc_vegp,                                   &
               soilbiogeochem_state_inst , soilbiogeochem_carbonstate_inst,         &
               cnveg_carbonflux_inst, cnveg_carbonstate_inst,                       &
               c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst,               &
               c14_soilbiogeochem_carbonstate_inst, &
               isotope='c14')
       end if

       call CStateUpdate3( num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
            cnveg_carbonflux_inst, cnveg_carbonstate_inst, soilbiogeochem_carbonstate_inst, &
            soilbiogeochem_carbonflux_inst)

       if ( use_c13 ) then
          call CStateUpdate3( num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
               c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst, c13_soilbiogeochem_carbonstate_inst, &
               c13_soilbiogeochem_carbonflux_inst)
       end if

       if ( use_c14 ) then
          call CStateUpdate3( num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
               c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst, c14_soilbiogeochem_carbonstate_inst, &
               c14_soilbiogeochem_carbonflux_inst)

          call C14Decay(bounds, num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
               c14_cnveg_carbonstate_inst, c14_soilbiogeochem_carbonstate_inst, &
               c14_cnveg_carbonflux_inst,  c14_soilbiogeochem_carbonflux_inst)
       end if
       call t_stopf('CNUpdate3')

    end if if_bgc_vegp3
    
    if ( use_nguardrail .and. num_bgc_vegp>0 ) then
       call t_startf('CNPrecisionControl')
       call CNPrecisionControl(bounds, num_bgc_vegp, filter_bgc_vegp, &
            cnveg_carbonstate_inst, c13_cnveg_carbonstate_inst, &
            c14_cnveg_carbonstate_inst, cnveg_nitrogenstate_inst)
       call t_stopf('CNPrecisionControl')
    end if

    end associate

  end subroutine CNDriverNoLeaching
  
  !-----------------------------------------------------------------------
  subroutine CNDriverLeaching(bounds, &
       num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
       num_actfirec, filter_actfirec, num_actfirep, filter_actfirep,&
       waterstatebulk_inst, waterfluxbulk_inst, &
       soilstate_inst, cnveg_state_inst, &
       cnveg_carbonflux_inst,cnveg_carbonstate_inst,soilbiogeochem_carbonstate_inst, &
       soilbiogeochem_carbonflux_inst,soilbiogeochem_state_inst, &
       cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, &
       soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst, &
       c13_cnveg_carbonstate_inst,c14_cnveg_carbonstate_inst, &
       c13_cnveg_carbonflux_inst,c14_cnveg_carbonflux_inst, &
       c13_soilbiogeochem_carbonstate_inst,c14_soilbiogeochem_carbonstate_inst,&
       c13_soilbiogeochem_carbonflux_inst,c14_soilbiogeochem_carbonflux_inst)
    !
    ! !DESCRIPTION:
    ! Update the nitrogen leaching rate as a function of soluble mineral N and total soil water outflow.
    ! Also update nitrogen state variables         
    !
    ! !USES:
    use SoilBiogeochemNLeachingMod, only: SoilBiogeochemNLeaching
    use CNNStateUpdate3Mod        , only: NStateUpdate3
    use CNNStateUpdate3Mod        , only: NStateUpdateLeaching
    use CNVegMatrixMod            , only: CNVegMatrix
    use CNSoilMatrixMod           , only: CNSoilMatrix
    use clm_time_manager          , only: is_first_step_of_this_run_segment,is_beg_curr_year,is_end_curr_year,get_curr_date
    use CNSharedParamsMod         , only: use_matrixcn
    use SoilBiogeochemDecompCascadeConType, only: use_soil_matrixcn
    !
    ! !ARGUMENTS:
    type(bounds_type)                       , intent(in)    :: bounds  
    integer                                 , intent(in)    :: num_bgc_soilc         ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_bgc_soilc(:)   ! filter for soil columns
    integer                                 , intent(in)    :: num_bgc_vegp         ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_bgc_vegp(:)   ! filter for soil patches
    integer                                 , intent(in)    :: num_actfirec         ! number of soil columns on fire in filter
    integer                                 , intent(in)    :: filter_actfirec(:)   ! filter for soil columns on fire
    integer                                 , intent(in)    :: num_actfirep         ! number of soil patches on fire in filter
    integer                                 , intent(in)    :: filter_actfirep(:)   ! filter for soil patches on fire
    type(waterstatebulk_type)                   , intent(in)    :: waterstatebulk_inst
    type(waterfluxbulk_type)                    , intent(inout)    :: waterfluxbulk_inst
    type(cnveg_state_type)                  , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonflux_type)             , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)            , intent(inout) :: cnveg_carbonstate_inst
    type(soilstate_type)                    , intent(inout) :: soilstate_inst
    type(soilbiogeochem_state_type)         , intent(inout) :: soilbiogeochem_state_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: soilbiogeochem_carbonstate_inst
    type(cnveg_nitrogenflux_type)           , intent(inout) :: cnveg_nitrogenflux_inst
    type(cnveg_nitrogenstate_type)          , intent(inout) :: cnveg_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(cnveg_carbonstate_type)            , intent(inout) :: c13_cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)            , intent(inout) :: c14_cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)             , intent(inout) :: c13_cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)             , intent(inout) :: c14_cnveg_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c13_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: c13_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c14_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: c14_soilbiogeochem_carbonflux_inst
    integer p,fp,yr,mon,day,sec
    !-----------------------------------------------------------------------
  
    ! Mineral nitrogen dynamics (deposition, fixation, leaching)
    
    call t_startf('SoilBiogeochemNLeaching')
    call SoilBiogeochemNLeaching(bounds, num_bgc_soilc, filter_bgc_soilc, &
         waterstatebulk_inst, waterfluxbulk_inst, soilbiogeochem_nitrogenstate_inst, &
         soilbiogeochem_nitrogenflux_inst)
    call NStateUpdateLeaching(num_bgc_soilc, filter_bgc_soilc, &
         soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst)
    call t_stopf('SoilBiogeochemNLeaching')

    ! Nitrogen state variable update, mortality fluxes.
    if(num_bgc_vegp>0)then
       call t_startf('NUpdate3')
       call NStateUpdate3(num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
            cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, &
            soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst)
       call t_stopf('NUpdate3')
    end if
    
    !--------------------------------------------------------------------------
    ! Solve the matrix solution and do the state update for matrix solution as
    ! part of that
    !--------------------------------------------------------------------------

    if ( use_matrixcn ) then
       call t_startf('CNVMatrix')
       call CNVegMatrix(bounds, num_bgc_vegp, filter_bgc_vegp(1:num_bgc_vegp), &
          num_actfirep, filter_actfirep, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, &
          cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, cnveg_state_inst,soilbiogeochem_nitrogenflux_inst, &
          c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst, c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst)
       call t_stopf('CNVMatrix')
    end if

    if(use_soil_matrixcn)then
       call t_startf('CNSoilMatrix')
       call CNSoilMatrix(bounds,num_bgc_soilc, filter_bgc_soilc(1:num_bgc_soilc), num_actfirec, filter_actfirec, &
       cnveg_carbonflux_inst,soilbiogeochem_carbonstate_inst, &
       soilbiogeochem_carbonflux_inst,soilbiogeochem_state_inst, &
       cnveg_nitrogenflux_inst, soilbiogeochem_nitrogenflux_inst, &
       soilbiogeochem_nitrogenstate_inst,c13_soilbiogeochem_carbonstate_inst,&
       c13_soilbiogeochem_carbonflux_inst,c14_soilbiogeochem_carbonstate_inst,&
       c14_soilbiogeochem_carbonflux_inst)  
    call t_stopf('CNSoilMatrix')
    end if
    
  end subroutine CNDriverLeaching

  !-----------------------------------------------------------------------
  subroutine CNDriverSummarizeStates(bounds, num_allc, filter_allc, &
       num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
       cnveg_carbonstate_inst, c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst, &
       cnveg_nitrogenstate_inst, &
       soilbiogeochem_carbonstate_inst, &
       c13_soilbiogeochem_carbonstate_inst, &
       c14_soilbiogeochem_carbonstate_inst, &
       soilbiogeochem_nitrogenstate_inst)
    !
    ! !DESCRIPTION:
    ! Call to all CN and SoilBiogeochem summary routines, for state variables
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type)                       , intent(in)    :: bounds  
    integer                                 , intent(in)    :: num_allc          ! number of columns in allc filter
    integer                                 , intent(in)    :: filter_allc(:)    ! filter for all active columns
    integer                                 , intent(in)    :: num_bgc_soilc         ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_bgc_soilc(:)   ! filter for soil columns
    integer                                 , intent(in)    :: num_bgc_vegp         ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_bgc_vegp(:)   ! filter for soil patches
    type(cnveg_carbonstate_type)            , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)            , intent(inout) :: c13_cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)            , intent(inout) :: c14_cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type)          , intent(inout) :: cnveg_nitrogenstate_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c13_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c14_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: begc,endc

    character(len=*), parameter :: subname = 'CNDriverSummarizeStates'
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc

    call t_startf('CNsum')

    ! ----------------------------------------------
    ! cnveg carbon/nitrogen state summary
    ! ----------------------------------------------
    call cnveg_carbonstate_inst%Summary(bounds, num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp)

    if ( use_c13 ) then
       call c13_cnveg_carbonstate_inst%Summary(bounds, num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp)
    end if

    if ( use_c14 ) then
       call c14_cnveg_carbonstate_inst%Summary(bounds, num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp)
    end if

    ! ----------------------------------------------
    ! soilbiogeochem carbon/nitrogen state summary
    ! RGK 02-23: soilbiogeochem summary now depends on
    !            cnveg summary, swapped call order
    ! ----------------------------------------------

    call soilbiogeochem_carbonstate_inst%summary(bounds, num_allc, filter_allc, &
         num_bgc_soilc, filter_bgc_soilc, cnveg_carbonstate_inst)
    if ( use_c13 ) then
       call c13_soilbiogeochem_carbonstate_inst%summary(bounds, num_allc, filter_allc, &
            num_bgc_soilc, filter_bgc_soilc, c13_cnveg_carbonstate_inst)
    end if
    if ( use_c14 ) then
       call c14_soilbiogeochem_carbonstate_inst%summary(bounds, num_allc, filter_allc, &
            num_bgc_soilc, filter_bgc_soilc, c14_cnveg_carbonstate_inst)
    end if
    
    
    ! RGK 02-23: This call will be moved to after cnveg nitr summary when we
    !            couple in FATES N
    

    call cnveg_nitrogenstate_inst%Summary(bounds, num_bgc_soilc, filter_bgc_soilc, &
         num_bgc_vegp, filter_bgc_vegp)

    call soilbiogeochem_nitrogenstate_inst%summary(bounds, num_allc, filter_allc, &
         num_bgc_soilc, filter_bgc_soilc, cnveg_nitrogenstate_inst)
    

    call t_stopf('CNsum')

  end subroutine CNDriverSummarizeStates

  !-----------------------------------------------------------------------
  subroutine CNDriverSummarizeFluxes(bounds, &
       num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
       cnveg_carbonflux_inst, c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst, &
       cnveg_nitrogenflux_inst, &
       c_products_inst, c13_products_inst, c14_products_inst, &
       soilbiogeochem_carbonflux_inst, &
       c13_soilbiogeochem_carbonflux_inst, &
       c14_soilbiogeochem_carbonflux_inst, &
       soilbiogeochem_carbonstate_inst, &
       c13_soilbiogeochem_carbonstate_inst, &
       c14_soilbiogeochem_carbonstate_inst, &
       soilbiogeochem_nitrogenstate_inst, &
       soilbiogeochem_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! Call to all CN and SoilBiogeochem summary routines, for state variables
    !
    ! !USES:
    use clm_varpar                        , only: ndecomp_cascade_transitions
    !
    ! !ARGUMENTS:
    type(bounds_type)                       , intent(in)    :: bounds  
    integer                                 , intent(in)    :: num_bgc_soilc         ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_bgc_soilc(:)   ! filter for soil columns
    integer                                 , intent(in)    :: num_bgc_vegp         ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_bgc_vegp(:)   ! filter for soil patches
    type(cnveg_carbonflux_type)             , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)             , intent(inout) :: c13_cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)             , intent(inout) :: c14_cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)           , intent(inout) :: cnveg_nitrogenflux_inst
    type(cn_products_type)                  , intent(in)    :: c_products_inst
    type(cn_products_type)                  , intent(in)    :: c13_products_inst
    type(cn_products_type)                  , intent(in)    :: c14_products_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: c13_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: c14_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(in)    :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonstate_type)   , intent(in)    :: c13_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonstate_type)   , intent(in)    :: c14_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(in)    :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: begp,endp
    integer :: begc,endc
    integer :: begg,endg

    character(len=*), parameter :: subname = 'CNDriverSummarizeFluxes'
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg = bounds%endg

    call t_startf('CNsum')

    ! ----------------------------------------------
    ! soilbiogeochem carbon/nitrogen flux summary
    ! ----------------------------------------------

    call soilbiogeochem_carbonflux_inst%Summary(bounds, num_bgc_soilc, filter_bgc_soilc, &
         num_bgc_vegp, filter_bgc_vegp, &
         soilbiogeochem_carbonflux_inst%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions), &
         soilbiogeochem_carbonstate_inst%cwdc_col(begc:endc), &
         soilbiogeochem_nitrogenstate_inst%cwdn_col(begc:endc), &
         leafc_to_litter_patch=cnveg_carbonflux_inst%leafc_to_litter_patch, &
         frootc_to_litter_patch=cnveg_carbonflux_inst%frootc_to_litter_patch)
    if ( use_c13 ) then
       call c13_soilbiogeochem_carbonflux_inst%Summary(bounds, num_bgc_soilc, filter_bgc_soilc, &
         num_bgc_vegp, filter_bgc_vegp, &
         c13_soilbiogeochem_carbonflux_inst%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions), &
         c13_soilbiogeochem_carbonstate_inst%cwdc_col(begc:endc), &
         soilbiogeochem_nitrogenstate_inst%cwdn_col(begc:endc), &
         leafc_to_litter_patch=c13_cnveg_carbonflux_inst%leafc_to_litter_patch, &
         frootc_to_litter_patch=c13_cnveg_carbonflux_inst%frootc_to_litter_patch)
    end if
    if ( use_c14 ) then
       call c14_soilbiogeochem_carbonflux_inst%Summary(bounds, num_bgc_soilc, filter_bgc_soilc, &
         num_bgc_vegp, filter_bgc_vegp, &
         c14_soilbiogeochem_carbonflux_inst%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions), &
         c14_soilbiogeochem_carbonstate_inst%cwdc_col(begc:endc), &
         soilbiogeochem_nitrogenstate_inst%cwdn_col(begc:endc), &
         leafc_to_litter_patch=c14_cnveg_carbonflux_inst%leafc_to_litter_patch, &
         frootc_to_litter_patch=c14_cnveg_carbonflux_inst%frootc_to_litter_patch)
    end if
    call soilbiogeochem_nitrogenflux_inst%Summary(bounds, num_bgc_soilc, filter_bgc_soilc)

    ! ----------------------------------------------
    ! cnveg carbon/nitrogen flux summary
    ! ----------------------------------------------

    if_bgc_vegp: if(num_bgc_vegp>0) then
       call t_startf('CNvegCflux_summary')
       call cnveg_carbonflux_inst%Summary(bounds, num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
            isotope='bulk', &
            soilbiogeochem_hr_col=soilbiogeochem_carbonflux_inst%hr_col(begc:endc), &
            soilbiogeochem_cwdhr_col=soilbiogeochem_carbonflux_inst%cwdhr_col(begc:endc), &
            soilbiogeochem_lithr_col=soilbiogeochem_carbonflux_inst%lithr_col(begc:endc), &  
            soilbiogeochem_decomp_cascade_ctransfer_col=&
            soilbiogeochem_carbonflux_inst%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions), &
            product_closs_grc=c_products_inst%product_loss_grc(begg:endg))
       
       if ( use_c13 ) then
          call c13_cnveg_carbonflux_inst%Summary(bounds, num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
               isotope='c13', &
               soilbiogeochem_hr_col=c13_soilbiogeochem_carbonflux_inst%hr_col(begc:endc), &
               soilbiogeochem_cwdhr_col=c13_soilbiogeochem_carbonflux_inst%cwdhr_col(begc:endc), &
               soilbiogeochem_lithr_col=c13_soilbiogeochem_carbonflux_inst%lithr_col(begc:endc), &  
               soilbiogeochem_decomp_cascade_ctransfer_col=&
               c13_soilbiogeochem_carbonflux_inst%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions), &
               product_closs_grc=c13_products_inst%product_loss_grc(begg:endg))
       end if
       
       if ( use_c14 ) then
          call c14_cnveg_carbonflux_inst%Summary(bounds, num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &
               isotope='c14', &
               soilbiogeochem_hr_col=c14_soilbiogeochem_carbonflux_inst%hr_col(begc:endc), &
               soilbiogeochem_cwdhr_col=c14_soilbiogeochem_carbonflux_inst%cwdhr_col(begc:endc), &
               soilbiogeochem_lithr_col=c14_soilbiogeochem_carbonflux_inst%lithr_col(begc:endc), &  
               soilbiogeochem_decomp_cascade_ctransfer_col=&
               c14_soilbiogeochem_carbonflux_inst%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions), &
               product_closs_grc=c14_products_inst%product_loss_grc(begg:endg))
       end if
       call t_stopf('CNvegCflux_summary')

       call cnveg_nitrogenflux_inst%Summary(bounds, num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp)
    end if if_bgc_vegp
    
    call t_stopf('CNsum')

  end subroutine CNDriverSummarizeFluxes

end  module CNDriverMod
