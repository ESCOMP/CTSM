module CNDriverMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Ecosystem dynamics: phenology, vegetation
  !
  ! !USES:
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use clm_varctl                      , only : use_c13, use_c14, use_fates, use_dynroot
  use dynSubgridControlMod            , only : get_do_harvest
  use decompMod                       , only : bounds_type
  use perf_mod                        , only : t_startf, t_stopf
  use clm_varctl                      , only : use_century_decomp, use_nitrif_denitrif, use_nguardrail
  use clm_varctl                      , only : use_crop
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
  use WaterstateType                  , only : waterstate_type
  use WaterfluxType                   , only : waterflux_type
  use atm2lndType                     , only : atm2lnd_type
  use SoilStateType                   , only : soilstate_type
  use TemperatureType                 , only : temperature_type 
  use PhotosynthesisMod               , only : photosyns_type
  use ch4Mod                          , only : ch4_type
  use EnergyFluxType                  , only : energyflux_type
  use SoilHydrologyType               , only : soilhydrology_type
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
    use CNFireMethodMod             , only : cnfire_method_type
    use SoilBiogeochemCompetitionMod, only : SoilBiogeochemCompetitionInit
    !
    ! !ARGUMENTS:
    type(bounds_type)                      , intent(in)    :: bounds      
    character(len=*)                       , intent(in)    :: NLFilename     ! Namelist filename
    class(cnfire_method_type)              , intent(inout) :: cnfire_method 
    !-----------------------------------------------------------------------
    call SoilBiogeochemCompetitionInit(bounds)
    call CNPhenologyInit(bounds)
    call cnfire_method%CNFireInit(bounds, NLFilename)
    
  end subroutine CNDriverInit

  !-----------------------------------------------------------------------
  subroutine CNDriverNoLeaching(bounds,                                                    &
       num_soilc, filter_soilc, num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb, &
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
       atm2lnd_inst, waterstate_inst, waterflux_inst,                                      &
       canopystate_inst, soilstate_inst, temperature_inst, crop_inst, ch4_inst,            &
       dgvs_inst, photosyns_inst, soilhydrology_inst, energyflux_inst,                     &
       nutrient_competition_method, cnfire_method, dribble_crophrv_xsmrpool_2atm)
    !
    ! !DESCRIPTION:
    ! The core CN code is executed here. Calculates fluxes for maintenance
    ! respiration, decomposition, allocation, phenology, and growth respiration.
    ! These routines happen on the radiation time step so that canopy structure
    ! stays synchronized with albedo calculations.
    !
    ! !USES:
    use clm_varpar                        , only: nlevgrnd, nlevdecomp_full 
    use clm_varpar                        , only: nlevdecomp, ndecomp_cascade_transitions, ndecomp_pools
    use subgridAveMod                     , only: p2c, p2c_2d
    use CropType                          , only: crop_type
    use CNNDynamicsMod                    , only: CNNDeposition,CNNFixation, CNNFert, CNSoyfix,CNFreeLivingFixation
    use CNMRespMod                        , only: CNMResp
    use CNFUNMod                          , only: CNFUNInit  !, CNFUN 
    use CNPhenologyMod                    , only: CNPhenology
    use CNGRespMod                        , only: CNGResp
    use CNFireMethodMod                   , only: cnfire_method_type
    use CNCIsoFluxMod                     , only: CIsoFlux1, CIsoFlux2, CIsoFlux2h, CIsoFlux3
    use CNC14DecayMod                     , only: C14Decay
    use CNCStateUpdate1Mod                , only: CStateUpdate1,CStateUpdate0
    use CNCStateUpdate2Mod                , only: CStateUpdate2, CStateUpdate2h
    use CNCStateUpdate3Mod                , only: CStateUpdate3
    use CNNStateUpdate1Mod                , only: NStateUpdate1
    use CNNStateUpdate2Mod                , only: NStateUpdate2, NStateUpdate2h
    use CNGapMortalityMod                 , only: CNGapMortality
    use CNSharedParamsMod                 , only: use_fun
    use dynHarvestMod                     , only: CNHarvest
    use SoilBiogeochemDecompCascadeBGCMod , only: decomp_rate_constants_bgc
    use SoilBiogeochemDecompCascadeCNMod  , only: decomp_rate_constants_cn
    use SoilBiogeochemCompetitionMod      , only: SoilBiogeochemCompetition
    use SoilBiogeochemDecompMod           , only: SoilBiogeochemDecomp
    use SoilBiogeochemLittVertTranspMod   , only: SoilBiogeochemLittVertTransp
    use SoilBiogeochemPotentialMod        , only: SoilBiogeochemPotential 
    use SoilBiogeochemVerticalProfileMod  , only: SoilBiogeochemVerticalProfile
    use SoilBiogeochemNitrifDenitrifMod   , only: SoilBiogeochemNitrifDenitrif
    use SoilBiogeochemNStateUpdate1Mod    , only: SoilBiogeochemNStateUpdate1
    use NutrientCompetitionMethodMod      , only: nutrient_competition_method_type
    use CNRootDynMod                      , only: CNRootDyn
    use CNPrecisionControlMod             , only: CNPrecisionControl
    !
    ! !ARGUMENTS:
    type(bounds_type)                       , intent(in)    :: bounds  
    integer                                 , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                                 , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    integer                                 , intent(in)    :: num_pcropp        ! number of prog. crop patches in filter
    integer                                 , intent(in)    :: filter_pcropp(:)  ! filter for prognostic crop patches
    logical                                 , intent(in)    :: doalb             ! true = surface albedo calculation time step
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
    type(atm2lnd_type)                      , intent(in)    :: atm2lnd_inst
    type(waterstate_type)                   , intent(in)    :: waterstate_inst
    type(waterflux_type)                    , intent(inout)    :: waterflux_inst
    type(canopystate_type)                  , intent(inout)    :: canopystate_inst
    type(soilstate_type)                    , intent(inout) :: soilstate_inst
    type(temperature_type)                  , intent(inout) :: temperature_inst
    type(crop_type)                         , intent(inout) :: crop_inst
    type(ch4_type)                          , intent(in)    :: ch4_inst
    type(dgvs_type)                         , intent(inout) :: dgvs_inst
    type(photosyns_type)                    , intent(in)    :: photosyns_inst
    type(soilhydrology_type)                , intent(in)    :: soilhydrology_inst
    type(energyflux_type)                   , intent(in)    :: energyflux_inst
    class(nutrient_competition_method_type) , intent(inout) :: nutrient_competition_method
    class(cnfire_method_type)               , intent(inout) :: cnfire_method
    logical                                 , intent(in)    :: dribble_crophrv_xsmrpool_2atm
    !
    ! !LOCAL VARIABLES:
    real(r8):: cn_decomp_pools(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
    real(r8):: p_decomp_cpool_loss(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_cascade_transitions) !potential C loss from one pool to another
    real(r8):: pmnf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_cascade_transitions) !potential mineral N flux, from one pool to another
    real(r8):: arepr(bounds%begp:bounds%endp) ! reproduction allocation coefficient (only used for use_crop)
    real(r8):: aroot(bounds%begp:bounds%endp) ! root allocation coefficient (only used for use_crop)
    integer :: begp,endp
    integer :: begc,endc

    integer :: dummy_to_make_pgi_happy
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    !real(r8) , intent(in)    :: rootfr_patch(bounds%begp:, 1:)          
    !integer  , intent(in)    :: altmax_lastyear_indx_col(bounds%begc:)  ! frost table depth (m)

    associate(                                                                    &
         crootfr_patch             => soilstate_inst%crootfr_patch              , & ! fraction of roots for carbon in each soil layer  (nlevgrnd)
         altmax_lastyear_indx_col  => canopystate_inst%altmax_lastyear_indx_col , & ! frost table depth (m)
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
    dummy_to_make_pgi_happy = ubound(filter_soilc, 1)
    call soilbiogeochem_carbonflux_inst%SetValues( &
         num_soilc, filter_soilc, 0._r8)
    if ( use_c13 ) then
       call c13_soilbiogeochem_carbonflux_inst%SetValues( &
            num_soilc, filter_soilc, 0._r8)
    end if
    if ( use_c14 ) then
       call c14_soilbiogeochem_carbonflux_inst%SetValues( &
            num_soilc, filter_soilc, 0._r8)
    end if

    call cnveg_carbonflux_inst%SetValues( &
         num_soilp, filter_soilp, 0._r8, &
         num_soilc, filter_soilc, 0._r8)
    if ( use_c13 ) then
       call c13_cnveg_carbonflux_inst%SetValues( &
            num_soilp, filter_soilp, 0._r8, &
            num_soilc, filter_soilc, 0._r8)
    end if
    if ( use_c14 ) then
       call c14_cnveg_carbonflux_inst%SetValues( &
            num_soilp, filter_soilp, 0._r8, &
            num_soilc, filter_soilc, 0._r8)
    end if

    call cnveg_nitrogenflux_inst%SetValues( &
         num_soilp, filter_soilp, 0._r8, &
         num_soilc, filter_soilc, 0._r8)

    call soilbiogeochem_nitrogenflux_inst%SetValues( &
         num_soilc, filter_soilc, 0._r8)

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
        call CNFreeLivingFixation( num_soilc, filter_soilc, &
             waterflux_inst, soilbiogeochem_nitrogenflux_inst)
        call t_stopf('CNFLivFixation')
    else
       call t_startf('CNFixation')
       call CNNFixation( num_soilc, filter_soilc, &
            cnveg_carbonflux_inst, soilbiogeochem_nitrogenflux_inst)
       call t_stopf('CNFixation')
    end if
  

    if (use_crop) then
       call CNNFert(bounds, num_soilc,filter_soilc, &
            cnveg_nitrogenflux_inst, soilbiogeochem_nitrogenflux_inst)

       if (.not. use_fun) then  ! if FUN is active, then soy fixation handled by FUN
          call  CNSoyfix (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               waterstate_inst, crop_inst, cnveg_state_inst, cnveg_nitrogenflux_inst , &
               soilbiogeochem_state_inst, soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)
       end if
    end if

    call t_startf('CNMResp')
    call CNMResp(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
         canopystate_inst, soilstate_inst, temperature_inst, photosyns_inst, &
         cnveg_carbonflux_inst, cnveg_nitrogenstate_inst)
    call t_stopf('CNMResp')

    !--------------------------------------------
    ! Soil Biogeochemistry
    !--------------------------------------------

    call t_startf('SoilBiogeochem')
    if (use_century_decomp) then
       call decomp_rate_constants_bgc(bounds, num_soilc, filter_soilc, &
            canopystate_inst, soilstate_inst, temperature_inst, ch4_inst, soilbiogeochem_carbonflux_inst)
    else
       call decomp_rate_constants_cn(bounds, num_soilc, filter_soilc, &
            canopystate_inst, soilstate_inst, temperature_inst, ch4_inst, soilbiogeochem_carbonflux_inst)
    end if

    ! calculate potential decomp rates and total immobilization demand (previously inlined in CNDecompAlloc)
    call SoilBiogeochemPotential (bounds, num_soilc, filter_soilc,                                                    &
         soilbiogeochem_state_inst, soilbiogeochem_carbonstate_inst, soilbiogeochem_carbonflux_inst,                  &
         soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst,                                         &
         cn_decomp_pools=cn_decomp_pools(begc:endc,1:nlevdecomp,1:ndecomp_pools), & 
         p_decomp_cpool_loss=p_decomp_cpool_loss(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions), &
         pmnf_decomp_cascade=pmnf_decomp_cascade(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions)) 

    ! calculate vertical profiles for distributing soil and litter C and N (previously subroutine decomp_vertprofiles called from CNDecompAlloc)
    call SoilBiogeochemVerticalProfile(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
         canopystate_inst, soilstate_inst,soilbiogeochem_state_inst)

    ! calculate nitrification and denitrification rates (previously subroutine nitrif_denitrif called from CNDecompAlloc)
    if (use_nitrif_denitrif) then 
       call SoilBiogeochemNitrifDenitrif(bounds, num_soilc, filter_soilc, &
            soilstate_inst, waterstate_inst, temperature_inst, ch4_inst, &
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
     if ( use_fun ) then
       call t_startf('CNPhenology_phase1')
       call CNPhenology (bounds, num_soilc, filter_soilc, num_soilp, &
            filter_soilp, num_pcropp, filter_pcropp, &
            doalb, waterstate_inst, temperature_inst, atm2lnd_inst, &
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

     call t_startf('calc_plant_nutrient_demand')
     call nutrient_competition_method%calc_plant_nutrient_demand ( &
         bounds, num_soilp, filter_soilp,                                 &
         photosyns_inst, crop_inst, canopystate_inst,                     &
         cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst, &
         c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst,            &
         cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst,               &
         soilbiogeochem_carbonflux_inst, soilbiogeochem_nitrogenstate_inst, &
         energyflux_inst, &
         aroot=aroot(begp:endp), arepr=arepr(begp:endp))

     ! get the column-averaged plant_ndemand (needed for following call to SoilBiogeochemCompetition)

     call p2c(bounds, num_soilc, filter_soilc,                    &
         cnveg_nitrogenflux_inst%plant_ndemand_patch(begp:endp), &
         soilbiogeochem_state_inst%plant_ndemand_col(begc:endc))
     call t_stopf('calc_plant_nutrient_demand')

     ! resolve plant/heterotroph competition for mineral N 
 
   
     call t_startf('soilbiogeochemcompetition')
     call SoilBiogeochemCompetition (bounds, num_soilc, filter_soilc,num_soilp, filter_soilp, waterstate_inst, &
                                     waterflux_inst,temperature_inst,soilstate_inst,cnveg_state_inst,          &
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
         bounds, num_soilp, filter_soilp, &
         cnveg_state_inst, crop_inst, canopystate_inst, &
         cnveg_carbonstate_inst, cnveg_carbonflux_inst, &
         c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst, &
         cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
         soilbiogeochem_nitrogenstate_inst, &
         aroot=aroot(begp:endp), &
         arepr=arepr(begp:endp), &
         fpg_col=soilbiogeochem_state_inst%fpg_col(begc:endc))
    call t_stopf('calc_plant_nutrient_competition')

    call t_stopf('CNDecompAlloc')

    !--------------------------------------------
    ! Calculate litter and soil decomposition rate
    !--------------------------------------------

    ! Calculation of actual immobilization and decomp rates, following
    ! resolution of plant/heterotroph  competition for mineral N (previously inlined in CNDecompAllocation in CNDecompMod)

    call t_startf('SoilBiogeochemDecomp')

    call SoilBiogeochemDecomp (bounds, num_soilc, filter_soilc,                                                       &
         soilbiogeochem_state_inst, soilbiogeochem_carbonstate_inst, soilbiogeochem_carbonflux_inst,                  &
         soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst,                                         &
         cn_decomp_pools=cn_decomp_pools(begc:endc,1:nlevdecomp,1:ndecomp_pools),                       & 
         p_decomp_cpool_loss=p_decomp_cpool_loss(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions), &
         pmnf_decomp_cascade=pmnf_decomp_cascade(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions)) 

    call t_stopf('SoilBiogeochemDecomp')

    !--------------------------------------------
    ! Phenology
    !--------------------------------------------

    ! CNphenology needs to be called after above calls, since it depends on current
    ! time-step fluxes to new growth on the lastlitterfall timestep in deciduous systems

    call t_startf('CNPhenology')

    if ( .not. use_fun ) then
       call CNPhenology (bounds, num_soilc, filter_soilc, num_soilp, &
            filter_soilp, num_pcropp, filter_pcropp, &
            doalb, waterstate_inst, temperature_inst, atm2lnd_inst, &
            crop_inst, canopystate_inst, soilstate_inst, dgvs_inst, &
            cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst, &
            cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
            c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst, &
            leaf_prof_patch=soilbiogeochem_state_inst%leaf_prof_patch(begp:endp,1:nlevdecomp_full), &
            froot_prof_patch=soilbiogeochem_state_inst%froot_prof_patch(begp:endp,1:nlevdecomp_full), &
            phase=1)
    end if
    call CNPhenology (bounds, num_soilc, filter_soilc, num_soilp, &
         filter_soilp, num_pcropp, filter_pcropp, &
         doalb, waterstate_inst, temperature_inst, atm2lnd_inst, &
         crop_inst, canopystate_inst, soilstate_inst, dgvs_inst, &
         cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst, &
         cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
         c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst, &
         leaf_prof_patch=soilbiogeochem_state_inst%leaf_prof_patch(begp:endp,1:nlevdecomp_full), &
         froot_prof_patch=soilbiogeochem_state_inst%froot_prof_patch(begp:endp,1:nlevdecomp_full), &
         phase=2)

    call t_stopf('CNPhenology')

    !--------------------------------------------
    ! Growth respiration
    !--------------------------------------------

    call t_startf('CNGResp')

    call CNGResp(num_soilp, filter_soilp,&
         cnveg_carbonflux_inst, canopystate_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst)  
         
    call t_stopf('CNGResp')

    !--------------------------------------------
    ! Dynamic Roots
    !--------------------------------------------

    if( use_dynroot ) then
        call t_startf('CNRootDyn')

        call CNRootDyn(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
             cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst,  &
             cnveg_state_inst, crop_inst,  soilstate_inst, soilbiogeochem_nitrogenstate_inst)

        call t_stopf('CNRootDyn')
     end if

    !--------------------------------------------
    ! CNUpdate0
    !--------------------------------------------

    call t_startf('CNUpdate0')

    call CStateUpdate0(num_soilp, filter_soilp, &
         cnveg_carbonflux_inst, cnveg_carbonstate_inst)

    if ( use_c13 ) then
       call CStateUpdate0(num_soilp, filter_soilp, &
            c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst)
    end if

    if ( use_c14 ) then
       call CStateUpdate0(num_soilp, filter_soilp, &
            c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst)
    end if

    call t_stopf('CNUpdate0')

    if ( use_nguardrail ) then
       call t_startf('CNPrecisionControl')
       call CNPrecisionControl(bounds, num_soilp, filter_soilp, &
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

       call CIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp,              &
            soilbiogeochem_state_inst,                                               &
            soilbiogeochem_carbonflux_inst,  soilbiogeochem_carbonstate_inst,        &
            cnveg_carbonflux_inst, cnveg_carbonstate_inst,                           &
            c13_soilbiogeochem_carbonflux_inst, c13_soilbiogeochem_carbonstate_inst, &
            c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst,                   &
            isotope='c13')
    end if
    if ( use_c14 ) then
       call CIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp,              &
            soilbiogeochem_state_inst,                                               &
            soilbiogeochem_carbonflux_inst,  soilbiogeochem_carbonstate_inst,        &
            cnveg_carbonflux_inst, cnveg_carbonstate_inst,                           &
            c14_soilbiogeochem_carbonflux_inst, c14_soilbiogeochem_carbonstate_inst, &
            c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst,                   &
            isotope='c14')
    end if

    ! Update all prognostic carbon state variables (except for gap-phase mortality and fire fluxes)
    call CStateUpdate1( num_soilc, filter_soilc, num_soilp, filter_soilp, &
         crop_inst, cnveg_carbonflux_inst, cnveg_carbonstate_inst, &
         soilbiogeochem_carbonflux_inst, dribble_crophrv_xsmrpool_2atm)
    if ( use_c13 ) then
       call CStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            crop_inst, c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst, &
            c13_soilbiogeochem_carbonflux_inst, dribble_crophrv_xsmrpool_2atm)
    end if
    if ( use_c14 ) then
       call CStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            crop_inst, c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst, &
            c14_soilbiogeochem_carbonflux_inst, dribble_crophrv_xsmrpool_2atm)
    end if

    ! Update all prognostic nitrogen state variables (except for gap-phase mortality and fire fluxes)
    call NStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
         cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)

    call t_stopf('CNUpdate1')

    if ( use_nguardrail ) then
       call t_startf('CNPrecisionControl')
       call CNPrecisionControl(bounds, num_soilp, filter_soilp, &
            cnveg_carbonstate_inst, c13_cnveg_carbonstate_inst, &
            c14_cnveg_carbonstate_inst, cnveg_nitrogenstate_inst)
       call t_stopf('CNPrecisionControl')
    end if

    call t_startf('SoilBiogeochemStateUpdate1')
    call SoilBiogeochemNStateUpdate1(num_soilc, filter_soilc,  &
         soilbiogeochem_state_inst, soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst)
    call t_stopf('SoilBiogeochemStateUpdate1')


    !--------------------------------------------
    ! Calculate vertical mixing of soil and litter pools
    !--------------------------------------------

    call t_startf('SoilBiogeochemLittVertTransp')

    call SoilBiogeochemLittVertTransp(bounds, num_soilc, filter_soilc,            &
         canopystate_inst, soilbiogeochem_state_inst,                             &
         soilbiogeochem_carbonstate_inst, soilbiogeochem_carbonflux_inst,         &
         c13_soilbiogeochem_carbonstate_inst, c13_soilbiogeochem_carbonflux_inst, &
         c14_soilbiogeochem_carbonstate_inst, c14_soilbiogeochem_carbonflux_inst, &
         soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)

    call t_stopf('SoilBiogeochemLittVertTransp')

    !--------------------------------------------
    ! Calculate the gap mortality carbon and nitrogen fluxes
    !--------------------------------------------

    call t_startf('CNGapMortality')

    call CNGapMortality (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp,                                &
         dgvs_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst,                                             &
         cnveg_carbonflux_inst, cnveg_nitrogenflux_inst,  canopystate_inst,                                       &   
         !cnveg_carbonflux_inst, cnveg_nitrogenflux_inst,                                                         &  
         leaf_prof_patch=soilbiogeochem_state_inst%leaf_prof_patch(begp:endp, 1:nlevdecomp_full),   &
         froot_prof_patch=soilbiogeochem_state_inst%froot_prof_patch(begp:endp, 1:nlevdecomp_full), & 
         croot_prof_patch=soilbiogeochem_state_inst%croot_prof_patch(begp:endp, 1:nlevdecomp_full), &
         stem_prof_patch=soilbiogeochem_state_inst%stem_prof_patch(begp:endp, 1:nlevdecomp_full))   

    call t_stopf('CNGapMortality')

    !--------------------------------------------
    ! Update2 (gap mortality)
    !--------------------------------------------

    call t_startf('CNUpdate2')

    ! Set the carbon isotopic fluxes for gap mortality
    if ( use_c13 ) then
       call CIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp,               &
            soilbiogeochem_state_inst, cnveg_carbonflux_inst, cnveg_carbonstate_inst, &
            iso_cnveg_carbonflux_inst=c13_cnveg_carbonflux_inst,                      &
            iso_cnveg_carbonstate_inst=c13_cnveg_carbonstate_inst,                    &
            isotope='c13')
    end if
    if ( use_c14 ) then
       call CIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp,               &
            soilbiogeochem_state_inst, cnveg_carbonflux_inst, cnveg_carbonstate_inst, &
            iso_cnveg_carbonflux_inst=c14_cnveg_carbonflux_inst,                      &
            iso_cnveg_carbonstate_inst=c14_cnveg_carbonstate_inst,                    &
            isotope='c14')
    end if

    ! Update all the prognostic carbon state variables affected by gap-phase mortality fluxes
    call CStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
         cnveg_carbonflux_inst, cnveg_carbonstate_inst, soilbiogeochem_carbonstate_inst)
    if ( use_c13 ) then
       call CStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst, c13_soilbiogeochem_carbonstate_inst)
    end if
    if ( use_c14 ) then
       call CStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst, c14_soilbiogeochem_carbonstate_inst)
    end if

    ! Update all the prognostic nitrogen state variables affected by gap-phase mortality fluxes
    call NStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
         cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, soilbiogeochem_nitrogenstate_inst)

    !--------------------------------------------
    ! Update2h (harvest)
    !--------------------------------------------

    ! Set harvest mortality routine 
    if (get_do_harvest()) then
       call CNHarvest(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            soilbiogeochem_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, &
            cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
    end if

    if ( use_c13 ) then
       call CIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp,   &
            soilbiogeochem_state_inst,                                     &
            cnveg_carbonflux_inst, cnveg_carbonstate_inst,                 &
            c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst,         &                         
            isotope='c13')
    end if
    if ( use_c14 ) then
       call CIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            soilbiogeochem_state_inst,                                     &
            cnveg_carbonflux_inst, cnveg_carbonstate_inst,                 &
            c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst,         &                         
            isotope='c14')
    end if

    call CStateUpdate2h( num_soilc, filter_soilc,  num_soilp, filter_soilp, &
         cnveg_carbonflux_inst, cnveg_carbonstate_inst, soilbiogeochem_carbonstate_inst)
    if ( use_c13 ) then
       call CStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst, c13_soilbiogeochem_carbonstate_inst)
    end if
    if ( use_c14 ) then
       call CStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst, c14_soilbiogeochem_carbonstate_inst)
    end if

    call NStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
         cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, soilbiogeochem_nitrogenstate_inst)
    call t_stopf('CNUpdate2')

    if ( use_nguardrail ) then
       call t_startf('CNPrecisionControl')
       call CNPrecisionControl(bounds, num_soilp, filter_soilp, &
            cnveg_carbonstate_inst, c13_cnveg_carbonstate_inst, &
            c14_cnveg_carbonstate_inst, cnveg_nitrogenstate_inst)
       call t_stopf('CNPrecisionControl')
    end if
    !--------------------------------------------
    ! Calculate loss fluxes from wood products pools
    ! and update product pool state variables
    !--------------------------------------------

    call t_startf('CNWoodProducts')
    call c_products_inst%UpdateProducts(bounds, &
         num_soilp, filter_soilp, &
         dwt_wood_product_gain_patch = cnveg_carbonflux_inst%dwt_wood_productc_gain_patch(begp:endp), &
         wood_harvest_patch = cnveg_carbonflux_inst%wood_harvestc_patch(begp:endp), &
         dwt_crop_product_gain_patch = cnveg_carbonflux_inst%dwt_crop_productc_gain_patch(begp:endp), &
         grain_to_cropprod_patch = cnveg_carbonflux_inst%grainc_to_cropprodc_patch(begp:endp))
    call t_stopf('CNWoodProducts')

    if (use_c13) then
       call c13_products_inst%UpdateProducts(bounds, &
            num_soilp, filter_soilp, &
            dwt_wood_product_gain_patch = c13_cnveg_carbonflux_inst%dwt_wood_productc_gain_patch(begp:endp), &
            wood_harvest_patch = c13_cnveg_carbonflux_inst%wood_harvestc_patch(begp:endp), &
            dwt_crop_product_gain_patch = c13_cnveg_carbonflux_inst%dwt_crop_productc_gain_patch(begp:endp), &
            grain_to_cropprod_patch = c13_cnveg_carbonflux_inst%grainc_to_cropprodc_patch(begp:endp))
    end if

    if (use_c14) then
       call c14_products_inst%UpdateProducts(bounds, &
            num_soilp, filter_soilp, &
            dwt_wood_product_gain_patch = c14_cnveg_carbonflux_inst%dwt_wood_productc_gain_patch(begp:endp), &
            wood_harvest_patch = c14_cnveg_carbonflux_inst%wood_harvestc_patch(begp:endp), &
            dwt_crop_product_gain_patch = c14_cnveg_carbonflux_inst%dwt_crop_productc_gain_patch(begp:endp), &
            grain_to_cropprod_patch = c14_cnveg_carbonflux_inst%grainc_to_cropprodc_patch(begp:endp))
    end if

    call n_products_inst%UpdateProducts(bounds, &
         num_soilp, filter_soilp, &
         dwt_wood_product_gain_patch = cnveg_nitrogenflux_inst%dwt_wood_productn_gain_patch(begp:endp), &
         wood_harvest_patch = cnveg_nitrogenflux_inst%wood_harvestn_patch(begp:endp), &
         dwt_crop_product_gain_patch = cnveg_nitrogenflux_inst%dwt_crop_productn_gain_patch(begp:endp), &
         grain_to_cropprod_patch = cnveg_nitrogenflux_inst%grainn_to_cropprodn_patch(begp:endp))

    !--------------------------------------------
    ! Calculate fire area and fluxes
    !--------------------------------------------

    call t_startf('CNFire')
    call cnfire_method%CNFireArea(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
         atm2lnd_inst, energyflux_inst, soilhydrology_inst, waterstate_inst, &
         cnveg_state_inst, cnveg_carbonstate_inst, &
         totlitc_col=soilbiogeochem_carbonstate_inst%totlitc_col(begc:endc), &
         decomp_cpools_vr_col=soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools), &
         t_soi17cm_col=temperature_inst%t_soi17cm_col(begc:endc))

    call cnfire_method%CNFireFluxes(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp,                                      &
         dgvs_inst, cnveg_state_inst,                                                                                              &
         cnveg_carbonstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst,                         &
         leaf_prof_patch=soilbiogeochem_state_inst%leaf_prof_patch(begp:endp, 1:nlevdecomp_full),                                  &
         froot_prof_patch=soilbiogeochem_state_inst%froot_prof_patch(begp:endp, 1:nlevdecomp_full),                                &
         croot_prof_patch=soilbiogeochem_state_inst%croot_prof_patch(begp:endp, 1:nlevdecomp_full),                                &
         stem_prof_patch=soilbiogeochem_state_inst%stem_prof_patch(begp:endp, 1:nlevdecomp_full),                                  &
         totsomc_col=soilbiogeochem_carbonstate_inst%totsomc_col(begc:endc),                                                       &
         decomp_cpools_vr_col=soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools),   &
         decomp_npools_vr_col=soilbiogeochem_nitrogenstate_inst%decomp_npools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools), &
         somc_fire_col=soilbiogeochem_carbonflux_inst%somc_fire_col(begc:endc))
    call t_stopf('CNFire')


    !--------------------------------------------
    ! Update3
    !--------------------------------------------

    call t_startf('CNUpdate3')
    if ( use_c13 ) then
       call CIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            soilbiogeochem_state_inst , soilbiogeochem_carbonstate_inst,         &
            cnveg_carbonflux_inst, cnveg_carbonstate_inst,                       &
            c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst,               &
            c13_soilbiogeochem_carbonstate_inst, &
            isotope='c13')
    end if
    if ( use_c14 ) then
       call CIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            soilbiogeochem_state_inst , soilbiogeochem_carbonstate_inst,         &
            cnveg_carbonflux_inst, cnveg_carbonstate_inst,                       &
            c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst,               &
            c14_soilbiogeochem_carbonstate_inst, &
            isotope='c14')
    end if

    call CStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
         cnveg_carbonflux_inst, cnveg_carbonstate_inst, soilbiogeochem_carbonstate_inst)

    if ( use_c13 ) then
       call CStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst, c13_soilbiogeochem_carbonstate_inst)
    end if

    if ( use_c14 ) then
       call CStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst, c14_soilbiogeochem_carbonstate_inst)

       call C14Decay(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            c14_cnveg_carbonstate_inst, c14_soilbiogeochem_carbonstate_inst)
    end if
    call t_stopf('CNUpdate3')

    if ( use_nguardrail ) then
       call t_startf('CNPrecisionControl')
       call CNPrecisionControl(bounds, num_soilp, filter_soilp, &
            cnveg_carbonstate_inst, c13_cnveg_carbonstate_inst, &
            c14_cnveg_carbonstate_inst, cnveg_nitrogenstate_inst)
       call t_stopf('CNPrecisionControl')
    end if

    end associate

  end subroutine CNDriverNoLeaching
  
  !-----------------------------------------------------------------------
  subroutine CNDriverLeaching(bounds, &
       num_soilc, filter_soilc, num_soilp, filter_soilp, &
       waterstate_inst, waterflux_inst, &
       cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, &
       soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst)
    !
    ! !DESCRIPTION:
    ! Update the nitrogen leaching rate as a function of soluble mineral N and total soil water outflow.
    ! Also update nitrogen state variables         
    !
    ! !USES:
    use SoilBiogeochemNLeachingMod, only: SoilBiogeochemNLeaching
    use CNNStateUpdate3Mod   , only: NStateUpdate3
    !
    ! !ARGUMENTS:
    type(bounds_type)                       , intent(in)    :: bounds  
    integer                                 , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                                 , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    type(waterstate_type)                   , intent(in)    :: waterstate_inst
    type(waterflux_type)                    , intent(inout)    :: waterflux_inst
    type(cnveg_nitrogenflux_type)           , intent(inout) :: cnveg_nitrogenflux_inst
    type(cnveg_nitrogenstate_type)          , intent(inout) :: cnveg_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    !-----------------------------------------------------------------------
  
    ! Mineral nitrogen dynamics (deposition, fixation, leaching)
    
    call t_startf('SoilBiogeochemNLeaching')
    call SoilBiogeochemNLeaching(bounds, num_soilc, filter_soilc, &
         waterstate_inst, waterflux_inst, soilbiogeochem_nitrogenstate_inst, &
         soilbiogeochem_nitrogenflux_inst)
    call t_stopf('SoilBiogeochemNLeaching')

    ! Nitrogen state variable update, mortality fluxes.

    call t_startf('NUpdate3')

    call NstateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
         cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, &
         soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst)

    call t_stopf('NUpdate3')

  end subroutine CNDriverLeaching

  !-----------------------------------------------------------------------
  subroutine CNDriverSummarizeStates(bounds, num_allc, filter_allc, &
       num_soilc, filter_soilc, num_soilp, filter_soilp, &
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
    integer                                 , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                                 , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_soilp(:)   ! filter for soil patches
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
    ! soilbiogeochem carbon/nitrogen state summary
    ! ----------------------------------------------

    call soilbiogeochem_carbonstate_inst%summary(bounds, num_allc, filter_allc)
    if ( use_c13 ) then
       call c13_soilbiogeochem_carbonstate_inst%summary(bounds, num_allc, filter_allc)
    end if
    if ( use_c14 ) then
       call c14_soilbiogeochem_carbonstate_inst%summary(bounds, num_allc, filter_allc)
    end if
    call soilbiogeochem_nitrogenstate_inst%summary(bounds, num_allc, filter_allc)

    ! ----------------------------------------------
    ! cnveg carbon/nitrogen state summary
    ! ----------------------------------------------

    call cnveg_carbonstate_inst%Summary(bounds, num_allc, filter_allc, &
         num_soilc, filter_soilc, num_soilp, filter_soilp, &
         soilbiogeochem_cwdc_col=soilbiogeochem_carbonstate_inst%cwdc_col(begc:endc), &
         soilbiogeochem_totlitc_col=soilbiogeochem_carbonstate_inst%totlitc_col(begc:endc), &
         soilbiogeochem_totsomc_col=soilbiogeochem_carbonstate_inst%totsomc_col(begc:endc), &
         soilbiogeochem_ctrunc_col=soilbiogeochem_carbonstate_inst%ctrunc_col(begc:endc))

    if ( use_c13 ) then
       call c13_cnveg_carbonstate_inst%Summary(bounds, num_allc, filter_allc, &
            num_soilc, filter_soilc, num_soilp, filter_soilp, &
            soilbiogeochem_cwdc_col=c13_soilbiogeochem_carbonstate_inst%cwdc_col(begc:endc), &
            soilbiogeochem_totlitc_col=c13_soilbiogeochem_carbonstate_inst%totlitc_col(begc:endc), &
            soilbiogeochem_totsomc_col=c13_soilbiogeochem_carbonstate_inst%totsomc_col(begc:endc), &
            soilbiogeochem_ctrunc_col=c13_soilbiogeochem_carbonstate_inst%ctrunc_col(begc:endc))
    end if

    if ( use_c14 ) then
       call c14_cnveg_carbonstate_inst%Summary(bounds, num_allc, filter_allc, &
            num_soilc, filter_soilc, num_soilp, filter_soilp, &
            soilbiogeochem_cwdc_col=c14_soilbiogeochem_carbonstate_inst%cwdc_col(begc:endc), &
            soilbiogeochem_totlitc_col=c14_soilbiogeochem_carbonstate_inst%totlitc_col(begc:endc), &
            soilbiogeochem_totsomc_col=c14_soilbiogeochem_carbonstate_inst%totsomc_col(begc:endc), &
            soilbiogeochem_ctrunc_col=c14_soilbiogeochem_carbonstate_inst%ctrunc_col(begc:endc))
    end if

    call cnveg_nitrogenstate_inst%Summary(bounds, num_allc, filter_allc, &
         num_soilc, filter_soilc, num_soilp, filter_soilp, &
         soilbiogeochem_nitrogenstate_inst)

    call t_stopf('CNsum')

  end subroutine CNDriverSummarizeStates

  !-----------------------------------------------------------------------
  subroutine CNDriverSummarizeFluxes(bounds, &
       num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnveg_carbonflux_inst, c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst, &
       cnveg_nitrogenflux_inst, &
       c_products_inst, c13_products_inst, c14_products_inst, &
       soilbiogeochem_carbonflux_inst, &
       c13_soilbiogeochem_carbonflux_inst, &
       c14_soilbiogeochem_carbonflux_inst, &
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
    integer                                 , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                                 , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_soilp(:)   ! filter for soil patches
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
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: begc,endc
    integer :: begg,endg

    character(len=*), parameter :: subname = 'CNDriverSummarizeFluxes'
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg = bounds%endg

    call t_startf('CNsum')

    ! ----------------------------------------------
    ! soilbiogeochem carbon/nitrogen flux summary
    ! ----------------------------------------------

    call soilbiogeochem_carbonflux_inst%Summary(bounds, num_soilc, filter_soilc)
    if ( use_c13 ) then
       call c13_soilbiogeochem_carbonflux_inst%Summary(bounds, num_soilc, filter_soilc)
    end if
    if ( use_c14 ) then
       call c14_soilbiogeochem_carbonflux_inst%Summary(bounds, num_soilc, filter_soilc)
    end if
    call soilbiogeochem_nitrogenflux_inst%Summary(bounds, num_soilc, filter_soilc)

    ! ----------------------------------------------
    ! cnveg carbon/nitrogen flux summary
    ! ----------------------------------------------

    call cnveg_carbonflux_inst%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
         isotope='bulk', &
         soilbiogeochem_hr_col=soilbiogeochem_carbonflux_inst%hr_col(begc:endc), &
         soilbiogeochem_lithr_col=soilbiogeochem_carbonflux_inst%lithr_col(begc:endc), &  
         soilbiogeochem_decomp_cascade_ctransfer_col=&
         soilbiogeochem_carbonflux_inst%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions), &
         product_closs_grc=c_products_inst%product_loss_grc(begg:endg))

    if ( use_c13 ) then
       call c13_cnveg_carbonflux_inst%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            isotope='c13', &
            soilbiogeochem_hr_col=c13_soilbiogeochem_carbonflux_inst%hr_col(begc:endc), &
            soilbiogeochem_lithr_col=c13_soilbiogeochem_carbonflux_inst%lithr_col(begc:endc), &  
            soilbiogeochem_decomp_cascade_ctransfer_col=&
            c13_soilbiogeochem_carbonflux_inst%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions), &
            product_closs_grc=c13_products_inst%product_loss_grc(begg:endg))
    end if

    if ( use_c14 ) then
       call c14_cnveg_carbonflux_inst%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            isotope='c14', &
            soilbiogeochem_hr_col=c14_soilbiogeochem_carbonflux_inst%hr_col(begc:endc), &
            soilbiogeochem_lithr_col=c14_soilbiogeochem_carbonflux_inst%lithr_col(begc:endc), &  
            soilbiogeochem_decomp_cascade_ctransfer_col=&
            c14_soilbiogeochem_carbonflux_inst%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions), &
            product_closs_grc=c14_products_inst%product_loss_grc(begg:endg))
    end if

    call cnveg_nitrogenflux_inst%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)

    call t_stopf('CNsum')

  end subroutine CNDriverSummarizeFluxes

end  module CNDriverMod
