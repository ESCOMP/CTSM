module ctsm_FatesBiogeochemDynMod

! This module creates a pathway to call the belowground biogeochemistry code as driven by the fates vegetation model 
! but bypassing the aboveground CN vegetation code.  It is modeled after the ctsm_CNDriverMod in its call sequence and 
! functionality.

  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use ctsm_VarCtl                      , only : use_c13, use_c14, use_fates
  use ctsm_Decomp                       , only : bounds_type
  use perf_mod                        , only : t_startf, t_stopf
  use shr_log_mod                     , only : errMsg => shr_log_errMsg
  use ctsm_AbortUtils                      , only : endrun
  use ctsm_VarCtl                      , only : use_century_decomp, use_nitrif_denitrif
  use ctsm_CNVegCarbonStateType	      , only : cnveg_carbonstate_type
  use ctsm_CNVegCarbonFluxType	      , only : cnveg_carbonflux_type
  use ctsm_SoilBiogeochemStateType         , only : soilbiogeochem_state_type
  use ctsm_SoilBiogeochemCarbonStateType   , only : soilbiogeochem_carbonstate_type
  use ctsm_SoilBiogeochemCarbonFluxType    , only : soilbiogeochem_carbonflux_type
  use ctsm_SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
  use ctsm_SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
  use ctsm_CanopyStateType                 , only : canopystate_type
  use ctsm_SoilStateType                   , only : soilstate_type
  use ctsm_SoilHydrologyType               , only : soilhydrology_type
  use ctsm_TemperatureType                 , only : temperature_type
  use ctsm_WaterFluxBulkType                   , only : waterfluxbulk_type
  use ctsm_ActiveLayer                  , only : active_layer_type
  use ctsm_Atm2LndType                     , only : atm2lnd_type
  use ctsm_SoilStateType                   , only : soilstate_type
  use ctsm_Methane                          , only : ch4_type


  ! public :: ctsm_FatesBiogeochemDynInit         ! BGC dynamics: initialization
  public :: ctsm_FatesBiogeochemDyn             ! BGC Dynamics
  public :: ctsm_FatesBiogeochemDynSummary      ! BGC dynamics: summary

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains


  !-----------------------------------------------------------------------
  subroutine ctsm_FatesBiogeochemDyn(bounds,                                                    &
       num_soilc, filter_soilc, num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb, &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst,                                      &
       soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst,                    &
       soilbiogeochem_state_inst,                                                          &
       soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst,                &
       c13_soilbiogeochem_carbonstate_inst, c13_soilbiogeochem_carbonflux_inst,            &
       c14_soilbiogeochem_carbonstate_inst, c14_soilbiogeochem_carbonflux_inst,            &
       active_layer_inst, atm2lnd_inst, waterfluxbulk_inst,                                &
       canopystate_inst, soilstate_inst, temperature_inst, crop_inst, ch4_inst)
    !
    ! !DESCRIPTION:

    !
    ! !USES:
    use ctsm_VarPar                        , only: nlevgrnd, nlevdecomp_full 
    use ctsm_VarPar                        , only: nlevdecomp, ndecomp_cascade_transitions, ndecomp_pools
    use ctsm_SubgridAve                     , only: p2c
    use ctsm_CropType                          , only: crop_type
    use ctsm_CNNDynamicsMod                    , only: CNNDeposition,CNNFixation, CNNFert, CNSoyfix
    use ctsm_CNMaintRespMod                        , only: ctsm_CNMaintResp
    use ctsm_CNPhenologyMod                    , only: ctsm_CNPhenology
    use ctsm_CNGrowthRespMod                        , only: ctsm_CNGrowthResp
    use ctsm_CNCIsoFluxes                     , only: CIsoFlux1, CIsoFlux2, CIsoFlux2h, CIsoFlux3
    use ctsm_CNC14DecayMod                     , only: C14Decay
    use ctsm_CNCStateUpdate1Mod                , only: CStateUpdate1,CStateUpdate0
    use ctsm_CNCStateUpdate2Mod                , only: CStateUpdate2, CStateUpdate2h
    use ctsm_CNCStateUpdate3Mod                , only: CStateUpdate3
    use ctsm_CNNStateUpdate1Mod                , only: NStateUpdate1
    use ctsm_CNNStateUpdate2Mod                , only: NStateUpdate2, NStateUpdate2h
    use ctsm_CNGapMortalityMod                 , only: ctsm_CNGapMortality
    use dynHarvestMod                     , only: CNHarvest
    use ctsm_SoilBiogeochemDecompCascadeBGC , only: decomp_rate_constants_bgc
    use ctsm_SoilBiogeochemDecompCascadeCN  , only: decomp_rate_constants_cn
    use ctsm_SoilBiogeochemCompetition      , only: SoilBiogeochemCompetition
    use ctsm_SoilBiogeochemDecomp           , only: SoilBiogeochemDecomp
    use ctsm_SoilBiogeochemLittVertTransp   , only: SoilBiogeochemLittVertTransp
    use ctsm_SoilBiogeochemPotential        , only: SoilBiogeochemPotential 
    use ctsm_SoilBiogeochemVerticalProfile  , only: SoilBiogeochemVerticalProfile
    use ctsm_SoilBiogeochemNitrifDenitrif   , only: SoilBiogeochemNitrifDenitrif
    use ctsm_SoilBiogeochemNStateUpdate1    , only: SoilBiogeochemNStateUpdate1
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
    type(cnveg_carbonflux_type)             , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)            , intent(inout) :: cnveg_carbonstate_inst
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
    type(waterfluxbulk_type)                    , intent(in)    :: waterfluxbulk_inst
    type(canopystate_type)                  , intent(in)    :: canopystate_inst
    type(soilstate_type)                    , intent(in)    :: soilstate_inst
    type(temperature_type)                  , intent(inout) :: temperature_inst
    type(crop_type)                         , intent(in)    :: crop_inst
    type(ch4_type)                          , intent(in)    :: ch4_inst
    !
    ! !LOCAL VARIABLES:
    real(r8):: cn_decomp_pools(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
    real(r8):: p_decomp_cpool_loss(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_cascade_transitions) !potential C loss from one pool to another
    real(r8):: pmnf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_cascade_transitions) !potential mineral N flux, from one pool to another
    real(r8):: arepr(bounds%begp:bounds%endp) ! reproduction allocation coefficient (only used for crop_prog)
    real(r8):: aroot(bounds%begp:bounds%endp) ! root allocation coefficient (only used for crop_prog)
    integer :: begp,endp
    integer :: begc,endc
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
    
    call t_startf('BGCZero')

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

    call t_stopf('BGCZero')

    ! --------------------------------------------------
    ! Nitrogen Deposition, Fixation and Respiration
    ! --------------------------------------------------

    ! call t_startf('CNDeposition')
    ! call CNNDeposition(bounds, &
    !      atm2lnd_inst, soilbiogeochem_nitrogenflux_inst)
    ! call t_stopf('CNDeposition')


    ! if (crop_prog) then
    !    call CNNFert(bounds, num_soilc,filter_soilc, &
    !         cnveg_nitrogenflux_inst, soilbiogeochem_nitrogenflux_inst)

    !    call  CNSoyfix (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
    !         waterstate_inst, crop_inst, cnveg_state_inst, cnveg_nitrogenflux_inst , &
    !         soilbiogeochem_state_inst, soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)
    ! end if

    !--------------------------------------------
    ! Soil Biogeochemistry
    !--------------------------------------------

    if (use_century_decomp) then
       call decomp_rate_constants_bgc(bounds, num_soilc, filter_soilc, &
            soilstate_inst, temperature_inst, ch4_inst, soilbiogeochem_carbonflux_inst)
    else
       call decomp_rate_constants_cn(bounds, num_soilc, filter_soilc, &
            soilstate_inst, temperature_inst, ch4_inst, soilbiogeochem_carbonflux_inst)
    end if

    ! calculate potential decomp rates and total immobilization demand (previously inlined in CNDecompAlloc)
    call SoilBiogeochemPotential (bounds, num_soilc, filter_soilc,                                                    &
         soilbiogeochem_state_inst, soilbiogeochem_carbonstate_inst, soilbiogeochem_carbonflux_inst,                  &
         soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst,                                         &
         cn_decomp_pools=cn_decomp_pools(begc:endc,1:nlevdecomp,1:ndecomp_pools), & 
         p_decomp_cpool_loss=p_decomp_cpool_loss(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions), &
         pmnf_decomp_cascade=pmnf_decomp_cascade(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions)) 


    !--------------------------------------------
    ! Resolve the competition between plants and soil heterotrophs 
    ! for available soil mineral N resource 
    !--------------------------------------------
    ! will add this back in when integrtating hte nutirent cycles


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
    ! Update1
    !--------------------------------------------

    call t_startf('BNGCUpdate1')


    ! Update all prognostic carbon state variables (except for gap-phase mortality and fire fluxes)
    call CStateUpdate1( num_soilc, filter_soilc, num_soilp, filter_soilp, &
         crop_inst, cnveg_carbonflux_inst, cnveg_carbonstate_inst, &
         soilbiogeochem_carbonflux_inst, dribble_crophrv_xsmrpool_2atm=.False.)

    call t_stopf('BNGCUpdate1')

    !--------------------------------------------
    ! Calculate vertical mixing of soil and litter pools
    !--------------------------------------------

    call t_startf('SoilBiogeochemLittVertTransp')

    call SoilBiogeochemLittVertTransp(bounds, num_soilc, filter_soilc,            &
         active_layer_inst, soilbiogeochem_state_inst,                            &
         soilbiogeochem_carbonstate_inst, soilbiogeochem_carbonflux_inst,         &
         c13_soilbiogeochem_carbonstate_inst, c13_soilbiogeochem_carbonflux_inst, &
         c14_soilbiogeochem_carbonstate_inst, c14_soilbiogeochem_carbonflux_inst, &
         soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)

    call t_stopf('SoilBiogeochemLittVertTransp')

    end associate

  end subroutine ctsm_FatesBiogeochemDyn


  !-----------------------------------------------------------------------
  subroutine ctsm_FatesBiogeochemDynSummary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst, &
       c13_soilbiogeochem_carbonflux_inst, c13_soilbiogeochem_carbonstate_inst, &
       c14_soilbiogeochem_carbonflux_inst, c14_soilbiogeochem_carbonstate_inst, &
       soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst, &
       clm_fates, nc)
    !
    ! !DESCRIPTION:
    ! Call to all CN and SoilBiogeochem summary routines
    ! also aggregate production and decomposition fluxes to whole-ecosystem balance fluxes
    !
    ! !USES:
    use ctsm_VarPar                        , only: ndecomp_cascade_transitions
    use ctsm_CNPrecisionControlMod             , only: ctsm_CNPrecisionControl
    use ctsm_SoilBiogeochemPrecisionControl , only: SoilBiogeochemPrecisionControl
    use ctsm_FatesInterface              , only: hlm_fates_interface_type
    !
    ! !ARGUMENTS:
    type(bounds_type)                       , intent(in)    :: bounds  
    integer                                 , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                                 , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: c13_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c13_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: c14_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c14_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(hlm_fates_interface_type)          , intent(inout) :: clm_fates
    integer                                 , intent(in)    :: nc  ! thread index
    !
    ! !LOCAL VARIABLES:
    integer :: begc,endc
    !-----------------------------------------------------------------------
  
    begc = bounds%begc; endc= bounds%endc

    ! Call to all summary routines

    call t_startf('BGCsum')

    ! Set controls on very low values in critical state variables 

    call SoilBiogeochemPrecisionControl(num_soilc, filter_soilc,  &
         soilbiogeochem_carbonstate_inst, c13_soilbiogeochem_carbonstate_inst, &
         c14_soilbiogeochem_carbonstate_inst,soilbiogeochem_nitrogenstate_inst)

    ! Note - all summary updates to cnveg_carbonstate_inst and cnveg_carbonflux_inst are done in 
    ! soilbiogeochem_carbonstate_inst%summary and CNVeg_carbonstate_inst%summary

    ! ----------------------------------------------
    ! soilbiogeochem carbon/nitrogen state summary
    ! ----------------------------------------------

    call soilbiogeochem_carbonstate_inst%summary(bounds, num_soilc, filter_soilc)
    if ( use_c13 ) then
       call c13_soilbiogeochem_carbonstate_inst%summary(bounds, num_soilc, filter_soilc)
    end if
    if ( use_c14 ) then
       call c14_soilbiogeochem_carbonstate_inst%summary(bounds, num_soilc, filter_soilc)
    end if
    ! call soilbiogeochem_nitrogenstate_inst%summary(bounds, num_soilc, filter_soilc)

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
    ! call soilbiogeochem_nitrogenflux_inst%Summary(bounds, num_soilc, filter_soilc)


    ! -----------------------------------------------------------------------------------
    ! fates veg carbon state and flux summary, Nitrogen (TBD) and Balance Checks
    ! -----------------------------------------------------------------------------------
    ! ----------------------------------------------
    ! fates veg nitrogen flux summary
    ! ----------------------------------------------
    ! ----------------------------------------------
    ! calculate balance checks on entire carbon cycle (FATES + BGC)
    ! ----------------------------------------------

    call clm_fates%wrap_bgc_summary(nc, soilbiogeochem_carbonflux_inst, &
                                        soilbiogeochem_carbonstate_inst)

    call t_stopf('BGCsum')

  end subroutine ctsm_FatesBiogeochemDynSummary

end  module ctsm_FatesBiogeochemDynMod
