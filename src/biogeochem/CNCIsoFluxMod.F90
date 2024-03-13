module CNCIsoFluxMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for carbon isotopic flux variable update, non-mortality fluxes.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varpar                         , only : ndecomp_cascade_transitions, nlevdecomp, ndecomp_pools
  use clm_varpar                         , only : i_litr_min, i_litr_max, i_met_lit
  use abortutils                         , only : endrun
  use pftconMod                          , only : pftcon
  use CNVegCarbonStateType               , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType                , only : cnveg_carbonflux_type
  use CropReprPoolsMod                   , only : nrepr, repr_grain_min, repr_grain_max, repr_structure_min, repr_structure_max
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use SoilBiogeochemStateType            , only : soilbiogeochem_state_type
  use SoilBiogeochemCarbonStateType      , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemCarbonFluxType       , only : soilbiogeochem_carbonflux_type
  use ColumnType                         , only : col                
  use PatchType                          , only : patch                
  use clm_varctl                         , only : use_crop
  use clm_varctl                         , only : use_grainproduct
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: CIsoFlux1
  public  :: CIsoFlux2
  public  :: CIsoFlux2h
  public  :: CIsoFlux2g
  public  :: CIsoFlux3
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: CNCIsoLitterToColumn
  private :: CNCIsoGapPftToColumn
  private :: CNCIsoHarvestPftToColumn
  private :: CNCIsoGrossUnrepPftToColumn
  private :: CIsoFluxCalc1d
  private :: CIsoFluxCalc2dFlux
  private :: CIsoFluxCalc2dBoth

  interface CIsoFluxCalc
     module procedure CIsoFluxCalc1d
     module procedure CIsoFluxCalc2dFlux
     module procedure CIsoFluxCalc2dBoth
  end interface CIsoFluxCalc

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp,         &
       soilbiogeochem_state_inst,                                                &
       soilbiogeochem_carbonflux_inst,  soilbiogeochem_carbonstate_inst,         &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst,                            &
       iso_soilbiogeochem_carbonflux_inst,  iso_soilbiogeochem_carbonstate_inst, &
       iso_cnveg_carbonflux_inst, iso_cnveg_carbonstate_inst,                    &
       isotope)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, set the carbon isotopic flux
    ! variables (except for gap-phase mortality and fire fluxes)
    !
    ! !ARGUMENTS:
    integer                               , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                               , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                               , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                               , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(soilbiogeochem_state_type)       , intent(in)    :: soilbiogeochem_state_inst
    type(soilbiogeochem_carbonflux_type)  , intent(in)    :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type) , intent(in)    :: soilbiogeochem_carbonstate_inst
    type(cnveg_carbonflux_type)           , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)          , intent(in)    :: cnveg_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)  , intent(inout) :: iso_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type) , intent(in)    :: iso_soilbiogeochem_carbonstate_inst
    type(cnveg_carbonflux_type)           , intent(inout) :: iso_cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)          , intent(in)    :: iso_cnveg_carbonstate_inst
    character(len=*)                      , intent(in)    :: isotope         ! 'c13' or 'c14'
    !
    ! !LOCAL VARIABLES:
    integer :: fp,l,fc,cc,j,k,p
    integer :: cdp 
    !-----------------------------------------------------------------------

    associate(                                                            &
         cascade_donor_pool    => decomp_cascade_con%cascade_donor_pool , & 
         soilbiogeochem_cs     => soilbiogeochem_carbonstate_inst       , &
         soilbiogeochem_cf     => soilbiogeochem_carbonflux_inst        , &
         cnveg_cf              => cnveg_carbonflux_inst                 , &
         cnveg_cs              => cnveg_carbonstate_inst                , &
         iso_cnveg_cf          => iso_cnveg_carbonflux_inst             , &
         iso_cnveg_cs          => iso_cnveg_carbonstate_inst            , &
         iso_soilbiogeochem_cs => iso_soilbiogeochem_carbonstate_inst   , &
         iso_soilbiogeochem_cf => iso_soilbiogeochem_carbonflux_inst      &
         )

      ! patch-level non-mortality fluxes
      
      ! Note: if the variables which are arguments to CIsoFluxCalc are ever changed to NOT be
      ! pointers, then the CIsoFluxCalc routine will need to be changed to declare the bounds
      ! of each argument, these bounds will need to be passed in, and - importantly for
      ! threading to work properly - the subroutine calls will need to be changed so that
      ! instead of 'call CIsoFluxCalc(foo, ...)' we have 'call CIsoFluxCalc(foo(begp:endp), ...)'.
     
      call CIsoFluxCalc(&
           iso_cnveg_cf%leafc_xfer_to_leafc_patch           , cnveg_cf%leafc_xfer_to_leafc_patch, &
           iso_cnveg_cs%leafc_xfer_patch                    , cnveg_cs%leafc_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%frootc_xfer_to_frootc_patch         , cnveg_cf%frootc_xfer_to_frootc_patch, &
           iso_cnveg_cs%frootc_xfer_patch                   , cnveg_cs%frootc_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%livestemc_xfer_to_livestemc_patch   , cnveg_cf%livestemc_xfer_to_livestemc_patch, &
           iso_cnveg_cs%livestemc_xfer_patch                , cnveg_cs%livestemc_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%deadstemc_xfer_to_deadstemc_patch   , cnveg_cf%deadstemc_xfer_to_deadstemc_patch, &
           iso_cnveg_cs%deadstemc_xfer_patch                , cnveg_cs%deadstemc_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%livecrootc_xfer_to_livecrootc_patch , cnveg_cf%livecrootc_xfer_to_livecrootc_patch, &
           iso_cnveg_cs%livecrootc_xfer_patch               , cnveg_cs%livecrootc_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%deadcrootc_xfer_to_deadcrootc_patch , cnveg_cf%deadcrootc_xfer_to_deadcrootc_patch, &
           iso_cnveg_cs%deadcrootc_xfer_patch               , cnveg_cs%deadcrootc_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%leafc_to_litter_patch               , cnveg_cf%leafc_to_litter_patch, &
           iso_cnveg_cs%leafc_patch                         , cnveg_cs%leafc_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%frootc_to_litter_patch              , cnveg_cf%frootc_to_litter_patch, &
           iso_cnveg_cs%frootc_patch                        , cnveg_cs%frootc_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%livestemc_to_deadstemc_patch        , cnveg_cf%livestemc_to_deadstemc_patch, &
           iso_cnveg_cs%livestemc_patch                     , cnveg_cs%livestemc_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%livecrootc_to_deadcrootc_patch      , cnveg_cf%livecrootc_to_deadcrootc_patch, &
           iso_cnveg_cs%livecrootc_patch                    , cnveg_cs%livecrootc_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%leaf_curmr_patch                    , cnveg_cf%leaf_curmr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%froot_curmr_patch                   , cnveg_cf%froot_curmr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%livestem_curmr_patch                , cnveg_cf%livestem_curmr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%livecroot_curmr_patch               , cnveg_cf%livecroot_curmr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%leaf_xsmr_patch                     , cnveg_cf%leaf_xsmr_patch, &
           iso_cnveg_cs%totvegc_patch                       , cnveg_cs%totvegc_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%froot_xsmr_patch                    , cnveg_cf%froot_xsmr_patch, &
           iso_cnveg_cs%totvegc_patch                       , cnveg_cs%totvegc_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%livestem_xsmr_patch                 , cnveg_cf%livestem_xsmr_patch, &
           iso_cnveg_cs%totvegc_patch                       , cnveg_cs%totvegc_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%livecroot_xsmr_patch                , cnveg_cf%livecroot_xsmr_patch, &
           iso_cnveg_cs%totvegc_patch                       , cnveg_cs%totvegc_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_xsmrpool_patch             , cnveg_cf%cpool_to_xsmrpool_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_leafc_patch                , cnveg_cf%cpool_to_leafc_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_leafc_storage_patch        , cnveg_cf%cpool_to_leafc_storage_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_frootc_patch               , cnveg_cf%cpool_to_frootc_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_frootc_storage_patch       , cnveg_cf%cpool_to_frootc_storage_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_livestemc_patch            , cnveg_cf%cpool_to_livestemc_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_livestemc_storage_patch    , cnveg_cf%cpool_to_livestemc_storage_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_deadstemc_patch            , cnveg_cf%cpool_to_deadstemc_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_deadstemc_storage_patch    , cnveg_cf%cpool_to_deadstemc_storage_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_livecrootc_patch           , cnveg_cf%cpool_to_livecrootc_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_livecrootc_storage_patch   , cnveg_cf%cpool_to_livecrootc_storage_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_deadcrootc_patch           , cnveg_cf%cpool_to_deadcrootc_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_deadcrootc_storage_patch   , cnveg_cf%cpool_to_deadcrootc_storage_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_leaf_gr_patch                 , cnveg_cf%cpool_leaf_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_froot_gr_patch                , cnveg_cf%cpool_froot_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_livestem_gr_patch             , cnveg_cf%cpool_livestem_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_deadstem_gr_patch             , cnveg_cf%cpool_deadstem_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_livecroot_gr_patch            , cnveg_cf%cpool_livecroot_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_deadcroot_gr_patch            , cnveg_cf%cpool_deadcroot_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_leaf_storage_gr_patch         , cnveg_cf%cpool_leaf_storage_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_froot_storage_gr_patch        , cnveg_cf%cpool_froot_storage_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_livestem_storage_gr_patch     , cnveg_cf%cpool_livestem_storage_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_deadstem_storage_gr_patch     , cnveg_cf%cpool_deadstem_storage_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_livecroot_storage_gr_patch    , cnveg_cf%cpool_livecroot_storage_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_deadcroot_storage_gr_patch    , cnveg_cf%cpool_deadcroot_storage_gr_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_gresp_storage_patch        , cnveg_cf%cpool_to_gresp_storage_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%transfer_leaf_gr_patch              , cnveg_cf%transfer_leaf_gr_patch, &
           iso_cnveg_cs%gresp_xfer_patch                    , cnveg_cs%gresp_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%transfer_froot_gr_patch             , cnveg_cf%transfer_froot_gr_patch, &
           iso_cnveg_cs%gresp_xfer_patch                    , cnveg_cs%gresp_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%transfer_livestem_gr_patch          , cnveg_cf%transfer_livestem_gr_patch, &
           iso_cnveg_cs%gresp_xfer_patch                    , cnveg_cs%gresp_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%transfer_deadstem_gr_patch          , cnveg_cf%transfer_deadstem_gr_patch, &
           iso_cnveg_cs%gresp_xfer_patch                    , cnveg_cs%gresp_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%transfer_livecroot_gr_patch         , cnveg_cf%transfer_livecroot_gr_patch, &
           iso_cnveg_cs%gresp_xfer_patch                    , cnveg_cs%gresp_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%transfer_deadcroot_gr_patch         , cnveg_cf%transfer_deadcroot_gr_patch, &
           iso_cnveg_cs%gresp_xfer_patch                    , cnveg_cs%gresp_xfer_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%leafc_storage_to_xfer_patch         , cnveg_cf%leafc_storage_to_xfer_patch, &
           iso_cnveg_cs%leafc_storage_patch                 , cnveg_cs%leafc_storage_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%frootc_storage_to_xfer_patch        , cnveg_cf%frootc_storage_to_xfer_patch, &
           iso_cnveg_cs%frootc_storage_patch                , cnveg_cs%frootc_storage_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%livestemc_storage_to_xfer_patch     , cnveg_cf%livestemc_storage_to_xfer_patch, &
           iso_cnveg_cs%livestemc_storage_patch             , cnveg_cs%livestemc_storage_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%deadstemc_storage_to_xfer_patch     , cnveg_cf%deadstemc_storage_to_xfer_patch, &
           iso_cnveg_cs%deadstemc_storage_patch             , cnveg_cs%deadstemc_storage_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%livecrootc_storage_to_xfer_patch    , cnveg_cf%livecrootc_storage_to_xfer_patch, &
           iso_cnveg_cs%livecrootc_storage_patch            , cnveg_cs%livecrootc_storage_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%deadcrootc_storage_to_xfer_patch    , cnveg_cf%deadcrootc_storage_to_xfer_patch, &
           iso_cnveg_cs%deadcrootc_storage_patch            , cnveg_cs%deadcrootc_storage_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gresp_storage_to_xfer_patch         , cnveg_cf%gresp_storage_to_xfer_patch, &
           iso_cnveg_cs%gresp_storage_patch                 , cnveg_cs%gresp_storage_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%soilc_change_patch                  , cnveg_cf%soilc_change_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      ! Note that cpool_to_resp_patch is a diagnostic flux and therefore this Iso flux calculation
      ! not strictly required. 
      call CIsoFluxCalc(&
           iso_cnveg_cf%cpool_to_resp_patch                 , cnveg_cf%cpool_to_resp_patch, &
           iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
           num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

      if ( use_crop )then
         call CIsoFluxCalc(&
              iso_cnveg_cf%reproductivec_xfer_to_reproductivec_patch, &
              cnveg_cf%reproductivec_xfer_to_reproductivec_patch, &
              iso_cnveg_cs%reproductivec_xfer_patch         , cnveg_cs%reproductivec_xfer_patch, &
              num_soilp                                     , filter_soilp, 1._r8, 0, isotope)

         call CIsoFluxCalc(&
              iso_cnveg_cf%repr_grainc_to_food_patch        , cnveg_cf%repr_grainc_to_food_patch, &
              iso_cnveg_cs%reproductivec_patch              , cnveg_cs%reproductivec_patch, &
              num_soilp                                     , filter_soilp, 1._r8, 0, isotope)

         call CIsoFluxCalc(&
              iso_cnveg_cf%leafc_to_biofuelc_patch          , cnveg_cf%leafc_to_biofuelc_patch, &
              iso_cnveg_cs%leafc_patch                      , cnveg_cs%leafc_patch, &
              num_soilp                                     , filter_soilp, 1._r8, 0, isotope)

         call CIsoFluxCalc(&
              iso_cnveg_cf%livestemc_to_biofuelc_patch      , cnveg_cf%livestemc_to_biofuelc_patch, &
              iso_cnveg_cs%livestemc_patch                  , cnveg_cs%livestemc_patch, &
              num_soilp                                     , filter_soilp, 1._r8, 0, isotope)

         call CIsoFluxCalc(&
              iso_cnveg_cf%leafc_to_removedresiduec_patch   , cnveg_cf%leafc_to_removedresiduec_patch, &
              iso_cnveg_cs%leafc_patch                      , cnveg_cs%leafc_patch, &
              num_soilp                                     , filter_soilp, 1._r8, 0, isotope)

         call CIsoFluxCalc(&
              iso_cnveg_cf%livestemc_to_removedresiduec_patch, cnveg_cf%livestemc_to_removedresiduec_patch, &
              iso_cnveg_cs%livestemc_patch                  , cnveg_cs%livestemc_patch, &
              num_soilp                                     , filter_soilp, 1._r8, 0, isotope)

         call CIsoFluxCalc(&
              iso_cnveg_cf%repr_grainc_to_seed_patch        , cnveg_cf%repr_grainc_to_seed_patch, &
              iso_cnveg_cs%reproductivec_patch              , cnveg_cs%reproductivec_patch, &
              num_soilp                                     , filter_soilp, 1._r8, 0, isotope)

         call CIsoFluxCalc(&
              iso_cnveg_cf%repr_structurec_to_cropprod_patch, cnveg_cf%repr_structurec_to_cropprod_patch, &
              iso_cnveg_cs%reproductivec_patch              , cnveg_cs%reproductivec_patch, &
              num_soilp                                     , filter_soilp, 1._r8, 0, isotope)

         call CIsoFluxCalc(&
              iso_cnveg_cf%repr_structurec_to_litter_patch  , cnveg_cf%repr_structurec_to_litter_patch, &
              iso_cnveg_cs%reproductivec_patch              , cnveg_cs%reproductivec_patch, &
              num_soilp                                     , filter_soilp, 1._r8, 0, isotope)

         call CIsoFluxCalc(&
              iso_cnveg_cf%crop_seedc_to_leaf_patch         , cnveg_cf%crop_seedc_to_leaf_patch, &
              iso_cnveg_cs%totvegc_patch                    , cnveg_cs%totvegc_patch, &
              num_soilp                                     , filter_soilp, 1._r8, 0, isotope)

         call CIsoFluxCalc(&
              iso_cnveg_cf%reproductive_curmr_patch         , cnveg_cf%reproductive_curmr_patch, &
              iso_cnveg_cs%cpool_patch                      , cnveg_cs%cpool_patch, &
              num_soilp                                     , filter_soilp, 1._r8, 0, isotope)

         call CIsoFluxCalc(&
              iso_cnveg_cf%reproductive_xsmr_patch          , cnveg_cf%reproductive_xsmr_patch, &
              iso_cnveg_cs%totvegc_patch                    , cnveg_cs%totvegc_patch, &
              num_soilp                                     , filter_soilp, 1._r8, 0, isotope)

         call CIsoFluxCalc(&
              iso_cnveg_cf%cpool_reproductive_gr_patch      , cnveg_cf%cpool_reproductive_gr_patch, &
              iso_cnveg_cs%cpool_patch                      , cnveg_cs%cpool_patch, &
              num_soilp                                     , filter_soilp, 1._r8, 0, isotope)

         call CIsoFluxCalc(&
              iso_cnveg_cf%cpool_to_reproductivec_patch     , cnveg_cf%cpool_to_reproductivec_patch, &
              iso_cnveg_cs%cpool_patch                      , cnveg_cs%cpool_patch, &
              num_soilp                                     , filter_soilp, 1._r8, 0, isotope)

         call CIsoFluxCalc(&
              iso_cnveg_cf%cpool_to_reproductivec_storage_patch, cnveg_cf%cpool_to_reproductivec_storage_patch, &
              iso_cnveg_cs%cpool_patch                         , cnveg_cs%cpool_patch, &
              num_soilp                                        , filter_soilp, 1._r8, 0, isotope)

         call CIsoFluxCalc(&
              iso_cnveg_cf%transfer_reproductive_gr_patch   , cnveg_cf%transfer_reproductive_gr_patch, &
              iso_cnveg_cs%gresp_xfer_patch                 , cnveg_cs%gresp_xfer_patch, &
              num_soilp                                     , filter_soilp, 1._r8, 0, isotope)

         call CIsoFluxCalc(&
              iso_cnveg_cf%cpool_reproductive_storage_gr_patch, cnveg_cf%cpool_reproductive_storage_gr_patch, &
              iso_cnveg_cs%cpool_patch                        , cnveg_cs%cpool_patch, &
              num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

         call CIsoFluxCalc(&
              iso_cnveg_cf%reproductivec_storage_to_xfer_patch, cnveg_cf%reproductivec_storage_to_xfer_patch, &
              iso_cnveg_cs%reproductivec_storage_patch        , cnveg_cs%reproductivec_storage_patch, &
              num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

         call CIsoFluxCalc(&
              iso_cnveg_cf%livestemc_to_litter_patch        , cnveg_cf%livestemc_to_litter_patch, &
              iso_cnveg_cs%livestemc_patch                  , cnveg_cs%livestemc_patch, &
              num_soilp                                     , filter_soilp, 1._r8, 0, isotope)

         do fp = 1,num_soilp
            p = filter_soilp(fp)
            iso_cnveg_cf%crop_harvestc_to_cropprodc_patch(p) = &
                 iso_cnveg_cf%leafc_to_biofuelc_patch(p) + &
                 iso_cnveg_cf%livestemc_to_biofuelc_patch(p) + &
                 iso_cnveg_cf%leafc_to_removedresiduec_patch(p) + &
                 iso_cnveg_cf%livestemc_to_removedresiduec_patch(p)
         end do

         if (use_grainproduct) then
            do k = repr_grain_min, repr_grain_max
               do fp = 1,num_soilp
                  p = filter_soilp(fp)
                  iso_cnveg_cf%crop_harvestc_to_cropprodc_patch(p) = &
                       iso_cnveg_cf%crop_harvestc_to_cropprodc_patch(p) + &
                       iso_cnveg_cf%repr_grainc_to_food_patch(p,k)
               end do
            end do
         end if

         do k = 1, nrepr
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               iso_cnveg_cf%reproductive_mr_patch(p,k) = &
                    iso_cnveg_cf%reproductive_xsmr_patch(p,k) + &
                    iso_cnveg_cf%reproductive_curmr_patch(p,k)
            end do
         end do

         do k = repr_structure_min, repr_structure_max
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               iso_cnveg_cf%crop_harvestc_to_cropprodc_patch(p) = &
                    iso_cnveg_cf%crop_harvestc_to_cropprodc_patch(p) + &
                    iso_cnveg_cf%repr_structurec_to_cropprod_patch(p,k)
            end do
         end do
      end if

      ! call routine to shift patch-level litterfall fluxes to column, for isotopes
      ! the non-isotope version of this routine is called in CNPhenologyMod.F90
      ! For later clean-up, it would be possible to generalize this function to operate on a single 
      ! patch-to-column flux.

      call CNCIsoLitterToColumn(num_soilp, filter_soilp, soilbiogeochem_state_inst, iso_cnveg_carbonflux_inst)

      ! column-level non-mortality fluxes

      do fc = 1,num_soilc
         cc = filter_soilc(fc)
         do j = 1, nlevdecomp
            do l = 1, ndecomp_cascade_transitions
               cdp = cascade_donor_pool(l)
               if ( soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,cdp) /= 0._r8) then
                  iso_soilbiogeochem_cf%decomp_cascade_hr_vr_col(cc,j,l)  =  &
                      soilbiogeochem_cf%decomp_cascade_hr_vr_col(cc,j,l) * &
                      (iso_soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,cdp) &
                         / soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,cdp)) * 1._r8
               else
                  iso_soilbiogeochem_cf%decomp_cascade_hr_vr_col(cc,j,l) = 0._r8
               end if
            end do
         end do
      end do

      do fc = 1,num_soilc
         cc = filter_soilc(fc)
         do j = 1, nlevdecomp
            do l = 1, ndecomp_cascade_transitions
               cdp = cascade_donor_pool(l)
               if ( soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,cdp) /= 0._r8) then
                  iso_soilbiogeochem_cf%decomp_cascade_ctransfer_vr_col(cc,j,l)  =  &
                      soilbiogeochem_cf%decomp_cascade_ctransfer_vr_col(cc,j,l) * &
                      (iso_soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,cdp) &
                      / soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,cdp)) * 1._r8
               else
                  iso_soilbiogeochem_cf%decomp_cascade_ctransfer_vr_col(cc,j,l) = 0._r8
               end if
            end do
         end do
      end do

    end associate

  end subroutine CIsoFlux1

  !-----------------------------------------------------------------------
  subroutine CIsoFlux2(num_soilp, filter_soilp, &
       soilbiogeochem_state_inst, &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst, &
       iso_cnveg_carbonflux_inst, iso_cnveg_carbonstate_inst, isotope)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, set the carbon isotopic fluxes for gap mortality
    !
    ! !ARGUMENTS:
    integer                         , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
    type(cnveg_carbonflux_type)     , intent(in)    :: cnveg_carbonflux_inst 
    type(cnveg_carbonstate_type)    , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: iso_cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)    , intent(in)    :: iso_cnveg_carbonstate_inst
    character(len=*)                , intent(in)    :: isotope         ! 'c13' or 'c14'

    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

    associate(                                               &
         cnveg_cf     => cnveg_carbonflux_inst     , &
         cnveg_cs     => cnveg_carbonstate_inst    , &
         iso_cnveg_cf => iso_cnveg_carbonflux_inst , &
         iso_cnveg_cs => iso_cnveg_carbonstate_inst        &
         )

      ! patch-level gap mortality fluxes
      call CIsoFluxCalc(&
           iso_cnveg_cf%m_leafc_to_litter_patch                          , cnveg_cf%m_leafc_to_litter_patch, &
           iso_cnveg_cs%leafc_patch                                      , cnveg_cs%leafc_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)
      call CIsoFluxCalc(&
           iso_cnveg_cf%m_leafc_storage_to_litter_patch                  , cnveg_cf%m_leafc_storage_to_litter_patch, &
           iso_cnveg_cs%leafc_storage_patch                              , cnveg_cs%leafc_storage_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_leafc_xfer_to_litter_patch                     , cnveg_cf%m_leafc_xfer_to_litter_patch, &
           iso_cnveg_cs%leafc_xfer_patch                                 , cnveg_cs%leafc_xfer_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_frootc_to_litter_patch                         , cnveg_cf%m_frootc_to_litter_patch, &
           iso_cnveg_cs%frootc_patch                                     , cnveg_cs%frootc_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_frootc_storage_to_litter_patch                 , cnveg_cf%m_frootc_storage_to_litter_patch, &
           iso_cnveg_cs%frootc_storage_patch                             , cnveg_cs%frootc_storage_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_frootc_xfer_to_litter_patch                    , cnveg_cf%m_frootc_xfer_to_litter_patch, &
           iso_cnveg_cs%frootc_xfer_patch                                , cnveg_cs%frootc_xfer_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livestemc_to_litter_patch                      , cnveg_cf%m_livestemc_to_litter_patch, &
           iso_cnveg_cs%livestemc_patch                                  , cnveg_cs%livestemc_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livestemc_storage_to_litter_patch              , cnveg_cf%m_livestemc_storage_to_litter_patch, &
           iso_cnveg_cs%livestemc_storage_patch                          , cnveg_cs%livestemc_storage_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livestemc_xfer_to_litter_patch                 , cnveg_cf%m_livestemc_xfer_to_litter_patch, &
           iso_cnveg_cs%livestemc_xfer_patch                             , cnveg_cs%livestemc_xfer_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadstemc_to_litter_patch                      , cnveg_cf%m_deadstemc_to_litter_patch, &
           iso_cnveg_cs%deadstemc_patch                                  , cnveg_cs%deadstemc_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadstemc_storage_to_litter_patch              , cnveg_cf%m_deadstemc_storage_to_litter_patch, &
           iso_cnveg_cs%deadstemc_storage_patch                          , cnveg_cs%deadstemc_storage_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadstemc_xfer_to_litter_patch                 , cnveg_cf%m_deadstemc_xfer_to_litter_patch, &
           iso_cnveg_cs%deadstemc_xfer_patch                             , cnveg_cs%deadstemc_xfer_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livecrootc_to_litter_patch                     , cnveg_cf%m_livecrootc_to_litter_patch, &
           iso_cnveg_cs%livecrootc_patch                                 , cnveg_cs%livecrootc_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livecrootc_storage_to_litter_patch             , cnveg_cf%m_livecrootc_storage_to_litter_patch, &
           iso_cnveg_cs%livecrootc_storage_patch                         , cnveg_cs%livecrootc_storage_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livecrootc_xfer_to_litter_patch                , cnveg_cf%m_livecrootc_xfer_to_litter_patch, &
           iso_cnveg_cs%livecrootc_xfer_patch                            , cnveg_cs%livecrootc_xfer_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadcrootc_to_litter_patch                     , cnveg_cf%m_deadcrootc_to_litter_patch, &
           iso_cnveg_cs%deadcrootc_patch                                 , cnveg_cs%deadcrootc_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadcrootc_storage_to_litter_patch             , cnveg_cf%m_deadcrootc_storage_to_litter_patch, &
           iso_cnveg_cs%deadcrootc_storage_patch                         , cnveg_cs%deadcrootc_storage_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadcrootc_xfer_to_litter_patch                , cnveg_cf%m_deadcrootc_xfer_to_litter_patch, &
           iso_cnveg_cs%deadcrootc_xfer_patch                            , cnveg_cs%deadcrootc_xfer_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_gresp_storage_to_litter_patch                  , cnveg_cf%m_gresp_storage_to_litter_patch, &
           iso_cnveg_cs%gresp_storage_patch                              , cnveg_cs%gresp_storage_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_gresp_xfer_to_litter_patch                     , cnveg_cf%m_gresp_xfer_to_litter_patch, &
           iso_cnveg_cs%gresp_xfer_patch                                 , cnveg_cs%gresp_xfer_patch, &
           num_soilp                                                     , filter_soilp, 1._r8, 0, isotope)

      ! call routine to shift patch-level gap mortality fluxes to column , for isotopes
      ! the non-isotope version of this routine is in CNGapMortalityMod.F90.

      call CNCIsoGapPftToColumn(num_soilp, filter_soilp, soilbiogeochem_state_inst, iso_cnveg_carbonflux_inst)

    end associate

  end subroutine CIsoFlux2

  !-----------------------------------------------------------------------
  subroutine CIsoFlux2h(num_soilp, filter_soilp,                             &
       soilbiogeochem_state_inst,                                            &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst,                        &
       iso_cnveg_carbonflux_inst, iso_cnveg_carbonstate_inst, isotope) 
    !
    ! !DESCRIPTION:
    ! set the carbon isotopic fluxes for harvest mortality
    !
    ! !ARGUMENTS:
    integer                           , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                           , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(soilbiogeochem_state_type)   , intent(in)    :: soilbiogeochem_state_inst
    type(cnveg_carbonflux_type)       , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)      , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)       , intent(inout) :: iso_cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)      , intent(in)    :: iso_cnveg_carbonstate_inst
    character(len=*)                  , intent(in)    :: isotope         ! 'c13' or 'c14'

    !-----------------------------------------------------------------------

    associate(                                               &
         cnveg_cf     => cnveg_carbonflux_inst           , &
         cnveg_cs     => cnveg_carbonstate_inst          , &
         iso_cnveg_cf => iso_cnveg_carbonflux_inst       , &
         iso_cnveg_cs => iso_cnveg_carbonstate_inst        &
         )

      ! patch-level gap mortality fluxes

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_leafc_to_litter_patch              , cnveg_cf%hrv_leafc_to_litter_patch, &
           iso_cnveg_cs%leafc_patch                            , cnveg_cs%leafc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_leafc_storage_to_litter_patch      , cnveg_cf%hrv_leafc_storage_to_litter_patch, &
           iso_cnveg_cs%leafc_storage_patch                    , cnveg_cs%leafc_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_leafc_xfer_to_litter_patch         , cnveg_cf%hrv_leafc_xfer_to_litter_patch, &
           iso_cnveg_cs%leafc_xfer_patch                       , cnveg_cs%leafc_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_frootc_to_litter_patch             , cnveg_cf%hrv_frootc_to_litter_patch, &
           iso_cnveg_cs%frootc_patch                           , cnveg_cs%frootc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_frootc_storage_to_litter_patch     , cnveg_cf%hrv_frootc_storage_to_litter_patch, &
           iso_cnveg_cs%frootc_storage_patch                   , cnveg_cs%frootc_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_frootc_xfer_to_litter_patch        , cnveg_cf%hrv_frootc_xfer_to_litter_patch, &
           iso_cnveg_cs%frootc_xfer_patch                      , cnveg_cs%frootc_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_livestemc_to_litter_patch          , cnveg_cf%hrv_livestemc_to_litter_patch, &
           iso_cnveg_cs%livestemc_patch                        , cnveg_cs%livestemc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_livestemc_storage_to_litter_patch  , cnveg_cf%hrv_livestemc_storage_to_litter_patch, &
           iso_cnveg_cs%livestemc_storage_patch                , cnveg_cs%livestemc_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_livestemc_xfer_to_litter_patch     , cnveg_cf%hrv_livestemc_xfer_to_litter_patch, &
           iso_cnveg_cs%livestemc_xfer_patch                   , cnveg_cs%livestemc_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%wood_harvestc_patch                    , cnveg_cf%wood_harvestc_patch, &
           iso_cnveg_cs%deadstemc_patch                        , cnveg_cs%deadstemc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_deadstemc_storage_to_litter_patch  , cnveg_cf%hrv_deadstemc_storage_to_litter_patch, &
           iso_cnveg_cs%deadstemc_storage_patch                , cnveg_cs%deadstemc_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_deadstemc_xfer_to_litter_patch     , cnveg_cf%hrv_deadstemc_xfer_to_litter_patch, &
           iso_cnveg_cs%deadstemc_xfer_patch                   , cnveg_cs%deadstemc_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_livecrootc_to_litter_patch         , cnveg_cf%hrv_livecrootc_to_litter_patch, &
           iso_cnveg_cs%livecrootc_patch                       , cnveg_cs%livecrootc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_livecrootc_storage_to_litter_patch , cnveg_cf%hrv_livecrootc_storage_to_litter_patch, &
           iso_cnveg_cs%livecrootc_storage_patch               , cnveg_cs%livecrootc_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_livecrootc_xfer_to_litter_patch    , cnveg_cf%hrv_livecrootc_xfer_to_litter_patch, &
           iso_cnveg_cs%livecrootc_xfer_patch                  , cnveg_cs%livecrootc_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_deadcrootc_to_litter_patch         , cnveg_cf%hrv_deadcrootc_to_litter_patch, &
           iso_cnveg_cs%deadcrootc_patch                       , cnveg_cs%deadcrootc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_deadcrootc_storage_to_litter_patch , cnveg_cf%hrv_deadcrootc_storage_to_litter_patch, &
           iso_cnveg_cs%deadcrootc_storage_patch               , cnveg_cs%deadcrootc_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_deadcrootc_xfer_to_litter_patch    , cnveg_cf%hrv_deadcrootc_xfer_to_litter_patch, &
           iso_cnveg_cs%deadcrootc_xfer_patch                  , cnveg_cs%deadcrootc_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_gresp_storage_to_litter_patch      , cnveg_cf%hrv_gresp_storage_to_litter_patch, &
           iso_cnveg_cs%gresp_storage_patch                    , cnveg_cs%gresp_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%hrv_gresp_xfer_to_litter_patch         , cnveg_cf%hrv_gresp_xfer_to_litter_patch, &
           iso_cnveg_cs%gresp_xfer_patch                       , cnveg_cs%gresp_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(& 
           iso_cnveg_cf%hrv_xsmrpool_to_atm_patch              , cnveg_cf%hrv_xsmrpool_to_atm_patch, &
           iso_cnveg_cs%totvegc_patch                          , cnveg_cs%totvegc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      ! call routine to shift patch-level gap mortality fluxes to column, 
      ! for isotopes the non-isotope version of this routine is in CNGapMortalityMod.F90.

      call CNCIsoHarvestPftToColumn(num_soilp, filter_soilp, soilbiogeochem_state_inst, iso_cnveg_carbonflux_inst)

    end associate

  end subroutine CIsoFlux2h

  !-----------------------------------------------------------------------
  subroutine CIsoFlux2g(num_soilp, filter_soilp,                             &
       soilbiogeochem_state_inst,                                            &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst,                        &
       iso_cnveg_carbonflux_inst, iso_cnveg_carbonstate_inst, isotope) 
    !
    ! !DESCRIPTION:
    ! set the carbon isotopic fluxes for gross unrepresented landcover change mortality
    !
    ! !ARGUMENTS:
    integer                           , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                           , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(soilbiogeochem_state_type)   , intent(in)    :: soilbiogeochem_state_inst
    type(cnveg_carbonflux_type)       , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)      , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)       , intent(inout) :: iso_cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)      , intent(in)    :: iso_cnveg_carbonstate_inst
    character(len=*)                  , intent(in)    :: isotope         ! 'c13' or 'c14'

    !-----------------------------------------------------------------------

    associate(                                               &
         cnveg_cf     => cnveg_carbonflux_inst           , &
         cnveg_cs     => cnveg_carbonstate_inst          , &
         iso_cnveg_cf => iso_cnveg_carbonflux_inst       , &
         iso_cnveg_cs => iso_cnveg_carbonstate_inst        &
         )

      ! patch-level gap mortality fluxes

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_leafc_to_litter_patch              , cnveg_cf%gru_leafc_to_litter_patch, &
           iso_cnveg_cs%leafc_patch                            , cnveg_cs%leafc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_leafc_storage_to_atm_patch         , cnveg_cf%gru_leafc_storage_to_atm_patch, &
           iso_cnveg_cs%leafc_storage_patch                    , cnveg_cs%leafc_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_leafc_xfer_to_atm_patch            , cnveg_cf%gru_leafc_xfer_to_atm_patch, &
           iso_cnveg_cs%leafc_xfer_patch                       , cnveg_cs%leafc_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_frootc_to_litter_patch             , cnveg_cf%gru_frootc_to_litter_patch, &
           iso_cnveg_cs%frootc_patch                           , cnveg_cs%frootc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_frootc_storage_to_atm_patch        , cnveg_cf%gru_frootc_storage_to_atm_patch, &
           iso_cnveg_cs%frootc_storage_patch                   , cnveg_cs%frootc_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_frootc_xfer_to_atm_patch           , cnveg_cf%gru_frootc_xfer_to_atm_patch, &
           iso_cnveg_cs%frootc_xfer_patch                      , cnveg_cs%frootc_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_livestemc_to_atm_patch             , cnveg_cf%gru_livestemc_to_atm_patch, &
           iso_cnveg_cs%livestemc_patch                        , cnveg_cs%livestemc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_livestemc_storage_to_atm_patch     , cnveg_cf%gru_livestemc_storage_to_atm_patch, &
           iso_cnveg_cs%livestemc_storage_patch                , cnveg_cs%livestemc_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_livestemc_xfer_to_atm_patch        , cnveg_cf%gru_livestemc_xfer_to_atm_patch, &
           iso_cnveg_cs%livestemc_xfer_patch                   , cnveg_cs%livestemc_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_deadstemc_to_atm_patch             , cnveg_cf%gru_deadstemc_to_atm_patch, &
           iso_cnveg_cs%deadstemc_patch                        , cnveg_cs%deadstemc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_wood_productc_gain_patch           , cnveg_cf%gru_wood_productc_gain_patch, &
           iso_cnveg_cs%deadstemc_patch                        , cnveg_cs%deadstemc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_deadstemc_storage_to_atm_patch     , cnveg_cf%gru_deadstemc_storage_to_atm_patch, &
           iso_cnveg_cs%deadstemc_storage_patch                , cnveg_cs%deadstemc_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_deadstemc_xfer_to_atm_patch        , cnveg_cf%gru_deadstemc_xfer_to_atm_patch, &
           iso_cnveg_cs%deadstemc_xfer_patch                   , cnveg_cs%deadstemc_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_livecrootc_to_litter_patch         , cnveg_cf%gru_livecrootc_to_litter_patch, &
           iso_cnveg_cs%livecrootc_patch                       , cnveg_cs%livecrootc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_livecrootc_storage_to_atm_patch    , cnveg_cf%gru_livecrootc_storage_to_atm_patch, &
           iso_cnveg_cs%livecrootc_storage_patch               , cnveg_cs%livecrootc_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_livecrootc_xfer_to_atm_patch       , cnveg_cf%gru_livecrootc_xfer_to_atm_patch, &
           iso_cnveg_cs%livecrootc_xfer_patch                  , cnveg_cs%livecrootc_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_deadcrootc_to_litter_patch         , cnveg_cf%gru_deadcrootc_to_litter_patch, &
           iso_cnveg_cs%deadcrootc_patch                       , cnveg_cs%deadcrootc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_deadcrootc_storage_to_atm_patch    , cnveg_cf%gru_deadcrootc_storage_to_atm_patch, &
           iso_cnveg_cs%deadcrootc_storage_patch               , cnveg_cs%deadcrootc_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_deadcrootc_xfer_to_atm_patch       , cnveg_cf%gru_deadcrootc_xfer_to_atm_patch, &
           iso_cnveg_cs%deadcrootc_xfer_patch                  , cnveg_cs%deadcrootc_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_gresp_storage_to_atm_patch         , cnveg_cf%gru_gresp_storage_to_atm_patch, &
           iso_cnveg_cs%gresp_storage_patch                    , cnveg_cs%gresp_storage_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%gru_gresp_xfer_to_atm_patch            , cnveg_cf%gru_gresp_xfer_to_atm_patch, &
           iso_cnveg_cs%gresp_xfer_patch                       , cnveg_cs%gresp_xfer_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(& 
           iso_cnveg_cf%gru_xsmrpool_to_atm_patch              , cnveg_cf%gru_xsmrpool_to_atm_patch, &
           iso_cnveg_cs%totvegc_patch                          , cnveg_cs%totvegc_patch, &
           num_soilp                                           , filter_soilp, 1._r8, 0, isotope)

      ! call routine to shift patch-level gap mortality fluxes to column, 
      ! for isotopes the non-isotope version of this routine is in CNGapMortalityMod.F90.

      call CNCIsoGrossUnrepPftToColumn(num_soilp, filter_soilp, soilbiogeochem_state_inst, iso_cnveg_carbonflux_inst)

    end associate

  end subroutine CIsoFlux2g

  !-----------------------------------------------------------------------
  subroutine CIsoFlux3(num_soilp, filter_soilp,                             &
       soilbiogeochem_state_inst , soilbiogeochem_carbonstate_inst,         &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst,                       &
       iso_cnveg_carbonflux_inst, iso_cnveg_carbonstate_inst,               &
       iso_soilbiogeochem_carbonstate_inst, isotope)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, set the carbon isotopic fluxes for fire mortality
    !
    ! !ARGUMENTS:
    integer                               , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                               , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(soilbiogeochem_state_type)       , intent(in)    :: soilbiogeochem_state_inst
    type(soilbiogeochem_carbonstate_type) , intent(in)    :: soilbiogeochem_carbonstate_inst
    type(cnveg_carbonflux_type)           , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)          , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)           , intent(inout) :: iso_cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)          , intent(in)    :: iso_cnveg_carbonstate_inst
    type(soilbiogeochem_carbonstate_type) , intent(in)    :: iso_soilbiogeochem_carbonstate_inst
    character(len=*)                      , intent(in)    :: isotope         ! 'c13' or 'c14'
    !
    ! !LOCAL VARIABLES:
    integer :: fp,pp,l,cc,j,i
    !-----------------------------------------------------------------------

    associate(                                                                 &
         ivt                                 => patch%itype                                                  , & ! Input:  [integer  (:)   ]  patch vegetation type                                
         wtcol                               => patch%wtcol                                                  , & ! Input:  [real(r8) (:)   ]  weight (relative to column) for this patch (0-1)    
         croot_prof                          => soilbiogeochem_state_inst%croot_prof_patch                   , & ! Input: [real(r8) (:,:) ]  (1/m) profile of coarse roots                          
         stem_prof                           => soilbiogeochem_state_inst%stem_prof_patch                    , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems                                 
         leaf_prof                           => soilbiogeochem_state_inst%leaf_prof_patch                    , & ! Input: [real(r8) (:,:) ]  (1/m) profile of leaves                          
         froot_prof                          => soilbiogeochem_state_inst%froot_prof_patch                   , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                                 
         soilbiogeochem_cs                   => soilbiogeochem_carbonstate_inst                              , &
         cnveg_cf                            => cnveg_carbonflux_inst                                        , &
         cnveg_cs                            => cnveg_carbonstate_inst                                       , &
         iso_cnveg_cf                        => iso_cnveg_carbonflux_inst                                    , &
         iso_cnveg_cs                        => iso_cnveg_carbonstate_inst                                   , &
         iso_soilbiogeochem_cs               => iso_soilbiogeochem_carbonstate_inst                          , &
         lf_f                                => pftcon%lf_f                                                  , & ! Input: [real(r8) (:,:)] leaf litter fractions
         fr_f                                => pftcon%fr_f                                                    & ! Input: [real(r8) (:,:)] fine root litter fractions
         )

      ! patch-level fire mortality fluxes

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_leafc_to_fire_patch              , cnveg_cf%m_leafc_to_fire_patch, &
           iso_cnveg_cs%leafc_patch                        , cnveg_cs%leafc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_leafc_storage_to_fire_patch      , cnveg_cf%m_leafc_storage_to_fire_patch, &
           iso_cnveg_cs%leafc_storage_patch                , cnveg_cs%leafc_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_leafc_xfer_to_fire_patch         , cnveg_cf%m_leafc_xfer_to_fire_patch, &
           iso_cnveg_cs%leafc_xfer_patch                   , cnveg_cs%leafc_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_frootc_to_fire_patch             , cnveg_cf%m_frootc_to_fire_patch, &
           iso_cnveg_cs%frootc_patch                       , cnveg_cs%frootc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_frootc_storage_to_fire_patch     , cnveg_cf%m_frootc_storage_to_fire_patch, &
           iso_cnveg_cs%frootc_storage_patch               , cnveg_cs%frootc_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_frootc_xfer_to_fire_patch        , cnveg_cf%m_frootc_xfer_to_fire_patch, &
           iso_cnveg_cs%frootc_xfer_patch                  , cnveg_cs%frootc_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livestemc_to_fire_patch          , cnveg_cf%m_livestemc_to_fire_patch, &
           iso_cnveg_cs%livestemc_patch                    , cnveg_cs%livestemc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livestemc_storage_to_fire_patch  , cnveg_cf%m_livestemc_storage_to_fire_patch, &
           iso_cnveg_cs%livestemc_storage_patch            , cnveg_cs%livestemc_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livestemc_xfer_to_fire_patch     , cnveg_cf%m_livestemc_xfer_to_fire_patch, &
           iso_cnveg_cs%livestemc_xfer_patch               , cnveg_cs%livestemc_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadstemc_to_fire_patch          , cnveg_cf%m_deadstemc_to_fire_patch, &
           iso_cnveg_cs%deadstemc_patch                    , cnveg_cs%deadstemc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadstemc_to_litter_fire_patch   , cnveg_cf%m_deadstemc_to_litter_fire_patch, &
           iso_cnveg_cs%deadstemc_patch                    , cnveg_cs%deadstemc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadstemc_storage_to_fire_patch  , cnveg_cf%m_deadstemc_storage_to_fire_patch, &
           iso_cnveg_cs%deadstemc_storage_patch            , cnveg_cs%deadstemc_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadstemc_xfer_to_fire_patch     , cnveg_cf%m_deadstemc_xfer_to_fire_patch, &
           iso_cnveg_cs%deadstemc_xfer_patch               , cnveg_cs%deadstemc_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livecrootc_to_fire_patch         , cnveg_cf%m_livecrootc_to_fire_patch, &
           iso_cnveg_cs%livecrootc_patch                   , cnveg_cs%livecrootc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livecrootc_storage_to_fire_patch , cnveg_cf%m_livecrootc_storage_to_fire_patch, &
           iso_cnveg_cs%livecrootc_storage_patch           , cnveg_cs%livecrootc_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livecrootc_xfer_to_fire_patch    , cnveg_cf%m_livecrootc_xfer_to_fire_patch, &
           iso_cnveg_cs%livecrootc_xfer_patch              , cnveg_cs%livecrootc_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadcrootc_to_fire_patch         , cnveg_cf%m_deadcrootc_to_fire_patch, &
           iso_cnveg_cs%deadcrootc_patch                   , cnveg_cs%deadcrootc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadcrootc_to_litter_fire_patch  , cnveg_cf%m_deadcrootc_to_litter_fire_patch, &
           iso_cnveg_cs%deadcrootc_patch                   , cnveg_cs%deadcrootc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadcrootc_storage_to_fire_patch , cnveg_cf%m_deadcrootc_storage_to_fire_patch, &
           iso_cnveg_cs%deadcrootc_storage_patch           , cnveg_cs%deadcrootc_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadcrootc_xfer_to_fire_patch    , cnveg_cf%m_deadcrootc_xfer_to_fire_patch, &
           iso_cnveg_cs%deadcrootc_xfer_patch              , cnveg_cs%deadcrootc_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_gresp_storage_to_fire_patch      , cnveg_cf%m_gresp_storage_to_fire_patch, &
           iso_cnveg_cs%gresp_storage_patch                , cnveg_cs%gresp_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_gresp_xfer_to_fire_patch         , cnveg_cf%m_gresp_xfer_to_fire_patch, &
           iso_cnveg_cs%gresp_xfer_patch                   , cnveg_cs%gresp_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_leafc_to_litter_fire_patch       , cnveg_cf%m_leafc_to_litter_fire_patch, &
           iso_cnveg_cs%leafc_patch                        , cnveg_cs%leafc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_leafc_storage_to_litter_fire_patch, cnveg_cf%m_leafc_storage_to_litter_fire_patch, &
           iso_cnveg_cs%leafc_storage_patch                , cnveg_cs%leafc_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_leafc_xfer_to_litter_fire_patch  , cnveg_cf%m_leafc_xfer_to_litter_fire_patch, &
           iso_cnveg_cs%leafc_xfer_patch                   , cnveg_cs%leafc_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livestemc_to_litter_fire_patch   , cnveg_cf%m_livestemc_to_litter_fire_patch, &
           iso_cnveg_cs%livestemc_patch                    , cnveg_cs%livestemc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livestemc_storage_to_litter_fire_patch, cnveg_cf%m_livestemc_storage_to_litter_fire_patch, &
           iso_cnveg_cs%livestemc_storage_patch            , cnveg_cs%livestemc_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livestemc_xfer_to_litter_fire_patch, cnveg_cf%m_livestemc_xfer_to_litter_fire_patch, &
           iso_cnveg_cs%livestemc_xfer_patch               , cnveg_cs%livestemc_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livestemc_to_deadstemc_fire_patch, cnveg_cf%m_livestemc_to_deadstemc_fire_patch, &
           iso_cnveg_cs%livestemc_patch                    , cnveg_cs%livestemc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadstemc_storage_to_litter_fire_patch, cnveg_cf%m_deadstemc_storage_to_litter_fire_patch, &
           iso_cnveg_cs%deadstemc_storage_patch            , cnveg_cs%deadstemc_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadstemc_xfer_to_litter_fire_patch, cnveg_cf%m_deadstemc_xfer_to_litter_fire_patch, &
           iso_cnveg_cs%deadstemc_xfer_patch               , cnveg_cs%deadstemc_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_frootc_to_litter_fire_patch      , cnveg_cf%m_frootc_to_litter_fire_patch, &
           iso_cnveg_cs%frootc_patch                       , cnveg_cs%frootc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_frootc_storage_to_litter_fire_patch, cnveg_cf%m_frootc_storage_to_litter_fire_patch, &
           iso_cnveg_cs%frootc_storage_patch               , cnveg_cs%frootc_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_frootc_xfer_to_litter_fire_patch , cnveg_cf%m_frootc_xfer_to_litter_fire_patch, &
           iso_cnveg_cs%frootc_xfer_patch                  , cnveg_cs%frootc_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livecrootc_to_litter_fire_patch  , cnveg_cf%m_livecrootc_to_litter_fire_patch, &
           iso_cnveg_cs%livecrootc_patch                   , cnveg_cs%livecrootc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livecrootc_storage_to_litter_fire_patch, cnveg_cf%m_livecrootc_storage_to_litter_fire_patch, &
           iso_cnveg_cs%livecrootc_storage_patch           , cnveg_cs%livecrootc_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livecrootc_xfer_to_litter_fire_patch, cnveg_cf%m_livecrootc_xfer_to_litter_fire_patch, &
           iso_cnveg_cs%livecrootc_xfer_patch              , cnveg_cs%livecrootc_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_livecrootc_to_deadcrootc_fire_patch, cnveg_cf%m_livecrootc_to_deadcrootc_fire_patch, &
           iso_cnveg_cs%livecrootc_patch                   , cnveg_cs%livecrootc_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadcrootc_storage_to_litter_fire_patch, cnveg_cf%m_deadcrootc_storage_to_litter_fire_patch, &
           iso_cnveg_cs%deadcrootc_storage_patch           , cnveg_cs%deadcrootc_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_deadcrootc_xfer_to_litter_fire_patch, cnveg_cf%m_deadcrootc_xfer_to_litter_fire_patch, &
           iso_cnveg_cs%deadcrootc_xfer_patch              , cnveg_cs%deadcrootc_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_gresp_storage_to_litter_fire_patch, cnveg_cf%m_gresp_storage_to_litter_fire_patch, &
           iso_cnveg_cs%gresp_storage_patch                , cnveg_cs%gresp_storage_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)

      call CIsoFluxCalc(&
           iso_cnveg_cf%m_gresp_xfer_to_litter_fire_patch  , cnveg_cf%m_gresp_xfer_to_litter_fire_patch, &
           iso_cnveg_cs%gresp_xfer_patch                   , cnveg_cs%gresp_xfer_patch, &
           num_soilp                                       , filter_soilp, 1._r8, 0, isotope)



      ! calculate the column-level flux of deadstem and deadcrootc to cwdc as the result of fire mortality.
      do fp = 1,num_soilp
         pp = filter_soilp(fp)
         cc = patch%column(pp)
         do j = 1, nlevdecomp
            iso_cnveg_cf%fire_mortality_c_to_cwdc_col(cc,j) = &
                 iso_cnveg_cf%fire_mortality_c_to_cwdc_col(cc,j) + &
                 (iso_cnveg_cf%m_deadstemc_to_litter_fire_patch(pp) + &
                 iso_cnveg_cf%m_livestemc_to_litter_fire_patch(pp)) * &
                 patch%wtcol(pp) * stem_prof(pp,j)
            iso_cnveg_cf%fire_mortality_c_to_cwdc_col(cc,j) = &
                 iso_cnveg_cf%fire_mortality_c_to_cwdc_col(cc,j) + &
                 (iso_cnveg_cf%m_deadcrootc_to_litter_fire_patch(pp) + &
                 iso_cnveg_cf%m_livecrootc_to_litter_fire_patch(pp)) * &
                 patch%wtcol(pp) * croot_prof(pp,j)
         end do

         do j = 1, nlevdecomp
            do l = 1, ndecomp_pools
               if ( soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,l) /= 0._r8) then
                  iso_cnveg_cf%m_decomp_cpools_to_fire_vr_col(cc,j,l)  =  &
                      cnveg_cf%m_decomp_cpools_to_fire_vr_col(cc,j,l) * &
                      (iso_soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,l) / &
                           soilbiogeochem_cs%decomp_cpools_vr_col(cc,j,l)) * 1._r8
               else
                  iso_cnveg_cf%m_decomp_cpools_to_fire_vr_col(cc,j,l) = 0._r8
               end if
            end do
         end do
      end do

      do fp = 1,num_soilp
         pp = filter_soilp(fp)
         cc = patch%column(pp)
         do j = 1, nlevdecomp
            iso_cnveg_cf%m_c_to_litr_fire_col(cc,j,i_met_lit) = &
                 iso_cnveg_cf%m_c_to_litr_fire_col(cc,j,i_met_lit) + &
                 ((iso_cnveg_cf%m_leafc_to_litter_fire_patch(pp) * lf_f(ivt(pp),i_met_lit) &
                 +iso_cnveg_cf%m_leafc_storage_to_litter_fire_patch(pp) + &
                 iso_cnveg_cf%m_leafc_xfer_to_litter_fire_patch(pp) + &
                 iso_cnveg_cf%m_gresp_storage_to_litter_fire_patch(pp) &
                 +iso_cnveg_cf%m_gresp_xfer_to_litter_fire_patch(pp))*leaf_prof(pp,j) + &
                 (iso_cnveg_cf%m_frootc_to_litter_fire_patch(pp) * fr_f(ivt(pp),i_met_lit) &
                 +iso_cnveg_cf%m_frootc_storage_to_litter_fire_patch(pp) + &
                 iso_cnveg_cf%m_frootc_xfer_to_litter_fire_patch(pp))*froot_prof(pp,j) &
                 +(iso_cnveg_cf%m_livestemc_storage_to_litter_fire_patch(pp) + &
                 iso_cnveg_cf%m_livestemc_xfer_to_litter_fire_patch(pp) &
                 +iso_cnveg_cf%m_deadstemc_storage_to_litter_fire_patch(pp) + &
                 iso_cnveg_cf%m_deadstemc_xfer_to_litter_fire_patch(pp))* stem_prof(pp,j)&
                 +(iso_cnveg_cf%m_livecrootc_storage_to_litter_fire_patch(pp) + &
                 iso_cnveg_cf%m_livecrootc_xfer_to_litter_fire_patch(pp) &
                 +iso_cnveg_cf%m_deadcrootc_storage_to_litter_fire_patch(pp) + &
                 iso_cnveg_cf%m_deadcrootc_xfer_to_litter_fire_patch(pp))* croot_prof(pp,j)) * patch%wtcol(pp)    

            ! Here metabolic litter is treated differently than other
            ! types of litter, so it remains outside this litter loop,
            ! in the line above
            do i = i_met_lit+1, i_litr_max
               iso_cnveg_cf%m_c_to_litr_fire_col(cc,j,i) = &
                    iso_cnveg_cf%m_c_to_litr_fire_col(cc,j,i) + &
                    (iso_cnveg_cf%m_leafc_to_litter_fire_patch(pp) * lf_f(ivt(pp),i) * leaf_prof(pp,j) + &
                    iso_cnveg_cf%m_frootc_to_litter_fire_patch(pp) * fr_f(ivt(pp),i) * froot_prof(pp,j)) * patch%wtcol(pp)
            end do
         end do
      end do

    end associate

  end subroutine CIsoFlux3

  !-----------------------------------------------------------------------
  subroutine CNCIsoLitterToColumn (num_soilp, filter_soilp, &
       soilbiogeochem_state_inst, iso_cnveg_carbonflux_inst)
    !
    ! !DESCRIPTION:
    ! called at the end of cn_phenology to gather all patch-level litterfall fluxes
    ! to the column level and assign them to the three litter pools
    !
    ! !USES:
!DML
    use pftconMod  , only : npcropmin
    use clm_varctl , only : use_grainproduct
!DML

    ! !ARGUMENTS:
    integer                         , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: iso_cnveg_carbonflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: fp,c,p,k,j,i
    !-----------------------------------------------------------------------

    associate(                                                                                     & 
         ivt                       => patch%itype                                             , & ! Input:  [integer  (:)   ]  patch vegetation type                                
         wtcol                     => patch%wtcol                                             , & ! Input:  [real(r8) (:)   ]  weight (relative to column) for this patch (0-1)    
         
         lf_f                      => pftcon%lf_f                                             , & ! Input:  [real(r8) (:,:) ]  leaf litter fractions
         fr_f                      => pftcon%fr_f                                             , & ! Input:  [real(r8) (:,:) ]  fine root litter fractions

         leaf_prof                 => soilbiogeochem_state_inst%leaf_prof_patch               , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
         froot_prof                => soilbiogeochem_state_inst%froot_prof_patch              , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots 
         
         leafc_to_litter           => iso_cnveg_carbonflux_inst%leafc_to_litter_patch         , & ! Input:  [real(r8) (:)   ]            
         frootc_to_litter          => iso_cnveg_carbonflux_inst%frootc_to_litter_patch        , & ! Input:  [real(r8) (:)   ] 
         livestemc_to_litter       => iso_cnveg_carbonflux_inst%livestemc_to_litter_patch     , & ! Input:  [real(r8) (:)   ]
         repr_grainc_to_food       => iso_cnveg_carbonflux_inst%repr_grainc_to_food_patch     , & ! Input:  [real(r8) (:,:) ]
         repr_structurec_to_litter => iso_cnveg_carbonflux_inst%repr_structurec_to_litter_patch,& ! Input:  [real(r8) (:,:) ]
         phenology_c_to_litr_c     => iso_cnveg_carbonflux_inst%phenology_c_to_litr_c_col       & ! InOut:  [real(r8) (:,:,:) ]  C fluxes associated with phenology (litterfall and crop) to litter pools (gC/m3/s)
         )

      do j = 1, nlevdecomp
         do fp = 1,num_soilp
            p = filter_soilp(fp)
            c = patch%column(p)

            do i = i_litr_min, i_litr_max
               phenology_c_to_litr_c(c,j,i) = &
                    phenology_c_to_litr_c(c,j,i) + &
                    ! leaf litter carbon fluxes
                    leafc_to_litter(p) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j) + &
                    ! fine root litter carbon fluxes
                    frootc_to_litter(p) * fr_f(ivt(p),i) * wtcol(p) * froot_prof(p,j)
            end do

!DML
            if (ivt(p) >= npcropmin) then ! add livestemc to litter
               ! stem litter carbon fluxes
               do i = i_litr_min, i_litr_max
                  phenology_c_to_litr_c(c,j,i) = &
                       phenology_c_to_litr_c(c,j,i) + &
                       livestemc_to_litter(p) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j)
               end do

               if (.not. use_grainproduct) then
                  ! grain litter carbon fluxes
                  do i = i_litr_min, i_litr_max
                     do k = repr_grain_min, repr_grain_max
                        phenology_c_to_litr_c(c,j,i) = &
                             phenology_c_to_litr_c(c,j,i) + &
                             repr_grainc_to_food(p,k) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j)
                     end do
                  end do
               end if

               ! reproductive structure litter carbon fluxes
               do i = i_litr_min, i_litr_max
                  do k = repr_structure_min, repr_structure_max
                     phenology_c_to_litr_c(c,j,i) = &
                          phenology_c_to_litr_c(c,j,i) + &
                          repr_structurec_to_litter(p,k) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j)
                  end do
               end do
            end if
            !DML
         end do
      end do

    end associate

   end subroutine CNCIsoLitterToColumn

   !-----------------------------------------------------------------------
   subroutine CNCIsoGapPftToColumn (num_soilp, filter_soilp, &
        soilbiogeochem_state_inst, iso_cnveg_carbonflux_inst)
     !
     ! !DESCRIPTION:
     ! gather all patch-level gap mortality fluxes
     ! to the column level and assign them to the three litter pools (+ cwd pool)
     !
     ! !ARGUMENTS:
     integer                         , intent(in)    :: num_soilp         ! number of soil patches in filter
     integer                         , intent(in)    :: filter_soilp(:)   ! soil patch filter
     type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
     type(cnveg_carbonflux_type)     , intent(inout) :: iso_cnveg_carbonflux_inst
     !
     ! !LOCAL VARIABLES:
     integer :: fp,c,p,j,i  ! indices
     !-----------------------------------------------------------------------

     associate(                                                                                             & 
          ivt                            => patch%itype                                                  , & ! Input:  [integer  (:)   ]  patch vegetation type                                
          wtcol                          => patch%wtcol                                                  , & ! Input:  [real(r8) (:)   ]  patch weight relative to column (0-1)               
          
          lf_f                           => pftcon%lf_f                                                , & ! Input:  leaf litter fractions
          fr_f                           => pftcon%fr_f                                                , & ! Input:  fine root litter fractions

          leaf_prof                      => soilbiogeochem_state_inst%leaf_prof_patch                  , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
          froot_prof                     => soilbiogeochem_state_inst%froot_prof_patch                 , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
          croot_prof                     => soilbiogeochem_state_inst%croot_prof_patch                 , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of coarse roots                   
          stem_prof                      => soilbiogeochem_state_inst%stem_prof_patch                  , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems                          
          
          m_leafc_to_litter              => iso_cnveg_carbonflux_inst%m_leafc_to_litter_patch              , & ! Input:  [real(r8) (:)   ]                                                    
          m_frootc_to_litter             => iso_cnveg_carbonflux_inst%m_frootc_to_litter_patch             , & ! Input:  [real(r8) (:)   ]                                                    
          m_livestemc_to_litter          => iso_cnveg_carbonflux_inst%m_livestemc_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
          m_deadstemc_to_litter          => iso_cnveg_carbonflux_inst%m_deadstemc_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
          m_livecrootc_to_litter         => iso_cnveg_carbonflux_inst%m_livecrootc_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          m_deadcrootc_to_litter         => iso_cnveg_carbonflux_inst%m_deadcrootc_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          m_leafc_storage_to_litter      => iso_cnveg_carbonflux_inst%m_leafc_storage_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
          m_frootc_storage_to_litter     => iso_cnveg_carbonflux_inst%m_frootc_storage_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
          m_livestemc_storage_to_litter  => iso_cnveg_carbonflux_inst%m_livestemc_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
          m_deadstemc_storage_to_litter  => iso_cnveg_carbonflux_inst%m_deadstemc_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
          m_livecrootc_storage_to_litter => iso_cnveg_carbonflux_inst%m_livecrootc_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
          m_deadcrootc_storage_to_litter => iso_cnveg_carbonflux_inst%m_deadcrootc_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
          m_gresp_storage_to_litter      => iso_cnveg_carbonflux_inst%m_gresp_storage_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
          m_leafc_xfer_to_litter         => iso_cnveg_carbonflux_inst%m_leafc_xfer_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          m_frootc_xfer_to_litter        => iso_cnveg_carbonflux_inst%m_frootc_xfer_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
          m_livestemc_xfer_to_litter     => iso_cnveg_carbonflux_inst%m_livestemc_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
          m_deadstemc_xfer_to_litter     => iso_cnveg_carbonflux_inst%m_deadstemc_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
          m_livecrootc_xfer_to_litter    => iso_cnveg_carbonflux_inst%m_livecrootc_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
          m_deadcrootc_xfer_to_litter    => iso_cnveg_carbonflux_inst%m_deadcrootc_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
          m_gresp_xfer_to_litter         => iso_cnveg_carbonflux_inst%m_gresp_xfer_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          
          gap_mortality_c_to_litr_c      => iso_cnveg_carbonflux_inst%gap_mortality_c_to_litr_c_col        , & ! InOut:  [real(r8) (:,:,:) ]  C fluxes associated with gap mortality to litter pools (gC/m3/s)
          gap_mortality_c_to_cwdc        => iso_cnveg_carbonflux_inst%gap_mortality_c_to_cwdc_col            & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with gap mortality to CWD pool (gC/m3/s)
          )
          
       do j = 1, nlevdecomp
          do fp = 1,num_soilp
             p = filter_soilp(fp)
             c = patch%column(p)

             do i = i_litr_min, i_litr_max
                ! leaf gap mortality carbon fluxes
                gap_mortality_c_to_litr_c(c,j,i) = &
                     gap_mortality_c_to_litr_c(c,j,i) + &
                     m_leafc_to_litter(p) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j)
                ! fine root gap mortality carbon fluxes
                gap_mortality_c_to_litr_c(c,j,i) = &
                     gap_mortality_c_to_litr_c(c,j,i) + &
                     m_frootc_to_litter(p) * fr_f(ivt(p),i) * wtcol(p) * froot_prof(p,j)
             end do

             ! wood gap mortality carbon fluxes
             gap_mortality_c_to_cwdc(c,j)  = gap_mortality_c_to_cwdc(c,j)  + &
                  m_livestemc_to_litter(p)  * wtcol(p) * stem_prof(p,j)
             gap_mortality_c_to_cwdc(c,j)  = gap_mortality_c_to_cwdc(c,j)  + &
                  m_deadstemc_to_litter(p)  * wtcol(p) * stem_prof(p,j)
             gap_mortality_c_to_cwdc(c,j) = gap_mortality_c_to_cwdc(c,j) + &
                  m_livecrootc_to_litter(p) * wtcol(p) * croot_prof(p,j)
             gap_mortality_c_to_cwdc(c,j) = gap_mortality_c_to_cwdc(c,j) + &
                  m_deadcrootc_to_litter(p) * wtcol(p) * croot_prof(p,j)

             ! Metabolic litter is treated differently than other types
             ! of litter, so it gets this additional line after the
             ! most recent loop over all litter types
             gap_mortality_c_to_litr_c(c,j,i_met_lit) = &
                  gap_mortality_c_to_litr_c(c,j,i_met_lit) + &
                  ! storage gap mortality carbon fluxes
                  m_leafc_storage_to_litter(p) * wtcol(p) * leaf_prof(p,j) + &
                  m_frootc_storage_to_litter(p) * wtcol(p) * froot_prof(p,j) + &
                  m_livestemc_storage_to_litter(p) * wtcol(p) * stem_prof(p,j) + &
                  m_deadstemc_storage_to_litter(p) * wtcol(p) * stem_prof(p,j) + &
                  m_livecrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
                  m_deadcrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
                  m_gresp_storage_to_litter(p) * wtcol(p) * leaf_prof(p,j) + &
                  ! transfer gap mortality carbon fluxes
                  m_leafc_xfer_to_litter(p) * wtcol(p) * leaf_prof(p,j) + &
                  m_frootc_xfer_to_litter(p) * wtcol(p) * froot_prof(p,j) + &
                  m_livestemc_xfer_to_litter(p) * wtcol(p) * stem_prof(p,j) + &
                  m_deadstemc_xfer_to_litter(p) * wtcol(p) * stem_prof(p,j) + &
                  m_livecrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
                  m_deadcrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
                  m_gresp_xfer_to_litter(p) * wtcol(p) * leaf_prof(p,j)

          end do
       end do


     end associate

   end subroutine CNCIsoGapPftToColumn

   !-----------------------------------------------------------------------
   subroutine CNCIsoHarvestPftToColumn (num_soilp, filter_soilp, &
        soilbiogeochem_state_inst, iso_cnveg_carbonflux_inst)
     !
     ! !DESCRIPTION:
     ! gather all patch-level harvest mortality fluxes
     ! to the column level and assign them to the litter, cwd, and wood product pools
     !
     ! !ARGUMENTS:
     integer                         , intent(in)    :: num_soilp         ! number of soil patches in filter
     integer                         , intent(in)    :: filter_soilp(:)   ! soil patch filter
     type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
     type(cnveg_carbonflux_type)     , intent(inout) :: iso_cnveg_carbonflux_inst
     !
     ! !LOCAL VARIABLES:
     integer :: fp,c,p,j,i  ! indices
     !-----------------------------------------------------------------------

     associate(                                                                                                  & 
          ivt                              => patch%itype                                                      , & ! Input:  [integer  (:)   ]  patch vegetation type                                
          wtcol                            => patch%wtcol                                                      , & ! Input:  [real(r8) (:)   ]  patch weight relative to column (0-1)               
          
          lf_f                             => pftcon%lf_f                                                      , & ! Input:  leaf litter fractions
          fr_f                             => pftcon%fr_f                                                      , & ! Input:  fine root litter fractions

          leaf_prof                        => soilbiogeochem_state_inst%leaf_prof_patch                        , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
          froot_prof                       => soilbiogeochem_state_inst%froot_prof_patch                       , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
          croot_prof                       => soilbiogeochem_state_inst%croot_prof_patch                       , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of coarse roots                   
          stem_prof                        => soilbiogeochem_state_inst%stem_prof_patch                        , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems                          
          
          hrv_leafc_to_litter              => iso_cnveg_carbonflux_inst%hrv_leafc_to_litter_patch              , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_frootc_to_litter             => iso_cnveg_carbonflux_inst%hrv_frootc_to_litter_patch             , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_livestemc_to_litter          => iso_cnveg_carbonflux_inst%hrv_livestemc_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
          pwood_harvestc                   => iso_cnveg_carbonflux_inst%wood_harvestc_patch                    , & ! Input:  [real(r8) (:)   ]
          hrv_livecrootc_to_litter         => iso_cnveg_carbonflux_inst%hrv_livecrootc_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_deadcrootc_to_litter         => iso_cnveg_carbonflux_inst%hrv_deadcrootc_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_leafc_storage_to_litter      => iso_cnveg_carbonflux_inst%hrv_leafc_storage_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_frootc_storage_to_litter     => iso_cnveg_carbonflux_inst%hrv_frootc_storage_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_livestemc_storage_to_litter  => iso_cnveg_carbonflux_inst%hrv_livestemc_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_deadstemc_storage_to_litter  => iso_cnveg_carbonflux_inst%hrv_deadstemc_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_livecrootc_storage_to_litter => iso_cnveg_carbonflux_inst%hrv_livecrootc_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_deadcrootc_storage_to_litter => iso_cnveg_carbonflux_inst%hrv_deadcrootc_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_gresp_storage_to_litter      => iso_cnveg_carbonflux_inst%hrv_gresp_storage_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_leafc_xfer_to_litter         => iso_cnveg_carbonflux_inst%hrv_leafc_xfer_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_frootc_xfer_to_litter        => iso_cnveg_carbonflux_inst%hrv_frootc_xfer_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_livestemc_xfer_to_litter     => iso_cnveg_carbonflux_inst%hrv_livestemc_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_deadstemc_xfer_to_litter     => iso_cnveg_carbonflux_inst%hrv_deadstemc_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_livecrootc_xfer_to_litter    => iso_cnveg_carbonflux_inst%hrv_livecrootc_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_deadcrootc_xfer_to_litter    => iso_cnveg_carbonflux_inst%hrv_deadcrootc_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
          hrv_gresp_xfer_to_litter         => iso_cnveg_carbonflux_inst%hrv_gresp_xfer_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          cwood_harvestc                   => iso_cnveg_carbonflux_inst%wood_harvestc_col                      , & ! Output:  [real(r8) (:)   ]
          harvest_c_to_litr_c              => iso_cnveg_carbonflux_inst%harvest_c_to_litr_c_col                , & ! Output: [real(r8) (:,:) ]  C fluxes associated with harvest to litter pools (gC/m3/s)
          harvest_c_to_cwdc                => iso_cnveg_carbonflux_inst%harvest_c_to_cwdc_col                    & ! Output: [real(r8) (:,:) ]  C fluxes associated with harvest to CWD pool (gC/m3/s)
          )

       do j = 1, nlevdecomp
          do fp = 1,num_soilp
             p = filter_soilp(fp)
             c = patch%column(p)
                
             do i = i_litr_min, i_litr_max
                ! leaf harvest mortality carbon fluxes
                harvest_c_to_litr_c(c,j,i) = &
                   harvest_c_to_litr_c(c,j,i) + &
                   hrv_leafc_to_litter(p) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j)

                ! fine root harvest mortality carbon fluxes
                harvest_c_to_litr_c(c,j,i) = &
                   harvest_c_to_litr_c(c,j,i) + &
                   hrv_frootc_to_litter(p) * fr_f(ivt(p),i) * wtcol(p) * froot_prof(p,j)
             end do

             ! wood harvest mortality carbon fluxes
             harvest_c_to_cwdc(c,j)  = harvest_c_to_cwdc(c,j)  + &
                  hrv_livestemc_to_litter(p)  * wtcol(p) * stem_prof(p,j)
             harvest_c_to_cwdc(c,j) = harvest_c_to_cwdc(c,j) + &
                  hrv_livecrootc_to_litter(p) * wtcol(p) * croot_prof(p,j)
             harvest_c_to_cwdc(c,j) = harvest_c_to_cwdc(c,j) + &
                  hrv_deadcrootc_to_litter(p) * wtcol(p) * croot_prof(p,j)

             ! Metabolic litter is treated differently than other types
             ! of litter, so it gets this additional line after the
             ! most recent loop over all litter types
             harvest_c_to_litr_c(c,j,i_met_lit) = &
                harvest_c_to_litr_c(c,j,i_met_lit) + &
                ! storage harvest mortality carbon fluxes
                hrv_leafc_storage_to_litter(p) * wtcol(p) * leaf_prof(p,j) + &
                hrv_frootc_storage_to_litter(p) * wtcol(p) * froot_prof(p,j) + &
                hrv_livestemc_storage_to_litter(p) * wtcol(p) * stem_prof(p,j) + &
                hrv_deadstemc_storage_to_litter(p) * wtcol(p) * stem_prof(p,j) + &
                hrv_livecrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
                hrv_deadcrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
                hrv_gresp_storage_to_litter(p) * wtcol(p) * leaf_prof(p,j) + &
                ! transfer harvest mortality carbon fluxes
                hrv_leafc_xfer_to_litter(p) * wtcol(p) * leaf_prof(p,j) + &
                hrv_frootc_xfer_to_litter(p) * wtcol(p) * froot_prof(p,j) + &
                hrv_livestemc_xfer_to_litter(p) * wtcol(p) * stem_prof(p,j) + &
                hrv_deadstemc_xfer_to_litter(p) * wtcol(p) * stem_prof(p,j) + &
                hrv_livecrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
                hrv_deadcrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j) + &
                hrv_gresp_xfer_to_litter(p) * wtcol(p) * leaf_prof(p,j)

          end do
       end do

       do fp = 1,num_soilp
          p = filter_soilp(fp)
          c = patch%column(p)
          cwood_harvestc(c) = cwood_harvestc(c) + &
               pwood_harvestc(p) * wtcol(p)
       end do

     end associate 

   end subroutine CNCIsoHarvestPftToColumn

   !-----------------------------------------------------------------------
   subroutine CNCIsoGrossUnrepPftToColumn (num_soilp, filter_soilp, &
        soilbiogeochem_state_inst, iso_cnveg_carbonflux_inst)
     !
     ! !DESCRIPTION:
     ! gather all patch-level gross unrepresented landcover change mortality fluxes
     ! to the column level and assign them to the litter, cwd, and wood product pools
     !
     ! !ARGUMENTS:
     integer                         , intent(in)    :: num_soilp         ! number of soil patches in filter
     integer                         , intent(in)    :: filter_soilp(:)   ! soil patch filter
     type(soilbiogeochem_state_type) , intent(in)    :: soilbiogeochem_state_inst
     type(cnveg_carbonflux_type)     , intent(inout) :: iso_cnveg_carbonflux_inst
     !
     ! !LOCAL VARIABLES:
     integer :: fp,c,p,j,i  ! indices
     !-----------------------------------------------------------------------

     associate(                                                                                                  & 
          ivt                              => patch%itype                                                      , & ! Input:  [integer  (:)   ]  patch vegetation type                                
          wtcol                            => patch%wtcol                                                      , & ! Input:  [real(r8) (:)   ]  patch weight relative to column (0-1)               
          
          lf_f                             => pftcon%lf_f                                                      , & ! Input:  leaf litter fractions
          fr_f                             => pftcon%fr_f                                                      , & ! Input:  fine root litter fractions
          
          leaf_prof                        => soilbiogeochem_state_inst%leaf_prof_patch                        , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
          froot_prof                       => soilbiogeochem_state_inst%froot_prof_patch                       , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
          croot_prof                       => soilbiogeochem_state_inst%croot_prof_patch                       , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of coarse roots                   
          stem_prof                        => soilbiogeochem_state_inst%stem_prof_patch                        , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems                          
          
          gru_leafc_to_litter              => iso_cnveg_carbonflux_inst%gru_leafc_to_litter_patch              , & ! Input:  [real(r8) (:)   ]                                                    
          gru_frootc_to_litter             => iso_cnveg_carbonflux_inst%gru_frootc_to_litter_patch             , & ! Input:  [real(r8) (:)   ]                                                    
          gru_livestemc_to_atm             => iso_cnveg_carbonflux_inst%gru_livestemc_to_atm_patch             , & ! Input:  [real(r8) (:)   ]                                                    
          gru_deadstemc_to_atm             => iso_cnveg_carbonflux_inst%gru_deadstemc_to_atm_patch             , & ! Input:  [real(r8) (:)   ]                                                    
          gru_wood_productc_gain           => iso_cnveg_carbonflux_inst%gru_wood_productc_gain_patch           , & ! Input:  [real(r8) (:)   ]
          gru_livecrootc_to_litter         => iso_cnveg_carbonflux_inst%gru_livecrootc_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          gru_deadcrootc_to_litter         => iso_cnveg_carbonflux_inst%gru_deadcrootc_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          gru_leafc_storage_to_atm         => iso_cnveg_carbonflux_inst%gru_leafc_storage_to_atm_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          gru_frootc_storage_to_atm        => iso_cnveg_carbonflux_inst%gru_frootc_storage_to_atm_patch        , & ! Input:  [real(r8) (:)   ]                                                    
          gru_livestemc_storage_to_atm     => iso_cnveg_carbonflux_inst%gru_livestemc_storage_to_atm_patch     , & ! Input:  [real(r8) (:)   ]                                                    
          gru_deadstemc_storage_to_atm     => iso_cnveg_carbonflux_inst%gru_deadstemc_storage_to_atm_patch     , & ! Input:  [real(r8) (:)   ]                                                    
          gru_livecrootc_storage_to_atm    => iso_cnveg_carbonflux_inst%gru_livecrootc_storage_to_atm_patch    , & ! Input:  [real(r8) (:)   ]                                                    
          gru_deadcrootc_storage_to_atm    => iso_cnveg_carbonflux_inst%gru_deadcrootc_storage_to_atm_patch    , & ! Input:  [real(r8) (:)   ]                                                    
          gru_gresp_storage_to_atm         => iso_cnveg_carbonflux_inst%gru_gresp_storage_to_atm_patch         , & ! Input:  [real(r8) (:)   ]                                                    
          gru_leafc_xfer_to_atm            => iso_cnveg_carbonflux_inst%gru_leafc_xfer_to_atm_patch            , & ! Input:  [real(r8) (:)   ]                                                    
          gru_frootc_xfer_to_atm           => iso_cnveg_carbonflux_inst%gru_frootc_xfer_to_atm_patch           , & ! Input:  [real(r8) (:)   ]                                                    
          gru_livestemc_xfer_to_atm        => iso_cnveg_carbonflux_inst%gru_livestemc_xfer_to_atm_patch        , & ! Input:  [real(r8) (:)   ]                                                    
          gru_deadstemc_xfer_to_atm        => iso_cnveg_carbonflux_inst%gru_deadstemc_xfer_to_atm_patch        , & ! Input:  [real(r8) (:)   ]                                                    
          gru_livecrootc_xfer_to_atm       => iso_cnveg_carbonflux_inst%gru_livecrootc_xfer_to_atm_patch       , & ! Input:  [real(r8) (:)   ]                                                    
          gru_deadcrootc_xfer_to_atm       => iso_cnveg_carbonflux_inst%gru_deadcrootc_xfer_to_atm_patch       , & ! Input:  [real(r8) (:)   ]                                                    
          gru_gresp_xfer_to_atm            => iso_cnveg_carbonflux_inst%gru_gresp_xfer_to_atm_patch            , & ! Input:  [real(r8) (:)   ]                                                    
          cwood_harvestc                   => iso_cnveg_carbonflux_inst%wood_harvestc_col                      , & ! Output:  [real(r8) (:)   ]
          gru_c_to_litr_c                  => iso_cnveg_carbonflux_inst%gru_c_to_litr_c_col                    , & ! Output: [real(r8) (:,:,:) ]  C fluxes associated with gross unrepresented landcover change to litter pools (gC/m3/s)
          gru_c_to_cwdc_c                  => iso_cnveg_carbonflux_inst%gru_c_to_cwdc_col                      , & ! Output: [real(r8) (:,:) ]  C fluxes associated with harvest to CWD pool (gC/m3/s)
          gru_wood_productc_gain_c         => iso_cnveg_carbonflux_inst%gru_wood_productc_gain_col               & ! Input:  [real(r8) (:)   ]
          )

       do j = 1, nlevdecomp
          do fp = 1,num_soilp
             p = filter_soilp(fp)
             c = patch%column(p)

             do i = i_litr_min, i_litr_max
                gru_c_to_litr_c(c,j,i) = &
                   gru_c_to_litr_c(c,j,i) + &
                   ! leaf gross unrepresented landcover change mortality carbon fluxes
                   gru_leafc_to_litter(p) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j) + &
                   ! fine root gross unrepresented landcover change mortality carbon fluxes
                   gru_frootc_to_litter(p) * fr_f(ivt(p),i) * wtcol(p) * froot_prof(p,j)
             end do

             ! coarse root gross unrepresented landcover change mortality carbon fluxes
             gru_c_to_cwdc_c(c,j) = gru_c_to_cwdc_c(c,j) + &
                  gru_livecrootc_to_litter(p) * wtcol(p) * croot_prof(p,j)
             gru_c_to_cwdc_c(c,j) = gru_c_to_cwdc_c(c,j) + &
                  gru_deadcrootc_to_litter(p) * wtcol(p) * croot_prof(p,j) 

          end do
       end do
   
       do fp = 1,num_soilp
          p = filter_soilp(fp)
          c = patch%column(p)

          ! wood gross unrepresented landcover change mortality carbon fluxes to product pools
          gru_wood_productc_gain_c(c)  = gru_wood_productc_gain_c(c)  + &
               gru_wood_productc_gain(p)  * wtcol(p)

       end do

     end associate 

   end subroutine CNCIsoGrossUnrepPftToColumn

   !-----------------------------------------------------------------------
   subroutine CIsoFluxCalc1d(&
        ciso_flux, ctot_flux, &
        ciso_state, ctot_state, &
        num, filter, frax_c13, diag, isotope)
     !
     ! !DESCRIPTION:
     ! On the radiation time step, set the carbon isotopic flux
     ! variables (except for gap-phase mortality and fire fluxes)
     !
     ! !ARGUMENTS:
     real(r8)         , intent(inout), pointer :: ciso_flux(:)  ! isoC flux
     real(r8)         , intent(in)   , pointer :: ctot_flux(:)  ! totC flux
     real(r8)         , intent(in)   , pointer :: ciso_state(:) ! isoC state, upstream pool
     real(r8)         , intent(in)   , pointer :: ctot_state(:) ! totC state, upstream pool
     real(r8)         , intent(in)             :: frax_c13      ! fractionation factor (1 = no fractionation) for C13
     integer          , intent(in)             :: num           ! number of filter members
     integer          , intent(in)             :: filter(:)     ! filter indices
     integer          , intent(in)             :: diag          ! 0=no diagnostics, 1=print diagnostics
     character(len=*) , intent(in)             :: isotope       ! 'c13' or 'c14'
     !
     ! ! LOCAL VARIABLES:
     integer  :: i,f     ! indices
     real(r8) :: temp
     real(r8) :: frax
     !-----------------------------------------------------------------------

     ! if C14, double the fractionation
     select case (isotope)
     case ('c14')
        frax = 1._r8 + (1._r8 - frax_c13) * 2._r8
     case ('c13')
        frax = frax_c13
     case default
        call endrun(msg='CNCIsoFluxMod: iso must be either c13 or c14'//errMsg(sourcefile, __LINE__))
     end select

     ! loop over the supplied filter
     do f = 1,num
        i = filter(f)
        if (ctot_state(i) /= 0._r8 .and. ciso_state(i) /= 0._r8) then
           ciso_flux(i) = ctot_flux(i) * (ciso_state(i)/ctot_state(i)) * frax
        else
           ciso_flux(i) = 0._r8
        end if

        if (diag == 1) then
           ! put diagnostic print statements here for isoC flux calculations
        end if
     end do

   end subroutine CIsoFluxCalc1d

   !-----------------------------------------------------------------------
   subroutine CIsoFluxCalc2dFlux(&
        ciso_flux, ctot_flux, &
        ciso_state, ctot_state, &
        num, filter, frax_c13, diag, isotope)
     !
     ! !DESCRIPTION:
     ! Wrapper to CIsoFluxCalc1d for just the flux being a 2-d variable
     !
     ! Loops over the second dimension of each flux variable to do a C Iso flux calc on each level
     !
     ! !ARGUMENTS:
     real(r8)         , intent(inout), pointer :: ciso_flux(:,:)  ! isoC flux
     real(r8)         , intent(in)   , pointer :: ctot_flux(:,:)  ! totC flux
     real(r8)         , intent(in)   , pointer :: ciso_state(:)   ! isoC state, upstream pool
     real(r8)         , intent(in)   , pointer :: ctot_state(:)   ! totC state, upstream pool
     real(r8)         , intent(in)             :: frax_c13        ! fractionation factor (1 = no fractionation) for C13
     integer          , intent(in)             :: num             ! number of filter members
     integer          , intent(in)             :: filter(:)       ! filter indices
     integer          , intent(in)             :: diag            ! 0=no diagnostics, 1=print diagnostics
     character(len=*) , intent(in)             :: isotope         ! 'c13' or 'c14'
     !
     ! !LOCAL VARIABLES:
     integer :: beg2d, end2d
     integer :: i
     real(r8), pointer :: ciso_flux_1d(:)
     real(r8), pointer :: ctot_flux_1d(:)

     character(len=*), parameter :: subname = 'CIsoFluxCalc2d'
     !-----------------------------------------------------------------------

     SHR_ASSERT_ALL_FL((lbound(ctot_flux) == lbound(ciso_flux)), sourcefile, __LINE__)
     SHR_ASSERT_ALL_FL((ubound(ctot_flux) == ubound(ciso_flux)), sourcefile, __LINE__)

     beg2d = lbound(ciso_flux, 2)
     end2d = ubound(ciso_flux, 2)

     do i = beg2d, end2d
        ciso_flux_1d => ciso_flux(:,i)
        ctot_flux_1d => ctot_flux(:,i)
        call CIsoFluxCalc1d(&
             ciso_flux  = ciso_flux_1d, &
             ctot_flux  = ctot_flux_1d, &
             ciso_state = ciso_state, &
             ctot_state = ctot_state, &
             num        = num, &
             filter     = filter, &
             frax_c13   = frax_c13, &
             diag       = diag, &
             isotope    = isotope)
     end do

   end subroutine CIsoFluxCalc2dFlux


   !-----------------------------------------------------------------------
   subroutine CIsoFluxCalc2dBoth(&
        ciso_flux, ctot_flux, &
        ciso_state, ctot_state, &
        num, filter, frax_c13, diag, isotope)
     !
     ! !DESCRIPTION:
     ! Wrapper to CIsoFluxCalc1d for both the flux and state being 2-d variables
     !
     ! Loops over the second dimension of each variable to do a C Iso flux calc on each level
     !
     ! !ARGUMENTS:
     real(r8)         , intent(inout), pointer :: ciso_flux(:,:)  ! isoC flux
     real(r8)         , intent(in)   , pointer :: ctot_flux(:,:)  ! totC flux
     real(r8)         , intent(in)   , pointer :: ciso_state(:,:) ! isoC state, upstream pool
     real(r8)         , intent(in)   , pointer :: ctot_state(:,:) ! totC state, upstream pool
     real(r8)         , intent(in)             :: frax_c13        ! fractionation factor (1 = no fractionation) for C13
     integer          , intent(in)             :: num             ! number of filter members
     integer          , intent(in)             :: filter(:)       ! filter indices
     integer          , intent(in)             :: diag            ! 0=no diagnostics, 1=print diagnostics
     character(len=*) , intent(in)             :: isotope         ! 'c13' or 'c14'
     !
     ! !LOCAL VARIABLES:
     integer :: beg2d, end2d
     integer :: i
     real(r8), pointer :: ciso_flux_1d(:)
     real(r8), pointer :: ctot_flux_1d(:)
     real(r8), pointer :: ciso_state_1d(:)
     real(r8), pointer :: ctot_state_1d(:)

     character(len=*), parameter :: subname = 'CIsoFluxCalc2d'
     !-----------------------------------------------------------------------

     SHR_ASSERT_ALL_FL((lbound(ctot_flux) == lbound(ciso_flux)), sourcefile, __LINE__)
     SHR_ASSERT_ALL_FL((ubound(ctot_flux) == ubound(ciso_flux)), sourcefile, __LINE__)
     ! Note that we do NOT compare the state and flux bounds: it is okay for the state
     ! variables to have wider bounds than the flux variables (e.g., the flux variables
     ! can apply only for the grain components, whereas the state variables apply over
     ! all reproductive components).
     SHR_ASSERT_ALL_FL((lbound(ctot_state) == lbound(ciso_state)), sourcefile, __LINE__)
     SHR_ASSERT_ALL_FL((ubound(ctot_state) == ubound(ciso_state)), sourcefile, __LINE__)

     beg2d = lbound(ciso_flux, 2)
     end2d = ubound(ciso_flux, 2)

     do i = beg2d, end2d
        ciso_flux_1d => ciso_flux(:,i)
        ctot_flux_1d => ctot_flux(:,i)
        ciso_state_1d => ciso_state(:,i)
        ctot_state_1d => ctot_state(:,i)
        call CIsoFluxCalc1d(&
             ciso_flux  = ciso_flux_1d, &
             ctot_flux  = ctot_flux_1d, &
             ciso_state = ciso_state_1d, &
             ctot_state = ctot_state_1d, &
             num        = num, &
             filter     = filter, &
             frax_c13   = frax_c13, &
             diag       = diag, &
             isotope    = isotope)
     end do

   end subroutine CIsoFluxCalc2dBoth

end module CNCIsoFluxMod
