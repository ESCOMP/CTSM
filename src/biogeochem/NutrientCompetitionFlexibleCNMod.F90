module NutrientCompetitionFlexibleCNMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! DESCRIPTION
  ! module contains different subroutines to do soil nutrient competition dynamics
  !
  ! FIXME(bja, 2015-08) This module was copied from
  ! NutrientCompetitionCLM45default then flexible cn modifications
  ! were added for the clm50 nitrogen science changes (r120). There is
  ! a significant amount of duplicate code between the two
  ! modules. They need to be reexamined and the common code pulled out
  ! into a common base class.
  !
  ! created by Jinyun Tang, Sep 8, 2014
  ! modified by Mariana Vertenstein, Nov 15, 2014
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use clm_time_manager    , only : get_step_size_real
  use decompMod           , only : bounds_type
  use LandunitType        , only : lun
  use ColumnType          , only : col
  use PatchType           , only : patch
  use pftconMod           , only : pftcon, npcropmin
  use NutrientCompetitionMethodMod, only : nutrient_competition_method_type
  use CropReprPoolsMod    , only : nrepr
  use CNPhenologyMod      , only : CropPhase
  use CropType            , only : cphase_leafemerge, cphase_grainfill
  use clm_varctl          , only : use_crop_agsys
  use CNSharedParamsMod   , only : use_matrixcn
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  public :: nutrient_competition_FlexibleCN_type
  !
  type, extends(nutrient_competition_method_type) :: nutrient_competition_FlexibleCN_type
     private
     real(r8), pointer :: actual_leafcn(:)                    ! leaf CN ratio used by flexible CN
     real(r8), pointer :: actual_storage_leafcn(:)            ! storage leaf CN ratio used by flexible CN
   contains
     ! public methocs
     procedure, public :: Init                                ! Initialization
     procedure, public :: calc_plant_nutrient_competition     ! calculate nutrient yield rate from competition
     procedure, public :: calc_plant_nutrient_demand          ! calculate plant nutrient demand
     !
     ! private methods
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: calc_plant_cn_alloc
     procedure, private :: calc_plant_nitrogen_demand
  end type nutrient_competition_FlexibleCN_type
  !
  interface nutrient_competition_FlexibleCN_type
     ! initialize a new nutrient_competition_FlexibleCN_type object
     module procedure constructor
  end interface nutrient_competition_FlexibleCN_type
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: calc_npool_to_components_flexiblecn  ! Calculate npool_to_* terms for a single patch using the FlexibleCN approach
  private :: calc_npool_to_components_agsys       ! Calculate npool_to_* terms for a single crop patch when using AgSys

  ! !PRIVATE DATA:
  logical,parameter :: matrixcheck_ph = .True.
  logical,parameter :: acc_ph = .False.

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  type(nutrient_competition_FlexibleCN_type) function constructor()
    !
    ! !DESCRIPTION:
    ! Creates an object of type nutrient_competition_FlexibleCN_type.
    ! For now, this is simply a place-holder.
  end function constructor

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize the class
    !
    class(nutrient_competition_FlexibleCN_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Allocate memory for the class data
    !
    ! !USES:
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    ! !ARGUMENTS:
    class(nutrient_competition_FlexibleCN_type) :: this
    type(bounds_type), intent(in) :: bounds

    allocate(this%actual_leafcn(bounds%begp:bounds%endp))         ; this%actual_leafcn(:)         = nan
    allocate(this%actual_storage_leafcn(bounds%begp:bounds%endp)) ; this%actual_storage_leafcn(:) = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Send data to history file
    !
    ! !USES:
    use histFileMod    , only : hist_addfld1d
    use clm_varcon     , only : spval
    !
    ! !ARGUMENTS:
    class(nutrient_competition_FlexibleCN_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp

    this%actual_leafcn(begp:endp) = spval
    call hist_addfld1d (fname='LEAFCN', units='gC/gN', &
         avgflag='A', long_name='Leaf CN ratio used for flexible CN', &
         ptr_patch=this%actual_leafcn )
    this%actual_storage_leafcn(begp:endp) = spval
    call hist_addfld1d (fname='LEAFCN_STORAGE', units='gC/gN', &
         avgflag='A', long_name='Storage Leaf CN ratio used for flexible CN', &
         ptr_patch=this%actual_storage_leafcn, default='inactive')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine calc_plant_nutrient_competition (this, &
       bounds, num_soilp, filter_soilp, &
       cnveg_state_inst, crop_inst, canopystate_inst, cnveg_carbonstate_inst, &
       cnveg_carbonflux_inst, &
       c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst, &
       cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
       soilbiogeochem_nitrogenstate_inst, fpg_col)
    !
    ! !USES:
    use CNVegStateType        , only : cnveg_state_type
    use CropType              , only : crop_type
    use CanopyStateType        , only : canopystate_type
    use CNVegCarbonStateType   , only : cnveg_carbonstate_type
    use CNVegCarbonFluxType   , only : cnveg_carbonflux_type
    use CNVegNitrogenStateType , only : cnveg_nitrogenstate_type
    use CNVegNitrogenFluxType , only : cnveg_nitrogenflux_type
    use SoilBiogeochemNitrogenStateType, only : soilbiogeochem_nitrogenstate_type
    !
    ! !ARGUMENTS:
    class(nutrient_competition_FlexibleCN_type) , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_state_type)          , intent(inout) :: cnveg_state_inst
    type(crop_type)                 , intent(in)    :: crop_inst
    type(canopystate_type)          , intent(in)    :: canopystate_inst
    type(cnveg_carbonstate_type)    , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: c13_cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: c14_cnveg_carbonflux_inst
    type(cnveg_nitrogenstate_type)  , intent(inout) :: cnveg_nitrogenstate_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type), intent(in)    :: soilbiogeochem_nitrogenstate_inst
    real(r8), intent(in)    :: fpg_col (bounds%begc:)

    call this%calc_plant_cn_alloc(bounds, num_soilp, filter_soilp,   &
         cnveg_state_inst, crop_inst, canopystate_inst, &
         cnveg_carbonstate_inst, cnveg_carbonflux_inst, c13_cnveg_carbonflux_inst, &
         c14_cnveg_carbonflux_inst, cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
         soilbiogeochem_nitrogenstate_inst, &
         fpg_col=fpg_col(bounds%begc:bounds%endc))

  end subroutine calc_plant_nutrient_competition

!-----------------------------------------------------------------------
  subroutine calc_plant_cn_alloc(this, bounds, num_soilp, filter_soilp,   &
       cnveg_state_inst, crop_inst, canopystate_inst, &
       cnveg_carbonstate_inst, cnveg_carbonflux_inst, c13_cnveg_carbonflux_inst, &
       c14_cnveg_carbonflux_inst, cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
       soilbiogeochem_nitrogenstate_inst, fpg_col)
    !
    ! !USES:
    use clm_varctl            , only : use_c13, use_c14, carbon_resp_opt
    use CNVegStateType        , only : cnveg_state_type
    use CropType              , only : crop_type
    use CanopyStateType        , only : canopystate_type
    use CNVegCarbonStateType   , only : cnveg_carbonstate_type
    use CNVegCarbonFluxType   , only : cnveg_carbonflux_type
    use CNVegNitrogenFluxType , only : cnveg_nitrogenflux_type
    use CNVegNitrogenStateType , only : cnveg_nitrogenstate_type
    use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
    use CNSharedParamsMod     , only : use_fun
    use CNPrecisionControlMod , only : n_min
    use clm_varcon            , only : spval
    !index for matrixcn
    use clm_varpar            , only : ileaf,ileaf_st,ileaf_xf,ifroot,ifroot_st,ifroot_xf,&
                                       ilivestem,ilivestem_st,ilivestem_xf,&
                                       ideadstem,ideadstem_st,ideadstem_xf,&
                                       ilivecroot,ilivecroot_st,ilivecroot_xf,&
                                       ideadcroot,ideadcroot_st,ideadcroot_xf,&
                                       igrain,igrain_st,igrain_xf,iretransn,ioutc,ioutn,nvegnpool
    use CNVegMatrixMod        , only : matrix_update_phn
    !
    ! !ARGUMENTS:
    class(nutrient_competition_FlexibleCN_type), intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_soilp        ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:)  ! filter for soil patches
    type(cnveg_state_type)          , intent(inout) :: cnveg_state_inst
    type(crop_type)                 , intent(in)    :: crop_inst
    type(canopystate_type)          , intent(in)    :: canopystate_inst
    type(cnveg_carbonstate_type)    , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: c13_cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: c14_cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    type(cnveg_nitrogenstate_type)  , intent(inout) :: cnveg_nitrogenstate_inst
    type(soilbiogeochem_nitrogenstate_type), intent(in) :: soilbiogeochem_nitrogenstate_inst
    real(r8)                        , intent(in)    :: fpg_col(bounds%begc:)
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,k              ! indices
    integer  :: fp                 ! lake filter patch index
    real(r8) :: f1,f2,f3,f4,g1,g2  ! allocation parameters
    real(r8) :: fcur               ! fraction of current psn displayed as growth
    real(r8) :: gresp_storage      ! temporary variable for growth resp to storage
    real(r8) :: nlc                ! temporary variable for total new leaf carbon allocation
    real(r8) :: f5(nrepr)          ! reproductive allocation parameters
    real(r8) :: dt                 ! model time step

    real(r8):: frootcn_storage_actual
    real(r8):: frootcn_actual
    real(r8):: livestemcn_storage_actual
    real(r8):: livestemcn_actual
    real(r8):: livecrootcn_storage_actual
    real(r8):: livecrootcn_actual
    real(r8):: leafcn_max
    real(r8):: frootcn_max
    real(r8):: livewdcn_max
    real(r8):: frac_resp
    real(r8):: npool_to_veg
    real(r8):: cpool_to_veg
    real(r8) :: tmp

    ! -----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(fpg_col) == (/bounds%endc/)) , sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(this%actual_storage_leafcn) >= (/bounds%endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((lbound(this%actual_storage_leafcn) <= (/bounds%begp/)), sourcefile, __LINE__)

    associate(                                                                                       &
         fpg                          => fpg_col                                                   , & ! Input:  [real(r8) (:)   ]  fraction of potential gpp (no units)

         ivt                          => patch%itype                                               , & ! Input:  [integer  (:) ]  patch vegetation type

         woody                        => pftcon%woody                                              , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)
         froot_leaf                   => pftcon%froot_leaf                                         , & ! Input:  allocation parameter: new fine root C per new leaf C (gC/gC)
         croot_stem                   => pftcon%croot_stem                                         , & ! Input:  allocation parameter: new coarse root C per new stem C (gC/gC)
         stem_leaf                    => pftcon%stem_leaf                                          , & ! Input:  allocation parameter: new stem c per new leaf C (gC/gC)
         flivewd                      => pftcon%flivewd                                            , & ! Input:  allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
         leafcn                       => pftcon%leafcn                                             , & ! Input:  leaf C:N (gC/gN)
         frootcn                      => pftcon%frootcn                                            , & ! Input:  fine root C:N (gC/gN)
         livewdcn                     => pftcon%livewdcn                                           , & ! Input:  live wood (phloem and ray parenchyma) C:N (gC/gN)
         fcur2                        => pftcon%fcur                                               , & ! Input:  allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
         grperc                       => pftcon%grperc                                             , & ! Input:  growth respiration parameter
         grpnow                       => pftcon%grpnow                                             , & ! Input:  growth respiration parameter
         evergreen                    => pftcon%evergreen                                          , & ! Input:  binary flag for evergreen leaf habit (0 or 1)

         croplive                     => crop_inst%croplive_patch                                  , & ! Input:  [logical  (:)   ]  flag, true if planted, not harvested

         peaklai                      => cnveg_state_inst%peaklai_patch                            , & ! Input:  [integer  (:)   ]  1: max allowed lai; 0: not at max
         aleaf                        => cnveg_state_inst%aleaf_patch                              , & ! Input:  [real(r8) (:)   ]  leaf allocation coefficient
         astem                        => cnveg_state_inst%astem_patch                              , & ! Input:  [real(r8) (:)   ]  stem allocation coefficient
         aroot                        => cnveg_state_inst%aroot_patch                              , & ! Input:  [real(r8) (:)   ]  root allocation coefficient
         arepr                        => cnveg_state_inst%arepr_patch                              , & ! Input:  [real(r8) (:,:) ]  reproductive allocation coefficient(s)
         ! aleaf_n, astem_n, aroot_n and arepr_n are also inputs when running with AgSys,
         ! but they cannot be associated here because these pointers may be unallocated
         ! if not running with AgSys
         c_allometry                  => cnveg_state_inst%c_allometry_patch                        , & ! Output: [real(r8) (:)   ]  C allocation index (DIM)

         annsum_npp                   => cnveg_carbonflux_inst%annsum_npp_patch                    , & ! Input:  [real(r8) (:)   ]  annual sum of NPP, for wood allocation
         availc                       => cnveg_carbonflux_inst%availc_patch                        , & ! Output: [real(r8) (:)   ]  C flux available for allocation (gC/m2/s)
         plant_calloc                 => cnveg_carbonflux_inst%plant_calloc_patch                  , & ! Output: [real(r8) (:)   ]  total allocated C flux (gC/m2/s)
         npp_growth                   => cnveg_carbonflux_inst%npp_growth_patch                    , & ! output:  [real(r8) (:) ] c for growth in fun. g/m2/s
         cpool_to_resp                => cnveg_carbonflux_inst%cpool_to_resp_patch                 , & ! output: [real(r8) (:)   ]
         cpool_to_leafc_resp          => cnveg_carbonflux_inst%cpool_to_leafc_resp_patch           , & ! Output: [real(r8) (:)   ]
         cpool_to_leafc_storage_resp  => cnveg_carbonflux_inst%cpool_to_leafc_storage_resp_patch   , & ! Output: [real(r8) (:)   ]
         cpool_to_frootc_resp         => cnveg_carbonflux_inst%cpool_to_frootc_resp_patch          , & ! Output: [real(r8) (:)   ]
         cpool_to_frootc_storage_resp => cnveg_carbonflux_inst%cpool_to_frootc_storage_resp_patch  , & ! Output: [real(r8) (:)   ]
         cpool_to_livecrootc_resp     => cnveg_carbonflux_inst%cpool_to_livecrootc_resp_patch      , & ! Output: [real(r8) (:)   ]
         cpool_to_livecrootc_storage_resp  => cnveg_carbonflux_inst%cpool_to_livecrootc_storage_resp_patch , & ! Output: [real(r8) (:)   ]
         cpool_to_livestemc_resp      => cnveg_carbonflux_inst%cpool_to_livestemc_resp_patch       , & ! Output: [real(r8) (:)   ]
         cpool_to_livestemc_storage_resp => cnveg_carbonflux_inst%cpool_to_livestemc_storage_resp_patch  , & ! Output: [real(r8) (:)   ]
         cpool_to_leafc               => cnveg_carbonflux_inst%cpool_to_leafc_patch                , & ! Output: [real(r8) (:)   ]
         cpool_to_leafc_storage       => cnveg_carbonflux_inst%cpool_to_leafc_storage_patch        , & ! Output: [real(r8) (:)   ]
         cpool_to_frootc              => cnveg_carbonflux_inst%cpool_to_frootc_patch               , & ! Output: [real(r8) (:)   ]
         cpool_to_frootc_storage      => cnveg_carbonflux_inst%cpool_to_frootc_storage_patch       , & ! Output: [real(r8) (:)   ]
         cpool_to_livestemc           => cnveg_carbonflux_inst%cpool_to_livestemc_patch            , & ! Output: [real(r8) (:)   ]
         cpool_to_livestemc_storage   => cnveg_carbonflux_inst%cpool_to_livestemc_storage_patch    , & ! Output: [real(r8) (:)   ]
         cpool_to_deadstemc           => cnveg_carbonflux_inst%cpool_to_deadstemc_patch            , & ! Output: [real(r8) (:)   ]
         cpool_to_deadstemc_storage   => cnveg_carbonflux_inst%cpool_to_deadstemc_storage_patch    , & ! Output: [real(r8) (:)   ]
         cpool_to_livecrootc          => cnveg_carbonflux_inst%cpool_to_livecrootc_patch           , & ! Output: [real(r8) (:)   ]
         cpool_to_livecrootc_storage  => cnveg_carbonflux_inst%cpool_to_livecrootc_storage_patch   , & ! Output: [real(r8) (:)   ]
         cpool_to_deadcrootc          => cnveg_carbonflux_inst%cpool_to_deadcrootc_patch           , & ! Output: [real(r8) (:)   ]
         cpool_to_deadcrootc_storage  => cnveg_carbonflux_inst%cpool_to_deadcrootc_storage_patch   , & ! Output: [real(r8) (:)   ]
         cpool_to_gresp_storage       => cnveg_carbonflux_inst%cpool_to_gresp_storage_patch        , & ! Output: [real(r8) (:)   ]  allocation to growth respiration storage (gC/m2/s)
         cpool_to_reproductivec              => cnveg_carbonflux_inst%cpool_to_reproductivec_patch               , & ! Output: [real(r8) (:,:)   ]  allocation to grain C (gC/m2/s)
         cpool_to_reproductivec_storage      => cnveg_carbonflux_inst%cpool_to_reproductivec_storage_patch       , & ! Output: [real(r8) (:,:)   ]  allocation to grain C storage (gC/m2/s)

         laisun                       => canopystate_inst%laisun_patch  , & ! Input:  [real(r8) (:)   ]  sunlit projected leaf area index
         laisha                       => canopystate_inst%laisha_patch  , & ! Input:  [real(r8) (:)   ]  shaded projected leaf area index
         smin_no3_vr                  => soilbiogeochem_nitrogenstate_inst%smin_no3_vr_col         , & ! Output: [real(r8) (:,:) ]  (gN/m3) soil mineral NO3
         leafn                        => cnveg_nitrogenstate_inst%leafn_patch                      , & ! Input:  [real(r8) (:)   ]  (gN/m2) leaf N
         leafn_storage                => cnveg_nitrogenstate_inst%leafn_storage_patch              , & ! Input:  [real(r8) (:)   ]  (gN/m2) leaf N
         npool                        => cnveg_nitrogenstate_inst%npool_patch                      , & ! Input:  [real(r8) (:)   ]  (gN/m2) temporary plant N pool
         plant_ndemand                => cnveg_nitrogenflux_inst%plant_ndemand_patch               , & ! Output: [real(r8) (:)   ]  N flux required to support initial GPP (gN/m2/s)
         plant_nalloc                 => cnveg_nitrogenflux_inst%plant_nalloc_patch                , & ! Output: [real(r8) (:)   ]  total allocated N flux (gN/m2/s)
         npool_to_reproductiven              => cnveg_nitrogenflux_inst%npool_to_reproductiven_patch             , & ! Output: [real(r8) (:,:)   ]  allocation to grain N (gN/m2/s)
         npool_to_reproductiven_storage      => cnveg_nitrogenflux_inst%npool_to_reproductiven_storage_patch     , & ! Output: [real(r8) (:,:)   ]  allocation to grain N storage (gN/m2/s)
         retransn_to_npool            => cnveg_nitrogenflux_inst%retransn_to_npool_patch           , & ! Output: [real(r8) (:)   ]  deployment of retranslocated N (gN/m2/s)
         retransn                     => cnveg_nitrogenstate_inst%retransn_patch           , & ! Input:  [real(r8) (:)   ]  (gN/m2) plant pool of retranslocated N
         sminn_to_npool               => cnveg_nitrogenflux_inst%sminn_to_npool_patch              , & ! Output: [real(r8) (:)   ]  deployment of soil mineral N uptake (gN/m2/s)
         npool_to_leafn               => cnveg_nitrogenflux_inst%npool_to_leafn_patch              , & ! Output: [real(r8) (:)   ]  allocation to leaf N (gN/m2/s)
         npool_to_leafn_storage       => cnveg_nitrogenflux_inst%npool_to_leafn_storage_patch      , & ! Output: [real(r8) (:)   ]  allocation to leaf N storage (gN/m2/s)
         npool_to_frootn              => cnveg_nitrogenflux_inst%npool_to_frootn_patch             , & ! Output: [real(r8) (:)   ]  allocation to fine root N (gN/m2/s)
         npool_to_frootn_storage      => cnveg_nitrogenflux_inst%npool_to_frootn_storage_patch     , & ! Output: [real(r8) (:)   ]  allocation to fine root N storage (gN/m2/s)
         npool_to_livestemn           => cnveg_nitrogenflux_inst%npool_to_livestemn_patch          , & ! Output: [real(r8) (:)   ]
         npool_to_livestemn_storage   => cnveg_nitrogenflux_inst%npool_to_livestemn_storage_patch  , & ! Output: [real(r8) (:)   ]
         npool_to_deadstemn           => cnveg_nitrogenflux_inst%npool_to_deadstemn_patch          , & ! Output: [real(r8) (:)   ]
         npool_to_deadstemn_storage   => cnveg_nitrogenflux_inst%npool_to_deadstemn_storage_patch  , & ! Output: [real(r8) (:)   ]
         npool_to_livecrootn          => cnveg_nitrogenflux_inst%npool_to_livecrootn_patch         , & ! Output: [real(r8) (:)   ]
         npool_to_livecrootn_storage  => cnveg_nitrogenflux_inst%npool_to_livecrootn_storage_patch , & ! Output: [real(r8) (:)   ]
         npool_to_deadcrootn          => cnveg_nitrogenflux_inst%npool_to_deadcrootn_patch         , & ! Output: [real(r8) (:)   ]
         npool_to_deadcrootn_storage  => cnveg_nitrogenflux_inst%npool_to_deadcrootn_storage_patch , & ! Output: [real(r8) (:)   ]
         Npassive                     => cnveg_nitrogenflux_inst%Npassive_patch                    , & ! Output:  [real(r8) (:) ]  Passive N uptake (gN/m2/s)
         Nfix                         => cnveg_nitrogenflux_inst%Nfix_patch                        , & ! Output:  [real(r8) (:) ]  Symbiotic BNF (gN/m2/s)
         Nactive                      => cnveg_nitrogenflux_inst%Nactive_patch                     , & ! Output:  [real(r8) (:) ] Mycorrhizal N uptake (gN/m2/s)
         Nnonmyc                      => cnveg_nitrogenflux_inst%Nnonmyc_patch                     , & ! Output:  [real(r8) (:) ] Non-mycorrhizal N uptake (gN/m2/s)
         Nam                          => cnveg_nitrogenflux_inst%Nam_patch                         , & ! Output:  [real(r8) (:) ]  AM uptake (gN/m2/s)
         Necm                         => cnveg_nitrogenflux_inst%Necm_patch                        , & ! Output:  [real(r8) (:) ]  ECM uptake (gN/m2/s)
         sminn_to_plant_fun           => cnveg_nitrogenflux_inst%sminn_to_plant_fun_patch          , & ! Output:  [real(r8) (:) ]  Total soil N uptake of FUN (gN/m2/s)

         iretransn_to_ileaf           => cnveg_nitrogenflux_inst%iretransn_to_ileaf_ph             , & ! Transfer index (from retranslocation pool to leaf pool)
         iretransn_to_ileafst         => cnveg_nitrogenflux_inst%iretransn_to_ileafst_ph           , & ! Transfer index (from retranslocation pool to leaf storage pool)
         iretransn_to_ifroot          => cnveg_nitrogenflux_inst%iretransn_to_ifroot_ph            , & ! Transfer index (from retranslocation pool to fine root pool)
         iretransn_to_ifrootst        => cnveg_nitrogenflux_inst%iretransn_to_ifrootst_ph          , & ! Transfer index (from retranslocation pool to fine root storage pool)
         iretransn_to_ilivestem       => cnveg_nitrogenflux_inst%iretransn_to_ilivestem_ph         , & ! Transfer index (from retranslocation pool to live stem pool)
         iretransn_to_ilivestemst     => cnveg_nitrogenflux_inst%iretransn_to_ilivestemst_ph       , & ! Transfer index (from retranslocation pool to live stem storage pool)
         iretransn_to_ideadstem       => cnveg_nitrogenflux_inst%iretransn_to_ideadstem_ph         , & ! Transfer index (from retranslocation pool to dead stem pool)
         iretransn_to_ideadstemst     => cnveg_nitrogenflux_inst%iretransn_to_ideadstemst_ph       , & ! Transfer index (from retranslocation pool to dead stem storage pool)
         iretransn_to_ilivecroot      => cnveg_nitrogenflux_inst%iretransn_to_ilivecroot_ph        , & ! Transfer index (from retranslocation pool to live coarse root pool)
         iretransn_to_ilivecrootst    => cnveg_nitrogenflux_inst%iretransn_to_ilivecrootst_ph      , & ! Transfer index (from retranslocation pool to live coarse root storage pool)
         iretransn_to_ideadcroot      => cnveg_nitrogenflux_inst%iretransn_to_ideadcroot_ph        , & ! Transfer index (from retranslocation pool to dead coarse root pool)
         iretransn_to_ideadcrootst    => cnveg_nitrogenflux_inst%iretransn_to_ideadcrootst_ph      , & ! Transfer index (from retranslocation pool to dead coarse root storage pool)
         iretransn_to_igrain          => cnveg_nitrogenflux_inst%iretransn_to_igrain_ph            , & ! Transfer index (from retranslocation pool to grain pool)
         iretransn_to_igrainst        => cnveg_nitrogenflux_inst%iretransn_to_igrainst_ph          , & ! Transfer index (from retranslocation pool to grain storage pool)
         iretransn_to_iout            => cnveg_nitrogenflux_inst%iretransn_to_iout_ph              , & ! Transfer index (from retranslocation pool to external)
         ileaf_to_iretransn           => cnveg_nitrogenflux_inst%ileaf_to_iretransn_ph             , & ! Transfer index (from leaf pool to retranslocation pools)
         ifroot_to_iretransn          => cnveg_nitrogenflux_inst%ifroot_to_iretransn_ph            , & ! Transfer index (from fine root pool to retranslocation pools)
         ilivestem_to_iretransn       => cnveg_nitrogenflux_inst%ilivestem_to_iretransn_ph           & ! Transfer index (from live stem pool to retranslocation pools)
         )

      ! set time steps
      dt = get_step_size_real()

      ! patch loop to distribute the available N between the competing patches
      ! on the basis of relative demand, and allocate C and N to new growth and storage

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = patch%column(p)

         ! set some local allocation variables
         f1 = froot_leaf(ivt(p))
         f2 = croot_stem(ivt(p))

         ! modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0,
         ! constrained so that it does not go lower than 0.2 (under negative annsum_npp)
         ! There was an error in this formula in previous version, where the coefficient
         ! was 0.004 instead of 0.0025.
         ! This variable allocation is only for trees. Shrubs have a constant
         ! allocation as specified in the pft-physiology file.  The value is also used
         ! as a trigger here: -1.0 means to use the dynamic allocation (trees).
         if (stem_leaf(ivt(p)) == -1._r8) then
            f3 = (2.7_r8/(1.0_r8+exp(-0.004_r8*(annsum_npp(p) - 300.0_r8)))) - 0.4_r8
         else
            f3 = stem_leaf(ivt(p))
         end if

         f4   = flivewd(ivt(p))
         g1   = grperc(ivt(p))
         g2   = grpnow(ivt(p))
         fcur = fcur2(ivt(p))

         if (evergreen(ivt(p)) == 1._r8) then
            fcur = 0.0_r8
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            if (croplive(p)) then
               f1 = aroot(p) / aleaf(p)
               f3 = astem(p) / aleaf(p)
               do k = 1, nrepr
                  f5(k) = arepr(p,k) / aleaf(p)
               end do
               g1 = 0.25_r8
            else
               f1 = 0._r8
               f3 = 0._r8
               do k = 1, nrepr
                  f5(k) = 0._r8
               end do
               g1 = 0.25_r8
            end if
         end if

         ! increase fcur linearly with ndays_active, until fcur reaches 1.0 at
         ! ndays_active = days/year.  This prevents the continued storage of C and N.
         ! turning off this correction (PET, 12/11/03), instead using bgtr in
         ! phenology algorithm.

         if(use_fun)then ! if we are using FUN, we get the N available from there.
           sminn_to_npool(p) = sminn_to_plant_fun(p)
         else ! no FUN. :( we get N available from the FPG calculation in soilbiogeochemistry competition.
           sminn_to_npool(p) = plant_ndemand(p) * fpg(c)
         endif

         plant_nalloc(p) = sminn_to_npool(p) + retransn_to_npool(p)
         if(use_matrixcn)then
            associate( &
              matrix_Ninput => cnveg_nitrogenflux_inst%matrix_Ninput_patch & ! N input of matrix
            )
            matrix_Ninput(p) =  sminn_to_npool(p)
            end associate
         end if

         if(use_fun)then
            plant_calloc(p)  = npp_growth(p)
            if(use_matrixcn)then
               cnveg_carbonflux_inst%matrix_Cinput_patch(p) = npp_growth(p)
            end if
         else
            plant_calloc(p)  = availc(p)
            if(use_matrixcn)then
               cnveg_carbonflux_inst%matrix_Cinput_patch(p) = availc(p)
            end if
         end if

         ! calculate the amount of new leaf C dictated by these allocation
         ! decisions, and calculate the daily fluxes of C and N to current
         ! growth and storage pools

         ! fcur is the proportion of this day's growth that is displayed now,
         ! the remainder going into storage for display next year through the
         ! transfer pools

         nlc = plant_calloc(p) / c_allometry(p)
         cpool_to_leafc(p)          = nlc * fcur
         cpool_to_leafc_storage(p)  = nlc * (1._r8 - fcur)
         cpool_to_frootc(p)         = nlc * f1 * fcur
         cpool_to_frootc_storage(p) = nlc * f1 * (1._r8 - fcur)
         if(use_matrixcn)then
            cpool_to_veg = cpool_to_leafc(p) + cpool_to_leafc_storage(p) &
                         + cpool_to_frootc(p) + cpool_to_frootc_storage(p)
         end if
         if (woody(ivt(p)) == 1._r8) then
            cpool_to_livestemc(p)          = nlc * f3 * f4 * fcur
            cpool_to_livestemc_storage(p)  = nlc * f3 * f4 * (1._r8 - fcur)
            cpool_to_deadstemc(p)          = nlc * f3 * (1._r8 - f4) * fcur
            cpool_to_deadstemc_storage(p)  = nlc * f3 * (1._r8 - f4) * (1._r8 - fcur)
            cpool_to_livecrootc(p)         = nlc * f2 * f3 * f4 * fcur
            cpool_to_livecrootc_storage(p) = nlc * f2 * f3 * f4 * (1._r8 - fcur)
            cpool_to_deadcrootc(p)         = nlc * f2 * f3 * (1._r8 - f4) * fcur
            cpool_to_deadcrootc_storage(p) = nlc * f2 * f3 * (1._r8 - f4) * (1._r8 - fcur)
            if(use_matrixcn)then
               cpool_to_veg = cpool_to_veg &
                            + cpool_to_livestemc(p)  + cpool_to_livestemc_storage(p) &
                            + cpool_to_deadstemc(p)  + cpool_to_deadstemc_storage(p) &
                            + cpool_to_livecrootc(p) + cpool_to_livecrootc_storage(p) &
                            + cpool_to_deadcrootc(p) + cpool_to_deadcrootc_storage(p)
            end if
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            cpool_to_livestemc(p)          = nlc * f3 * f4 * fcur
            cpool_to_livestemc_storage(p)  = nlc * f3 * f4 * (1._r8 - fcur)
            cpool_to_deadstemc(p)          = nlc * f3 * (1._r8 - f4) * fcur
            cpool_to_deadstemc_storage(p)  = nlc * f3 * (1._r8 - f4) * (1._r8 - fcur)
            cpool_to_livecrootc(p)         = nlc * f2 * f3 * f4 * fcur
            cpool_to_livecrootc_storage(p) = nlc * f2 * f3 * f4 * (1._r8 - fcur)
            cpool_to_deadcrootc(p)         = nlc * f2 * f3 * (1._r8 - f4) * fcur
            cpool_to_deadcrootc_storage(p) = nlc * f2 * f3 * (1._r8 - f4) * (1._r8 - fcur)
            do k = 1, nrepr
               cpool_to_reproductivec(p,k)         = nlc * f5(k) * fcur
               cpool_to_reproductivec_storage(p,k) = nlc * f5(k) * (1._r8 -fcur)
            end do
            if(use_matrixcn)then
               cpool_to_veg = cpool_to_veg &
                            + cpool_to_livestemc(p)  + cpool_to_livestemc_storage(p) &
                            + cpool_to_deadstemc(p)  + cpool_to_deadstemc_storage(p) &
                            + cpool_to_livecrootc(p) + cpool_to_livecrootc_storage(p) &
                            + cpool_to_deadcrootc(p) + cpool_to_deadcrootc_storage(p)
               do k = 1, nrepr
                  cpool_to_veg = cpool_to_veg &
                               + cpool_to_reproductivec(p,k) + cpool_to_reproductivec_storage(p,k)
               end do
            end if
         end if

         if (use_matrixcn) then
            associate( &
               matrix_Cinput => cnveg_carbonflux_inst%matrix_Cinput_patch, & ! C input of matrix 
               matrix_alloc  => cnveg_carbonflux_inst%matrix_alloc_patch   & ! B-matrix for carbon allocation
            )
            matrix_Cinput(p) = cpool_to_veg
            if(cpool_to_veg .ne. 0)then
               matrix_alloc(p,ileaf)     = cpool_to_leafc(p)  / cpool_to_veg
               matrix_alloc(p,ileaf_st)  = cpool_to_leafc_storage(p)  / cpool_to_veg
               matrix_alloc(p,ifroot)    = cpool_to_frootc(p) / cpool_to_veg
               matrix_alloc(p,ifroot_st) = cpool_to_frootc_storage(p) / cpool_to_veg
            end if

            if (woody(ivt(p)) == 1._r8) then
               if(cpool_to_veg .ne. 0)then
                  matrix_alloc(p,ilivestem)     = cpool_to_livestemc(p)  / cpool_to_veg
                  matrix_alloc(p,ilivestem_st)  = cpool_to_livestemc_storage(p)  / cpool_to_veg
                  matrix_alloc(p,ideadstem)     = cpool_to_deadstemc(p)  / cpool_to_veg
                  matrix_alloc(p,ideadstem_st)  = cpool_to_deadstemc_storage(p)  / cpool_to_veg
                  matrix_alloc(p,ilivecroot)    = cpool_to_livecrootc(p)  / cpool_to_veg
                  matrix_alloc(p,ilivecroot_st) = cpool_to_livecrootc_storage(p)  / cpool_to_veg
                  matrix_alloc(p,ideadcroot)    = cpool_to_deadcrootc(p)  / cpool_to_veg
                  matrix_alloc(p,ideadcroot_st) = cpool_to_deadcrootc_storage(p)  / cpool_to_veg
               end if
            end if
            if (ivt(p) >= npcropmin) then ! skip 2 generic crops
               if(cpool_to_veg .ne. 0)then
                  matrix_alloc(p,ilivestem)     = cpool_to_livestemc(p)  / cpool_to_veg
                  matrix_alloc(p,ilivestem_st)  = cpool_to_livestemc_storage(p)  / cpool_to_veg
                  matrix_alloc(p,ideadstem)     = cpool_to_deadstemc(p)  / cpool_to_veg
                  matrix_alloc(p,ideadstem_st)  = cpool_to_deadstemc_storage(p)  / cpool_to_veg
                  matrix_alloc(p,ilivecroot)    = cpool_to_livecrootc(p) / cpool_to_veg
                  matrix_alloc(p,ilivecroot_st) = cpool_to_livecrootc_storage(p) / cpool_to_veg
                  matrix_alloc(p,ideadcroot)    = cpool_to_deadcrootc(p) / cpool_to_veg
                  matrix_alloc(p,ideadcroot_st) = cpool_to_deadcrootc_storage(p) / cpool_to_veg
                  matrix_alloc(p,igrain)        = 0.0_r8
                  matrix_alloc(p,igrain_st)     = 0.0_r8
                  do k = 1, nrepr
                     matrix_alloc(p,igrain)        = matrix_alloc(p,igrain) + cpool_to_reproductivec(p,k) / cpool_to_veg
                     matrix_alloc(p,igrain_st)     = matrix_alloc(p,igrain_st) + cpool_to_reproductivec_storage(p,k) / cpool_to_veg
                  end do
              end if
           end if
           end associate
         end if !use_matrixcn

         ! Calculate the amount of carbon that needs to go into growth
         ! respiration storage to satisfy all of the storage growth demands.
         ! Allows for the fraction of growth respiration that is released at the
         ! time of fixation, versus the remaining fraction that is stored for
         ! release at the time of display. Note that all the growth respiration
         ! fluxes that get released on a given timestep are calculated in growth_resp(),
         ! but that the storage of C for growth resp during display of transferred
         ! growth is assigned here.

         gresp_storage = cpool_to_leafc_storage(p) + cpool_to_frootc_storage(p)
         if (woody(ivt(p)) == 1._r8) then
            gresp_storage = gresp_storage + cpool_to_livestemc_storage(p)
            gresp_storage = gresp_storage + cpool_to_deadstemc_storage(p)

            gresp_storage = gresp_storage + cpool_to_livecrootc_storage(p)
            gresp_storage = gresp_storage + cpool_to_deadcrootc_storage(p)
         end if
         if (ivt(p) >= npcropmin) then     ! skip 2 generic crops
            gresp_storage = gresp_storage + cpool_to_livestemc_storage(p)
            do k = 1, nrepr
               gresp_storage = gresp_storage + cpool_to_reproductivec_storage(p,k)
            end do
         end if
         cpool_to_gresp_storage(p) = gresp_storage * g1 * (1._r8 - g2)

         if (use_crop_agsys .and. ivt(p) >= npcropmin) then
            call calc_npool_to_components_agsys( &
                 ! Inputs
                 npool = npool(p), &
                 fcur = fcur, &
                 f4 = f4, &
                 ! The following inputs cannot appear in the associate statement at the
                 ! top of the subroutine because these pointers may be unallocated if not
                 ! running with AgSys:
                 aleaf_n = cnveg_state_inst%aleaf_n_patch(p), &
                 astem_n = cnveg_state_inst%astem_n_patch(p), &
                 aroot_n = cnveg_state_inst%aroot_n_patch(p), &
                 arepr_n = cnveg_state_inst%arepr_n_patch(p,:), &

                 ! Outputs
                 npool_to_leafn = npool_to_leafn(p), &
                 npool_to_leafn_storage = npool_to_leafn_storage(p), &
                 npool_to_frootn = npool_to_frootn(p), &
                 npool_to_frootn_storage = npool_to_frootn_storage(p), &
                 npool_to_livestemn = npool_to_livestemn(p), &
                 npool_to_livestemn_storage = npool_to_livestemn_storage(p), &
                 npool_to_deadstemn = npool_to_deadstemn(p), &
                 npool_to_deadstemn_storage = npool_to_deadstemn_storage(p), &
                 npool_to_livecrootn = npool_to_livecrootn(p), &
                 npool_to_livecrootn_storage = npool_to_livecrootn_storage(p), &
                 npool_to_deadcrootn = npool_to_deadcrootn(p), &
                 npool_to_deadcrootn_storage = npool_to_deadcrootn_storage(p), &
                 npool_to_reproductiven = npool_to_reproductiven(p,:), &
                 npool_to_reproductiven_storage = npool_to_reproductiven_storage(p,:))
         else
            call calc_npool_to_components_flexiblecn( &
                 ! Inputs
                 npool = npool(p), &
                 ivt = ivt(p), &
                 nlc = nlc, &
                 fcur = fcur, &
                 f1 = f1, &
                 f2 = f2, &
                 f3 = f3, &
                 f4 = f4, &
                 f5 = f5, &

                 ! Outputs
                 npool_to_leafn = npool_to_leafn(p), &
                 npool_to_leafn_storage = npool_to_leafn_storage(p), &
                 npool_to_frootn = npool_to_frootn(p), &
                 npool_to_frootn_storage = npool_to_frootn_storage(p), &
                 npool_to_livestemn = npool_to_livestemn(p), &
                 npool_to_livestemn_storage = npool_to_livestemn_storage(p), &
                 npool_to_deadstemn = npool_to_deadstemn(p), &
                 npool_to_deadstemn_storage = npool_to_deadstemn_storage(p), &
                 npool_to_livecrootn = npool_to_livecrootn(p), &
                 npool_to_livecrootn_storage = npool_to_livecrootn_storage(p), &
                 npool_to_deadcrootn = npool_to_deadcrootn(p), &
                 npool_to_deadcrootn_storage = npool_to_deadcrootn_storage(p), &
                 npool_to_reproductiven = npool_to_reproductiven(p,:), &
                 npool_to_reproductiven_storage = npool_to_reproductiven_storage(p,:))
         end if

         cpool_to_resp(p) = 0.0_r8
         cpool_to_leafc_resp(p) = 0.0_r8
         cpool_to_leafc_storage_resp(p) = 0.0_r8
         cpool_to_frootc_resp(p) = 0.0_r8
         cpool_to_frootc_storage_resp(p) = 0.0_r8
         cpool_to_livecrootc_resp(p) = 0.0_r8
         cpool_to_livecrootc_storage_resp(p) = 0.0_r8
         cpool_to_livestemc_resp(p) = 0.0_r8
         cpool_to_livestemc_storage_resp(p) = 0.0_r8

         if ( laisun(p)+laisha(p) > 0.0_r8 ) then
            if (cnveg_nitrogenstate_inst%leafn_storage_patch(p) == 0.0_r8 ) then
               ! to avoid division by zero, and also to make actual_leafncn(p) a very large number if leafn(p) is zero
               this%actual_storage_leafcn(p) = spval
            else
               ! leaf CN ratio
               this%actual_storage_leafcn(p) = cnveg_carbonstate_inst%leafc_storage_patch(p)  &
                    / cnveg_nitrogenstate_inst%leafn_storage_patch(p)
            end if
         end if

         if (carbon_resp_opt == 1 .AND. laisun(p)+laisha(p) > 0.0_r8) then
            ! computing carbon to nitrogen ratio of different plant parts


            if (cnveg_nitrogenstate_inst%frootn_storage_patch(p) == 0.0_r8) then
               ! to avoid division by zero, and also to make frootcn_actual(p) a very large number if frootc(p) is zero
               frootcn_actual = cnveg_carbonstate_inst%frootc_storage_patch(p) / n_min
            else
               ! fine root CN ratio
               frootcn_actual = cnveg_carbonstate_inst%frootc_storage_patch(p) / cnveg_nitrogenstate_inst%frootn_storage_patch(p)
            end if

            if (woody(ivt(p)) == 1._r8) then

               if (cnveg_nitrogenstate_inst%livestemn_storage_patch(p) == 0.0_r8) then
                  ! to avoid division by zero, and also to make livestemcn_actual(p) a very large number if livestemc(p) is zero
                  livestemcn_actual = cnveg_carbonstate_inst%livestemc_storage_patch(p) / n_min
               else
                  ! live stem CN ratio
                  livestemcn_actual =  cnveg_carbonstate_inst%livestemc_storage_patch(p) / &
                       cnveg_nitrogenstate_inst%livestemn_storage_patch(p)
               end if

               if (cnveg_nitrogenstate_inst%livecrootn_storage_patch(p) == 0.0_r8) then
                  ! to avoid division by zero, and also to make livecrootcn_actual(p) a very large number if livecrootc(p) is zero
                  livecrootcn_actual = cnveg_carbonstate_inst%livecrootc_storage_patch(p) / n_min
               else
                  ! live coarse root CN ratio
                  livecrootcn_actual = cnveg_carbonstate_inst%livecrootc_storage_patch(p) / &
                       cnveg_nitrogenstate_inst%livecrootn_storage_patch(p)
               end if
            end if

            if (ivt(p) >= npcropmin) then ! skip 2 generic crops

               if (cnveg_nitrogenstate_inst%livestemn_storage_patch(p) == 0.0_r8) then
                  ! to avoid division by zero, and also to make livestemcn_actual(p) a very large number if livestemc(p) is zero
                  livestemcn_actual = cnveg_carbonstate_inst%livestemc_storage_patch(p) / n_min
               else
                  ! live stem CN ratio
                  livestemcn_actual =  cnveg_carbonstate_inst%livestemc_storage_patch(p) / &
                       cnveg_nitrogenstate_inst%livestemn_storage_patch(p)
               end if

               if (cnveg_nitrogenstate_inst%livecrootn_storage_patch(p) == 0.0_r8) then
                  ! to avoid division by zero, and also to make livecrootcn_actual(p) a very large number if livecrootc(p) is zero
                  livecrootcn_actual = cnveg_carbonstate_inst%livecrootc_storage_patch(p) / n_min
               else
                  ! live coarse root CN ratio
                  livecrootcn_actual = cnveg_carbonstate_inst%livecrootc_storage_patch(p) / &
                       cnveg_nitrogenstate_inst%livecrootn_storage_patch(p)
               end if
            end if

            leafcn_max = leafcn(ivt(p)) + 15.0_r8
            frootcn_max = frootcn(ivt(p)) + 15.0_r8

            ! Note that for high CN ratio stress the plant part does not retranslocate nitrogen as the plant part will need the N
            ! if high leaf CN ratio (i.e., high leaf C compared to N) then turnover extra C
            if (this%actual_storage_leafcn(p) > leafcn_max) then

               frac_resp =  (this%actual_storage_leafcn(p) - leafcn_max) / 10.0_r8
               frac_resp = min(1.0_r8, max(0.0_r8, frac_resp))

               cpool_to_leafc_resp(p)          = frac_resp * cpool_to_leafc(p)
               cpool_to_leafc_storage_resp(p)  = frac_resp * cpool_to_leafc_storage(p)

               !cpool_to_leafc(p) = cpool_to_leafc(p) - cpool_to_leafc_resp(p)
               !cpool_to_leafc_storage(p) = cpool_to_leafc_storage(p) - cpool_to_leafc_storage_resp(p)

            end if

            ! if high fine root CN ratio (i.e., high fine root C compared to N) then turnover extra C
            if (frootcn_actual > frootcn_max) then

               frac_resp =  (frootcn_actual - frootcn_max) / 10.0_r8
               frac_resp = min(1.0_r8, max(0.0_r8, frac_resp))

               cpool_to_frootc_resp(p)         = frac_resp * cpool_to_frootc(p)
               cpool_to_frootc_storage_resp(p) = frac_resp * cpool_to_frootc_storage(p)

               !cpool_to_frootc(p) = cpool_to_frootc(p) - cpool_to_frootc_resp(p)
               !cpool_to_frootc_storage(p) = cpool_to_frootc_storage(p) - cpool_to_frootc_storage_resp(p)

            end if

            if (woody(ivt(p)) == 1._r8) then

               livewdcn_max = livewdcn(ivt(p)) + 15.0_r8

               ! if high coarse root CN ratio (i.e., high coarse root C compared to N) then turnover extra C
               if (livecrootcn_actual > livewdcn_max) then

                  frac_resp =  (livecrootcn_actual - livewdcn_max) / 10.0_r8
                  frac_resp = min(1.0_r8, max(0.0_r8, frac_resp))

                  cpool_to_livecrootc_resp(p)         = frac_resp * cpool_to_livecrootc(p)
                  cpool_to_livecrootc_storage_resp(p) = frac_resp * cpool_to_livecrootc_storage(p)

                  !cpool_to_livecrootc(p) = cpool_to_livecrootc(p) - cpool_to_livecrootc_resp(p)
                  !cpool_to_livecrootc_storage(p) = cpool_to_livecrootc_storage(p) - cpool_to_livecrootc_storage_resp(p)

               end if

               ! if high stem CN ratio (i.e., high stem C compared to N) then turnover extra C
               if (livestemcn_actual > livewdcn_max) then

                  frac_resp =  (livestemcn_actual - livewdcn_max) / 10.0_r8
                  frac_resp = min(1.0_r8, max(0.0_r8, frac_resp))

                  cpool_to_livestemc_resp(p)          = frac_resp * cpool_to_livestemc(p)
                  cpool_to_livestemc_storage_resp(p)  = frac_resp * cpool_to_livestemc_storage(p)

                  !cpool_to_livestemc(p) = cpool_to_livestemc(p) - cpool_to_livestemc_resp(p)
                  !cpool_to_livestemc_storage(p) = cpool_to_livestemc_storage(p) - cpool_to_livestemc_storage_resp(p)

               end if

            end if

            if (ivt(p) >= npcropmin) then ! skip 2 generic crops

               livewdcn_max = livewdcn(ivt(p)) + 15.0_r8

               ! if high coarse root CN ratio (i.e., high coarse root C compared to N) then turnover extra C
               if (livecrootcn_actual > livewdcn_max) then

                  frac_resp =  (livecrootcn_actual - livewdcn_max) / 10.0_r8
                  frac_resp = min(1.0_r8, max(0.0_r8, frac_resp))

                  cpool_to_livecrootc_resp(p)         = frac_resp * cpool_to_livecrootc(p)
                  cpool_to_livecrootc_storage_resp(p) = frac_resp * cpool_to_livecrootc_storage(p)

                  !cpool_to_livecrootc(p) = cpool_to_livecrootc(p) - cpool_to_livecrootc_resp(p)
                  !cpool_to_livecrootc_storage(p) = cpool_to_livecrootc_storage(p) - cpool_to_livecrootc_storage_resp(p)

               end if

               ! if high stem CN ratio (i.e., high stem C compared to N) then turnover extra C
               if (livestemcn_actual > livewdcn_max) then

                  frac_resp =  (livestemcn_actual - livewdcn_max) / 10.0_r8
                  frac_resp = min(1.0_r8, max(0.0_r8, frac_resp))

                  cpool_to_livestemc_resp(p)          = frac_resp * cpool_to_livestemc(p)
                  cpool_to_livestemc_storage_resp(p)  = frac_resp * cpool_to_livestemc_storage(p)

                  !cpool_to_livestemc(p) = cpool_to_livestemc(p) - cpool_to_livestemc_resp(p)
                  !cpool_to_livestemc_storage(p) = cpool_to_livestemc_storage(p) - cpool_to_livestemc_storage_resp(p)

               end if

            end if

            cpool_to_resp(p) = cpool_to_leafc_resp(p) + cpool_to_leafc_storage_resp(p)  + cpool_to_frootc_resp(p) + &
                 cpool_to_frootc_storage_resp(p) + cpool_to_livecrootc_resp(p) + cpool_to_livecrootc_storage_resp(p) + &
                 cpool_to_livestemc_resp(p) + cpool_to_livestemc_storage_resp(p)

            if(use_matrixcn)then
               cnveg_carbonflux_inst%matrix_Cinput_patch(p) = cnveg_carbonflux_inst%matrix_Cinput_patch(p) - cpool_to_resp(p)
            end if

         end if   ! end of if (carbon_resp_opt == 1 .AND. laisun(p)+laisha(p) > 0.0_r8) then

         !if (cnveg_nitrogenstate_inst%leafn_storage_patch(p) < n_min .or. laisun(p)+laisha(p) <= 0.0_r8) then
         !! to make output on history missing value
         !this%actual_storage_leafcn(p) = spval
         !end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if(use_matrixcn)then
           associate( &
              matrix_Ninput     => cnveg_nitrogenflux_inst%matrix_Ninput_patch,  & ! N input of matrix
              matrix_nalloc     => cnveg_nitrogenflux_inst%matrix_nalloc_patch,  & ! B-matrix for nitrogen allocation
              psnsun_to_cpool   => cnveg_carbonflux_inst%psnsun_to_cpool_patch,  & ! 
              psnshade_to_cpool => cnveg_carbonflux_inst%psnshade_to_cpool_patch & ! 
           )
           if(use_c13 .and. psnsun_to_cpool(p)+psnshade_to_cpool(p).ne. 0._r8)then
               associate( &
                   matrix_C13input => cnveg_carbonflux_inst%matrix_C13input_patch & ! C13 input of matrix
               )
               matrix_C13input(p) = plant_calloc(p) * &
                                ((c13_cnveg_carbonflux_inst%psnsun_to_cpool_patch(p)+ c13_cnveg_carbonflux_inst%psnshade_to_cpool_patch(p))/ &
                                (psnsun_to_cpool(p)+psnshade_to_cpool(p)))
               end associate
           end if
           if(use_c14 .and. psnsun_to_cpool(p)+psnshade_to_cpool(p).ne. 0._r8)then
               associate( &
                  matrix_C14input => cnveg_carbonflux_inst%matrix_C14input_patch & ! C14 input of matrix
               )
               matrix_C14input(p) = plant_calloc(p) * &
                                ((c14_cnveg_carbonflux_inst%psnsun_to_cpool_patch(p)+ c14_cnveg_carbonflux_inst%psnshade_to_cpool_patch(p))/ &
                                (psnsun_to_cpool(p)+psnshade_to_cpool(p)))
               end associate
           end if
           npool_to_veg = npool_to_leafn(p) + npool_to_leafn_storage(p) &
                         + npool_to_frootn(p) + npool_to_frootn_storage(p) &
                         + npool_to_livestemn(p) + npool_to_livestemn_storage(p) &
                         + npool_to_deadstemn(p) + npool_to_deadstemn_storage(p) &
                         + npool_to_livecrootn(p) + npool_to_livecrootn_storage(p)  &
                         + npool_to_deadcrootn(p) + npool_to_deadcrootn_storage(p)   
           if (ivt(p) >= npcropmin)then
               npool_to_veg = npool_to_veg + npool_to_reproductiven(p,1) + npool_to_reproductiven_storage(p,1)
           end if
           if(npool_to_veg .ne. 0._r8)then
               matrix_nalloc(p,ileaf         ) = npool_to_leafn(p)              / npool_to_veg
               matrix_nalloc(p,ileaf_st      ) = npool_to_leafn_storage(p)      / npool_to_veg
               matrix_nalloc(p,ifroot        ) = npool_to_frootn(p)             / npool_to_veg
               matrix_nalloc(p,ifroot_st     ) = npool_to_frootn_storage(p)     / npool_to_veg
               matrix_nalloc(p,ilivestem     ) = npool_to_livestemn(p)          / npool_to_veg
               matrix_nalloc(p,ilivestem_st  ) = npool_to_livestemn_storage(p)  / npool_to_veg
               matrix_nalloc(p,ideadstem     ) = npool_to_deadstemn(p)          / npool_to_veg
               matrix_nalloc(p,ideadstem_st  ) = npool_to_deadstemn_storage(p)  / npool_to_veg
               matrix_nalloc(p,ilivecroot    ) = npool_to_livecrootn(p)         / npool_to_veg
               matrix_nalloc(p,ilivecroot_st ) = npool_to_livecrootn_storage(p) / npool_to_veg
               matrix_nalloc(p,ideadcroot    ) = npool_to_deadcrootn(p)         / npool_to_veg
               matrix_nalloc(p,ideadcroot_st ) = npool_to_deadcrootn_storage(p) / npool_to_veg
               if (ivt(p) >= npcropmin)then
                  matrix_nalloc(p,igrain     ) = npool_to_reproductiven(p,1)             / npool_to_veg 
                  matrix_nalloc(p,igrain_st  ) = npool_to_reproductiven_storage(p,1)     / npool_to_veg 
               end if
               matrix_Ninput(p) = npool_to_veg - retransn_to_npool(p)
           else
               if(retransn(p) .ne. 0._r8)then
                  retransn_to_npool(p) = retransn(p) * matrix_update_phn(p,iretransn_to_iout,retransn_to_npool(p)/retransn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,.True.)
               end if
           end if

           if(retransn(p) .ne. 0._r8)then
               tmp = matrix_update_phn(p,iretransn_to_ileaf             ,matrix_nalloc(p,ileaf )         * retransn_to_npool(p) / retransn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,.True.)
               tmp = matrix_update_phn(p,iretransn_to_ileafst           ,matrix_nalloc(p,ileaf_st )      * retransn_to_npool(p) / retransn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,.True.)
               tmp = matrix_update_phn(p,iretransn_to_ifroot            ,matrix_nalloc(p,ifroot )        * retransn_to_npool(p) / retransn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,.True.)
               tmp = matrix_update_phn(p,iretransn_to_ifrootst          ,matrix_nalloc(p,ifroot_st )     * retransn_to_npool(p) / retransn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,.True.)
               tmp = matrix_update_phn(p,iretransn_to_ilivestem         ,matrix_nalloc(p,ilivestem )     * retransn_to_npool(p) / retransn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,.True.)
               tmp = matrix_update_phn(p,iretransn_to_ilivestemst       ,matrix_nalloc(p,ilivestem_st )  * retransn_to_npool(p) / retransn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,.True.)
               tmp = matrix_update_phn(p,iretransn_to_ideadstem         ,matrix_nalloc(p,ideadstem )     * retransn_to_npool(p) / retransn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,.True.)
               tmp = matrix_update_phn(p,iretransn_to_ideadstemst       ,matrix_nalloc(p,ideadstem_st )  * retransn_to_npool(p) / retransn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,.True.)
               tmp = matrix_update_phn(p,iretransn_to_ilivecroot        ,matrix_nalloc(p,ilivecroot )    * retransn_to_npool(p) / retransn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,.True.)
               tmp = matrix_update_phn(p,iretransn_to_ilivecrootst      ,matrix_nalloc(p,ilivecroot_st ) * retransn_to_npool(p) / retransn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,.True.)
               tmp = matrix_update_phn(p,iretransn_to_ideadcroot        ,matrix_nalloc(p,ideadcroot )    * retransn_to_npool(p) / retransn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,.True.)
               tmp = matrix_update_phn(p,iretransn_to_ideadcrootst      ,matrix_nalloc(p,ideadcroot_st ) * retransn_to_npool(p) / retransn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,.True.)
               if(ivt(p) >= npcropmin)then
                 tmp = matrix_update_phn(p,iretransn_to_igrain   ,matrix_nalloc(p,igrain    ) * retransn_to_npool(p) / retransn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,.True.)
                 tmp = matrix_update_phn(p,iretransn_to_igrainst ,matrix_nalloc(p,igrain_st ) * retransn_to_npool(p) / retransn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,.True.)
               end if
           end if
           end associate
         end if !end use_matrixcn  
      end do ! end patch loop

    end associate

  end subroutine calc_plant_cn_alloc

  !-----------------------------------------------------------------------
  subroutine calc_npool_to_components_flexiblecn( &
       npool, ivt, nlc, fcur, f1, f2, f3, f4, f5, &
       npool_to_leafn, npool_to_leafn_storage, &
       npool_to_frootn, npool_to_frootn_storage, &
       npool_to_livestemn, npool_to_livestemn_storage, &
       npool_to_deadstemn, npool_to_deadstemn_storage, &
       npool_to_livecrootn, npool_to_livecrootn_storage, &
       npool_to_deadcrootn, npool_to_deadcrootn_storage, &
       npool_to_reproductiven, npool_to_reproductiven_storage)
    !
    ! !DESCRIPTION:
    ! Calculate npool_to_* terms for a single patch using the FlexibleCN approach
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: npool ! temporary plant N pool (gN/m2)
    integer, intent(in) :: ivt ! vegetation type
    real(r8), intent(in) :: nlc ! new leaf carbon allocation (gC/m2/s)
    real(r8), intent(in) :: fcur ! fraction of current psn displayed as growth
    real(r8), intent(in) :: f1 ! C allocation parameter - fine_root:leaf ratio
    real(r8), intent(in) :: f2 ! C allocation parameter - coarse_root:stem ratio
    real(r8), intent(in) :: f3 ! C allocation parameter - stem:leaf ratio
    real(r8), intent(in) :: f4 ! C allocation parameter - fraction of new wood that is live
    real(r8), intent(in) :: f5(:) ! C allocation parameter - repr:leaf ratio for each crop reproductive pool

    ! Each of the following output variables is in units of gN/m2/s; they are
    ! intent(inout) because some may remain unchanged in some circumstances.
    real(r8), intent(inout) :: npool_to_leafn
    real(r8), intent(inout) :: npool_to_leafn_storage
    real(r8), intent(inout) :: npool_to_frootn
    real(r8), intent(inout) :: npool_to_frootn_storage
    real(r8), intent(inout) :: npool_to_livestemn
    real(r8), intent(inout) :: npool_to_livestemn_storage
    real(r8), intent(inout) :: npool_to_deadstemn
    real(r8), intent(inout) :: npool_to_deadstemn_storage
    real(r8), intent(inout) :: npool_to_livecrootn
    real(r8), intent(inout) :: npool_to_livecrootn_storage
    real(r8), intent(inout) :: npool_to_deadcrootn
    real(r8), intent(inout) :: npool_to_deadcrootn_storage
    real(r8), intent(inout) :: npool_to_reproductiven(:)
    real(r8), intent(inout) :: npool_to_reproductiven_storage(:)

    !
    ! !LOCAL VARIABLES:
    real(r8) :: cnl,cnfr,cnlw,cndw ! C:N ratios for leaf, fine root, and wood
    real(r8) :: cng                ! C:N ratio for grain (= cnlw for now; slevis)
    real(r8) :: dt                 ! model time step
    integer  :: k

    real(r8) :: npool_to_reproductiven_demand_tot
    real(r8) :: npool_to_reproductiven_storage_demand_tot
    real(r8) :: npool_to_leafn_demand
    real(r8) :: npool_to_leafn_storage_demand
    real(r8) :: npool_to_frootn_demand
    real(r8) :: npool_to_frootn_storage_demand
    real(r8) :: npool_to_livestemn_demand
    real(r8) :: npool_to_livestemn_storage_demand
    real(r8) :: npool_to_livecrootn_demand
    real(r8) :: npool_to_livecrootn_storage_demand
    real(r8) :: npool_to_deadstemn_demand
    real(r8) :: npool_to_deadstemn_storage_demand
    real(r8) :: npool_to_deadcrootn_demand
    real(r8) :: npool_to_deadcrootn_storage_demand
    real(r8) :: npool_to_reproductiven_demand(nrepr)
    real(r8) :: npool_to_reproductiven_storage_demand(nrepr)
    real(r8) :: total_plant_Ndemand
    real(r8) :: frNdemand_npool_to_leafn
    real(r8) :: frNdemand_npool_to_leafn_storage
    real(r8) :: frNdemand_npool_to_frootn
    real(r8) :: frNdemand_npool_to_frootn_storage
    real(r8) :: frNdemand_npool_to_livestemn
    real(r8) :: frNdemand_npool_to_livestemn_storage
    real(r8) :: frNdemand_npool_to_deadstemn
    real(r8) :: frNdemand_npool_to_deadstemn_storage
    real(r8) :: frNdemand_npool_to_livecrootn
    real(r8) :: frNdemand_npool_to_livecrootn_storage
    real(r8) :: frNdemand_npool_to_deadcrootn
    real(r8) :: frNdemand_npool_to_deadcrootn_storage
    real(r8) :: frNdemand_npool_to_reproductiven(nrepr)
    real(r8) :: frNdemand_npool_to_reproductiven_storage(nrepr)

    character(len=*), parameter :: subname = 'calc_npool_to_components_flexiblecn'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(f5) == [nrepr]), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(npool_to_reproductiven) == [nrepr]), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(npool_to_reproductiven_storage) == [nrepr]), sourcefile, __LINE__)

    associate( &
         woody    => pftcon%woody    , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)
         leafcn   => pftcon%leafcn   , & ! Input:  leaf C:N (gC/gN)
         frootcn  => pftcon%frootcn  , & ! Input:  fine root C:N (gC/gN)
         livewdcn => pftcon%livewdcn , & ! Input:  live wood (phloem and ray parenchyma) C:N (gC/gN)
         deadwdcn => pftcon%deadwdcn , & ! Input:  dead wood (xylem and heartwood) C:N (gC/gN)
         graincn  => pftcon%graincn    & ! Input:  grain C:N (gC/gN)
         )

    dt = get_step_size_real()

    cnl  = leafcn(ivt)
    cnfr = frootcn(ivt)
    cnlw = livewdcn(ivt)
    cndw = deadwdcn(ivt)

    ! computing 1.) fractional N demand and 2.) N allocation after uptake for different plant parts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! computing nitrogen demand for different pools based on carbon allocated and CN ratio
    npool_to_leafn_demand          = (nlc / cnl) * fcur
    npool_to_leafn_storage_demand  = (nlc / cnl) * (1._r8 - fcur)
    npool_to_frootn_demand         = (nlc * f1 / cnfr) * fcur
    npool_to_frootn_storage_demand = (nlc * f1 / cnfr) * (1._r8 - fcur)
    if (woody(ivt) == 1._r8) then

       npool_to_livestemn_demand          = (nlc * f3 * f4 / cnlw) * fcur
       npool_to_livestemn_storage_demand  = (nlc * f3 * f4 / cnlw) * (1._r8 - fcur)
       npool_to_deadstemn_demand          = (nlc * f3 * (1._r8 - f4) / cndw) * fcur
       npool_to_deadstemn_storage_demand  = (nlc * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
       npool_to_livecrootn_demand         = (nlc * f2 * f3 * f4 / cnlw) * fcur
       npool_to_livecrootn_storage_demand = (nlc * f2 * f3 * f4 / cnlw) * (1._r8 - fcur)
       npool_to_deadcrootn_demand         = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * fcur
       npool_to_deadcrootn_storage_demand = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
    end if
    if (ivt >= npcropmin) then ! skip 2 generic crops

       cng = graincn(ivt)
       npool_to_livestemn_demand          = (nlc * f3 * f4 / cnlw) * fcur
       npool_to_livestemn_storage_demand  = (nlc * f3 * f4 / cnlw) * (1._r8 - fcur)
       npool_to_deadstemn_demand          = (nlc * f3 * (1._r8 - f4) / cndw) * fcur
       npool_to_deadstemn_storage_demand  = (nlc * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
       npool_to_livecrootn_demand         = (nlc * f2 * f3 * f4 / cnlw) * fcur
       npool_to_livecrootn_storage_demand = (nlc * f2 * f3 * f4 / cnlw) * (1._r8 - fcur)
       npool_to_deadcrootn_demand         = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * fcur
       npool_to_deadcrootn_storage_demand = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
       do k = 1, nrepr
          npool_to_reproductiven_demand(k)         = (nlc * f5(k) / cng) * fcur
          npool_to_reproductiven_storage_demand(k) = (nlc * f5(k) / cng) * (1._r8 -fcur)
       end do
    end if


    ! computing 1.) fractional N demand for different plant parts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    total_plant_Ndemand = npool_to_leafn_demand + npool_to_leafn_storage_demand + &
         npool_to_frootn_demand + npool_to_frootn_storage_demand

    if (woody(ivt) == 1._r8) then

       total_plant_Ndemand = npool_to_leafn_demand + npool_to_leafn_storage_demand + &
            npool_to_frootn_demand + npool_to_frootn_storage_demand + &
            npool_to_livestemn_demand + npool_to_livestemn_storage_demand + npool_to_deadstemn_demand + &
            npool_to_deadstemn_storage_demand  + &
            npool_to_livecrootn_demand + npool_to_livecrootn_storage_demand + npool_to_deadcrootn_demand + &
            npool_to_deadcrootn_storage_demand

    end if
    if (ivt >= npcropmin) then ! skip 2 generic crops

       npool_to_reproductiven_demand_tot = 0._r8
       npool_to_reproductiven_storage_demand_tot = 0._r8
       do k = 1, nrepr
          npool_to_reproductiven_demand_tot = npool_to_reproductiven_demand_tot + &
               npool_to_reproductiven_demand(k)
          npool_to_reproductiven_storage_demand_tot = npool_to_reproductiven_storage_demand_tot + &
               npool_to_reproductiven_storage_demand(k)
       end do

       total_plant_Ndemand = npool_to_leafn_demand + npool_to_leafn_storage_demand + &
            npool_to_frootn_demand + npool_to_frootn_storage_demand + &
            npool_to_livestemn_demand + npool_to_livestemn_storage_demand + npool_to_deadstemn_demand + &
            npool_to_deadstemn_storage_demand  + &
            npool_to_livecrootn_demand + npool_to_livecrootn_storage_demand + npool_to_deadcrootn_demand + &
            npool_to_deadcrootn_storage_demand + &
            npool_to_reproductiven_demand_tot + npool_to_reproductiven_storage_demand_tot

    end if

    if (total_plant_Ndemand == 0.0_r8) then    ! removing division by zero

       frNdemand_npool_to_leafn = 0.0_r8
       frNdemand_npool_to_leafn_storage = 0.0_r8
       frNdemand_npool_to_frootn = 0.0_r8
       frNdemand_npool_to_frootn_storage = 0.0_r8
       if (woody(ivt) == 1._r8) then

          frNdemand_npool_to_livestemn = 0.0_r8
          frNdemand_npool_to_livestemn_storage  = 0.0_r8
          frNdemand_npool_to_deadstemn = 0.0_r8
          frNdemand_npool_to_deadstemn_storage = 0.0_r8
          frNdemand_npool_to_livecrootn = 0.0_r8
          frNdemand_npool_to_livecrootn_storage = 0.0_r8
          frNdemand_npool_to_deadcrootn = 0.0_r8
          frNdemand_npool_to_deadcrootn_storage = 0.0_r8
       end if
       if (ivt >= npcropmin) then ! skip 2 generic crops

          frNdemand_npool_to_livestemn = 0.0_r8
          frNdemand_npool_to_livestemn_storage = 0.0_r8
          frNdemand_npool_to_deadstemn = 0.0_r8
          frNdemand_npool_to_deadstemn_storage = 0.0_r8
          frNdemand_npool_to_livecrootn = 0.0_r8
          frNdemand_npool_to_livecrootn_storage = 0.0_r8
          frNdemand_npool_to_deadcrootn = 0.0_r8
          frNdemand_npool_to_deadcrootn_storage = 0.0_r8
          do k = 1, nrepr
             frNdemand_npool_to_reproductiven(k) = 0.0_r8
             frNdemand_npool_to_reproductiven_storage(k) = 0.0_r8
          end do
       end if

    else

       frNdemand_npool_to_leafn = npool_to_leafn_demand / total_plant_Ndemand
       frNdemand_npool_to_leafn_storage = npool_to_leafn_storage_demand / total_plant_Ndemand
       frNdemand_npool_to_frootn = npool_to_frootn_demand / total_plant_Ndemand
       frNdemand_npool_to_frootn_storage = npool_to_frootn_storage_demand / total_plant_Ndemand
       if (woody(ivt) == 1._r8) then

          frNdemand_npool_to_livestemn = npool_to_livestemn_demand / total_plant_Ndemand
          frNdemand_npool_to_livestemn_storage  = npool_to_livestemn_storage_demand / total_plant_Ndemand
          frNdemand_npool_to_deadstemn = npool_to_deadstemn_demand / total_plant_Ndemand
          frNdemand_npool_to_deadstemn_storage = npool_to_deadstemn_storage_demand / total_plant_Ndemand
          frNdemand_npool_to_livecrootn = npool_to_livecrootn_demand / total_plant_Ndemand
          frNdemand_npool_to_livecrootn_storage = npool_to_livecrootn_storage_demand / total_plant_Ndemand
          frNdemand_npool_to_deadcrootn = npool_to_deadcrootn_demand / total_plant_Ndemand
          frNdemand_npool_to_deadcrootn_storage = npool_to_deadcrootn_storage_demand / total_plant_Ndemand
       end if
       if (ivt >= npcropmin) then ! skip 2 generic crops

          frNdemand_npool_to_livestemn = npool_to_livestemn_demand / total_plant_Ndemand
          frNdemand_npool_to_livestemn_storage = npool_to_livestemn_storage_demand / total_plant_Ndemand
          frNdemand_npool_to_deadstemn = npool_to_deadstemn_demand / total_plant_Ndemand
          frNdemand_npool_to_deadstemn_storage = npool_to_deadstemn_storage_demand / total_plant_Ndemand
          frNdemand_npool_to_livecrootn = npool_to_livecrootn_demand / total_plant_Ndemand
          frNdemand_npool_to_livecrootn_storage = npool_to_livecrootn_storage_demand / total_plant_Ndemand
          frNdemand_npool_to_deadcrootn = npool_to_deadcrootn_demand / total_plant_Ndemand
          frNdemand_npool_to_deadcrootn_storage = npool_to_deadcrootn_storage_demand / total_plant_Ndemand
          do k = 1, nrepr
             frNdemand_npool_to_reproductiven(k) = &
                  npool_to_reproductiven_demand(k) / total_plant_Ndemand
             frNdemand_npool_to_reproductiven_storage(k) = &
                  npool_to_reproductiven_storage_demand(k) / total_plant_Ndemand
          end do
       end if

    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! computing N allocation for different plant parts
    ! allocating allocation to different plant parts in proportion to the fractional demand
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    npool_to_leafn = frNdemand_npool_to_leafn * npool / dt
    npool_to_leafn_storage  = frNdemand_npool_to_leafn_storage * npool / dt
    npool_to_frootn = frNdemand_npool_to_frootn * npool / dt
    npool_to_frootn_storage = frNdemand_npool_to_frootn_storage * npool / dt
    if (woody(ivt) == 1._r8) then
       npool_to_livestemn = frNdemand_npool_to_livestemn * npool / dt
       npool_to_livestemn_storage = frNdemand_npool_to_livestemn_storage * npool / dt
       npool_to_deadstemn = frNdemand_npool_to_deadstemn * npool / dt
       npool_to_deadstemn_storage = frNdemand_npool_to_deadstemn_storage * npool / dt
       npool_to_livecrootn = frNdemand_npool_to_livecrootn * npool / dt
       npool_to_livecrootn_storage = frNdemand_npool_to_livecrootn_storage * npool / dt
       npool_to_deadcrootn = frNdemand_npool_to_deadcrootn * npool / dt
       npool_to_deadcrootn_storage = frNdemand_npool_to_deadcrootn_storage * npool / dt
    end if
    if (ivt >= npcropmin) then ! skip 2 generic crops
       npool_to_livestemn = frNdemand_npool_to_livestemn * npool / dt
       npool_to_livestemn_storage = frNdemand_npool_to_livestemn_storage * npool / dt
       npool_to_deadstemn = frNdemand_npool_to_deadstemn * npool / dt
       npool_to_deadstemn_storage = frNdemand_npool_to_deadstemn_storage * npool / dt
       npool_to_livecrootn = frNdemand_npool_to_livecrootn * npool / dt
       npool_to_livecrootn_storage = frNdemand_npool_to_livecrootn_storage * npool / dt
       npool_to_deadcrootn = frNdemand_npool_to_deadcrootn * npool / dt
       npool_to_deadcrootn_storage = frNdemand_npool_to_deadcrootn_storage * npool / dt
       do k = 1, nrepr
          npool_to_reproductiven(k) = frNdemand_npool_to_reproductiven(k) * npool / dt
          npool_to_reproductiven_storage(k) = frNdemand_npool_to_reproductiven_storage(k) * npool / dt
       end do
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end associate
  end subroutine calc_npool_to_components_flexiblecn

  !-----------------------------------------------------------------------
  subroutine calc_npool_to_components_agsys( &
       npool, fcur, f4, aleaf_n, astem_n, aroot_n, arepr_n, &
       npool_to_leafn, npool_to_leafn_storage, &
       npool_to_frootn, npool_to_frootn_storage, &
       npool_to_livestemn, npool_to_livestemn_storage, &
       npool_to_deadstemn, npool_to_deadstemn_storage, &
       npool_to_livecrootn, npool_to_livecrootn_storage, &
       npool_to_deadcrootn, npool_to_deadcrootn_storage, &
       npool_to_reproductiven, npool_to_reproductiven_storage)
    !
    ! !DESCRIPTION:
    ! Calculate npool_to_* terms for a single crop patch when using AgSys
    !
    ! Note that this assumes that there is no allocation to coarse roots
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: npool      ! temporary plant N pool (gN/m2)
    real(r8), intent(in) :: fcur       ! fraction of current psn displayed as growth
    real(r8), intent(in) :: f4         ! C allocation parameter - fraction of new wood that is live
    real(r8), intent(in) :: aleaf_n    ! leaf allocation coefficient for N
    real(r8), intent(in) :: astem_n    ! stem allocation coefficient for N
    real(r8), intent(in) :: aroot_n    ! root allocation coefficient for N
    real(r8), intent(in) :: arepr_n(:) ! reproductive allocation coefficient(s) for N

    ! Each of the following output variables is in units of gN/m2/s
    real(r8), intent(out) :: npool_to_leafn
    real(r8), intent(out) :: npool_to_leafn_storage
    real(r8), intent(out) :: npool_to_frootn
    real(r8), intent(out) :: npool_to_frootn_storage
    real(r8), intent(out) :: npool_to_livestemn
    real(r8), intent(out) :: npool_to_livestemn_storage
    real(r8), intent(out) :: npool_to_deadstemn
    real(r8), intent(out) :: npool_to_deadstemn_storage
    real(r8), intent(out) :: npool_to_livecrootn
    real(r8), intent(out) :: npool_to_livecrootn_storage
    real(r8), intent(out) :: npool_to_deadcrootn
    real(r8), intent(out) :: npool_to_deadcrootn_storage
    real(r8), intent(out) :: npool_to_reproductiven(:)
    real(r8), intent(out) :: npool_to_reproductiven_storage(:)

    !
    ! !LOCAL VARIABLES:
    real(r8) :: dt                 ! model time step
    integer  :: k

    character(len=*), parameter :: subname = 'calc_npool_to_components_agsys'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(arepr_n) == [nrepr]), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(npool_to_reproductiven) == [nrepr]), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(npool_to_reproductiven_storage) == [nrepr]), sourcefile, __LINE__)

    dt = get_step_size_real()

    npool_to_leafn         = aleaf_n * fcur           * npool / dt
    npool_to_leafn_storage = aleaf_n * (1._r8 - fcur) * npool / dt

    npool_to_frootn         = aroot_n * fcur           * npool / dt
    npool_to_frootn_storage = aroot_n * (1._r8 - fcur) * npool / dt

    npool_to_livestemn         = astem_n * f4           * fcur           * npool / dt
    npool_to_livestemn_storage = astem_n * f4           * (1._r8 - fcur) * npool / dt
    npool_to_deadstemn         = astem_n * (1._r8 - f4) * fcur           * npool / dt
    npool_to_deadstemn_storage = astem_n * (1._r8 - f4) * (1._r8 - fcur) * npool / dt

    ! Assume no allocation to coarse roots for crops. If there *were* allocation to coarse
    ! roots (via a non-zero croot_stem), we would have bigger issues with consistency
    ! between AgSys's desired allocation and the actual C/N allocation: the way this
    ! allocation is currently formulated, non-zero allocation to coarse roots implies that
    ! things like aleaf and arepr aren't truly the fractional allocation to leaves and
    ! reproductive organs: the actual allocation fractions end up being reduced somewhat
    ! via a normalizing factor that differs from 1.
    !
    ! It's not really necessary to explicitly set these to 0 every time step, but we do
    ! it to try make it obvious that this will need to change if we ever want to have
    ! non-zero allocation to coarse roots for crops.
    npool_to_livecrootn         = 0._r8
    npool_to_livecrootn_storage = 0._r8
    npool_to_deadcrootn         = 0._r8
    npool_to_deadcrootn_storage = 0._r8

    do k = 1, nrepr
       npool_to_reproductiven(k)         = arepr_n(k) * fcur           * npool / dt
       npool_to_reproductiven_storage(k) = arepr_n(k) * (1._r8 - fcur) * npool / dt
    end do

  end subroutine calc_npool_to_components_agsys

! -----------------------------------------------------------------------
  subroutine calc_plant_nutrient_demand(this, bounds,                          &
       num_p, filter_p, call_is_for_pcrop,                                     &
       crop_inst, canopystate_inst,                                            &
       cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst,        &
       cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
       soilbiogeochem_carbonflux_inst, soilbiogeochem_nitrogenstate_inst, &
       energyflux_inst)
    !
    ! !USES:
    use CanopyStateType        , only : canopystate_type
    use CropType               , only : crop_type
    use CNVegStateType         , only : cnveg_state_type
    use CNVegCarbonStateType   , only : cnveg_carbonstate_type
    use CNVegNitrogenStateType , only : cnveg_nitrogenstate_type
    use CNVegCarbonFluxType    , only : cnveg_carbonflux_type
    use CNVegNitrogenFluxType  , only : cnveg_nitrogenflux_type
    use SoilBiogeochemCarbonFluxType, only : soilbiogeochem_carbonflux_type
    use SoilBiogeochemNitrogenStateType, only : soilbiogeochem_nitrogenstate_type
    use EnergyFluxType         , only : energyflux_type     !
    ! !ARGUMENTS:
    class(nutrient_competition_FlexibleCN_type), intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds

    ! This subroutine is meant to be called separately for non-prognostic-crop points and
    ! prognostic-crop points. (The reason for this is so that the call for prognostic-crop
    ! points can be skipped when a separate crop model is calculating these variables.) In
    ! the call for non-prognostic-crop points, this filter should be the soilnopcropp
    ! filter and call_is_for_pcrop should be false; in the call for prognostic-crop
    ! points, this filter should be the pcropp filter and call_is_for_pcrop should be
    ! true.
    integer                         , intent(in)    :: num_p        ! number of patches in filter
    integer                         , intent(in)    :: filter_p(:)  ! patch filter
    logical                         , intent(in)    :: call_is_for_pcrop

    type(crop_type)                 , intent(in)    :: crop_inst
    type(canopystate_type)          , intent(in)    :: canopystate_inst
    type(cnveg_state_type)          , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)    , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)     , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_nitrogenstate_type)  , intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_carbonflux_type)   , intent(in) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_nitrogenstate_type), intent(in)    :: soilbiogeochem_nitrogenstate_inst
    type(energyflux_type)           , intent(in)    :: energyflux_inst
    !-----------------------------------------------------------------------

    call this%calc_plant_nitrogen_demand(bounds,                           &
       num_p, filter_p, call_is_for_pcrop,                                 &
       crop_inst, canopystate_inst,                                        &
       cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst,    &
       cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst,                  &
       soilbiogeochem_carbonflux_inst, soilbiogeochem_nitrogenstate_inst, &
       energyflux_inst)

  end subroutine calc_plant_nutrient_demand

  !-----------------------------------------------------------------------
  subroutine calc_plant_nitrogen_demand(this, bounds,                           &
       num_p, filter_p, call_is_for_pcrop,                                      &
       crop_inst, canopystate_inst,                                             &
       cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst,         &
       cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
       soilbiogeochem_carbonflux_inst, soilbiogeochem_nitrogenstate_inst, &
       energyflux_inst)
    !
    ! !DESCRIPTION:
    ! Sets the following output variables that are used elsewhere:
    ! - plant_ndemand
    ! - retransn_to_npool
    ! - leafn_to_retransn
    ! - frootn_to_retransn
    ! - livestemn_to_retransn
    !
    ! !USES:
    use pftconMod              , only : ntmp_soybean, nirrig_tmp_soybean
    use pftconMod              , only : ntrp_soybean, nirrig_trp_soybean
    use clm_varcon             , only : dzsoi_decomp
    use clm_varpar             , only : nlevdecomp
    use CanopyStateType        , only : canopystate_type
    use CropType               , only : crop_type
    use CNVegStateType         , only : cnveg_state_type
    use CNVegCarbonStateType   , only : cnveg_carbonstate_type
    use CNVegCarbonFluxType    , only : cnveg_carbonflux_type
    use CNVegNitrogenStateType , only : cnveg_nitrogenstate_type
    use CNVegNitrogenFluxType  , only : cnveg_nitrogenflux_type
    use SoilBiogeochemCarbonFluxType, only : soilbiogeochem_carbonflux_type
    use SoilBiogeochemNitrogenStateType, only : soilbiogeochem_nitrogenstate_type
    use EnergyFluxType         , only : energyflux_type     !
    use CNSharedParamsMod      , only : use_fun
    use CNPrecisionControlMod  , only : n_min
    use clm_varcon             , only : spval
    use clm_varpar             , only : ileaf,ileaf_st,ileaf_xf,ifroot,ifroot_st,ifroot_xf,&
                                       ilivestem,ilivestem_st,ilivestem_xf,&
                                       ideadstem,ideadstem_st,ideadstem_xf,&
                                       ilivecroot,ilivecroot_st,ilivecroot_xf,&
                                       ideadcroot,ideadcroot_st,ideadcroot_xf,&
                                       igrain,igrain_st,igrain_xf,iretransn,ioutc,ioutn
    use CNVegMatrixMod         , only : matrix_update_phn
    ! !ARGUMENTS:
    class(nutrient_competition_FlexibleCN_type), intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds

    ! This subroutine is meant to be called separately for non-prognostic-crop points and
    ! prognostic-crop points. (The reason for this is so that the call for prognostic-crop
    ! points can be skipped when a separate crop model is calculating these variables.) In
    ! the call for non-prognostic-crop points, this filter should be the soilnopcropp
    ! filter and call_is_for_pcrop should be false; in the call for prognostic-crop
    ! points, this filter should be the pcropp filter and call_is_for_pcrop should be
    ! true.
    integer                         , intent(in)    :: num_p        ! number of patches in filter
    integer                         , intent(in)    :: filter_p(:)  ! patch filter
    logical                         , intent(in)    :: call_is_for_pcrop

    type(crop_type)                 , intent(in)    :: crop_inst
    type(canopystate_type)          , intent(in)    :: canopystate_inst
    type(cnveg_state_type)          , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)    , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)     , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_nitrogenstate_type)  , intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_carbonflux_type)   , intent(in) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_nitrogenstate_type), intent(in) :: soilbiogeochem_nitrogenstate_inst
    type(energyflux_type)           , intent(in)    :: energyflux_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c, p, j                                    ! indices
    integer  :: fp                                         ! lake filter patch index
    real(r8) :: t1                                         ! temporary variable
    real(r8) :: dt                                         ! model time step
    real(r8) :: f_N              (bounds%begp:bounds%endp)
    real(r8) :: Kmin
    real(r8) :: leafcn_max
    real(r8) :: leafcn_min
    real(r8) :: nscalar
    real(r8) :: sminn_total
    real(r8) :: substrate_term
    real(r8) :: temp_scalar
    real(r8) :: Vmax_N
    real(r8) :: crop_phase       (bounds%begp:bounds%endp)
    real(r8) :: allocation_leaf  (bounds%begp:bounds%endp)
    real(r8) :: allocation_stem  (bounds%begp:bounds%endp)
    real(r8) :: allocation_froot (bounds%begp:bounds%endp)
    real(r8) :: tmp

    character(len=*), parameter :: subname = "calc_plant_nitrogen_demand"
    ! -----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(this%actual_leafcn) >= (/bounds%endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((lbound(this%actual_leafcn) <= (/bounds%begp/)), sourcefile, __LINE__)

    associate(                                                                        &
         ivt                   => patch%itype                                        ,  & ! Input:  [integer  (:) ]  patch vegetation type

         leafcn                => pftcon%leafcn                                    ,  & ! Input:  leaf C:N (gC/gN)
         fleafcn               => pftcon%fleafcn                                    , & ! Input:  leaf c:n during organ fill
         ffrootcn              => pftcon%ffrootcn                                   , & ! Input:  froot c:n during organ fill
         fstemcn               => pftcon%fstemcn                                    , & ! Input:  stem c:n during organ fill
         astemf                => pftcon%astemf                                     , & ! Input:  parameter used below
         season_decid          => pftcon%season_decid                               , & ! Input:  binary flag for seasonal-deciduous leaf habit (0 or 1)
         stress_decid          => pftcon%stress_decid                               , & ! Input:  binary flag for stress-deciduous leaf habit (0 or 1)

         laisun                => canopystate_inst%laisun_patch                     , & ! Input:  [real(r8) (:)   ]  sunlit projected leaf area index
         laisha                => canopystate_inst%laisha_patch                     , & ! Input:  [real(r8) (:)   ]  shaded projected leaf area index

         croplive              => crop_inst%croplive_patch                          , & ! Input:  [logical  (:)   ]  flag, true if planted, not harvested

         astem                 => cnveg_state_inst%astem_patch                      , & ! Input: [real(r8) (:)   ]  stem allocation coefficient
         c_allometry           => cnveg_state_inst%c_allometry_patch                , & ! Input: [real(r8) (:)   ]  C allocation index (DIM)
         n_allometry           => cnveg_state_inst%n_allometry_patch                , & ! Input: [real(r8) (:)   ]  N allocation index (DIM)
         annsum_potential_gpp  => cnveg_state_inst%annsum_potential_gpp_patch       , & ! Input:  [real(r8) (:)   ]  annual sum of potential GPP
         annmax_retransn       => cnveg_state_inst%annmax_retransn_patch            , & ! Input:  [real(r8) (:)   ]  annual max of retranslocated N pool
         grain_flag            => cnveg_state_inst%grain_flag_patch                 , & ! Output: [real(r8) (:)   ]  1: grain fill stage; 0: not
         tempsum_potential_gpp => cnveg_state_inst%tempsum_potential_gpp_patch      , & ! Output: [real(r8) (:)   ]  temporary annual sum of potential GPP
         tempmax_retransn      => cnveg_state_inst%tempmax_retransn_patch           , & ! Output: [real(r8) (:)   ]  temporary annual max of retranslocated N pool (gN/m2)

         leafc                 => cnveg_carbonstate_inst%leafc_patch                , & ! Input:  [real(r8) (:)   ]
         frootc                => cnveg_carbonstate_inst%frootc_patch               , & ! Input:  [real(r8) (:)   ]
         livestemc             => cnveg_carbonstate_inst%livestemc_patch            , & ! Input:  [real(r8) (:)   ]
         livecrootc            => cnveg_carbonstate_inst%livecrootc_patch           , & ! Input:  [real(r8) (:)   ]
         retransn              => cnveg_nitrogenstate_inst%retransn_patch           , & ! Input:  [real(r8) (:)   ]  (gN/m2) plant pool of retranslocated N

         gpp                   => cnveg_carbonflux_inst%gpp_before_downreg_patch    , & ! Input:  [real(r8) (:)   ]  GPP flux before downregulation (gC/m2/s)
         availc                => cnveg_carbonflux_inst%availc_patch                , & ! Input:  [real(r8) (:)   ]  C flux available for allocation (gC/m2/s)

         leafn                 => cnveg_nitrogenstate_inst%leafn_patch              , & ! Input:  [real(r8) (:)   ]  (gN/m2) leaf N
         plant_ndemand         => cnveg_nitrogenflux_inst%plant_ndemand_patch       , & ! Output: [real(r8) (:)   ]  N flux required to support initial GPP (gN/m2/s)
         avail_retransn        => cnveg_nitrogenflux_inst%avail_retransn_patch      , & ! Output: [real(r8) (:)   ]  N flux available from retranslocation pool (gN/m2/s)
         retransn_to_npool     => cnveg_nitrogenflux_inst%retransn_to_npool_patch   , & ! Output: [real(r8) (:)   ]  deployment of retranslocated N (gN/m2/s)
         leafn_to_retransn     => cnveg_nitrogenflux_inst%leafn_to_retransn_patch   , & ! Output: [real(r8) (:)   ]
         frootn_to_retransn    => cnveg_nitrogenflux_inst%frootn_to_retransn_patch  , & ! Output: [real(r8) (:)   ]
         livestemn_to_retransn => cnveg_nitrogenflux_inst%livestemn_to_retransn_patch,& ! Output: [real(r8) (:)   ]
         livestemn             => cnveg_nitrogenstate_inst%livestemn_patch          , & ! Input:  [real(r8) (:)   ]  (gN/m2) livestem N
         frootn                => cnveg_nitrogenstate_inst%frootn_patch             , & ! Input:  [real(r8) (:)   ]  (gN/m2) fine root N
         sminn_vr              => soilbiogeochem_nitrogenstate_inst%sminn_vr_col    , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral N
         t_scalar              => soilbiogeochem_carbonflux_inst%t_scalar_col       , & ! Input:  [real(r8) (:,:) ]  soil temperature scalar for decomp
         ileaf_to_iretransn_phn     => cnveg_nitrogenflux_inst%ileaf_to_iretransn_ph, &
         ifroot_to_iretransn_phn    => cnveg_nitrogenflux_inst%ifroot_to_iretransn_ph, &
         ilivestem_to_iretransn_phn => cnveg_nitrogenflux_inst%ilivestem_to_iretransn_ph &
         )

      ! set time steps
      dt = get_step_size_real()

      ! loop over patches to assess the total plant N demand
      do fp = 1, num_p
         p = filter_p(fp)

         ! when we have "if (leafn(p) == 0.0_r8)" below then we
         ! have floating overflow (out of floating point range)
         ! error in "actual_leafcn(p) = leafc(p) / leafn(p)"
         if (leafn(p) < n_min ) then
            ! to avoid division by zero, and to set leafcn to missing value for history files
            this%actual_leafcn(p) = spval
         else
            ! leaf CN ratio
            this%actual_leafcn(p) = leafc(p)  / leafn(p)
         end if

         leafcn_min = leafcn(ivt(p)) - 10.0_r8
         leafcn_max = leafcn(ivt(p)) + 10.0_r8

         this%actual_leafcn(p) = max( this%actual_leafcn(p), leafcn_min-0.0001_r8 )
         this%actual_leafcn(p) = min( this%actual_leafcn(p), leafcn_max )

         nscalar = (this%actual_leafcn(p) - leafcn_min ) / (leafcn_max - leafcn_min)  ! Nitrogen scaler factor
         nscalar = min( max(0.0_r8, nscalar), 1.0_r8 )

         c = patch%column(p)
         sminn_total = 0.0_r8
         do j = 1, nlevdecomp
            sminn_total = sminn_total + sminn_vr(c,j) * dzsoi_decomp(j)
         end do
         Kmin = 1.0_r8
         substrate_term = sminn_total / (sminn_total + Kmin)

         c = patch%column(p)
         temp_scalar=t_scalar(c,1)
         temp_scalar = min( max(0.0_r8, temp_scalar), 1.0_r8 )

         if(use_fun)then ! in FUN, plant_ndemand is just used as a maximum draw on soil N pools.
             plant_ndemand(p) = availc(p)*(n_allometry(p)/c_allometry(p))
         else !FUN
            if (laisun(p)+laisha(p) > 0.0_r8) then
               Vmax_N = 2.7E-8_r8
               plant_ndemand(p) =  Vmax_N  * frootc(p) * substrate_term * temp_scalar * nscalar
            else
               plant_ndemand(p) = 0.0_r8
            end if

            if (this%actual_leafcn(p) < leafcn_min )then
               plant_ndemand(p) = 0.0_r8
            end if
         end if  !FUN

         !if (leafn(p) < n_min ) then
            !! to set leafcn to missing value for history files
            !this%actual_leafcn(p) = spval
         !end if

         ! retranslocated N deployment depends on seasonal cycle of potential GPP
         ! (requires one year run to accumulate demand)

         tempsum_potential_gpp(p) = tempsum_potential_gpp(p) + gpp(p)

         ! Adding the following line to carry max retransn info to CN Annual Update
         tempmax_retransn(p) = max(tempmax_retransn(p),retransn(p))
      end do

      if (call_is_for_pcrop) then
         call CropPhase(bounds, num_p, filter_p, crop_inst, cnveg_state_inst, &
              crop_phase = crop_phase(bounds%begp:bounds%endp))

         do fp = 1, num_p
            p = filter_p(fp)

            if (croplive(p)) then
               if (crop_phase(p) == cphase_leafemerge) then
                  grain_flag(p) = 0._r8 ! setting to 0 while in phase 2
               else if (crop_phase(p) == cphase_grainfill) then
                  !Beth's retranslocation of leafn, stemn, rootn to organ
                  !Filter excess plant N to retransn pool for organ N
                  !Only do one time then hold grain_flag till onset next season

                  ! slevis: Will astem ever = astemf exactly?
                  ! Beth's response: ...looks like astem can equal astemf under the right circumstances.
                  !It might be worth a rewrite to capture what I was trying to do, but the retranslocation for
                  !corn and wheat begins at the beginning of the grain fill stage, but for soybean I was holding it
                  !until after the leaf and stem decline were complete. Looking at how astem is calculated, once the
                  !stem decline is near complete, astem should (usually) be set to astemf. The reason for holding off
                  !on soybean is that the retranslocation scheme begins at the beginning of the grain phase, when the
                  !leaf and stem are still growing, but declining. Since carbon is still getting allocated and now
                  !there is more nitrogen available, the nitrogen can be diverted from grain. For corn and wheat
                  !the impact was probably enough to boost productivity, but for soybean the nitrogen was better off
                  !fulfilling the grain fill. It seems that if the peak lai is reached for soybean though that this
                  !would be bypassed altogether, not the intended outcome. I checked several of my output files and
                  !they all seemed to be going through the retranslocation loop for soybean - good news.

                  if (astem(p) == astemf(ivt(p)) .or. &
                       (ivt(p) /= ntmp_soybean .and. ivt(p) /= nirrig_tmp_soybean .and.&
                       ivt(p) /= ntrp_soybean .and. ivt(p) /= nirrig_trp_soybean)) then
                     if (grain_flag(p) == 0._r8) then
                        t1 = 1 / dt
                        leafn_to_retransn(p) = t1 * max(leafn(p)- (leafc(p) / fleafcn(ivt(p))),0._r8)
                        livestemn_to_retransn(p) = t1 * max(livestemn(p) - (livestemc(p) / fstemcn(ivt(p))),0._r8)
                        frootn_to_retransn(p) = 0._r8
                        if (ffrootcn(ivt(p)) > 0._r8) then
                           frootn_to_retransn(p) = t1 * max(frootn(p) - (frootc(p) / ffrootcn(ivt(p))),0._r8)
                        end if
                        grain_flag(p) = 1._r8
                        if(use_matrixcn)then
                           if(leafn(p) .ne. 0._r8)then
                               leafn_to_retransn(p) = leafn(p) * matrix_update_phn(p,ileaf_to_iretransn_phn,leafn_to_retransn(p) / leafn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                           end if
                           if(frootn(p) .ne. 0._r8)then
                               frootn_to_retransn(p) = frootn(p) * matrix_update_phn(p,ifroot_to_iretransn_phn,frootn_to_retransn(p) / frootn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                           end if
                           if(livestemn(p) .ne. 0._r8)then
                               livestemn_to_retransn(p) = livestemn(p) * matrix_update_phn(p,ilivestem_to_iretransn_phn,livestemn_to_retransn(p) / livestemn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                           end if
                        end if

                     end if
                  end if
               end if
            end if
         end do
      end if

      ! Beth's code: crops pull from retransn pool only during grain fill;
      !              retransn pool has N from leaves, stems, and roots for
      !              retranslocation
      if (call_is_for_pcrop) then
         do fp = 1, num_p
            p = filter_p(fp)

            if (grain_flag(p) == 1._r8) then
               avail_retransn(p) = plant_ndemand(p)
            else
               avail_retransn(p) = 0.0_r8
            end if
         end do
      else
         do fp = 1, num_p
            p = filter_p(fp)

            if (annsum_potential_gpp(p) > 0._r8) then
               avail_retransn(p) = (annmax_retransn(p)/2._r8)*(gpp(p)/annsum_potential_gpp(p))/dt
            else
               avail_retransn(p) = 0.0_r8
            end if
         end do
      end if

      do fp = 1, num_p
         p = filter_p(fp)

         ! make sure available retrans N doesn't exceed storage
         avail_retransn(p) = min(avail_retransn(p), retransn(p)/dt)

         ! modify plant N demand according to the availability of
         ! retranslocated N
         ! take from retransn pool at most the flux required to meet
         ! plant ndemand

         if (plant_ndemand(p) > avail_retransn(p)) then
            retransn_to_npool(p) = avail_retransn(p)
         else
            retransn_to_npool(p) = plant_ndemand(p)
         end if

         if ( .not. use_fun ) then
            plant_ndemand(p) = plant_ndemand(p) - retransn_to_npool(p)
         end if

      end do

    end associate

  end subroutine calc_plant_nitrogen_demand

end module NutrientCompetitionFlexibleCNMod
