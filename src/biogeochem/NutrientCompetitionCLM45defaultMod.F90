module NutrientCompetitionCLM45defaultMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! DESCRIPTION
  ! module contains different subroutines to do soil nutrient competition dynamics
  !
  ! created by Jinyun Tang, Sep 8, 2014
  ! modified by Mariana Vertenstein, Nov 15, 2014
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use decompMod           , only : bounds_type
  use LandunitType        , only : lun
  use ColumnType          , only : col
  use PatchType           , only : patch
  use NutrientCompetitionMethodMod, only : nutrient_competition_method_type
  use CropReprPoolsMod    , only : nrepr
  use CNPhenologyMod      , only : CropPhase
  use CropType            , only : cphase_leafemerge, cphase_grainfill
  use clm_varctl          , only : iulog
  use abortutils          , only : endrun
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  public :: nutrient_competition_clm45default_type
  !
  type, extends(nutrient_competition_method_type) :: nutrient_competition_clm45default_type
     private
   contains
     ! public methocs
     procedure, public :: init                                ! Initialize the class
     procedure, public :: calc_plant_nutrient_competition     ! calculate nutrient yield rate from competition
     procedure, public :: calc_plant_nutrient_demand          ! calculate plant nutrient demand
     !
     ! private methods
     procedure, private:: calc_plant_cn_alloc
     procedure, private:: calc_plant_nitrogen_demand
  end type nutrient_competition_clm45default_type
  !
  interface nutrient_competition_clm45default_type
     ! initialize a new nutrient_competition_clm45default_type object
     module procedure constructor
  end interface nutrient_competition_clm45default_type
  !

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  type(nutrient_competition_clm45default_type) function constructor()
    !
    ! !DESCRIPTION:
    ! Creates an object of type nutrient_competition_clm45default_type.
    ! For now, this is simply a place-holder.

  end function constructor

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize the class (currently empty for this version)
    !
    class(nutrient_competition_clm45default_type) :: this
    type(bounds_type), intent(in) :: bounds

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine calc_plant_nutrient_competition (this,                   &
          bounds, num_soilp, filter_soilp,                            &
          cnveg_state_inst, crop_inst, canopystate_inst, cnveg_carbonstate_inst, &
          cnveg_carbonflux_inst,                                      &
          c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst,       &
          cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst,          &
          soilbiogeochem_nitrogenstate_inst, fpg_col)
    !
    ! !USES:
    use CNVegStateType        , only : cnveg_state_type
    use CropType              , only : crop_type
    use CanopyStateType        , only : canopystate_type
    use CNVegCarbonStateType  , only : cnveg_carbonstate_type
    use CNVegCarbonFluxType   , only : cnveg_carbonflux_type
    use CNVegNitrogenStateType, only : cnveg_nitrogenstate_type
    use CNVegNitrogenFluxType , only : cnveg_nitrogenflux_type
    use SoilBiogeochemNitrogenStateType, only : soilbiogeochem_nitrogenstate_type
    !
    ! !ARGUMENTS:
    class(nutrient_competition_clm45default_type), intent(inout) :: this
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
    type(cnveg_nitrogenstate_type)  , intent(inout) :: cnveg_nitrogenstate_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type), intent(in)    :: soilbiogeochem_nitrogenstate_inst
    real(r8)                        , intent(in)    :: fpg_col(bounds%begc:)

    call this%calc_plant_cn_alloc (bounds, num_soilp, filter_soilp,        &
         cnveg_state_inst, crop_inst, canopystate_inst, &
         cnveg_carbonstate_inst, cnveg_carbonflux_inst, c13_cnveg_carbonflux_inst, &
         c14_cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
         fpg_col=fpg_col(bounds%begc:bounds%endc))

  end subroutine calc_plant_nutrient_competition

  !-----------------------------------------------------------------------
  subroutine calc_plant_cn_alloc (this, bounds, num_soilp, filter_soilp,   &
       cnveg_state_inst, crop_inst, canopystate_inst, &
       cnveg_carbonstate_inst, cnveg_carbonflux_inst, c13_cnveg_carbonflux_inst, &
       c14_cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, fpg_col)
    !
    ! !USES:
    use pftconMod             , only : pftcon, npcropmin
    use clm_varctl            , only : use_c13, use_c14
    use CNVegStateType        , only : cnveg_state_type
    use CropType              , only : crop_type
    use CanopyStateType        , only : canopystate_type
    use CNVegCarbonStateType   , only : cnveg_carbonstate_type
    use CNVegCarbonFluxType   , only : cnveg_carbonflux_type
    use CNVegNitrogenFluxType , only : cnveg_nitrogenflux_type
    use CNSharedParamsMod     , only : use_fun
    use shr_infnan_mod        , only : shr_infnan_isnan

    !
    ! !ARGUMENTS:
    class(nutrient_competition_clm45default_type), intent(in) :: this
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
    real(r8)                        , intent(in)    :: fpg_col(bounds%begc:)
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,l,j,k          ! indices
    integer :: fp                 ! lake filter patch index
    real(r8):: f1,f2,f3,f4,g1,g2  ! allocation parameters
    real(r8):: cnl,cnfr,cnlw,cndw ! C:N ratios for leaf, fine root, and wood
    real(r8):: fcur               ! fraction of current psn displayed as growth
    real(r8):: gresp_storage      ! temporary variable for growth resp to storage
    real(r8):: nlc                ! temporary variable for total new leaf carbon allocation
    real(r8):: f5(nrepr)          ! reproductive allocation parameters
    real(r8):: cng                ! C:N ratio for grain (= cnlw for now; slevis)
    real(r8):: fsmn(bounds%begp:bounds%endp)  ! A emperate variable for adjusting FUN uptakes
   !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(fpg_col) == (/bounds%endc/)), sourcefile, __LINE__)

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
         deadwdcn                     => pftcon%deadwdcn                                           , & ! Input:  dead wood (xylem and heartwood) C:N (gC/gN)
         fcur2                        => pftcon%fcur                                               , & ! Input:  allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
         graincn                      => pftcon%graincn                                            , & ! Input:  grain C:N (gC/gN)
         grperc                       => pftcon%grperc                                             , & ! Input:  growth respiration parameter
         grpnow                       => pftcon%grpnow                                             , & ! Input:  growth respiration parameter

         croplive                     => crop_inst%croplive_patch                                  , & ! Input:  [logical  (:)   ]  flag, true if planted, not harvested

         peaklai                      => cnveg_state_inst%peaklai_patch                            , & ! Input:  [integer  (:)   ]  1: max allowed lai; 0: not at max
         aleaf                        => cnveg_state_inst%aleaf_patch                              , & ! Input:  [real(r8) (:)   ]  leaf allocation coefficient
         astem                        => cnveg_state_inst%astem_patch                              , & ! Input:  [real(r8) (:)   ]  stem allocation coefficient
         aroot                        => cnveg_state_inst%aroot_patch                              , & ! Input:  [real(r8) (:)   ]  root allocation coefficient
         arepr                        => cnveg_state_inst%arepr_patch                              , & ! Input:  [real(r8) (:,:) ]  reproductive allocation coefficient(s)
         c_allometry                  => cnveg_state_inst%c_allometry_patch                        , & ! Output: [real(r8) (:)   ]  C allocation index (DIM)
         n_allometry                  => cnveg_state_inst%n_allometry_patch                        , & ! Output: [real(r8) (:)   ]  N allocation index (DIM)
         downreg                      => cnveg_state_inst%downreg_patch                            , & ! Output: [real(r8) (:)   ]  fractional reduction in GPP due to N limitation (DIM)

         annsum_npp                   => cnveg_carbonflux_inst%annsum_npp_patch                    , & ! Input:  [real(r8) (:)   ]  annual sum of NPP, for wood allocation
         gpp                          => cnveg_carbonflux_inst%gpp_before_downreg_patch            , & ! Output: [real(r8) (:)   ]  GPP flux before downregulation (gC/m2/s)
         availc                       => cnveg_carbonflux_inst%availc_patch                        , & ! Output: [real(r8) (:)   ]  C flux available for allocation (gC/m2/s)
         excess_cflux                 => cnveg_carbonflux_inst%excess_cflux_patch                  , & ! Output: [real(r8) (:)   ]  C flux not allocated due to downregulation (gC/m2/s)
         plant_calloc                 => cnveg_carbonflux_inst%plant_calloc_patch                  , & ! Output: [real(r8) (:)   ]  total allocated C flux (gC/m2/s)
         psnsun_to_cpool              => cnveg_carbonflux_inst%psnsun_to_cpool_patch               , & ! Output: [real(r8) (:)   ]
         psnshade_to_cpool            => cnveg_carbonflux_inst%psnshade_to_cpool_patch             , & ! Output: [real(r8) (:)   ]
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

         plant_ndemand                => cnveg_nitrogenflux_inst%plant_ndemand_patch               , & ! Output: [real(r8) (:)   ]  N flux required to support initial GPP (gN/m2/s)
         plant_nalloc                 => cnveg_nitrogenflux_inst%plant_nalloc_patch                , & ! Output: [real(r8) (:)   ]  total allocated N flux (gN/m2/s)
         npool_to_reproductiven              => cnveg_nitrogenflux_inst%npool_to_reproductiven_patch             , & ! Output: [real(r8) (:,:)   ]  allocation to grain N (gN/m2/s)
         npool_to_reproductiven_storage      => cnveg_nitrogenflux_inst%npool_to_reproductiven_storage_patch     , & ! Output: [real(r8) (:,:)   ]  allocation to grain N storage (gN/m2/s)
         retransn_to_npool            => cnveg_nitrogenflux_inst%retransn_to_npool_patch           , & ! Output: [real(r8) (:)   ]  deployment of retranslocated N (gN/m2/s)
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
         Nactive                      => cnveg_nitrogenflux_inst%Nactive_patch                     , & ! Output:  [real(r8) (:) ]  Mycorrhizal N uptake (gN/m2/s)
         Nnonmyc                      => cnveg_nitrogenflux_inst%Nnonmyc_patch                     , & ! Output:  [real(r8) (:) ]  Non-mycorrhizal N uptake (gN/m2/s)
         Nam                          => cnveg_nitrogenflux_inst%Nam_patch                         , & ! Output:  [real(r8) (:) ]  AM uptake (gN/m2/s)
         Necm                         => cnveg_nitrogenflux_inst%Necm_patch                        , & ! Output:  [real(r8) (:) ]  ECM uptake (gN/m2/s)
         sminn_to_plant_fun           => cnveg_nitrogenflux_inst%sminn_to_plant_fun_patch            & ! Output:  [real(r8) (:) ]  Total N uptake of FUN (gN/m2/s)
         )

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
         cnl  = leafcn(ivt(p))
         cnfr = frootcn(ivt(p))
         cnlw = livewdcn(ivt(p))
         cndw = deadwdcn(ivt(p))
         fcur = fcur2(ivt(p))

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
           if (croplive(p).and.(.not.shr_infnan_isnan(aleaf(p)))) then
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

         if(use_fun)then ! if we are using FUN, we get the N available from there.
            sminn_to_npool(p) = sminn_to_plant_fun(p)
         else ! no FUN. :( we get N available from the FPG calculation in soilbiogeochemistry competition.
            sminn_to_npool(p) = plant_ndemand(p) * fpg(c)
         endif

         plant_nalloc(p) = sminn_to_npool(p) + retransn_to_npool(p)
         plant_calloc(p) = plant_nalloc(p) * (c_allometry(p)/n_allometry(p))

         if(.not.use_fun)then  !ORIGINAL CLM(CN) downregulation code.
            excess_cflux(p) = availc(p) - plant_calloc(p)
	    ! reduce gpp fluxes due to N limitation
	    if (gpp(p) > 0.0_r8) then
	       downreg(p) = excess_cflux(p)/gpp(p)

	       psnsun_to_cpool(p)   = psnsun_to_cpool(p)  *(1._r8 - downreg(p))
	       psnshade_to_cpool(p) = psnshade_to_cpool(p)*(1._r8 - downreg(p))

	       if ( use_c13 ) then
	          c13_cnveg_carbonflux_inst%psnsun_to_cpool_patch(p)   = &
	               c13_cnveg_carbonflux_inst%psnsun_to_cpool_patch(p)  *(1._r8 - downreg(p))
	          c13_cnveg_carbonflux_inst%psnshade_to_cpool_patch(p) = &
	               c13_cnveg_carbonflux_inst%psnshade_to_cpool_patch(p)*(1._r8 - downreg(p))
	       endif
	       if ( use_c14 ) then
	          c14_cnveg_carbonflux_inst%psnsun_to_cpool_patch(p)   = &
	               c14_cnveg_carbonflux_inst%psnsun_to_cpool_patch(p)  *(1._r8 - downreg(p))
	          c14_cnveg_carbonflux_inst%psnshade_to_cpool_patch(p) = &
	               c14_cnveg_carbonflux_inst%psnshade_to_cpool_patch(p)*(1._r8 - downreg(p))
	       endif
	    end if

	 end if !use_fun

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
         if (woody(ivt(p)) == 1._r8) then
            cpool_to_livestemc(p)          = nlc * f3 * f4 * fcur
            cpool_to_livestemc_storage(p)  = nlc * f3 * f4 * (1._r8 - fcur)
            cpool_to_deadstemc(p)          = nlc * f3 * (1._r8 - f4) * fcur
            cpool_to_deadstemc_storage(p)  = nlc * f3 * (1._r8 - f4) * (1._r8 - fcur)
            cpool_to_livecrootc(p)         = nlc * f2 * f3 * f4 * fcur
            cpool_to_livecrootc_storage(p) = nlc * f2 * f3 * f4 * (1._r8 - fcur)
            cpool_to_deadcrootc(p)         = nlc * f2 * f3 * (1._r8 - f4) * fcur
            cpool_to_deadcrootc_storage(p) = nlc * f2 * f3 * (1._r8 - f4) * (1._r8 - fcur)
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
         end if

         ! corresponding N fluxes
         npool_to_leafn(p)          = (nlc / cnl) * fcur
         npool_to_leafn_storage(p)  = (nlc / cnl) * (1._r8 - fcur)
         npool_to_frootn(p)         = (nlc * f1 / cnfr) * fcur
         npool_to_frootn_storage(p) = (nlc * f1 / cnfr) * (1._r8 - fcur)
         if (woody(ivt(p)) == 1._r8) then
            npool_to_livestemn(p)          = (nlc * f3 * f4 / cnlw) * fcur
            npool_to_livestemn_storage(p)  = (nlc * f3 * f4 / cnlw) * (1._r8 - fcur)
            npool_to_deadstemn(p)          = (nlc * f3 * (1._r8 - f4) / cndw) * fcur
            npool_to_deadstemn_storage(p)  = (nlc * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
            npool_to_livecrootn(p)         = (nlc * f2 * f3 * f4 / cnlw) * fcur
            npool_to_livecrootn_storage(p) = (nlc * f2 * f3 * f4 / cnlw) * (1._r8 - fcur)
            npool_to_deadcrootn(p)         = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * fcur
            npool_to_deadcrootn_storage(p) = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            cng = graincn(ivt(p))
            npool_to_livestemn(p)          = (nlc * f3 * f4 / cnlw) * fcur
            npool_to_livestemn_storage(p)  = (nlc * f3 * f4 / cnlw) * (1._r8 - fcur)
            npool_to_deadstemn(p)          = (nlc * f3 * (1._r8 - f4) / cndw) * fcur
            npool_to_deadstemn_storage(p)  = (nlc * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
            npool_to_livecrootn(p)         = (nlc * f2 * f3 * f4 / cnlw) * fcur
            npool_to_livecrootn_storage(p) = (nlc * f2 * f3 * f4 / cnlw) * (1._r8 - fcur)
            npool_to_deadcrootn(p)         = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * fcur
            npool_to_deadcrootn_storage(p) = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
            do k = 1, nrepr
               npool_to_reproductiven(p,k)         = (nlc * f5(k) / cng) * fcur
               npool_to_reproductiven_storage(p,k) = (nlc * f5(k) / cng) * (1._r8 -fcur)
            end do
         end if

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
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            gresp_storage = gresp_storage + cpool_to_livestemc_storage(p)
            do k = 1, nrepr
               gresp_storage = gresp_storage + cpool_to_reproductivec_storage(p,k)
            end do
         end if
         cpool_to_gresp_storage(p) = gresp_storage * g1 * (1._r8 - g2)

      end do ! end patch loop

    end associate

  end subroutine calc_plant_cn_alloc

  !-----------------------------------------------------------------------
  subroutine calc_plant_nutrient_demand(this, bounds,                          &
       num_p, filter_p, call_is_for_pcrop,                                     &
       crop_inst, canopystate_inst,                                            &
       cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst,        &
       cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst,                      &
       soilbiogeochem_carbonflux_inst, soilbiogeochem_nitrogenstate_inst,      &
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
    use EnergyFluxType         , only : energyflux_type
    !
    ! !ARGUMENTS:
    class(nutrient_competition_clm45default_type), intent(inout) :: this
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
    type(canopystate_type)          , intent(in)    :: canopystate_inst ! unused in this version
    type(cnveg_state_type)          , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)    , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)     , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_nitrogenstate_type)  , intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_carbonflux_type)   , intent(in) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_nitrogenstate_type), intent(in) :: soilbiogeochem_nitrogenstate_inst
    type(energyflux_type)           , intent(in)    :: energyflux_inst
    !-----------------------------------------------------------------------

    call this%calc_plant_nitrogen_demand(bounds,                           &
       num_p, filter_p, call_is_for_pcrop,                                 &
       crop_inst,                                                          &
       cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst,    &
       cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst)

  end subroutine calc_plant_nutrient_demand

  !-----------------------------------------------------------------------
  subroutine calc_plant_nitrogen_demand(this, bounds,                           &
       num_p, filter_p, call_is_for_pcrop,                                      &
       crop_inst,                                                               &
       cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst,         &
       cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst)
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
    use pftconMod              , only : npcropmin, pftcon
    use pftconMod              , only : ntmp_soybean, nirrig_tmp_soybean
    use pftconMod              , only : ntrp_soybean, nirrig_trp_soybean
    use clm_time_manager       , only : get_step_size_real
    use CropType               , only : crop_type
    use CNVegStateType         , only : cnveg_state_type
    use CNVegCarbonStateType   , only : cnveg_carbonstate_type
    use CNVegNitrogenStateType , only : cnveg_nitrogenstate_type
    use CNVegCarbonFluxType    , only : cnveg_carbonflux_type
    use CNVegNitrogenFluxType  , only : cnveg_nitrogenflux_type
    use CNSharedParamsMod      , only : use_fun
    !
    ! !ARGUMENTS:
    class(nutrient_competition_clm45default_type), intent(in) :: this
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
    type(cnveg_state_type)          , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)    , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)     , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_nitrogenstate_type)  , intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: p,l,j              ! indices
    integer :: fp                 ! lake filter patch index
    real(r8):: t1                 ! temporary variable
    real(r8):: dt                 ! model time step
    real(r8):: crop_phase(bounds%begp:bounds%endp)

    character(len=*), parameter :: subname = "calc_plant_nitrogen_demand"
    !-----------------------------------------------------------------------

    associate(                                                                        &
         ivt                   => patch%itype                                        ,  & ! Input:  [integer  (:) ]  patch vegetation type

         leafcn                => pftcon%leafcn                                    ,  & ! Input:  leaf C:N (gC/gN)
         frootcn               => pftcon%frootcn                                   ,  & ! Input:  fine root C:N (gC/gN)
         livewdcn              => pftcon%livewdcn                                  , & ! Input:  live wood (phloem and ray parenchyma) C:N (gC/gN)
         graincn               => pftcon%graincn                                    , & ! Input:  grain C:N (gC/gN)
         fleafcn               => pftcon%fleafcn                                    , & ! Input:  leaf c:n during organ fill
         ffrootcn              => pftcon%ffrootcn                                   , & ! Input:  froot c:n during organ fill
         fstemcn               => pftcon%fstemcn                                    , & ! Input:  stem c:n during organ fill
         astemf                => pftcon%astemf                                     , & ! Input:  parameter used below
         season_decid          => pftcon%season_decid                               , & ! Input:  binary flag for seasonal-deciduous leaf habit (0 or 1)
         stress_decid          => pftcon%stress_decid                               , & ! Input:  binary flag for stress-deciduous leaf habit (0 or 1)


         croplive              => crop_inst%croplive_patch                          , & ! Input:  [logical  (:)   ]  flag, true if planted, not harvested

         astem                 => cnveg_state_inst%astem_patch                      , & ! Input: [real(r8) (:)   ]  stem allocation coefficient
         c_allometry           => cnveg_state_inst%c_allometry_patch                , & ! Input:  [real(r8) (:)   ]  C allocation index (DIM)
         n_allometry           => cnveg_state_inst%n_allometry_patch                , & ! Input:  [real(r8) (:)   ]  N allocation index (DIM)
         annsum_potential_gpp  => cnveg_state_inst%annsum_potential_gpp_patch       , & ! Input:  [real(r8) (:)   ]  annual sum of potential GPP
         annmax_retransn       => cnveg_state_inst%annmax_retransn_patch            , & ! Input:  [real(r8) (:)   ]  annual max of retranslocated N pool
         grain_flag            => cnveg_state_inst%grain_flag_patch                 , & ! Output: [real(r8) (:)   ]  1: grain fill stage; 0: not
         tempsum_potential_gpp => cnveg_state_inst%tempsum_potential_gpp_patch      , & ! Output: [real(r8) (:)   ]  temporary annual sum of potential GPP
         tempmax_retransn      => cnveg_state_inst%tempmax_retransn_patch           , & ! Output: [real(r8) (:)   ]  temporary annual max of retranslocated N pool (gN/m2)

         leafc                 => cnveg_carbonstate_inst%leafc_patch                , & ! Input:  [real(r8) (:)   ]
         frootc                => cnveg_carbonstate_inst%frootc_patch               , & ! Input:  [real(r8) (:)   ]
         livestemc             => cnveg_carbonstate_inst%livestemc_patch            , & ! Input:  [real(r8) (:)   ]

         retransn              => cnveg_nitrogenstate_inst%retransn_patch           , & ! Input:  [real(r8) (:)   ]  (gN/m2) plant pool of retranslocated N

         gpp                   => cnveg_carbonflux_inst%gpp_before_downreg_patch    , & ! Input:  [real(r8) (:)   ]  GPP flux before downregulation (gC/m2/s)
         availc                => cnveg_carbonflux_inst%availc_patch                , & ! Input:  [real(r8) (:)   ]  C flux available for allocation (gC/m2/s)

         plant_ndemand         => cnveg_nitrogenflux_inst%plant_ndemand_patch       , & ! Output: [real(r8) (:)   ]  N flux required to support initial GPP (gN/m2/s)
         avail_retransn        => cnveg_nitrogenflux_inst%avail_retransn_patch      , & ! Output: [real(r8) (:)   ]  N flux available from retranslocation pool (gN/m2/s)
         retransn_to_npool     => cnveg_nitrogenflux_inst%retransn_to_npool_patch   , & ! Output: [real(r8) (:)   ]  deployment of retranslocated N (gN/m2/s)
         leafn_to_retransn     => cnveg_nitrogenflux_inst%leafn_to_retransn_patch   , & ! Output: [real(r8) (:)   ]
         frootn_to_retransn    => cnveg_nitrogenflux_inst%frootn_to_retransn_patch  , & ! Output: [real(r8) (:)   ]
         livestemn_to_retransn => cnveg_nitrogenflux_inst%livestemn_to_retransn_patch & ! Output: [real(r8) (:)   ]
         )

      ! set time steps
      dt = get_step_size_real()

      ! loop over patches to assess the total plant N demand
      do fp = 1, num_p
         p = filter_p(fp)

         plant_ndemand(p) = availc(p)*(n_allometry(p)/c_allometry(p))

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
                     if (grain_flag(p) == 0._r8)then
                        if(.not.use_fun) then
                           t1 = 1 / dt
                           leafn_to_retransn(p) = t1 * ((leafc(p) / leafcn(ivt(p))) - (leafc(p) / &
                                fleafcn(ivt(p))))
                           livestemn_to_retransn(p) = t1 * ((livestemc(p) / livewdcn(ivt(p))) - (livestemc(p) / &
                                fstemcn(ivt(p))))
                           frootn_to_retransn(p) = 0._r8
                           if (ffrootcn(ivt(p)) > 0._r8) then
                              frootn_to_retransn(p) = t1 * ((frootc(p) / frootcn(ivt(p))) - (frootc(p) / &
                                   ffrootcn(ivt(p))))
                           end if
                        else !leafn retrans flux is handled in phenology
                           frootn_to_retransn(p) = 0._r8
                           livestemn_to_retransn(p)=0.0_r8
                        end if !fun
                        grain_flag(p) = 1._r8

                     end if
                  end if
               end if
            end if
         end do
      end if

      ! Beth's code: crops pull from retransn pool only during grain fill;
      !              retransn pool has N from leaves, stems, and roots for
      !              retranslocation
      if (.not. use_fun) then
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
            else
               if (season_decid(ivt(p)) == 1._r8.or.stress_decid(ivt(p))==1._r8) then
                  plant_ndemand(p) = plant_ndemand(p) - retransn_to_npool(p)
               end if
            end if
         end do

      end if  !use_fun

      end associate

    end subroutine calc_plant_nitrogen_demand

end module NutrientCompetitionCLM45defaultMod
