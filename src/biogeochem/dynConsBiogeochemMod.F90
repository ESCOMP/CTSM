module dynConsBiogeochemMod

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Handle conservation of biogeochemical quantities (C & N) with dynamic land cover.
  !
  ! !USES:
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use shr_log_mod                  , only : errMsg => shr_log_errMsg
  use decompMod                    , only : bounds_type
  use abortutils                   , only : endrun
  use clm_varctl                   , only : iulog, use_c13, use_c14, use_lch4
  use pftconMod                    , only : pftcon
  use CanopyStateType              , only : canopystate_type
  use PhotosynthesisMod            , only : photosyns_type
  use CNVegStateType               , only : cnveg_state_type
  use CNVegCarbonStateType         , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType          , only : cnveg_carbonflux_type
  use CNVegNitrogenStateType       , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType        , only : cnveg_nitrogenflux_type
  use SoilBiogeochemStateType      , only : soilBiogeochem_state_type
  use SoilBiogeochemCarbonFluxType , only : soilBiogeochem_carbonflux_type
  use SoilBiogeochemCarbonStateType, only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemNitrogenStateType, only : soilbiogeochem_nitrogenstate_type
  use ch4Mod                       , only : ch4_type
  use LandunitType                 , only : lun                
  use ColumnType                   , only : col                
  use PatchType                    , only : patch                
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  !
  public :: dyn_cnbal_patch
  public :: dyn_cnbal_col

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dyn_cnbal_patch(bounds, &
       num_soilp_with_inactive, filter_soilp_with_inactive, &
       prior_weights, patch_state_updater, &
       canopystate_inst, photosyns_inst, cnveg_state_inst,                                &
       cnveg_carbonstate_inst, c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst,    &
       cnveg_carbonflux_inst, c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst,       &
       cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, soilbiogeochem_carbonflux_inst, &
       soilbiogeochem_state_inst)
    !
    ! !DESCRIPTION:
    ! Modify patch-level state and flux variables to maintain carbon and nitrogen balance with
    ! dynamic patch-weights.
    !
    ! !USES:
    use shr_const_mod      , only : SHR_CONST_PDB
    use landunit_varcon    , only : istsoil, istcrop
    use clm_varpar         , only : nlevdecomp
    use clm_varcon         , only : c13ratio, c14ratio, c3_r2, c4_r2
    use clm_time_manager   , only : get_step_size
    use dynPriorWeightsMod , only : prior_weights_type
    use dynPatchStateUpdaterMod, only : patch_state_updater_type
    !
    ! !ARGUMENTS:
    type(bounds_type)                    , intent(in)    :: bounds        
    integer                              , intent(in)    :: num_soilp_with_inactive ! number of points in filter
    integer                              , intent(in)    :: filter_soilp_with_inactive(:) ! soil patch filter that includes inactive points
    type(prior_weights_type)             , intent(in)    :: prior_weights ! weights prior to the subgrid weight updates
    type(patch_state_updater_type)       , intent(in)    :: patch_state_updater
    type(canopystate_type)               , intent(inout) :: canopystate_inst
    type(photosyns_type)                 , intent(inout) :: photosyns_inst
    type(cnveg_state_type)               , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)         , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)         , intent(inout) :: c13_cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)         , intent(inout) :: c14_cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)          , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)          , intent(inout) :: c13_cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)          , intent(inout) :: c14_cnveg_carbonflux_inst
    type(cnveg_nitrogenstate_type)       , intent(inout) :: cnveg_nitrogenstate_inst
    type(cnveg_nitrogenflux_type)        , intent(inout) :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_carbonflux_type) , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_state_type)      , intent(in)    :: soilbiogeochem_state_inst
    !
    ! !LOCAL VARIABLES:
    integer                       :: p,c,l,g,j                     ! indices
    integer                       :: begp, endp
    integer                       :: ier                           ! error code
    real(r8)                      :: dt                            ! land model time step (sec)
    real(r8), allocatable         :: dwt(:)                        ! change in patch weight (relative to gridcell)
    real(r8), allocatable         :: dwt_leafc_seed(:)             ! patch-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_leafn_seed(:)             ! patch-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_deadstemc_seed(:)         ! patch-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_deadstemn_seed(:)         ! patch-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_frootc_to_litter(:)       ! patch-level mass loss due to weight shift (expressed per unit COLUMN area)
    real(r8), allocatable         :: dwt_livecrootc_to_litter(:)   ! patch-level mass loss due to weight shift (expressed per unit COLUMN area)
    real(r8), allocatable         :: dwt_deadcrootc_to_litter(:)   ! patch-level mass loss due to weight shift (expressed per unit COLUMN area)
    real(r8), allocatable, target :: dwt_frootn_to_litter(:)       ! patch-level mass loss due to weight shift (expressed per unit COLUMN area)
    real(r8), allocatable, target :: dwt_livecrootn_to_litter(:)   ! patch-level mass loss due to weight shift (expressed per unit COLUMN area)
    real(r8), allocatable, target :: dwt_deadcrootn_to_litter(:)   ! patch-level mass loss due to weight shift (expressed per unit COLUMN area)
    real(r8), allocatable         :: conv_cflux(:)                 ! patch-level mass loss due to weight shift (expressed per unit GRIDCELL area)
    real(r8), allocatable         :: wood_product_cflux(:)         ! patch-level mass loss due to weight shift (expressed per unit GRIDCELL area)
    real(r8), allocatable         :: crop_product_cflux(:)         ! patch-level mass loss due to weight shift (expressed per unit GRIDCELL area)
    real(r8), allocatable, target :: conv_nflux(:)                 ! patch-level mass loss due to weight shift (expressed per unit GRIDCELL area)
    real(r8), allocatable         :: wood_product_nflux(:)         ! patch-level mass loss due to weight shift (expressed per unit GRIDCELL area)
    real(r8), allocatable         :: crop_product_nflux(:)         ! patch-level mass loss due to weight shift (expressed per unit GRIDCELL area)
    character(len=32)             :: subname='dyn_cbal'            ! subroutine name
    !! C13
    real(r8), allocatable         :: dwt_leafc13_seed(:)           ! patch-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_deadstemc13_seed(:)       ! patch-level mass gain due to seeding of new area
    real(r8), allocatable, target :: dwt_frootc13_to_litter(:)     ! patch-level mass loss due to weight shift (expressed per unit COLUMN area)
    real(r8), allocatable, target :: dwt_livecrootc13_to_litter(:) ! patch-level mass loss due to weight shift (expressed per unit COLUMN area)
    real(r8), allocatable, target :: dwt_deadcrootc13_to_litter(:) ! patch-level mass loss due to weight shift (expressed per unit COLUMN area)
    real(r8), allocatable, target :: conv_c13flux(:)               ! patch-level mass loss due to weight shift (expressed per unit GRIDCELL area)
    real(r8), allocatable         :: wood_product_c13flux(:)       ! patch-level mass loss due to weight shift (expressed per unit GRIDCELL area)
    real(r8), allocatable         :: crop_product_c13flux(:)       ! patch-level mass loss due to weight shift (expressed per unit GRIDCELL area)
    !! C14
    real(r8), allocatable         :: dwt_leafc14_seed(:)           ! patch-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_deadstemc14_seed(:)       ! patch-level mass gain due to seeding of new area
    real(r8), allocatable, target :: dwt_frootc14_to_litter(:)     ! patch-level mass loss due to weight shift (expressed per unit COLUMN area)
    real(r8), allocatable, target :: dwt_livecrootc14_to_litter(:) ! patch-level mass loss due to weight shift (expressed per unit COLUMN area)
    real(r8), allocatable, target :: dwt_deadcrootc14_to_litter(:) ! patch-level mass loss due to weight shift (expressed per unit COLUMN area)
    real(r8), allocatable, target :: conv_c14flux(:)               ! patch-level mass loss due to weight shift (expressed per unit GRIDCELL area)
    real(r8), allocatable         :: wood_product_c14flux(:)       ! patch-level mass loss due to weight shift (expressed per unit GRIDCELL area)
    real(r8), allocatable         :: crop_product_c14flux(:)       ! patch-level mass loss due to weight shift (expressed per unit GRIDCELL area)

    logical  :: patch_initiating(bounds%begp:bounds%endp)

    ! amounts to add to growing patches
    real(r8), parameter :: leafc_seed = 1._r8
    real(r8), parameter :: deadstemc_seed = 0.1_r8

    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    ! Allocate patch-level mass loss arrays
    allocate(dwt(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt'
          call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    allocate(dwt_leafc_seed(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafc_seed'
          call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    allocate(dwt_leafn_seed(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafn_seed'
          call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    allocate(dwt_deadstemc_seed(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemc_seed'
          call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    allocate(dwt_deadstemn_seed(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemn_seed'
          call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    allocate(dwt_frootc_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_frootc_to_litter'
          call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    allocate(dwt_livecrootc_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_livecrootc_to_litter'
          call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    allocate(dwt_deadcrootc_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadcrootc_to_litter'
          call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    allocate(dwt_frootn_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_frootn_to_litter'
          call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    allocate(dwt_livecrootn_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_livecrootn_to_litter'
          call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    allocate(dwt_deadcrootn_to_litter(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadcrootn_to_litter'
          call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    allocate(conv_cflux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for conv_cflux'
          call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    allocate(wood_product_cflux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for wood_product_cflux'
          call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    allocate(crop_product_cflux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for crop_product_cflux'
          call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    allocate(conv_nflux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for conv_nflux'
          call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    allocate(wood_product_nflux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for wood_product_nflux'
          call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    allocate(crop_product_nflux(begp:endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for crop_product_nflux'
          call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    if ( use_c13 ) then
       allocate(dwt_leafc13_seed(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafc13_seed'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       allocate(dwt_deadstemc13_seed(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemc13_seed'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       allocate(dwt_frootc13_to_litter(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_frootc13_to_litter'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       allocate(dwt_livecrootc13_to_litter(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_livecrootc13_to_litter'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       allocate(dwt_deadcrootc13_to_litter(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadcrootc13_to_litter'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       allocate(conv_c13flux(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for conv_c13flux'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       allocate(wood_product_c13flux(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for wood_product_c13flux'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       allocate(crop_product_c13flux(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for crop_product_c13flux'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    endif
    if ( use_c14 ) then
       allocate(dwt_leafc14_seed(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafc14_seed'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       allocate(dwt_deadstemc14_seed(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemc14_seed'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       allocate(dwt_frootc14_to_litter(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_frootc14_to_litter'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       allocate(dwt_livecrootc14_to_litter(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_livecrootc14_to_litter'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       allocate(dwt_deadcrootc14_to_litter(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadcrootc14_to_litter'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       allocate(conv_c14flux(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for conv_c14flux'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       allocate(wood_product_c14flux(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for wood_product_c14flux'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       allocate(crop_product_c14flux(begp:endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for crop_product_c14flux'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    endif
    
    ! Get time step
    dt = real( get_step_size(), r8 )

    patch_initiating = patch_state_updater%patch_initiating(bounds)

    do p = begp,endp
       c = patch%column(p)
       ! initialize all the patch-level local flux arrays
       dwt(p) = 0._r8
       dwt_leafc_seed(p) = 0._r8
       dwt_leafn_seed(p) = 0._r8
       dwt_deadstemc_seed(p) = 0._r8
       dwt_deadstemn_seed(p) = 0._r8
       dwt_frootc_to_litter(p) = 0._r8
       dwt_livecrootc_to_litter(p) = 0._r8
       dwt_deadcrootc_to_litter(p) = 0._r8
       dwt_frootn_to_litter(p) = 0._r8
       dwt_livecrootn_to_litter(p) = 0._r8
       dwt_deadcrootn_to_litter(p) = 0._r8
       conv_cflux(p) = 0._r8
       wood_product_cflux(p) = 0._r8
       crop_product_cflux(p) = 0._r8
       conv_nflux(p) = 0._r8
       wood_product_nflux(p) = 0._r8
       crop_product_nflux(p) = 0._r8
       
       if ( use_c13 ) then
          dwt_leafc13_seed(p) = 0._r8
          dwt_deadstemc13_seed(p) = 0._r8
          dwt_frootc13_to_litter(p) = 0._r8
          dwt_livecrootc13_to_litter(p) = 0._r8
          dwt_deadcrootc13_to_litter(p) = 0._r8
          conv_c13flux(p) = 0._r8
          wood_product_c13flux(p) = 0._r8
          crop_product_c13flux(p) = 0._r8
       endif
       
       if ( use_c14 ) then
          dwt_leafc14_seed(p) = 0._r8
          dwt_deadstemc14_seed(p) = 0._r8
          dwt_frootc14_to_litter(p) = 0._r8
          dwt_livecrootc14_to_litter(p) = 0._r8
          dwt_deadcrootc14_to_litter(p) = 0._r8
          conv_c14flux(p) = 0._r8
          wood_product_c14flux(p) = 0._r8
          crop_product_c14flux(p) = 0._r8
       endif
       
       l = patch%landunit(p)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          
          ! calculate the change in weight for the timestep
          dwt(p) = patch%wtgcell(p)-prior_weights%pwtgcell(p)

          ! Identify patches that are initiating on this timestep and set all the
          ! necessary state and flux variables

          ! TODO(wjs, 2016-06-01) It seems like this code should be moved to a new
          ! subroutine that is solely responsible for initializing newly-active patches

          ! NOTE(wjs, 2016-06-01) We could replace the check of patch_initiating with a
          ! check of something like (patch%active(p) .and. .not.
          ! prior_weights%pactive(p)). That would mean that 0-weight active patches would
          ! remain in their spunup state when they began to grow. But I think that's a bad
          ! idea, because it means that the evolution of the system depends on which
          ! 0-weight patches we choose to keep active. It seems better to reinitialize
          ! patches as soon as they grow to > 0 area, as is done here.

          if (patch_initiating(p)) then

             canopystate_inst%laisun_patch(p) = 0._r8
             canopystate_inst%laisha_patch(p) = 0._r8

             cnveg_state_inst%dormant_flag_patch(p)          = 1._r8
             cnveg_state_inst%days_active_patch(p)           = 0._r8
             cnveg_state_inst%onset_flag_patch(p)            = 0._r8
             cnveg_state_inst%onset_counter_patch(p)         = 0._r8
             cnveg_state_inst%onset_gddflag_patch(p)         = 0._r8
             cnveg_state_inst%onset_fdd_patch(p)             = 0._r8
             cnveg_state_inst%onset_gdd_patch(p)             = 0._r8
             cnveg_state_inst%onset_swi_patch(p)             = 0._r8
             cnveg_state_inst%offset_flag_patch(p)           = 0._r8
             cnveg_state_inst%offset_counter_patch(p)        = 0._r8
             cnveg_state_inst%offset_fdd_patch(p)            = 0._r8
             cnveg_state_inst%offset_swi_patch(p)            = 0._r8
             cnveg_state_inst%lgsf_patch(p)                  = 0._r8
             cnveg_state_inst%bglfr_patch(p)                 = 0._r8
             cnveg_state_inst%bgtr_patch(p)                  = 0._r8
             cnveg_state_inst%annavg_t2m_patch(p)            = cnveg_state_inst%annavg_t2m_col(c)
             cnveg_state_inst%tempavg_t2m_patch(p)           = 0._r8
             cnveg_state_inst%c_allometry_patch(p)           = 0._r8
             cnveg_state_inst%n_allometry_patch(p)           = 0._r8
             cnveg_state_inst%tempsum_potential_gpp_patch(p) = 0._r8
             cnveg_state_inst%annsum_potential_gpp_patch(p)  = 0._r8
             cnveg_state_inst%tempmax_retransn_patch(p)      = 0._r8
             cnveg_state_inst%annmax_retransn_patch(p)       = 0._r8
             cnveg_state_inst%downreg_patch(p)               = 0._r8

             cnveg_carbonflux_inst%xsmrpool_recover_patch(p)      = 0._r8
             cnveg_carbonflux_inst%plant_calloc_patch(p)          = 0._r8
             cnveg_carbonflux_inst%excess_cflux_patch(p)          = 0._r8
             cnveg_carbonflux_inst%prev_leafc_to_litter_patch(p)  = 0._r8
             cnveg_carbonflux_inst%prev_frootc_to_litter_patch(p) = 0._r8
             cnveg_carbonflux_inst%availc_patch(p)                = 0._r8
             cnveg_carbonflux_inst%gpp_before_downreg_patch(p)    = 0._r8

             cnveg_carbonflux_inst%tempsum_npp_patch(p)       = 0._r8
             cnveg_carbonflux_inst%annsum_npp_patch(p)        = 0._r8

             cnveg_nitrogenflux_inst%plant_ndemand_patch(p)       = 0._r8
             cnveg_nitrogenflux_inst%avail_retransn_patch(p)      = 0._r8
             cnveg_nitrogenflux_inst%plant_nalloc_patch(p)        = 0._r8

             if ( use_c13 ) then
                c13_cnveg_carbonflux_inst%xsmrpool_c13ratio_patch(p) = c13ratio
             end if

             call photosyns_inst%NewPatchinit(p)

          end if  ! end initialization of new patch
       end if     ! is soil
    end do        ! patch loop

    ! Determine annually-smoothed (dribbled) change in weight
    call CNveg_state_inst%dwt_dribbler_patch%set_curr_delta(bounds, &
         dwt(bounds%begp:bounds%endp))
    call CNveg_state_inst%dwt_dribbler_patch%get_dribbled_delta(bounds, &
         CNveg_state_inst%dwt_smoothed_patch(bounds%begp:bounds%endp))

    ! Adjust patch variables and compute associated fluxes for changing patch areas

    call cnveg_carbonstate_inst%DynamicPatchAdjustments(bounds, &
         num_soilp_with_inactive, filter_soilp_with_inactive, &
         patch_state_updater, &
         leafc_seed = leafc_seed, &
         deadstemc_seed = deadstemc_seed, &
         conv_cflux = conv_cflux(begp:endp), &
         wood_product_cflux = wood_product_cflux(begp:endp), &
         crop_product_cflux = crop_product_cflux(begp:endp), &
         dwt_frootc_to_litter = dwt_frootc_to_litter(begp:endp), &
         dwt_livecrootc_to_litter = dwt_livecrootc_to_litter(begp:endp), &
         dwt_deadcrootc_to_litter = dwt_deadcrootc_to_litter(begp:endp), &
         dwt_leafc_seed = dwt_leafc_seed(begp:endp), &
         dwt_deadstemc_seed = dwt_deadstemc_seed(begp:endp))

    ! These fluxes are computed as negative quantities, but are expected to be positive,
    ! so flip the signs
    do p = begp, endp
       dwt_frootc_to_litter(p) = -1._r8 * dwt_frootc_to_litter(p)
       dwt_livecrootc_to_litter(p) = -1._r8 * dwt_livecrootc_to_litter(p)
       dwt_deadcrootc_to_litter(p) = -1._r8 * dwt_deadcrootc_to_litter(p)
    end do

    if (use_c13) then
       call c13_cnveg_carbonstate_inst%DynamicPatchAdjustments(bounds, &
            num_soilp_with_inactive, filter_soilp_with_inactive, &
            patch_state_updater, &
            leafc_seed = leafc_seed, &
            deadstemc_seed = deadstemc_seed, &
            conv_cflux = conv_c13flux(begp:endp), &
            wood_product_cflux = wood_product_c13flux(begp:endp), &
            crop_product_cflux = crop_product_c13flux(begp:endp), &
            dwt_frootc_to_litter = dwt_frootc13_to_litter(begp:endp), &
            dwt_livecrootc_to_litter = dwt_livecrootc13_to_litter(begp:endp), &
            dwt_deadcrootc_to_litter = dwt_deadcrootc13_to_litter(begp:endp), &
            dwt_leafc_seed = dwt_leafc13_seed(begp:endp), &
            dwt_deadstemc_seed = dwt_deadstemc13_seed(begp:endp))

       ! These fluxes are computed as negative quantities, but are expected to be positive,
       ! so flip the signs
       do p = begp,endp
          dwt_frootc13_to_litter(p) = -1._r8 * dwt_frootc13_to_litter(p)
          dwt_livecrootc13_to_litter(p) = -1._r8 * dwt_livecrootc13_to_litter(p)
          dwt_deadcrootc13_to_litter(p) = -1._r8 * dwt_deadcrootc13_to_litter(p)
       end do

    end if

    if (use_c14) then
       call c14_cnveg_carbonstate_inst%DynamicPatchAdjustments(bounds, &
            num_soilp_with_inactive, filter_soilp_with_inactive, &
            patch_state_updater, &
            leafc_seed = leafc_seed, &
            deadstemc_seed = deadstemc_seed, &
            conv_cflux = conv_c14flux(begp:endp), &
            wood_product_cflux = wood_product_c14flux(begp:endp), &
            crop_product_cflux = crop_product_c14flux(begp:endp), &
            dwt_frootc_to_litter = dwt_frootc14_to_litter(begp:endp), &
            dwt_livecrootc_to_litter = dwt_livecrootc14_to_litter(begp:endp), &
            dwt_deadcrootc_to_litter = dwt_deadcrootc14_to_litter(begp:endp), &
            dwt_leafc_seed = dwt_leafc14_seed(begp:endp), &
            dwt_deadstemc_seed = dwt_deadstemc14_seed(begp:endp))

       ! These fluxes are computed as negative quantities, but are expected to be positive,
       ! so flip the signs
       do p = begp,endp
          dwt_frootc14_to_litter(p) = -1._r8 * dwt_frootc14_to_litter(p)
          dwt_livecrootc14_to_litter(p) = -1._r8 * dwt_livecrootc14_to_litter(p)
          dwt_deadcrootc14_to_litter(p) = -1._r8 * dwt_deadcrootc14_to_litter(p)
       end do

    end if

    call cnveg_nitrogenstate_inst%DynamicPatchAdjustments(bounds, &
         num_soilp_with_inactive, filter_soilp_with_inactive, &
         patch_state_updater, &
         leafc_seed = leafc_seed, &
         deadstemc_seed = deadstemc_seed, &
         conv_nflux = conv_nflux(begp:endp), &
         wood_product_nflux = wood_product_nflux(begp:endp), &
         crop_product_nflux = crop_product_nflux(begp:endp), &
         dwt_frootn_to_litter = dwt_frootn_to_litter(begp:endp), &
         dwt_livecrootn_to_litter = dwt_livecrootn_to_litter(begp:endp), &
         dwt_deadcrootn_to_litter = dwt_deadcrootn_to_litter(begp:endp), &
         dwt_leafn_seed = dwt_leafn_seed(begp:endp), &
         dwt_deadstemn_seed = dwt_deadstemn_seed(begp:endp))

    ! These fluxes are computed as negative quantities, but are expected to be positive,
    ! so flip the signs
    do p = begp,endp
       dwt_frootn_to_litter(p) = -1._r8 * dwt_frootn_to_litter(p)
       dwt_livecrootn_to_litter(p) = -1._r8 * dwt_livecrootn_to_litter(p)
       dwt_deadcrootn_to_litter(p) = -1._r8 * dwt_deadcrootn_to_litter(p)
    end do

    ! calculate column-level seeding fluxes
    do p = bounds%begp, bounds%endp
       g = patch%gridcell(p)

       ! C fluxes
       cnveg_carbonflux_inst%dwt_seedc_to_leaf_patch(p) = dwt_leafc_seed(p)/dt
       cnveg_carbonflux_inst%dwt_seedc_to_leaf_grc(g) = &
            cnveg_carbonflux_inst%dwt_seedc_to_leaf_grc(g) + &
            cnveg_carbonflux_inst%dwt_seedc_to_leaf_patch(p)

       cnveg_carbonflux_inst%dwt_seedc_to_deadstem_patch(p) = dwt_deadstemc_seed(p)/dt
       cnveg_carbonflux_inst%dwt_seedc_to_deadstem_grc(g) = &
            cnveg_carbonflux_inst%dwt_seedc_to_deadstem_grc(g) + &
            cnveg_carbonflux_inst%dwt_seedc_to_deadstem_patch(p)

       if ( use_c13 ) then
          c13_cnveg_carbonflux_inst%dwt_seedc_to_leaf_patch(p) = dwt_leafc13_seed(p)/dt
          c13_cnveg_carbonflux_inst%dwt_seedc_to_leaf_grc(g) = &
               c13_cnveg_carbonflux_inst%dwt_seedc_to_leaf_grc(g) + &
               c13_cnveg_carbonflux_inst%dwt_seedc_to_leaf_patch(p)

          c13_cnveg_carbonflux_inst%dwt_seedc_to_deadstem_patch(p) = dwt_deadstemc13_seed(p)/dt
          c13_cnveg_carbonflux_inst%dwt_seedc_to_deadstem_grc(g) = &
               c13_cnveg_carbonflux_inst%dwt_seedc_to_deadstem_grc(g) + &
               c13_cnveg_carbonflux_inst%dwt_seedc_to_deadstem_patch(p)
       endif

       if ( use_c14 ) then	
          c14_cnveg_carbonflux_inst%dwt_seedc_to_leaf_patch(p) = dwt_leafc14_seed(p)/dt
          c14_cnveg_carbonflux_inst%dwt_seedc_to_leaf_grc(g) = &
               c14_cnveg_carbonflux_inst%dwt_seedc_to_leaf_grc(g) + &
               c14_cnveg_carbonflux_inst%dwt_seedc_to_leaf_patch(p)

          c14_cnveg_carbonflux_inst%dwt_seedc_to_deadstem_patch(p) = dwt_deadstemc14_seed(p)/dt
          c14_cnveg_carbonflux_inst%dwt_seedc_to_deadstem_grc(g) = &
               c14_cnveg_carbonflux_inst%dwt_seedc_to_deadstem_grc(g) + &
               c14_cnveg_carbonflux_inst%dwt_seedc_to_deadstem_patch(p)
       endif

       ! N fluxes
       cnveg_nitrogenflux_inst%dwt_seedn_to_leaf_patch(p) = dwt_leafn_seed(p)/dt
       cnveg_nitrogenflux_inst%dwt_seedn_to_leaf_grc(g) = &
            cnveg_nitrogenflux_inst%dwt_seedn_to_leaf_grc(g) + &
            cnveg_nitrogenflux_inst%dwt_seedn_to_leaf_patch(p)

       cnveg_nitrogenflux_inst%dwt_seedn_to_deadstem_patch(p) = dwt_deadstemn_seed(p)/dt
       cnveg_nitrogenflux_inst%dwt_seedn_to_deadstem_grc(g) = &
            cnveg_nitrogenflux_inst%dwt_seedn_to_deadstem_grc(g) + &
            cnveg_nitrogenflux_inst%dwt_seedn_to_deadstem_patch(p)

    end do
    

    ! calculate patch-to-column slash fluxes into litter and CWD pools
    do p = bounds%begp, bounds%endp
       c = patch%column(p)

       ! fine and coarse root to litter and CWD slash carbon fluxes
       cnveg_carbonflux_inst%dwt_slash_cflux_col(c) = &
               cnveg_carbonflux_inst%dwt_slash_cflux_col(c) + &
               dwt_frootc_to_litter(p)/dt + &
               dwt_livecrootc_to_litter(p)/dt + &
               dwt_deadcrootc_to_litter(p)/dt

       if ( use_c13 ) then
          c13_cnveg_carbonflux_inst%dwt_slash_cflux_col(c) = &
                  c13_cnveg_carbonflux_inst%dwt_slash_cflux_col(c) + &
                  dwt_frootc13_to_litter(p)/dt + &
                  dwt_livecrootc13_to_litter(p)/dt + &
                  dwt_deadcrootc13_to_litter(p)/dt
       endif

       if ( use_c14 ) then
          c14_cnveg_carbonflux_inst%dwt_slash_cflux_col(c) = &
                  c14_cnveg_carbonflux_inst%dwt_slash_cflux_col(c) + &
                  dwt_frootc14_to_litter(p)/dt + &
                  dwt_livecrootc14_to_litter(p)/dt + &
                  dwt_deadcrootc14_to_litter(p)/dt
       endif

    end do

    
    ! calculate patch-to-column for fluxes into litter and CWD pools
    do j = 1, nlevdecomp
       do p = bounds%begp, bounds%endp
          c = patch%column(p)

          ! fine root litter carbon fluxes
          cnveg_carbonflux_inst%dwt_frootc_to_litr_met_c_col(c,j) = &
               cnveg_carbonflux_inst%dwt_frootc_to_litr_met_c_col(c,j) + &
               (dwt_frootc_to_litter(p)*pftcon%fr_flab(patch%itype(p)))/dt &
               * soilbiogeochem_state_inst%froot_prof_patch(p,j)

          cnveg_carbonflux_inst%dwt_frootc_to_litr_cel_c_col(c,j) = &
               cnveg_carbonflux_inst%dwt_frootc_to_litr_cel_c_col(c,j) + &
               (dwt_frootc_to_litter(p)*pftcon%fr_fcel(patch%itype(p)))/dt &
               * soilbiogeochem_state_inst%froot_prof_patch(p,j)

          cnveg_carbonflux_inst%dwt_frootc_to_litr_lig_c_col(c,j) = &
               cnveg_carbonflux_inst%dwt_frootc_to_litr_lig_c_col(c,j) + &
               (dwt_frootc_to_litter(p)*pftcon%fr_flig(patch%itype(p)))/dt &
               * soilbiogeochem_state_inst%froot_prof_patch(p,j)


          ! fine root litter nitrogen fluxes
          cnveg_nitrogenflux_inst%dwt_frootn_to_litr_met_n_col(c,j) = &
               cnveg_nitrogenflux_inst%dwt_frootn_to_litr_met_n_col(c,j) + &
               (dwt_frootn_to_litter(p)*pftcon%fr_flab(patch%itype(p)))/dt &
               * soilbiogeochem_state_inst%froot_prof_patch(p,j)
          cnveg_nitrogenflux_inst%dwt_frootn_to_litr_cel_n_col(c,j) = &
               cnveg_nitrogenflux_inst%dwt_frootn_to_litr_cel_n_col(c,j) + &
               (dwt_frootn_to_litter(p)*pftcon%fr_fcel(patch%itype(p)))/dt &
               * soilbiogeochem_state_inst%froot_prof_patch(p,j)

          cnveg_nitrogenflux_inst%dwt_frootn_to_litr_lig_n_col(c,j) = &
               cnveg_nitrogenflux_inst%dwt_frootn_to_litr_lig_n_col(c,j) + &
               (dwt_frootn_to_litter(p)*pftcon%fr_flig(patch%itype(p)))/dt &
               * soilbiogeochem_state_inst%froot_prof_patch(p,j)

          ! livecroot fluxes to cwd
          cnveg_carbonflux_inst%dwt_livecrootc_to_cwdc_col(c,j) = &
               cnveg_carbonflux_inst%dwt_livecrootc_to_cwdc_col(c,j) + &
               (dwt_livecrootc_to_litter(p))/dt * soilbiogeochem_state_inst%croot_prof_patch(p,j)

          cnveg_nitrogenflux_inst%dwt_livecrootn_to_cwdn_col(c,j) = &
               cnveg_nitrogenflux_inst%dwt_livecrootn_to_cwdn_col(c,j) + &
               (dwt_livecrootn_to_litter(p))/dt * soilbiogeochem_state_inst%croot_prof_patch(p,j)

          ! deadcroot fluxes to cwd
          cnveg_carbonflux_inst%dwt_deadcrootc_to_cwdc_col(c,j) = &
               cnveg_carbonflux_inst%dwt_deadcrootc_to_cwdc_col(c,j) + &
               (dwt_deadcrootc_to_litter(p))/dt * soilbiogeochem_state_inst%croot_prof_patch(p,j)

          cnveg_nitrogenflux_inst%dwt_deadcrootn_to_cwdn_col(c,j) = &
               cnveg_nitrogenflux_inst%dwt_deadcrootn_to_cwdn_col(c,j) + &
               (dwt_deadcrootn_to_litter(p))/dt * soilbiogeochem_state_inst%croot_prof_patch(p,j)

          if ( use_c13 ) then
             ! C13 fine root litter fluxes
             c13_cnveg_carbonflux_inst%dwt_frootc_to_litr_met_c_col(c,j) = &
                  c13_cnveg_carbonflux_inst%dwt_frootc_to_litr_met_c_col(c,j) + &
                  (dwt_frootc13_to_litter(p)*pftcon%fr_flab(patch%itype(p)))/dt &
                  * soilbiogeochem_state_inst%froot_prof_patch(p,j)

             c13_cnveg_carbonflux_inst%dwt_frootc_to_litr_cel_c_col(c,j) = &
                  c13_cnveg_carbonflux_inst%dwt_frootc_to_litr_cel_c_col(c,j) + &
                  (dwt_frootc13_to_litter(p)*pftcon%fr_fcel(patch%itype(p)))/dt &
                  * soilbiogeochem_state_inst%froot_prof_patch(p,j)

             c13_cnveg_carbonflux_inst%dwt_frootc_to_litr_lig_c_col(c,j) = &
                  c13_cnveg_carbonflux_inst%dwt_frootc_to_litr_lig_c_col(c,j) + &
                  (dwt_frootc13_to_litter(p)*pftcon%fr_flig(patch%itype(p)))/dt &
                  * soilbiogeochem_state_inst%froot_prof_patch(p,j)

             ! livecroot fluxes to cwd
             c13_cnveg_carbonflux_inst%dwt_livecrootc_to_cwdc_col(c,j) = &
                  c13_cnveg_carbonflux_inst%dwt_livecrootc_to_cwdc_col(c,j) + &
                  (dwt_livecrootc13_to_litter(p))/dt * soilbiogeochem_state_inst%croot_prof_patch(p,j)

             ! deadcroot fluxes to cwd
             c13_cnveg_carbonflux_inst%dwt_deadcrootc_to_cwdc_col(c,j) = &
                  c13_cnveg_carbonflux_inst%dwt_deadcrootc_to_cwdc_col(c,j) + &
                  (dwt_deadcrootc13_to_litter(p))/dt * soilbiogeochem_state_inst%croot_prof_patch(p,j)

          endif

          if ( use_c14 ) then                   
             ! C14 fine root litter fluxes
             c14_cnveg_carbonflux_inst%dwt_frootc_to_litr_met_c_col(c,j) = &
                  c14_cnveg_carbonflux_inst%dwt_frootc_to_litr_met_c_col(c,j) + &
                  (dwt_frootc14_to_litter(p)*pftcon%fr_flab(patch%itype(p)))/dt &
                  * soilbiogeochem_state_inst%froot_prof_patch(p,j)

             c14_cnveg_carbonflux_inst%dwt_frootc_to_litr_cel_c_col(c,j) = &
                  c14_cnveg_carbonflux_inst%dwt_frootc_to_litr_cel_c_col(c,j) + &
                  (dwt_frootc14_to_litter(p)*pftcon%fr_fcel(patch%itype(p)))/dt &
                  * soilbiogeochem_state_inst%froot_prof_patch(p,j)

             c14_cnveg_carbonflux_inst%dwt_frootc_to_litr_lig_c_col(c,j) = &
                  c14_cnveg_carbonflux_inst%dwt_frootc_to_litr_lig_c_col(c,j) + &
                  (dwt_frootc14_to_litter(p)*pftcon%fr_flig(patch%itype(p)))/dt &
                  * soilbiogeochem_state_inst%froot_prof_patch(p,j)

             ! livecroot fluxes to cwd
             c14_cnveg_carbonflux_inst%dwt_livecrootc_to_cwdc_col(c,j) = &
                  c14_cnveg_carbonflux_inst%dwt_livecrootc_to_cwdc_col(c,j) + &
                  (dwt_livecrootc14_to_litter(p))/dt * soilbiogeochem_state_inst%croot_prof_patch(p,j)

             ! deadcroot fluxes to cwd
             c14_cnveg_carbonflux_inst%dwt_deadcrootc_to_cwdc_col(c,j) = &
                  c14_cnveg_carbonflux_inst%dwt_deadcrootc_to_cwdc_col(c,j) + &
                  (dwt_deadcrootc14_to_litter(p))/dt * soilbiogeochem_state_inst%croot_prof_patch(p,j)
          endif

       end do
    end do

    ! Store fluxes into product pools. Note that the temporary conv_cflux, wood_product_cflux
    ! (and similar) fluxes are accumulated as negative values, but the values stored in
    ! carbonflux_inst and nitrogenflux_inst are positive values.
    do p = begp, endp
       cnveg_carbonflux_inst%dwt_wood_productc_gain_patch(p) = -wood_product_cflux(p)/dt
       cnveg_carbonflux_inst%dwt_crop_productc_gain_patch(p) = -crop_product_cflux(p)/dt
       if (use_c13) then
          c13_cnveg_carbonflux_inst%dwt_wood_productc_gain_patch(p) = -wood_product_c13flux(p)/dt
          c13_cnveg_carbonflux_inst%dwt_crop_productc_gain_patch(p) = -crop_product_c13flux(p)/dt
       end if
       if (use_c14) then
          c14_cnveg_carbonflux_inst%dwt_wood_productc_gain_patch(p) = -wood_product_c14flux(p)/dt
          c14_cnveg_carbonflux_inst%dwt_crop_productc_gain_patch(p) = -crop_product_c14flux(p)/dt
       end if
       cnveg_nitrogenflux_inst%dwt_wood_productn_gain_patch(p) = -wood_product_nflux(p)/dt
       cnveg_nitrogenflux_inst%dwt_crop_productn_gain_patch(p) = -crop_product_nflux(p)/dt
    end do

    ! Set column-level conversion fluxes

    do p = bounds%begp, bounds%endp
       g = patch%gridcell(p)

       ! Note that patch-level fluxes are stored per unit GRIDCELL area - thus, we don't
       ! need to multiply by the patch's gridcell weight when translating patch-level
       ! fluxes into gridcell-level fluxes.
       
       cnveg_carbonflux_inst%dwt_conv_cflux_patch(p) = -conv_cflux(p)/dt
       cnveg_carbonflux_inst%dwt_conv_cflux_grc(g) = &
            cnveg_carbonflux_inst%dwt_conv_cflux_grc(g) + &
            cnveg_carbonflux_inst%dwt_conv_cflux_patch(p)

       if ( use_c13 ) then
          ! C13 column-level flux updates
          c13_cnveg_carbonflux_inst%dwt_conv_cflux_patch(p) = -conv_c13flux(p)/dt
          c13_cnveg_carbonflux_inst%dwt_conv_cflux_grc(g) = &
               c13_cnveg_carbonflux_inst%dwt_conv_cflux_grc(g) + &
               c13_cnveg_carbonflux_inst%dwt_conv_cflux_patch(p)
       endif

       if ( use_c14 ) then
          ! C14 column-level flux updates
          c14_cnveg_carbonflux_inst%dwt_conv_cflux_patch(p) = -conv_c14flux(p)/dt
          c14_cnveg_carbonflux_inst%dwt_conv_cflux_grc(g) = &
               c14_cnveg_carbonflux_inst%dwt_conv_cflux_grc(g) + &
               c14_cnveg_carbonflux_inst%dwt_conv_cflux_patch(p)
       endif

       cnveg_nitrogenflux_inst%dwt_conv_nflux_patch(p) = -conv_nflux(p)/dt
       cnveg_nitrogenflux_inst%dwt_conv_nflux_grc(g) = &
            cnveg_nitrogenflux_inst%dwt_conv_nflux_grc(g) + &
            cnveg_nitrogenflux_inst%dwt_conv_nflux_patch(p)

    end do

    ! Deallocate patch-level flux arrays
    deallocate(dwt)
    deallocate(dwt_leafc_seed)
    deallocate(dwt_leafn_seed)
    deallocate(dwt_deadstemc_seed)
    deallocate(dwt_deadstemn_seed)
    deallocate(dwt_frootc_to_litter)
    deallocate(dwt_livecrootc_to_litter)
    deallocate(dwt_deadcrootc_to_litter)
    deallocate(dwt_frootn_to_litter)
    deallocate(dwt_livecrootn_to_litter)
    deallocate(dwt_deadcrootn_to_litter)
    deallocate(conv_cflux)
    deallocate(wood_product_cflux)
    deallocate(crop_product_cflux)
    deallocate(conv_nflux)
    deallocate(wood_product_nflux)
    deallocate(crop_product_nflux)
             
    if ( use_c13 ) then
       deallocate(dwt_leafc13_seed)
       deallocate(dwt_deadstemc13_seed)
       deallocate(dwt_frootc13_to_litter)
       deallocate(dwt_livecrootc13_to_litter)
       deallocate(dwt_deadcrootc13_to_litter)
       deallocate(conv_c13flux)
       deallocate(wood_product_c13flux)
       deallocate(crop_product_c13flux)
    endif
             
    if ( use_c14 ) then
       deallocate(dwt_leafc14_seed)
       deallocate(dwt_deadstemc14_seed)
       deallocate(dwt_frootc14_to_litter)
       deallocate(dwt_livecrootc14_to_litter)
       deallocate(dwt_deadcrootc14_to_litter)
       deallocate(conv_c14flux)
       deallocate(wood_product_c14flux)
       deallocate(crop_product_c14flux)
    endif
    
   end subroutine dyn_cnbal_patch

   !-----------------------------------------------------------------------
   subroutine dyn_cnbal_col(bounds, clump_index, column_state_updater, &
        soilbiogeochem_carbonstate_inst, c13_soilbiogeochem_carbonstate_inst, &
        c14_soilbiogeochem_carbonstate_inst, soilbiogeochem_nitrogenstate_inst, &
        ch4_inst)
     !
     ! !DESCRIPTION:
     ! Modify column-level state variables to maintain carbon and nitrogen balance with
     ! dynamic column weights.
     !
     ! !USES:
     use dynColumnStateUpdaterMod, only : column_state_updater_type
     !
     ! !ARGUMENTS:
     type(bounds_type)                       , intent(in)    :: bounds        

     ! Index of clump on which we're currently operating. Note that this implies that this
     ! routine must be called from within a clump loop.
     integer                                 , intent(in)    :: clump_index

     type(column_state_updater_type)         , intent(in)    :: column_state_updater
     type(soilbiogeochem_carbonstate_type)   , intent(inout) :: soilbiogeochem_carbonstate_inst
     type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c13_soilbiogeochem_carbonstate_inst
     type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c14_soilbiogeochem_carbonstate_inst
     type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
     type(ch4_type)                          , intent(inout) :: ch4_inst
     !
     ! !LOCAL VARIABLES:

     character(len=*), parameter :: subname = 'dyn_cnbal_col'
     !-----------------------------------------------------------------------

     call soilbiogeochem_carbonstate_inst%DynamicColumnAdjustments(bounds, clump_index, &
          column_state_updater)
     if (use_c13) then
        call c13_soilbiogeochem_carbonstate_inst%DynamicColumnAdjustments(bounds, clump_index, &
             column_state_updater)
     end if
     if (use_c14) then
        call c14_soilbiogeochem_carbonstate_inst%DynamicColumnAdjustments(bounds, clump_index, &
             column_state_updater)
     end if
     
     call soilbiogeochem_nitrogenstate_inst%DynamicColumnAdjustments(bounds, clump_index, &
          column_state_updater)

     if (use_lch4) then
        call ch4_inst%DynamicColumnAdjustments(bounds, clump_index, column_state_updater)
     end if

   end subroutine dyn_cnbal_col


end module dynConsBiogeochemMod
