module CNCStateUpdate1Mod

  !-----------------------------------------------------------------------
  ! Module for carbon state variable update, non-mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varpar                         , only : ndecomp_cascade_transitions, nlevdecomp
  use clm_time_manager                   , only : get_step_size_real
  use clm_varpar                         , only : i_litr_min, i_litr_max, i_cwd
  use clm_varpar                         , only : i_met_lit, i_str_lit, i_phys_som, i_chem_som
  use pftconMod                          , only : npcropmin, nc3crop, pftcon
  use abortutils                         , only : endrun
  use decompMod                          , only : bounds_type
  use CNVegCarbonStateType               , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType                , only : cnveg_carbonflux_type
  use CropType                           , only : crop_type
  use CropReprPoolsMod                   , only : nrepr, repr_grain_min, repr_grain_max
  use CropReprPoolsMod                   , only : repr_structure_min, repr_structure_max
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con, decomp_method, mimics_decomp, use_soil_matrixcn
  use SoilBiogeochemCarbonFluxType       , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemCarbonStateType      , only : soilbiogeochem_carbonstate_type
  use PatchType                          , only : patch
  use CNSharedParamsMod                  , only : use_matrixcn
  use CLMFatesInterfaceMod               , only : hlm_fates_interface_type
  use ColumnType                         , only : col
  
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CStateUpdateDynPatch
  public :: CStateUpdate0
  public :: CStateUpdate1
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CStateUpdateDynPatch(bounds, num_soilc_with_inactive, filter_soilc_with_inactive, &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst, soilbiogeochem_carbonstate_inst)
    !
    ! !DESCRIPTION:
    ! Update carbon states based on fluxes from dyn_cnbal_patch
    ! This routine is not called with FATES active.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds      
    integer, intent(in) :: num_soilc_with_inactive       ! number of columns in soil filter
    integer, intent(in) :: filter_soilc_with_inactive(:) ! soil column filter that includes inactive points
    type(cnveg_carbonflux_type)           , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)          , intent(inout) :: cnveg_carbonstate_inst
    type(soilbiogeochem_carbonstate_type) , intent(inout) :: soilbiogeochem_carbonstate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c   ! column index
    integer  :: fc  ! column filter index
    integer  :: g   ! gridcell index
    integer  :: j   ! level index
    integer  :: i   ! litter pool index
    real(r8) :: dt  ! time step (seconds)

    character(len=*), parameter :: subname = 'CStateUpdateDynPatch'
    !-----------------------------------------------------------------------

    associate( &
         cf_veg => cnveg_carbonflux_inst  , &
         cs_veg => cnveg_carbonstate_inst , &
         cs_soil => soilbiogeochem_carbonstate_inst &
         )

    dt = get_step_size_real()

    do j = 1,nlevdecomp
       do fc = 1, num_soilc_with_inactive
          c = filter_soilc_with_inactive(fc)
          do i = i_litr_min, i_litr_max
             cs_soil%decomp_cpools_vr_col(c,j,i) = &
                  cs_soil%decomp_cpools_vr_col(c,j,i) + &
                  cf_veg%dwt_frootc_to_litr_c_col(c,j,i) * dt
          end do
          cs_soil%decomp_cpools_vr_col(c,j,i_cwd) = cs_soil%decomp_cpools_vr_col(c,j,i_cwd) + &
               ( cf_veg%dwt_livecrootc_to_cwdc_col(c,j) + cf_veg%dwt_deadcrootc_to_cwdc_col(c,j) ) * dt
       end do
    end do

    do g = bounds%begg, bounds%endg
       cs_veg%seedc_grc(g) = cs_veg%seedc_grc(g) - cf_veg%dwt_seedc_to_leaf_grc(g) * dt
       cs_veg%seedc_grc(g) = cs_veg%seedc_grc(g) - cf_veg%dwt_seedc_to_deadstem_grc(g) * dt
    end do

    end associate

  end subroutine CStateUpdateDynPatch

  !-----------------------------------------------------------------------
  subroutine CStateUpdate0(num_soilp, filter_soilp, &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update cpool carbon state
    !
    ! !ARGUMENTS:
    integer                      , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                      , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_carbonflux_type)  , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type) , intent(inout) :: cnveg_carbonstate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: p  ! indices
    integer :: fp ! lake filter indices
    real(r8):: dt ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                             & 
         cf_veg => cnveg_carbonflux_inst , &
         cs_veg => cnveg_carbonstate_inst  &
         )

      ! set time steps
      dt = get_step_size_real()

      ! gross photosynthesis fluxes
      do fp = 1,num_soilp
         p = filter_soilp(fp)
         cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) + cf_veg%psnsun_to_cpool_patch(p)*dt
         cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) + cf_veg%psnshade_to_cpool_patch(p)*dt
      end do

     
    end associate

  end subroutine CStateUpdate0

  !-----------------------------------------------------------------------
  subroutine CStateUpdate1( num_soilc, filter_soilc, num_soilp, filter_soilp, &
       crop_inst, cnveg_carbonflux_inst, cnveg_carbonstate_inst, &
       soilbiogeochem_carbonflux_inst, dribble_crophrv_xsmrpool_2atm, &
       clm_fates, clump_index)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic carbon state
    ! variables (except for gap-phase mortality and fire fluxes)
    !
    use clm_varctl    , only : carbon_resp_opt
    use CNVegMatrixMod, only : matrix_update_phc
    ! !ARGUMENTS:
    integer                              , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                              , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                              , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                              , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(crop_type)                      , intent(in)    :: crop_inst
    type(cnveg_carbonflux_type)          , intent(inout) :: cnveg_carbonflux_inst ! See note below for xsmrpool_to_atm_patch
    type(cnveg_carbonstate_type)         , intent(inout) :: cnveg_carbonstate_inst
    type(soilbiogeochem_carbonflux_type) , intent(inout) :: soilbiogeochem_carbonflux_inst
    logical                              , intent(in)    :: dribble_crophrv_xsmrpool_2atm
    type(hlm_fates_interface_type)       , intent(inout) :: clm_fates
    integer                              , intent(in)    :: clump_index
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,k,l,i  ! indices
    integer  :: fp,fc     ! filter indices
    real(r8) :: dt        ! radiation time step (seconds)
    real(r8) :: check_cpool
    real(r8) :: cpool_delta
    real(r8), parameter :: kprod05 = 1.44e-7_r8  ! decay constant for 0.5-year product pool (1/s) (lose ~90% over a half year)
    !-----------------------------------------------------------------------

    associate(                                                               & 
         ivt                   => patch%itype                                , & ! Input:  [integer  (:)     ]  patch vegetation type                                

         mimics_fi             => pftcon%mimics_fi                         , & ! Input: MIMICS parameter fi
         woody                 => pftcon%woody                             , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)

         cascade_donor_pool    => decomp_cascade_con%cascade_donor_pool    , & ! Input:  [integer  (:)     ]  which pool is C taken from for a given decomposition step
         cascade_receiver_pool => decomp_cascade_con%cascade_receiver_pool , & ! Input:  [integer  (:)     ]  which pool is C added to for a given decomposition step

         harvdate              => crop_inst%harvdate_patch                 , & ! Input:  [integer  (:)     ]  harvest date                                       

         cf_veg                => cnveg_carbonflux_inst                    , & ! Output:
         cs_veg                => cnveg_carbonstate_inst                   , & ! Output:
         cf_soil               => soilbiogeochem_carbonflux_inst             & ! Output:
         )

      ! set time steps
      dt = get_step_size_real()

      ! Below is the input into the soil biogeochemistry model

      fc_loop: do fc = 1,num_soilc
         c = filter_soilc(fc)

         fates_if: if( col%is_fates(c) ) then

            ! If this is a fates column, then we ask fates for the
            ! litter fluxes, the following routine simply copies
            ! prepared litter c flux boundary conditions into
            ! cf_soil%decomp_cpools_sourcesink_col
            
            call clm_fates%UpdateCLitterfluxes(cf_soil,clump_index,c)
            
         else
            do j = 1,nlevdecomp
               !
               ! State update without the matrix solution
               !
               if (.not. use_soil_matrixcn) then
                  ! phenology and dynamic land cover fluxes
                  if (decomp_method == mimics_decomp) then
                     do i = i_litr_min, i_litr_max  ! in MIMICS these are 1 and 2
                        cf_soil%decomp_cpools_sourcesink_col(c,j,i) = (1 - mimics_fi(i)) * &
                           cf_veg%phenology_c_to_litr_c_col(c,j,i) * dt
                     end do
                     cf_soil%decomp_cpools_sourcesink_col(c,j,i_phys_som) = mimics_fi(1) * &
                        cf_veg%phenology_c_to_litr_c_col(c,j,i_met_lit) * dt
                     cf_soil%decomp_cpools_sourcesink_col(c,j,i_chem_som) = mimics_fi(2) * &
                        cf_veg%phenology_c_to_litr_c_col(c,j,i_str_lit) * dt
                  else
                     do i = i_litr_min, i_litr_max
                        cf_soil%decomp_cpools_sourcesink_col(c,j,i) = &
                             cf_veg%phenology_c_to_litr_c_col(c,j,i) * dt
                     end do
                  end if

                  ! NOTE(wjs, 2017-01-02) This used to be set to a non-zero value, but the
                  ! terms have been moved to CStateUpdateDynPatch. I think this is zeroed every
                  ! time step, but to be safe, I'm explicitly setting it to zero here.
                  cf_soil%decomp_cpools_sourcesink_col(c,j,i_cwd) = 0._r8

               else
                  !
                  ! For the matrix solution the actual state update comes after the matrix
                  ! multiply in SoilMatrix, but the matrix needs to be setup with
                  ! the equivalent of above. Those changes can be here or in the
                  ! native subroutines dealing with that field
                  !
                  ! phenology and dynamic land cover fluxes
                  do i = i_litr_min, i_litr_max
                     cf_soil%matrix_Cinput%V(c,j+(i-1)*nlevdecomp) = &
                          cf_soil%matrix_Cinput%V(c,j+(i-1)*nlevdecomp) + cf_veg%phenology_c_to_litr_c_col(c,j,i) *dt
                  end do
               end if
            end do

         end if fates_if

         do j = 1,nlevdecomp
            !
            ! State update without the matrix solution
            !
            if (.not. use_soil_matrixcn) then
               ! litter and SOM HR fluxes
               do k = 1, ndecomp_cascade_transitions
                  cf_soil%decomp_cpools_sourcesink_col(c,j,cascade_donor_pool(k)) = &
                       cf_soil%decomp_cpools_sourcesink_col(c,j,cascade_donor_pool(k)) &
                       - ( cf_soil%decomp_cascade_hr_vr_col(c,j,k) + cf_soil%decomp_cascade_ctransfer_vr_col(c,j,k)) * dt
                  if ( cascade_receiver_pool(k) /= 0 ) then  ! skip terminal transitions
                     cf_soil%decomp_cpools_sourcesink_col(c,j,cascade_receiver_pool(k)) = &
                          cf_soil%decomp_cpools_sourcesink_col(c,j,cascade_receiver_pool(k)) &
                          + cf_soil%decomp_cascade_ctransfer_vr_col(c,j,k) * dt
                  end if
               end do
            end if
         end do
            
      end do fc_loop

      soilpatch_loop: do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = patch%column(p)

         ! phenology: transfer growth fluxes
         ! TODO slevis: improve indentation
        if(.not. use_matrixcn)then
           ! NOTE: Any changes that go here MUST be applied to the matrix
           ! version as well
           cs_veg%leafc_patch(p)           = cs_veg%leafc_patch(p)       + cf_veg%leafc_xfer_to_leafc_patch(p)*dt
           cs_veg%leafc_xfer_patch(p)      = cs_veg%leafc_xfer_patch(p)  - cf_veg%leafc_xfer_to_leafc_patch(p)*dt
           cs_veg%frootc_patch(p)          = cs_veg%frootc_patch(p)      + cf_veg%frootc_xfer_to_frootc_patch(p)*dt
           cs_veg%frootc_xfer_patch(p)     = cs_veg%frootc_xfer_patch(p) - cf_veg%frootc_xfer_to_frootc_patch(p)*dt
           if (woody(ivt(p)) == 1._r8) then
              cs_veg%livestemc_patch(p)       = cs_veg%livestemc_patch(p)       + cf_veg%livestemc_xfer_to_livestemc_patch(p)*dt
              cs_veg%livestemc_xfer_patch(p)  = cs_veg%livestemc_xfer_patch(p)  - cf_veg%livestemc_xfer_to_livestemc_patch(p)*dt
              cs_veg%deadstemc_patch(p)       = cs_veg%deadstemc_patch(p)       + cf_veg%deadstemc_xfer_to_deadstemc_patch(p)*dt
              cs_veg%deadstemc_xfer_patch(p)  = cs_veg%deadstemc_xfer_patch(p)  - cf_veg%deadstemc_xfer_to_deadstemc_patch(p)*dt
              cs_veg%livecrootc_patch(p)      = cs_veg%livecrootc_patch(p)      + cf_veg%livecrootc_xfer_to_livecrootc_patch(p)*dt
              cs_veg%livecrootc_xfer_patch(p) = cs_veg%livecrootc_xfer_patch(p) - cf_veg%livecrootc_xfer_to_livecrootc_patch(p)*dt
              cs_veg%deadcrootc_patch(p)      = cs_veg%deadcrootc_patch(p)      + cf_veg%deadcrootc_xfer_to_deadcrootc_patch(p)*dt
              cs_veg%deadcrootc_xfer_patch(p) = cs_veg%deadcrootc_xfer_patch(p) - cf_veg%deadcrootc_xfer_to_deadcrootc_patch(p)*dt
           end if
           if (ivt(p) >= npcropmin) then ! skip 2 generic crops
             ! lines here for consistency; the transfer terms are zero
              cs_veg%livestemc_patch(p)       = cs_veg%livestemc_patch(p)      + cf_veg%livestemc_xfer_to_livestemc_patch(p)*dt
              cs_veg%livestemc_xfer_patch(p)  = cs_veg%livestemc_xfer_patch(p) - cf_veg%livestemc_xfer_to_livestemc_patch(p)*dt
              do k = 1, nrepr
                 cs_veg%reproductivec_patch(p,k) = cs_veg%reproductivec_patch(p,k) &
                      + cf_veg%reproductivec_xfer_to_reproductivec_patch(p,k)*dt
                 cs_veg%reproductivec_xfer_patch(p,k) = cs_veg%reproductivec_xfer_patch(p,k) &
                      - cf_veg%reproductivec_xfer_to_reproductivec_patch(p,k)*dt
              end do
           end if

         ! phenology: litterfall fluxes
           cs_veg%leafc_patch(p) = cs_veg%leafc_patch(p) - cf_veg%leafc_to_litter_patch(p)*dt
           cs_veg%frootc_patch(p) = cs_veg%frootc_patch(p) - cf_veg%frootc_to_litter_patch(p)*dt
         
         ! livewood turnover fluxes
           if (woody(ivt(p)) == 1._r8) then
              cs_veg%livestemc_patch(p)  = cs_veg%livestemc_patch(p)  - cf_veg%livestemc_to_deadstemc_patch(p)*dt
              cs_veg%deadstemc_patch(p)  = cs_veg%deadstemc_patch(p)  + cf_veg%livestemc_to_deadstemc_patch(p)*dt
              cs_veg%livecrootc_patch(p) = cs_veg%livecrootc_patch(p) - cf_veg%livecrootc_to_deadcrootc_patch(p)*dt
              cs_veg%deadcrootc_patch(p) = cs_veg%deadcrootc_patch(p) + cf_veg%livecrootc_to_deadcrootc_patch(p)*dt
           end if
           if (ivt(p) >= npcropmin) then ! skip 2 generic crops
              cs_veg%livestemc_patch(p)  = cs_veg%livestemc_patch(p)  - cf_veg%livestemc_to_litter_patch(p)*dt
              cs_veg%livestemc_patch(p)  = cs_veg%livestemc_patch(p)  - &
                   (cf_veg%livestemc_to_biofuelc_patch(p) + cf_veg%livestemc_to_removedresiduec_patch(p))*dt
              cs_veg%leafc_patch(p)      = cs_veg%leafc_patch(p)      - &
                   (cf_veg%leafc_to_biofuelc_patch(p) + cf_veg%leafc_to_removedresiduec_patch(p))*dt
              cs_veg%cropseedc_deficit_patch(p) = cs_veg%cropseedc_deficit_patch(p) &
                   - cf_veg%crop_seedc_to_leaf_patch(p) * dt
              do k = repr_grain_min, repr_grain_max
                 cs_veg%reproductivec_patch(p,k)   = cs_veg%reproductivec_patch(p,k) &
                      - (cf_veg%repr_grainc_to_food_patch(p,k) + cf_veg%repr_grainc_to_seed_patch(p,k))*dt
                 cs_veg%cropseedc_deficit_patch(p) = cs_veg%cropseedc_deficit_patch(p) &
                      + cf_veg%repr_grainc_to_seed_patch(p,k) * dt
              end do
              do k = repr_structure_min, repr_structure_max
                 cs_veg%reproductivec_patch(p,k) = cs_veg%reproductivec_patch(p,k) &
                      - (cf_veg%repr_structurec_to_cropprod_patch(p,k) + cf_veg%repr_structurec_to_litter_patch(p,k))*dt
              end do
           end if
        else
           ! NOTE: Changes for above that apply for matrix code are in CNPhenology EBK (11/26/2019)

           ! This part below MUST match exactly the code for the non-matrix part
           ! above!
           if (ivt(p) >= npcropmin) then
              cs_veg%cropseedc_deficit_patch(p) = cs_veg%cropseedc_deficit_patch(p) &
                   - cf_veg%crop_seedc_to_leaf_patch(p) * dt
              do k = repr_grain_min, repr_grain_max
                 cs_veg%cropseedc_deficit_patch(p) = cs_veg%cropseedc_deficit_patch(p) &
                      + cf_veg%repr_grainc_to_seed_patch(p,k) * dt
              end do
           end if
        end if !not use_matrixcn

        check_cpool = cs_veg%cpool_patch(p)- cf_veg%psnsun_to_cpool_patch(p)*dt-cf_veg%psnshade_to_cpool_patch(p)*dt
        cpool_delta  =  cs_veg%cpool_patch(p) 

           ! maintenance respiration fluxes from cpool

           cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%cpool_to_xsmrpool_patch(p)*dt
           cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%leaf_curmr_patch(p)*dt
           cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%froot_curmr_patch(p)*dt
           If (woody(ivt(p)) == 1._r8) then
              cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%livestem_curmr_patch(p)*dt
              cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%livecroot_curmr_patch(p)*dt
           end if
           if (ivt(p) >= npcropmin) then ! skip 2 generic crops
              cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%livestem_curmr_patch(p)*dt
              do k = 1, nrepr
                 cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%reproductive_curmr_patch(p,k)*dt
              end do
           end if
         
         
           cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) -  cf_veg%cpool_to_resp_patch(p)*dt

           !RF Add in the carbon spent on uptake respiration. 
           cs_veg%cpool_patch(p)= cs_veg%cpool_patch(p) - cf_veg%soilc_change_patch(p)*dt
           
           ! maintenance respiration fluxes from xsmrpool
           cs_veg%xsmrpool_patch(p) = cs_veg%xsmrpool_patch(p) + cf_veg%cpool_to_xsmrpool_patch(p)*dt
           cs_veg%xsmrpool_patch(p) = cs_veg%xsmrpool_patch(p) - cf_veg%leaf_xsmr_patch(p)*dt
           cs_veg%xsmrpool_patch(p) = cs_veg%xsmrpool_patch(p) - cf_veg%froot_xsmr_patch(p)*dt
           if (woody(ivt(p)) == 1._r8) then
              cs_veg%xsmrpool_patch(p) = cs_veg%xsmrpool_patch(p) - cf_veg%livestem_xsmr_patch(p)*dt
              cs_veg%xsmrpool_patch(p) = cs_veg%xsmrpool_patch(p) - cf_veg%livecroot_xsmr_patch(p)*dt
           end if

           ! allocation fluxes
           if (carbon_resp_opt == 1) then
              cf_veg%cpool_to_leafc_patch(p) = cf_veg%cpool_to_leafc_patch(p) - cf_veg%cpool_to_leafc_resp_patch(p)
              cf_veg%cpool_to_leafc_storage_patch(p) = cf_veg%cpool_to_leafc_storage_patch(p) - &
                   cf_veg%cpool_to_leafc_storage_resp_patch(p)
              cf_veg%cpool_to_frootc_patch(p) = cf_veg%cpool_to_frootc_patch(p) - cf_veg%cpool_to_frootc_resp_patch(p)
              cf_veg%cpool_to_frootc_storage_patch(p) = cf_veg%cpool_to_frootc_storage_patch(p) &
                 - cf_veg%cpool_to_frootc_storage_resp_patch(p)
         end if
         cs_veg%cpool_patch(p)             = cs_veg%cpool_patch(p)          - cf_veg%cpool_to_leafc_patch(p)*dt
         cs_veg%cpool_patch(p)             = cs_veg%cpool_patch(p)          - cf_veg%cpool_to_leafc_storage_patch(p)*dt
         cs_veg%cpool_patch(p)             = cs_veg%cpool_patch(p)          - cf_veg%cpool_to_frootc_patch(p)*dt
         cs_veg%cpool_patch(p)             = cs_veg%cpool_patch(p)          - cf_veg%cpool_to_frootc_storage_patch(p)*dt 
        if(.not. use_matrixcn) then
           cs_veg%leafc_patch(p)           = cs_veg%leafc_patch(p)          + cf_veg%cpool_to_leafc_patch(p)*dt
           cs_veg%leafc_storage_patch(p)   = cs_veg%leafc_storage_patch(p)  + cf_veg%cpool_to_leafc_storage_patch(p)*dt
           cs_veg%frootc_patch(p)          = cs_veg%frootc_patch(p)         + cf_veg%cpool_to_frootc_patch(p)*dt
           cs_veg%frootc_storage_patch(p)  = cs_veg%frootc_storage_patch(p) + cf_veg%cpool_to_frootc_storage_patch(p)*dt
         else
           ! NOTE: The equivalent changes for matrix code are in CNPhenology EBK (11/26/2019)
         end if !not use_matrixcn
         if (woody(ivt(p)) == 1._r8) then
            if (carbon_resp_opt == 1) then
               cf_veg%cpool_to_livecrootc_patch(p) = cf_veg%cpool_to_livecrootc_patch(p) - cf_veg%cpool_to_livecrootc_resp_patch(p)
               cf_veg%cpool_to_livecrootc_storage_patch(p) = cf_veg%cpool_to_livecrootc_storage_patch(p) - &
                    cf_veg%cpool_to_livecrootc_storage_resp_patch(p)
               cf_veg%cpool_to_livestemc_patch(p) = cf_veg%cpool_to_livestemc_patch(p) - cf_veg%cpool_to_livestemc_resp_patch(p)
               cf_veg%cpool_to_livestemc_storage_patch(p) = cf_veg%cpool_to_livestemc_storage_patch(p) - &
                 cf_veg%cpool_to_livestemc_storage_resp_patch(p)
            end if
            cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_livestemc_patch(p)*dt
            cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_livestemc_storage_patch(p)*dt
            cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_deadstemc_patch(p)*dt
            cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_deadstemc_storage_patch(p)*dt
            cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_livecrootc_patch(p)*dt
            cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_livecrootc_storage_patch(p)*dt
            cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_deadcrootc_patch(p)*dt
            cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_deadcrootc_storage_patch(p)*dt 
            if(.not. use_matrixcn)then
               cs_veg%livestemc_patch(p)          = cs_veg%livestemc_patch(p)          + cf_veg%cpool_to_livestemc_patch(p)*dt
               cs_veg%livestemc_storage_patch(p)  = cs_veg%livestemc_storage_patch(p)  + cf_veg%cpool_to_livestemc_storage_patch(p)*dt
               cs_veg%deadstemc_patch(p)          = cs_veg%deadstemc_patch(p)          + cf_veg%cpool_to_deadstemc_patch(p)*dt
               cs_veg%deadstemc_storage_patch(p)  = cs_veg%deadstemc_storage_patch(p)  + cf_veg%cpool_to_deadstemc_storage_patch(p)*dt
               cs_veg%livecrootc_patch(p)         = cs_veg%livecrootc_patch(p)         + cf_veg%cpool_to_livecrootc_patch(p)*dt
               cs_veg%livecrootc_storage_patch(p) = cs_veg%livecrootc_storage_patch(p) + cf_veg%cpool_to_livecrootc_storage_patch(p)*dt
               cs_veg%deadcrootc_patch(p)         = cs_veg%deadcrootc_patch(p)         + cf_veg%cpool_to_deadcrootc_patch(p)*dt
               cs_veg%deadcrootc_storage_patch(p) = cs_veg%deadcrootc_storage_patch(p) + cf_veg%cpool_to_deadcrootc_storage_patch(p)*dt
            else
               ! NOTE: The equivalent changes for matrix code are in CNPhenology EBK (11/26/2019)
            end if !not use_matrixcn
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            if (carbon_resp_opt == 1) then
               cf_veg%cpool_to_livestemc_patch(p) = cf_veg%cpool_to_livestemc_patch(p) - cf_veg%cpool_to_livestemc_resp_patch(p)
               cf_veg%cpool_to_livestemc_storage_patch(p) = cf_veg%cpool_to_livestemc_storage_patch(p) - &
                 cf_veg%cpool_to_livestemc_storage_resp_patch(p)
            end if
            cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_livestemc_patch(p)*dt
            cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_livestemc_storage_patch(p)*dt
            do k = 1, nrepr
                cs_veg%cpool_patch(p)          = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_reproductivec_patch(p,k)*dt
                cs_veg%cpool_patch(p)          = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_reproductivec_storage_patch(p,k)*dt
            end do
            if(.not. use_matrixcn)then
               cs_veg%livestemc_patch(p)          = cs_veg%livestemc_patch(p)          + cf_veg%cpool_to_livestemc_patch(p)*dt
               cs_veg%livestemc_storage_patch(p)  = cs_veg%livestemc_storage_patch(p)  + cf_veg%cpool_to_livestemc_storage_patch(p)*dt
               do k = 1, nrepr
                  cs_veg%reproductivec_patch(p,k)         = cs_veg%reproductivec_patch(p,k) &
                                                          + cf_veg%cpool_to_reproductivec_patch(p,k)*dt
                  cs_veg%reproductivec_storage_patch(p,k) = cs_veg%reproductivec_storage_patch(p,k) &
                                                          + cf_veg%cpool_to_reproductivec_storage_patch(p,k)*dt
               end do
            else
               ! NOTE: The equivalent changes for matrix code are in CNPhenology EBK (11/26/2019)
            end if  !not use_matrixcn
         end if

         ! growth respiration fluxes for current growth
         cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%cpool_leaf_gr_patch(p)*dt
         cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%cpool_froot_gr_patch(p)*dt

         if (woody(ivt(p)) == 1._r8) then
            cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%cpool_livestem_gr_patch(p)*dt
            cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%cpool_deadstem_gr_patch(p)*dt
            cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%cpool_livecroot_gr_patch(p)*dt
            cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%cpool_deadcroot_gr_patch(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%cpool_livestem_gr_patch(p)*dt
            do k = 1, nrepr
               cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%cpool_reproductive_gr_patch(p,k)*dt
            end do
         end if

         ! growth respiration for transfer growth
         cs_veg%gresp_xfer_patch(p) = cs_veg%gresp_xfer_patch(p) - cf_veg%transfer_leaf_gr_patch(p)*dt
         cs_veg%gresp_xfer_patch(p) = cs_veg%gresp_xfer_patch(p) - cf_veg%transfer_froot_gr_patch(p)*dt
         if (woody(ivt(p)) == 1._r8) then
            cs_veg%gresp_xfer_patch(p) = cs_veg%gresp_xfer_patch(p) - cf_veg%transfer_livestem_gr_patch(p)*dt
            cs_veg%gresp_xfer_patch(p) = cs_veg%gresp_xfer_patch(p) - cf_veg%transfer_deadstem_gr_patch(p)*dt
            cs_veg%gresp_xfer_patch(p) = cs_veg%gresp_xfer_patch(p) - cf_veg%transfer_livecroot_gr_patch(p)*dt
            cs_veg%gresp_xfer_patch(p) = cs_veg%gresp_xfer_patch(p) - cf_veg%transfer_deadcroot_gr_patch(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            cs_veg%gresp_xfer_patch(p) = cs_veg%gresp_xfer_patch(p) - cf_veg%transfer_livestem_gr_patch(p)*dt
            do k = 1, nrepr
               cs_veg%gresp_xfer_patch(p) = cs_veg%gresp_xfer_patch(p) - cf_veg%transfer_reproductive_gr_patch(p,k)*dt
            end do
         end if

         ! growth respiration at time of storage
         cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%cpool_leaf_storage_gr_patch(p)*dt
         cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%cpool_froot_storage_gr_patch(p)*dt

         if (woody(ivt(p)) == 1._r8) then
            cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%cpool_livestem_storage_gr_patch(p)*dt
            cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%cpool_deadstem_storage_gr_patch(p)*dt
            cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%cpool_livecroot_storage_gr_patch(p)*dt
            cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%cpool_deadcroot_storage_gr_patch(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%cpool_livestem_storage_gr_patch(p)*dt

            do k = 1, nrepr
               cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%cpool_reproductive_storage_gr_patch(p,k)*dt
            end do

         end if

         ! growth respiration stored for release during transfer growth
         cs_veg%cpool_patch(p)         = cs_veg%cpool_patch(p)         - cf_veg%cpool_to_gresp_storage_patch(p)*dt
         cs_veg%gresp_storage_patch(p) = cs_veg%gresp_storage_patch(p) + cf_veg%cpool_to_gresp_storage_patch(p)*dt

         ! move storage pools into transfer pools
        if(.not. use_matrixcn)then
          cs_veg%leafc_storage_patch(p)  = cs_veg%leafc_storage_patch(p)  - cf_veg%leafc_storage_to_xfer_patch(p)*dt
          cs_veg%leafc_xfer_patch(p)     = cs_veg%leafc_xfer_patch(p)     + cf_veg%leafc_storage_to_xfer_patch(p)*dt
          cs_veg%frootc_storage_patch(p) = cs_veg%frootc_storage_patch(p) - cf_veg%frootc_storage_to_xfer_patch(p)*dt
          cs_veg%frootc_xfer_patch(p)    = cs_veg%frootc_xfer_patch(p)    + cf_veg%frootc_storage_to_xfer_patch(p)*dt
         else
           ! NOTE: The equivalent changes for matrix code are in CNPhenology EBK (11/26/2019)
         end if !not use_matrixcn
         if (woody(ivt(p)) == 1._r8) then
            cs_veg%gresp_storage_patch(p)      = cs_veg%gresp_storage_patch(p)     - cf_veg%gresp_storage_to_xfer_patch(p)*dt
            cs_veg%gresp_xfer_patch(p)         = cs_veg%gresp_xfer_patch(p)        + cf_veg%gresp_storage_to_xfer_patch(p)*dt
            if(.not. use_matrixcn)then
               cs_veg%livestemc_storage_patch(p)  = cs_veg%livestemc_storage_patch(p) - cf_veg%livestemc_storage_to_xfer_patch(p)*dt
               cs_veg%livestemc_xfer_patch(p)     = cs_veg%livestemc_xfer_patch(p)    + cf_veg%livestemc_storage_to_xfer_patch(p)*dt
               cs_veg%deadstemc_storage_patch(p)  = cs_veg%deadstemc_storage_patch(p) - cf_veg%deadstemc_storage_to_xfer_patch(p)*dt
               cs_veg%deadstemc_xfer_patch(p)     = cs_veg%deadstemc_xfer_patch(p)    + cf_veg%deadstemc_storage_to_xfer_patch(p)*dt
               cs_veg%livecrootc_storage_patch(p) = cs_veg%livecrootc_storage_patch(p)- cf_veg%livecrootc_storage_to_xfer_patch(p)*dt
               cs_veg%livecrootc_xfer_patch(p)    = cs_veg%livecrootc_xfer_patch(p)   + cf_veg%livecrootc_storage_to_xfer_patch(p)*dt
               cs_veg%deadcrootc_storage_patch(p) = cs_veg%deadcrootc_storage_patch(p)- cf_veg%deadcrootc_storage_to_xfer_patch(p)*dt
               cs_veg%deadcrootc_xfer_patch(p)    = cs_veg%deadcrootc_xfer_patch(p)   + cf_veg%deadcrootc_storage_to_xfer_patch(p)*dt
            else
               ! NOTE: The equivalent changes for matrix code are in CNPhenology EBK (11/26/2019)
            end if !not use_matrixcn
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ! lines here for consistency; the transfer terms are zero
            if(.not. use_matrixcn)then
               ! lines here for consistency; the transfer terms are zero
               cs_veg%livestemc_storage_patch(p)  = cs_veg%livestemc_storage_patch(p) - cf_veg%livestemc_storage_to_xfer_patch(p)*dt
               cs_veg%livestemc_xfer_patch(p)     = cs_veg%livestemc_xfer_patch(p)    + cf_veg%livestemc_storage_to_xfer_patch(p)*dt
               do k = 1, nrepr
                  cs_veg%reproductivec_storage_patch(p,k) = cs_veg%reproductivec_storage_patch(p,k) &
                       - cf_veg%reproductivec_storage_to_xfer_patch(p,k)*dt
                  cs_veg%reproductivec_xfer_patch(p,k) = cs_veg%reproductivec_xfer_patch(p,k) &
                       + cf_veg%reproductivec_storage_to_xfer_patch(p,k)*dt
               end do
            else
               ! NOTE: The equivalent changes for matrix code are in CNPhenology EBK (11/26/2019)
            end if !not use_matrixcn
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            cs_veg%xsmrpool_patch(p) = cs_veg%xsmrpool_patch(p) - cf_veg%livestem_xsmr_patch(p)*dt
            do k = 1, nrepr
               cs_veg%xsmrpool_patch(p) = cs_veg%xsmrpool_patch(p) - cf_veg%reproductive_xsmr_patch(p,k)*dt
            end do
            if (harvdate(p) < 999) then ! beginning at harvest, send to atm
               ! TODO (mv, 11-02-2014) the following lines are why the cf_veg is
               ! an intent(inout)
               ! fluxes should not be updated in this module - not sure where
               ! this belongs
               ! DML (06-20-2017) While debugging crop isotope code, found that cpool_patch and frootc_patch 
               ! could occasionally be very small but nonzero numbers after crop harvest, which persists 
               ! through to next planting and for reasons that could not 100%
               ! isolate, caused C12/C13 ratios to occasionally go out of
               ! bounds. Zeroing out these small pools and putting them into the flux to the
               ! atmosphere solved many of the crop isotope problems

               if ( .not. dribble_crophrv_xsmrpool_2atm ) then
                  cf_veg%xsmrpool_to_atm_patch(p) = cf_veg%xsmrpool_to_atm_patch(p) + cs_veg%xsmrpool_patch(p)/dt
                  cf_veg%xsmrpool_to_atm_patch(p) = cf_veg%xsmrpool_to_atm_patch(p) + cs_veg%cpool_patch(p)/dt
                  if(.not. use_matrixcn)then
                     cf_veg%xsmrpool_to_atm_patch(p) = cf_veg%xsmrpool_to_atm_patch(p) + cs_veg%frootc_patch(p)/dt
                  else
                     cf_veg%xsmrpool_to_atm_patch(p) = cf_veg%xsmrpool_to_atm_patch(p) &
                       + cs_veg%frootc_patch(p) * matrix_update_phc(p,cf_veg%ifroot_to_iout_ph,1._r8/dt,dt,cnveg_carbonflux_inst,.true.,.true.)
                  end if
                  ! Save xsmrpool, cpool, frootc to loss state variable for
                  ! dribbling
               else
                  ! EBK: 10/08/2020 this could potentially change answers by
                  ! roundoff relative to the baseline (becuase frootc isn't
                  ! alsto subtracted here)
                  cs_veg%xsmrpool_loss_patch(p) = cs_veg%xsmrpool_loss_patch(p) + &
                                                  cs_veg%xsmrpool_patch(p) + &
                                                  cs_veg%cpool_patch(p)
                  if(.not. use_matrixcn)then
                     cs_veg%xsmrpool_loss_patch(p) = cs_veg%xsmrpool_loss_patch(p) + cs_veg%frootc_patch(p)
                  else
                     cs_veg%xsmrpool_loss_patch(p) = cs_veg%xsmrpool_loss_patch(p) &
                       + cs_veg%frootc_patch(p) * matrix_update_phc(p,cf_veg%ifroot_to_iout_ph,1._r8/dt,dt,cnveg_carbonflux_inst,.true.,.true.)
                  end if
               end if
               if (.not. use_matrixcn) then
                  cs_veg%frootc_patch(p)          = 0._r8
               end if
               cs_veg%xsmrpool_patch(p)        = 0._r8
               cs_veg%cpool_patch(p)           = 0._r8
            end if
            ! Slowly release xsmrpool to atmosphere
            if ( dribble_crophrv_xsmrpool_2atm ) then
               ! calculate flux of xsmrpool loss to atm
               cf_veg%xsmrpool_to_atm_patch(p) = cs_veg%xsmrpool_loss_patch(p) * kprod05

               ! update xsmrpool loss state
               cs_veg%xsmrpool_loss_patch(p) = cs_veg%xsmrpool_loss_patch(p) - cf_veg%xsmrpool_to_atm_patch(p) * dt
            end if
         end if
      end do soilpatch_loop ! end of patch loop
    
    end associate
  
  end subroutine CStateUpdate1

end module CNCStateUpdate1Mod
