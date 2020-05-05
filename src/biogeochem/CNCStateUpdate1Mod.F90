module CNCStateUpdate1Mod

  !-----------------------------------------------------------------------
  ! Module for carbon state variable update, non-mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varpar                         , only : ndecomp_cascade_transitions, nlevdecomp
  use clm_time_manager                   , only : get_step_size, get_step_size_real
  use clm_varpar                         , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use pftconMod                          , only : npcropmin, nc3crop, pftcon
  use abortutils                         , only : endrun
  use decompMod                          , only : bounds_type
  use CNVegCarbonStateType               , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType                , only : cnveg_carbonflux_type
  use CropType                           , only : crop_type
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use SoilBiogeochemCarbonFluxType       , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemCarbonStateType      , only : soilbiogeochem_carbonstate_type
  use PatchType                          , only : patch
  use clm_varctl                         , only : use_fates, use_cn, iulog
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
    real(r8) :: dt  ! time step (seconds)

    character(len=*), parameter :: subname = 'CStateUpdateDynPatch'
    !-----------------------------------------------------------------------

    associate( &
         cf_veg => cnveg_carbonflux_inst  , &
         cs_veg => cnveg_carbonstate_inst , &
         cs_soil => soilbiogeochem_carbonstate_inst &
         )

    dt = get_step_size_real()

    if (.not. use_fates) then
       do j = 1,nlevdecomp
          do fc = 1, num_soilc_with_inactive
             c = filter_soilc_with_inactive(fc)
             cs_soil%decomp_cpools_vr_col(c,j,i_met_lit) = cs_soil%decomp_cpools_vr_col(c,j,i_met_lit) + &
                  cf_veg%dwt_frootc_to_litr_met_c_col(c,j) * dt
             cs_soil%decomp_cpools_vr_col(c,j,i_cel_lit) = cs_soil%decomp_cpools_vr_col(c,j,i_cel_lit) + &
                  cf_veg%dwt_frootc_to_litr_cel_c_col(c,j) * dt
             cs_soil%decomp_cpools_vr_col(c,j,i_lig_lit) = cs_soil%decomp_cpools_vr_col(c,j,i_lig_lit) + &
                  cf_veg%dwt_frootc_to_litr_lig_c_col(c,j) * dt
             cs_soil%decomp_cpools_vr_col(c,j,i_cwd) = cs_soil%decomp_cpools_vr_col(c,j,i_cwd) + &
                  ( cf_veg%dwt_livecrootc_to_cwdc_col(c,j) + cf_veg%dwt_deadcrootc_to_cwdc_col(c,j) ) * dt
          end do
       end do

       do g = bounds%begg, bounds%endg
          cs_veg%seedc_grc(g) = cs_veg%seedc_grc(g) - cf_veg%dwt_seedc_to_leaf_grc(g) * dt
          cs_veg%seedc_grc(g) = cs_veg%seedc_grc(g) - cf_veg%dwt_seedc_to_deadstem_grc(g) * dt
       end do

    end if

    ! TODO(wjs, 2017-01-02) Do we need to move some of the FATES fluxes into here (from
    ! CStateUpdate1) if use_fates is true? Specifically, some portion or all of the fluxes
    ! from these updates in CStateUpdate1:
    ! cf_soil%decomp_cpools_sourcesink_col(c,j,i_met_lit) = cf_soil%FATES_c_to_litr_lab_c_col(c,j) * dt
    ! cf_soil%decomp_cpools_sourcesink_col(c,j,i_cel_lit) = cf_soil%FATES_c_to_litr_cel_c_col(c,j) * dt
    ! cf_soil%decomp_cpools_sourcesink_col(c,j,i_lig_lit) = cf_soil%FATES_c_to_litr_lig_c_col(c,j) * dt

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
      dt = real( get_step_size(), r8 )



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
       soilbiogeochem_carbonflux_inst, dribble_crophrv_xsmrpool_2atm)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic carbon state
    ! variables (except for gap-phase mortality and fire fluxes)
    !
    use clm_varctl, only : carbon_resp_opt
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
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,k,l ! indices
    integer  :: fp,fc     ! filter indices
    real(r8) :: dt        ! radiation time step (seconds)
    real(r8) :: check_cpool
    real(r8) :: cpool_delta
    real(r8), parameter :: kprod05 = 1.44e-7  ! decay constant for 0.5-year product pool (1/s) (lose ~90% over a half year)
    !-----------------------------------------------------------------------

    associate(                                                               & 
         ivt                   => patch%itype                                , & ! Input:  [integer  (:)     ]  patch vegetation type                                

         woody                 => pftcon%woody                             , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)

         cascade_donor_pool    => decomp_cascade_con%cascade_donor_pool    , & ! Input:  [integer  (:)     ]  which pool is C taken from for a given decomposition step
         cascade_receiver_pool => decomp_cascade_con%cascade_receiver_pool , & ! Input:  [integer  (:)     ]  which pool is C added to for a given decomposition step

         harvdate              => crop_inst%harvdate_patch                 , & ! Input:  [integer  (:)     ]  harvest date                                       

         cf_veg                => cnveg_carbonflux_inst                    , & ! Output:
         cs_veg                => cnveg_carbonstate_inst                   , & ! Output:
         cf_soil               => soilbiogeochem_carbonflux_inst             & ! Output:
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      ! Below is the input into the soil biogeochemistry model

      ! plant to litter fluxes
      if (.not. use_fates) then    
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               ! phenology and dynamic land cover fluxes
               cf_soil%decomp_cpools_sourcesink_col(c,j,i_met_lit) = &
                    cf_veg%phenology_c_to_litr_met_c_col(c,j) *dt
               cf_soil%decomp_cpools_sourcesink_col(c,j,i_cel_lit) = &
                    cf_veg%phenology_c_to_litr_cel_c_col(c,j) *dt
               cf_soil%decomp_cpools_sourcesink_col(c,j,i_lig_lit) = &
                    cf_veg%phenology_c_to_litr_lig_c_col(c,j) *dt
   
               ! NOTE(wjs, 2017-01-02) This used to be set to a non-zero value, but the
               ! terms have been moved to CStateUpdateDynPatch. I think this is zeroed every
               ! time step, but to be safe, I'm explicitly setting it to zero here.
               cf_soil%decomp_cpools_sourcesink_col(c,j,i_cwd) = 0._r8
            end do
         end do
      else  !use_fates
         ! here add all fates litterfall and CWD breakdown to litter fluxes
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               ! TODO(wjs, 2017-01-02) Should some portion or all of the following fluxes
               ! be moved to the updates in CStateUpdateDynPatch?
               cf_soil%decomp_cpools_sourcesink_col(c,j,i_met_lit) = cf_soil%FATES_c_to_litr_lab_c_col(c,j) * dt
               cf_soil%decomp_cpools_sourcesink_col(c,j,i_cel_lit) = cf_soil%FATES_c_to_litr_cel_c_col(c,j) * dt
               cf_soil%decomp_cpools_sourcesink_col(c,j,i_lig_lit) = cf_soil%FATES_c_to_litr_lig_c_col(c,j) * dt
            end do
         end do
      endif
         
      ! litter and SOM HR fluxes
      do k = 1, ndecomp_cascade_transitions
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               cf_soil%decomp_cpools_sourcesink_col(c,j,cascade_donor_pool(k)) = &
                    cf_soil%decomp_cpools_sourcesink_col(c,j,cascade_donor_pool(k)) &
                    - ( cf_soil%decomp_cascade_hr_vr_col(c,j,k) + cf_soil%decomp_cascade_ctransfer_vr_col(c,j,k)) *dt
            end do
         end do
      end do
      do k = 1, ndecomp_cascade_transitions
         if ( cascade_receiver_pool(k) /= 0 ) then  ! skip terminal transitions
            do j = 1,nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  cf_soil%decomp_cpools_sourcesink_col(c,j,cascade_receiver_pool(k)) = &
                       cf_soil%decomp_cpools_sourcesink_col(c,j,cascade_receiver_pool(k)) &
                       + cf_soil%decomp_cascade_ctransfer_vr_col(c,j,k)*dt
               end do
            end do
         end if
      end do

      if (.not. use_fates) then    
        do fp = 1,num_soilp
           p = filter_soilp(fp)
           c = patch%column(p)
  
           ! phenology: transfer growth fluxes
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
              cs_veg%grainc_patch(p)          = cs_veg%grainc_patch(p)         + cf_veg%grainc_xfer_to_grainc_patch(p)*dt
              cs_veg%grainc_xfer_patch(p)     = cs_veg%grainc_xfer_patch(p)    - cf_veg%grainc_xfer_to_grainc_patch(p)*dt
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
              cs_veg%grainc_patch(p)     = cs_veg%grainc_patch(p) &
                   - (cf_veg%grainc_to_food_patch(p) + cf_veg%grainc_to_seed_patch(p))*dt
              cs_veg%cropseedc_deficit_patch(p) = cs_veg%cropseedc_deficit_patch(p) &
                   - cf_veg%crop_seedc_to_leaf_patch(p) * dt &
                   + cf_veg%grainc_to_seed_patch(p) * dt
           end if
         
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
              cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%grain_curmr_patch(p)*dt
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
           cs_veg%cpool_patch(p)           = cs_veg%cpool_patch(p)          - cf_veg%cpool_to_leafc_patch(p)*dt
           cs_veg%leafc_patch(p)           = cs_veg%leafc_patch(p)          + cf_veg%cpool_to_leafc_patch(p)*dt
           cs_veg%cpool_patch(p)           = cs_veg%cpool_patch(p)          - cf_veg%cpool_to_leafc_storage_patch(p)*dt
           cs_veg%leafc_storage_patch(p)   = cs_veg%leafc_storage_patch(p)  + cf_veg%cpool_to_leafc_storage_patch(p)*dt
           cs_veg%cpool_patch(p)           = cs_veg%cpool_patch(p)          - cf_veg%cpool_to_frootc_patch(p)*dt
           cs_veg%frootc_patch(p)          = cs_veg%frootc_patch(p)         + cf_veg%cpool_to_frootc_patch(p)*dt
           cs_veg%cpool_patch(p)           = cs_veg%cpool_patch(p)          - cf_veg%cpool_to_frootc_storage_patch(p)*dt

           cs_veg%frootc_storage_patch(p)  = cs_veg%frootc_storage_patch(p) + cf_veg%cpool_to_frootc_storage_patch(p)*dt
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
              cs_veg%livestemc_patch(p)          = cs_veg%livestemc_patch(p)          + cf_veg%cpool_to_livestemc_patch(p)*dt
              cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_livestemc_storage_patch(p)*dt
              cs_veg%livestemc_storage_patch(p)  = cs_veg%livestemc_storage_patch(p)  + cf_veg%cpool_to_livestemc_storage_patch(p)*dt
              cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_deadstemc_patch(p)*dt
              cs_veg%deadstemc_patch(p)          = cs_veg%deadstemc_patch(p)          + cf_veg%cpool_to_deadstemc_patch(p)*dt
              cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_deadstemc_storage_patch(p)*dt
              cs_veg%deadstemc_storage_patch(p)  = cs_veg%deadstemc_storage_patch(p)  + cf_veg%cpool_to_deadstemc_storage_patch(p)*dt
              cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_livecrootc_patch(p)*dt
              cs_veg%livecrootc_patch(p)         = cs_veg%livecrootc_patch(p)         + cf_veg%cpool_to_livecrootc_patch(p)*dt
              cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_livecrootc_storage_patch(p)*dt
              cs_veg%livecrootc_storage_patch(p) = cs_veg%livecrootc_storage_patch(p) + cf_veg%cpool_to_livecrootc_storage_patch(p)*dt
              cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_deadcrootc_patch(p)*dt
              cs_veg%deadcrootc_patch(p)         = cs_veg%deadcrootc_patch(p)         + cf_veg%cpool_to_deadcrootc_patch(p)*dt
              cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_deadcrootc_storage_patch(p)*dt
              cs_veg%deadcrootc_storage_patch(p) = cs_veg%deadcrootc_storage_patch(p) + cf_veg%cpool_to_deadcrootc_storage_patch(p)*dt
           end if
           if (ivt(p) >= npcropmin) then ! skip 2 generic crops
              if (carbon_resp_opt == 1) then
                 cf_veg%cpool_to_livestemc_patch(p) = cf_veg%cpool_to_livestemc_patch(p) - cf_veg%cpool_to_livestemc_resp_patch(p)
                 cf_veg%cpool_to_livestemc_storage_patch(p) = cf_veg%cpool_to_livestemc_storage_patch(p) - &
                 cf_veg%cpool_to_livestemc_storage_resp_patch(p)
              end if
              cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_livestemc_patch(p)*dt
              cs_veg%livestemc_patch(p)          = cs_veg%livestemc_patch(p)          + cf_veg%cpool_to_livestemc_patch(p)*dt
              cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_livestemc_storage_patch(p)*dt
              cs_veg%livestemc_storage_patch(p)  = cs_veg%livestemc_storage_patch(p)  + cf_veg%cpool_to_livestemc_storage_patch(p)*dt
              cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_grainc_patch(p)*dt
              cs_veg%grainc_patch(p)             = cs_veg%grainc_patch(p)             + cf_veg%cpool_to_grainc_patch(p)*dt
              cs_veg%cpool_patch(p)              = cs_veg%cpool_patch(p)              - cf_veg%cpool_to_grainc_storage_patch(p)*dt
              cs_veg%grainc_storage_patch(p)     = cs_veg%grainc_storage_patch(p)     + cf_veg%cpool_to_grainc_storage_patch(p)*dt
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
              cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%cpool_grain_gr_patch(p)*dt
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
              cs_veg%gresp_xfer_patch(p) = cs_veg%gresp_xfer_patch(p) - cf_veg%transfer_grain_gr_patch(p)*dt
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
   
              cs_veg%cpool_patch(p) = cs_veg%cpool_patch(p) - cf_veg%cpool_grain_storage_gr_patch(p)*dt

           end if

           ! growth respiration stored for release during transfer growth
           cs_veg%cpool_patch(p)         = cs_veg%cpool_patch(p)         - cf_veg%cpool_to_gresp_storage_patch(p)*dt
           cs_veg%gresp_storage_patch(p) = cs_veg%gresp_storage_patch(p) + cf_veg%cpool_to_gresp_storage_patch(p)*dt

           ! move storage pools into transfer pools
           cs_veg%leafc_storage_patch(p)  = cs_veg%leafc_storage_patch(p)  - cf_veg%leafc_storage_to_xfer_patch(p)*dt
           cs_veg%leafc_xfer_patch(p)     = cs_veg%leafc_xfer_patch(p)     + cf_veg%leafc_storage_to_xfer_patch(p)*dt
           cs_veg%frootc_storage_patch(p) = cs_veg%frootc_storage_patch(p) - cf_veg%frootc_storage_to_xfer_patch(p)*dt
           cs_veg%frootc_xfer_patch(p)    = cs_veg%frootc_xfer_patch(p)    + cf_veg%frootc_storage_to_xfer_patch(p)*dt
           if (woody(ivt(p)) == 1._r8) then
              cs_veg%livestemc_storage_patch(p)  = cs_veg%livestemc_storage_patch(p) - cf_veg%livestemc_storage_to_xfer_patch(p)*dt
              cs_veg%livestemc_xfer_patch(p)     = cs_veg%livestemc_xfer_patch(p)    + cf_veg%livestemc_storage_to_xfer_patch(p)*dt
              cs_veg%deadstemc_storage_patch(p)  = cs_veg%deadstemc_storage_patch(p) - cf_veg%deadstemc_storage_to_xfer_patch(p)*dt
              cs_veg%deadstemc_xfer_patch(p)     = cs_veg%deadstemc_xfer_patch(p)    + cf_veg%deadstemc_storage_to_xfer_patch(p)*dt
              cs_veg%livecrootc_storage_patch(p) = cs_veg%livecrootc_storage_patch(p)- cf_veg%livecrootc_storage_to_xfer_patch(p)*dt
              cs_veg%livecrootc_xfer_patch(p)    = cs_veg%livecrootc_xfer_patch(p)   + cf_veg%livecrootc_storage_to_xfer_patch(p)*dt
              cs_veg%deadcrootc_storage_patch(p) = cs_veg%deadcrootc_storage_patch(p)- cf_veg%deadcrootc_storage_to_xfer_patch(p)*dt
              cs_veg%deadcrootc_xfer_patch(p)    = cs_veg%deadcrootc_xfer_patch(p)   + cf_veg%deadcrootc_storage_to_xfer_patch(p)*dt
              cs_veg%gresp_storage_patch(p)      = cs_veg%gresp_storage_patch(p)     - cf_veg%gresp_storage_to_xfer_patch(p)*dt
              cs_veg%gresp_xfer_patch(p)         = cs_veg%gresp_xfer_patch(p)        + cf_veg%gresp_storage_to_xfer_patch(p)*dt
           end if
           if (ivt(p) >= npcropmin) then ! skip 2 generic crops
              ! lines here for consistency; the transfer terms are zero
              cs_veg%livestemc_storage_patch(p)  = cs_veg%livestemc_storage_patch(p) - cf_veg%livestemc_storage_to_xfer_patch(p)*dt
              cs_veg%livestemc_xfer_patch(p)     = cs_veg%livestemc_xfer_patch(p)    + cf_veg%livestemc_storage_to_xfer_patch(p)*dt
              cs_veg%grainc_storage_patch(p)     = cs_veg%grainc_storage_patch(p)    - cf_veg%grainc_storage_to_xfer_patch(p)*dt
              cs_veg%grainc_xfer_patch(p)        = cs_veg%grainc_xfer_patch(p)       + cf_veg%grainc_storage_to_xfer_patch(p)*dt
           end if

           if (ivt(p) >= npcropmin) then ! skip 2 generic crops
              cs_veg%xsmrpool_patch(p) = cs_veg%xsmrpool_patch(p) - cf_veg%livestem_xsmr_patch(p)*dt
              cs_veg%xsmrpool_patch(p) = cs_veg%xsmrpool_patch(p) - cf_veg%grain_xsmr_patch(p)*dt
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
  
                 ! Instantly release XSMRPOOL to atmosphere
                 if ( .not. dribble_crophrv_xsmrpool_2atm ) then
                    cf_veg%xsmrpool_to_atm_patch(p) = cf_veg%xsmrpool_to_atm_patch(p) + cs_veg%xsmrpool_patch(p)/dt
                    cf_veg%xsmrpool_to_atm_patch(p) = cf_veg%xsmrpool_to_atm_patch(p) + cs_veg%cpool_patch(p)/dt
                    cf_veg%xsmrpool_to_atm_patch(p) = cf_veg%xsmrpool_to_atm_patch(p) + cs_veg%frootc_patch(p)/dt
                 ! Save xsmrpool, cpool, frootc to loss state variable for dribbling
                 else
                    cs_veg%xsmrpool_loss_patch(p) = cs_veg%xsmrpool_loss_patch(p) + &
                                                    cs_veg%xsmrpool_patch(p) + &
                                                    cs_veg%cpool_patch(p) + &
                                                    cs_veg%frootc_patch(p)
                 end if
                 cs_veg%xsmrpool_patch(p)        = 0._r8
                 cs_veg%cpool_patch(p)           = 0._r8
                 cs_veg%frootc_patch(p)          = 0._r8
              end if

              ! Slowly release xsmrpool to atmosphere
              if ( dribble_crophrv_xsmrpool_2atm ) then
                 ! calculate flux of xsmrpool loss to atm
                 cf_veg%xsmrpool_to_atm_patch(p) = cs_veg%xsmrpool_loss_patch(p) * kprod05
  
                 ! update xsmrpool loss state
                 cs_veg%xsmrpool_loss_patch(p) = cs_veg%xsmrpool_loss_patch(p) - cf_veg%xsmrpool_to_atm_patch(p) * dt
              end if
  
           end if

         
        end do ! end of patch loop
      end if   ! end of NOT fates
    
    end associate
  
  end subroutine CStateUpdate1

end module CNCStateUpdate1Mod
