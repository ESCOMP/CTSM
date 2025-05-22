   module CNCStateUpdate2Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for carbon state variable update, mortality fluxes.
  ! When the matrix solution is being used (use_matrixcn and use_soil_matrixcn)
  ! only some state updates are done here, the other state updates happen
  ! after the matrix is solved in VegMatrix and SoilMatrix.
  !
  ! !USES:
  use shr_kind_mod                   , only : r8 => shr_kind_r8
  use shr_log_mod                    , only : errMsg => shr_log_errMsg
  use abortutils                     , only : endrun
  use clm_time_manager               , only : get_step_size_real
  use clm_varpar                     , only : i_litr_min, i_litr_max, nlevdecomp, i_cwd
  use CNvegCarbonStateType           , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType            , only : cnveg_carbonflux_type
  use SoilBiogeochemCarbonStatetype  , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemCarbonFluxtype   , only : soilbiogeochem_carbonflux_type
  use CNSharedParamsMod              , only : use_matrixcn
  use SoilBiogeochemDecompCascadeConType , only : use_soil_matrixcn
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CStateUpdate2
  public:: CStateUpdate2h
  public:: CStateUpdate2g
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst, soilbiogeochem_carbonstate_inst, &
       soilbiogeochem_carbonflux_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic carbon state
    ! variables affected by gap-phase mortality fluxes
    !
    ! !ARGUMENTS:
    integer                                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_carbonflux_type)            , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)           , intent(inout) :: cnveg_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)   , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)  , intent(inout) :: soilbiogeochem_carbonstate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,i  ! indices
    integer  :: fp,fc  ! lake filter indices
    real(r8) :: dt     ! radiation time step (seconds)
    !-----------------------------------------------------------------------
    
    associate(                                      & 
         cf_veg => cnveg_carbonflux_inst          , &
         cs_veg => cnveg_carbonstate_inst         , &

         cf_soil => soilbiogeochem_carbonflux_inst, &
         cs_soil => soilbiogeochem_carbonstate_inst  &
         )

      ! set time steps
      dt = get_step_size_real()

      ! column level carbon fluxes from gap-phase mortality
      do j = 1,nlevdecomp
         ! column loop
         do fc = 1,num_soilc
            c = filter_soilc(fc)

            ! column gap mortality fluxes

            !
            ! State update without the matrix solution
            !
            if (.not. use_soil_matrixcn)then
               do i = i_litr_min, i_litr_max
                  cs_soil%decomp_cpools_vr_col(c,j,i) = &
                    cs_soil%decomp_cpools_vr_col(c,j,i) + &
                    cf_veg%gap_mortality_c_to_litr_c_col(c,j,i) * dt
               end do
               ! Currently i_cwd .ne. i_litr_max + 1 if .not. fates and
               !           i_cwd = 0 if fates, so not including in the i-loop
               cs_soil%decomp_cpools_vr_col(c,j,i_cwd) = &
                 cs_soil%decomp_cpools_vr_col(c,j,i_cwd) + cf_veg%gap_mortality_c_to_cwdc_col(c,j) * dt
            !
            ! For the matrix solution the actual state update comes after the matrix
            ! multiply in SoilMatrix, but the matrix needs to be setup with
            ! the equivalent of above. Those changes can be here or in the
            ! native subroutines dealing with that field
            !
            else
              ! Match above for soil-matrix
               do i = i_litr_min, i_litr_max
                  cf_soil%matrix_Cinput%V(c,j+(i-1)*nlevdecomp) = &
                    cf_soil%matrix_Cinput%V(c,j+(i-1)*nlevdecomp) + cf_veg%gap_mortality_c_to_litr_c_col(c,j,i) * dt
               end do
               ! Currently i_cwd .ne. i_litr_max + 1 if .not. fates and
               !           i_cwd = 0 if fates, so not including in the i-loop
               cf_soil%matrix_Cinput%V(c,j+(i_cwd-1)*nlevdecomp) =     &
                 cf_soil%matrix_Cinput%V(c,j+(i_cwd-1)*nlevdecomp)     + cf_veg%gap_mortality_c_to_cwdc_col(c,j) * dt
            end if !soil_matrix
         end do
      end do

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         cs_veg%gresp_storage_patch(p) = cs_veg%gresp_storage_patch(p)           &
              - cf_veg%m_gresp_storage_to_litter_patch(p) * dt
         cs_veg%gresp_xfer_patch(p) = cs_veg%gresp_xfer_patch(p)                 &
              - cf_veg%m_gresp_xfer_to_litter_patch(p) * dt
         !
         ! State update without the matrix solution
         !
         if(.not.  use_matrixcn)then
            ! patch-level carbon fluxes from gap-phase mortality
            ! displayed pools
            cs_veg%leafc_patch(p) = cs_veg%leafc_patch(p)                           &
              - cf_veg%m_leafc_to_litter_patch(p) * dt
            cs_veg%frootc_patch(p) = cs_veg%frootc_patch(p)                         &
              - cf_veg%m_frootc_to_litter_patch(p) * dt
            cs_veg%livestemc_patch(p) = cs_veg%livestemc_patch(p)                   &
              - cf_veg%m_livestemc_to_litter_patch(p) * dt
            cs_veg%deadstemc_patch(p) = cs_veg%deadstemc_patch(p)                   &
              - cf_veg%m_deadstemc_to_litter_patch(p) * dt
            cs_veg%livecrootc_patch(p) = cs_veg%livecrootc_patch(p)                 &
              - cf_veg%m_livecrootc_to_litter_patch(p) * dt
            cs_veg%deadcrootc_patch(p) = cs_veg%deadcrootc_patch(p)                 &
              - cf_veg%m_deadcrootc_to_litter_patch(p) * dt

            ! storage pools
            cs_veg%leafc_storage_patch(p) = cs_veg%leafc_storage_patch(p)           &
              - cf_veg%m_leafc_storage_to_litter_patch(p) * dt
            cs_veg%frootc_storage_patch(p) = cs_veg%frootc_storage_patch(p)         &
              - cf_veg%m_frootc_storage_to_litter_patch(p) * dt
            cs_veg%livestemc_storage_patch(p) = cs_veg%livestemc_storage_patch(p)   &
              - cf_veg%m_livestemc_storage_to_litter_patch(p) * dt
            cs_veg%deadstemc_storage_patch(p) = cs_veg%deadstemc_storage_patch(p)   &
              - cf_veg%m_deadstemc_storage_to_litter_patch(p) * dt
            cs_veg%livecrootc_storage_patch(p) = cs_veg%livecrootc_storage_patch(p) &
              - cf_veg%m_livecrootc_storage_to_litter_patch(p) * dt
            cs_veg%deadcrootc_storage_patch(p) = cs_veg%deadcrootc_storage_patch(p) &
              - cf_veg%m_deadcrootc_storage_to_litter_patch(p) * dt

            ! transfer pools
            cs_veg%leafc_xfer_patch(p) = cs_veg%leafc_xfer_patch(p)                 &
              - cf_veg%m_leafc_xfer_to_litter_patch(p) * dt
            cs_veg%frootc_xfer_patch(p) = cs_veg%frootc_xfer_patch(p)               &
              - cf_veg%m_frootc_xfer_to_litter_patch(p) * dt
            cs_veg%livestemc_xfer_patch(p) = cs_veg%livestemc_xfer_patch(p)         &
              - cf_veg%m_livestemc_xfer_to_litter_patch(p) * dt
            cs_veg%deadstemc_xfer_patch(p) = cs_veg%deadstemc_xfer_patch(p)         &
              - cf_veg%m_deadstemc_xfer_to_litter_patch(p) * dt
            cs_veg%livecrootc_xfer_patch(p) = cs_veg%livecrootc_xfer_patch(p)       &
              - cf_veg%m_livecrootc_xfer_to_litter_patch(p) * dt
            cs_veg%deadcrootc_xfer_patch(p) = cs_veg%deadcrootc_xfer_patch(p)       &
              - cf_veg%m_deadcrootc_xfer_to_litter_patch(p) * dt
         !
         ! For the matrix solution the actual state update comes after the matrix
         ! multiply in VegMatrix, but the matrix needs to be setup with
         ! the equivalent of above. Those changes can be here or in the
         ! native subroutines dealing with that field
         !
         else
            ! NOTE: The matrix version of this is in CNGapMortality (EBK 11/25/2019)
         end if !not use_matrixcn
      end do ! end of patch loop

    end associate

  end subroutine CStateUpdate2

  !-----------------------------------------------------------------------
  subroutine CStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst, soilbiogeochem_carbonstate_inst, &
       soilbiogeochem_carbonflux_inst)
    !
    ! !DESCRIPTION:
    ! Update all the prognostic carbon state
    ! variables affected by harvest mortality fluxes
    !
    ! !ARGUMENTS:
    integer                                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_carbonflux_type)            , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)           , intent(inout) :: cnveg_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)   , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)  , intent(inout) :: soilbiogeochem_carbonstate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,k,l,i  ! indices
    integer :: fp,fc     ! lake filter indices
    real(r8):: dt        ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                                     & 
         cf_veg => cnveg_carbonflux_inst         , &
         cs_veg => cnveg_carbonstate_inst        , &
         cf_soil => soilbiogeochem_carbonflux_inst, & 
         cs_soil => soilbiogeochem_carbonstate_inst &
         )
     
      ! set time steps
      dt = get_step_size_real()

      ! column level carbon fluxes from harvest mortality
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)

            ! column harvest fluxes

            !
            ! State update without the matrix solution
            !
            if (.not. use_soil_matrixcn)then
               do i = i_litr_min, i_litr_max
                  cs_soil%decomp_cpools_vr_col(c,j,i) = &
                    cs_soil%decomp_cpools_vr_col(c,j,i) + &
                    cf_veg%harvest_c_to_litr_c_col(c,j,i) * dt
               end do
               ! Currently i_cwd .ne. i_litr_max + 1 if .not. fates and
               !           i_cwd = 0 if fates, so not including in the i-loop
               cs_soil%decomp_cpools_vr_col(c,j,i_cwd) = &
                    cs_soil%decomp_cpools_vr_col(c,j,i_cwd) + cf_veg%harvest_c_to_cwdc_col(c,j)  * dt
            !
            ! For the matrix solution the actual state update comes after the matrix
            ! multiply in SoilMatrix, but the matrix needs to be setup with
            ! the equivalent of above. Those changes can be here or in the
            ! native subroutines dealing with that field
            !
            else
               ! Match above for matrix method
               do i = i_litr_min, i_litr_max
                  cf_soil%matrix_Cinput%V(c,j+(i-1)*nlevdecomp) = &
                    cf_soil%matrix_Cinput%V(c,j+(i-1)*nlevdecomp) + cf_veg%harvest_c_to_litr_c_col(c,j,i) * dt
               end do
               ! Currently i_cwd .ne. i_litr_max + 1 if .not. fates and
               !           i_cwd = 0 if fates, so not including in the i-loop
               cf_soil%matrix_Cinput%V(c,j+(i_cwd-1)*nlevdecomp) = &
                 cf_soil%matrix_Cinput%V(c,j+(i_cwd-1)*nlevdecomp) + cf_veg%harvest_c_to_cwdc_col(c,j) * dt
            end if

            ! wood to product pools - states updated in CNProducts
         end do
      end do

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! xsmrpool
         cs_veg%xsmrpool_patch(p) = cs_veg%xsmrpool_patch(p)                     &
              - cf_veg%hrv_xsmrpool_to_atm_patch(p) * dt

         ! patch-level carbon fluxes from harvest mortality
         ! storage pools
         cs_veg%gresp_storage_patch(p) = cs_veg%gresp_storage_patch(p)           &
              - cf_veg%hrv_gresp_storage_to_litter_patch(p) * dt

         ! transfer pools
         cs_veg%gresp_xfer_patch(p) = cs_veg%gresp_xfer_patch(p)                 &
              - cf_veg%hrv_gresp_xfer_to_litter_patch(p) * dt


         !
         ! State update without the matrix solution
         !
         if(.not. use_matrixcn)then
            ! displayed pools
            cs_veg%leafc_patch(p) = cs_veg%leafc_patch(p)                           &
              - cf_veg%hrv_leafc_to_litter_patch(p) * dt
            cs_veg%frootc_patch(p) = cs_veg%frootc_patch(p)                         &
              - cf_veg%hrv_frootc_to_litter_patch(p) * dt
            cs_veg%livestemc_patch(p) = cs_veg%livestemc_patch(p)                   &
              - cf_veg%hrv_livestemc_to_litter_patch(p) * dt
            cs_veg%deadstemc_patch(p) = cs_veg%deadstemc_patch(p)                   &
              - cf_veg%wood_harvestc_patch(p) * dt
            cs_veg%livecrootc_patch(p) = cs_veg%livecrootc_patch(p)                 &
              - cf_veg%hrv_livecrootc_to_litter_patch(p) * dt
            cs_veg%deadcrootc_patch(p) = cs_veg%deadcrootc_patch(p)                 &
              - cf_veg%hrv_deadcrootc_to_litter_patch(p) * dt

            ! storage pools
            cs_veg%leafc_storage_patch(p) = cs_veg%leafc_storage_patch(p)           &
              - cf_veg%hrv_leafc_storage_to_litter_patch(p) * dt
            cs_veg%frootc_storage_patch(p) = cs_veg%frootc_storage_patch(p)         &
              - cf_veg%hrv_frootc_storage_to_litter_patch(p) * dt
            cs_veg%livestemc_storage_patch(p) = cs_veg%livestemc_storage_patch(p)   &
              - cf_veg%hrv_livestemc_storage_to_litter_patch(p) * dt
            cs_veg%deadstemc_storage_patch(p) = cs_veg%deadstemc_storage_patch(p)   &
              - cf_veg%hrv_deadstemc_storage_to_litter_patch(p) * dt
            cs_veg%livecrootc_storage_patch(p) = cs_veg%livecrootc_storage_patch(p) &
              - cf_veg%hrv_livecrootc_storage_to_litter_patch(p) * dt
            cs_veg%deadcrootc_storage_patch(p) = cs_veg%deadcrootc_storage_patch(p) &
              - cf_veg%hrv_deadcrootc_storage_to_litter_patch(p) * dt

            ! transfer pools
            cs_veg%leafc_xfer_patch(p) = cs_veg%leafc_xfer_patch(p)                 &
              - cf_veg%hrv_leafc_xfer_to_litter_patch(p) * dt
            cs_veg%frootc_xfer_patch(p) = cs_veg%frootc_xfer_patch(p)               &
              - cf_veg%hrv_frootc_xfer_to_litter_patch(p) * dt
            cs_veg%livestemc_xfer_patch(p) = cs_veg%livestemc_xfer_patch(p)         &
              - cf_veg%hrv_livestemc_xfer_to_litter_patch(p) * dt
            cs_veg%deadstemc_xfer_patch(p) = cs_veg%deadstemc_xfer_patch(p)         &
              - cf_veg%hrv_deadstemc_xfer_to_litter_patch(p) * dt
            cs_veg%livecrootc_xfer_patch(p) = cs_veg%livecrootc_xfer_patch(p)       &
              - cf_veg%hrv_livecrootc_xfer_to_litter_patch(p) * dt
            cs_veg%deadcrootc_xfer_patch(p) = cs_veg%deadcrootc_xfer_patch(p)       &
              - cf_veg%hrv_deadcrootc_xfer_to_litter_patch(p) * dt
         !
         ! For the matrix solution the actual state update comes after the matrix
         ! multiply in VegMatrix, but the matrix needs to be setup with
         ! the equivalent of above. Those changes can be here or in the
         ! native subroutines dealing with that field
         !
         else
            ! NOTE: The matrix equivalent of the above is in CNHarvest (EBK 11/25/2019)
         end if

      end do ! end of patch loop

    end associate

  end subroutine CStateUpdate2h

  !-----------------------------------------------------------------------
  subroutine CStateUpdate2g(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst, &
       soilbiogeochem_carbonstate_inst, soilbiogeochem_carbonflux_inst)
    !
    ! !DESCRIPTION:
    ! Update all the prognostic carbon state
    ! variables affected by gross unrepresented landcover change mortality fluxes
    !
    ! !ARGUMENTS:
    integer                                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_carbonflux_type)            , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)           , intent(inout) :: cnveg_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)   , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)  , intent(inout) :: soilbiogeochem_carbonstate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,k,l,i  ! indices
    integer :: fp,fc     ! lake filter indices
    real(r8):: dt        ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                                     & 
         cf_veg => cnveg_carbonflux_inst         , &
         cs_veg => cnveg_carbonstate_inst        , &
         cf_soil => soilbiogeochem_carbonflux_inst, &
         cs_soil => soilbiogeochem_carbonstate_inst &
         )
     
      ! set time steps
      dt = get_step_size_real()

      ! column level carbon fluxes from gross unrepresented landcover change mortality
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)

            if (.not. use_soil_matrixcn)then
               ! column gross unrepresented landcover change fluxes
               do i = i_litr_min, i_litr_max
                  cs_soil%decomp_cpools_vr_col(c,j,i) = &
                     cs_soil%decomp_cpools_vr_col(c,j,i) + cf_veg%gru_c_to_litr_c_col(c,j,i) * dt
               end do
               cs_soil%decomp_cpools_vr_col(c,j,i_cwd) = &
                    cs_soil%decomp_cpools_vr_col(c,j,i_cwd) + cf_veg%gru_c_to_cwdc_col(c,j) * dt

               ! wood to product pools - states updated in CNProducts
            else
              ! Match above for soil-matrix
               do i = i_litr_min, i_litr_max
                  cf_soil%matrix_Cinput%V(c,j+(i-1)*nlevdecomp) = &
                    cf_soil%matrix_Cinput%V(c,j+(i-1)*nlevdecomp) + cf_veg%gru_c_to_litr_c_col(c,j,i) * dt
               end do
               ! Currently i_cwd .ne. i_litr_max + 1 if .not. fates and
               !           i_cwd = 0 if fates, so not including in the i-loop
               cf_soil%matrix_Cinput%V(c,j+(i_cwd-1)*nlevdecomp) =     &
                 cf_soil%matrix_Cinput%V(c,j+(i_cwd-1)*nlevdecomp) + cf_veg%gru_c_to_cwdc_col(c,j) * dt
            end if !soil_matrix
         end do
      end do

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! xsmrpool
         cs_veg%xsmrpool_patch(p) = cs_veg%xsmrpool_patch(p) &
              - cf_veg%gru_xsmrpool_to_atm_patch(p) * dt
         ! gresp storage pool
         cs_veg%gresp_storage_patch(p) = cs_veg%gresp_storage_patch(p) &
              - cf_veg%gru_gresp_storage_to_atm_patch(p) * dt
         ! gresp transfer pool
         cs_veg%gresp_xfer_patch(p) = cs_veg%gresp_xfer_patch(p) &
              - cf_veg%gru_gresp_xfer_to_atm_patch(p) * dt

         ! patch-level carbon fluxes from gross unrepresented landcover change mortality
         !
         ! State update without the matrix solution
         !
         if(.not. use_matrixcn)then
            ! displayed pools
            cs_veg%leafc_patch(p) = cs_veg%leafc_patch(p)                           &
                 - cf_veg%gru_leafc_to_litter_patch(p) * dt
            cs_veg%frootc_patch(p) = cs_veg%frootc_patch(p)                         &
                 - cf_veg%gru_frootc_to_litter_patch(p) * dt
            cs_veg%livestemc_patch(p) = cs_veg%livestemc_patch(p)                   &
                 - cf_veg%gru_livestemc_to_atm_patch(p) * dt
            cs_veg%deadstemc_patch(p) = cs_veg%deadstemc_patch(p)                   &
                 - cf_veg%gru_deadstemc_to_atm_patch(p) * dt
            cs_veg%deadstemc_patch(p) = cs_veg%deadstemc_patch(p)                   &
                 - cf_veg%gru_wood_productc_gain_patch(p) * dt
            cs_veg%livecrootc_patch(p) = cs_veg%livecrootc_patch(p)                 &
                 - cf_veg%gru_livecrootc_to_litter_patch(p) * dt
            cs_veg%deadcrootc_patch(p) = cs_veg%deadcrootc_patch(p)                 &
                 - cf_veg%gru_deadcrootc_to_litter_patch(p) * dt

            ! storage pools
            cs_veg%leafc_storage_patch(p) = cs_veg%leafc_storage_patch(p)           &
                 - cf_veg%gru_leafc_storage_to_atm_patch(p) * dt
            cs_veg%frootc_storage_patch(p) = cs_veg%frootc_storage_patch(p)         &
                 - cf_veg%gru_frootc_storage_to_atm_patch(p) * dt
            cs_veg%livestemc_storage_patch(p) = cs_veg%livestemc_storage_patch(p)   &
                 - cf_veg%gru_livestemc_storage_to_atm_patch(p) * dt
            cs_veg%deadstemc_storage_patch(p) = cs_veg%deadstemc_storage_patch(p)   &
                 - cf_veg%gru_deadstemc_storage_to_atm_patch(p) * dt
            cs_veg%livecrootc_storage_patch(p) = cs_veg%livecrootc_storage_patch(p) &
                 - cf_veg%gru_livecrootc_storage_to_atm_patch(p) * dt
            cs_veg%deadcrootc_storage_patch(p) = cs_veg%deadcrootc_storage_patch(p) &
                 - cf_veg%gru_deadcrootc_storage_to_atm_patch(p) * dt

            ! transfer pools
            cs_veg%leafc_xfer_patch(p) = cs_veg%leafc_xfer_patch(p)                 &
                 - cf_veg%gru_leafc_xfer_to_atm_patch(p) * dt
            cs_veg%frootc_xfer_patch(p) = cs_veg%frootc_xfer_patch(p)               &
                 - cf_veg%gru_frootc_xfer_to_atm_patch(p) * dt
            cs_veg%livestemc_xfer_patch(p) = cs_veg%livestemc_xfer_patch(p)         &
                 - cf_veg%gru_livestemc_xfer_to_atm_patch(p) * dt
            cs_veg%deadstemc_xfer_patch(p) = cs_veg%deadstemc_xfer_patch(p)         &
                 - cf_veg%gru_deadstemc_xfer_to_atm_patch(p) * dt
            cs_veg%livecrootc_xfer_patch(p) = cs_veg%livecrootc_xfer_patch(p)       &
                 - cf_veg%gru_livecrootc_xfer_to_atm_patch(p) * dt
            cs_veg%deadcrootc_xfer_patch(p) = cs_veg%deadcrootc_xfer_patch(p)       &
                 - cf_veg%gru_deadcrootc_xfer_to_atm_patch(p) * dt

         else
            ! NB (slevis) The matrix equivalent of the above is in
            ! dynGrossUnrepMod::CNGrossUnrep*
         end if
      end do ! end of patch loop

    end associate

  end subroutine CStateUpdate2g

end module CNCStateUpdate2Mod
