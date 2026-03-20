module CNNStateUpdate2Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for nitrogen state variable update, mortality fluxes.
  ! When the matrix solution is being used (use_matrixcn and use_soil_matrixcn)
  ! only some state updates are done here, the other state updates happen
  ! after the matrix is solved in VegMatrix and SoilMatrix.
  !
  ! !USES:
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use clm_time_manager                , only : get_step_size_real
  use clm_varpar                      , only : nlevsoi, nlevdecomp
  use clm_varpar                      , only : i_litr_min, i_litr_max, i_cwd
  use clm_varctl                      , only : iulog
  use SoilBiogeochemDecompCascadeConType, only : use_soil_matrixcn
  use CNSharedParamsMod               , only : use_matrixcn
  use CNVegNitrogenStateType          , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType           , only : cnveg_nitrogenflux_type
  use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: NStateUpdate2
  public:: NStateUpdate2h
  public:: NStateUpdate2g
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine NStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, soilbiogeochem_nitrogenstate_inst, &
       soilbiogeochem_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic nitrogen state
    ! variables affected by gap-phase mortality fluxes
    ! NOTE - associate statements have been removed where there are
    ! no science equations. This increases readability and maintainability
    !
    ! !ARGUMENTS:
    integer                                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                                 , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_nitrogenflux_type)           , intent(in)    :: cnveg_nitrogenflux_inst
    type(cnveg_nitrogenstate_type)          , intent(inout) :: cnveg_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,l,i  ! indices
    integer  :: fp,fc   ! lake filter indices
    real(r8) :: dt      ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                                        & 
         nf_veg => cnveg_nitrogenflux_inst          , &
         ns_veg => cnveg_nitrogenstate_inst         , &
         nf_soil => soilbiogeochem_nitrogenflux_inst, &
         ns_soil => soilbiogeochem_nitrogenstate_inst  &
         )

      ! set time steps
      dt = get_step_size_real()

      ! column-level nitrogen fluxes from gap-phase mortality

      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)

            !
            ! State update without the matrix solution
            !
            if (.not. use_soil_matrixcn)then 
               do i = i_litr_min, i_litr_max
                  ns_soil%decomp_npools_vr_col(c,j,i) = &
                     ns_soil%decomp_npools_vr_col(c,j,i) + &
                     nf_veg%gap_mortality_n_to_litr_n_col(c,j,i) * dt
               end do
               ! Currently i_cwd .ne. i_litr_max + 1 if .not. fates and
               !           i_cwd = 0 if fates, so not including in the i-loop
               ns_soil%decomp_npools_vr_col(c,j,i_cwd)     = &
                 ns_soil%decomp_npools_vr_col(c,j,i_cwd)     + nf_veg%gap_mortality_n_to_cwdn_col(c,j)       * dt
            !
            ! For the matrix solution the actual state update comes after the matrix
            ! multiply in SoilMatrix, but the matrix needs to be setup with
            ! the equivalent of above. Those changes can be here or in the
            ! native subroutines dealing with that field
            !
            else
               ! Do above for the matrix solution
               do i = i_litr_min, i_litr_max
                  nf_soil%matrix_Ninput%V(c,j+(i-1)*nlevdecomp) = &
                    nf_soil%matrix_Ninput%V(c,j+(i-1)*nlevdecomp) + nf_veg%gap_mortality_n_to_litr_n_col(c,j,i) * dt
               end do
               ! Currently i_cwd .ne. i_litr_max + 1 if .not. fates and
               !           i_cwd = 0 if fates, so not including in the i-loop
               nf_soil%matrix_Ninput%V(c,j+(i_cwd-1)*nlevdecomp)     = &
                 nf_soil%matrix_Ninput%V(c,j+(i_cwd-1)*nlevdecomp)     + nf_veg%gap_mortality_n_to_cwdn_col(c,j)       * dt
            end if !not use_soil_matrix
         end do
      end do

      ! patch -level nitrogen fluxes from gap-phase mortality

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         !
         ! State update without the matrix solution
         !
         if(.not.  use_matrixcn)then
            ! displayed pools
            ns_veg%leafn_patch(p) =  ns_veg%leafn_patch(p)                           &
              - nf_veg%m_leafn_to_litter_patch(p) * dt
            ns_veg%frootn_patch(p) =  ns_veg%frootn_patch(p)                         &
              - nf_veg%m_frootn_to_litter_patch(p) * dt
            ns_veg%livestemn_patch(p) =  ns_veg%livestemn_patch(p)                   &
              - nf_veg%m_livestemn_to_litter_patch(p) * dt
            ns_veg%deadstemn_patch(p) =  ns_veg%deadstemn_patch(p)                   &
              - nf_veg%m_deadstemn_to_litter_patch(p) * dt
            ns_veg%livecrootn_patch(p) =  ns_veg%livecrootn_patch(p)                 &
              - nf_veg%m_livecrootn_to_litter_patch(p) * dt
            ns_veg%deadcrootn_patch(p) =  ns_veg%deadcrootn_patch(p)                 &
              - nf_veg%m_deadcrootn_to_litter_patch(p) * dt
            ns_veg%retransn_patch(p) =  ns_veg%retransn_patch(p)                     &
              - nf_veg%m_retransn_to_litter_patch(p) * dt

            ! storage pools
            ns_veg%leafn_storage_patch(p) =  ns_veg%leafn_storage_patch(p)           &
              - nf_veg%m_leafn_storage_to_litter_patch(p) * dt
            ns_veg%frootn_storage_patch(p) =  ns_veg%frootn_storage_patch(p)         &
              - nf_veg%m_frootn_storage_to_litter_patch(p) * dt
            ns_veg%livestemn_storage_patch(p) =  ns_veg%livestemn_storage_patch(p)   &
              - nf_veg%m_livestemn_storage_to_litter_patch(p) * dt
            ns_veg%deadstemn_storage_patch(p) =  ns_veg%deadstemn_storage_patch(p)   &
              - nf_veg%m_deadstemn_storage_to_litter_patch(p) * dt
            ns_veg%livecrootn_storage_patch(p) =  ns_veg%livecrootn_storage_patch(p) &
              - nf_veg%m_livecrootn_storage_to_litter_patch(p) * dt
            ns_veg%deadcrootn_storage_patch(p) =  ns_veg%deadcrootn_storage_patch(p) &
              - nf_veg%m_deadcrootn_storage_to_litter_patch(p) * dt

            ! transfer pools
            ns_veg%leafn_xfer_patch(p) =  ns_veg%leafn_xfer_patch(p)                 &
              - nf_veg%m_leafn_xfer_to_litter_patch(p) * dt
            ns_veg%frootn_xfer_patch(p) =  ns_veg%frootn_xfer_patch(p)               &
              - nf_veg%m_frootn_xfer_to_litter_patch(p) * dt
            ns_veg%livestemn_xfer_patch(p) =  ns_veg%livestemn_xfer_patch(p)         &
              - nf_veg%m_livestemn_xfer_to_litter_patch(p) * dt
            ns_veg%deadstemn_xfer_patch(p) =  ns_veg%deadstemn_xfer_patch(p)         &
              - nf_veg%m_deadstemn_xfer_to_litter_patch(p) * dt
            ns_veg%livecrootn_xfer_patch(p) =  ns_veg%livecrootn_xfer_patch(p)       &
              - nf_veg%m_livecrootn_xfer_to_litter_patch(p) * dt
            ns_veg%deadcrootn_xfer_patch(p) =  ns_veg%deadcrootn_xfer_patch(p)       &
              - nf_veg%m_deadcrootn_xfer_to_litter_patch(p) * dt

         !
         ! For the matrix solution the actual state update comes after the matrix
         ! multiply in VegMatrix, but the matrix needs to be setup with
         ! the equivalent of above. Those changes can be here or in the
         ! native subroutines dealing with that field
         !
         else
            ! NOTE: The equivalent changes for matrix code are in dynHarvest::CNHarvest EBK (11/26/2019)
         end if !not use_matrixcn
      end do

    end associate

  end subroutine NStateUpdate2

  !-----------------------------------------------------------------------
  subroutine NStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, soilbiogeochem_nitrogenstate_inst, &
       soilbiogeochem_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! Update all the prognostic nitrogen state
    ! variables affected by harvest mortality fluxes
    ! NOTE - associate statements have been removed where there are
    ! no science equations. This increases readability and maintainability
    !
    ! !ARGUMENTS:
    integer                                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                                 , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_nitrogenflux_type)           , intent(in)    :: cnveg_nitrogenflux_inst
    type(cnveg_nitrogenstate_type)          , intent(inout) :: cnveg_nitrogenstate_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l,i  ! indices
    integer :: fp,fc   ! lake filter indices
    real(r8):: dt      ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                                 & 
         nf_veg  => cnveg_nitrogenflux_inst  , &
         ns_veg  => cnveg_nitrogenstate_inst , &
         nf_soil => soilbiogeochem_nitrogenflux_inst ,  &
         ns_soil => soilbiogeochem_nitrogenstate_inst   &
         )

      ! set time steps
      dt = get_step_size_real()

      ! column-level nitrogen fluxes from harvest mortality

      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            !
            ! State update without the matrix solution
            !
            if (.not. use_soil_matrixcn)then
               do i = i_litr_min, i_litr_max
                  ns_soil%decomp_npools_vr_col(c,j,i) = &
                     ns_soil%decomp_npools_vr_col(c,j,i) + &
                     nf_veg%harvest_n_to_litr_n_col(c,j,i) * dt
               end do
               ! Currently i_cwd .ne. i_litr_max + 1 if .not. fates and
               !           i_cwd = 0 if fates, so not including in the i-loop
               ns_soil%decomp_npools_vr_col(c,j,i_cwd)     = &
                 ns_soil%decomp_npools_vr_col(c,j,i_cwd)     + nf_veg%harvest_n_to_cwdn_col(c,j)       * dt
            !
            ! For the matrix solution the actual state update comes after the matrix
            ! multiply in SoilMatrix, but the matrix needs to be setup with
            ! the equivalent of above. Those changes can be here or in the
            ! native subroutines dealing with that field
            !
            else
               ! Do above for the matrix solution
               do i = i_litr_min, i_litr_max
                  nf_soil%matrix_Ninput%V(c,j+(i-1)*nlevdecomp) = &
                    nf_soil%matrix_Ninput%V(c,j+(i-1)*nlevdecomp) + nf_veg%harvest_n_to_litr_n_col(c,j,i) * dt
               end do
               ! Currently i_cwd .ne. i_litr_max + 1 if .not. fates and
               !           i_cwd = 0 if fates, so not including in the i-loop
               nf_soil%matrix_Ninput%V(c,j+(i_cwd-1)*nlevdecomp) = &
                 nf_soil%matrix_Ninput%V(c,j+(i_cwd-1)*nlevdecomp) + nf_veg%harvest_n_to_cwdn_col(c,j) * dt
            end if !not use_soil_matrixcn
         end do
      end do

      ! patch-level nitrogen fluxes from harvest mortality

      do fp = 1,num_soilp
         p = filter_soilp(fp)


         !
         ! State update without the matrix solution
         !
         if(.not. use_matrixcn)then
            ! displayed pools
            ns_veg%leafn_patch(p) = ns_veg%leafn_patch(p)                           &
              - nf_veg%hrv_leafn_to_litter_patch(p) * dt
            ns_veg%frootn_patch(p) = ns_veg%frootn_patch(p)                         &
              - nf_veg%hrv_frootn_to_litter_patch(p) * dt
            ns_veg%livestemn_patch(p) = ns_veg%livestemn_patch(p)                   &
              - nf_veg%hrv_livestemn_to_litter_patch(p) * dt
            ns_veg%deadstemn_patch(p) = ns_veg%deadstemn_patch(p)                   &
              - nf_veg%wood_harvestn_patch(p) * dt
            ns_veg%livecrootn_patch(p) = ns_veg%livecrootn_patch(p)                 &
              - nf_veg%hrv_livecrootn_to_litter_patch(p) * dt
            ns_veg%deadcrootn_patch(p) = ns_veg%deadcrootn_patch(p)                 &
              - nf_veg%hrv_deadcrootn_to_litter_patch(p) * dt
            ns_veg%retransn_patch(p) = ns_veg%retransn_patch(p)                     &
              - nf_veg%hrv_retransn_to_litter_patch(p) * dt

            ! storage pools
            ns_veg%leafn_storage_patch(p) = ns_veg%leafn_storage_patch(p)           &
              - nf_veg%hrv_leafn_storage_to_litter_patch(p) * dt
            ns_veg%frootn_storage_patch(p) = ns_veg%frootn_storage_patch(p)         &
              - nf_veg%hrv_frootn_storage_to_litter_patch(p) * dt
            ns_veg%livestemn_storage_patch(p) = ns_veg%livestemn_storage_patch(p)   &
              - nf_veg%hrv_livestemn_storage_to_litter_patch(p) * dt
            ns_veg%deadstemn_storage_patch(p) = ns_veg%deadstemn_storage_patch(p)   &
              - nf_veg%hrv_deadstemn_storage_to_litter_patch(p) * dt
            ns_veg%livecrootn_storage_patch(p) = ns_veg%livecrootn_storage_patch(p) &
              - nf_veg%hrv_livecrootn_storage_to_litter_patch(p) * dt
            ns_veg%deadcrootn_storage_patch(p) = ns_veg%deadcrootn_storage_patch(p) &
              - nf_veg%hrv_deadcrootn_storage_to_litter_patch(p) * dt

            ! transfer pools
            ns_veg%leafn_xfer_patch(p) = ns_veg%leafn_xfer_patch(p)                 &
              - nf_veg%hrv_leafn_xfer_to_litter_patch(p) *dt
            ns_veg%frootn_xfer_patch(p) = ns_veg%frootn_xfer_patch(p)               &
              - nf_veg%hrv_frootn_xfer_to_litter_patch(p) *dt
            ns_veg%livestemn_xfer_patch(p) = ns_veg%livestemn_xfer_patch(p)         &
              - nf_veg%hrv_livestemn_xfer_to_litter_patch(p) *dt
            ns_veg%deadstemn_xfer_patch(p) = ns_veg%deadstemn_xfer_patch(p)         &
              - nf_veg%hrv_deadstemn_xfer_to_litter_patch(p) *dt
            ns_veg%livecrootn_xfer_patch(p) = ns_veg%livecrootn_xfer_patch(p)       &
              - nf_veg%hrv_livecrootn_xfer_to_litter_patch(p) *dt
            ns_veg%deadcrootn_xfer_patch(p) = ns_veg%deadcrootn_xfer_patch(p)       &
              - nf_veg%hrv_deadcrootn_xfer_to_litter_patch(p) *dt
         !
         ! For the matrix solution the actual state update comes after the matrix
         ! multiply in VegMatrix, but the matrix needs to be setup with
         ! the equivalent of above. Those changes can be here or in the
         ! native subroutines dealing with that field
         !
         else
            ! NOTE: The equivalent changes for matrix code are in dynHarvest::CNHarvest EBK (11/26/2019)
         end if !not use_matrixcn

      end do

    end associate

  end subroutine NStateUpdate2h

  !-----------------------------------------------------------------------
  subroutine NStateUpdate2g(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, &
       soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! Update all the prognostic nitrogen state
    ! variables affected by gross unrepresented landcover change mortality fluxes
    ! NOTE - associate statements have been removed where there are
    ! no science equations. This increases readability and maintainability
    !
    ! !ARGUMENTS:
    integer                                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                                 , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_nitrogenflux_type)           , intent(in)    :: cnveg_nitrogenflux_inst
    type(cnveg_nitrogenstate_type)          , intent(inout) :: cnveg_nitrogenstate_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l,i  ! indices
    integer :: fp,fc   ! lake filter indices
    real(r8):: dt      ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                                 & 
         nf_veg  => cnveg_nitrogenflux_inst  , &
         ns_veg  => cnveg_nitrogenstate_inst , &
         nf_soil => soilbiogeochem_nitrogenflux_inst, &
         ns_soil => soilbiogeochem_nitrogenstate_inst   &
         )

      ! set time steps
      dt = get_step_size_real()

      ! column-level nitrogen fluxes from gross unrepresented landcover change mortality

      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            !
            ! State update without the matrix solution
            !
            if (.not. use_soil_matrixcn)then
               do i = i_litr_min, i_litr_max
                  ns_soil%decomp_npools_vr_col(c,j,i) = &
                     ns_soil%decomp_npools_vr_col(c,j,i) + nf_veg%gru_n_to_litr_n_col(c,j,i) * dt
               end do
               ! Currently i_cwd .ne. i_litr_max + 1 if .not. fates and
               !           i_cwd = 0 if fates, so not including in the i-loop
               ns_soil%decomp_npools_vr_col(c,j,i_cwd) = &
                    ns_soil%decomp_npools_vr_col(c,j,i_cwd) + nf_veg%gru_n_to_cwdn_col(c,j) * dt
            else
               ! Do above for the matrix solution
               do i = i_litr_min, i_litr_max
                  nf_soil%matrix_Ninput%V(c,j+(i-1)*nlevdecomp) = &
                    nf_soil%matrix_Ninput%V(c,j+(i-1)*nlevdecomp) + nf_veg%gru_n_to_litr_n_col(c,j,i) * dt
               end do
               ! Currently i_cwd .ne. i_litr_max + 1 if .not. fates and
               !           i_cwd = 0 if fates, so not including in the i-loop
               nf_soil%matrix_Ninput%V(c,j+(i_cwd-1)*nlevdecomp) = &
                 nf_soil%matrix_Ninput%V(c,j+(i_cwd-1)*nlevdecomp) + nf_veg%gru_n_to_cwdn_col(c,j) * dt
            end if !not use_soil_matrixcn
         end do
      end do

      ! patch-level nitrogen fluxes from gross unrepresented landcover change mortality

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         !
         ! State update without the matrix solution
         !
         if(.not. use_matrixcn)then
            ! displayed pools
            ns_veg%leafn_patch(p) = ns_veg%leafn_patch(p)                           &
                 - nf_veg%gru_leafn_to_litter_patch(p) * dt
            ns_veg%frootn_patch(p) = ns_veg%frootn_patch(p)                         &
                 - nf_veg%gru_frootn_to_litter_patch(p) * dt
            ns_veg%livestemn_patch(p) = ns_veg%livestemn_patch(p)                   &
                 - nf_veg%gru_livestemn_to_atm_patch(p) * dt
            ns_veg%deadstemn_patch(p) = ns_veg%deadstemn_patch(p)                   &
                 - nf_veg%gru_deadstemn_to_atm_patch(p) * dt
            ns_veg%deadstemn_patch(p) = ns_veg%deadstemn_patch(p)                   &
                 - nf_veg%gru_wood_productn_gain_patch(p) * dt
            ns_veg%livecrootn_patch(p) = ns_veg%livecrootn_patch(p)                 &
                 - nf_veg%gru_livecrootn_to_litter_patch(p) * dt
            ns_veg%deadcrootn_patch(p) = ns_veg%deadcrootn_patch(p)                 &
                 - nf_veg%gru_deadcrootn_to_litter_patch(p) * dt
            ns_veg%retransn_patch(p) = ns_veg%retransn_patch(p)                     &
                 - nf_veg%gru_retransn_to_litter_patch(p) * dt

            ! storage pools
            ns_veg%leafn_storage_patch(p) = ns_veg%leafn_storage_patch(p)           &
                 - nf_veg%gru_leafn_storage_to_atm_patch(p) * dt
            ns_veg%frootn_storage_patch(p) = ns_veg%frootn_storage_patch(p)         &
                 - nf_veg%gru_frootn_storage_to_atm_patch(p) * dt
            ns_veg%livestemn_storage_patch(p) = ns_veg%livestemn_storage_patch(p)   &
                 - nf_veg%gru_livestemn_storage_to_atm_patch(p) * dt
            ns_veg%deadstemn_storage_patch(p) = ns_veg%deadstemn_storage_patch(p)   &
                 - nf_veg%gru_deadstemn_storage_to_atm_patch(p) * dt
            ns_veg%livecrootn_storage_patch(p) = ns_veg%livecrootn_storage_patch(p) &
                 - nf_veg%gru_livecrootn_storage_to_atm_patch(p) * dt
            ns_veg%deadcrootn_storage_patch(p) = ns_veg%deadcrootn_storage_patch(p) &
                 - nf_veg%gru_deadcrootn_storage_to_atm_patch(p) * dt

            ! transfer pools
            ns_veg%leafn_xfer_patch(p) = ns_veg%leafn_xfer_patch(p)                 &
                 - nf_veg%gru_leafn_xfer_to_atm_patch(p) *dt
            ns_veg%frootn_xfer_patch(p) = ns_veg%frootn_xfer_patch(p)               &
                 - nf_veg%gru_frootn_xfer_to_atm_patch(p) *dt
            ns_veg%livestemn_xfer_patch(p) = ns_veg%livestemn_xfer_patch(p)         &
                 - nf_veg%gru_livestemn_xfer_to_atm_patch(p) *dt
            ns_veg%deadstemn_xfer_patch(p) = ns_veg%deadstemn_xfer_patch(p)         &
                 - nf_veg%gru_deadstemn_xfer_to_atm_patch(p) *dt
            ns_veg%livecrootn_xfer_patch(p) = ns_veg%livecrootn_xfer_patch(p)       &
                 - nf_veg%gru_livecrootn_xfer_to_atm_patch(p) *dt
            ns_veg%deadcrootn_xfer_patch(p) = ns_veg%deadcrootn_xfer_patch(p)       &
                 - nf_veg%gru_deadcrootn_xfer_to_atm_patch(p) *dt
         else
            ! NB (slevis) The equivalent changes for matrix code are in
            ! dynGrossUnrepMod::CNGrossUnrep*
         end if !not use_matrixcn

      end do

    end associate

  end subroutine NStateUpdate2g

end module CNNStateUpdate2Mod
