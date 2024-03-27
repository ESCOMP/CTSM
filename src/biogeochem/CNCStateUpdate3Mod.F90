module CNCStateUpdate3Mod

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
  use clm_varpar                     , only : nlevdecomp, ndecomp_pools, i_cwd, i_litr_min, i_litr_max
  use CNVegCarbonStateType           , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType            , only : cnveg_carbonflux_type
  use SoilBiogeochemCarbonStateType  , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemCarbonFluxType   , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemDecompCascadeConType , only : use_soil_matrixcn
  use CNSharedParamsMod              , only : use_matrixcn
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CStateUpdate3
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst, soilbiogeochem_carbonstate_inst,&
       soilbiogeochem_carbonflux_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic carbon state
    ! variables affected by fire fluxes
    !
    ! !ARGUMENTS:
    integer                                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_carbonflux_type)            , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)           , intent(inout) :: cnveg_carbonstate_inst
    type(soilbiogeochem_carbonstate_type)  , intent(inout) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)   , intent(inout) :: soilbiogeochem_carbonflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l,k,i  ! indices
    integer :: fp,fc     ! lake filter indices
    real(r8):: dt        ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                                      & 
         cf_veg => cnveg_carbonflux_inst ,          & ! Input
         cs_veg => cnveg_carbonstate_inst,          & ! Output
         cf_soil => soilbiogeochem_carbonflux_inst, & ! Output
         cs_soil => soilbiogeochem_carbonstate_inst & ! Output
         )

      ! set time steps
      dt = get_step_size_real()

      ! column level carbon fluxes from fire
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            ! patch-level wood to column-level CWD (uncombusted wood)

            !
            ! State update without the matrix solution
            !
            if (.not. use_soil_matrixcn) then
               cs_soil%decomp_cpools_vr_col(c,j,i_cwd) = cs_soil%decomp_cpools_vr_col(c,j,i_cwd) + &
                 cf_veg%fire_mortality_c_to_cwdc_col(c,j) * dt

               ! patch-level wood to column-level litter (uncombusted wood)
               do i = i_litr_min, i_litr_max
                  cs_soil%decomp_cpools_vr_col(c,j,i) = &
                     cs_soil%decomp_cpools_vr_col(c,j,i) + &
                     cf_veg%m_c_to_litr_fire_col(c,j,i) * dt
               end do
            !
            ! For the matrix solution the actual state update comes after the matrix
            ! multiply in SoilMatrix, but the matrix needs to be setup with
            ! the equivalent of above. Those changes can be here or in the
            ! native subroutines dealing with that field
            !
            else
               ! Match above for matrix terms
               ! patch-level wood to column-level CWD (uncombusted wood)
               cf_soil%matrix_Cinput%V(c,j+(i_cwd-1)*nlevdecomp) = cf_soil%matrix_Cinput%V(c,j+(i_cwd-1)*nlevdecomp) + &
                 cf_veg%fire_mortality_c_to_cwdc_col(c,j) * dt
               ! patch-level wood to column-level litter (uncombusted wood)
               do i = i_litr_min, i_litr_max
                  cf_soil%matrix_Cinput%V(c,j+(i-1)*nlevdecomp) = cf_soil%matrix_Cinput%V(c,j+(i-1)*nlevdecomp) + &
                    cf_veg%m_c_to_litr_fire_col(c,j,i)* dt
               end do
            end if
         end do
      end do

      ! litter and CWD losses to fire

      !
      ! State update without the matrix solution
      !
      if(.not. use_soil_matrixcn)then
         do l = 1, ndecomp_pools
            do j = 1, nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  cs_soil%decomp_cpools_vr_col(c,j,l) = cs_soil%decomp_cpools_vr_col(c,j,l) - &
                    cf_veg%m_decomp_cpools_to_fire_vr_col(c,j,l) * dt
               end do
            end do
         end do
      end if

      ! patch-level carbon fluxes from fire
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         cs_veg%gresp_storage_patch(p) = cs_veg%gresp_storage_patch(p) -           &
              cf_veg%m_gresp_storage_to_fire_patch(p) * dt
         cs_veg%gresp_storage_patch(p) = cs_veg%gresp_storage_patch(p) -           &
              cf_veg%m_gresp_storage_to_litter_fire_patch(p) * dt
         cs_veg%gresp_xfer_patch(p) = cs_veg%gresp_xfer_patch(p) -                 &
              cf_veg%m_gresp_xfer_to_fire_patch(p) * dt
         cs_veg%gresp_xfer_patch(p) = cs_veg%gresp_xfer_patch(p) -                 &
              cf_veg%m_gresp_xfer_to_litter_fire_patch(p) * dt  
         !
         ! State update without the matrix solution
         !
         if(.not. use_matrixcn)then 
            ! displayed pools
            cs_veg%leafc_patch(p) = cs_veg%leafc_patch(p) -                           &
              cf_veg%m_leafc_to_fire_patch(p) * dt
            cs_veg%leafc_patch(p) = cs_veg%leafc_patch(p) -                           &
              cf_veg%m_leafc_to_litter_fire_patch(p) * dt
            cs_veg%frootc_patch(p) = cs_veg%frootc_patch(p) -                         &
              cf_veg%m_frootc_to_fire_patch(p) * dt
            cs_veg%frootc_patch(p) = cs_veg%frootc_patch(p) -                         &
              cf_veg%m_frootc_to_litter_fire_patch(p) * dt
            cs_veg%livestemc_patch(p) = cs_veg%livestemc_patch(p) -                   &
              cf_veg%m_livestemc_to_fire_patch(p) * dt
            cs_veg%livestemc_patch(p) = cs_veg%livestemc_patch(p) -                   &
              cf_veg%m_livestemc_to_litter_fire_patch(p) * dt  -                   &
              cf_veg%m_livestemc_to_deadstemc_fire_patch(p) * dt
            cs_veg%deadstemc_patch(p) = cs_veg%deadstemc_patch(p) -                   &
              cf_veg%m_deadstemc_to_fire_patch(p) * dt
            cs_veg%deadstemc_patch(p) = cs_veg%deadstemc_patch(p) -                   &
              cf_veg%m_deadstemc_to_litter_fire_patch(p) * dt  +                   &
              cf_veg%m_livestemc_to_deadstemc_fire_patch(p) * dt
            cs_veg%livecrootc_patch(p) = cs_veg%livecrootc_patch(p) -                 &
              cf_veg%m_livecrootc_to_fire_patch(p) * dt
            cs_veg%livecrootc_patch(p) = cs_veg%livecrootc_patch(p) -                 &
              cf_veg%m_livecrootc_to_litter_fire_patch(p) * dt   -                 &
              cf_veg%m_livecrootc_to_deadcrootc_fire_patch(p) * dt
            cs_veg%deadcrootc_patch(p) = cs_veg%deadcrootc_patch(p) -                 &
              cf_veg%m_deadcrootc_to_fire_patch(p) * dt
            cs_veg%deadcrootc_patch(p) = cs_veg%deadcrootc_patch(p) -                 &
              cf_veg%m_deadcrootc_to_litter_fire_patch(p)* dt    +                 &
              cf_veg%m_livecrootc_to_deadcrootc_fire_patch(p) * dt

            ! storage pools
            cs_veg%leafc_storage_patch(p) = cs_veg%leafc_storage_patch(p) -           &
              cf_veg%m_leafc_storage_to_fire_patch(p) * dt
            cs_veg%leafc_storage_patch(p) = cs_veg%leafc_storage_patch(p) -           &
              cf_veg%m_leafc_storage_to_litter_fire_patch(p) * dt
            cs_veg%frootc_storage_patch(p) = cs_veg%frootc_storage_patch(p) -         &
              cf_veg%m_frootc_storage_to_fire_patch(p) * dt
            cs_veg%frootc_storage_patch(p) = cs_veg%frootc_storage_patch(p) -         &
              cf_veg%m_frootc_storage_to_litter_fire_patch(p) * dt
            cs_veg%livestemc_storage_patch(p) = cs_veg%livestemc_storage_patch(p) -   &
              cf_veg%m_livestemc_storage_to_fire_patch(p) * dt
            cs_veg%livestemc_storage_patch(p) = cs_veg%livestemc_storage_patch(p) -   &
              cf_veg%m_livestemc_storage_to_litter_fire_patch(p) * dt
            cs_veg%deadstemc_storage_patch(p) = cs_veg%deadstemc_storage_patch(p) -   &
              cf_veg%m_deadstemc_storage_to_fire_patch(p) * dt
            cs_veg%deadstemc_storage_patch(p) = cs_veg%deadstemc_storage_patch(p) -   &
              cf_veg%m_deadstemc_storage_to_litter_fire_patch(p) * dt
            cs_veg%livecrootc_storage_patch(p) = cs_veg%livecrootc_storage_patch(p) - &
              cf_veg%m_livecrootc_storage_to_fire_patch(p) * dt
            cs_veg%livecrootc_storage_patch(p) = cs_veg%livecrootc_storage_patch(p) - &
              cf_veg%m_livecrootc_storage_to_litter_fire_patch(p)* dt
            cs_veg%deadcrootc_storage_patch(p) = cs_veg%deadcrootc_storage_patch(p) - &
              cf_veg%m_deadcrootc_storage_to_fire_patch(p) * dt
            cs_veg%deadcrootc_storage_patch(p) = cs_veg%deadcrootc_storage_patch(p) - &
              cf_veg%m_deadcrootc_storage_to_litter_fire_patch(p)* dt

            ! transfer pools
            cs_veg%leafc_xfer_patch(p) = cs_veg%leafc_xfer_patch(p) -                 &
              cf_veg%m_leafc_xfer_to_fire_patch(p) * dt
            cs_veg%leafc_xfer_patch(p) = cs_veg%leafc_xfer_patch(p) -                 &
              cf_veg%m_leafc_xfer_to_litter_fire_patch(p) * dt
            cs_veg%frootc_xfer_patch(p) = cs_veg%frootc_xfer_patch(p) -               &
              cf_veg%m_frootc_xfer_to_fire_patch(p) * dt
            cs_veg%frootc_xfer_patch(p) = cs_veg%frootc_xfer_patch(p) -               &
              cf_veg%m_frootc_xfer_to_litter_fire_patch(p) * dt
            cs_veg%livestemc_xfer_patch(p) = cs_veg%livestemc_xfer_patch(p) -         &
              cf_veg%m_livestemc_xfer_to_fire_patch(p) * dt
            cs_veg%livestemc_xfer_patch(p) = cs_veg%livestemc_xfer_patch(p) -         &
              cf_veg%m_livestemc_xfer_to_litter_fire_patch(p) * dt
            cs_veg%deadstemc_xfer_patch(p) = cs_veg%deadstemc_xfer_patch(p) -         &
              cf_veg%m_deadstemc_xfer_to_fire_patch(p) * dt
            cs_veg%deadstemc_xfer_patch(p) = cs_veg%deadstemc_xfer_patch(p) -         &
              cf_veg%m_deadstemc_xfer_to_litter_fire_patch(p) * dt
            cs_veg%livecrootc_xfer_patch(p) = cs_veg%livecrootc_xfer_patch(p) -       &
              cf_veg%m_livecrootc_xfer_to_fire_patch(p) * dt
            cs_veg%livecrootc_xfer_patch(p) = cs_veg%livecrootc_xfer_patch(p) -       &
              cf_veg%m_livecrootc_xfer_to_litter_fire_patch(p)* dt
            cs_veg%deadcrootc_xfer_patch(p) = cs_veg%deadcrootc_xfer_patch(p) -       &
              cf_veg%m_deadcrootc_xfer_to_fire_patch(p) * dt
            cs_veg%deadcrootc_xfer_patch(p) = cs_veg%deadcrootc_xfer_patch(p) -       &
              cf_veg%m_deadcrootc_xfer_to_litter_fire_patch(p)* dt
         !
         ! For the matrix solution the actual state update comes after the matrix
         ! multiply in VegMatrix, but the matrix needs to be setup with
         ! the equivalent of above. Those changes can be here or in the
         ! native subroutines dealing with that field
         !
         else
            ! NOTE: The equivalent changes for matrix code are in CNFireBase and CNFireLi2014 codes EBK (11/26/2019)
         end if !not use_matrixcn
      end do ! end of patch loop

    end associate

  end subroutine CStateUpdate3

end module CNCStateUpdate3Mod
