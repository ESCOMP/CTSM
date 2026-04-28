module SoilBiogeochemCompetition_mod

  !-----------------------------------------------------------------------
  ! Standalone extraction of SoilBiogeochemCompetition from
  ! src/soilbiogeochem/SoilBiogeochemCompetitionMod.F90.
  !
  ! Stage 4: shr_kind_mod dependency removed; r8 is now defined locally
  ! via selected_real_kind(12), which matches CTSM's shr_kind_r8 in
  ! practice. The module has no non-intrinsic dependencies.
  !-----------------------------------------------------------------------

  implicit none
  private

  integer, parameter, public :: r8 = selected_real_kind(12)

  public :: SoilBiogeochemCompetition

contains

  !-----------------------------------------------------------------------
   subroutine SoilBiogeochemCompetition( &
       ! sizes / index ranges
       begc, endc, nlevdecomp, ndecomp_cascade_transitions, &
       num_bgc_soilc, filter_bgc_soilc,                     &
       ! scalar config (was module state / runtime flags / params_inst)
       dt, bdnr,                                            &
       use_nitrif_denitrif, carbon_only,                    &
       decomp_method, mimics_decomp, i_cop_mic, i_oli_mic,  &
       compet_plant_no3, compet_plant_nh4,                  &
       compet_decomp_no3, compet_decomp_nh4,                &
       compet_denit, compet_nit,                            &
       ! 1D arrays
       dzsoi_decomp, cascade_receiver_pool, landunit,       &
       ! state-block fields (was soilbiogeochem_state_inst%*_col)
       fpg, fpi, fpi_vr, nfixation_prof, plant_ndemand,     &
       ! n-state fields
       sminn_vr, smin_nh4_vr, smin_no3_vr,                  &
       ! c-flux fields
       c_overflow_vr,                                       &
       ! n-flux fields
       pot_f_nit_vr, pot_f_denit_vr, f_nit_vr, f_denit_vr,  &
       potential_immob, actual_immob, sminn_to_plant,       &
       sminn_to_denit_excess_vr,                            &
       actual_immob_no3_vr, actual_immob_nh4_vr,            &
       smin_no3_to_plant_vr, smin_nh4_to_plant_vr,          &
       n2_n2o_ratio_denit_vr, f_n2o_denit_vr, f_n2o_nit_vr, &
       supplement_to_sminn_vr, sminn_to_plant_vr,           &
       potential_immob_vr, actual_immob_vr,                 &
       ! 3D arrays
       pmnf_decomp_cascade, p_decomp_cn_gain)
    !
    ! !ARGUMENTS:
    integer , intent(in) :: begc, endc                                  ! column index range (was bounds%begc:bounds%endc)
    integer , intent(in) :: nlevdecomp                                  ! number of biogeochemically active soil layers
    integer , intent(in) :: ndecomp_cascade_transitions                 ! number of decomposition cascade transitions
    integer , intent(in) :: num_bgc_soilc                               ! number of soil columns in filter
    integer , intent(in) :: filter_bgc_soilc(:)                         ! filter for soil columns
    real(r8), intent(in) :: dt                                          ! decomp timestep (seconds)
    real(r8), intent(in) :: bdnr                                        ! bulk denitrification rate (1/s)
    logical , intent(in) :: use_nitrif_denitrif                         ! true => use nitrif/denitrif branch
    logical , intent(in) :: carbon_only                                 ! true => carbon-only mode (was allocate_carbon_only())
    integer , intent(in) :: decomp_method                               ! type of decomposition method
    integer , intent(in) :: mimics_decomp                               ! id value of MIMICS decomposition method
    integer , intent(in) :: i_cop_mic                                   ! copiotrophic microbial pool index
    integer , intent(in) :: i_oli_mic                                   ! oligotrophic microbial pool index
    real(r8), intent(in) :: compet_plant_no3                            ! relative competitiveness of plants for NO3
    real(r8), intent(in) :: compet_plant_nh4                            ! relative competitiveness of plants for NH4
    real(r8), intent(in) :: compet_decomp_no3                           ! relative competitiveness of immobilizers for NO3
    real(r8), intent(in) :: compet_decomp_nh4                           ! relative competitiveness of immobilizers for NH4
    real(r8), intent(in) :: compet_denit                                ! relative competitiveness of denitrifiers for NO3
    real(r8), intent(in) :: compet_nit                                  ! relative competitiveness of nitrifiers for NH4
    real(r8), intent(in) :: dzsoi_decomp(:)                             ! per-layer thickness for biogeochemical layers
    integer , pointer    :: cascade_receiver_pool(:)                    ! which pool is C added to for a given decomposition step
    integer , pointer    :: landunit(:)                                 ! landunit index per column (was col%landunit)
    real(r8), pointer    :: fpg(:)                                      ! fraction of potential gpp
    real(r8), pointer    :: fpi(:)                                      ! fraction of potential immobilization
    real(r8), pointer    :: fpi_vr(:,:)                                 ! fraction of potential immobilization (per layer)
    real(r8), pointer    :: nfixation_prof(:,:)
    real(r8), pointer    :: plant_ndemand(:)                            ! column-level plant N demand
    real(r8), pointer    :: sminn_vr(:,:)                               ! (gN/m3) soil mineral N
    real(r8), pointer    :: smin_nh4_vr(:,:)                            ! (gN/m3) soil mineral NH4
    real(r8), pointer    :: smin_no3_vr(:,:)                            ! (gN/m3) soil mineral NO3
    real(r8), pointer    :: c_overflow_vr(:,:,:)                        ! (gC/m3/s) C rejected by microbes that cannot process it
    real(r8), pointer    :: pot_f_nit_vr(:,:)                           ! (gN/m3/s) potential soil nitrification flux
    real(r8), pointer    :: pot_f_denit_vr(:,:)                         ! (gN/m3/s) potential soil denitrification flux
    real(r8), pointer    :: f_nit_vr(:,:)                               ! (gN/m3/s) soil nitrification flux
    real(r8), pointer    :: f_denit_vr(:,:)                             ! (gN/m3/s) soil denitrification flux
    real(r8), pointer    :: potential_immob(:)
    real(r8), pointer    :: actual_immob(:)
    real(r8), pointer    :: sminn_to_plant(:)
    real(r8), pointer    :: sminn_to_denit_excess_vr(:,:)
    real(r8), pointer    :: actual_immob_no3_vr(:,:)
    real(r8), pointer    :: actual_immob_nh4_vr(:,:)
    real(r8), pointer    :: smin_no3_to_plant_vr(:,:)
    real(r8), pointer    :: smin_nh4_to_plant_vr(:,:)
    real(r8), pointer    :: n2_n2o_ratio_denit_vr(:,:)                  ! ratio of N2 to N2O production by denitrification [gN/gN]
    real(r8), pointer    :: f_n2o_denit_vr(:,:)                         ! flux of N2O from denitrification [gN/m3/s]
    real(r8), pointer    :: f_n2o_nit_vr(:,:)                           ! flux of N2O from nitrification [gN/m3/s]
    real(r8), pointer    :: supplement_to_sminn_vr(:,:)
    real(r8), pointer    :: sminn_to_plant_vr(:,:)
    real(r8), pointer    :: potential_immob_vr(:,:)
    real(r8), pointer    :: actual_immob_vr(:,:)
    real(r8), intent(in) :: pmnf_decomp_cascade(begc:,1:,1:)            ! potential mineral N flux from one pool to another (gN/m3/s)
    real(r8), intent(in) :: p_decomp_cn_gain(begc:,1:,1:)               ! C:N ratio of the flux gained by the receiver pool
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: nitrif_n2o_loss_frac = 6.e-4_r8              ! fraction of N lost as N2O in nitrification (Li et al., 2000)
    integer  :: c,p,l,pi,j,k                                            ! indices
    integer  :: fc                                                      ! filter column index
    real(r8) :: amnf_immob_vr                                           ! actual mineral N flux from immobilization (gN/m3/s)
    real(r8) :: n_deficit_vr                                            ! microbial N deficit, vertically resolved (gN/m3/s)
    real(r8) :: fpi_no3_vr(begc:endc,1:nlevdecomp)                      ! fraction of potential immobilization supplied by no3
    real(r8) :: fpi_nh4_vr(begc:endc,1:nlevdecomp)                      ! fraction of potential immobilization supplied by nh4
    real(r8) :: sum_nh4_demand(begc:endc,1:nlevdecomp)
    real(r8) :: sum_nh4_demand_scaled(begc:endc,1:nlevdecomp)
    real(r8) :: sum_no3_demand(begc:endc,1:nlevdecomp)
    real(r8) :: sum_no3_demand_scaled(begc:endc,1:nlevdecomp)
    real(r8) :: sum_ndemand_vr(begc:endc, 1:nlevdecomp)                 ! total column N demand (gN/m3/s) at a given level
    real(r8) :: nuptake_prof(begc:endc, 1:nlevdecomp)
    real(r8) :: sminn_tot(begc:endc)
    integer  :: nlimit(begc:endc,0:nlevdecomp)                          ! flag for N limitation
    integer  :: nlimit_no3(begc:endc,0:nlevdecomp)                      ! flag for NO3 limitation
    integer  :: nlimit_nh4(begc:endc,0:nlevdecomp)                      ! flag for NH4 limitation
    real(r8) :: residual_sminn_vr(begc:endc, 1:nlevdecomp)
    real(r8) :: residual_sminn(begc:endc)
    real(r8) :: residual_smin_nh4_vr(begc:endc, 1:nlevdecomp)
    real(r8) :: residual_smin_no3_vr(begc:endc, 1:nlevdecomp)
    real(r8) :: residual_smin_nh4(begc:endc)
    real(r8) :: residual_smin_no3(begc:endc)
    real(r8) :: residual_plant_ndemand(begc:endc)
    !-----------------------------------------------------------------------

      ! column loops to resolve plant/heterotroph competition for mineral N

      if_nitrif: if (.not. use_nitrif_denitrif) then

         ! init sminn_tot
         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            sminn_tot(c) = 0.
         end do

         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               sminn_tot(c) = sminn_tot(c) + sminn_vr(c,j) * dzsoi_decomp(j)
            end do
         end do

         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               if (sminn_tot(c)  >  0.) then
                  nuptake_prof(c,j) = sminn_vr(c,j) / sminn_tot(c)
               else
                  nuptake_prof(c,j) = nfixation_prof(c,j)
               endif
            end do
         end do

         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               sum_ndemand_vr(c,j) = plant_ndemand(c) * nuptake_prof(c,j) + potential_immob_vr(c,j)
            end do
         end do

         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               l = landunit(c)
               if (sum_ndemand_vr(c,j)*dt < sminn_vr(c,j)) then

                  ! N availability is not limiting immobilization or plant
                  ! uptake, and both can proceed at their potential rates
                  nlimit(c,j) = 0
                  fpi_vr(c,j) = 1.0_r8
                  actual_immob_vr(c,j) = potential_immob_vr(c,j)
                  sminn_to_plant_vr(c,j) = plant_ndemand(c) * nuptake_prof(c,j)
               else if ( carbon_only ) then !.or. &
                  ! this code block controls the addition of N to sminn pool
                  ! to eliminate any N limitation, when Carbon_Only is set.  This lets the
                  ! model behave essentially as a carbon-only model, but with the
                  ! benefit of keeping track of the N additions needed to
                  ! eliminate N limitations, so there is still a diagnostic quantity
                  ! that describes the degree of N limitation at steady-state.

                  nlimit(c,j) = 1
                  fpi_vr(c,j) = 1.0_r8
                  actual_immob_vr(c,j) = potential_immob_vr(c,j)
                  sminn_to_plant_vr(c,j) =  plant_ndemand(c) * nuptake_prof(c,j)
                  supplement_to_sminn_vr(c,j) = sum_ndemand_vr(c,j) - (sminn_vr(c,j)/dt)
               else
                  ! N availability can not satisfy the sum of immobilization and
                  ! plant growth demands, so these two demands compete for available
                  ! soil mineral N resource.

                  nlimit(c,j) = 1
                  if (sum_ndemand_vr(c,j) > 0.0_r8) then
                     actual_immob_vr(c,j) = (sminn_vr(c,j)/dt)*(potential_immob_vr(c,j) / sum_ndemand_vr(c,j))
                  else
                     actual_immob_vr(c,j) = 0.0_r8
                  end if

                  if (potential_immob_vr(c,j) > 0.0_r8) then
                     fpi_vr(c,j) = actual_immob_vr(c,j) / potential_immob_vr(c,j)
                  else
                     fpi_vr(c,j) = 0.0_r8
                  end if

                  sminn_to_plant_vr(c,j) = (sminn_vr(c,j)/dt) - actual_immob_vr(c,j)
               end if
            end do
         end do

         ! sum up N fluxes to plant
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               sminn_to_plant(c) = sminn_to_plant(c) + sminn_to_plant_vr(c,j) * dzsoi_decomp(j)
            end do
         end do

         ! give plants a second pass to see if there is any mineral N left over with which to satisfy residual N demand.
         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            residual_sminn(c) = 0._r8
         end do

         ! sum up total N left over after initial plant and immobilization fluxes
         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            residual_plant_ndemand(c) = plant_ndemand(c) - sminn_to_plant(c)
         end do
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               if (residual_plant_ndemand(c)  >  0._r8 ) then
                  if (nlimit(c,j) .eq. 0) then
                     residual_sminn_vr(c,j) = max(sminn_vr(c,j) - (actual_immob_vr(c,j) + sminn_to_plant_vr(c,j) ) * dt, 0._r8)
                     residual_sminn(c) = residual_sminn(c) + residual_sminn_vr(c,j) * dzsoi_decomp(j)
                  else
                     residual_sminn_vr(c,j)  = 0._r8
                  endif
               endif
            end do
         end do

         ! distribute residual N to plants
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               if ( residual_plant_ndemand(c)  >  0._r8 .and. residual_sminn(c)  >  0._r8 .and. nlimit(c,j) .eq. 0) then
                  sminn_to_plant_vr(c,j) = sminn_to_plant_vr(c,j) + residual_sminn_vr(c,j) * &
                       min(( residual_plant_ndemand(c) *  dt ) / residual_sminn(c), 1._r8) / dt
               endif
            end do
         end do

         ! re-sum up N fluxes to plant
         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            sminn_to_plant(c) = 0._r8
         end do
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               sminn_to_plant(c) = sminn_to_plant(c) + sminn_to_plant_vr(c,j) * dzsoi_decomp(j)
               sum_ndemand_vr(c,j) = potential_immob_vr(c,j) + sminn_to_plant_vr(c,j)
            end do
         end do

         ! under conditions of excess N, some proportion is assumed to
         ! be lost to denitrification, in addition to the constant
         ! proportion lost in the decomposition pathways
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               if ((sminn_to_plant_vr(c,j) + actual_immob_vr(c,j))*dt < sminn_vr(c,j)) then
                  sminn_to_denit_excess_vr(c,j) = max(bdnr*((sminn_vr(c,j)/dt) - sum_ndemand_vr(c,j)),0._r8)
               else
                  sminn_to_denit_excess_vr(c,j) = 0._r8
               endif
            end do
         end do

         ! sum up N fluxes to immobilization
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               actual_immob(c) = actual_immob(c) + actual_immob_vr(c,j) * dzsoi_decomp(j)
               potential_immob(c) = potential_immob(c) + potential_immob_vr(c,j) * dzsoi_decomp(j)
            end do
         end do

         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            ! calculate the fraction of potential growth that can be
            ! acheived with the N available to plants
            if (plant_ndemand(c) > 0.0_r8) then
               fpg(c) = sminn_to_plant(c) / plant_ndemand(c)
            else
               fpg(c) = 1.0_r8
            end if

            ! calculate the fraction of immobilization realized (for diagnostic purposes)
            if (potential_immob(c) > 0.0_r8) then
               fpi(c) = actual_immob(c) / potential_immob(c)
            else
               fpi(c) = 1.0_r8
            end if
         end do

      else  !----------NITRIF_DENITRIF-------------!

         ! column loops to resolve plant/heterotroph/nitrifier/denitrifier competition for mineral N

         ! init total mineral N pools
         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            sminn_tot(c) = 0.
         end do

         ! sum up total mineral N pools
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               sminn_tot(c) = sminn_tot(c) + (smin_no3_vr(c,j) + smin_nh4_vr(c,j)) * dzsoi_decomp(j)
            end do
         end do

         ! define N uptake profile for initial vertical distribution of plant N uptake, assuming plant seeks N from where it is most abundant
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               if (sminn_tot(c)  >  0.) then
                  nuptake_prof(c,j) = sminn_vr(c,j) / sminn_tot(c)
               else
                  nuptake_prof(c,j) = nfixation_prof(c,j)
               endif
            end do
         end do

         ! main column/vertical loop
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               l = landunit(c)

               !  first compete for nh4
               sum_nh4_demand(c,j) = plant_ndemand(c) * nuptake_prof(c,j) + potential_immob_vr(c,j) + pot_f_nit_vr(c,j)
               sum_nh4_demand_scaled(c,j) = plant_ndemand(c)* nuptake_prof(c,j) * compet_plant_nh4 + &
                    potential_immob_vr(c,j)*compet_decomp_nh4 + pot_f_nit_vr(c,j)*compet_nit

               if (sum_nh4_demand(c,j)*dt < smin_nh4_vr(c,j)) then

                  ! NH4 availability is not limiting immobilization or plant
                  ! uptake, and all can proceed at their potential rates
                  nlimit_nh4(c,j) = 0
                  fpi_nh4_vr(c,j) = 1.0_r8
                  actual_immob_nh4_vr(c,j) = potential_immob_vr(c,j)
                  !RF added new term.

                  f_nit_vr(c,j) = pot_f_nit_vr(c,j)

                  smin_nh4_to_plant_vr(c,j) = plant_ndemand(c) * nuptake_prof(c,j)

               else

                  ! NH4 availability can not satisfy the sum of immobilization, nitrification, and
                  ! plant growth demands, so these three demands compete for available
                  ! soil mineral NH4 resource.
                  nlimit_nh4(c,j) = 1
                  if (sum_nh4_demand(c,j) > 0.0_r8) then
                  ! RF microbes compete based on the hypothesised plant demand.
                     actual_immob_nh4_vr(c,j) = min((smin_nh4_vr(c,j)/dt)*(potential_immob_vr(c,j)* &
                          compet_decomp_nh4 / sum_nh4_demand_scaled(c,j)), potential_immob_vr(c,j))

                     f_nit_vr(c,j) =  min((smin_nh4_vr(c,j)/dt)*(pot_f_nit_vr(c,j)*compet_nit / &
                          sum_nh4_demand_scaled(c,j)), pot_f_nit_vr(c,j))

                     smin_nh4_to_plant_vr(c,j) = min((smin_nh4_vr(c,j)/dt)*(plant_ndemand(c)* &
                      nuptake_prof(c,j)*compet_plant_nh4 / sum_nh4_demand_scaled(c,j)), plant_ndemand(c)*nuptake_prof(c,j))

                  else
                     actual_immob_nh4_vr(c,j) = 0.0_r8
                     smin_nh4_to_plant_vr(c,j) = 0.0_r8
                     f_nit_vr(c,j) = 0.0_r8
                  end if

                  if (potential_immob_vr(c,j) > 0.0_r8) then
                     fpi_nh4_vr(c,j) = actual_immob_nh4_vr(c,j) / potential_immob_vr(c,j)
                  else
                     fpi_nh4_vr(c,j) = 0.0_r8
                  end if

               end if

               if (decomp_method == mimics_decomp) then
                  ! turn off fpi for MIMICS and only lets plants
                  ! take up available mineral nitrogen.
                  ! TODO slevis: -ve or tiny sminn_vr could cause problems
                  fpi_nh4_vr(c,j) = 1.0_r8
                  actual_immob_nh4_vr(c,j) = potential_immob_vr(c,j)
               end if

               sum_no3_demand(c,j) = (plant_ndemand(c)*nuptake_prof(c,j)-smin_nh4_to_plant_vr(c,j)) + &
              (potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j)) + pot_f_denit_vr(c,j)
               sum_no3_demand_scaled(c,j) = (plant_ndemand(c)*nuptake_prof(c,j) &
                                             -smin_nh4_to_plant_vr(c,j))*compet_plant_no3 + &
              (potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j))*compet_decomp_no3 + pot_f_denit_vr(c,j)*compet_denit

               if (sum_no3_demand(c,j)*dt < smin_no3_vr(c,j)) then

                  ! NO3 availability is not limiting immobilization or plant
                  ! uptake, and all can proceed at their potential rates
                  nlimit_no3(c,j) = 0
                  fpi_no3_vr(c,j) = 1.0_r8 -  fpi_nh4_vr(c,j)
                  actual_immob_no3_vr(c,j) = (potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j))

                  f_denit_vr(c,j) = pot_f_denit_vr(c,j)

                  smin_no3_to_plant_vr(c,j) = (plant_ndemand(c)*nuptake_prof(c,j)-smin_nh4_to_plant_vr(c,j))

               else

                  ! NO3 availability can not satisfy the sum of immobilization, denitrification, and
                  ! plant growth demands, so these three demands compete for available
                  ! soil mineral NO3 resource.
                  nlimit_no3(c,j) = 1

                  if (sum_no3_demand(c,j) > 0.0_r8) then
                     actual_immob_no3_vr(c,j) = min((smin_no3_vr(c,j)/dt)*((potential_immob_vr(c,j)- &
                     actual_immob_nh4_vr(c,j))*compet_decomp_no3 / sum_no3_demand_scaled(c,j)), &
                               potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j))

                     smin_no3_to_plant_vr(c,j) = min((smin_no3_vr(c,j)/dt)*((plant_ndemand(c)* &
                               nuptake_prof(c,j)-smin_nh4_to_plant_vr(c,j))*compet_plant_no3 / sum_no3_demand_scaled(c,j)), &
                               plant_ndemand(c)*nuptake_prof(c,j)-smin_nh4_to_plant_vr(c,j))

                     f_denit_vr(c,j) = min((smin_no3_vr(c,j)/dt)*(pot_f_denit_vr(c,j)*compet_denit / &
                               sum_no3_demand_scaled(c,j)), pot_f_denit_vr(c,j))

                  else ! no no3 demand. no uptake fluxes.
                     actual_immob_no3_vr(c,j) = 0.0_r8
                     smin_no3_to_plant_vr(c,j) = 0.0_r8
                     f_denit_vr(c,j) = 0.0_r8

                  end if !any no3 demand?




                  if (potential_immob_vr(c,j) > 0.0_r8) then
                     fpi_no3_vr(c,j) = actual_immob_no3_vr(c,j) / potential_immob_vr(c,j)
                  else
                     fpi_no3_vr(c,j) = 0.0_r8
                  end if

               end if

               if (decomp_method == mimics_decomp) then
                  ! turn off fpi for MIMICS and only lets plants
                  ! take up available mineral nitrogen.
                  ! TODO slevis: -ve or tiny sminn_vr could cause problems
                  fpi_no3_vr(c,j) = 1.0_r8 - fpi_nh4_vr(c,j)  ! => 0
                  actual_immob_no3_vr(c,j) = potential_immob_vr(c,j) - &
                                             actual_immob_nh4_vr(c,j)  ! => 0
               end if

               ! n2o emissions: n2o from nitr is const fraction, n2o from denitr is calculated in nitrif_denitrif
               f_n2o_nit_vr(c,j) = f_nit_vr(c,j) * nitrif_n2o_loss_frac
               f_n2o_denit_vr(c,j) = f_denit_vr(c,j) / (1._r8 + n2_n2o_ratio_denit_vr(c,j))


               ! this code block controls the addition of N to sminn pool
               ! to eliminate any N limitation, when Carbon_Only is set.  This lets the
               ! model behave essentially as a carbon-only model, but with the
               ! benefit of keeping track of the N additions needed to
               ! eliminate N limitations, so there is still a diagnostic quantity
               ! that describes the degree of N limitation at steady-state.

               if ( carbon_only ) then !.or. &
                  if ( fpi_no3_vr(c,j) + fpi_nh4_vr(c,j) < 1._r8 ) then
                     fpi_nh4_vr(c,j) = 1.0_r8 - fpi_no3_vr(c,j)
                     supplement_to_sminn_vr(c,j) = (potential_immob_vr(c,j) &
                                                  - actual_immob_no3_vr(c,j)) - actual_immob_nh4_vr(c,j)
                     ! update to new values that satisfy demand
                     actual_immob_nh4_vr(c,j) = potential_immob_vr(c,j) -  actual_immob_no3_vr(c,j)
                  end if
                  if ( smin_no3_to_plant_vr(c,j) + smin_nh4_to_plant_vr(c,j) < plant_ndemand(c)*nuptake_prof(c,j) ) then
                     supplement_to_sminn_vr(c,j) = supplement_to_sminn_vr(c,j) + &
                          (plant_ndemand(c)*nuptake_prof(c,j) - smin_no3_to_plant_vr(c,j)) - smin_nh4_to_plant_vr(c,j)  ! use old values
                     smin_nh4_to_plant_vr(c,j) = plant_ndemand(c)*nuptake_prof(c,j) - smin_no3_to_plant_vr(c,j)
                  end if
                  sminn_to_plant_vr(c,j) = smin_no3_to_plant_vr(c,j) + smin_nh4_to_plant_vr(c,j)
               end if

               ! sum up no3 and nh4 fluxes
               fpi_vr(c,j) = fpi_no3_vr(c,j) + fpi_nh4_vr(c,j)
               sminn_to_plant_vr(c,j) = smin_no3_to_plant_vr(c,j) + smin_nh4_to_plant_vr(c,j)
               actual_immob_vr(c,j) = actual_immob_no3_vr(c,j) + actual_immob_nh4_vr(c,j)
            end do
         end do

         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            ! sum up N fluxes to plant after initial competition
            sminn_to_plant(c) = 0._r8
         end do
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               sminn_to_plant(c) = sminn_to_plant(c) + sminn_to_plant_vr(c,j) * dzsoi_decomp(j)
            end do
         end do

         if (decomp_method == mimics_decomp) then
            do j = 1, nlevdecomp
               do fc=1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)

                  do k = 1, ndecomp_cascade_transitions
                     if (cascade_receiver_pool(k) == i_cop_mic .or. &
                         cascade_receiver_pool(k) == i_oli_mic) then
                        sum_ndemand_vr(c,j) = sum_no3_demand_scaled(c,j) + &
                                              sum_nh4_demand_scaled(c,j)
                        if (pmnf_decomp_cascade(c,j,k) > 0.0_r8 .and. &
                            sum_ndemand_vr(c,j) > 0.0_r8) then
                           amnf_immob_vr = (sminn_vr(c,j) / dt) * &
                                           (pmnf_decomp_cascade(c,j,k) / &
                                            sum_ndemand_vr(c,j))
                           n_deficit_vr = pmnf_decomp_cascade(c,j,k) - &
                                          amnf_immob_vr
                           c_overflow_vr(c,j,k) = &
                              n_deficit_vr * p_decomp_cn_gain(c,j,cascade_receiver_pool(k))
                        else  ! not pmnf and sum_ndemand > 0
                           c_overflow_vr(c,j,k) = 0.0_r8
                        end if
                     else  ! not microbes receiving
                        c_overflow_vr(c,j,k) = 0.0_r8
                     end if
                  end do
               end do
            end do
         else  ! not mimics_decomp
            c_overflow_vr(:,:,:) = 0.0_r8
         end if

         ! give plants a second pass to see if there is any mineral N left over with which to satisfy residual N demand.
         ! first take frm nh4 pool; then take from no3 pool
         do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               residual_plant_ndemand(c) = plant_ndemand(c) - sminn_to_plant(c)
               residual_smin_nh4(c) = 0._r8
            end do
            do j = 1, nlevdecomp
               do fc=1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)
                  if (residual_plant_ndemand(c)  >  0._r8 ) then
                     if (nlimit_nh4(c,j) .eq. 0) then
                        residual_smin_nh4_vr(c,j) = max(smin_nh4_vr(c,j) - (actual_immob_nh4_vr(c,j) + &
                                                    smin_nh4_to_plant_vr(c,j) + f_nit_vr(c,j) ) * dt, 0._r8)

                        residual_smin_nh4(c) = residual_smin_nh4(c) + residual_smin_nh4_vr(c,j) * dzsoi_decomp(j)
                     else
                        residual_smin_nh4_vr(c,j)  = 0._r8
                     endif

                     if ( residual_smin_nh4(c) > 0._r8 .and. nlimit_nh4(c,j) .eq. 0 ) then
                        smin_nh4_to_plant_vr(c,j) = smin_nh4_to_plant_vr(c,j) + residual_smin_nh4_vr(c,j) * &
                             min(( residual_plant_ndemand(c) *  dt ) / residual_smin_nh4(c), 1._r8) / dt
                     endif
                  end if
               end do
            end do

            ! re-sum up N fluxes to plant after second pass for nh4
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               sminn_to_plant(c) = 0._r8
            end do
            do j = 1, nlevdecomp
               do fc=1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)
                  sminn_to_plant_vr(c,j) = smin_nh4_to_plant_vr(c,j) + smin_no3_to_plant_vr(c,j)
                  sminn_to_plant(c) = sminn_to_plant(c) + (sminn_to_plant_vr(c,j)) * dzsoi_decomp(j)
               end do
            end do

            !
            ! and now do second pass for no3
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               residual_plant_ndemand(c) = plant_ndemand(c) - sminn_to_plant(c)
               residual_smin_no3(c) = 0._r8
            end do

            do j = 1, nlevdecomp
               do fc=1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)
                  if (residual_plant_ndemand(c) > 0._r8 ) then
                     if (nlimit_no3(c,j) .eq. 0) then
                       residual_smin_no3_vr(c,j) = max(smin_no3_vr(c,j) - (actual_immob_no3_vr(c,j) + &
                                                   smin_no3_to_plant_vr(c,j) + f_denit_vr(c,j) ) * dt, 0._r8)
                        residual_smin_no3(c) = residual_smin_no3(c) + residual_smin_no3_vr(c,j) * dzsoi_decomp(j)
                     else
                        residual_smin_no3_vr(c,j)  = 0._r8
                     endif

                     if ( residual_smin_no3(c) > 0._r8 .and. nlimit_no3(c,j) .eq. 0) then
                        smin_no3_to_plant_vr(c,j) = smin_no3_to_plant_vr(c,j) + residual_smin_no3_vr(c,j) * &
                             min(( residual_plant_ndemand(c) *  dt ) / residual_smin_no3(c), 1._r8) / dt
                     endif
                  endif
               end do
            end do

            ! re-sum up N fluxes to plant after second passes of both no3 and nh4
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               sminn_to_plant(c) = 0._r8
            end do
            do j = 1, nlevdecomp
               do fc=1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)
                  sminn_to_plant_vr(c,j) = smin_nh4_to_plant_vr(c,j) + smin_no3_to_plant_vr(c,j)
                  sminn_to_plant(c) = sminn_to_plant(c) + (sminn_to_plant_vr(c,j)) * dzsoi_decomp(j)
               end do
            end do

         ! sum up N fluxes to immobilization
         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            actual_immob(c) = 0._r8
            potential_immob(c) = 0._r8
         end do
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               actual_immob(c) = actual_immob(c) + actual_immob_vr(c,j) * dzsoi_decomp(j)
               potential_immob(c) = potential_immob(c) + potential_immob_vr(c,j) * dzsoi_decomp(j)
            end do
         end do




         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            ! calculate the fraction of potential growth that can be
            ! acheived with the N available to plants
            ! calculate the fraction of immobilization realized (for diagnostic purposes)
            if (plant_ndemand(c) > 0.0_r8) then
               fpg(c) = sminn_to_plant(c) / plant_ndemand(c)
            else
               fpg(c) = 1._r8
            end if

            if (potential_immob(c) > 0.0_r8) then
               fpi(c) = actual_immob(c) / potential_immob(c)
            else
               fpi(c) = 1._r8
            end if
         end do ! end of column loops

      end if if_nitrif  !end of if_not_use_nitrif_denitrif

  end subroutine SoilBiogeochemCompetition

end module SoilBiogeochemCompetition_mod
