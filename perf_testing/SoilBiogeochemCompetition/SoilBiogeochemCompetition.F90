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
    use perf_timers_mod, only : perf_timer_start, perf_timer_stop
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

         ! Initialize sminn_tot
         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            sminn_tot(c) = 0.
         end do

         ! Get total soil mineral N
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               sminn_tot(c) = sminn_tot(c) + sminn_vr(c,j) * dzsoi_decomp(j)
            end do
         end do

         ! Get N uptake profile (fraction of plant uptake coming from each soil layer)
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

         ! Get total column N demand from each soil layer
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               sum_ndemand_vr(c,j) = plant_ndemand(c) * nuptake_prof(c,j) + potential_immob_vr(c,j)
            end do
         end do

         ! Get actual plant N uptake from each soil layer
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

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Give plants a second pass to see if there is any mineral N left over
         ! with which to satisfy residual N demand.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! sum up total N left over after initial plant and immobilization fluxes
         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            residual_sminn(c) = 0._r8
         end do
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

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Done with second pass
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

         ! Hoisted !$acc data region for read-only inputs that are constant
         ! across all kernels in the canonical path. Inner kernels reference
         ! these as 'present' so transfers happen once per call to this
         ! routine (i.e. once per timestep in the real model), not once per
         ! kernel. The copyin list grows as more kernels get OpenACC-ified.
         !$acc data copyin(smin_nh4_vr, smin_no3_vr,             &
         !$acc&            dzsoi_decomp, filter_bgc_soilc,       &
         !$acc&            sminn_vr, nfixation_prof)

         ! column loops to resolve plant/heterotroph/nitrifier/denitrifier competition for mineral N

         ! init total mineral N pools
         call perf_timer_start('init_sminn_tot')
         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            sminn_tot(c) = 0.
         end do
         call perf_timer_stop('init_sminn_tot')

         ! Inner data region scoped to the two kernels that use the
         ! routine-local automatics sminn_tot and nuptake_prof. sminn_tot
         ! was just initialized to zero on the host (init_sminn_tot
         ! above); copyin brings those zeros to device. accum_sminn_tot
         ! then accumulates into it on device, and compute_nuptake_prof
         ! reads it on device — no host round-trip between the two
         ! kernels. nuptake_prof is written on device and copied out at
         ! region end so the host-side main_competition (still on CPU)
         ! sees the final values.
         !$acc data copyin(sminn_tot) copyout(nuptake_prof)

         ! sum up total mineral N pools.
         ! GPU/multicore (_OPENACC): parallelize over fc, serialize j inside
         ! each thread (sminn_tot(c) is accumulated across j for each c —
         ! keep that reduction serial per-thread). CPU-serial: original loop
         ! order (j outer, fc inner) is more cache-friendly because
         ! smin_no3_vr(c,j) etc. are column-major. Body and end-do's are
         ! shared; only the loop opening differs.
         call perf_timer_start('accum_sminn_tot')
#ifdef _OPENACC
         !$acc parallel loop default(present)
         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            do j = 1, nlevdecomp
#else
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
#endif
               call accum_sminn_tot(sminn_tot(c), smin_no3_vr(c,j), smin_nh4_vr(c,j), dzsoi_decomp(j))
            end do
         end do
         call perf_timer_stop('accum_sminn_tot')

         ! define N uptake profile for initial vertical distribution of plant N uptake, assuming plant seeks N from where it is most abundant.
         ! Each (c,j) writes to its own nuptake_prof(c,j); no reduction —
         ! safe to parallelize both loops together via collapse(2).
         call perf_timer_start('compute_nuptake_prof')
         !$acc parallel loop collapse(2) default(present)
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               call compute_nuptake_prof(nuptake_prof(c,j), sminn_tot(c), sminn_vr(c,j), nfixation_prof(c,j))
            end do
         end do
         call perf_timer_stop('compute_nuptake_prof')

         !$acc end data

         ! main column/vertical loop
         call perf_timer_start('main_competition')
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               l = landunit(c)

               !  first compete for nh4
               call compete_nh4( &
                    sum_nh4_demand(c,j), sum_nh4_demand_scaled(c,j), nlimit_nh4(c,j), &
                    fpi_nh4_vr(c,j), actual_immob_nh4_vr(c,j), &
                    f_nit_vr(c,j), smin_nh4_to_plant_vr(c,j), &
                    plant_ndemand(c), nuptake_prof(c,j), &
                    potential_immob_vr(c,j), pot_f_nit_vr(c,j), smin_nh4_vr(c,j), &
                    dt, compet_plant_nh4, compet_decomp_nh4, compet_nit, &
                    decomp_method, mimics_decomp)

               ! then compete for no3
               call compete_no3( &
                    sum_no3_demand(c,j), sum_no3_demand_scaled(c,j), nlimit_no3(c,j), &
                    fpi_no3_vr(c,j), actual_immob_no3_vr(c,j), &
                    f_denit_vr(c,j), smin_no3_to_plant_vr(c,j), &
                    plant_ndemand(c), nuptake_prof(c,j), &
                    smin_nh4_to_plant_vr(c,j), actual_immob_nh4_vr(c,j), fpi_nh4_vr(c,j), &
                    potential_immob_vr(c,j), pot_f_denit_vr(c,j), smin_no3_vr(c,j), &
                    dt, compet_plant_no3, compet_decomp_no3, compet_denit, &
                    decomp_method, mimics_decomp)

               ! n2o emissions: n2o from nitr is const fraction, n2o from denitr is calculated in nitrif_denitrif
               call compute_n2o_emissions( &
                    f_n2o_nit_vr(c,j), f_n2o_denit_vr(c,j), &
                    f_nit_vr(c,j), f_denit_vr(c,j), n2_n2o_ratio_denit_vr(c,j), &
                    nitrif_n2o_loss_frac)


               ! this code block controls the addition of N to sminn pool
               ! to eliminate any N limitation, when Carbon_Only is set.  This lets the
               ! model behave essentially as a carbon-only model, but with the
               ! benefit of keeping track of the N additions needed to
               ! eliminate N limitations, so there is still a diagnostic quantity
               ! that describes the degree of N limitation at steady-state.
               call apply_carbon_only_adjustment( &
                    fpi_nh4_vr(c,j), supplement_to_sminn_vr(c,j), &
                    actual_immob_nh4_vr(c,j), smin_nh4_to_plant_vr(c,j), &
                    sminn_to_plant_vr(c,j), &
                    fpi_no3_vr(c,j), actual_immob_no3_vr(c,j), &
                    smin_no3_to_plant_vr(c,j), &
                    potential_immob_vr(c,j), plant_ndemand(c), nuptake_prof(c,j), &
                    carbon_only)

               ! sum up no3 and nh4 fluxes
               call compute_competition_summary( &
                    fpi_vr(c,j), sminn_to_plant_vr(c,j), actual_immob_vr(c,j), &
                    fpi_no3_vr(c,j), fpi_nh4_vr(c,j), &
                    smin_no3_to_plant_vr(c,j), smin_nh4_to_plant_vr(c,j), &
                    actual_immob_no3_vr(c,j), actual_immob_nh4_vr(c,j))
            end do
         end do
         call perf_timer_stop('main_competition')

         ! sum up N fluxes to plant after initial competition
         call perf_timer_start('sum_sminn_to_plant')
         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            sminn_to_plant(c) = 0._r8
         end do
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               call accum_dz_weighted(sminn_to_plant(c), sminn_to_plant_vr(c,j), dzsoi_decomp(j))
            end do
         end do
         call perf_timer_stop('sum_sminn_to_plant')

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
         call perf_timer_start('residual_uptake_nh4')
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
                        residual_smin_nh4_vr(c,j) = compute_residual_smin_vr( &
                             smin_nh4_vr(c,j), actual_immob_nh4_vr(c,j), smin_nh4_to_plant_vr(c,j), f_nit_vr(c,j), dt)

                        residual_smin_nh4(c) = residual_smin_nh4(c) + residual_smin_nh4_vr(c,j) * dzsoi_decomp(j)
                     else
                        residual_smin_nh4_vr(c,j)  = 0._r8
                     endif

                     if ( residual_smin_nh4(c) > 0._r8 .and. nlimit_nh4(c,j) .eq. 0 ) then
                        smin_nh4_to_plant_vr(c,j) = distribute_residual_to_plant( &
                             smin_nh4_to_plant_vr(c,j), residual_smin_nh4_vr(c,j), residual_plant_ndemand(c), residual_smin_nh4(c), dt)
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
            call perf_timer_stop('residual_uptake_nh4')

            !
            ! and now do second pass for no3
            call perf_timer_start('residual_uptake_no3')
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
                       residual_smin_no3_vr(c,j) = compute_residual_smin_vr( &
                             smin_no3_vr(c,j), actual_immob_no3_vr(c,j), smin_no3_to_plant_vr(c,j), f_denit_vr(c,j), dt)
                        residual_smin_no3(c) = residual_smin_no3(c) + residual_smin_no3_vr(c,j) * dzsoi_decomp(j)
                     else
                        residual_smin_no3_vr(c,j)  = 0._r8
                     endif

                     if ( residual_smin_no3(c) > 0._r8 .and. nlimit_no3(c,j) .eq. 0) then
                        smin_no3_to_plant_vr(c,j) = distribute_residual_to_plant( &
                             smin_no3_to_plant_vr(c,j), residual_smin_no3_vr(c,j), residual_plant_ndemand(c), residual_smin_no3(c), dt)
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
            call perf_timer_stop('residual_uptake_no3')

         ! sum up N fluxes to immobilization
         call perf_timer_start('sum_immobilization')
         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            actual_immob(c) = 0._r8
            potential_immob(c) = 0._r8
         end do
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               call accum_dz_weighted(actual_immob(c),    actual_immob_vr(c,j),    dzsoi_decomp(j))
               call accum_dz_weighted(potential_immob(c), potential_immob_vr(c,j), dzsoi_decomp(j))
            end do
         end do
         call perf_timer_stop('sum_immobilization')



         call perf_timer_start('compute_fpg_fpi')
         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            ! calculate the fraction of potential growth that can be
            ! acheived with the N available to plants
            ! calculate the fraction of immobilization realized (for diagnostic purposes)
            fpg(c) = compute_fraction_or_one(sminn_to_plant(c), plant_ndemand(c))
            fpi(c) = compute_fraction_or_one(actual_immob(c),   potential_immob(c))
         end do ! end of column loops
         call perf_timer_stop('compute_fpg_fpi')

         !$acc end data

      end if if_nitrif  !end of if_not_use_nitrif_denitrif

  end subroutine SoilBiogeochemCompetition

  !-----------------------------------------------------------------------
  pure subroutine accum_sminn_tot(sminn_tot, smin_no3_vr, smin_nh4_vr, dzsoi_decomp)
    !$acc routine seq
    real(r8), intent(inout) :: sminn_tot
    real(r8), intent(in)    :: smin_no3_vr, smin_nh4_vr, dzsoi_decomp
    sminn_tot = sminn_tot + (smin_no3_vr + smin_nh4_vr) * dzsoi_decomp
  end subroutine accum_sminn_tot

  !-----------------------------------------------------------------------
  pure subroutine compute_nuptake_prof(nuptake_prof, sminn_tot, sminn_vr, nfixation_prof)
    !$acc routine seq
    real(r8), intent(out) :: nuptake_prof
    real(r8), intent(in)  :: sminn_tot, sminn_vr, nfixation_prof
    if (sminn_tot  >  0.) then
       nuptake_prof = sminn_vr / sminn_tot
    else
       nuptake_prof = nfixation_prof
    endif
  end subroutine compute_nuptake_prof

  !-----------------------------------------------------------------------
  pure subroutine compete_nh4( &
       sum_nh4_demand, sum_nh4_demand_scaled, nlimit_nh4, &
       fpi_nh4_vr, actual_immob_nh4_vr, &
       f_nit_vr, smin_nh4_to_plant_vr, &
       plant_ndemand, nuptake_prof, &
       potential_immob_vr, pot_f_nit_vr, smin_nh4_vr, &
       dt, compet_plant_nh4, compet_decomp_nh4, compet_nit, &
       decomp_method, mimics_decomp)
    real(r8), intent(out) :: sum_nh4_demand, sum_nh4_demand_scaled
    integer , intent(out) :: nlimit_nh4
    real(r8), intent(out) :: fpi_nh4_vr, actual_immob_nh4_vr
    real(r8), intent(out) :: f_nit_vr, smin_nh4_to_plant_vr
    real(r8), intent(in)  :: plant_ndemand, nuptake_prof
    real(r8), intent(in)  :: potential_immob_vr, pot_f_nit_vr, smin_nh4_vr
    real(r8), intent(in)  :: dt, compet_plant_nh4, compet_decomp_nh4, compet_nit
    integer , intent(in)  :: decomp_method, mimics_decomp

               sum_nh4_demand = plant_ndemand * nuptake_prof + potential_immob_vr + pot_f_nit_vr
               sum_nh4_demand_scaled = plant_ndemand* nuptake_prof * compet_plant_nh4 + &
                    potential_immob_vr*compet_decomp_nh4 + pot_f_nit_vr*compet_nit

               if (sum_nh4_demand*dt < smin_nh4_vr) then

                  ! NH4 availability is not limiting immobilization or plant
                  ! uptake, and all can proceed at their potential rates
                  nlimit_nh4 = 0
                  fpi_nh4_vr = 1.0_r8
                  actual_immob_nh4_vr = potential_immob_vr
                  !RF added new term.

                  f_nit_vr = pot_f_nit_vr

                  smin_nh4_to_plant_vr = plant_ndemand * nuptake_prof

               else

                  ! NH4 availability can not satisfy the sum of immobilization, nitrification, and
                  ! plant growth demands, so these three demands compete for available
                  ! soil mineral NH4 resource.
                  nlimit_nh4 = 1
                  if (sum_nh4_demand > 0.0_r8) then
                  ! RF microbes compete based on the hypothesised plant demand.
                     actual_immob_nh4_vr = min((smin_nh4_vr/dt)*(potential_immob_vr* &
                          compet_decomp_nh4 / sum_nh4_demand_scaled), potential_immob_vr)

                     f_nit_vr =  min((smin_nh4_vr/dt)*(pot_f_nit_vr*compet_nit / &
                          sum_nh4_demand_scaled), pot_f_nit_vr)

                     smin_nh4_to_plant_vr = min((smin_nh4_vr/dt)*(plant_ndemand* &
                      nuptake_prof*compet_plant_nh4 / sum_nh4_demand_scaled), plant_ndemand*nuptake_prof)

                  else
                     actual_immob_nh4_vr = 0.0_r8
                     smin_nh4_to_plant_vr = 0.0_r8
                     f_nit_vr = 0.0_r8
                  end if

                  if (potential_immob_vr > 0.0_r8) then
                     fpi_nh4_vr = actual_immob_nh4_vr / potential_immob_vr
                  else
                     fpi_nh4_vr = 0.0_r8
                  end if

               end if

               if (decomp_method == mimics_decomp) then
                  ! turn off fpi for MIMICS and only lets plants
                  ! take up available mineral nitrogen.
                  ! TODO slevis: -ve or tiny sminn_vr could cause problems
                  fpi_nh4_vr = 1.0_r8
                  actual_immob_nh4_vr = potential_immob_vr
               end if
  end subroutine compete_nh4

  !-----------------------------------------------------------------------
  pure subroutine compete_no3( &
       sum_no3_demand, sum_no3_demand_scaled, nlimit_no3, &
       fpi_no3_vr, actual_immob_no3_vr, &
       f_denit_vr, smin_no3_to_plant_vr, &
       plant_ndemand, nuptake_prof, &
       smin_nh4_to_plant_vr, actual_immob_nh4_vr, fpi_nh4_vr, &
       potential_immob_vr, pot_f_denit_vr, smin_no3_vr, &
       dt, compet_plant_no3, compet_decomp_no3, compet_denit, &
       decomp_method, mimics_decomp)
    real(r8), intent(out) :: sum_no3_demand, sum_no3_demand_scaled
    integer , intent(out) :: nlimit_no3
    real(r8), intent(out) :: fpi_no3_vr, actual_immob_no3_vr
    real(r8), intent(out) :: f_denit_vr, smin_no3_to_plant_vr
    real(r8), intent(in)  :: plant_ndemand, nuptake_prof
    real(r8), intent(in)  :: smin_nh4_to_plant_vr, actual_immob_nh4_vr, fpi_nh4_vr
    real(r8), intent(in)  :: potential_immob_vr, pot_f_denit_vr, smin_no3_vr
    real(r8), intent(in)  :: dt, compet_plant_no3, compet_decomp_no3, compet_denit
    integer , intent(in)  :: decomp_method, mimics_decomp

               sum_no3_demand = (plant_ndemand*nuptake_prof-smin_nh4_to_plant_vr) + &
              (potential_immob_vr-actual_immob_nh4_vr) + pot_f_denit_vr
               sum_no3_demand_scaled = (plant_ndemand*nuptake_prof &
                                             -smin_nh4_to_plant_vr)*compet_plant_no3 + &
              (potential_immob_vr-actual_immob_nh4_vr)*compet_decomp_no3 + pot_f_denit_vr*compet_denit

               if (sum_no3_demand*dt < smin_no3_vr) then

                  ! NO3 availability is not limiting immobilization or plant
                  ! uptake, and all can proceed at their potential rates
                  nlimit_no3 = 0
                  fpi_no3_vr = 1.0_r8 -  fpi_nh4_vr
                  actual_immob_no3_vr = (potential_immob_vr-actual_immob_nh4_vr)

                  f_denit_vr = pot_f_denit_vr

                  smin_no3_to_plant_vr = (plant_ndemand*nuptake_prof-smin_nh4_to_plant_vr)

               else

                  ! NO3 availability can not satisfy the sum of immobilization, denitrification, and
                  ! plant growth demands, so these three demands compete for available
                  ! soil mineral NO3 resource.
                  nlimit_no3 = 1

                  if (sum_no3_demand > 0.0_r8) then
                     actual_immob_no3_vr = min((smin_no3_vr/dt)*((potential_immob_vr- &
                     actual_immob_nh4_vr)*compet_decomp_no3 / sum_no3_demand_scaled), &
                               potential_immob_vr-actual_immob_nh4_vr)

                     smin_no3_to_plant_vr = min((smin_no3_vr/dt)*((plant_ndemand* &
                               nuptake_prof-smin_nh4_to_plant_vr)*compet_plant_no3 / sum_no3_demand_scaled), &
                               plant_ndemand*nuptake_prof-smin_nh4_to_plant_vr)

                     f_denit_vr = min((smin_no3_vr/dt)*(pot_f_denit_vr*compet_denit / &
                               sum_no3_demand_scaled), pot_f_denit_vr)

                  else ! no no3 demand. no uptake fluxes.
                     actual_immob_no3_vr = 0.0_r8
                     smin_no3_to_plant_vr = 0.0_r8
                     f_denit_vr = 0.0_r8

                  end if !any no3 demand?




                  if (potential_immob_vr > 0.0_r8) then
                     fpi_no3_vr = actual_immob_no3_vr / potential_immob_vr
                  else
                     fpi_no3_vr = 0.0_r8
                  end if

               end if

               if (decomp_method == mimics_decomp) then
                  ! turn off fpi for MIMICS and only lets plants
                  ! take up available mineral nitrogen.
                  ! TODO slevis: -ve or tiny sminn_vr could cause problems
                  fpi_no3_vr = 1.0_r8 - fpi_nh4_vr  ! => 0
                  actual_immob_no3_vr = potential_immob_vr - &
                                             actual_immob_nh4_vr  ! => 0
               end if
  end subroutine compete_no3

  !-----------------------------------------------------------------------
  pure subroutine compute_n2o_emissions( &
       f_n2o_nit_vr, f_n2o_denit_vr, &
       f_nit_vr, f_denit_vr, n2_n2o_ratio_denit_vr, &
       nitrif_n2o_loss_frac)
    real(r8), intent(out) :: f_n2o_nit_vr, f_n2o_denit_vr
    real(r8), intent(in)  :: f_nit_vr, f_denit_vr, n2_n2o_ratio_denit_vr
    real(r8), intent(in)  :: nitrif_n2o_loss_frac
               f_n2o_nit_vr = f_nit_vr * nitrif_n2o_loss_frac
               f_n2o_denit_vr = f_denit_vr / (1._r8 + n2_n2o_ratio_denit_vr)
  end subroutine compute_n2o_emissions

  !-----------------------------------------------------------------------
  pure subroutine apply_carbon_only_adjustment( &
       fpi_nh4_vr, supplement_to_sminn_vr, &
       actual_immob_nh4_vr, smin_nh4_to_plant_vr, &
       sminn_to_plant_vr, &
       fpi_no3_vr, actual_immob_no3_vr, &
       smin_no3_to_plant_vr, &
       potential_immob_vr, plant_ndemand, nuptake_prof, &
       carbon_only)
    real(r8), intent(inout) :: fpi_nh4_vr, supplement_to_sminn_vr
    real(r8), intent(inout) :: actual_immob_nh4_vr, smin_nh4_to_plant_vr
    real(r8), intent(inout) :: sminn_to_plant_vr
    real(r8), intent(in)    :: fpi_no3_vr, actual_immob_no3_vr
    real(r8), intent(in)    :: smin_no3_to_plant_vr
    real(r8), intent(in)    :: potential_immob_vr, plant_ndemand, nuptake_prof
    logical , intent(in)    :: carbon_only

               if ( carbon_only ) then !.or. &
                  if ( fpi_no3_vr + fpi_nh4_vr < 1._r8 ) then
                     fpi_nh4_vr = 1.0_r8 - fpi_no3_vr
                     supplement_to_sminn_vr = (potential_immob_vr &
                                                  - actual_immob_no3_vr) - actual_immob_nh4_vr
                     ! update to new values that satisfy demand
                     actual_immob_nh4_vr = potential_immob_vr -  actual_immob_no3_vr
                  end if
                  if ( smin_no3_to_plant_vr + smin_nh4_to_plant_vr < plant_ndemand*nuptake_prof ) then
                     supplement_to_sminn_vr = supplement_to_sminn_vr + &
                          (plant_ndemand*nuptake_prof - smin_no3_to_plant_vr) - smin_nh4_to_plant_vr  ! use old values
                     smin_nh4_to_plant_vr = plant_ndemand*nuptake_prof - smin_no3_to_plant_vr
                  end if
                  sminn_to_plant_vr = smin_no3_to_plant_vr + smin_nh4_to_plant_vr
               end if
  end subroutine apply_carbon_only_adjustment

  !-----------------------------------------------------------------------
  pure subroutine compute_competition_summary( &
       fpi_vr, sminn_to_plant_vr, actual_immob_vr, &
       fpi_no3_vr, fpi_nh4_vr, &
       smin_no3_to_plant_vr, smin_nh4_to_plant_vr, &
       actual_immob_no3_vr, actual_immob_nh4_vr)
    real(r8), intent(out) :: fpi_vr, sminn_to_plant_vr, actual_immob_vr
    real(r8), intent(in)  :: fpi_no3_vr, fpi_nh4_vr
    real(r8), intent(in)  :: smin_no3_to_plant_vr, smin_nh4_to_plant_vr
    real(r8), intent(in)  :: actual_immob_no3_vr, actual_immob_nh4_vr
               fpi_vr = fpi_no3_vr + fpi_nh4_vr
               sminn_to_plant_vr = smin_no3_to_plant_vr + smin_nh4_to_plant_vr
               actual_immob_vr = actual_immob_no3_vr + actual_immob_nh4_vr
  end subroutine compute_competition_summary

  !-----------------------------------------------------------------------
  ! Generic per-layer dzsoi-weighted accumulation: column_total += value_vr * dz.
  ! Used to vertically integrate sminn_to_plant, actual_immob, potential_immob.
  pure subroutine accum_dz_weighted(column_total, value_vr, dzsoi_decomp)
    real(r8), intent(inout) :: column_total
    real(r8), intent(in)    :: value_vr, dzsoi_decomp
    column_total = column_total + value_vr * dzsoi_decomp
  end subroutine accum_dz_weighted

  !-----------------------------------------------------------------------
  ! Per-layer leftover mineral N after first-pass demands (used for both
  ! NH4 and NO3). f_loss is f_nit_vr for NH4, f_denit_vr for NO3.
  pure function compute_residual_smin_vr( &
       smin_vr, actual_immob_vr, smin_to_plant_vr, f_loss_vr, dt) result(residual_smin_vr)
    real(r8) :: residual_smin_vr
    real(r8), intent(in) :: smin_vr, actual_immob_vr, smin_to_plant_vr, f_loss_vr, dt
    residual_smin_vr = max(smin_vr - (actual_immob_vr + smin_to_plant_vr + f_loss_vr ) * dt, 0._r8)
  end function compute_residual_smin_vr

  !-----------------------------------------------------------------------
  ! Distribute layer-wise residual N to satisfy residual plant demand
  ! (used for both NH4 and NO3).
  pure function distribute_residual_to_plant( &
       smin_to_plant_vr, residual_smin_vr, residual_plant_ndemand, residual_smin, dt) result(smin_to_plant_vr_new)
    real(r8) :: smin_to_plant_vr_new
    real(r8), intent(in) :: smin_to_plant_vr, residual_smin_vr, residual_plant_ndemand, residual_smin, dt
    smin_to_plant_vr_new = smin_to_plant_vr + residual_smin_vr * &
         min(( residual_plant_ndemand *  dt ) / residual_smin, 1._r8) / dt
  end function distribute_residual_to_plant

  !-----------------------------------------------------------------------
  ! Defensive fraction: numerator/denominator if denominator > 0, else 1.
  ! Used for fpg (sminn_to_plant / plant_ndemand) and fpi (actual_immob /
  ! potential_immob) — both naturally return 1 when there's no demand.
  pure function compute_fraction_or_one(numerator, denominator) result(frac)
    real(r8) :: frac
    real(r8), intent(in) :: numerator, denominator
    if (denominator > 0.0_r8) then
       frac = numerator / denominator
    else
       frac = 1._r8
    end if
  end function compute_fraction_or_one

end module SoilBiogeochemCompetition_mod
