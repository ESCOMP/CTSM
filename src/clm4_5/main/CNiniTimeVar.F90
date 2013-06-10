!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNiniTimeVar
!
! !INTERFACE:
subroutine CNiniTimeVar()

#ifdef CN
!
! !DESCRIPTION:
! Initializes time varying variables used only in
! coupled carbon-nitrogen mode (CN):
!
! !USES:
   use clmtype
   use clm_atmlnd  , only: clm_a2l
   use shr_kind_mod, only: r8 => shr_kind_r8
   use clm_varcon  , only: istsoil, zsoi
   use clm_varpar  , only: nlevgrnd, nlevsoi, nlevdecomp, ndecomp_pools, nlevdecomp_full
   use clm_varcon  , only: istcrop, c13ratio, c14ratio
   use clm_varctl  , only: use_c13, use_c14
   use pftvarcon   , only: noveg
   use pftvarcon   , only: npcropmin
   use decompMod   , only: get_proc_bounds
   use surfrdMod   , only: crop_prog
!
! !ARGUMENTS:
   implicit none
!
! !CALLED FROM:
! subroutine iniTimeVar in file iniTimeVar.F90
!
! !REVISION HISTORY:
! 10/21/03: Created by Peter Thornton
! F. Li and S. Levis (11/06/12)
!
! local pointers to implicit in arguments
!
   real(r8), pointer :: evergreen(:) ! binary flag for evergreen leaf habit (0 or 1)
   real(r8), pointer :: woody(:)     ! binary flag for woody lifeform (1=woody, 0=not woody)
   real(r8), pointer :: leafcn(:)    ! leaf C:N (gC/gN)
   real(r8), pointer :: deadwdcn(:)  ! dead wood (xylem and heartwood) C:N (gC/gN)
   integer , pointer :: ivt(:)       ! pft vegetation type
   logical , pointer :: lakpoi(:)    ! true => landunit is a lake point
   integer , pointer :: plandunit(:) ! landunit index associated with each pft
   integer , pointer :: clandunit(:) ! landunit index associated with each column
   integer , pointer :: itypelun(:)  ! landunit type
!
! local pointers to implicit out arguments
!
   real(r8), pointer :: forc_hgt_u_pft(:)    !observational height of wind at pft-level [m]
   real(r8), pointer :: annsum_counter(:) ! seconds since last annual accumulator turnover
   real(r8), pointer :: cannsum_npp(:)    ! annual sum of NPP, averaged from pft-level (gC/m2/yr)
   real(r8), pointer :: cannavg_t2m(:)    !annual average of 2m air temperature, averaged from pft-level (K)
   real(r8), pointer :: sminn(:)              ! (gN/m2) soil mineral N
   real(r8), pointer :: cwdc(:)               ! (gC/m2) coarse woody debris C
   real(r8), pointer :: cwdn(:)               ! (gN/m2) coarse woody debris N
   real(r8), pointer :: decomp_cpools(:,:)         ! (gC/m2)  decomposing (litter, cwd, soil) c pools
   real(r8), pointer :: decomp_cpools_vr(:,:,:)    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   real(r8), pointer :: decomp_cpools_1m(:,:)           ! (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter
   real(r8), pointer :: decomp_npools(:,:)         ! (gC/m2)  decomposing (litter, cwd, soil) N pools
   real(r8), pointer :: decomp_npools_vr(:,:,:)    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
   real(r8), pointer :: decomp_npools_1m(:,:)           ! (gN/m2)  diagnostic: decomposing (litter, cwd, soil) N pools to 1 meter
   real(r8), pointer :: sminn_vr(:,:)              ! (gN/m3) soil mineral N
   real(r8), pointer :: col_ctrunc(:)              ! (gC/m2) column-level sink for C truncation (diagnostic)
   real(r8), pointer :: col_ctrunc_vr(:,:)         ! (gC/m3) column-level sink for C truncation (prognostic)
   real(r8), pointer :: col_ntrunc_vr(:,:)         ! (gN/m3) column-level sink for N truncation
   real(r8), pointer :: nfixation_prof(:,:)     ! (1/m) profile for N fixation additions
   real(r8), pointer :: ndep_prof(:,:)     ! (1/m) profile for N fixation additions

   real(r8), pointer :: fpi_vr(:,:)
   real(r8), pointer :: alt(:)
   real(r8), pointer :: altmax(:)
   real(r8), pointer :: altmax_lastyear(:)
   integer, pointer :: alt_indx(:)
   integer, pointer :: altmax_indx(:)
   integer, pointer :: altmax_lastyear_indx(:)
   real(r8), pointer :: som_adv_coef(:,:)
   real(r8), pointer :: som_diffus_coef(:,:)
#ifdef NITRIF_DENITRIF
   real(r8), pointer :: smin_nh4_vr(:,:)           ! (gN/m3) soil mineral NH4 pool
   real(r8), pointer :: smin_no3_vr(:,:)           ! (gN/m3) soil mineral NO3 pool
   real(r8), pointer :: smin_nh4(:)           ! (gN/m2) soil mineral NH4 pool
   real(r8), pointer :: smin_no3(:)           ! (gN/m2) soil mineral NO3 pool
#endif

   real(r8), pointer :: leafc(:)              ! (gC/m2) leaf C
   real(r8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real(r8), pointer :: grainc(:)             ! (gC/m2) grain C
   real(r8), pointer :: grainc_storage(:)     ! (gC/m2) grain C storage
   real(r8), pointer :: grainc_xfer(:)        ! (gC/m2) grain C transfer
   real(r8), pointer :: frootc(:)             ! (gC/m2) fine root C
   real(r8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real(r8), pointer :: livestemc(:)          ! (gC/m2) live stem C
   real(r8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
   real(r8), pointer :: deadstemc(:)          ! (gC/m2) dead stem C
   real(r8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real(r8), pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real(r8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real(r8), pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
   real(r8), pointer :: cpool(:)              ! (gC/m2) temporary photosynthate C pool
   real(r8), pointer :: xsmrpool(:)           ! (gC/m2) abstract C pool to meet excess MR demand
   real(r8), pointer :: leafn(:)              ! (gN/m2) leaf N
   real(r8), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real(r8), pointer :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
   real(r8), pointer :: grainn(:)             ! (gN/m2) grain N
   real(r8), pointer :: grainn_storage(:)     ! (gN/m2) grain N storage
   real(r8), pointer :: grainn_xfer(:)        ! (gN/m2) grain N transfer
   real(r8), pointer :: frootn(:)             ! (gN/m2) fine root N
   real(r8), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real(r8), pointer :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
   real(r8), pointer :: livestemn(:)          ! (gN/m2) live stem N
   real(r8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real(r8), pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
   real(r8), pointer :: deadstemn(:)          ! (gN/m2) dead stem N
   real(r8), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real(r8), pointer :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
   real(r8), pointer :: livecrootn(:)         ! (gN/m2) live coarse root N
   real(r8), pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real(r8), pointer :: livecrootn_xfer(:)    ! (gN/m2) live coarse root N transfer
   real(r8), pointer :: deadcrootn(:)         ! (gN/m2) dead coarse root N
   real(r8), pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   real(r8), pointer :: deadcrootn_xfer(:)    ! (gN/m2) dead coarse root N transfer
   real(r8), pointer :: retransn(:)           ! (gN/m2) plant pool of retranslocated N
   real(r8), pointer :: npool(:)              ! (gN/m2) temporary plant N pool
   real(r8), pointer :: psnsun(:)             ! sunlit leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: psnsha(:)             ! shaded leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: c13_psnsun(:)             ! sunlit leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: c13_psnsha(:)             ! shaded leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: c14_psnsun(:)             ! sunlit leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: c14_psnsha(:)             ! shaded leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: laisun(:)             ! sunlit projected leaf area index
   real(r8), pointer :: laisha(:)             ! shaded projected leaf area index
   real(r8), pointer :: dormant_flag(:)       ! dormancy flag
   real(r8), pointer :: days_active(:)        ! number of days since last dormancy
   real(r8), pointer :: onset_flag(:)         ! onset flag
   real(r8), pointer :: onset_counter(:)      ! onset days counter
   real(r8), pointer :: onset_gddflag(:)      ! onset flag for growing degree day sum
   real(r8), pointer :: onset_fdd(:)          ! onset freezing degree days counter
   real(r8), pointer :: onset_gdd(:)          ! onset growing degree days
   real(r8), pointer :: onset_swi(:)          ! onset soil water index
   real(r8), pointer :: offset_flag(:)        ! offset flag
   real(r8), pointer :: offset_counter(:)     ! offset days counter
   real(r8), pointer :: offset_fdd(:)         ! offset freezing degree days counter
   real(r8), pointer :: offset_swi(:)         ! offset soil water index
   real(r8), pointer :: lgsf(:)               ! long growing season factor [0-1]
   real(r8), pointer :: bglfr(:)              ! background litterfall rate (1/s)
   real(r8), pointer :: bgtr(:)               ! background transfer rate (1/s)
   real(r8), pointer :: dayl(:)               ! daylength (seconds)
   real(r8), pointer :: prev_dayl(:)          ! daylength from previous timestep (seconds)
   real(r8), pointer :: annavg_t2m(:)         ! annual average 2m air temperature (K)
   real(r8), pointer :: tempavg_t2m(:)        ! temporary average 2m air temperature (K)
   real(r8), pointer :: gpp(:)                ! GPP flux before downregulation (gC/m2/s)
   real(r8), pointer :: availc(:)             ! C flux available for allocation (gC/m2/s)
   real(r8), pointer :: xsmrpool_recover(:)   ! C flux assigned to recovery of negative cpool (gC/m2/s)
   real(r8), pointer :: xsmrpool_c13ratio(:)  ! C flux assigned to recovery of negative cpool (gC/m2/s)
   real(r8), pointer :: alloc_pnow(:)         ! fraction of current allocation to display as new growth (DIM)
   real(r8), pointer :: c_allometry(:)        ! C allocation index (DIM)
   real(r8), pointer :: n_allometry(:)        ! N allocation index (DIM)
   real(r8), pointer :: plant_ndemand(:)      ! N flux required to support initial GPP (gN/m2/s)
   real(r8), pointer :: tempsum_potential_gpp(:) ! temporary annual sum of plant_ndemand
   real(r8), pointer :: annsum_potential_gpp(:)  ! annual sum of plant_ndemand
   real(r8), pointer :: tempmax_retransn(:)   ! temporary max of retranslocated N pool (gN/m2)
   real(r8), pointer :: annmax_retransn(:)    ! annual max of retranslocated N pool (gN/m2)
   real(r8), pointer :: avail_retransn(:)     ! N flux available from retranslocation pool (gN/m2/s)
   real(r8), pointer :: plant_nalloc(:)       ! total allocated N flux (gN/m2/s)
   real(r8), pointer :: plant_calloc(:)       ! total allocated C flux (gC/m2/s)
   real(r8), pointer :: excess_cflux(:)       ! C flux not allocated due to downregulation (gC/m2/s)
   real(r8), pointer :: downreg(:)            ! fractional reduction in GPP due to N limitation (DIM)
   real(r8), pointer :: tempsum_npp(:)        ! temporary annual sum of NPP
   real(r8), pointer :: annsum_npp(:)         ! annual sum of NPP
#if (defined CNDV)
   real(r8), pointer :: tempsum_litfall(:)    ! temporary annual sum of litfall
   real(r8), pointer :: annsum_litfall(:)     ! annual sum of litfall
#endif
   real(r8), pointer :: rc13_canair(:)        !C13O2/C12O2 in canopy air
   real(r8), pointer :: rc13_psnsun(:)        !C13O2/C12O2 in sunlit canopy psn flux
   real(r8), pointer :: rc13_psnsha(:)        !C13O2/C12O2 in shaded canopy psn flux
   real(r8), pointer :: alphapsnsun(:)        !sunlit 13c fractionation ([])
   real(r8), pointer :: alphapsnsha(:)        !shaded 13c fractionation ([])
   real(r8), pointer :: qflx_drain(:)         ! sub-surface runoff (mm H2O /s)
   real(r8), pointer :: qflx_surf(:)          ! surface runoff (mm H2O /s)
   real(r8), pointer :: qflx_irrig(:)         !irrigation flux (mm H2O/s)

   ! fire-related variables changed by F. Li and S. Levis
   real(r8), pointer :: wf(:)                 ! soil moisture in top 0.05 m
   real(r8), pointer :: wf2(:) 
   real(r8), pointer :: nfire(:)              ! fire counts/km2/timestep
   real(r8), pointer :: baf_crop(:)          ! burned area fraction in crop
   real(r8), pointer :: baf_peatf(:)         ! burned area fraction in peatland
   real(r8), pointer :: fbac(:)          
   real(r8), pointer :: fbac1(:)          
   real(r8), pointer :: farea_burned(:)       ! timestep fractional area burned (proportion)

   real(r8), pointer :: totcolc(:)            ! (gC/m2) total column carbon, incl veg and cpool
   real(r8), pointer :: totecosysc(:)         ! (gC/m2) total ecosystem carbon, incl veg but excl cpool
   real(r8), pointer :: totlitc(:)            ! (gC/m2) total litter carbon
   real(r8), pointer :: totsomc(:)            ! (gC/m2) total soil organic matter carbon
   real(r8), pointer :: totlitc_1m(:)         ! (gC/m2) total litter carbon to 1 meter
   real(r8), pointer :: totsomc_1m(:)         ! (gC/m2) total soil organic matter carbon to 1 meter

   real(r8), pointer :: woodc(:)              ! (gC/m2) pft-level wood C
   real(r8), pointer :: totcoln(:)            ! (gN/m2) total column nitrogen, incl veg
   real(r8), pointer :: totecosysn(:)         ! (gN/m2) total ecosystem nitrogen, incl veg
   real(r8), pointer :: totlitn(:)            ! (gN/m2) total litter nitrogen
   real(r8), pointer :: totsomn(:)            ! (gN/m2) total soil organic matter nitrogen
   real(r8), pointer :: totlitn_1m(:)         ! (gN/m2) total litter nitrogen to 1 meter
   real(r8), pointer :: totsomn_1m(:)         ! (gN/m2) total soil organic matter nitrogen to 1 meter
   real(r8), pointer :: dispvegc(:)           ! (gC/m2) displayed veg carbon, excluding storage and cpool
   real(r8), pointer :: pft_ctrunc(:)         ! (gC/m2) pft-level sink for C truncation
   real(r8), pointer :: storvegc(:)           ! (gC/m2) stored vegetation carbon, excluding cpool
   real(r8), pointer :: totpftc(:)            ! (gC/m2) total pft-level carbon, including cpool
   real(r8), pointer :: totvegc(:)            ! (gC/m2) total vegetation carbon, excluding cpool
   real(r8), pointer :: prev_frootc_to_litter(:)!previous timestep froot C litterfall flux (gC/m2/s)
   real(r8), pointer :: prev_leafc_to_litter(:) !previous timestep leaf C litterfall flux (gC/m2/s)
   real(r8), pointer :: dispvegn(:)           ! (gN/m2) displayed veg nitrogen, excluding storage
   real(r8), pointer :: pft_ntrunc(:)         ! (gN/m2) pft-level sink for N truncation
   real(r8), pointer :: storvegn(:)           ! (gN/m2) stored vegetation nitrogen
   real(r8), pointer :: totpftn(:)            ! (gN/m2) total pft-level nitrogen
   real(r8), pointer :: totvegn(:)            ! (gN/m2) total vegetation nitrogen
   !!! C13
   real(r8), pointer :: cwdc13(:)               ! (gC/m2) coarse woody debris C
   real(r8), pointer :: decomp_c13pools(:,:)              ! (gC/m2)  decomposing (litter, cwd, soil) c pools
   real(r8), pointer :: decomp_c13pools_vr(:,:,:)         ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   real(r8), pointer :: decomp_c13pools_1m(:,:)           ! (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter
   real(r8), pointer :: c13_col_ctrunc_vr(:,:)            ! (gC/m3) C truncation term
   real(r8), pointer :: leafc13(:)              ! (gC/m2) leaf C
   real(r8), pointer :: leafc13_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: leafc13_xfer(:)         ! (gC/m2) leaf C transfer
   real(r8), pointer :: frootc13(:)             ! (gC/m2) fine root C
   real(r8), pointer :: frootc13_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: frootc13_xfer(:)        ! (gC/m2) fine root C transfer
   real(r8), pointer :: livestemc13(:)          ! (gC/m2) live stem C
   real(r8), pointer :: livestemc13_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: livestemc13_xfer(:)     ! (gC/m2) live stem C transfer
   real(r8), pointer :: deadstemc13(:)          ! (gC/m2) dead stem C
   real(r8), pointer :: deadstemc13_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: deadstemc13_xfer(:)     ! (gC/m2) dead stem C transfer
   real(r8), pointer :: livecrootc13(:)         ! (gC/m2) live coarse root C
   real(r8), pointer :: livecrootc13_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: livecrootc13_xfer(:)    ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: deadcrootc13(:)         ! (gC/m2) dead coarse root C
   real(r8), pointer :: deadcrootc13_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: deadcrootc13_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: c13_gresp_storage(:)    ! (gC/m2) growth respiration storage
   real(r8), pointer :: c13_gresp_xfer(:)       ! (gC/m2) growth respiration transfer
   real(r8), pointer :: c13pool(:)              ! (gC/m2) temporary photosynthate C pool
   real(r8), pointer :: c13xsmrpool(:)          ! (gC/m2) temporary photosynthate C pool
   real(r8), pointer :: c13_pft_ctrunc(:)       ! (gC/m2) C truncation term
   real(r8), pointer :: totvegc13(:)            ! (gC/m2) total vegetation carbon, excluding cpool

   !!! C14
   real(r8), pointer :: cwdc14(:)               ! (gC/m2) coarse woody debris C
   real(r8), pointer :: decomp_c14pools(:,:)              ! (gC/m2)  decomposing (litter, cwd, soil) c pools
   real(r8), pointer :: decomp_c14pools_vr(:,:,:)         ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   real(r8), pointer :: decomp_c14pools_1m(:,:)           ! (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter
   real(r8), pointer :: c14_col_ctrunc_vr(:,:)            ! (gC/m3) C truncation term
   real(r8), pointer :: leafc14(:)              ! (gC/m2) leaf C
   real(r8), pointer :: leafc14_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: leafc14_xfer(:)         ! (gC/m2) leaf C transfer
   real(r8), pointer :: frootc14(:)             ! (gC/m2) fine root C
   real(r8), pointer :: frootc14_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: frootc14_xfer(:)        ! (gC/m2) fine root C transfer
   real(r8), pointer :: livestemc14(:)          ! (gC/m2) live stem C
   real(r8), pointer :: livestemc14_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: livestemc14_xfer(:)     ! (gC/m2) live stem C transfer
   real(r8), pointer :: deadstemc14(:)          ! (gC/m2) dead stem C
   real(r8), pointer :: deadstemc14_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: deadstemc14_xfer(:)     ! (gC/m2) dead stem C transfer
   real(r8), pointer :: livecrootc14(:)         ! (gC/m2) live coarse root C
   real(r8), pointer :: livecrootc14_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: livecrootc14_xfer(:)    ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: deadcrootc14(:)         ! (gC/m2) dead coarse root C
   real(r8), pointer :: deadcrootc14_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: deadcrootc14_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: c14_gresp_storage(:)    ! (gC/m2) growth respiration storage
   real(r8), pointer :: c14_gresp_xfer(:)       ! (gC/m2) growth respiration transfer
   real(r8), pointer :: c14pool(:)              ! (gC/m2) temporary photosynthate C pool
   real(r8), pointer :: c14xsmrpool(:)          ! (gC/m2) temporary photosynthate C pool
   real(r8), pointer :: c14_pft_ctrunc(:)       ! (gC/m2) C truncation term
   real(r8), pointer :: totvegc14(:)            ! (gC/m2) total vegetation carbon, excluding cpool
   real(r8), pointer :: rc14_atm(:)             !C14O2/C12O2 in atmosphere

   ! dynamic landuse variables
   real(r8), pointer :: seedc(:)              ! (gC/m2) column-level pool for seeding new PFTs
   real(r8), pointer :: prod10c(:)            ! (gC/m2) wood product C pool, 10-year lifespan
   real(r8), pointer :: prod100c(:)           ! (gC/m2) wood product C pool, 100-year lifespan
   real(r8), pointer :: totprodc(:)           ! (gC/m2) total wood product C

   !!! C13
   real(r8), pointer :: seedc13(:)              ! (gC/m2) column-level pool for seeding new PFTs
   real(r8), pointer :: prod10c13(:)          ! (gC/m2) wood product C13 pool, 10-year lifespan
   real(r8), pointer :: prod100c13(:)         ! (gC/m2) wood product C13 pool, 100-year lifespan
   real(r8), pointer :: totprodc13(:)         ! (gC/m2) total wood product C13
   !!! C14
   real(r8), pointer :: seedc14(:)              ! (gC/m2) column-level pool for seeding new PFTs
   real(r8), pointer :: prod10c14(:)          ! (gC/m2) wood product C14 pool, 10-year lifespan
   real(r8), pointer :: prod100c14(:)         ! (gC/m2) wood product C14 pool, 100-year lifespan
   real(r8), pointer :: totprodc14(:)         ! (gC/m2) total wood product C14

   real(r8), pointer :: seedn(:)              ! (gN/m2) column-level pool for seeding new PFTs
   real(r8), pointer :: prod10n(:)            ! (gN/m2) wood product N pool, 10-year lifespan
   real(r8), pointer :: prod100n(:)           ! (gN/m2) wood product N pool, 100-year lifespan
   real(r8), pointer :: totprodn(:)           ! (gN/m2) total wood product N
   real(r8), pointer :: initial_cn_ratio(:)               ! c:n ratio for initialization of pools
   real(r8), pointer :: initial_stock(:)                  ! initial concentration for seeding at spinup

   ! crop
   real(r8), pointer :: fert_counter(:)
   real(r8), pointer :: fert(:)
   real(r8), pointer :: soyfixn(:)
   real(r8), pointer :: grain_flag(:)
!
! !LOCAL VARIABLES:
   integer :: g,l,c,p,j,k      ! indices
   integer :: begp, endp   ! per-clump/proc beginning and ending pft indices
   integer :: begc, endc   ! per-clump/proc beginning and ending column indices
   integer :: begl, endl   ! per-clump/proc beginning and ending landunit indices
   integer :: begg, endg   ! per-clump/proc gridcell ending gridcell indices

!EOP
!-----------------------------------------------------------------------

    ! assign local pointers at the gridcell level

    ! assign local pointers at the landunit level
    lakpoi                         =>lun%lakpoi
    itypelun                       => lun%itype

    ! assign local pointers at the column level
    clandunit                      =>col%landunit
    annsum_counter                 => cps%annsum_counter
    cannsum_npp                    => cps%cannsum_npp
    cannavg_t2m                    => cps%cannavg_t2m

     !fire related variables changed by F. Li and S. Levis
    wf                             => cps%wf
    wf2                            => cps%wf2
    nfire                          => cps%nfire
    baf_crop                       => cps%baf_crop
    baf_peatf                      => cps%baf_peatf
    fbac                           => cps%fbac
    fbac1                          => cps%fbac1
    farea_burned                   => cps%farea_burned

    qflx_drain                     => cwf%qflx_drain
    qflx_surf                      => cwf%qflx_surf
    decomp_cpools                  => ccs%decomp_cpools
    decomp_cpools_1m               => ccs%decomp_cpools_1m
    decomp_cpools_vr               => ccs%decomp_cpools_vr
    decomp_npools                  => cns%decomp_npools
    decomp_npools_vr               => cns%decomp_npools_vr
    decomp_npools_1m               => cns%decomp_npools_1m
    nfixation_prof                    => cps%nfixation_prof
    ndep_prof                         => cps%ndep_prof
    qflx_irrig                     => cwf%qflx_irrig
    
    ! dynamic landuse variables
    seedc                          => ccs%seedc
    prod10c                        => ccs%prod10c
    prod100c                       => ccs%prod100c
    totprodc                       => ccs%totprodc
    seedn                          => cns%seedn
    prod10n                        => cns%prod10n
    prod100n                       => cns%prod100n
    totprodn                       => cns%totprodn
    sminn                          => cns%sminn
    col_ctrunc                        => ccs%col_ctrunc    
    sminn_vr                          => cns%sminn_vr
    col_ctrunc_vr                     => ccs%col_ctrunc_vr
    col_ntrunc_vr                     => cns%col_ntrunc_vr
    
    fpi_vr                                 => cps%fpi_vr
    alt                                    => cps%alt
    altmax                                 => cps%altmax
    altmax_lastyear                        => cps%altmax_lastyear
    som_adv_coef                           => cps%som_adv_coef
    som_diffus_coef                        => cps%som_diffus_coef
    alt_indx                               => cps%alt_indx
    altmax_indx                            => cps%altmax_indx
    altmax_lastyear_indx                   => cps%altmax_lastyear_indx
#ifdef NITRIF_DENITRIF
    smin_nh4_vr                       => cns%smin_nh4_vr
    smin_no3_vr                       => cns%smin_no3_vr
    smin_nh4                       => cns%smin_nh4
    smin_no3                       => cns%smin_no3
#endif
    
    totcolc                        => ccs%totcolc
    cwdc                           => ccs%cwdc
    totecosysc                     => ccs%totecosysc
    totlitc                        => ccs%totlitc
    totsomc                        => ccs%totsomc
    totlitc_1m                     => ccs%totlitc_1m
    totsomc_1m                     => ccs%totsomc_1m
    
    totcoln                        => cns%totcoln
    cwdn                           => cns%cwdn
    totecosysn                     => cns%totecosysn
    totlitn                        => cns%totlitn
    totsomn                        => cns%totsomn
    totlitn_1m                     => cns%totlitn_1m
    totsomn_1m                     => cns%totsomn_1m
    if ( use_c13 ) then
       seedc13                          => cc13s%seedc
       prod10c13                        => cc13s%prod10c
       prod100c13                       => cc13s%prod100c
       totprodc13                       => cc13s%totprodc
       cwdc13                           => cc13s%cwdc
       decomp_c13pools                  => cc13s%decomp_cpools
       decomp_c13pools_vr               => cc13s%decomp_cpools_vr
       c13_col_ctrunc_vr                => cc13s%col_ctrunc_vr
       decomp_c13pools_1m               => cc13s%decomp_cpools_1m
       c13_psnsun                       => pc13f%psnsun
       c13_psnsha                       => pc13f%psnsha
       xsmrpool_c13ratio                => pepv%xsmrpool_c13ratio
       alphapsnsun                      => pps%alphapsnsun
       alphapsnsha                      => pps%alphapsnsha
       leafc13                          => pc13s%leafc
       leafc13_storage                  => pc13s%leafc_storage
       leafc13_xfer                     => pc13s%leafc_xfer
       frootc13                         => pc13s%frootc
       frootc13_storage                 => pc13s%frootc_storage
       frootc13_xfer                    => pc13s%frootc_xfer
       livestemc13                      => pc13s%livestemc
       livestemc13_storage              => pc13s%livestemc_storage
       livestemc13_xfer                 => pc13s%livestemc_xfer
       deadstemc13                      => pc13s%deadstemc
       deadstemc13_storage              => pc13s%deadstemc_storage
       deadstemc13_xfer                 => pc13s%deadstemc_xfer
       livecrootc13                     => pc13s%livecrootc
       livecrootc13_storage             => pc13s%livecrootc_storage
       livecrootc13_xfer                => pc13s%livecrootc_xfer
       deadcrootc13                     => pc13s%deadcrootc
       deadcrootc13_storage             => pc13s%deadcrootc_storage
       deadcrootc13_xfer                => pc13s%deadcrootc_xfer
       c13_gresp_storage                => pc13s%gresp_storage
       c13_gresp_xfer                   => pc13s%gresp_xfer
       c13pool                          => pc13s%cpool
       c13xsmrpool                      => pc13s%xsmrpool
       c13_pft_ctrunc                   => pc13s%pft_ctrunc
       totvegc13                        => pc13s%totvegc
       rc13_canair                      => pepv%rc13_canair
       rc13_psnsun                      => pepv%rc13_psnsun
       rc13_psnsha                      => pepv%rc13_psnsha
    endif
    if ( use_c14 ) then
       seedc14                          => cc14s%seedc
       prod10c14                        => cc14s%prod10c
       prod100c14                       => cc14s%prod100c
       totprodc14                       => cc14s%totprodc
       cwdc14                           => cc14s%cwdc
       decomp_c14pools                  => cc14s%decomp_cpools
       decomp_c14pools_vr               => cc14s%decomp_cpools_vr
       c14_col_ctrunc_vr                => cc14s%col_ctrunc_vr
       decomp_c14pools_1m               => cc14s%decomp_cpools_1m
       c14_psnsun                       => pc14f%psnsun
       c14_psnsha                       => pc14f%psnsha
       leafc14                          => pc14s%leafc
       leafc14_storage                  => pc14s%leafc_storage
       leafc14_xfer                     => pc14s%leafc_xfer
       frootc14                         => pc14s%frootc
       frootc14_storage                 => pc14s%frootc_storage
       frootc14_xfer                    => pc14s%frootc_xfer
       livestemc14                      => pc14s%livestemc
       livestemc14_storage              => pc14s%livestemc_storage
       livestemc14_xfer                 => pc14s%livestemc_xfer
       deadstemc14                      => pc14s%deadstemc
       deadstemc14_storage              => pc14s%deadstemc_storage
       deadstemc14_xfer                 => pc14s%deadstemc_xfer
       livecrootc14                     => pc14s%livecrootc
       livecrootc14_storage             => pc14s%livecrootc_storage
       livecrootc14_xfer                => pc14s%livecrootc_xfer
       deadcrootc14                     => pc14s%deadcrootc
       deadcrootc14_storage             => pc14s%deadcrootc_storage
       deadcrootc14_xfer                => pc14s%deadcrootc_xfer
       c14_gresp_storage                => pc14s%gresp_storage
       c14_gresp_xfer                   => pc14s%gresp_xfer
       c14pool                          => pc14s%cpool
       c14xsmrpool                      => pc14s%xsmrpool
       c14_pft_ctrunc                   => pc14s%pft_ctrunc
       totvegc14                        => pc14s%totvegc
       rc14_atm                         => pepv%rc14_atm
    endif
    ! crop
    soyfixn                          => pnf%soyfixn
    fert                             => pnf%fert
    fert_counter                     => pepv%fert_counter
    grain_flag                       => pepv%grain_flag
    ! assign local pointers at the pft level
    ivt                            =>pft%itype
    plandunit                      =>pft%landunit
    leafc                          => pcs%leafc
    leafc_storage                  => pcs%leafc_storage
    leafc_xfer                     => pcs%leafc_xfer
    grainc                         => pcs%grainc
    grainc_storage                 => pcs%grainc_storage
    grainc_xfer                    => pcs%grainc_xfer
    frootc                         => pcs%frootc
    frootc_storage                 => pcs%frootc_storage
    frootc_xfer                    => pcs%frootc_xfer
    livestemc                      => pcs%livestemc
    livestemc_storage              => pcs%livestemc_storage
    livestemc_xfer                 => pcs%livestemc_xfer
    deadstemc                      => pcs%deadstemc
    deadstemc_storage              => pcs%deadstemc_storage
    deadstemc_xfer                 => pcs%deadstemc_xfer
    livecrootc                     => pcs%livecrootc
    livecrootc_storage             => pcs%livecrootc_storage
    livecrootc_xfer                => pcs%livecrootc_xfer
    deadcrootc                     => pcs%deadcrootc
    deadcrootc_storage             => pcs%deadcrootc_storage
    deadcrootc_xfer                => pcs%deadcrootc_xfer
    gresp_storage                  => pcs%gresp_storage
    gresp_xfer                     => pcs%gresp_xfer
    cpool                          => pcs%cpool
    xsmrpool                       => pcs%xsmrpool
    forc_hgt_u_pft                 => pps%forc_hgt_u_pft
    woodc                          => pcs%woodc
    leafn                          => pns%leafn
    leafn_storage                  => pns%leafn_storage
    leafn_xfer                     => pns%leafn_xfer
    grainn                         => pns%grainn
    grainn_storage                 => pns%grainn_storage
    grainn_xfer                    => pns%grainn_xfer
    frootn                         => pns%frootn
    frootn_storage                 => pns%frootn_storage
    frootn_xfer                    => pns%frootn_xfer
    livestemn                      => pns%livestemn
    livestemn_storage              => pns%livestemn_storage
    livestemn_xfer                 => pns%livestemn_xfer
    deadstemn                      => pns%deadstemn
    deadstemn_storage              => pns%deadstemn_storage
    deadstemn_xfer                 => pns%deadstemn_xfer
    livecrootn                     => pns%livecrootn
    livecrootn_storage             => pns%livecrootn_storage
    livecrootn_xfer                => pns%livecrootn_xfer
    deadcrootn                     => pns%deadcrootn
    deadcrootn_storage             => pns%deadcrootn_storage
    deadcrootn_xfer                => pns%deadcrootn_xfer
    retransn                       => pns%retransn
    npool                          => pns%npool
    psnsun                         => pcf%psnsun
    psnsha                         => pcf%psnsha
    laisun                         => pps%laisun
    laisha                         => pps%laisha
    dormant_flag                   => pepv%dormant_flag
    days_active                    => pepv%days_active
    onset_flag                     => pepv%onset_flag
    onset_counter                  => pepv%onset_counter
    onset_gddflag                  => pepv%onset_gddflag
    onset_fdd                      => pepv%onset_fdd
    onset_gdd                      => pepv%onset_gdd
    onset_swi                      => pepv%onset_swi
    offset_flag                    => pepv%offset_flag
    offset_counter                 => pepv%offset_counter
    offset_fdd                     => pepv%offset_fdd
    offset_swi                     => pepv%offset_swi
    lgsf                           => pepv%lgsf
    bglfr                          => pepv%bglfr
    bgtr                           => pepv%bgtr
    dayl                           => pepv%dayl
    prev_dayl                      => pepv%prev_dayl
    annavg_t2m                     => pepv%annavg_t2m
    tempavg_t2m                    => pepv%tempavg_t2m
    gpp                            => pepv%gpp
    availc                         => pepv%availc
    xsmrpool_recover                  => pepv%xsmrpool_recover
    alloc_pnow                     => pepv%alloc_pnow
    c_allometry                    => pepv%c_allometry
    n_allometry                    => pepv%n_allometry
    plant_ndemand                  => pepv%plant_ndemand
    tempsum_potential_gpp          => pepv%tempsum_potential_gpp
    annsum_potential_gpp           => pepv%annsum_potential_gpp
    tempmax_retransn               => pepv%tempmax_retransn
    annmax_retransn                => pepv%annmax_retransn
    avail_retransn                 => pepv%avail_retransn
    plant_nalloc                   => pepv%plant_nalloc
    plant_calloc                   => pepv%plant_calloc
    excess_cflux                   => pepv%excess_cflux
    downreg                        => pepv%downreg
    tempsum_npp                    => pepv%tempsum_npp
    annsum_npp                     => pepv%annsum_npp
#if (defined CNDV)
    tempsum_litfall                => pepv%tempsum_litfall
    annsum_litfall                 => pepv%annsum_litfall
#endif
    dispvegc                       => pcs%dispvegc
    pft_ctrunc                     => pcs%pft_ctrunc
    storvegc                       => pcs%storvegc
    totpftc                        => pcs%totpftc
    totvegc                        => pcs%totvegc
    prev_frootc_to_litter          => pepv%prev_frootc_to_litter
    prev_leafc_to_litter           => pepv%prev_leafc_to_litter
    dispvegn                       => pns%dispvegn
    pft_ntrunc                     => pns%pft_ntrunc
    storvegn                       => pns%storvegn
    totpftn                        => pns%totpftn
    totvegn                        => pns%totvegn
    
    ! assign local pointers for ecophysiological constants
    evergreen                      => pftcon%evergreen
    woody                          => pftcon%woody
    leafcn                         => pftcon%leafcn
    deadwdcn                       => pftcon%deadwdcn
    ! decomposoition parameters
   initial_cn_ratio                        => decomp_cascade_con%initial_cn_ratio
   initial_stock                           => decomp_cascade_con%initial_stock

   ! Determine subgrid bounds on this processor
   call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

   ! Added 5/4/04, PET: initialize forc_hgt_u (gridcell-level),
   ! since this is not initialized before first call to CNVegStructUpdate,
   ! and it is required to set the upper bound for canopy top height.
   ! Changed 3/21/08, KO: still needed but don't have sufficient information 
   ! to set this properly (e.g., pft-level displacement height and roughness 
   ! length). So leave at 30m.
   do p = begp, endp
      forc_hgt_u_pft(p) = 30._r8
   end do

   ! initialize column-level variables
   do c = begc, endc
      l = clandunit(c)
      if (itypelun(l) == istsoil .or. itypelun(l) == istcrop) then

         ! column physical state variables
         annsum_counter(c) = 0._r8
         cannsum_npp(c)    = 0._r8
         cannavg_t2m(c)    = 280._r8
         ! fire related variables changed by F. Li and S. Levis
         wf(c) = 1.0_r8  ! it needs to be non zero so the first time step has no fires
         wf2(c) = 1.0_r8
         nfire(c) = 0._r8
         baf_crop(c) = 0._r8
         baf_peatf(c) = 0._r8
         fbac(c) = 0._r8
         fbac1(c) = 0._r8
         farea_burned(c) = 0._r8

         
         ! needed for CNNLeaching
         qflx_drain(c) = 0._r8
         qflx_surf(c) = 0._r8
         qflx_irrig(c) = 0._r8

         ! column carbon state variable initialization
         do j = 1, nlevdecomp
            do k = 1, ndecomp_pools
               if (zsoi(j) .lt. 0.3 ) then  !! only initialize upper soil column
                  decomp_cpools_vr(c,j,k) = initial_stock(k)
               else
                  decomp_cpools_vr(c,j,k) = 0._r8
               endif
            end do
            col_ctrunc_vr(c,j) = 0._r8
         end do
         if ( nlevdecomp .gt. 1 ) then
            do j = nlevdecomp+1, nlevdecomp_full
               do k = 1, ndecomp_pools
                  decomp_cpools_vr(c,j,k) = 0._r8
               end do
               col_ctrunc_vr(c,j) = 0._r8
            end do
         end if
         do k = 1, ndecomp_pools
            decomp_cpools(c,k) = initial_stock(k)
            decomp_cpools_1m(c,k) = initial_stock(k)
         end do

         do j = 1, nlevdecomp_full
            ! initialize fpi_vr so that levels below nlevsoi are not nans
            fpi_vr(c,j) = 0._r8
            som_adv_coef(c,j) = 0._r8
            som_diffus_coef(c,j) = 0._r8
            ! here initialize the profiles for converting to vertically resolved carbon pools
            nfixation_prof(c,j) = 0._r8
            ndep_prof(c,j) = 0._r8
         end do

         ! and define alt variables to be zero
         alt(c) = 0._r8
         altmax(c) = 0._r8
         altmax_lastyear(c) = 0._r8
         alt_indx(c) = 0
         altmax_indx(c) = 0
         altmax_lastyear_indx = 0

         cwdc(c)       = 0._r8
         col_ctrunc(c) = 0._r8
         totlitc(c)    = 0._r8
         totsomc(c)    = 0._r8
         totlitc_1m(c) = 0._r8
         totsomc_1m(c) = 0._r8
         totecosysc(c) = 0._r8
         totcolc(c)    = 0._r8

         if ( use_c13 ) then
            do j = 1, nlevdecomp
               do k = 1, ndecomp_pools
                  decomp_c13pools_vr(c,j,k) = decomp_cpools_vr(c,j,k) * c13ratio
               end do
               c13_col_ctrunc_vr(c,j) = col_ctrunc_vr(c,j) * c13ratio
            end do
            if ( nlevdecomp .gt. 1 ) then
               do j = nlevdecomp+1, nlevdecomp_full
                  do k = 1, ndecomp_pools
                     decomp_c13pools_vr(c,j,k) = 0._r8
                  end do
                  c13_col_ctrunc_vr(c,j) = 0._r8
               end do
            end if
            cwdc13(c) = cwdc(c) * c13ratio
            do k = 1, ndecomp_pools
               decomp_c13pools(c,k) = decomp_cpools(c,k) * c13ratio
               decomp_c13pools_1m(c,k) = decomp_cpools_1m(c,k) * c13ratio
            end do
         endif

         if ( use_c14 ) then
            do j = 1, nlevdecomp
               do k = 1, ndecomp_pools
                  decomp_c14pools_vr(c,j,k) = decomp_cpools_vr(c,j,k) * c14ratio
               end do
               c14_col_ctrunc_vr(c,j) = col_ctrunc_vr(c,j) * c14ratio
            end do
            if ( nlevdecomp .gt. 1 ) then
               do j = nlevdecomp+1, nlevdecomp_full
                  do k = 1, ndecomp_pools
                     decomp_c14pools_vr(c,j,k) = 0._r8
                  end do
                  c14_col_ctrunc_vr(c,j) = 0._r8
               end do
            end if
            cwdc14(c) = cwdc(c) * c14ratio
            do k = 1, ndecomp_pools
               decomp_c14pools(c,k) = decomp_cpools(c,k) * c14ratio
               decomp_c14pools_1m(c,k) = decomp_cpools_1m(c,k) * c14ratio
            end do
         endif

         ! column nitrogen state variables
         sminn(c) = 0._r8
         do j = 1, nlevdecomp
            do k = 1, ndecomp_pools
               decomp_npools_vr(c,j,k) = decomp_cpools_vr(c,j,k) / initial_cn_ratio(k)
            end do
            sminn_vr(c,j) = 0._r8
            col_ntrunc_vr(c,j) = 0._r8
         end do
         if ( nlevdecomp .gt. 1 ) then
            do j = nlevdecomp+1, nlevdecomp_full
               do k = 1, ndecomp_pools
                  decomp_npools_vr(c,j,k) = 0._r8
               end do
               sminn_vr(c,j) = 0._r8
               col_ntrunc_vr(c,j) = 0._r8
            end do
         end if
         do k = 1, ndecomp_pools
            decomp_npools(c,k) = decomp_cpools(c,k) / initial_cn_ratio(k)
            decomp_npools_1m(c,k) = decomp_cpools_1m(c,k) / initial_cn_ratio(k)
         end do

#ifdef NITRIF_DENITRIF
         do j = 1, nlevdecomp_full
            smin_nh4_vr(c,j) = 0._r8
            smin_no3_vr(c,j) = 0._r8
         end do
         smin_nh4(c) = 0._r8
         smin_no3(c) = 0._r8
#endif
         totlitn(c)    = 0._r8
         totsomn(c)    = 0._r8
         totlitn_1m(c) = 0._r8
         totsomn_1m(c) = 0._r8
         totecosysn(c) = 0._r8
         totcoln(c)    = 0._r8
         cwdn(c)       = 0._r8

	 ! dynamic landcover state variables
     seedc(c)  = 0._r8
	 prod10c(c)    = 0._r8
	 prod100c(c)   = 0._r8
	 totprodc(c)   = 0._r8

         if ( use_c13 ) then
            seedc13(c)    = 0._r8
            prod10c13(c)  = 0._r8
            prod100c13(c) = 0._r8
            totprodc13(c) = 0._r8
         endif

         if ( use_c14 ) then
            seedc14(c)    = 0._r8
            prod10c14(c)  = 0._r8
            prod100c14(c) = 0._r8
            totprodc14(c) = 0._r8
         endif

	 seedn(c)      = 0._r8
	 prod10n(c)    = 0._r8
	 prod100n(c)   = 0._r8
	 totprodn(c)   = 0._r8
	 
	 ! also initialize dynamic landcover fluxes so that they have
	 ! real values on first timestep, prior to calling pftdyn_cnbal
	 ccf%dwt_seedc_to_leaf(c) = 0._r8
	 ccf%dwt_seedc_to_deadstem(c) = 0._r8
	 ccf%dwt_conv_cflux(c) = 0._r8
        ccf%lf_conv_cflux(c) = 0._r8
	 ccf%dwt_prod10c_gain(c) = 0._r8
	 ccf%prod10c_loss(c) = 0._r8
	 ccf%dwt_prod100c_gain(c) = 0._r8
	 ccf%prod100c_loss(c) = 0._r8
         do j = 1, nlevdecomp_full
            ccf%dwt_frootc_to_litr_met_c(c,j) = 0._r8
            ccf%dwt_frootc_to_litr_cel_c(c,j) = 0._r8
            ccf%dwt_frootc_to_litr_lig_c(c,j) = 0._r8
            ccf%dwt_livecrootc_to_cwdc(c,j) = 0._r8
            ccf%dwt_deadcrootc_to_cwdc(c,j) = 0._r8
          end do
	 ccf%dwt_closs(c) = 0._r8

         if ( use_c13 ) then
            cc13f%dwt_seedc_to_leaf(c) = 0._r8
            cc13f%dwt_seedc_to_deadstem(c) = 0._r8
            cc13f%dwt_conv_cflux(c) = 0._r8
            cc13f%dwt_prod10c_gain(c) = 0._r8
            cc13f%prod10c_loss(c) = 0._r8
            cc13f%dwt_prod100c_gain(c) = 0._r8
            cc13f%prod100c_loss(c) = 0._r8
            do j = 1, nlevdecomp_full
               cc13f%dwt_frootc_to_litr_met_c(c,j) = 0._r8
               cc13f%dwt_frootc_to_litr_cel_c(c,j) = 0._r8
               cc13f%dwt_frootc_to_litr_lig_c(c,j) = 0._r8
               cc13f%dwt_livecrootc_to_cwdc(c,j) = 0._r8
               cc13f%dwt_deadcrootc_to_cwdc(c,j) = 0._r8
            end do
            cc13f%dwt_closs(c) = 0._r8
         endif         
         
         if ( use_c14 ) then
            cc14f%dwt_seedc_to_leaf(c) = 0._r8
            cc14f%dwt_seedc_to_deadstem(c) = 0._r8
            cc14f%dwt_conv_cflux(c) = 0._r8
            cc14f%dwt_prod10c_gain(c) = 0._r8
            cc14f%prod10c_loss(c) = 0._r8
            cc14f%dwt_prod100c_gain(c) = 0._r8
            cc14f%prod100c_loss(c) = 0._r8
            do j = 1, nlevdecomp_full
               cc14f%dwt_frootc_to_litr_met_c(c,j) = 0._r8
               cc14f%dwt_frootc_to_litr_cel_c(c,j) = 0._r8
               cc14f%dwt_frootc_to_litr_lig_c(c,j) = 0._r8
               cc14f%dwt_livecrootc_to_cwdc(c,j) = 0._r8
               cc14f%dwt_deadcrootc_to_cwdc(c,j) = 0._r8
            end do
            cc14f%dwt_closs(c) = 0._r8
         endif

	 cnf%dwt_seedn_to_leaf(c) = 0._r8
	 cnf%dwt_seedn_to_deadstem(c) = 0._r8
	 cnf%dwt_conv_nflux(c) = 0._r8
	 cnf%dwt_prod10n_gain(c) = 0._r8
	 cnf%prod10n_loss(c) = 0._r8
	 cnf%dwt_prod100n_gain(c) = 0._r8
	 cnf%prod100n_loss(c) = 0._r8
         do j = 1, nlevdecomp_full
            cnf%dwt_frootn_to_litr_met_n(c,j) = 0._r8
            cnf%dwt_frootn_to_litr_cel_n(c,j) = 0._r8
            cnf%dwt_frootn_to_litr_lig_n(c,j) = 0._r8
            cnf%dwt_livecrootn_to_cwdn(c,j) = 0._r8
            cnf%dwt_deadcrootn_to_cwdn(c,j) = 0._r8
         end do
	 cnf%dwt_nloss(c) = 0._r8
      end if
   end do

   ! initialize pft-level variables
   do p = begp, endp
      l = plandunit(p)
      if (itypelun(l) == istsoil .or. itypelun(l) == istcrop) then

         ! carbon state variables
         if (ivt(p) == noveg) then
            leafc(p) = 0._r8
            leafc_storage(p) = 0._r8
         else
            if (evergreen(ivt(p)) == 1._r8) then
               leafc(p) = 1._r8
               leafc_storage(p) = 0._r8
            else if (ivt(p) >= npcropmin) then ! prognostic crop types
               leafc(p) = 0._r8
               leafc_storage(p) = 0._r8
            else
               leafc(p) = 0._r8
               leafc_storage(p) = 1._r8
            end if
         end if
         leafc_xfer(p) = 0._r8
         if ( crop_prog )then
            grainc(p) = 0._r8
            grainc_storage(p) = 0._r8
            grainc_xfer(p) = 0._r8
            fert(p) = 0._r8
            soyfixn(p) = 0._r8
         end if
         frootc(p) = 0._r8
         frootc_storage(p) = 0._r8
         frootc_xfer(p) = 0._r8
         livestemc(p) = 0._r8
         livestemc_storage(p) = 0._r8
         livestemc_xfer(p) = 0._r8

         ! tree types need to be initialized with some stem mass so that
         ! roughness length is not zero in canopy flux calculation

         if (woody(ivt(p)) == 1._r8) then
            deadstemc(p) = 0.1_r8
         else
            deadstemc(p) = 0._r8
         end if

         deadstemc_storage(p) = 0._r8
         deadstemc_xfer(p) = 0._r8
         livecrootc(p) = 0._r8
         livecrootc_storage(p) = 0._r8
         livecrootc_xfer(p) = 0._r8
         deadcrootc(p) = 0._r8
         deadcrootc_storage(p) = 0._r8
         deadcrootc_xfer(p) = 0._r8
         gresp_storage(p) = 0._r8
         gresp_xfer(p) = 0._r8
         cpool(p) = 0._r8
         xsmrpool(p) = 0._r8
         pft_ctrunc(p) = 0._r8
         dispvegc(p) = 0._r8
         storvegc(p) = 0._r8
         totpftc(p)  = 0._r8
         ! calculate totvegc explicitly so that it is available for the isotope 
         ! code on the first time step.
         totvegc(p)  = leafc(p) + leafc_storage(p) + leafc_xfer(p) + frootc(p) +  &
            frootc_storage(p) + frootc_xfer(p) + livestemc(p) + livestemc_storage(p) +  &
            livestemc_xfer(p) + deadstemc(p) + deadstemc_storage(p) + deadstemc_xfer(p) +  &
            livecrootc(p) + livecrootc_storage(p) + livecrootc_xfer(p) + deadcrootc(p) +  &
            deadcrootc_storage(p) + deadcrootc_xfer(p) + gresp_storage(p) +  &
            gresp_xfer(p) + cpool(p)

         if ( crop_prog )then
            totvegc(p) = totvegc(p) + grainc(p) + grainc_storage(p) + grainc_xfer(p)
         end if

         woodc(p)    = 0._r8


         if ( use_c13 ) then
            leafc13(p)               = leafc(p)               * c13ratio
            leafc13_storage(p)       = leafc_storage(p)       * c13ratio
            leafc13_xfer(p)          = leafc_xfer(p)          * c13ratio
            frootc13(p)              = frootc(p)              * c13ratio
            frootc13_storage(p)      = frootc_storage(p)      * c13ratio
            frootc13_xfer(p)         = frootc_xfer(p)         * c13ratio
            livestemc13(p)           = livestemc(p)           * c13ratio
            livestemc13_storage(p)   = livestemc_storage(p)   * c13ratio
            livestemc13_xfer(p)      = livestemc_xfer(p)      * c13ratio
            deadstemc13(p)           = deadstemc(p)           * c13ratio
            deadstemc13_storage(p)   = deadstemc_storage(p)   * c13ratio
            deadstemc13_xfer(p)      = deadstemc_xfer(p)      * c13ratio
            livecrootc13(p)          = livecrootc(p)          * c13ratio
            livecrootc13_storage(p)  = livecrootc_storage(p)  * c13ratio
            livecrootc13_xfer(p)     = livecrootc_xfer(p)     * c13ratio
            deadcrootc13(p)          = deadcrootc(p)          * c13ratio
            deadcrootc13_storage(p)  = deadcrootc_storage(p)  * c13ratio
            deadcrootc13_xfer(p)     = deadcrootc_xfer(p)     * c13ratio
            c13_gresp_storage(p)     = gresp_storage(p)       * c13ratio
            c13_gresp_xfer(p)        = gresp_xfer(p)          * c13ratio
            c13pool(p)               = cpool(p)               * c13ratio
            c13xsmrpool(p)           = xsmrpool(p)            * c13ratio
            c13_pft_ctrunc(p)        = pft_ctrunc(p)          * c13ratio
            
            ! calculate totvegc explicitly so that it is available for the isotope 
            ! code on the first time step.
            totvegc13(p)  = leafc13(p) + leafc13_storage(p) + leafc13_xfer(p) + frootc13(p) +  &
                 frootc13_storage(p) + frootc13_xfer(p) + livestemc13(p) + livestemc13_storage(p) +  &
                 livestemc13_xfer(p) + deadstemc13(p) + deadstemc13_storage(p) + deadstemc13_xfer(p) +  &
                 livecrootc13(p) + livecrootc13_storage(p) + livecrootc13_xfer(p) + deadcrootc13(p) +  &
                 deadcrootc13_storage(p) + deadcrootc13_xfer(p) + c13_gresp_storage(p) +  &
                 c13_gresp_xfer(p) + c13pool(p)
         endif
            
         if ( use_c14 ) then
            leafc14(p)               = leafc(p)               * c14ratio
            leafc14_storage(p)       = leafc_storage(p)       * c14ratio
            leafc14_xfer(p)          = leafc_xfer(p)          * c14ratio
            frootc14(p)              = frootc(p)              * c14ratio
            frootc14_storage(p)      = frootc_storage(p)      * c14ratio
            frootc14_xfer(p)         = frootc_xfer(p)         * c14ratio
            livestemc14(p)           = livestemc(p)           * c14ratio
            livestemc14_storage(p)   = livestemc_storage(p)   * c14ratio
            livestemc14_xfer(p)      = livestemc_xfer(p)      * c14ratio
            deadstemc14(p)           = deadstemc(p)           * c14ratio
            deadstemc14_storage(p)   = deadstemc_storage(p)   * c14ratio
            deadstemc14_xfer(p)      = deadstemc_xfer(p)      * c14ratio
            livecrootc14(p)          = livecrootc(p)          * c14ratio
            livecrootc14_storage(p)  = livecrootc_storage(p)  * c14ratio
            livecrootc14_xfer(p)     = livecrootc_xfer(p)     * c14ratio
            deadcrootc14(p)          = deadcrootc(p)          * c14ratio
            deadcrootc14_storage(p)  = deadcrootc_storage(p)  * c14ratio
            deadcrootc14_xfer(p)     = deadcrootc_xfer(p)     * c14ratio
            c14_gresp_storage(p)     = gresp_storage(p)       * c14ratio
            c14_gresp_xfer(p)        = gresp_xfer(p)          * c14ratio
            c14pool(p)               = cpool(p)               * c14ratio
            c14xsmrpool(p)           = xsmrpool(p)            * c14ratio
            c14_pft_ctrunc(p)        = pft_ctrunc(p)          * c14ratio
            
            ! calculate totvegc explicitly so that it is available for the isotope 
            ! code on the first time step.
            totvegc14(p)  = leafc14(p) + leafc14_storage(p) + leafc14_xfer(p) + frootc14(p) +  &
                 frootc14_storage(p) + frootc14_xfer(p) + livestemc14(p) + livestemc14_storage(p) +  &
                 livestemc14_xfer(p) + deadstemc14(p) + deadstemc14_storage(p) + deadstemc14_xfer(p) +  &
                 livecrootc14(p) + livecrootc14_storage(p) + livecrootc14_xfer(p) + deadcrootc14(p) +  &
                 deadcrootc14_storage(p) + deadcrootc14_xfer(p) + c14_gresp_storage(p) +  &
                 c14_gresp_xfer(p) + c14pool(p)
            
            rc14_atm(p) = c14ratio
         endif
                                
         ! nitrogen state variables
         if (ivt(p) == noveg) then
            leafn(p) = 0._r8
            leafn_storage(p) = 0._r8
         else
            leafn(p) = leafc(p) / leafcn(ivt(p))
            leafn_storage(p) = leafc_storage(p) / leafcn(ivt(p))
         end if

         leafn_xfer(p) = 0._r8
         if ( crop_prog )then
            grainn(p) = 0._r8
            grainn_storage(p) = 0._r8
            grainn_xfer(p) = 0._r8
         end if
         frootn(p) = 0._r8
         frootn_storage(p) = 0._r8
         frootn_xfer(p) = 0._r8
         livestemn(p) = 0._r8
         livestemn_storage(p) = 0._r8
         livestemn_xfer(p) = 0._r8

         ! tree types need to be initialized with some stem mass so that
         ! roughness length is not zero in canopy flux calculation

         if (woody(ivt(p)) == 1._r8) then
            deadstemn(p) = deadstemc(p) / deadwdcn(ivt(p))
         else
            deadstemn(p) = 0._r8
         end if

         deadstemn_storage(p) = 0._r8
         deadstemn_xfer(p) = 0._r8
         livecrootn(p) = 0._r8
         livecrootn_storage(p) = 0._r8
         livecrootn_xfer(p) = 0._r8
         deadcrootn(p) = 0._r8
         deadcrootn_storage(p) = 0._r8
         deadcrootn_xfer(p) = 0._r8
         retransn(p) = 0._r8
         npool(p) = 0._r8
         pft_ntrunc(p) = 0._r8
         dispvegn(p) = 0._r8
         storvegn(p) = 0._r8
         totvegn(p)  = 0._r8
         totpftn(p)  = 0._r8

         ! initialization for psnsun and psnsha required for
         ! proper arbitrary initialization of allocation routine
         ! in initial ecosysdyn call

         psnsun(p) = 0._r8
         psnsha(p) = 0._r8

         if ( use_c13 ) then
            c13_psnsun(p) = 0._r8
            c13_psnsha(p) = 0._r8
         endif
         
         if ( use_c14 ) then
            c14_psnsun(p) = 0._r8
            c14_psnsha(p) = 0._r8
         endif

         laisun(p) = 0._r8
         laisha(p) = 0._r8

         ! ecophysiological variables
         ! phenology variables
         dormant_flag(p) = 1._r8
         days_active(p) = 0._r8
         onset_flag(p) = 0._r8
         onset_counter(p) = 0._r8
         onset_gddflag(p) = 0._r8
         onset_fdd(p) = 0._r8
         onset_gdd(p) = 0._r8
         onset_swi(p) = 0.0_r8
         offset_flag(p) = 0._r8
         offset_counter(p) = 0._r8
         offset_fdd(p) = 0._r8
         offset_swi(p) = 0._r8
         lgsf(p) = 0._r8
         bglfr(p) = 0._r8
         bgtr(p) = 0._r8
         annavg_t2m(p) = 280._r8
         tempavg_t2m(p) = 0._r8
         fert_counter(p) = 0._r8
         grain_flag(p) = 0._r8

         ! non-phenology variables
         gpp(p) = 0._r8
         availc(p) = 0._r8
         xsmrpool_recover(p) = 0._r8
         alloc_pnow(p) = 1._r8
         c_allometry(p) = 0._r8
         n_allometry(p) = 0._r8
         plant_ndemand(p) = 0._r8
         tempsum_potential_gpp(p) = 0._r8
         annsum_potential_gpp(p) = 0._r8
         tempmax_retransn(p) = 0._r8
         annmax_retransn(p) = 0._r8
         avail_retransn(p) = 0._r8
         plant_nalloc(p) = 0._r8
         plant_calloc(p) = 0._r8
         excess_cflux(p) = 0._r8
         downreg(p) = 0._r8
         prev_leafc_to_litter(p) = 0._r8
         prev_frootc_to_litter(p) = 0._r8
         tempsum_npp(p) = 0._r8
         annsum_npp(p) = 0._r8
#if (defined CNDV)
         tempsum_litfall(p) = 0._r8
         annsum_litfall(p) = 0._r8
#endif

         if ( use_c13 ) then
            xsmrpool_c13ratio(p) = c13ratio
            rc13_canair(p) = 0._r8
            rc13_psnsun(p) = 0._r8
            rc13_psnsha(p) = 0._r8
            alphapsnsun(p) = 0._r8
            alphapsnsha(p) = 0._r8
         endif
         
      end if   ! end of if-istsoil block
   end do   ! end of loop over pfts  
#endif

end subroutine CNiniTimeVar
