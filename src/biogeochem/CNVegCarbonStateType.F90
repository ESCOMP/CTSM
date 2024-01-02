module CNVegCarbonStateType

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_const_mod  , only : SHR_CONST_PDB
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use pftconMod      , only : noveg, npcropmin, pftcon, nc3crop, nc3irrig
  use clm_varcon     , only : spval, c3_r2, c4_r2, c14ratio
  use clm_varctl     , only : iulog, use_cndv, use_crop
  use CNSharedParamsMod, only : use_matrixcn
  use decompMod      , only : bounds_type
  use abortutils     , only : endrun
  use spmdMod        , only : masterproc 
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use PatchType      , only : patch
  use CNSpeciesMod   , only : species_from_string, CN_SPECIES_C12
  use dynPatchStateUpdaterMod, only : patch_state_updater_type
  use CNVegComputeSeedMod, only : ComputeSeedAmounts
  use CropReprPoolsMod   , only : nrepr, get_repr_hist_fname, get_repr_rest_fname, get_repr_longname
  ! 
  ! !PUBLIC TYPES:
  implicit none
  private
  !

  type, public :: cnveg_carbonstate_type

     integer :: species  ! c12, c13, c14

     real(r8), pointer :: reproductivec_patch               (:,:) ! (gC/m2) reproductive (e.g., grain) C (crop model)
     real(r8), pointer :: reproductivec_storage_patch       (:,:) ! (gC/m2) reproductive (e.g., grain) C storage (crop model)
     real(r8), pointer :: reproductivec_xfer_patch          (:,:) ! (gC/m2) reproductive (e.g., grain) C transfer (crop model)
     real(r8), pointer :: matrix_cap_reproc_patch             (:) ! (gC/m2) Capacity of grain C
     real(r8), pointer :: matrix_cap_reproc_storage_patch     (:) ! (gC/m2) Capacity of grain storage C
     real(r8), pointer :: matrix_cap_reproc_xfer_patch        (:) ! (gC/m2) Capacity of grain transfer C
     real(r8), pointer :: leafc_patch                         (:) ! (gC/m2) leaf C
     real(r8), pointer :: leafc_storage_patch                 (:) ! (gC/m2) leaf C storage
     real(r8), pointer :: leafc_xfer_patch                    (:) ! (gC/m2) leaf C transfer
     real(r8), pointer :: matrix_cap_leafc_patch              (:) ! (gC/m2) Capacity of leaf C
     real(r8), pointer :: matrix_cap_leafc_storage_patch      (:) ! (gC/m2) Capacity of leaf C storage
     real(r8), pointer :: matrix_cap_leafc_xfer_patch         (:) ! (gC/m2) Capacity of leaf C transfer
     real(r8), pointer :: leafc_storage_xfer_acc_patch        (:) ! (gC/m2) Accmulated leaf C transfer
     real(r8), pointer :: storage_cdemand_patch               (:) ! (gC/m2)       C use from the C storage pool 
     real(r8), pointer :: frootc_patch                        (:) ! (gC/m2) fine root C
     real(r8), pointer :: frootc_storage_patch                (:) ! (gC/m2) fine root C storage
     real(r8), pointer :: frootc_xfer_patch                   (:) ! (gC/m2) fine root C transfer
     real(r8), pointer :: matrix_cap_frootc_patch             (:) ! (gC/m2) Capacity of fine root C
     real(r8), pointer :: matrix_cap_frootc_storage_patch     (:) ! (gC/m2) Capacity of fine root C storage
     real(r8), pointer :: matrix_cap_frootc_xfer_patch        (:) ! (gC/m2) Capacity of fine root C transfer
     real(r8), pointer :: livestemc_patch                     (:) ! (gC/m2) live stem C
     real(r8), pointer :: livestemc_storage_patch             (:) ! (gC/m2) live stem C storage
     real(r8), pointer :: livestemc_xfer_patch                (:) ! (gC/m2) live stem C transfer
     real(r8), pointer :: matrix_cap_livestemc_patch          (:) ! (gC/m2) Capacity of live stem C
     real(r8), pointer :: matrix_cap_livestemc_storage_patch  (:) ! (gC/m2) Capacity of live stem C storage
     real(r8), pointer :: matrix_cap_livestemc_xfer_patch     (:) ! (gC/m2) Capacity of live stem C transfer
     real(r8), pointer :: deadstemc_patch                     (:) ! (gC/m2) dead stem C
     real(r8), pointer :: deadstemc_storage_patch             (:) ! (gC/m2) dead stem C storage
     real(r8), pointer :: deadstemc_xfer_patch                (:) ! (gC/m2) dead stem C transfer
     real(r8), pointer :: matrix_cap_deadstemc_patch          (:) ! (gC/m2) Capacity of dead stem C
     real(r8), pointer :: matrix_cap_deadstemc_storage_patch  (:) ! (gC/m2) Capacity of dead stem C storage
     real(r8), pointer :: matrix_cap_deadstemc_xfer_patch     (:) ! (gC/m2) Capacity of dead stem C transfer
     real(r8), pointer :: livecrootc_patch                    (:) ! (gC/m2) live coarse root C
     real(r8), pointer :: livecrootc_storage_patch            (:) ! (gC/m2) live coarse root C storage
     real(r8), pointer :: livecrootc_xfer_patch               (:) ! (gC/m2) live coarse root C transfer
     real(r8), pointer :: matrix_cap_livecrootc_patch         (:) ! (gC/m2) Capacity of live coarse root C
     real(r8), pointer :: matrix_cap_livecrootc_storage_patch (:) ! (gC/m2) Capacity of live coarse root C storage
     real(r8), pointer :: matrix_cap_livecrootc_xfer_patch    (:) ! (gC/m2) Capacity of live coarse root C transfer
     real(r8), pointer :: deadcrootc_patch                    (:) ! (gC/m2) dead coarse root C
     real(r8), pointer :: deadcrootc_storage_patch            (:) ! (gC/m2) dead coarse root C storage
     real(r8), pointer :: deadcrootc_xfer_patch               (:) ! (gC/m2) dead coarse root C transfer
     real(r8), pointer :: matrix_cap_deadcrootc_patch         (:) ! (gC/m2) Capacity of dead coarse root C
     real(r8), pointer :: matrix_cap_deadcrootc_storage_patch (:) ! (gC/m2) Capacity of dead coarse root C storage
     real(r8), pointer :: matrix_cap_deadcrootc_xfer_patch    (:) ! (gC/m2) Capacity of dead coarse root C transfer
     real(r8), pointer :: gresp_storage_patch                 (:) ! (gC/m2) growth respiration storage
     real(r8), pointer :: gresp_xfer_patch                    (:) ! (gC/m2) growth respiration transfer
     real(r8), pointer :: cpool_patch                         (:) ! (gC/m2) temporary photosynthate C pool
     real(r8), pointer :: xsmrpool_patch                      (:) ! (gC/m2) abstract C pool to meet excess MR demand
     real(r8), pointer :: xsmrpool_loss_patch                 (:) ! (gC/m2) abstract C pool to meet excess MR demand loss
     real(r8), pointer :: ctrunc_patch                        (:) ! (gC/m2) patch-level sink for C truncation
     real(r8), pointer :: woodc_patch                         (:) ! (gC/m2) wood C
     real(r8), pointer :: leafcmax_patch                      (:) ! (gC/m2) ann max leaf C
     real(r8), pointer :: rootc_col                           (:) ! (gC/m2) root carbon at column level (fire)
     real(r8), pointer :: leafc_col                           (:) ! (gC/m2) column-level leafc (fire)
     real(r8), pointer :: deadstemc_col                       (:) ! (gC/m2) column-level deadstemc (fire)
     real(r8), pointer :: fuelc_col                           (:) ! fuel load outside cropland
     real(r8), pointer :: fuelc_crop_col                      (:) ! fuel load for cropland
     real(r8), pointer :: cropseedc_deficit_patch             (:) ! (gC/m2) pool for seeding new crop growth; this is a NEGATIVE term, indicating the amount of seed usage that needs to be repaid
! initial pool size of year for matrix
     real(r8), pointer :: leafc0_patch                        (:) ! (gC/m2) Initial value of leaf C for SASU
     real(r8), pointer :: leafc0_storage_patch                (:) ! (gC/m2) Initial value of leaf C storage for SASU
     real(r8), pointer :: leafc0_xfer_patch                   (:) ! (gC/m2) Initial value of leaf C transfer for SASU
     real(r8), pointer :: frootc0_patch                       (:) ! (gC/m2) Initial value of fine root C for SASU
     real(r8), pointer :: frootc0_storage_patch               (:) ! (gC/m2) Initial value of fine root C storage for SASU
     real(r8), pointer :: frootc0_xfer_patch                  (:) ! (gC/m2) Initial value of fine root C transfer for SASU
     real(r8), pointer :: livestemc0_patch                    (:) ! (gC/m2) Initial value of live stem C for SASU
     real(r8), pointer :: livestemc0_storage_patch            (:) ! (gC/m2) Initial value of live stem C storage for SASU
     real(r8), pointer :: livestemc0_xfer_patch               (:) ! (gC/m2) Initial value of live stem C transfer for SASU
     real(r8), pointer :: deadstemc0_patch                    (:) ! (gC/m2) Initial value of dead stem C for SASU
     real(r8), pointer :: deadstemc0_storage_patch            (:) ! (gC/m2) Initial value of dead stem C storage for SASU
     real(r8), pointer :: deadstemc0_xfer_patch               (:) ! (gC/m2) Initial value of dead stem C transfer for SASU
     real(r8), pointer :: livecrootc0_patch                   (:) ! (gC/m2) Initial value of live coarse root C for SASU
     real(r8), pointer :: livecrootc0_storage_patch           (:) ! (gC/m2) Initial value of live coarse root C storage for SASU
     real(r8), pointer :: livecrootc0_xfer_patch              (:) ! (gC/m2) Initial value of live coarse root C transfer for SASU
     real(r8), pointer :: deadcrootc0_patch                   (:) ! (gC/m2) Initial value of dead coarse root C for SASU
     real(r8), pointer :: deadcrootc0_storage_patch           (:) ! (gC/m2) Initial value of dead coarse root C storage for SASU
     real(r8), pointer :: deadcrootc0_xfer_patch              (:) ! (gC/m2) Initial value of dead coarse root C transfer for SASU
     real(r8), pointer :: reproc0_patch                       (:) ! (gC/m2) Initial value of fine grain C for SASU
     real(r8), pointer :: reproc0_storage_patch               (:) ! (gC/m2) Initial value of fine grain C storage for SASU
     real(r8), pointer :: reproc0_xfer_patch                  (:) ! (gC/m2) Initial value of fine grain C transfer for SASU

     ! pools for dynamic landcover
     real(r8), pointer :: seedc_grc                           (:) ! (gC/m2) gridcell-level pool for seeding new PFTs via dynamic landcover

     ! summary (diagnostic) state variables, not involved in mass balance
     real(r8), pointer :: dispvegc_patch                      (:) ! (gC/m2) displayed veg carbon, excluding storage and cpool
     real(r8), pointer :: storvegc_patch                      (:) ! (gC/m2) stored vegetation carbon, excluding cpool
     
     logical, private  :: dribble_crophrv_xsmrpool_2atm           ! Flag to indicate if should harvest xsmrpool to the atmosphere
                                                                  ! it originates and is defined in CNVegetationFacade.F90
     
     ! Total C pools
     real(r8), pointer :: totc_patch                          (:) ! (gC/m2) total patch-level carbon, including cpool
     real(r8), pointer :: totvegc_patch                       (:) ! (gC/m2) total vegetation carbon, excluding cpool
     real(r8), pointer :: totvegc_col                         (:) ! (gC/m2) total vegetation carbon, excluding cpool averaged to column (p2c)      
     real(r8), pointer :: totc_p2c_col                        (:) ! (gC/m2) totc_patch averaged to col

! Accumulation variables are accumulated for a whole year. They are used for matrix spinup and calculation of diagnostic variables
     real(r8), pointer :: matrix_calloc_leaf_acc_patch        (:) ! (gC/m2/year) Input C allocated to leaf during this year 
     real(r8), pointer :: matrix_calloc_leafst_acc_patch      (:) ! (gC/m2/year) Input C allocated to leaf storage during this year
     real(r8), pointer :: matrix_calloc_froot_acc_patch       (:) ! (gC/m2/year) Input C allocated to fine root during this year
     real(r8), pointer :: matrix_calloc_frootst_acc_patch     (:) ! (gC/m2/year) Input C allocated to fine root storage during this year
     real(r8), pointer :: matrix_calloc_livestem_acc_patch    (:) ! (gC/m2/year) Input C allocated to live stem during this year
     real(r8), pointer :: matrix_calloc_livestemst_acc_patch  (:) ! (gC/m2/year) Input C allocated to live stem storage during this year
     real(r8), pointer :: matrix_calloc_deadstem_acc_patch    (:) ! (gC/m2/year) Input C allocated to dead stem during this year
     real(r8), pointer :: matrix_calloc_deadstemst_acc_patch  (:) ! (gC/m2/year) Input C allocated to dead stem storage during this year
     real(r8), pointer :: matrix_calloc_livecroot_acc_patch   (:) ! (gC/m2/year) Input C allocated to live coarse root during this year
     real(r8), pointer :: matrix_calloc_livecrootst_acc_patch (:) ! (gC/m2/year) Input C allocated to live coarse root storage during this year
     real(r8), pointer :: matrix_calloc_deadcroot_acc_patch   (:) ! (gC/m2/year) Input C allocated to dead coarse root during this year
     real(r8), pointer :: matrix_calloc_deadcrootst_acc_patch (:) ! (gC/m2/year) Input C allocated to dead coarse root storage during this year
     real(r8), pointer :: matrix_calloc_grain_acc_patch       (:) ! (gC/m2/year) Input C allocated to grain during this year
     real(r8), pointer :: matrix_calloc_grainst_acc_patch     (:) ! (gC/m2/year) Input C allocated to grain storage during this year

     real(r8), pointer :: matrix_ctransfer_leafst_to_leafxf_acc_patch           (:) ! (gC/m2/year) C transfer from leaf storage to leaf transfer pool during this year
     real(r8), pointer :: matrix_ctransfer_leafxf_to_leaf_acc_patch             (:) ! (gC/m2/year) C transfer from leaf transfer to leaf pool during this year
     real(r8), pointer :: matrix_ctransfer_frootst_to_frootxf_acc_patch         (:) ! (gC/m2/year) C transfer from fine root storage to fine root transfer pool during this year
     real(r8), pointer :: matrix_ctransfer_frootxf_to_froot_acc_patch           (:) ! (gC/m2/year) C transfer from fine root transfer to fine root pool during this year
     real(r8), pointer :: matrix_ctransfer_livestemst_to_livestemxf_acc_patch   (:) ! (gC/m2/year) C transfer from live stem storage to live stem transfer pool during this year
     real(r8), pointer :: matrix_ctransfer_livestemxf_to_livestem_acc_patch     (:) ! (gC/m2/year) C transfer from live stem transfer to live stem pool during this year
     real(r8), pointer :: matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch   (:) ! (gC/m2/year) C transfer from dead stem storage to dead stem transfer pool during this year
     real(r8), pointer :: matrix_ctransfer_deadstemxf_to_deadstem_acc_patch     (:) ! (gC/m2/year) C transfer from dead stem transfer to dead stem pool during this year
     real(r8), pointer :: matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch (:) ! (gC/m2/year) C transfer from live coarse root storage to live coarse root transfer pool during this year
     real(r8), pointer :: matrix_ctransfer_livecrootxf_to_livecroot_acc_patch   (:) ! (gC/m2/year) C transfer from live coarse root transfer to live coarse root pool during this year
     real(r8), pointer :: matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch (:) ! (gC/m2/year) C transfer from dead coarse root storage to dead coarse root transfer pool during this year
     real(r8), pointer :: matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch   (:) ! (gC/m2/year) C transfer from dead coarse root transfer to dead coarse root pool during this year
     real(r8), pointer :: matrix_ctransfer_grainst_to_grainxf_acc_patch         (:) ! (gC/m2/year) C transfer from grain storage to grain transfer pool during this year
     real(r8), pointer :: matrix_ctransfer_grainxf_to_grain_acc_patch           (:) ! (gC/m2/year) C transfer from grain transfer to grain pool during this year
     real(r8), pointer :: matrix_ctransfer_livestem_to_deadstem_acc_patch       (:) ! (gC/m2/year) C transfer from live stem to dead stem pool during this year
     real(r8), pointer :: matrix_ctransfer_livecroot_to_deadcroot_acc_patch     (:) ! (gC/m2/year) C transfer from live coarse root to dead coarse root pool during this year

     real(r8), pointer :: matrix_cturnover_leaf_acc_patch             (:) ! (gC/m2/year) C turnover from leaf
     real(r8), pointer :: matrix_cturnover_leafst_acc_patch           (:) ! (gC/m2/year) C turnover from leaf storage
     real(r8), pointer :: matrix_cturnover_leafxf_acc_patch           (:) ! (gC/m2/year) C turnover from leaf transfer
     real(r8), pointer :: matrix_cturnover_froot_acc_patch            (:) ! (gC/m2/year) C turnover from fine root
     real(r8), pointer :: matrix_cturnover_frootst_acc_patch          (:) ! (gC/m2/year) C turnover from fine root storage
     real(r8), pointer :: matrix_cturnover_frootxf_acc_patch          (:) ! (gC/m2/year) C turnover from fine root transfer
     real(r8), pointer :: matrix_cturnover_livestem_acc_patch         (:) ! (gC/m2/year) C turnover from live stem
     real(r8), pointer :: matrix_cturnover_livestemst_acc_patch       (:) ! (gC/m2/year) C turnover from live stem storage
     real(r8), pointer :: matrix_cturnover_livestemxf_acc_patch       (:) ! (gC/m2/year) C turnover from live stem transfer
     real(r8), pointer :: matrix_cturnover_deadstem_acc_patch         (:) ! (gC/m2/year) C turnover from dead stem
     real(r8), pointer :: matrix_cturnover_deadstemst_acc_patch       (:) ! (gC/m2/year) C turnover from dead stem storage
     real(r8), pointer :: matrix_cturnover_deadstemxf_acc_patch       (:) ! (gC/m2/year) C turnover from dead stem transfer
     real(r8), pointer :: matrix_cturnover_livecroot_acc_patch        (:) ! (gC/m2/year) C turnover from live coarse root
     real(r8), pointer :: matrix_cturnover_livecrootst_acc_patch      (:) ! (gC/m2/year) C turnover from live coarse root storage
     real(r8), pointer :: matrix_cturnover_livecrootxf_acc_patch      (:) ! (gC/m2/year) C turnover from live coarse root transfer
     real(r8), pointer :: matrix_cturnover_deadcroot_acc_patch        (:) ! (gC/m2/year) C turnover from dead coarse root
     real(r8), pointer :: matrix_cturnover_deadcrootst_acc_patch      (:) ! (gC/m2/year) C turnover from dead coarse root storage
     real(r8), pointer :: matrix_cturnover_deadcrootxf_acc_patch      (:) ! (gC/m2/year) C turnover from dead coarse root transfer
     real(r8), pointer :: matrix_cturnover_grain_acc_patch            (:) ! (gC/m2/year) C turnover from grain 
     real(r8), pointer :: matrix_cturnover_grainst_acc_patch          (:) ! (gC/m2/year) C turnover from grain storage
     real(r8), pointer :: matrix_cturnover_grainxf_acc_patch          (:) ! (gC/m2/year) C turnover from grain transfer

     real(r8), pointer :: grainc_SASUsave_patch               (:) ! (gC/m2) grain C (crop model)
     real(r8), pointer :: grainc_storage_SASUsave_patch       (:) ! (gC/m2) grain C storage (crop model)
     real(r8), pointer :: leafc_SASUsave_patch                (:) ! (gC/m2) leaf C
     real(r8), pointer :: leafc_storage_SASUsave_patch        (:) ! (gC/m2) leaf C storage
     real(r8), pointer :: leafc_xfer_SASUsave_patch           (:) ! (gC/m2) leaf C transfer
     real(r8), pointer :: frootc_SASUsave_patch               (:) ! (gC/m2) fine root C
     real(r8), pointer :: frootc_storage_SASUsave_patch       (:) ! (gC/m2) fine root C storage
     real(r8), pointer :: frootc_xfer_SASUsave_patch          (:) ! (gC/m2) fine root C transfer
     real(r8), pointer :: livestemc_SASUsave_patch            (:) ! (gC/m2) live stem C
     real(r8), pointer :: livestemc_storage_SASUsave_patch    (:) ! (gC/m2) live stem C storage
     real(r8), pointer :: livestemc_xfer_SASUsave_patch       (:) ! (gC/m2) live stem C transfer
     real(r8), pointer :: deadstemc_SASUsave_patch            (:) ! (gC/m2) dead stem C
     real(r8), pointer :: deadstemc_storage_SASUsave_patch    (:) ! (gC/m2) dead stem C storage
     real(r8), pointer :: deadstemc_xfer_SASUsave_patch       (:) ! (gC/m2) dead stem C transfer
     real(r8), pointer :: livecrootc_SASUsave_patch           (:) ! (gC/m2) live coarse root C
     real(r8), pointer :: livecrootc_storage_SASUsave_patch   (:) ! (gC/m2) live coarse root C storage
     real(r8), pointer :: livecrootc_xfer_SASUsave_patch      (:) ! (gC/m2) live coarse root C transfer
     real(r8), pointer :: deadcrootc_SASUsave_patch           (:) ! (gC/m2) dead coarse root C
     real(r8), pointer :: deadcrootc_storage_SASUsave_patch   (:) ! (gC/m2) dead coarse root C storage
     real(r8), pointer :: deadcrootc_xfer_SASUsave_patch      (:) ! (gC/m2) dead coarse root C transfer

   contains

     procedure , public  :: Init   
     procedure , public  :: SetValues
     procedure , public  :: ZeroDwt
     procedure , public  :: Restart
     procedure , public  :: Summary => Summary_carbonstate
     procedure , public  :: DynamicPatchAdjustments   ! adjust state variables when patch areas change
     
     procedure , private :: InitAllocate    ! Allocate arrays
     procedure , private :: InitReadNML     ! Read in namelist
     procedure , private :: InitHistory     ! Initialize history
     procedure , private :: InitCold        ! Initialize arrays for a cold-start

  end type cnveg_carbonstate_type

  real(r8), public  :: spinup_factor_deadwood = 1.0_r8        ! Spinup factor used for this simulation
  real(r8), public  :: spinup_factor_AD       = 10.0_r8       ! Spinup factor used when in Accelerated Decomposition mode

  ! !PRIVATE DATA:

  type, private :: cnvegcarbonstate_const_type
      ! !PRIVATE MEMBER DATA:
      real(r8) :: initial_vegC = 20._r8    ! Initial vegetation carbon for leafc/frootc and storage
  end type
  type(cnvegcarbonstate_const_type), private :: cnvegcstate_const    ! Constants used here
  character(len=*), parameter :: sourcefile = &
       __FILE__

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, carbon_type, ratio, NLFilename, &
                  dribble_crophrv_xsmrpool_2atm, alloc_full_veg, c12_cnveg_carbonstate_inst)

    class(cnveg_carbonstate_type)                       :: this
    type(bounds_type)            , intent(in)           :: bounds  
    real(r8)                     , intent(in)           :: ratio
    character(len=*)             , intent(in)           :: carbon_type                ! Carbon isotope type C12, C13 or C1
    character(len=*)             , intent(in)           :: NLFilename                 ! Namelist filename
    logical                      , intent(in)           :: dribble_crophrv_xsmrpool_2atm
    logical                      , intent(in)           :: alloc_full_veg               ! total number of bgc patches (non-fates)
    type(cnveg_carbonstate_type) , intent(in), optional :: c12_cnveg_carbonstate_inst ! cnveg_carbonstate for C12 (if C13 or C14)
    !-----------------------------------------------------------------------

    this%species = species_from_string(carbon_type)

    this%dribble_crophrv_xsmrpool_2atm = dribble_crophrv_xsmrpool_2atm

    call this%InitAllocate ( bounds, alloc_full_veg)
    if(alloc_full_veg)then
       call this%InitReadNML  ( NLFilename )
       call this%InitHistory ( bounds, carbon_type)
       if (present(c12_cnveg_carbonstate_inst)) then
          call this%InitCold  ( bounds, ratio, carbon_type, c12_cnveg_carbonstate_inst )
       else
          call this%InitCold  ( bounds, ratio, carbon_type )
       end if
    end if
    
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitReadNML(this, NLFilename)
    !
    ! !DESCRIPTION:
    ! Read the namelist for CNVegCarbonState
    !
    !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    !
    ! !ARGUMENTS:
    class(cnveg_carbonstate_type)                       :: this
    character(len=*)             , intent(in)           :: NLFilename                 ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'InitReadNML'
    character(len=*), parameter :: nmlname = 'cnvegcarbonstate'   ! MUST match what is in namelist below
    !-----------------------------------------------------------------------
    real(r8) :: initial_vegC
    namelist /cnvegcarbonstate/ initial_vegC

    initial_vegC = cnvegcstate_const%initial_vegC

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=cnvegcarbonstate, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (initial_vegC            , mpicom)

    cnvegcstate_const%initial_vegC = initial_vegC

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=cnvegcarbonstate)    ! Name here MUST be the same as in nmlname above!
       write(iulog,*) ' '
    end if

    !-----------------------------------------------------------------------

  end subroutine InitReadNML

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds, alloc_full_veg)
    !
    ! !ARGUMENTS:
    class (cnveg_carbonstate_type) :: this
    type(bounds_type), intent(in) :: bounds
    logical,intent(in)            :: alloc_full_veg ! Total number of bgc patches on the proc (non_fates)
    !
    ! !LOCAL VARIABLES:
    integer           :: begp,endp
    integer           :: begc,endc
    integer           :: begg,endg
    !------------------------------------------------------------------------
    
    if(alloc_full_veg)then
       begp = bounds%begp; endp = bounds%endp
       begc = bounds%begc; endc = bounds%endc
       begg = bounds%begg; endg = bounds%endg
    else
       begp = 0;endp=0
       begc = 0;endc=0
       begg = 0;endg=0
    end if

    allocate(this%leafc_patch                            (begp:endp)) ; this%leafc_patch                        (:) = nan
    allocate(this%leafc_storage_patch                    (begp:endp)) ; this%leafc_storage_patch                (:) = nan
    allocate(this%leafc_xfer_patch                       (begp:endp)) ; this%leafc_xfer_patch                   (:) = nan
    if(use_matrixcn)then
       allocate(this%matrix_cap_leafc_patch              (begp:endp)) ; this%matrix_cap_leafc_patch             (:) = nan
       allocate(this%matrix_cap_leafc_storage_patch      (begp:endp)) ; this%matrix_cap_leafc_storage_patch     (:) = nan
       allocate(this%matrix_cap_leafc_xfer_patch         (begp:endp)) ; this%matrix_cap_leafc_xfer_patch        (:) = nan
    end if
    allocate(this%leafc_storage_xfer_acc_patch           (begp:endp)) ; this%leafc_storage_xfer_acc_patch       (:) = nan
    allocate(this%storage_cdemand_patch                  (begp:endp)) ; this%storage_cdemand_patch              (:) = nan
    allocate(this%frootc_patch                           (begp:endp)) ; this%frootc_patch                       (:) = nan
    allocate(this%frootc_storage_patch                   (begp:endp)) ; this%frootc_storage_patch               (:) = nan
    allocate(this%frootc_xfer_patch                      (begp:endp)) ; this%frootc_xfer_patch                  (:) = nan
    if(use_matrixcn)then
       allocate(this%matrix_cap_frootc_patch             (begp:endp)) ; this%matrix_cap_frootc_patch            (:) = nan
       allocate(this%matrix_cap_frootc_storage_patch     (begp:endp)) ; this%matrix_cap_frootc_storage_patch    (:) = nan
       allocate(this%matrix_cap_frootc_xfer_patch        (begp:endp)) ; this%matrix_cap_frootc_xfer_patch       (:) = nan
    end if
    allocate(this%livestemc_patch                        (begp:endp)) ; this%livestemc_patch                    (:) = nan
    allocate(this%livestemc_storage_patch                (begp:endp)) ; this%livestemc_storage_patch            (:) = nan
    allocate(this%livestemc_xfer_patch                   (begp:endp)) ; this%livestemc_xfer_patch               (:) = nan
    if(use_matrixcn)then
       allocate(this%matrix_cap_livestemc_patch          (begp:endp)) ; this%matrix_cap_livestemc_patch         (:) = nan
       allocate(this%matrix_cap_livestemc_storage_patch  (begp:endp)) ; this%matrix_cap_livestemc_storage_patch (:) = nan
       allocate(this%matrix_cap_livestemc_xfer_patch     (begp:endp)) ; this%matrix_cap_livestemc_xfer_patch    (:) = nan
    end if
    allocate(this%deadstemc_patch                        (begp:endp)) ; this%deadstemc_patch                    (:) = nan
    allocate(this%deadstemc_storage_patch                (begp:endp)) ; this%deadstemc_storage_patch            (:) = nan
    allocate(this%deadstemc_xfer_patch                   (begp:endp)) ; this%deadstemc_xfer_patch               (:) = nan
    if(use_matrixcn)then
       allocate(this%matrix_cap_deadstemc_patch          (begp:endp)) ; this%matrix_cap_deadstemc_patch         (:) = nan
       allocate(this%matrix_cap_deadstemc_storage_patch  (begp:endp)) ; this%matrix_cap_deadstemc_storage_patch (:) = nan
       allocate(this%matrix_cap_deadstemc_xfer_patch     (begp:endp)) ; this%matrix_cap_deadstemc_xfer_patch    (:) = nan
    end if
    allocate(this%livecrootc_patch                       (begp:endp)) ; this%livecrootc_patch                   (:) = nan
    allocate(this%livecrootc_storage_patch               (begp:endp)) ; this%livecrootc_storage_patch           (:) = nan
    allocate(this%livecrootc_xfer_patch                  (begp:endp)) ; this%livecrootc_xfer_patch              (:) = nan
    if(use_matrixcn)then
       allocate(this%matrix_cap_livecrootc_patch         (begp:endp)) ; this%matrix_cap_livecrootc_patch        (:) = nan
       allocate(this%matrix_cap_livecrootc_storage_patch (begp:endp)) ; this%matrix_cap_livecrootc_storage_patch(:) = nan
       allocate(this%matrix_cap_livecrootc_xfer_patch    (begp:endp)) ; this%matrix_cap_livecrootc_xfer_patch   (:) = nan
    end if
    allocate(this%deadcrootc_patch                       (begp:endp)) ; this%deadcrootc_patch                   (:) = nan
    allocate(this%deadcrootc_storage_patch               (begp:endp)) ; this%deadcrootc_storage_patch           (:) = nan
    allocate(this%deadcrootc_xfer_patch                  (begp:endp)) ; this%deadcrootc_xfer_patch              (:) = nan
    if(use_matrixcn)then
       allocate(this%matrix_cap_deadcrootc_patch         (begp:endp)) ; this%matrix_cap_deadcrootc_patch        (:) = nan
       allocate(this%matrix_cap_deadcrootc_storage_patch (begp:endp)) ; this%matrix_cap_deadcrootc_storage_patch(:) = nan
       allocate(this%matrix_cap_deadcrootc_xfer_patch    (begp:endp)) ; this%matrix_cap_deadcrootc_xfer_patch   (:) = nan
    end if
    allocate(this%gresp_storage_patch                    (begp:endp)) ; this%gresp_storage_patch                (:) = nan
    allocate(this%gresp_xfer_patch                       (begp:endp)) ; this%gresp_xfer_patch                   (:) = nan
    allocate(this%cpool_patch                            (begp:endp)) ; this%cpool_patch                        (:) = nan
    allocate(this%xsmrpool_patch                         (begp:endp)) ; this%xsmrpool_patch                     (:) = nan
    allocate(this%xsmrpool_loss_patch                    (begp:endp)) ; this%xsmrpool_loss_patch                (:) = nan
    allocate(this%ctrunc_patch                           (begp:endp)) ; this%ctrunc_patch                       (:) = nan
    allocate(this%dispvegc_patch                         (begp:endp)) ; this%dispvegc_patch                     (:) = nan
    allocate(this%storvegc_patch                         (begp:endp)) ; this%storvegc_patch                     (:) = nan
    allocate(this%leafcmax_patch                         (begp:endp)) ; this%leafcmax_patch                     (:) = nan
    allocate(this%totc_patch                             (begp:endp))  ; this%totc_patch                        (:) = nan
    allocate(this%reproductivec_patch             (begp:endp, nrepr)) ; this%reproductivec_patch               (:,:) = nan
    allocate(this%reproductivec_storage_patch     (begp:endp, nrepr)) ; this%reproductivec_storage_patch       (:,:) = nan
    allocate(this%reproductivec_xfer_patch        (begp:endp, nrepr)) ; this%reproductivec_xfer_patch          (:,:) = nan
    if(use_matrixcn)then
       allocate(this%matrix_cap_reproc_patch             (begp:endp)) ; this%matrix_cap_reproc_patch            (:) = nan
       allocate(this%matrix_cap_reproc_storage_patch     (begp:endp)) ; this%matrix_cap_reproc_storage_patch    (:) = nan
       allocate(this%matrix_cap_reproc_xfer_patch        (begp:endp)) ; this%matrix_cap_reproc_xfer_patch       (:) = nan
    end if
    allocate(this%woodc_patch                            (begp:endp)) ; this%woodc_patch                        (:) = nan     
!initial pool size of year for matrix
    if(use_matrixcn)then
       allocate(this%leafc0_patch                        (begp:endp)) ; this%leafc0_patch                       (:) = nan
       allocate(this%leafc0_storage_patch                (begp:endp)) ; this%leafc0_storage_patch               (:) = nan
       allocate(this%leafc0_xfer_patch                   (begp:endp)) ; this%leafc0_xfer_patch                  (:) = nan
       allocate(this%frootc0_patch                       (begp:endp)) ; this%frootc0_patch                      (:) = nan
       allocate(this%frootc0_storage_patch               (begp:endp)) ; this%frootc0_storage_patch              (:) = nan
       allocate(this%frootc0_xfer_patch                  (begp:endp)) ; this%frootc0_xfer_patch                 (:) = nan
       allocate(this%livestemc0_patch                    (begp:endp)) ; this%livestemc0_patch                   (:) = nan
       allocate(this%livestemc0_storage_patch            (begp:endp)) ; this%livestemc0_storage_patch           (:) = nan
       allocate(this%livestemc0_xfer_patch               (begp:endp)) ; this%livestemc0_xfer_patch              (:) = nan
       allocate(this%deadstemc0_patch                    (begp:endp)) ; this%deadstemc0_patch                   (:) = nan
       allocate(this%deadstemc0_storage_patch            (begp:endp)) ; this%deadstemc0_storage_patch           (:) = nan
       allocate(this%deadstemc0_xfer_patch               (begp:endp)) ; this%deadstemc0_xfer_patch              (:) = nan
       allocate(this%livecrootc0_patch                   (begp:endp)) ; this%livecrootc0_patch                  (:) = nan
       allocate(this%livecrootc0_storage_patch           (begp:endp)) ; this%livecrootc0_storage_patch          (:) = nan
       allocate(this%livecrootc0_xfer_patch              (begp:endp)) ; this%livecrootc0_xfer_patch             (:) = nan
       allocate(this%deadcrootc0_patch                   (begp:endp)) ; this%deadcrootc0_patch                  (:) = nan
       allocate(this%deadcrootc0_storage_patch           (begp:endp)) ; this%deadcrootc0_storage_patch          (:) = nan
       allocate(this%deadcrootc0_xfer_patch              (begp:endp)) ; this%deadcrootc0_xfer_patch             (:) = nan
       allocate(this%reproc0_patch                       (begp:endp)) ; this%reproc0_patch                      (:) = nan
       allocate(this%reproc0_storage_patch               (begp:endp)) ; this%reproc0_storage_patch              (:) = nan
       allocate(this%reproc0_xfer_patch                  (begp:endp)) ; this%reproc0_xfer_patch                 (:) = nan
 
       allocate(this%leafc_SASUsave_patch                (begp:endp)) ; this%leafc_SASUsave_patch               (:) = nan
       allocate(this%leafc_storage_SASUsave_patch        (begp:endp)) ; this%leafc_storage_SASUsave_patch       (:) = nan
       allocate(this%leafc_xfer_SASUsave_patch           (begp:endp)) ; this%leafc_xfer_SASUsave_patch          (:) = nan
       allocate(this%frootc_SASUsave_patch               (begp:endp)) ; this%frootc_SASUsave_patch              (:) = nan
       allocate(this%frootc_storage_SASUsave_patch       (begp:endp)) ; this%frootc_storage_SASUsave_patch      (:) = nan
       allocate(this%frootc_xfer_SASUsave_patch          (begp:endp)) ; this%frootc_xfer_SASUsave_patch         (:) = nan
       allocate(this%livestemc_SASUsave_patch            (begp:endp)) ; this%livestemc_SASUsave_patch           (:) = nan
       allocate(this%livestemc_storage_SASUsave_patch    (begp:endp)) ; this%livestemc_storage_SASUsave_patch   (:) = nan
       allocate(this%livestemc_xfer_SASUsave_patch       (begp:endp)) ; this%livestemc_xfer_SASUsave_patch      (:) = nan
       allocate(this%deadstemc_SASUsave_patch            (begp:endp)) ; this%deadstemc_SASUsave_patch           (:) = nan
       allocate(this%deadstemc_storage_SASUsave_patch    (begp:endp)) ; this%deadstemc_storage_SASUsave_patch   (:) = nan
       allocate(this%deadstemc_xfer_SASUsave_patch       (begp:endp)) ; this%deadstemc_xfer_SASUsave_patch      (:) = nan
       allocate(this%livecrootc_SASUsave_patch           (begp:endp)) ; this%livecrootc_SASUsave_patch          (:) = nan
       allocate(this%livecrootc_storage_SASUsave_patch   (begp:endp)) ; this%livecrootc_storage_SASUsave_patch  (:) = nan
       allocate(this%livecrootc_xfer_SASUsave_patch      (begp:endp)) ; this%livecrootc_xfer_SASUsave_patch     (:) = nan
       allocate(this%deadcrootc_SASUsave_patch           (begp:endp)) ; this%deadcrootc_SASUsave_patch          (:) = nan
       allocate(this%deadcrootc_storage_SASUsave_patch   (begp:endp)) ; this%deadcrootc_storage_SASUsave_patch  (:) = nan
       allocate(this%deadcrootc_xfer_SASUsave_patch      (begp:endp)) ; this%deadcrootc_xfer_SASUsave_patch     (:) = nan
       allocate(this%grainc_SASUsave_patch               (begp:endp)) ; this%grainc_SASUsave_patch              (:) = nan
       allocate(this%grainc_storage_SASUsave_patch       (begp:endp)) ; this%grainc_storage_SASUsave_patch      (:) = nan

       allocate(this%matrix_calloc_leaf_acc_patch        (begp:endp)); this%matrix_calloc_leaf_acc_patch        (:) = nan
       allocate(this%matrix_calloc_leafst_acc_patch      (begp:endp)); this%matrix_calloc_leafst_acc_patch      (:) = nan
       allocate(this%matrix_calloc_froot_acc_patch       (begp:endp)); this%matrix_calloc_froot_acc_patch       (:) = nan
       allocate(this%matrix_calloc_frootst_acc_patch     (begp:endp)); this%matrix_calloc_frootst_acc_patch     (:) = nan
       allocate(this%matrix_calloc_livestem_acc_patch    (begp:endp)); this%matrix_calloc_livestem_acc_patch    (:) = nan
       allocate(this%matrix_calloc_livestemst_acc_patch  (begp:endp)); this%matrix_calloc_livestemst_acc_patch  (:) = nan
       allocate(this%matrix_calloc_deadstem_acc_patch    (begp:endp)); this%matrix_calloc_deadstem_acc_patch    (:) = nan
       allocate(this%matrix_calloc_deadstemst_acc_patch  (begp:endp)); this%matrix_calloc_deadstemst_acc_patch  (:) = nan
       allocate(this%matrix_calloc_livecroot_acc_patch   (begp:endp)); this%matrix_calloc_livecroot_acc_patch   (:) = nan
       allocate(this%matrix_calloc_livecrootst_acc_patch (begp:endp)); this%matrix_calloc_livecrootst_acc_patch (:) = nan
       allocate(this%matrix_calloc_deadcroot_acc_patch   (begp:endp)); this%matrix_calloc_deadcroot_acc_patch   (:) = nan
       allocate(this%matrix_calloc_deadcrootst_acc_patch (begp:endp)); this%matrix_calloc_deadcrootst_acc_patch (:) = nan
       allocate(this%matrix_calloc_grain_acc_patch       (begp:endp)); this%matrix_calloc_grain_acc_patch       (:) = nan
       allocate(this%matrix_calloc_grainst_acc_patch     (begp:endp)); this%matrix_calloc_grainst_acc_patch     (:) = nan

       allocate(this%matrix_ctransfer_leafst_to_leafxf_acc_patch           (begp:endp))
       this%matrix_ctransfer_leafst_to_leafxf_acc_patch                    (:) = nan
       allocate(this%matrix_ctransfer_leafxf_to_leaf_acc_patch             (begp:endp))
       this%matrix_ctransfer_leafxf_to_leaf_acc_patch                      (:) = nan
       allocate(this%matrix_ctransfer_frootst_to_frootxf_acc_patch         (begp:endp))
       this%matrix_ctransfer_frootst_to_frootxf_acc_patch                  (:) = nan
       allocate(this%matrix_ctransfer_frootxf_to_froot_acc_patch           (begp:endp))
       this%matrix_ctransfer_frootxf_to_froot_acc_patch                    (:) = nan
       allocate(this%matrix_ctransfer_livestemst_to_livestemxf_acc_patch   (begp:endp))
       this%matrix_ctransfer_livestemst_to_livestemxf_acc_patch            (:) = nan
       allocate(this%matrix_ctransfer_livestemxf_to_livestem_acc_patch     (begp:endp))
       this%matrix_ctransfer_livestemxf_to_livestem_acc_patch              (:) = nan
       allocate(this%matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch   (begp:endp))
       this%matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch            (:) = nan
       allocate(this%matrix_ctransfer_deadstemxf_to_deadstem_acc_patch     (begp:endp))
       this%matrix_ctransfer_deadstemxf_to_deadstem_acc_patch              (:) = nan
       allocate(this%matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch (begp:endp))
       this%matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch          (:) = nan
       allocate(this%matrix_ctransfer_livecrootxf_to_livecroot_acc_patch   (begp:endp))
       this%matrix_ctransfer_livecrootxf_to_livecroot_acc_patch            (:) = nan
       allocate(this%matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch (begp:endp))
       this%matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch          (:) = nan
       allocate(this%matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch   (begp:endp))
       this%matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch            (:) = nan
       allocate(this%matrix_ctransfer_grainst_to_grainxf_acc_patch         (begp:endp))
       this%matrix_ctransfer_grainst_to_grainxf_acc_patch                  (:) = nan
       allocate(this%matrix_ctransfer_grainxf_to_grain_acc_patch           (begp:endp))
       this%matrix_ctransfer_grainxf_to_grain_acc_patch                    (:) = nan
       allocate(this%matrix_ctransfer_livestem_to_deadstem_acc_patch       (begp:endp))
       this%matrix_ctransfer_livestem_to_deadstem_acc_patch                (:) = nan
       allocate(this%matrix_ctransfer_livecroot_to_deadcroot_acc_patch     (begp:endp))
       this%matrix_ctransfer_livecroot_to_deadcroot_acc_patch              (:) = nan

       allocate(this%matrix_cturnover_leaf_acc_patch        (begp:endp)) ; this%matrix_cturnover_leaf_acc_patch        (:) = nan
       allocate(this%matrix_cturnover_leafst_acc_patch      (begp:endp)) ; this%matrix_cturnover_leafst_acc_patch      (:) = nan
       allocate(this%matrix_cturnover_leafxf_acc_patch      (begp:endp)) ; this%matrix_cturnover_leafxf_acc_patch      (:) = nan
       allocate(this%matrix_cturnover_froot_acc_patch       (begp:endp)) ; this%matrix_cturnover_froot_acc_patch       (:) = nan
       allocate(this%matrix_cturnover_frootst_acc_patch     (begp:endp)) ; this%matrix_cturnover_frootst_acc_patch     (:) = nan
       allocate(this%matrix_cturnover_frootxf_acc_patch     (begp:endp)) ; this%matrix_cturnover_frootxf_acc_patch     (:) = nan
       allocate(this%matrix_cturnover_livestem_acc_patch    (begp:endp)) ; this%matrix_cturnover_livestem_acc_patch    (:) = nan
       allocate(this%matrix_cturnover_livestemst_acc_patch  (begp:endp)) ; this%matrix_cturnover_livestemst_acc_patch  (:) = nan
       allocate(this%matrix_cturnover_livestemxf_acc_patch  (begp:endp)) ; this%matrix_cturnover_livestemxf_acc_patch  (:) = nan
       allocate(this%matrix_cturnover_deadstem_acc_patch    (begp:endp)) ; this%matrix_cturnover_deadstem_acc_patch    (:) = nan
       allocate(this%matrix_cturnover_deadstemst_acc_patch  (begp:endp)) ; this%matrix_cturnover_deadstemst_acc_patch  (:) = nan
       allocate(this%matrix_cturnover_deadstemxf_acc_patch  (begp:endp)) ; this%matrix_cturnover_deadstemxf_acc_patch  (:) = nan
       allocate(this%matrix_cturnover_livecroot_acc_patch   (begp:endp)) ; this%matrix_cturnover_livecroot_acc_patch   (:) = nan
       allocate(this%matrix_cturnover_livecrootst_acc_patch (begp:endp)) ; this%matrix_cturnover_livecrootst_acc_patch (:) = nan
       allocate(this%matrix_cturnover_livecrootxf_acc_patch (begp:endp)) ; this%matrix_cturnover_livecrootxf_acc_patch (:) = nan
       allocate(this%matrix_cturnover_deadcroot_acc_patch   (begp:endp)) ; this%matrix_cturnover_deadcroot_acc_patch   (:) = nan
       allocate(this%matrix_cturnover_deadcrootst_acc_patch (begp:endp)) ; this%matrix_cturnover_deadcrootst_acc_patch (:) = nan
       allocate(this%matrix_cturnover_deadcrootxf_acc_patch (begp:endp)) ; this%matrix_cturnover_deadcrootxf_acc_patch (:) = nan
       allocate(this%matrix_cturnover_grain_acc_patch       (begp:endp)) ; this%matrix_cturnover_grain_acc_patch       (:) = nan
       allocate(this%matrix_cturnover_grainst_acc_patch     (begp:endp)) ; this%matrix_cturnover_grainst_acc_patch     (:) = nan
       allocate(this%matrix_cturnover_grainxf_acc_patch     (begp:endp)) ; this%matrix_cturnover_grainxf_acc_patch     (:) = nan
    end if

    allocate(this%cropseedc_deficit_patch  (begp:endp)) ; this%cropseedc_deficit_patch  (:) = nan
    allocate(this%seedc_grc                (begg:endg)) ; this%seedc_grc                (:) = nan
    allocate(this%rootc_col                (begc:endc)) ; this%rootc_col                (:) = nan
    allocate(this%leafc_col                (begc:endc)) ; this%leafc_col                (:) = nan
    allocate(this%deadstemc_col            (begc:endc)) ; this%deadstemc_col            (:) = nan
    allocate(this%fuelc_col                (begc:endc)) ; this%fuelc_col                (:) = nan
    allocate(this%fuelc_crop_col           (begc:endc)) ; this%fuelc_crop_col           (:) = nan

    allocate(this%totvegc_patch            (begp:endp)) ; this%totvegc_patch            (:) = nan
    allocate(this%totvegc_col              (begc:endc)) ; this%totvegc_col              (:) = nan

    allocate(this%totc_p2c_col             (begc:endc)) ; this%totc_p2c_col             (:) = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds, carbon_type)
    !
    ! !DESCRIPTION:
    ! add history fields for all CN variables, always set as default='inactive'
    !
    ! !USES:
    use clm_varctl , only : use_c13, use_c14
    use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp 
    !
    ! !ARGUMENTS:
    class (cnveg_carbonstate_type) :: this
    type(bounds_type)         , intent(in) :: bounds 
    character(len=*)          , intent(in) :: carbon_type ! one of ['c12', c13','c14']
    !
    ! !LOCAL VARIABLES:
    integer           :: k,l,ii,jj 
    character(10)     :: active
    integer           :: begp,endp
    integer           :: begc,endc
    integer           :: begg,endg 
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    real(r8), pointer :: data2dptr(:,:) ! temp. pointer for slicing larger arrays
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    !-------------------------------
    ! C12 state variables
    !-------------------------------

    if (carbon_type == 'c12') then

       if (use_crop) then
          this%reproductivec_patch(begp:endp,:) = spval
          do k = 1, nrepr
             data1dptr => this%reproductivec_patch(:,k)
             call hist_addfld1d ( &
                  ! e.g., GRAINC
                  fname=get_repr_hist_fname(k)//'C', &
                  units='gC/m^2', &
                  avgflag='A', &
                  long_name=get_repr_longname(k)//' C (does not equal yield)', &
                  ptr_patch=data1dptr)
          end do

          this%cropseedc_deficit_patch(begp:endp) = spval
          call hist_addfld1d (fname='CROPSEEDC_DEFICIT', units='gC/m^2', &
               avgflag='A', long_name='C used for crop seed that needs to be repaid', &
               ptr_patch=this%cropseedc_deficit_patch)

          this%xsmrpool_loss_patch(begp:endp) = spval
          call hist_addfld1d (fname='XSMRPOOL_LOSS', units='gC/m^2', &
               avgflag='A', long_name='temporary photosynthate C pool loss', &
               ptr_patch=this%xsmrpool_loss_patch, default='inactive')
       end if
       
       this%woodc_patch(begp:endp) = spval
       call hist_addfld1d (fname='WOODC', units='gC/m^2', &
            avgflag='A', long_name='wood C', &
            ptr_patch=this%woodc_patch)

       this%leafc_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC', units='gC/m^2', &
            avgflag='A', long_name='leaf C', &
            ptr_patch=this%leafc_patch)

       this%leafc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='leaf C storage', &
            ptr_patch=this%leafc_storage_patch, default='inactive')    

       this%leafc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_XFER', units='gC/m^2', &
            avgflag='A', long_name='leaf C transfer', &
            ptr_patch=this%leafc_xfer_patch, default='inactive')    

       if(use_matrixcn)then
          this%matrix_cap_leafc_patch(begp:endp) = spval
          call hist_addfld1d (fname='LEAFC_CAP', units='gC/m^2', &
               avgflag='I', long_name='leaf C capacity', &
               ptr_patch=this%matrix_cap_leafc_patch)

          this%matrix_cap_leafc_storage_patch(begp:endp) = spval
          call hist_addfld1d (fname='LEAFC_STORAGE_CAP', units='gC/m^2', &
               avgflag='I', long_name='leaf C storage capacity', &
               ptr_patch=this%matrix_cap_leafc_storage_patch, default='inactive')    

          this%matrix_cap_leafc_xfer_patch(begp:endp) = spval
          call hist_addfld1d (fname='LEAFC_XFER_CAP', units='gC/m^2', &
               avgflag='I', long_name='leaf C transfer capacity', &
               ptr_patch=this%matrix_cap_leafc_xfer_patch, default='inactive')    
       end if

       this%leafc_storage_xfer_acc_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_STORAGE_XFER_ACC', units='gC/m^2', &
            avgflag='A', long_name='Accumulated leaf C transfer', &
            ptr_patch=this%leafc_storage_xfer_acc_patch, default='inactive')

       this%storage_cdemand_patch(begp:endp) = spval
       call hist_addfld1d (fname='STORAGE_CDEMAND', units='gC/m^2', &
            avgflag='A', long_name='C use from the C storage pool', &
            ptr_patch=this%storage_cdemand_patch, default='inactive')

       this%frootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC', units='gC/m^2', &
            avgflag='A', long_name='fine root C', &
            ptr_patch=this%frootc_patch)

       this%frootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='fine root C storage', &
            ptr_patch=this%frootc_storage_patch, default='inactive')   

       this%frootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC_XFER', units='gC/m^2', &
            avgflag='A', long_name='fine root C transfer', &
            ptr_patch=this%frootc_xfer_patch, default='inactive')    

       if(use_matrixcn)then
          this%matrix_cap_frootc_patch(begp:endp) = spval
          call hist_addfld1d (fname='FROOTC_CAP', units='gC/m^2', &
               avgflag='I', long_name='fine root C capacity', &
               ptr_patch=this%matrix_cap_frootc_patch)

          this%matrix_cap_frootc_storage_patch(begp:endp) = spval
          call hist_addfld1d (fname='FROOTC_STORAGE_CAP', units='gC/m^2', &
               avgflag='I', long_name='fine root C storage capacity', &
               ptr_patch=this%matrix_cap_frootc_storage_patch, default='inactive')   

          this%matrix_cap_frootc_xfer_patch(begp:endp) = spval
          call hist_addfld1d (fname='FROOTC_XFER_CAP', units='gC/m^2', &
               avgflag='I', long_name='fine root C transfer capacity', &
               ptr_patch=this%matrix_cap_frootc_xfer_patch, default='inactive')    
       end if

       this%livestemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC', units='gC/m^2', &
            avgflag='A', long_name='live stem C', &
            ptr_patch=this%livestemc_patch)

       this%livestemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='live stem C storage', &
            ptr_patch=this%livestemc_storage_patch, default='inactive')    

       this%livestemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC_XFER', units='gC/m^2', &
            avgflag='A', long_name='live stem C transfer', &
            ptr_patch=this%livestemc_xfer_patch, default='inactive')     

       if(use_matrixcn)then
          this%matrix_cap_livestemc_patch(begp:endp) = spval
          call hist_addfld1d (fname='LIVESTEMC_CAP', units='gC/m^2', &
               avgflag='I', long_name='live stem C capacity', &
               ptr_patch=this%matrix_cap_livestemc_patch)

          this%matrix_cap_livestemc_storage_patch(begp:endp) = spval
          call hist_addfld1d (fname='LIVESTEMC_STORAGE_CAP', units='gC/m^2', &
               avgflag='I', long_name='live stem C storage capcity', &
               ptr_patch=this%matrix_cap_livestemc_storage_patch, default='inactive')    

          this%matrix_cap_livestemc_xfer_patch(begp:endp) = spval
          call hist_addfld1d (fname='LIVESTEMC_XFER_CAP', units='gC/m^2', &
               avgflag='I', long_name='live stem C transfer capacity', &
               ptr_patch=this%matrix_cap_livestemc_xfer_patch, default='inactive')     
       end if

       this%deadstemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC', units='gC/m^2', &
            avgflag='A', long_name='dead stem C', &
            ptr_patch=this%deadstemc_patch)

       this%deadstemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='dead stem C storage', &
            ptr_patch=this%deadstemc_storage_patch, default='inactive')    

       this%deadstemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC_XFER', units='gC/m^2', &
            avgflag='A', long_name='dead stem C transfer', &
            ptr_patch=this%deadstemc_xfer_patch, default='inactive')    

       if(use_matrixcn)then
          this%matrix_cap_deadstemc_patch(begp:endp) = spval
          call hist_addfld1d (fname='DEADSTEMC_CAP', units='gC/m^2', &
               avgflag='I', long_name='dead stem C capacity', &
               ptr_patch=this%matrix_cap_deadstemc_patch)

          this%matrix_cap_deadstemc_storage_patch(begp:endp) = spval
          call hist_addfld1d (fname='DEADSTEMC_STORAGE_CAP', units='gC/m^2', &
               avgflag='I', long_name='dead stem C storage capacity', &
               ptr_patch=this%matrix_cap_deadstemc_storage_patch, default='inactive')    

          this%matrix_cap_deadstemc_xfer_patch(begp:endp) = spval
          call hist_addfld1d (fname='DEADSTEMC_XFER_CAP', units='gC/m^2', &
               avgflag='I', long_name='dead stem C transfer capacity', &
               ptr_patch=this%matrix_cap_deadstemc_xfer_patch, default='inactive')    
       end if

       this%livecrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC', units='gC/m^2', &
            avgflag='A', long_name='live coarse root C', &
            ptr_patch=this%livecrootc_patch)

       this%livecrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='live coarse root C storage', &
            ptr_patch=this%livecrootc_storage_patch, default='inactive')     

       this%livecrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC_XFER', units='gC/m^2', &
            avgflag='A', long_name='live coarse root C transfer', &
            ptr_patch=this%livecrootc_xfer_patch, default='inactive')    

       if(use_matrixcn)then
          this%matrix_cap_livecrootc_patch(begp:endp) = spval
          call hist_addfld1d (fname='LIVECROOTC_CAP', units='gC/m^2', &
               avgflag='I', long_name='live coarse root C capacity', &
               ptr_patch=this%matrix_cap_livecrootc_patch)

          this%matrix_cap_livecrootc_storage_patch(begp:endp) = spval
          call hist_addfld1d (fname='LIVECROOTC_STORAGE_CAP', units='gC/m^2', &
               avgflag='I', long_name='live coarse root C storage capacity', &
            ptr_patch=this%matrix_cap_livecrootc_storage_patch, default='inactive')     

          this%matrix_cap_livecrootc_xfer_patch(begp:endp) = spval
          call hist_addfld1d (fname='LIVECROOTC_XFER_CAP', units='gC/m^2', &
               avgflag='I', long_name='live coarse root C transfer capacity', &
               ptr_patch=this%matrix_cap_livecrootc_xfer_patch, default='inactive')    
       end if

       this%deadcrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC', units='gC/m^2', &
            avgflag='A', long_name='dead coarse root C', &
            ptr_patch=this%deadcrootc_patch)

       this%deadcrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='dead coarse root C storage', &
            ptr_patch=this%deadcrootc_storage_patch, default='inactive')   

       this%deadcrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC_XFER', units='gC/m^2', &
            avgflag='A', long_name='dead coarse root C transfer', &
            ptr_patch=this%deadcrootc_xfer_patch, default='inactive')   

       if(use_matrixcn)then
          this%matrix_cap_deadcrootc_patch(begp:endp) = spval
          call hist_addfld1d (fname='DEADCROOTC_CAP', units='gC/m^2', &
               avgflag='I', long_name='dead coarse root C capacity', &
               ptr_patch=this%matrix_cap_deadcrootc_patch)

          this%matrix_cap_deadcrootc_storage_patch(begp:endp) = spval
          call hist_addfld1d (fname='DEADCROOTC_STORAGE_CAP', units='gC/m^2', &
               avgflag='I', long_name='dead coarse root C storage capacity', &
               ptr_patch=this%matrix_cap_deadcrootc_storage_patch, default='inactive')   

          this%matrix_cap_deadcrootc_xfer_patch(begp:endp) = spval
          call hist_addfld1d (fname='DEADCROOTC_XFER_CAP', units='gC/m^2', &
               avgflag='I', long_name='dead coarse root C transfer capacity', &
               ptr_patch=this%matrix_cap_deadcrootc_xfer_patch, default='inactive')   
       end if

       this%gresp_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='GRESP_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='growth respiration storage', &
            ptr_patch=this%gresp_storage_patch, default='inactive')    

       this%gresp_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='GRESP_XFER', units='gC/m^2', &
            avgflag='A', long_name='growth respiration transfer', &
            ptr_patch=this%gresp_xfer_patch, default='inactive')     

       this%cpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL', units='gC/m^2', &
            avgflag='A', long_name='temporary photosynthate C pool', &
            ptr_patch=this%cpool_patch)

       this%xsmrpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='XSMRPOOL', units='gC/m^2', &
            avgflag='A', long_name='temporary photosynthate C pool', &
            ptr_patch=this%xsmrpool_patch)

       this%ctrunc_patch(begp:endp) = spval
       call hist_addfld1d (fname='PFT_CTRUNC', units='gC/m^2', &
            avgflag='A', long_name='patch-level sink for C truncation', &
            ptr_patch=this%ctrunc_patch, default='inactive')

       this%dispvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='DISPVEGC', units='gC/m^2', &
            avgflag='A', long_name='displayed veg carbon, excluding storage and cpool', &
            ptr_patch=this%dispvegc_patch)

       this%storvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='STORVEGC', units='gC/m^2', &
            avgflag='A', long_name='stored vegetation carbon, excluding cpool', &
            ptr_patch=this%storvegc_patch)

       this%totvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='TOTVEGC', units='gC/m^2', &
            avgflag='A', long_name='total vegetation carbon, excluding cpool', &
            ptr_patch=this%totvegc_patch)

       this%totc_patch(begp:endp) = spval
       call hist_addfld1d (fname='TOTPFTC', units='gC/m^2', &
            avgflag='A', long_name='total patch-level carbon, including cpool', &
            ptr_patch=this%totc_patch)

       this%seedc_grc(begg:endg) = spval
       call hist_addfld1d (fname='SEEDC', units='gC/m^2', &
            avgflag='A', long_name='pool for seeding new PFTs via dynamic landcover', &
            ptr_gcell=this%seedc_grc)

       this%fuelc_col(begc:endc) = spval
       call hist_addfld1d (fname='FUELC', units='gC/m^2', &
            avgflag='A', long_name='fuel load', &
            ptr_col=this%fuelc_col)

    end if

    !-------------------------------
    ! C13 state variables 
    !-------------------------------

    if ( carbon_type == 'c13' ) then

       this%leafc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC', units='gC13/m^2', &
            avgflag='A', long_name='C13 leaf C', &
            ptr_patch=this%leafc_patch, default='inactive')

       this%leafc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 leaf C storage', &
            ptr_patch=this%leafc_storage_patch, default='inactive')

       this%leafc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 leaf C transfer', &
            ptr_patch=this%leafc_xfer_patch, default='inactive')

       if(use_matrixcn)then
          this%matrix_cap_leafc_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_LEAFC_CAP', units='gC13/m^2', &
               avgflag='I', long_name='C13 leaf C capacity', &
               ptr_patch=this%matrix_cap_leafc_patch)

          this%matrix_cap_leafc_storage_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_LEAFC_STORAGE_CAP', units='gC13/m^2', &
               avgflag='I', long_name='C13 leaf C storage capacity', &
               ptr_patch=this%matrix_cap_leafc_storage_patch)!, default='inactive')    

          this%matrix_cap_leafc_xfer_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_LEAFC_XFER_CAP', units='gC13/m^2', &
               avgflag='I', long_name='C13 leaf C transfer capacity', &
               ptr_patch=this%matrix_cap_leafc_xfer_patch)!, default='inactive')    
       end if

       this%leafc_storage_xfer_acc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC_STORAGE_XFER_ACC', units='gC13/m^2', &
            avgflag='A', long_name='Accumulated C13 leaf C transfer', &
            ptr_patch=this%leafc_storage_xfer_acc_patch, default='inactive')

       this%frootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC', units='gC13/m^2', &
            avgflag='A', long_name='C13 fine root C', &
            ptr_patch=this%frootc_patch, default='inactive')

       this%frootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 fine root C storage', &
            ptr_patch=this%frootc_storage_patch, default='inactive')

       this%frootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 fine root C transfer', &
            ptr_patch=this%frootc_xfer_patch, default='inactive')

       if(use_matrixcn)then
          this%matrix_cap_frootc_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_FROOTC_CAP', units='gC13/m^2', &
               avgflag='I', long_name='C13 fine root C capacity', &
               ptr_patch=this%matrix_cap_frootc_patch)

          this%matrix_cap_frootc_storage_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_FROOTC_STORAGE_CAP', units='gC13/m^2', &
               avgflag='I', long_name='C13 fine root C storage capacity', &
               ptr_patch=this%matrix_cap_frootc_storage_patch)!, default='inactive')   

          this%matrix_cap_frootc_xfer_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_FROOTC_XFER_CAP', units='gC13/m^2', &
               avgflag='I', long_name='C13 fine root C transfer capacity', &
               ptr_patch=this%matrix_cap_frootc_xfer_patch)!, default='inactive')    
       end if

       this%livestemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC', units='gC13/m^2', &
            avgflag='A', long_name='C13 live stem C', &
            ptr_patch=this%livestemc_patch, default='inactive')

       this%livestemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 live stem C storage', &
            ptr_patch=this%livestemc_storage_patch, default='inactive')

       this%livestemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 live stem C transfer', &
            ptr_patch=this%livestemc_xfer_patch, default='inactive')

       if(use_matrixcn)then
          this%matrix_cap_livestemc_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_LIVESTEMC_CAP', units='gC13/m^2', &
               avgflag='I', long_name='C13 live stem C capacity', &
               ptr_patch=this%matrix_cap_livestemc_patch)

          this%matrix_cap_livestemc_storage_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_LIVESTEMC_STORAGE_CAP', units='gC13/m^2', &
               avgflag='I', long_name='C13 live stem C storage capcity', &
               ptr_patch=this%matrix_cap_livestemc_storage_patch)!, default='inactive')    

          this%matrix_cap_livestemc_xfer_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_LIVESTEMC_XFER_CAP', units='gC13/m^2', &
               avgflag='I', long_name='C13 live stem C transfer capacity', &
               ptr_patch=this%matrix_cap_livestemc_xfer_patch)!, default='inactive')     
       end if

       this%deadstemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADSTEMC', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead stem C', &
            ptr_patch=this%deadstemc_patch, default='inactive')

       this%deadstemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADSTEMC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead stem C storage', &
            ptr_patch=this%deadstemc_storage_patch, default='inactive')

       this%deadstemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADSTEMC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead stem C transfer', &
            ptr_patch=this%deadstemc_xfer_patch, default='inactive')

       if(use_matrixcn)then
          this%matrix_cap_deadstemc_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_DEADSTEMC_CAP', units='gC13/m^2', &
               avgflag='I', long_name='C13 dead stem C capacity', &
               ptr_patch=this%matrix_cap_deadstemc_patch)

          this%matrix_cap_deadstemc_storage_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_DEADSTEMC_STORAGE_CAP', units='gC13/m^2', &
               avgflag='I', long_name='C13 dead stem C storage capacity', &
               ptr_patch=this%matrix_cap_deadstemc_storage_patch)!, default='inactive')    

          this%matrix_cap_deadstemc_xfer_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_DEADSTEMC_XFER_CAP', units='gC13/m^2', &
               avgflag='I', long_name='C13 dead stem C transfer capacity', &
               ptr_patch=this%matrix_cap_deadstemc_xfer_patch)!, default='inactive')    
       end if

       this%livecrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC', units='gC13/m^2', &
            avgflag='A', long_name='C13 live coarse root C', &
            ptr_patch=this%livecrootc_patch, default='inactive')

       this%livecrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 live coarse root C storage', &
            ptr_patch=this%livecrootc_storage_patch, default='inactive')

       this%livecrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 live coarse root C transfer', &
            ptr_patch=this%livecrootc_xfer_patch, default='inactive')

       if(use_matrixcn)then
          this%matrix_cap_livecrootc_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_LIVECROOTC_CAP', units='gC13/m^2', &
               avgflag='I', long_name='C13 live coarse root C capacity', &
               ptr_patch=this%matrix_cap_livecrootc_patch)

          this%matrix_cap_livecrootc_storage_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_LIVECROOTC_STORAGE_CAP', units='gC13/m^2', &
               avgflag='I', long_name='C13 live coarse root C storage capacity', &
            ptr_patch=this%matrix_cap_livecrootc_storage_patch)!, default='inactive')     

          this%matrix_cap_livecrootc_xfer_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_LIVECROOTC_XFER_CAP', units='gC13/m^2', &
               avgflag='I', long_name='C13 live coarse root C transfer capacity', &
               ptr_patch=this%matrix_cap_livecrootc_xfer_patch)!, default='inactive')    
       end if

       this%deadcrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADCROOTC', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead coarse root C', &
            ptr_patch=this%deadcrootc_patch, default='inactive')

       this%deadcrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADCROOTC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead coarse root C storage', &
            ptr_patch=this%deadcrootc_storage_patch,  default='inactive')

       this%deadcrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADCROOTC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead coarse root C transfer', &
            ptr_patch=this%deadcrootc_xfer_patch, default='inactive')

       if(use_matrixcn)then
          this%matrix_cap_deadcrootc_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_DEADCROOTC_CAP', units='gC13/m^2', &
               avgflag='I', long_name='C13 dead coarse root C capacity', &
               ptr_patch=this%matrix_cap_deadcrootc_patch)

          this%matrix_cap_deadcrootc_storage_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_DEADCROOTC_STORAGE_CAP', units='gC13/m^2', &
               avgflag='I', long_name='C13 dead coarse root C storage capacity', &
               ptr_patch=this%matrix_cap_deadcrootc_storage_patch)!, default='inactive')   

          this%matrix_cap_deadcrootc_xfer_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_DEADCROOTC_XFER_CAP', units='gC13/m^2', &
               avgflag='I', long_name='C13 dead coarse root C transfer capacity', &
               ptr_patch=this%matrix_cap_deadcrootc_xfer_patch)!, default='inactive')   
       end if

       this%gresp_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_GRESP_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 growth respiration storage', &
            ptr_patch=this%gresp_storage_patch, default='inactive')

       this%gresp_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_GRESP_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 growth respiration transfer', &
            ptr_patch=this%gresp_xfer_patch, default='inactive')

       this%cpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL', units='gC13/m^2', &
            avgflag='A', long_name='C13 temporary photosynthate C pool', &
            ptr_patch=this%cpool_patch, default='inactive')

       this%xsmrpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_XSMRPOOL', units='gC13/m^2', &
            avgflag='A', long_name='C13 temporary photosynthate C pool', &
            ptr_patch=this%xsmrpool_patch, default='inactive')

       this%ctrunc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_PFT_CTRUNC', units='gC13/m^2', &
            avgflag='A', long_name='C13 patch-level sink for C truncation', &
            ptr_patch=this%ctrunc_patch, default='inactive')

       this%dispvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DISPVEGC', units='gC13/m^2', &
            avgflag='A', long_name='C13 displayed veg carbon, excluding storage and cpool', &
            ptr_patch=this%dispvegc_patch, default='inactive')

       this%storvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_STORVEGC', units='gC13/m^2', &
            avgflag='A', long_name='C13 stored vegetation carbon, excluding cpool', &
            ptr_patch=this%storvegc_patch, default='inactive')

       this%totvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_TOTVEGC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total vegetation carbon, excluding cpool', &
            ptr_patch=this%totvegc_patch)

       this%totc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_TOTPFTC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total patch-level carbon, including cpool', &
            ptr_patch=this%totc_patch, default='inactive')

       this%seedc_grc(begg:endg) = spval
       call hist_addfld1d (fname='C13_SEEDC', units='gC13/m^2', &
            avgflag='A', long_name='C13 pool for seeding new PFTs via dynamic landcover', &
            ptr_gcell=this%seedc_grc, default='inactive')

       if (use_crop) then
          this%reproductivec_patch(begp:endp,:) = spval
          do k = 1, nrepr
             data1dptr => this%reproductivec_patch(:,k)
             call hist_addfld1d ( &
                  ! e.g., C13_GRAINC
                  fname='C13_'//get_repr_hist_fname(k)//'C', &
                  units='gC/m^2', &
                  avgflag='A', &
                  long_name='C13 '//get_repr_longname(k)//' C (does not equal yield)', &
                  ptr_patch=data1dptr, default='inactive')
          end do

          this%cropseedc_deficit_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_CROPSEEDC_DEFICIT', units='gC/m^2', &
               avgflag='A', long_name='C13 C used for crop seed that needs to be repaid', &
               ptr_patch=this%cropseedc_deficit_patch, default='inactive')

          this%xsmrpool_loss_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_XSMRPOOL_LOSS', units='gC13/m^2', &
               avgflag='A', long_name='C13 temporary photosynthate C pool loss', &
               ptr_patch=this%xsmrpool_loss_patch, default='inactive')
       end if


    endif

    !-------------------------------
    ! C14 state variables 
    !-------------------------------

    if ( carbon_type == 'c14') then

       this%leafc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC', units='gC14/m^2', &
            avgflag='A', long_name='C14 leaf C', &
            ptr_patch=this%leafc_patch, default='inactive')

       this%leafc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 leaf C storage', &
            ptr_patch=this%leafc_storage_patch, default='inactive')

       this%leafc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 leaf C transfer', &
            ptr_patch=this%leafc_xfer_patch, default='inactive')

        this%leafc_storage_xfer_acc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_LEAFC_STORAGE_XFER_ACC', units='gC14/m^2', &
             avgflag='A', long_name='Accumulated C14 leaf C transfer', &
             ptr_patch=this%leafc_storage_xfer_acc_patch, default='inactive')

       if(use_matrixcn)then
          this%matrix_cap_leafc_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_LEAFC_CAP', units='gC14/m^2', &
               avgflag='I', long_name='C14 leaf C capacity', &
               ptr_patch=this%matrix_cap_leafc_patch)

          this%matrix_cap_leafc_storage_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_LEAFC_STORAGE_CAP', units='gC14/m^2', &
               avgflag='I', long_name='C14 leaf C storage capacity', &
               ptr_patch=this%matrix_cap_leafc_storage_patch)!, default='inactive')    

          this%matrix_cap_leafc_xfer_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_LEAFC_XFER_CAP', units='gC14/m^2', &
               avgflag='I', long_name='C14 leaf C transfer capacity', &
               ptr_patch=this%matrix_cap_leafc_xfer_patch)!, default='inactive')    
       end if

       this%frootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC', units='gC14/m^2', &
            avgflag='A', long_name='C14 fine root C', &
            ptr_patch=this%frootc_patch, default='inactive')

       this%frootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 fine root C storage', &
            ptr_patch=this%frootc_storage_patch, default='inactive')

       this%frootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 fine root C transfer', &
            ptr_patch=this%frootc_xfer_patch, default='inactive')

       if(use_matrixcn)then
          this%matrix_cap_frootc_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_FROOTC_CAP', units='gC14/m^2', &
               avgflag='I', long_name='C14 fine root C capacity', &
               ptr_patch=this%matrix_cap_frootc_patch)

          this%matrix_cap_frootc_storage_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_FROOTC_STORAGE_CAP', units='gC14/m^2', &
               avgflag='I', long_name='C14 fine root C storage capacity', &
               ptr_patch=this%matrix_cap_frootc_storage_patch)!, default='inactive')   

          this%matrix_cap_frootc_xfer_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_FROOTC_XFER_CAP', units='gC14/m^2', &
               avgflag='I', long_name='C14 fine root C transfer capacity', &
               ptr_patch=this%matrix_cap_frootc_xfer_patch)!, default='inactive')    
       end if

       this%livestemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC', units='gC14/m^2', &
            avgflag='A', long_name='C14 live stem C', &
            ptr_patch=this%livestemc_patch, default='inactive')

       this%livestemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 live stem C storage', &
            ptr_patch=this%livestemc_storage_patch, default='inactive')

       this%livestemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 live stem C transfer', &
            ptr_patch=this%livestemc_xfer_patch, default='inactive')

       if(use_matrixcn)then
          this%matrix_cap_livestemc_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_LIVESTEMC_CAP', units='gC14/m^2', &
               avgflag='I', long_name='C14 live stem C capacity', &
               ptr_patch=this%matrix_cap_livestemc_patch)

          this%matrix_cap_livestemc_storage_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_LIVESTEMC_STORAGE_CAP', units='gC14/m^2', &
               avgflag='I', long_name='C14 live stem C storage capcity', &
               ptr_patch=this%matrix_cap_livestemc_storage_patch)!, default='inactive')    

          this%matrix_cap_livestemc_xfer_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_LIVESTEMC_XFER_CAP', units='gC14/m^2', &
               avgflag='I', long_name='C14 live stem C transfer capacity', &
               ptr_patch=this%matrix_cap_livestemc_xfer_patch)!, default='inactive')     
       end if

       this%deadstemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADSTEMC', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead stem C', &
            ptr_patch=this%deadstemc_patch, default='inactive')

       this%deadstemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADSTEMC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead stem C storage', &
            ptr_patch=this%deadstemc_storage_patch, default='inactive')

       this%deadstemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADSTEMC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead stem C transfer', &
            ptr_patch=this%deadstemc_xfer_patch, default='inactive')

       if(use_matrixcn)then
          this%matrix_cap_deadstemc_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_DEADSTEMC_CAP', units='gC14/m^2', &
               avgflag='I', long_name='C14 dead stem C capacity', &
               ptr_patch=this%matrix_cap_deadstemc_patch)

          this%matrix_cap_deadstemc_storage_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_DEADSTEMC_STORAGE_CAP', units='gC14/m^2', &
               avgflag='I', long_name='C14 dead stem C storage capacity', &
               ptr_patch=this%matrix_cap_deadstemc_storage_patch)!, default='inactive')    

          this%matrix_cap_deadstemc_xfer_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_DEADSTEMC_XFER_CAP', units='gC14/m^2', &
               avgflag='I', long_name='C14 dead stem C transfer capacity', &
               ptr_patch=this%matrix_cap_deadstemc_xfer_patch)!, default='inactive')    
       end if

       this%livecrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC', units='gC14/m^2', &
            avgflag='A', long_name='C14 live coarse root C', &
            ptr_patch=this%livecrootc_patch, default='inactive')

       this%livecrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 live coarse root C storage', &
            ptr_patch=this%livecrootc_storage_patch, default='inactive')

       this%livecrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 live coarse root C transfer', &
            ptr_patch=this%livecrootc_xfer_patch, default='inactive')

       if(use_matrixcn)then
          this%matrix_cap_livecrootc_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_LIVECROOTC_CAP', units='gC14/m^2', &
               avgflag='I', long_name='C14 live coarse root C capacity', &
               ptr_patch=this%matrix_cap_livecrootc_patch)

          this%matrix_cap_livecrootc_storage_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_LIVECROOTC_STORAGE_CAP', units='gC14/m^2', &
               avgflag='I', long_name='C14 live coarse root C storage capacity', &
            ptr_patch=this%matrix_cap_livecrootc_storage_patch)!, default='inactive')     

          this%matrix_cap_livecrootc_xfer_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_LIVECROOTC_XFER_CAP', units='gC14/m^2', &
               avgflag='I', long_name='C14 live coarse root C transfer capacity', &
               ptr_patch=this%matrix_cap_livecrootc_xfer_patch)!, default='inactive')    
       end if

       this%deadcrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADCROOTC', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead coarse root C', &
            ptr_patch=this%deadcrootc_patch, default='inactive')

       this%deadcrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADCROOTC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead coarse root C storage', &
            ptr_patch=this%deadcrootc_storage_patch,  default='inactive')

       this%deadcrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADCROOTC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead coarse root C transfer', &
            ptr_patch=this%deadcrootc_xfer_patch, default='inactive')

       if(use_matrixcn)then
          this%matrix_cap_deadcrootc_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_DEADCROOTC_CAP', units='gC14/m^2', &
               avgflag='I', long_name='C14 dead coarse root C capacity', &
               ptr_patch=this%matrix_cap_deadcrootc_patch)

          this%matrix_cap_deadcrootc_storage_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_DEADCROOTC_STORAGE_CAP', units='gC14/m^2', &
               avgflag='I', long_name='C14 dead coarse root C storage capacity', &
               ptr_patch=this%matrix_cap_deadcrootc_storage_patch)!, default='inactive')   

          this%matrix_cap_deadcrootc_xfer_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_DEADCROOTC_XFER_CAP', units='gC14/m^2', &
               avgflag='I', long_name='C14 dead coarse root C transfer capacity', &
               ptr_patch=this%matrix_cap_deadcrootc_xfer_patch)!, default='inactive')   
       end if

       this%gresp_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_GRESP_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 growth respiration storage', &
            ptr_patch=this%gresp_storage_patch, default='inactive')

       this%gresp_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_GRESP_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 growth respiration transfer', &
            ptr_patch=this%gresp_xfer_patch, default='inactive')

       this%cpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL', units='gC14/m^2', &
            avgflag='A', long_name='C14 temporary photosynthate C pool', &
            ptr_patch=this%cpool_patch, default='inactive')

       this%xsmrpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_XSMRPOOL', units='gC14/m^2', &
            avgflag='A', long_name='C14 temporary photosynthate C pool', &
            ptr_patch=this%xsmrpool_patch, default='inactive')

       this%ctrunc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_PFT_CTRUNC', units='gC14/m^2', &
            avgflag='A', long_name='C14 patch-level sink for C truncation', &
            ptr_patch=this%ctrunc_patch, default='inactive')

       this%dispvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DISPVEGC', units='gC14/m^2', &
            avgflag='A', long_name='C14 displayed veg carbon, excluding storage and cpool', &
            ptr_patch=this%dispvegc_patch, default='inactive')

       this%storvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_STORVEGC', units='gC14/m^2', &
            avgflag='A', long_name='C14 stored vegetation carbon, excluding cpool', &
            ptr_patch=this%storvegc_patch, default='inactive')

       this%totvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_TOTVEGC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total vegetation carbon, excluding cpool', &
            ptr_patch=this%totvegc_patch)

       this%totc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_TOTPFTC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total patch-level carbon, including cpool', &
            ptr_patch=this%totc_patch, default='inactive')

       this%seedc_grc(begg:endg) = spval
       call hist_addfld1d (fname='C14_SEEDC', units='gC14/m^2', &
            avgflag='A', long_name='C14 pool for seeding new PFTs via dynamic landcover', &
            ptr_gcell=this%seedc_grc, default='inactive')

       if (use_crop) then
          this%reproductivec_patch(begp:endp,:) = spval
          do k = 1, nrepr
             data1dptr => this%reproductivec_patch(:,k)
             call hist_addfld1d ( &
                  ! e.g., C14_GRAINC
                  fname='C14_'//get_repr_hist_fname(k)//'C', units='gC/m^2', &
                  avgflag='A', &
                  long_name='C14 '//get_repr_longname(k)//' C (does not equal yield)', &
                  ptr_patch=data1dptr, default='inactive')
          end do

          this%cropseedc_deficit_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_CROPSEEDC_DEFICIT', units='gC/m^2', &
               avgflag='A', long_name='C14 C used for crop seed that needs to be repaid', &
               ptr_patch=this%cropseedc_deficit_patch, default='inactive')

          this%xsmrpool_loss_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_XSMRPOOL_LOSS', units='gC14/m^2', &
               avgflag='A', long_name='C14 temporary photosynthate C pool loss', &
               ptr_patch=this%xsmrpool_loss_patch, default='inactive')
       end if


    endif

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, ratio, carbon_type, c12_cnveg_carbonstate_inst)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    ! !USES:
    use landunit_varcon  , only : istsoil, istcrop
    use clm_time_manager , only : is_restart, get_nstep
    use clm_varctl, only : MM_Nuptake_opt, spinup_state
    !
    ! !ARGUMENTS:
    class(cnveg_carbonstate_type)                       :: this 
    type(bounds_type)            , intent(in)           :: bounds  
    real(r8)                     , intent(in)           :: ratio              ! Standard isotope ratio
    character(len=*)             , intent(in)           :: carbon_type        ! 'c12' or 'c13' or 'c14'
    type(cnveg_carbonstate_type) , optional, intent(in) :: c12_cnveg_carbonstate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,l,g,j,k,i
    integer  :: fc                                       ! filter index
    integer  :: num_special_col                          ! number of good values in special_col filter
    integer  :: num_special_patch                        ! number of good values in special_patch filter
    integer  :: special_col(bounds%endc-bounds%begc+1)   ! special landunit filter - columns
    integer  :: special_patch(bounds%endp-bounds%begp+1) ! special landunit filter - patches
    !-----------------------------------------------------------------------

    if (carbon_type == 'c13' .or. carbon_type == 'c14') then
       if (.not. present(c12_cnveg_carbonstate_inst)) then
          call endrun(msg=' ERROR: for C13 or C14 must pass in c12_cnveg_carbonstate_inst as argument' //&
               errMsg(sourcefile, __LINE__))
       end if
    else
       if ( spinup_state == 2 ) spinup_factor_deadwood = spinup_factor_AD
    end if

    ! Set column filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    ! Set patch filters

    num_special_patch = 0
    do p = bounds%begp,bounds%endp
       l = patch%landunit(p)
       if (lun%ifspecial(l)) then
          num_special_patch = num_special_patch + 1
          special_patch(num_special_patch) = p
       end if
    end do

    !-----------------------------------------------
    ! initialize patch-level carbon state variables
    !-----------------------------------------------

    do p = bounds%begp,bounds%endp

       this%leafcmax_patch(p) = 0._r8

       l = patch%landunit(p)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

          if (patch%itype(p) == noveg) then
             this%leafc_patch(p)                        = 0._r8
             this%leafc_storage_patch(p)                = 0._r8
             this%frootc_patch(p)                       = 0._r8            
             this%frootc_storage_patch(p)               = 0._r8    
             if(use_matrixcn)then
                this%matrix_cap_leafc_patch(p)          = 0._r8
                this%matrix_cap_leafc_storage_patch(p)  = 0._r8
                this%matrix_cap_frootc_patch(p)         = 0._r8            
                this%matrix_cap_frootc_storage_patch(p) = 0._r8    
             end if
          else
             if (pftcon%evergreen(patch%itype(p)) == 1._r8) then
                this%leafc_patch(p)                        = cnvegcstate_const%initial_vegC * ratio     
                this%leafc_storage_patch(p)                = 0._r8
                this%frootc_patch(p)                       = cnvegcstate_const%initial_vegC * ratio           
                this%frootc_storage_patch(p)               = 0._r8    
                if(use_matrixcn)then
                   this%matrix_cap_leafc_patch(p)          = cnvegcstate_const%initial_vegC * ratio     
                   this%matrix_cap_leafc_storage_patch(p)  = 0._r8
                   this%matrix_cap_frootc_patch(p)         = cnvegcstate_const%initial_vegC * ratio           
                   this%matrix_cap_frootc_storage_patch(p) = 0._r8    
                end if
             else if (patch%itype(p) >= npcropmin) then ! prognostic crop types
                this%leafc_patch(p)                        = 0._r8
                this%leafc_storage_patch(p)                = 0._r8
                this%frootc_patch(p)                       = 0._r8            
                this%frootc_storage_patch(p)               = 0._r8    
                if(use_matrixcn)then
                   this%matrix_cap_leafc_patch(p)          = 0._r8
                   this%matrix_cap_leafc_storage_patch(p)  = 0._r8
                   this%matrix_cap_frootc_patch(p)         = 0._r8            
                   this%matrix_cap_frootc_storage_patch(p) = 0._r8    
                end if
             else
                this%leafc_patch(p)                        = 0._r8
                this%leafc_storage_patch(p)                = cnvegcstate_const%initial_vegC * ratio   
                this%frootc_patch(p)                       = 0._r8            
                this%frootc_storage_patch(p)               = cnvegcstate_const%initial_vegC * ratio   
                if(use_matrixcn)then
                   this%matrix_cap_leafc_patch(p)          = 0._r8
                   this%matrix_cap_leafc_storage_patch(p)  = cnvegcstate_const%initial_vegC * ratio   
                   this%matrix_cap_frootc_patch(p)         = 0._r8            
                   this%matrix_cap_frootc_storage_patch(p) = cnvegcstate_const%initial_vegC * ratio   
                end if
             end if
          end if
          this%leafc_xfer_patch(p)                         = 0._r8
          if(use_matrixcn)then
             this%matrix_cap_leafc_xfer_patch(p)           = 0._r8
          end if
          this%leafc_storage_xfer_acc_patch(p)             = 0._r8
          this%storage_cdemand_patch(p)                    = 0._r8

          if (MM_Nuptake_opt .eqv. .false.) then  ! if not running in floating CN ratio option 
             this%frootc_patch(p)                          = 0._r8 
             this%frootc_storage_patch(p)                  = 0._r8 
             if(use_matrixcn)then
                this%matrix_cap_frootc_patch(p)            = 0._r8 
                this%matrix_cap_frootc_storage_patch(p)    = 0._r8 
             end if
          end if     
          this%frootc_xfer_patch(p)                        = 0._r8 
          if(use_matrixcn)then
             this%matrix_cap_frootc_xfer_patch(p)          = 0._r8 
          end if

          this%livestemc_patch(p)                          = 0._r8 
          this%livestemc_storage_patch(p)                  = 0._r8 
          this%livestemc_xfer_patch(p)                     = 0._r8 
          if(use_matrixcn)then
             this%matrix_cap_livestemc_patch(p)            = 0._r8 
             this%matrix_cap_livestemc_storage_patch(p)    = 0._r8 
             this%matrix_cap_livestemc_xfer_patch(p)       = 0._r8 
          end if

          if (pftcon%woody(patch%itype(p)) == 1._r8) then
             this%deadstemc_patch(p)                       = 0.1_r8 * ratio
             if(use_matrixcn)then
                this%matrix_cap_deadstemc_patch(p)         = 0.1_r8 * ratio
             end if
          else
             this%deadstemc_patch(p)                       = 0._r8 
             if(use_matrixcn)then
                this%matrix_cap_deadstemc_patch(p)         = 0._r8 
             end if
          end if
          this%deadstemc_storage_patch(p)                  = 0._r8 
          this%deadstemc_xfer_patch(p)                     = 0._r8 
          if(use_matrixcn)then
             this%matrix_cap_deadstemc_storage_patch(p)    = 0._r8 
             this%matrix_cap_deadstemc_xfer_patch(p)       = 0._r8 
          end if

          this%livecrootc_patch(p)                         = 0._r8 
          this%livecrootc_storage_patch(p)                 = 0._r8 
          this%livecrootc_xfer_patch(p)                    = 0._r8 
          if(use_matrixcn)then
             this%matrix_cap_livecrootc_patch(p)           = 0._r8 
             this%matrix_cap_livecrootc_storage_patch(p)   = 0._r8 
             this%matrix_cap_livecrootc_xfer_patch(p)      = 0._r8 
          end if

          this%deadcrootc_patch(p)                         = 0._r8 
          this%deadcrootc_storage_patch(p)                 = 0._r8 
          this%deadcrootc_xfer_patch(p)                    = 0._r8 
          if(use_matrixcn)then
             this%matrix_cap_deadcrootc_patch(p)           = 0._r8 
             this%matrix_cap_deadcrootc_storage_patch(p)   = 0._r8 
             this%matrix_cap_deadcrootc_xfer_patch(p)      = 0._r8 
          end if

          this%gresp_storage_patch(p)      = 0._r8 
          this%gresp_xfer_patch(p)         = 0._r8 

          this%cpool_patch(p)              = 0._r8 
          this%xsmrpool_patch(p)           = 0._r8 
          this%ctrunc_patch(p)             = 0._r8 
          this%dispvegc_patch(p)           = 0._r8 
          this%storvegc_patch(p)           = 0._r8 
          this%woodc_patch(p)              = 0._r8
          this%totc_patch(p)               = 0._r8 
!!!!initial pool size for matrix
          if(use_matrixcn)then
             this%leafc0_patch(p)              = 1.e-30_r8
             this%leafc0_storage_patch(p)      = 1.e-30_r8
             this%leafc0_xfer_patch(p)         = 1.e-30_r8
             this%frootc0_patch(p)             = 1.e-30_r8            
             this%frootc0_storage_patch(p)     = 1.e-30_r8  
             this%frootc0_xfer_patch(p)        = 1.e-30_r8 

             this%livestemc0_patch(p)          = 1.e-30_r8 
             this%livestemc0_storage_patch(p)  = 1.e-30_r8 
             this%livestemc0_xfer_patch(p)     = 1.e-30_r8
             this%deadstemc0_patch(p)          = 1.e-30_r8
             this%deadstemc0_storage_patch(p)  = 1.e-30_r8 
             this%deadstemc0_xfer_patch(p)     = 1.e-30_r8 

             this%livecrootc0_patch(p)         = 1.e-30_r8 
             this%livecrootc0_storage_patch(p) = 1.e-30_r8 
             this%livecrootc0_xfer_patch(p)    = 1.e-30_r8 

             this%deadcrootc0_patch(p)         = 1.e-30_r8 
             this%deadcrootc0_storage_patch(p) = 1.e-30_r8 
             this%deadcrootc0_xfer_patch(p)    = 1.e-30_r8

             this%reproc0_patch(p)             = 1.e-30_r8
             this%reproc0_storage_patch(p)     = 1.e-30_r8
             this%reproc0_xfer_patch(p)        = 1.e-30_r8

             this%leafc_SASUsave_patch(p)              = 0._r8
             this%leafc_storage_SASUsave_patch(p)      = 0._r8
             this%leafc_xfer_SASUsave_patch(p)         = 0._r8
             this%frootc_SASUsave_patch(p)             = 0._r8
             this%frootc_storage_SASUsave_patch(p)     = 0._r8
             this%frootc_xfer_SASUsave_patch(p)        = 0._r8
             this%livestemc_SASUsave_patch(p)          = 0._r8
             this%livestemc_storage_SASUsave_patch(p)  = 0._r8
             this%livestemc_xfer_SASUsave_patch(p)     = 0._r8
             this%deadstemc_SASUsave_patch(p)          = 0._r8
             this%deadstemc_storage_SASUsave_patch(p)  = 0._r8
             this%deadstemc_xfer_SASUsave_patch(p)     = 0._r8
             this%livecrootc_SASUsave_patch(p)         = 0._r8
             this%livecrootc_storage_SASUsave_patch(p) = 0._r8
             this%livecrootc_xfer_SASUsave_patch(p)    = 0._r8
             this%deadcrootc_SASUsave_patch(p)         = 0._r8
             this%deadcrootc_storage_SASUsave_patch(p) = 0._r8
             this%deadcrootc_xfer_SASUsave_patch(p)    = 0._r8
             this%grainc_SASUsave_patch(p)             = 0._r8
             this%grainc_storage_SASUsave_patch(p)     = 0._r8

             this%matrix_calloc_leaf_acc_patch(p)                           = 0._r8
             this%matrix_calloc_leafst_acc_patch(p)                         = 0._r8
             this%matrix_calloc_froot_acc_patch(p)                          = 0._r8
             this%matrix_calloc_frootst_acc_patch(p)                        = 0._r8
             this%matrix_calloc_livestem_acc_patch(p)                       = 0._r8
             this%matrix_calloc_livestemst_acc_patch(p)                     = 0._r8
             this%matrix_calloc_deadstem_acc_patch(p)                       = 0._r8
             this%matrix_calloc_deadstemst_acc_patch(p)                     = 0._r8
             this%matrix_calloc_livecroot_acc_patch(p)                      = 0._r8
             this%matrix_calloc_livecrootst_acc_patch(p)                    = 0._r8
             this%matrix_calloc_deadcroot_acc_patch(p)                      = 0._r8
             this%matrix_calloc_deadcrootst_acc_patch(p)                    = 0._r8

             this%matrix_ctransfer_leafst_to_leafxf_acc_patch(p)            = 0._r8
             this%matrix_ctransfer_leafxf_to_leaf_acc_patch(p)              = 0._r8
             this%matrix_ctransfer_frootst_to_frootxf_acc_patch(p)          = 0._r8
             this%matrix_ctransfer_frootxf_to_froot_acc_patch(p)            = 0._r8
             this%matrix_ctransfer_livestemst_to_livestemxf_acc_patch(p)    = 0._r8
             this%matrix_ctransfer_livestemxf_to_livestem_acc_patch(p)      = 0._r8
             this%matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch(p)    = 0._r8
             this%matrix_ctransfer_deadstemxf_to_deadstem_acc_patch(p)      = 0._r8
             this%matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch(p)  = 0._r8
             this%matrix_ctransfer_livecrootxf_to_livecroot_acc_patch(p)    = 0._r8
             this%matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch(p)  = 0._r8
             this%matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch(p)    = 0._r8
             this%matrix_ctransfer_livestem_to_deadstem_acc_patch(p)        = 0._r8
             this%matrix_ctransfer_livecroot_to_deadcroot_acc_patch(p)      = 0._r8

             this%matrix_cturnover_leaf_acc_patch(p)                        = 0._r8
             this%matrix_cturnover_leafst_acc_patch(p)                      = 0._r8
             this%matrix_cturnover_leafxf_acc_patch(p)                      = 0._r8
             this%matrix_cturnover_froot_acc_patch(p)                       = 0._r8
             this%matrix_cturnover_frootst_acc_patch(p)                     = 0._r8
             this%matrix_cturnover_frootxf_acc_patch(p)                     = 0._r8
             this%matrix_cturnover_livestem_acc_patch(p)                    = 0._r8
             this%matrix_cturnover_livestemst_acc_patch(p)                  = 0._r8
             this%matrix_cturnover_livestemxf_acc_patch(p)                  = 0._r8
             this%matrix_cturnover_deadstem_acc_patch(p)                    = 0._r8
             this%matrix_cturnover_deadstemst_acc_patch(p)                  = 0._r8
             this%matrix_cturnover_deadstemxf_acc_patch(p)                  = 0._r8
             this%matrix_cturnover_livecroot_acc_patch(p)                   = 0._r8
             this%matrix_cturnover_livecrootst_acc_patch(p)                 = 0._r8
             this%matrix_cturnover_livecrootxf_acc_patch(p)                 = 0._r8
             this%matrix_cturnover_deadcroot_acc_patch(p)                   = 0._r8
             this%matrix_cturnover_deadcrootst_acc_patch(p)                 = 0._r8
             this%matrix_cturnover_deadcrootxf_acc_patch(p)                 = 0._r8
          end if
  

          if ( use_crop )then
             this%reproductivec_patch(p,:)                                  = 0._r8
             this%reproductivec_storage_patch(p,:)                          = 0._r8
             this%reproductivec_xfer_patch(p,:)                             = 0._r8
             this%cropseedc_deficit_patch(p)                                = 0._r8
             this%xsmrpool_loss_patch(p)                                    = 0._r8 
             if(use_matrixcn)then
                this%matrix_cap_reproc_patch(p)                             = 0._r8            
                this%matrix_cap_reproc_storage_patch(p)                     = 0._r8    
                this%matrix_cap_reproc_xfer_patch(p)                        = 0._r8    
                ! I think these need to change as well...
                this%matrix_calloc_grain_acc_patch(p)                       = 0._r8            
                this%matrix_calloc_grainst_acc_patch(p)                     = 0._r8    
                this%matrix_ctransfer_grainst_to_grainxf_acc_patch(p)       = 0._r8
                this%matrix_ctransfer_grainxf_to_grain_acc_patch(p)         = 0._r8
                this%matrix_cturnover_grain_acc_patch(p)                    = 0._r8
                this%matrix_cturnover_grainst_acc_patch(p)                  = 0._r8
                this%matrix_cturnover_grainxf_acc_patch(p)                  = 0._r8
             end if
          end if

       endif

    end do

    ! -----------------------------------------------
    ! initialize column-level variables
    ! -----------------------------------------------

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
!          this%totgrainc_col(c)  = 0._r8

          ! total carbon pools
          this%totc_p2c_col(c)   = 0._r8
          
       end if
    end do


    do g = bounds%begg, bounds%endg
       this%seedc_grc(g) = 0._r8
    end do

    ! initialize fields for special filters

    call this%SetValues (&
         num_patch=num_special_patch, filter_patch=special_patch, value_patch=0._r8, &
         num_column=num_special_col, filter_column=special_col, value_column=0._r8)

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart ( this,  bounds, ncid, flag, carbon_type, reseed_dead_plants, &
                       c12_cnveg_carbonstate_inst, filter_reseed_patch, &
                       num_reseed_patch, spinup_factor4deadwood )
    !
    ! !DESCRIPTION: 
    ! Read/write CN restart data for carbon state
    !
    ! !USES:
    use shr_infnan_mod   , only : isnan => shr_infnan_isnan, nan => shr_infnan_nan, assignment(=)
    use clm_varcon       , only : c13ratio, c14ratio
    use clm_varctl       , only : spinup_state, use_cndv, MM_Nuptake_opt
    use clm_varctl       , only : spinup_state, use_cndv, MM_Nuptake_opt
    use clm_time_manager , only : is_restart
    use landunit_varcon  , only : istsoil, istcrop 
    use spmdMod          , only : mpicom
    use shr_mpi_mod      , only : shr_mpi_sum
    use restUtilMod
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class (cnveg_carbonstate_type)                               :: this
    type(bounds_type)                     , intent(in)           :: bounds 
    type(file_desc_t)                     , intent(inout)        :: ncid   ! netcdf id
    character(len=*)                      , intent(in)           :: flag   !'read' or 'write'
    character(len=*)                      , intent(in)           :: carbon_type ! 'c12' or 'c13' or 'c14'
    logical                               , intent(in)           :: reseed_dead_plants
    type (cnveg_carbonstate_type)         , intent(in), optional :: c12_cnveg_carbonstate_inst 
    integer                               , intent(out), optional :: filter_reseed_patch(:)
    integer                               , intent(out), optional :: num_reseed_patch
    real(r8)                              , intent(out), optional :: spinup_factor4deadwood
    !
    ! !LOCAL VARIABLES:
    integer            :: i,j,k,l,c,p
    real(r8)           :: ratio
    character(len=256) :: varname   ! temporary
    logical            :: readvar
    integer            :: idata
    logical            :: exit_spinup  = .false.
    logical            :: enter_spinup = .false.
    ! spinup state as read from restart file, for determining whether to enter or exit spinup mode.
    integer            :: restart_file_spinup_state
    integer            :: total_num_reseed_patch      ! Total number of patches to reseed across all processors
    real(r8), pointer  :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    real(r8), parameter:: totvegcthresh = 1.0_r8      ! Total vegetation carbon threshold to reseed dead vegetation

    logical            :: missing_ciso    ! whether C isotope fields are missing from the input file, despite the run containing C isotopes

    !------------------------------------------------------------------------

    missing_ciso = .false.

    if (carbon_type == 'c13' .or. carbon_type == 'c14') then
       if (.not. present(c12_cnveg_carbonstate_inst)) then
          call endrun(msg=' ERROR: for C14 must pass in c12_cnveg_carbonstate_inst as argument' //&
               errMsg(sourcefile, __LINE__))
       end if
    end if
    if (carbon_type == 'c12') then
       ratio = 1._r8
    else if (carbon_type == 'c13') then
       ratio = c13ratio
    else if (carbon_type == 'c14') then
       ratio = c14ratio
    end if

    if ( (      present(num_reseed_patch) .and. .not. present(filter_reseed_patch)) &
    .or. (.not. present(num_reseed_patch) .and.       present(filter_reseed_patch) ) )then
       call endrun(msg=' ERROR: filter_reseed_patch and num_reseed_patch both need to be entered ' //&
       errMsg(sourcefile, __LINE__))
    end if
    if ( present(num_reseed_patch) )then
       num_reseed_patch = 0
       filter_reseed_patch(:) = -1
    end if

    !--------------------------------
    ! patch carbon state variables (c12)
    !--------------------------------

    if (carbon_type == 'c12') then
       call restartvar(ncid=ncid, flag=flag, varname='leafc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='leafc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_xfer_patch) 

       if(use_matrixcn)then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_cap', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_leafc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_cap', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_leafc_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer_cap', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_leafc_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='leafc0', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc0_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='leafc0_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc0_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='leafc0_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc0_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_leaf_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_leaf_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_leafst_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_leafst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_leafst_to_leafxf_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_leafst_to_leafxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_leafxf_to_leaf_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_leafxf_to_leaf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_leaf_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_leaf_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_leafst_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_leafst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctrunover_leafxf_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_leafxf_acc_patch) 

       end if

       call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_xfer_acc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_storage_xfer_acc_patch)
 
       call restartvar(ncid=ncid, flag=flag, varname='storage_cdemand', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%storage_cdemand_patch)

       call restartvar(ncid=ncid, flag=flag, varname='frootc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='frootc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_xfer_patch) 

       if(use_matrixcn)then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_cap', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_frootc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='frootc_storage_cap', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_frootc_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer_cap', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_frootc_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='frootc0', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc0_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='frootc0_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc0_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='frootc0_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc0_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_froot_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_froot_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_frootst_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_frootst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_frootst_to_frootxf_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_frootst_to_frootxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_frootxf_to_froot_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_frootxf_to_froot_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_froot_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_froot_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_frootst_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_frootst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_frootxf_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_frootxf_acc_patch) 
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livestemc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_xfer_patch) 

       if(use_matrixcn)then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_cap', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_livestemc_patch) 
  
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage_cap', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_livestemc_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer_cap', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_livestemc_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livestemc0', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc0_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livestemc0_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc0_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livestemc0_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc0_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_livestem_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_livestem_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_livestemst_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_livestemst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_livestemst_to_livestemxf_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_livestemst_to_livestemxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_livestemxf_to_livestem_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_livestemxf_to_livestem_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_livestem_to_deadstem_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_livestem_to_deadstem_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_livestem_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_livestem_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_livestemst_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_livestemst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_livestemxf_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_livestemxf_acc_patch) 
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_xfer_patch) 

       if(use_matrixcn)then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_cap', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_deadstemc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage_cap', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_deadstemc_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer_cap', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_deadstemc_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc0', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc0_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc0_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc0_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc0_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc0_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_deadstem_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_deadstem_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_deadstemst_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_deadstemst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_deadstemst_to_deadstemxf_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_deadstemxf_to_deadstem_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_deadstemxf_to_deadstem_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_deadstem_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_deadstem_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_deadstemst_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_deadstemst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_deadstemxf_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_deadstemxf_acc_patch) 
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_xfer_patch) 

       if(use_matrixcn)then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_cap', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_livecrootc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage_cap', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_livecrootc_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer_cap', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_livecrootc_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc0', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc0_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc0_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc0_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc0_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc0_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_livecroot_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_livecroot_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_livecrootst_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_livecrootst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_livecrootst_to_livecrootxf_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_livecrootxf_to_livecroot_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_livecrootxf_to_livecroot_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_livecroot_to_deadcroot_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_livecroot_to_deadcroot_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_livecroot_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_livecroot_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_livecrootst_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_livecrootst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_livecrootxf_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_livecrootxf_acc_patch) 
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_xfer_patch) 

       if(use_matrixcn)then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_cap', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_deadcrootc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage_cap', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_deadcrootc_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer_cap', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_deadcrootc_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc0', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc0_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc0_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc0_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc0_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc0_xfer_patch) 
!
          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_deadcroot_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_deadcroot_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_deadcrootst_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_deadcrootst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_deadcrootst_to_deadcrootxf_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_deadcrootxf_to_deadcroot_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_deadcroot_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_deadcroot_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_deadcrootst_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_deadcrootst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_deadcrootxf_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_deadcrootxf_acc_patch) 
       end if

       if(use_matrixcn .and. use_crop)then
          call restartvar(ncid=ncid, flag=flag,  varname='reproc0', xtype=ncd_double,  &
               dim1name='pft',    long_name='initial grain C', units='gC/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%reproc0_patch)
 
          call restartvar(ncid=ncid, flag=flag,  varname='reproc0_storage', xtype=ncd_double,  &
               dim1name='pft',    long_name='initial grain C storage', units='gC/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%reproc0_storage_patch)

          call restartvar(ncid=ncid, flag=flag,  varname='reproc0_xfer', xtype=ncd_double,  &
               dim1name='pft',    long_name='initial grain C transfer', units='gC/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%reproc0_xfer_patch)

          call restartvar(ncid=ncid, flag=flag,  varname='matrix_calloc_grain_acc', xtype=ncd_double,  &
               dim1name='pft',    long_name='C accumulated allocation to grain', units='gC/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_grain_acc_patch)
 
          call restartvar(ncid=ncid, flag=flag,  varname='matrix_calloc_grainst_acc', xtype=ncd_double,  &
               dim1name='pft',    long_name='C accumulated allocation to grain storage', units='gC/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_grainst_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_grainst_to_grainxf_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_grainst_to_grainxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_grainxf_to_grain_acc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_grainxf_to_grain_acc_patch)

          call restartvar(ncid=ncid, flag=flag,  varname='matrix_cturnover_grain_acc', xtype=ncd_double,  &
               dim1name='pft',    long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_grain_acc_patch)
 
          call restartvar(ncid=ncid, flag=flag,  varname='matrix_cturnover_grainst_acc', xtype=ncd_double,  &
               dim1name='pft',    long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_grainst_acc_patch)
 
          call restartvar(ncid=ncid, flag=flag,  varname='matrix_cturnover_grainxf_acc', xtype=ncd_double,  &
               dim1name='pft',    long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_grainxf_acc_patch)
       end if

       call restartvar(ncid=ncid, flag=flag, varname='gresp_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%gresp_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%gresp_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='cpool', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%cpool_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='xsmrpool', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_patch) 

       if (use_crop) then
          call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_loss', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_loss_patch) 
          if (flag == 'read' .and. (.not. readvar) ) then
              this%xsmrpool_loss_patch(bounds%begp:bounds%endp) = 0._r8
          end if
       end if

       call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%ctrunc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='leafcmax', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafcmax_patch)

       if (flag == 'read') then
          call restartvar(ncid=ncid, flag=flag, varname='spinup_state', xtype=ncd_int, &
            long_name='Spinup state of the model that wrote this restart file: ' &
            // ' 0 = normal model mode, 1 = AD spinup, 2 = AAD spinup', units='', &
            interpinic_flag='copy', readvar=readvar,  data=idata)

          if (readvar) then
             restart_file_spinup_state = idata
          else
             restart_file_spinup_state = spinup_state
             if ( masterproc ) then
                write(iulog,*) ' CNRest: WARNING!  Restart file does not contain info ' &
                      // ' on spinup state used to generate the restart file. '
                 write(iulog,*) '   Assuming the same as current setting: ', spinup_state
             end if
          end if
       end if

       if (flag == 'read' .and. spinup_state /= restart_file_spinup_state .and. .not. use_cndv) then
          if ( masterproc ) write(iulog, *) 'exit_spinup ',exit_spinup,' restart_file_spinup_state ',restart_file_spinup_state
          if (spinup_state <= 1 .and. restart_file_spinup_state == 2 ) then
             if ( masterproc ) write(iulog,*) ' CNRest: taking Dead wood C pools out of AD spinup mode'
             exit_spinup = .true.
             if ( masterproc ) write(iulog, *) 'Multiplying stemc and crootc by ', spinup_factor_AD, ' for exit spinup'
             do i = bounds%begp,bounds%endp
                this%deadstemc_patch(i) = this%deadstemc_patch(i) * spinup_factor_AD
                this%deadcrootc_patch(i) = this%deadcrootc_patch(i) * spinup_factor_AD
             end do
          else if (spinup_state == 2 .and. restart_file_spinup_state <= 1 )then
             if (spinup_state == 2 .and. restart_file_spinup_state <= 1 )then
                if ( masterproc ) write(iulog,*) ' CNRest: taking Dead wood C pools into AD spinup mode'
                enter_spinup = .true.
                if ( masterproc ) write(iulog, *) 'Dividing stemc and crootc by ', spinup_factor_AD, 'for enter spinup '
                do i = bounds%begp,bounds%endp
                   this%deadstemc_patch(i) = this%deadstemc_patch(i) / spinup_factor_AD
                   this%deadcrootc_patch(i) = this%deadcrootc_patch(i) / spinup_factor_AD
                end do
             end if
          end if
       end if

       call restartvar(ncid=ncid, flag=flag, varname='totvegc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totvegc_patch) 
       ! totvegc_col needed for resetting soil carbon stocks during AD spinup exit
       call restartvar(ncid=ncid, flag=flag, varname='totvegc_col', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totvegc_col)


       if (  flag == 'read' .and. (enter_spinup .or. (reseed_dead_plants .and. .not. is_restart())) .and. .not. use_cndv) then
             if ( masterproc ) write(iulog, *) 'Reseeding dead plants for CNVegCarbonState'
             ! If a pft is dead or near-dead (indicated by totvegc <= totvegcthresh) then we reseed that
             ! pft according to the cold start protocol in the InitCold subroutine.
             ! Thus, the variable totvegc is required to be read before here
             ! so that if it is zero for a given pft, the pft can be reseeded.
             do i = bounds%begp,bounds%endp
                if (this%totvegc_patch(i) .le. totvegcthresh) then
                   !-----------------------------------------------
                   ! initialize patch-level carbon state variables
                   !-----------------------------------------------

                   this%leafcmax_patch(i) = 0._r8

                   l = patch%landunit(i)
                   if (lun%itype(l) == istsoil  .or. patch%itype(i) == nc3crop .or. patch%itype(i) == nc3irrig)then
                      if ( present(num_reseed_patch) ) then
                         num_reseed_patch = num_reseed_patch + 1
                         filter_reseed_patch(num_reseed_patch) = i
                      end if

                      if (patch%itype(i) == noveg) then
                         this%leafc_patch(i)                           = 0._r8
                         this%leafc_storage_patch(i)                   = 0._r8
                         this%frootc_patch(i)                          = 0._r8            
                         this%frootc_storage_patch(i)                  = 0._r8    
                         if(use_matrixcn)then
                            this%matrix_cap_leafc_patch(i)             = 0._r8
                            this%matrix_cap_leafc_storage_patch(i)     = 0._r8
                            this%matrix_cap_frootc_patch(i)            = 0._r8            
                            this%matrix_cap_frootc_storage_patch(i)    = 0._r8    
                         end if
                      else
                         if (pftcon%evergreen(patch%itype(i)) == 1._r8) then
                            this%leafc_patch(i)                        = cnvegcstate_const%initial_vegC * ratio     
                            this%leafc_storage_patch(i)                = 0._r8
                            this%frootc_patch(i)                       = cnvegcstate_const%initial_vegC * ratio           
                            this%frootc_storage_patch(i)               = 0._r8    
                            if(use_matrixcn)then
                               this%matrix_cap_leafc_patch(i)          = cnvegcstate_const%initial_vegC * ratio     
                               this%matrix_cap_leafc_storage_patch(i)  = 0._r8
                               this%matrix_cap_frootc_patch(i)         = cnvegcstate_const%initial_vegC * ratio           
                               this%matrix_cap_frootc_storage_patch(i) = 0._r8    
                            end if
                         else
                            this%leafc_patch(i)                        = 0._r8
                            this%leafc_storage_patch(i)                = cnvegcstate_const%initial_vegC * ratio   
                            this%frootc_patch(i)                       = 0._r8            
                            this%frootc_storage_patch(i)               = cnvegcstate_const%initial_vegC * ratio   
                            if(use_matrixcn)then
                               this%matrix_cap_leafc_patch(i)          = 0._r8
                               this%matrix_cap_leafc_storage_patch(i)  = cnvegcstate_const%initial_vegC * ratio   
                               this%matrix_cap_frootc_patch(i)         = 0._r8            
                               this%matrix_cap_frootc_storage_patch(i) = cnvegcstate_const%initial_vegC * ratio   
                            end if
                         end if
                      end if
                      this%leafc_xfer_patch(i)                         = 0._r8
                      if(use_matrixcn)then
                         this%matrix_cap_leafc_xfer_patch(i)           = 0._r8
                      end if
                      this%leafc_storage_xfer_acc_patch(i)             = 0._r8
                      this%storage_cdemand_patch(i)                    = 0._r8

                      if (MM_Nuptake_opt .eqv. .false.) then  ! if not running in floating CN ratio option 
                         this%frootc_patch(i)                          = 0._r8 
                         this%frootc_storage_patch(i)                  = 0._r8 
                         if(use_matrixcn)then
                            this%matrix_cap_frootc_patch(i)            = 0._r8 
                            this%matrix_cap_frootc_storage_patch(i)    = 0._r8 
                         end if
                      end if     
                      this%frootc_xfer_patch(i)                        = 0._r8 

                      this%livestemc_patch(i)                          = 0._r8 
                      this%livestemc_storage_patch(i)                  = 0._r8 
                      this%livestemc_xfer_patch(i)                     = 0._r8 
                      if(use_matrixcn)then
                         this%matrix_cap_frootc_xfer_patch(i)          = 0._r8 
                         this%matrix_cap_livestemc_patch(i)            = 0._r8 
                         this%matrix_cap_livestemc_storage_patch(i)    = 0._r8 
                         this%matrix_cap_livestemc_xfer_patch(i)       = 0._r8 
                      end if

                      if (pftcon%woody(patch%itype(i)) == 1._r8) then
                         this%deadstemc_patch(i)                       = 0.1_r8 * ratio
                         if(use_matrixcn)then
                            this%matrix_cap_deadstemc_patch(i)         = 0.1_r8 * ratio
                         end if
                      else
                         this%deadstemc_patch(i)                       = 0._r8 
                         if(use_matrixcn)then
                            this%matrix_cap_deadstemc_patch(i)         = 0._r8 
                         end if
                      end if
                      this%deadstemc_storage_patch(i)                  = 0._r8 
                      this%deadstemc_xfer_patch(i)                     = 0._r8 
                      if(use_matrixcn)then
                         this%matrix_cap_deadstemc_storage_patch(i)    = 0._r8 
                         this%matrix_cap_deadstemc_xfer_patch(i)       = 0._r8 
                      end if

                      this%livecrootc_patch(i)                         = 0._r8 
                      this%livecrootc_storage_patch(i)                 = 0._r8 
                      this%livecrootc_xfer_patch(i)                    = 0._r8 

                      this%deadcrootc_patch(i)                         = 0._r8 
                      this%deadcrootc_storage_patch(i)                 = 0._r8 
                      this%deadcrootc_xfer_patch(i)                    = 0._r8 

                      if(use_matrixcn)then
                         this%matrix_cap_livecrootc_patch(i)           = 0._r8 
                         this%matrix_cap_livecrootc_storage_patch(i)   = 0._r8 
                         this%matrix_cap_livecrootc_xfer_patch(i)      = 0._r8 
 
                         this%matrix_cap_deadcrootc_patch(i)           = 0._r8 
                         this%matrix_cap_deadcrootc_storage_patch(i)   = 0._r8 
                         this%matrix_cap_deadcrootc_xfer_patch(i)      = 0._r8 
                      end if

                      this%gresp_storage_patch(i)                      = 0._r8 
                      this%gresp_xfer_patch(i)                         = 0._r8 

                      this%cpool_patch(i)                              = 0._r8 
                      this%xsmrpool_patch(i)                           = 0._r8 
                      this%ctrunc_patch(i)                             = 0._r8 
                      this%dispvegc_patch(i)                           = 0._r8 
                      this%storvegc_patch(i)                           = 0._r8 
                      this%woodc_patch(i)                              = 0._r8
                      this%totc_patch(i)                               = 0._r8 

                      if ( use_crop )then
                         this%reproductivec_patch(i,:)         = 0._r8
                         this%reproductivec_storage_patch(i,:) = 0._r8
                         this%reproductivec_xfer_patch(i,:)    = 0._r8
                         if(use_matrixcn)then
                            this%reproc0_patch(i)                      = 0._r8 
                            this%reproc0_storage_patch(i)              = 0._r8 
                            this%reproc0_xfer_patch(i)                 = 0._r8 
                            this%matrix_cap_reproc_patch(i)            = 0._r8 
                            this%matrix_cap_reproc_storage_patch(i)    = 0._r8 
                            this%matrix_cap_reproc_xfer_patch(i)       = 0._r8 
                         end if
                         this%cropseedc_deficit_patch(i)               = 0._r8
                         this%xsmrpool_loss_patch(i)                   = 0._r8 
                      end if

                      ! calculate totvegc explicitly so that it is available for the isotope 
                      ! code on the first time step.

                      this%totvegc_patch(i) = &
                           this%leafc_patch(i)              + &
                           this%leafc_storage_patch(i)      + &
                           this%leafc_xfer_patch(i)         + &
                           this%frootc_patch(i)             + &
                           this%frootc_storage_patch(i)     + &
                           this%frootc_xfer_patch(i)        + &
                           this%livestemc_patch(i)          + &
                           this%livestemc_storage_patch(i)  + &
                           this%livestemc_xfer_patch(i)     + &
                           this%deadstemc_patch(i)          + &
                           this%deadstemc_storage_patch(i)  + &
                           this%deadstemc_xfer_patch(i)     + &
                           this%livecrootc_patch(i)         + &
                           this%livecrootc_storage_patch(i) + &
                           this%livecrootc_xfer_patch(i)    + &
                           this%deadcrootc_patch(i)         + &
                           this%deadcrootc_storage_patch(i) + &
                           this%deadcrootc_xfer_patch(i)    + &
                           this%gresp_storage_patch(i)      + &
                           this%gresp_xfer_patch(i)         + &
                           this%cpool_patch(i)

                      if ( use_crop )then
                         do k = 1, nrepr
                            this%totvegc_patch(i) =         &
                                 this%totvegc_patch(i)    + &
                                 this%reproductivec_patch(i,k)         + &
                                 this%reproductivec_storage_patch(i,k) + &
                                 this%reproductivec_xfer_patch(i,k)
                         end do
                      end if

                   endif
                end if
             end do
             if ( present(num_reseed_patch) ) then
                call shr_mpi_sum( num_reseed_patch, total_num_reseed_patch, mpicom )
                if ( masterproc ) write(iulog,*) 'Total num_reseed, over all tasks = ', total_num_reseed_patch
             end if
       end if

    end if

    !--------------------------------
    ! C13 patch carbon state variables 
    !--------------------------------

    if ( carbon_type == 'c13')  then

       call restartvar(ncid=ncid, flag=flag, varname='totvegc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totvegc_patch) 
       if (flag=='read' .and. .not. readvar) then
          ! BUG(kwo, 2023-10-19, ESCOMP/ctsm#2119) There is a bug that causes incorrect values for
          ! C isotopes if running from a case without C isotopes (an initial file) to a case with C
          ! isotopes (https://github.com/ESCOMP/ctsm/issues/2119). Here we check if the user
          ! is doing this and abort if they are. This particular check is covering the case
          ! when use_init_interp=.false.  There is a similar check (but for the purpose of working around
          ! a different bug) in initInterp.F90. This check below should be removed if bug #2119 is ever
          ! fully resolved (i.e., we decide we need to support going from a case without C isotopes (an 
          ! initial file) to a case with C isotopes), and replaced by the logic shown below for .e.g, 
          ! totvegc_col_13, where totvegc_col_13 is initialized with atmospheric c13 values.
          ! We arbitrarily check totvegc_13 (we could pick any c13 restart field).
          if (masterproc) then
             write(iulog,*) 'Cannot initialize from a run without c13 to a run with c13,'
             write(iulog,*) 'due to <https://github.com/ESCOMP/ctsm/issues/2119>.'
             write(iulog,*) 'Either use an input initial conditions file with c13 information,'
             write(iulog,*) 'or re-spinup from cold start.'
          end if
          missing_ciso = .true.
       end if

    end if

    if ( carbon_type == 'c14')  then
       call restartvar(ncid=ncid, flag=flag, varname='totvegc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totvegc_patch) 
       if (flag=='read' .and. .not. readvar) then
          ! BUG(kwo, 2023-10-19, ESCOMP/ctsm#2119) There is a bug that causes incorrect values for
          ! C isotopes if running from a case without C isotopes (an initial file) to a case with C
          ! isotopes (https://github.com/ESCOMP/ctsm/issues/2119). Here we check if the user
          ! is doing this and abort if they are. This particular check is covering the case
          ! when use_init_interp=.false.  There is a similar check (but for the purpose of working around
          ! a different bug) in initInterp.F90. This check below should be removed if bug #2119 is ever
          ! fully resolved (i.e., we decide we need to support going from a case without C isotopes (an 
          ! initial file) to a case with C isotopes), and replaced by the logic shown below for .e.g, 
          ! totvegc_col_14, where totvegc_col_l4 is initialized with atmospheric c14 values.
          ! We arbitrarily check totvegc_14 (we could pick any c14 restart field).
          if (masterproc) then
             write(iulog,*) 'Cannot interpolate from a run without c14 to a run with c14,'
             write(iulog,*) 'due to <https://github.com/ESCOMP/ctsm/issues/2119>.'
             write(iulog,*) 'Either use an input initial conditions file with c14 information,'
             write(iulog,*) 'or re-spinup from cold start.'
          end if
          missing_ciso = .true.
       endif
    end if

    if (missing_ciso) then
       call endrun(msg='Cannot initialize from a run without c13/c14 to a run with c13/c14', &
            additional_msg=errMsg(sourcefile, __LINE__))
    end if

    if ( carbon_type == 'c13')  then

       call restartvar(ncid=ncid, flag=flag, varname='totvegc_col_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totvegc_col)
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing cnveg_carbonstate_inst%totvegc with atmospheric c13 value'
          do i = bounds%begc,bounds%endc
             if (this%totvegc_col(i) /= spval .and. &
                 .not. isnan(this%totvegc_col(i)) ) then
                this%totvegc_col(i) = c12_cnveg_carbonstate_inst%totvegc_col(i) * c13ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='leafc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_patch)
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%leafc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%leafc_patch(i) = c12_cnveg_carbonstate_inst%leafc_patch(i) * c3_r2
             else
                this%leafc_patch(i) = c12_cnveg_carbonstate_inst%leafc_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%leafc_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%leafc_storage_patch(i) = c12_cnveg_carbonstate_inst%leafc_storage_patch(i) * c3_r2
             else
                this%leafc_storage_patch(i) = c12_cnveg_carbonstate_inst%leafc_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%leafc_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%leafc_xfer_patch(i) = c12_cnveg_carbonstate_inst%leafc_xfer_patch(i) * c3_r2
             else
                this%leafc_xfer_patch(i) = c12_cnveg_carbonstate_inst%leafc_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       if(use_matrixcn)then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_cap_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_leafc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_cap_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_leafc_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer_cap_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_leafc_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='leafc0_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc0_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='leafc0_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc0_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='leafc0_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc0_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_leaf_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_leaf_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_leafst_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_leafst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_leafst_to_leafxf_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_leafst_to_leafxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_leafxf_to_leaf_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_leafxf_to_leaf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_leaf_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_leaf_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_leafst_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_leafst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctrunover_leafxf_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_leafxf_acc_patch) 

       end if

       call restartvar(ncid=ncid, flag=flag, varname='frootc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%frootc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%frootc_patch(i) = c12_cnveg_carbonstate_inst%frootc_patch(i) * c3_r2
             else
                this%frootc_patch(i) = c12_cnveg_carbonstate_inst%frootc_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='frootc_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%frootc_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%frootc_storage_patch(i) = c12_cnveg_carbonstate_inst%frootc_storage_patch(i) * c3_r2
             else
                this%frootc_storage_patch(i) = c12_cnveg_carbonstate_inst%frootc_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%frootc_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%frootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%frootc_xfer_patch(i) * c3_r2
             else
                this%frootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%frootc_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       if(use_matrixcn)then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_cap_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_frootc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='frootc_storage_cap_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_frootc_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer_cap_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_frootc_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='frootc0_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc0_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='frootc0_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc0_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='frootc0_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc0_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_froot_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_froot_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_frootst_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_frootst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_frootst_to_frootxf_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_frootst_to_frootxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_frootxf_to_froot_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_frootxf_to_froot_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_froot_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_froot_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_frootst_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_frootst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_frootxf_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_frootxf_acc_patch) 
       end if


       call restartvar(ncid=ncid, flag=flag, varname='livestemc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livestemc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%livestemc_patch(i) = c12_cnveg_carbonstate_inst%livestemc_patch(i) * c3_r2
             else
                this%livestemc_patch(i) = c12_cnveg_carbonstate_inst%livestemc_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livestemc_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%livestemc_storage_patch(i) = c12_cnveg_carbonstate_inst%livestemc_storage_patch(i) * c3_r2
             else
                this%livestemc_storage_patch(i) = c12_cnveg_carbonstate_inst%livestemc_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livestemc_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%livestemc_xfer_patch(i) = c12_cnveg_carbonstate_inst%livestemc_xfer_patch(i) * c3_r2
             else
                this%livestemc_xfer_patch(i) = c12_cnveg_carbonstate_inst%livestemc_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       if(use_matrixcn)then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_cap_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_livestemc_patch) 
  
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage_cap_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_livestemc_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer_cap_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_livestemc_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livestemc0_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc0_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livestemc0_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc0_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livestemc0_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc0_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_livestem_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_livestem_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_livestemst_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_livestemst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_livestemst_to_livestemxf_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_livestemst_to_livestemxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_livestemxf_to_livestem_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_livestemxf_to_livestem_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_livestem_to_deadstem_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_livestem_to_deadstem_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_livestem_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_livestem_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_livestemst_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_livestemst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_livestemxf_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_livestemxf_acc_patch) 
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadstemc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%deadstemc_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_patch(i) * c3_r2
             else
                this%deadstemc_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadstemc_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%deadstemc_storage_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_storage_patch(i) * c3_r2
             else
                this%deadstemc_storage_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadstemc_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%deadstemc_xfer_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_xfer_patch(i) * c3_r2
             else
                this%deadstemc_xfer_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       if(use_matrixcn)then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_cap_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_deadstemc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage_cap_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_deadstemc_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer_cap_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_deadstemc_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc0_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc0_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc0_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc0_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc0_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc0_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_deadstem_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_deadstem_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_deadstemst_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_deadstemst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_deadstemst_to_deadstemxf_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_deadstemxf_to_deadstem_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_deadstemxf_to_deadstem_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_deadstem_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_deadstem_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_deadstemst_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_deadstemst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_deadstemxf_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_deadstemxf_acc_patch) 
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livecrootc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%livecrootc_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_patch(i) * c3_r2
             else
                this%livecrootc_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livecrootc_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%livecrootc_storage_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_storage_patch(i) * c3_r2
             else
                this%livecrootc_storage_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livecrootc_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%livecrootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_xfer_patch(i) * c3_r2
             else
                this%livecrootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       if(use_matrixcn)then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_cap_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_livecrootc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage_cap_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_livecrootc_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer_cap_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_livecrootc_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc0_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc0_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc0_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc0_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc0_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc0_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_livecroot_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_livecroot_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_livecrootst_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_livecrootst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_livecrootst_to_livecrootxf_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_livecrootxf_to_livecroot_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_livecrootxf_to_livecroot_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_livecroot_to_deadcroot_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_livecroot_to_deadcroot_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_livecroot_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_livecroot_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_livecrootst_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_livecrootst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_livecrootxf_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_livecrootxf_acc_patch) 
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadcrootc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%deadcrootc_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_patch(i) * c3_r2
             else
                this%deadcrootc_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadcrootc_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%deadcrootc_storage_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_storage_patch(i) * c3_r2
             else
                this%deadcrootc_storage_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadcrootc_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%deadcrootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_xfer_patch(i) * c3_r2
             else
                this%deadcrootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       if(use_matrixcn)then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_cap_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_deadcrootc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage_cap_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_deadcrootc_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer_cap_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_deadcrootc_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc0_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc0_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc0_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc0_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc0_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc0_xfer_patch) 
!
          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_deadcroot_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_deadcroot_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_deadcrootst_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_deadcrootst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_deadcrootxf_to_deadcroot_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_deadcroot_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_deadcroot_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_deadcrootst_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_deadcrootst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_deadcrootxf_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_deadcrootxf_acc_patch) 
       end if

       if(use_matrixcn .and. use_crop)then
          call restartvar(ncid=ncid, flag=flag,  varname='reproc0_13', xtype=ncd_double,  &
               dim1name='pft',    long_name='initial grain C13', units='gC13/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%reproc0_patch)
 
          call restartvar(ncid=ncid, flag=flag,  varname='reproc0_storage_13', xtype=ncd_double,  &
               dim1name='pft',    long_name='initial grain C13 storage', units='gC13/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%reproc0_storage_patch)

          call restartvar(ncid=ncid, flag=flag,  varname='reproc0_xfer_13', xtype=ncd_double,  &
               dim1name='pft',    long_name='initial grain C13 transfer', units='gC13/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%reproc0_xfer_patch)

          call restartvar(ncid=ncid, flag=flag,  varname='matrix_calloc_grain_acc_13', xtype=ncd_double,  &
               dim1name='pft',    long_name='C13 accumulated allocation to grain', units='gC13/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_grain_acc_patch)
 
          call restartvar(ncid=ncid, flag=flag,  varname='matrix_calloc_grainst_acc_13', xtype=ncd_double,  &
               dim1name='pft',    long_name='C13 accumulated allocation to grain storage', units='gC13/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_grainst_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_grainst_to_grainxf_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_grainst_to_grainxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_grainxf_to_grain_acc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_grainxf_to_grain_acc_patch)

          call restartvar(ncid=ncid, flag=flag,  varname='matrix_cturnover_grain_acc_13', xtype=ncd_double,  &
               dim1name='pft',    long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_grain_acc_patch)
 
          call restartvar(ncid=ncid, flag=flag,  varname='matrix_cturnover_grainst_acc_13', xtype=ncd_double,  &
               dim1name='pft',    long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_grainst_acc_patch)
 
          call restartvar(ncid=ncid, flag=flag,  varname='matrix_cturnover_grainxf_acc_13', xtype=ncd_double,  &
               dim1name='pft',    long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_grainxf_acc_patch)
       end if

       call restartvar(ncid=ncid, flag=flag, varname='gresp_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%gresp_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%gresp_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%gresp_storage_patch(i) = c12_cnveg_carbonstate_inst%gresp_storage_patch(i) * c3_r2
             else
                this%gresp_storage_patch(i) = c12_cnveg_carbonstate_inst%gresp_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%gresp_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%gresp_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%gresp_xfer_patch(i) = c12_cnveg_carbonstate_inst%gresp_xfer_patch(i) * c3_r2
             else
                this%gresp_xfer_patch(i) = c12_cnveg_carbonstate_inst%gresp_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='cpool_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%cpool_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%cpool with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%cpool_patch(i) = c12_cnveg_carbonstate_inst%cpool_patch(i) * c3_r2
             else
                this%cpool_patch(i) = c12_cnveg_carbonstate_inst%cpool_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%xsmrpool with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%xsmrpool_patch(i) = c12_cnveg_carbonstate_inst%xsmrpool_patch(i) * c3_r2
             else
                this%xsmrpool_patch(i) = c12_cnveg_carbonstate_inst%xsmrpool_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_loss_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_loss_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%xsmrpool_loss with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%xsmrpool_loss_patch(i) = c12_cnveg_carbonstate_inst%xsmrpool_loss_patch(i) * c3_r2
             else
                this%xsmrpool_loss_patch(i) = c12_cnveg_carbonstate_inst%xsmrpool_loss_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%ctrunc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%ctrunc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%ctrunc_patch(i) = c12_cnveg_carbonstate_inst%ctrunc_patch(i) * c3_r2
             else
                this%ctrunc_patch(i) = c12_cnveg_carbonstate_inst%ctrunc_patch(i) * c4_r2
             endif
          end do
       end if

    end if

    !--------------------------------
    ! C14 patch carbon state variables 
    !--------------------------------

    if ( carbon_type == 'c14')  then

       call restartvar(ncid=ncid, flag=flag, varname='totvegc_col_14', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totvegc_col)
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing cnveg_carbonstate_inst%totvegc with atmospheric c14 value'
          do i = bounds%begc,bounds%endc
             if (this%totvegc_col(i) /= spval .and. &
                 .not. isnan(this%totvegc_col(i)) ) then
                this%totvegc_col(i) = c12_cnveg_carbonstate_inst%totvegc_col(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='leafc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%leafc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%leafc_patch(i) /= spval .and. &
                  .not. isnan(this%leafc_patch(i)) ) then
                this%leafc_patch(i) = c12_cnveg_carbonstate_inst%leafc_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%leafc_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%leafc_storage_patch(i) /= spval .and. &
                  .not. isnan(this%leafc_storage_patch(i)) ) then
                this%leafc_storage_patch(i) = c12_cnveg_carbonstate_inst%leafc_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer_14', xtype=ncd_double,  &
            dim1name='pft',    long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%leafc_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%leafc_xfer_patch(i) /= spval .and. .not. isnan(this%leafc_xfer_patch(i)) ) then
                this%leafc_xfer_patch(i) = c12_cnveg_carbonstate_inst%leafc_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       if(use_matrixcn)then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_cap_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_leafc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_cap_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_leafc_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer_cap_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_leafc_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='leafc0_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc0_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='leafc0_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc0_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='leafc0_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc0_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_leaf_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_leaf_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_leafst_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_leafst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_leafst_to_leafxf_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_leafst_to_leafxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_leafxf_to_leaf_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_leafxf_to_leaf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_leaf_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_leaf_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_leafst_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_leafst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctrunover_leafxf_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_leafxf_acc_patch) 

       end if

       call restartvar(ncid=ncid, flag=flag, varname='frootc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%frootc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%frootc_patch(i) /= spval .and. &
                  .not. isnan(this%frootc_patch(i)) ) then
                this%frootc_patch(i) = c12_cnveg_carbonstate_inst%frootc_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='frootc_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%frootc_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%frootc_storage_patch(i) /= spval .and. &
                  .not. isnan(this%frootc_storage_patch(i)) ) then
                this%frootc_storage_patch(i) = c12_cnveg_carbonstate_inst%frootc_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%frootc_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%frootc_xfer_patch(i) /= spval .and. &
                  .not. isnan(this%frootc_xfer_patch(i)) ) then
                this%frootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%frootc_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       if(use_matrixcn)then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_cap_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_frootc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='frootc_storage_cap_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_frootc_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer_cap_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_frootc_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='frootc0_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc0_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='frootc0_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc0_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='frootc0_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc0_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_froot_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_froot_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_frootst_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_frootst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_frootst_to_frootxf_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_frootst_to_frootxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_frootxf_to_froot_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_frootxf_to_froot_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_froot_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_froot_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_frootst_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_frootst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_frootxf_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_frootxf_acc_patch) 
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livestemc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%livestemc_patch(i) /= spval .and. .not. isnan(this%livestemc_patch(i)) ) then
                this%livestemc_patch(i) = c12_cnveg_carbonstate_inst%livestemc_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livestemc_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%livestemc_storage_patch(i) /= spval .and. .not. isnan(this%livestemc_storage_patch(i)) ) then
                this%livestemc_storage_patch(i) = c12_cnveg_carbonstate_inst%livestemc_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livestemc_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%livestemc_xfer_patch(i) /= spval .and. .not. isnan(this%livestemc_xfer_patch(i)) ) then
                this%livestemc_xfer_patch(i) = c12_cnveg_carbonstate_inst%livestemc_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       if(use_matrixcn)then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_cap_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_livestemc_patch) 
  
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage_cap_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_livestemc_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer_cap_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_livestemc_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livestemc0_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc0_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livestemc0_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc0_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livestemc0_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc0_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_livestem_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_livestem_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_livestemst_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_livestemst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_livestemst_to_livestemxf_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_livestemst_to_livestemxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_livestemxf_to_livestem_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_livestemxf_to_livestem_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_livestem_to_deadstem_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_livestem_to_deadstem_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_livestem_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_livestem_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_livestemst_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_livestemst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_livestemxf_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_livestemxf_acc_patch) 
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadstemc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%deadstemc_patch(i) /= spval .and. .not. isnan(this%deadstemc_patch(i)) ) then
                this%deadstemc_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadstemc_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%deadstemc_storage_patch(i) /= spval .and. .not. isnan(this%deadstemc_storage_patch(i)) ) then
                this%deadstemc_storage_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadstemc_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%deadstemc_xfer_patch(i) /= spval .and. .not. isnan(this%deadstemc_xfer_patch(i)) ) then
                this%deadstemc_xfer_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       if(use_matrixcn)then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_cap_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_deadstemc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage_cap_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_deadstemc_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer_cap_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_deadstemc_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc0_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc0_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc0_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc0_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadstemc0_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc0_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_deadstem_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_deadstem_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_deadstemst_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_deadstemst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_deadstemst_to_deadstemxf_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_deadstemxf_to_deadstem_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_deadstemxf_to_deadstem_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_deadstem_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_deadstem_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_deadstemst_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_deadstemst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_deadstemxf_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_deadstemxf_acc_patch) 
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livecrootc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%livecrootc_patch(i) /= spval .and. .not. isnan(this%livecrootc_patch(i)) ) then
                this%livecrootc_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livecrootc_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%livecrootc_storage_patch(i) /= spval .and. .not. isnan(this%livecrootc_storage_patch(i)) ) then
                this%livecrootc_storage_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%livecrootc_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%livecrootc_xfer_patch(i) /= spval .and. .not. isnan(this%livecrootc_xfer_patch(i)) ) then
                this%livecrootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       if(use_matrixcn)then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_cap_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_livecrootc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage_cap_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_livecrootc_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer_cap_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_livecrootc_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc0_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc0_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc0_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc0_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='livecrootc0_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc0_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_livecroot_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_livecroot_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_livecrootst_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_livecrootst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_livecrootst_to_livecrootxf_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_livecrootxf_to_livecroot_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_livecrootxf_to_livecroot_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_livecroot_to_deadcroot_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_livecroot_to_deadcroot_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_livecroot_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_livecroot_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_livecrootst_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_livecrootst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_livecrootxf_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_livecrootxf_acc_patch) 
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadcrootc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%deadcrootc_patch(i) /= spval .and. .not. isnan(this%deadcrootc_patch(i)) ) then
                this%deadcrootc_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadcrootc_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%deadcrootc_storage_patch(i) /= spval .and. .not. isnan(this%deadcrootc_storage_patch(i)) ) then
                this%deadcrootc_storage_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%deadcrootc_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%deadcrootc_xfer_patch(i) /= spval .and. .not. isnan(this%deadcrootc_xfer_patch(i)) ) then
                this%deadcrootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       if(use_matrixcn)then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_cap_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_deadcrootc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage_cap_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_deadcrootc_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer_cap_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cap_deadcrootc_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc0_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc0_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc0_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc0_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc0_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc0_xfer_patch) 
!
          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_deadcroot_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_deadcroot_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_calloc_deadcrootst_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_deadcrootst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_deadcrootxf_to_deadcroot_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_deadcroot_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_deadcroot_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_deadcrootst_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_deadcrootst_acc_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='matrix_cturnover_deadcrootxf_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_deadcrootxf_acc_patch) 
       end if

       if(use_matrixcn .and. use_crop)then
          call restartvar(ncid=ncid, flag=flag,  varname='reproc0_14', xtype=ncd_double,  &
               dim1name='pft',    long_name='initial grain C14', units='gC14/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%reproc0_patch)
 
          call restartvar(ncid=ncid, flag=flag,  varname='reproc0_storage_14', xtype=ncd_double,  &
               dim1name='pft',    long_name='initial grain C14 storage', units='gC14/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%reproc0_storage_patch)

          call restartvar(ncid=ncid, flag=flag,  varname='reproc0_xfer_14', xtype=ncd_double,  &
               dim1name='pft',    long_name='initial grain C14 transfer', units='gC14/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%reproc0_xfer_patch)

          call restartvar(ncid=ncid, flag=flag,  varname='matrix_calloc_grain_acc_14', xtype=ncd_double,  &
               dim1name='pft',    long_name='C14 accumulated allocation to grain', units='gC14/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_grain_acc_patch)
 
          call restartvar(ncid=ncid, flag=flag,  varname='matrix_calloc_grainst_acc_14', xtype=ncd_double,  &
               dim1name='pft',    long_name='C14 accumulated allocation to grain storage', units='gC14/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_calloc_grainst_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_grainst_to_grainxf_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_grainst_to_grainxf_acc_patch)

          call restartvar(ncid=ncid, flag=flag, varname='matrix_ctransfer_grainxf_to_grain_acc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_ctransfer_grainxf_to_grain_acc_patch)

          call restartvar(ncid=ncid, flag=flag,  varname='matrix_cturnover_grain_acc_14', xtype=ncd_double,  &
               dim1name='pft',    long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_grain_acc_patch)
 
          call restartvar(ncid=ncid, flag=flag,  varname='matrix_cturnover_grainst_acc_14', xtype=ncd_double,  &
               dim1name='pft',    long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_grainst_acc_patch)
 
          call restartvar(ncid=ncid, flag=flag,  varname='matrix_cturnover_grainxf_acc_14', xtype=ncd_double,  &
               dim1name='pft',    long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%matrix_cturnover_grainxf_acc_patch)
       end if

       call restartvar(ncid=ncid, flag=flag, varname='gresp_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%gresp_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%gresp_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%gresp_storage_patch(i) /= spval .and. .not. isnan(this%gresp_storage_patch(i)) ) then
                this%gresp_storage_patch(i) = c12_cnveg_carbonstate_inst%gresp_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%gresp_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%gresp_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%gresp_xfer_patch(i) /= spval .and. .not. isnan(this%gresp_xfer_patch(i)) ) then
                this%gresp_xfer_patch(i) = c12_cnveg_carbonstate_inst%gresp_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='cpool_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%cpool_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%cpool_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%cpool_patch(i) /= spval .and. .not. isnan(this%cpool_patch(i)) ) then
                this%cpool_patch(i) = c12_cnveg_carbonstate_inst%cpool_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%xsmrpool_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%xsmrpool_patch(i) /= spval .and. .not. isnan(this%xsmrpool_patch(i)) ) then
                this%xsmrpool_patch(i) = c12_cnveg_carbonstate_inst%xsmrpool_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_loss_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_loss_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%xsmrpool_loss_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%xsmrpool_loss_patch(i) /= spval .and. .not. isnan(this%xsmrpool_loss_patch(i)) ) then
                this%xsmrpool_loss_patch(i) = c12_cnveg_carbonstate_inst%xsmrpool_loss_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%ctrunc_patch) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%ctrunc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%ctrunc_patch(i) /= spval .and. .not. isnan(this%ctrunc_patch(i)) ) then
                this%ctrunc_patch(i) = c12_cnveg_carbonstate_inst%ctrunc_patch(i) * c14ratio
             endif
          end do
       end if

    end if

    !--------------------------------
    ! patch prognostic crop variables
    !--------------------------------

    if (use_crop) then
       if (carbon_type == 'c12') then
          do k = 1, nrepr
             data1dptr => this%reproductivec_patch(:,k)
             ! e.g., reproductivec
             varname = get_repr_rest_fname(k)//'c'
             call restartvar(ncid=ncid, flag=flag,  varname=varname, &
                  xtype=ncd_double,  &
                  dim1name='pft', &
                  long_name=get_repr_longname(k)//' C', &
                  units='gC/m2', &
                  interpinic_flag='interp', readvar=readvar, data=data1dptr)
          end do

          do k = 1, nrepr
             data1dptr => this%reproductivec_storage_patch(:,k)
             ! e.g., reproductivec_storage
             varname = get_repr_rest_fname(k)//'c_storage'
             call restartvar(ncid=ncid, flag=flag,  varname=varname, &
                  xtype=ncd_double,  &
                  dim1name='pft', &
                  long_name=get_repr_longname(k)//' C storage', &
                  units='gC/m2', &
                  interpinic_flag='interp', readvar=readvar, data=data1dptr)
          end do

          do k = 1, nrepr
             data1dptr => this%reproductivec_xfer_patch(:,k)
             ! e.g., reproductivec_xfer
             varname = get_repr_rest_fname(k)//'c_xfer'
             call restartvar(ncid=ncid, flag=flag,  varname=varname, &
                  xtype=ncd_double,  &
                  dim1name='pft', &
                  long_name=get_repr_longname(k)//' C transfer', &
                  units='gC/m2', &
                  interpinic_flag='interp', readvar=readvar, data=data1dptr)
          end do

          call restartvar(ncid=ncid, flag=flag, varname='cropseedc_deficit', xtype=ncd_double,  &
               dim1name='pft', long_name='pool for seeding new crop growth', units='gC/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%cropseedc_deficit_patch)
       end if

       if (carbon_type == 'c13') then
          do k = 1, nrepr
             data1dptr => this%reproductivec_patch(:,k)
             ! e.g., reprocuctive_13
             varname = get_repr_rest_fname(k)//'c_13'
             call restartvar(ncid=ncid, flag=flag, varname=varname, &
                  xtype=ncd_double,  &
                  dim1name='pft', &
                  long_name='c13 '//get_repr_longname(k)//' C', &
                  units='gC13/m2', &
                  interpinic_flag='interp', readvar=readvar, data=data1dptr)
             if (flag=='read' .and. .not. readvar) then
                call set_missing_from_template( &
                     my_var = data1dptr, &
                     template_var = c12_cnveg_carbonstate_inst%reproductivec_patch(:,k), &
                     multiplier = c3_r2)
             end if
          end do

          do k = 1, nrepr
             data1dptr => this%reproductivec_storage_patch(:,k)
             ! e.g., reproductivec_13_storage
             varname = get_repr_rest_fname(k)//'c_13_storage'
             call restartvar(ncid=ncid, flag=flag, varname=varname, &
                  xtype=ncd_double,  &
                  dim1name='pft', &
                  long_name='c13 '//get_repr_longname(k)//' C storage', &
                  units='gC13/m2', &
                  interpinic_flag='interp', readvar=readvar, data=data1dptr)
             if (flag=='read' .and. .not. readvar) then
                call set_missing_from_template( &
                     my_var = data1dptr, &
                     template_var = c12_cnveg_carbonstate_inst%reproductivec_storage_patch(:,k), &
                     multiplier = c3_r2)
             end if
          end do

          do k = 1, nrepr
             data1dptr => this%reproductivec_xfer_patch(:,k)
             ! e.g., reproductivec_13_xfer
             varname = get_repr_rest_fname(k)//'c_13_xfer'
             call restartvar(ncid=ncid, flag=flag, varname=varname, &
                  xtype=ncd_double,  &
                  dim1name='pft', &
                  long_name='c13 '//get_repr_longname(k)//' C transfer', &
                  units='gC13/m2', &
                  interpinic_flag='interp', readvar=readvar, data=data1dptr)
             if (flag=='read' .and. .not. readvar) then
                call set_missing_from_template( &
                     my_var = data1dptr, &
                     template_var = c12_cnveg_carbonstate_inst%reproductivec_xfer_patch(:,k), &
                     multiplier = c3_r2)
             end if
          end do

          call restartvar(ncid=ncid, flag=flag, varname='cropseedc_13_deficit', xtype=ncd_double,  &
               dim1name='pft', long_name='pool for seeding new crop growth', units='gC13/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%cropseedc_deficit_patch)
          if (flag=='read' .and. .not. readvar) then
             call set_missing_from_template( &
                  my_var = this%cropseedc_deficit_patch, &
                  template_var = c12_cnveg_carbonstate_inst%cropseedc_deficit_patch, &
                  multiplier = c3_r2)
          end if
       end if

       if ( carbon_type == 'c14' ) then

          do k = 1, nrepr
             data1dptr => this%reproductivec_patch(:,k)
             ! e.g., reproductivec_14
             varname = get_repr_rest_fname(k)//'c_14'
             call restartvar(ncid=ncid, flag=flag, varname=varname, &
                  xtype=ncd_double,  &
                  dim1name='pft', &
                  long_name='c14 '//get_repr_longname(k)//' C', &
                  units='gC14/m2', &
                  interpinic_flag='interp', readvar=readvar, data=data1dptr)
             if (flag=='read' .and. .not. readvar) then
                call set_missing_from_template( &
                     my_var = data1dptr, &
                     template_var = c12_cnveg_carbonstate_inst%reproductivec_patch(:,k), &
                     multiplier = c3_r2)
             end if
          end do

          do k = 1, nrepr
             data1dptr => this%reproductivec_storage_patch(:,k)
             ! e.g., reproductivec_14_storage
             varname = get_repr_rest_fname(k)//'c_14_storage'
             call restartvar(ncid=ncid, flag=flag, varname=varname, &
                  xtype=ncd_double,  &
                  dim1name='pft', &
                  long_name='c14 '//get_repr_longname(k)//' C storage', &
                  units='gC14/m2', &
                  interpinic_flag='interp', readvar=readvar, data=data1dptr)
             if (flag=='read' .and. .not. readvar) then
                call set_missing_from_template( &
                     my_var = data1dptr, &
                     template_var = c12_cnveg_carbonstate_inst%reproductivec_storage_patch(:,k), &
                     multiplier = c3_r2)
             end if
          end do

          do k = 1, nrepr
             data1dptr => this%reproductivec_xfer_patch(:,k)
             ! e.g., reproductivec_14_xfer
             varname = get_repr_rest_fname(k)//'c_14_xfer'
             call restartvar(ncid=ncid, flag=flag, varname=varname, &
                  xtype=ncd_double,  &
                  dim1name='pft', &
                  long_name='c14 '//get_repr_longname(k)//' C transfer', &
                  units='gC14/m2', &
                  interpinic_flag='interp', readvar=readvar, data=data1dptr)
             if (flag=='read' .and. .not. readvar) then
                call set_missing_from_template( &
                     my_var = data1dptr, &
                     template_var = c12_cnveg_carbonstate_inst%reproductivec_xfer_patch(:,k), &
                     multiplier = c3_r2)
             end if
          end do

          call restartvar(ncid=ncid, flag=flag, varname='cropseedc_14_deficit', xtype=ncd_double,  &
               dim1name='pft', long_name='pool for seeding new crop growth', units='gC14/m2', &
               interpinic_flag='interp', readvar=readvar, data=this%cropseedc_deficit_patch)
          if (flag=='read' .and. .not. readvar) then
             if ( masterproc ) write(iulog,*) 'initializing this%cropseedc_deficit_patch with atmospheric c14 value'
             call set_missing_from_template( &
                  my_var = this%cropseedc_deficit_patch, &
                  template_var = c12_cnveg_carbonstate_inst%cropseedc_deficit_patch, &
                  multiplier = c14ratio)
          end if
       end if
    end if

    !--------------------------------
    ! gridcell carbon state variables
    !--------------------------------

    if (carbon_type == 'c12') then
       ! BACKWARDS_COMPATIBILITY(wjs, 2017-01-12) Naming this with a _g suffix in order
       ! to distinguish it from the old column-level seedc restart variable
       call restartvar(ncid=ncid, flag=flag, varname='seedc_g', xtype=ncd_double,  &
            dim1name='gridcell', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%seedc_grc) 
    end if

    !--------------------------------
    ! C13 gridcell carbon state variables
    !--------------------------------

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='seedc_13_g', xtype=ncd_double,  &
            dim1name='gridcell', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%seedc_grc) 
       if (flag=='read' .and. .not. readvar) then
          call set_missing_from_template( &
               my_var = this%seedc_grc, &
               template_var = c12_cnveg_carbonstate_inst%seedc_grc, &
               multiplier = c3_r2)
       end if
    end if

    !--------------------------------
    ! C14 column carbon state variables
    !--------------------------------

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='seedc_14_g', xtype=ncd_double,  &
            dim1name='gridcell', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%seedc_grc) 
       if (flag=='read' .and. .not. readvar) then
          if ( masterproc ) write(iulog,*) 'initializing this%seedc_grc with atmospheric c14 value'
          call set_missing_from_template( &
               my_var = this%seedc_grc, &
               template_var = c12_cnveg_carbonstate_inst%seedc_grc, &
               multiplier = c14ratio)
       end if
    end if

    ! Output spinup factor for deadwood (dead stem and dead course root)
    if ( present(spinup_factor4deadwood) ) spinup_factor4deadwood = spinup_factor_AD

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine SetValues ( this, &
       num_patch, filter_patch, value_patch, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set carbon state variables
    !
    ! !ARGUMENTS:
    class (cnveg_carbonstate_type) :: this
    integer , intent(in) :: num_patch
    integer , intent(in) :: filter_patch(:)
    real(r8), intent(in) :: value_patch
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i,j,k,l     ! loop index
    !------------------------------------------------------------------------

    do fi = 1,num_patch
       i  = filter_patch(fi)
       this%leafc_patch(i)              = value_patch
       this%leafc_storage_patch(i)      = value_patch
       this%leafc_xfer_patch(i)         = value_patch
       this%leafc_storage_xfer_acc_patch(i) = value_patch
       this%storage_cdemand_patch(i)        = value_patch        
       this%frootc_patch(i)             = value_patch
       this%frootc_storage_patch(i)     = value_patch
       this%frootc_xfer_patch(i)        = value_patch
       this%livestemc_patch(i)          = value_patch
       this%livestemc_storage_patch(i)  = value_patch
       this%livestemc_xfer_patch(i)     = value_patch
       this%deadstemc_patch(i)          = value_patch
       this%deadstemc_storage_patch(i)  = value_patch
       this%deadstemc_xfer_patch(i)     = value_patch
       this%livecrootc_patch(i)         = value_patch
       this%livecrootc_storage_patch(i) = value_patch
       this%livecrootc_xfer_patch(i)    = value_patch
       this%deadcrootc_patch(i)         = value_patch
       this%deadcrootc_storage_patch(i) = value_patch
       this%deadcrootc_xfer_patch(i)    = value_patch
       if(use_matrixcn)then
          this%matrix_cap_leafc_patch(i)              = value_patch
          this%matrix_cap_leafc_storage_patch(i)      = value_patch
          this%matrix_cap_leafc_xfer_patch(i)         = value_patch
          this%matrix_cap_frootc_patch(i)             = value_patch
          this%matrix_cap_frootc_storage_patch(i)     = value_patch
          this%matrix_cap_frootc_xfer_patch(i)        = value_patch
          this%matrix_cap_livestemc_patch(i)          = value_patch
          this%matrix_cap_livestemc_storage_patch(i)  = value_patch
          this%matrix_cap_livestemc_xfer_patch(i)     = value_patch
          this%matrix_cap_deadstemc_patch(i)          = value_patch
          this%matrix_cap_deadstemc_storage_patch(i)  = value_patch
          this%matrix_cap_deadstemc_xfer_patch(i)     = value_patch
          this%matrix_cap_livecrootc_patch(i)         = value_patch
          this%matrix_cap_livecrootc_storage_patch(i) = value_patch
          this%matrix_cap_livecrootc_xfer_patch(i)    = value_patch
          this%matrix_cap_deadcrootc_patch(i)         = value_patch
          this%matrix_cap_deadcrootc_storage_patch(i) = value_patch
          this%matrix_cap_deadcrootc_xfer_patch(i)    = value_patch

          this%leafc0_patch(i)              = value_patch
          this%leafc0_storage_patch(i)      = value_patch
          this%leafc0_xfer_patch(i)         = value_patch   
          this%frootc0_patch(i)             = value_patch
          this%frootc0_storage_patch(i)     = value_patch
          this%frootc0_xfer_patch(i)        = value_patch
          this%livestemc0_patch(i)          = value_patch
          this%livestemc0_storage_patch(i)  = value_patch
          this%livestemc0_xfer_patch(i)     = value_patch
          this%deadstemc0_patch(i)          = value_patch
          this%deadstemc0_storage_patch(i)  = value_patch
          this%deadstemc0_xfer_patch(i)     = value_patch
          this%livecrootc0_patch(i)         = value_patch
          this%livecrootc0_storage_patch(i) = value_patch
          this%livecrootc0_xfer_patch(i)    = value_patch
          this%deadcrootc0_patch(i)         = value_patch
          this%deadcrootc0_storage_patch(i) = value_patch
          this%deadcrootc0_xfer_patch(i)    = value_patch
          this%reproc0_patch(i)             = value_patch
          this%reproc0_storage_patch(i)     = value_patch
          this%reproc0_xfer_patch(i)        = value_patch
!!!!matrix
          this%matrix_calloc_leaf_acc_patch(i)        =  value_patch
          this%matrix_calloc_leafst_acc_patch(i)      =  value_patch
          this%matrix_calloc_froot_acc_patch(i)       =  value_patch
          this%matrix_calloc_frootst_acc_patch(i)     =  value_patch
          this%matrix_calloc_livestem_acc_patch(i)    =  value_patch
          this%matrix_calloc_livestemst_acc_patch(i)  =  value_patch
          this%matrix_calloc_deadstem_acc_patch(i)    =  value_patch
          this%matrix_calloc_deadstemst_acc_patch(i)  =  value_patch
          this%matrix_calloc_livecroot_acc_patch(i)   =  value_patch
          this%matrix_calloc_livecrootst_acc_patch(i) =  value_patch
          this%matrix_calloc_deadcroot_acc_patch(i)   =  value_patch
          this%matrix_calloc_deadcrootst_acc_patch(i) =  value_patch

          this%matrix_ctransfer_leafst_to_leafxf_acc_patch           (i) = value_patch
          this%matrix_ctransfer_leafxf_to_leaf_acc_patch             (i) = value_patch
          this%matrix_ctransfer_frootst_to_frootxf_acc_patch         (i) = value_patch
          this%matrix_ctransfer_frootxf_to_froot_acc_patch           (i) = value_patch
          this%matrix_ctransfer_livestemst_to_livestemxf_acc_patch   (i) = value_patch
          this%matrix_ctransfer_livestemxf_to_livestem_acc_patch     (i) = value_patch
          this%matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch   (i) = value_patch
          this%matrix_ctransfer_deadstemxf_to_deadstem_acc_patch     (i) = value_patch
          this%matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch (i) = value_patch
          this%matrix_ctransfer_livecrootxf_to_livecroot_acc_patch   (i) = value_patch
          this%matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch (i) = value_patch
          this%matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch   (i) = value_patch
          this%matrix_ctransfer_livestem_to_deadstem_acc_patch       (i) = value_patch
          this%matrix_ctransfer_livecroot_to_deadcroot_acc_patch     (i) = value_patch

          this%matrix_cturnover_leaf_acc_patch(i)        = value_patch
          this%matrix_cturnover_leafst_acc_patch(i)      = value_patch
          this%matrix_cturnover_leafxf_acc_patch(i)      = value_patch   
          this%matrix_cturnover_froot_acc_patch(i)       = value_patch
          this%matrix_cturnover_frootst_acc_patch(i)     = value_patch
          this%matrix_cturnover_frootxf_acc_patch(i)     = value_patch   
          this%matrix_cturnover_livestem_acc_patch(i)    = value_patch
          this%matrix_cturnover_livestemst_acc_patch(i)  = value_patch
          this%matrix_cturnover_livestemxf_acc_patch(i)  = value_patch   
          this%matrix_cturnover_deadstem_acc_patch(i)    = value_patch
          this%matrix_cturnover_deadstemst_acc_patch(i)  = value_patch
          this%matrix_cturnover_deadstemxf_acc_patch(i)  = value_patch   
          this%matrix_cturnover_livecroot_acc_patch(i)   = value_patch
          this%matrix_cturnover_livecrootst_acc_patch(i) = value_patch
          this%matrix_cturnover_livecrootxf_acc_patch(i) = value_patch   
          this%matrix_cturnover_deadcroot_acc_patch(i)   = value_patch
          this%matrix_cturnover_deadcrootst_acc_patch(i) = value_patch
          this%matrix_cturnover_deadcrootxf_acc_patch(i) = value_patch   
       end if
       this%gresp_storage_patch(i)      = value_patch
       this%gresp_xfer_patch(i)         = value_patch
       this%cpool_patch(i)              = value_patch
       this%xsmrpool_patch(i)           = value_patch
       this%ctrunc_patch(i)             = value_patch
       this%dispvegc_patch(i)           = value_patch
       this%storvegc_patch(i)           = value_patch
       this%woodc_patch(i)              = value_patch
       this%totvegc_patch(i)            = value_patch
       this%totc_patch(i)               = value_patch
       if ( use_crop ) then
          if(use_matrixcn)then
             this%matrix_calloc_grain_acc_patch(i)                  = value_patch
             this%matrix_calloc_grainst_acc_patch(i)                = value_patch
             this%matrix_ctransfer_grainst_to_grainxf_acc_patch (i) = value_patch
             this%matrix_ctransfer_grainxf_to_grain_acc_patch   (i) = value_patch
             this%matrix_cturnover_grain_acc_patch(i)               = value_patch
             this%matrix_cturnover_grainst_acc_patch(i)             = value_patch
             this%matrix_cturnover_grainxf_acc_patch(i)             = value_patch
          end if
          this%cropseedc_deficit_patch(i)  = value_patch
          this%xsmrpool_loss_patch(i)   = value_patch
       end if
    end do

    if (use_crop) then
       do k = 1, nrepr
          do fi = 1,num_patch
             i  = filter_patch(fi)
             this%reproductivec_patch(i,k)          = value_patch
             this%reproductivec_storage_patch(i,k)  = value_patch
             this%reproductivec_xfer_patch(i,k)     = value_patch
          end do
       end do
       if(use_matrixcn)then
          do fi = 1,num_column
             i  = filter_column(fi)
             this%matrix_cap_reproc_patch(i)           = value_patch
             this%matrix_cap_reproc_storage_patch(i)   = value_patch
             this%matrix_cap_reproc_xfer_patch(i)      = value_patch
           end do
       end if
    end if

    do fi = 1,num_column
       i  = filter_column(fi)
       this%rootc_col(i)                = value_column
       this%leafc_col(i)                = value_column
       this%deadstemc_col(i)            = value_column
       this%fuelc_col(i)                = value_column
       this%fuelc_crop_col(i)           = value_column
       this%totvegc_col(i)              = value_column
       this%totc_p2c_col(i)             = value_column
    end do

  end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine ZeroDwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(cnveg_carbonstate_type) :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: p          ! indices
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       this%dispvegc_patch(p)   = 0._r8
       this%storvegc_patch(p)   = 0._r8
       this%totc_patch(p)       = 0._r8
    end do

  end subroutine ZeroDwt

  !-----------------------------------------------------------------------
  subroutine Summary_carbonstate(this, bounds, num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp)

    !
    ! !USES:
    use subgridAveMod, only : p2c
    use clm_time_manager , only : get_nstep

    !
    ! !DESCRIPTION:
    ! Perform patch and column-level carbon summary calculations
    !
    ! !ARGUMENTS:
    class(cnveg_carbonstate_type)  :: this
    type(bounds_type) , intent(in) :: bounds          
    integer           , intent(in) :: num_bgc_soilc       ! number of bgc soil columns in filter
    integer           , intent(in) :: filter_bgc_soilc(:) ! filter for bgc soil columns
    integer           , intent(in) :: num_bgc_vegp       ! number of soil patches in filter
    integer           , intent(in) :: filter_bgc_vegp(:) ! filter for soil patches

    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,k,l       ! indices
    integer  :: fp,fc           ! lake filter indices
    !-----------------------------------------------------------------------

    ! calculate patch -level summary of carbon state
    do fp = 1,num_bgc_vegp
       p = filter_bgc_vegp(fp)

       ! displayed vegetation carbon, excluding storage and cpool (DISPVEGC)
       this%dispvegc_patch(p) =        &
            this%leafc_patch(p)      + &
            this%frootc_patch(p)     + &
            this%livestemc_patch(p)  + &
            this%deadstemc_patch(p)  + &
            this%livecrootc_patch(p) + &
            this%deadcrootc_patch(p)

       ! stored vegetation carbon, excluding cpool (STORVEGC)
       this%storvegc_patch(p) =                &
            this%cpool_patch(p)              + &
            this%leafc_storage_patch(p)      + &
            this%frootc_storage_patch(p)     + &
            this%livestemc_storage_patch(p)  + &
            this%deadstemc_storage_patch(p)  + &
            this%livecrootc_storage_patch(p) + &
            this%deadcrootc_storage_patch(p) + &
            this%leafc_xfer_patch(p)         + &
            this%frootc_xfer_patch(p)        + &
            this%livestemc_xfer_patch(p)     + &
            this%deadstemc_xfer_patch(p)     + &
            this%livecrootc_xfer_patch(p)    + &
            this%deadcrootc_xfer_patch(p)    + &
            this%gresp_storage_patch(p)      + &
            this%gresp_xfer_patch(p)

       if ( use_crop .and. patch%itype(p) >= npcropmin )then
          do k = 1, nrepr
             this%storvegc_patch(p) =            &
                  this%storvegc_patch(p)       + &
                  this%reproductivec_storage_patch(p,k) + &
                  this%reproductivec_xfer_patch(p,k)

             this%dispvegc_patch(p) =            &
                  this%dispvegc_patch(p)       + &
                  this%reproductivec_patch(p,k)
          end do
       end if

       ! total vegetation carbon, excluding cpool (TOTVEGC)
       this%totvegc_patch(p) = &
            this%dispvegc_patch(p) + &
            this%storvegc_patch(p)

       ! total patch-level carbon, including xsmrpool, ctrunc
       this%totc_patch(p) = &
            this%totvegc_patch(p) + &
            this%xsmrpool_patch(p) + &
            this%ctrunc_patch(p)

       if (use_crop) then 
          this%totc_patch(p) = this%totc_patch(p) + this%cropseedc_deficit_patch(p) + &
               this%xsmrpool_loss_patch(p)
       end if

       ! (WOODC) - wood C
       this%woodc_patch(p) = &
            this%deadstemc_patch(p)    + &
            this%livestemc_patch(p)    + &
            this%deadcrootc_patch(p)   + &
            this%livecrootc_patch(p)

    end do

    ! --------------------------------------------
    ! column level summary
    ! --------------------------------------------
    if(num_bgc_vegp>0)then
       call p2c(bounds, num_bgc_soilc, filter_bgc_soilc, &
            this%totvegc_patch(bounds%begp:bounds%endp), &
            this%totvegc_col(bounds%begc:bounds%endc))
       
       call p2c(bounds, num_bgc_soilc, filter_bgc_soilc, &
            this%totc_patch(bounds%begp:bounds%endp), &
            this%totc_p2c_col(bounds%begc:bounds%endc))
    end if
       
  end subroutine Summary_carbonstate

  !-----------------------------------------------------------------------
  subroutine DynamicPatchAdjustments(this, bounds, &
       num_soilp_with_inactive, filter_soilp_with_inactive, &
       patch_state_updater, &
       leafc_seed, deadstemc_seed, &
       conv_cflux, wood_product_cflux, crop_product_cflux, &
       dwt_frootc_to_litter, &
       dwt_livecrootc_to_litter, &
       dwt_deadcrootc_to_litter, &
       dwt_leafc_seed, &
       dwt_deadstemc_seed)
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cnveg_carbonstate_type)   , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_soilp_with_inactive ! number of points in filter
    integer                         , intent(in)    :: filter_soilp_with_inactive(:) ! soil patch filter that includes inactive points
    type(patch_state_updater_type)  , intent(in)    :: patch_state_updater
    real(r8)                        , intent(in)    :: leafc_seed  ! seed amount for leaf C
    real(r8)                        , intent(in)    :: deadstemc_seed ! seed amount for deadstem C
    real(r8)                        , intent(inout) :: conv_cflux( bounds%begp: )  ! patch-level conversion C flux to atm (expressed per unit GRIDCELL area)
    real(r8)                        , intent(inout) :: wood_product_cflux( bounds%begp: ) ! patch-level product C flux (expressed per unit GRIDCELL area)
    real(r8)                        , intent(inout) :: crop_product_cflux( bounds%begp: ) ! patch-level crop product C flux (expressed per unit GRIDCELL area)
    real(r8)                        , intent(inout) :: dwt_frootc_to_litter( bounds%begp: ) ! patch-level fine root C to litter (expressed per unit COLUMN area)
    real(r8)                        , intent(inout) :: dwt_livecrootc_to_litter( bounds%begp: ) ! patch-level live coarse root C to litter (expressed per unit COLUMN area)
    real(r8)                        , intent(inout) :: dwt_deadcrootc_to_litter( bounds%begp: ) ! patch-level live coarse root C to litter (expressed per unit COLUMN area)
    real(r8)                        , intent(inout) :: dwt_leafc_seed( bounds%begp: ) ! patch-level mass gain due to seeding of new area: leaf C (expressed per unit GRIDCELL area)
    real(r8)                        , intent(inout) :: dwt_deadstemc_seed( bounds%begp: ) ! patch-level mass gain due to seeding of new area: deadstem C (expressed per unit GRIDCELL area)
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: k

    logical  :: old_weight_was_zero(bounds%begp:bounds%endp)
    logical  :: patch_grew(bounds%begp:bounds%endp)

    ! The following are only set for growing patches:
    real(r8) :: seed_leafc_patch(bounds%begp:bounds%endp)
    real(r8) :: seed_leafc_storage_patch(bounds%begp:bounds%endp)
    real(r8) :: seed_leafc_xfer_patch(bounds%begp:bounds%endp)
    real(r8) :: seed_deadstemc_patch(bounds%begp:bounds%endp)

    character(len=*), parameter :: subname = 'DynamicPatchAdjustments'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    SHR_ASSERT_ALL_FL((ubound(conv_cflux) == (/endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(wood_product_cflux) == (/endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(crop_product_cflux) == (/endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(dwt_frootc_to_litter) == (/endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(dwt_livecrootc_to_litter) == (/endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(dwt_deadcrootc_to_litter) == (/endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(dwt_leafc_seed) == (/endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(dwt_deadstemc_seed) == (/endp/)), sourcefile, __LINE__)

    old_weight_was_zero = patch_state_updater%old_weight_was_zero(bounds)
    patch_grew = patch_state_updater%patch_grew(bounds)

    call ComputeSeedAmounts(bounds, &
         num_soilp_with_inactive, filter_soilp_with_inactive, &
         species = this%species, &
         leafc_seed = leafc_seed, &
         deadstemc_seed = deadstemc_seed, &
         leaf_patch = this%leafc_patch(begp:endp), &
         leaf_storage_patch = this%leafc_storage_patch(begp:endp), &
         leaf_xfer_patch = this%leafc_xfer_patch(begp:endp), &

         ! Calculations only needed for patches that grew:
         compute_here_patch = patch_grew(begp:endp), &

         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero(begp:endp), &

         seed_leaf_patch = seed_leafc_patch(begp:endp), &
         seed_leaf_storage_patch = seed_leafc_storage_patch(begp:endp), &
         seed_leaf_xfer_patch = seed_leafc_xfer_patch(begp:endp), &
         seed_deadstem_patch = seed_deadstemc_patch(begp:endp))

    call update_patch_state( &
         var = this%leafc_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp), &
         seed = seed_leafc_patch(begp:endp), &
         seed_addition = dwt_leafc_seed(begp:endp))

    call update_patch_state( &
         var = this%leafc_storage_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp), &
         seed = seed_leafc_storage_patch(begp:endp), &
         seed_addition = dwt_leafc_seed(begp:endp))

    call update_patch_state( &
         var = this%leafc_xfer_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp), &
         seed = seed_leafc_xfer_patch(begp:endp), &
         seed_addition = dwt_leafc_seed(begp:endp))

    call update_patch_state( &
         var = this%frootc_patch(begp:endp), &
         flux_out_col_area = dwt_frootc_to_litter(begp:endp))

    call update_patch_state( &
         var = this%frootc_storage_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%frootc_xfer_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%livestemc_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%livestemc_storage_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%livestemc_xfer_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call patch_state_updater%update_patch_state_partition_flux_by_type(bounds, &
         num_soilp_with_inactive, filter_soilp_with_inactive, &
         flux1_fraction_by_pft_type = pftcon%pconv, &
         var = this%deadstemc_patch(begp:endp), &
         flux1_out = conv_cflux(begp:endp), &
         flux2_out = wood_product_cflux(begp:endp), &
         seed = seed_deadstemc_patch(begp:endp), &
         seed_addition = dwt_deadstemc_seed(begp:endp))

    call update_patch_state( &
         var = this%deadstemc_storage_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%deadstemc_xfer_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%livecrootc_patch(begp:endp), &
         flux_out_col_area = dwt_livecrootc_to_litter(begp:endp))

    call update_patch_state( &
         var = this%livecrootc_storage_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%livecrootc_xfer_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%deadcrootc_patch(begp:endp), &
         flux_out_col_area = dwt_deadcrootc_to_litter(begp:endp))

    call update_patch_state( &
         var = this%deadcrootc_storage_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%deadcrootc_xfer_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%gresp_storage_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%gresp_xfer_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%cpool_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%xsmrpool_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    call update_patch_state( &
         var = this%ctrunc_patch(begp:endp), &
         flux_out_grc_area = conv_cflux(begp:endp))

    if (use_crop) then
       do k = 1, nrepr
          call update_patch_state( &
               var = this%reproductivec_patch(begp:endp, k), &
               flux_out_grc_area = crop_product_cflux(begp:endp))
       end do

       do k = 1, nrepr
          call update_patch_state( &
               var = this%reproductivec_storage_patch(begp:endp, k), &
               flux_out_grc_area = conv_cflux(begp:endp))
       end do

       do k = 1, nrepr
          call update_patch_state( &
               var = this%reproductivec_xfer_patch(begp:endp, k), &
               flux_out_grc_area = conv_cflux(begp:endp))
       end do

       ! This is a negative pool. So any deficit that we haven't repaid gets sucked out
       ! of the atmosphere.
       call update_patch_state( &
            var = this%cropseedc_deficit_patch(begp:endp), &
            flux_out_grc_area = conv_cflux(begp:endp))

       call update_patch_state( &
            var = this%xsmrpool_loss_patch(begp:endp), &
            flux_out_grc_area = conv_cflux(begp:endp))
    end if

  contains
    subroutine update_patch_state(var, flux_out_col_area, flux_out_grc_area, &
         seed, seed_addition)
      ! Wraps call to update_patch_state, in order to remove duplication
      real(r8), intent(inout) :: var( bounds%begp: )
      real(r8), intent(inout), optional :: flux_out_col_area( bounds%begp: )
      real(r8), intent(inout), optional :: flux_out_grc_area( bounds%begp: )
      real(r8), intent(in), optional :: seed( bounds%begp: )
      real(r8), intent(inout), optional :: seed_addition( bounds%begp: )

      call patch_state_updater%update_patch_state(bounds, &
         num_soilp_with_inactive, filter_soilp_with_inactive, &
         var = var, &
         flux_out_col_area = flux_out_col_area, &
         flux_out_grc_area = flux_out_grc_area, &
         seed = seed, &
         seed_addition = seed_addition)
    end subroutine update_patch_state

  end subroutine DynamicPatchAdjustments

end module CNVegCarbonStateType
