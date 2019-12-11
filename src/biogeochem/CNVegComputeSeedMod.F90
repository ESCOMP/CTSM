module CNVegComputeSeedMod

  !-----------------------------------------------------------------------
  ! Module to compute seed amounts for new patch areas
  !
  ! !USES:
#include "shr_assert.h"

  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use pftconMod      , only : pftcon, noveg
  use clm_varcon     , only : c3_r2, c4_r2, c14ratio
  use clm_varctl     , only : iulog
  use PatchType      , only : patch
  use abortutils     , only : endrun
  use CNSpeciesMod   , only : CN_SPECIES_C12, CN_SPECIES_C13, CN_SPECIES_C14, CN_SPECIES_N
  !
  ! !PUBLIC ROUTINES:
  implicit none
  private

  public :: ComputeSeedAmounts

  ! !PRIVATE ROUTINES:

  private :: SpeciesTypeMultiplier
  private :: LeafProportions  ! compute leaf proportions (leaf, storage and xfer)

  ! !PRIVATE DATA:

  integer, parameter :: COMPONENT_LEAF = 1
  integer, parameter :: COMPONENT_DEADWOOD = 2

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine ComputeSeedAmounts(bounds, &
       num_soilp_with_inactive, filter_soilp_with_inactive, &
       species, &
       leafc_seed, deadstemc_seed, &
       leaf_patch, leaf_storage_patch, leaf_xfer_patch, &
       compute_here_patch, ignore_current_state_patch, &
       seed_leaf_patch, seed_leaf_storage_patch, seed_leaf_xfer_patch, &
       seed_deadstem_patch)
    !
    ! !DESCRIPTION:
    ! Compute seed amounts for patches that increase in area, for various variables, for
    ! the given species (c12, c13, c14 or n).
    !
    ! The output variables are only set for patches inside the filter, where
    ! compute_here_patch is true; for other patches, they remain at their original values.
    !
    ! Note that, regardless of the species, leafc_seed and deadstemc_seed are specified
    ! in terms of gC/m2; these amounts are converted to the amount of the given species
    ! here.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type)              , intent(in) :: bounds
    integer                        , intent(in) :: num_soilp_with_inactive       ! number of points in filter
    integer                        , intent(in) :: filter_soilp_with_inactive(:) ! soil patch filter that includes inactive points
    integer                        , intent(in) :: species                       ! which C/N species we're operating on; should be one of the values in CNSpeciesMod
    real(r8)                       , intent(in) :: leafc_seed                    ! seed amount for leaf C
    real(r8)                       , intent(in) :: deadstemc_seed                ! seed amount for deadstem C
    real(r8)                       , intent(in) :: leaf_patch( bounds%begp: )   ! current leaf C or N content (g/m2)
    real(r8)                       , intent(in) :: leaf_storage_patch( bounds%begp: ) ! current leaf C or N storage content (g/m2)
    real(r8)                       , intent(in) :: leaf_xfer_patch( bounds%begp: )    ! current leaf C or N xfer content (g/m2)

    ! whether to compute outputs for each patch
    logical, intent(in) :: compute_here_patch( bounds%begp: )

    ! If ignore_current_state is true, then use default leaf proportions rather than
    ! proportions based on current state.
    logical, intent(in) :: ignore_current_state_patch( bounds%begp: )

    real(r8), intent(inout) :: seed_leaf_patch( bounds%begp: ) ! seed amount for leaf itself for this species (g/m2)
    real(r8), intent(inout) :: seed_leaf_storage_patch( bounds%begp: ) ! seed amount for leaf storage for this species (g/m2)
    real(r8), intent(inout) :: seed_leaf_xfer_patch( bounds%begp: ) ! seed amount for leaf xfer for this species (g/m2)
    real(r8), intent(inout) :: seed_deadstem_patch( bounds%begp: ) ! seed amount for deadstem for this species (g/m2)
    !
    ! !LOCAL VARIABLES:
    integer :: fp, p
    integer :: begp, endp
    real(r8) :: my_leaf_seed
    real(r8) :: my_deadstem_seed
    integer  :: pft_type
    real(r8) :: pleaf
    real(r8) :: pstor
    real(r8) :: pxfer

    character(len=*), parameter :: subname = 'ComputeSeedAmounts'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    SHR_ASSERT_ALL_FL((ubound(leaf_patch) == (/endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(leaf_storage_patch) == (/endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(leaf_xfer_patch) == (/endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(compute_here_patch) == (/endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(ignore_current_state_patch) == (/endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(seed_leaf_patch) == (/endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(seed_leaf_storage_patch) == (/endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(seed_leaf_xfer_patch) == (/endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(seed_deadstem_patch) == (/endp/)), sourcefile, __LINE__)


    do fp = 1, num_soilp_with_inactive
       p = filter_soilp_with_inactive(fp)

       if (compute_here_patch(p)) then

          my_leaf_seed = 0._r8
          my_deadstem_seed = 0._r8

          pft_type = patch%itype(p)

          call LeafProportions( &
               ignore_current_state = ignore_current_state_patch(p), &
               pft_type = pft_type, &
               leaf = leaf_patch(p), &
               leaf_storage = leaf_storage_patch(p), &
               leaf_xfer = leaf_xfer_patch(p), &
               pleaf = pleaf, &
               pstorage = pstor, &
               pxfer = pxfer)

          if (pft_type /= noveg) then
             my_leaf_seed = leafc_seed * &
                  SpeciesTypeMultiplier(species, pft_type, COMPONENT_LEAF)
             if (pftcon%woody(pft_type) == 1._r8) then
                my_deadstem_seed = deadstemc_seed * &
                     SpeciesTypeMultiplier(species, pft_type, COMPONENT_DEADWOOD)
             end if
          end if

          seed_leaf_patch(p) = my_leaf_seed * pleaf
          seed_leaf_storage_patch(p) = my_leaf_seed * pstor
          seed_leaf_xfer_patch(p) = my_leaf_seed * pxfer
          seed_deadstem_patch(p) = my_deadstem_seed
       end if

    end do

  end subroutine ComputeSeedAmounts


  !-----------------------------------------------------------------------
  function SpeciesTypeMultiplier(species, pft_type, component) result(multiplier)
    !
    ! !DESCRIPTION:
    ! Returns a multiplier based on the species type. This multiplier is
    ! meant to be applied to some state variable expressed in terms of g C, translating
    ! this value into an appropriate value for c13, c14 or n.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8) :: multiplier  ! function result
    integer, intent(in) :: species ! which C/N species we're operating on; should be one of the values in CNSpeciesMod
    integer, intent(in) :: pft_type
    integer, intent(in) :: component ! which plant component; should be one of the COMPONENT_* parameters defined in this module
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'SpeciesTypeMultiplier'
    !-----------------------------------------------------------------------

    select case (species)
    case (CN_SPECIES_C12)
       multiplier = 1._r8

    case (CN_SPECIES_C13)
       if (pftcon%c3psn(pft_type) == 1._r8) then
          multiplier = c3_r2
       else
          multiplier = c4_r2
       end if

    case (CN_SPECIES_C14)
       ! 14c state is initialized assuming initial "modern" 14C of 1.e-12
       multiplier = c14ratio

    case (CN_SPECIES_N)
       select case (component)
       case (COMPONENT_LEAF)
          multiplier = 1._r8 / pftcon%leafcn(pft_type)
       case (COMPONENT_DEADWOOD)
          multiplier = 1._r8 / pftcon%deadwdcn(pft_type)
       case default
          write(iulog,*) subname//' ERROR: unknown component: ', component
          call endrun(subname//': unknown component')
       end select

    case default
       write(iulog,*) subname//' ERROR: unknown species: ', species
       call endrun(subname//': unknown species')
    end select

  end function SpeciesTypeMultiplier


  !-----------------------------------------------------------------------
  subroutine LeafProportions(ignore_current_state, &
       pft_type, &
       leaf, leaf_storage, leaf_xfer, &
       pleaf, pstorage, pxfer)
    !
    ! !DESCRIPTION:
    ! Compute leaf proportions (leaf, storage and xfer)
    !
    ! If ignore_current_state is true, then use default proportions rather than
    ! proportions based on current state. (Also use default proportions if total leaf mass
    ! is 0 for this patch.)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    logical, intent(in) :: ignore_current_state  ! see comment above
    integer , intent(in) :: pft_type
    real(r8), intent(in) :: leaf  ! g/m2 leaf C or N
    real(r8), intent(in) :: leaf_storage ! g/m2 leaf C or N storage
    real(r8), intent(in) :: leaf_xfer ! g/m2 leaf C or N transfer

    real(r8), intent(out) :: pleaf  ! proportion in leaf itself
    real(r8), intent(out) :: pstorage ! proportion in leaf storage
    real(r8), intent(out) :: pxfer ! proportion in leaf xfer
    !
    ! !LOCAL VARIABLES:
    real(r8) :: tot_leaf

    character(len=*), parameter :: subname = 'LeafProportions'
    !-----------------------------------------------------------------------

    tot_leaf = leaf + leaf_storage + leaf_xfer
    pleaf = 0._r8
    pstorage = 0._r8
    pxfer = 0._r8

    if (tot_leaf == 0._r8 .or. ignore_current_state) then
       if (pftcon%evergreen(pft_type) == 1._r8) then
          pleaf = 1._r8
       else
          pstorage = 1._r8
       end if
    else
       pleaf = leaf/tot_leaf
       pstorage = leaf_storage/tot_leaf
       pxfer = leaf_xfer/tot_leaf
    end if

  end subroutine LeafProportions

end module CNVegComputeSeedMod
