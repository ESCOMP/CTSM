module dynSubgridDriverMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! High-level routines for dynamic subgrid areas (prescribed transient PFTs, CNDV, and
  ! dynamic landunits).
  !
  ! !USES:
  use clmtype
  use dynPriorWeightsMod , only : prior_weights_type
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save

  public :: dynSubgrid_init             ! initialize transient land cover
  public :: dynSubgrid_driver           ! top-level driver for transient land cover

  !
  ! !PRIVATE TYPES:
  type(prior_weights_type) :: prior_weights ! saved weights from before the subgrid weight updates
  !---------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dynSubgrid_init(bounds)
    !
    ! !DESCRIPTION:
    ! Determine initial subgrid weights for prescribed transient PFTs, CNDV, and/or
    ! dynamic landunits. Note that these weights will be overwritten in a restart run.
    !
    ! This should be called from initialization. 
    !
    ! Note that dynpft_init / dynpft_interp need to be called from outside any loops over
    ! clumps - so this routine needs to be called from outside any loops over clumps.
    !
    ! !USES:
    use clm_varctl    , only : fpftdyn, use_cndv
    use decompMod     , only : bounds_type, BOUNDS_LEVEL_PROC
    use dynpftFileMod , only : dynpft_init
    use dynHarvestMod , only : dynHarvest_init
    use dynCNDVMod    , only : dynCNDV_init
    use reweightMod   , only : compute_higher_order_weights
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! processor-level bounds
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'dynSubgrid_init'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    prior_weights = prior_weights_type(bounds)

    ! Initialize stuff for prescribed transient PFTs
    if (fpftdyn /= ' ') then
       call dynpft_init(bounds)
    end if

    ! Initialize stuff for harvest (currently shares the fpftdyn file)
    if (fpftdyn /= ' ') then
       call dynHarvest_init(bounds)
    end if

    if (use_cndv) then
       call dynCNDV_init(bounds)
    end if

    call compute_higher_order_weights(bounds)
    
  end subroutine dynSubgrid_init

  !-----------------------------------------------------------------------
  subroutine dynSubgrid_driver(bounds_proc)
    !
    ! !DESCRIPTION:
    ! Update subgrid weights for prescribed transient PFTs, CNDV, and/or dynamic
    ! landunits. Also do related adjustments (water & energy, carbon & nitrogen).
    !
    ! This should be called every time step in CLM's run loop.
    !
    ! Note that this routine operates partly at the proc-level (outside an OMP region),
    ! and partly at the clump level (inside OMP regions). Thus, this must be called from
    ! OUTSIDE any loops over clumps in the driver.
    !
    ! !USES:
    use clm_varctl           , only : fpftdyn, use_cndv, use_cn
    use decompMod            , only : bounds_type, get_proc_clumps, get_clump_bounds, &
                                      BOUNDS_LEVEL_PROC
    use dynLandunitAreaMod   , only : update_landunit_weights
    use dynConsBiogeophysMod , only : dyn_hwcontent_init, dyn_hwcontent_final
    use dynConsBiogeochemMod , only : dyn_cnbal_pft
    use dynpftFileMod        , only : dynpft_interp
    use dynHarvestMod        , only : dynHarvest_interp
    use dynCNDVMod           , only : dynCNDV_interp
    use reweightMod          , only : compute_higher_order_weights, reweightWrapup
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds_proc  ! processor-level bounds
    !
    ! !LOCAL VARIABLES:
    integer           :: nclumps      ! number of clumps on this processor
    integer           :: nc           ! clump index
    type(bounds_type) :: bounds_clump ! clump-level bounds

    character(len=*), parameter :: subname = 'dynSubgrid_driver'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds_proc%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    nclumps = get_proc_clumps()
    
    ! ==========================================================================
    ! Do initialization, prior to land cover change
    ! ==========================================================================

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1, nclumps
       call get_clump_bounds(nc, bounds_clump)

       call dyn_hwcontent_init(bounds_clump)
       call prior_weights%set_prior_weights(bounds_clump)
    end do
    !$OMP END PARALLEL DO

    ! ==========================================================================
    ! Do land cover change that requires I/O, and thus must be outside a threaded region
    ! ==========================================================================

    if (fpftdyn /= ' ') then
       call dynpft_interp(bounds_proc)
    end if

    if (fpftdyn /= ' ') then
       call dynHarvest_interp(bounds_proc)
    end if

    ! ==========================================================================
    ! Do everything else related to land cover change
    ! ==========================================================================

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1, nclumps
       call get_clump_bounds(nc, bounds_clump)

       if (use_cndv) then
          call dynCNDV_interp(bounds_clump)
       end if

       ! Everything following this point in this loop only needs to be called if we have
       ! actually changed some weights in this time step. However, it doesn't hurt for
       ! these things to be called all the time (except for a minor performance hit), so
       ! I'm leaving them outside any conditionals for simplicity

       call update_landunit_weights(bounds_clump)

       call compute_higher_order_weights(bounds_clump)

       call reweightWrapup(bounds_clump)

       call dyn_hwcontent_final(bounds_clump)

       if (use_cn) then
          call dyn_cnbal_pft(bounds_clump, prior_weights)
       end if

    end do
    !$OMP END PARALLEL DO

  end subroutine dynSubgrid_driver

end module dynSubgridDriverMod
