module dynSubgridDriverMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! High-level routines for dynamic subgrid areas (prescribed transient Patches, CNDV, and
  ! dynamic landunits).
  !
  ! !USES:
  use decompMod                    , only : bounds_type, BOUNDS_LEVEL_PROC, BOUNDS_LEVEL_CLUMP
  use decompMod                    , only : get_proc_clumps, get_clump_bounds
  use dynSubgridControlMod         , only : get_flanduse_timeseries
  use dynSubgridControlMod         , only : get_do_transient_pfts, get_do_transient_crops
  use dynSubgridControlMod         , only : get_do_harvest
  use dynPriorWeightsMod           , only : prior_weights_type
  use dynPatchStateUpdaterMod      , only : patch_state_updater_type
  use dynColumnStateUpdaterMod     , only : column_state_updater_type
  use dynpftFileMod                , only : dynpft_init, dynpft_interp
  use dyncropFileMod               , only : dyncrop_init, dyncrop_interp
  use dynHarvestMod                , only : dynHarvest_init, dynHarvest_interp
  use dynLandunitAreaMod           , only : update_landunit_weights
  use subgridWeightsMod            , only : compute_higher_order_weights, set_subgrid_diagnostic_fields
  use reweightMod                  , only : reweight_wrapup
  use glcBehaviorMod               , only : glc_behavior_type
  use UrbanParamsType              , only : urbanparams_type
  use CanopyStateType              , only : canopystate_type
  use CNVegetationFacade           , only : cn_vegetation_type
  use SoilBiogeochemStateType      , only : soilBiogeochem_state_type
  use SoilBiogeochemCarbonFluxType , only : soilBiogeochem_carbonflux_type
  use SoilBiogeochemCarbonStateType, only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemNitrogenStateType, only : soilbiogeochem_nitrogenstate_type
  use ch4Mod,                        only : ch4_type
  use EnergyFluxType               , only : energyflux_type
  use PhotosynthesisMod            , only : photosyns_type
  use SoilHydrologyType            , only : soilhydrology_type  
  use SoilStateType                , only : soilstate_type
  use WaterFluxBulkType                , only : waterfluxbulk_type
  use WaterStateBulkType               , only : waterstatebulk_type
  use WaterDiagnosticBulkType               , only : waterdiagnosticbulk_type
  use WaterBalanceType               , only : waterbalance_type
  use TemperatureType              , only : temperature_type
  use CropType                     , only : crop_type
  use glc2lndMod                   , only : glc2lnd_type
  use filterMod                    , only : filter, filter_inactive_and_active
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  !
  public :: dynSubgrid_init                  ! initialize transient land cover
  public :: dynSubgrid_driver                ! top-level driver for transient land cover
  public :: dynSubgrid_wrapup_weight_changes ! reconcile various variables after subgrid weights change

  !
  ! !PRIVATE TYPES:

  ! saved weights from before the subgrid weight updates
  type(prior_weights_type), target :: prior_weights

  ! object used to update patch-level states after subgrid weight updates
  type(patch_state_updater_type), target :: patch_state_updater

  ! object used to update column-level states after subgrid weight updates
  type(column_state_updater_type), target :: column_state_updater
  !---------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dynSubgrid_init(bounds_proc, glc_behavior, crop_inst)
    !
    ! !DESCRIPTION:
    ! Initialize objects needed for prescribed transient PFTs, CNDV, and/or dynamic
    ! landunits.
    !
    ! Also sets initial subgrid weight for aspects prescribed from file (transient PFTs
    ! and transient crops). These initial weights will be overwritten in a restart run,
    ! or any other run that starts up from initial conditions (except startup runs that
    ! use init_interp).
    !
    ! This should be called from initialization, after dynSubgridControl is initialized. 
    !
    ! Note that dynpft_init needs to be called from outside any loops over clumps - so
    ! this routine needs to be called from outside any loops over clumps.
    !
    ! !ARGUMENTS:
    type(bounds_type)       , intent(in)    :: bounds_proc ! processor-level bounds
    type(glc_behavior_type) , intent(in)    :: glc_behavior
    type(crop_type)         , intent(inout) :: crop_inst
    !
    ! !LOCAL VARIABLES:
    integer           :: nclumps      ! number of clumps on this processor
    integer           :: nc           ! clump index
    type(bounds_type) :: bounds_clump ! clump-level bounds
    character(len=*), parameter :: subname = 'dynSubgrid_init'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds_proc%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    nclumps = get_proc_clumps()

    prior_weights = prior_weights_type(bounds_proc)
    patch_state_updater = patch_state_updater_type(bounds_proc)
    column_state_updater = column_state_updater_type(bounds_proc, nclumps)

    ! Initialize stuff for prescribed transient Patches
    if (get_do_transient_pfts()) then
       call dynpft_init(bounds_proc, dynpft_filename=get_flanduse_timeseries())
    end if

    ! Initialize stuff for prescribed transient crops
    if (get_do_transient_crops()) then
       call dyncrop_init(bounds_proc, dyncrop_filename=get_flanduse_timeseries())
    end if

    ! Initialize stuff for harvest. Note that, currently, the harvest data are on the
    ! flanduse_timeseries file. However, this could theoretically be changed so that the
    ! harvest data were separated from the pftdyn data, allowing them to differ in the
    ! years over which they apply.
    if (get_do_harvest()) then
       call dynHarvest_init(bounds_proc, harvest_filename=get_flanduse_timeseries())
    end if

    ! ------------------------------------------------------------------------
    ! Set initial subgrid weights for aspects that are read from file. This is relevant
    ! for cold start and use_init_interp-based initialization.
    ! ------------------------------------------------------------------------

    if (get_do_transient_pfts()) then
       call dynpft_interp(bounds_proc)
    end if

    if (get_do_transient_crops()) then
       call dyncrop_interp(bounds_proc, crop_inst)
    end if

    ! (We don't bother calling dynHarvest_interp, because the harvest information isn't
    ! needed until the run loop. Harvest has nothing to do with subgrid weights, and in
    ! some respects doesn't even really belong in this module at all.)

    ! The following is only needed if there were actually weight changes due to the above
    ! interp calls, but it doesn't hurt to always run this code:

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1, nclumps
       call get_clump_bounds(nc, bounds_clump)

       call dynSubgrid_wrapup_weight_changes(bounds_clump, glc_behavior)
    end do
    !$OMP END PARALLEL DO

  end subroutine dynSubgrid_init

  !-----------------------------------------------------------------------
  subroutine dynSubgrid_driver(bounds_proc,                                            &
       urbanparams_inst, soilstate_inst, soilhydrology_inst,           &
       waterstatebulk_inst, waterdiagnosticbulk_inst, waterbalancebulk_inst, &
       waterfluxbulk_inst, temperature_inst, energyflux_inst,             &
       canopystate_inst, photosyns_inst, crop_inst, glc2lnd_inst, bgc_vegetation_inst,          &
       soilbiogeochem_state_inst, soilbiogeochem_carbonstate_inst, &
       c13_soilbiogeochem_carbonstate_inst, c14_soilbiogeochem_carbonstate_inst,       &
       soilbiogeochem_nitrogenstate_inst, soilbiogeochem_carbonflux_inst, ch4_inst, &
       glc_behavior)
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
    use clm_varctl           , only : use_cn, use_fates
    use dynInitColumnsMod    , only : initialize_new_columns
    use dynConsBiogeophysMod , only : dyn_hwcontent_init, dyn_hwcontent_final
    use dynEDMod             , only : dyn_ED
    !
    ! !ARGUMENTS:
    type(bounds_type)                    , intent(in)    :: bounds_proc  ! processor-level bounds
    type(urbanparams_type)               , intent(in)    :: urbanparams_inst
    type(soilstate_type)                 , intent(in)    :: soilstate_inst
    type(soilhydrology_type)             , intent(inout) :: soilhydrology_inst
    type(waterstatebulk_type)                , intent(inout) :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)                , intent(inout) :: waterdiagnosticbulk_inst
    type(waterbalance_type)                , intent(inout) :: waterbalancebulk_inst
    type(waterfluxbulk_type)                 , intent(inout) :: waterfluxbulk_inst
    type(temperature_type)               , intent(inout) :: temperature_inst
    type(energyflux_type)                , intent(inout) :: energyflux_inst
    type(canopystate_type)               , intent(inout) :: canopystate_inst
    type(photosyns_type)                 , intent(inout) :: photosyns_inst
    type(crop_type)                      , intent(inout) :: crop_inst
    type(glc2lnd_type)                   , intent(inout) :: glc2lnd_inst
    type(cn_vegetation_type)             , intent(inout) :: bgc_vegetation_inst
    type(soilbiogeochem_state_type)      , intent(in)    :: soilbiogeochem_state_inst
    type(soilbiogeochem_carbonstate_type), intent(inout) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonstate_type), intent(inout) :: c13_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonstate_type), intent(inout) :: c14_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_nitrogenstate_type), intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_carbonflux_type) , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(ch4_type)                       , intent(inout) :: ch4_inst
    type(glc_behavior_type)              , intent(in)    :: glc_behavior
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

       call dyn_hwcontent_init(bounds_clump, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            filter(nc)%num_lakec, filter(nc)%lakec, &
            urbanparams_inst, soilstate_inst, soilhydrology_inst, &
            waterstatebulk_inst, waterdiagnosticbulk_inst, waterbalancebulk_inst, &
            waterfluxbulk_inst, temperature_inst, energyflux_inst)

       call prior_weights%set_prior_weights(bounds_clump)
       call patch_state_updater%set_old_weights(bounds_clump)
       call column_state_updater%set_old_weights(bounds_clump)
    end do
    !$OMP END PARALLEL DO

    ! ==========================================================================
    ! Do land cover change that requires I/O, and thus must be outside a threaded region
    ! ==========================================================================

    if (get_do_transient_pfts()) then
       call dynpft_interp(bounds_proc)
    end if

    if (get_do_transient_crops()) then
       call dyncrop_interp(bounds_proc,crop_inst)
    end if

    if (get_do_harvest()) then
       call dynHarvest_interp(bounds_proc)
    end if

    ! ==========================================================================
    ! Do land cover change that does not require I/O
    ! ==========================================================================

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1, nclumps
       call get_clump_bounds(nc, bounds_clump)

       call bgc_vegetation_inst%UpdateSubgridWeights(bounds_clump)
       
       if (use_fates) then
          call dyn_ED(bounds_clump)
       end if

       call glc2lnd_inst%update_glc2lnd_fracs( &
            bounds = bounds_clump)

       ! ========================================================================
       ! Do wrapup stuff after land cover change
       !
       ! Everything following this point in this loop only needs to be called if we have
       ! actually changed some weights in this time step. However, it doesn't do any harm
       ! (other than a small performance hit) to call this stuff all the time, so we do so
       ! for simplicity and safety.
       !
       ! NOTE(wjs, 2017-02-24) I'm not positive that the above paragraph is 100% true. It
       ! is at least *mostly* true, but there may be some subtleties, like resetting of
       ! some variables, that are needed even in (some) time steps where we haven't
       ! changed weights.
       ! ========================================================================

       call dynSubgrid_wrapup_weight_changes(bounds_clump, glc_behavior)

       call patch_state_updater%set_new_weights(bounds_clump)
       call column_state_updater%set_new_weights(bounds_clump, nc)

       call set_subgrid_diagnostic_fields(bounds_clump)

       call initialize_new_columns(bounds_clump, &
            prior_weights%cactive(bounds_clump%begc:bounds_clump%endc), &
            temperature_inst, waterstatebulk_inst, soilhydrology_inst)

       call dyn_hwcontent_final(bounds_clump, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            filter(nc)%num_lakec, filter(nc)%lakec, &
            urbanparams_inst, soilstate_inst, soilhydrology_inst, &
            waterstatebulk_inst, waterdiagnosticbulk_inst, waterbalancebulk_inst, &
            waterfluxbulk_inst, temperature_inst, energyflux_inst)

       if (use_cn) then
          call bgc_vegetation_inst%DynamicAreaConservation(bounds_clump, nc, &
               filter_inactive_and_active(nc)%num_soilp, filter_inactive_and_active(nc)%soilp, &
               filter_inactive_and_active(nc)%num_soilc, filter_inactive_and_active(nc)%soilc, &
               prior_weights, patch_state_updater, column_state_updater, &
               canopystate_inst, photosyns_inst, &
               soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst, &
               c13_soilbiogeochem_carbonstate_inst, c14_soilbiogeochem_carbonstate_inst, &
               soilbiogeochem_nitrogenstate_inst, ch4_inst, soilbiogeochem_state_inst)
       end if

    end do
    !$OMP END PARALLEL DO

  end subroutine dynSubgrid_driver

  !-----------------------------------------------------------------------
  subroutine dynSubgrid_wrapup_weight_changes(bounds_clump, glc_behavior)
    !
    ! !DESCRIPTION:
    ! Reconcile various variables after subgrid weights change
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type)       , intent(in) :: bounds_clump ! clump-level bounds
    type(glc_behavior_type) , intent(in) :: glc_behavior
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'dynSubgrid_wrapup_weight_changes'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds_clump%level == BOUNDS_LEVEL_CLUMP, subname // ': argument must be CLUMP-level bounds')

    call update_landunit_weights(bounds_clump)

    call compute_higher_order_weights(bounds_clump)

    ! Here: filters are re-made
    !
    ! This call requires clump-level bounds, which is why we need to ensure that the
    ! argument to this routine is clump-level bounds
    call reweight_wrapup(bounds_clump, glc_behavior)

  end subroutine dynSubgrid_wrapup_weight_changes


end module dynSubgridDriverMod
