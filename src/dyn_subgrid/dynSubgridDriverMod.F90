module dynSubgridDriverMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! High-level routines for dynamic subgrid areas (prescribed transient Patches, CNDV, and
  ! dynamic landunits).
  !
  ! !USES:
  use dynSubgridControlMod         , only : get_flanduse_timeseries
  use dynSubgridControlMod         , only : get_do_transient_pfts, get_do_transient_crops
  use dynSubgridControlMod         , only : get_do_harvest
  use dynPriorWeightsMod           , only : prior_weights_type
  use dynPatchStateUpdaterMod      , only : patch_state_updater_type
  use dynColumnStateUpdaterMod     , only : column_state_updater_type
  use UrbanParamsType              , only : urbanparams_type
  use CanopyStateType              , only : canopystate_type
  use CNVegetationFacade           , only : cn_vegetation_type
  use SoilBiogeochemStateType      , only : soilBiogeochem_state_type
  use SoilBiogeochemCarbonFluxType , only : soilBiogeochem_carbonflux_type
  use SoilBiogeochemCarbonStateType, only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemNitrogenStateType, only : soilbiogeochem_nitrogenstate_type
  use ch4Mod,                        only : ch4_type
  use EnergyFluxType               , only : energyflux_type
  use LakeStateType                , only : lakestate_type
  use PhotosynthesisMod            , only : photosyns_type
  use SoilHydrologyType            , only : soilhydrology_type  
  use SoilStateType                , only : soilstate_type
  use WaterfluxType                , only : waterflux_type
  use WaterstateType               , only : waterstate_type
  use TemperatureType              , only : temperature_type
  use glc2lndMod                   , only : glc2lnd_type
  use filterMod                    , only : filter_inactive_and_active
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  !
  public :: dynSubgrid_init             ! initialize transient land cover
  public :: dynSubgrid_driver           ! top-level driver for transient land cover
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
  subroutine dynSubgrid_init(bounds)
    !
    ! !DESCRIPTION:
    ! Initialize objects needed for prescribed transient PFTs, CNDV, and/or dynamic
    ! landunits. 
    !
    ! This should be called from initialization, after dynSubgridControl is initialized. 
    !
    ! Note that no subgrid weights are updated here - so in initialization, weights will
    ! stay at the values read from the surface dataset, then there will be a potentially
    ! large adjustment of weights in the first time step. The reason why no weights are
    ! updated here is that prognostic subgrid weight information (e.g., glacier cover
    ! from CISM) is not available in initialization; thus, for consistency, we also avoid
    ! applying prescribed transient weights in initialization.
    !
    ! However, the above note is only relevant for a cold start run: For a restart run or
    ! a run with initial conditions, subgrid weights will be read from the restart file
    ! after this routine is called. This is also relevant for a run that uses
    ! interpolated initial conditions (use_init_interp = .true.), since subgrid weights
    ! are not interpolated, so start at their cold start values.
    !
    ! Note that dynpft_init needs to be called from outside any loops over clumps - so
    ! this routine needs to be called from outside any loops over clumps.
    !
    ! !USES:
    use clm_varctl        , only : use_cndv
    use decompMod         , only : bounds_type, BOUNDS_LEVEL_PROC
    use dynpftFileMod     , only : dynpft_init
    use dyncropFileMod    , only : dyncrop_init
    use dynHarvestMod     , only : dynHarvest_init
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds     ! processor-level bounds
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'dynSubgrid_init'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    prior_weights = prior_weights_type(bounds)
    patch_state_updater = patch_state_updater_type(bounds)
    column_state_updater = column_state_updater_type(bounds)

    ! Initialize stuff for prescribed transient Patches
    if (get_do_transient_pfts()) then
       call dynpft_init(bounds, dynpft_filename=get_flanduse_timeseries())
    end if

    ! Initialize stuff for prescribed transient crops
    if (get_do_transient_crops()) then
       call dyncrop_init(bounds, dyncrop_filename=get_flanduse_timeseries())
    end if

    ! Initialize stuff for harvest. Note that, currently, the harvest data are on the
    ! flanduse_timeseries file. However, this could theoretically be changed so that the
    ! harvest data were separated from the pftdyn data, allowing them to differ in the
    ! years over which they apply.
    if (get_do_harvest()) then
       call dynHarvest_init(bounds, harvest_filename=get_flanduse_timeseries())
    end if

  end subroutine dynSubgrid_init

  !-----------------------------------------------------------------------
  subroutine dynSubgrid_driver(bounds_proc,                                            &
       urbanparams_inst, soilstate_inst, soilhydrology_inst, lakestate_inst,           &
       waterstate_inst, waterflux_inst, temperature_inst, energyflux_inst,             &
       canopystate_inst, photosyns_inst, glc2lnd_inst, bgc_vegetation_inst,          &
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
    use clm_time_manager     , only : is_first_step
    use clm_varctl           , only : is_cold_start
    use clm_varctl           , only : use_cndv, use_cn, create_glacier_mec_landunit, use_ed
    use decompMod            , only : bounds_type, get_proc_clumps, get_clump_bounds
    use decompMod            , only : BOUNDS_LEVEL_PROC
    use dynLandunitAreaMod   , only : update_landunit_weights
    use dynInitColumnsMod    , only : initialize_new_columns
    use dynConsBiogeophysMod , only : dyn_hwcontent_init, dyn_hwcontent_final
    use dynpftFileMod        , only : dynpft_interp
    use dynCropFileMod       , only : dyncrop_interp
    use dynHarvestMod        , only : dynHarvest_interp
    use dynEDMod             , only : dyn_ED
    use reweightMod          , only : reweight_wrapup
    use subgridWeightsMod    , only : compute_higher_order_weights, set_subgrid_diagnostic_fields
    use glcBehaviorMod       , only : glc_behavior_type
    !
    ! !ARGUMENTS:
    type(bounds_type)                    , intent(in)    :: bounds_proc  ! processor-level bounds
    type(urbanparams_type)               , intent(in)    :: urbanparams_inst
    type(soilstate_type)                 , intent(in)    :: soilstate_inst
    type(soilhydrology_type)             , intent(in)    :: soilhydrology_inst
    type(lakestate_type)                 , intent(in)    :: lakestate_inst
    type(waterstate_type)                , intent(inout) :: waterstate_inst
    type(waterflux_type)                 , intent(inout) :: waterflux_inst
    type(temperature_type)               , intent(inout) :: temperature_inst
    type(energyflux_type)                , intent(inout) :: energyflux_inst
    type(canopystate_type)               , intent(inout) :: canopystate_inst
    type(photosyns_type)                 , intent(inout) :: photosyns_inst
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
    logical           :: first_step_cold_start ! true if this is the first step since cold start

    ! These are used if this is the first step of a cold start
    type(prior_weights_type), target :: new_weights
    type(patch_state_updater_type), target :: patch_state_updater_new_weights
    type(column_state_updater_type), target :: column_state_updater_new_weights

    ! These point to the appropriate prior_weights, patch_state_updater and
    ! column_state_updater instances, depending on whether it's a cold start
    type(prior_weights_type), pointer :: my_prior_weights
    type(patch_state_updater_type), pointer :: my_patch_state_updater
    type(column_state_updater_type), pointer :: my_column_state_updater

    character(len=*), parameter :: subname = 'dynSubgrid_driver'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds_proc%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    nclumps = get_proc_clumps()
    first_step_cold_start = (is_first_step() .and. is_cold_start)

    ! ==========================================================================
    ! Do initialization, prior to land cover change
    ! ==========================================================================

    if (first_step_cold_start) then
       ! These objects need to be constructed outside a loop over clumps
       new_weights = prior_weights_type(bounds_proc)
       patch_state_updater_new_weights = patch_state_updater_type(bounds_proc)
       column_state_updater_new_weights = column_state_updater_type(bounds_proc)
    end if

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1, nclumps
       call get_clump_bounds(nc, bounds_clump)

       call dyn_hwcontent_init(bounds_clump, &
            urbanparams_inst, soilstate_inst, soilhydrology_inst, lakestate_inst, &
            waterstate_inst, waterflux_inst, temperature_inst, energyflux_inst)

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
       call dyncrop_interp(bounds_proc)
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
       
       if (use_ed) then
          call dyn_ED(bounds_clump)
       end if

       if (create_glacier_mec_landunit) then
          call glc2lnd_inst%update_glc2lnd_non_topo( &
               bounds = bounds_clump, &
               glc_behavior = glc_behavior)
       end if

       ! ========================================================================
       ! Do wrapup stuff after land cover change
       !
       ! Everything following this point in this loop only needs to be called if we have
       ! actually changed some weights in this time step. This is also required in the
       ! first time step of the run to update filters to reflect state of CISM
       ! (particularly mask that is past through coupler).
       !
       ! However, it doesn't do any harm (other than a small performance hit) to call
       ! this stuff all the time, so we do so for simplicity and safety.
       ! ========================================================================

       call update_landunit_weights(bounds_clump)

       call compute_higher_order_weights(bounds_clump)

       ! Here: filters are re-made
       call reweight_wrapup(bounds_clump, glc_behavior)

       call patch_state_updater%set_new_weights(bounds_clump)
       call column_state_updater%set_new_weights(bounds_clump)

       call set_subgrid_diagnostic_fields(bounds_clump)

       call initialize_new_columns(bounds_clump, &
            prior_weights%cactive(bounds_clump%begc:bounds_clump%endc), &
            temperature_inst)

       if (first_step_cold_start) then
          ! In the first timestep of a cold start run, we want to avoid doing any state
          ! adjustments. This is because we expect big transients in the first time step,
          ! since transient subgrid weights aren't updated in initialization. To
          ! accomplish this, we set up a prior_weights object and a column_state_updater
          ! object that say that there were no weight updates in this time step - i.e.,
          ! that the old weights are the same as the new weights.
          call new_weights%set_prior_weights(bounds_clump)
          call patch_state_updater_new_weights%set_old_weights(bounds_clump)
          call patch_state_updater_new_weights%set_new_weights(bounds_clump)
          call column_state_updater_new_weights%set_old_weights(bounds_clump)
          call column_state_updater_new_weights%set_new_weights(bounds_clump)
          my_prior_weights => new_weights
          my_patch_state_updater => patch_state_updater_new_weights
          my_column_state_updater => column_state_updater_new_weights
       else
          my_prior_weights => prior_weights
          my_patch_state_updater => patch_state_updater
          my_column_state_updater => column_state_updater
       end if

       call dyn_hwcontent_final(bounds_clump, first_step_cold_start, &
            urbanparams_inst, soilstate_inst, soilhydrology_inst, lakestate_inst, &
            waterstate_inst, waterflux_inst, temperature_inst, energyflux_inst)

       if (use_cn) then
          call bgc_vegetation_inst%DynamicAreaConservation(bounds_clump, &
               filter_inactive_and_active(nc)%num_soilp, filter_inactive_and_active(nc)%soilp, &
               my_prior_weights, my_patch_state_updater, my_column_state_updater, &
               canopystate_inst, photosyns_inst, &
               soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst, &
               c13_soilbiogeochem_carbonstate_inst, c14_soilbiogeochem_carbonstate_inst, &
               soilbiogeochem_nitrogenstate_inst, ch4_inst, soilbiogeochem_state_inst)
       end if

    end do
    !$OMP END PARALLEL DO

  end subroutine dynSubgrid_driver

end module dynSubgridDriverMod
