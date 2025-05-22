module clm_driver

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module provides the main CLM driver physics calling sequence.  Most
  ! computations occurs over ``clumps'' of gridcells (and associated subgrid
  ! scale entities) assigned to each MPI process. Computation is further
  ! parallelized by looping over clumps on each process using shared memory OpenMP.
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use clm_varctl             , only : iulog, use_fates, use_fates_sp, use_fates_bgc
  use clm_varctl             , only : use_cn, use_lch4, use_noio, use_c13, use_c14
  use CNSharedParamsMod      , only : use_matrixcn
  use clm_varctl             , only : use_crop, irrigate, ndep_from_cpl
  use clm_varctl             , only : use_soil_moisture_streams
  use clm_varctl             , only : use_cropcal_streams
  use clm_time_manager       , only : get_nstep, is_beg_curr_day, is_beg_curr_year
  use clm_time_manager       , only : get_prev_date, is_first_step
  use clm_varpar             , only : nlevsno, nlevgrnd
  use clm_varorb             , only : obliqr
  use spmdMod                , only : masterproc, mpicom
  use decompMod              , only : get_proc_clumps, get_clump_bounds, get_proc_bounds, bounds_type
  use filterMod              , only : filter, filter_inactive_and_active
  use filterMod              , only : setExposedvegpFilter
  use histFileMod            , only : hist_update_hbuf, hist_htapes_wrapup
  use restFileMod            , only : restFile_write, restFile_filename
  use abortutils             , only : endrun
  !
  use dynSubgridDriverMod    , only : dynSubgrid_driver, dynSubgrid_wrapup_weight_changes
  use BalanceCheckMod        , only : WaterGridcellBalance, BeginWaterColumnBalance, BalanceCheck
  !
  use BiogeophysPreFluxCalcsMod  , only : BiogeophysPreFluxCalcs
  use SurfaceHumidityMod     , only : CalculateSurfaceHumidity
  use UrbanTimeVarType       , only : urbantv_type
  use SoilTemperatureMod     , only : SoilTemperature
  use LakeTemperatureMod     , only : LakeTemperature
  !
  use BareGroundFluxesMod    , only : BareGroundFluxes
  use CanopyFluxesMod        , only : CanopyFluxes
  use SoilFluxesMod          , only : SoilFluxes ! (formerly Biogeophysics2Mod)
  use UrbanFluxesMod         , only : UrbanFluxes
  use LakeFluxesMod          , only : LakeFluxes
  !
  use HydrologyNoDrainageMod , only : CalcAndWithdrawIrrigationFluxes, HandleNewSnow, HydrologyNoDrainage ! (formerly Hydrology2Mod)
  use HydrologyDrainageMod   , only : HydrologyDrainage   ! (formerly Hydrology2Mod)
  use CanopyHydrologyMod     , only : CanopyInterceptionAndThroughfall
  use SurfaceWaterMod        , only : UpdateFracH2oSfc
  use LakeHydrologyMod       , only : LakeHydrology
  use SoilWaterMovementMod   , only : use_aquifer_layer
  !
  use AerosolMod             , only : AerosolMasses
  use SnowSnicarMod          , only : SnowAge_grain
  use SurfaceAlbedoMod       , only : SurfaceAlbedo
  use SurfaceAlbedoMod       , only : UpdateZenithAngles
  use UrbanAlbedoMod         , only : UrbanAlbedo
  !
  use SurfaceRadiationMod    , only : SurfaceRadiation, CanopySunShadeFracs
  use UrbanRadiationMod      , only : UrbanRadiation
  !
  use SoilBiogeochemVerticalProfileMod   , only : SoilBiogeochemVerticalProfile
  use SatellitePhenologyMod  , only : CalcSatellitePhenologyTimeInterp, interpMonthlyVeg, UpdateSatellitePhenologyCanopy
  use ndepStreamMod          , only : ndep_interp
  use cropcalStreamMod       , only : cropcal_advance, cropcal_interp
  use ch4Mod                 , only : ch4, ch4_init_gridcell_balance_check, ch4_init_column_balance_check
  use VOCEmissionMod         , only : VOCEmission
  !
  use filterMod              , only : setFilters
  !
  use atm2lndMod             , only : downscale_forcings, set_atm2lnd_water_tracers
  use lnd2atmMod             , only : lnd2atm
  use lnd2glcMod             , only : lnd2glc_type
  !
  use shr_drydep_mod         , only : n_drydep
  use DryDepVelocity         , only : depvel_compute
  !
  use DaylengthMod           , only : UpdateDaylength
  use perf_mod
  !
  use GridcellType           , only : grc
  use LandunitType           , only : lun
  use ColumnType             , only : col
  use PatchType              , only : patch
  use clm_instMod
  use SoilMoistureStreamMod  , only : PrescribedSoilMoistureInterp, PrescribedSoilMoistureAdvance
  use SoilBiogeochemDecompCascadeConType , only : no_soil_decomp, decomp_method
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: clm_drv            ! Main clm driver
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: clm_drv_patch2col
  private :: clm_drv_init      ! Initialization of variables needed from previous timestep
  private :: write_diagnostic  ! Write diagnostic information to log file

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine clm_drv(doalb, nextsw_cday, declinp1, declin, rstwr, nlend, rdate, rof_prognostic)
    !
    ! !DESCRIPTION:
    !
    ! First phase of the clm driver calling the clm physics. An outline of
    ! the calling tree is given in the description of this module.
    !
    ! !USES:
    use clm_time_manager      , only : get_curr_date
    use clm_varctl            , only : use_lai_streams, fates_spitfire_mode
    use clm_varctl            , only : fates_seeddisp_cadence
    use laiStreamMod          , only : lai_advance
    use FATESFireFactoryMod   , only : scalar_lightning
    use FatesInterfaceTypesMod, only : fates_dispersal_cadence_none
    !
    ! !ARGUMENTS:
    implicit none
    logical ,        intent(in) :: doalb       ! true if time for surface albedo calc
    real(r8),        intent(in) :: nextsw_cday ! calendar day for nstep+1
    real(r8),        intent(in) :: declinp1    ! declination angle for next time step
    real(r8),        intent(in) :: declin      ! declination angle for current time step
    logical,         intent(in) :: rstwr       ! true => write restart file this step
    logical,         intent(in) :: nlend       ! true => end of run on this step
    character(len=*),intent(in) :: rdate       ! restart file time stamp for name

    ! Whether we're running with a prognostic ROF component. This shouldn't change from
    ! timestep to timestep, but we pass it into the driver loop because it isn't available
    ! in initialization.
    logical,         intent(in) :: rof_prognostic  ! whether we're running with a prognostic ROF component
    !
    ! !LOCAL VARIABLES:
    integer              :: nstep                   ! time step number
    integer              :: nc, c, p, l, g          ! indices
    integer              :: nclumps                 ! number of clumps on this processor
    integer              :: yr                      ! year (0, ...)
    integer              :: mon                     ! month (1, ..., 12)
    integer              :: day                     ! day of month (1, ..., 31)
    integer              :: sec                     ! seconds of the day
    character(len=256)   :: filer                   ! restart file name
    integer              :: ier                     ! error code
    logical              :: need_glacier_initialization ! true if we need to initialize glacier areas in this time step
    type(bounds_type)    :: bounds_clump
    type(bounds_type)    :: bounds_proc

    ! COMPILER_BUG(wjs, 2016-02-24, pgi 15.10) These temporary allocatable arrays are
    ! needed to work around pgi compiler bugs, as noted below
    real(r8), allocatable :: downreg_patch(:)
    real(r8), allocatable :: leafn_patch(:)
    real(r8), allocatable :: agnpp_patch(:)
    real(r8), allocatable :: bgnpp_patch(:)
    real(r8), allocatable :: annsum_npp_patch(:)
    real(r8), allocatable :: rr_patch(:)
    real(r8), allocatable :: net_carbon_exchange_grc(:)
    real(r8), allocatable :: froot_carbon(:)
    real(r8), allocatable :: croot_carbon(:)

    ! COMPILER_BUG(wjs, 2014-11-29, pgi 14.7) Workaround for internal compiler error with
    ! pgi 14.7 ('normalize_forall_array: non-conformable'), which appears in the call to
    ! CalcIrrigationNeeded. Simply declaring this variable makes the ICE go away.
    real(r8), allocatable :: dummy1_to_make_pgi_happy(:)
    !-----------------------------------------------------------------------

    ! Determine processor bounds and clumps for this processor

    call get_proc_bounds(bounds_proc)
    nclumps = get_proc_clumps()

    ! ========================================================================
    ! In the first time step of a startup or hybrid run, we want to update CLM's glacier
    ! areas to match those given by GLC. This is because, in initialization, we do not yet
    ! know GLC's glacier areas, so CLM's glacier areas are based on the surface dataset
    ! (for a cold start or init_interp run) or the initial conditions file (in a
    ! non-init_interp, non-cold start run) - which may not match GLC's glacier areas for
    ! this configuration. (Coupling fields from GLC aren't received until the run loop.)
    ! Thus, CLM will see a potentially large, fictitious glacier area change in the first
    ! time step. We don't want this fictitious area change to result in any state or flux
    ! adjustments. Thus, we apply this area change here, at the start of the driver loop,
    ! so that in dynSubgrid_driver, it will look like there is no glacier area change in
    ! the first time step. (See
    ! https://github.com/ESCOMP/ctsm/issues/340#issuecomment-410483131 for more
    ! discussion on this.)
    !
    ! This needs to happen very early in the run loop, before any balance checks are
    ! initialized, because - by design - this doesn't conserve mass at the grid cell
    ! level. (The whole point of this code block is that we adjust areas without doing
    ! the typical state or flux adjustments that need to accompany those area changes for
    ! conservation.)
    !
    ! This accomplishes approximately the same effect that we would get if we were able to
    ! update glacier areas in initialization. The one difference - and minor, theoretical
    ! problem - that could arise from this start-of-run-loop update is: If the first time
    ! step of the CESM run loop looked like: (1) GLC runs and updates glacier area (i.e.,
    ! glacier area changes in the first time step compared with what was set in
    ! initialization); (2) coupler passes new glacier area to CLM; (3) CLM runs. Then the
    ! code here would mean that the true change in glacier area between initialization and
    ! the first time step would be ignored as far as state and flux adjustments are
    ! concerned. But this is unlikely to be an issue in practice: Currently GLC doesn't
    ! update this frequently, and even if it did, the change in glacier area in a single
    ! time step would typically be very small.
    !
    ! If we are ever able to change the CESM initialization sequence so that GLC fields
    ! are passed to CLM in initialization, then this code block can be removed.
    ! ========================================================================

    need_glacier_initialization = is_first_step()
    
    if (need_glacier_initialization) then
       !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
       do nc = 1, nclumps
          call get_clump_bounds(nc, bounds_clump)

          call glc2lnd_inst%update_glc2lnd_fracs( &
               bounds = bounds_clump)

          call dynSubgrid_wrapup_weight_changes(bounds_clump, glc_behavior)

       end do
       !$OMP END PARALLEL DO
    end if

    ! ============================================================================
    ! Specified phenology
    ! Done in SP mode, FATES-SP mode and also when dry-deposition is active
    ! ============================================================================
    
    if (use_cn) then
       ! For dry-deposition need to call CLMSP so that mlaidiff is obtained
       ! NOTE: This is also true of FATES below
       if ( n_drydep > 0 ) then
          call t_startf('interpMonthlyVeg')
          call interpMonthlyVeg(bounds_proc, canopystate_inst)
          call t_stopf('interpMonthlyVeg')
       endif

    elseif(use_fates) then

       ! For FATES-Specified phenology mode interpolate the weights for
       ! time-interpolation of monthly vegetation data (as in SP mode below)
       ! Also for FATES with dry-deposition as above need to call CLMSP so that mlaidiff is obtained
       !if ( use_fates_sp .or. (n_drydep > 0 ) ) then    ! Replace with this when we have dry-deposition working
       ! For now don't allow for dry-deposition because of issues in #1044 EBK Jun/17/2022
       if ( use_fates_sp ) then
          call t_startf('interpMonthlyVeg')
          call interpMonthlyVeg(bounds_proc, canopystate_inst)
          call t_stopf('interpMonthlyVeg')
       end if

    else

       ! Determine weights for time interpolation of monthly vegetation data.
       ! This also determines whether it is time to read new monthly vegetation and
       ! obtain updated leaf area index [mlai1,mlai2], stem area index [msai1,msai2],
       ! vegetation top [mhvt1,mhvt2] and vegetation bottom [mhvb1,mhvb2]. The
       ! weights obtained here are used in subroutine SatellitePhenology to obtain time
       ! interpolated values.
       ! This is also done for FATES-SP mode above
       if ( doalb .or. ( n_drydep > 0 ) )then
          call t_startf('interpMonthlyVeg')
          call interpMonthlyVeg(bounds_proc, canopystate_inst)
          call t_stopf('interpMonthlyVeg')
       end if

    end if
    
    ! ==================================================================================
    ! Determine decomp vertical profiles
    !
    ! These routines (alt_calc & decomp_vertprofiles) need to be called before
    ! pftdyn_cnbal, and it appears that they need to be called before pftdyn_interp and
    ! the associated filter updates, too (otherwise we get a carbon balance error)
    ! ==================================================================================

    !$OMP PARALLEL DO PRIVATE (nc,bounds_clump)
    do nc = 1,nclumps
       call get_clump_bounds(nc, bounds_clump)

       ! BUG(wjs, 2014-12-15, bugz 2107) Because of the placement of the following
       ! routines (alt_calc and SoilBiogeochemVerticalProfile) in the driver sequence -
       ! they are called very early in each timestep, before weights are adjusted and
       ! filters are updated - it may be necessary for these routines to compute values
       ! over inactive as well as active points (since some inactive points may soon
       ! become active) - so that's what is done now. Currently, it seems to be okay to do
       ! this, because the variables computed here seem to only depend on quantities that
       ! are valid over inactive as well as active points.

       call t_startf("decomp_vert")
       call active_layer_inst%alt_calc(filter_inactive_and_active(nc)%num_soilc, filter_inactive_and_active(nc)%soilc, &
            temperature_inst)

       ! Filter bgc_soilc operates on all non-sp soil columns
       ! Filter bgc_vegp  operates on all non-fates, non-sp patches (use_cn) on soil
       if ((use_cn .or. use_fates_bgc) .and. decomp_method /= no_soil_decomp) then
          call SoilBiogeochemVerticalProfile(bounds_clump                                       , &
               filter_inactive_and_active(nc)%num_bgc_soilc, filter_inactive_and_active(nc)%bgc_soilc   , &
               filter_inactive_and_active(nc)%num_bgc_vegp, filter_inactive_and_active(nc)%bgc_vegp   , &
               active_layer_inst, soilstate_inst, soilbiogeochem_state_inst)
       end if

       call t_stopf("decomp_vert")
    end do
    !$OMP END PARALLEL DO

    ! ============================================================================
    ! Initialize the mass balance checks for carbon and nitrogen, and zero fluxes for
    ! transient land cover
    ! ============================================================================

    if (use_cn) then
       !$OMP PARALLEL DO PRIVATE (nc,bounds_clump)
       do nc = 1,nclumps
          call get_clump_bounds(nc, bounds_clump)

          call t_startf('cninit')

          call bgc_vegetation_inst%InitEachTimeStep(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc)

          call t_stopf('cninit')
       end do
       !$OMP END PARALLEL DO
    end if

    !$OMP PARALLEL DO PRIVATE (nc,bounds_clump)
    do nc = 1,nclumps
       call get_clump_bounds(nc, bounds_clump)

       call t_startf('begcnbal_grc')
       if (use_cn .or. use_fates_bgc) then
          ! Initialize gridcell-level balance check
          call bgc_vegetation_inst%InitGridcellBalance(bounds_clump, &
               filter(nc)%num_allc, filter(nc)%allc, &
               filter(nc)%num_bgc_soilc, filter(nc)%bgc_soilc, &
               filter(nc)%num_bgc_vegp, filter(nc)%bgc_vegp, &
               soilbiogeochem_carbonstate_inst, &
               c13_soilbiogeochem_carbonstate_inst, &
               c14_soilbiogeochem_carbonstate_inst, &
               soilbiogeochem_nitrogenstate_inst)
       end if
       if (use_lch4) then
          call ch4_init_gridcell_balance_check(bounds_clump, &
               filter(nc)%num_nolakec, filter(nc)%nolakec, &
               filter(nc)%num_lakec, filter(nc)%lakec, &
               ch4_inst)
       end if
       call t_stopf('begcnbal_grc')

       call t_startf('begwbal')
       call WaterGridcellBalance(bounds_clump, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            filter(nc)%num_lakec, filter(nc)%lakec, &
            water_inst, lakestate_inst, &
            use_aquifer_layer = use_aquifer_layer(), flag = 'begwb')
       call t_stopf('begwbal')
    end do
    !$OMP END PARALLEL DO

    ! ============================================================================
    ! Update subgrid weights with dynamic landcover (prescribed transient patches,
    ! CNDV, and or dynamic landunits), and do related adjustments. Note that this
    ! call needs to happen outside loops over nclumps.
    ! ============================================================================

    call t_startf('dyn_subgrid')
    call dynSubgrid_driver(bounds_proc,                                               &
         urbanparams_inst, soilstate_inst, water_inst,                       &
         temperature_inst, energyflux_inst, lakestate_inst, &
         canopystate_inst, photosyns_inst, crop_inst, glc2lnd_inst, bgc_vegetation_inst, &
         soilbiogeochem_state_inst, soilbiogeochem_carbonstate_inst,                  &
         c13_soilbiogeochem_carbonstate_inst, c14_soilbiogeochem_carbonstate_inst,    &
         soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst, &
         soilbiogeochem_carbonflux_inst, ch4_inst, glc_behavior)
    call t_stopf('dyn_subgrid')

    ! ============================================================================
    ! If soil moisture is prescribed from data streams set it here
    ! NOTE: This call needs to happen outside loops over nclumps (as streams are not threadsafe).
    ! ============================================================================
    if (use_soil_moisture_streams) then
       call t_startf('prescribed_sm')
       call PrescribedSoilMoistureAdvance( bounds_proc )
       call t_stopf('prescribed_sm')
    endif
    ! ============================================================================
    ! Initialize the column-level mass balance checks for water, carbon & nitrogen.
    !
    ! For water: Currently, I believe this needs to be done after weights are updated for
    ! prescribed transient patches or CNDV, because column-level water is not generally
    ! conserved when weights change (instead the difference is put in the grid cell-level
    ! terms, qflx_liq_dynbal, etc.). Grid cell-level balance
    ! checks ensure that the grid cell-level water is conserved by considering
    ! qflx_liq_dynbal and calling WaterGridcellBalance
    ! before the weight updates.
    !
    ! For carbon & nitrogen: This needs to be done after dynSubgrid_driver, because the
    ! changes due to dynamic area adjustments can break column-level conservation
    ! ============================================================================

    !$OMP PARALLEL DO PRIVATE (nc,bounds_clump)
    do nc = 1,nclumps
       call get_clump_bounds(nc, bounds_clump)

       if (use_soil_moisture_streams) then
          call t_startf('prescribed_sm')
          call PrescribedSoilMoistureInterp(bounds_clump, soilstate_inst, &
               water_inst%waterstatebulk_inst)
          call t_stopf('prescribed_sm')
       endif
       call t_startf('begwbal')
       call BeginWaterColumnBalance(bounds_clump,             &
            filter(nc)%num_nolakec, filter(nc)%nolakec,       &
            filter(nc)%num_lakec, filter(nc)%lakec,           &
            water_inst, soilhydrology_inst, lakestate_inst, &
            use_aquifer_layer = use_aquifer_layer())

       call t_stopf('begwbal')

       call t_startf('begcnbal_col')
       if (use_cn .or. use_fates_bgc) then
          ! Initialize column-level balance check
          call bgc_vegetation_inst%InitColumnBalance(bounds_clump, &
               filter(nc)%num_allc, filter(nc)%allc, &
               filter(nc)%num_bgc_soilc, filter(nc)%bgc_soilc, &
               filter(nc)%num_bgc_vegp, filter(nc)%bgc_vegp, &
               soilbiogeochem_carbonstate_inst, &
               c13_soilbiogeochem_carbonstate_inst, &
               c14_soilbiogeochem_carbonstate_inst, &
               soilbiogeochem_nitrogenstate_inst)
       end if

       if (use_lch4) then
          call ch4_init_column_balance_check(bounds_clump, &
               filter(nc)%num_nolakec, filter(nc)%nolakec, &
               filter(nc)%num_lakec, filter(nc)%lakec, &
               ch4_inst)
       end if
       call t_stopf('begcnbal_col')

    end do
    !$OMP END PARALLEL DO

    ! ============================================================================
    ! Update dynamic N deposition field, on albedo timestep
    ! currently being done outside clumps loop, but no reason why it couldn't be
    ! re-written to go inside.
    ! ============================================================================

    if (use_cn) then ! .or. use_fates_bgc) then (ndep with fates will be added soon)
       if (.not. ndep_from_cpl) then
          call ndep_interp(bounds_proc, atm2lnd_inst)
       end if
    end if
    
    if(use_cn) then
       call t_startf('bgc_interp')
       call bgc_vegetation_inst%InterpFileInputs(bounds_proc)
       call t_stopf('bgc_interp')
       ! fates_spitfire_mode is assigned an integer value in the namelist
       ! see bld/namelist_files/namelist_definition_clm4_5.xml for details
    else if (fates_spitfire_mode > scalar_lightning) then
       call clm_fates%InterpFileInputs(bounds_proc)
    end if

    ! Get time varying urban data
    call urbantv_inst%urbantv_interp(bounds_proc)

    ! When LAI streams are being used
    ! NOTE: This call needs to happen outside loops over nclumps (as streams are not threadsafe)
    if (doalb .and. use_lai_streams) then
       call lai_advance(bounds_proc)
    endif

    ! When crop calendar streams are being used
    ! NOTE: This call needs to happen outside loops over nclumps (as streams are not threadsafe)
    if (use_crop .and. use_cropcal_streams .and. is_beg_curr_year()) then
      call cropcal_advance( bounds_proc )
    end if

    ! ============================================================================
    ! Initialize variables from previous time step, downscale atm forcings, and
    ! Determine canopy interception and precipitation onto ground surface.
    ! Determine the fraction of foliage covered by water and the fraction
    ! of foliage that is dry and transpiring. Initialize snow layer if the
    ! snow accumulation exceeds 10 mm.
    ! ============================================================================

    !$OMP PARALLEL DO PRIVATE (nc,l,c, bounds_clump, downreg_patch, leafn_patch, agnpp_patch, bgnpp_patch, annsum_npp_patch, rr_patch, froot_carbon, croot_carbon)
    do nc = 1,nclumps
       call get_clump_bounds(nc, bounds_clump)

       call t_startf('drvinit')

       call UpdateDaylength(bounds_clump, declin=declin, obliquity=obliqr)

       ! Initialze variables needed for new driver time step
       call clm_drv_init(bounds_clump, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            filter(nc)%num_nolakep, filter(nc)%nolakep, &
            filter(nc)%num_soilp  , filter(nc)%soilp,   &
            canopystate_inst, water_inst%waterstatebulk_inst, &
            water_inst%waterdiagnosticbulk_inst, &
            energyflux_inst)

       call topo_inst%UpdateTopo(bounds_clump, &
            filter(nc)%num_icec, filter(nc)%icec, &
            glc2lnd_inst, glc_behavior, &
            atm_topo = atm2lnd_inst%forc_topo_grc(bounds_clump%begg:bounds_clump%endg))

       call downscale_forcings(bounds_clump, &
            topo_inst, atm2lnd_inst, surfalb_inst, water_inst%wateratm2lndbulk_inst, &
            eflx_sh_precip_conversion = energyflux_inst%eflx_sh_precip_conversion_col(bounds_clump%begc:bounds_clump%endc))

       call set_atm2lnd_water_tracers(bounds_clump, &
            filter(nc)%num_allc, filter(nc)%allc, &
            water_inst)

       if (water_inst%DoConsistencyCheck()) then
          call t_startf("tracer_consistency_check")
          call water_inst%TracerConsistencyCheck(bounds_clump, 'after downscale_forcings')
          call t_stopf("tracer_consistency_check")
       end if

       ! Update filters that depend on variables set in clm_drv_init

       call setExposedvegpFilter(bounds_clump, &
            canopystate_inst%frac_veg_nosno_patch(bounds_clump%begp:bounds_clump%endp))

       call t_stopf('drvinit')

       if (irrigate) then

          call t_startf('irrigationwithdraw')

          call CalcAndWithdrawIrrigationFluxes( &
               bounds = bounds_clump, &
               num_soilc = filter(nc)%num_soilc, &
               filter_soilc = filter(nc)%soilc, &
               num_soilp = filter(nc)%num_soilp, &
               filter_soilp = filter(nc)%soilp, &
               soilhydrology_inst = soilhydrology_inst, &
               soilstate_inst = soilstate_inst, &
               irrigation_inst = irrigation_inst, &
               water_inst = water_inst)

          call t_stopf('irrigationwithdraw')

       end if

       if (water_inst%DoConsistencyCheck()) then
          call t_startf("tracer_consistency_check")
          call water_inst%TracerConsistencyCheck(bounds_clump, 'after CalcAndWithdrawIrrigationFluxes')
          call t_stopf("tracer_consistency_check")
       end if

       ! ============================================================================
       ! First Stage of Hydrology
       ! (1) water storage of intercepted precipitation
       ! (2) direct throughfall and canopy drainage of precipitation
       ! (3) fraction of foliage covered by water and the fraction is dry and transpiring
       ! (4) snow layer initialization if the snow accumulation exceeds 10 mm.
       ! ============================================================================

       call t_startf('hydro1')

       call CanopyInterceptionAndThroughfall(bounds_clump, &
            filter(nc)%num_soilp, filter(nc)%soilp, &
            filter(nc)%num_nolakep, filter(nc)%nolakep, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            patch, col, canopystate_inst, atm2lnd_inst, water_inst)

       call HandleNewSnow(bounds_clump, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            scf_method, &
            atm2lnd_inst, temperature_inst, &
            aerosol_inst, water_inst)

       ! update surface water fraction (this may modify frac_sno)
       call UpdateFracH2oSfc(bounds_clump, &
            filter(nc)%num_soilc, filter(nc)%soilc, &
            water_inst)

       if (water_inst%DoConsistencyCheck()) then
          call t_startf("tracer_consistency_check")
          call water_inst%TracerConsistencyCheck(bounds_clump, 'after first stage of hydrology')
          call t_stopf("tracer_consistency_check")
       end if

       call t_stopf('hydro1')

       ! ============================================================================
       ! Surface Radiation
       ! ============================================================================

       call t_startf('surfrad')

       ! Surface Radiation primarily for non-urban columns

       ! Most of the surface radiation calculations are agnostic to the forest-model
       ! but the calculations of the fractions of sunlit and shaded canopies
       ! are specific, calculate them first.
       ! The nourbanp filter is set in dySubgrid_driver (earlier in this call)
       ! over the patch index range defined by bounds_clump%begp:bounds_proc%endp

       if(use_fates) then
          call clm_fates%wrap_sunfrac(nc,atm2lnd_inst, canopystate_inst)
       else
          call CanopySunShadeFracs(filter(nc)%nourbanp,filter(nc)%num_nourbanp,     &
                                   atm2lnd_inst, surfalb_inst, canopystate_inst,    &
                                   solarabs_inst)
       end if



       call SurfaceRadiation(bounds_clump,                                 &
            filter(nc)%num_nourbanp, filter(nc)%nourbanp,                  &
            filter(nc)%num_urbanp, filter(nc)%urbanp,                      &
            filter(nc)%num_urbanc, filter(nc)%urbanc,                      &
            atm2lnd_inst, water_inst%waterdiagnosticbulk_inst, &
            canopystate_inst, surfalb_inst, &
            solarabs_inst, surfrad_inst)

       ! Surface Radiation for only urban columns

       call UrbanRadiation(bounds_clump,                                       &
            filter(nc)%num_nourbanl, filter(nc)%nourbanl,                      &
            filter(nc)%num_urbanl, filter(nc)%urbanl,                          &
            filter(nc)%num_urbanc, filter(nc)%urbanc,                          &
            filter(nc)%num_urbanp, filter(nc)%urbanp,                          &
            atm2lnd_inst, water_inst%waterdiagnosticbulk_inst, &
            temperature_inst, urbanparams_inst, &
            solarabs_inst, surfalb_inst, energyflux_inst)

       call t_stopf('surfrad')

       ! ============================================================================
       ! Determine leaf temperature and surface fluxes based on ground
       ! temperature from previous time step.
       ! ============================================================================

       call t_startf('bgp1')

       call BiogeophysPreFluxCalcs(bounds_clump,                                   &
            filter(nc)%num_nolakec, filter(nc)%nolakec,                       &
            filter(nc)%num_nolakep, filter(nc)%nolakep,                       &
            filter(nc)%num_urbanc, filter(nc)%urbanc,                         &
            clm_fates,                                                        &
            atm2lnd_inst, canopystate_inst, energyflux_inst, frictionvel_inst, &
            soilstate_inst, temperature_inst, &
            water_inst%wateratm2lndbulk_inst, water_inst%waterdiagnosticbulk_inst, &
            water_inst%waterstatebulk_inst, water_inst%waterfluxbulk_inst)

       call ozone_inst%CalcOzoneStress(bounds_clump, &
            filter(nc)%num_exposedvegp, filter(nc)%exposedvegp, &
            filter(nc)%num_noexposedvegp, filter(nc)%noexposedvegp)

       ! TODO(wjs, 2019-10-02) I'd like to keep moving this down until it is below
       ! LakeFluxes... I'll probably leave it in place there.
       if (water_inst%DoConsistencyCheck()) then
          call t_startf("tracer_consistency_check")
          call water_inst%TracerConsistencyCheck(bounds_clump, 'after BiogeophysPreFluxCalcs')
          call t_stopf("tracer_consistency_check")
       end if

       call CalculateSurfaceHumidity(bounds_clump,                                   &
            filter(nc)%num_nolakec, filter(nc)%nolakec,                       &
            atm2lnd_inst, temperature_inst, &
            water_inst%waterstatebulk_inst, water_inst%wateratm2lndbulk_inst, &
            soilstate_inst, water_inst%waterdiagnosticbulk_inst)
       call t_stopf('bgp1')

       ! ============================================================================
       ! Determine fluxes
       ! ============================================================================

       call t_startf('bgp_fluxes')
       call t_startf('bgflux')

       ! Bareground fluxes for all patches except lakes and urban landunits

       call BareGroundFluxes(bounds_clump,                                 &
            filter(nc)%num_noexposedvegp, filter(nc)%noexposedvegp,          &
            atm2lnd_inst, soilstate_inst,                &
            frictionvel_inst, ch4_inst, energyflux_inst, temperature_inst, &
            water_inst%waterfluxbulk_inst, water_inst%waterstatebulk_inst, &
            water_inst%waterdiagnosticbulk_inst, water_inst%wateratm2lndbulk_inst, &
            photosyns_inst, humanindex_inst, canopystate_inst)
       call t_stopf('bgflux')

       ! non-bareground fluxes for all patches except lakes and urban landunits
       ! Calculate canopy temperature, latent and sensible fluxes from the canopy,
       ! and leaf water change by evapotranspiration

       call t_startf('canflux')

       ! COMPILER_BUG(wjs, 2016-02-24, pgi 15.10) In principle, we should be able to make
       ! these function calls inline in the CanopyFluxes argument list. However, with pgi
       ! 15.10, that results in the dummy arguments having the wrong size (I suspect size
       ! 0, based on similar pgi compiler bugs that we have run into before). Also note
       ! that I don't have explicit bounds on the left-hand-side of these assignments:
       ! excluding these explicit bounds seemed to be needed to get around other compiler
       ! bugs.
       allocate(downreg_patch(bounds_clump%begp:bounds_clump%endp))
       allocate(leafn_patch(bounds_clump%begp:bounds_clump%endp))
       downreg_patch = bgc_vegetation_inst%get_downreg_patch(bounds_clump)
       leafn_patch = bgc_vegetation_inst%get_leafn_patch(bounds_clump)

       allocate(froot_carbon(bounds_clump%begp:bounds_clump%endp))
       allocate(croot_carbon(bounds_clump%begp:bounds_clump%endp))
       froot_carbon = bgc_vegetation_inst%get_froot_carbon_patch( &
            bounds_clump, canopystate_inst%tlai_patch(bounds_clump%begp:bounds_clump%endp))
       croot_carbon = bgc_vegetation_inst%get_croot_carbon_patch( &
            bounds_clump, canopystate_inst%tlai_patch(bounds_clump%begp:bounds_clump%endp))

       call CanopyFluxes(bounds_clump,                                                      &
            filter(nc)%num_exposedvegp, filter(nc)%exposedvegp,                             &
            clm_fates,nc,                                                                   &
            active_layer_inst, atm2lnd_inst, canopystate_inst,                              &
            energyflux_inst, frictionvel_inst, soilstate_inst, solarabs_inst, surfalb_inst, &
            temperature_inst, water_inst%waterfluxbulk_inst, water_inst%waterstatebulk_inst, &
            water_inst%waterdiagnosticbulk_inst, water_inst%wateratm2lndbulk_inst,          &
            ch4_inst, ozone_inst, photosyns_inst, &
            humanindex_inst, soil_water_retention_curve, &
            downreg_patch = downreg_patch(bounds_clump%begp:bounds_clump%endp), &
            leafn_patch = leafn_patch(bounds_clump%begp:bounds_clump%endp), &
            froot_carbon = froot_carbon(bounds_clump%begp:bounds_clump%endp), &
            croot_carbon = croot_carbon(bounds_clump%begp:bounds_clump%endp))
       deallocate(downreg_patch, leafn_patch, froot_carbon, croot_carbon)
       call t_stopf('canflux')

       ! Fluxes for all urban landunits

       call t_startf('uflux')
       call UrbanFluxes(bounds_clump,                                         &
            filter(nc)%num_nourbanl, filter(nc)%nourbanl,                     &
            filter(nc)%num_urbanl, filter(nc)%urbanl,                         &
            filter(nc)%num_urbanc, filter(nc)%urbanc,                         &
            filter(nc)%num_urbanp, filter(nc)%urbanp,                         &
            atm2lnd_inst, urbanparams_inst, soilstate_inst, temperature_inst,   &
            water_inst%waterstatebulk_inst, water_inst%waterdiagnosticbulk_inst, &
            frictionvel_inst, energyflux_inst, water_inst%waterfluxbulk_inst, &
            water_inst%wateratm2lndbulk_inst, humanindex_inst)
       call t_stopf('uflux')

       ! Fluxes for all lake landunits

       call t_startf('bgplake')
       call LakeFluxes(bounds_clump,                                         &
            filter(nc)%num_lakec, filter(nc)%lakec,                          &
            filter(nc)%num_lakep, filter(nc)%lakep,                          &
            atm2lnd_inst, solarabs_inst, frictionvel_inst, temperature_inst, &
            energyflux_inst, water_inst%waterstatebulk_inst,                 &
            water_inst%waterdiagnosticbulk_inst,                             &
            water_inst%waterfluxbulk_inst, water_inst%wateratm2lndbulk_inst, &
            lakestate_inst,&
            humanindex_inst)
       call t_stopf('bgplake')

       call frictionvel_inst%SetActualRoughnessLengths( &
            bounds = bounds_clump, &
            num_exposedvegp = filter(nc)%num_exposedvegp, &
            filter_exposedvegp = filter(nc)%exposedvegp, &
            num_noexposedvegp = filter(nc)%num_noexposedvegp, &
            filter_noexposedvegp = filter(nc)%noexposedvegp, &
            num_urbanp = filter(nc)%num_urbanp, &
            filter_urbanp = filter(nc)%urbanp, &
            num_lakep = filter(nc)%num_lakep, &
            filter_lakep = filter(nc)%lakep)

       call t_stopf('bgp_fluxes')

       if (irrigate) then

          ! ============================================================================
          ! Determine irrigation needed for future time steps
          ! ============================================================================

          ! NOTE(wjs, 2016-09-08) The placement of this call in the driver is historical: it
          ! used to be that it had to come after btran was computed. Now it no longer depends
          ! on btran, so it could be moved earlier in the driver loop - possibly even
          ! immediately before ApplyIrrigation, which would be a more clear place to put it.

          call t_startf('irrigationneeded')
          call irrigation_inst%CalcIrrigationNeeded( &
               bounds             = bounds_clump, &
               num_exposedvegp    = filter(nc)%num_exposedvegp, &
               filter_exposedvegp = filter(nc)%exposedvegp, &
               elai               = canopystate_inst%elai_patch(bounds_clump%begp:bounds_clump%endp), &
               t_soisno           = temperature_inst%t_soisno_col(bounds_clump%begc:bounds_clump%endc  , 1:nlevgrnd), &
               eff_porosity       = soilstate_inst%eff_porosity_col(bounds_clump%begc:bounds_clump%endc, 1:nlevgrnd), &
               h2osoi_liq         = water_inst%waterstatebulk_inst%h2osoi_liq_col&
               (bounds_clump%begc:bounds_clump%endc , 1:nlevgrnd), &
               volr               = water_inst%wateratm2lndbulk_inst%volrmch_grc(bounds_clump%begg:bounds_clump%endg), &
               rof_prognostic     = rof_prognostic)
          call t_stopf('irrigationneeded')

       end if

       ! ============================================================================
       ! DUST and VOC emissions
       ! ============================================================================

       call t_startf('bgc')

       ! Dust mobilization (C. Zender's modified codes)
       call dust_emis_inst%DustEmission(bounds_clump,                                       &
            filter(nc)%num_nolakep, filter(nc)%nolakep,                      &
            atm2lnd_inst, soilstate_inst, canopystate_inst, &
            water_inst%waterstatebulk_inst, water_inst%waterdiagnosticbulk_inst, &
            frictionvel_inst)

       ! Dust dry deposition (C. Zender's modified codes)
       call dust_emis_inst%DustDryDep(bounds_clump, &
            atm2lnd_inst, frictionvel_inst)

       ! VOC emission (A. Guenther's MEGAN (2006) model; Wang et al., 2022, 2024a, 2024b)
       call VOCEmission(bounds_clump,                                         &
               filter(nc)%num_soilp, filter(nc)%soilp,                           &
               atm2lnd_inst, canopystate_inst, photosyns_inst, temperature_inst, &
               vocemis_inst, energyflux_inst)

       call t_stopf('bgc')

       ! ============================================================================
       ! Determine temperatures
       ! ============================================================================

       ! Set lake temperature

       call t_startf('lakeTemp')
       call LakeTemperature(bounds_clump,                                             &
            filter(nc)%num_lakec, filter(nc)%lakec,                                   &
            filter(nc)%num_lakep, filter(nc)%lakep,                                   &
            solarabs_inst, soilstate_inst, water_inst%waterstatebulk_inst, &
            water_inst%waterdiagnosticbulk_inst, water_inst%waterfluxbulk_inst, ch4_inst, &
            energyflux_inst, temperature_inst, lakestate_inst)
       call t_stopf('lakeTemp')

       ! Set soil/snow temperatures including ground temperature

       call t_startf('soiltemperature')
       call SoilTemperature(bounds_clump,                                                      &
            filter(nc)%num_urbanl  , filter(nc)%urbanl,                                        &
            filter(nc)%num_urbanc  , filter(nc)%urbanc,                                        &
            filter(nc)%num_nolakep , filter(nc)%nolakep,                                       &
            filter(nc)%num_nolakec , filter(nc)%nolakec,                                       &
            atm2lnd_inst, urbanparams_inst, canopystate_inst, water_inst%waterstatebulk_inst, &
            water_inst%waterdiagnosticbulk_inst, water_inst%waterfluxbulk_inst, &
            solarabs_inst, soilstate_inst, energyflux_inst,  temperature_inst, urbantv_inst)

       ! The following is called immediately after SoilTemperature so that melted ice is
       ! converted back to solid ice as soon as possible
       call glacier_smb_inst%HandleIceMelt(bounds_clump, &
            filter(nc)%num_do_smb_c, filter(nc)%do_smb_c, &
            water_inst%waterstatebulk_inst, water_inst%waterfluxbulk_inst)
       call t_stopf('soiltemperature')

       ! ============================================================================
       ! update surface fluxes for new ground temperature.
       ! ============================================================================

       call t_startf('bgp2')
       call SoilFluxes(bounds_clump,                                                          &
            filter(nc)%num_urbanl,  filter(nc)%urbanl,                                        &
            filter(nc)%num_urbanp,  filter(nc)%urbanp,                                        &
            filter(nc)%num_nolakec, filter(nc)%nolakec,                                       &
            filter(nc)%num_nolakep, filter(nc)%nolakep,                                       &
            atm2lnd_inst, solarabs_inst, temperature_inst, canopystate_inst, &
            water_inst%waterstatebulk_inst, water_inst%waterdiagnosticbulk_inst, &
            energyflux_inst, water_inst%waterfluxbulk_inst)
       call t_stopf('bgp2')

       ! ============================================================================
       ! Perform averaging from patch level to column level
       ! ============================================================================

       call t_startf('patch2col')
       call clm_drv_patch2col(bounds_clump, &
            filter(nc)%num_allc, filter(nc)%allc, filter(nc)%num_nolakec, filter(nc)%nolakec, &
            energyflux_inst, water_inst%waterfluxbulk_inst)
       call t_stopf('patch2col')

       ! ============================================================================
       ! Vertical (column) soil and surface hydrology
       ! ============================================================================

       ! Note that filter_snowc and filter_nosnowc are returned by
       ! LakeHydrology after the new snow filter is built

       call t_startf('hydro_without_drainage')

       call HydrologyNoDrainage(bounds_clump,                                &
            filter(nc)%num_nolakec, filter(nc)%nolakec,                      &
            filter(nc)%num_hydrologyc, filter(nc)%hydrologyc,                &
            filter(nc)%num_urbanc, filter(nc)%urbanc,                        &
            filter(nc)%num_snowc, filter(nc)%snowc,                          &
            filter(nc)%num_nosnowc, filter(nc)%nosnowc,                      &
            clm_fates,                                                         &
            atm2lnd_inst, soilstate_inst, energyflux_inst, temperature_inst,   &
            water_inst, soilhydrology_inst, &
            saturated_excess_runoff_inst, &
            infiltration_excess_runoff_inst, &
            aerosol_inst, canopystate_inst, scf_method, soil_water_retention_curve, topo_inst)

       ! The following needs to be done after HydrologyNoDrainage (because it needs
       ! waterfluxbulk_inst%qflx_snwcp_ice_col), but before HydrologyDrainage (because
       ! HydrologyDrainage calls glacier_smb_inst%AdjustRunoffTerms, which depends on
       ! ComputeSurfaceMassBalance having already been called).
       call glacier_smb_inst%ComputeSurfaceMassBalance(bounds_clump, &
            filter(nc)%num_allc, filter(nc)%allc, &
            filter(nc)%num_do_smb_c, filter(nc)%do_smb_c, &
            glc2lnd_inst, water_inst%waterstatebulk_inst, water_inst%waterfluxbulk_inst)

       !  Calculate column-integrated aerosol masses, and
       !  mass concentrations for radiative calculations and output
       !  (based on new snow level state, after SnowFilter is rebuilt.
       !  NEEDS TO BE AFTER SnowFiler is rebuilt, otherwise there
       !  can be zero snow layers but an active column in filter)

       call AerosolMasses( bounds_clump,                                   &
            num_on=filter(nc)%num_snowc, filter_on=filter(nc)%snowc,       &
            num_off=filter(nc)%num_nosnowc, filter_off=filter(nc)%nosnowc, &
            waterfluxbulk_inst = water_inst%waterfluxbulk_inst,            &
            waterstatebulk_inst = water_inst%waterstatebulk_inst,          &
            waterdiagnosticbulk_inst = water_inst%waterdiagnosticbulk_inst, &
            aerosol_inst=aerosol_inst)

       call t_stopf('hydro_without_drainage')

       ! ============================================================================
       ! Lake hydrology
       ! ============================================================================

       ! Note that filter_lakesnowc and filter_lakenosnowc are returned by
       ! LakeHydrology after the new snow filter is built

       call t_startf('hylake')
       call LakeHydrology(bounds_clump,                                                      &
            filter(nc)%num_lakec, filter(nc)%lakec,                                          &
            filter(nc)%num_lakep, filter(nc)%lakep,                                          &
            filter(nc)%num_lakesnowc, filter(nc)%lakesnowc,                                  &
            filter(nc)%num_lakenosnowc, filter(nc)%lakenosnowc,                              &
            scf_method, water_inst, &
            atm2lnd_inst, temperature_inst, soilstate_inst, &
            energyflux_inst, aerosol_inst, lakestate_inst, topo_inst)

       !  Calculate column-integrated aerosol masses, and
       !  mass concentrations for radiative calculations and output
       !  (based on new snow level state, after SnowFilter is rebuilt.
       !  NEEDS TO BE AFTER SnowFiler is rebuilt, otherwise there
       !  can be zero snow layers but an active column in filter)

       call AerosolMasses(bounds_clump,                                            &
            num_on=filter(nc)%num_lakesnowc, filter_on=filter(nc)%lakesnowc,       &
            num_off=filter(nc)%num_lakenosnowc, filter_off=filter(nc)%lakenosnowc, &
            waterfluxbulk_inst = water_inst%waterfluxbulk_inst,                    &
            waterstatebulk_inst = water_inst%waterstatebulk_inst,                  &
            waterdiagnosticbulk_inst = water_inst%waterdiagnosticbulk_inst,        &
            aerosol_inst=aerosol_inst)

       ! Must be done here because must use a snow filter for lake columns

       call SnowAge_grain(bounds_clump,                         &
            filter(nc)%num_lakesnowc, filter(nc)%lakesnowc,     &
            filter(nc)%num_lakenosnowc, filter(nc)%lakenosnowc, &
            water_inst%waterfluxbulk_inst, water_inst%waterstatebulk_inst, &
            water_inst%waterdiagnosticbulk_inst, temperature_inst, &
            atm2lnd_inst)

       call t_stopf('hylake')

       ! ============================================================================
       ! ! Fraction of soil covered by snow (Z.-L. Yang U. Texas)
       ! ============================================================================
       call t_startf('snow_init')

       do c = bounds_clump%begc,bounds_clump%endc
          l = col%landunit(c)
          if (lun%urbpoi(l)) then
             ! Urban landunit use Bonan 1996 (LSM Technical Note)
             water_inst%waterdiagnosticbulk_inst%frac_sno_col(c) = &
                  min( water_inst%waterdiagnosticbulk_inst%snow_depth_col(c)/0.05_r8, 1._r8)
          end if
       end do

       ! ============================================================================
       ! Snow aging routine based on Flanner and Zender (2006), Linking snowpack
       ! microphysics and albedo evolution, JGR, and Brun (1989), Investigation of
       ! wet-snow metamorphism in respect of liquid-water content, Ann. Glaciol.
       ! ============================================================================
       ! Note the snow filters here do not include lakes
       ! TODO: move this up

       call SnowAge_grain(bounds_clump,                 &
            filter(nc)%num_snowc, filter(nc)%snowc,     &
            filter(nc)%num_nosnowc, filter(nc)%nosnowc, &
            water_inst%waterfluxbulk_inst, water_inst%waterstatebulk_inst, &
            water_inst%waterdiagnosticbulk_inst, temperature_inst, &
            atm2lnd_inst)
       call t_stopf('snow_init')

       ! ============================================================================
       ! Ecosystem dynamics: Uses CN, CNDV, or static parameterizations
       ! ============================================================================

       ! FIX(SPM,032414)  push these checks into the routines below and/or make this consistent.

             ! fully prognostic canopy structure and C-N biogeochemistry
             ! - CNDV defined: prognostic biogeography; else prescribed
       ! - crop model:  crop algorithms called from within CNDriver

       ! Filter bgc_soilc operates on all non-sp soil columns
       ! Filter bgc_vegp  operates on all non-fates, non-sp patches (use_cn) on soil

       if(use_cn .or. use_fates_bgc)then
          call t_startf('ecosysdyn')
          call bgc_vegetation_inst%EcosystemDynamicsPreDrainage(bounds_clump,            &
               filter(nc)%num_bgc_soilc, filter(nc)%bgc_soilc,                       &
               filter(nc)%num_bgc_vegp, filter(nc)%bgc_vegp,                       &
               filter(nc)%num_actfirec, filter(nc)%actfirec,                 &
               filter(nc)%num_actfirep, filter(nc)%actfirep,                 &
               filter(nc)%num_pcropp, filter(nc)%pcropp, &
               filter(nc)%num_soilnopcropp, filter(nc)%soilnopcropp, &
               filter(nc)%num_exposedvegp, filter(nc)%exposedvegp, &
               filter(nc)%num_noexposedvegp, filter(nc)%noexposedvegp, &
               soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst,         &
               c13_soilbiogeochem_carbonflux_inst, c13_soilbiogeochem_carbonstate_inst, &
               c14_soilbiogeochem_carbonflux_inst, c14_soilbiogeochem_carbonstate_inst, &
               soilbiogeochem_state_inst,                                               &
               soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst,     &
               active_layer_inst, clm_fates, &
               atm2lnd_inst, water_inst%waterstatebulk_inst, &
               water_inst%waterdiagnosticbulk_inst, water_inst%waterfluxbulk_inst,      &
               water_inst%wateratm2lndbulk_inst, canopystate_inst, soilstate_inst, temperature_inst, &
               soil_water_retention_curve, crop_inst, ch4_inst, &
               photosyns_inst, saturated_excess_runoff_inst, energyflux_inst,          &
               nutrient_competition_method, fireemis_inst)
          call t_stopf('ecosysdyn')
       end if

       ! Prescribed biogeography - prescribed canopy structure, some prognostic carbon fluxes

       if (((.not. use_cn) .and. (.not. use_fates) .and. (doalb))) then
          call t_startf('SatellitePhenology')
          call CalcSatellitePhenologyTimeInterp(bounds_clump, filter(nc)%num_nolakep, filter(nc)%nolakep, &
               canopystate_inst)
          call UpdateSatellitePhenologyCanopy(bounds_clump, filter(nc)%num_nolakep, filter(nc)%nolakep, &
               water_inst%waterdiagnosticbulk_inst, canopystate_inst)
          call t_stopf('SatellitePhenology')
       end if

       if (use_fates_sp.and.doalb) then
          call t_startf('SatellitePhenology')

          ! FATES satellite phenology mode needs to include all active and inactive patch-level soil
          ! filters due to the translation between the hlm pfts and the fates pfts.
          ! E.g. in FATES, an active PFT vector of 1, 0, 0, 0, 1, 0, 1, 0 would be mapped into
          ! the host land model as 1, 1, 1, 0, 0, 0, 0.  As such, the 'active' filter would only
          ! use the first three points, which would incorrectly represent the interpolated values.
          call CalcSatellitePhenologyTimeInterp(bounds_clump, &
               filter_inactive_and_active(nc)%num_soilp, filter_inactive_and_active(nc)%soilp, &
               canopystate_inst)
          call t_stopf('SatellitePhenology')

       end if

       ! Dry Deposition of chemical tracers (Wesely (1998) parameterizaion)

       call t_startf('depvel')
       call depvel_compute(bounds_clump, &
            atm2lnd_inst, canopystate_inst, water_inst%waterstatebulk_inst, &
            water_inst%waterdiagnosticbulk_inst, water_inst%wateratm2lndbulk_inst, &
            frictionvel_inst, photosyns_inst, drydepvel_inst)
       call t_stopf('depvel')

       if (use_crop .and. use_cropcal_streams .and. is_beg_curr_year()) then
          ! ============================================================================
          ! Update crop calendars
          ! ============================================================================
          call cropcal_interp(bounds_clump, filter_inactive_and_active(nc)%num_pcropp, &
               filter_inactive_and_active(nc)%pcropp, .false., crop_inst)
       end if

       ! ============================================================================
       ! Calculate soil/snow hydrology with drainage (subsurface runoff)
       ! ============================================================================

       call t_startf('hydro2_drainage')

       call HydrologyDrainage(bounds_clump,                   &
            filter(nc)%num_nolakec, filter(nc)%nolakec,       &
            filter(nc)%num_hydrologyc, filter(nc)%hydrologyc, &
            filter(nc)%num_urbanc, filter(nc)%urbanc,         &
            filter(nc)%num_do_smb_c, filter(nc)%do_smb_c,     &
            glc2lnd_inst, temperature_inst,                   &
            soilhydrology_inst, soilstate_inst, water_inst%waterstatebulk_inst, &
            water_inst%waterdiagnosticbulk_inst, water_inst%waterbalancebulk_inst, &
            water_inst%waterfluxbulk_inst, water_inst%wateratm2lndbulk_inst, &
            glacier_smb_inst)

       call t_stopf('hydro2_drainage')

       
       if (use_cn .or. use_fates_bgc) then
          call t_startf('EcosysDynPostDrainage')
          call bgc_vegetation_inst%EcosystemDynamicsPostDrainage(bounds_clump, &
               filter(nc)%num_allc, filter(nc)%allc, &
               filter(nc)%num_bgc_soilc, filter(nc)%bgc_soilc, &
               filter(nc)%num_bgc_vegp, filter(nc)%bgc_vegp, &
               filter(nc)%num_actfirec, filter(nc)%actfirec,                 &
               filter(nc)%num_actfirep, filter(nc)%actfirep,                 &
               doalb, crop_inst, &
               soilstate_inst, soilbiogeochem_state_inst, &
               water_inst%waterstatebulk_inst, water_inst%waterdiagnosticbulk_inst, &
               water_inst%waterfluxbulk_inst, frictionvel_inst, canopystate_inst, &
               soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst, &
               c13_soilbiogeochem_carbonflux_inst, c13_soilbiogeochem_carbonstate_inst, &
               c14_soilbiogeochem_carbonflux_inst, c14_soilbiogeochem_carbonstate_inst, &
               soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst)
          call t_stopf('EcosysDynPostDrainage')
       end if

       if ( use_fates) then

          ! FATES has its own running mean functions, such as 24hr
          ! vegetation temperature and exponential moving averages
          ! for leaf photosynthetic acclimation temperature. These
          ! moving averages are updated here
          call clm_fates%WrapUpdateFatesRmean(nc,temperature_inst)


          call clm_fates%wrap_update_hifrq_hist(bounds_clump, &
               soilbiogeochem_carbonflux_inst, &
               soilbiogeochem_carbonstate_inst)


          if( is_beg_curr_day() ) then

             ! --------------------------------------------------------------------------
             ! This is the main call to FATES dynamics
             ! --------------------------------------------------------------------------

             call clm_fates%dynamics_driv( nc, bounds_clump,                        &
                  atm2lnd_inst, soilstate_inst, temperature_inst, active_layer_inst, &
                  water_inst%waterstatebulk_inst, water_inst%waterdiagnosticbulk_inst, &
                  water_inst%wateratm2lndbulk_inst, canopystate_inst, soilbiogeochem_carbonflux_inst, &
                  frictionvel_inst, soil_water_retention_curve)

             ! TODO(wjs, 2016-04-01) I think this setFilters call should be replaced by a
             ! call to reweight_wrapup, if it's needed at all.
             call setFilters( bounds_clump, glc_behavior )

          end if



       end if ! use_fates branch

       ! ============================================================================
       ! Create summaries of water diagnostic terms
       ! ============================================================================

       call water_inst%Summary(bounds_clump, &
            filter(nc)%num_soilp, filter(nc)%soilp, &
            filter(nc)%num_allc, filter(nc)%allc, &
            filter(nc)%num_nolakec, filter(nc)%nolakec)

       ! ============================================================================
       ! Check the carbon and nitrogen balance
       ! ============================================================================

       if(use_cn .or. use_fates_bgc)then
          call t_startf('cnbalchk')
          call bgc_vegetation_inst%BalanceCheck( &
               bounds_clump, filter(nc)%num_bgc_soilc, filter(nc)%bgc_soilc, &
               soilbiogeochem_carbonflux_inst, & 
               soilbiogeochem_nitrogenflux_inst, &
               soilbiogeochem_carbonstate_inst, & 
               soilbiogeochem_nitrogenstate_inst, &
               atm2lnd_inst, clm_fates )
          call t_stopf('cnbalchk')
       end if
       
       ! Calculation of methane fluxes

       if (use_lch4) then
          call t_startf('ch4')
          ! COMPILER_BUG(wjs, 2016-02-24, pgi 15.10) In principle, we should be able to
          ! make these function calls inline in the CanopyFluxes argument list. However,
          ! with pgi 15.10, that results in the dummy arguments having the wrong size (I
          ! suspect size 0, based on similar pgi compiler bugs that we have run into
          ! before). Also note that I don't have explicit bounds on the left-hand-side of
          ! these assignments: excluding these explicit bounds seemed to be needed to get
          ! around other compiler bugs.
          allocate(agnpp_patch(bounds_clump%begp:bounds_clump%endp))
          allocate(bgnpp_patch(bounds_clump%begp:bounds_clump%endp))
          allocate(annsum_npp_patch(bounds_clump%begp:bounds_clump%endp))
          allocate(rr_patch(bounds_clump%begp:bounds_clump%endp))
          agnpp_patch = bgc_vegetation_inst%get_agnpp_patch(bounds_clump)
          bgnpp_patch = bgc_vegetation_inst%get_bgnpp_patch(bounds_clump)
          annsum_npp_patch = bgc_vegetation_inst%get_annsum_npp_patch(bounds_clump)
          rr_patch = bgc_vegetation_inst%get_root_respiration_patch(bounds_clump)

          call ch4 (bounds_clump,                                                                  &
               filter(nc)%num_soilc, filter(nc)%soilc,                                             &
               filter(nc)%num_lakec, filter(nc)%lakec,                                             &
               filter(nc)%num_nolakec, filter(nc)%nolakec,                                         &
               filter(nc)%num_soilp, filter(nc)%soilp,                                             &
               atm2lnd_inst, lakestate_inst, canopystate_inst, soilstate_inst, soilhydrology_inst, &
               temperature_inst, energyflux_inst, water_inst%waterstatebulk_inst, &
               water_inst%waterdiagnosticbulk_inst, water_inst%waterfluxbulk_inst,                 &
               soilbiogeochem_carbonflux_inst,                          &
               soilbiogeochem_nitrogenflux_inst, ch4_inst, lnd2atm_inst, &
               clm_fates, &
               agnpp = agnpp_patch(bounds_clump%begp:bounds_clump%endp), &
               bgnpp = bgnpp_patch(bounds_clump%begp:bounds_clump%endp), &
               annsum_npp = annsum_npp_patch(bounds_clump%begp:bounds_clump%endp), &
               rr = rr_patch(bounds_clump%begp:bounds_clump%endp))
          deallocate(agnpp_patch, bgnpp_patch, annsum_npp_patch, rr_patch)
          call t_stopf('ch4')
       end if

       ! ============================================================================
       ! Determine albedos for next time step
       ! ============================================================================
       
       if (.not.doalb) then


          if(use_fates) then
             ! During branch runs and non continue_run restarts, the doalb flag
             ! does not trigger correctly for fates runs (and non-fates?), and thus
             ! the zenith angles are not calculated and ready when radiation scattering
             ! needs to occur.
             call UpdateZenithAngles(bounds_clump, surfalb_inst, nextsw_cday, declinp1)
             call clm_fates%wrap_canopy_radiation(bounds_clump, nc, &
                  water_inst%waterdiagnosticbulk_inst%fcansno_patch(bounds_clump%begp:bounds_clump%endp), &
                  surfalb_inst)
          end if
          
       else

          ! Albedos for non-urban columns
          call t_startf('surfalb')
          call SurfaceAlbedo(bounds_clump,                      &
               nc,                                              &
               filter_inactive_and_active(nc)%num_nourbanc,     &
               filter_inactive_and_active(nc)%nourbanc,         &
               filter_inactive_and_active(nc)%num_nourbanp,     &
               filter_inactive_and_active(nc)%nourbanp,         &
               filter_inactive_and_active(nc)%num_urbanc,       &
               filter_inactive_and_active(nc)%urbanc,           &
               filter_inactive_and_active(nc)%num_urbanp,       &
               filter_inactive_and_active(nc)%urbanp,           &
               nextsw_cday, declinp1,                           &
               clm_fates,                                       &
               aerosol_inst, canopystate_inst, &
               water_inst%waterstatebulk_inst, &
               water_inst%waterdiagnosticbulk_inst, &
               lakestate_inst, temperature_inst, surfalb_inst)

          call t_stopf('surfalb')

          ! Albedos for urban columns
          if (filter_inactive_and_active(nc)%num_urbanl > 0) then
             call t_startf('urbalb')
             call UrbanAlbedo(bounds_clump,                  &
                  filter_inactive_and_active(nc)%num_urbanl, &
                  filter_inactive_and_active(nc)%urbanl,     &
                  filter_inactive_and_active(nc)%num_urbanc, &
                  filter_inactive_and_active(nc)%urbanc,     &
                  filter_inactive_and_active(nc)%num_urbanp, &
                  filter_inactive_and_active(nc)%urbanp,     &
                  water_inst%waterstatebulk_inst, &
                  water_inst%waterdiagnosticbulk_inst, &
                  urbanparams_inst,         &
                  solarabs_inst, surfalb_inst)
             call t_stopf('urbalb')
          end if

       end if

    end do
    !$OMP END PARALLEL DO


    ! Pass fates seed dispersal information to neighboring gridcells across
    ! all MPI tasks.  Note that WrapGlobalSeedDispersal calls an MPI collective routine
    ! and as such WrapGlobalSeedDispersal should be called outside of OMP threaded loop regions
    if (use_fates) then
       if (fates_seeddisp_cadence /= fates_dispersal_cadence_none) call clm_fates%WrapGlobalSeedDispersal()
    end if

    ! ============================================================================
    ! Determine gridcell averaged properties to send to atm
    ! ============================================================================

    ! BUG(wjs, 2019-07-12, ESCOMP/ctsm#762) Remove this block once above code is fully tracerized
    ! We need this call here so that tracer values are consistent with bulk in lnd2atm;
    ! without it, we get an error in CheckSnowConsistency.
    if (water_inst%DoConsistencyCheck()) then
       !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
       do nc = 1,nclumps
          call get_clump_bounds(nc, bounds_clump)
          call water_inst%ResetCheckedTracers(bounds_clump)
          call water_inst%TracerConsistencyCheck(bounds_clump, 'after pre-lnd2atm reset')
       end do
       !$OMP END PARALLEL DO
    end if

    call t_startf('lnd2atm')
    ! COMPILER_BUG(wjs, 2016-02-24, pgi 15.10) In principle, we should be able to make
    ! this function call inline in the CanopyFluxes argument list. However, with pgi
    ! 15.10, that results in the dummy argument having the wrong size (I suspect size 0,
    ! based on similar pgi compiler bugs that we have run into before). Also note that I
    ! don't have explicit bounds on the left-hand-side of this assignment: excluding these
    ! explicit bounds seemed to be needed to get around other compiler bugs.
    allocate(net_carbon_exchange_grc(bounds_proc%begg:bounds_proc%endg))
    net_carbon_exchange_grc = bgc_vegetation_inst%get_net_carbon_exchange_grc(bounds_proc)

    call lnd2atm(bounds_proc,                                            &
         atm2lnd_inst, surfalb_inst, temperature_inst, frictionvel_inst, &
         water_inst, &
         energyflux_inst, solarabs_inst, drydepvel_inst,       &
         vocemis_inst, fireemis_inst, dust_emis_inst, ch4_inst, glc_behavior, &
         lnd2atm_inst, &
         net_carbon_exchange_grc = net_carbon_exchange_grc(bounds_proc%begg:bounds_proc%endg))
    deallocate(net_carbon_exchange_grc)
    call t_stopf('lnd2atm')

    ! ============================================================================
    ! Determine gridcell averaged properties to send to glc
    ! ============================================================================

    call t_startf('lnd2glc')
    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1,nclumps
       call get_clump_bounds(nc, bounds_clump)
       call lnd2glc_inst%update_lnd2glc(bounds_clump,       &
            filter(nc)%num_do_smb_c, filter(nc)%do_smb_c,   &
            temperature_inst, water_inst%waterfluxbulk_inst, topo_inst,    &
            init=.false.)
    end do
    !$OMP END PARALLEL DO
    call t_stopf('lnd2glc')

    ! ==========================================================================
    ! Check the energy and water balance
    ! ==========================================================================

    call t_startf('balchk')
    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1,nclumps
       call get_clump_bounds(nc, bounds_clump)
       call WaterGridcellBalance(bounds_clump, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            filter(nc)%num_lakec, filter(nc)%lakec, &
            water_inst, lakestate_inst, &
            use_aquifer_layer = use_aquifer_layer(), flag = 'endwb')
       call BalanceCheck(bounds_clump, &
            filter(nc)%num_allc, filter(nc)%allc, &
            atm2lnd_inst, solarabs_inst, water_inst%waterfluxbulk_inst, &
            water_inst%waterstatebulk_inst, water_inst%waterdiagnosticbulk_inst, &
            water_inst%waterbalancebulk_inst, water_inst%wateratm2lndbulk_inst, &
            water_inst%waterlnd2atmbulk_inst, surfalb_inst, energyflux_inst, &
            canopystate_inst)
    end do
    !$OMP END PARALLEL DO
    call t_stopf('balchk')

    ! ============================================================================
    ! Write global average diagnostics to standard output
    ! ============================================================================

    nstep = get_nstep()

    ! ============================================================================
    ! Update accumulators
    ! ============================================================================

    ! FIX(SPM,032414) double check why this isn't called for ED
    ! FIX(SPM, 082814) - in the fates branch RF and I commented out the if(.not.
    ! use_fates) then statement ... double check if this is required and why

    call t_startf('accum')

    call atm2lnd_inst%UpdateAccVars(bounds_proc)

    call temperature_inst%UpdateAccVars(bounds_proc, crop_inst)

    call canopystate_inst%UpdateAccVars(bounds_proc)

    call water_inst%UpdateAccVars(bounds_proc)

    call energyflux_inst%UpdateAccVars(bounds_proc)

    ! COMPILER_BUG(wjs, 2014-11-30, pgi 14.7) For pgi 14.7 to be happy when
    ! compiling this threaded, I needed to change the dummy arguments to be
    ! pointers, and get rid of the explicit bounds in the subroutine call.
    ! call bgc_vegetation_inst%UpdateAccVars(bounds_proc, &
    !      t_a10_patch=temperature_inst%t_a10_patch(bounds_proc%begp:bounds_proc%endp), &
    !      t_ref2m_patch=temperature_inst%t_ref2m_patch(bounds_proc%begp:bounds_proc%endp))
    call bgc_vegetation_inst%UpdateAccVars(bounds_proc, &
         t_a10_patch=temperature_inst%t_a10_patch, &
         t_ref2m_patch=temperature_inst%t_ref2m_patch)

    if (use_crop) then
       call crop_inst%CropUpdateAccVars(bounds_proc, &
            temperature_inst%t_ref2m_patch, temperature_inst%t_soisno_col)
    end if

    if(use_fates) then
       call clm_fates%UpdateAccVars(bounds_proc)
    end if

    call t_stopf('accum')

    ! ============================================================================
    ! Update history buffer
    ! ============================================================================

    call t_startf('hbuf')
    call hist_update_hbuf(bounds_proc)
    call t_stopf('hbuf')

    ! ============================================================================
    ! Call dv (dynamic vegetation)
    ! ============================================================================

    if (use_cn) then
       call t_startf('endTSdynveg')
       nclumps = get_proc_clumps()
       !$OMP PARALLEL DO PRIVATE (nc,bounds_clump)
       do nc = 1,nclumps
          call get_clump_bounds(nc, bounds_clump)
          call bgc_vegetation_inst%EndOfTimeStepVegDynamics(bounds_clump, &
               filter(nc)%num_natvegp, filter(nc)%natvegp, &
               atm2lnd_inst, water_inst%wateratm2lndbulk_inst)
       end do
       !$OMP END PARALLEL DO
       call t_stopf('endTSdynveg')
    end if

    ! ============================================================================
    ! History/Restart output
    ! ============================================================================

    if (.not. use_noio) then

       call t_startf('clm_drv_io')

       ! Create history and write history tapes if appropriate
       call t_startf('clm_drv_io_htapes')

       call hist_htapes_wrapup( rstwr, nlend, bounds_proc,                      &
            soilstate_inst%watsat_col(bounds_proc%begc:bounds_proc%endc, 1:),   &
            soilstate_inst%sucsat_col(bounds_proc%begc:bounds_proc%endc, 1:),   &
            soilstate_inst%bsw_col(bounds_proc%begc:bounds_proc%endc, 1:),      &
            soilstate_inst%hksat_col(bounds_proc%begc:bounds_proc%endc, 1:),    &
            soilstate_inst%cellsand_col(bounds_proc%begc:bounds_proc%endc, 1:), &
            soilstate_inst%cellclay_col(bounds_proc%begc:bounds_proc%endc, 1:))

       call t_stopf('clm_drv_io_htapes')

       if (use_cn) then
          call bgc_vegetation_inst%WriteHistory(bounds_proc)
       end if

       ! Write restart/initial files if appropriate
       if (rstwr) then
          call t_startf('clm_drv_io_wrest')
          filer = restFile_filename(rdate=rdate)

          call restFile_write( bounds_proc, filer, &
               writing_finidat_interp_dest_file=.false., rdate=rdate )

          call t_stopf('clm_drv_io_wrest')
       end if
       call t_stopf('clm_drv_io')

    end if

    ! BUG(wjs, 2019-07-12, ESCOMP/ctsm#762) Remove this block once above code is fully tracerized
    ! (I'm not sure if we need this in addition to the call before lnd2atm, but it
    ! doesn't hurt to have this extra call to ResetCheckedTracers for now.)
    if (water_inst%DoConsistencyCheck()) then
       !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
       do nc = 1,nclumps
          call get_clump_bounds(nc, bounds_clump)
          call water_inst%ResetCheckedTracers(bounds_clump)
          call water_inst%TracerConsistencyCheck(bounds_clump, 'after end-of-loop reset')
       end do
       !$OMP END PARALLEL DO
    end if

  end subroutine clm_drv

  !-----------------------------------------------------------------------
  subroutine clm_drv_init(bounds, &
       num_nolakec, filter_nolakec, &
       num_nolakep, filter_nolakep, &
       num_soilp  , filter_soilp, &
       canopystate_inst, waterstatebulk_inst, &
       waterdiagnosticbulk_inst, energyflux_inst)
    !
    ! !DESCRIPTION:
    ! Initialization of clm driver variables needed from previous timestep
    !
    ! !USES:
    use shr_kind_mod       , only : r8 => shr_kind_r8
    use shr_infnan_mod     , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar         , only : nlevsno
    use CanopyStateType    , only : canopystate_type
    use WaterBalanceType   , only : waterbalance_type
    use WaterDiagnosticBulkType     , only : waterdiagnosticbulk_type
    use WaterFluxBulkType      , only : waterfluxbulk_type
    use WaterStateBulkType , only : waterstatebulk_type
    use EnergyFluxType     , only : energyflux_type
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)    :: bounds
    integer               , intent(in)    :: num_nolakec       ! number of non-lake points in column filter
    integer               , intent(in)    :: filter_nolakec(:) ! column filter for non-lake points
    integer               , intent(in)    :: num_nolakep       ! number of non-lake points in patch filter
    integer               , intent(in)    :: filter_nolakep(:) ! patch filter for non-lake points
    integer               , intent(in)    :: num_soilp         ! number of soil points in patch filter
    integer               , intent(in)    :: filter_soilp(:)   ! patch filter for soil points
    type(canopystate_type), intent(inout) :: canopystate_inst
    type(waterstatebulk_type) , intent(inout) :: waterstatebulk_inst
    type(waterdiagnosticbulk_type) , intent(inout) :: waterdiagnosticbulk_inst
    type(energyflux_type) , intent(inout) :: energyflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c, p, f, j              ! indices
    integer :: fp, fc                  ! filter indices
    !-----------------------------------------------------------------------

    associate(                                                             &
         snl                => col%snl                                   , & ! Input:  [integer  (:)   ]  number of snow layers

         h2osoi_ice         => waterstatebulk_inst%h2osoi_ice_col            , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)
         h2osoi_liq         => waterstatebulk_inst%h2osoi_liq_col            , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)
         frac_iceold        => waterdiagnosticbulk_inst%frac_iceold_col           , & ! Output: [real(r8) (:,:) ]  fraction of ice relative to the tot water
         elai               => canopystate_inst%elai_patch               , & ! Input:  [real(r8) (:)   ]  one-sided leaf area index with burying by snow
         esai               => canopystate_inst%esai_patch               , & ! Input:  [real(r8) (:)   ]  one-sided stem area index with burying by snow
         frac_veg_nosno     => canopystate_inst%frac_veg_nosno_patch     , & ! Output: [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]
         frac_veg_nosno_alb => canopystate_inst%frac_veg_nosno_alb_patch , & ! Output: [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]

         eflx_bot           => energyflux_inst%eflx_bot_col              , & ! Output: [real(r8) (:)   ]  heat flux from beneath soil/ice column (W/m**2)

         cisun_z            => photosyns_inst%cisun_z_patch              , & ! Output: [real(r8) (:)   ]  intracellular sunlit leaf CO2 (Pa)
         cisha_z            => photosyns_inst%cisha_z_patch                & ! Output: [real(r8) (:)   ]  intracellular shaded leaf CO2 (Pa)
         )

      ! Initialize intracellular CO2 (Pa) parameters each timestep for use in VOCEmission
      do p = bounds%begp,bounds%endp
         cisun_z(p,:) = -999._r8
         cisha_z(p,:) = -999._r8
      end do

      do c = bounds%begc,bounds%endc
         ! Reset flux from beneath soil/ice column
         eflx_bot(c)  = 0._r8
      end do

      ! Initialize fraction of vegetation not covered by snow

      do p = bounds%begp,bounds%endp
         if (patch%active(p)) then
            frac_veg_nosno(p) = frac_veg_nosno_alb(p)
         else
            frac_veg_nosno(p) = 0._r8
         end if
      end do

      ! Initialize set of previous time-step variables
      ! Ice fraction of snow at previous time step

      do j = -nlevsno+1,0
         do f = 1, num_nolakec
            c = filter_nolakec(f)
            if (j >= snl(c) + 1) then
               frac_iceold(c,j) = h2osoi_ice(c,j)/(h2osoi_liq(c,j)+h2osoi_ice(c,j))
            end if
         end do
      end do

    end associate

  end subroutine clm_drv_init

  !-----------------------------------------------------------------------
  subroutine clm_drv_patch2col (bounds, &
       num_allc, filter_allc, num_nolakec, filter_nolakec, &
       energyflux_inst, waterfluxbulk_inst)
    !
    ! !DESCRIPTION:
    ! Averages over all patches for variables defined over both soil and lake to provide
    ! the column-level averages of flux variables defined at the patch level.
    !
    ! !USES:
    use WaterFluxBulkType  , only : waterfluxbulk_type
    use EnergyFluxType , only : energyflux_type
    use subgridAveMod  , only : p2c
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)    :: bounds
    integer               , intent(in)    :: num_allc          ! number of column points in allc filter
    integer               , intent(in)    :: filter_allc(:)    ! column filter for all active points
    integer               , intent(in)    :: num_nolakec       ! number of column non-lake points in column filter
    integer               , intent(in)    :: filter_nolakec(:) ! column filter for non-lake points
    type(waterfluxbulk_type)  , intent(inout) :: waterfluxbulk_inst
    type(energyflux_type) , intent(inout) :: energyflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c,fc              ! indices
    ! -----------------------------------------------------------------

    ! Note: lake points are excluded from many of the following
    ! averages. For some fields, this is because the field doesn't
    ! apply over lakes. However, for many others, this is because the
    ! field is computed in LakeHydrologyMod, which is called after
    ! this routine; thus, for lakes, the column-level values of these
    ! fields are explicitly set in LakeHydrologyMod. (The fields that
    ! are included here for lakes are computed elsewhere, e.g., in
    ! LakeFluxesMod.)

    ! Averaging for patch evaporative flux variables

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterfluxbulk_inst%qflx_ev_snow_patch(bounds%begp:bounds%endp), &
         waterfluxbulk_inst%qflx_ev_snow_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterfluxbulk_inst%qflx_ev_soil_patch(bounds%begp:bounds%endp), &
         waterfluxbulk_inst%qflx_ev_soil_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterfluxbulk_inst%qflx_ev_h2osfc_patch(bounds%begp:bounds%endp), &
         waterfluxbulk_inst%qflx_ev_h2osfc_col(bounds%begc:bounds%endc))

    ! Averaging for patch water flux variables

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterfluxbulk_inst%qflx_evap_soi_patch(bounds%begp:bounds%endp), &
         waterfluxbulk_inst%qflx_evap_soi_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterfluxbulk_inst%qflx_evap_tot_patch(bounds%begp:bounds%endp), &
         waterfluxbulk_inst%qflx_evap_tot_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterfluxbulk_inst%qflx_tran_veg_patch(bounds%begp:bounds%endp), &
         waterfluxbulk_inst%qflx_tran_veg_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterfluxbulk_inst%qflx_liqevap_from_top_layer_patch(bounds%begp:bounds%endp), &
         waterfluxbulk_inst%qflx_liqevap_from_top_layer_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_allc, filter_allc, &
         waterfluxbulk_inst%qflx_evap_soi_patch(bounds%begp:bounds%endp), &
         waterfluxbulk_inst%qflx_evap_soi_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterfluxbulk_inst%qflx_liqdew_to_top_layer_patch(bounds%begp:bounds%endp), &
         waterfluxbulk_inst%qflx_liqdew_to_top_layer_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterfluxbulk_inst%qflx_solidevap_from_top_layer_patch(bounds%begp:bounds%endp), &
         waterfluxbulk_inst%qflx_solidevap_from_top_layer_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterfluxbulk_inst%qflx_soliddew_to_top_layer_patch(bounds%begp:bounds%endp), &
         waterfluxbulk_inst%qflx_soliddew_to_top_layer_col(bounds%begc:bounds%endc))

  end subroutine clm_drv_patch2col

  !------------------------------------------------------------------------
  subroutine write_diagnostic (bounds, nstep, lnd2atm_inst)
    !
    ! !DESCRIPTION:
    ! Write diagnostic surface temperature output each timestep.  Written to
    ! be fast but not bit-for-bit because order of summations can change each
    ! timestep.
    !
    ! !USES:
    use decompMod  , only : get_proc_global
    use spmdMod    , only : masterproc, npes, MPI_REAL8
    use spmdMod    , only : MPI_STATUS_SIZE, mpicom, MPI_SUM
    use shr_sys_mod, only : shr_sys_flush
    use abortutils , only : endrun
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use lnd2atmType, only : lnd2atm_type
    !
    ! !ARGUMENTS:
    type(bounds_type)  , intent(in) :: bounds
    integer            , intent(in) :: nstep      !model time step
    type(lnd2atm_type) , intent(in) :: lnd2atm_inst
    !
    ! !REVISION HISTORY:
    ! Created by Mariana Vertenstein
    !
    ! !LOCAL VARIABLES:
    integer :: p                       ! loop index
    integer :: numg                    ! total number of gridcells across all processors
    integer :: ier                     ! error status
    real(r8):: psum                    ! partial sum of ts
    real(r8):: tsum                    ! sum of ts
    real(r8):: tsxyav                  ! average ts for diagnostic output
    integer :: status(MPI_STATUS_SIZE) ! mpi status
    !------------------------------------------------------------------------

    call get_proc_global(ng=numg)

    if (masterproc) then
       write(iulog,*)'clm: completed timestep ',nstep
       call shr_sys_flush(iulog)
    end if

  end subroutine write_diagnostic

end module clm_driver
