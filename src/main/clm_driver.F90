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
  use clm_varctl             , only : wrtdia, iulog, use_fates
  use clm_varctl             , only : use_cn, use_lch4, use_noio, use_c13, use_c14
  use clm_varctl             , only : use_crop, irrigate, ndep_from_cpl
  use clm_time_manager       , only : get_nstep, is_beg_curr_day
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
  use BalanceCheckMod        , only : BeginWaterBalance, BalanceCheck
  !
  use CanopyTemperatureMod   , only : CanopyTemperature ! (formerly Biogeophysics1Mod)
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
  use HydrologyNoDrainageMod , only : CalcAndWithdrawIrrigationFluxes, HydrologyNoDrainage ! (formerly Hydrology2Mod)
  use HydrologyDrainageMod   , only : HydrologyDrainage   ! (formerly Hydrology2Mod)
  use CanopyHydrologyMod     , only : CanopyHydrology     ! (formerly Hydrology1Mod)
  use LakeHydrologyMod       , only : LakeHydrology
  use SoilWaterMovementMod   , only : use_aquifer_layer
  !
  use AerosolMod             , only : AerosolMasses  
  use SnowSnicarMod          , only : SnowAge_grain
  use SurfaceAlbedoMod       , only : SurfaceAlbedo
  use UrbanAlbedoMod         , only : UrbanAlbedo
  !
  use SurfaceRadiationMod    , only : SurfaceRadiation, CanopySunShadeFracs
  use UrbanRadiationMod      , only : UrbanRadiation
  !
  use SoilBiogeochemVerticalProfileMod   , only : SoilBiogeochemVerticalProfile
  use SatellitePhenologyMod  , only : SatellitePhenology, interpMonthlyVeg
  use ndepStreamMod          , only : ndep_interp
  use ActiveLayerMod         , only : alt_calc
  use ch4Mod                 , only : ch4, ch4_init_balance_check
  use DUSTMod                , only : DustDryDep, DustEmission
  use VOCEmissionMod         , only : VOCEmission
  !
  use filterMod              , only : setFilters
  !
  use atm2lndMod             , only : downscale_forcings, set_atm2lnd_water_tracers
  use lnd2atmMod             , only : lnd2atm
  use lnd2glcMod             , only : lnd2glc_type
  !
  use seq_drydep_mod         , only : n_drydep, drydep_method, DD_XLND
  use DryDepVelocity         , only : depvel_compute
  !
  use DaylengthMod           , only : UpdateDaylength
  use perf_mod
  !
  use clm_instMod            , only : nutrient_competition_method
  use GridcellType           , only : grc                
  use LandunitType           , only : lun                
  use ColumnType             , only : col                
  use PatchType              , only : patch                
  use clm_instMod
  use clm_instMod            , only : soil_water_retention_curve
  use EDBGCDynMod            , only : EDBGCDyn, EDBGCDynSummary
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
    use clm_time_manager, only : get_curr_date    
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
    ! ============================================================================

    if (use_cn) then 
       ! For dry-deposition need to call CLMSP so that mlaidiff is obtained
       if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
          call t_startf('interpMonthlyVeg')
          call interpMonthlyVeg(bounds_proc, canopystate_inst)
          call t_stopf('interpMonthlyVeg')
       endif

    else
       ! Determine weights for time interpolation of monthly vegetation data.
       ! This also determines whether it is time to read new monthly vegetation and
       ! obtain updated leaf area index [mlai1,mlai2], stem area index [msai1,msai2],
       ! vegetation top [mhvt1,mhvt2] and vegetation bottom [mhvb1,mhvb2]. The
       ! weights obtained here are used in subroutine SatellitePhenology to obtain time
       ! interpolated values.
       if (doalb .or. ( n_drydep > 0 .and. drydep_method == DD_XLND )) then
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
       call alt_calc(filter_inactive_and_active(nc)%num_soilc, filter_inactive_and_active(nc)%soilc, &
            temperature_inst, canopystate_inst) 

       if (use_cn) then
          call SoilBiogeochemVerticalProfile(bounds_clump                                       , &
               filter_inactive_and_active(nc)%num_soilc, filter_inactive_and_active(nc)%soilc   , &
               filter_inactive_and_active(nc)%num_soilp, filter_inactive_and_active(nc)%soilp   , &
               canopystate_inst, soilstate_inst, soilbiogeochem_state_inst)
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

    ! ============================================================================
    ! Update subgrid weights with dynamic landcover (prescribed transient patches,
    ! CNDV, and or dynamic landunits), and do related adjustments. Note that this
    ! call needs to happen outside loops over nclumps.
    ! ============================================================================

    call t_startf('dyn_subgrid')
    call dynSubgrid_driver(bounds_proc,                                               &
         urbanparams_inst, soilstate_inst, water_inst,                       &
         temperature_inst, energyflux_inst,          &
         canopystate_inst, photosyns_inst, crop_inst, glc2lnd_inst, bgc_vegetation_inst, &
         soilbiogeochem_state_inst, soilbiogeochem_carbonstate_inst,                  &
         c13_soilbiogeochem_carbonstate_inst, c14_soilbiogeochem_carbonstate_inst,    &
         soilbiogeochem_nitrogenstate_inst, soilbiogeochem_carbonflux_inst, ch4_inst, &
         glc_behavior)
    call t_stopf('dyn_subgrid')

    ! ============================================================================
    ! Initialize the column-level mass balance checks for water, carbon & nitrogen.
    !
    ! For water: Currently, I believe this needs to be done after weights are updated for
    ! prescribed transient patches or CNDV, because column-level water is not generally
    ! conserved when weights change (instead the difference is put in the grid cell-level
    ! terms, qflx_liq_dynbal, etc.). In the future, we may want to change the balance
    ! checks to ensure that the grid cell-level water is conserved, considering
    ! qflx_liq_dynbal; in this case, the call to BeginWaterBalance should be moved to
    ! before the weight updates.
    !
    ! For carbon & nitrogen: This needs to be done after dynSubgrid_driver, because the
    ! changes due to dynamic area adjustments can break column-level conservation
    ! ============================================================================

    !$OMP PARALLEL DO PRIVATE (nc,bounds_clump)
    do nc = 1,nclumps
       call get_clump_bounds(nc, bounds_clump)

       call t_startf('begwbal')
       call BeginWaterBalance(bounds_clump,                   &
            filter(nc)%num_nolakec, filter(nc)%nolakec,       &
            filter(nc)%num_lakec, filter(nc)%lakec,           &
            water_inst, soilhydrology_inst, &
            use_aquifer_layer = use_aquifer_layer())

       call t_stopf('begwbal')

       call t_startf('begcnbal_col')
       if (use_cn) then
          call bgc_vegetation_inst%InitColumnBalance(bounds_clump, &
               filter(nc)%num_allc, filter(nc)%allc, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               filter(nc)%num_soilp, filter(nc)%soilp, &
               soilbiogeochem_carbonstate_inst, &
               c13_soilbiogeochem_carbonstate_inst, &
               c14_soilbiogeochem_carbonstate_inst, &
               soilbiogeochem_nitrogenstate_inst)
       end if

       if (use_lch4) then
          call ch4_init_balance_check(bounds_clump, &
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

    if (use_cn) then
       call t_startf('bgc_interp')
       if (.not. ndep_from_cpl) then
          call ndep_interp(bounds_proc, atm2lnd_inst)
       end if
       call bgc_vegetation_inst%InterpFileInputs(bounds_proc)
       call t_stopf('bgc_interp')
    end if

    ! Get time varying urban data
    call urbantv_inst%urbantv_interp(bounds_proc)

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
            filter(nc)%num_icemecc, filter(nc)%icemecc, &
            glc2lnd_inst, glc_behavior, &
            atm_topo = atm2lnd_inst%forc_topo_grc(bounds_clump%begg:bounds_clump%endg))

       call downscale_forcings(bounds_clump, &
            topo_inst, atm2lnd_inst, water_inst%wateratm2lndbulk_inst, &
            eflx_sh_precip_conversion = energyflux_inst%eflx_sh_precip_conversion_col(bounds_clump%begc:bounds_clump%endc))

       call set_atm2lnd_water_tracers(bounds_clump, &
            filter(nc)%num_allc, filter(nc)%allc, &
            water_inst)

       if (water_inst%DoConsistencyCheck()) then
          ! BUG(wjs, 2018-09-05, ESCOMP/ctsm#498) Eventually do tracer consistency checks
          ! every time step
          if (get_nstep() == 0) then
             call t_startf("tracer_consistency_check")
             call water_inst%TracerConsistencyCheck(bounds_clump, 'after downscale_forcings')
             call t_stopf("tracer_consistency_check")
          end if
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
          ! BUG(wjs, 2018-09-05, ESCOMP/ctsm#498) Eventually do tracer consistency checks
          ! every time step
          if (get_nstep() == 0) then
             call t_startf("tracer_consistency_check")
             call water_inst%TracerConsistencyCheck(bounds_clump, 'after CalcAndWithdrawIrrigationFluxes')
             call t_stopf("tracer_consistency_check")
          end if
       end if

       ! ============================================================================
       ! Canopy Hydrology
       ! (1) water storage of intercepted precipitation
       ! (2) direct throughfall and canopy drainage of precipitation
       ! (3) fraction of foliage covered by water and the fraction is dry and transpiring
       ! (4) snow layer initialization if the snow accumulation exceeds 10 mm.
       ! ============================================================================

       call t_startf('canhydro')
       call CanopyHydrology(bounds_clump, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            filter(nc)%num_nolakep, filter(nc)%nolakep, &
            atm2lnd_inst, canopystate_inst, temperature_inst, &
            aerosol_inst, water_inst%waterstatebulk_inst, &
            water_inst%waterdiagnosticbulk_inst, &
            water_inst%waterfluxbulk_inst, &
            water_inst%wateratm2lndbulk_inst)
       call t_stopf('canhydro')

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
       call CanopyTemperature(bounds_clump,                                   &
            filter(nc)%num_nolakec, filter(nc)%nolakec,                       &
            filter(nc)%num_nolakep, filter(nc)%nolakep,                       &
            clm_fates,                                                        &
            atm2lnd_inst, canopystate_inst, soilstate_inst, frictionvel_inst, &
            water_inst%waterstatebulk_inst, water_inst%waterdiagnosticbulk_inst, &
            water_inst%waterfluxbulk_inst, water_inst%wateratm2lndbulk_inst,  &
            energyflux_inst, temperature_inst)
       call t_stopf('bgp1')

       ! ============================================================================
       ! Determine fluxes
       ! ============================================================================

       call t_startf('bgflux')

       ! Bareground fluxes for all patches except lakes and urban landunits

       call BareGroundFluxes(bounds_clump,                                 &
            filter(nc)%num_noexposedvegp, filter(nc)%noexposedvegp,          &
            atm2lnd_inst, soilstate_inst,                &
            frictionvel_inst, ch4_inst, energyflux_inst, temperature_inst, &
            water_inst%waterfluxbulk_inst, water_inst%waterstatebulk_inst, &
            water_inst%waterdiagnosticbulk_inst, water_inst%wateratm2lndbulk_inst, &
            photosyns_inst, humanindex_inst)
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
            atm2lnd_inst, canopystate_inst,                                                 &
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
       call DustEmission(bounds_clump,                                       &
            filter(nc)%num_nolakep, filter(nc)%nolakep,                      &
            atm2lnd_inst, soilstate_inst, canopystate_inst, &
            water_inst%waterstatebulk_inst, water_inst%waterdiagnosticbulk_inst, &
            frictionvel_inst, dust_inst)

       ! Dust dry deposition (C. Zender's modified codes)
       call DustDryDep(bounds_clump, &
            atm2lnd_inst, frictionvel_inst, dust_inst)

       ! VOC emission (A. Guenther's MEGAN (2006) model)
       call VOCEmission(bounds_clump,                                         &
               filter(nc)%num_soilp, filter(nc)%soilp,                           &
               atm2lnd_inst, canopystate_inst, photosyns_inst, temperature_inst, &
               vocemis_inst)

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
            water_inst%waterfluxbulk_inst, water_inst%waterstatebulk_inst, &
            water_inst%waterdiagnosticbulk_inst, soilhydrology_inst, &
            saturated_excess_runoff_inst, &
            infiltration_excess_runoff_inst, &
            aerosol_inst, canopystate_inst, soil_water_retention_curve, topo_inst)

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
            atm2lnd_inst, temperature_inst, soilstate_inst, water_inst%waterstatebulk_inst, &
            water_inst%waterdiagnosticbulk_inst, water_inst%waterbalancebulk_inst, &
            water_inst%waterfluxbulk_inst, water_inst%wateratm2lndbulk_inst, &
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
             
       if (use_cn) then 
          call t_startf('ecosysdyn')
          call bgc_vegetation_inst%EcosystemDynamicsPreDrainage(bounds_clump,            &
                  filter(nc)%num_soilc, filter(nc)%soilc,                       &
                  filter(nc)%num_soilp, filter(nc)%soilp,                       &
                  filter(nc)%num_pcropp, filter(nc)%pcropp, doalb,              &
               soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst,         &
               c13_soilbiogeochem_carbonflux_inst, c13_soilbiogeochem_carbonstate_inst, &
               c14_soilbiogeochem_carbonflux_inst, c14_soilbiogeochem_carbonstate_inst, &
               soilbiogeochem_state_inst,                                               &
               soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst,     &
               atm2lnd_inst, water_inst%waterstatebulk_inst, &
               water_inst%waterdiagnosticbulk_inst, water_inst%waterfluxbulk_inst,      &
               water_inst%wateratm2lndbulk_inst, canopystate_inst, soilstate_inst, temperature_inst, crop_inst, ch4_inst, &
               photosyns_inst, saturated_excess_runoff_inst, energyflux_inst,          &
               nutrient_competition_method, fireemis_inst)

          call t_stopf('ecosysdyn')

       end if

                ! Prescribed biogeography - prescribed canopy structure, some prognostic carbon fluxes

       if ((.not. use_cn) .and. (.not. use_fates) .and. (doalb)) then 
          call t_startf('SatellitePhenology')
          call SatellitePhenology(bounds_clump, filter(nc)%num_nolakep, filter(nc)%nolakep, &
               water_inst%waterdiagnosticbulk_inst, canopystate_inst)
          call t_stopf('SatellitePhenology')
       end if

       ! Dry Deposition of chemical tracers (Wesely (1998) parameterizaion)
          
       call t_startf('depvel')
       call depvel_compute(bounds_clump, &
            atm2lnd_inst, canopystate_inst, water_inst%waterstatebulk_inst, &
            water_inst%waterdiagnosticbulk_inst, water_inst%wateratm2lndbulk_inst, &
            frictionvel_inst, photosyns_inst, drydepvel_inst)
       call t_stopf('depvel')

       ! ============================================================================
       ! Calculate soil/snow hydrology with drainage (subsurface runoff)
       ! ============================================================================

       call t_startf('hydro2_drainage')

       call HydrologyDrainage(bounds_clump,                   &
            filter(nc)%num_nolakec, filter(nc)%nolakec,       &
            filter(nc)%num_hydrologyc, filter(nc)%hydrologyc, &
            filter(nc)%num_urbanc, filter(nc)%urbanc,         &                 
            filter(nc)%num_do_smb_c, filter(nc)%do_smb_c,     &                
            atm2lnd_inst, glc2lnd_inst, temperature_inst,     &
            soilhydrology_inst, soilstate_inst, water_inst%waterstatebulk_inst, &
            water_inst%waterdiagnosticbulk_inst, water_inst%waterbalancebulk_inst, &
            water_inst%waterfluxbulk_inst, water_inst%wateratm2lndbulk_inst, &
            glacier_smb_inst)

       call t_stopf('hydro2_drainage')     

       if (use_cn) then

          call t_startf('EcosysDynPostDrainage')     
          call bgc_vegetation_inst%EcosystemDynamicsPostDrainage(bounds_clump, &
               filter(nc)%num_allc, filter(nc)%allc, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               filter(nc)%num_soilp, filter(nc)%soilp, &
               doalb, crop_inst, &
               water_inst%waterstatebulk_inst, water_inst%waterdiagnosticbulk_inst, &
               water_inst%waterfluxbulk_inst, frictionvel_inst, canopystate_inst, &
               soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst, &
               c13_soilbiogeochem_carbonflux_inst, c13_soilbiogeochem_carbonstate_inst, &
               c14_soilbiogeochem_carbonflux_inst, c14_soilbiogeochem_carbonstate_inst, &
               soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst)
          call t_stopf('EcosysDynPostDrainage')     

       end if

       if ( use_fates  .and. is_beg_curr_day() ) then ! run fates at the start of each day
          
          if ( masterproc ) then
             write(iulog,*)  'clm: calling FATES model ', get_nstep()
          end if

          call clm_fates%dynamics_driv( nc, bounds_clump,                        &
               atm2lnd_inst, soilstate_inst, temperature_inst,                   &
               water_inst%waterstatebulk_inst, water_inst%waterdiagnosticbulk_inst, &
               water_inst%wateratm2lndbulk_inst, canopystate_inst, soilbiogeochem_carbonflux_inst,&
               frictionvel_inst)
          
          ! TODO(wjs, 2016-04-01) I think this setFilters call should be replaced by a
          ! call to reweight_wrapup, if it's needed at all.
          call setFilters( bounds_clump, glc_behavior )
          
       end if ! use_fates branch
       
       
       if ( use_fates ) then
          
          call EDBGCDyn(bounds_clump,                                                              &
               filter(nc)%num_soilc, filter(nc)%soilc,                                             &
               filter(nc)%num_soilp, filter(nc)%soilp,                                             &
               filter(nc)%num_pcropp, filter(nc)%pcropp, doalb,                                    &
               bgc_vegetation_inst%cnveg_carbonflux_inst, &
               bgc_vegetation_inst%cnveg_carbonstate_inst, &
               soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst,                    &
               soilbiogeochem_state_inst,                                                          &
               soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst,                &
               c13_soilbiogeochem_carbonstate_inst, c13_soilbiogeochem_carbonflux_inst,            &
               c14_soilbiogeochem_carbonstate_inst, c14_soilbiogeochem_carbonflux_inst,            &
               atm2lnd_inst, water_inst%waterfluxbulk_inst,                                      &
               canopystate_inst, soilstate_inst, temperature_inst, crop_inst, ch4_inst)

          call EDBGCDynSummary(bounds_clump,                                             &
                filter(nc)%num_soilc, filter(nc)%soilc,                                  &
                filter(nc)%num_soilp, filter(nc)%soilp,                                  &
                soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst,         &
                c13_soilbiogeochem_carbonflux_inst, c13_soilbiogeochem_carbonstate_inst, &
                c14_soilbiogeochem_carbonflux_inst, c14_soilbiogeochem_carbonstate_inst, &
                soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst,     &
                clm_fates, nc)
       end if



       ! ============================================================================
       ! Check the energy and water balance and also carbon and nitrogen balance
       ! ============================================================================

       call t_startf('balchk')
       call BalanceCheck(bounds_clump, &
            atm2lnd_inst, solarabs_inst, water_inst%waterfluxbulk_inst, &
            water_inst%waterstatebulk_inst, water_inst%waterdiagnosticbulk_inst, &
            water_inst%waterbalancebulk_inst, water_inst%wateratm2lndbulk_inst, &
            surfalb_inst, energyflux_inst, canopystate_inst)
       call t_stopf('balchk')

       ! ============================================================================
       ! Check the carbon and nitrogen balance
       ! ============================================================================

       if (use_cn) then
          call t_startf('cnbalchk')
          call bgc_vegetation_inst%BalanceCheck( &
               bounds_clump, filter(nc)%num_soilc, filter(nc)%soilc, &
               soilbiogeochem_carbonflux_inst, soilbiogeochem_nitrogenflux_inst)
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

       if (doalb) then

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

    ! ============================================================================
    ! Determine gridcell averaged properties to send to atm
    ! ============================================================================

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
         vocemis_inst, fireemis_inst, dust_inst, ch4_inst, glc_behavior, &
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

    ! ============================================================================
    ! Write global average diagnostics to standard output
    ! ============================================================================

    nstep = get_nstep()
    if (wrtdia) call mpi_barrier(mpicom,ier)
    call t_startf('wrtdiag')
    call write_diagnostic(bounds_proc, wrtdia, nstep, lnd2atm_inst)
    call t_stopf('wrtdiag')

    ! ============================================================================
    ! Update accumulators
    ! ============================================================================

    ! FIX(SPM,032414) double check why this isn't called for ED
    ! FIX(SPM, 082814) - in the fates branch RF and I commented out the if(.not.
    ! use_fates) then statement ... double check if this is required and why

    if (nstep > 0) then
       call t_startf('accum')

       call atm2lnd_inst%UpdateAccVars(bounds_proc)

       call temperature_inst%UpdateAccVars(bounds_proc)

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

       call t_stopf('accum')
    end if

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

       call hist_htapes_wrapup( rstwr, nlend, bounds_proc,                    &
            soilstate_inst%watsat_col(bounds_proc%begc:bounds_proc%endc, 1:), &
            soilstate_inst%sucsat_col(bounds_proc%begc:bounds_proc%endc, 1:), &
            soilstate_inst%bsw_col(bounds_proc%begc:bounds_proc%endc, 1:),    &
            soilstate_inst%hksat_col(bounds_proc%begc:bounds_proc%endc, 1:))

       call t_stopf('clm_drv_io_htapes')

       if (use_cn) then
          call bgc_vegetation_inst%WriteHistory(bounds_proc)
       end if

       ! Write restart/initial files if appropriate
       if (rstwr) then
          call t_startf('clm_drv_io_wrest')
          filer = restFile_filename(rdate=rdate)

          call restFile_write( bounds_proc, filer, rdate=rdate )

          call t_stopf('clm_drv_io_wrest')
       end if
       call t_stopf('clm_drv_io')

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
         waterfluxbulk_inst%qflx_rain_grnd_patch(bounds%begp:bounds%endp), &
         waterfluxbulk_inst%qflx_rain_grnd_col(bounds%begc:bounds%endc))
    
    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterfluxbulk_inst%qflx_snow_grnd_patch(bounds%begp:bounds%endp), &
         waterfluxbulk_inst%qflx_snow_grnd_col(bounds%begc:bounds%endc))
    
    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterfluxbulk_inst%qflx_tran_veg_patch(bounds%begp:bounds%endp), &
         waterfluxbulk_inst%qflx_tran_veg_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterfluxbulk_inst%qflx_evap_grnd_patch(bounds%begp:bounds%endp), &
         waterfluxbulk_inst%qflx_evap_grnd_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_allc, filter_allc, &
         waterfluxbulk_inst%qflx_evap_soi_patch(bounds%begp:bounds%endp), &
         waterfluxbulk_inst%qflx_evap_soi_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterfluxbulk_inst%qflx_prec_grnd_patch(bounds%begp:bounds%endp), &
         waterfluxbulk_inst%qflx_prec_grnd_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterfluxbulk_inst%qflx_dew_grnd_patch(bounds%begp:bounds%endp), &
         waterfluxbulk_inst%qflx_dew_grnd_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterfluxbulk_inst%qflx_sub_snow_patch(bounds%begp:bounds%endp), &
         waterfluxbulk_inst%qflx_sub_snow_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterfluxbulk_inst%qflx_dew_snow_patch(bounds%begp:bounds%endp), &
         waterfluxbulk_inst%qflx_dew_snow_col(bounds%begc:bounds%endc))

  end subroutine clm_drv_patch2col

  !------------------------------------------------------------------------
  subroutine write_diagnostic (bounds, wrtdia, nstep, lnd2atm_inst)
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
    logical            , intent(in) :: wrtdia     !true => write diagnostic
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

    if (wrtdia) then

       call t_barrierf('sync_write_diag', mpicom)
       psum = sum(lnd2atm_inst%t_rad_grc(bounds%begg:bounds%endg))
       call mpi_reduce(psum, tsum, 1, MPI_REAL8, MPI_SUM, 0, mpicom, ier)
       if (ier/=0) then
          write(iulog,*) 'write_diagnostic: Error in mpi_reduce()'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       if (masterproc) then
          tsxyav = tsum / numg
          write(iulog,1000) nstep, tsxyav
          call shr_sys_flush(iulog)
       end if

    else

       if (masterproc) then
          write(iulog,*)'clm: completed timestep ',nstep
          call shr_sys_flush(iulog)
       end if

    endif

1000 format (1x,'nstep = ',i10,'   TS = ',f21.15)

  end subroutine write_diagnostic

end module clm_driver
