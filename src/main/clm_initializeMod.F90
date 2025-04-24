module clm_initializeMod

  !-----------------------------------------------------------------------
  ! Performs land model initialization
  !-----------------------------------------------------------------------

  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_sys_mod           , only : shr_sys_flush
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use spmdMod               , only : masterproc, mpicom
  use decompMod             , only : bounds_type, get_proc_bounds, get_proc_clumps, get_clump_bounds
  use abortutils            , only : endrun
  use clm_varctl            , only : nsrest, nsrStartup, nsrContinue, nsrBranch
  use clm_varctl            , only : use_fates_sp, use_fates_bgc, use_fates
  use clm_varctl            , only : is_cold_start
  use clm_varctl            , only : iulog
  use clm_varctl            , only : use_lch4, use_cn, use_cndv, use_c13, use_c14, nhillslope
  use clm_varctl            , only : use_soil_moisture_streams
  use clm_instur            , only : wt_lunit, urban_valid, wt_nat_patch, wt_cft, fert_cft
  use clm_instur            , only : irrig_method, wt_glc_mec, topo_glc_mec, pct_lake_max, pct_urban_max, ncolumns_hillslope
  use perf_mod              , only : t_startf, t_stopf
  use readParamsMod         , only : readParameters
  use ncdio_pio             , only : file_desc_t
  use GridcellType          , only : grc           ! instance
  use LandunitType          , only : lun           ! instance
  use ColumnType            , only : col           ! instance
  use PatchType             , only : patch         ! instance
  use reweightMod           , only : reweight_wrapup
  use filterMod             , only : allocFilters, filter, filter_inactive_and_active
  use CLMFatesInterfaceMod  , only : CLMFatesGlobals1,CLMFatesGlobals2
  use CLMFatesInterfaceMod  , only : CLMFatesTimesteps
  use dynSubgridControlMod  , only : dynSubgridControl_init, get_reset_dynbal_baselines
  use SelfTestDriver        , only : self_test_driver
  use SoilMoistureStreamMod , only : PrescribedSoilMoistureInit
  use clm_instMod
  !
  implicit none
  private  ! By default everything is private
  !
  public :: initialize1  ! Phase one initialization
  public :: initialize2  ! Phase two initialization

  integer :: actual_numcft  ! numcft from sfc dataset
  integer :: actual_nlevurb ! nlevurb from sfc dataset
  integer :: actual_numpft  ! numpft from sfc dataset

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

  subroutine initialize1(dtime)
    !
    ! !DESCRIPTION:
    ! CLM initialization first phase
    !
    ! !USES:
    use clm_varpar           , only: clm_varpar_init
    use clm_varcon           , only: clm_varcon_init
    use landunit_varcon      , only: landunit_varcon_init
    use clm_varctl           , only: fsurdat, version
    use surfrdMod            , only: surfrd_get_num_patches, surfrd_get_nlevurb, surfrd_compat_check
    use controlMod           , only: control_init, control_print, NLFilename
    use ncdio_pio            , only: ncd_pio_init
    use initGridCellsMod     , only: initGridCells
    use UrbanParamsType      , only: IsSimpleBuildTemp
    use dynSubgridControlMod , only: dynSubgridControl_init
    use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_par_init
    use CropReprPoolsMod     , only: crop_repr_pools_init
    use HillslopeHydrologyMod, only: hillslope_properties_init
    !
    ! !ARGUMENTS
    integer, intent(in) :: dtime    ! model time step (seconds)
    !
    ! !LOCAL VARIABLES:
    integer           :: ier                     ! error status
    integer           :: i,j,n,k,c,l,g           ! indices
    integer           :: nl                      ! gdc and glo lnd indices
    integer           :: ns, ni, nj              ! global grid sizes
    integer           :: begg, endg              ! processor bounds
    type(bounds_type) :: bounds_proc
    type(bounds_type) :: bounds_clump
    integer           :: nclumps                 ! number of clumps on this processor
    integer           :: nc                      ! clump index
    integer           :: actual_maxsoil_patches  ! value from surface dataset
    integer ,pointer  :: amask(:)                ! global land mask
    character(len=32) :: subname = 'initialize1' ! subroutine name
    !-----------------------------------------------------------------------

    call t_startf('clm_init1')

    ! Initialize run control variables, timestep

    if ( masterproc )then
       write(iulog,*) trim(version)
       write(iulog,*)
       write(iulog,*) 'Attempting to initialize the land model .....'
       write(iulog,*)
       call shr_sys_flush(iulog)
    endif

    call control_init(dtime)
    call ncd_pio_init()
    call surfrd_compat_check(fsurdat)
    call surfrd_get_num_patches(fsurdat, actual_maxsoil_patches, actual_numpft, actual_numcft)
    call surfrd_get_nlevurb(fsurdat, actual_nlevurb)

    ! If fates is on, we override actual_maxsoil_patches. FATES dictates the
    ! number of patches per column.  We still use numcft from the surface
    ! file though...
    if(use_fates) then
       call CLMFatesGlobals1(actual_numpft, actual_numcft, actual_maxsoil_patches)
    end if

    call clm_varpar_init(actual_maxsoil_patches, actual_numpft, actual_numcft, actual_nlevurb)
    call decomp_cascade_par_init( NLFilename )
    call clm_varcon_init( IsSimpleBuildTemp() )
    call landunit_varcon_init()
    if (masterproc) call control_print()
    call dynSubgridControl_init(NLFilename)
    call crop_repr_pools_init()
    call hillslope_properties_init(NLFilename)

    call t_stopf('clm_init1')

  end subroutine initialize1

  !-----------------------------------------------------------------------
  subroutine initialize2(ni,nj, currtime)
    !
    ! !DESCRIPTION:
    ! CLM initialization second phase
    !
    ! !USES:
    use ESMF                          , only : ESMF_Time
    use clm_varcon                    , only : spval
    use clm_varpar                    , only : natpft_lb, natpft_ub, cft_lb, cft_ub, maxpatch_glc
    use clm_varpar                    , only : surfpft_lb, surfpft_ub
    use clm_varpar                    , only : nlevsno
    use clm_varpar                    , only : natpft_size,cft_size
    use clm_varctl                    , only : fsurdat, hillslope_file
    use clm_varctl                    , only : finidat, finidat_interp_source, finidat_interp_dest
    use clm_varctl                    , only : use_cn, use_fates, use_fates_luh
    use clm_varctl                    , only : use_crop, ndep_from_cpl, fates_spitfire_mode
    use clm_varctl                    , only : use_hillslope
    use clm_varorb                    , only : eccen, mvelpp, lambm0, obliqr
    use clm_varctl                    , only : use_cropcal_streams
    use landunit_varcon               , only : landunit_varcon_init, max_lunit, numurbl
    use pftconMod                     , only : pftcon
    use decompInitMod                 , only : decompInit_clumps, decompInit_glcp
    use domainMod                     , only : domain_check, ldomain, domain_init
    use surfrdMod                     , only : surfrd_get_data
    use controlMod                    , only : NLFilename, fluh_timeseries
    use initGridCellsMod              , only : initGridCells
    use ch4varcon                     , only : ch4conrd
    use UrbanParamsType               , only : UrbanInput, IsSimpleBuildTemp
    use shr_orb_mod                   , only : shr_orb_decl
    use shr_drydep_mod                , only : n_drydep
    use accumulMod                    , only : print_accum_fields
    use clm_time_manager              , only : get_step_size_real, get_curr_calday
    use clm_time_manager              , only : get_curr_date, get_nstep, advance_timestep
    use clm_time_manager              , only : timemgr_init, timemgr_restart_io, timemgr_restart, is_restart
    use CIsoAtmTimeseriesMod          , only : C14_init_BombSpike, use_c14_bombspike, C13_init_TimeSeries, use_c13_timeseries
    use DaylengthMod                  , only : InitDaylength
    use dynSubgridDriverMod           , only : dynSubgrid_init
    use dynConsBiogeophysMod          , only : dyn_hwcontent_set_baselines
    use fileutils                     , only : getfil
    use initInterpMod                 , only : initInterp
    use subgridWeightsMod             , only : init_subgrid_weights_mod
    use histFileMod                   , only : hist_htapes_build, htapes_fieldlist, hist_printflds
    use histFileMod                   , only : hist_addfld1d, hist_addfld2d, no_snow_normal
    use restFileMod                   , only : restFile_getfile, restFile_open, restFile_close
    use restFileMod                   , only : restFile_read, restFile_write
    use ndepStreamMod                 , only : ndep_init, ndep_interp
    use cropcalStreamMod              , only : cropcal_init, cropcal_interp, cropcal_advance
    use LakeCon                       , only : LakeConInit
    use SatellitePhenologyMod         , only : SatellitePhenologyInit, readAnnualVegetation, interpMonthlyVeg
    use SatellitePhenologyMod         , only : CalcSatellitePhenologyTimeInterp
    use SnowSnicarMod                 , only : SnowAge_init, SnowOptics_init
    use lnd2atmMod                    , only : lnd2atm_minimal
    use controlMod                    , only : NLFilename, check_missing_initdata_status
    use clm_instMod                   , only : clm_fates
    use BalanceCheckMod               , only : BalanceCheckInit
    use CNSharedParamsMod             , only : CNParamsSetSoilDepth
    use NutrientCompetitionFactoryMod , only : create_nutrient_competition_method
    use FATESFireFactoryMod           , only : scalar_lightning
    use dynFATESLandUseChangeMod      , only : dynFatesLandUseInit
    use HillslopeHydrologyMod         , only : InitHillslope
    !
    ! !ARGUMENTS
    integer, intent(in) :: ni, nj         ! global grid sizes
    type(ESMF_Time), intent(in) :: currtime
    !
    ! !LOCAL VARIABLES:
    integer            :: c,g,i,j,k,l,n,p ! indices
    integer            :: yr              ! current year (0, ...)
    integer            :: mon             ! current month (1 -> 12)
    integer            :: day             ! current day (1 -> 31)
    integer            :: ncsec           ! current time of day [seconds]
    character(len=256) :: fnamer          ! name of netcdf restart file
    character(len=256) :: pnamer          ! full pathname of netcdf restart file
    character(len=256) :: locfn           ! local file name
    type(file_desc_t)  :: ncid            ! netcdf id
    real(r8)           :: dtime           ! time step increment (sec)
    integer            :: nstep           ! model time step
    real(r8)           :: calday          ! calendar day for nstep
    real(r8)           :: caldaym1        ! calendar day for nstep-1
    real(r8)           :: declin          ! solar declination angle in radians for nstep
    real(r8)           :: declinm1        ! solar declination angle in radians for nstep-1
    real(r8)           :: eccf            ! earth orbit eccentricity factor
    type(bounds_type)  :: bounds_proc     ! processor bounds
    type(bounds_type)  :: bounds_clump    ! clump bounds
    integer            :: nclumps         ! number of clumps on this processor
    integer            :: nc              ! clump index
    logical            :: reset_dynbal_baselines_all_columns
    logical            :: reset_dynbal_baselines_lake_columns
    integer            :: begg, endg
    integer            :: iun
    integer            :: klen
    integer            :: ioe
    integer            :: ier
    logical            :: lexists
    real(r8), pointer  :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
    character(len=32)  :: subname = 'initialize2' ! subroutine name
    !-----------------------------------------------------------------------

    call t_startf('clm_init2')

    ! Get processor bounds for gridcells
    call get_proc_bounds(bounds_proc)
    begg = bounds_proc%begg; endg = bounds_proc%endg

    ! Initialize glc behavior
    call glc_behavior%Init(begg, endg, NLFilename)

    ! Initialize urban model input (initialize urbinp data structure)
    ! This needs to be called BEFORE the call to surfrd_get_data since
    ! that will call surfrd_get_special which in turn calls check_urban
    call UrbanInput(begg, endg, mode='initialize')

    ! Allocate surface grid dynamic memory (just gridcell bounds dependent)
    allocate (wt_lunit     (begg:endg, max_lunit           ))
    allocate (urban_valid  (begg:endg                      ))
    allocate (wt_cft       (begg:endg, cft_lb:cft_ub       ))
    allocate (fert_cft     (begg:endg, cft_lb:cft_ub       ))
    allocate (irrig_method (begg:endg, cft_lb:cft_ub       ))
    allocate (wt_glc_mec   (begg:endg, maxpatch_glc     ))
    allocate (topo_glc_mec (begg:endg, maxpatch_glc     ))
    allocate (pct_lake_max (begg:endg                      ))
    allocate (pct_urban_max(begg:endg, numurbl             ))
    if (use_hillslope) then
       allocate (ncolumns_hillslope  (begg:endg            ))
    endif
    allocate (wt_nat_patch (begg:endg, surfpft_lb:surfpft_ub ))

    ! Read list of Patches and their corresponding parameter values
    ! Independent of model resolution, Needs to stay before surfrd_get_data
    call pftcon%Init()

    ! Read surface dataset and set up subgrid weight arrays
    call surfrd_get_data(begg, endg, ldomain, fsurdat, hillslope_file, actual_numcft)

    if(use_fates) then

       ! Ask Fates to evaluate its own dimensioning needs.
       ! This determines the total amount of space it requires in its largest
       ! dimension.  We are currently calling that the "cohort" dimension, but
       ! it is really a utility dimension that captures the models largest
       ! size need.
       ! Sets:
       !   fates_maxElementsPerPatch
       !   fates_maxElementsPerSite (where a site is roughly equivalent to a column)
       ! (Note: fates_maxELementsPerSite is the critical variable used by CLM
       ! to allocate space)
       ! This also sets up various global constants in FATES
       ! ------------------------------------------------------------------------

       call CLMFatesGlobals2()

    end if

    ! Determine decomposition of subgrid scale landunits, columns, patches
    call decompInit_clumps(ni, nj, glc_behavior)

    ! *** Get ALL processor bounds - for gridcells, landunit, columns and patches ***
    call get_proc_bounds(bounds_proc)

    ! Allocate memory for subgrid data structures
    ! This is needed here BEFORE the following call to initGridcells
    ! Note that the assumption is made that none of the subgrid initialization
    ! can depend on other elements of the subgrid in the calls below
    call grc%Init  (bounds_proc%begg, bounds_proc%endg)
    call lun%Init  (bounds_proc%begl, bounds_proc%endl)
    call col%Init  (bounds_proc%begc, bounds_proc%endc)
    call patch%Init(bounds_proc%begp, bounds_proc%endp)

    ! Build hierarchy and topological info for derived types
    ! This is needed here for the following call to decompInit_glcp
    nclumps = get_proc_clumps()
    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1, nclumps
       call get_clump_bounds(nc, bounds_clump)
       call initGridCells(bounds_clump, glc_behavior)
    end do
    !$OMP END PARALLEL DO

    ! Set global seg maps for gridcells, landlunits, columns and patches
    call decompInit_glcp(ni, nj, glc_behavior)

    if (use_hillslope) then
       ! Initialize hillslope properties
       call InitHillslope(bounds_proc, hillslope_file)
    endif

    ! Set filters
    call allocFilters()

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1, nclumps
       call get_clump_bounds(nc, bounds_clump)
       call reweight_wrapup(bounds_clump, glc_behavior)
    end do
    !$OMP END PARALLEL DO

    ! Set CH4 Model Parameters from namelist.
    ! Need to do before initTimeConst so that it knows whether to
    ! look for several optional parameters on surfdata file.
    if (use_lch4) then
       call ch4conrd()
    end if

    ! Run any requested self-tests
    call self_test_driver(bounds_proc)

    ! Deallocate surface grid dynamic memory for variables that aren't needed elsewhere.
    ! Some things are kept until the end of initialize2; urban_valid is kept through the
    ! end of the run for error checking, pct_urban_max is kept through the end of the run
    ! for reweighting in subgridWeights.
    deallocate (wt_lunit, wt_cft, wt_glc_mec, pct_lake_max)
    if (use_hillslope)  deallocate (ncolumns_hillslope)

    ! Determine processor bounds and clumps for this processor
    call get_proc_bounds(bounds_proc)
    nclumps = get_proc_clumps()

    ! Read in parameters files
    call clm_instReadNML( NLFilename )
    allocate(nutrient_competition_method, &
         source=create_nutrient_competition_method(bounds_proc))
    call readParameters(photosyns_inst)

    
    ! Initialize time manager
    if (nsrest == nsrStartup) then
       call timemgr_init()
    else
       call timemgr_init(curr_date_in=currtime)
       call restFile_getfile(file=fnamer, path=pnamer)
       call restFile_open( flag='read', file=fnamer, ncid=ncid )
       call timemgr_restart_io( ncid=ncid, flag='read' )
       call restFile_close( ncid=ncid )
       call timemgr_restart()
    end if

    ! Pass model timestep info to FATES
    if (use_fates) call CLMFatesTimesteps()

    ! Initialize daylength from the previous time step (needed so prev_dayl can be set correctly)
    call t_startf('init_orbd')
    calday = get_curr_calday(reuse_day_365_for_day_366=.true.)
    call shr_orb_decl( calday, eccen, mvelpp, lambm0, obliqr, declin, eccf )
    dtime = get_step_size_real()
    caldaym1 = get_curr_calday(offset=-int(dtime), reuse_day_365_for_day_366=.true.)
    call shr_orb_decl( caldaym1, eccen, mvelpp, lambm0, obliqr, declinm1, eccf )
    call t_stopf('init_orbd')
    call InitDaylength(bounds_proc, declin=declin, declinm1=declinm1, obliquity=obliqr)

    ! Initialize Balance checking (after time-manager)
    call BalanceCheckInit()

    ! History file variables
    if (use_cn) then
       call hist_addfld1d (fname='DAYL',  units='s', &
            avgflag='A', long_name='daylength', &
            ptr_gcell=grc%dayl, default='inactive')

       call hist_addfld1d (fname='PREV_DAYL', units='s', &
            avgflag='A', long_name='daylength from previous timestep', &
            ptr_gcell=grc%prev_dayl, default='inactive')
    end if

    ! Initialize component data structures
    ! Note: new logic is in place that sets all the history fields to spval so
    ! there is no guesswork in the initialization to nans of the allocated variables
    ! First put in history calls for subgrid data structures - these cannot appear in the
    ! module for the subgrid data definition due to circular dependencies that are introduced

    data2dptr => col%dz(:,-nlevsno+1:0)
    col%dz(bounds_proc%begc:bounds_proc%endc,:) = spval
    call hist_addfld2d (fname='SNO_Z', units='m', type2d='levsno',  &
         avgflag='A', long_name='Snow layer thicknesses', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    call hist_addfld2d (fname='SNO_Z_ICE', units='m', type2d='levsno',  &
         avgflag='A', long_name='Snow layer thicknesses (ice landunits only)', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, &
         l2g_scale_type='ice', default='inactive')

    col%zii(bounds_proc%begc:bounds_proc%endc) = spval
    call hist_addfld1d (fname='ZII', units='m', &
         avgflag='A', long_name='convective boundary height', &
         ptr_col=col%zii, default='inactive')

    ! Initialize instances of all derived types as well as time constant variables
    call clm_instInit(bounds_proc)

    call CNParamsSetSoilDepth()
    ! Initialize SNICAR optical and aging parameters
    call SnowOptics_init( ) ! SNICAR optical parameters:
    call SnowAge_init( )    ! SNICAR aging   parameters:

    ! Print history field info to standard out
    call hist_printflds()

    ! Initializate dynamic subgrid weights (for prescribed transient Patches, CNDV
    ! and/or dynamic landunits); note that these will be overwritten in a restart run
    call t_startf('init_dyn_subgrid')
    call init_subgrid_weights_mod(bounds_proc)
    call dynSubgrid_init(bounds_proc, glc_behavior, crop_inst)
    call t_stopf('init_dyn_subgrid')

    ! Initialize fates LUH2 usage
    if (use_fates_luh) then
       call dynFatesLandUseInit(bounds_proc, fluh_timeseries)
    end if

    ! Initialize baseline water and energy states needed for dynamic subgrid operation
    ! This will be overwritten by the restart file, but needs to be done for a cold start
    ! case.
    ! BACKWARDS_COMPATIBILITY(wjs, 2019-03-05) dyn_hwcontent_set_baselines is called again
    ! later in initialization if reset_dynbal_baselines is set. I think we could just have
    ! a single call in that location (adding some logic to also do the call if we're doing
    ! a cold start) once we can assume that all finidat files have the necessary restart
    ! fields on them. But for now, having the extra call here handles the case where the
    ! relevant restart fields are missing from an old finidat file.
    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1,nclumps
       call get_clump_bounds(nc, bounds_clump)

       call dyn_hwcontent_set_baselines(bounds_clump, &
            filter_inactive_and_active(nc)%num_icec, &
            filter_inactive_and_active(nc)%icec, &
            filter_inactive_and_active(nc)%num_lakec, &
            filter_inactive_and_active(nc)%lakec, &
            urbanparams_inst, soilstate_inst, lakestate_inst, water_inst, temperature_inst, &
            reset_all_baselines = .true., &
            ! reset_lake_baselines is irrelevant since reset_all_baselines is true
            reset_lake_baselines = .false.)
    end do
    !$OMP END PARALLEL DO

    ! Initialize modules (after time-manager initialization in most cases)
    if (use_cn .or. use_fates) then
       call bgc_vegetation_inst%Init2(bounds_proc, NLFilename)
    end if

    if (use_cn) then

       ! NOTE(wjs, 2016-02-23) Maybe the rest of the body of this conditional should also
       ! be moved into bgc_vegetation_inst%Init2
       if (n_drydep > 0) then
          ! Must do this also when drydeposition is used so that estimates of monthly
          ! differences in LAI can be computed
          ! Also do this for FATES see below
          call SatellitePhenologyInit(bounds_proc)
       end if
       if ( use_c14 .and. use_c14_bombspike ) then
          call C14_init_BombSpike()
       end if
       if ( use_c13 .and. use_c13_timeseries ) then
          call C13_init_TimeSeries()
       end if

    else ! FATES OR Satellite phenology

       ! For SP FATES-SP Initialize SP
       ! Also for FATES with Dry-Deposition on as well (see above)
       !if(use_fates_sp .or. (.not.use_cn) .or. (n_drydep > 0) )then  !  Replace with this when we have dry-deposition working
       ! For now don't allow for dry-deposition because of issues in #1044 EBK Jun/17/2022
       if( use_fates_sp .or. .not. use_fates )then
          call SatellitePhenologyInit(bounds_proc)
       end if

       ! fates_spitfire_mode is assigned an integer value in the namelist
       ! see bld/namelist_files/namelist_definition_clm4_5.xml for details
       if(use_fates .and. (fates_spitfire_mode > scalar_lightning)) then
          call clm_fates%Init2(bounds_proc, NLFilename)
       end if
    end if
    if (use_soil_moisture_streams) then
       call PrescribedSoilMoistureInit(bounds_proc)
    endif

    ! On restart only - process the history namelist.
    ! Later the namelist from the restart file will be used.  This allows basic
    ! checking to make sure you didn't try to change the history namelist on restart.
    if (nsrest == nsrContinue ) then
       call htapes_fieldlist()
    end if

    ! Read restart/initial info
    is_cold_start = .false.
    reset_dynbal_baselines_lake_columns = .false.
    if (nsrest == nsrStartup) then
       if (finidat == ' ') then
          if (finidat_interp_source == ' ') then
             is_cold_start = .true.
             if (masterproc) then
                write(iulog,'(a)')'Using cold start initial conditions '
             end if
          else
             if (masterproc) then
                write(iulog,'(a)')'Interpolating initial conditions from '//trim(finidat_interp_source)
                write(iulog,'(a)')'Creating new initial conditions file '//trim(finidat_interp_dest)
             end if
          end if
       else
          if (trim(finidat) == trim(finidat_interp_dest)) then
             ! Check to see if status file for finidat exists
             call check_missing_initdata_status(finidat_interp_dest)
          end if
          if (masterproc) then
             write(iulog,'(a)')'Reading initial conditions from file '//trim(finidat)
          end if
          call getfil( finidat, fnamer, 0 )
          call restFile_read(bounds_proc, fnamer, glc_behavior, &
               reset_dynbal_baselines_lake_columns = reset_dynbal_baselines_lake_columns)
       end if
    else if ((nsrest == nsrContinue) .or. (nsrest == nsrBranch)) then
       if (masterproc) then
          write(iulog,'(a)')'Reading restart file '//trim(fnamer)
       end if
       call restFile_read(bounds_proc, fnamer, glc_behavior, &
            reset_dynbal_baselines_lake_columns = reset_dynbal_baselines_lake_columns)
    end if

    ! If appropriate, create interpolated initial conditions
    if (nsrest == nsrStartup .and. finidat_interp_source /= ' ') then

       ! Check that finidat is not cold start - abort if it is
       if (finidat /= ' ') then
          call endrun(msg='ERROR clm_initializeMod: '//&
               'finidat and finidat_interp_source cannot both be non-blank')
       end if

       ! Determine name if finidat_interp_dest status file
       klen = len_trim(finidat_interp_dest) - 3 ! remove the .nc
       locfn = finidat_interp_dest(1:klen)//'.status'

       ! Remove file if it already exists
       if (masterproc) then
          inquire(file=trim(locfn), exist=lexists)
          if (lexists) then
             open(unit=9876, file=locfn, status='old', iostat=ioe)
             if (ioe == 0) then
                close(9876, status='delete')
             end if
          end if
       end if
       call mpi_barrier(mpicom,ier)

       ! Create new template file using cold start
       call restFile_write(bounds_proc, finidat_interp_dest, writing_finidat_interp_dest_file=.true.)

       ! Interpolate finidat onto new template file
       call getfil( finidat_interp_source, fnamer,  0 )
       call initInterp(filei=fnamer, fileo=finidat_interp_dest, bounds=bounds_proc, &
            glc_behavior=glc_behavior)

       ! Read new interpolated conditions file back in
       call restFile_read(bounds_proc, finidat_interp_dest, glc_behavior, &
            reset_dynbal_baselines_lake_columns = reset_dynbal_baselines_lake_columns)

       ! Reset finidat to now be finidat_interp_dest
       ! (to be compatible with routines still using finidat)
       finidat = trim(finidat_interp_dest)

       ! Write out finidat status flag
       call mpi_barrier(mpicom,ier)
       if (masterproc) then
          open (newunit=iun, file=locfn, status='unknown',  iostat=ioe)
          if (ioe /= 0) then
             call endrun(msg='ERROR failed to open file '//trim(locfn))
          end if
          write(iun,'(a)')'Successfully wrote out '//trim(locfn)
          close(iun)
          write(iulog,'(a)')' Successfully wrote finidat status file '//trim(locfn)
       end if
    end if

    ! If requested, reset dynbal baselines
    ! This needs to happen after reading the restart file (including after reading the
    ! interpolated restart file, if applicable).
    reset_dynbal_baselines_all_columns = get_reset_dynbal_baselines()
    if (nsrest == nsrBranch) then
       if (reset_dynbal_baselines_all_columns) then
          call endrun(msg='ERROR clm_initializeMod: '//&
               'Cannot set reset_dynbal_baselines in a branch run')
       end if
    else if (nsrest == nsrContinue) then
       ! It's okay for the reset_dynbal_baselines flag to remain set in a continue
       ! run, but we'll ignore it. (This way, the user doesn't have to change their
       ! namelist file for the continue run.)
       reset_dynbal_baselines_all_columns = .false.
    end if
    ! Note that we will still honor reset_dynbal_baselines_lake_columns even in a branch
    ! or continue run: even in these runs, we want to reset those baselines if they are
    ! wrong on the restart file.

    if (masterproc) then
       if (reset_dynbal_baselines_all_columns) then
          write(iulog,*) ' '
          write(iulog,*) 'Resetting dynbal baselines for all columns'
          write(iulog,*) ' '
       else if (reset_dynbal_baselines_lake_columns) then
          write(iulog,*) ' '
          write(iulog,*) 'Resetting dynbal baselines for lake columns'
          write(iulog,*) ' '
       end if
    end if

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1,nclumps
       call get_clump_bounds(nc, bounds_clump)

       call dyn_hwcontent_set_baselines(bounds_clump, &
            filter_inactive_and_active(nc)%num_icec, &
            filter_inactive_and_active(nc)%icec, &
            filter_inactive_and_active(nc)%num_lakec, &
            filter_inactive_and_active(nc)%lakec, &
            urbanparams_inst, soilstate_inst, lakestate_inst, &
            water_inst, temperature_inst, &
            reset_all_baselines = reset_dynbal_baselines_all_columns, &
            reset_lake_baselines = reset_dynbal_baselines_lake_columns)
    end do
    !$OMP END PARALLEL DO

    ! Initialize nitrogen deposition
    if (use_cn ) then !.or. use_fates_bgc) then (ndep with fates will be added soon RGK)
       call t_startf('init_ndep')
       if (.not. ndep_from_cpl) then
          call ndep_init(bounds_proc, NLFilename)
          call ndep_interp(bounds_proc, atm2lnd_inst)
       end if
       call t_stopf('init_ndep')
    end if

    ! Initialize crop calendars
    if (use_crop) then
      call t_startf('init_cropcal')
      call cropcal_init(bounds_proc)
      if (use_cropcal_streams) then
        call cropcal_advance( bounds_proc )
        !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
        do nc = 1,nclumps
           call get_clump_bounds(nc, bounds_clump)
           call cropcal_interp(bounds_clump, filter_inactive_and_active(nc)%num_pcropp, &
                filter_inactive_and_active(nc)%pcropp, .true., crop_inst)
        end do
        !$OMP END PARALLEL DO
      end if
      call t_stopf('init_cropcal')
    end if

    ! Initialize active history fields.
    ! This is only done if not a restart run. If a restart run, then this
    ! information has already been obtained from the restart data read above.
    ! Note that routine hist_htapes_build needs time manager information,
    ! so this call must be made after the restart information has been read.
    if (nsrest /= nsrContinue) then
       call hist_htapes_build()
    end if

    ! Initialize variables that are associated with accumulated fields.
    ! The following is called for both initial and restart runs and must
    ! must be called after the restart file is read
    call atm2lnd_inst%initAccVars(bounds_proc)
    call temperature_inst%initAccVars(bounds_proc)
    call water_inst%initAccVars(bounds_proc)
    call energyflux_inst%initAccVars(bounds_proc)
    call canopystate_inst%initAccVars(bounds_proc)
    call bgc_vegetation_inst%initAccVars(bounds_proc)
    if (use_crop) then
       call crop_inst%initAccVars(bounds_proc)
    end if

    if ( use_fates )then
       call clm_fates%initAccVars(bounds_proc)
    end if

    ! Read monthly vegetation
    ! Even if CN or FATES is on, and dry-deposition is active, read CLMSP annual vegetation
    ! to get estimates of monthly LAI
    if ( n_drydep > 0 ) then
       call readAnnualVegetation(bounds_proc, canopystate_inst)
       ! Call interpMonthlyVeg for dry-deposition so that mlaidiff will be calculated
       ! This needs to be done even if FATES, CN or CNDV is on!
       call interpMonthlyVeg(bounds_proc, canopystate_inst)
    ! If fates has satellite phenology enabled, get the monthly veg values
    ! prior to the first call to SatellitePhenology()
    elseif ( use_fates_sp ) then
       call interpMonthlyVeg(bounds_proc, canopystate_inst)
    end if
    
    ! Determine gridcell averaged properties to send to atm
    if (nsrest == nsrStartup) then
       call t_startf('init_map2gc')
       call lnd2atm_minimal(bounds_proc, &
            water_inst, surfalb_inst, energyflux_inst, lnd2atm_inst)
       call t_stopf('init_map2gc')
    end if

    ! Initialize sno export state to send to glc
    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1,nclumps
       call get_clump_bounds(nc, bounds_clump)

       call t_startf('init_lnd2glc')
       call lnd2glc_inst%update_lnd2glc(bounds_clump,       &
            filter(nc)%num_do_smb_c, filter(nc)%do_smb_c,   &
            temperature_inst, water_inst%waterfluxbulk_inst, topo_inst, &
            init=.true.)
       call t_stopf('init_lnd2glc')
    end do
    !$OMP END PARALLEL DO

    ! Deallocate wt_nat_patch
    ! wt_nat_patch was allocated in initialize1, but needed to be kept around through
    ! initialize2 for some consistency checking; now it can be deallocated
    deallocate(wt_nat_patch)

    ! Initialise the fates model state structure
    if ( use_fates .and. .not. (is_restart() .or. nsrest .eq. nsrBranch) .and. finidat == ' ') then
       ! If fates is using satellite phenology mode, make sure to call the SatellitePhenology
       ! procedure prior to init_coldstart which will eventually call leaf_area_profile
       if ( use_fates_sp ) then
          !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
          do nc = 1,nclumps
             call get_clump_bounds(nc, bounds_clump)
             ! FATES satellite phenology mode needs to include all active and inactive patch-level soil
             ! filters due to the translation between the hlm pfts and the fates pfts.
             ! E.g. in FATES, an active PFT vector of 1, 0, 0, 0, 1, 0, 1, 0 would be mapped into
             ! the host land model as 1, 1, 1, 0, 0, 0, 0.  As such, the 'active' filter would only
             ! use the first three points, which would incorrectly represent the interpolated values.
             call CalcSatellitePhenologyTimeInterp(bounds_clump, &
                  filter_inactive_and_active(nc)%num_soilp, filter_inactive_and_active(nc)%soilp, &
                  canopystate_inst)

          end do
          !$OMP END PARALLEL DO
       end if
       
       call clm_fates%init_coldstart(water_inst%waterstatebulk_inst, &
            water_inst%waterdiagnosticbulk_inst, canopystate_inst, &
            soilstate_inst, soilbiogeochem_carbonflux_inst)
    end if
    
    ! topo_glc_mec was allocated in initialize1, but needed to be kept around through
    ! initialize2 because it is used to initialize other variables; now it can be deallocated
    deallocate(topo_glc_mec, fert_cft, irrig_method)

    ! Write log output for end of initialization
    call t_startf('init_wlog')
    if (masterproc) then
       write(iulog,*) 'Successfully initialized the land model'
       if (nsrest == nsrStartup) then
          write(iulog,*) 'begin initial run at: '
       else
          write(iulog,*) 'begin continuation run at:'
       end if
       call get_curr_date(yr, mon, day, ncsec)
       write(iulog,*) '   nstep= ',get_nstep(), ' year= ',yr,' month= ',mon,&
            ' day= ',day,' seconds= ',ncsec
       write(iulog,*)
       write(iulog,'(72a1)') ("*",i=1,60)
       write(iulog,*)
    endif
    call t_stopf('init_wlog')

    if (water_inst%DoConsistencyCheck()) then
       !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
       do nc = 1,nclumps
          call get_clump_bounds(nc, bounds_clump)
          call water_inst%TracerConsistencyCheck(bounds_clump, 'end of initialization')
       end do
       !$OMP END PARALLEL DO
    end if

    call t_stopf('clm_init2')

  end subroutine initialize2

end module clm_initializeMod
