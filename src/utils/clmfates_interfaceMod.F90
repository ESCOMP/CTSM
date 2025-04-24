module CLMFatesInterfaceMod

   ! -------------------------------------------------------------------------------------
   ! This module contains various functions and definitions to aid in the
   ! coupling of the FATES library/API with the CLM/ALM/ATS/etc model driver.
   ! All connections between the two models should occur in this file alone.
   !
   ! This is also the only location where CLM code is allowed to see FATES memory
   ! structures.
   ! The routines here, that call FATES library routines, cannot pass most types defined
   ! by the driving land model (HLM), only native type arrays (int,real,log, etc), implementations
   ! of fates abstract classes, and references into fates boundary condition structures.
   !
   ! Note that CLM/ALM does use Shared Memory Parallelism (SMP), where processes such as
   ! the update of state variables are forked.  However, IO is not assumed to be
   ! threadsafe and therefore memory spaces reserved for IO must be continuous vectors,
   ! and moreover they must be pushed/pulled from history IO for each individual
   ! bounds_proc memory space as a unit.
   !
   ! Therefore, the state variables in the clm_fates communicator is vectorized by
   ! threadcount, and the IO communication arrays are not.
   !
   !
   ! Conventions:
   ! keep line widths within 90 spaces
   ! HLM acronym = Host Land Model
   !
   ! -------------------------------------------------------------------------------------

   !  use ed_driver_interface, only:

   ! Used CLM Modules
#include "shr_assert.h"
   use PatchType         , only : patch
   use shr_kind_mod      , only : r8 => shr_kind_r8
   use decompMod         , only : bounds_type, subgrid_level_column
   use WaterStateBulkType    , only : waterstatebulk_type
   use WaterDiagnosticBulkType    , only : waterdiagnosticbulk_type
   use WaterFluxBulkType     , only : waterfluxbulk_type
   use Wateratm2lndBulkType     , only : wateratm2lndbulk_type
   use ActiveLayerMod    , only : active_layer_type
   use CanopyStateType   , only : canopystate_type
   use TemperatureType   , only : temperature_type
   use EnergyFluxType    , only : energyflux_type
   use SoilStateType     , only : soilstate_type
   use CNProductsMod     , only : cn_products_type
   use clm_varctl        , only : iulog
   use clm_varctl        , only : fates_parteh_mode
   use PRTGenericMod     , only : prt_cnp_flex_allom_hyp
   use clm_varctl        , only : use_fates
   use clm_varctl        , only : fates_spitfire_mode
   use clm_varctl        , only : use_fates_tree_damage
   use clm_varctl        , only : use_fates_planthydro
   use clm_varctl        , only : use_fates_cohort_age_tracking
   use clm_varctl        , only : use_fates_daylength_factor
   use clm_varctl        , only : fates_photosynth_acclimation
   use clm_varctl        , only : use_fates_ed_st3
   use clm_varctl        , only : use_fates_ed_prescribed_phys
   use clm_varctl        , only : fates_harvest_mode
   use clm_varctl        , only : fates_stomatal_model
   use clm_varctl        , only : fates_stomatal_assimilation
   use clm_varctl        , only : fates_leafresp_model
   use clm_varctl        , only : fates_cstarvation_model
   use clm_varctl        , only : fates_regeneration_model
   use clm_varctl        , only : fates_hydro_solver
   use clm_varctl        , only : fates_radiation_model
   use clm_varctl        , only : fates_electron_transport_model
   use clm_varctl        , only : use_fates_inventory_init
   use clm_varctl        , only : use_fates_fixed_biogeog
   use clm_varctl        , only : use_fates_nocomp
   use clm_varctl        , only : use_fates_sp
   use clm_varctl        , only : use_fates_luh
   use clm_varctl        , only : use_fates_lupft
   use clm_varctl        , only : use_fates_potentialveg
   use clm_varctl        , only : flandusepftdat
   use clm_varctl        , only : fates_seeddisp_cadence
   use clm_varctl        , only : fates_inventory_ctrl_filename
   use clm_varctl        , only : use_nitrif_denitrif
   use clm_varctl        , only : use_lch4
   use clm_varctl        , only : fates_history_dimlevel
   use clm_varctl        , only : nsrest, nsrBranch
   use clm_varcon        , only : tfrz
   use clm_varcon        , only : spval
   use clm_varcon        , only : denice
   use clm_varcon        , only : ispval
   use clm_varcon        , only : sum_to_1_tol
   use clm_varpar        , only : surfpft_lb,surfpft_ub
   use clm_varpar        , only : numrad
   use clm_varpar        , only : ivis
   use clm_varpar        , only : inir
   use clm_varpar        , only : nlevdecomp
   use clm_varpar        , only : nlevdecomp_full
   use clm_varpar        , only : nlevsoi
   use PhotosynthesisMod , only : photosyns_type
   use atm2lndType       , only : atm2lnd_type
   use SurfaceAlbedoType , only : surfalb_type
   use SolarAbsorbedType , only : solarabs_type
   use SoilBiogeochemCarbonFluxType, only :  soilbiogeochem_carbonflux_type
   use SoilBiogeochemCarbonStateType, only : soilbiogeochem_carbonstate_type
   use SoilBiogeochemNitrogenFluxType, only :  soilbiogeochem_nitrogenflux_type
   use SoilBiogeochemNitrogenStateType, only : soilbiogeochem_nitrogenstate_type
   use FrictionVelocityMod  , only : frictionvel_type
   use clm_time_manager  , only : is_restart, is_first_restart_step
   use ncdio_pio         , only : file_desc_t, ncd_int, ncd_double
   use restUtilMod,        only : restartvar
   use clm_time_manager  , only : get_curr_days_per_year, &
                                  get_curr_date,     &
                                  get_ref_date,      &
                                  timemgr_datediff,  &
                                  is_beg_curr_day,   &
                                  get_step_size_real,&
                                  get_nstep
   use spmdMod           , only : masterproc
   use decompMod         , only : get_proc_bounds,   &
                                  get_proc_clumps,   &
                                  get_proc_global,   &
                                  get_clump_bounds
   use SoilBiogeochemDecompCascadeConType , only : mimics_decomp, decomp_method
   use SoilBiogeochemDecompCascadeConType , only : no_soil_decomp, century_decomp
   use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
   use GridCellType      , only : grc
   use ColumnType        , only : col
   use LandunitType      , only : lun
   use landunit_varcon   , only : istsoil
   use abortutils        , only : endrun
   use shr_log_mod       , only : errMsg => shr_log_errMsg
   use clm_varcon        , only : dzsoi_decomp
   use FuncPedotransferMod, only: get_ipedof
   use CLMFatesParamInterfaceMod, only: fates_param_reader_ctsm_impl
!   use SoilWaterPlantSinkMod, only : Compute_EffecRootFrac_And_VertTranSink_Default

   ! Used FATES Modules
   use FatesInterfaceMod , only : fates_interface_type
   use FatesInterfaceMod, only : FatesInterfaceInit, FatesReportParameters
   use FatesInterfaceMod, only : SetFatesGlobalElements1
   use FatesInterfaceMod, only : SetFatesGlobalElements2
   use FatesInterfaceMod     , only : allocate_bcin
   use FatesInterfaceMod     , only : allocate_bcout
   use FatesInterfaceMod     , only : allocate_bcpconst
   use FatesInterfaceMod     , only : set_bcpconst
   use FatesInterfaceMod     , only : zero_bcs
   use FatesInterfaceMod     , only : SetFatesTime
   use FatesInterfaceMod     , only : set_fates_ctrlparms
   use FatesInterfaceMod     , only : UpdateFatesRMeansTStep
   use FatesInterfaceMod     , only : InitTimeAveragingGlobals
   
   use FatesParametersInterface, only : fates_param_reader_type
   use FatesParametersInterface, only : fates_parameters_type

   use FatesInterfaceMod     , only : DetermineGridCellNeighbors
   use FatesHistoryInterfaceMod, only : fates_hist
   use FatesRestartInterfaceMod, only : fates_restart_interface_type

   use FatesPatchMod         , only : fates_patch_type
   use PRTGenericMod         , only : num_elements
   use FatesInterfaceTypesMod, only : hlm_stepsize
   use FatesInterfaceTypesMod, only : fates_maxPatchesPerSite
   use FatesInterfaceTypesMod, only : hlm_num_luh2_states
   use FatesInterfaceTypesMod, only : fates_dispersal_cadence_none
   use EDMainMod             , only : ed_ecosystem_dynamics
   use EDMainMod             , only : ed_update_site
   use EDInitMod             , only : zero_site
   use EDInitMod             , only : init_site_vars
   use EDInitMod             , only : init_patches
   use EDInitMod             , only : set_site_properties
   use EDPftVarcon           , only : EDpftvarcon_inst
   use FatesRadiationDriveMod, only : FatesSunShadeFracs, FatesNormalizedCanopyRadiation
   use EDBtranMod            , only : btran_ed, &
                                      get_active_suction_layers
   use EDCanopyStructureMod  , only : canopy_summarization, update_hlm_dynamics
   use EDCanopyStructureMod  , only : UpdateFatesAvgSnowDepth
   use FatesPlantRespPhotosynthMod, only : FatesPlantRespPhotosynthDrive
   use EDAccumulateFluxesMod , only : AccumulateFluxes_ED
   use FatesSoilBGCFluxMod    , only : FluxIntoLitterPools
   use FatesSoilBGCFluxMod    , only : UnPackNutrientAquisitionBCs
   use FatesPlantHydraulicsMod, only : hydraulics_drive
   use FatesPlantHydraulicsMod, only : HydrSiteColdStart
   use FatesPlantHydraulicsMod, only : InitHydrSites
   use FatesPlantHydraulicsMod, only : RestartHydrStates
   use FATESFireBase          , only : fates_fire_base_type
   use FATESFireFactoryMod    , only : no_fire, scalar_lightning, successful_ignitions,&
                                       anthro_ignitions, anthro_suppression
   use dynHarvestMod          , only : num_harvest_inst, harvest_varnames
   use dynHarvestMod          , only : harvest_units, mass_units, unitless_units
   use dynHarvestMod          , only : dynHarvest_interp_resolve_harvesttypes
   use FatesConstantsMod      , only : hlm_harvest_area_fraction
   use FatesConstantsMod      , only : hlm_harvest_carbon
   use FatesDispersalMod      , only : lneighbors, dispersal_type, IsItDispersalTime

   use perf_mod               , only : t_startf, t_stopf

   use dynFATESLandUseChangeMod, only : num_landuse_transition_vars, num_landuse_state_vars
   use dynFATESLandUseChangeMod, only : landuse_transitions, landuse_states
   use dynFATESLandUseChangeMod, only : landuse_transition_varnames, landuse_state_varnames
   use dynFATESLandUseChangeMod, only : num_landuse_harvest_vars
   use dynFATESLandUseChangeMod, only : fates_harvest_no_logging
   use dynFATESLandUseChangeMod, only : fates_harvest_clmlanduse
   use dynFATESLandUseChangeMod, only : fates_harvest_luh_area
   use dynFATESLandUseChangeMod, only : fates_harvest_luh_mass
   use dynFATESLandUseChangeMod, only : landuse_harvest
   use dynFATESLandUseChangeMod, only : landuse_harvest_units
   use dynFATESLandUseChangeMod, only : landuse_harvest_varnames

   implicit none

   type, public :: f2hmap_type

      ! This is the associated column index of each FATES site
      integer, allocatable :: fcolumn (:)

      ! This is the associated site index of any HLM columns
      ! This vector may be sparse, and non-sites have index 0
      integer, allocatable :: hsites  (:)

   end type f2hmap_type


   type, public :: hlm_fates_interface_type

      !      private


      ! See above for descriptions of the sub-types populated
      ! by thread.  This type is somewhat self-explanatory, in that it simply
      ! breaks up memory and process by thread.  Each thread will have its
      ! own list of sites, and boundary conditions for those sites

      type(fates_interface_type), allocatable :: fates (:)

      ! This memory structure is used to map fates sites
      ! into the host model.  Currently, the FATES site
      ! and its column number matching are its only members

      type(f2hmap_type), allocatable  :: f2hmap(:)

      ! fates_restart is the inteface calss for restarting the model
      type(fates_restart_interface_type) :: fates_restart

      ! fates_fire_data_method determines the fire data passed from HLM to FATES
      class(fates_fire_base_type), allocatable :: fates_fire_data_method

      ! Type structure that holds allocatable arrays for mpi-based seed dispersal
      type(dispersal_type) :: fates_seed

   contains

      procedure, public :: init
      procedure, public :: check_hlm_active
      procedure, public :: restart
      procedure, public :: init_coldstart
      procedure, public :: dynamics_driv
      procedure, public :: wrap_sunfrac
      procedure, public :: wrap_btran
      procedure, public :: wrap_photosynthesis
      procedure, public :: wrap_accumulatefluxes
      procedure, public :: prep_canopyfluxes
      procedure, public :: wrap_canopy_radiation
      procedure, public :: wrap_update_hifrq_hist
      procedure, public :: TransferZ0mDisp
      procedure, public :: InterpFileInputs  ! Interpolate inputs from files
      procedure, public :: Init2  ! Initialization after determining subgrid weights
      procedure, public :: InitAccBuffer ! Initialize any accumulation buffers
      procedure, public :: InitAccVars   ! Initialize any accumulation variables
      procedure, public :: UpdateAccVars ! Update any accumulation variables
      procedure, private :: init_history_io
      procedure, private :: wrap_update_hlmfates_dyn
      procedure, private :: init_soil_depths
      procedure, public  :: ComputeRootSoilFlux
      procedure, public  :: wrap_hydraulics_drive
      procedure, public  :: WrapUpdateFatesRmean
      procedure, public  :: wrap_WoodProducts
      procedure, public  :: WrapGlobalSeedDispersal
      procedure, public  :: WrapUpdateFatesSeedInOut
      procedure, public  :: UpdateCLitterFluxes
      procedure, public  :: UpdateNLitterFluxes
   end type hlm_fates_interface_type

   ! hlm_bounds_to_fates_bounds is not currently called outside the interface.
   ! Although there may be good reasons to, I privatized it so that the next
   ! developer will at least question its usage (RGK)
   private :: hlm_bounds_to_fates_bounds

   ! The GetAndSetTime function is used to get the current time from the CLM
   ! time procedures and then set to the fates global time variables during restart,
   ! init_coldstart, and dynamics_driv function calls
   private :: GetAndSetTime

   logical :: debug  = .false.
   logical, private, allocatable :: copy_fates_var(:)  ! .true. -> copy variable from FATES to CTSM in clmfates_interface
                                                       ! .false. -> do not copy in clmfates_interface and use value already in memory

   character(len=*), parameter, private :: sourcefile = &
        __FILE__

   integer, parameter :: num_landuse_pft_vars = 4

   public  :: CLMFatesGlobals1
   public  :: CLMFatesGlobals2
   public  :: CrossRefHistoryFields

 contains

   subroutine CLMFatesGlobals1(surf_numpft,surf_numcft,maxsoil_patches)

     ! --------------------------------------------------------------------------------
     ! This is the first call to fates
     ! We open the fates parameter file. And use that and some info on
     ! namelist variables to determine how many patches need to be allocated
     ! in CTSM
     ! --------------------------------------------------------------------------------

     integer,intent(in)                             :: surf_numpft
     integer,intent(in)                             :: surf_numcft
     integer,intent(out)                            :: maxsoil_patches
     integer                                        :: pass_use_fixed_biogeog
     integer                                        :: pass_use_nocomp
     integer                                        :: pass_use_sp
     integer                                        :: pass_masterproc
     integer                                        :: pass_use_luh2
     logical                                        :: verbose_output
     type(fates_param_reader_ctsm_impl)             :: var_reader
     
     call t_startf('fates_globals1')

     if (use_fates) then

        verbose_output = .false.
        call FatesInterfaceInit(iulog, verbose_output)

        ! Force FATES parameters that are recieve type, to the unset value
        call set_fates_ctrlparms('flush_to_unset')

        ! Send parameters individually
        
        if(use_fates_fixed_biogeog)then
           pass_use_fixed_biogeog = 1
        else
           pass_use_fixed_biogeog = 0
        end if
        call set_fates_ctrlparms('use_fixed_biogeog',ival=pass_use_fixed_biogeog)
        
        if(use_fates_nocomp)then
           pass_use_nocomp = 1
        else
           pass_use_nocomp = 0
        end if
        call set_fates_ctrlparms('use_nocomp',ival=pass_use_nocomp)
        
        if(use_fates_sp)then
           pass_use_sp = 1
        else
           pass_use_sp = 0
        end if
        call set_fates_ctrlparms('use_sp',ival=pass_use_sp)
        
        if(masterproc)then
           pass_masterproc = 1
        else
           pass_masterproc = 0
        end if
        call set_fates_ctrlparms('masterproc',ival=pass_masterproc)

        ! FATES landuse modes
        if(use_fates_luh) then
           pass_use_luh2 = 1
        else
           pass_use_luh2 = 0
        end if
        call set_fates_ctrlparms('use_luh2',ival=pass_use_luh2)
        
     end if


     ! The following call reads in the parameter file
     ! and then uses that to determine the number of patches
     ! FATES requires. We pass that to CLM here
     ! so that it can perform some of its allocations.
     ! During init 2, we will perform more size checks
     ! and allocations on the FATES side, which require
     ! some allocations from CLM (like soil layering)

     call SetFatesGlobalElements1(use_fates,surf_numpft,surf_numcft,var_reader)

     maxsoil_patches = fates_maxPatchesPerSite
     
     call t_stopf('fates_globals1')

     return
   end subroutine CLMFatesGlobals1

   ! ===================================================================================

   subroutine CLMFatesGlobals2()

     ! --------------------------------------------------------------------------------
     ! This is one of the first calls to fates
     ! Used for setting dimensions.  This MUST
     ! be called after NL variables are specified and
     ! after the FATES parameter file has been read in
     ! Aside from setting global dimension info, which
     ! is used in the history file, we also transfer
     ! over the NL variables to FATES global settings.
     ! --------------------------------------------------------------------------------
       
     integer                                        :: pass_vertsoilc
     integer                                        :: pass_ch4
     integer                                        :: pass_spitfire
     integer                                        :: pass_ed_st3
     integer                                        :: pass_num_lu_harvest_cats
     integer                                        :: pass_lu_harvest
     integer                                        :: pass_logging
     integer                                        :: pass_ed_prescribed_phys
     integer                                        :: pass_planthydro
     integer                                        :: pass_inventory_init
     integer                                        :: pass_is_restart
     integer                                        :: pass_cohort_age_tracking
     integer                                        :: pass_tree_damage
     integer                                        :: pass_use_potentialveg
     integer                                        :: pass_num_luh_states
     integer                                        :: pass_num_luh_transitions
     integer                                        :: pass_photosynth_acclimation_switch
     integer                                        :: pass_daylength_factor_switch
     integer                                        :: pass_stomatal_model
     integer                                        :: pass_stomatal_assimilation
     integer                                        :: pass_leafresp_model
     integer                                        :: pass_cstarvation_model
     integer                                        :: pass_regeneration_model
     integer                                        :: pass_hydro_solver
     integer                                        :: pass_radiation_model
     integer                                        :: pass_electron_transport_model

     call t_startf('fates_globals2')

     if (use_fates) then

        

        ! Send parameters individually
        call set_fates_ctrlparms('num_sw_bbands',ival=numrad)
        call set_fates_ctrlparms('vis_sw_index',ival=ivis)
        call set_fates_ctrlparms('nir_sw_index',ival=inir)

        call set_fates_ctrlparms('num_lev_soil',ival=nlevsoi)
        call set_fates_ctrlparms('hlm_name',cval='CLM')
        call set_fates_ctrlparms('hio_ignore_val',rval=spval)
        call set_fates_ctrlparms('soilwater_ipedof',ival=get_ipedof(0))
        
        call set_fates_ctrlparms('parteh_mode',ival=fates_parteh_mode)
        call set_fates_ctrlparms('seeddisp_cadence',ival=fates_seeddisp_cadence)


        call set_fates_ctrlparms('hist_hifrq_dimlevel',ival=fates_history_dimlevel(1))
        call set_fates_ctrlparms('hist_dynam_dimlevel',ival=fates_history_dimlevel(2))
        
        ! CTSM-FATES is not fully coupled (yet)
        ! So lets tell fates to use the RD competition mechanism
        ! which has fewer boundary conditions (simpler)
        call set_fates_ctrlparms('nu_com',cval='RD')

        if (decomp_method == mimics_decomp) then
           call set_fates_ctrlparms('decomp_method',cval='MIMICS')
        elseif(decomp_method == century_decomp ) then
           call set_fates_ctrlparms('decomp_method',cval='CENTURY')
        elseif(decomp_method == no_soil_decomp ) then
           call set_fates_ctrlparms('decomp_method',cval='NONE')
        end if

        if(use_fates_tree_damage)then
           pass_tree_damage = 1
        else
           pass_tree_damage = 0
        end if
        call set_fates_ctrlparms('use_tree_damage',ival=pass_tree_damage)
        
        ! These may be in a non-limiting status (ie when supplements)
        ! are added, but they are always allocated and cycled non-the less
        ! FATES may want to interact differently with other models
        ! that don't even have these arrays allocated.
        ! FATES also checks that if NO3 is cycled in ELM, then
        ! any plant affinity parameters are checked.

        if(use_nitrif_denitrif) then
           call set_fates_ctrlparms('nitrogen_spec',ival=1)
        else
           call set_fates_ctrlparms('nitrogen_spec',ival=2)
        end if

        ! Phosphorus is not tracked in CLM
        call set_fates_ctrlparms('phosphorus_spec',ival=0)


        call set_fates_ctrlparms('spitfire_mode',ival=fates_spitfire_mode)
        call set_fates_ctrlparms('sf_nofire_def',ival=no_fire)
        call set_fates_ctrlparms('sf_scalar_lightning_def',ival=scalar_lightning)
        call set_fates_ctrlparms('sf_successful_ignitions_def',ival=successful_ignitions)
        call set_fates_ctrlparms('sf_anthro_ignitions_def',ival=anthro_ignitions)

        ! This has no variable on the FATES side yet (RGK)
        !call set_fates_ctrlparms('sf_anthro_suppression_def',ival=anthro_suppression)

        if(is_restart() .or. nsrest .eq. nsrBranch) then
           pass_is_restart = 1
        else
           pass_is_restart = 0
        end if
        call set_fates_ctrlparms('is_restart',ival=pass_is_restart)

        if(use_lch4) then
           pass_ch4 = 1
        else
           pass_ch4 = 0
        end if
        call set_fates_ctrlparms('use_ch4',ival=pass_ch4)

        ! use_vertsoilc: Carbon soil layer profile is assumed to be on all the time now
        pass_vertsoilc = 1
        call set_fates_ctrlparms('use_vertsoilc',ival=pass_vertsoilc)

        if(use_fates_ed_st3) then
           pass_ed_st3 = 1
        else
           pass_ed_st3 = 0
        end if
        call set_fates_ctrlparms('use_ed_st3',ival=pass_ed_st3)

        if(use_fates_ed_prescribed_phys) then
           pass_ed_prescribed_phys = 1
        else
           pass_ed_prescribed_phys = 0
        end if
        call set_fates_ctrlparms('use_ed_prescribed_phys',ival=pass_ed_prescribed_phys)

        if(use_fates_planthydro) then
           pass_planthydro = 1
        else
           pass_planthydro = 0
        end if
        call set_fates_ctrlparms('use_planthydro',ival=pass_planthydro)

        if(use_fates_cohort_age_tracking) then
           pass_cohort_age_tracking = 1
        else
           pass_cohort_age_tracking = 0
        end if
        call set_fates_ctrlparms('use_cohort_age_tracking',ival=pass_cohort_age_tracking)

        if(trim(fates_photosynth_acclimation) == 'kumarathunge2019') then
           pass_photosynth_acclimation_switch = 1
        else if(trim(fates_photosynth_acclimation) == 'nonacclimating') then
           pass_photosynth_acclimation_switch = 0
        end if
        call set_fates_ctrlparms('photosynth_acclimation',ival=pass_photosynth_acclimation_switch)

        if(use_fates_daylength_factor) then
           pass_daylength_factor_switch = 1
        else
           pass_daylength_factor_switch = 0
        end if
        call set_fates_ctrlparms('use_daylength_factor_switch',ival=pass_daylength_factor_switch)

        if (trim(fates_stomatal_model) == 'ballberry1987') then
           pass_stomatal_model = 1
        else if (trim(fates_stomatal_model) == 'medlyn2011') then
           pass_stomatal_model = 2
        end if
        call set_fates_ctrlparms('stomatal_model',ival=pass_stomatal_model)

        if (trim(fates_stomatal_assimilation) == 'net') then
           pass_stomatal_assimilation = 1
        else if (trim(fates_stomatal_assimilation) == 'gross') then
           pass_stomatal_assimilation = 2
        end if
        call set_fates_ctrlparms('stomatal_assim_model',ival=pass_stomatal_assimilation)

        if (trim(fates_leafresp_model) == 'ryan1991') then
           pass_leafresp_model = 1
        else if (trim(fates_leafresp_model) == 'atkin2017') then
           pass_leafresp_model = 2
        end if
        call set_fates_ctrlparms('maintresp_leaf_model',ival=pass_leafresp_model)

        if (trim(fates_cstarvation_model) == 'linear') then
           pass_cstarvation_model = 1
        else if (trim(fates_cstarvation_model) == 'expontential') then
           pass_cstarvation_model = 2
        end if
        call set_fates_ctrlparms('mort_cstarvation_model',ival=pass_cstarvation_model)

        if (trim(fates_regeneration_model) == 'default') then
           pass_regeneration_model = 1
        else if (trim(fates_regeneration_model) == 'trs') then
           pass_regeneration_model = 2
        else if (trim(fates_regeneration_model) == 'trs_no_seed_dyn') then
           pass_regeneration_model = 3
        end if
        call set_fates_ctrlparms('regeneration_model',ival=pass_regeneration_model)

        if (trim(fates_hydro_solver) == '1D_Taylor') then
           pass_hydro_solver = 1
        else if (trim(fates_hydro_solver) == '2D_Picard') then
           pass_hydro_solver = 2
        else if (trim(fates_hydro_solver) == '2D_Newton') then
           pass_hydro_solver = 3
        end if
        call set_fates_ctrlparms('hydr_solver',ival=pass_hydro_solver)

        if (trim(fates_radiation_model) == 'norman') then
           pass_radiation_model = 1
        else if (trim(fates_radiation_model) == 'twostream') then
           pass_radiation_model = 2
        end if
        call set_fates_ctrlparms('radiation_model',ival=pass_radiation_model)

        if (trim(fates_electron_transport_model) == 'FvCB1980') then
           pass_radiation_model = 1
        else if (trim(fates_electron_transport_model) == 'JohnsonBerry2021') then
           pass_radiation_model = 2
        end if
        call set_fates_ctrlparms('electron_transport_model',ival=pass_electron_transport_model)

        ! FATES logging and harvest modes
        pass_logging = 0
        pass_lu_harvest = 0
        pass_num_lu_harvest_cats = 0
        if (trim(fates_harvest_mode) /= fates_harvest_no_logging) then
           pass_logging = 1 ! Time driven logging, without landuse harvest
           ! CLM landuse timeseries driven harvest rates
           if (trim(fates_harvest_mode) == fates_harvest_clmlanduse) then
              pass_num_lu_harvest_cats = num_harvest_inst
              pass_lu_harvest = 1

           ! LUH2 landuse timeseries driven  harvest rates
           else if (trim(fates_harvest_mode)== fates_harvest_luh_area .or. &
                    trim(fates_harvest_mode)== fates_harvest_luh_mass) then
              pass_lu_harvest = 1
              pass_num_lu_harvest_cats = num_landuse_harvest_vars
           end if
        end if

        call set_fates_ctrlparms('use_lu_harvest',ival=pass_lu_harvest)
        call set_fates_ctrlparms('num_lu_harvest_cats',ival=pass_num_lu_harvest_cats)
        call set_fates_ctrlparms('use_logging',ival=pass_logging)

        ! FATES landuse modes
        if(use_fates_luh) then
           pass_num_luh_states = num_landuse_state_vars
           pass_num_luh_transitions = num_landuse_transition_vars
        else
           pass_num_luh_states = 0
           pass_num_luh_transitions = 0
        end if
        call set_fates_ctrlparms('num_luh2_states',ival=pass_num_luh_states)
        call set_fates_ctrlparms('num_luh2_transitions',ival=pass_num_luh_transitions)

        if ( use_fates_potentialveg ) then
           pass_use_potentialveg = 1
        else
           pass_use_potentialveg = 0
        end if
        call set_fates_ctrlparms('use_fates_potentialveg',ival=pass_use_potentialveg)

        if(use_fates_inventory_init) then
           pass_inventory_init = 1
        else
           pass_inventory_init = 0
        end if
        call set_fates_ctrlparms('use_inventory_init',ival=pass_inventory_init)

        call set_fates_ctrlparms('inventory_ctrl_file',cval=fates_inventory_ctrl_filename)


        ! Check through FATES parameters to see if all have been set
        call set_fates_ctrlparms('check_allset')

     end if

     ! This determines the total amount of space it requires in its largest
     ! dimension.  We are currently calling that the "cohort" dimension, but
     ! it is really a utility dimension that captures the models largest
     ! size need.
     ! Sets:
     ! fates_maxElementsPerPatch
     ! num_elements
     ! fates_maxElementsPerSite (where a site is roughly equivalent to a column)
     ! (Note: this needs to be called when use_fates=.false. as well, becuase
     ! it will return some nominal dimension sizes of 1

     call SetFatesGlobalElements2(use_fates)

     call t_stopf('fates_globals2')

     return
   end subroutine CLMFatesGlobals2

   ! ===================================================================================

   subroutine CrossRefHistoryFields

     ! This routine only needs to be called on the masterproc.
     ! Here we cross reference the CLM history master
     ! list and make sure that all fields that start
     ! with fates have been allocated. If it has
     ! not, then we give a more constructive error
     ! message than what is possible in PIO. The user
     ! most likely needs to increase the history density
     ! level

     use histFileMod, only: getname
     use histFileMod, only: hist_fincl1,hist_fincl2,hist_fincl3,hist_fincl4
     use histFileMod, only: hist_fincl5,hist_fincl6,hist_fincl7,hist_fincl8
     use histFileMod, only: hist_fincl9,hist_fincl10
     use histFileMod, only: max_tapes, max_flds, max_namlen

     integer :: t     ! iterator index for history tapes
     integer :: f     ! iterator index for registered history field names
     integer :: nh    ! iterator index for fates registered history
     logical :: is_fates_field ! Does this start with FATES_ ?
     logical :: found ! if true, than the history field is either
                      ! not part of the fates set, or was found in
                      ! the fates set
     character(len=64) :: fincl_name
     ! This is a copy of the public in histFileMod, copied
     ! here because it isn't filled at the time of this call
     character(len=max_namlen+2) :: fincl(max_flds,max_tapes)
     
     fincl(:,1)  = hist_fincl1(:)
     fincl(:,2)  = hist_fincl2(:)
     fincl(:,3)  = hist_fincl3(:)
     fincl(:,4)  = hist_fincl4(:)
     fincl(:,5)  = hist_fincl5(:)
     fincl(:,6)  = hist_fincl6(:)
     fincl(:,7)  = hist_fincl7(:)
     fincl(:,8)  = hist_fincl8(:)
     fincl(:,9)  = hist_fincl9(:)
     fincl(:,10) = hist_fincl10(:)
     
     do t = 1,max_tapes

        f = 1
        search_fields: do while (f < max_flds .and. fincl(f,t) /= ' ')
           
           fincl_name = getname(fincl(f,t))
           is_fates_field = fincl_name(1:6)=='FATES_'
           
           if(is_fates_field) then
              found = .false.
              do_fates_hist: do nh = 1,fates_hist%num_history_vars()
                 if(trim(fates_hist%hvars(nh)%vname) == &
                      trim(fincl_name)) then
                    found=.true.
                    exit do_fates_hist
                 end if
              end do do_fates_hist
              
              if(.not.found)then
                 write(iulog,*) 'the history field: ',trim(fincl_name)
                 write(iulog,*) 'was requested in the namelist, but was'
                 write(iulog,*) 'not found in the list of fates_hist%hvars.'
                 write(iulog,*) 'Most likely, this is because this history variable'
                 write(iulog,*) 'was specified in the user namelist, but the user'
                 write(iulog,*) 'specified a FATES history output dimension level'
                 write(iulog,*) 'that does not contain that variable in its valid set.'
                 write(iulog,*) 'You may have to increase the namelist setting: fates_history_dimlevel'
                 write(iulog,*) 'current fates_history_dimlevel: ',fates_history_dimlevel(:)
                 !uncomment if you want to list all fates history variables in registry
                 !do_fates_hist2: do nh = 1,fates_hist%num_history_vars()
                 !   write(iulog,*) trim(fates_hist%hvars(nh)%vname)
                 !end do do_fates_hist2
                 call endrun(msg=errMsg(sourcefile, __LINE__))
              end if
           end if
           f = f + 1
        end do search_fields
        
     end do
   end subroutine CrossRefHistoryFields

   
   ! ===================================================================================
  
   subroutine CLMFatesTimesteps()
     
     hlm_stepsize = get_step_size_real()

     call InitTimeAveragingGlobals()

     return
   end subroutine CLMFatesTimesteps
   
   
   ! ====================================================================================

   subroutine init(this, bounds_proc, flandusepftdat)

      ! ---------------------------------------------------------------------------------
      ! This initializes the hlm_fates_interface_type
      !
      ! sites is the root of the fates state hierarchy (instantaneous info on
      ! the state of the ecosystem).  As such, it governs the connection points between
      ! the host (which also dictates its allocation) and its patch structures.
      !
      ! sites may associate with different scales in different models. In
      ! CLM, it is being designed to relate to column scale.
      !
      ! This global may become relegated to this module.
      !
      ! Note: CLM/ALM currently wants sites to be allocated even if ed
      ! is not turned on
      ! ---------------------------------------------------------------------------------

      use spmdMod,                only : npes
      use decompMod,              only : procinfo
      use FatesInterfaceTypesMod, only : numpft_fates => numpft
      use FatesParameterDerivedMod, only : param_derived
      use subgridMod, only :  natveg_patch_exists
      use clm_instur       , only : wt_nat_patch
      use FATESFireFactoryMod , only: create_fates_fire_data_method

      ! Input Arguments
      class(hlm_fates_interface_type), intent(inout) :: this
      type(bounds_type),intent(in)                   :: bounds_proc
      character(len=*), intent(in)                   :: flandusepftdat

      ! local variables
      integer                                        :: nclumps   ! Number of threads
      integer                                        :: nc        ! thread index
      integer                                        :: s         ! FATES site index
      integer                                        :: c         ! HLM column index
      integer                                        :: l         ! HLM LU index
      integer                                        :: g         ! HLM grid index
      integer                                        :: m         ! HLM PFT index
      integer                                        :: ft       ! FATES PFT index
      integer                                        :: pi,pf
      integer, allocatable                           :: collist (:)
      type(bounds_type)                              :: bounds_clump
      integer                                        :: nmaxcol
      integer                                        :: ndecomp
      integer                                        :: numg

      real(r8), allocatable :: landuse_pft_map(:,:,:)
      real(r8), allocatable :: landuse_bareground(:)

      ! Initialize the FATES communicators with the HLM
      ! This involves to stages
      ! 1) allocate the vectors
      ! 2) add the history variables defined in clm_inst to the history machinery


      call t_startf('fates_init')

      ! Parameter Routines
      call param_derived%Init( numpft_fates )


      ! Initialize dispersal
      if (fates_seeddisp_cadence /= fates_dispersal_cadence_none) then

         ! Initialize fates global seed dispersal array for all nodes
         call get_proc_global(ng=numg)
         call this%fates_seed%init(npes,numg,procinfo%ncells,numpft_fates)

         ! Initialize the array of nearest neighbors for fates-driven grid cell communications
         ! This must be called after surfrd_get_data and decompInit_lnd
         call DetermineGridCellNeighbors(lneighbors,this%fates_seed,numg)
      end if

      nclumps = get_proc_clumps()
      allocate(this%fates(nclumps))
      allocate(this%f2hmap(nclumps))

      if(debug)then
         write(iulog,*) 'clm_fates%init():  allocating for ',nclumps,' threads'
      end if

      ! Retrieve the landuse x pft static data if the file is present
      if (use_fates_fixed_biogeog .and. use_fates_luh) then
         call GetLandusePFTData(bounds_proc, flandusepftdat, landuse_pft_map, landuse_bareground)
      end if

      nclumps = get_proc_clumps()

      allocate(copy_fates_var(bounds_proc%begc:bounds_proc%endc))
      copy_fates_var(:) = .false.

      !$OMP PARALLEL DO PRIVATE (nc,bounds_clump,nmaxcol,s,c,l,g,collist,pi,pf,ft)
      do nc = 1,nclumps
         call get_clump_bounds(nc, bounds_clump)

         nmaxcol = bounds_clump%endc - bounds_clump%begc + 1

         allocate(collist(1:nmaxcol))

         ! Allocate the mapping that points columns to FATES sites, 0 is NA
         allocate(this%f2hmap(nc)%hsites(bounds_clump%begc:bounds_clump%endc))

         ! Initialize all columns with a zero index, which indicates no FATES site
         this%f2hmap(nc)%hsites(:) = 0

         s = 0
         do c = bounds_clump%begc,bounds_clump%endc
            l = col%landunit(c)

            ! These are the key constraints that determine if this column
            ! will have a FATES site associated with it

            ! INTERF-TODO: WE HAVE NOT FILTERED OUT FATES SITES ON INACTIVE COLUMNS.. YET
            ! NEED A RUN-TIME ROUTINE THAT CLEARS AND REWRITES THE SITE LIST

            if (lun%itype(l) == istsoil ) then
               s = s + 1
               collist(s) = c
               this%f2hmap(nc)%hsites(c) = s
               col%is_fates(c) = .true.
               if(debug)then
                  write(iulog,*) 'clm_fates%init(): thread',nc,': found column',c,'with lu',l
                  write(iulog,*) 'LU type:', lun%itype(l)
               end if
            endif

         enddo

         if(debug)then
            write(iulog,*) 'clm_fates%init(): thread',nc,': allocated ',s,' sites'
         end if

         ! Allocate vectors that match FATES sites with HLM columns
         ! RGK: Sites and fcolumns are forced as args during clm_driv() as of 6/4/2016
         ! We may have to give these a dummy allocation of 1, which should
         ! not be a problem since we always iterate on nsites.

         allocate(this%f2hmap(nc)%fcolumn(s))

         ! Assign the h2hmap indexing
         this%f2hmap(nc)%fcolumn(1:s)         =  collist(1:s)

         ! Deallocate the temporary arrays
         deallocate(collist)

         ! Set the number of FATES sites
         this%fates(nc)%nsites = s

         ! Allocate the FATES sites
         allocate (this%fates(nc)%sites(this%fates(nc)%nsites))

         ! Allocate the FATES boundary arrays (in)
         allocate(this%fates(nc)%bc_in(this%fates(nc)%nsites))

         ! Allocate the FATES boundary arrays (out)
         allocate(this%fates(nc)%bc_out(this%fates(nc)%nsites))


         ! Parameter Constants defined by FATES, but used in ELM
         ! Note that FATES has its parameters defined, so we can also set the values
         call allocate_bcpconst(this%fates(nc)%bc_pconst,nlevdecomp)

         ! This also needs
         call set_bcpconst(this%fates(nc)%bc_pconst,nlevdecomp)


         ! Allocate and Initialize the Boundary Condition Arrays
         ! These are staticaly allocated at maximums, so
         ! No information about the patch or cohort structure is needed at this step


         do s = 1, this%fates(nc)%nsites

            c = this%f2hmap(nc)%fcolumn(s)

            this%fates(nc)%sites(s)%h_gid = c

            ndecomp = col%nbedrock(c)

            call allocate_bcin(this%fates(nc)%bc_in(s),col%nbedrock(c),ndecomp, &
                               num_harvest_inst, num_landuse_state_vars, num_landuse_transition_vars, &
                               surfpft_lb,surfpft_ub)
            call allocate_bcout(this%fates(nc)%bc_out(s),col%nbedrock(c),ndecomp)
            call zero_bcs(this%fates(nc),s)

            ! Pass any grid-cell derived attributes to the site
            ! ---------------------------------------------------------------------------
            c = this%f2hmap(nc)%fcolumn(s)
            g = col%gridcell(c)
            this%fates(nc)%sites(s)%lat = grc%latdeg(g)
            this%fates(nc)%sites(s)%lon = grc%londeg(g)

            ! Transfer the landuse x pft data to fates via bc_in if file is given
            if (use_fates_fixed_biogeog) then
               if (use_fates_luh) then
                  this%fates(nc)%bc_in(s)%pft_areafrac_lu(:,1:num_landuse_pft_vars) = landuse_pft_map(g,:,1:num_landuse_pft_vars)
                  this%fates(nc)%bc_in(s)%baregroundfrac = landuse_bareground(g)
               else
                  ! initialize static layers for reduced complexity FATES versions from HLM
                  ! maybe make this into a subroutine of it's own later.
                  this%fates(nc)%bc_in(s)%pft_areafrac(:)=0._r8
                  do m = surfpft_lb,surfpft_ub
                     ft = m - surfpft_lb
                     this%fates(nc)%bc_in(s)%pft_areafrac(ft)=wt_nat_patch(g,m)
                  end do

                  if (abs(sum(this%fates(nc)%bc_in(s)%pft_areafrac(surfpft_lb:surfpft_ub)) - 1.0_r8) >    sum_to_1_tol) then
                     write(iulog,*) 'pft_area error in interfc ', s, sum(this%fates(nc)%bc_in(s)   %pft_areafrac(:)) - 1.0_r8
                     call endrun(msg=errMsg(sourcefile, __LINE__))
                  end if
               end if
            end if
          end do !site

        ! Initialize site-level static quantities dictated by the HLM
        ! currently ground layering depth
         call this%init_soil_depths(nc)

         if (use_fates_planthydro) then
            call InitHydrSites(this%fates(nc)%sites,this%fates(nc)%bc_in)
         end if


         if( this%fates(nc)%nsites == 0 ) then
            write(iulog,*) 'Clump ',nc,' had no valid FATES sites'
            write(iulog,*) 'This will likely cause problems until code is improved'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if


         ! Set patch itypes on natural veg columns to nonsense
         ! This will force a crash if the model outside of FATES tries to think
         ! of the patch as a PFT.

         do s = 1, this%fates(nc)%nsites
            c = this%f2hmap(nc)%fcolumn(s)
            pi = col%patchi(c)+1
            pf = col%patchf(c)
!            patch%itype(pi:pf) = ispval
            patch%is_fates(pi:pf) = .true.
         end do

      end do
      !$OMP END PARALLEL DO

      call this%init_history_io(bounds_proc)

      ! Report Fates Parameters (debug flag in lower level routines)
      call FatesReportParameters(masterproc)

      ! Fire data to send to FATES
      call create_fates_fire_data_method( this%fates_fire_data_method )

      ! deallocate the local landuse x pft array
      if (use_fates_fixed_biogeog .and. use_fates_luh) then
         deallocate(landuse_pft_map)
         deallocate(landuse_bareground)
      end if

      call t_stopf('fates_init')

    end subroutine init

    ! ===================================================================================

    subroutine check_hlm_active(this, nc, bounds_clump)

      ! ---------------------------------------------------------------------------------
      ! This subroutine is not currently used.  It is just a utility that may come
      ! in handy when we have dynamic sites in FATES
      ! ---------------------------------------------------------------------------------

      class(hlm_fates_interface_type), intent(inout) :: this
      integer                                        :: nc
      type(bounds_type),intent(in)                   :: bounds_clump

      ! local variables
      integer :: c

      call t_startf('fates_check_hlm_active')

      do c = bounds_clump%begc,bounds_clump%endc

         ! FATES ACTIVE BUT HLM IS NOT
         if(this%f2hmap(nc)%hsites(c)>0 .and. .not.col%active(c)) then

            write(iulog,*) 'INACTIVE COLUMN WITH ACTIVE FATES SITE'
            write(iulog,*) 'c = ',c
            call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, &
                 msg=errMsg(sourcefile, __LINE__))

         elseif (this%f2hmap(nc)%hsites(c)==0 .and. col%active(c)) then

            write(iulog,*) 'ACTIVE COLUMN WITH INACTIVE FATES SITE'
            write(iulog,*) 'c = ',c
            call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, &
                 msg=errMsg(sourcefile, __LINE__))
         end if
      end do

      call t_stopf('fates_check_hlm_active')

   end subroutine check_hlm_active

   ! ------------------------------------------------------------------------------------

   subroutine dynamics_driv(this, nc, bounds_clump,      &
         atm2lnd_inst, soilstate_inst, temperature_inst, active_layer_inst, &
         waterstatebulk_inst, waterdiagnosticbulk_inst, wateratm2lndbulk_inst, &
         canopystate_inst, soilbiogeochem_carbonflux_inst, frictionvel_inst, &
         soil_water_retention_curve)

      ! This wrapper is called daily from clm_driver
      ! This wrapper calls ed_driver, which is the daily dynamics component of FATES
      ! ed_driver is not a hlm_fates_inst_type procedure because we need an extra step
      ! to process array bounding information

      ! !USES
      use FATESFireFactoryMod, only: scalar_lightning, anthro_ignitions, anthro_suppression
      use subgridMod, only :  natveg_patch_exists

      ! !ARGUMENTS:

      class(hlm_fates_interface_type), intent(inout) :: this
      type(bounds_type),intent(in)                   :: bounds_clump
      type(atm2lnd_type)      , intent(in)           :: atm2lnd_inst
      type(soilstate_type)    , intent(in)           :: soilstate_inst
      type(temperature_type)  , intent(in)           :: temperature_inst
      type(active_layer_type) , intent(in)           :: active_layer_inst
      integer                 , intent(in)           :: nc
      type(waterstatebulk_type)   , intent(inout)        :: waterstatebulk_inst
      type(waterdiagnosticbulk_type)   , intent(inout)        :: waterdiagnosticbulk_inst
      type(wateratm2lndbulk_type)   , intent(inout)        :: wateratm2lndbulk_inst
      type(canopystate_type)  , intent(inout)        :: canopystate_inst
      type(soilbiogeochem_carbonflux_type), intent(inout) :: soilbiogeochem_carbonflux_inst
      type(frictionvel_type)  , intent(inout)        :: frictionvel_inst
      class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve

      ! !LOCAL VARIABLES:
      integer  :: s                        ! site index
      integer  :: g                        ! grid-cell index (HLM)
      integer  :: c                        ! column index (HLM)
      integer  :: j                        ! Soil layer index
      integer  :: ifp                      ! patch index ft
      integer  :: p                        ! HLM patch index
      integer  :: nlevsoil                 ! number of soil layers at the site
      integer  :: nld_si                   ! site specific number of decomposition layers
      integer  :: ft                       ! plant functional type
      real(r8), pointer :: lnfm24(:)       ! 24-hour averaged lightning data
      real(r8), pointer :: gdp_lf_col(:)          ! gdp data
      integer  :: ier
      integer  :: begg,endg
      real(r8) :: harvest_rates(bounds_clump%begg:bounds_clump%endg,num_harvest_inst)
      real(r8) :: s_node, smp_node         ! local for relative water content and potential
      logical  :: after_start_of_harvest_ts
      integer  :: iharv
      !-----------------------------------------------------------------------

      ! ---------------------------------------------------------------------------------
      ! Part I.
      ! Prepare input boundary conditions for FATES dynamics
      ! Note that timing information is the same across all sites, this may
      ! seem redundant, but it is possible that we may have asynchronous site simulations
      ! one day.  The cost of holding site level boundary conditions is minimal
      ! and it keeps all the boundaries in one location
      ! ---------------------------------------------------------------------------------


      call t_startf('fates_dynamics_daily_driver')

      begg = bounds_clump%begg; endg = bounds_clump%endg

      ! Set the FATES global time and date variables
      call GetAndSetTime

      ! Get harvest rates for CLM landuse timeseries driven rates
      if (trim(fates_harvest_mode) == fates_harvest_clmlanduse) then
         call dynHarvest_interp_resolve_harvesttypes(bounds_clump, &
              harvest_rates=harvest_rates(begg:endg,1:num_harvest_inst), &
              after_start_of_harvest_ts=after_start_of_harvest_ts)
      endif

      if (fates_spitfire_mode > scalar_lightning) then
         allocate(lnfm24(bounds_clump%begg:bounds_clump%endg), stat=ier)
         if (ier /= 0) then
            call endrun(msg="allocation error for lnfm24"//&
                 errmsg(sourcefile, __LINE__))
         endif
         lnfm24 = this%fates_fire_data_method%GetLight24()
      end if
      
      if (fates_spitfire_mode == anthro_suppression) then
         allocate(gdp_lf_col(bounds_clump%begc:bounds_clump%endc), stat=ier)
         if (ier /= 0) then
            call endrun(msg="allocation error for gdp"//&
                 errmsg(sourcefile, __LINE__))
         endif
         gdp_lf_col = this%fates_fire_data_method%GetGDP()
      end if

      do s=1,this%fates(nc)%nsites
         c = this%f2hmap(nc)%fcolumn(s)
         g = col%gridcell(c)

         if (fates_spitfire_mode > scalar_lightning) then
            do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno
               
               this%fates(nc)%bc_in(s)%lightning24(ifp) = lnfm24(g) * 24._r8  ! #/km2/hr to #/km2/day
               
               if (fates_spitfire_mode .ge. anthro_ignitions) then
                  this%fates(nc)%bc_in(s)%pop_density(ifp) = this%fates_fire_data_method%forc_hdm(g)
               end if

            end do ! ifp

            if (fates_spitfire_mode == anthro_suppression) then
               ! Placeholder for future fates use of gdp - comment out before integration
               !this%fates(nc)%bc_in(s)%gdp = gdp_lf_col(c) ! k US$/capita(g)
            end if
         end if

         nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil

         ! Decomposition fluxes
         if ( decomp_method /= no_soil_decomp )then
            this%fates(nc)%bc_in(s)%w_scalar_sisl(1:nlevsoil) = soilbiogeochem_carbonflux_inst%w_scalar_col(c,1:nlevsoil)
            this%fates(nc)%bc_in(s)%t_scalar_sisl(1:nlevsoil) = soilbiogeochem_carbonflux_inst%t_scalar_col(c,1:nlevsoil)
         else
            this%fates(nc)%bc_in(s)%w_scalar_sisl(1:nlevsoil) = 0.0_r8
            this%fates(nc)%bc_in(s)%t_scalar_sisl(1:nlevsoil) = 0.0_r8
         end if

         ! Soil water
         this%fates(nc)%bc_in(s)%h2o_liqvol_sl(1:nlevsoil)  = &
               waterstatebulk_inst%h2osoi_vol_col(c,1:nlevsoil)

         this%fates(nc)%bc_in(s)%max_rooting_depth_index_col = &
              min(nlevsoil, active_layer_inst%altmax_lastyear_indx_col(c))

         nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil
         do j = 1,nlevsoil
            this%fates(nc)%bc_in(s)%tempk_sl(j) = temperature_inst%t_soisno_col(c,j)
         end do

         call get_active_suction_layers(this%fates(nc)%nsites, &
             this%fates(nc)%sites,  &
             this%fates(nc)%bc_in,  &
             this%fates(nc)%bc_out)

         do j = 1,nlevsoil
            if(this%fates(nc)%bc_out(s)%active_suction_sl(j)) then
               s_node = max(waterstatebulk_inst%h2osoi_vol_col(c,j)/soilstate_inst%eff_porosity_col(c,j) ,0.01_r8)
               call soil_water_retention_curve%soil_suction(c,j,s_node, soilstate_inst, smp_node)
               this%fates(nc)%bc_in(s)%smp_sl(j)           = smp_node
            end if
         end do


         do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno !for vegetated patches
            ! Mapping between  IFP space (1,2,3) and HLM P space (looping by IFP)
            p = ifp+col%patchi(c)

            this%fates(nc)%bc_in(s)%precip24_pa(ifp) = &
                  wateratm2lndbulk_inst%prec24_patch(p)

            this%fates(nc)%bc_in(s)%relhumid24_pa(ifp) = &
                  wateratm2lndbulk_inst%rh24_patch(p)

            this%fates(nc)%bc_in(s)%wind24_pa(ifp) = &
                  atm2lnd_inst%wind24_patch(p)

         end do

         ! Here we use the same logic as the pft_areafrac initialization to get an array with values for each pft
         ! in FATES.
         ! N.B. Fow now these are fixed values pending HLM updates.
         if(use_fates_sp)then
           do ft = surfpft_lb,surfpft_ub
               ! here we are mapping from P space in the HLM to FT space in the sp_input arrays.
               p = ft + col%patchi(c) ! for an FT of 1 we want to use
               this%fates(nc)%bc_in(s)%hlm_sp_tlai(ft) = canopystate_inst%tlai_input_patch(p)
               this%fates(nc)%bc_in(s)%hlm_sp_tsai(ft) = canopystate_inst%tsai_input_patch(p)
               this%fates(nc)%bc_in(s)%hlm_sp_htop(ft) = canopystate_inst%htop_input_patch(p)
               if(canopystate_inst%htop_input_patch(p).lt.1.0e-20)then ! zero htop causes inifinite/nans. This is
                 this%fates(nc)%bc_in(s)%hlm_sp_htop(ft) = 0.01_r8
               endif
           end do ! p
         end if ! SP

         if(use_fates_planthydro)then
            this%fates(nc)%bc_in(s)%hksat_sisl(1:nlevsoil)  = soilstate_inst%hksat_col(c,1:nlevsoil)
            this%fates(nc)%bc_in(s)%watsat_sisl(1:nlevsoil) = soilstate_inst%watsat_col(c,1:nlevsoil)
            this%fates(nc)%bc_in(s)%watres_sisl(1:nlevsoil) = soilstate_inst%watres_col(c,1:nlevsoil)
            this%fates(nc)%bc_in(s)%sucsat_sisl(1:nlevsoil) = soilstate_inst%sucsat_col(c,1:nlevsoil)
            this%fates(nc)%bc_in(s)%bsw_sisl(1:nlevsoil)    = soilstate_inst%bsw_col(c,1:nlevsoil)
            this%fates(nc)%bc_in(s)%h2o_liq_sisl(1:nlevsoil) =  waterstatebulk_inst%h2osoi_liq_col(c,1:nlevsoil)
         end if

         ! get the harvest data, which is by gridcell
         ! for now there is one veg column per gridcell, so store all harvest data in each site
         ! this will eventually change
         ! today's hlm harvest flag needs to be set no matter what
         if (trim(fates_harvest_mode) == fates_harvest_clmlanduse) then
            if (after_start_of_harvest_ts) then
               this%fates(nc)%bc_in(s)%hlm_harvest_rates(1:num_harvest_inst) = harvest_rates(g,1:num_harvest_inst)
            else
               this%fates(nc)%bc_in(s)%hlm_harvest_rates(1:num_harvest_inst) = 0._r8
            end if
            this%fates(nc)%bc_in(s)%hlm_harvest_catnames(1:num_harvest_inst) = harvest_varnames(1:num_harvest_inst)

            ! also pass the units that the harvest rates are specified in
            if (trim(harvest_units) == trim(unitless_units)) then
               this%fates(nc)%bc_in(s)%hlm_harvest_units = hlm_harvest_area_fraction
            else if (trim(harvest_units) == trim(mass_units)) then
               this%fates(nc)%bc_in(s)%hlm_harvest_units = hlm_harvest_carbon
            else
               write(iulog,*) 'units field not one of the specified options.'
               write(iulog,*) harvest_units
               call endrun(msg=errMsg(sourcefile, __LINE__))
           end if

         else if (trim(fates_harvest_mode) == fates_harvest_luh_area .or. &
                  trim(fates_harvest_mode) == fates_harvest_luh_mass) then
              this%fates(nc)%bc_in(s)%hlm_harvest_rates = landuse_harvest(:,g)
              this%fates(nc)%bc_in(s)%hlm_harvest_catnames = landuse_harvest_varnames
              this%fates(nc)%bc_in(s)%hlm_harvest_units = landuse_harvest_units
         endif

         if (use_fates_luh) then
               this%fates(nc)%bc_in(s)%hlm_luh_states = landuse_states(:,g)
               this%fates(nc)%bc_in(s)%hlm_luh_state_names = landuse_state_varnames
               this%fates(nc)%bc_in(s)%hlm_luh_transitions = landuse_transitions(:,g)
               this%fates(nc)%bc_in(s)%hlm_luh_transition_names = landuse_transition_varnames
         end if

      end do

      ! Nutrient uptake fluxes have been accumulating with each short
      ! timestep, here, we unload them from the boundary condition
      ! structures into the cohort structures.
      call UnPackNutrientAquisitionBCs(this%fates(nc)%sites, this%fates(nc)%bc_in)

      ! Distribute any seeds from neighboring gridcells into the current gridcell
      ! Global seed availability array populated by WrapGlobalSeedDispersal call
      if (fates_seeddisp_cadence /= fates_dispersal_cadence_none) then
         call this%WrapUpdateFatesSeedInOut(bounds_clump)
      end if

      ! ---------------------------------------------------------------------------------
      ! Part II: Call the FATES model now that input boundary conditions have been
      ! provided.
      ! ---------------------------------------------------------------------------------

      do s = 1,this%fates(nc)%nsites

            call ed_ecosystem_dynamics(this%fates(nc)%sites(s),    &
                  this%fates(nc)%bc_in(s), &
                  this%fates(nc)%bc_out(s))

            call ed_update_site(this%fates(nc)%sites(s), &
                  this%fates(nc)%bc_in(s), &
                  this%fates(nc)%bc_out(s), &
                  is_restarting = .false.)
      enddo



      ! ---------------------------------------------------------------------------------
      ! Part III.2 (continued).
      ! Update diagnostics of the FATES ecosystem structure that are used in the HLM.
      ! ---------------------------------------------------------------------------------
      call this%wrap_update_hlmfates_dyn(nc,               &
                                         bounds_clump,     &
                                         waterdiagnosticbulk_inst,  &
                                         canopystate_inst, &
                                         soilbiogeochem_carbonflux_inst, &
                                         .false.)

      ! ---------------------------------------------------------------------------------
      ! Part IV:
      ! Update history IO fields that depend on ecosystem dynamics
      ! ---------------------------------------------------------------------------------
      call fates_hist%update_history_dyn( nc,                    &
                                          this%fates(nc)%nsites, &
                                          this%fates(nc)%sites,  &
                                          this%fates(nc)%bc_in )

      call t_stopf('fates_dynamics_daily_driver')

      return
   end subroutine dynamics_driv

   ! ===============================================================================

   subroutine UpdateNLitterFluxes(this,soilbiogeochem_nitrogenflux_inst,ci,c)

     use clm_varpar, only : i_met_lit

     class(hlm_fates_interface_type), intent(inout)       :: this
     type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: soilbiogeochem_nitrogenflux_inst
     integer                        , intent(in)          :: ci         ! clump index
     integer                        , intent(in)          :: c          ! column index

     integer  :: s                        ! site index
     real(r8) :: dtime
     integer  :: i_lig_lit, i_cel_lit     ! indices for lignan and cellulose

     dtime = get_step_size_real()
     s = this%f2hmap(ci)%hsites(c)
     
     associate(nf_soil => soilbiogeochem_nitrogenflux_inst)

       nf_soil%decomp_npools_sourcesink_col(c,:,:) = 0._r8
       
       if ( .not. use_fates_sp ) then

          ! (gC/m3/timestep)
          !nf_soil%decomp_npools_sourcesink_col(c,1:nlevdecomp,i_met_lit) = &
          !     nf_soil%decomp_npools_sourcesink_col(c,1:nlevdecomp,i_met_lit) + &
          !     this%fates(ci)%bc_out(s)%litt_flux_lab_n_si(1:nlevdecomp)*dtime

          ! Used for mass balance checking (gC/m2/s)
          !nf_soil%fates_litter_flux(c) = sum(this%fates(ci)%bc_out(s)%litt_flux_lab_n_si(1:nlevdecomp) * &
          !                                   this%fates(ci)%bc_in(s)%dz_decomp_sisl(1:nlevdecomp))
          
          i_cel_lit = i_met_lit + 1
          
          !nf_soil%decomp_npools_sourcesink_col(c,1:nlevdecomp,i_cel_lit) = &
          !     nf_soil%decomp_npools_sourcesink_col(c,1:nlevdecomp,i_cel_lit) + &
          !     this%fates(ci)%bc_out(s)%litt_flux_cel_n_si(1:nlevdecomp)*dtime

          !nf_soil%fates_litter_flux(c) = nf_soil%fates_litter_flux(c) + &
          !     sum(this%fates(ci)%bc_out(s)%litt_flux_cel_n_si(1:nlevdecomp) * &
          !         this%fates(ci)%bc_in(s)%dz_decomp_sisl(1:nlevdecomp))

          if (decomp_method == mimics_decomp) then
             ! Mimics has a structural pool, which is cellulose and lignan
             i_lig_lit = i_cel_lit
          elseif(decomp_method == century_decomp ) then
             ! CENTURY has a separate lignan pool from cellulose
             i_lig_lit = i_cel_lit + 1
          end if
        
          !nf_soil%decomp_npools_sourcesink_col(c,1:nlevdecomp,i_lig_lit) = &
          !     nf_soil%decomp_npools_sourcesink_col(c,1:nlevdecomp,i_lig_lit) + &
          !     this%fates(ci)%bc_out(s)%litt_flux_lig_n_si(1:nlevdecomp)*dtime
          
          !nf_soil%fates_litter_flux(c) = nf_soil%fates_litter_flux(c) + &
          !     sum(this%fates(ci)%bc_out(s)%litt_flux_lig_n_si(1:nlevdecomp) * &
          !         this%fates(ci)%bc_in(s)%dz_decomp_sisl(1:nlevdecomp))

          nf_soil%fates_litter_flux = 0._r8
          
       else

          ! In SP mode their is no mass flux between the two 
          nf_soil%fates_litter_flux = 0._r8
          
       end if

     end associate
     
     return
   end subroutine UpdateNLitterFluxes

   ! ===========================================================
   
   subroutine UpdateCLitterFluxes(this,soilbiogeochem_carbonflux_inst,ci,c)

     use clm_varpar, only : i_met_lit

     class(hlm_fates_interface_type), intent(inout)       :: this
     type(soilbiogeochem_carbonflux_type) , intent(inout) :: soilbiogeochem_carbonflux_inst
     integer                        , intent(in)          :: ci         ! clump index
     integer                        , intent(in)          :: c          ! column index
     
     integer  :: s                        ! site index
     real(r8) :: dtime
     integer  :: i_lig_lit, i_cel_lit     ! indices for lignan and cellulose
     
     dtime = get_step_size_real()
     s = this%f2hmap(ci)%hsites(c)

     associate(cf_soil => soilbiogeochem_carbonflux_inst)

       ! This is zeroed in CNDriverNoLeaching -> soilbiogeochem_carbonflux_inst%SetValues()
       ! Which is called prior to this call, which is later in the CNDriverNoLeaching()
       ! routine.
       ! cf_soil%decomp_cpools_sourcesink_col(c,:,:) = 0._r8
       
       if ( .not. use_fates_sp ) then


          call FluxIntoLitterPools(this%fates(ci)%sites(s), &
                                   this%fates(ci)%bc_in(s), &
                                   this%fates(ci)%bc_out(s))

          ! (gC/m3/timestep)
          cf_soil%decomp_cpools_sourcesink_col(c,1:nlevdecomp,i_met_lit) = &
               cf_soil%decomp_cpools_sourcesink_col(c,1:nlevdecomp,i_met_lit) + &
               this%fates(ci)%bc_out(s)%litt_flux_lab_c_si(1:nlevdecomp)*dtime

          ! Used for mass balance checking (gC/m2/s)
          cf_soil%fates_litter_flux(c) = sum(this%fates(ci)%bc_out(s)%litt_flux_lab_c_si(1:nlevdecomp) * &
                                             this%fates(ci)%bc_in(s)%dz_decomp_sisl(1:nlevdecomp))
          
          i_cel_lit = i_met_lit + 1
          
          cf_soil%decomp_cpools_sourcesink_col(c,1:nlevdecomp,i_cel_lit) = &
               cf_soil%decomp_cpools_sourcesink_col(c,1:nlevdecomp,i_cel_lit) + &
               this%fates(ci)%bc_out(s)%litt_flux_cel_c_si(1:nlevdecomp)*dtime

          cf_soil%fates_litter_flux(c) = cf_soil%fates_litter_flux(c) + &
               sum(this%fates(ci)%bc_out(s)%litt_flux_cel_c_si(1:nlevdecomp) * &
                   this%fates(ci)%bc_in(s)%dz_decomp_sisl(1:nlevdecomp))

          if (decomp_method == mimics_decomp) then
             ! Mimics has a structural pool, which is cellulose and lignan
             i_lig_lit = i_cel_lit
          elseif(decomp_method == century_decomp ) then
             ! CENTURY has a separate lignan pool from cellulose
             i_lig_lit = i_cel_lit + 1
          end if
        
          cf_soil%decomp_cpools_sourcesink_col(c,1:nlevdecomp,i_lig_lit) = &
               cf_soil%decomp_cpools_sourcesink_col(c,1:nlevdecomp,i_lig_lit) + &
               this%fates(ci)%bc_out(s)%litt_flux_lig_c_si(1:nlevdecomp)*dtime
          
          cf_soil%fates_litter_flux(c) = cf_soil%fates_litter_flux(c) + &
               sum(this%fates(ci)%bc_out(s)%litt_flux_lig_c_si(1:nlevdecomp) * &
                   this%fates(ci)%bc_in(s)%dz_decomp_sisl(1:nlevdecomp))
          
       else
          ! In SP mode their is no mass flux between the two 
          
          cf_soil%fates_litter_flux = 0._r8
       end if
          
     end associate

     return
   end subroutine UpdateCLitterFluxes

   ! ===================================================================================
   
   subroutine wrap_update_hlmfates_dyn(this, nc, bounds_clump,      &
        waterdiagnosticbulk_inst, canopystate_inst, &
        soilbiogeochem_carbonflux_inst, is_initing_from_restart)

      ! ---------------------------------------------------------------------------------
      ! This routine handles the updating of vegetation canopy diagnostics, (such as lai)
      ! that either requires HLM boundary conditions (like snow accumulation) or
      ! provides boundary conditions (such as vegetation fractional coverage)
      ! ---------------------------------------------------------------------------------

     class(hlm_fates_interface_type), intent(inout) :: this
     type(bounds_type),intent(in)                   :: bounds_clump
     integer                 , intent(in)           :: nc
     type(waterdiagnosticbulk_type)   , intent(inout)        :: waterdiagnosticbulk_inst
     type(canopystate_type)  , intent(inout)        :: canopystate_inst
     type(soilbiogeochem_carbonflux_type), intent(inout) :: soilbiogeochem_carbonflux_inst
                   

     ! is this being called during a read from restart sequence (if so then use the restarted fates
     ! snow depth variable rather than the CLM variable).
     logical                 , intent(in)           :: is_initing_from_restart

     integer :: npatch  ! number of patches in each site
     integer :: ifp     ! index FATES patch
     integer :: p       ! HLM patch index
     integer :: s       ! site index
     integer :: c       ! column index
     integer :: g       ! grid cell

     logical :: dispersal_flag ! local flag to pass to the inside of the site loop
     real(r8) :: areacheck
     call t_startf('fates_wrap_update_hlmfates_dyn')

     associate(                                &
         tlai => canopystate_inst%tlai_patch , &
         elai => canopystate_inst%elai_patch , &
         tsai => canopystate_inst%tsai_patch , &
         esai => canopystate_inst%esai_patch , &
         htop => canopystate_inst%htop_patch , &
         hbot => canopystate_inst%hbot_patch , &
         z0m  => canopystate_inst%z0m_patch  , & ! Output: [real(r8) (:)   ] momentum roughness length (m)
         displa => canopystate_inst%displa_patch, &
         dleaf_patch => canopystate_inst%dleaf_patch, &
         snow_depth => waterdiagnosticbulk_inst%snow_depth_col, &
         frac_sno_eff => waterdiagnosticbulk_inst%frac_sno_eff_col, &
         frac_veg_nosno_alb => canopystate_inst%frac_veg_nosno_alb_patch)


       ! Process input boundary conditions to FATES
       ! --------------------------------------------------------------------------------
       do s=1,this%fates(nc)%nsites
          c = this%f2hmap(nc)%fcolumn(s)
          this%fates(nc)%bc_in(s)%snow_depth_si   = snow_depth(c)
          this%fates(nc)%bc_in(s)%frac_sno_eff_si = frac_sno_eff(c)
       end do

       ! Only update the fates internal snow burial if this is not a restart
       if (.not. is_initing_from_restart) then
          call UpdateFatesAvgSnowDepth(this%fates(nc)%sites,this%fates(nc)%bc_in)
       end if

       ! Canopy diagnostics for FATES
       call canopy_summarization(this%fates(nc)%nsites, &
            this%fates(nc)%sites,  &
            this%fates(nc)%bc_in)            

       ! Canopy diagnostic outputs for HLM
       call update_hlm_dynamics(this%fates(nc)%nsites, &
            this%fates(nc)%sites,  &
            this%f2hmap(nc)%fcolumn, &
            this%fates(nc)%bc_out )

       !------------------------------------------------------------------------
       ! FATES calculation of ligninNratio
       !------------------------------------------------------------------------
       ! If it's the first timestep of a restart & copy_fates_var(c) = .false.
       ! (this will happen in the first timestep of any restart)
       ! then skip this variable because a more accurate value was obtained
       ! from the restart file.
       ! Note 1. I had hoped is_first_restart_step() alone would suffice, but
       ! is_first_restart_step() remained true for the whole first day in my
       ! test. Hence I added the check for copy_fates_var(c) which changes to
       ! .true. the first time that we come through this code in a restart.
       ! Note 2. copy_fates_var is a column-level array to make thread-safe.
       ! Note 3. This if statement is a workaround for restart tests to pass.
       ! Ideally, the fates variable would have the correct value at restart,
       ! but rgknox explains that accomplishing this is more complex given the
       ! current FATES treatment of other litter fluxes passed to the CTSM.
       !------------------------------------------------------------------------
       if ( decomp_method /= no_soil_decomp )then
          do s = 1, this%fates(nc)%nsites
             c = this%f2hmap(nc)%fcolumn(s)
             if (is_first_restart_step() .and. .not. copy_fates_var(c)) then
                copy_fates_var(c) = .true.
             else
                soilbiogeochem_carbonflux_inst%litr_lig_c_to_n_col(c) = &
                  this%fates(nc)%bc_out(s)%litt_flux_ligc_per_n
             end if
          end do
       end if

       !---------------------------------------------------------------------------------
       ! CHANGING STORED WATER DURING PLANT DYNAMICS IS NOT FULLY IMPLEMENTED
       ! LEAVING AS A PLACE-HOLDER FOR NOW.
       !       ! Diagnose water storage in canopy if hydraulics is on
       !       ! This updates the internal value and the bc_out value.
       !       ! If hydraulics is off, it returns 0 storage
       if ( use_fates_planthydro ) then
          do s = 1, this%fates(nc)%nsites
             c = this%f2hmap(nc)%fcolumn(s)
             waterdiagnosticbulk_inst%total_plant_stored_h2o_col(c) = &
                  this%fates(nc)%bc_out(s)%plant_stored_h2o_si
          end do
       end if
       !---------------------------------------------------------------------------------

       ! Convert FATES dynamics into HLM usable information
       ! Initialize weighting variables (note FATES is the only HLM module
       ! that uses "is_veg" and "is_bareground".  The entire purpose of these
       ! variables is to inform patch%wtcol(p).  wt_ed is imposed on wtcol,
       ! but only for FATES columns.

       patch%is_veg(bounds_clump%begp:bounds_clump%endp)        = .false.
       patch%is_bareground(bounds_clump%begp:bounds_clump%endp) = .false.
       patch%wt_ed(bounds_clump%begp:bounds_clump%endp)         = 0.0_r8

       if (fates_seeddisp_cadence /= fates_dispersal_cadence_none) then
          ! Zero the outgoing_local seed values prior to populating with the most recent seed update
          this%fates_seed%outgoing_local(:,:) = 0._r8

          ! check the dispersal time once outside of the loop and set a flag to pass in
          dispersal_flag = .false.
          if (IsItDispersalTime()) dispersal_flag = .true.
       end if

      do s = 1,this%fates(nc)%nsites

          c = this%f2hmap(nc)%fcolumn(s)
          g = col%gridcell(c)

          ! Accumulate seeds from sites to the gridcell local outgoing buffer
          if (fates_seeddisp_cadence /= fates_dispersal_cadence_none) then
             if (dispersal_flag) this%fates_seed%outgoing_local(:,g) = this%fates(nc)%sites(s)%seed_out(:)
          end if

          ! Other modules may have AI's we only flush values
          ! that are on the naturally vegetated columns
          elai(col%patchi(c):col%patchf(c)) = 0.0_r8
          esai(col%patchi(c):col%patchf(c)) = 0.0_r8
          hbot(col%patchi(c):col%patchf(c)) = 0.0_r8
          
          tlai(col%patchi(c):col%patchf(c)) = 0.0_r8
          tsai(col%patchi(c):col%patchf(c)) = 0.0_r8
          htop(col%patchi(c):col%patchf(c)) = 0.0_r8

          ! FATES does not dictate bare-ground so turbulent
          ! variables are not over-written.
          z0m(col%patchi(c)+1:col%patchf(c)) = 0.0_r8
          displa(col%patchi(c)+1:col%patchf(c)) = 0.0_r8
          dleaf_patch(col%patchi(c)+1:col%patchf(c)) = 0.0_r8

          frac_veg_nosno_alb(col%patchi(c):col%patchf(c)) = 0.0_r8

          ! Set the bareground patch indicator
          patch%is_bareground(col%patchi(c)) = .true.
          npatch = this%fates(nc)%sites(s)%youngest_patch%patchno

          ! Precision errors on the canopy_fraction_pa sum, even small (e-12)
          ! do exist, and can create potentially negetive bare-soil fractions
          ! (ie -1e-12 or smaller). Even though this is effectively zero,
          ! it can generate weird logic scenarios in the ctsm/elm code, so we
          ! protext it here with a lower bound of 0.0_r8.

          patch%wt_ed(col%patchi(c)) = max(0.0_r8, &

               1.0_r8-sum(this%fates(nc)%bc_out(s)%canopy_fraction_pa(1:npatch)))

          patch%sp_pftorder_index(col%patchi(c)) = 0 !bg is the 0th patch in the SP FATES structure

          areacheck = patch%wt_ed(col%patchi(c)) ! this is where we start the areachecking

          do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno
             ! for the vegetated patches

             p = ifp+col%patchi(c)

             ! bc_out(s)%canopy_fraction_pa(ifp) is the area fraction
             ! the site's total ground area that is occupied by the
             ! area footprint of the current patch's vegetation canopy

             patch%is_veg(p) = .true.
             patch%wt_ed(p)  = this%fates(nc)%bc_out(s)%canopy_fraction_pa(ifp)
             areacheck = areacheck + patch%wt_ed(p)

             elai(p) = this%fates(nc)%bc_out(s)%elai_pa(ifp)
             esai(p) = this%fates(nc)%bc_out(s)%esai_pa(ifp)
             hbot(p) = this%fates(nc)%bc_out(s)%hbot_pa(ifp)

             tlai(p) = this%fates(nc)%bc_out(s)%tlai_pa(ifp)
             tsai(p) = this%fates(nc)%bc_out(s)%tsai_pa(ifp)
             htop(p) = this%fates(nc)%bc_out(s)%htop_pa(ifp)
             
             if(use_fates_sp.and.abs(tlai(p) - &
                                 this%fates(nc)%bc_out(s)%tlai_pa(ifp)).gt.1e-09)then
               write(iulog,*) 'fates lai not like hlm lai',tlai(p),this%fates(nc)%bc_out(s)%tlai_pa(ifp),ifp
             endif

             frac_veg_nosno_alb(p) = this%fates(nc)%bc_out(s)%frac_veg_nosno_alb_pa(ifp)

             ! Note that while we pass the following values at this point
             ! we have to send the same values after each time-step because
             ! the HLM keeps changing the value and re-setting, so we
             ! re-send instead of re-set. See clm_fates%TransferZ0mDisp()
             z0m(p)    = this%fates(nc)%bc_out(s)%z0m_pa(ifp)
             displa(p) = this%fates(nc)%bc_out(s)%displa_pa(ifp)
             dleaf_patch(p) = this%fates(nc)%bc_out(s)%dleaf_pa(ifp)
          end do ! veg pach

          if(abs(areacheck - 1.0_r8).gt.1.e-9_r8)then
            write(iulog,*) 'area wrong in interface',areacheck - 1.0_r8
            call endrun(msg=errMsg(sourcefile, __LINE__))
          endif

       end do
     end associate

     call t_stopf('fates_wrap_update_hlmfates_dyn')

   end subroutine wrap_update_hlmfates_dyn

   ! ====================================================================================

   subroutine restart( this, bounds_proc, ncid, flag, waterdiagnosticbulk_inst, &
        waterstatebulk_inst, canopystate_inst, soilstate_inst, &
        active_layer_inst, soilbiogeochem_carbonflux_inst, &
        soilbiogeochem_nitrogenflux_inst)

      ! ---------------------------------------------------------------------------------
      ! The ability to restart the model is handled through three different types of calls
      ! "Define" the variables in the restart file, we "read" those variables into memory
      ! or "write" data into the file from memory.  This subroutine accomodates all three
      ! of those modes through the "flag" argument.  FATES as an external model also
      ! requires an initialization step, where we set-up the dimensions, allocate and
      ! flush the memory space that is used to transfer data in and out of the file.  This
      ! Only occurs once, where as the define step occurs every time a file is opened.
      !
      ! Note: waterstate_inst and canopystate_inst are arguments only because following
      ! the reading of variables, it is necessary to update diagnostics of the canopy
      ! throug the interface call clm_fates%wrap_update_hlmfates_dyn() which requires
      ! this information from the HLM.
      ! ---------------------------------------------------------------------------------


     use FatesConstantsMod, only : fates_long_string_length
     use FatesIODimensionsMod, only: fates_bounds_type
     use FatesIOVariableKindMod, only : site_r8, site_int, cohort_r8, cohort_int
     use EDMainMod, only :        ed_update_site
     use FatesInterfaceTypesMod, only:  fates_maxElementsPerSite

      ! Arguments

      class(hlm_fates_interface_type), intent(inout) :: this
      type(bounds_type)              , intent(in)    :: bounds_proc
      type(file_desc_t)              , intent(inout) :: ncid    ! netcdf id
      character(len=*)               , intent(in)    :: flag
      type(waterdiagnosticbulk_type) , intent(inout) :: waterdiagnosticbulk_inst
      type(waterstatebulk_type)      , intent(inout) :: waterstatebulk_inst
      type(canopystate_type)         , intent(inout) :: canopystate_inst
      type(soilstate_type)           , intent(inout) :: soilstate_inst
      type(active_layer_type)        , intent(in)    :: active_layer_inst
      type(soilbiogeochem_carbonflux_type), intent(inout) :: soilbiogeochem_carbonflux_inst
      type(soilbiogeochem_nitrogenflux_type), intent(inout) :: soilbiogeochem_nitrogenflux_inst
      
      ! Locals
      type(bounds_type) :: bounds_clump
      integer           :: nc
      integer           :: nclumps
      type(fates_bounds_type) :: fates_bounds
      type(fates_bounds_type) :: fates_clump
      integer                 :: c   ! HLM column index
      integer                 :: s   ! Fates site index
      integer                 :: g   ! grid-cell index
      integer                 :: p   ! HLM patch index
      integer                 :: ft  ! plant functional type
      integer                 :: dk_index
      integer                 :: nlevsoil
      character(len=fates_long_string_length) :: ioname
      integer                 :: nvar
      integer                 :: ivar
      logical                 :: readvar
      logical, save           :: initialized = .false.

     call t_startf('fates_restart')

      nclumps = get_proc_clumps()

      ! ---------------------------------------------------------------------------------
      ! note (rgk: 11-2016) The history and restart intialization process assumes
      ! that the number of site/columns active is a static entity.  Thus
      ! we only allocate the mapping tables for the column/sites we start with.
      ! If/when we start having dynamic column/sites (for reasons uknown as of yet)
      ! we will need to re-evaluate the allocation of the mapping tables so they
      ! can be unallocated,reallocated and set every time a new column/site is spawned
      ! ---------------------------------------------------------------------------------

      ! ---------------------------------------------------------------------------------
      ! Only initialize the FATES restart structures the first time it is called
      ! Note that the allocations involved with initialization are static.
      ! This is because the array spaces for IO span the entire column, patch and cohort
      ! range on the proc.
      ! With DYNAMIC LANDUNITS or SPAWNING NEW OR CULLING OLD SITES:
      ! we will in that case have to de-allocate, reallocate and then re-set the mapping
      ! tables:  this%fates_restart%restart_map(nc)
      ! I think that is it...
      ! ---------------------------------------------------------------------------------

      ! Set the FATES global time and date variables
      call GetAndSetTime

      if(.not.initialized) then

         initialized=.true.

         ! ------------------------------------------------------------------------------
         ! PART I: Set FATES DIMENSIONING INFORMATION
         ! ------------------------------------------------------------------------------

         call hlm_bounds_to_fates_bounds(bounds_proc, fates_bounds)

         call this%fates_restart%Init(nclumps, fates_bounds)

         ! Define the bounds on the first dimension for each thread
         !$OMP PARALLEL DO PRIVATE (nc,bounds_clump,fates_clump)
         do nc = 1,nclumps
            call get_clump_bounds(nc, bounds_clump)

            ! thread bounds for patch
            call hlm_bounds_to_fates_bounds(bounds_clump, fates_clump)
            call this%fates_restart%SetThreadBoundsEach(nc, fates_clump)
         end do
         !$OMP END PARALLEL DO

         !$OMP PARALLEL DO PRIVATE (nc,s,c,g,bounds_clump)
         do nc = 1,nclumps

            call get_clump_bounds(nc, bounds_clump)
            allocate(this%fates_restart%restart_map(nc)%site_index(this%fates(nc)%nsites))
            allocate(this%fates_restart%restart_map(nc)%cohort1_index(this%fates(nc)%nsites))
            do s=1,this%fates(nc)%nsites
               c = this%f2hmap(nc)%fcolumn(s)
               this%fates_restart%restart_map(nc)%site_index(s)   = c
               g = col%gridcell(c)
               this%fates_restart%restart_map(nc)%cohort1_index(s) = (g-1)*fates_maxElementsPerSite + 1
            end do

         end do
         !$OMP END PARALLEL DO

         ! ------------------------------------------------------------------------------------
         ! PART II: USE THE JUST DEFINED DIMENSIONS TO ASSEMBLE THE VALID IO TYPES
         ! INTERF-TODO: THESE CAN ALL BE EMBEDDED INTO A SUBROUTINE IN HISTORYIOMOD
         ! ------------------------------------------------------------------------------------
         call this%fates_restart%assemble_restart_output_types()


         ! ------------------------------------------------------------------------------------
         ! PART III: DEFINE THE LIST OF OUTPUT VARIABLE OBJECTS, AND REGISTER THEM WITH THE
         ! HLM ACCORDING TO THEIR TYPES
         ! ------------------------------------------------------------------------------------
         call this%fates_restart%initialize_restart_vars()

      end if

      ! ---------------------------------------------------------------------------------
      ! If we are writing, we must loop through our linked list structures and transfer the
      ! information in the linked lists (FATES state memory) to the output vectors.
      ! ---------------------------------------------------------------------------------

      if(flag=='write')then
         !$OMP PARALLEL DO PRIVATE (nc)
         do nc = 1, nclumps
            if (this%fates(nc)%nsites>0) then
               call this%fates_restart%set_restart_vectors(nc,this%fates(nc)%nsites, &
                                                           this%fates(nc)%sites)
            end if
         end do
         !$OMP END PARALLEL DO
      end if

      ! ---------------------------------------------------------------------------------
      ! In all cases, iterate through the list of variable objects
      ! and either define, write or read to the NC buffer
      ! This seems strange, but keep in mind that the call to restartvar()
      ! has a different function in all three cases.
      ! ---------------------------------------------------------------------------------

      nvar = this%fates_restart%num_restart_vars()
      do ivar = 1, nvar

         associate( vname => this%fates_restart%rvars(ivar)%vname, &
              vunits      => this%fates_restart%rvars(ivar)%units,   &
              vlong       => this%fates_restart%rvars(ivar)%long )

           dk_index = this%fates_restart%rvars(ivar)%dim_kinds_index
           ioname = trim(this%fates_restart%dim_kinds(dk_index)%name)

           select case(trim(ioname))
           case(cohort_r8)

              call restartvar(ncid=ncid, flag=flag, varname=trim(vname), &
                    xtype=ncd_double,dim1name=trim('cohort'),long_name=trim(vlong), &
                    units=trim(vunits),interpinic_flag='interp', &
                    data=this%fates_restart%rvars(ivar)%r81d,readvar=readvar)

           case(site_r8)

              call restartvar(ncid=ncid, flag=flag, varname=trim(vname), &
                    xtype=ncd_double,dim1name=trim('column'),long_name=trim(vlong), &
                    units=trim(vunits),interpinic_flag='interp', &
                    data=this%fates_restart%rvars(ivar)%r81d,readvar=readvar)

           case(cohort_int)

              call restartvar(ncid=ncid, flag=flag, varname=trim(vname), &
                    xtype=ncd_int,dim1name=trim('cohort'),long_name=trim(vlong), &
                    units=trim(vunits),interpinic_flag='interp', &
                    data=this%fates_restart%rvars(ivar)%int1d,readvar=readvar)

           case(site_int)

              call restartvar(ncid=ncid, flag=flag, varname=trim(vname), &
                    xtype=ncd_int,dim1name=trim('column'),long_name=trim(vlong), &
                    units=trim(vunits),interpinic_flag='interp', &
                    data=this%fates_restart%rvars(ivar)%int1d,readvar=readvar)

           case default
              write(iulog,*) 'A FATES iotype was created that was not registerred'
              write(iulog,*) 'in CLM.:',trim(ioname)
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end select

         end associate
      end do

      ! ---------------------------------------------------------------------------------
      ! If we are in a read mode, then we have just populated the sparse vectors
      ! in the IO object list. The data in these vectors needs to be transferred
      ! to the linked lists to populate the state memory.
      ! ---------------------------------------------------------------------------------

      if(flag=='read')then

         !$OMP PARALLEL DO PRIVATE (nc,bounds_clump,s)
         do nc = 1, nclumps
            if (this%fates(nc)%nsites>0) then

               call get_clump_bounds(nc, bounds_clump)

               ! ------------------------------------------------------------------------
               ! Convert newly read-in vectors into the FATES namelist state variables
               ! ------------------------------------------------------------------------
               call this%fates_restart%create_patchcohort_structure(nc, &
                    this%fates(nc)%nsites, this%fates(nc)%sites, this%fates(nc)%bc_in, &
                    this%fates(nc)%bc_out)

               call this%fates_restart%get_restart_vectors(nc, this%fates(nc)%nsites, &
                    this%fates(nc)%sites )

               ! I think ed_update_site and update_hlmfates_dyn are doing some similar
               ! update type stuff, should consolidate (rgk 11-2016)
               do s = 1,this%fates(nc)%nsites

                  c = this%f2hmap(nc)%fcolumn(s)

                  this%fates(nc)%bc_in(s)%max_rooting_depth_index_col = &
                       min(this%fates(nc)%bc_in(s)%nlevsoil, active_layer_inst%altmax_lastyear_indx_col(c))

                  ! When restarting the model, this subroutine has several
                  ! procedures that are incremental or don't need to be performed for 
                  ! during the restart sequence. For the prior, we don't want the restarted
                  ! run to call these routines more than would had been called during
                  ! a continuous simulation period, as it would change results. So 
                  ! we pass in the "is_restarting=.true." flag so we can bypass those procedures

                  call ed_update_site( this%fates(nc)%sites(s), &
                        this%fates(nc)%bc_in(s), &
                        this%fates(nc)%bc_out(s), &
                        is_restarting = .true. )

               end do

               if(use_fates_sp)then
                  do s = 1,this%fates(nc)%nsites
                     c = this%f2hmap(nc)%fcolumn(s)
                     do ft = surfpft_lb,surfpft_ub  !set of pfts in HLM
                        ! here we are mapping from P space in the HLM to FT space in the sp_input arrays.
                        p = ft + col%patchi(c) ! for an FT of 1 we want to use
                        this%fates(nc)%bc_in(s)%hlm_sp_tlai(ft) = canopystate_inst%tlai_input_patch(p)
                        this%fates(nc)%bc_in(s)%hlm_sp_tsai(ft) = canopystate_inst%tsai_input_patch(p)
                        this%fates(nc)%bc_in(s)%hlm_sp_htop(ft) = canopystate_inst%htop_input_patch(p)
                        if(canopystate_inst%htop_input_patch(p).lt.1.0e-20)then ! zero htop causes inifinite/nans. This is
                           this%fates(nc)%bc_in(s)%hlm_sp_htop(ft) = 0.01_r8
                        endif
                     end do ! p
                  end do ! c
                end if ! SP

               ! ------------------------------------------------------------------------
               ! Re-populate all the hydraulics variables that are dependent
               ! on the key hydro state variables and plant carbon/geometry
               ! ------------------------------------------------------------------------
               if (use_fates_planthydro) then

                  do s = 1,this%fates(nc)%nsites
                     c = this%f2hmap(nc)%fcolumn(s)
                     nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil
                     this%fates(nc)%bc_in(s)%hksat_sisl(1:nlevsoil) = &
                          soilstate_inst%hksat_col(c,1:nlevsoil)

                     this%fates(nc)%bc_in(s)%watsat_sisl(1:nlevsoil) = &
                          soilstate_inst%watsat_col(c,1:nlevsoil)

                     this%fates(nc)%bc_in(s)%watres_sisl(1:nlevsoil) = &
                          soilstate_inst%watres_col(c,1:nlevsoil)

                     this%fates(nc)%bc_in(s)%sucsat_sisl(1:nlevsoil) = &
                          soilstate_inst%sucsat_col(c,1:nlevsoil)

                     this%fates(nc)%bc_in(s)%bsw_sisl(1:nlevsoil) = &
                          soilstate_inst%bsw_col(c,1:nlevsoil)

                     this%fates(nc)%bc_in(s)%h2o_liq_sisl(1:nlevsoil) = &
                          waterstatebulk_inst%h2osoi_liq_col(c,1:nlevsoil)
                  end do


                  call RestartHydrStates(this%fates(nc)%sites,  &
                                         this%fates(nc)%nsites, &
                                         this%fates(nc)%bc_in,  &
                                         this%fates(nc)%bc_out)
               end if

               ! ------------------------------------------------------------------------
               ! Update diagnostics of FATES ecosystem structure used in HLM.
               ! ------------------------------------------------------------------------
               call this%wrap_update_hlmfates_dyn(nc,bounds_clump, &
                     waterdiagnosticbulk_inst,canopystate_inst, &
                     soilbiogeochem_carbonflux_inst, .true.)

               ! ------------------------------------------------------------------------
               ! Update the 3D patch level radiation absorption fractions
               ! ------------------------------------------------------------------------
               call this%fates_restart%update_3dpatch_radiation(this%fates(nc)%nsites, &
                                                                this%fates(nc)%sites, &
                                                                this%fates(nc)%bc_out)


               ! ------------------------------------------------------------------------
               ! Flush FATES history variables.
               ! The flushing process sets all columns outside of FATES perview to
               ! the ignore value. This only needs to be done once, because FATES will
               ! not overwright values outside the columns that it is in charge of.
               ! ------------------------------------------------------------------------

               call fates_hist%flush_all_hvars(nc)    
               
               call fates_hist%update_history_dyn( nc,                     &
                                                   this%fates(nc)%nsites,  &
                                                   this%fates(nc)%sites,   &
                                                   this%fates(nc)%bc_in)


            end if
         end do
         !$OMP END PARALLEL DO

         ! Disperse seeds
         if (fates_seeddisp_cadence /= fates_dispersal_cadence_none) then
            call this%WrapGlobalSeedDispersal(is_restart_flag=.true.)
         end if

      end if

     call t_stopf('fates_restart')

      return
   end subroutine restart

   !=====================================================================================

   subroutine init_coldstart(this, waterstatebulk_inst, waterdiagnosticbulk_inst, &
        canopystate_inst, soilstate_inst, soilbiogeochem_carbonflux_inst)


     ! Arguments
     class(hlm_fates_interface_type), intent(inout) :: this
     type(waterstatebulk_type)          , intent(inout) :: waterstatebulk_inst
     type(waterdiagnosticbulk_type)          , intent(inout) :: waterdiagnosticbulk_inst
     type(canopystate_type)         , intent(inout) :: canopystate_inst
     type(soilstate_type)           , intent(inout) :: soilstate_inst
     type(soilbiogeochem_carbonflux_type), intent(inout) :: soilbiogeochem_carbonflux_inst

     ! locals
     integer                                        :: nclumps
     integer                                        :: nc
     type(bounds_type)                              :: bounds_clump
     ! locals
     real(r8) :: vol_ice
     real(r8) :: eff_porosity
     integer :: nlevsoil  ! Number of soil layers at each site
     integer :: i, j
     integer :: s
     integer :: c, g
     integer :: p   ! HLM patch index
     integer :: ft  ! plant functional type


     call t_startf('fates_initcoldstart')

     ! Set the FATES global time and date variables
     call GetAndSetTime

     nclumps = get_proc_clumps()

     !$OMP PARALLEL DO PRIVATE (nc,bounds_clump,s,c,j,vol_ice,eff_porosity)
     do nc = 1, nclumps

        if ( this%fates(nc)%nsites>0 ) then

           call get_clump_bounds(nc, bounds_clump)

           do s = 1,this%fates(nc)%nsites
              call init_site_vars(this%fates(nc)%sites(s), &
                                  this%fates(nc)%bc_in(s), &
                                  this%fates(nc)%bc_out(s) )
              call zero_site(this%fates(nc)%sites(s))
           end do

           call set_site_properties(this%fates(nc)%nsites, &
                                    this%fates(nc)%sites,  &
                                    this%fates(nc)%bc_in)


           ! ----------------------------------------------------------------------------
           ! Initialize satellite phenology values if turned on
           ! ----------------------------------------------------------------------------
            if(use_fates_sp)then
               do s = 1,this%fates(nc)%nsites
                  c = this%f2hmap(nc)%fcolumn(s)
                  do ft = surfpft_lb,surfpft_ub
                     ! here we are mapping from P space in the HLM to FT space in the sp_input arrays.
                     p = ft + col%patchi(c) ! for an FT of 1 we want to use
                     this%fates(nc)%bc_in(s)%hlm_sp_tlai(ft) = canopystate_inst%tlai_input_patch(p)
                     this%fates(nc)%bc_in(s)%hlm_sp_tsai(ft) = canopystate_inst%tsai_input_patch(p)
                     this%fates(nc)%bc_in(s)%hlm_sp_htop(ft) = canopystate_inst%htop_input_patch(p)
                     if(canopystate_inst%htop_input_patch(p).lt.1.0e-20)then ! zero htop causes inifinite/nans. This is
                        this%fates(nc)%bc_in(s)%hlm_sp_htop(ft) = 0.01_r8
                     endif
                  end do ! p
               end do ! c
            end if ! SP

           ! ----------------------------------------------------------------------------
           ! Initialize Hydraulics Code if turned on
           ! Called prior to init_patches(). Site level rhizosphere shells must
           ! be set prior to cohort initialization.
           ! ----------------------------------------------------------------------------

           if (use_fates_planthydro) then

              do s = 1,this%fates(nc)%nsites
                 c = this%f2hmap(nc)%fcolumn(s)

                 nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil

                 this%fates(nc)%bc_in(s)%watsat_sisl(1:nlevsoil) = &
                      soilstate_inst%watsat_col(c,1:nlevsoil)

                 this%fates(nc)%bc_in(s)%sucsat_sisl(1:nlevsoil) = &
                      soilstate_inst%sucsat_col(c,1:nlevsoil)

                 this%fates(nc)%bc_in(s)%bsw_sisl(1:nlevsoil) = &
                      soilstate_inst%bsw_col(c,1:nlevsoil)

                 this%fates(nc)%bc_in(s)%h2o_liq_sisl(1:nlevsoil) = &
                      waterstatebulk_inst%h2osoi_liq_col(c,1:nlevsoil)

                 this%fates(nc)%bc_in(s)%hksat_sisl(1:nlevsoil) = &
                       soilstate_inst%hksat_col(c,1:nlevsoil)

                 do j = 1, nlevsoil
                    vol_ice = min(soilstate_inst%watsat_col(c,j), &
                          waterstatebulk_inst%h2osoi_ice_col(c,j)/(col%dz(c,j)*denice))
                    eff_porosity = max(0.01_r8,soilstate_inst%watsat_col(c,j)-vol_ice)
                    this%fates(nc)%bc_in(s)%eff_porosity_sl(j) = eff_porosity
                 end do

              end do

              call HydrSiteColdStart(this%fates(nc)%sites,this%fates(nc)%bc_in)
           end if

           do s = 1,this%fates(nc)%nsites
              c = this%f2hmap(nc)%fcolumn(s)
              g = col%gridcell(c)

              if (use_fates_luh) then
                    this%fates(nc)%bc_in(s)%hlm_luh_states = landuse_states(:,g)
                    this%fates(nc)%bc_in(s)%hlm_luh_state_names = landuse_state_varnames
                    this%fates(nc)%bc_in(s)%hlm_luh_transitions = landuse_transitions(:,g)
                    this%fates(nc)%bc_in(s)%hlm_luh_transition_names = landuse_transition_varnames

                    if (trim(fates_harvest_mode) == fates_harvest_luh_area .or. &
                        trim(fates_harvest_mode) == fates_harvest_luh_mass) then
                       this%fates(nc)%bc_in(s)%hlm_harvest_rates = landuse_harvest(:,g)
                       this%fates(nc)%bc_in(s)%hlm_harvest_catnames = landuse_harvest_varnames
                       this%fates(nc)%bc_in(s)%hlm_harvest_units = landuse_harvest_units
                    end if
              end if
           end do

           ! Newly initialized patches need a starting temperature
           call init_patches(this%fates(nc)%nsites, this%fates(nc)%sites, &
                             this%fates(nc)%bc_in)

           do s = 1,this%fates(nc)%nsites

              c = this%f2hmap(nc)%fcolumn(s)

              ! CLM has not calculated the maximum rooting depth
              ! from last year, it won't until the beginning of the
              ! time-step loop. Therefore, we just initialize fluxes
              ! into the litter pool in a trivial way prior to timestepping
              this%fates(nc)%bc_in(s)%max_rooting_depth_index_col = this%fates(nc)%bc_in(s)%nlevsoil

              call ed_update_site(this%fates(nc)%sites(s), &
                    this%fates(nc)%bc_in(s), &
                    this%fates(nc)%bc_out(s), &
                    is_restarting = .false.)

           end do

           ! ------------------------------------------------------------------------
           ! Update diagnostics of FATES ecosystem structure used in HLM.
           ! ------------------------------------------------------------------------
           call this%wrap_update_hlmfates_dyn(nc,bounds_clump, &
                waterdiagnosticbulk_inst,canopystate_inst, &
                soilbiogeochem_carbonflux_inst, .false.)

           ! ------------------------------------------------------------------------
           ! Flush and zero FATES history variables.
           ! The flushing process sets all columns outside of FATES perview to
           ! the ignore value. This only needs to be done once, because FATES will
           ! not overwright values outside the columns that it is in charge of.
           ! We also start off by setting all values on FATES columns to zero.
           ! ------------------------------------------------------------------------

           call fates_hist%flush_all_hvars(nc)
           
           call fates_hist%update_history_dyn( nc,                     &
                this%fates(nc)%nsites,                                  &
                this%fates(nc)%sites,                                   &
                this%fates(nc)%bc_in) 

        end if
     end do
     !$OMP END PARALLEL DO

     call t_stopf('fates_initcoldstart')

   end subroutine init_coldstart

   ! ======================================================================================

   subroutine wrap_sunfrac(this,nc,atm2lnd_inst,canopystate_inst)

      ! ---------------------------------------------------------------------------------
      ! This interface function is a wrapper call on ED_SunShadeFracs. The only
      ! returned variable is a patch vector, fsun_patch, which describes the fraction
      ! of the canopy that is exposed to sun.
      ! ---------------------------------------------------------------------------------

      ! Input Arguments
      class(hlm_fates_interface_type), intent(inout) :: this

      integer, intent(in)                  :: nc

      ! direct and diffuse downwelling radiation (W/m2)
      type(atm2lnd_type),intent(in)        :: atm2lnd_inst

      ! Input/Output Arguments to CLM
      type(canopystate_type),intent(inout) :: canopystate_inst

      ! Local Variables
      integer  :: p                           ! global index of the host patch
      integer  :: g                           ! global index of the host gridcell
      integer  :: c                           ! global index of the host column

      integer  :: s                           ! FATES site index
      integer  :: ifp                         ! FATEs patch index
                                              ! this is the order increment of patch
                                              ! on the site

      type(fates_patch_type), pointer :: cpatch  ! c"urrent" patch  INTERF-TODO: SHOULD
                                              ! BE HIDDEN AS A FATES PRIVATE

     call t_startf('fates_wrapsunfrac')

      associate( forc_solad => atm2lnd_inst%forc_solad_not_downscaled_grc, &
                 forc_solai => atm2lnd_inst%forc_solai_grc, &
                 fsun       => canopystate_inst%fsun_patch, &
                 laisun     => canopystate_inst%laisun_patch, &
                 laisha     => canopystate_inst%laisha_patch )

        ! -------------------------------------------------------------------------------
        ! Convert input BC's
        ! The sun-shade calculations are performed only on FATES patches
        ! -------------------------------------------------------------------------------

        do s = 1, this%fates(nc)%nsites
           c = this%f2hmap(nc)%fcolumn(s)
           g = col%gridcell(c)

           do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno
              this%fates(nc)%bc_in(s)%solad_parb(ifp,:) = forc_solad(g,:)
              this%fates(nc)%bc_in(s)%solai_parb(ifp,:) = forc_solai(g,:)
           end do
        end do


        ! -------------------------------------------------------------------------------
        ! Call FATES public function to calculate internal sun/shade structures
        ! as well as total patch sun/shade fraction output boundary condition
        ! -------------------------------------------------------------------------------

        call FatesSunShadeFracs(this%fates(nc)%nsites, &
             this%fates(nc)%sites,  &
             this%fates(nc)%bc_in,  &
             this%fates(nc)%bc_out)

        ! -------------------------------------------------------------------------------
        ! Transfer the FATES output boundary condition for canopy sun/shade fraction
        ! to the HLM
        ! -------------------------------------------------------------------------------

        do s = 1, this%fates(nc)%nsites
           c = this%f2hmap(nc)%fcolumn(s)
           do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno

              p = ifp+col%patchi(c)

              fsun(p)   = this%fates(nc)%bc_out(s)%fsun_pa(ifp)
              laisun(p) = this%fates(nc)%bc_out(s)%laisun_pa(ifp)
              laisha(p) = this%fates(nc)%bc_out(s)%laisha_pa(ifp)
           end do
        end do

      end associate

     call t_stopf('fates_wrapsunfrac')

   end subroutine wrap_sunfrac

   ! ===================================================================================

   subroutine prep_canopyfluxes(this, nc, fn, filterp, photosyns_inst)

     ! ----------------------------------------------------------------------
     ! the main function for calculating photosynthesis is called within a
     ! loop based on convergence.  Some intitializations, including
     ! canopy resistance must be intitialized before the loop
     ! ----------------------------------------------------------------------

     ! Arguments
     class(hlm_fates_interface_type), intent(inout) :: this
     integer, intent(in)                            :: nc
     integer, intent(in)                            :: fn
     integer, intent(in)                            :: filterp(fn)
     type(photosyns_type),intent(inout)             :: photosyns_inst
     ! locals
     integer                                        :: f,p,c,s
     ! parameters
     integer,parameter                              :: rsmax0 = 2.e4_r8

     call t_startf('fates_prepcanfluxes')

     do s = 1, this%fates(nc)%nsites
        ! filter flag == 1 means that this patch has not been called for photosynthesis
        this%fates(nc)%bc_in(s)%filter_photo_pa(:) = 1

        ! set transpiration input boundary condition to zero. The exposed
        ! vegetation filter may not even call every patch.
        if (use_fates_planthydro) then
            this%fates(nc)%bc_in(s)%qflx_transp_pa(:) = 0._r8
        end if

     end do

     call t_stopf('fates_prepcanfluxes')

  end subroutine prep_canopyfluxes

   ! ====================================================================================

   subroutine wrap_btran(this,nc,fn,filterc,soilstate_inst, &
                         waterdiagnosticbulk_inst, temperature_inst, energyflux_inst,  &
                         soil_water_retention_curve)

      ! ---------------------------------------------------------------------------------
      ! This subroutine calculates btran for FATES, this will be an input boundary
      ! condition for FATES photosynthesis/transpiration.
      !
      ! This subroutine also calculates rootr
      !
      ! ---------------------------------------------------------------------------------

      use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type

      ! Arguments
      class(hlm_fates_interface_type), intent(inout) :: this
      integer                , intent(in)            :: nc
      integer                , intent(in)            :: fn
      integer                , intent(in)            :: filterc(fn) ! This is a list of
                                                                        ! columns with exposed veg
      type(soilstate_type)   , intent(inout)         :: soilstate_inst
      type(waterdiagnosticbulk_type)  , intent(in)            :: waterdiagnosticbulk_inst
      type(temperature_type) , intent(in)            :: temperature_inst
      type(energyflux_type)  , intent(inout)         :: energyflux_inst
      class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve

      ! local variables
      real(r8) :: smp_node ! Soil suction potential, negative, [mm]
      real(r8) :: s_node
      integer  :: s
      integer  :: c
      integer  :: j
      integer  :: ifp
      integer  :: p
      integer  :: nlevsoil

     call t_startf('fates_wrapbtran')

      associate(&
         sucsat      => soilstate_inst%sucsat_col           , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
         watsat      => soilstate_inst%watsat_col           , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         bsw         => soilstate_inst%bsw_col              , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
         eff_porosity => soilstate_inst%eff_porosity_col    , & ! Input:  [real(r8) (:,:) ]  effective porosity = porosity - vol_ice
         t_soisno    => temperature_inst%t_soisno_col       , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)
         h2osoi_liqvol => waterdiagnosticbulk_inst%h2osoi_liqvol_col , & ! Input: [real(r8) (:,:) ]  liquid volumetric moisture, will be used for BeTR
         btran       => energyflux_inst%btran_patch         , & ! Output: [real(r8) (:)   ]  transpiration wetness factor (0 to 1)
         rresis      => energyflux_inst%rresis_patch        , & ! Output: [real(r8) (:,:) ]  root resistance by layer (0-1)  (nlevgrnd)
         rootr       => soilstate_inst%rootr_patch          & ! Output: [real(r8) (:,:) ]  Fraction of water uptake in each layer
         )

        ! -------------------------------------------------------------------------------
        ! Convert input BC's
        ! Critical step: a filter is being passed in that dictates which columns have
        ! exposed vegetation (above snow).  This is necessary, because various hydrologic
        ! variables like h2osoi_liqvol are not calculated and will have uninitialized
        ! values outside this list.
        !
        ! bc_in(s)%filter_btran      (this is in, but is also used in this subroutine)
        !
        ! We also filter a second time within this list by determining which soil layers
        ! have conditions for active uptake based on soil moisture and temperature. This
        ! must be determined by FATES (science stuff).  But the list of layers and patches
        ! needs to be passed back to the interface, because it then needs to request
        ! suction on these layers via CLM/ALM functions.  We cannot wide-swath calculate
        ! this on all layers, because values with no moisture or low temps will generate
        ! unstable values and cause sigtraps.
        ! -------------------------------------------------------------------------------

        do s = 1, this%fates(nc)%nsites
           c = this%f2hmap(nc)%fcolumn(s)
           nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil

           ! Check to see if this column is in the exposed veg filter
           if( any(filterc==c) )then

              this%fates(nc)%bc_in(s)%filter_btran = .true.
              do j = 1,nlevsoil
                 this%fates(nc)%bc_in(s)%tempk_sl(j)         = t_soisno(c,j)
                 this%fates(nc)%bc_in(s)%h2o_liqvol_sl(j)    = h2osoi_liqvol(c,j)
                 this%fates(nc)%bc_in(s)%eff_porosity_sl(j)  = eff_porosity(c,j)
                 this%fates(nc)%bc_in(s)%watsat_sl(j)        = watsat(c,j)
              end do

           else
              this%fates(nc)%bc_in(s)%filter_btran = .false.
              this%fates(nc)%bc_in(s)%tempk_sl(:)         = -999._r8
              this%fates(nc)%bc_in(s)%h2o_liqvol_sl(:)    = -999._r8
              this%fates(nc)%bc_in(s)%eff_porosity_sl(:)  = -999._r8
              this%fates(nc)%bc_in(s)%watsat_sl(:)        = -999._r8
           end if

        end do

        ! -------------------------------------------------------------------------------
        ! This function evaluates the ground layer to determine if
        ! root water uptake can happen, and soil suction should even
        ! be calculated.  We ask FATES for a boundary condition output
        ! logical because we don't want science calculations in the interface
        ! yet... hydrology (suction calculation) is provided by the host
        ! so we need fates to tell us where to calculate suction
        ! but not calculate it itself. Yeah, complicated, but thats life.
        ! -------------------------------------------------------------------------------
        call get_active_suction_layers(this%fates(nc)%nsites, &
             this%fates(nc)%sites,  &
             this%fates(nc)%bc_in,  &
             this%fates(nc)%bc_out)

        ! Now that the active layers of water uptake have been decided by fates
        ! Calculate the suction that is passed back to fates
        ! Note that the filter_btran is unioned with active_suction_sl

        do s = 1, this%fates(nc)%nsites
           c = this%f2hmap(nc)%fcolumn(s)
           nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil
           do j = 1,nlevsoil
              if(this%fates(nc)%bc_out(s)%active_suction_sl(j)) then
                 s_node = max(h2osoi_liqvol(c,j)/eff_porosity(c,j),0.01_r8)
                 call soil_water_retention_curve%soil_suction(c,j,s_node, soilstate_inst, smp_node)
                 this%fates(nc)%bc_in(s)%smp_sl(j)           = smp_node
              end if
           end do
        end do

        ! -------------------------------------------------------------------------------
        ! Suction and active uptake layers calculated, lets calculate uptake (btran)
        ! This will calculate internals, as well as output boundary conditions:
        ! btran, rootr
        ! -------------------------------------------------------------------------------

        call btran_ed(this%fates(nc)%nsites, &
             this%fates(nc)%sites,  &
             this%fates(nc)%bc_in,  &
             this%fates(nc)%bc_out)

        ! -------------------------------------------------------------------------------
        ! Convert output BC's
        ! For CLM/ALM this wrapper provides return variables that should
        ! be similar to that of calc_root_moist_stress().  However,
        ! CLM/ALM-FATES simulations will no make use of rresis or btran
        ! outside of FATES. We do not have code in place to calculate or
        ! rresis right now, so we force to bad.  We have btran calculated so we
        ! pass it in case people want diagnostics.  rootr is actually the only
        ! variable that will be used, as it is needed to help distribute the
        ! the transpiration sink to the appropriate layers. (RGK)
        ! -------------------------------------------------------------------------------

        do s = 1, this%fates(nc)%nsites
           nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil
           c = this%f2hmap(nc)%fcolumn(s)
           do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno

              p = ifp+col%patchi(c)

              do j = 1,nlevsoil

                 rresis(p,j) = -999.9  ! We do not calculate this correctly
                 ! it should not thought of as valid output until we decide to.
                 rootr(p,j)  = this%fates(nc)%bc_out(s)%rootr_pasl(ifp,j)
                 btran(p)    = this%fates(nc)%bc_out(s)%btran_pa(ifp)

              end do
           end do
        end do
      end associate

     call t_stopf('fates_wrapbtran')

   end subroutine wrap_btran

   ! ====================================================================================

   subroutine wrap_photosynthesis(this, nc, bounds, fn, filterp, &
         esat_tv, eair, oair, cair, rb, dayl_factor,             &
         atm2lnd_inst, temperature_inst, canopystate_inst, photosyns_inst)

    use shr_log_mod       , only : errMsg => shr_log_errMsg
    use abortutils        , only : endrun
    use decompMod         , only : bounds_type
    use clm_varcon        , only : tfrz, namep
    use clm_varctl        , only : iulog
    use PatchType         , only : patch
    use quadraticMod      , only : quadratic
    use EDtypesMod        , only : ed_site_type
    use FatesPatchMod,      only : fates_patch_type
    use FatesCohortMod    , only : fates_cohort_type

    !
    ! !ARGUMENTS:
    class(hlm_fates_interface_type), intent(inout) :: this
    integer                , intent(in)            :: nc                          ! clump index
    type(bounds_type)      , intent(in)            :: bounds
    integer                , intent(in)            :: fn                          ! size of pft filter
    integer                , intent(in)            :: filterp(fn)                 ! pft filter
    real(r8)               , intent(in)            :: esat_tv(bounds%begp: )      ! saturation vapor pressure at t_veg (Pa)
    real(r8)               , intent(in)            :: eair( bounds%begp: )        ! vapor pressure of canopy air (Pa)
    real(r8)               , intent(in)            :: oair( bounds%begp: )        ! Atmospheric O2 partial pressure (Pa)
    real(r8)               , intent(in)            :: cair( bounds%begp: )        ! Atmospheric CO2 partial pressure (Pa)
    real(r8)               , intent(in)            :: rb( bounds%begp: )          ! boundary layer resistance (s/m)
    real(r8)               , intent(in)            :: dayl_factor( bounds%begp: ) ! scalar (0-1) for daylength
    type(atm2lnd_type)     , intent(in)            :: atm2lnd_inst
    type(temperature_type) , intent(in)            :: temperature_inst
    type(canopystate_type) , intent(inout)         :: canopystate_inst
    type(photosyns_type)   , intent(inout)         :: photosyns_inst

    integer                                        :: nlevsoil  ! number of soil layers in this site
    integer                                        :: s,c,p,ifp,j,icp
    real(r8)                                       :: dtime
    call t_startf('fates_psn')

    associate(&
          t_soisno  => temperature_inst%t_soisno_col , &
          t_veg     => temperature_inst%t_veg_patch  , &
          tgcm      => temperature_inst%thm_patch    , &
          forc_pbot => atm2lnd_inst%forc_pbot_downscaled_col, &
          rssun     => photosyns_inst%rssun_patch  , &
          rssha     => photosyns_inst%rssha_patch,   &
          psnsun    => photosyns_inst%psnsun_patch,  &
          psnsha    => photosyns_inst%psnsha_patch)

      do s = 1, this%fates(nc)%nsites

         c = this%f2hmap(nc)%fcolumn(s)

         nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil

         do j = 1,nlevsoil
            this%fates(nc)%bc_in(s)%t_soisno_sl(j)   = t_soisno(c,j)  ! soil temperature (Kelvin)
        end do
         this%fates(nc)%bc_in(s)%forc_pbot           = forc_pbot(c)   ! atmospheric pressure (Pa)

         do ifp = 1,this%fates(nc)%sites(s)%youngest_patch%patchno
            p = ifp+col%patchi(c)
            ! Check to see if this patch is in the filter
            ! Note that this filter is most likely changing size, and getting smaller
            ! and smaller as more patch have converged on solution

            if( any(filterp==p) )then

               ! This filter is flushed to 1 before the canopyflux stability iterator
               ! It is set to status 2 if it is an active patch within the iterative loop
               ! After photosynthesis is called, it is upgraded to 3 if it was called.
               ! After all iterations we can evaluate which patches have a final flag
               ! of 3 to check if we missed any.
               this%fates(nc)%bc_in(s)%filter_photo_pa(ifp) = 2

               this%fates(nc)%bc_in(s)%dayl_factor_pa(ifp) = dayl_factor(p) ! scalar (0-1) for daylength
               this%fates(nc)%bc_in(s)%esat_tv_pa(ifp)     = esat_tv(p)     ! saturation vapor pressure at t_veg (Pa)
               this%fates(nc)%bc_in(s)%eair_pa(ifp)        = eair(p)        ! vapor pressure of canopy air (Pa)
               this%fates(nc)%bc_in(s)%oair_pa(ifp)        = oair(p)        ! Atmospheric O2 partial pressure (Pa)
               this%fates(nc)%bc_in(s)%cair_pa(ifp)        = cair(p)        ! Atmospheric CO2 partial pressure (Pa)
               this%fates(nc)%bc_in(s)%rb_pa(ifp)          = rb(p)          ! boundary layer resistance (s/m)
               this%fates(nc)%bc_in(s)%t_veg_pa(ifp)       = t_veg(p)       ! vegetation temperature (Kelvin)
               this%fates(nc)%bc_in(s)%tgcm_pa(ifp)        = tgcm(p)        ! air temperature at agcm reference height (kelvin)
            end if
         end do
      end do

      dtime = get_step_size_real()

      ! Call photosynthesis

      call FatesPlantRespPhotosynthDrive (this%fates(nc)%nsites, &
                                 this%fates(nc)%sites,  &
                                this%fates(nc)%bc_in,  &
                                this%fates(nc)%bc_out, &
                                dtime)

      ! Perform a double check to see if all patches on naturally vegetated columns
      ! were activated for photosynthesis
      ! ---------------------------------------------------------------------------------
      do icp = 1,fn
         p = filterp(icp)
         c = patch%column(p)
         s = this%f2hmap(nc)%hsites(c)
         ! do if structure here and only pass natveg columns
         ifp = p-col%patchi(c)

         if(this%fates(nc)%bc_in(s)%filter_photo_pa(ifp) /= 2)then
            write(iulog,*) 'Not all patches on the natveg column in the photosynthesis'
            write(iulog,*) 'filter ran photosynthesis s p icp ifp ilter',s,p,icp,ifp
            call endrun(msg=errMsg(sourcefile, __LINE__))
         else
            this%fates(nc)%bc_in(s)%filter_photo_pa(ifp) = 3
            rssun(p) = this%fates(nc)%bc_out(s)%rssun_pa(ifp)
            rssha(p) = this%fates(nc)%bc_out(s)%rssha_pa(ifp)

            ! These fields are marked with a bad-value flag
            photosyns_inst%psnsun_patch(p)   = spval
            photosyns_inst%psnsha_patch(p)   = spval
         end if
      end do

    end associate

    call t_stopf('fates_psn')

 end subroutine wrap_photosynthesis

 ! ======================================================================================

 subroutine wrap_accumulatefluxes(this, nc, fn, filterp)

   ! !ARGUMENTS:
   class(hlm_fates_interface_type), intent(inout) :: this
   integer                , intent(in)            :: nc                   ! clump index
   integer                , intent(in)            :: fn                   ! size of pft filter
   integer                , intent(in)            :: filterp(fn)          ! pft filter

   integer                                        :: s,c,p,ifp,icp
   real(r8)                                       :: dtime

   call t_startf('fates_wrapaccfluxes')

    ! Run a check on the filter
    do icp = 1,fn
       p = filterp(icp)
       c = patch%column(p)
       s = this%f2hmap(nc)%hsites(c)
       ifp = p-col%patchi(c)
       if(this%fates(nc)%bc_in(s)%filter_photo_pa(ifp) /= 3)then
            write(iulog,*) 'Not all patches on the natveg column in the canopys'
            write(iulog,*) 'filter ran canopy fluxes s p icp ifp ilter',s,p,icp,ifp
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end do


    dtime = get_step_size_real()
    call  AccumulateFluxes_ED(this%fates(nc)%nsites,  &
                               this%fates(nc)%sites, &
                               this%fates(nc)%bc_in,  &
                               this%fates(nc)%bc_out, &
                               dtime)

   call t_stopf('fates_wrapaccfluxes')

 end subroutine wrap_accumulatefluxes

 ! ======================================================================================

 subroutine wrap_WoodProducts(this, bounds_clump, num_soilc, filter_soilc, &
                              c_products_inst, n_products_inst)

   ! !ARGUMENTS:
   class(hlm_fates_interface_type), intent(inout) :: this
   type(bounds_type)              , intent(in)    :: bounds_clump
   integer                        , intent(in)    :: num_soilc          ! size of column filter
   integer                        , intent(in)    :: filter_soilc(:)    ! column filter
   type(cn_products_type)         , intent(inout) :: c_products_inst
   type(cn_products_type)         , intent(inout) :: n_products_inst
   
   ! Locals
   integer                                        :: s,c,g,fc
   integer                                        :: ci                 ! Clump index
      
   ci = bounds_clump%clump_index
   
   ! Loop over columns
   do fc = 1, num_soilc
      
      c = filter_soilc(fc)
      g = col%gridcell(c)
      s = this%f2hmap(ci)%hsites(c)
     
      ! Shijie: Pass harvested wood products to CLM product pools
      c_products_inst%hrv_deadstem_to_prod10_grc(g) = &
           c_products_inst%hrv_deadstem_to_prod10_grc(g) + &
           this%fates(ci)%bc_out(s)%hrv_deadstemc_to_prod10c
      
      c_products_inst%hrv_deadstem_to_prod100_grc(g) = &
           c_products_inst%hrv_deadstem_to_prod100_grc(g) + &
           this%fates(ci)%bc_out(s)%hrv_deadstemc_to_prod100c

      ! If N cycling is on
      if(fates_parteh_mode == prt_cnp_flex_allom_hyp ) then
         
         !n_products_inst%hrv_deadstem_to_prod10_grc(g) = &
         !     n_products_inst%hrv_deadstem_to_prod10_grc(g) + &
         !     this%fates(ci)%bc_out(s)%hrv_deadstemc_to_prod10c
         
         !n_products_inst%hrv_deadstem_to_prod100_grc(g) = &
         !     n_products_inst%hrv_deadstem_to_prod100_grc(g) + &
         !     this%fates(ci)%bc_out(s)%hrv_deadstemc_to_prod100c
         
      end if

          
   end do
    
   return
 end subroutine wrap_WoodProducts
 
 ! ======================================================================================

 subroutine wrap_canopy_radiation(this, bounds_clump, nc, fcansno, surfalb_inst)


    ! Pass boundary conditions (zenith angle, ground albedo and snow-cover)
    ! to FATES, and have fates call normalized radiation scattering.
    ! These methods are normalized because it is the zenith of the next
    ! timestep, and because of atmospheric model coupling, we need to provide
    ! an albedo. This normalized solution will be scaled by the actual
    ! downwelling radiation during the wrap sunshade fraction call

   
    ! Arguments
    class(hlm_fates_interface_type), intent(inout) :: this
    type(bounds_type),  intent(in)             :: bounds_clump
    integer            , intent(in)            :: nc ! clump index
    real(r8)           , intent(in)            :: fcansno( bounds_clump%begp: )
    type(surfalb_type) , intent(inout)         :: surfalb_inst
    
    ! locals
    integer                                    :: s,c,p,ifp,g

    call t_startf('fates_wrapcanopyradiation')

    associate(&
         coszen_col   =>    surfalb_inst%coszen_col         , & !in
         albgrd_col   =>    surfalb_inst%albgrd_col         , & !in
         albgri_col   =>    surfalb_inst%albgri_col         , & !in
         albd         =>    surfalb_inst%albd_patch         , & !out
         albi         =>    surfalb_inst%albi_patch         , & !out
         fabd         =>    surfalb_inst%fabd_patch         , & !out
         fabi         =>    surfalb_inst%fabi_patch         , & !out
         ftdd         =>    surfalb_inst%ftdd_patch         , & !out
         ftid         =>    surfalb_inst%ftid_patch         , & !out
         ftii         =>    surfalb_inst%ftii_patch)            !out


    do s = 1, this%fates(nc)%nsites

       c = this%f2hmap(nc)%fcolumn(s)
       g = col%gridcell(c)

       this%fates(nc)%bc_in(s)%coszen = coszen_col(c)

       ! FATES only does short-wave radiation scattering when
       ! the zenith angle is positive. In the future, FATES
       ! may perform calculations on diffuse radiation during
       ! dawn and dusk when the sun is below the horizon, but not yet

       do ifp = 1,this%fates(nc)%sites(s)%youngest_patch%patchno
          p = ifp+col%patchi(c)
          this%fates(nc)%bc_in(s)%fcansno_pa(ifp) = fcansno(p)
       end do

       
       if(coszen_col(c) > 0._r8) then

          this%fates(nc)%bc_in(s)%albgr_dir_rb(:) = albgrd_col(c,:)
          this%fates(nc)%bc_in(s)%albgr_dif_rb(:) = albgri_col(c,:)
       else

          ! This will ensure a crash in FATES if it tries
          this%fates(nc)%bc_in(s)%albgr_dir_rb(:) = spval
          this%fates(nc)%bc_in(s)%albgr_dif_rb(:) = spval
          
       end if

    end do

    call FatesNormalizedCanopyRadiation( &
            this%fates(nc)%sites, &
            this%fates(nc)%bc_in,  &
            this%fates(nc)%bc_out)

    ! Pass FATES BC's back to HLM
    ! -----------------------------------------------------------------------------------
    do s = 1, this%fates(nc)%nsites

       c = this%f2hmap(nc)%fcolumn(s)
       
       do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno
          
          p = ifp+col%patchi(c)

          albd(p,:) = this%fates(nc)%bc_out(s)%albd_parb(ifp,:)
          albi(p,:) = this%fates(nc)%bc_out(s)%albi_parb(ifp,:)
          fabd(p,:) = this%fates(nc)%bc_out(s)%fabd_parb(ifp,:)
          fabi(p,:) = this%fates(nc)%bc_out(s)%fabi_parb(ifp,:)
          ftdd(p,:) = this%fates(nc)%bc_out(s)%ftdd_parb(ifp,:)
          ftid(p,:) = this%fates(nc)%bc_out(s)%ftid_parb(ifp,:)
          ftii(p,:) = this%fates(nc)%bc_out(s)%ftii_parb(ifp,:)
          
       end do
    end do
          
  end associate

  call t_stopf('fates_wrapcanopyradiation')

 end subroutine wrap_canopy_radiation

 ! ======================================================================================

 subroutine WrapGlobalSeedDispersal(this,is_restart_flag)

   ! Call mpi procedure to provide the global seed output distribution array to every gridcell.
   ! This could be conducted with a more sophisticated halo-type structure or distributed graph.

   use decompMod,              only : procinfo
   use spmdMod,                only : MPI_REAL8, MPI_SUM, mpicom
   use FatesDispersalMod,      only : lneighbors, neighbor_type, dispersal_type
   use FatesInterfaceTypesMod, only : numpft_fates => numpft

   ! Arguments
   class(hlm_fates_interface_type), intent(inout) :: this
   logical, optional                              :: is_restart_flag 

   ! Local
   integer :: numg ! total number of gridcells across all processors
   integer :: ier  ! error code
   integer :: g    ! gridcell index

   logical :: set_restart_flag ! local logical variable to pass to IsItDispersalTime
                               ! if optional is_restart_flag is true
#ifdef _OPENMP
   logical, external :: omp_in_parallel
#endif

   type (neighbor_type),  pointer :: neighbor

   ! Check to see if we are not in a threaded region.  Fail the run if this returns true.
#ifdef _OPENMP
   if (omp_in_parallel()) then
      call endrun(msg='clmfates interface error: MPI routine called within threaded region'//&
           errMsg(sourcefile, __LINE__))
   end if
#endif

   ! This should only be run once per day
   if(is_beg_curr_day()) then

      ! If WrapGlobalSeedDispersal is being called at the end a fates restart call,
      ! pass .false. to the set_dispersed_flag to avoid updating the 
      ! global dispersal date
      set_restart_flag = .true.
      if (present(is_restart_flag)) then
         if (is_restart_flag) set_restart_flag = .false.
      end if

      call t_startf('fates-seed-mpi_reduce')

      if (IsItDispersalTime(setdispersedflag=set_restart_flag)) then

         ! Re-initialize incoming seed buffer for this time step
         this%fates_seed%incoming_global(:,:) = 0._r8
         this%fates_seed%outgoing_global(:,:) = 0._r8

         ! Distribute outgoing seed data across all MPI tasks
         ! This method of grid cell communications is inefficient in that it creates global information
         ! across all processes.  A future update should communicate to the minimum set of neighbors.
         call MPI_Allgatherv(this%fates_seed%outgoing_local, procinfo%ncells*numpft_fates, MPI_REAL8, &
                             this%fates_seed%outgoing_global, this%fates_seed%ncells_array*numpft_fates, this%fates_seed%begg_array*numpft_fates, &
                             MPI_REAL8, mpicom, ier)
         if (ier /= 0) then
            call endrun(msg='clmfates interface error: MPI_Allgatherv failed'//&
                 errMsg(sourcefile, __LINE__))
         end if

         ! zero outgoing local for all gridcells outside threaded region now that we've passed them out
         this%fates_seed%outgoing_local(:,:) = 0._r8

         ! Calculate the current gridcell incoming seed for each gridcell index
         ! This should be conducted outside of a threaded region to provide access to
         ! the neighbor%gindex which might not be available in via the clumped index
         call get_proc_global(ng=numg)
         do g = 1, numg
            neighbor => lneighbors(g)%first_neighbor
            do while (associated(neighbor))
               ! This also applies the same neighborhood distribution scheme to all pfts
               ! This needs to have a per pft density probability value
               this%fates_seed%incoming_global(:,g) = this%fates_seed%incoming_global(:,g) + &
                                                    this%fates_seed%outgoing_global(:,neighbor%gindex) * &
                                                    neighbor%density_prob(:) / lneighbors(g)%neighbor_count
               neighbor => neighbor%next_neighbor
            end do
         end do
      end if
   end if
  
   call t_stopf('fates-seed-mpi_reduce')

 end subroutine WrapGlobalSeedDispersal

 ! ======================================================================================

 subroutine WrapUpdateFatesSeedInOut(this,bounds_clump)

    ! This subroutine passes the globally dispersed seed via WrapGlobalSeedDispersal, incoming_global
    ! to the fates local process seed_in site object.  It also resets the fates seed_out
    ! in preparation for fates to update the seeds being dispersed out.

    ! Arguments
    class(hlm_fates_interface_type), intent(inout) :: this
    type(bounds_type),  intent(in)                 :: bounds_clump

    integer  :: g                           ! global index of the host gridcell
    integer  :: c                           ! global index of the host column
    integer  :: s                           ! FATES site index
    integer  :: nc                          ! clump index

    call t_startf('fates-seed-disperse')

    nc = bounds_clump%clump_index

    ! Check that it is the beginning of the current dispersal time step
    if (IsItDispersalTime()) then
       do s = 1, this%fates(nc)%nsites
          c = this%f2hmap(nc)%fcolumn(s)
          g = col%gridcell(c)

             ! assuming equal area for all sites, seed_id_global in [kg/grid/day], seed_in in [kg/site/day]
             this%fates(nc)%sites(s)%seed_in(:) = this%fates_seed%incoming_global(:,g)
             this%fates(nc)%sites(s)%seed_out(:) = 0._r8  ! reset seed_out
       end do
    else
       ! if it is not the dispersing time, pass in zero
       this%fates(nc)%sites(s)%seed_in(:) = 0._r8
    end if

    call t_stopf('fates-seed-disperse')

 end subroutine WrapUpdateFatesSeedInOut

 ! ======================================================================================

 subroutine wrap_update_hifrq_hist(this, bounds_clump, &
                                   soilbiogeochem_carbonflux_inst,     &
                                   soilbiogeochem_carbonstate_inst)

    ! Arguments
    class(hlm_fates_interface_type), intent(inout)    :: this
    type(bounds_type), intent(in)                     :: bounds_clump
    type(soilbiogeochem_carbonflux_type), intent(in)  :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type), intent(in) :: soilbiogeochem_carbonstate_inst

    ! locals
    real(r8) :: dtime
    integer  :: s, c, nc

    call t_startf('fates_wrap_update_hifrq_hist')

    associate(&
        hr            => soilbiogeochem_carbonflux_inst%hr_col,       & ! (gC/m2/s) total heterotrophic respiration
        totsomc       => soilbiogeochem_carbonstate_inst%totsomc_col, & ! (gC/m2) total soil organic matter carbon
        totlitc       => soilbiogeochem_carbonstate_inst%totlitc_col)   ! (gC/m2) total litter carbon in BGC pools

      nc = bounds_clump%clump_index

      if ( decomp_method /= no_soil_decomp )then
         ! Summarize Net Fluxes
         do s = 1, this%fates(nc)%nsites
            c = this%f2hmap(nc)%fcolumn(s)
            this%fates(nc)%bc_in(s)%tot_het_resp = hr(c)
            this%fates(nc)%bc_in(s)%tot_somc     = totsomc(c)
            this%fates(nc)%bc_in(s)%tot_litc     = totlitc(c)
         end do
      else
         ! Summarize Net Fluxes
         do s = 1, this%fates(nc)%nsites
            this%fates(nc)%bc_in(s)%tot_het_resp = 0.0_r8
            this%fates(nc)%bc_in(s)%tot_somc     = 0.0_r8
            this%fates(nc)%bc_in(s)%tot_litc     = 0.0_r8
         end do
      end if

      dtime = get_step_size_real()

      ! Update history variables that track these variables
      call fates_hist%update_history_hifrq(nc, &
            this%fates(nc)%nsites,  &
            this%fates(nc)%sites,   &
            this%fates(nc)%bc_in,   &
            this%fates(nc)%bc_out,  &
            dtime)
      
    end associate

    call t_stopf('fates_wrap_update_hifrq_hist')

  end subroutine wrap_update_hifrq_hist

 ! ======================================================================================


 subroutine TransferZ0mDisp(this, bounds_clump, z0m_patch, displa_patch)

    ! Arguments
    class(hlm_fates_interface_type), intent(in) :: this
    type(bounds_type),intent(in)                :: bounds_clump
    real(r8), intent(inout) :: z0m_patch( bounds_clump%begp: )  ! patch momentum roughness length (m)
    real(r8), intent(inout) :: displa_patch( bounds_clump%begp: )  ! patch displacement height (m)

    ! Locals
    integer :: ci   ! Current clump index
    integer :: s    ! Site index
    integer :: c    ! Column index
    integer :: ifp  ! Fates patch index
    integer :: p    ! CLM patch index

    call t_startf('fates_transferz0disp')

    SHR_ASSERT_FL((ubound(z0m_patch, 1) == bounds_clump%endp), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(displa_patch, 1) == bounds_clump%endp), sourcefile, __LINE__)

    ci = bounds_clump%clump_index

    do s = 1, this%fates(ci)%nsites
       c = this%f2hmap(ci)%fcolumn(s)

       z0m_patch(col%patchi(c)+1:col%patchf(c)) = 0.0_r8
       displa_patch(col%patchi(c)+1:col%patchf(c)) = 0.0_r8

         do ifp = 1, this%fates(ci)%sites(s)%youngest_patch%patchno
          p = ifp+col%patchi(c)
          z0m_patch(p) = this%fates(ci)%bc_out(s)%z0m_pa(ifp)
          displa_patch(p) = this%fates(ci)%bc_out(s)%displa_pa(ifp)
       end do
    end do

    call t_stopf('fates_transferz0disp')

    return
 end subroutine TransferZ0mDisp

 !-----------------------------------------------------------------------

  subroutine InterpFileInputs(this, bounds)
    !
    ! !DESCRIPTION:
    ! Interpolate inputs from files
    !
    ! NOTE(wjs, 2016-02-23) Stuff done here could probably be done at the end of
    ! InitEachTimeStep, rather than in this separate routine, except for the
    ! fact that
    ! (currently) this Interp stuff is done with proc bounds rather thna clump
    ! bounds. I
    ! think that is needed so that you don't update a given stream multiple
    ! times. If we
    ! rework the handling of threading / clumps so that there is a separate
    ! object for
    ! each clump, then I think this problem would disappear - at which point we
    ! could
    ! remove this Interp routine, moving its body to the end of
    ! InitEachTimeStep.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(hlm_fates_interface_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InterpFileInputs'
    !-----------------------------------------------------------------------

    call t_startf('fates_interpfileinputs')

    call this%fates_fire_data_method%FireInterp(bounds)

    call t_stopf('fates_interpfileinputs')

  end subroutine InterpFileInputs

  !-----------------------------------------------------------------------
  subroutine Init2(this, bounds, NLFilename)
    !
    ! !DESCRIPTION:
    ! Initialization after subgrid weights are determined
    !
    ! This copy should only be called if use_fates is .true.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(hlm_fates_interface_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    character(len=*), intent(in) :: NLFilename ! namelist filename
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Init2'
    !-----------------------------------------------------------------------

    call t_startf('fates_init2')

    call this%fates_fire_data_method%FireInit(bounds, NLFilename)

    call t_stopf('fates_init2')

  end subroutine Init2

  !-----------------------------------------------------------------------
  subroutine InitAccBuffer(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialized any accumulation buffers needed for FATES
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(hlm_fates_interface_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitAccBuffer'
    !-----------------------------------------------------------------------

    call t_startf('fates_initaccbuff')

    call this%fates_fire_data_method%InitAccBuffer( bounds )

    call t_stopf('fates_initaccbuff')

  end subroutine InitAccBuffer


  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialized any accumulation variables needed for FATES
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(hlm_fates_interface_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitAccVars'
    !-----------------------------------------------------------------------

    call t_startf('fates_initaccvars')

    call this%fates_fire_data_method%InitAccVars( bounds )

    call t_stopf('fates_initaccvars')

  end subroutine InitAccVars

  !-----------------------------------------------------------------------

    subroutine UpdateAccVars(this, bounds_proc)
   
    !
    ! !DESCRIPTION:
    ! Update any accumulation variables needed for FATES
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(hlm_fates_interface_type), intent(inout) :: this
    type(bounds_type), intent(in)                  :: bounds_proc
    !
  
    character(len=*), parameter :: subname = 'UpdateAccVars'
    !-----------------------------------------------------------------------

    call t_startf('fates_updateaccvars')

    call this%fates_fire_data_method%UpdateAccVars( bounds_proc )
    
    call t_stopf('fates_updateaccvars')

  end subroutine UpdateAccVars
  
 ! ======================================================================================

  subroutine WrapUpdateFatesRmean(this, nc, temperature_inst)

    class(hlm_fates_interface_type), intent(inout) :: this
    type(temperature_type), intent(in)             :: temperature_inst

    ! !LOCAL VARIABLES:
    integer                     :: nc,s,c,p,ifp  ! indices and loop counters

    do s = 1, this%fates(nc)%nsites
       c = this%f2hmap(nc)%fcolumn(s)
       do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno
          p = ifp+col%patchi(c)
          this%fates(nc)%bc_in(s)%t_veg_pa(ifp) = temperature_inst%t_veg_patch(p)
       end do
    end do

    call UpdateFatesRMeansTStep(this%fates(nc)%sites,this%fates(nc)%bc_in,this%fates(nc)%bc_out)
    
  end subroutine WrapUpdateFatesRmean
  
 ! ======================================================================================
  
 subroutine init_history_io(this,bounds_proc)

   use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp

   use FatesConstantsMod, only : fates_short_string_length, fates_long_string_length
   use FatesIOVariableKindMod, only : site_r8, site_soil_r8, site_size_pft_r8
   use FatesIOVariableKindMod, only : site_size_r8, site_pft_r8, site_age_r8
   use FatesIOVariableKindMod, only : site_coage_r8, site_coage_pft_r8
   use FatesIOVariableKindMod, only : site_fuel_r8, site_cwdsc_r8, site_scag_r8
   use FatesIOVariableKindMod, only : site_scagpft_r8, site_agepft_r8
   use FatesIOVariableKindMod, only : site_can_r8, site_cnlf_r8, site_cnlfpft_r8
   use FatesIOVariableKindMod, only : site_height_r8, site_elem_r8, site_elpft_r8
   use FatesIOVariableKindMod, only : site_elcwd_r8, site_elage_r8, site_agefuel_r8
   use FatesIOVariableKindMod, only : site_cdpf_r8, site_cdsc_r8, site_clscpf_r8
   use FatesIOVariableKindMod, only : site_landuse_r8, site_lulu_r8
   use FatesIODimensionsMod, only : fates_bounds_type


   ! Arguments
   class(hlm_fates_interface_type), intent(inout) :: this
   type(bounds_type),intent(in)                   :: bounds_proc  ! Currently "proc"


   ! Locals
   type(bounds_type)                              :: bounds_clump
   integer :: nvar  ! number of IO variables found
   integer :: ivar  ! variable index 1:nvar
   integer :: nc    ! thread counter 1:nclumps
   integer :: nclumps ! number of threads on this proc
   integer :: s     ! FATES site index
   integer :: c     ! ALM/CLM column index
   character(len=fates_short_string_length) :: dim2name
   character(len=fates_long_string_length) :: ioname
   integer :: d_index, dk_index

   type(fates_bounds_type) :: fates_bounds
   type(fates_bounds_type) :: fates_clump

    call t_startf('fates_inithistoryio')

   ! This routine initializes the types of output variables
   ! not the variables themselves, just the types
   ! ---------------------------------------------------------------------------------

   nclumps = get_proc_clumps()

   ! ------------------------------------------------------------------------------------
   ! PART I: Set FATES DIMENSIONING INFORMATION
   !
   ! -------------------------------------------------------------------------------
   ! Those who wish add variables that require new dimensions, please
   ! see FATES: FatesHistoryInterfaceMod.F90.  Dimension types are defined at the top of the
   ! module, and a new explicitly named instance of that type should be created.
   ! With this new dimension, a new output type/kind can contain that dimension.
   ! A new type/kind can be added to the dim_kinds structure, which defines its members
   ! in created in init_dim_kinds_maps().  Make sure to increase the size of fates_num_dim_kinds.
   ! A type/kind of output is defined by the data type (ie r8,int,..)
   ! and the dimensions.  Keep in mind that 3D variables (or 4D if you include time)
   ! are not really supported in CLM/ALM right now.  There are ways around this
   ! limitations by creating combined dimensions, for instance the size+pft dimension
   ! "scpf"
   ! ------------------------------------------------------------------------------------

   call hlm_bounds_to_fates_bounds(bounds_proc, fates_bounds)

   call fates_hist%Init(nclumps, fates_bounds)

   ! Define the bounds on the first dimension for each thread
   !$OMP PARALLEL DO PRIVATE (nc,bounds_clump,fates_clump)
   do nc = 1,nclumps

      call get_clump_bounds(nc, bounds_clump)

      ! thread bounds for patch
      call hlm_bounds_to_fates_bounds(bounds_clump, fates_clump)
      call fates_hist%SetThreadBoundsEach(nc, fates_clump)
   end do
   !$OMP END PARALLEL DO

   ! ------------------------------------------------------------------------------------
   ! PART II: USE THE JUST DEFINED DIMENSIONS TO ASSEMBLE THE VALID IO TYPES
   ! INTERF-TODO: THESE CAN ALL BE EMBEDDED INTO A SUBROUTINE IN HISTORYIOMOD
   ! ------------------------------------------------------------------------------------
   call fates_hist%assemble_history_output_types()

   ! ------------------------------------------------------------------------------------
   ! PART III: DEFINE THE LIST OF OUTPUT VARIABLE OBJECTS, AND REGISTER THEM WITH THE
   ! HLM ACCORDING TO THEIR TYPES
   ! ------------------------------------------------------------------------------------
   call fates_hist%initialize_history_vars()
   nvar = fates_hist%num_history_vars()

   call CrossRefHistoryFields()
   
   do ivar = 1, nvar

      associate( vname    => fates_hist%hvars(ivar)%vname, &
                 vunits   => fates_hist%hvars(ivar)%units,   &
                 vlong    => fates_hist%hvars(ivar)%long, &
                 vdefault => fates_hist%hvars(ivar)%use_default, &
                 vavgflag => fates_hist%hvars(ivar)%avgflag)

        dk_index = fates_hist%hvars(ivar)%dim_kinds_index
        ioname = trim(fates_hist%dim_kinds(dk_index)%name)

        select case(trim(ioname))
        case(site_r8)
           call hist_addfld1d(fname=trim(vname),units=trim(vunits),         &
                              avgflag=trim(vavgflag),long_name=trim(vlong), &
                              ptr_col=fates_hist%hvars(ivar)%r81d,      &
                              default=trim(vdefault),                       &
                              set_lake=0._r8,set_urb=0._r8)

        case(site_soil_r8, site_size_pft_r8, site_size_r8, site_pft_r8, &
             site_age_r8, site_height_r8, site_coage_r8,site_coage_pft_r8, &
             site_fuel_r8, site_cwdsc_r8, site_clscpf_r8, &
             site_can_r8,site_cnlf_r8, site_cnlfpft_r8, site_scag_r8, &
             site_scagpft_r8, site_agepft_r8, site_elem_r8, site_elpft_r8, &
             site_elcwd_r8, site_elage_r8, site_agefuel_r8, &
             site_cdsc_r8, site_cdpf_r8, &
             site_landuse_r8, site_lulu_r8)


           d_index = fates_hist%dim_kinds(dk_index)%dim2_index
           dim2name = fates_hist%dim_bounds(d_index)%name
           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         &
                              type2d=trim(dim2name),                        &
                              avgflag=trim(vavgflag),long_name=trim(vlong), &
                              ptr_col=fates_hist%hvars(ivar)%r82d,    &
                              default=trim(vdefault))


        case default
           write(iulog,*) 'A FATES iotype was created that was not registerred'
           write(iulog,*) 'in CLM.:',trim(ioname)
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end select

      end associate

   end do

   call t_stopf('fates_inithistoryio')

 end subroutine init_history_io

 ! ======================================================================================

 subroutine init_soil_depths(this, nc)

    ! Input Arguments
    class(hlm_fates_interface_type), intent(inout) :: this
    integer,intent(in)                             :: nc   ! Clump

    ! Locals
    integer :: s  ! site index
    integer :: c  ! column index
    integer :: j  ! Depth index
    integer :: nlevsoil
    integer :: nlevdecomp

    call t_startf('fates_initsoildepths')

    do s = 1, this%fates(nc)%nsites
       c = this%f2hmap(nc)%fcolumn(s)
       nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil
       nlevdecomp = this%fates(nc)%bc_in(s)%nlevdecomp

       this%fates(nc)%bc_in(s)%zi_sisl(0:nlevsoil)    = col%zi(c,0:nlevsoil)
       this%fates(nc)%bc_in(s)%dz_sisl(1:nlevsoil)    = col%dz(c,1:nlevsoil)
       this%fates(nc)%bc_in(s)%z_sisl(1:nlevsoil)     = col%z(c,1:nlevsoil)
       this%fates(nc)%bc_in(s)%dz_decomp_sisl(1:nlevdecomp) = &
             dzsoi_decomp(1:nlevdecomp)

       do j=1,nlevsoil
          this%fates(nc)%bc_in(s)%decomp_id(j) = j
          ! Check to make sure that dz = dz_decomp_sisl when vertical soil dynamics
          ! are active
          if(abs(this%fates(nc)%bc_in(s)%dz_decomp_sisl(j)-this%fates(nc)%bc_in(s)%dz_sisl(j))>1.e-10_r8)then
             write(iulog,*) 'when vertical soil decomp dynamics are on'
             write(iulog,*) 'fates assumes that the decomposition depths equal the soil depths'
             write(iulog,*) 'layer: ',j
             write(iulog,*) 'dz_decomp_sisl(j): ',this%fates(nc)%bc_in(s)%dz_decomp_sisl(j)
             write(iulog,*) 'dz_sisl(j): ',this%fates(nc)%bc_in(s)%dz_sisl(j)
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       end do

    end do

    call t_stopf('fates_initsoildepths')

    return
 end subroutine init_soil_depths

 ! ======================================================================================

 subroutine ComputeRootSoilFlux(this, bounds_clump, num_filterc, filterc, &
       soilstate_inst, waterfluxbulk_inst)

    class(hlm_fates_interface_type), intent(inout) :: this
    type(bounds_type),intent(in)                   :: bounds_clump
    integer,intent(in)                             :: num_filterc
    integer,intent(in)                             :: filterc(num_filterc)
    type(soilstate_type), intent(inout)            :: soilstate_inst
    type(waterfluxbulk_type), intent(inout)            :: waterfluxbulk_inst

    ! locals
    integer :: s
    integer :: c
    integer :: l
    integer :: nc
    integer :: num_filter_fates
    integer :: nlevsoil


    if( .not. use_fates_planthydro ) return

    call t_startf('fates_rootsoilflux')

    nc = bounds_clump%clump_index

    ! Perform a check that the number of columns submitted to fates for
    ! root water sink is the same that was expected in the hydrology filter
    num_filter_fates = 0
    do s = 1,num_filterc
       l = col%landunit(filterc(s))
       if (lun%itype(l) == istsoil ) then
          num_filter_fates = num_filter_fates + 1
       end if
    end do

    if(num_filter_fates /= this%fates(nc)%nsites )then
       write(iulog,*) 'The HLM list of natural veg columns during root water transfer'
       write(iulog,*) 'is not the same size as the fates site list?'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    do s = 1, this%fates(nc)%nsites
       c = this%f2hmap(nc)%fcolumn(s)
       nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil

       ! This is the water removed from the soil layers by roots (or added)
       waterfluxbulk_inst%qflx_rootsoi_col(c,1:nlevsoil) = this%fates(nc)%bc_out(s)%qflx_soil2root_sisl(1:nlevsoil)

       ! This is the total amount of water transferred to surface runoff
       ! (this is generated potentially from supersaturating soils
       ! (currently this is unnecessary)
       ! waterflux_inst%qflx_drain_vr_col(c,1:nlevsoil) = this%fates(nc)%bc_out(s)%qflx_ro_sisl(1:nlevsoil)


    end do

    call t_stopf('fates_rootsoilflux')

 end subroutine ComputeRootSoilFlux

 ! ======================================================================================

 subroutine wrap_hydraulics_drive(this, bounds_clump, nc, &
                                 soilstate_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, waterfluxbulk_inst, &
                                 fn, filterp, solarabs_inst, energyflux_inst)


   class(hlm_fates_interface_type), intent(inout) :: this
   type(bounds_type),intent(in)                   :: bounds_clump
   integer,intent(in)                             :: nc
   integer, intent(in)                            :: fn
   integer, intent(in)                            :: filterp(fn)
   type(soilstate_type)    , intent(inout)        :: soilstate_inst
   type(waterstatebulk_type)   , intent(inout)        :: waterstatebulk_inst
   type(waterdiagnosticbulk_type)   , intent(inout)        :: waterdiagnosticbulk_inst
   type(waterfluxbulk_type)    , intent(inout)        :: waterfluxbulk_inst
   type(solarabs_type)     , intent(inout)        :: solarabs_inst
   type(energyflux_type)   , intent(inout)        :: energyflux_inst

    ! locals
   integer :: s
   integer :: c
   integer :: j
   integer :: ifp
   integer :: p
   integer :: f
   integer :: nlevsoil
   integer :: icp
   real(r8) :: dtime

   if ( .not.use_fates_planthydro ) return

   call t_startf('fates_wraphydrodriv')

   dtime = get_step_size_real()

   ! Prepare Input Boundary Conditions
   ! ------------------------------------------------------------------------------------

   do s = 1, this%fates(nc)%nsites

      c = this%f2hmap(nc)%fcolumn(s)

      nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil

      this%fates(nc)%bc_in(s)%smpmin_si                 = &
            soilstate_inst%smpmin_col(c)
      this%fates(nc)%bc_in(s)%watsat_sisl(1:nlevsoil)    = &
            soilstate_inst%watsat_col(c,1:nlevsoil)
      this%fates(nc)%bc_in(s)%watres_sisl(1:nlevsoil)    = &
            soilstate_inst%watres_col(c,1:nlevsoil)
      this%fates(nc)%bc_in(s)%sucsat_sisl(1:nlevsoil)     = &
            soilstate_inst%sucsat_col(c,1:nlevsoil)
      this%fates(nc)%bc_in(s)%bsw_sisl(1:nlevsoil)        = &
            soilstate_inst%bsw_col(c,1:nlevsoil)
      this%fates(nc)%bc_in(s)%h2o_liq_sisl(1:nlevsoil)    = &
            waterstatebulk_inst%h2osoi_liq_col(c,1:nlevsoil)
      this%fates(nc)%bc_in(s)%eff_porosity_sl(1:nlevsoil) = &
            soilstate_inst%eff_porosity_col(c,1:nlevsoil)

        do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno
         p = ifp+col%patchi(c)

         this%fates(nc)%bc_in(s)%swrad_net_pa(ifp) = solarabs_inst%fsa_patch(p)
         this%fates(nc)%bc_in(s)%lwrad_net_pa(ifp) = energyflux_inst%eflx_lwrad_net_patch(p)
        end do
   end do

   ! The exposed vegetation filter "filterp" dictates which patches
   ! had their transpiration updated during canopy_fluxes(). Patches
   ! not in the filter had been zero'd during prep_canopyfluxes().

   do f = 1,fn
      p = filterp(f)
      c = patch%column(p)
      s = this%f2hmap(nc)%hsites(c)
      ifp = p - col%patchi(c)

      this%fates(nc)%bc_in(s)%qflx_transp_pa(ifp) = waterfluxbulk_inst%qflx_tran_veg_patch(p)
   end do

   ! Call Fates Hydraulics
   ! ------------------------------------------------------------------------------------


   call hydraulics_drive(this%fates(nc)%nsites, &
            this%fates(nc)%sites,  &
            this%fates(nc)%bc_in,  &
            this%fates(nc)%bc_out, &
            dtime)

   ! Prepare Output Boundary Conditions
   ! ------------------------------------------------------------------------------------

   do s = 1, this%fates(nc)%nsites
      c = this%f2hmap(nc)%fcolumn(s)
      waterdiagnosticbulk_inst%total_plant_stored_h2o_col(c) = &
            this%fates(nc)%bc_out(s)%plant_stored_h2o_si
   end do

   call fates_hist%update_history_hydraulics(nc, &
         this%fates(nc)%nsites, &
         this%fates(nc)%sites, &
         this%fates(nc)%bc_in, &
         dtime)

   call t_stopf('fates_wraphydrodriv')

   return
 end subroutine wrap_hydraulics_drive

 ! ======================================================================================

 subroutine hlm_bounds_to_fates_bounds(hlm, fates)

   use FatesIODimensionsMod,   only : fates_bounds_type
   use FatesInterfaceTypesMod, only : nlevsclass, nlevage, nlevcoage
   use FatesInterfaceTypesMod, only : nlevheight
   use FatesInterfaceTypesMod, only : nlevdamage
   use FatesFuelClassesMod,    only : num_fuel_classes
   use FatesLitterMod,         only : ncwd
   use EDParamsMod,            only : nlevleaf, nclmax
   use FatesInterfaceTypesMod, only : numpft_fates => numpft
   use FatesConstantsMod, only : n_landuse_cats


   type(bounds_type), intent(in) :: hlm
   type(fates_bounds_type), intent(out) :: fates

   call t_startf('fates_hlm2fatesbnds')

   fates%cohort_begin = hlm%begcohort
   fates%cohort_end = hlm%endcohort

   fates%column_begin = hlm%begc
   fates%column_end = hlm%endc
   
   fates%soil_begin = 1
   fates%soil_end = nlevsoi

   fates%sizepft_class_begin = 1
   fates%sizepft_class_end = nlevsclass * numpft_fates

   fates%size_class_begin = 1
   fates%size_class_end = nlevsclass

   fates%coagepf_class_begin = 1
   fates%coagepf_class_end = nlevcoage * numpft_fates

   fates%coage_class_begin = 1
   fates%coage_class_end = nlevcoage

   fates%pft_class_begin = 1
   fates%pft_class_end = numpft_fates

   fates%age_class_begin = 1
   fates%age_class_end = nlevage

   fates%height_begin = 1
   fates%height_end = nlevheight

   fates%sizeage_class_begin = 1
   fates%sizeage_class_end   = nlevsclass * nlevage

   fates%agepft_class_begin = 1
   fates%agepft_class_end   = nlevage * numpft_fates

   fates%sizeagepft_class_begin = 1
   fates%sizeagepft_class_end   = nlevsclass * nlevage * numpft_fates

   fates%fuel_begin = 1
   fates%fuel_end = num_fuel_classes

   fates%cwdsc_begin = 1
   fates%cwdsc_end = ncwd

   fates%can_begin = 1
   fates%can_end = nclmax

   fates%cnlf_begin = 1
   fates%cnlf_end = nlevleaf * nclmax

   fates%cnlfpft_begin = 1
   fates%cnlfpft_end = nlevleaf * nclmax * numpft_fates

   fates%elem_begin = 1
   fates%elem_end   = num_elements

   fates%elpft_begin = 1
   fates%elpft_end   = num_elements * numpft_fates

   fates%elcwd_begin = 1
   fates%elcwd_end   = num_elements * ncwd

   fates%elage_begin = 1
   fates%elage_end   = num_elements * nlevage

   fates%agefuel_begin = 1
   fates%agefuel_end   = nlevage * num_fuel_classes

   fates%cdpf_begin = 1
   fates%cdpf_end = nlevdamage * numpft_fates * nlevsclass

   fates%cdsc_begin = 1
   fates%cdsc_end = nlevdamage * nlevsclass

   fates%cdam_begin = 1
   fates%cdam_end = nlevdamage

   fates%clscpf_begin = 1
   fates%clscpf_end = numpft_fates * nlevsclass * nclmax

   fates%landuse_begin = 1
   fates%landuse_end   = n_landuse_cats

   fates%lulu_begin = 1
   fates%lulu_end   = n_landuse_cats * n_landuse_cats
   
   call t_stopf('fates_hlm2fatesbnds')

 end subroutine hlm_bounds_to_fates_bounds

 ! ======================================================================================

 subroutine GetAndSetTime()

   ! CLM MODULES
   use clm_time_manager  , only : get_curr_days_per_year, &
                                  get_curr_date,     &
                                  get_ref_date,      &
                                  timemgr_datediff

   ! FATES MODULES
   use FatesInterfaceMod     , only : SetFatesTime

   ! LOCAL VARIABLES
   integer  :: yr                       ! year (0, ...)
   integer  :: mon                      ! month (1, ..., 12)
   integer  :: day                      ! day of month (1, ..., 31)
   integer  :: sec                      ! seconds of the day
   integer  :: current_year
   integer  :: current_month
   integer  :: current_day
   integer  :: current_tod
   integer  :: current_date
   integer  :: jan01_curr_year
   integer  :: reference_date
   integer  :: days_per_year
   real(r8) :: model_day
   real(r8) :: day_of_year

   call t_startf('fates_getandsettime')

   ! Get the current date and determine the set the start of the current year
   call get_curr_date(current_year,current_month,current_day,current_tod)
   current_date = current_year*10000 + current_month*100 + current_day
   jan01_curr_year = current_year*10000 + 100 + 1

   ! Get the reference date components and compute the date
   call get_ref_date(yr, mon, day, sec)
   reference_date = yr*10000 + mon*100 + day

   ! Get the defined number of days per year
   days_per_year = get_curr_days_per_year()

   ! Determine the model day
   call timemgr_datediff(reference_date, sec, current_date, current_tod, model_day)

   ! Determine the current DOY
   call timemgr_datediff(jan01_curr_year,0,current_date,sec,day_of_year)

   ! Set the FATES global time variables
   call SetFatesTime(current_year, current_month, &
                     current_day, current_tod, &
                     current_date, reference_date, &
                     model_day, floor(day_of_year), &
                     days_per_year, 1.0_r8/dble(days_per_year))

   call t_stopf('fates_getandsettime')

 end subroutine GetAndSetTime

 ! ======================================================================================

 subroutine GetLandusePFTData(bounds, landuse_pft_file, landuse_pft_map, landuse_bareground)

   ! !DESCRIPTION:
   ! Read in static landuse x pft file

   ! !USES:
   use fileutils , only : getfil
   use ncdio_pio , only : file_desc_t, ncd_io, ncd_inqdlen
   use ncdio_pio , only : ncd_pio_openfile, ncd_pio_closefile
   use decompMod , only : BOUNDS_LEVEL_PROC
   use clm_varcon, only : grlnd
   use FatesConstantsMod, only : fates_unset_r8


   ! !ARGUMENTS:
   type(bounds_type), intent(in)        :: bounds            ! proc-level bounds
   character(len=*) , intent(in)        :: landuse_pft_file  ! name of file containing static landuse x pft information
   real(r8), allocatable, intent(inout) :: landuse_pft_map(:,:,:)
   real(r8), allocatable, intent(inout) :: landuse_bareground(:)

   ! !LOCAL VARIABLES
   integer            :: varnum                    ! variable number
   integer            :: dimid, dimlen             ! dimension id number and length
   integer            :: ier                       ! error id
   character(len=256) :: locfn                     ! local file name
   type(file_desc_t)  :: ncid                      ! netcdf id
   real(r8), pointer  :: arraylocal(:,:)           ! local array for reading fraction data
   real(r8), pointer  :: arraylocal_bareground(:)  ! local array for reading bareground data
   logical            :: readvar                   ! true => variable is on dataset
   !character(len=16), parameter :: grlnd  = 'lndgrid'      ! name of lndgrid

   integer, parameter :: dim_landuse_pft = 14

   ! Land use name arrays
   character(len=10), parameter  :: landuse_pft_map_varnames(num_landuse_pft_vars) = &
                    [character(len=10)  :: 'frac_primr','frac_secnd','frac_pastr','frac_range'] !need to move 'frac_surf' to a different variable

   character(len=*), parameter :: subname = 'GetLandusePFTData'

   !-----------------------------------------------------------------------

   ! Check to see if the landuse file name has been provided
   ! Note: getfile checks this as well
   if (masterproc) then
      write(iulog,*) 'Attempting to read landuse x pft data .....'
      if (landuse_pft_file == ' ') then
         write(iulog,*)'landuse_pft_file must be specified'
         call endrun(msg=errMsg(__FILE__, __LINE__))
      endif
   endif

   ! Initialize the landuse x pft arrays and initialize to unset
    allocate(landuse_pft_map(bounds%begg:bounds%endg,dim_landuse_pft,num_landuse_pft_vars),stat=ier)
    if (ier /= 0) then
       call endrun(msg=' allocation error for landuse_pft_map'//errMsg(__FILE__, __LINE__))
    end if
    landuse_pft_map = fates_unset_r8

    allocate(landuse_bareground(bounds%begg:bounds%endg),stat=ier)
    if (ier /= 0) then
       call endrun(msg=' allocation error for landuse_bareground'//errMsg(__FILE__, __LINE__))
    end if
    landuse_bareground = fates_unset_r8


   ! Get the local filename and open the file
   call getfil(landuse_pft_file, locfn, 0)
   call ncd_pio_openfile (ncid, trim(locfn), 0)

   ! Check that natpft dimension on the file matches the target array dimensions
   call ncd_inqdlen(ncid, dimid, dimlen, 'natpft')
   if (dimlen /= dim_landuse_pft) then
      write(iulog,*) 'natpft dimensions on the landuse x pft file do not match target array size'
      call endrun(msg=errMsg(__FILE__, __LINE__))
   end if

   ! Allocate a temporary array since ncdio expects a pointer
   allocate(arraylocal(bounds%begg:bounds%endg,dim_landuse_pft))
   allocate(arraylocal_bareground(bounds%begg:bounds%endg))

   ! Read the landuse x pft data from file
   do varnum = 1, num_landuse_pft_vars
      call ncd_io(ncid=ncid, varname=landuse_pft_map_varnames(varnum), flag='read', &
                  data=arraylocal, dim1name=grlnd, readvar=readvar)
      if (.not. readvar) &
         call endrun(msg='ERROR: '//trim(landuse_pft_map_varnames(varnum))// &
                         ' NOT on landuse x pft file'//errMsg(__FILE__, __LINE__))
      landuse_pft_map(bounds%begg:bounds%endg,:,varnum) = arraylocal(bounds%begg:bounds%endg,:)
   end do

   ! Read the bareground data from file.  This is per gridcell only.
   call ncd_io(ncid=ncid, varname='frac_brgnd', flag='read', &
               data=arraylocal_bareground, dim1name=grlnd, readvar=readvar)
   if (.not. readvar) call endrun(msg='ERROR: frac_brgnd NOT on landuse x pft file'//errMsg(__FILE__, __LINE__))
   landuse_bareground(bounds%begg:bounds%endg) = arraylocal_bareground(bounds%begg:bounds%endg)

   ! Deallocate the temporary local array point and close the file
   deallocate(arraylocal)
   deallocate(arraylocal_bareground)
   call ncd_pio_closefile(ncid)

 end subroutine GetLandusePFTData


 !-----------------------------------------------------------------------

end module CLMFatesInterfaceMod
