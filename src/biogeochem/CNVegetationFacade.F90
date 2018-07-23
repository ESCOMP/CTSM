module CNVegetationFacade

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Facade for the CN Vegetation subsystem.
  !
  ! (A "facade", in software engineering terms, is a unified interface to a set of
  ! interfaces in a subsystem. The facade defines a higher-level interface that makes the
  ! subsystem easier to use.)
  !
  ! NOTE(wjs, 2016-02-19) I envision that we will introduce an abstract base class
  ! (VegBase). Then both CNVeg and EDVeg will extend VegBase. The rest of the CLM code can
  ! then have an instance of VegBase, which depending on the run, can be either a CNVeg or
  ! EDVeg instance.
  !
  ! In addition, we probably want an implementation when running without CN or fates - i.e.,
  ! an SPVeg inst. This would provide implementations for get_leafn_patch,
  ! get_downreg_patch, etc., so that we don't need to handle the non-cn case here (note
  ! that, currently, we return NaN for most of these getters, because these arrays are
  ! invalid and shouldn't be used when running in SP mode). Also, in its EcosystemDynamics
  ! routine, it would call SatellitePhenology (but note that the desired interface for
  ! EcosystemDynamics would be quite different... could just pass everything needed by any
  ! model, and ignore unneeded arguments). Then we can get rid of comments in this module
  ! like, "only call if use_cn is true", as well as use_cn conditionals in this module.
  !
  ! NOTE(wjs, 2016-02-23) Currently, SatellitePhenology is called even when running with
  ! CN, for the sake of dry deposition. This seems weird to me, and my gut feeling -
  ! without understanding it well - is that this should be rewritten to depend on LAI from
  ! CN rather than from satellite phenology. Until that is done, the separation between SP
  ! and other Veg modes will be messier.
  !
  ! NOTE(wjs, 2016-02-23) Currently, this class coordinates calls to soil BGC routines as
  ! well as veg BGC routines (even though it doesn't contain any soil BGC types). This is
  ! because CNDriver coordinates both the veg & soil BGC. We should probably split up
  ! CNDriver so that there is a cleaner separation between veg BGC and soil BGC, to allow
  ! easier swapping of (for example) CN and ED. At that point, this class could
  ! coordinate just the calls to veg BGC routines, with a similar facade class
  ! coordinating the calls to soil BGC routines.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use shr_infnan_mod                  , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod                     , only : errMsg => shr_log_errMsg
  use perf_mod                        , only : t_startf, t_stopf
  use decompMod                       , only : bounds_type
  use clm_varctl                      , only : iulog, use_cn, use_cndv, use_c13, use_c14
  use abortutils                      , only : endrun
  use spmdMod                         , only : masterproc
  use clm_time_manager                , only : get_curr_date, get_ref_date
  use clm_time_manager                , only : get_nstep, is_end_curr_year, is_first_step
  use CNBalanceCheckMod               , only : cn_balance_type
  use CNVegStateType                  , only : cnveg_state_type
  use CNVegCarbonFluxType             , only : cnveg_carbonflux_type
  use CNVegCarbonStateType            , only : cnveg_carbonstate_type
  use CNVegNitrogenFluxType           , only : cnveg_nitrogenflux_type
  use CNVegNitrogenStateType          , only : cnveg_nitrogenstate_type
  use CNFireMethodMod                 , only : cnfire_method_type
  use CNProductsMod                   , only : cn_products_type
  use NutrientCompetitionMethodMod    , only : nutrient_competition_method_type
  use SpeciesIsotopeType              , only : species_isotope_type
  use SpeciesNonIsotopeType           , only : species_non_isotope_type
  use CanopyStateType                 , only : canopystate_type
  use PhotosynthesisMod               , only : photosyns_type
  use atm2lndType                     , only : atm2lnd_type
  use WaterStateBulkType                  , only : waterstatebulk_type
  use WaterDiagnosticBulkType                  , only : waterdiagnosticbulk_type
  use WaterFluxBulkType                   , only : waterfluxbulk_type
  use SoilStateType                   , only : soilstate_type
  use TemperatureType                 , only : temperature_type 
  use CropType                        , only : crop_type
  use ch4Mod                          , only : ch4_type
  use CNDVType                        , only : dgvs_type
  use CNDVDriverMod                   , only : CNDVDriver, CNDVHIST
  use EnergyFluxType                  , only : energyflux_type
  use SaturatedExcessRunoffMod        , only : saturated_excess_runoff_type
  use FrictionVelocityMod             , only : frictionvel_type
  use SoilBiogeochemStateType         , only : soilBiogeochem_state_type
  use SoilBiogeochemCarbonStateType   , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemCarbonFluxType    , only : soilBiogeochem_carbonflux_type
  use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
  use CNFireEmissionsMod              , only : fireemis_type, CNFireEmisUpdate
  use CNDriverMod                     , only : CNDriverInit
  use CNDriverMod                     , only : CNDriverSummarizeStates, CNDriverSummarizeFluxes
  use CNDriverMod                     , only : CNDriverNoLeaching, CNDriverLeaching
  use CNCStateUpdate1Mod              , only : CStateUpdateDynPatch
  use CNNStateUpdate1Mod              , only : NStateUpdateDynPatch
  use CNVegStructUpdateMod            , only : CNVegStructUpdate
  use CNAnnualUpdateMod               , only : CNAnnualUpdate
  use dynConsBiogeochemMod            , only : dyn_cnbal_patch, dyn_cnbal_col
  use dynCNDVMod                      , only : dynCNDV_init, dynCNDV_interp
  use CNPrecisionControlMod           , only: CNPrecisionControl
  use SoilBiogeochemPrecisionControlMod , only: SoilBiogeochemPrecisionControl
  !
  implicit none
  private

  ! !PUBLIC TYPES:

  type, public :: cn_vegetation_type
     ! FIXME(bja, 2016-06) These need to be public for use when fates is
     ! turned on. Should either be moved out of here or create some ED
     ! version of the facade....
     type(cnveg_state_type)         :: cnveg_state_inst
     type(cnveg_carbonstate_type)   :: cnveg_carbonstate_inst
     type(cnveg_carbonflux_type)    :: cnveg_carbonflux_inst

     !X!private

     type(cnveg_carbonstate_type)   :: c13_cnveg_carbonstate_inst
     type(cnveg_carbonstate_type)   :: c14_cnveg_carbonstate_inst
     type(cnveg_carbonflux_type)    :: c13_cnveg_carbonflux_inst
     type(cnveg_carbonflux_type)    :: c14_cnveg_carbonflux_inst
     type(cnveg_nitrogenstate_type) :: cnveg_nitrogenstate_inst
     type(cnveg_nitrogenflux_type)  :: cnveg_nitrogenflux_inst

     type(cn_products_type)         :: c_products_inst
     type(cn_products_type)         :: c13_products_inst
     type(cn_products_type)         :: c14_products_inst
     type(cn_products_type)         :: n_products_inst

     type(cn_balance_type)          :: cn_balance_inst
     class(cnfire_method_type), allocatable :: cnfire_method
     type(dgvs_type)                :: dgvs_inst

     ! Control variables
     logical, private :: reseed_dead_plants    ! Flag to indicate if should reseed dead plants when starting up the model

     ! TODO(wjs, 2016-02-19) Evaluate whether some other variables should be moved in
     ! here. Whether they should be moved in depends on how tightly they are tied in with
     ! the other CN Vegetation stuff. A question to ask is: Is this module used when
     ! running with SP or ED? If so, then it should probably remain outside of CNVeg.
     !
     ! From the clm_instMod section on "CN vegetation types":
     ! - nutrient_competition_method
     !   - I'm pretty sure this should be moved into here; it's just a little messy to do
     !     so, because of how it's initialized (specifically, the call to readParameters
     !     in clm_initializeMod).
     !
     ! From the clm_instMod section on "general biogeochem types":
     ! - ch4_inst
     !   - probably not: really seems to belong in soilbiogeochem
     ! - crop_inst
     ! - dust_inst
     ! - vocemis_inst
     ! - fireemis_inst
     ! - drydepvel_inst
     
   contains
     procedure, public :: Init
     procedure, public :: InitAccBuffer
     procedure, public :: InitAccVars
     procedure, public :: UpdateAccVars
     procedure, public :: Restart

     procedure, public :: Init2                         ! Do initialization in initialize phase, after subgrid weights are determined
     procedure, public :: InitEachTimeStep              ! Do initializations at the start of each time step
     procedure, public :: InterpFileInputs              ! Interpolate inputs from files
     procedure, public :: UpdateSubgridWeights          ! Update subgrid weights if running with prognostic patch weights
     procedure, public :: DynamicAreaConservation       ! Conserve C & N with updates in subgrid weights
     procedure, public :: InitColumnBalance             ! Set the starting point for col-level balance checks
     procedure, public :: EcosystemDynamicsPreDrainage  ! Do the main science that needs to be done before hydrology-drainage
     procedure, public :: EcosystemDynamicsPostDrainage ! Do the main science that needs to be done after hydrology-drainage
     procedure, public :: BalanceCheck                  ! Check the carbon and nitrogen balance
     procedure, public :: EndOfTimeStepVegDynamics      ! Do vegetation dynamics that should be done at the end of each time step
     procedure, public :: WriteHistory                  ! Do any history writes that are specific to veg dynamics

     procedure, public :: get_net_carbon_exchange_grc   ! Get gridcell-level net carbon exchange array
     procedure, public :: get_leafn_patch               ! Get patch-level leaf nitrogen array
     procedure, public :: get_downreg_patch             ! Get patch-level downregulation array
     procedure, public :: get_root_respiration_patch    ! Get patch-level root respiration array
     procedure, public :: get_annsum_npp_patch          ! Get patch-level annual sum NPP array
     procedure, public :: get_agnpp_patch               ! Get patch-level aboveground NPP array
     procedure, public :: get_bgnpp_patch               ! Get patch-level belowground NPP array
     procedure, public :: get_froot_carbon_patch        ! Get patch-level fine root carbon array
     procedure, public :: get_croot_carbon_patch        ! Get patch-level coarse root carbon array
     procedure, public :: get_totvegc_col               ! Get column-level total vegetation carbon array

     procedure, private :: CNReadNML                    ! Read in the CN general namelist
  end type cn_vegetation_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds, NLFilename)
    !
    ! !DESCRIPTION:
    ! Initialize a CNVeg object.
    !
    ! Should be called regardless of whether use_cn is true
    !
    ! !USES:
    use CNFireFactoryMod , only : create_cnfire_method
    use clm_varcon       , only : c13ratio, c14ratio
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    character(len=*) , intent(in)    :: NLFilename ! namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp

    character(len=*), parameter :: subname = 'Init'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    ! Note - always initialize the memory for cnveg_state_inst (used in biogeophys/)
    call this%cnveg_state_inst%Init(bounds)
    
    if (use_cn) then

       ! Read in the general CN namelist
       call this%CNReadNML( NLFilename )    ! MUST be called first as passes down control information to others

       call this%cnveg_carbonstate_inst%Init(bounds, carbon_type='c12', ratio=1._r8, NLFilename=NLFilename)
       if (use_c13) then
          call this%c13_cnveg_carbonstate_inst%Init(bounds, carbon_type='c13', ratio=c13ratio, &
               NLFilename=NLFilename, c12_cnveg_carbonstate_inst=this%cnveg_carbonstate_inst)
       end if
       if (use_c14) then
          call this%c14_cnveg_carbonstate_inst%Init(bounds, carbon_type='c14', ratio=c14ratio, &
               NLFilename=NLFilename, c12_cnveg_carbonstate_inst=this%cnveg_carbonstate_inst)
       end if
       call this%cnveg_carbonflux_inst%Init(bounds, carbon_type='c12')
       if (use_c13) then
          call this%c13_cnveg_carbonflux_inst%Init(bounds, carbon_type='c13')
       end if
       if (use_c14) then
          call this%c14_cnveg_carbonflux_inst%Init(bounds, carbon_type='c14')
       end if
       call this%cnveg_nitrogenstate_inst%Init(bounds,                   &
            this%cnveg_carbonstate_inst%leafc_patch(begp:endp),          &
            this%cnveg_carbonstate_inst%leafc_storage_patch(begp:endp),  &
            this%cnveg_carbonstate_inst%frootc_patch(begp:endp),         &
            this%cnveg_carbonstate_inst%frootc_storage_patch(begp:endp), &
            this%cnveg_carbonstate_inst%deadstemc_patch(begp:endp) )
       call this%cnveg_nitrogenflux_inst%Init(bounds) 

       call this%c_products_inst%Init(bounds, species_non_isotope_type('C'))
       if (use_c13) then
          call this%c13_products_inst%Init(bounds, species_isotope_type('C', '13'))
       end if
       if (use_c14) then
          call this%c14_products_inst%Init(bounds, species_isotope_type('C', '14'))
       end if
       call this%n_products_inst%Init(bounds, species_non_isotope_type('N'))

       call this%cn_balance_inst%Init(bounds)

       ! Initialize the memory for the dgvs_inst data structure regardless of whether
       ! use_cndv is true so that it can be used in associate statements (nag compiler
       ! complains otherwise)
       call this%dgvs_inst%Init(bounds)
    end if

    allocate(this%cnfire_method, &
         source=create_cnfire_method(NLFilename))

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine CNReadNML( this, NLFilename )
    !
    ! !DESCRIPTION:
    ! Read in the general CN control namelist
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(inout) :: this
    character(len=*)         , intent(in)    :: NLFilename                 ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'CNReadNML'
    character(len=*), parameter :: nmlname = 'cn_general'   ! MUST match what is in namelist below
    !-----------------------------------------------------------------------
    logical :: reseed_dead_plants
    namelist /cn_general/ reseed_dead_plants

    reseed_dead_plants = this%reseed_dead_plants

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=cn_general, iostat=ierr)   ! Namelist name here MUST be the same as in nmlname above!
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (reseed_dead_plants      , mpicom)

    this%reseed_dead_plants = reseed_dead_plants

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=cn_general)    ! Name here MUST be the same as in nmlname above!
       write(iulog,*) ' '
    end if

    !-----------------------------------------------------------------------

  end subroutine CNReadNML


  !-----------------------------------------------------------------------
  subroutine InitAccBuffer(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for types contained here 
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitAccBuffer'
    !-----------------------------------------------------------------------

    if (use_cndv) then
       call this%dgvs_inst%InitAccBuffer(bounds)
    end if

  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize variables that are associated with accumulated fields
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitAccVars'
    !-----------------------------------------------------------------------

    if (use_cndv) then
       call this%dgvs_inst%initAccVars(bounds)
    end if

  end subroutine InitAccVars

  !-----------------------------------------------------------------------
  subroutine UpdateAccVars(this, bounds, t_a10_patch, t_ref2m_patch)
    !
    ! !DESCRIPTION:
    ! Update accumulated variables
    !
    ! Should be called every time step
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    ! NOTE(wjs, 2016-02-23) These need to be pointers to agree with the interface of
    ! UpdateAccVars in CNDVType (they are pointers there as a workaround for a compiler
    ! bug).
    real(r8), pointer , intent(in)   :: t_a10_patch(:)      ! 10-day running mean of the 2 m temperature (K)
    real(r8), pointer , intent(in)   :: t_ref2m_patch(:)    ! 2 m height surface air temperature (K)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'UpdateAccVars'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(t_a10_patch) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(t_ref2m_patch) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))

    if (use_cndv) then
       call this%dgvs_inst%UpdateAccVars(bounds, &
            t_a10_patch = t_a10_patch, &
            t_ref2m_patch = t_ref2m_patch)
    end if

  end subroutine UpdateAccVars


  !-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Handle restart (read / write) for CNVeg
    !
    ! Should be called regardless of whether use_cn is true
    !
    ! !USES:
    use ncdio_pio, only : file_desc_t
    use clm_varcon, only : c3_r2, c14ratio
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    integer  :: reseed_patch(bounds%endp-bounds%begp+1)
    integer  :: num_reseed_patch
    !
    ! !LOCAL VARIABLES:

    integer :: begp, endp

    character(len=*), parameter :: subname = 'Restart'
    !-----------------------------------------------------------------------

    if (use_cn) then
       begp = bounds%begp
       endp = bounds%endp
       call this%cnveg_carbonstate_inst%restart(bounds, ncid, flag=flag, carbon_type='c12', &
               reseed_dead_plants=this%reseed_dead_plants, filter_reseed_patch=reseed_patch, &
               num_reseed_patch=num_reseed_patch )
       if ( flag /= 'read' .and. num_reseed_patch /= 0 )then
          call endrun(msg="ERROR num_reseed should be zero and is not"//errmsg(sourcefile, __LINE__))
       end if
       if (use_c13) then
          call this%c13_cnveg_carbonstate_inst%restart(bounds, ncid, flag=flag, carbon_type='c13', &
               reseed_dead_plants=this%reseed_dead_plants, c12_cnveg_carbonstate_inst=this%cnveg_carbonstate_inst)
       end if
       if (use_c14) then
          call this%c14_cnveg_carbonstate_inst%restart(bounds, ncid, flag=flag, carbon_type='c14', &
               reseed_dead_plants=this%reseed_dead_plants, c12_cnveg_carbonstate_inst=this%cnveg_carbonstate_inst)
       end if

       call this%cnveg_carbonflux_inst%restart(bounds, ncid, flag=flag, carbon_type='c12')
       if (use_c13) then
          call this%c13_cnveg_carbonflux_inst%restart(bounds, ncid, flag=flag, carbon_type='c13')
       end if
       if (use_c14) then
          call this%c14_cnveg_carbonflux_inst%restart(bounds, ncid, flag=flag, carbon_type='c14')
       end if

       call this%cnveg_nitrogenstate_inst%restart(bounds, ncid, flag=flag,  &
            leafc_patch=this%cnveg_carbonstate_inst%leafc_patch(begp:endp),         &
            leafc_storage_patch=this%cnveg_carbonstate_inst%leafc_storage_patch(begp:endp), &
            frootc_patch=this%cnveg_carbonstate_inst%frootc_patch(begp:endp), &
            frootc_storage_patch=this%cnveg_carbonstate_inst%frootc_storage_patch(begp:endp), &
            deadstemc_patch=this%cnveg_carbonstate_inst%deadstemc_patch(begp:endp), &
            filter_reseed_patch=reseed_patch, num_reseed_patch=num_reseed_patch)
       call this%cnveg_nitrogenflux_inst%restart(bounds, ncid, flag=flag)
       call this%cnveg_state_inst%restart(bounds, ncid, flag=flag, &
            cnveg_carbonstate=this%cnveg_carbonstate_inst, &
            cnveg_nitrogenstate=this%cnveg_nitrogenstate_inst, &
            filter_reseed_patch=reseed_patch, num_reseed_patch=num_reseed_patch)

       call this%c_products_inst%restart(bounds, ncid, flag)
       if (use_c13) then
          call this%c13_products_inst%restart(bounds, ncid, flag, &
               template_for_missing_fields = this%c_products_inst, &
               template_multiplier = c3_r2)
       end if
       if (use_c14) then
          call this%c14_products_inst%restart(bounds, ncid, flag, &
               template_for_missing_fields = this%c_products_inst, &
               template_multiplier = c14ratio)
       end if
       call this%n_products_inst%restart(bounds, ncid, flag)

    end if

    if (use_cndv) then
       call this%dgvs_inst%Restart(bounds, ncid, flag=flag)
    end if

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine Init2(this, bounds, NLFilename)
    !
    ! !DESCRIPTION:
    ! Do initialization that is needed in the initialize phase, after subgrid weights are
    ! determined
    !
    ! Should only be called if use_cn is true
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type) , intent(inout) :: this
    type(bounds_type) , intent(in)    :: bounds
    character(len=*)  , intent(in)    :: NLFilename ! namelist filename
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Init2'
    !-----------------------------------------------------------------------

    call CNDriverInit(bounds, NLFilename, this%cnfire_method)

    if (use_cndv) then
       call dynCNDV_init(bounds, this%dgvs_inst)
    end if

  end subroutine Init2


  !-----------------------------------------------------------------------
  subroutine InitEachTimeStep(this, bounds, num_soilc, filter_soilc)
    !
    ! !DESCRIPTION:
    ! Do initializations that need to be done at the start of every time step
    !
    ! This includes zeroing fluxes
    !
    ! Should only be called if use_cn is true
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type) , intent(inout) :: this
    type(bounds_type) , intent(in)    :: bounds
    integer           , intent(in)    :: num_soilc       ! number of soil columns filter
    integer           , intent(in)    :: filter_soilc(:) ! filter for soil columns
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitEachTimeStep'
    !-----------------------------------------------------------------------

    call this%cnveg_carbonflux_inst%ZeroDWT(bounds)
    if (use_c13) then
       call this%c13_cnveg_carbonflux_inst%ZeroDWT(bounds)
    end if
    if (use_c14) then
       call this%c14_cnveg_carbonflux_inst%ZeroDWT(bounds)
    end if
    call this%cnveg_nitrogenflux_inst%ZeroDWT(bounds)
    call this%cnveg_carbonstate_inst%ZeroDWT(bounds)
    call this%cnveg_nitrogenstate_inst%ZeroDWT(bounds)

  end subroutine InitEachTimeStep

  !-----------------------------------------------------------------------
  subroutine InterpFileInputs(this, bounds)
    !
    ! !DESCRIPTION:
    ! Interpolate inputs from files
    !
    ! NOTE(wjs, 2016-02-23) Stuff done here could probably be done at the end of
    ! InitEachTimeStep, rather than in this separate routine, except for the fact that
    ! (currently) this Interp stuff is done with proc bounds rather thna clump bounds. I
    ! think that is needed so that you don't update a given stream multiple times. If we
    ! rework the handling of threading / clumps so that there is a separate object for
    ! each clump, then I think this problem would disappear - at which point we could
    ! remove this Interp routine, moving its body to the end of InitEachTimeStep.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type) , intent(inout) :: this
    type(bounds_type) , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InterpFileInputs'
    !-----------------------------------------------------------------------

    call this%cnfire_method%CNFireInterp(bounds)

  end subroutine InterpFileInputs


  !-----------------------------------------------------------------------
  subroutine UpdateSubgridWeights(this, bounds)
    !
    ! !DESCRIPTION:
    ! Update subgrid weights if running with prognostic patch weights
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type) , intent(inout) :: this
    type(bounds_type) , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'UpdateSubgridWeights'
    !-----------------------------------------------------------------------

    if (use_cndv) then
       call dynCNDV_interp(bounds, this%dgvs_inst)
    end if

  end subroutine UpdateSubgridWeights


  !-----------------------------------------------------------------------
  subroutine DynamicAreaConservation(this, bounds, clump_index, &
       num_soilp_with_inactive, filter_soilp_with_inactive, &
       num_soilc_with_inactive, filter_soilc_with_inactive, &
       prior_weights, patch_state_updater, column_state_updater, &
       canopystate_inst, photosyns_inst, &
       soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst, &
       c13_soilbiogeochem_carbonstate_inst, c14_soilbiogeochem_carbonstate_inst, &
       soilbiogeochem_nitrogenstate_inst, ch4_inst, soilbiogeochem_state_inst)
    !
    ! !DESCRIPTION:
    ! Conserve C & N with updates in subgrid weights
    !
    ! Should only be called if use_cn is true
    !
    ! !USES:
    use dynPriorWeightsMod      , only : prior_weights_type
    use dynPatchStateUpdaterMod, only : patch_state_updater_type
    use dynColumnStateUpdaterMod, only : column_state_updater_type
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(inout) :: this
    type(bounds_type)                       , intent(in)    :: bounds        

    ! Index of clump on which we're currently operating. Note that this implies that this
    ! routine must be called from within a clump loop.
    integer                                 , intent(in)    :: clump_index

    integer                                 , intent(in)    :: num_soilp_with_inactive ! number of points in filter_soilp_with_inactive
    integer                                 , intent(in)    :: filter_soilp_with_inactive(:) ! soil patch filter that includes inactive points
    integer                                 , intent(in)    :: num_soilc_with_inactive ! number of points in filter_soilc_with_inactive
    integer                                 , intent(in)    :: filter_soilc_with_inactive(:) ! soil column filter that includes inactive points
    type(prior_weights_type)                , intent(in)    :: prior_weights         ! weights prior to the subgrid weight updates
    type(patch_state_updater_type)          , intent(in)    :: patch_state_updater
    type(column_state_updater_type)         , intent(in)    :: column_state_updater
    type(canopystate_type)                  , intent(inout) :: canopystate_inst
    type(photosyns_type)                    , intent(inout) :: photosyns_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c13_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c14_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(ch4_type)                          , intent(inout) :: ch4_inst
    type(soilbiogeochem_state_type)         , intent(in)    :: soilbiogeochem_state_inst
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'DynamicAreaConservation'
    !-----------------------------------------------------------------------

    call t_startf('dyn_cnbal_patch')
    call dyn_cnbal_patch(bounds, &
         num_soilp_with_inactive, filter_soilp_with_inactive, &
         prior_weights, patch_state_updater, &
         canopystate_inst, photosyns_inst, &
         this%cnveg_state_inst, &
         this%cnveg_carbonstate_inst, this%c13_cnveg_carbonstate_inst, this%c14_cnveg_carbonstate_inst, &
         this%cnveg_carbonflux_inst, this%c13_cnveg_carbonflux_inst, this%c14_cnveg_carbonflux_inst, &
         this%cnveg_nitrogenstate_inst, this%cnveg_nitrogenflux_inst, &
         soilbiogeochem_carbonflux_inst, soilbiogeochem_state_inst)
    call t_stopf('dyn_cnbal_patch')

    ! It is important to update column-level state variables based on the fluxes
    ! generated by dyn_cnbal_patch (which handles the change in aboveground / patch-level
    ! C/N due to shrinking patches), before calling dyn_cnbal_col (which handles the
    ! change in belowground / column-level C/N due to changing column areas). This way,
    ! any aboveground biomass which is sent to litter or soil due to shrinking patch
    ! areas is accounted for by the column-level conservation. This is important if
    ! column weights on the grid cell are changing at the same time as patch weights on
    ! the grid cell (which will typically be the case when columns change in area).
    !
    ! The filters here need to include inactive points as well as active points so that
    ! we correctly update column states in columns that have just shrunk to 0 area -
    ! since those column states are still important in the following dyn_cnbal_col.
    call t_startf('CNUpdateDynPatch')
    call CStateUpdateDynPatch(bounds, num_soilc_with_inactive, filter_soilc_with_inactive, &
         this%cnveg_carbonflux_inst, this%cnveg_carbonstate_inst, &
         soilbiogeochem_carbonstate_inst)
    if (use_c13) then
       call CStateUpdateDynPatch(bounds, num_soilc_with_inactive, filter_soilc_with_inactive, &
            this%c13_cnveg_carbonflux_inst, this%c13_cnveg_carbonstate_inst, &
            soilbiogeochem_carbonstate_inst)
    end if
    if (use_c14) then
       call CStateUpdateDynPatch(bounds, num_soilc_with_inactive, filter_soilc_with_inactive, &
            this%c14_cnveg_carbonflux_inst, this%c14_cnveg_carbonstate_inst, &
            soilbiogeochem_carbonstate_inst)
    end if
    call NStateUpdateDynPatch(bounds, num_soilc_with_inactive, filter_soilc_with_inactive, &
         this%cnveg_nitrogenflux_inst, this%cnveg_nitrogenstate_inst, &
         soilbiogeochem_nitrogenstate_inst)
    call t_stopf('CNUpdateDynPatch')

    call t_startf('dyn_cnbal_col')
    call dyn_cnbal_col(bounds, clump_index, column_state_updater, &
         soilbiogeochem_carbonstate_inst, c13_soilbiogeochem_carbonstate_inst, &
         c14_soilbiogeochem_carbonstate_inst, soilbiogeochem_nitrogenstate_inst, &
         ch4_inst)
    call t_stopf('dyn_cnbal_col')

  end subroutine DynamicAreaConservation

  !-----------------------------------------------------------------------
  subroutine InitColumnBalance(this, bounds, num_allc, filter_allc, &
       num_soilc, filter_soilc, num_soilp, filter_soilp, &
       soilbiogeochem_carbonstate_inst, &
       c13_soilbiogeochem_carbonstate_inst, &
       c14_soilbiogeochem_carbonstate_inst, &
       soilbiogeochem_nitrogenstate_inst)
    !
    ! !DESCRIPTION:
    ! Set the starting point for column-level balance checks.
    !
    ! This should be called after DynamicAreaConservation, since the changes made by
    ! DynamicAreaConservation can break column-level conservation checks.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type)               , intent(inout) :: this
    type(bounds_type)                       , intent(in)    :: bounds  
    integer                                 , intent(in)    :: num_allc          ! number of columns in allc filter
    integer                                 , intent(in)    :: filter_allc(:)    ! filter for all active columns
    integer                                 , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                                 , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c13_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c14_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitColumnBalance'
    !-----------------------------------------------------------------------

    call CNDriverSummarizeStates(bounds, &
         num_allc, filter_allc, &
         num_soilc, filter_soilc, &
         num_soilp, filter_soilp, &
         this%cnveg_carbonstate_inst, &
         this%c13_cnveg_carbonstate_inst, &
         this%c14_cnveg_carbonstate_inst, &
         this%cnveg_nitrogenstate_inst, &
         soilbiogeochem_carbonstate_inst, &
         c13_soilbiogeochem_carbonstate_inst, &
         c14_soilbiogeochem_carbonstate_inst, &
         soilbiogeochem_nitrogenstate_inst)

    call this%cn_balance_inst%BeginCNBalance( &
         bounds, num_soilc, filter_soilc, &
         this%cnveg_carbonstate_inst, this%cnveg_nitrogenstate_inst)

  end subroutine InitColumnBalance


  !-----------------------------------------------------------------------
  subroutine EcosystemDynamicsPreDrainage(this, bounds, &
       num_soilc, filter_soilc, &
       num_soilp, filter_soilp, &
       num_pcropp, filter_pcropp, &
       doalb, &
       soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst,         &
       c13_soilbiogeochem_carbonflux_inst, c13_soilbiogeochem_carbonstate_inst, &
       c14_soilbiogeochem_carbonflux_inst, c14_soilbiogeochem_carbonstate_inst, &
       soilbiogeochem_state_inst,                                               &
       soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst,     &
       atm2lnd_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, waterfluxbulk_inst,                           &
       canopystate_inst, soilstate_inst, temperature_inst, crop_inst, ch4_inst, &
       photosyns_inst, saturated_excess_runoff_inst, energyflux_inst,          &
       nutrient_competition_method, fireemis_inst)
    !
    ! !DESCRIPTION:
    ! Do the main science for CN vegetation that needs to be done before hydrology-drainage
    !
    ! Should only be called if use_cn is true
    !
    ! !USES:

    !
    ! !ARGUMENTS:
    class(cn_vegetation_type)               , intent(inout) :: this
    type(bounds_type)                       , intent(in)    :: bounds  
    integer                                 , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                                 , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    integer                                 , intent(in)    :: num_pcropp        ! number of prog. crop patches in filter
    integer                                 , intent(in)    :: filter_pcropp(:)  ! filter for prognostic crop patches
    logical                                 , intent(in)    :: doalb             ! true = surface albedo calculation time step
    type(soilbiogeochem_state_type)         , intent(inout) :: soilbiogeochem_state_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: c13_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c13_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: c14_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c14_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(atm2lnd_type)                      , intent(in)    :: atm2lnd_inst 
    type(waterstatebulk_type)                   , intent(in)    :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)                   , intent(in)    :: waterdiagnosticbulk_inst
    type(waterfluxbulk_type)                    , intent(inout) :: waterfluxbulk_inst
    type(canopystate_type)                  , intent(inout) :: canopystate_inst
    type(soilstate_type)                    , intent(inout) :: soilstate_inst
    type(temperature_type)                  , intent(inout) :: temperature_inst
    type(crop_type)                         , intent(inout) :: crop_inst
    type(ch4_type)                          , intent(in)    :: ch4_inst
    type(photosyns_type)                    , intent(in)    :: photosyns_inst
    type(saturated_excess_runoff_type)      , intent(in)    :: saturated_excess_runoff_inst
    type(energyflux_type)                   , intent(in)    :: energyflux_inst
    class(nutrient_competition_method_type) , intent(inout) :: nutrient_competition_method
    type(fireemis_type)                     , intent(inout) :: fireemis_inst
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'EcosystemDynamicsPreDrainage'
    !-----------------------------------------------------------------------

    call crop_inst%CropIncrementYear(num_pcropp, filter_pcropp)

    call CNDriverNoLeaching(bounds,                                         &
         num_soilc, filter_soilc,                       &
         num_soilp, filter_soilp,                       &
         num_pcropp, filter_pcropp, doalb,              &
         this%cnveg_state_inst,                                                        &
         this%cnveg_carbonflux_inst, this%cnveg_carbonstate_inst,                           &
         this%c13_cnveg_carbonflux_inst, this%c13_cnveg_carbonstate_inst,                   &
         this%c14_cnveg_carbonflux_inst, this%c14_cnveg_carbonstate_inst,                   &
         this%cnveg_nitrogenflux_inst, this%cnveg_nitrogenstate_inst,                       &
         this%c_products_inst, this%c13_products_inst, this%c14_products_inst,    &
         this%n_products_inst,                                                    &
         soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst,         &
         c13_soilbiogeochem_carbonflux_inst, c13_soilbiogeochem_carbonstate_inst, &
         c14_soilbiogeochem_carbonflux_inst, c14_soilbiogeochem_carbonstate_inst, &
         soilbiogeochem_state_inst,                                               &
         soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst,     &
         atm2lnd_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, waterfluxbulk_inst,                           &
         canopystate_inst, soilstate_inst, temperature_inst, crop_inst, ch4_inst, &
         this%dgvs_inst, photosyns_inst, saturated_excess_runoff_inst, energyflux_inst,          &
         nutrient_competition_method, this%cnfire_method)

    ! fire carbon emissions 
    call CNFireEmisUpdate(bounds, num_soilp, filter_soilp, &
         this%cnveg_carbonflux_inst, this%cnveg_carbonstate_inst, fireemis_inst )

    call CNAnnualUpdate(bounds,            &
         num_soilc, filter_soilc, &
         num_soilp, filter_soilp, &
         this%cnveg_state_inst, this%cnveg_carbonflux_inst)

  end subroutine EcosystemDynamicsPreDrainage

  !-----------------------------------------------------------------------
  subroutine EcosystemDynamicsPostDrainage(this, bounds, num_allc, filter_allc, &
       num_soilc, filter_soilc, num_soilp, filter_soilp, doalb, crop_inst, &
       waterstatebulk_inst, waterdiagnosticbulk_inst, waterfluxbulk_inst, frictionvel_inst, canopystate_inst, &
       soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst, &
       c13_soilbiogeochem_carbonflux_inst, c13_soilbiogeochem_carbonstate_inst, &
       c14_soilbiogeochem_carbonflux_inst, c14_soilbiogeochem_carbonstate_inst, &
       soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst)
    !
    ! !DESCRIPTION:
    ! Do the main science for CN vegetation that needs to be done after hydrology-drainage
    !
    ! Should only be called if use_cn is true
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type)               , intent(inout) :: this
    type(bounds_type)                       , intent(in)    :: bounds  
    integer                                 , intent(in)    :: num_allc          ! number of columns in allc filter
    integer                                 , intent(in)    :: filter_allc(:)    ! filter for all active columns
    integer                                 , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                                 , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    logical                                 , intent(in)    :: doalb             ! true = surface albedo calculation time step
    type(crop_type)                         , intent(in)    :: crop_inst
    type(waterstatebulk_type)                   , intent(in)    :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)                   , intent(in)    :: waterdiagnosticbulk_inst
    type(waterfluxbulk_type)                    , intent(inout) :: waterfluxbulk_inst
    type(frictionvel_type)                  , intent(in)    :: frictionvel_inst
    type(canopystate_type)                  , intent(inout) :: canopystate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: c13_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c13_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: c14_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c14_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'EcosystemDynamicsPostDrainage'
    !-----------------------------------------------------------------------

    ! Update the nitrogen leaching rate as a function of soluble mineral N 
    ! and total soil water outflow.
    
    call CNDriverLeaching(bounds, &
         num_soilc, filter_soilc, &
         num_soilp, filter_soilp, &
         waterstatebulk_inst, waterfluxbulk_inst, &
         this%cnveg_nitrogenflux_inst, this%cnveg_nitrogenstate_inst, &
         soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst)

    ! Set controls on very low values in critical state variables 

    call t_startf('CNPrecisionControl')
    call CNPrecisionControl(bounds, num_soilp, filter_soilp, &
         this%cnveg_carbonstate_inst, this%c13_cnveg_carbonstate_inst, &
         this%c14_cnveg_carbonstate_inst, this%cnveg_nitrogenstate_inst)
    call t_stopf('CNPrecisionControl')

    call t_startf('SoilBiogeochemPrecisionControl')
    call SoilBiogeochemPrecisionControl(num_soilc, filter_soilc,  &
         soilbiogeochem_carbonstate_inst, c13_soilbiogeochem_carbonstate_inst, &
         c14_soilbiogeochem_carbonstate_inst,soilbiogeochem_nitrogenstate_inst)
    call t_stopf('SoilBiogeochemPrecisionControl')

    ! Call to all CN summary routines

    call  CNDriverSummarizeStates(bounds, &
         num_allc, filter_allc, &
         num_soilc, filter_soilc, &
         num_soilp, filter_soilp, &
         this%cnveg_carbonstate_inst, &
         this%c13_cnveg_carbonstate_inst, &
         this%c14_cnveg_carbonstate_inst, &
         this%cnveg_nitrogenstate_inst, &
         soilbiogeochem_carbonstate_inst, &
         c13_soilbiogeochem_carbonstate_inst, &
         c14_soilbiogeochem_carbonstate_inst, &
         soilbiogeochem_nitrogenstate_inst)

    call  CNDriverSummarizeFluxes(bounds, &
         num_soilc, filter_soilc, &
         num_soilp, filter_soilp, &
         this%cnveg_carbonflux_inst, &
         this%c13_cnveg_carbonflux_inst, &
         this%c14_cnveg_carbonflux_inst, &
         this%cnveg_nitrogenflux_inst, &
         this%c_products_inst, this%c13_products_inst, this%c14_products_inst, &
         soilbiogeochem_carbonflux_inst, &
         c13_soilbiogeochem_carbonflux_inst, &
         c14_soilbiogeochem_carbonflux_inst, &
         soilbiogeochem_nitrogenflux_inst)

    ! On the radiation time step, use C state variables to calculate
    ! vegetation structure (LAI, SAI, height)

    if (doalb) then   
       call CNVegStructUpdate(num_soilp, filter_soilp, &
            waterdiagnosticbulk_inst, frictionvel_inst, this%dgvs_inst, this%cnveg_state_inst, &
            crop_inst, this%cnveg_carbonstate_inst, canopystate_inst)
    end if

  end subroutine EcosystemDynamicsPostDrainage

  !-----------------------------------------------------------------------
  subroutine BalanceCheck(this, bounds, num_soilc, filter_soilc, &
       soilbiogeochem_carbonflux_inst, soilbiogeochem_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! Check the carbon and nitrogen balance
    !
    ! Should only be called if use_cn is true
    !
    ! !USES:
    use clm_time_manager   , only : get_nstep_since_startup_or_lastDA_restart_or_pause
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type)               , intent(inout) :: this
    type(bounds_type)                       , intent(in)    :: bounds  
    integer                                 , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer              :: DA_nstep                   ! time step number

    character(len=*), parameter :: subname = 'BalanceCheck'
    !-----------------------------------------------------------------------

    DA_nstep = get_nstep_since_startup_or_lastDA_restart_or_pause()
    if (DA_nstep < 2 )then
       if (masterproc) then
          write(iulog,*) '--WARNING-- skipping CN balance check for first timesteps after startup or data assimilation'
       end if
    else

       call this%cn_balance_inst%CBalanceCheck( &
            bounds, num_soilc, filter_soilc, &
            soilbiogeochem_carbonflux_inst, &
            this%cnveg_carbonflux_inst, this%cnveg_carbonstate_inst)

       call this%cn_balance_inst%NBalanceCheck( &
            bounds, num_soilc, filter_soilc, &
            soilbiogeochem_nitrogenflux_inst, &
            this%cnveg_nitrogenflux_inst, this%cnveg_nitrogenstate_inst)

    end if

  end subroutine BalanceCheck

  !-----------------------------------------------------------------------
  subroutine EndOfTimeStepVegDynamics(this, bounds, num_natvegp, filter_natvegp, &
       atm2lnd_inst)
    !
    ! !DESCRIPTION:
    ! Do vegetation dynamics that should be done at the end of each time step
    !
    ! Should only be called if use_cn is true
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(inout) :: this
    type(bounds_type)  , intent(in)    :: bounds                  
    integer            , intent(inout) :: num_natvegp       ! number of naturally-vegetated patches in filter
    integer            , intent(inout) :: filter_natvegp(:) ! filter for naturally-vegetated patches
    type(atm2lnd_type) , intent(inout) :: atm2lnd_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: nstep  ! time step number
    integer  :: yr     ! year (0, ...)
    integer  :: mon    ! month (1, ..., 12)
    integer  :: day    ! day of month (1, ..., 31)
    integer  :: sec    ! seconds of the day
    integer  :: ncdate ! current date
    integer  :: nbdate ! base date (reference date)
    integer  :: kyr    ! thousand years, equals 2 at end of first year

    character(len=*), parameter :: subname = 'EndOfTimeStepVegDynamics'
    !-----------------------------------------------------------------------

    if (use_cndv) then
       ! Call dv (dynamic vegetation) at last time step of year

       call t_startf('d2dgvm')
       if (is_end_curr_year() .and. .not. is_first_step())  then

          ! Get date info.  kyr is used in lpj().  At end of first year, kyr = 2.
          call get_curr_date(yr, mon, day, sec)
          ncdate = yr*10000 + mon*100 + day
          call get_ref_date(yr, mon, day, sec)
          nbdate = yr*10000 + mon*100 + day
          kyr = ncdate/10000 - nbdate/10000 + 1

          if (masterproc) then
             nstep = get_nstep()
             write(iulog,*) 'End of year. CNDV called now: ncdate=', &
                  ncdate,' nbdate=',nbdate,' kyr=',kyr,' nstep=', nstep
          end if

          call CNDVDriver(bounds, &
               num_natvegp, filter_natvegp, kyr,  &
               atm2lnd_inst, &
               this%cnveg_carbonflux_inst, this%cnveg_carbonstate_inst, this%dgvs_inst)
       end if
       call t_stopf('d2dgvm')
    end if

  end subroutine EndOfTimeStepVegDynamics

  !-----------------------------------------------------------------------
  subroutine WriteHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Do any history writes that are specific to vegetation dynamics
    !
    ! NOTE(wjs, 2016-02-23) This could probably be combined with
    ! EndOfTimeStepVegDynamics, except for the fact that (currently) history writes are
    ! done with proc bounds rather than clump bounds. If that were changed, then the body
    ! of this could be moved into EndOfTimeStepVegDynamics, inside a "if (.not.
    ! use_noio)" conditional.
    !
    ! Should only be called if use_cn is true
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(in) :: this
    type(bounds_type)  , intent(in) :: bounds                  
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'WriteHistory'
    !-----------------------------------------------------------------------

    ! Write to CNDV history buffer if appropriate
    if (use_cndv) then
       if (is_end_curr_year() .and. .not. is_first_step())  then
          call t_startf('clm_drv_io_hdgvm')
          call CNDVHist( bounds, this%dgvs_inst )
          if (masterproc) write(iulog,*) 'Annual CNDV calculations are complete'
          call t_stopf('clm_drv_io_hdgvm')
       end if
    end if

  end subroutine WriteHistory


  !-----------------------------------------------------------------------
  function get_net_carbon_exchange_grc(this, bounds) result(net_carbon_exchange_grc)
    !
    ! !DESCRIPTION:
    ! Get gridcell-level net carbon exchange array
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8) :: net_carbon_exchange_grc(bounds%begg:bounds%endg)  ! function result: net carbon exchange between land and atmosphere, includes fire, landuse, harvest and hrv_xsmrpool flux, positive for source (gC/m2/s)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_net_carbon_exchange_grc'
    !-----------------------------------------------------------------------

    if (use_cn) then
       net_carbon_exchange_grc(bounds%begg:bounds%endg) = &
            -this%cnveg_carbonflux_inst%nbp_grc(bounds%begg:bounds%endg)
    else
       net_carbon_exchange_grc(bounds%begg:bounds%endg) = 0._r8
    end if

  end function get_net_carbon_exchange_grc


  !-----------------------------------------------------------------------
  function get_leafn_patch(this, bounds) result(leafn_patch)
    !
    ! !DESCRIPTION:
    ! Get patch-level leaf nitrogen array
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8) :: leafn_patch(bounds%begp:bounds%endp)  ! function result: leaf N (gN/m2)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_leafn_patch'
    !-----------------------------------------------------------------------

    if (use_cn) then
       leafn_patch(bounds%begp:bounds%endp) = &
            this%cnveg_nitrogenstate_inst%leafn_patch(bounds%begp:bounds%endp)
    else
       leafn_patch(bounds%begp:bounds%endp) = nan
    end if

  end function get_leafn_patch

  !-----------------------------------------------------------------------
  function get_downreg_patch(this, bounds) result(downreg_patch)
    !
    ! !DESCRIPTION:
    ! Get patch-level downregulation array
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8) :: downreg_patch(bounds%begp:bounds%endp)  ! function result: fractional reduction in GPP due to N limitation (dimensionless)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_downreg_patch'
    !-----------------------------------------------------------------------

    if (use_cn) then
       downreg_patch(bounds%begp:bounds%endp) = &
            this%cnveg_state_inst%downreg_patch(bounds%begp:bounds%endp)
    else
       downreg_patch(bounds%begp:bounds%endp) = nan
    end if

  end function get_downreg_patch

  !-----------------------------------------------------------------------
  function get_root_respiration_patch(this, bounds) result(root_respiration_patch)
    !
    ! !DESCRIPTION:
    ! Get patch-level root respiration array
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8) :: root_respiration_patch(bounds%begp:bounds%endp)  ! function result: root respiration (fine root MR + total root GR) (gC/m2/s)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_root_respiration_patch'
    !-----------------------------------------------------------------------

    if (use_cn) then
       root_respiration_patch(bounds%begp:bounds%endp) = &
            this%cnveg_carbonflux_inst%rr_patch(bounds%begp:bounds%endp)
    else
       root_respiration_patch(bounds%begp:bounds%endp) = nan
    end if

  end function get_root_respiration_patch

  ! TODO(wjs, 2016-02-19) annsum_npp, agnpp and bgnpp are all needed for the estimation
  ! of tillers in ch4Mod. Rather than providing getters for these three things so that
  ! ch4Mod can estimate tillers, it would probably be better if the tiller estimation
  ! algorithm was moved into some CNVeg-specific module, and then tillers could be
  ! queried directly.

  !-----------------------------------------------------------------------
  function get_annsum_npp_patch(this, bounds) result(annsum_npp_patch)
    !
    ! !DESCRIPTION:
    ! Get patch-level annual sum NPP array
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8) :: annsum_npp_patch(bounds%begp:bounds%endp)  ! function result: annual sum NPP (gC/m2/yr)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_annsum_npp_patch'
    !-----------------------------------------------------------------------

    if (use_cn) then
       annsum_npp_patch(bounds%begp:bounds%endp) = &
            this%cnveg_carbonflux_inst%annsum_npp_patch(bounds%begp:bounds%endp)
    else
       annsum_npp_patch(bounds%begp:bounds%endp) = nan
    end if

  end function get_annsum_npp_patch

  !-----------------------------------------------------------------------
  function get_agnpp_patch(this, bounds) result(agnpp_patch)
    !
    ! !DESCRIPTION:
    ! Get patch-level aboveground NPP array
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8) :: agnpp_patch(bounds%begp:bounds%endp)  ! function result: aboveground NPP (gC/m2/s)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_agnpp_patch'
    !-----------------------------------------------------------------------

    if (use_cn) then
       agnpp_patch(bounds%begp:bounds%endp) = &
            this%cnveg_carbonflux_inst%agnpp_patch(bounds%begp:bounds%endp)
    else
       agnpp_patch(bounds%begp:bounds%endp) = nan
    end if

  end function get_agnpp_patch

  !-----------------------------------------------------------------------
  function get_bgnpp_patch(this, bounds) result(bgnpp_patch)
    !
    ! !DESCRIPTION:
    ! Get patch-level belowground NPP array
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8) :: bgnpp_patch(bounds%begp:bounds%endp)  ! function result: belowground NPP (gC/m2/s)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_bgnpp_patch'
    !-----------------------------------------------------------------------

    if (use_cn) then
       bgnpp_patch(bounds%begp:bounds%endp) = &
            this%cnveg_carbonflux_inst%bgnpp_patch(bounds%begp:bounds%endp)
    else
       bgnpp_patch(bounds%begp:bounds%endp) = nan
    end if

  end function get_bgnpp_patch

  !-----------------------------------------------------------------------
  function get_froot_carbon_patch(this, bounds, tlai) result(froot_carbon_patch)
    !
    ! !DESCRIPTION:
    ! Get patch-level fine root carbon array
    !
    ! !USES:
    use pftconMod           , only : pftcon
    use PatchType           , only : patch
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8)         , intent(in) :: tlai( bounds%begp: )
    real(r8) :: froot_carbon_patch(bounds%begp:bounds%endp)  ! function result: (gC/m2)
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'get_froot_carbon_patch'
    integer :: p
    !-----------------------------------------------------------------------

    if (use_cn) then
       froot_carbon_patch(bounds%begp:bounds%endp) = &
            this%cnveg_carbonstate_inst%frootc_patch(bounds%begp:bounds%endp)
    else
! To get leaf biomass:
! bleaf = LAI / slatop
! g/m2 =  m2/m2  / m2/g 
! To get root biomass: 
! broot = bleaf * froot_leaf(ivt(p))
! g/m2 = g/m2 * g/g
       do p=bounds%begp, bounds%endp
          if (pftcon%slatop(patch%itype(p)) > 0._r8) then 
             froot_carbon_patch(p) = tlai(p) &
                  / pftcon%slatop(patch%itype(p)) &
                  *pftcon%froot_leaf(patch%itype(p))
          else
             froot_carbon_patch(p) = 0._r8
          endif
       enddo
    end if

  end function get_froot_carbon_patch

  !-----------------------------------------------------------------------
  function get_croot_carbon_patch(this, bounds, tlai) result(croot_carbon_patch)
    !
    ! !DESCRIPTION:
    ! Get patch-level live coarse root carbon array
    !
    ! !USES:
    use pftconMod           , only : pftcon
    use PatchType           , only : patch
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8)         , intent(in) :: tlai( bounds%begp: )
    real(r8) :: croot_carbon_patch(bounds%begp:bounds%endp)  ! function result: (gC/m2)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_croot_carbon_patch'
    integer :: p
    !-----------------------------------------------------------------------

    if (use_cn) then
       croot_carbon_patch(bounds%begp:bounds%endp) = &
            this%cnveg_carbonstate_inst%livecrootc_patch(bounds%begp:bounds%endp)
    else
! To get leaf biomass:
! bleaf = LAI / slatop
! g/m2 =  m2/m2  / m2/g 
! To get root biomass: 
! broot = bleaf * froot_leaf(ivt(p))
! g/m2 = g/m2 * g/g
       do p=bounds%begp, bounds%endp
          if (pftcon%slatop(patch%itype(p)) > 0._r8) then 
             croot_carbon_patch(p) = tlai(p) &
                  / pftcon%slatop(patch%itype(p)) &
                  *pftcon%stem_leaf(patch%itype(p)) &
                  *pftcon%croot_stem(patch%itype(p))
          else
             croot_carbon_patch(p) = 0._r8
          endif
       enddo
    end if

  end function get_croot_carbon_patch

  !-----------------------------------------------------------------------
  function get_totvegc_col(this, bounds) result(totvegc_col)
    !
    ! !DESCRIPTION:
    ! Get column-level total vegetation carbon array
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8) :: totvegc_col(bounds%begc:bounds%endc)  ! function result: (gC/m2)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_totvegc_col'
    !-----------------------------------------------------------------------

    if (use_cn) then
       totvegc_col(bounds%begc:bounds%endc) = &
            this%cnveg_carbonstate_inst%totvegc_col(bounds%begc:bounds%endc)
    else
       totvegc_col(bounds%begc:bounds%endc) = nan
    end if

  end function get_totvegc_col


end module CNVegetationFacade
