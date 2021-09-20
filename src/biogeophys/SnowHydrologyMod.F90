! -*- mode: f90; indent-tabs-mode: nil; f90-do-indent:3; f90-if-indent:3; f90-type-indent:3; f90-program-indent:2; f90-associate-indent:0; f90-continuation-indent:5  -*-
module SnowHydrologyMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate snow hydrology.
  ! - Using as input aerosol deposition from atmosphere model calculate
  !   aerosol fluxes and masses in each layer - need for surface albedo calculation
  ! - Change of snow mass and the snow water onto soil
  ! - Change in snow layer thickness due to compaction
  ! - Combine snow layers less than a min thickness
  ! - Subdivide snow layers if they exceed maximum thickness
  ! - Construct snow/no-snow filters

  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : bounds_type, subgrid_level_column
  use abortutils      , only : endrun
  use column_varcon   , only : icol_roof, icol_sunwall, icol_shadewall
  use clm_varpar      , only : nlevsno, nlevsoi, nlevgrnd, nlevmaxurbgrnd
  use clm_varctl      , only : iulog, use_subgrid_fluxes
  use clm_varcon      , only : h2osno_max, hfus, denh2o, denice, rpi, spval, tfrz
  use clm_varcon      , only : cpice, cpliq
  use atm2lndType     , only : atm2lnd_type
  use AerosolMod      , only : aerosol_type, AerosolFluxes
  use TemperatureType , only : temperature_type
  use WaterType       , only : water_type
  use WaterFluxBulkType   , only : waterfluxbulk_type
  use WaterStateBulkType  , only : waterstatebulk_type
  use WaterDiagnosticBulkType  , only : waterdiagnosticbulk_type
  use SnowCoverFractionBaseMod, only : snow_cover_fraction_base_type
  use LandunitType    , only : landunit_type, lun
  use TopoMod, only : topo_type
  use ColumnType      , only : column_type, col
  use landunit_varcon , only : istsoil, istdlak, istsoil, istwet, istice, istcrop
  use clm_time_manager, only : get_step_size_real, get_nstep
  use filterColMod    , only : filter_col_type, col_filter_from_filter_and_logical_array
  use LakeCon         , only : lsadz
  use NumericsMod     , only : truncate_small_values_one_lev
  use WaterTracerUtils, only : CalcTracerFromBulk, CalcTracerFromBulkMasked
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SnowHydrologyClean         ! Deallocate variables for unit testing
  public :: SnowHydrology_readnl       ! Read namelist
  public :: UpdateQuantitiesForNewSnow ! Update various snow-related quantities to account for new snow
  public :: RemoveSnowFromThawedWetlands ! Remove snow from thawed wetlands
  public :: InitializeExplicitSnowPack ! Initialize an explicit snow pack in columns where this is warranted based on snow depth
  public :: SnowWater                  ! Change of snow mass and the snow water onto soil
  public :: SnowCompaction             ! Change in snow layer thickness due to compaction
  public :: CombineSnowLayers          ! Combine snow layers less than a min thickness
  public :: DivideSnowLayers           ! Subdivide snow layers if they exceed maximum thickness
  public :: ZeroEmptySnowLayers        ! Set empty snow layers to zero
  public :: InitSnowLayers             ! Initialize cold-start snow layer thickness
  public :: BuildSnowFilter            ! Construct snow/no-snow filters
  public :: SnowCapping                ! Remove snow mass for capped columns
  public :: NewSnowBulkDensity         ! Compute bulk density of any newly-fallen snow
  public :: readParams                 ! Read in parameters on parameter file

  ! The following are public just for the sake of unit testing:
  public :: SnowCappingExcess          ! Determine the excess snow that needs to be capped
  public :: SnowHydrologySetControlForTesting ! Set some of the control settings

  type, private :: params_type
      real(r8) :: wimp                  ! Water impremeable if porosity less than wimp (unitless)
      real(r8) :: ssi                   ! Irreducible water saturation of snow (unitless)
      real(r8) :: drift_gs              ! Wind drift compaction / grain size (fixed value for now) (unitless)
      real(r8) :: eta0_anderson         ! Viscosity coefficent from Anderson1976 (kg*s/m2)
      real(r8) :: eta0_vionnet          ! Viscosity coefficent from Vionnet2012 (kg*s/m2)
      real(r8) :: wind_snowcompact_fact ! Reference wind above which fresh snow density is (substantially) increased (m/s)
      real(r8) :: rho_max               ! Wind drift compaction / maximum density (kg/m3)
      real(r8) :: tau_ref               ! Wind drift compaction / reference time (48*3600) (s)
      real(r8) :: scvng_fct_mlt_sf      ! Scaling factor modifying scavenging factors for BC, OC, and dust species inclusion in meltwater (-)
      real(r8) :: ceta                  ! Overburden compaction constant (kg/m3)
  end type params_type
  type(params_type), private ::  params_inst

  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: BulkDiag_NewSnowDiagnostics ! Update various snow-related diagnostic quantities to account for new snow
  private :: UpdateState_AddNewSnow      ! Update h2osno_no_layers or h2osoi_ice based on new snow
  private :: BuildFilter_ThawedWetlandThinSnowpack ! Build a column-level filter of thawed wetland columns with a thin snowpack
  private :: UpdateState_RemoveSnowFromThawedWetlands ! For bulk or one tracer: remove snow from thawed wetlands, for state variables
  private :: Bulk_RemoveSnowFromThawedWetlands ! Remove snow from thawed wetlands, for bulk-only quantities
  private :: BuildFilter_SnowpackInitialized ! Build a column-level filter of columns where an explicit snow pack needs to be initialized
  private :: UpdateState_InitializeSnowPack ! For bulk or one tracer: initialize water state variables for columns in which an explicit snow pack is being newly initialized
  private :: Bulk_InitializeSnowPack ! Initialize an explicit snow pack in columns where this is warranted based on snow depth, for bulk-only quantities
  private :: UpdateState_TopLayerFluxes ! Update top layer of snow pack with various fluxes into and out of the top layer
  private :: BulkFlux_SnowPercolation ! Calculate liquid percolation through the snow pack, for bulk water
  private :: TracerFlux_SnowPercolation ! Calculate liquid percolation through the snow pack, for one tracer
  private :: UpdateState_SnowPercolation ! Update h2osoi_liq for snow percolation, for bulk or one tracer
  private :: CalcAndApplyAerosolFluxes ! Calculate and apply fluxes of aerosols through the snow pack
  private :: PostPercolation_AdjustLayerThicknesses ! Adjust layer thickness for any water+ice content changes after percolation through the snow pack
  private :: BulkDiag_SnowWaterAccumulatedSnow ! Update int_snow, and reset accumulated snow when no snow present
  private :: SumFlux_AddSnowPercolation ! Calculate summed fluxes accounting for qflx_snow_percolation and similar fluxes
  private :: InitFlux_SnowCapping ! Initialize snow capping fluxes to 0
  private :: BulkFlux_SnowCappingFluxes ! Calculate snow capping fluxes and related terms for bulk water
  private :: TracerFlux_SnowCappingFluxes ! Calculate snow capping fluxes and related terms for one tracer
  private :: UpdateState_RemoveSnowCappingFluxes ! Remove snow capping fluxes from h2osoi_ice and h2osoi_liq
  private :: SnowCappingUpdateDzAndAerosols ! Following snow capping, adjust dz and aerosol masses in bottom snow layer
  private :: Combo                  ! Returns the combined variables: dz, t, wliq, wice.
  private :: MassWeightedSnowRadius ! Mass weighted snow grain size
  !
  ! !PUBLIC DATA MEMBERS:
  !  Aerosol species indices:
  !  1= hydrophillic black carbon
  !  2= hydrophobic black carbon
  !  3= hydrophilic organic carbon
  !  4= hydrophobic organic carbon
  !  5= dust species 1
  !  6= dust species 2
  !  7= dust species 3
  !  8= dust species 4
  !
  real(r8), public, parameter :: scvng_fct_mlt_bcphi = 0.20_r8 ! scavenging factor for hydrophillic BC inclusion in meltwater [frc]
  real(r8), public, parameter :: scvng_fct_mlt_bcpho = 0.03_r8 ! scavenging factor for hydrophobic BC inclusion in meltwater  [frc]
  real(r8), public, parameter :: scvng_fct_mlt_ocphi = 0.20_r8 ! scavenging factor for hydrophillic OC inclusion in meltwater [frc]
  real(r8), public, parameter :: scvng_fct_mlt_ocpho = 0.03_r8 ! scavenging factor for hydrophobic OC inclusion in meltwater  [frc]
  real(r8), public, parameter :: scvng_fct_mlt_dst1  = 0.02_r8 ! scavenging factor for dust species 1 inclusion in meltwater  [frc]
  real(r8), public, parameter :: scvng_fct_mlt_dst2  = 0.02_r8 ! scavenging factor for dust species 2 inclusion in meltwater  [frc]
  real(r8), public, parameter :: scvng_fct_mlt_dst3  = 0.01_r8 ! scavenging factor for dust species 3 inclusion in meltwater  [frc]
  real(r8), public, parameter :: scvng_fct_mlt_dst4  = 0.01_r8 ! scavenging factor for dust species 4 inclusion in meltwater  [frc]

  ! The following are public for the sake of unit testing
  integer, parameter, public :: LoTmpDnsSlater2017            = 2    ! For temperature below -15C use equation from Slater 2017
  integer, parameter, public :: LoTmpDnsTruncatedAnderson1976 = 1    ! Truncate low temp. snow density from the Anderson-1976 version at -15C

  !
  ! !PRIVATE DATA MEMBERS:

  integer, parameter :: OverburdenCompactionMethodAnderson1976 = 1
  integer, parameter :: OverburdenCompactionMethodVionnet2012  = 2

  ! If true, the density of new snow depends on wind speed, and there is also
  ! wind-dependent snow compaction
  logical  :: wind_dependent_snow_density                      ! If snow density depends on wind or not
  real(r8) :: snow_dzmin_1, snow_dzmax_l_1, snow_dzmax_u_1  ! namelist-defined top snow layer information
  real(r8) :: snow_dzmin_2, snow_dzmax_l_2, snow_dzmax_u_2  ! namelist-defined 2nd snow layer information
  real(r8), private, allocatable :: dzmin(:)  !min snow thickness of layer
  real(r8), private, allocatable :: dzmax_l(:)  !max snow thickness of layer when no layers beneath
  real(r8), private, allocatable :: dzmax_u(:)  !max snow thickness of layer when layers beneath

  integer  :: overburden_compaction_method = -1
  integer  :: new_snow_density            = LoTmpDnsSlater2017 ! Snow density type
  real(r8) :: upplim_destruct_metamorph   = 100.0_r8           ! Upper Limit on Destructive Metamorphism Compaction [kg/m3]
  real(r8) :: overburden_compress_Tfactor = 0.08_r8            ! snow compaction overburden exponential factor (1/K)

  ! ------------------------------------------------------------------------
  ! Parameters controlling the resetting of the snow pack
  ! ------------------------------------------------------------------------

  logical  :: reset_snow      = .false.  ! If set to true, we reset the non-glc snow pack, based on the following parameters
  logical  :: reset_snow_glc  = .false.  ! If set to true, we reset the glc snow pack, based reset_snow_glc_ela

  ! Default for reset_snow_glc_ela implies that snow will be reset for all glacier columns if reset_snow_glc = .true.
  real(r8) :: reset_snow_glc_ela = 1.e9_r8  ! equilibrium line altitude (m); snow is reset for glacier columns below this elevation if reset_snow_glc = .true. (ignored if reset_snow_glc = .false.)

  ! The following are public simply to support unit testing

  ! 35 mm was chosen by Raymond Sellevold, based on finding the location in Greenland
  ! with the least amount of snowfall from Sept. 1 (roughly the end of the melt season)
  ! and Jan. 1 (when we typically start simulations). This location with the least amount
  ! of snowfall had an average of 35 mm snow fall over this 4-month period.
  real(r8), parameter, public :: reset_snow_h2osno = 35._r8  ! mm SWE to reset the snow pack to

  ! We scale the number of reset time steps with the number of snow layers, since we can
  ! remove up to one layer per time step. In the absence of snow accumulation, we might
  ! be able to get away with 1 reset time step per layer. However, we specify a larger
  ! number to be more robust.
  real(r8), parameter, public :: reset_snow_timesteps_per_layer = 4

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SnowHydrology_readnl( NLFilename)
    !
    ! !DESCRIPTION:
    ! Read the namelist for SnowHydrology
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=64) :: snow_overburden_compaction_method
    character(len=25) :: lotmp_snowdensity_method

    character(len=*), parameter :: subname = 'SnowHydrology_readnl'
    !-----------------------------------------------------------------------

    namelist /clm_snowhydrology_inparm/ &
         wind_dependent_snow_density, snow_overburden_compaction_method, &
         lotmp_snowdensity_method, upplim_destruct_metamorph, &
         overburden_compress_Tfactor, &
         reset_snow, reset_snow_glc, reset_snow_glc_ela, &
         snow_dzmin_1, snow_dzmax_l_1, snow_dzmax_u_1, &
         snow_dzmin_2, snow_dzmax_l_2, snow_dzmax_u_2

    ! Initialize options to default values, in case they are not specified in the namelist
    wind_dependent_snow_density = .false.
    snow_overburden_compaction_method = ' '
    snow_dzmin_1 = nan
    snow_dzmin_2 = nan
    snow_dzmax_l_1 = nan
    snow_dzmax_l_2 = nan
    snow_dzmax_u_1 = nan
    snow_dzmax_u_2 = nan

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in clm_SnowHydrology_inparm  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, 'clm_SnowHydrology_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, clm_snowhydrology_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading clm_snowhydrology_inparm namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR finding clm_snowhydrology_inparm namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (wind_dependent_snow_density, mpicom)
    call shr_mpi_bcast (snow_overburden_compaction_method, mpicom)
    call shr_mpi_bcast (lotmp_snowdensity_method   , mpicom)
    call shr_mpi_bcast (upplim_destruct_metamorph  , mpicom)
    call shr_mpi_bcast (overburden_compress_Tfactor, mpicom)
    call shr_mpi_bcast (reset_snow                 , mpicom)
    call shr_mpi_bcast (reset_snow_glc             , mpicom)
    call shr_mpi_bcast (reset_snow_glc_ela         , mpicom)
    call shr_mpi_bcast (snow_dzmin_1, mpicom)
    call shr_mpi_bcast (snow_dzmin_2, mpicom)
    call shr_mpi_bcast (snow_dzmax_l_1, mpicom)
    call shr_mpi_bcast (snow_dzmax_l_2, mpicom)
    call shr_mpi_bcast (snow_dzmax_u_1, mpicom)
    call shr_mpi_bcast (snow_dzmax_u_2, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'SnowHydrology settings:'
       write(iulog,nml=clm_snowhydrology_inparm)
       write(iulog,*) ' '
    end if

    if (      trim(lotmp_snowdensity_method) == 'Slater2017' ) then
       new_snow_density = LoTmpDnsSlater2017
    else if ( trim(lotmp_snowdensity_method) == 'TruncatedAnderson1976' ) then
       new_snow_density = LoTmpDnsTruncatedAnderson1976
    else
       call endrun(msg="ERROR bad lotmp_snowdensity_method name"//errmsg(sourcefile, __LINE__))
    end if

    if (trim(snow_overburden_compaction_method) == 'Anderson1976') then
       overburden_compaction_method = OverburdenCompactionMethodAnderson1976
    else if (trim(snow_overburden_compaction_method) == 'Vionnet2012') then
       overburden_compaction_method = OverburdenCompactionMethodVionnet2012
    else
       call endrun(msg="ERROR bad snow_overburden_compaction_method name"// &
            errMsg(sourcefile, __LINE__))
    end if

  end subroutine SnowHydrology_readnl

  !----------------------------------------------------------------------------
  subroutine readParams( ncid )
    !
    ! !USES:
    use ncdio_pio, only: file_desc_t
    use paramUtilMod, only: readNcdioScalar
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'readParams_SnowHydrology'
    !--------------------------------------------------------------------

    ! Water impremeable if porosity less than wimp (unitless)
    call readNcdioScalar(ncid, 'wimp', subname, params_inst%wimp)
    ! Irreducible water saturation of snow (unitless)
    call readNcdioScalar(ncid, 'ssi', subname, params_inst%ssi)
    ! Wind drift compaction / grain size (fixed value for now) (unitless)
    call readNcdioScalar(ncid, 'drift_gs', subname, params_inst%drift_gs)
    ! Viscosity coefficent from Anderson1976 (kg*s/m2)
    call readNcdioScalar(ncid, 'eta0_anderson', subname, params_inst%eta0_anderson)
    ! Viscosity coefficent from Vionnet2012 (kg*s/m2)
    call readNcdioScalar(ncid, 'eta0_vionnet', subname, params_inst%eta0_vionnet)
    ! Reference wind above which fresh snow density is (substantially) increased (m/s)
    call readNcdioScalar(ncid, 'wind_snowcompact_fact', subname, params_inst%wind_snowcompact_fact)
    ! Wind drift compaction / maximum density (kg/m3)
    call readNcdioScalar(ncid, 'rho_max', subname, params_inst%rho_max)
    ! Wind drift compaction / reference time (48*3600) (s)
    call readNcdioScalar(ncid, 'tau_ref', subname, params_inst%tau_ref)
    ! Scaling factor modifying scavenging factors for BC, OC, and dust species inclusion in meltwater (-)
    call readNcdioScalar(ncid, 'scvng_fct_mlt_sf', subname, params_inst%scvng_fct_mlt_sf)
    ! Overburden compaction constant (kg/m3)
    call readNcdioScalar(ncid, 'ceta', subname, params_inst%ceta)

  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine UpdateQuantitiesForNewSnow(bounds, num_c, filter_c, &
       scf_method, atm2lnd_inst, water_inst)
    !
    ! !DESCRIPTION:
    ! Update various snow-related quantities to account for new snow
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_c          ! number of column points in column filter
    integer                , intent(in)    :: filter_c(:)    ! column filter
    class(snow_cover_fraction_base_type), intent(in) :: scf_method
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(water_type)       , intent(inout) :: water_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: i     ! index of water tracer or bulk
    real(r8) :: dtime ! land model time step (sec)
    real(r8) :: bifall(bounds%begc:bounds%endc)                   ! bulk density of newly fallen dry snow [kg/m3]
    real(r8) :: h2osno_total(bounds%begc:bounds%endc)             ! total snow water (mm H2O)

    character(len=*), parameter :: subname = 'UpdateQuantitiesForNewSnow'
    !-----------------------------------------------------------------------

    associate( &
         begc => bounds%begc, &
         endc => bounds%endc, &

         b_waterflux_inst       => water_inst%waterfluxbulk_inst, &
         b_waterstate_inst      => water_inst%waterstatebulk_inst, &
         b_waterdiagnostic_inst => water_inst%waterdiagnosticbulk_inst &
         )

    ! Get time step
    dtime = get_step_size_real()

    call NewSnowBulkDensity(bounds, num_c, filter_c, &
         atm2lnd_inst, bifall(bounds%begc:bounds%endc))

    call b_waterstate_inst%CalculateTotalH2osno(bounds, num_c, filter_c, &
         caller = 'HandleNewSnow', &
         h2osno_total = h2osno_total(bounds%begc:bounds%endc))

    call BulkDiag_NewSnowDiagnostics(bounds, num_c, filter_c, &
         ! Inputs
         scf_method          = scf_method, &
         dtime               = dtime, &
         lun_itype_col       = col%lun_itype(begc:endc), &
         urbpoi              = col%urbpoi(begc:endc), &
         snl                 = col%snl(begc:endc), &
         bifall              = bifall(begc:endc), &
         h2osno_total        = h2osno_total(begc:endc), &
         h2osoi_ice          = b_waterstate_inst%h2osoi_ice_col(begc:endc,:), &
         h2osoi_liq          = b_waterstate_inst%h2osoi_liq_col(begc:endc,:), &
         qflx_snow_grnd      = b_waterflux_inst%qflx_snow_grnd_col(begc:endc), &
         qflx_snow_drain     = b_waterflux_inst%qflx_snow_drain_col(begc:endc), &
         ! Outputs
         dz                  = col%dz(begc:endc,:), &
         int_snow            = b_waterstate_inst%int_snow_col(begc:endc), &
         swe_old             = b_waterdiagnostic_inst%swe_old_col(begc:endc,:), &
         frac_sno            = b_waterdiagnostic_inst%frac_sno_col(begc:endc), &
         frac_sno_eff        = b_waterdiagnostic_inst%frac_sno_eff_col(begc:endc), &
         snow_depth          = b_waterdiagnostic_inst%snow_depth_col(begc:endc))

    do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
       associate(w => water_inst%bulk_and_tracers(i))
       call UpdateState_AddNewSnow(bounds, num_c, filter_c, &
            ! Inputs
            dtime            = dtime, &
            snl              = col%snl(begc:endc), &
            qflx_snow_grnd   = w%waterflux_inst%qflx_snow_grnd_col(begc:endc), &
            ! Outputs
            h2osno_no_layers = w%waterstate_inst%h2osno_no_layers_col(begc:endc), &
            h2osoi_ice       = w%waterstate_inst%h2osoi_ice_col(begc:endc,:))
       end associate
    end do

    end associate

  end subroutine UpdateQuantitiesForNewSnow

  !-----------------------------------------------------------------------
  subroutine BulkDiag_NewSnowDiagnostics(bounds, num_c, filter_c, &
       scf_method, &
       dtime, lun_itype_col, urbpoi, snl, bifall, h2osno_total, h2osoi_ice, h2osoi_liq, &
       qflx_snow_grnd, qflx_snow_drain, &
       dz, int_snow, swe_old, frac_sno, frac_sno_eff, snow_depth)
    !
    ! !DESCRIPTION:
    ! Update various snow-related diagnostic quantities to account for new snow
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_c
    integer, intent(in) :: filter_c(:)

    class(snow_cover_fraction_base_type), intent(in) :: scf_method
    real(r8)                  , intent(in)    :: dtime                           ! land model time step (sec)
    integer                   , intent(in)    :: lun_itype_col( bounds%begc: )   ! landunit type for each column
    logical                   , intent(in)    :: urbpoi( bounds%begc: )          ! true=>urban point
    integer                   , intent(in)    :: snl( bounds%begc: )             ! negative number of snow layers
    real(r8)                  , intent(in)    :: bifall( bounds%begc: )          ! bulk density of newly fallen dry snow [kg/m3]
    real(r8)                  , intent(in)    :: h2osno_total( bounds%begc: )    ! total snow water (mm H2O)
    real(r8)                  , intent(in)    :: h2osoi_ice( bounds%begc: , -nlevsno+1: ) ! ice lens (kg/m2)
    real(r8)                  , intent(in)    :: h2osoi_liq( bounds%begc: , -nlevsno+1: ) ! liquid water (kg/m2)
    real(r8)                  , intent(in)    :: qflx_snow_grnd( bounds%begc: )  ! snow on ground after interception (mm H2O/s)
    real(r8)                  , intent(in)    :: qflx_snow_drain( bounds%begc: ) ! drainage from snow pack from previous time step (mm H2O/s)
    real(r8)                  , intent(inout) :: dz( bounds%begc: , -nlevsno+1: ) ! layer depth (m)
    real(r8)                  , intent(inout) :: int_snow( bounds%begc: )        ! integrated snowfall (mm H2O)
    real(r8)                  , intent(inout) :: swe_old( bounds%begc:, -nlevsno+1: )         ! snow water before update (mm H2O)
    real(r8)                  , intent(inout) :: frac_sno( bounds%begc: )        ! fraction of ground covered by snow (0 to 1)
    real(r8)                  , intent(inout) :: frac_sno_eff( bounds%begc: )    ! eff. fraction of ground covered by snow (0 to 1)
    real(r8)                  , intent(inout) :: snow_depth( bounds%begc: )      ! snow height (m)
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c
    integer  :: j
    real(r8) :: dz_snowf                              ! layer thickness rate change due to precipitation [mm/s]
    real(r8) :: temp_snow_depth(bounds%begc:bounds%endc) ! snow depth prior to updating [mm]
    real(r8) :: snowmelt(bounds%begc:bounds%endc)
    real(r8) :: newsnow(bounds%begc:bounds%endc)
    real(r8) :: z_avg                                 ! grid cell average snow depth

    character(len=*), parameter :: subname = 'BulkDiag_NewSnowDiagnostics'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(lun_itype_col, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(urbpoi, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(snl, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(bifall, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osno_total, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osoi_ice, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osoi_liq, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_snow_grnd, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_snow_drain, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(dz, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(int_snow, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(swe_old) == [bounds%endc, 0]), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(frac_sno, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(frac_sno_eff, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(snow_depth, 1) == bounds%endc), sourcefile, __LINE__)

    associate( &
         begc => bounds%begc, &
         endc => bounds%endc  &
         )

    do fc = 1, num_c
       c = filter_c(fc)
       
       ! Use Alta relationship, Anderson(1976); LaChapelle(1961),
       ! U.S.Department of Agriculture Forest Service, Project F,
       ! Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.

       ! set temporary variables prior to updating
       temp_snow_depth(c) = snow_depth(c)
       ! save initial snow content
       do j= -nlevsno+1,snl(c)
          swe_old(c,j) = 0.0_r8
       end do
       do j= snl(c)+1,0
          swe_old(c,j)=h2osoi_liq(c,j)+h2osoi_ice(c,j)
       end do

       ! all snow falls on ground, no snow on h2osfc (note that qflx_snow_h2osfc is
       ! currently set to 0 always in CanopyHydrologyMod)
       newsnow(c) = qflx_snow_grnd(c) * dtime

       ! update int_snow
       int_snow(c) = max(int_snow(c),h2osno_total(c)) !h2osno_total could be larger due to frost

       ! snowmelt from previous time step * dtime
       snowmelt(c) = qflx_snow_drain(c) * dtime
    end do

    call scf_method%UpdateSnowDepthAndFrac(bounds, num_c, filter_c, &
         ! Inputs
         lun_itype_col = lun_itype_col(begc:endc), &
         urbpoi        = urbpoi(begc:endc), &
         h2osno_total  = h2osno_total(begc:endc), &
         snowmelt      = snowmelt(begc:endc), &
         int_snow      = int_snow(begc:endc), &
         newsnow       = newsnow(begc:endc), &
         bifall        = bifall(begc:endc), &
         ! Outputs
         snow_depth    = snow_depth(begc:endc), &
         frac_sno      = frac_sno(begc:endc), &
         frac_sno_eff  = frac_sno_eff(begc:endc))

    do fc = 1, num_c
       c = filter_c(fc)
       if (h2osno_total(c) == 0._r8 .and. newsnow(c) > 0._r8) then
          ! Snow pack started at 0, then adding new snow; reset int_snow prior to
          ! adding newsnow
          int_snow(c) = 0._r8
       end if
    end do

    call scf_method%AddNewsnowToIntsnow(bounds, num_c, filter_c, &
         ! Inputs
         newsnow      = newsnow(begc:endc), &
         h2osno_total = h2osno_total(begc:endc), &
         frac_sno     = frac_sno(begc:endc), &
         ! Outputs
         int_snow     = int_snow(begc:endc))

    do fc = 1, num_c
       c = filter_c(fc)

       ! update change in snow depth
       if (snl(c) < 0) then
          dz_snowf = (snow_depth(c) - temp_snow_depth(c)) / dtime
          dz(c,snl(c)+1) = dz(c,snl(c)+1)+dz_snowf*dtime
       end if

    end do

    end associate

  end subroutine BulkDiag_NewSnowDiagnostics

  !-----------------------------------------------------------------------
  subroutine UpdateState_AddNewSnow(bounds, num_c, filter_c, &
       dtime, snl, qflx_snow_grnd, &
       h2osno_no_layers, h2osoi_ice)
    !
    ! !DESCRIPTION:
    ! Update h2osno_no_layers or h2osoi_ice based on new snow
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_c
    integer, intent(in) :: filter_c(:)

    real(r8) , intent(in)    :: dtime                                    ! land model time step (sec)
    integer  , intent(in)    :: snl( bounds%begc: )                      ! negative number of snow layers
    real(r8) , intent(in)    :: qflx_snow_grnd( bounds%begc: )           ! snow on ground after interception (mm H2O/s)
    real(r8) , intent(inout) :: h2osno_no_layers( bounds%begc: )         ! snow that is not resolved into layers (kg/m2)
    real(r8) , intent(inout) :: h2osoi_ice( bounds%begc: , -nlevsno+1: ) ! ice lens (kg/m2)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c

    character(len=*), parameter :: subname = 'UpdateState_AddNewSnow'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(snl, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_snow_grnd, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osno_no_layers, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osoi_ice, 1) == bounds%endc), sourcefile, __LINE__)

    do fc = 1, num_c
       c = filter_c(fc)

       if (snl(c) == 0) then
          h2osno_no_layers(c) = h2osno_no_layers(c) + (qflx_snow_grnd(c) * dtime)
       else
          ! The change of ice partial density of surface node due to precipitation. Only
          ! ice part of snowfall is added here, the liquid part will be added later.
          h2osoi_ice(c,snl(c)+1) = h2osoi_ice(c,snl(c)+1) + (qflx_snow_grnd(c) * dtime)
       end if
    end do

  end subroutine UpdateState_AddNewSnow

  !-----------------------------------------------------------------------
  subroutine RemoveSnowFromThawedWetlands(bounds, num_nolakec, filter_nolakec, &
       temperature_inst, water_inst)
    !
    ! !DESCRIPTION:
    ! Remove snow from thawed wetlands
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_nolakec
    integer, intent(in) :: filter_nolakec(:)

    type(temperature_type) , intent(in)    :: temperature_inst
    type(water_type)       , intent(inout) :: water_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: i     ! index of water tracer or bulk
    type(filter_col_type) :: thawed_wetland_thin_snowpack_filterc ! column filter: thawed wetland columns with a thin (no-layer) snow pack

    character(len=*), parameter :: subname = 'RemoveSnowFromThawedWetlands'
    !-----------------------------------------------------------------------

    associate( &
         begc => bounds%begc, &
         endc => bounds%endc, &

         b_waterdiagnostic_inst => water_inst%waterdiagnosticbulk_inst &
         )

    ! BUG(wjs, 2019-06-05, ESCOMP/ctsm#735) It seems like the intended behavior here is to
    ! zero out the entire snow pack. Currently, however, this only zeros out a very thin
    ! (zero-layer) snow pack. For now, I'm adding an snl==0 conditional to make this
    ! behavior explicit; long-term, we'd like to change this to actually zero out the
    ! whole snow pack. (At that time, this snow removal should probably be done from a
    ! more appropriate point of the driver loop, as noted in comments in issue #735.)

    call BuildFilter_ThawedWetlandThinSnowpack(bounds, num_nolakec, filter_nolakec, &
         ! Inputs
         t_grnd        = temperature_inst%t_grnd_col(begc:endc), &
         lun_itype_col = col%lun_itype(begc:endc), &
         snl           = col%snl(begc:endc), &
         ! Outputs
         thawed_wetland_thin_snowpack_filterc = thawed_wetland_thin_snowpack_filterc)

    do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
       associate(w => water_inst%bulk_and_tracers(i))
       call UpdateState_RemoveSnowFromThawedWetlands(bounds, thawed_wetland_thin_snowpack_filterc, &
            ! Outputs
            h2osno_no_layers = w%waterstate_inst%h2osno_no_layers_col(begc:endc))
       end associate
    end do

    call Bulk_RemoveSnowFromThawedWetlands(bounds, thawed_wetland_thin_snowpack_filterc, &
         ! Outputs
         snow_depth = b_waterdiagnostic_inst%snow_depth_col(begc:endc))

    end associate

  end subroutine RemoveSnowFromThawedWetlands

  !-----------------------------------------------------------------------
  subroutine BuildFilter_ThawedWetlandThinSnowpack(bounds, num_nolakec, filter_nolakec, &
       t_grnd, lun_itype_col, snl, thawed_wetland_thin_snowpack_filterc)
    !
    ! !DESCRIPTION:
    ! Build a column-level filter of thawed wetland columns with a thin snowpack (which
    ! we will subsequently remove)
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_nolakec
    integer, intent(in) :: filter_nolakec(:)

    real(r8)              , intent(in)  :: t_grnd( bounds%begc: )               ! ground temperature (Kelvin)
    integer               , intent(in)  :: lun_itype_col( bounds%begc: )        ! landunit type for each column
    integer               , intent(in)  :: snl( bounds%begc: )                  ! negative number of snow layers
    type(filter_col_type) , intent(out) :: thawed_wetland_thin_snowpack_filterc ! column filter: thawed wetland columns with a thin (no-layer) snow pack
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c
    logical :: thawed_wetland_thin_snowpack(bounds%begc:bounds%endc)

    character(len=*), parameter :: subname = 'BuildFilter_ThawedWetlandThinSnowpack'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(t_grnd, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(lun_itype_col, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(snl, 1) == bounds%endc), sourcefile, __LINE__)

    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       thawed_wetland_thin_snowpack(c) = .false.
       if (lun_itype_col(c) == istwet .and. t_grnd(c) > tfrz .and. snl(c) == 0) then
          thawed_wetland_thin_snowpack(c) = .true.
       end if
    end do

    thawed_wetland_thin_snowpack_filterc = col_filter_from_filter_and_logical_array( &
         bounds = bounds, &
         num_orig = num_nolakec, &
         filter_orig = filter_nolakec, &
         logical_col = thawed_wetland_thin_snowpack(bounds%begc:bounds%endc))

  end subroutine BuildFilter_ThawedWetlandThinSnowpack

  !-----------------------------------------------------------------------
  subroutine UpdateState_RemoveSnowFromThawedWetlands(bounds, &
       thawed_wetland_thin_snowpack_filterc, &
       h2osno_no_layers)
    !
    ! !DESCRIPTION:
    ! For bulk or one tracer: remove snow from thawed wetlands
    !
    ! This routine applies to state variables that apply to both bulk and tracers
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    type(filter_col_type), intent(in) :: thawed_wetland_thin_snowpack_filterc ! column filter: thawed wetland columns with a thin (no-layer) snow pack

    real(r8) , intent(inout) :: h2osno_no_layers( bounds%begc: ) ! snow that is not resolved into layers (kg/m2)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c

    character(len=*), parameter :: subname = 'UpdateState_RemoveSnowFromThawedWetlands'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(h2osno_no_layers, 1) == bounds%endc), sourcefile, __LINE__)

    do fc = 1, thawed_wetland_thin_snowpack_filterc%num
       c = thawed_wetland_thin_snowpack_filterc%indices(fc)

       h2osno_no_layers(c) = 0._r8
    end do

  end subroutine UpdateState_RemoveSnowFromThawedWetlands

  !-----------------------------------------------------------------------
  subroutine Bulk_RemoveSnowFromThawedWetlands(bounds, &
       thawed_wetland_thin_snowpack_filterc, &
       snow_depth)
    !
    ! !DESCRIPTION:
    ! Remove snow from thawed wetlands
    !
    ! This routine just operates on bulk-only quantities; state variables are handled in
    ! UpdateState_RemoveSnowFromThawedWetlands.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    type(filter_col_type), intent(in) :: thawed_wetland_thin_snowpack_filterc ! column filter: thawed wetland columns with a thin (no-layer) snow pack

    real(r8)            , intent(inout) :: snow_depth( bounds%begc: ) ! snow height (m)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c

    character(len=*), parameter :: subname = 'Bulk_RemoveSnowFromThawedWetlands'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(snow_depth, 1) == bounds%endc), sourcefile, __LINE__)

    do fc = 1, thawed_wetland_thin_snowpack_filterc%num
       c = thawed_wetland_thin_snowpack_filterc%indices(fc)

       snow_depth(c) = 0._r8
    end do

  end subroutine Bulk_RemoveSnowFromThawedWetlands

  !-----------------------------------------------------------------------
  subroutine InitializeExplicitSnowPack(bounds, num_c, filter_c, &
       atm2lnd_inst, temperature_inst, aerosol_inst, water_inst)
    !
    ! !DESCRIPTION:
    ! Initialize an explicit snow pack in columns where this is warranted based on snow
    ! depth
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_c          ! number of column points in column filter
    integer                , intent(in)    :: filter_c(:)    ! column filter
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(temperature_type) , intent(inout) :: temperature_inst
    type(aerosol_type)     , intent(inout) :: aerosol_inst
    type(water_type)       , intent(inout) :: water_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: i     ! index of water tracer or bulk
    type(filter_col_type) :: snowpack_initialized_filterc         ! column filter: columns where an explicit snow pack is initialized

    character(len=*), parameter :: subname = 'InitializeExplicitSnowPack'
    !-----------------------------------------------------------------------

    associate( &
         begc => bounds%begc, &
         endc => bounds%endc, &

         b_waterflux_inst       => water_inst%waterfluxbulk_inst, &
         b_waterdiagnostic_inst => water_inst%waterdiagnosticbulk_inst  &
         )

    call BuildFilter_SnowpackInitialized(bounds, num_c, filter_c, &
         ! Inputs
         snl                  = col%snl(begc:endc), &
         lun_itype_col        = col%lun_itype(begc:endc), &
         frac_sno_eff         = b_waterdiagnostic_inst%frac_sno_eff_col(begc:endc), &
         snow_depth           = b_waterdiagnostic_inst%snow_depth_col(begc:endc), &
         qflx_snow_grnd       = b_waterflux_inst%qflx_snow_grnd_col(begc:endc), &
         ! Outputs
         snowpack_initialized_filterc = snowpack_initialized_filterc)

    do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
       associate(w => water_inst%bulk_and_tracers(i))
       call UpdateState_InitializeSnowPack(bounds, snowpack_initialized_filterc, &
            ! Outputs
            h2osno_no_layers     = w%waterstate_inst%h2osno_no_layers_col(begc:endc), &
            h2osoi_ice           = w%waterstate_inst%h2osoi_ice_col(begc:endc,:), &
            h2osoi_liq           = w%waterstate_inst%h2osoi_liq_col(begc:endc,:))
       end associate
    end do

    call Bulk_InitializeSnowPack(bounds, snowpack_initialized_filterc, &
         ! Inputs
         forc_t      = atm2lnd_inst%forc_t_downscaled_col(begc:endc), &
         snow_depth  = b_waterdiagnostic_inst%snow_depth_col(begc:endc), &
         ! Outputs
         snl         = col%snl(begc:endc), &
         zi          = col%zi(begc:endc,:), &
         dz          = col%dz(begc:endc,:), &
         z           = col%z(begc:endc,:), &
         t_soisno    = temperature_inst%t_soisno_col(begc:endc,:), &
         frac_iceold = b_waterdiagnostic_inst%frac_iceold_col(begc:endc,:))

    ! intitialize SNICAR variables for fresh snow:
    call aerosol_inst%ResetFilter( &
         num_c    = snowpack_initialized_filterc%num, &
         filter_c = snowpack_initialized_filterc%indices)
    call b_waterdiagnostic_inst%ResetBulkFilter( &
         num_c    = snowpack_initialized_filterc%num, &
         filter_c = snowpack_initialized_filterc%indices)

    end associate

  end subroutine InitializeExplicitSnowPack

  !-----------------------------------------------------------------------
  subroutine BuildFilter_SnowpackInitialized(bounds, num_c, filter_c, &
       snl, lun_itype_col, frac_sno_eff, snow_depth, qflx_snow_grnd, &
       snowpack_initialized_filterc)
    !
    ! !DESCRIPTION:
    ! Build a column-level filter of columns where an explicit snow pack needs to be initialized
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_c
    integer, intent(in) :: filter_c(:)

    integer               , intent(in)  :: snl( bounds%begc: )            ! negative number of snow layers
    integer               , intent(in)  :: lun_itype_col( bounds%begc: )  ! landunit type for each column
    real(r8)              , intent(in)  :: frac_sno_eff( bounds%begc: )   ! fraction of ground covered by snow (0 to 1)
    real(r8)              , intent(in)  :: snow_depth( bounds%begc: )     ! snow height (m)
    real(r8)              , intent(in)  :: qflx_snow_grnd( bounds%begc: ) ! snow on ground after interception (mm H2O/s)
    type(filter_col_type) , intent(out) :: snowpack_initialized_filterc   ! column filter: columns where an explicit snow pack is initialized
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c
    logical :: snowpack_initialized(bounds%begc:bounds%endc)

    character(len=*), parameter :: subname = 'BuildFilter_SnowpackInitialized'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(snl, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(lun_itype_col, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(frac_sno_eff, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(snow_depth, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_snow_grnd, 1) == bounds%endc), sourcefile, __LINE__)

    do fc = 1, num_c
       c = filter_c(fc)

       if (lun_itype_col(c) == istdlak) then
          ! NOTE(wjs, 2019-08-23) Note two differences from the standard case: (1)
          ! addition of lsadz (see notes in LakeCon.F90, where this is defined, for
          ! details); (2) inclusion of a qflx_snow_grnd(c) > 0 criteria. I'm not sure why
          ! (2) is needed here and not for other landunits, but I'm keeping it here to
          ! maintain answers as before.
          snowpack_initialized(c) = ( &
               snl(c) == 0 .and. &
               frac_sno_eff(c)*snow_depth(c) >= (dzmin(1) + lsadz) .and. &
               qflx_snow_grnd(c) > 0.0_r8)
       else
          snowpack_initialized(c) = ( &
               snl(c) == 0 .and. &
               frac_sno_eff(c)*snow_depth(c) >= dzmin(1))
       end if

    end do

    snowpack_initialized_filterc = col_filter_from_filter_and_logical_array( &
         bounds = bounds, &
         num_orig = num_c, &
         filter_orig = filter_c, &
         logical_col = snowpack_initialized(bounds%begc:bounds%endc))

  end subroutine BuildFilter_SnowpackInitialized

  !-----------------------------------------------------------------------
  subroutine UpdateState_InitializeSnowPack(bounds, snowpack_initialized_filterc, &
       h2osno_no_layers, h2osoi_ice, h2osoi_liq)
    !
    ! !DESCRIPTION:
    ! For bulk or one tracer: initialize water state variables for columns in which an
    ! explicit snow pack is being newly initialized.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    type(filter_col_type), intent(in) :: snowpack_initialized_filterc ! column filter: columns where an explicit snow pack is initialized

    real(r8) , intent(inout) :: h2osno_no_layers( bounds%begc: )         ! snow that is not resolved into layers (kg/m2)
    real(r8) , intent(inout) :: h2osoi_ice( bounds%begc: , -nlevsno+1: ) ! ice lens (kg/m2)
    real(r8) , intent(inout) :: h2osoi_liq( bounds%begc: , -nlevsno+1: ) ! liquid water (kg/m2)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c

    character(len=*), parameter :: subname = 'UpdateState_InitializeSnowPack'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(h2osno_no_layers, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osoi_ice, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osoi_liq, 1) == bounds%endc), sourcefile, __LINE__)

    do fc = 1, snowpack_initialized_filterc%num
       c = snowpack_initialized_filterc%indices(fc)

       h2osoi_ice(c,0) = h2osno_no_layers(c)
       h2osoi_liq(c,0) = 0._r8
       h2osno_no_layers(c) = 0._r8
    end do

  end subroutine UpdateState_InitializeSnowPack

  !-----------------------------------------------------------------------
  subroutine Bulk_InitializeSnowPack(bounds, snowpack_initialized_filterc, &
       forc_t, snow_depth, snl, zi, dz, z, t_soisno, frac_iceold)
    !
    ! !DESCRIPTION:
    ! Initialize an explicit snow pack in columns where this is warranted based on snow depth
    !
    ! This routine just operates on bulk-only quantities; state variables are handled in
    ! UpdateState_InitializeSnowPack.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    type(filter_col_type), intent(in) :: snowpack_initialized_filterc ! column filter: columns where an explicit snow pack is initialized

    real(r8)                       , intent(in)    :: forc_t( bounds%begc: )                    ! atmospheric temperature (Kelvin)
    real(r8)                       , intent(in)    :: snow_depth( bounds%begc: )                ! snow height (m)
    integer                        , intent(inout) :: snl( bounds%begc: )                       ! negative number of snow layers
    real(r8)                       , intent(inout) :: zi( bounds%begc: , -nlevsno+0: )          ! interface level below a "z" level (m)
    real(r8)                       , intent(inout) :: dz( bounds%begc: , -nlevsno+1: )          ! layer thickness (m)
    real(r8)                       , intent(inout) :: z( bounds%begc: , -nlevsno+1: )           ! layer depth (m)
    real(r8)                       , intent(inout) :: t_soisno( bounds%begc: , -nlevsno+1: )    ! soil temperature (Kelvin)
    real(r8)                       , intent(inout) :: frac_iceold( bounds%begc: , -nlevsno+1: ) ! fraction of ice relative to the tot water
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c

    character(len=*), parameter :: subname = 'Bulk_InitializeSnowPack'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(forc_t, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(snow_depth, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(snl, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(zi) == [bounds%endc, nlevmaxurbgrnd]), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(dz) == [bounds%endc, nlevmaxurbgrnd]), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(z) == [bounds%endc, nlevmaxurbgrnd]), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(t_soisno) == [bounds%endc, nlevmaxurbgrnd]), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(frac_iceold) == [bounds%endc, nlevgrnd]), sourcefile, __LINE__)

    do fc = 1, snowpack_initialized_filterc%num
       c = snowpack_initialized_filterc%indices(fc)

       snl(c) = -1
       dz(c,0) = snow_depth(c)
       z(c,0) = -0.5_r8*dz(c,0)
       zi(c,-1) = -dz(c,0)
       ! Currently, the water temperature for the precipitation is simply set
       ! as the surface air temperature
       t_soisno(c,0) = min(tfrz, forc_t(c))

       ! This value of frac_iceold makes sense together with the state initialization:
       ! h2osoi_ice is non-zero, while h2osoi_liq is zero.
       frac_iceold(c,0) = 1._r8
    end do

  end subroutine Bulk_InitializeSnowPack

  !-----------------------------------------------------------------------
  subroutine SnowWater(bounds, &
       num_snowc, filter_snowc, num_nosnowc, filter_nosnowc, &
       atm2lnd_inst, aerosol_inst, water_inst)
    !
    ! !DESCRIPTION:
    ! Evaluate the change of snow mass and the snow water onto soil.
    ! Water flow within snow is computed by an explicit and non-physical
    ! based scheme, which permits a part of liquid water over the holding
    ! capacity (a tentative value is used, i.e. equal to 0.033*porosity) to
    ! percolate into the underlying layer.  Except for cases where the
    ! porosity of one of the two neighboring layers is less than 0.05, zero
    ! flow is assumed. The water flow out of the bottom of the snow pack will
    ! participate as the input of the soil water and runoff.  This subroutine
    ! uses a filter for columns containing snow which must be constructed prior
    ! to being called.
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)    :: bounds
    integer               , intent(in)    :: num_snowc         ! number of snow points in column filter
    integer               , intent(in)    :: filter_snowc(:)   ! column filter for snow points
    integer               , intent(in)    :: num_nosnowc       ! number of non-snow points in column filter
    integer               , intent(in)    :: filter_nosnowc(:) ! column filter for non-snow points
    type(atm2lnd_type)    , intent(in)    :: atm2lnd_inst
    type(aerosol_type)    , intent(inout) :: aerosol_inst
    type(water_type)      , intent(inout) :: water_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: i                                                  ! index of water tracer or bulk
    integer  :: g                                                  ! gridcell loop index
    integer  :: c, j, fc, l                                        ! do loop/array indices
    real(r8) :: dtime                                              ! land model time step (sec)
    !-----------------------------------------------------------------------

    associate( &
         begc => bounds%begc, &
         endc => bounds%endc, &

         b_waterflux_inst       => water_inst%waterfluxbulk_inst, &
         b_waterstate_inst      => water_inst%waterstatebulk_inst, &
         b_waterdiagnostic_inst => water_inst%waterdiagnosticbulk_inst &
         )

    ! Determine model time step

    dtime = get_step_size_real()

    do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
       associate(w => water_inst%bulk_and_tracers(i))
       call UpdateState_TopLayerFluxes(bounds, num_snowc, filter_snowc, &
            ! Inputs
            name           = water_inst%GetBulkOrTracerName(i), &
            dtime          = dtime, &
            snl            = col%snl(begc:endc), &
            frac_sno_eff   = b_waterdiagnostic_inst%frac_sno_eff_col(begc:endc), &
            qflx_soliddew_to_top_layer    = w%waterflux_inst%qflx_soliddew_to_top_layer_col(begc:endc), &
            qflx_solidevap_from_top_layer = w%waterflux_inst%qflx_solidevap_from_top_layer_col(begc:endc), &
            qflx_liq_grnd                 = w%waterflux_inst%qflx_liq_grnd_col(begc:endc), &
            qflx_liqdew_to_top_layer      = w%waterflux_inst%qflx_liqdew_to_top_layer_col(begc:endc), &
            qflx_liqevap_from_top_layer   = w%waterflux_inst%qflx_liqevap_from_top_layer_col(begc:endc), &
            ! Outputs
            h2osoi_ice     = w%waterstate_inst%h2osoi_ice_col(begc:endc,:), &
            h2osoi_liq     = w%waterstate_inst%h2osoi_liq_col(begc:endc,:))
       end associate
    end do

    call BulkFlux_SnowPercolation(bounds, num_snowc, filter_snowc, &
         ! Inputs
         dtime                 = dtime, &
         snl                   = col%snl(begc:endc), &
         dz                    = col%dz(begc:endc,:), &
         frac_sno_eff          = b_waterdiagnostic_inst%frac_sno_eff_col(begc:endc), &
         h2osoi_ice            = b_waterstate_inst%h2osoi_ice_col(begc:endc,:), &
         h2osoi_liq            = b_waterstate_inst%h2osoi_liq_col(begc:endc,:), &
         ! Outputs
         qflx_snow_percolation = b_waterflux_inst%qflx_snow_percolation_col(begc:endc,:))

    do i = water_inst%tracers_beg, water_inst%tracers_end
       associate(w => water_inst%bulk_and_tracers(i))
       call TracerFlux_SnowPercolation(bounds, num_snowc, filter_snowc, &
            ! Inputs
            snl                        = col%snl(begc:endc), &
            bulk_h2osoi_liq            = b_waterstate_inst%h2osoi_liq_col(begc:endc,:), &
            bulk_qflx_snow_percolation = b_waterflux_inst%qflx_snow_percolation_col(begc:endc,:), &
            trac_h2osoi_liq            = w%waterstate_inst%h2osoi_liq_col(begc:endc,:), &
            ! Outputs
            trac_qflx_snow_percolation = w%waterflux_inst%qflx_snow_percolation_col(begc:endc,:))
       end associate
    end do

    do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
       associate(w => water_inst%bulk_and_tracers(i))
       call UpdateState_SnowPercolation(bounds, num_snowc, filter_snowc, &
            ! Inputs
            dtime                 = dtime, &
            snl                   = col%snl(begc:endc), &
            qflx_snow_percolation = w%waterflux_inst%qflx_snow_percolation_col(begc:endc,:), &
            ! Outputs
            h2osoi_liq            = w%waterstate_inst%h2osoi_liq_col(begc:endc,:))
       end associate
    end do

    call CalcAndApplyAerosolFluxes(bounds, num_snowc, filter_snowc, &
         ! Inputs
         dtime                 = dtime, &
         snl                   = col%snl(begc:endc), &
         h2osoi_ice            = b_waterstate_inst%h2osoi_ice_col(begc:endc,:), &
         h2osoi_liq            = b_waterstate_inst%h2osoi_liq_col(begc:endc,:), &
         qflx_snow_percolation = b_waterflux_inst%qflx_snow_percolation_col(begc:endc,:), &
         atm2lnd_inst          = atm2lnd_inst, &
         ! Outputs
         aerosol_inst          = aerosol_inst)

    call PostPercolation_AdjustLayerThicknesses(bounds, num_snowc, filter_snowc, &
         ! Inputs
         snl        = col%snl(begc:endc), &
         h2osoi_ice = b_waterstate_inst%h2osoi_ice_col(begc:endc,:), &
         h2osoi_liq = b_waterstate_inst%h2osoi_liq_col(begc:endc,:), &
         ! Outputs
         dz         = col%dz(begc:endc,:))

    call BulkDiag_SnowWaterAccumulatedSnow(bounds, &
         num_snowc, filter_snowc, num_nosnowc, filter_nosnowc, &
         ! Inputs
         dtime            = dtime, &
         frac_sno_eff     = b_waterdiagnostic_inst%frac_sno_eff_col(begc:endc), &
         qflx_soliddew_to_top_layer = b_waterflux_inst%qflx_soliddew_to_top_layer_col(begc:endc), &
         qflx_liqdew_to_top_layer   = b_waterflux_inst%qflx_liqdew_to_top_layer_col(begc:endc), &
         qflx_liq_grnd              = b_waterflux_inst%qflx_liq_grnd_col(begc:endc), &
         h2osno_no_layers           = b_waterstate_inst%h2osno_no_layers_col(begc:endc), &
         ! Outputs
         int_snow         = b_waterstate_inst%int_snow_col(begc:endc), &
         frac_sno         = b_waterdiagnostic_inst%frac_sno_col(begc:endc), &
         snow_depth       = b_waterdiagnostic_inst%snow_depth_col(begc:endc))
         
    do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
       associate(w => water_inst%bulk_and_tracers(i))
       call SumFlux_AddSnowPercolation(bounds, &
            num_snowc, filter_snowc, num_nosnowc, filter_nosnowc, &
            ! Inputs
            frac_sno_eff                 = b_waterdiagnostic_inst%frac_sno_eff_col(begc:endc), &
            qflx_snow_percolation_bottom = w%waterflux_inst%qflx_snow_percolation_col(begc:endc, 0), &
            qflx_liq_grnd                = w%waterflux_inst%qflx_liq_grnd_col(begc:endc), &
            qflx_snomelt                 = w%waterflux_inst%qflx_snomelt_col(begc:endc), &
            ! Outputs
            qflx_snow_drain              = w%waterflux_inst%qflx_snow_drain_col(begc:endc), &
            qflx_rain_plus_snomelt       = w%waterflux_inst%qflx_rain_plus_snomelt_col(begc:endc))
       end associate
    end do

    end associate
  end subroutine SnowWater

  !-----------------------------------------------------------------------
  subroutine UpdateState_TopLayerFluxes(bounds, num_snowc, filter_snowc, &
       name, dtime, snl, frac_sno_eff, &
       qflx_soliddew_to_top_layer, qflx_solidevap_from_top_layer, qflx_liq_grnd, &
       qflx_liqdew_to_top_layer, qflx_liqevap_from_top_layer, h2osoi_ice, h2osoi_liq)
    !
    ! !DESCRIPTION:
    ! Update top layer of snow pack with various fluxes into and out of the top layer,
    ! for bulk or one tracer.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_snowc
    integer, intent(in) :: filter_snowc(:)

    character(len=*) , intent(in) :: name                           ! Name of bulk or this tracer (for output in case there's an error)
    real(r8)         , intent(in) :: dtime                          ! land model time step (sec)
    integer          , intent(in) :: snl( bounds%begc: )            ! negative number of snow layers
    real(r8)         , intent(in) :: frac_sno_eff( bounds%begc: )   ! eff. fraction of ground covered by snow (0 to 1)
    real(r8)         , intent(in) :: qflx_liqdew_to_top_layer( bounds%begc: ) ! rate of liquid water deposited on top soil or snow layer (dew) (mm H2O /s)
    real(r8)         , intent(in) :: qflx_solidevap_from_top_layer( bounds%begc: ) ! rate of ice evaporated from top soil or snow layer (sublimation) (mm H2O /s)
    real(r8)         , intent(in) :: qflx_liq_grnd( bounds%begc: )  ! liquid on ground after interception (mm H2O/s)
    real(r8)         , intent(in) :: qflx_soliddew_to_top_layer( bounds%begc: ) ! rate of solid water deposited on top soil or snow layer (frost) (mm H2O /s)
    real(r8)         , intent(in) :: qflx_liqevap_from_top_layer( bounds%begc: ) ! rate of liquid water evaporated from top soil or snow layer (mm H2O/s)

    real(r8) , intent(inout) :: h2osoi_ice( bounds%begc: , -nlevsno+1: ) ! ice lens (kg/m2)
    real(r8) , intent(inout) :: h2osoi_liq( bounds%begc: , -nlevsno+1: ) ! liquid water (kg/m2)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c
    integer  :: lev_top(bounds%begc:bounds%endc)                   ! index of the top snow level
    real(r8) :: h2osoi_ice_top_orig(bounds%begc:bounds%endc)       ! h2osoi_ice in top snow layer before state update
    real(r8) :: h2osoi_liq_top_orig(bounds%begc:bounds%endc)       ! h2osoi_liq in top snow layer before state update

    character(len=*), parameter :: subname = 'UpdateState_TopLayerFluxes'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(snl, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(frac_sno_eff, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_soliddew_to_top_layer, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_solidevap_from_top_layer, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_liq_grnd, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_liqdew_to_top_layer, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_liqevap_from_top_layer, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osoi_ice, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osoi_liq, 1) == bounds%endc), sourcefile, __LINE__)

    ! Renew the mass of ice lens (h2osoi_ice) and liquid (h2osoi_liq) in the
    ! surface snow layer resulting from sublimation (frost) / evaporation (condense)

    do fc = 1,num_snowc
       c = filter_snowc(fc)

       lev_top(c) = snl(c)+1
       h2osoi_ice_top_orig(c) = h2osoi_ice(c,lev_top(c))
       h2osoi_liq_top_orig(c) = h2osoi_liq(c,lev_top(c))

       h2osoi_ice(c,lev_top(c)) = h2osoi_ice(c,lev_top(c)) &
            + frac_sno_eff(c) * (qflx_soliddew_to_top_layer(c) &
            - qflx_solidevap_from_top_layer(c)) * dtime

       h2osoi_liq(c,lev_top(c)) = h2osoi_liq(c,lev_top(c)) +  &
            frac_sno_eff(c) * (qflx_liq_grnd(c) + qflx_liqdew_to_top_layer(c) &
            - qflx_liqevap_from_top_layer(c)) * dtime
    end do

    ! If states were supposed to go to 0 but instead ended up near-0 (positive or
    ! negative), truncate to exactly 0.
    call truncate_small_values_one_lev( &
         num_f = num_snowc, &
         filter_f = filter_snowc, &
         lb = bounds%begc, &
         ub = bounds%endc, &
         lev_lb = -nlevsno+1, &
         lev = lev_top(bounds%begc:bounds%endc), &
         data_baseline = h2osoi_ice_top_orig(bounds%begc:bounds%endc), &
         data = h2osoi_ice(bounds%begc:bounds%endc, :), &
         custom_rel_epsilon = 1.e-12_r8)
    call truncate_small_values_one_lev( &
         num_f = num_snowc, &
         filter_f = filter_snowc, &
         lb = bounds%begc, &
         ub = bounds%endc, &
         lev_lb = -nlevsno+1, &
         lev = lev_top(bounds%begc:bounds%endc), &
         data_baseline = h2osoi_liq_top_orig(bounds%begc:bounds%endc), &
         data = h2osoi_liq(bounds%begc:bounds%endc, :), &
         custom_rel_epsilon = 1.e-12_r8)

    ! Make sure that we don't have any negative residuals - i.e., that we didn't try to
    ! remove more ice or liquid than was initially present.
    do fc = 1, num_snowc
       c = filter_snowc(fc)

       if (h2osoi_ice(c,lev_top(c)) < 0._r8) then
          write(iulog,*) "ERROR: In UpdateState_TopLayerFluxes, h2osoi_ice has gone significantly negative"
          write(iulog,*) "Bulk/tracer name = ", name
          write(iulog,*) "c, lev_top(c) = ", c, lev_top(c)
          write(iulog,*) "h2osoi_ice_top_orig = ", h2osoi_ice_top_orig(c)
          write(iulog,*) "h2osoi_ice          = ", h2osoi_ice(c,lev_top(c))
          write(iulog,*) "frac_sno_eff        = ", frac_sno_eff(c)
          write(iulog,*) "qflx_soliddew_to_top_layer*dtime = ", qflx_soliddew_to_top_layer(c)*dtime
          write(iulog,*) "qflx_solidevap_from_top_layer*dtime = ", qflx_solidevap_from_top_layer(c)*dtime
          call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, &
               msg="In UpdateState_TopLayerFluxes, h2osoi_ice has gone significantly negative")
       end if

       if (h2osoi_liq(c,lev_top(c)) < 0._r8) then
          write(iulog,*) "ERROR: In UpdateState_TopLayerFluxes, h2osoi_liq has gone significantly negative"
          write(iulog,*) "Bulk/tracer name = ", name
          write(iulog,*) "c, lev_top(c) = ", c, lev_top(c)
          write(iulog,*) "h2osoi_liq_top_orig  = ", h2osoi_liq_top_orig(c)
          write(iulog,*) "h2osoi_liq           = ", h2osoi_liq(c,lev_top(c))
          write(iulog,*) "frac_sno_eff         = ", frac_sno_eff(c)
          write(iulog,*) "qflx_liq_grnd*dtime  = ", qflx_liq_grnd(c)*dtime
          write(iulog,*) "qflx_liqdew_to_top_layer*dtime  = ", qflx_liqdew_to_top_layer(c)*dtime
          write(iulog,*) "qflx_liqevap_from_top_layer*dtime = ", qflx_liqevap_from_top_layer(c)*dtime
          call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, &
               msg="In UpdateState_TopLayerFluxes, h2osoi_liq has gone significantly negative")
       end if

    end do

  end subroutine UpdateState_TopLayerFluxes

  !-----------------------------------------------------------------------
  subroutine BulkFlux_SnowPercolation(bounds, num_snowc, filter_snowc, &
       dtime, snl, dz, frac_sno_eff, h2osoi_ice, h2osoi_liq, &
       qflx_snow_percolation)
    !
    ! !DESCRIPTION:
    ! Calculate liquid percolation through the snow pack, for bulk water
    !
    ! qflx_snow_percolation(c,j) gives the percolation out of the bottom of layer j, into
    ! the top of layer j+1.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_snowc
    integer, intent(in) :: filter_snowc(:)

    real(r8) , intent(in)    :: dtime                                    ! land model time step (sec)
    integer  , intent(in)    :: snl( bounds%begc: )                      ! negative number of snow layers
    real(r8) , intent(in)    :: dz( bounds%begc: , -nlevsno+1: )         ! layer depth (m)
    real(r8) , intent(in)    :: frac_sno_eff( bounds%begc: )             ! eff. fraction of ground covered by snow (0 to 1)
    real(r8) , intent(in)    :: h2osoi_ice( bounds%begc: , -nlevsno+1: ) ! ice lens (kg/m2)
    real(r8) , intent(in)    :: h2osoi_liq( bounds%begc: , -nlevsno+1: ) ! liquid water (kg/m2)

    real(r8) , intent(inout) :: qflx_snow_percolation( bounds%begc: , -nlevsno+1: ) ! liquid percolation out of the bottom of snow layer j (mm H2O /s)
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c
    integer  :: j
    real(r8) :: vol_liq(bounds%begc:bounds%endc,-nlevsno+1:0)      ! partial volume of liquid water in layer
    real(r8) :: vol_ice(bounds%begc:bounds%endc,-nlevsno+1:0)      ! partial volume of ice lens in layer
    real(r8) :: eff_porosity(bounds%begc:bounds%endc,-nlevsno+1:0) ! effective porosity = porosity - vol_ice

    character(len=*), parameter :: subname = 'BulkFlux_SnowPercolation'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(snl, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(dz, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(frac_sno_eff, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osoi_ice, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osoi_liq, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_snow_percolation, 1) == bounds%endc), sourcefile, __LINE__)

    ! Porosity and partial volume

    do j = -nlevsno+1, 0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             ! need to scale dz by frac_sno to convert to grid cell average depth
             vol_ice(c,j)      = min(1._r8, h2osoi_ice(c,j)/(dz(c,j)*frac_sno_eff(c)*denice))
             eff_porosity(c,j) = 1._r8 - vol_ice(c,j)
             vol_liq(c,j)      = min(eff_porosity(c,j),h2osoi_liq(c,j)/(dz(c,j)*frac_sno_eff(c)*denh2o))
          end if
       end do
    end do

    ! Calculate qflx_snow_percolation from each layer; only valid for j >= snl(c)+1
    !
    ! Capillary forces within snow are usually two or more orders of magnitude
    ! less than those of gravity. Only gravity terms are considered.
    ! The general expression for water flow is "K * ss**3", however,
    ! no effective parameterization for "K".  Thus, a very simple consideration
    ! (not physically based) is introduced:
    ! when the liquid water of layer exceeds the layer's holding
    ! capacity, the excess meltwater adds to the underlying neighbor layer.

    do j = -nlevsno+1, 0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             if (j <= -1) then
                ! No runoff over snow surface, just ponding on surface
                if (eff_porosity(c,j) < params_inst%wimp .OR. eff_porosity(c,j+1) < params_inst%wimp) then
                   qflx_snow_percolation(c,j) = 0._r8
                else
                   ! dz must be scaled by frac_sno to obtain gridcell average value
                   qflx_snow_percolation(c,j) = max(0._r8,(vol_liq(c,j) &
                        - params_inst%ssi*eff_porosity(c,j))*dz(c,j)*frac_sno_eff(c))
                   qflx_snow_percolation(c,j) = min(qflx_snow_percolation(c,j),(1._r8-vol_ice(c,j+1) &
                        - vol_liq(c,j+1))*dz(c,j+1)*frac_sno_eff(c))
                end if
             else
                qflx_snow_percolation(c,j) = max(0._r8,(vol_liq(c,j) &
                     - params_inst%ssi*eff_porosity(c,j))*dz(c,j)*frac_sno_eff(c))
             end if
             qflx_snow_percolation(c,j) = (qflx_snow_percolation(c,j)*1000._r8)/dtime
          end if
       end do
    end do

  end subroutine BulkFlux_SnowPercolation

  !-----------------------------------------------------------------------
  subroutine TracerFlux_SnowPercolation(bounds, num_snowc, filter_snowc, &
       snl, bulk_h2osoi_liq, bulk_qflx_snow_percolation, trac_h2osoi_liq, &
       trac_qflx_snow_percolation)
    !
    ! !DESCRIPTION:
    ! Calculate liquid percolation through the snow pack, for one tracer
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_snowc
    integer, intent(in) :: filter_snowc(:)

    ! For description of arguments, see comments in BulkFlux_SnowPercolation. Here,
    ! bulk_* variables refer to bulk water and trac_* variables refer to the given water
    ! tracer.
    integer  , intent(in)    :: snl( bounds%begc: )
    real(r8) , intent(in)    :: bulk_h2osoi_liq( bounds%begc: , -nlevsno+1: )
    real(r8) , intent(in)    :: bulk_qflx_snow_percolation( bounds%begc: , -nlevsno+1: )
    real(r8) , intent(in)    :: trac_h2osoi_liq( bounds%begc: , -nlevsno+1: )
    real(r8) , intent(inout) :: trac_qflx_snow_percolation( bounds%begc: , -nlevsno+1: )
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c
    integer :: j
    logical :: snow_layer_exists(bounds%begc:bounds%endc)  ! Whether the current snow layer exists in each column

    character(len=*), parameter :: subname = 'TracerFlux_SnowPercolation'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(snl, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(bulk_h2osoi_liq, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(bulk_qflx_snow_percolation, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(trac_h2osoi_liq, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(trac_qflx_snow_percolation, 1) == bounds%endc), sourcefile, __LINE__)

    associate( &
         begc => bounds%begc, &
         endc => bounds%endc  &
         )

    do j = -nlevsno+1, 0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             snow_layer_exists(c) = .true.
          else
             snow_layer_exists(c) = .false.
          end if
       end do

       call CalcTracerFromBulkMasked( &
            subgrid_level = subgrid_level_column, &
            lb            = begc, &
            num_pts       = num_snowc, &
            filter_pts    = filter_snowc, &
            mask_array    = snow_layer_exists(begc:endc), &
            bulk_source   = bulk_h2osoi_liq(begc:endc, j), &
            bulk_val      = bulk_qflx_snow_percolation(begc:endc, j), &
            tracer_source = trac_h2osoi_liq(begc:endc, j), &
            tracer_val    = trac_qflx_snow_percolation(begc:endc, j))

    end do

    end associate

  end subroutine TracerFlux_SnowPercolation

  !-----------------------------------------------------------------------
  subroutine UpdateState_SnowPercolation(bounds, num_snowc, filter_snowc, &
       dtime, snl, qflx_snow_percolation, h2osoi_liq)
    !
    ! !DESCRIPTION:
    ! Update h2osoi_liq for snow percolation, for bulk or one tracer
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_snowc
    integer, intent(in) :: filter_snowc(:)

    real(r8) , intent(in)    :: dtime                                    ! land model time step (sec)
    integer  , intent(in)    :: snl( bounds%begc: )                      ! negative number of snow layers
    real(r8) , intent(in)    :: qflx_snow_percolation( bounds%begc: , -nlevsno+1: ) ! liquid percolation out of the bottom of snow layer j (mm H2O /s)
    real(r8) , intent(inout) :: h2osoi_liq( bounds%begc: , -nlevsno+1: ) ! liquid water (kg/m2)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c
    integer :: j

    character(len=*), parameter :: subname = 'UpdateState_SnowPercolation'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(snl, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_snow_percolation, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osoi_liq, 1) == bounds%endc), sourcefile, __LINE__)

    do j = -nlevsno+1, 0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then

             ! For layers below the top layer, add percolation from layer above
             if (j >= snl(c)+2) then
                h2osoi_liq(c,j) = h2osoi_liq(c,j) + qflx_snow_percolation(c,j-1) * dtime
             end if

             ! Subtract percolation out of this layer
             h2osoi_liq(c,j) = h2osoi_liq(c,j) - qflx_snow_percolation(c,j) * dtime
          end if
       end do
    end do

  end subroutine UpdateState_SnowPercolation

  !-----------------------------------------------------------------------
  subroutine CalcAndApplyAerosolFluxes(bounds, num_snowc, filter_snowc, &
       dtime, snl, h2osoi_ice, h2osoi_liq, qflx_snow_percolation, atm2lnd_inst, &
       aerosol_inst)
    !
    ! !DESCRIPTION:
    ! Calculate and apply fluxes of aerosols through the snow pack
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_snowc
    integer, intent(in) :: filter_snowc(:)

    real(r8)           , intent(in)    :: dtime                                               ! land model time step (sec)
    integer            , intent(in)    :: snl( bounds%begc: )                                 ! negative number of snow layers
    real(r8)           , intent(in)    :: h2osoi_ice( bounds%begc: , -nlevsno+1: )            ! ice lens (kg/m2)
    real(r8)           , intent(in)    :: h2osoi_liq( bounds%begc: , -nlevsno+1: )            ! liquid water (kg/m2)
    real(r8)           , intent(in)    :: qflx_snow_percolation( bounds%begc: , -nlevsno+1: ) ! liquid percolation out of the bottom of snow layer j (mm H2O /s)
    type(atm2lnd_type) , intent(in)    :: atm2lnd_inst
    type(aerosol_type) , intent(inout) :: aerosol_inst
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c
    integer :: j
    real(r8) :: mss_liqice                                         ! mass of liquid+ice in a layer
    real(r8) :: qin_bc_phi  (bounds%begc:bounds%endc)              ! flux of hydrophilic BC into   layer [kg/s]
    real(r8) :: qout_bc_phi (bounds%begc:bounds%endc)              ! flux of hydrophilic BC out of layer [kg/s]
    real(r8) :: qin_bc_pho  (bounds%begc:bounds%endc)              ! flux of hydrophobic BC into   layer [kg/s]
    real(r8) :: qout_bc_pho (bounds%begc:bounds%endc)              ! flux of hydrophobic BC out of layer [kg/s]
    real(r8) :: qin_oc_phi  (bounds%begc:bounds%endc)              ! flux of hydrophilic OC into   layer [kg/s]
    real(r8) :: qout_oc_phi (bounds%begc:bounds%endc)              ! flux of hydrophilic OC out of layer [kg/s]
    real(r8) :: qin_oc_pho  (bounds%begc:bounds%endc)              ! flux of hydrophobic OC into   layer [kg/s]
    real(r8) :: qout_oc_pho (bounds%begc:bounds%endc)              ! flux of hydrophobic OC out of layer [kg/s]
    real(r8) :: qin_dst1    (bounds%begc:bounds%endc)              ! flux of dust species 1 into   layer [kg/s]
    real(r8) :: qout_dst1   (bounds%begc:bounds%endc)              ! flux of dust species 1 out of layer [kg/s]
    real(r8) :: qin_dst2    (bounds%begc:bounds%endc)              ! flux of dust species 2 into   layer [kg/s]
    real(r8) :: qout_dst2   (bounds%begc:bounds%endc)              ! flux of dust species 2 out of layer [kg/s]
    real(r8) :: qin_dst3    (bounds%begc:bounds%endc)              ! flux of dust species 3 into   layer [kg/s]
    real(r8) :: qout_dst3   (bounds%begc:bounds%endc)              ! flux of dust species 3 out of layer [kg/s]
    real(r8) :: qin_dst4    (bounds%begc:bounds%endc)              ! flux of dust species 4 into   layer [kg/s]
    real(r8) :: qout_dst4   (bounds%begc:bounds%endc)              ! flux of dust species 4 out of layer [kg/s]

    character(len=*), parameter :: subname = 'CalcAndApplyAerosolFluxes'
    !-----------------------------------------------------------------------

    associate( &
         mss_bcphi      => aerosol_inst%mss_bcphi_col        , & ! Output: [real(r8) (:,:) ] hydrophillic BC mass in snow (col,lyr) [kg]
         mss_bcpho      => aerosol_inst%mss_bcpho_col        , & ! Output: [real(r8) (:,:) ] hydrophobic  BC mass in snow (col,lyr) [kg]
         mss_ocphi      => aerosol_inst%mss_ocphi_col        , & ! Output: [real(r8) (:,:) ] hydrophillic OC mass in snow (col,lyr) [kg]
         mss_ocpho      => aerosol_inst%mss_ocpho_col        , & ! Output: [real(r8) (:,:) ] hydrophobic  OC mass in snow (col,lyr) [kg]
         mss_dst1       => aerosol_inst%mss_dst1_col         , & ! Output: [real(r8) (:,:) ] mass of dust species 1 in snow (col,lyr) [kg]
         mss_dst2       => aerosol_inst%mss_dst2_col         , & ! Output: [real(r8) (:,:) ] mass of dust species 2 in snow (col,lyr) [kg]
         mss_dst3       => aerosol_inst%mss_dst3_col         , & ! Output: [real(r8) (:,:) ] mass of dust species 3 in snow (col,lyr) [kg]
         mss_dst4       => aerosol_inst%mss_dst4_col           & ! Output: [real(r8) (:,:) ] mass of dust species 4 in snow (col,lyr) [kg]
         )

    ! Compute aerosol fluxes through snowpack:
    ! 1) compute aerosol mass in each layer
    ! 2) add aerosol mass flux from above layer to mass of this layer
    ! 3) qout_xxx is mass flux of aerosol species xxx out bottom of
    !    layer in water flow, proportional to (current) concentration
    !    of aerosol in layer multiplied by a scavenging ratio.
    ! 4) update mass of aerosol in top layer, accordingly
    ! 5) update mass concentration of aerosol accordingly

    do fc = 1, num_snowc
       c = filter_snowc(fc)

       qin_bc_phi (c) = 0._r8
       qin_bc_pho (c) = 0._r8
       qin_oc_phi (c) = 0._r8
       qin_oc_pho (c) = 0._r8
       qin_dst1   (c) = 0._r8
       qin_dst2   (c) = 0._r8
       qin_dst3   (c) = 0._r8
       qin_dst4   (c) = 0._r8
    end do

    do j = -nlevsno+1, 0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then

             mss_bcphi(c,j) = mss_bcphi(c,j) + qin_bc_phi(c) * dtime
             mss_bcpho(c,j) = mss_bcpho(c,j) + qin_bc_pho(c) * dtime
             mss_ocphi(c,j) = mss_ocphi(c,j) + qin_oc_phi(c) * dtime
             mss_ocpho(c,j) = mss_ocpho(c,j) + qin_oc_pho(c) * dtime

             mss_dst1(c,j)  = mss_dst1(c,j) + qin_dst1(c) * dtime
             mss_dst2(c,j)  = mss_dst2(c,j) + qin_dst2(c) * dtime
             mss_dst3(c,j)  = mss_dst3(c,j) + qin_dst3(c) * dtime
             mss_dst4(c,j)  = mss_dst4(c,j) + qin_dst4(c) * dtime

             ! mass of ice+water: in extremely rare circumstances, this can
             ! be zero, even though there is a snow layer defined. In
             ! this case, set the mass to a very small value to
             ! prevent division by zero.

             mss_liqice = h2osoi_liq(c,j)+h2osoi_ice(c,j)
             if (mss_liqice < 1E-30_r8) then
                mss_liqice = 1E-30_r8
             endif

             ! BCPHI:
             ! 1. flux with meltwater:
             qout_bc_phi(c) = qflx_snow_percolation(c,j)*params_inst%scvng_fct_mlt_sf* &
                              scvng_fct_mlt_bcphi*(mss_bcphi(c,j)/mss_liqice)
             if (qout_bc_phi(c)*dtime > mss_bcphi(c,j)) then
                qout_bc_phi(c) = mss_bcphi(c,j)/dtime
                mss_bcphi(c,j) = 0._r8
             else
                mss_bcphi(c,j) = mss_bcphi(c,j) - qout_bc_phi(c) * dtime
             end if
             qin_bc_phi(c) = qout_bc_phi(c)

             ! BCPHO:
             ! 1. flux with meltwater:
             qout_bc_pho(c) = qflx_snow_percolation(c,j)*params_inst%scvng_fct_mlt_sf* &
                              scvng_fct_mlt_bcpho*(mss_bcpho(c,j)/mss_liqice)
             if (qout_bc_pho(c)*dtime > mss_bcpho(c,j)) then
                qout_bc_pho(c) = mss_bcpho(c,j)/dtime
                mss_bcpho(c,j) = 0._r8
             else
                mss_bcpho(c,j) = mss_bcpho(c,j) - qout_bc_pho(c) * dtime
             end if
             qin_bc_pho(c) = qout_bc_pho(c)

             ! OCPHI:
             ! 1. flux with meltwater:
             qout_oc_phi(c) = qflx_snow_percolation(c,j)*params_inst%scvng_fct_mlt_sf* &
                              scvng_fct_mlt_ocphi*(mss_ocphi(c,j)/mss_liqice)
             if (qout_oc_phi(c)*dtime > mss_ocphi(c,j)) then
                qout_oc_phi(c) = mss_ocphi(c,j)/dtime
                mss_ocphi(c,j) = 0._r8
             else
                mss_ocphi(c,j) = mss_ocphi(c,j) - qout_oc_phi(c) * dtime
             end if
             qin_oc_phi(c) = qout_oc_phi(c)

             ! OCPHO:
             ! 1. flux with meltwater:
             qout_oc_pho(c) = qflx_snow_percolation(c,j)*params_inst%scvng_fct_mlt_sf* &
                              scvng_fct_mlt_ocpho*(mss_ocpho(c,j)/mss_liqice)
             if (qout_oc_pho(c)*dtime > mss_ocpho(c,j)) then
                qout_oc_pho(c) = mss_ocpho(c,j)/dtime
                mss_ocpho(c,j) = 0._r8
             else
                mss_ocpho(c,j) = mss_ocpho(c,j) - qout_oc_pho(c) * dtime
             end if
             qin_oc_pho(c) = qout_oc_pho(c)

             ! DUST 1:
             ! 1. flux with meltwater:
             qout_dst1(c) = qflx_snow_percolation(c,j)*params_inst%scvng_fct_mlt_sf* &
                            scvng_fct_mlt_dst1*(mss_dst1(c,j)/mss_liqice)
             if (qout_dst1(c)*dtime > mss_dst1(c,j)) then
                qout_dst1(c) = mss_dst1(c,j)/dtime
                mss_dst1(c,j) = 0._r8
             else
                mss_dst1(c,j) = mss_dst1(c,j) - qout_dst1(c) * dtime
             end if
             qin_dst1(c) = qout_dst1(c)

             ! DUST 2:
             ! 1. flux with meltwater:
             qout_dst2(c) = qflx_snow_percolation(c,j)*params_inst%scvng_fct_mlt_sf* &
                            scvng_fct_mlt_dst2*(mss_dst2(c,j)/mss_liqice)
             if (qout_dst2(c)*dtime > mss_dst2(c,j)) then
                qout_dst2(c) = mss_dst2(c,j)/dtime
                mss_dst2(c,j) = 0._r8
             else
                mss_dst2(c,j) = mss_dst2(c,j) - qout_dst2(c) * dtime
             end if
             qin_dst2(c) = qout_dst2(c)

             ! DUST 3:
             ! 1. flux with meltwater:
             qout_dst3(c) = qflx_snow_percolation(c,j)*params_inst%scvng_fct_mlt_sf* &
                            scvng_fct_mlt_dst3*(mss_dst3(c,j)/mss_liqice)
             if (qout_dst3(c)*dtime > mss_dst3(c,j)) then
                qout_dst3(c) = mss_dst3(c,j)/dtime
                mss_dst3(c,j) = 0._r8
             else
                mss_dst3(c,j) = mss_dst3(c,j) - qout_dst3(c) * dtime
             end if
             qin_dst3(c) = qout_dst3(c)

             ! DUST 4:
             ! 1. flux with meltwater:
             qout_dst4(c) = qflx_snow_percolation(c,j)*params_inst%scvng_fct_mlt_sf* &
                            scvng_fct_mlt_dst4*(mss_dst4(c,j)/mss_liqice)
             if (qout_dst4(c)*dtime > mss_dst4(c,j)) then
                qout_dst4(c) = mss_dst4(c,j)/dtime
                mss_dst4(c,j) = 0._r8
             else
                mss_dst4(c,j) = mss_dst4(c,j) - qout_dst4(c) * dtime
             end if
             qin_dst4(c) = qout_dst4(c)

          end if
       end do
    end do

    ! Compute aerosol fluxes through snowpack and aerosol deposition fluxes into top layer

    call AerosolFluxes(bounds, num_snowc, filter_snowc, &
         atm2lnd_inst, aerosol_inst)

    end associate

  end subroutine CalcAndApplyAerosolFluxes

  !-----------------------------------------------------------------------
  subroutine PostPercolation_AdjustLayerThicknesses(bounds, num_snowc, filter_snowc, &
       snl, h2osoi_ice, h2osoi_liq, dz)
    !
    ! !DESCRIPTION:
    ! Adjust layer thickness for any water+ice content changes after percolation through
    ! the snow pack
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_snowc
    integer, intent(in) :: filter_snowc(:)

    integer  , intent(in)    :: snl( bounds%begc: )                      ! negative number of snow layers
    real(r8) , intent(in)    :: h2osoi_ice( bounds%begc: , -nlevsno+1: ) ! ice lens (kg/m2)
    real(r8) , intent(in)    :: h2osoi_liq( bounds%begc: , -nlevsno+1: ) ! liquid water (kg/m2)
    real(r8) , intent(inout) :: dz( bounds%begc: , -nlevsno+1: )         ! layer depth (m)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c
    integer :: j

    character(len=*), parameter :: subname = 'PostPercolation_AdjustLayerThicknesses'
    !-----------------------------------------------------------------------

    ! Adjust layer thickness for any water+ice content changes in excess of previous
    ! layer thickness. Strictly speaking, only necessary for top snow layer, but doing
    ! it for all snow layers will catch problems with older initial files.
    ! Layer interfaces (zi) and node depths (z) do not need adjustment here because they
    ! are adjusted in CombineSnowLayers and are not used up to that point.

    do j = -nlevsno+1, 0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             dz(c,j) = max(dz(c,j),h2osoi_liq(c,j)/denh2o + h2osoi_ice(c,j)/denice)
          end if
       end do
    end do

  end subroutine PostPercolation_AdjustLayerThicknesses

  !-----------------------------------------------------------------------
  subroutine BulkDiag_SnowWaterAccumulatedSnow(bounds, &
       num_snowc, filter_snowc, num_nosnowc, filter_nosnowc, &
       dtime, frac_sno_eff, qflx_soliddew_to_top_layer, qflx_liqdew_to_top_layer, &
       qflx_liq_grnd, h2osno_no_layers, int_snow, frac_sno, snow_depth)
    !
    ! !DESCRIPTION:
    ! Update int_snow, and reset accumulated snow when no snow present
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_snowc
    integer, intent(in) :: filter_snowc(:)
    integer, intent(in) :: num_nosnowc
    integer, intent(in) :: filter_nosnowc(:)

    real(r8) , intent(in)    :: dtime                            ! land model time step (sec)
    real(r8) , intent(in)    :: frac_sno_eff( bounds%begc: )     ! eff. fraction of ground covered by snow (0 to 1)
    real(r8) , intent(in)    :: qflx_soliddew_to_top_layer( bounds%begc: ) ! rate of solid water deposited on top soil or snow layer (frost) (mm H2O /s)
    real(r8) , intent(in)    :: qflx_liqdew_to_top_layer( bounds%begc: ) ! rate of liquid water deposited on top soil or snow layer (dew) (mm H2O /s)
    real(r8) , intent(in)    :: qflx_liq_grnd( bounds%begc: )    ! liquid on ground after interception (mm H2O/s)
    real(r8) , intent(in)    :: h2osno_no_layers( bounds%begc: ) ! snow that is not resolved into layers (kg/m2)

    real(r8) , intent(inout) :: int_snow( bounds%begc: )         ! integrated snowfall (mm H2O)
    real(r8) , intent(inout) :: frac_sno( bounds%begc: )         ! fraction of ground covered by snow (0 to 1)
    real(r8) , intent(inout) :: snow_depth( bounds%begc: )       ! snow height (m)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c

    character(len=*), parameter :: subname = 'BulkDiag_SnowWaterAccumulatedSnow'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(frac_sno_eff, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_soliddew_to_top_layer, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_liqdew_to_top_layer, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_liq_grnd, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osno_no_layers, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(int_snow, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(frac_sno, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(snow_depth, 1) == bounds%endc), sourcefile, __LINE__)

    do fc = 1, num_snowc
       c = filter_snowc(fc)

       int_snow(c) = int_snow(c) + frac_sno_eff(c) &
            * (qflx_soliddew_to_top_layer(c) + qflx_liqdew_to_top_layer(c) &
            + qflx_liq_grnd(c)) * dtime
    end do

    do fc = 1, num_nosnowc
       c = filter_nosnowc(fc)

       ! reset accumulated snow when no snow present
       if (h2osno_no_layers(c) <= 0._r8) then
          int_snow(c) = 0._r8
          frac_sno(c) = 0._r8
          snow_depth(c) = 0._r8
       end if
    end do

  end subroutine BulkDiag_SnowWaterAccumulatedSnow

  !-----------------------------------------------------------------------
  subroutine SumFlux_AddSnowPercolation(bounds, &
       num_snowc, filter_snowc, num_nosnowc, filter_nosnowc, &
       frac_sno_eff, qflx_snow_percolation_bottom, qflx_liq_grnd, qflx_snomelt, &
       qflx_snow_drain, qflx_rain_plus_snomelt)
    !
    ! !DESCRIPTION:
    ! Calculate summed fluxes accounting for qflx_snow_percolation and similar fluxes
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_snowc
    integer, intent(in) :: filter_snowc(:)
    integer, intent(in) :: num_nosnowc
    integer, intent(in) :: filter_nosnowc(:)

    real(r8) , intent(in)    :: frac_sno_eff( bounds%begc: )                 ! eff. fraction of ground covered by snow (0 to 1)
    real(r8) , intent(in)    :: qflx_snow_percolation_bottom( bounds%begc: ) ! liquid percolation out of the bottom of the snow pack (mm H2O /s)
    real(r8) , intent(in)    :: qflx_liq_grnd( bounds%begc: )                ! liquid on ground after interception (mm H2O/s)
    real(r8) , intent(in)    :: qflx_snomelt( bounds%begc: )                 ! snow melt (mm H2O /s)

    real(r8) , intent(inout) :: qflx_snow_drain( bounds%begc: )              ! drainage from snow pack from previous time step (mm H2O/s)
    real(r8) , intent(inout) :: qflx_rain_plus_snomelt( bounds%begc: )       ! rain plus snow melt falling on the soil (mm/s)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c

    character(len=*), parameter :: subname = 'SumFlux_AddSnowPercolation'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(frac_sno_eff, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_snow_percolation_bottom, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_liq_grnd, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_snomelt, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_snow_drain, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_rain_plus_snomelt, 1) == bounds%endc), sourcefile, __LINE__)

    do fc = 1, num_snowc
       c = filter_snowc(fc)

       qflx_snow_drain(c) = qflx_snow_drain(c) + qflx_snow_percolation_bottom(c)
       qflx_rain_plus_snomelt(c) = qflx_snow_percolation_bottom(c) &
            + (1.0_r8 - frac_sno_eff(c)) * qflx_liq_grnd(c)
    end do

    do fc = 1, num_nosnowc
       c = filter_nosnowc(fc)

       qflx_snow_drain(c) = qflx_snomelt(c)
       qflx_rain_plus_snomelt(c) = qflx_liq_grnd(c) + qflx_snomelt(c)
    end do

  end subroutine SumFlux_AddSnowPercolation

  !-----------------------------------------------------------------------
  subroutine SnowCompaction(bounds, num_snowc, filter_snowc, &
       scf_method, &
       temperature_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, atm2lnd_inst)
    !
    ! !DESCRIPTION:
    ! Determine the change in snow layer thickness due to compaction and
    ! settling.
    ! Three metamorphisms of changing snow characteristics are implemented,
    ! i.e., destructive, overburden, and melt. The treatments of the former
    ! two are from SNTHERM.89 and SNTHERM.99 (1991, 1999). The contribution
    ! due to melt metamorphism is simply taken as a ratio of snow ice
    ! fraction after the melting versus before the melting.
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in) :: bounds
    integer                , intent(in) :: num_snowc       ! number of column snow points in column filter
    integer                , intent(in) :: filter_snowc(:) ! column filter for snow points
    class(snow_cover_fraction_base_type), intent(in) :: scf_method
    type(temperature_type) , intent(in) :: temperature_inst
    type(waterstatebulk_type)  , intent(in) :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)  , intent(in) :: waterdiagnosticbulk_inst
    type(atm2lnd_type)     , intent(in) :: atm2lnd_inst
    !
    ! !LOCAL VARIABLES:
    integer :: j, l, c, fc                      ! indices
    integer :: g                                ! gridcell index
    real(r8):: dtime                            ! land model time step (sec)
    ! parameters
    real(r8), parameter :: c3 = 2.777e-6_r8     ! [1/s]
    real(r8), parameter :: c4 = 0.04_r8         ! [1/K]
    real(r8), parameter :: c5 = 2.0_r8          !
    !
    real(r8) :: burden(bounds%begc:bounds%endc)  ! pressure of overlying snow [kg/m2]
    real(r8) :: zpseudo(bounds%begc:bounds%endc) ! wind drift compaction / pseudo depth (only valid if wind_dependent_snow_density is .true.)
    logical  :: mobile(bounds%begc:bounds%endc)  ! current snow layer is mobile, i.e. susceptible to wind drift (only valid if wind_dependent_snow_density is .true.)
    real(r8) :: ddz1   ! Rate of settling of snowpack due to destructive metamorphism.
    real(r8) :: ddz2   ! Rate of compaction of snowpack due to overburden.
    real(r8) :: ddz3   ! Rate of compaction of snowpack due to melt [1/s]
    real(r8) :: dexpf  ! expf=exp(-c4*(273.15-t_soisno)).
    real(r8) :: fi     ! Fraction of ice relative to the total water content at current time step
    real(r8) :: td     ! t_soisno - tfrz [K]
    real(r8) :: pdzdtc ! Nodal rate of change in fractional-thickness due to compaction [fraction/s]
    real(r8) :: void   ! void (1 - vol_ice - vol_liq)
    real(r8) :: wx     ! water mass (ice+liquid) [kg/m2]
    real(r8) :: bi     ! partial density of ice [kg/m3]
    real(r8) :: wsum   ! snowpack total water mass (ice+liquid) [kg/m2]
    real(r8) :: fsno_melt
    real(r8) :: ddz4   ! Rate of compaction of snowpack due to wind drift.
    !-----------------------------------------------------------------------

    associate( &
         snl          => col%snl                          , & ! Input:  [integer (:)    ] number of snow layers
         lakpoi       => lun%lakpoi                       , & ! Input:  [logical  (:)   ] true => landunit is a lake point
         urbpoi       => lun%urbpoi                       , & ! Input:  [logical  (:)   ] true => landunit is an urban point
         forc_wind    => atm2lnd_inst%forc_wind_grc       , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed (m/s)

         t_soisno     => temperature_inst%t_soisno_col    , & ! Input:  [real(r8) (:,:) ] soil temperature (Kelvin)
         imelt        => temperature_inst%imelt_col       , & ! Input:  [integer (:,:)  ] flag for melting (=1), freezing (=2), Not=0

         frac_sno     => waterdiagnosticbulk_inst%frac_sno_eff_col , & ! Input:  [real(r8) (:)   ] snow covered fraction
         frac_h2osfc  => waterdiagnosticbulk_inst%frac_h2osfc_col  , & ! Input:  [real(r8) (:)   ] fraction of ground covered by surface water (0 to 1)
         swe_old      => waterdiagnosticbulk_inst%swe_old_col      , & ! Input:  [real(r8) (:,:) ] initial swe values
         int_snow     => waterstatebulk_inst%int_snow_col     , & ! Input:  [real(r8) (:)   ] integrated snowfall [mm]
         frac_iceold  => waterdiagnosticbulk_inst%frac_iceold_col  , & ! Input:  [real(r8) (:,:) ] fraction of ice relative to the tot water
         h2osoi_ice   => waterstatebulk_inst%h2osoi_ice_col   , & ! Input:  [real(r8) (:,:) ] ice lens (kg/m2)
         h2osoi_liq   => waterstatebulk_inst%h2osoi_liq_col   , & ! Input:  [real(r8) (:,:) ] liquid water (kg/m2)

         dz           => col%dz                             & ! Output: [real(r8) (: ,:) ] layer depth (m)
    )

    ! Get time step

    dtime = get_step_size_real()

    ! Begin calculation - note that the following column loops are only invoked if snl(c) < 0

    do fc = 1, num_snowc
       c = filter_snowc(fc)
       
       burden(c)  = 0._r8
       zpseudo(c) = 0._r8
       mobile(c)  = .true.
    end do

    do j = -nlevsno+1, 0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          g = col%gridcell(c)
          if (j >= snl(c)+1) then

             wx = (h2osoi_ice(c,j) + h2osoi_liq(c,j))
             void = 1._r8 - (h2osoi_ice(c,j)/denice + h2osoi_liq(c,j)/denh2o)&
                  /(frac_sno(c) * dz(c,j))

             ! Allow compaction only for non-saturated node and higher ice lens node.
             if (void > 0.001_r8 .and. h2osoi_ice(c,j) > .1_r8) then

                bi = h2osoi_ice(c,j) / (frac_sno(c) * dz(c,j))
                fi = h2osoi_ice(c,j) / wx
                td = tfrz-t_soisno(c,j)
                dexpf = exp(-c4*td)

                ! Settling as a result of destructive metamorphism

                ddz1 = -c3*dexpf
                if (bi > upplim_destruct_metamorph) ddz1 = ddz1*exp(-46.0e-3_r8*(bi-upplim_destruct_metamorph))

                ! Liquid water term

                if (h2osoi_liq(c,j) > 0.01_r8*dz(c,j)*frac_sno(c)) ddz1=ddz1*c5

                select case (overburden_compaction_method)
                case (OverburdenCompactionMethodAnderson1976)
                   ddz2 = OverburdenCompactionAnderson1976( &
                        burden = burden(c), &
                        wx = wx, &
                        td = td, &
                        bi = bi)

                case (OverburdenCompactionMethodVionnet2012)
                   ddz2 = OverburdenCompactionVionnet2012( &
                        h2osoi_liq = h2osoi_liq(c,j), &
                        dz = dz(c,j), &
                        burden = burden(c), &
                        wx = wx, &
                        td = td, &
                        bi = bi)

                case default
                   call endrun(msg="Unknown overburden_compaction_method")
                end select

                ! Compaction occurring during melt

                if (imelt(c,j) == 1) then
                   l = col%landunit(c)
                   ! For consistency with other uses of use_subgrid_fluxes, we apply this
                   ! code over all landunits other than lake and urban. (Elsewhere, the
                   ! uses of use_subgrid_fluxes are in a nolake filter, and check
                   ! .not. urbpoi.)
                   if(use_subgrid_fluxes .and. (.not. lakpoi(l) .and. .not. urbpoi(l))) then
                      ! first term is delta mass over mass
                      ddz3 = max(0._r8,min(1._r8,(swe_old(c,j) - wx)/wx))

                      ! 2nd term is delta fsno over fsno, allowing for negative values for ddz3
                      if((swe_old(c,j) - wx) > 0._r8) then
                         wsum = sum(h2osoi_liq(c,snl(c)+1:0)+h2osoi_ice(c,snl(c)+1:0))
                         fsno_melt = scf_method%FracSnowDuringMelt( &
                              c            = c, &
                              h2osno_total = wsum, &
                              int_snow     = int_snow(c))
                         
                         ! Ensure sum of snow and surface water fractions are <= 1 after update
                         !
                         ! Note that there is a similar adjustment in subroutine
                         ! FracH2oSfc (related to frac_sno); these two should be kept in
                         ! sync (e.g., if a 3rd fraction is ever added in one place, it
                         ! needs to be added in the other place, too).
                         if ((fsno_melt + frac_h2osfc(c)) > 1._r8) then
                            fsno_melt = 1._r8 - frac_h2osfc(c)
                         end if

                         ddz3 = ddz3 - max(0._r8,(fsno_melt - frac_sno(c))/frac_sno(c))
                      endif
                      ddz3 = -1._r8/dtime * ddz3
                   else
                      ddz3 = - 1._r8/dtime * max(0._r8,(frac_iceold(c,j) - fi)/frac_iceold(c,j))
                   endif
                else
                   ddz3 = 0._r8
                end if

                ! Compaction occurring due to wind drift
                if (wind_dependent_snow_density) then
                   call WindDriftCompaction( &
                        bi = bi, &
                        forc_wind = forc_wind(g), &
                        dz = dz(c,j), &
                        zpseudo = zpseudo(c), &
                        mobile = mobile(c), &
                        compaction_rate = ddz4)
                else
                   ddz4 = 0.0_r8
                end if

                ! Time rate of fractional change in dz (units of s-1)
                pdzdtc = ddz1 + ddz2 + ddz3 + ddz4

                ! The change in dz due to compaction
                ! Limit compaction to be no greater than fully saturated layer thickness
                dz(c,j) = max(dz(c,j) * (1._r8+pdzdtc*dtime),(h2osoi_ice(c,j)/denice+ h2osoi_liq(c,j)/denh2o)/frac_sno(c))

             else
                ! saturated node is immobile
                !
                ! This is only needed if wind_dependent_snow_density is true, but it's
                ! simplest just to update mobile always
                mobile(c) = .false.
             end if

             ! Pressure of overlying snow

             burden(c) = burden(c) + wx

          end if
       end do
    end do

    end associate
  end subroutine SnowCompaction

  !-----------------------------------------------------------------------
  subroutine CombineSnowLayers(bounds, num_snowc, filter_snowc, &
        aerosol_inst, temperature_inst, water_inst)
    !
    ! !DESCRIPTION:
    ! Combine snow layers that are less than a minimum thickness or mass
    ! If the snow element thickness or mass is less than a prescribed minimum,
    ! then it is combined with a neighboring element.  The subroutine
    ! clm\_combo.f90 then executes the combination of mass and energy.
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(inout) :: num_snowc       ! number of column snow points in column filter
    integer                , intent(inout) :: filter_snowc(:) ! column filter for snow points
    type(aerosol_type)     , intent(inout) :: aerosol_inst
    type(temperature_type) , intent(inout) :: temperature_inst
    type(water_type)       , intent(inout) :: water_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c, fc                            ! column indices
    integer :: i,k                              ! loop indices
    integer :: j,l                              ! node indices
    integer :: wi                               ! index of water tracer or bulk
    integer :: msn_old(bounds%begc:bounds%endc) ! number of top snow layer
    integer :: mssi(bounds%begc:bounds%endc)    ! node index
    integer :: neibor                           ! adjacent node selected for combination
    real(r8):: h2osno_total(bounds%begc:bounds%endc) ! total snow water (mm H2O)
    real(r8):: zwice(water_inst%bulk_and_tracers_beg:water_inst%bulk_and_tracers_end, bounds%begc:bounds%endc)  ! total ice mass in snow, for bulk and each tracer
    real(r8):: zwliq(water_inst%bulk_and_tracers_beg:water_inst%bulk_and_tracers_end, bounds%begc:bounds%endc)  ! total liquid water in snow, for bulk and each tracer
    real(r8):: dzminloc(nlevsno)  ! minimum thickness of snow layer (local)
    real(r8):: dtime                            !land model time step (sec)

    !-----------------------------------------------------------------------

    ! In contrast to most routines, this one operates on a mix of bulk-only quantities
    ! and bulk-and-tracer quantities. Variable names like h2osoi_liq_bulk refer to bulk
    ! quantities. Where bulk-and-tracer quantities are referenced, they are referred to
    ! like w%waterstate_inst%h2osoi_liq_col.

    associate( &
         b_waterstate_inst      => water_inst%waterstatebulk_inst, &
         b_waterdiagnostic_inst => water_inst%waterdiagnosticbulk_inst &
         )

    associate( &
         ltype            => lun%itype                           , & ! Input:  [integer  (:)   ] landunit type
         urbpoi           => lun%urbpoi                          , & ! Input:  [logical  (:)   ] true => landunit is an urban point

         t_soisno         => temperature_inst%t_soisno_col       , & ! Output: [real(r8) (:,:) ] soil temperature (Kelvin)

         mss_bcphi        => aerosol_inst%mss_bcphi_col          , & ! Output: [real(r8) (:,:) ] hydrophilic BC mass in snow (col,lyr) [kg]
         mss_bcpho        => aerosol_inst%mss_bcpho_col          , & ! Output: [real(r8) (:,:) ] hydrophobic BC mass in snow (col,lyr) [kg]
         mss_ocphi        => aerosol_inst%mss_ocphi_col          , & ! Output: [real(r8) (:,:) ] hydrophilic OC mass in snow (col,lyr) [kg]
         mss_ocpho        => aerosol_inst%mss_ocpho_col          , & ! Output: [real(r8) (:,:) ] hydrophobic OC mass in snow (col,lyr) [kg]
         mss_dst1         => aerosol_inst%mss_dst1_col           , & ! Output: [real(r8) (:,:) ] dust species 1 mass in snow (col,lyr) [kg]
         mss_dst2         => aerosol_inst%mss_dst2_col           , & ! Output: [real(r8) (:,:) ] dust species 2 mass in snow (col,lyr) [kg]
         mss_dst3         => aerosol_inst%mss_dst3_col           , & ! Output: [real(r8) (:,:) ] dust species 3 mass in snow (col,lyr) [kg]
         mss_dst4         => aerosol_inst%mss_dst4_col           , & ! Output: [real(r8) (:,:) ] dust species 4 mass in snow (col,lyr) [kg]

         frac_sno         => b_waterdiagnostic_inst%frac_sno_col        , & ! Input:  [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)
         frac_sno_eff     => b_waterdiagnostic_inst%frac_sno_eff_col    , & ! Input:  [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)
         snow_depth       => b_waterdiagnostic_inst%snow_depth_col      , & ! Output: [real(r8) (:)   ] snow height (m)
         int_snow         => b_waterstate_inst%int_snow_col        , & ! Output:  [real(r8) (:)   ] integrated snowfall [mm]
         snw_rds          => b_waterdiagnostic_inst%snw_rds_col         , & ! Output: [real(r8) (:,:) ] effective snow grain radius (col,lyr) [microns, m^-6]

         ! The following associates, with suffix _bulk, refer to bulk water. This is to
         ! distinguish them from references in this routine like
         ! w%waterstate_inst%h2osoi_ice_col, which refer to the current bulk or tracer.
         h2osno_no_layers_bulk => b_waterstate_inst%h2osno_no_layers_col, & ! Output: [real(r8) (:)   ]  snow that is not resolved into layers (kg/m2)
         h2osoi_ice_bulk  => b_waterstate_inst%h2osoi_ice_col      , & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)
         h2osoi_liq_bulk  => b_waterstate_inst%h2osoi_liq_col      , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)

         snl              => col%snl                             , & ! Output: [integer  (:)   ] number of snow layers
         dz               => col%dz                              , & ! Output: [real(r8) (:,:) ] layer depth (m)
         zi               => col%zi                              , & ! Output: [real(r8) (:,:) ] interface level below a "z" level (m)
         z                => col%z                                 & ! Output: [real(r8) (:,:) ] layer thickness (m)
    )

    ! Determine model time step

    dtime = get_step_size_real()

    ! Check the mass of ice lens of snow, when the total is less than a small value,
    ! combine it with the underlying neighbor.

    dzminloc(:) = dzmin(:) ! dzmin will stay constant between timesteps

    ! Add lsadz to dzmin for lakes
    ! Determine whether called from LakeHydrology
    ! Note: this assumes that this function is called separately with the lake-snow and non-lake-snow filters.
    if (num_snowc > 0) then
       c = filter_snowc(1)
       l = col%landunit(c)
       if (ltype(l) == istdlak) then ! Called from LakeHydrology
          dzminloc(:) = dzmin(:) + lsadz
       end if
    end if

    do fc = 1, num_snowc
       c = filter_snowc(fc)

       msn_old(c) = snl(c)

       do wi = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
          associate(w => water_inst%bulk_and_tracers(wi))

          w%waterflux_inst%qflx_sl_top_soil_col(c) = 0._r8

          end associate
       end do
    end do

    ! The following loop is NOT VECTORIZED

    do fc = 1, num_snowc
       c = filter_snowc(fc)
       l = col%landunit(c)
       do j = msn_old(c)+1,0
          ! use 0.01 to avoid runaway ice buildup
          if (h2osoi_ice_bulk(c,j) <= .01_r8) then
             if (j < 0 .or. (ltype(l) == istsoil .or. urbpoi(l) .or. ltype(l) == istcrop)) then
                ! Note that, for landunits other than soil, crop and urban, the above
                ! conditional prevents us from trying to transfer the bottom snow layer's
                ! water content to the soil, since there is no soil to receive ti.

                do wi = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
                   associate(w => water_inst%bulk_and_tracers(wi))

                   w%waterstate_inst%h2osoi_liq_col(c,j+1) = &
                        w%waterstate_inst%h2osoi_liq_col(c,j+1) + &
                        w%waterstate_inst%h2osoi_liq_col(c,j)

                   w%waterstate_inst%h2osoi_ice_col(c,j+1) = &
                        w%waterstate_inst%h2osoi_ice_col(c,j+1) + &
                        w%waterstate_inst%h2osoi_ice_col(c,j)

                   end associate
                end do
             end if

             if (j < 0) then
                dz(c,j+1) = dz(c,j+1) + dz(c,j)

                mss_bcphi(c,j+1) = mss_bcphi(c,j+1)  + mss_bcphi(c,j)
                mss_bcpho(c,j+1) = mss_bcpho(c,j+1)  + mss_bcpho(c,j)
                mss_ocphi(c,j+1) = mss_ocphi(c,j+1)  + mss_ocphi(c,j)
                mss_ocpho(c,j+1) = mss_ocpho(c,j+1)  + mss_ocpho(c,j)
                mss_dst1(c,j+1)  = mss_dst1(c,j+1)   + mss_dst1(c,j)
                mss_dst2(c,j+1)  = mss_dst2(c,j+1)   + mss_dst2(c,j)
                mss_dst3(c,j+1)  = mss_dst3(c,j+1)   + mss_dst3(c,j)
                mss_dst4(c,j+1)  = mss_dst4(c,j+1)   + mss_dst4(c,j)

                ! NOTE: Temperature, and similarly snw_rds, of the
                ! underlying snow layer are NOT adjusted in this case.
                ! Because the layer being eliminated has a small mass,
                ! this should not make a large difference, but it
                ! would be more thorough to do so.
             end if

             if (j == 0) then

                do wi = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
                   associate(w => water_inst%bulk_and_tracers(wi))

                   w%waterflux_inst%qflx_sl_top_soil_col(c) = &
                        (w%waterstate_inst%h2osoi_liq_col(c,j) + w%waterstate_inst%h2osoi_ice_col(c,j))/dtime

                   end associate
                end do
             end if

             ! shift all elements above this down one.
             if (j > snl(c)+1 .and. snl(c) < -1) then
                do i = j, snl(c)+2, -1
                   do wi = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
                      associate(w => water_inst%bulk_and_tracers(wi))

                      w%waterstate_inst%h2osoi_liq_col(c,i) = w%waterstate_inst%h2osoi_liq_col(c,i-1)
                      w%waterstate_inst%h2osoi_ice_col(c,i) = w%waterstate_inst%h2osoi_ice_col(c,i-1)

                      end associate
                   end do

                   t_soisno(c,i)   = t_soisno(c,i-1)

                   mss_bcphi(c,i)   = mss_bcphi(c,i-1)
                   mss_bcpho(c,i)   = mss_bcpho(c,i-1)
                   mss_ocphi(c,i)   = mss_ocphi(c,i-1)
                   mss_ocpho(c,i)   = mss_ocpho(c,i-1)
                   mss_dst1(c,i)    = mss_dst1(c,i-1)
                   mss_dst2(c,i)    = mss_dst2(c,i-1)
                   mss_dst3(c,i)    = mss_dst3(c,i-1)
                   mss_dst4(c,i)    = mss_dst4(c,i-1)
                   snw_rds(c,i)     = snw_rds(c,i-1)

                   dz(c,i)         = dz(c,i-1)
                end do
             end if
             snl(c) = snl(c) + 1
          end if
       end do
    end do

    do fc = 1, num_snowc
       c = filter_snowc(fc)
       snow_depth(c) = 0._r8
       ! See note in the following loop regarding why we are setting h2osno_total inline
       ! rather than relying on CalculateTotalH2osno.
       h2osno_total(c) = 0._r8

       do wi = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
          zwice(wi,c)  = 0._r8
          zwliq(wi,c)  = 0._r8
       end do

    end do

    do j = -nlevsno+1,0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             do wi = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
                associate(w => water_inst%bulk_and_tracers(wi))

                zwice(wi,c)  = zwice(wi,c) + w%waterstate_inst%h2osoi_ice_col(c,j)
                zwliq(wi,c)  = zwliq(wi,c) + w%waterstate_inst%h2osoi_liq_col(c,j)

                end associate
             end do

             snow_depth(c) = snow_depth(c) + dz(c,j)
             ! We generally compute h2osno_total with CalculateTotalH2osno. Here we
             ! calculate it inline for two reasons: (1) we're calculating other related
             ! variables here anyway; and (2) the consistency checks invoked by
             ! CalculateTotalH2osno in debug mode will sometimes fail here because we
             ! haven't yet zeroed out layers that just disappeared. Because of (2), if we
             ! wanted to use that routine to calculate h2osno_total here, we would need
             ! to first call ZeroEmptySnowLayers.
             h2osno_total(c) = h2osno_total(c) + h2osoi_ice_bulk(c,j) + h2osoi_liq_bulk(c,j)
          end if
       end do
    end do

    ! Check the snow depth - all snow gone
    ! The liquid water assumes ponding on soil surface.

    do fc = 1, num_snowc
       c = filter_snowc(fc)
       l = col%landunit(c)
       if (snow_depth(c) > 0._r8) then
          if ((ltype(l) == istdlak .and. snow_depth(c) < dzmin(1) + lsadz ) .or. &
               ((ltype(l) /= istdlak) .and. ((frac_sno_eff(c)*snow_depth(c) < dzmin(1))  &
               .or. (h2osno_total(c)/(frac_sno_eff(c)*snow_depth(c)) < 50._r8)))) then


             do wi = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
                associate(w => water_inst%bulk_and_tracers(wi))

                ! The explicit snow pack is disappearing. Transfer ice to
                ! h2osno_no_layers and (for soil landunits) transfer liquid water from
                ! snow pack to layer 1 (soil).

                w%waterstate_inst%h2osno_no_layers_col(c) = zwice(wi,c)
                if (ltype(l) == istsoil .or. urbpoi(l) .or. ltype(l) == istcrop) then
                   w%waterstate_inst%h2osoi_liq_col(c,1) = &
                        w%waterstate_inst%h2osoi_liq_col(c,1) + zwliq(wi,c)
                end if

                end associate
             end do

             snl(c) = 0
             h2osno_total(c) = h2osno_no_layers_bulk(c)

             mss_bcphi(c,:) = 0._r8
             mss_bcpho(c,:) = 0._r8
             mss_ocphi(c,:) = 0._r8
             mss_ocpho(c,:) = 0._r8
             mss_dst1(c,:)  = 0._r8
             mss_dst2(c,:)  = 0._r8
             mss_dst3(c,:)  = 0._r8
             mss_dst4(c,:)  = 0._r8

             if (h2osno_no_layers_bulk(c) <= 0._r8) then
                snow_depth(c) = 0._r8
             end if

          endif
       end if
       if (h2osno_total(c) <= 0._r8) then
          snow_depth(c) = 0._r8
          frac_sno(c) = 0._r8
          frac_sno_eff(c) = 0._r8
          int_snow(c) = 0._r8
       endif
    end do

    ! Check the snow depth - snow layers combined
    ! The following loop IS NOT VECTORIZED

    do fc = 1, num_snowc
       c = filter_snowc(fc)

       ! Two or more layers

       if (snl(c) < -1) then

          msn_old(c) = snl(c)
          mssi(c) = 1

          do i = msn_old(c)+1,0
             if ((frac_sno_eff(c)*dz(c,i) < dzminloc(mssi(c))) .or. &
                  ((h2osoi_ice_bulk(c,i) + h2osoi_liq_bulk(c,i))/(frac_sno_eff(c)*dz(c,i)) < 50._r8)) then
                if (i == snl(c)+1) then
                   ! If top node is removed, combine with bottom neighbor.
                   neibor = i + 1
                else if (i == 0) then
                   ! If the bottom neighbor is not snow, combine with the top neighbor.
                   neibor = i - 1
                else
                   ! If none of the above special cases apply, combine with the thinnest neighbor
                   neibor = i + 1
                   if ((dz(c,i-1)+dz(c,i)) < (dz(c,i+1)+dz(c,i))) neibor = i-1
                end if

                ! Node l and j are combined and stored as node j.
                if (neibor > i) then
                   j = neibor
                   l = i
                else
                   j = i
                   l = neibor
                end if

                ! this should be included in 'Combo' for consistency,
                ! but functionally it is the same to do it here
                mss_bcphi(c,j)=mss_bcphi(c,j)+mss_bcphi(c,l)
                mss_bcpho(c,j)=mss_bcpho(c,j)+mss_bcpho(c,l)
                mss_ocphi(c,j)=mss_ocphi(c,j)+mss_ocphi(c,l)
                mss_ocpho(c,j)=mss_ocpho(c,j)+mss_ocpho(c,l)
                mss_dst1(c,j)=mss_dst1(c,j)+mss_dst1(c,l)
                mss_dst2(c,j)=mss_dst2(c,j)+mss_dst2(c,l)
                mss_dst3(c,j)=mss_dst3(c,j)+mss_dst3(c,l)
                mss_dst4(c,j)=mss_dst4(c,j)+mss_dst4(c,l)

                ! mass-weighted combination of effective grain size:
                snw_rds(c,j) = (snw_rds(c,j)*(h2osoi_liq_bulk(c,j)+h2osoi_ice_bulk(c,j)) + &
                     snw_rds(c,l)*(h2osoi_liq_bulk(c,l)+h2osoi_ice_bulk(c,l))) / &
                     (h2osoi_liq_bulk(c,j)+h2osoi_ice_bulk(c,j)+h2osoi_liq_bulk(c,l)+h2osoi_ice_bulk(c,l))

                call Combo (dz(c,j), h2osoi_liq_bulk(c,j), h2osoi_ice_bulk(c,j), &
                     t_soisno(c,j), dz(c,l), h2osoi_liq_bulk(c,l), h2osoi_ice_bulk(c,l), t_soisno(c,l) )

                ! Bulk already combined in Combo; here we just need to loop over tracers
                ! and do a similar combination for them.
                do wi = water_inst%tracers_beg, water_inst%tracers_end
                   associate(w => water_inst%bulk_and_tracers(wi))

                   w%waterstate_inst%h2osoi_ice_col(c,j) = &
                        w%waterstate_inst%h2osoi_ice_col(c,j) + w%waterstate_inst%h2osoi_ice_col(c,l)
                   w%waterstate_inst%h2osoi_liq_col(c,j) = &
                        w%waterstate_inst%h2osoi_liq_col(c,j) + w%waterstate_inst%h2osoi_liq_col(c,l)

                   end associate
                end do

                ! Now shift all elements above this down one.
                if (j-1 > snl(c)+1) then

                   do k = j-1, snl(c)+2, -1
                      do wi = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
                         associate(w => water_inst%bulk_and_tracers(wi))

                         w%waterstate_inst%h2osoi_ice_col(c,k) = w%waterstate_inst%h2osoi_ice_col(c,k-1)
                         w%waterstate_inst%h2osoi_liq_col(c,k) = w%waterstate_inst%h2osoi_liq_col(c,k-1)

                         end associate
                      end do

                      t_soisno(c,k) = t_soisno(c,k-1)

                      mss_bcphi(c,k) = mss_bcphi(c,k-1)
                      mss_bcpho(c,k) = mss_bcpho(c,k-1)
                      mss_ocphi(c,k) = mss_ocphi(c,k-1)
                      mss_ocpho(c,k) = mss_ocpho(c,k-1)
                      mss_dst1(c,k)  = mss_dst1(c,k-1)
                      mss_dst2(c,k)  = mss_dst2(c,k-1)
                      mss_dst3(c,k)  = mss_dst3(c,k-1)
                      mss_dst4(c,k)  = mss_dst4(c,k-1)
                      snw_rds(c,k)   = snw_rds(c,k-1)

                      dz(c,k) = dz(c,k-1)
                   end do
                end if

                ! Decrease the number of snow layers
                snl(c) = snl(c) + 1
                if (snl(c) >= -1) EXIT

             else

                ! The layer thickness is greater than the prescribed minimum value
                mssi(c) = mssi(c) + 1

             end if
          end do

       end if

    end do

    ! Reset the node depth and the depth of layer interface

    do j = 0, -nlevsno+1, -1
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c) + 1) then
             z(c,j) = zi(c,j) - 0.5_r8*dz(c,j)
             zi(c,j-1) = zi(c,j) - dz(c,j)
          end if
       end do
    end do

    end associate
    end associate
  end subroutine CombineSnowLayers

  !-----------------------------------------------------------------------
  subroutine DivideSnowLayers(bounds, num_snowc, filter_snowc, &
        aerosol_inst, temperature_inst, water_inst, is_lake)
    !
    ! !DESCRIPTION:
    ! Subdivides snow layers if they exceed their prescribed maximum thickness.
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_snowc       ! number of column snow points in column filter
    integer                , intent(in)    :: filter_snowc(:) ! column filter for snow points
    type(aerosol_type)     , intent(inout) :: aerosol_inst
    type(temperature_type) , intent(inout) :: temperature_inst
    type(water_type)       , intent(inout) :: water_inst
    logical                , intent(in)    :: is_lake  !TODO - this should be examined and removed in the future
    !
    ! !LOCAL VARIABLES:
    integer  :: j, c, fc, k                              ! indices
    integer  :: wi                                       ! index of water tracer or bulk
    integer  :: i_bulk                                   ! index of bulk water
    real(r8) :: drr                                      ! thickness of the combined [m]
    integer  :: msno                                     ! number of snow layer 1 (top) to msno (bottom)
    real(r8) :: dzsno(bounds%begc:bounds%endc,nlevsno)   ! Snow layer thickness [m]
    real(r8) :: swice(water_inst%bulk_and_tracers_beg:water_inst%bulk_and_tracers_end, bounds%begc:bounds%endc, nlevsno) ! Partial volume of ice, for bulk and each tracer [m3/m3]
    real(r8) :: swliq(water_inst%bulk_and_tracers_beg:water_inst%bulk_and_tracers_end, bounds%begc:bounds%endc, nlevsno) ! Partial volume of liquid water, for bulk and each tracer [m3/m3]
    real(r8) :: tsno(bounds%begc:bounds%endc ,nlevsno)   ! Nodel temperature [K]
    real(r8) :: zwice(water_inst%bulk_and_tracers_beg:water_inst%bulk_and_tracers_end) ! temporary
    real(r8) :: zwliq(water_inst%bulk_and_tracers_beg:water_inst%bulk_and_tracers_end) ! temporary
    real(r8) :: propor                                   ! temporary
    real(r8) :: dtdz                                     ! temporary
    ! temporary variables mimicking the structure of other layer division variables
    real(r8) :: mbc_phi(bounds%begc:bounds%endc,nlevsno) ! mass of BC in each snow layer
    real(r8) :: zmbc_phi                                 ! temporary
    real(r8) :: mbc_pho(bounds%begc:bounds%endc,nlevsno) ! mass of BC in each snow layer
    real(r8) :: zmbc_pho                                 ! temporary
    real(r8) :: moc_phi(bounds%begc:bounds%endc,nlevsno) ! mass of OC in each snow layer
    real(r8) :: zmoc_phi                                 ! temporary
    real(r8) :: moc_pho(bounds%begc:bounds%endc,nlevsno) ! mass of OC in each snow layer
    real(r8) :: zmoc_pho                                 ! temporary
    real(r8) :: mdst1(bounds%begc:bounds%endc,nlevsno)   ! mass of dust 1 in each snow layer
    real(r8) :: zmdst1                                   ! temporary
    real(r8) :: mdst2(bounds%begc:bounds%endc,nlevsno)   ! mass of dust 2 in each snow layer
    real(r8) :: zmdst2                                   ! temporary
    real(r8) :: mdst3(bounds%begc:bounds%endc,nlevsno)   ! mass of dust 3 in each snow layer
    real(r8) :: zmdst3                                   ! temporary
    real(r8) :: mdst4(bounds%begc:bounds%endc,nlevsno)   ! mass of dust 4 in each snow layer
    real(r8) :: zmdst4                                   ! temporary
    real(r8) :: rds(bounds%begc:bounds%endc,nlevsno)
    ! Variables for consistency check
    real(r8) :: dztot(bounds%begc:bounds%endc)
    real(r8) :: snwicetot(water_inst%bulk_and_tracers_beg:water_inst%bulk_and_tracers_end, bounds%begc:bounds%endc)
    real(r8) :: snwliqtot(water_inst%bulk_and_tracers_beg:water_inst%bulk_and_tracers_end, bounds%begc:bounds%endc)
    real(r8) :: offset ! temporary
    !-----------------------------------------------------------------------

    ! In contrast to most routines, this one operates on a mix of bulk-only quantities and
    ! bulk-and-tracer quantities. Where bulk-and-tracer quantities are referenced, they
    ! are referred to like w%waterstate_inst%h2osoi_liq_col.

    associate( &
         b_waterdiagnostic_inst => water_inst%waterdiagnosticbulk_inst, &
         i_bulk                 => water_inst%i_bulk &
         )

    associate( &
         t_soisno   => temperature_inst%t_soisno_col    , & ! Output: [real(r8) (:,:) ] soil temperature (Kelvin)

         frac_sno   => b_waterdiagnostic_inst%frac_sno_eff_col , & ! Output: [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)
         snw_rds    => b_waterdiagnostic_inst%snw_rds_col      , & ! Output: [real(r8) (:,:) ] effective snow grain radius (col,lyr) [microns, m^-6]

         mss_bcphi  => aerosol_inst%mss_bcphi_col       , & ! Output: [real(r8) (:,:) ] hydrophilic BC mass in snow (col,lyr) [kg]
         mss_bcpho  => aerosol_inst%mss_bcpho_col       , & ! Output: [real(r8) (:,:) ] hydrophobic BC mass in snow (col,lyr) [kg]
         mss_ocphi  => aerosol_inst%mss_ocphi_col       , & ! Output: [real(r8) (:,:) ] hydrophilic OC mass in snow (col,lyr) [kg]
         mss_ocpho  => aerosol_inst%mss_ocpho_col       , & ! Output: [real(r8) (:,:) ] hydrophobic OC mass in snow (col,lyr) [kg]
         mss_dst1   => aerosol_inst%mss_dst1_col        , & ! Output: [real(r8) (:,:) ] dust species 1 mass in snow (col,lyr) [kg]
         mss_dst2   => aerosol_inst%mss_dst2_col        , & ! Output: [real(r8) (:,:) ] dust species 2 mass in snow (col,lyr) [kg]
         mss_dst3   => aerosol_inst%mss_dst3_col        , & ! Output: [real(r8) (:,:) ] dust species 3 mass in snow (col,lyr) [kg]
         mss_dst4   => aerosol_inst%mss_dst4_col        , & ! Output: [real(r8) (:,:) ] dust species 4 mass in snow (col,lyr) [kg]

         snl        => col%snl                          , & ! Output: [integer  (:)   ] number of snow layers
         dz         => col%dz                           , & ! Output: [real(r8) (:,:) ] layer depth (m)
         zi         => col%zi                           , & ! Output: [real(r8) (:,:) ] interface level below a "z" level (m)
         z          => col%z                              & ! Output: [real(r8) (:,:) ] layer thickness (m)
    )

    if ( is_lake ) then
       ! Initialize for consistency check
       do j = -nlevsno+1,0
          do fc = 1, num_snowc
             c = filter_snowc(fc)

             if (j == -nlevsno+1) then
                dztot(c) = 0._r8

                do wi = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
                   snwicetot(wi,c) = 0._r8
                   snwliqtot(wi,c) = 0._r8
                end do
             end if

             if (j >= snl(c)+1) then
                dztot(c) = dztot(c) + dz(c,j)

                do wi = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
                   associate(w => water_inst%bulk_and_tracers(wi))

                   snwicetot(wi,c) = snwicetot(wi,c) + w%waterstate_inst%h2osoi_ice_col(c,j)
                   snwliqtot(wi,c) = snwliqtot(wi,c) + w%waterstate_inst%h2osoi_liq_col(c,j)

                   end associate
                end do
             end if
          end do
       end do
    end if

    ! Begin calculation - note that the following column loops are only invoked
    ! for snow-covered columns

    do j = 1,nlevsno
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j <= abs(snl(c))) then
             if (is_lake) then
                dzsno(c,j) = dz(c,j+snl(c))
             else
                dzsno(c,j) = frac_sno(c)*dz(c,j+snl(c))
             end if

             do wi = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
                associate(w => water_inst%bulk_and_tracers(wi))

                swice(wi,c,j) = w%waterstate_inst%h2osoi_ice_col(c,j+snl(c))
                swliq(wi,c,j) = w%waterstate_inst%h2osoi_liq_col(c,j+snl(c))

                end associate
             end do

             tsno(c,j)  = t_soisno(c,j+snl(c))

             mbc_phi(c,j) = mss_bcphi(c,j+snl(c))
             mbc_pho(c,j) = mss_bcpho(c,j+snl(c))
             moc_phi(c,j) = mss_ocphi(c,j+snl(c))
             moc_pho(c,j) = mss_ocpho(c,j+snl(c))
             mdst1(c,j)   = mss_dst1(c,j+snl(c))
             mdst2(c,j)   = mss_dst2(c,j+snl(c))
             mdst3(c,j)   = mss_dst3(c,j+snl(c))
             mdst4(c,j)   = mss_dst4(c,j+snl(c))
             rds(c,j)     = snw_rds(c,j+snl(c))
          end if
       end do
    end do

    loop_snowcolumns: do fc = 1, num_snowc
       c = filter_snowc(fc)

       msno = abs(snl(c))

       ! Now traverse layers from top to bottom in a dynamic way, as the total
       ! number of layers (msno) may increase during the loop.
       ! Impose k < nlevsno; the special case 'k == nlevsno' is not relevant,
       ! as it is neither allowed to subdivide nor does it have layers below.
       k = 1
       loop_layers: do while( k <= msno .and. k < nlevsno )

          ! Current layer is bottom layer
          if (k == msno) then

             if (is_lake) then
                offset = 2._r8 * lsadz
             else
                offset = 0._r8
             end if

             if (dzsno(c,k) > dzmax_l(k) + offset) then
                ! Subdivide layer into two layers with equal thickness, water
                ! content, ice content and temperature
                msno = msno + 1
                dzsno(c,k)     = dzsno(c,k) / 2.0_r8
                dzsno(c,k+1)   = dzsno(c,k)

                do wi = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
                   swice(wi,c,k)     = swice(wi,c,k) / 2.0_r8
                   swice(wi,c,k+1)   = swice(wi,c,k)
                   swliq(wi,c,k)     = swliq(wi,c,k) / 2.0_r8
                   swliq(wi,c,k+1)   = swliq(wi,c,k)
                end do

                if (k == 1) then
                   ! special case
                   tsno(c,k+1)    = tsno(c,k)
                else
                   ! use temperature gradient
                   dtdz           = (tsno(c,k-1) - tsno(c,k))/((dzsno(c,k-1)+2*dzsno(c,k))/2.0_r8)
                   tsno(c,k+1) = tsno(c,k) - dtdz*dzsno(c,k)/2.0_r8
                   if (tsno(c,k+1) >= tfrz) then
                      tsno(c,k+1)  = tsno(c,k)
                   else
                      tsno(c,k) = tsno(c,k) + dtdz*dzsno(c,k)/2.0_r8
                   endif
                end if

                mbc_phi(c,k)   = mbc_phi(c,k) / 2.0_r8
                mbc_phi(c,k+1) = mbc_phi(c,k)
                mbc_pho(c,k)   = mbc_pho(c,k) / 2.0_r8
                mbc_pho(c,k+1) = mbc_pho(c,k)
                moc_phi(c,k)   = moc_phi(c,k) / 2.0_r8
                moc_phi(c,k+1) = moc_phi(c,k)
                moc_pho(c,k)   = moc_pho(c,k) / 2.0_r8
                moc_pho(c,k+1) = moc_pho(c,k)
                mdst1(c,k)     = mdst1(c,k) / 2.0_r8
                mdst1(c,k+1)   = mdst1(c,k)
                mdst2(c,k)     = mdst2(c,k) / 2.0_r8
                mdst2(c,k+1)   = mdst2(c,k)
                mdst3(c,k)     = mdst3(c,k) / 2.0_r8
                mdst3(c,k+1)   = mdst3(c,k)
                mdst4(c,k)     = mdst4(c,k) / 2.0_r8
                mdst4(c,k+1)   = mdst4(c,k)

                rds(c,k+1)     = rds(c,k)
             end if
          end if

          ! There are layers below (note this is not exclusive with previous
          ! if-statement, since msno may have increased in the previous if-statement)
          if (k < msno) then

             if (is_lake) then
                offset = lsadz
             else
                offset = 0._r8
             end if

             if (dzsno(c,k) > dzmax_u(k) + offset ) then
                ! Only dump excess snow to underlying layer in a conservative fashion.
                ! Other quantities will depend on the height of the excess snow: a ratio is used for this.
                drr      = dzsno(c,k) - dzmax_u(k) - offset

                propor   = drr/dzsno(c,k)
                do wi = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
                   zwice(wi) = propor*swice(wi,c,k)
                   zwliq(wi) = propor*swliq(wi,c,k)
                end do
                zmbc_phi = propor*mbc_phi(c,k)
                zmbc_pho = propor*mbc_pho(c,k)
                zmoc_phi = propor*moc_phi(c,k)
                zmoc_pho = propor*moc_pho(c,k)
                zmdst1   = propor*mdst1(c,k)
                zmdst2   = propor*mdst2(c,k)
                zmdst3   = propor*mdst3(c,k)
                zmdst4   = propor*mdst4(c,k)

                propor         = (dzmax_u(k)+offset)/dzsno(c,k)
                do wi = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
                   swice(wi,c,k) = propor*swice(wi,c,k)
                   swliq(wi,c,k) = propor*swliq(wi,c,k)
                end do
                mbc_phi(c,k)   = propor*mbc_phi(c,k)
                mbc_pho(c,k)   = propor*mbc_pho(c,k)
                moc_phi(c,k)   = propor*moc_phi(c,k)
                moc_pho(c,k)   = propor*moc_pho(c,k)
                mdst1(c,k)     = propor*mdst1(c,k)
                mdst2(c,k)     = propor*mdst2(c,k)
                mdst3(c,k)     = propor*mdst3(c,k)
                mdst4(c,k)     = propor*mdst4(c,k)

                ! Set depth layer k to maximum allowed value
                dzsno(c,k)  = dzmax_u(k)  + offset

                mbc_phi(c,k+1) = mbc_phi(c,k+1)+zmbc_phi  ! (combo)
                mbc_pho(c,k+1) = mbc_pho(c,k+1)+zmbc_pho  ! (combo)
                moc_phi(c,k+1) = moc_phi(c,k+1)+zmoc_phi  ! (combo)
                moc_pho(c,k+1) = moc_pho(c,k+1)+zmoc_pho  ! (combo)
                mdst1(c,k+1)   = mdst1(c,k+1)+zmdst1  ! (combo)
                mdst2(c,k+1)   = mdst2(c,k+1)+zmdst2  ! (combo)
                mdst3(c,k+1)   = mdst3(c,k+1)+zmdst3  ! (combo)
                mdst4(c,k+1)   = mdst4(c,k+1)+zmdst4  ! (combo)

                ! Mass-weighted combination of radius
                rds(c,k+1) = MassWeightedSnowRadius( rds(c,k), rds(c,k+1), &
                     (swliq(i_bulk,c,k+1)+swice(i_bulk,c,k+1)), (zwliq(i_bulk)+zwice(i_bulk)) )

                call Combo (dzsno(c,k+1), swliq(i_bulk,c,k+1), swice(i_bulk,c,k+1), tsno(c,k+1), drr, &
                     zwliq(i_bulk), zwice(i_bulk), tsno(c,k))

                ! Bulk already combined in Combo; here we just need to loop over tracers
                ! and do a similar combination for them
                do wi = water_inst%tracers_beg, water_inst%tracers_end
                   swliq(wi,c,k+1) = swliq(wi,c,k+1) + zwliq(wi)
                   swice(wi,c,k+1) = swice(wi,c,k+1) + zwice(wi)
                end do
             end if
          end if
          k = k+1
       end do loop_layers

       snl(c) = -msno

    end do loop_snowcolumns

    do j = -nlevsno+1,0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             if (is_lake) then
                dz(c,j) = dzsno(c,j-snl(c))
             else
                dz(c,j) = dzsno(c,j-snl(c))/frac_sno(c)
             end if

             do wi = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
                associate(w => water_inst%bulk_and_tracers(wi))

                w%waterstate_inst%h2osoi_ice_col(c,j) = swice(wi,c,j-snl(c))
                w%waterstate_inst%h2osoi_liq_col(c,j) = swliq(wi,c,j-snl(c))

                end associate
             end do

             t_soisno(c,j)   = tsno(c,j-snl(c))
             mss_bcphi(c,j)   = mbc_phi(c,j-snl(c))
             mss_bcpho(c,j)   = mbc_pho(c,j-snl(c))
             mss_ocphi(c,j)   = moc_phi(c,j-snl(c))
             mss_ocpho(c,j)   = moc_pho(c,j-snl(c))
             mss_dst1(c,j)    = mdst1(c,j-snl(c))
             mss_dst2(c,j)    = mdst2(c,j-snl(c))
             mss_dst3(c,j)    = mdst3(c,j-snl(c))
             mss_dst4(c,j)    = mdst4(c,j-snl(c))
             snw_rds(c,j)     = rds(c,j-snl(c))

          end if
       end do
    end do

    ! Consistency check
    if (is_lake) then
       do j = -nlevsno + 1, 0
          do fc = 1, num_snowc
             c = filter_snowc(fc)

             if (j >= snl(c)+1) then
                dztot(c) = dztot(c) - dz(c,j)

                do wi = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
                   associate(w => water_inst%bulk_and_tracers(wi))

                   snwicetot(wi,c) = snwicetot(wi,c) - w%waterstate_inst%h2osoi_ice_col(c,j)
                   snwliqtot(wi,c) = snwliqtot(wi,c) - w%waterstate_inst%h2osoi_liq_col(c,j)

                   end associate
                end do
             end if

             if (j == 0) then
                if ( abs(dztot(c)) > 1.e-10_r8) then
                   write(iulog,*)'Inconsistency in SnowDivision_Lake! c, remainders', &
                        'dztot = ',c,dztot(c)
                   call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, msg=errmsg(sourcefile, __LINE__))
                end if

                do wi = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
                   if ( abs(snwicetot(wi,c)) > 1.e-7_r8 .or. abs(snwliqtot(wi,c)) > 1.e-7_r8 ) then
                      write(iulog,*)'Inconsistency in SnowDivision_Lake! wi, c, remainders', &
                           'snwicetot, snwliqtot = ',wi,c,snwicetot(wi,c),snwliqtot(wi,c)
                      call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, msg=errmsg(sourcefile, __LINE__))
                   end if
                end do
             end if
          end do
       end do
    end if

    ! Reset the node depth and the depth of layer interface

    do j = 0, -nlevsno+1, -1
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             z(c,j)    = zi(c,j) - 0.5_r8*dz(c,j)
             zi(c,j-1) = zi(c,j) - dz(c,j)
          end if
       end do
    end do

    end associate
    end associate
  end subroutine DivideSnowLayers

  !-----------------------------------------------------------------------
  subroutine ZeroEmptySnowLayers(bounds, num_snowc, filter_snowc, &
       col, water_inst, temperature_inst)
    !
    ! !DESCRIPTION:
    ! Set empty snow layers to zero
    !
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds
    integer                   , intent(in)    :: num_snowc       ! number of column snow points in column filter
    integer                   , intent(in)    :: filter_snowc(:) ! column filter for snow points
    type(column_type)         , intent(inout) :: col
    type(water_type)          , intent(inout) :: water_inst
    type(temperature_type)    , intent(inout) :: temperature_inst
    !
    ! !LOCAL VARIABLES:
    integer :: j
    integer :: fc, c
    integer :: wi                                       ! index of water tracer or bulk

    character(len=*), parameter :: subname = 'ZeroEmptySnowLayers'
    !-----------------------------------------------------------------------

    ! In contrast to most routines, this one operates on a mix of bulk-only quantities and
    ! bulk-and-tracer quantities. Where bulk-and-tracer quantities are referenced, they
    ! are referred to like w%waterstate_inst%h2osoi_liq_col.

    associate( &
         snl                => col%snl                                , & ! Input:  [integer  (:)   ]  number of snow layers
         z                  => col%z                                  , & ! Output: [real(r8) (:,:) ]  layer depth  (m)
         dz                 => col%dz                                 , & ! Output: [real(r8) (:,:) ]  layer thickness depth (m)
         zi                 => col%zi                                 , & ! Output: [real(r8) (:,:) ]  interface depth (m)
         t_soisno           => temperature_inst%t_soisno_col            & ! Output: [real(r8) (:,:) ]  soil temperature (Kelvin)
         )

    do j = -nlevsno+1,0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j <= snl(c) .and. snl(c) > -nlevsno) then
             do wi = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
                associate(w => water_inst%bulk_and_tracers(wi))
                w%waterstate_inst%h2osoi_ice_col(c,j) = 0._r8
                w%waterstate_inst%h2osoi_liq_col(c,j) = 0._r8
                end associate
             end do
             t_soisno(c,j)  = 0._r8
             dz(c,j)    = 0._r8
             z(c,j)     = 0._r8
             zi(c,j-1)  = 0._r8
          end if
       end do
    end do

    end associate

  end subroutine ZeroEmptySnowLayers

  !-----------------------------------------------------------------------
  subroutine InitSnowLayers (bounds, snow_depth)
    !
    ! !DESCRIPTION:
    ! Initialize snow layer depth from specified total depth.
    !
    use spmdMod, only : masterproc
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    real(r8)               , intent(in)    :: snow_depth(bounds%begc:)
    !
    !
    ! LOCAL VARAIBLES:
    integer               :: c,l,j              ! indices
    real(r8)              :: minbound, maxbound ! helper variables
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(snow_depth)  == (/bounds%endc/)), sourcefile, __LINE__)

    associate( &
         snl => col%snl,   & ! Output: [integer (:)    ]  number of snow layers
         dz  => col%dz,    & ! Output: [real(r8) (:,:) ]  layer thickness (m)  (-nlevsno+1:nlevmaxurbgrnd)
         z   => col%z,     & ! Output: [real(r8) (:,:) ]  layer depth (m) (-nlevsno+1:nlevmaxurbgrnd)
         zi  => col%zi     & ! Output: [real(r8) (:,:) ]  interface level below a "z" level (m) (-nlevsno+0:nlevmaxurbgrnd)
    )

    allocate(dzmin(1:nlevsno))
    allocate(dzmax_l(1:nlevsno))
    allocate(dzmax_u(1:nlevsno))

    ! These three variables determine the vertical structure of the snow pack:
    ! dzmin: minimum snow thickness of layer
    ! dzmax_l: maximum snow thickness of layer when no layers beneath
    ! dzmax_u: maximum snow thickness of layer when layers beneath
    dzmin(1) = snow_dzmin_1  ! default or user-defined value from namelist
    dzmax_l(1) = snow_dzmax_l_1  ! same comment
    dzmax_u(1) = snow_dzmax_u_1  ! same comment
    dzmin(2) = snow_dzmin_2  ! default or user-defined value from namelist
    dzmax_l(2) = snow_dzmax_l_2  ! same comment
    dzmax_u(2) = snow_dzmax_u_2  ! same comment
    do j = 3, nlevsno
       dzmin(j) = dzmax_u(j-1) * 0.5_r8
       dzmax_u(j) = 2._r8 * dzmax_u(j-1) + 0.01_r8
       dzmax_l(j) = dzmax_u(j) + dzmax_l(j-1)
       if (j == nlevsno) then
          dzmax_u(j) = huge(1._r8)
          dzmax_l(j) = huge(1._r8)
       end if
    end do

    ! Error check loops
    do j = 2, nlevsno
       if (dzmin(j) <= dzmin(j-1)) then
          write(iulog,*) 'ERROR at snow layer j =', j, ' because dzmin(j) =', dzmin(j), ' and dzmin(j-1) =', dzmin(j-1)
          call endrun(msg="ERROR dzmin(j) cannot be <= dzmin(j-1)"// &
               errMsg(sourcefile, __LINE__))
       end if
       if (dzmax_u(j) <= dzmax_u(j-1)) then
          write(iulog,*) 'ERROR at snow layer j =', j, ' because dzmax_u(j) =', dzmax_u(j), ' and dzmax_u(j-1) =', dzmax_u(j-1)
          call endrun(msg="ERROR dzmax_u(j) cannot be <= dzmax_u(j-1)"// &
               errMsg(sourcefile, __LINE__))
       end if
       if (dzmax_l(j) <= dzmax_l(j-1)) then
          write(iulog,*) 'ERROR at snow layer j =', j, ' because dzmax_l(j) =', dzmax_l(j), ' and dzmax_l(j-1) =', dzmax_l(j-1)
          call endrun(msg="ERROR dzmax_l(j) cannot be <= dzmax_l(j-1)"// &
               errMsg(sourcefile, __LINE__))
       end if
    end do
    do j = 1, nlevsno
       if (dzmin(j) >= dzmax_u(j)) then
          write(iulog,*) 'ERROR at snow layer j =', j, ' because dzmin(j) =', dzmin(j), ' and dzmax_u(j) =', dzmax_u(j)
          call endrun(msg="ERROR dzmin(j) cannot be >= dzmax_u(j)"// &
               errMsg(sourcefile, __LINE__))
       end if
    end do
    do j = 1, nlevsno-1
       if (dzmax_u(j) >= dzmax_l(j)) then
          write(iulog,*) 'ERROR at snow layer j =', j, ' because dzmax_u(j) =', dzmax_u(j), ' and dzmax_l(j) =', dzmax_l(j)
          call endrun(msg="ERROR dzmax_u(j) cannot be >= dzmax_l(j)"// &
               errMsg(sourcefile, __LINE__))
       end if
    end do

    if (masterproc) then
       write(iulog,*) 'dzmin =', dzmin
       write(iulog,*) 'dzmax_l =', dzmax_l
       write(iulog,*) 'dzmax_u =', dzmax_u
    end if

    loop_columns: do c = bounds%begc,bounds%endc
       l = col%landunit(c)

       dz(c,-nlevsno+1: 0) = spval
       z (c,-nlevsno+1: 0) = spval
       zi(c,-nlevsno  :-1) = spval

       ! Special case: lake
       if (lun%lakpoi(l)) then
          snl(c)              = 0
          dz(c,-nlevsno+1:0)  = 0._r8
          z(c,-nlevsno+1:0)   = 0._r8
          zi(c,-nlevsno+0:0)  = 0._r8
          cycle
       end if

       ! LvK 9-JUN-2015: in CanopyHydrologyMod , snow_depth is scaled with frac_sno
       ! Here we do not apply scaling to snow_depth, so inconsistent? TODO

       ! Special case: too little snow for snowpack existence
       if (snow_depth(c) < dzmin(1)) then
          snl(c)              = 0
          dz(c,-nlevsno+1:0)  = 0._r8
          z(c,-nlevsno+1:0)   = 0._r8
          zi(c,-nlevsno+0:0)  = 0._r8
          cycle
       end if

       ! There has to be at least one snow layer
       snl(c)   = -1
       minbound = dzmin(1)
       maxbound = dzmax_l(1)

       if (snow_depth(c) >= minbound .and. snow_depth(c) <= maxbound) then
          ! Special case: single layer
          dz(c,0) = snow_depth(c)

       else
          ! Search for appropriate number of layers (snl) by increasing the number
          ! the number of layers and check for matching bounds.
          snl(c) = snl(c) - 1
          minbound = maxbound
          maxbound = sum(dzmax_u(1:-snl(c)))

          do while(snow_depth(c) > maxbound .and. -snl(c) < nlevsno )
             snl(c) = snl(c) - 1
             minbound = maxbound
             maxbound = sum(dzmax_u(1:-snl(c)))
          end do

          ! Set thickness of all layers except bottom two
          do j = 1, -snl(c)-2
             dz(c,j+snl(c))  = dzmax_u(j)
          enddo

          ! Determine whether the two bottom layers should be equal in size,
          ! or not. The rule here is: always create equal size when possible.
          if (snow_depth(c) <= sum(dzmax_u(1:-snl(c)-2)) + 2 * dzmax_u(-snl(c)-1)) then
             dz(c,-1) = (snow_depth(c) - sum(dzmax_u(1:-snl(c)-2))) / 2._r8
             dz(c,0)  = dz(c,-1)
          else
             dz(c,-1) = dzmax_u(-snl(c)-1)
             dz(c,0)  = snow_depth(c) - sum(dzmax_u(1:-snl(c)-1))
          endif
       endif

       ! Initialize the node depth and the depth of layer interface
       do j = 0, snl(c)+1, -1
          z(c,j)    = zi(c,j) - 0.5_r8*dz(c,j)
          zi(c,j-1) = zi(c,j) - dz(c,j)
       end do

    end do loop_columns

    end associate
  end subroutine InitSnowLayers

  !-----------------------------------------------------------------------
  subroutine SnowCapping(bounds, num_initc, filter_initc, num_snowc, filter_snowc, &
       topo_inst, aerosol_inst, water_inst )
    !
    ! !DESCRIPTION:
    ! Removes mass from bottom snow layer for columns that exceed the maximum snow depth.
    ! This routine is called twice: once for non-lake columns and once for lake columns. 
    ! The initialization of the snow capping fluxes should only be done ONCE for each group,
    ! therefore they are a passed as an extra argument (filter_initc). 
    ! Density and temperature of the layer are conserved (density needs some work, temperature is a state
    ! variable)
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_initc       ! number of column points that need to be initialized
    integer                , intent(in)    :: filter_initc(:) ! column filter for points that need to be initialized
    integer                , intent(in)    :: num_snowc       ! number of column snow points in column filter
    integer                , intent(in)    :: filter_snowc(:) ! column filter for snow points
    class(topo_type)       , intent(in)    :: topo_inst
    type(aerosol_type)     , intent(inout) :: aerosol_inst
    type(water_type)       , intent(inout) :: water_inst
    !
    ! !LOCAL VARIABLES:
    integer    :: i                                ! index of water tracer or bulk
    real(r8)   :: dtime                            ! land model time step (sec)
    real(r8)   :: h2osno_total(bounds%begc:bounds%endc)  ! total snow water (mm H2O)
    real(r8)   :: rho_orig_bottom(bounds%begc:bounds%endc) ! partial density of ice in bottom snow layer, before updates (not scaled with frac_sno) [kg/m3]
    real(r8)   :: frac_adjust(bounds%begc:bounds%endc) ! fraction of mass remaining after capping
    type(filter_col_type) :: snow_capping_filterc ! column filter: columns undergoing snow capping

    !-----------------------------------------------------------------------

    associate( &
         begc => bounds%begc, &
         endc => bounds%endc, &

         b_waterflux_inst  => water_inst%waterfluxbulk_inst, &
         b_waterstate_inst => water_inst%waterstatebulk_inst &
         )

    ! Determine model time step
    dtime = get_step_size_real()

    do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
       associate(w => water_inst%bulk_and_tracers(i))
       call InitFlux_SnowCapping(bounds, num_initc, filter_initc, &
            ! Outputs
            qflx_snwcp_ice           = w%waterflux_inst%qflx_snwcp_ice_col(begc:endc), &
            qflx_snwcp_liq           = w%waterflux_inst%qflx_snwcp_liq_col(begc:endc), &
            qflx_snwcp_discarded_ice = w%waterflux_inst%qflx_snwcp_discarded_ice_col(begc:endc), &
            qflx_snwcp_discarded_liq = w%waterflux_inst%qflx_snwcp_discarded_liq_col(begc:endc))
       end associate
    end do

    call b_waterstate_inst%CalculateTotalH2osno(bounds, num_snowc, filter_snowc, &
         caller = 'SnowCapping', &
         h2osno_total = h2osno_total(begc:endc))

    call BulkFlux_SnowCappingFluxes(bounds, num_snowc, filter_snowc, &
         ! Inputs
         dtime                    = dtime, &
         dz_bottom                = col%dz(begc:endc, 0), &
         topo                     = topo_inst%topo_col(begc:endc), &
         h2osno_total             = h2osno_total(begc:endc), &
         h2osoi_ice_bottom        = b_waterstate_inst%h2osoi_ice_col(begc:endc, 0), &
         h2osoi_liq_bottom        = b_waterstate_inst%h2osoi_liq_col(begc:endc, 0), &
         ! Outputs
         snow_capping_filterc     = snow_capping_filterc, &
         rho_orig_bottom          = rho_orig_bottom(begc:endc), &
         frac_adjust              = frac_adjust(begc:endc), &
         qflx_snwcp_ice           = b_waterflux_inst%qflx_snwcp_ice_col(begc:endc), &
         qflx_snwcp_liq           = b_waterflux_inst%qflx_snwcp_liq_col(begc:endc), &
         qflx_snwcp_discarded_ice = b_waterflux_inst%qflx_snwcp_discarded_ice_col(begc:endc), &
         qflx_snwcp_discarded_liq = b_waterflux_inst%qflx_snwcp_discarded_liq_col(begc:endc))

    do i = water_inst%tracers_beg, water_inst%tracers_end
       associate(w => water_inst%bulk_and_tracers(i))
       call TracerFlux_SnowCappingFluxes(bounds, snow_capping_filterc, &
            ! Inputs
            bulk_h2osoi_ice_bottom        = b_waterstate_inst%h2osoi_ice_col(begc:endc, 0), &
            bulk_h2osoi_liq_bottom        = b_waterstate_inst%h2osoi_liq_col(begc:endc, 0), &
            bulk_qflx_snwcp_ice           = b_waterflux_inst%qflx_snwcp_ice_col(begc:endc), &
            bulk_qflx_snwcp_liq           = b_waterflux_inst%qflx_snwcp_liq_col(begc:endc), &
            bulk_qflx_snwcp_discarded_ice = b_waterflux_inst%qflx_snwcp_discarded_ice_col(begc:endc), &
            bulk_qflx_snwcp_discarded_liq = b_waterflux_inst%qflx_snwcp_discarded_liq_col(begc:endc), &
            trac_h2osoi_ice_bottom        = w%waterstate_inst%h2osoi_ice_col(begc:endc, 0), &
            trac_h2osoi_liq_bottom        = w%waterstate_inst%h2osoi_liq_col(begc:endc, 0), &
            ! Outputs
            trac_qflx_snwcp_ice           = w%waterflux_inst%qflx_snwcp_ice_col(begc:endc), &
            trac_qflx_snwcp_liq           = w%waterflux_inst%qflx_snwcp_liq_col(begc:endc), &
            trac_qflx_snwcp_discarded_ice = w%waterflux_inst%qflx_snwcp_discarded_ice_col(begc:endc), &
            trac_qflx_snwcp_discarded_liq = w%waterflux_inst%qflx_snwcp_discarded_liq_col(begc:endc))
       end associate
    end do

    do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
       associate(w => water_inst%bulk_and_tracers(i))
       call UpdateState_RemoveSnowCappingFluxes(bounds, snow_capping_filterc, &
            ! Inputs
            dtime                    = dtime, &
            qflx_snwcp_ice           = w%waterflux_inst%qflx_snwcp_ice_col(begc:endc), &
            qflx_snwcp_liq           = w%waterflux_inst%qflx_snwcp_liq_col(begc:endc), &
            qflx_snwcp_discarded_ice = w%waterflux_inst%qflx_snwcp_discarded_ice_col(begc:endc), &
            qflx_snwcp_discarded_liq = w%waterflux_inst%qflx_snwcp_discarded_liq_col(begc:endc), &
            ! Outputs
            h2osoi_ice_bottom        = w%waterstate_inst%h2osoi_ice_col(begc:endc, 0), &
            h2osoi_liq_bottom        = w%waterstate_inst%h2osoi_liq_col(begc:endc, 0))
       end associate
    end do

    call SnowCappingUpdateDzAndAerosols(bounds, snow_capping_filterc, &
         ! Inputs
         rho_orig_bottom   = rho_orig_bottom(begc:endc), &
         h2osoi_ice_bottom = b_waterstate_inst%h2osoi_ice_col(begc:endc, 0), &
         frac_adjust       = frac_adjust(begc:endc), &
         ! Outputs
         dz_bottom         = col%dz(begc:endc, 0), &
         mss_bcphi_bottom  = aerosol_inst%mss_bcphi_col(begc:endc, 0), &
         mss_bcpho_bottom  = aerosol_inst%mss_bcpho_col(begc:endc, 0), &
         mss_ocphi_bottom  = aerosol_inst%mss_ocphi_col(begc:endc, 0), &
         mss_ocpho_bottom  = aerosol_inst%mss_ocpho_col(begc:endc, 0), &
         mss_dst1_bottom   = aerosol_inst%mss_dst1_col(begc:endc, 0), &
         mss_dst2_bottom   = aerosol_inst%mss_dst2_col(begc:endc, 0), &
         mss_dst3_bottom   = aerosol_inst%mss_dst3_col(begc:endc, 0), &
         mss_dst4_bottom   = aerosol_inst%mss_dst4_col(begc:endc, 0))

    end associate
  end subroutine SnowCapping

  !-----------------------------------------------------------------------
  subroutine InitFlux_SnowCapping(bounds, num_initc, filter_initc, &
       qflx_snwcp_ice, qflx_snwcp_liq, qflx_snwcp_discarded_ice, qflx_snwcp_discarded_liq)
    !
    ! !DESCRIPTION:
    ! Initialize snow capping fluxes to 0
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_initc       ! number of column points that need to be initialized
    integer                , intent(in)    :: filter_initc(:) ! column filter for points that need to be initialized

    real(r8), intent(inout) :: qflx_snwcp_ice( bounds%begc: ) ! excess solid h2o due to snow capping (outgoing) (mm H2O /s)
    real(r8), intent(inout) :: qflx_snwcp_liq( bounds%begc: ) ! excess liquid h2o due to snow capping (outgoing) (mm H2O /s)
    real(r8), intent(inout) :: qflx_snwcp_discarded_ice( bounds%begc: ) ! excess solid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s)
    real(r8), intent(inout) :: qflx_snwcp_discarded_liq( bounds%begc: ) ! excess liquid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c

    character(len=*), parameter :: subname = 'InitFlux_SnowCapping'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(qflx_snwcp_ice, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_snwcp_liq, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_snwcp_discarded_ice, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_snwcp_discarded_liq, 1) == bounds%endc), sourcefile, __LINE__)

    do fc = 1, num_initc
       c = filter_initc(fc)
       qflx_snwcp_ice(c) = 0.0_r8
       qflx_snwcp_liq(c) = 0.0_r8
       qflx_snwcp_discarded_ice(c) = 0.0_r8
       qflx_snwcp_discarded_liq(c) = 0.0_r8
    end do

  end subroutine InitFlux_SnowCapping

  !-----------------------------------------------------------------------
  subroutine BulkFlux_SnowCappingFluxes(bounds, num_snowc, filter_snowc, &
       dtime, dz_bottom, topo, h2osno_total, h2osoi_ice_bottom, h2osoi_liq_bottom, &
       snow_capping_filterc, rho_orig_bottom, frac_adjust, &
       qflx_snwcp_ice, qflx_snwcp_liq, qflx_snwcp_discarded_ice, qflx_snwcp_discarded_liq)
    !
    ! !DESCRIPTION:
    ! Calculate snow capping fluxes and related terms for bulk water
    !
    ! The output arrays are set within the points given by snow_capping_filterc (which is
    ! a subset of the snowc filter); elsewhere, they are left at their original values
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds
    integer           , intent(in) :: num_snowc       ! number of column snow points in column filter
    integer           , intent(in) :: filter_snowc(:) ! column filter for snow points

    real(r8) , intent(in) :: dtime                             ! land model time step (sec)
    real(r8) , intent(in) :: dz_bottom( bounds%begc: )         ! layer depth of bottom snow layer (m)
    real(r8) , intent(in) :: topo( bounds%begc: )              ! column surface height (m)
    real(r8) , intent(in) :: h2osno_total( bounds%begc: )      ! total snow water (mm H2O)
    real(r8) , intent(in) :: h2osoi_ice_bottom( bounds%begc: ) ! ice lens in bottom snow layer (kg/m2)
    real(r8) , intent(in) :: h2osoi_liq_bottom( bounds%begc: ) ! liquid water in bottom snow layer (kg/m2)

    type(filter_col_type) , intent(out)   :: snow_capping_filterc                     ! column filter: columns undergoing snow capping
    real(r8)              , intent(inout) :: rho_orig_bottom( bounds%begc: )          ! partial density of ice in bottom snow layer, before updates (not scaled with frac_sno) (kg/m3)
    real(r8)              , intent(inout) :: frac_adjust( bounds%begc: )              ! fraction of mass remaining after capping
    real(r8)              , intent(inout) :: qflx_snwcp_ice( bounds%begc: )           ! excess solid h2o due to snow capping (outgoing) (mm H2O /s)
    real(r8)              , intent(inout) :: qflx_snwcp_liq( bounds%begc: )           ! excess liquid h2o due to snow capping (outgoing) (mm H2O /s)
    real(r8)              , intent(inout) :: qflx_snwcp_discarded_ice( bounds%begc: ) ! excess solid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s)
    real(r8)              , intent(inout) :: qflx_snwcp_discarded_liq( bounds%begc: ) ! excess liquid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s)
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c
    real(r8) :: mss_snwcp_tot                               ! total snow capping mass [kg/m2] 
    real(r8) :: mss_snow_bottom_lyr                         ! total snow mass (ice+liquid) in bottom layer [kg/m2]
    real(r8) :: snwcp_flux_ice                              ! snow capping flux (ice) [kg/m2]
    real(r8) :: snwcp_flux_liq                              ! snow capping flux (liquid) [kg/m2]
    real(r8) :: icefrac                                     ! fraction of ice mass w.r.t. total mass [unitless]
    real(r8) :: h2osno_excess(bounds%begc:bounds%endc)      ! excess snow that needs to be capped [mm H2O]
    logical  :: apply_runoff(bounds%begc:bounds%endc)       ! for columns with capping, whether the capping flux should be sent to runoff
    logical  :: column_has_capping(bounds%begc:bounds%endc) ! true for columns with snow capping

    ! Always keep at least this fraction of the bottom snow layer when doing snow capping
    ! This needs to be slightly greater than 0 to avoid roundoff problems
    real(r8), parameter :: min_snow_to_keep = 1.e-3_r8  ! fraction of bottom snow layer to keep with capping

    character(len=*), parameter :: subname = 'BulkFlux_SnowCappingFluxes'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(dz_bottom, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(topo, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osno_total, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osoi_ice_bottom, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osoi_liq_bottom, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(rho_orig_bottom, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(frac_adjust, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_snwcp_ice, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_snwcp_liq, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_snwcp_discarded_ice, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_snwcp_discarded_liq, 1) == bounds%endc), sourcefile, __LINE__)

    call SnowCappingExcess(bounds, num_snowc, filter_snowc, &
         h2osno = h2osno_total(bounds%begc:bounds%endc), &
         topo = topo(bounds%begc:bounds%endc), &
         h2osno_excess = h2osno_excess(bounds%begc:bounds%endc), &
         apply_runoff = apply_runoff(bounds%begc:bounds%endc))

    do fc = 1, num_snowc
       c = filter_snowc(fc)
       if (h2osno_excess(c) > 0._r8) then
          column_has_capping(c) = .true.
       else
          column_has_capping(c) = .false.
       end if
    end do

    snow_capping_filterc = col_filter_from_filter_and_logical_array( &
         bounds = bounds, &
         num_orig = num_snowc, &
         filter_orig = filter_snowc, &
         logical_col = column_has_capping(bounds%begc:bounds%endc))

    do fc = 1, snow_capping_filterc%num
       c = snow_capping_filterc%indices(fc)

       rho_orig_bottom(c) = h2osoi_ice_bottom(c) / dz_bottom(c) ! ice only

       mss_snow_bottom_lyr = h2osoi_ice_bottom(c) + h2osoi_liq_bottom(c) 
       mss_snwcp_tot = min(h2osno_excess(c), mss_snow_bottom_lyr * (1._r8 - min_snow_to_keep)) ! Can't remove more mass than available

       ! Ratio of snow/liquid in bottom layer determines partitioning of runoff fluxes
       icefrac = h2osoi_ice_bottom(c) / mss_snow_bottom_lyr
       snwcp_flux_ice = mss_snwcp_tot/dtime * icefrac
       snwcp_flux_liq = mss_snwcp_tot/dtime * (1._r8 - icefrac)
       if (apply_runoff(c)) then
          qflx_snwcp_ice(c) = snwcp_flux_ice
          qflx_snwcp_liq(c) = snwcp_flux_liq
       else
          qflx_snwcp_discarded_ice(c) = snwcp_flux_ice
          qflx_snwcp_discarded_liq(c) = snwcp_flux_liq
       end if

       frac_adjust(c) = (mss_snow_bottom_lyr - mss_snwcp_tot) / mss_snow_bottom_lyr

    end do

  end subroutine BulkFlux_SnowCappingFluxes

  !-----------------------------------------------------------------------
  subroutine SnowCappingExcess(bounds, num_snowc, filter_snowc, &
       h2osno, topo, h2osno_excess, apply_runoff)
    !
    ! !DESCRIPTION:
    ! Determine the amount of excess snow that needs to be capped
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer  , intent(in)  :: num_snowc                     ! number of column snow points in column filter
    integer  , intent(in)  :: filter_snowc(:)               ! column filter for snow points
    real(r8) , intent(in)  :: h2osno( bounds%begc: )        ! snow water (mm H2O)
    real(r8) , intent(in)  :: topo( bounds%begc: )          ! column surface height (m)
    real(r8) , intent(out) :: h2osno_excess( bounds%begc: ) ! excess snow that needs to be capped (mm H2O)
    logical  , intent(out) :: apply_runoff( bounds%begc: )  ! whether capped snow should be sent to runoff; only valid where h2osno_excess > 0
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c, l
    integer :: reset_snow_timesteps
    logical :: is_reset_snow_active  ! whether snow resetting is active in this time step for at least some points

    character(len=*), parameter :: subname = 'SnowCappingExcess'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(h2osno) == (/bounds%endc/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(topo) == (/bounds%endc/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(h2osno_excess) == (/bounds%endc/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(apply_runoff) == (/bounds%endc/)), sourcefile, __LINE__)

    do fc = 1, num_snowc
       c = filter_snowc(fc)
       h2osno_excess(c) = 0._r8
       if (h2osno(c) > h2osno_max) then
          h2osno_excess(c) = h2osno(c) - h2osno_max
          apply_runoff(c) = .true.
       end if
    end do

    ! Implement snow resetting (i.e., resetting points that have h2osno greater than some
    ! value) by applying the snow capping scheme to a value that's smaller than
    ! h2osno_max, but NOT sending the resulting capping flux to the coupler. It is easier
    ! to implement the resetting this way than to try to manually reset the snow pack,
    ! because there are so many snow pack variables that need to be kept consistent if
    ! doing this resetting manually. Note that we need to continue to apply the resetting
    ! for some number of time steps, because we can remove at most one snow layer per
    ! time step.

    ! It is important that this snow resetting comes after the standard check (h2osno(c) >
    ! h2osno_max), so that we override any standard capping.
    is_reset_snow_active = .false.
    if (reset_snow .or. reset_snow_glc) then
       reset_snow_timesteps = reset_snow_timesteps_per_layer * nlevsno
       if (get_nstep() <= reset_snow_timesteps) then
          is_reset_snow_active = .true.
       end if
    end if

    if (is_reset_snow_active) then
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          l = col%landunit(c)
          if ((lun%itype(l) /= istice) .and. &
               reset_snow .and. &
               (h2osno(c) > reset_snow_h2osno)) then
             h2osno_excess(c) = h2osno(c) - reset_snow_h2osno
             apply_runoff(c) = .false.
          else if ((lun%itype(l) == istice) .and. &
               reset_snow_glc .and. &
               (h2osno(c) > reset_snow_h2osno) .and. &
               (topo(c) <= reset_snow_glc_ela)) then
             h2osno_excess(c) = h2osno(c) - reset_snow_h2osno
             apply_runoff(c) = .false.
          end if
       end do
    end if

  end subroutine SnowCappingExcess

  !-----------------------------------------------------------------------
  subroutine TracerFlux_SnowCappingFluxes(bounds, snow_capping_filterc, &
       bulk_h2osoi_ice_bottom, bulk_h2osoi_liq_bottom, &
       bulk_qflx_snwcp_ice, bulk_qflx_snwcp_liq, bulk_qflx_snwcp_discarded_ice, bulk_qflx_snwcp_discarded_liq, &
       trac_h2osoi_ice_bottom, trac_h2osoi_liq_bottom, &
       trac_qflx_snwcp_ice, trac_qflx_snwcp_liq, trac_qflx_snwcp_discarded_ice, trac_qflx_snwcp_discarded_liq)
    !
    ! !DESCRIPTION:
    ! Calculate snow capping fluxes and related terms for one tracer
    !
    ! !ARGUMENTS:
    !
    type(bounds_type)     , intent(in) :: bounds
    type(filter_col_type) , intent(in) :: snow_capping_filterc ! column filter: columns undergoing snow capping
    
    ! For description of arguments, see comments in BulkFlux_SnowCappingFluxes. Here,
    ! bulk_* variables refer to bulk water and trac_* variables refer to the given water
    ! tracer.
    real(r8), intent(in) :: bulk_h2osoi_ice_bottom( bounds%begc: )
    real(r8), intent(in) :: bulk_h2osoi_liq_bottom( bounds%begc: )
    real(r8), intent(in) :: bulk_qflx_snwcp_ice( bounds%begc: )
    real(r8), intent(in) :: bulk_qflx_snwcp_liq( bounds%begc: )
    real(r8), intent(in) :: bulk_qflx_snwcp_discarded_ice( bounds%begc: )
    real(r8), intent(in) :: bulk_qflx_snwcp_discarded_liq( bounds%begc: )
    real(r8), intent(in) :: trac_h2osoi_ice_bottom( bounds%begc: )
    real(r8), intent(in) :: trac_h2osoi_liq_bottom( bounds%begc: )
    real(r8), intent(inout) :: trac_qflx_snwcp_ice( bounds%begc: )
    real(r8), intent(inout) :: trac_qflx_snwcp_liq( bounds%begc: )
    real(r8), intent(inout) :: trac_qflx_snwcp_discarded_ice( bounds%begc: )
    real(r8), intent(inout) :: trac_qflx_snwcp_discarded_liq( bounds%begc: )

    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'TracerFlux_SnowCappingFluxes'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(bulk_h2osoi_ice_bottom, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(bulk_h2osoi_liq_bottom, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(bulk_qflx_snwcp_ice, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(bulk_qflx_snwcp_liq, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(bulk_qflx_snwcp_discarded_ice, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(bulk_qflx_snwcp_discarded_liq, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(trac_h2osoi_ice_bottom, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(trac_h2osoi_liq_bottom, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(trac_qflx_snwcp_ice, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(trac_qflx_snwcp_liq, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(trac_qflx_snwcp_discarded_ice, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(trac_qflx_snwcp_discarded_liq, 1) == bounds%endc), sourcefile, __LINE__)

    associate( &
         begc => bounds%begc, &
         endc => bounds%endc  &
         )

    call CalcTracerFromBulk( &
         subgrid_level = subgrid_level_column, &
         lb            = begc, &
         num_pts       = snow_capping_filterc%num, &
         filter_pts    = snow_capping_filterc%indices, &
         bulk_source   = bulk_h2osoi_ice_bottom(begc:endc), &
         bulk_val      = bulk_qflx_snwcp_ice(begc:endc), &
         tracer_source = trac_h2osoi_ice_bottom(begc:endc), &
         tracer_val    = trac_qflx_snwcp_ice(begc:endc))

    call CalcTracerFromBulk( &
         subgrid_level = subgrid_level_column, &
         lb            = begc, &
         num_pts       = snow_capping_filterc%num, &
         filter_pts    = snow_capping_filterc%indices, &
         bulk_source   = bulk_h2osoi_liq_bottom(begc:endc), &
         bulk_val      = bulk_qflx_snwcp_liq(begc:endc), &
         tracer_source = trac_h2osoi_liq_bottom(begc:endc), &
         tracer_val    = trac_qflx_snwcp_liq(begc:endc))

    call CalcTracerFromBulk( &
         subgrid_level = subgrid_level_column, &
         lb            = begc, &
         num_pts       = snow_capping_filterc%num, &
         filter_pts    = snow_capping_filterc%indices, &
         bulk_source   = bulk_h2osoi_ice_bottom(begc:endc), &
         bulk_val      = bulk_qflx_snwcp_discarded_ice(begc:endc), &
         tracer_source = trac_h2osoi_ice_bottom(begc:endc), &
         tracer_val    = trac_qflx_snwcp_discarded_ice(begc:endc))

    call CalcTracerFromBulk( &
         subgrid_level = subgrid_level_column, &
         lb            = begc, &
         num_pts       = snow_capping_filterc%num, &
         filter_pts    = snow_capping_filterc%indices, &
         bulk_source   = bulk_h2osoi_liq_bottom(begc:endc), &
         bulk_val      = bulk_qflx_snwcp_discarded_liq(begc:endc), &
         tracer_source = trac_h2osoi_liq_bottom(begc:endc), &
         tracer_val    = trac_qflx_snwcp_discarded_liq(begc:endc))

    end associate

  end subroutine TracerFlux_SnowCappingFluxes

  !-----------------------------------------------------------------------
  subroutine UpdateState_RemoveSnowCappingFluxes(bounds, snow_capping_filterc, &
       dtime, qflx_snwcp_ice, qflx_snwcp_liq, qflx_snwcp_discarded_ice, qflx_snwcp_discarded_liq, &
       h2osoi_ice_bottom, h2osoi_liq_bottom)
    !
    ! !DESCRIPTION:
    ! Remove snow capping fluxes from h2osoi_ice and h2osoi_liq
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in) :: bounds
    type(filter_col_type) , intent(in) :: snow_capping_filterc ! column filter: columns undergoing snow capping

    real(r8) , intent(in)    :: dtime                                    ! land model time step (sec)
    real(r8) , intent(in)    :: qflx_snwcp_ice( bounds%begc: )           ! excess solid h2o due to snow capping (outgoing) (mm H2O /s)
    real(r8) , intent(in)    :: qflx_snwcp_liq( bounds%begc: )           ! excess liquid h2o due to snow capping (outgoing) (mm H2O /s)
    real(r8) , intent(in)    :: qflx_snwcp_discarded_ice( bounds%begc: ) ! excess solid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s)
    real(r8) , intent(in)    :: qflx_snwcp_discarded_liq( bounds%begc: ) ! excess liquid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s)
    real(r8) , intent(inout) :: h2osoi_ice_bottom( bounds%begc: )        ! ice lens in bottom snow layer (kg/m2)
    real(r8) , intent(inout) :: h2osoi_liq_bottom( bounds%begc: )        ! liquid water in bottom snow layer (kg/m2)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c

    character(len=*), parameter :: subname = 'UpdateState_RemoveSnowCappingFluxes'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(qflx_snwcp_ice, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_snwcp_liq, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_snwcp_discarded_ice, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(qflx_snwcp_discarded_liq, 1) == bounds%endc), sourcefile, __LINE__)

    do fc = 1, snow_capping_filterc%num
       c = snow_capping_filterc%indices(fc)

       h2osoi_ice_bottom(c) = h2osoi_ice_bottom(c) - (qflx_snwcp_ice(c) + qflx_snwcp_discarded_ice(c))*dtime
       h2osoi_liq_bottom(c) = h2osoi_liq_bottom(c) - (qflx_snwcp_liq(c) + qflx_snwcp_discarded_liq(c))*dtime

       ! Check that water capacity is still positive
       if (h2osoi_ice_bottom(c) < 0._r8 .or. h2osoi_liq_bottom(c) < 0._r8 ) then
          write(iulog,*)'ERROR: capping procedure failed (negative mass remaining) c = ',c
          write(iulog,*)'h2osoi_ice_bottom = ', h2osoi_ice_bottom(c), ' h2osoi_liq_bottom = ', h2osoi_liq_bottom(c)
          call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, msg=errmsg(sourcefile, __LINE__))
       end if

    end do

  end subroutine UpdateState_RemoveSnowCappingFluxes

  !-----------------------------------------------------------------------
  subroutine SnowCappingUpdateDzAndAerosols(bounds, snow_capping_filterc, &
       rho_orig_bottom, h2osoi_ice_bottom, frac_adjust, dz_bottom, &
       mss_bcphi_bottom, mss_bcpho_bottom, mss_ocphi_bottom, mss_ocpho_bottom, &
       mss_dst1_bottom, mss_dst2_bottom, mss_dst3_bottom, mss_dst4_bottom)
    !
    ! !DESCRIPTION:
    ! Following snow capping, adjust dz and aerosol masses in bottom snow layer
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in) :: bounds
    type(filter_col_type) , intent(in) :: snow_capping_filterc ! column filter: columns undergoing snow capping

    real(r8) , intent(in)    :: rho_orig_bottom( bounds%begc: )   ! partial density of ice in bottom snow layer, before updates (not scaled with frac_sno) (kg/m3)
    real(r8) , intent(in)    :: h2osoi_ice_bottom( bounds%begc: ) ! ice lens in bottom snow layer (kg/m2)
    real(r8) , intent(in)    :: frac_adjust( bounds%begc: )       ! fraction of mass remaining after capping
    real(r8) , intent(inout) :: dz_bottom( bounds%begc: )         ! layer depth of bottom snow layer (m)
    real(r8) , intent(inout) :: mss_bcphi_bottom( bounds%begc: )  ! hydrophilic BC mass in snow, bottom layer (kg)
    real(r8) , intent(inout) :: mss_bcpho_bottom( bounds%begc: )  ! hydrophobic BC mass in snow, bottom layer (kg)
    real(r8) , intent(inout) :: mss_ocphi_bottom( bounds%begc: )  ! hydrophilic OC mass in snow, bottom layer (kg)
    real(r8) , intent(inout) :: mss_ocpho_bottom( bounds%begc: )  ! hydrophobic OC mass in snow, bottom layer (kg)
    real(r8) , intent(inout) :: mss_dst1_bottom( bounds%begc: )   ! dust species 1 mass in snow, bottom layer (kg)
    real(r8) , intent(inout) :: mss_dst2_bottom( bounds%begc: )   ! dust species 2 mass in snow, bottom layer (kg)
    real(r8) , intent(inout) :: mss_dst3_bottom( bounds%begc: )   ! dust species 3 mass in snow, bottom layer (kg)
    real(r8) , intent(inout) :: mss_dst4_bottom( bounds%begc: )   ! dust species 4 mass in snow, bottom layer (kg)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c

    character(len=*), parameter :: subname = 'SnowCappingUpdateDzAndAerosols'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(rho_orig_bottom, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osoi_ice_bottom, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(frac_adjust, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(dz_bottom, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(mss_bcphi_bottom, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(mss_bcpho_bottom, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(mss_ocphi_bottom, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(mss_ocpho_bottom, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(mss_dst1_bottom, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(mss_dst2_bottom, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(mss_dst3_bottom, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(mss_dst4_bottom, 1) == bounds%endc), sourcefile, __LINE__)

    do fc = 1, snow_capping_filterc%num
       c = snow_capping_filterc%indices(fc)

       ! Scale dz such that ice density (or: pore space) is conserved
       !
       ! Avoid scaling dz for very low ice densities. This can occur, in principle, if
       ! the layer is mostly liquid water. Furthermore, this check is critical in the
       ! unlikely event that rho is 0, which can happen if the layer is entirely liquid
       ! water.
       if (rho_orig_bottom(c) > 1.0_r8) then
          dz_bottom(c) = h2osoi_ice_bottom(c) / rho_orig_bottom(c)
       end if

       ! Correct the bottom layer aerosol mass to account for snow capping.
       ! This approach conserves the aerosol mass concentration but not aerosol mass.
       mss_bcphi_bottom(c) = mss_bcphi_bottom(c) * frac_adjust(c)
       mss_bcpho_bottom(c) = mss_bcpho_bottom(c) * frac_adjust(c)
       mss_ocphi_bottom(c) = mss_ocphi_bottom(c) * frac_adjust(c)
       mss_ocpho_bottom(c) = mss_ocpho_bottom(c) * frac_adjust(c)
       mss_dst1_bottom(c)  = mss_dst1_bottom(c) * frac_adjust(c)
       mss_dst2_bottom(c)  = mss_dst2_bottom(c) * frac_adjust(c)
       mss_dst3_bottom(c)  = mss_dst3_bottom(c) * frac_adjust(c)
       mss_dst4_bottom(c)  = mss_dst4_bottom(c) * frac_adjust(c)

    end do

  end subroutine SnowCappingUpdateDzAndAerosols


  !-----------------------------------------------------------------------
  subroutine NewSnowBulkDensity(bounds, num_c, filter_c, atm2lnd_inst, bifall)
    !
    ! !DESCRIPTION:
    ! Compute the bulk density of any newly-fallen snow.
    !
    ! The return value is placed in bifall. Only columns within the given filter are set:
    ! all other columns remain at their original values.
    !
    ! !ARGUMENTS:
    type(bounds_type)  , intent(in)    :: bounds
    integer            , intent(in)    :: num_c                ! number of columns in filterc
    integer            , intent(in)    :: filter_c(:)          ! column-level filter to operate on
    type(atm2lnd_type) , intent(in)    :: atm2lnd_inst
    real(r8)           , intent(inout) :: bifall(bounds%begc:) ! bulk density of newly fallen dry snow [kg/m3]
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c, g
    real(r8) :: t_for_bifall_degC  ! temperature to use in bifall equation (deg C)

    character(len=*), parameter :: subname = 'NewSnowBulkDensity'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(bifall) == (/bounds%endc/)), sourcefile, __LINE__)

    associate( &
         forc_t      => atm2lnd_inst%forc_t_downscaled_col , & ! Input:  [real(r8) (:)   ]  atmospheric temperature (Kelvin)        
         forc_wind   => atm2lnd_inst%forc_wind_grc           & ! Input:  [real(r8) (:)   ]  atmospheric wind speed (m/s)
         )

    do fc = 1, num_c
       c = filter_c(fc)
       g = col%gridcell(c)

       if (forc_t(c) > tfrz + 2._r8) then
          bifall(c) = 50._r8 + 1.7_r8*(17.0_r8)**1.5_r8
       else if (forc_t(c) > tfrz - 15._r8) then
          bifall(c) = 50._r8 + 1.7_r8*(forc_t(c) - tfrz + 15._r8)**1.5_r8
       else if ( new_snow_density == LoTmpDnsTruncatedAnderson1976 ) then
          bifall(c) = 50._r8
       else if (new_snow_density == LoTmpDnsSlater2017) then 
          ! Andrew Slater: A temp of about -15C gives the nicest
          ! "blower" powder, but as you get colder the flake size decreases so
          ! density goes up. e.g. the smaller snow crystals from the Arctic and Antarctic
          ! winters
          if (forc_t(c) > tfrz - 57.55_r8) then
             t_for_bifall_degC = (forc_t(c)-tfrz)
          else
             ! Below -57.55 deg C, the following function starts to decrease with
             ! decreasing temperatures. Limit the function to avoid this turning over.
             t_for_bifall_degC = -57.55_r8
          end if
          bifall(c) = -(50._r8/15._r8 + 0.0333_r8*15_r8)*t_for_bifall_degC - 0.0333_r8*t_for_bifall_degC**2
       end if

       if (wind_dependent_snow_density .and. forc_wind(g) > 0.1_r8 ) then
       ! Density offset for wind-driven compaction, initial ideas based on Liston et. al (2007) J. Glaciology,
       ! 53(181), 241-255. Modified for a continuous wind impact and slightly more sensitive
       ! to wind - Andrew Slater, 2016
          bifall(c) = bifall(c) + (266.861_r8 * ((1._r8 + &
                      TANH(forc_wind(g)/params_inst%wind_snowcompact_fact))/2._r8)**8.8_r8)
       end if

    end do

    end associate

  end subroutine NewSnowBulkDensity

  !-----------------------------------------------------------------------
  pure function OverburdenCompactionAnderson1976(burden, wx, td, bi) &
       result(compaction_rate)
    !
    ! !DESCRIPTION:
    ! Compute snow overburden compaction for a single column and level using the Anderson
    ! 1976 formula
    !
    ! From Anderson 1976: A point energy and mass balance model of a snow cover, NOAA
    ! Technical Report NWS 19
    !
    ! !ARGUMENTS:
    real(r8) :: compaction_rate ! function result
    real(r8) , intent(in) :: burden ! pressure of overlying snow in this column [kg/m2]
    real(r8) , intent(in) :: wx     ! water mass (ice+liquid) [kg/m2]
    real(r8) , intent(in) :: td     ! t_soisno - tfrz [K]
    real(r8) , intent(in) :: bi     ! partial density of ice [kg/m3]
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: c2 = 23.e-3_r8       ! [m3/kg]

    character(len=*), parameter :: subname = 'OverburdenCompactionAnderson1976'
    !-----------------------------------------------------------------------

    compaction_rate = -(burden+wx/2._r8)*exp(-overburden_compress_Tfactor*td - c2*bi)/params_inst%eta0_anderson

  end function OverburdenCompactionAnderson1976

  !-----------------------------------------------------------------------
  function OverburdenCompactionVionnet2012(h2osoi_liq, dz, burden, wx, td, bi) &
       result(compaction_rate)
    !
    ! !DESCRIPTION:
    ! Compute snow overburden compaction for a single column and level using the Vionnet
    ! et al. 2012 formula
    !
    ! From Vionnet V et al. 2012, "The detailed snowpack scheme Crocus and its
    ! implementation in SURFEX v7.2", Geosci. Model Dev. 5, 773–791.
    !
    ! Preconditions (required to avoid divide by 0):
    ! - dz > 0
    ! - bi > 0
    !
    ! !ARGUMENTS:
    real(r8) :: compaction_rate ! function result
    real(r8) , intent(in) :: h2osoi_liq ! liquid water in this column and level [kg/m2]
    real(r8) , intent(in) :: dz         ! layer depth for this column and level [m]
    real(r8) , intent(in) :: burden     ! pressure of overlying snow in this column [kg/m2]
    real(r8) , intent(in) :: wx         ! water mass (ice+liquid) [kg/m2]
    real(r8) , intent(in) :: td         ! t_soisno - tfrz [K]
    real(r8) , intent(in) :: bi         ! partial density of ice [kg/m3]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: f1, f2                          ! overburden compaction modifiers to viscosity
    real(r8) :: eta                             ! Viscosity

    real(r8), parameter :: aeta = 0.1_r8        ! overburden compaction constant [1/K]
    real(r8), parameter :: beta = 0.023_r8      ! overburden compaction constant [m3/kg]

    character(len=*), parameter :: subname = 'OverburdenCompactionVionnet2012'
    !-----------------------------------------------------------------------

    f1 = 1._r8 / (1._r8 + 60._r8 * h2osoi_liq / (denh2o * dz))
    f2 = 4.0_r8 ! currently fixed to maximum value, holds in absence of angular grains
    eta = f1*f2*(bi/params_inst%ceta)*exp(aeta*td + beta*bi)*params_inst%eta0_vionnet
    compaction_rate = -(burden+wx/2._r8) / eta

  end function OverburdenCompactionVionnet2012

  !-----------------------------------------------------------------------
  subroutine WindDriftCompaction(bi, forc_wind, dz, &
       zpseudo, mobile, compaction_rate)
    !
    ! !DESCRIPTION:
    !
    ! Compute wind drift compaction for a single column and level.
    !
    ! Also updates zpseudo and mobile for this column. However, zpseudo remains unchanged
    ! if mobile is already false or becomes false within this subroutine.
    !
    ! The structure of the updates done here for zpseudo and mobile requires that this
    ! subroutine be called first for the top layer of snow, then for the 2nd layer down,
    ! etc. - and finally for the bottom layer. Before beginning the loops over layers,
    ! mobile should be initialized to .true. and zpseudo should be initialized to 0.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8) , intent(in)    :: bi              ! partial density of ice [kg/m3]
    real(r8) , intent(in)    :: forc_wind       ! atmospheric wind speed [m/s]
    real(r8) , intent(in)    :: dz              ! layer depth for this column and level [m]
    real(r8) , intent(inout) :: zpseudo         ! wind drift compaction / pseudo depth for this column at this layer
    logical  , intent(inout) :: mobile          ! whether this snow column is still mobile at this layer (i.e., susceptible to wind drift)
    real(r8) , intent(out)   :: compaction_rate ! rate of compaction of snowpack due to wind drift, for the current column and layer
    !
    ! !LOCAL VARIABLES:
    real(r8) :: Frho        ! Mobility density factor [-]
    real(r8) :: MO          ! Mobility index [-]
    real(r8) :: SI          ! Driftability index [-]
    real(r8) :: gamma_drift ! Scaling factor for wind drift time scale [-]
    real(r8) :: tau_inverse ! Inverse of the effective time scale [1/s]

    real(r8), parameter :: rho_min = 50._r8      ! wind drift compaction / minimum density [kg/m3]
    real(r8), parameter :: drift_sph = 1.0_r8    ! wind drift compaction / sphericity

    character(len=*), parameter :: subname = 'WindDriftCompaction'
    !-----------------------------------------------------------------------

    if (mobile) then
       Frho = 1.25_r8 - 0.0042_r8*(max(rho_min, bi)-rho_min)
       ! assuming dendricity = 0, sphericity = 1, grain size = 0.35 mm Non-dendritic snow
       MO = 0.34_r8 * (-0.583_r8*params_inst%drift_gs - 0.833_r8*drift_sph + 0.833_r8) + 0.66_r8*Frho
       SI = -2.868_r8 * exp(-0.085_r8*forc_wind) + 1._r8 + MO

       if (SI > 0.0_r8) then
          SI = min(SI, 3.25_r8)
          ! Increase zpseudo (wind drift / pseudo depth) to the middle of
          ! the pseudo-node for the sake of the following calculation
          zpseudo = zpseudo + 0.5_r8 * dz * (3.25_r8 - SI)
          gamma_drift = SI*exp(-zpseudo/0.1_r8)
          tau_inverse = gamma_drift / params_inst%tau_ref
          compaction_rate = -max(0.0_r8, params_inst%rho_max-bi) * tau_inverse
          ! Further increase zpseudo to the bottom of the pseudo-node for
          ! the sake of calculations done on the underlying layer (i.e.,
          ! the next time through the j loop).
          zpseudo = zpseudo + 0.5_r8 * dz * (3.25_r8 - SI)
       else  ! SI <= 0
          mobile = .false.
          compaction_rate = 0._r8
       end if
    else  ! .not. mobile
       compaction_rate = 0._r8
    end if

  end subroutine WindDriftCompaction

  !-----------------------------------------------------------------------
  subroutine Combo(dz,  wliq,  wice, t, dz2, wliq2, wice2, t2)
    !
    ! !DESCRIPTION:
    ! Combines two elements and returns the following combined
    ! variables: dz, t, wliq, wice.
    ! The combined temperature is based on the equation:
    ! the sum of the enthalpies of the two elements =
    ! that of the combined element.
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)    :: dz2   ! nodal thickness of 2 elements being combined [m]
    real(r8), intent(in)    :: wliq2 ! liquid water of element 2 [kg/m2]
    real(r8), intent(in)    :: wice2 ! ice of element 2 [kg/m2]
    real(r8), intent(in)    :: t2    ! nodal temperature of element 2 [K]
    real(r8), intent(inout) :: dz    ! nodal thickness of 1 elements being combined [m]
    real(r8), intent(inout) :: wliq  ! liquid water of element 1
    real(r8), intent(inout) :: wice  ! ice of element 1 [kg/m2]
    real(r8), intent(inout) :: t     ! nodel temperature of element 1 [K]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dzc   ! Total thickness of nodes 1 and 2 (dzc=dz+dz2).
    real(r8) :: wliqc ! Combined liquid water [kg/m2]
    real(r8) :: wicec ! Combined ice [kg/m2]
    real(r8) :: tc    ! Combined node temperature [K]
    real(r8) :: h     ! enthalpy of element 1 [J/m2]
    real(r8) :: h2    ! enthalpy of element 2 [J/m2]
    real(r8) :: hc    ! temporary
    !-----------------------------------------------------------------------

    dzc = dz+dz2
    wicec = (wice+wice2)
    wliqc = (wliq+wliq2)
    h = (cpice*wice+cpliq*wliq) * (t-tfrz)+hfus*wliq
    h2= (cpice*wice2+cpliq*wliq2) * (t2-tfrz)+hfus*wliq2

    hc = h + h2
    tc = tfrz + (hc - hfus*wliqc) / (cpice*wicec + cpliq*wliqc)

    dz = dzc
    wice = wicec
    wliq = wliqc
    t = tc

  end subroutine Combo

  !-----------------------------------------------------------------------
  function MassWeightedSnowRadius( rds1, rds2, swtot, zwtot ) result(mass_weighted_snowradius)
    !
    ! !DESCRIPTION:
    ! Calculate the mass weighted snow radius when two layers are combined
    !
    ! !USES:
    use AerosolMod   , only : snw_rds_min
    use SnowSnicarMod, only : snw_rds_max
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(IN) :: rds1         ! Layer 1 radius
    real(r8), intent(IN) :: rds2         ! Layer 2 radius
    real(r8), intent(IN) :: swtot        ! snow water total layer 2
    real(r8), intent(IN) :: zwtot        ! snow water total layer 1
    real(r8) :: mass_weighted_snowradius ! resulting bounded mass weighted snow radius

    SHR_ASSERT_FL( (swtot+zwtot > 0.0_r8), sourcefile, __LINE__)
    mass_weighted_snowradius = (rds2*swtot + rds1*zwtot)/(swtot+zwtot)

    if (      mass_weighted_snowradius > snw_rds_max ) then
       mass_weighted_snowradius = snw_rds_max
    else if ( mass_weighted_snowradius < snw_rds_min ) then
       mass_weighted_snowradius = snw_rds_min
    end if
  end function MassWeightedSnowRadius

  !-----------------------------------------------------------------------
  subroutine BuildSnowFilter(bounds, num_nolakec, filter_nolakec, &
       num_snowc, filter_snowc, num_nosnowc, filter_nosnowc)
    !
    ! !DESCRIPTION:
    ! Constructs snow filter for use in vectorized loops for snow hydrology.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)  :: bounds
    integer           , intent(in)  :: num_nolakec       ! number of column non-lake points in column filter
    integer           , intent(in)  :: filter_nolakec(:) ! column filter for non-lake points
    integer           , intent(out) :: num_snowc         ! number of column snow points in column filter
    integer           , intent(out) :: filter_snowc(:)   ! column filter for snow points
    integer           , intent(out) :: num_nosnowc       ! number of column non-snow points in column filter
    integer           , intent(out) :: filter_nosnowc(:) ! column filter for non-snow points
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c
    !-----------------------------------------------------------------------

    ! Build snow/no-snow filters for other subroutines

    num_snowc = 0
    num_nosnowc = 0
    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       if (col%snl(c) < 0) then
          num_snowc = num_snowc + 1
          filter_snowc(num_snowc) = c
       else
          num_nosnowc = num_nosnowc + 1
          filter_nosnowc(num_nosnowc) = c
       end if
    end do
  end subroutine BuildSnowFilter

  subroutine SnowHydrologySetControlForTesting( set_winddep_snowdensity, set_new_snow_density, &
       set_reset_snow, set_reset_snow_glc, set_reset_snow_glc_ela)
    !
    ! !DESCRIPTION:
    ! Sets some of the control settings for SnowHydrologyMod
    ! NOTE: THIS IS JUST HERE AS AN INTERFACE FOR UNIT TESTING.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    logical, intent(in), optional :: set_winddep_snowdensity  ! Set wind dependent snow density
    integer, intent(in), optional :: set_new_snow_density     ! snow density method
    logical, intent(in), optional :: set_reset_snow           ! whether to reset the snow pack, non-glacier points
    logical, intent(in), optional :: set_reset_snow_glc       ! whether to reset the snow pack, glacier points
    real(r8), intent(in), optional :: set_reset_snow_glc_ela  ! elevation below which to reset the snow pack if set_reset_snow_glc is true (m)
    !-----------------------------------------------------------------------
    if (present(set_winddep_snowdensity)) then
       wind_dependent_snow_density = set_winddep_snowdensity
    end if
    if (present(set_new_snow_density)) then
       new_snow_density            = set_new_snow_density
    end if
    if (present(set_reset_snow)) then
       reset_snow = set_reset_snow
    end if
    if (present(set_reset_snow_glc)) then
       reset_snow_glc = set_reset_snow_glc
    end if
    if (present(set_reset_snow_glc_ela)) then
       reset_snow_glc_ela = set_reset_snow_glc_ela
    end if
    snow_dzmin_1 = 0.010_r8  ! The same default values specified in...
    snow_dzmin_2 = 0.015_r8  ! /bld/namelist_files/namelist_defaults_ctsm.xml
    snow_dzmax_l_1 = 0.03_r8  ! and used when alternate values do not
    snow_dzmax_l_2 = 0.07_r8  ! get set in
    snow_dzmax_u_1 = 0.02_r8  ! user_nl_clm
    snow_dzmax_u_2 = 0.05_r8

    params_inst%wind_snowcompact_fact = 5.0_r8

  end subroutine SnowHydrologySetControlForTesting

  !-----------------------------------------------------------------------
  subroutine SnowHydrologyClean()
    !
    ! !DESCRIPTION:
    ! Deallocate memory
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'Clean'
    !-----------------------------------------------------------------------

     deallocate(dzmin)
     deallocate(dzmax_l)
     deallocate(dzmax_u)

  end subroutine SnowHydrologyClean

end module SnowHydrologyMod
