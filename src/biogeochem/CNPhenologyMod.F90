module CNPhenologyMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !MODULE: CNPhenologyMod
  !
  ! !DESCRIPTION:
  ! Module holding routines used in phenology model for coupled carbon
  ! nitrogen code.
  !
  ! !USES:
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use shr_log_mod                     , only : errMsg => shr_log_errMsg
  use shr_sys_mod                     , only : shr_sys_flush
  use decompMod                       , only : bounds_type
  use clm_varpar                      , only : ileaf,ileaf_st,ileaf_xf,ifroot,ifroot_st,ifroot_xf,&
                                               ilivestem,ilivestem_st,ilivestem_xf,&
                                               ideadstem,ideadstem_st,ideadstem_xf,&
                                               ilivecroot,ilivecroot_st,ilivecroot_xf,&
                                               ideadcroot,ideadcroot_st,ideadcroot_xf,&
                                               igrain,igrain_st,igrain_xf,iretransn,ioutc,ioutn
  use clm_varpar                      , only : maxveg, nlevdecomp_full, mxsowings, mxharvests
  use clm_varpar                      , only : i_litr_min, i_litr_max
  use clm_varctl                      , only : iulog, use_cndv
  use CNSharedParamsMod               , only : use_matrixcn
  use clm_varctl                      , only : for_testing_no_crop_seed_replenishment
  use clm_varcon                      , only : tfrz
  use abortutils                      , only : endrun
  use CanopyStateType                 , only : canopystate_type
  use CNDVType                        , only : dgvs_type
  use CNVegstateType                  , only : cnveg_state_type
  use CNVegCarbonStateType            , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType             , only : cnveg_carbonflux_type
  use CNVegnitrogenstateType          , only : cnveg_nitrogenstate_type
  use CNVegnitrogenfluxType           , only : cnveg_nitrogenflux_type
  use CropType                        , only : crop_type
  use CropType                        , only : cphase_planted, cphase_leafemerge
  use CropType                        , only : cphase_grainfill, cphase_harvest
  use pftconMod                       , only : pftcon
  use SoilStateType                   , only : soilstate_type
  use TemperatureType                 , only : temperature_type
  use WaterDiagnosticBulkType         , only : waterdiagnosticbulk_type
  use Wateratm2lndBulkType            , only : wateratm2lndbulk_type
  use initVerticalMod                 , only : find_soil_layer_containing_depth
  use CropReprPoolsMod                , only : nrepr, repr_grain_min, repr_grain_max, repr_structure_min, repr_structure_max
  use ColumnType                      , only : col
  use GridcellType                    , only : grc                
  use PatchType                       , only : patch   
  use atm2lndType                     , only : atm2lnd_type             
  use CNVegMatrixMod                  , only : matrix_update_phc, matrix_update_phn
  use CNVegMatrixMod                  , only : matrix_update_gmc, matrix_update_gmn
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readParams           ! Read parameters
  public :: CNPhenologyreadNML   ! Read namelist
  public :: CNPhenologyInit      ! Initialization
  public :: CNPhenology          ! Update
  public :: CropPhase            ! Get the current phase of each crop patch
  public :: DaysPastPlanting     ! Get how many days it's been since crop was planted

  ! !PUBLIC for unit testing
  public :: CNPhenologySetNML         ! Set the namelist setttings explicitly for unit tests
  public :: CNPhenologySetParams      ! Set the parameters explicitly for unit tests
  public :: SeasonalDecidOnset        ! Logical function to determine is seasonal decidious onset should be triggered
  public :: SeasonalCriticalDaylength ! Critical day length needed for Seasonal decidious offset
  public :: get_swindow
  public :: was_sown_in_this_window

  ! !PRIVITE MEMBER FIUNCTIONS:
  private :: CNPhenologyClimate             ! Get climatological everages to figure out triggers for Phenology
  private :: CNEvergreenPhenology           ! Phenology for evergreen plants
  private :: CNSeasonDecidPhenology         ! Phenology for seasonal decidious platnts
  private :: CNStressDecidPhenology         ! Phenology for stress deciidous plants
  private :: CropPhenology                  ! Phenology for crops
  private :: CropPhenologyInit              ! Initialize phenology for crops
  private :: vernalization                  ! Vernalization (overwinterring) of crops
  private :: CNOnsetGrowth                  ! Leaf Onset growth
  private :: CNOffsetLitterfall             ! Leaf Offset litter fall
  private :: CNBackgroundLitterfall         ! Background litter fall
  private :: CNLivewoodTurnover             ! Liver wood turnover to deadwood
  private :: CNCropHarvestToProductPools    ! Move crop harvest to product pools
  private :: CNLitterToColumn               ! Move litter ofrom patch to column level
  !
  ! !PRIVATE DATA MEMBERS:
  type, private :: params_type
     real(r8) :: crit_dayl             ! critical day length for senescence (sec) 
                                       ! (11 hrs [39,600 sec] from White 2001)
     real(r8) :: crit_dayl_at_high_lat ! critical day length for senescence at high latitudes (sec) 
                                       ! (in Eitel 2019 this was 54000 [15 hrs])
     real(r8) :: crit_dayl_lat_slope   ! Slope of time for critical day length with latitude (sec/deg) 
                                       ! (Birch et. all 2021 it was 720 see line below)
                                       ! 15hr-11hr/(65N-45N)=linear slope = 720 min/latitude (Birch et. al 2021)
     real(r8) :: ndays_off             ! number of days to complete leaf offset
     real(r8) :: fstor2tran            ! fraction of storage to move to transfer for each onset
     real(r8) :: crit_onset_fdd        ! critical number of freezing days to set gdd counter
     real(r8) :: crit_onset_swi        ! critical number of days > soilpsi_on for onset
     real(r8) :: soilpsi_on            ! critical soil water potential for leaf onset
     real(r8) :: crit_offset_fdd       ! critical number of freezing days to initiate offset
     real(r8) :: crit_offset_swi       ! critical number of water stress days to initiate offset
     real(r8) :: soilpsi_off           ! critical soil water potential for leaf offset
     real(r8) :: lwtop                 ! live wood turnover proportion (annual fraction)
     real(r8) :: phenology_soil_depth  ! soil depth used for measuring states for phenology triggers
     real(r8) :: snow5d_thresh_for_onset ! 5-day snow depth threshold for leaf onset
  end type params_type

  type(params_type) :: params_inst

  real(r8) :: dt                            ! time step delta t (seconds)
  real(r8) :: fracday                       ! dtime as a fraction of day
  real(r8) :: crit_dayl                     ! critical daylength for offset (seconds)
  real(r8) :: ndays_off                     ! number of days to complete offset
  real(r8) :: fstor2tran                    ! fraction of storage to move to transfer on each onset
  real(r8) :: crit_onset_fdd                ! critical number of freezing days
  real(r8) :: crit_onset_swi                ! water stress days for offset trigger
  real(r8) :: soilpsi_on                    ! water potential for onset trigger (MPa)
  real(r8) :: crit_offset_fdd               ! critical number of freezing degree days to trigger offset
  real(r8) :: crit_offset_swi               ! water stress days for offset trigger
  real(r8) :: soilpsi_off                   ! water potential for offset trigger (MPa)
  real(r8) :: lwtop                         ! live wood turnover proportion (annual fraction)
  integer  :: phenology_soil_layer          ! soil layer used for measuring states for phenology triggers

  ! CropPhenology variables and constants
  real(r8) :: p1d, p1v                      ! photoperiod factor constants for crop vernalization
  real(r8) :: hti                           ! cold hardening index threshold for vernalization
  real(r8) :: tbase                         ! base temperature for vernalization

  integer, parameter :: NOT_Planted   = 999 ! If not planted   yet in year
  integer, parameter :: NOT_Harvested = 999 ! If not harvested yet in year
  integer, parameter :: inNH       = 1      ! Northern Hemisphere
  integer, parameter :: inSH       = 2      ! Southern Hemisphere
  integer, pointer   :: inhemi(:)           ! Hemisphere that patch is in 

  integer, allocatable :: minplantjday(:,:) ! minimum planting julian day
  integer, allocatable :: maxplantjday(:,:) ! maximum planting julian day
  integer              :: jdayyrstart(inSH) ! julian day of start of year

  logical,parameter :: matrixcheck_ph = .True.                   ! Matrix check
  logical,parameter :: acc_ph = .False.                          ! Another matrix check

  real(r8), private :: initial_seed_at_planting        = 3._r8   ! Initial seed at planting

  real(r8)         :: min_gddmaturity = 1._r8     ! Weird things can happen if gddmaturity is tiny
  logical,  public :: generate_crop_gdds = .false. ! If true, harvest the day before next sowing
  logical,  public :: use_mxmat = .true.           ! If true, ignore crop maximum growing season length

  ! For use with adapt_cropcal_rx_cultivar_gdds .true.
  real(r8), parameter :: min_gdd20_baseline = 0._r8  ! If gdd20_baseline_patch is ≤ this, do not consider baseline.

  ! Constants for seasonal decidious leaf onset and offset
  logical,  private :: onset_thresh_depends_on_veg     = .false. ! If onset threshold depends on vegetation type
  integer,  public, parameter :: critical_daylight_constant           = 1
  integer,  public, parameter :: critical_daylight_depends_on_lat     = critical_daylight_constant + 1
  integer,  public, parameter :: critical_daylight_depends_on_veg     = critical_daylight_depends_on_lat + 1
  integer,  public, parameter :: critical_daylight_depends_on_latnveg = critical_daylight_depends_on_veg + 1
  integer,  private :: critical_daylight_method = critical_daylight_constant
  ! For determining leaf offset latitude that's considered high latitude (see Eitel 2019)
  real(r8), parameter :: critical_offset_high_lat         = 65._r8     ! Start of what's considered high latitude (degrees)

  real(r8), parameter :: HARVEST_REASON_MATURE = 1._r8
  real(r8), parameter :: HARVEST_REASON_MAXSEASLENGTH = 2._r8
  real(r8), parameter :: HARVEST_REASON_SOWNBADDEC31 = 3._r8
  real(r8), parameter :: HARVEST_REASON_SOWTODAY = 4._r8
  real(r8), parameter :: HARVEST_REASON_SOWTOMORROW = 5._r8
  real(r8), parameter :: HARVEST_REASON_IDOPTOMORROW = 6._r8
  real(r8), parameter :: HARVEST_REASON_VERNFREEZEKILL = 7._r8

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNPhenologyReadNML( NLFilename )
    !
    ! !DESCRIPTION:
    ! Read the namelist for CNPhenology
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=25) :: min_critical_dayl_method   ! Method to determine critical day length for onset
    character(len=*), parameter :: subname = 'CNPhenologyReadNML'
    character(len=*), parameter :: nmlname = 'cnphenology'
    !-----------------------------------------------------------------------
    namelist /cnphenology/ initial_seed_at_planting, onset_thresh_depends_on_veg, &
                           min_critical_dayl_method, generate_crop_gdds, &
                           use_mxmat

    ! Initialize options to default values, in case they are not specified in
    ! the namelist

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=cnphenology, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (initial_seed_at_planting,    mpicom)
    call shr_mpi_bcast (onset_thresh_depends_on_veg, mpicom)
    call shr_mpi_bcast (min_critical_dayl_method,     mpicom)
    call shr_mpi_bcast (generate_crop_gdds,          mpicom)
    call shr_mpi_bcast (use_mxmat,                   mpicom)

    if (      min_critical_dayl_method == "DependsOnLat"       )then
       critical_daylight_method = critical_daylight_depends_on_lat
    else if ( min_critical_dayl_method == "DependsOnVeg"       )then
       critical_daylight_method = critical_daylight_depends_on_veg
    else if ( min_critical_dayl_method == "DependsOnLatAndVeg" )then
       critical_daylight_method = critical_daylight_depends_on_latnveg
    else if ( min_critical_dayl_method == "Constant"           )then
       critical_daylight_method = critical_daylight_constant
    else
       call endrun(msg="ERROR min_critical_dayl_method is NOT set to a valid value"//errmsg(sourcefile, __LINE__))
    end if

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=cnphenology)
       write(iulog,*) ' '
    end if


    !-----------------------------------------------------------------------
    
  end subroutine CNPhenologyReadNML

  !-----------------------------------------------------------------------
  subroutine CNPhenologySetNML( input_onset_thresh_depends_on_veg, input_critical_daylight_method )
    !
    ! !DESCRIPTION:
    ! Set the namelist items for unit-testing
    !
    logical, intent(in) :: input_onset_thresh_depends_on_veg
    integer, intent(in) :: input_critical_daylight_method

    onset_thresh_depends_on_veg = input_onset_thresh_depends_on_veg
    critical_daylight_method = input_critical_daylight_method
    if ( (critical_daylight_method < critical_daylight_constant) .or. &
         (critical_daylight_method > critical_daylight_depends_on_latnveg) )then
       call endrun(msg="ERROR critical_daylight_method "//errmsg(sourcefile, __LINE__))
    end if
  end subroutine CNPhenologySetNML
  
  !-----------------------------------------------------------------------
  subroutine CNPhenologySetParams( )
    !
    ! !DESCRIPTION:
    ! Set the parameters for unit-testing
    !
    params_inst%crit_dayl             = 39200._r8     ! Seconds
    params_inst%crit_dayl_at_high_lat = 54000._r8     ! Seconds
    params_inst%crit_dayl_lat_slope   = 720._r8       ! Seconds / degree
    params_inst%ndays_off             = 30._r8        ! Days
    params_inst%fstor2tran            = 0.5           ! Fraction
    params_inst%crit_onset_fdd        = 15._r8        ! Days
    params_inst%crit_onset_swi        = 15._r8        ! Days
    params_inst%soilpsi_on            = -0.6_r8       ! MPa
    params_inst%crit_offset_fdd       = 15._r8        ! Days
    params_inst%crit_offset_swi       = 15._r8        ! Days
    params_inst%soilpsi_off           = -0.8          ! MPa
    params_inst%lwtop                 = 0.7_r8        ! Fraction
    params_inst%phenology_soil_depth  = 0.08_r8       ! m
    params_inst%snow5d_thresh_for_onset = 0.2_r8      ! m
  end subroutine CNPhenologySetParams
  
  !-----------------------------------------------------------------------
  subroutine readParams ( ncid )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ncdio_pio    , only: file_desc_t
    use paramUtilMod , only : readNcdioScalar

    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter  :: subname = 'readParams_CNPhenology'
    !-----------------------------------------------------------------------

    call readNcdioScalar(ncid, 'crit_dayl', subname, params_inst%crit_dayl)
    call readNcdioScalar(ncid, 'crit_dayl_at_high_lat', subname, params_inst%crit_dayl_at_high_lat)
    call readNcdioScalar(ncid, 'crit_dayl_lat_slope', subname, params_inst%crit_dayl_lat_slope)
    call readNcdioScalar(ncid, 'ndays_off', subname, params_inst%ndays_off)
    call readNcdioScalar(ncid, 'fstor2tran', subname, params_inst%fstor2tran)
    call readNcdioScalar(ncid, 'crit_onset_fdd', subname, params_inst%crit_onset_fdd)
    call readNcdioScalar(ncid, 'crit_onset_swi', subname, params_inst%crit_onset_swi)
    call readNcdioScalar(ncid, 'soilpsi_on', subname, params_inst%soilpsi_on)
    call readNcdioScalar(ncid, 'crit_offset_fdd', subname, params_inst%crit_offset_fdd)
    call readNcdioScalar(ncid, 'crit_offset_swi', subname, params_inst%crit_offset_swi)
    call readNcdioScalar(ncid, 'soilpsi_off', subname, params_inst%soilpsi_off)
    call readNcdioScalar(ncid, 'lwtop_ann', subname, params_inst%lwtop)
    call readNcdioScalar(ncid, 'phenology_soil_depth', subname, params_inst%phenology_soil_depth)
    call readNcdioScalar(ncid, 'snow5d_thresh_for_onset', subname, params_inst%snow5d_thresh_for_onset)
    

  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine CNPhenology (bounds, num_soilc, filter_soilc, num_soilp, &
       filter_soilp, num_pcropp, filter_pcropp, &
       waterdiagnosticbulk_inst, wateratm2lndbulk_inst, temperature_inst, atm2lnd_inst, crop_inst, &
       canopystate_inst, soilstate_inst, dgvs_inst, &
       cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst,    &
       cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
       c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst, &
       leaf_prof_patch, froot_prof_patch, phase)
    ! !USES:
    use clm_time_manager , only: is_first_step
    use CNSharedParamsMod, only: use_fun
    !
    ! !DESCRIPTION:
    ! Dynamic phenology routine for coupled carbon-nitrogen code (CN)
    ! 1. grass phenology
    !
    ! !ARGUMENTS:
    type(bounds_type)              , intent(in)    :: bounds
    integer                        , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                        , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                        , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                        , intent(in)    :: filter_soilp(:) ! filter for soil patches
    integer                        , intent(in)    :: num_pcropp      ! number of prog. crop patches in filter
    integer                        , intent(in)    :: filter_pcropp(:)! filter for prognostic crop patches
    type(waterdiagnosticbulk_type)          , intent(in)    :: waterdiagnosticbulk_inst
    type(wateratm2lndbulk_type)          , intent(in)    :: wateratm2lndbulk_inst
    type(temperature_type)         , intent(inout) :: temperature_inst
    type(atm2lnd_type)             , intent(in)    :: atm2lnd_inst
    type(crop_type)                , intent(inout) :: crop_inst
    type(canopystate_type)         , intent(in)    :: canopystate_inst
    type(soilstate_type)           , intent(in)    :: soilstate_inst
    type(dgvs_type)                , intent(inout) :: dgvs_inst
    type(cnveg_state_type)         , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenstate_type) , intent(inout) :: cnveg_nitrogenstate_inst
    type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: c13_cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: c14_cnveg_carbonstate_inst
    real(r8)                       , intent(in)    :: leaf_prof_patch(bounds%begp:,1:)
    real(r8)                       , intent(in)    :: froot_prof_patch(bounds%begp:,1:)
    integer                        , intent(in)    :: phase
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(leaf_prof_patch)   == (/bounds%endp,nlevdecomp_full/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(froot_prof_patch)  == (/bounds%endp,nlevdecomp_full/)), sourcefile, __LINE__)

    ! each of the following phenology type routines includes a filter
    ! to operate only on the relevant patches


    if ( phase == 1 ) then
       call CNPhenologyClimate(num_soilp, filter_soilp, &
            temperature_inst, cnveg_state_inst, crop_inst)
   
       call CNEvergreenPhenology(num_soilp, filter_soilp, &
            cnveg_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)

       call CNSeasonDecidPhenology(num_soilp, filter_soilp, &
            temperature_inst, waterdiagnosticbulk_inst, cnveg_state_inst, dgvs_inst, &
            cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)

       call CNStressDecidPhenology(num_soilp, filter_soilp,   &
            soilstate_inst, temperature_inst, atm2lnd_inst, wateratm2lndbulk_inst, cnveg_state_inst, &
            cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)

       if (num_pcropp > 0) then
          call CropPhenology(num_pcropp, filter_pcropp, &
               waterdiagnosticbulk_inst, temperature_inst, crop_inst, canopystate_inst, cnveg_state_inst, &
               cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
               c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst)
       end if
    else if ( phase == 2 ) then
       ! the same onset and offset routines are called regardless of
       ! phenology type - they depend only on onset_flag, offset_flag, bglfr, and bgtr

       call CNOnsetGrowth(num_soilp, filter_soilp, &
            cnveg_state_inst, &
            cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)

       call CNOffsetLitterfall(num_soilp, filter_soilp, &
            cnveg_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
            crop_inst)

       call CNBackgroundLitterfall(num_soilp, filter_soilp, &
            cnveg_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
   
       call CNLivewoodTurnover(num_soilp, filter_soilp, &
            cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)

       call CNCropHarvestToProductPools(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, &
            cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)

       ! gather all patch-level litterfall fluxes to the column for litter C and N inputs

       call CNLitterToColumn(bounds, num_soilp, filter_soilp, &
            cnveg_state_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
            leaf_prof_patch(bounds%begp:bounds%endp,1:nlevdecomp_full), & 
            froot_prof_patch(bounds%begp:bounds%endp,1:nlevdecomp_full))
    else
       call endrun( 'bad phase' )
    end if

  end subroutine CNPhenology

  !-----------------------------------------------------------------------
  subroutine CNPhenologyInit(bounds)
    !
    ! !DESCRIPTION:
    ! Initialization of CNPhenology. Must be called after time-manager is
    ! initialized, and after pftcon file is read in.
    !
    ! !USES:
    use clm_time_manager, only: get_step_size_real
    use clm_varctl      , only: use_crop, use_cropcal_rx_swindows
    use clm_varcon      , only: secspday
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  
    !------------------------------------------------------------------------

    !
    ! Get time-step and what fraction of a day it is
    !
    dt      = get_step_size_real()
    fracday = dt/secspday

    ! set constants for CNSeasonDecidPhenology 
    ! (critical daylength from Biome-BGC, v4.1.2)
    crit_dayl=params_inst%crit_dayl

    ! Set constants for CNSeasonDecidPhenology and CNStressDecidPhenology
    ndays_off=params_inst%ndays_off

    ! set transfer parameters
    fstor2tran=params_inst%fstor2tran

    call find_soil_layer_containing_depth( &
         depth = params_inst%phenology_soil_depth, &
         layer = phenology_soil_layer)

    ! -----------------------------------------
    ! Constants for CNStressDecidPhenology
    ! -----------------------------------------

    ! onset parameters
    crit_onset_fdd=params_inst%crit_onset_fdd
    ! critical onset gdd now being calculated as a function of annual
    ! average 2m temp.
    ! crit_onset_gdd = 150.0 ! c3 grass value
    ! crit_onset_gdd = 1000.0   ! c4 grass value
    crit_onset_swi=params_inst%crit_onset_swi
    soilpsi_on=params_inst%soilpsi_on

    ! offset parameters
    crit_offset_fdd=params_inst%crit_offset_fdd
    crit_offset_swi=params_inst%crit_offset_swi
    soilpsi_off=params_inst%soilpsi_off    

    ! -----------------------------------------
    ! Constants for CNLivewoodTurnover
    ! -----------------------------------------

    ! set the global parameter for livewood turnover rate
    ! define as an annual fraction (0.7), and convert to fraction per second
    lwtop=params_inst%lwtop/31536000.0_r8 !annual fraction converted to per second

    ! -----------------------------------------
    ! Call any subroutine specific initialization routines
    ! -----------------------------------------

    if ( use_crop ) call CropPhenologyInit(bounds)

    ! Error checking for parameters
    if ( (critical_daylight_method == critical_daylight_depends_on_lat) .or. &
         (critical_daylight_method == critical_daylight_depends_on_veg) .or.  &
         (critical_daylight_method == critical_daylight_depends_on_latnveg) )then
        if ( params_inst%crit_dayl_at_high_lat < params_inst%crit_dayl )then
            call endrun(msg="ERROR crit_dayl_at_high_lat should be higher than crit_dayl on paramsfile"//errmsg(sourcefile, __LINE__))
        end if
        if ( params_inst%crit_dayl_at_high_lat >= secspday )then
            call endrun(msg="ERROR crit_dayl_at_high_lat can NOT be higher than seconds in a day"//errmsg(sourcefile, __LINE__))
        end if
        if ( params_inst%crit_dayl >= secspday )then
            call endrun(msg="ERROR crit_dayl can NOT be higher than seconds in a day"//errmsg(sourcefile, __LINE__))
        end if
    end if
    if ( (critical_daylight_method == critical_daylight_depends_on_lat) .or. &
         (critical_daylight_method == critical_daylight_depends_on_latnveg) )then
        if ( params_inst%crit_dayl_lat_slope <= 0.0_r8 )then
            call endrun(msg="ERROR crit_dayl_lat_slope can not be negative or zero"//errmsg(sourcefile, __LINE__))
        end if
        if ( params_inst%crit_dayl_lat_slope >= (secspday - params_inst%crit_dayl_at_high_lat)/ &
                                                 (90._r8 - critical_offset_high_lat) )then
            call endrun(msg="ERROR crit_dayl_lat_slope cannot allow crit_dayl longer than a day"//errmsg(sourcefile, __LINE__))
        end if
    end if

  end subroutine CNPhenologyInit

  !-----------------------------------------------------------------------
  subroutine CNPhenologyClimate (num_soilp, filter_soilp, &
       temperature_inst, cnveg_state_inst, crop_inst)
    !
    ! !DESCRIPTION:
    ! For coupled carbon-nitrogen code (CN).
    !
    ! !USES:
    use clm_time_manager , only : get_curr_days_per_year
    use clm_time_manager , only : get_curr_date, is_first_step
    !
    ! !ARGUMENTS:
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(temperature_type) , intent(inout) :: temperature_inst
    type(cnveg_state_type) , intent(inout) :: cnveg_state_inst
    type(crop_type)        , intent(inout) :: crop_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: p       ! indices
    integer  :: fp      ! filter patch index
    real(r8) :: dayspyr ! days per year (days)
    !-----------------------------------------------------------------------

    associate(                                                & 
         nyrs_crop_active => crop_inst%nyrs_crop_active_patch,   & ! InOut:  [integer (:)  ]  number of years this crop patch has been active
         
         t_ref2m        => temperature_inst%t_ref2m_patch ,   & ! Input:  [real(r8) (:) ]  2m air temperature (K)
         
         tempavg_t2m    => cnveg_state_inst%tempavg_t2m_patch & ! Output: [real(r8) (:) ]  temp. avg 2m air temperature (K)                  
         )

      ! set time steps
      
      dayspyr = get_curr_days_per_year()

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         tempavg_t2m(p) = tempavg_t2m(p) + t_ref2m(p) * (fracday/dayspyr)
      end do

    end associate

  end subroutine CNPhenologyClimate

  !-----------------------------------------------------------------------
  subroutine CNEvergreenPhenology (num_soilp, filter_soilp , &
       cnveg_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)   
       ! cnveg_state_inst) 
    !
    ! !DESCRIPTION:
    ! For coupled carbon-nitrogen code (CN).
    !
    ! !USES:
    use clm_varcon       , only : secspday
    use clm_time_manager , only : get_average_days_per_year
    use clm_varctl       , only : CN_evergreen_phenology_opt   
    !
    ! !ARGUMENTS:
    integer           , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer           , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_state_type), intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: cnveg_carbonstate_inst    
    type(cnveg_nitrogenstate_type) , intent(inout) :: cnveg_nitrogenstate_inst  
    type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst     
    type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst   
    !
    ! !LOCAL VARIABLES:
    real(r8):: avg_dayspyr                ! Average days per year
    integer :: p                          ! indices
    integer :: fp                         ! filter patch index
    
    real(r8):: tranr 				      
    real(r8):: t1                         ! temporary variable 
    !-----------------------------------------------------------------------

    associate(                                        & 
         ivt        => patch%itype                    , & ! Input:  [integer  (:) ]  patch vegetation type                                

         evergreen  => pftcon%evergreen             , & ! Input:  binary flag for evergreen leaf habit (0 or 1)     
         leaf_long  => pftcon%leaf_long             , & ! Input:  leaf longevity (yrs)  
         
         woody                               =>    pftcon%woody                                                         , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)     
         
   		 leafc_storage                       =>    cnveg_carbonstate_inst%leafc_storage_patch                           , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C storage                             
   		 frootc_storage                      =>    cnveg_carbonstate_inst%frootc_storage_patch                          , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C storage                         
   		 livestemc_storage                   =>    cnveg_carbonstate_inst%livestemc_storage_patch                       , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C storage                         
   		 deadstemc_storage                   =>    cnveg_carbonstate_inst%deadstemc_storage_patch                       , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C storage                         
   		 livecrootc_storage                  =>    cnveg_carbonstate_inst%livecrootc_storage_patch                      , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C storage                  
   		 deadcrootc_storage                  =>    cnveg_carbonstate_inst%deadcrootc_storage_patch                      , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C storage                  
   		 gresp_storage                       =>    cnveg_carbonstate_inst%gresp_storage_patch                           , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration storage   
   		 leafc_xfer                          =>    cnveg_carbonstate_inst%leafc_xfer_patch                              , & ! InOut:  [real(r8) (:)]  (gC/m2) leaf C transfer                            
   		 frootc_xfer                         =>    cnveg_carbonstate_inst%frootc_xfer_patch                             , & ! InOut:  [real(r8) (:)]  (gC/m2) fine root C transfer                       
   		 livestemc_xfer                      =>    cnveg_carbonstate_inst%livestemc_xfer_patch                          , & ! InOut:  [real(r8) (:)]  (gC/m2) live stem C transfer                       
   		 deadstemc_xfer                      =>    cnveg_carbonstate_inst%deadstemc_xfer_patch                          , & ! InOut:  [real(r8) (:)]  (gC/m2) dead stem C transfer                       
   		 livecrootc_xfer                     =>    cnveg_carbonstate_inst%livecrootc_xfer_patch                         , & ! InOut:  [real(r8) (:)]  (gC/m2) live coarse root C transfer                
   		 deadcrootc_xfer                     =>    cnveg_carbonstate_inst%deadcrootc_xfer_patch                         , & ! InOut:  [real(r8) (:)]  (gC/m2) dead coarse root C transfer   
   		                
   		 leafn_storage                       =>    cnveg_nitrogenstate_inst%leafn_storage_patch                         , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N storage                              
   		 frootn_storage                      =>    cnveg_nitrogenstate_inst%frootn_storage_patch                        , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N storage                         
   		 livestemn_storage                   =>    cnveg_nitrogenstate_inst%livestemn_storage_patch                     , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N storage                         
   		 deadstemn_storage                   =>    cnveg_nitrogenstate_inst%deadstemn_storage_patch                     , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N storage                         
   		 livecrootn_storage                  =>    cnveg_nitrogenstate_inst%livecrootn_storage_patch                    , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N storage                  
   		 deadcrootn_storage                  =>    cnveg_nitrogenstate_inst%deadcrootn_storage_patch                    , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N storage               
   		 leafn_xfer                          =>    cnveg_nitrogenstate_inst%leafn_xfer_patch                            , & ! InOut:  [real(r8) (:)]  (gN/m2) leaf N transfer                            
   		 frootn_xfer                         =>    cnveg_nitrogenstate_inst%frootn_xfer_patch                           , & ! InOut:  [real(r8) (:)]  (gN/m2) fine root N transfer                       
   		 livestemn_xfer                      =>    cnveg_nitrogenstate_inst%livestemn_xfer_patch                        , & ! InOut:  [real(r8) (:)]  (gN/m2) live stem N transfer                       
   		 deadstemn_xfer                      =>    cnveg_nitrogenstate_inst%deadstemn_xfer_patch                        , & ! InOut:  [real(r8) (:)]  (gN/m2) dead stem N transfer                       
   		 livecrootn_xfer                     =>    cnveg_nitrogenstate_inst%livecrootn_xfer_patch                       , & ! InOut:  [real(r8) (:)]  (gN/m2) live coarse root N transfer                
   		 deadcrootn_xfer                     =>    cnveg_nitrogenstate_inst%deadcrootn_xfer_patch                       , & ! InOut:  [real(r8) (:)]  (gN/m2) dead coarse root N transfer     
   		                 
   		 leafc_storage_to_xfer               =>    cnveg_carbonflux_inst%leafc_storage_to_xfer_patch                    , & ! InOut:  [real(r8) (:)]                                                     
   		 frootc_storage_to_xfer              =>    cnveg_carbonflux_inst%frootc_storage_to_xfer_patch                   , & ! InOut:  [real(r8) (:)]                                                     
   		 livestemc_storage_to_xfer           =>    cnveg_carbonflux_inst%livestemc_storage_to_xfer_patch                , & ! InOut:  [real(r8) (:)]                                                     
   		 deadstemc_storage_to_xfer           =>    cnveg_carbonflux_inst%deadstemc_storage_to_xfer_patch                , & ! InOut:  [real(r8) (:)]                                                     
   		 livecrootc_storage_to_xfer          =>    cnveg_carbonflux_inst%livecrootc_storage_to_xfer_patch               , & ! InOut:  [real(r8) (:)]                                                     
   		 deadcrootc_storage_to_xfer          =>    cnveg_carbonflux_inst%deadcrootc_storage_to_xfer_patch               , & ! InOut:  [real(r8) (:)]                                                     
   		 gresp_storage_to_xfer               =>    cnveg_carbonflux_inst%gresp_storage_to_xfer_patch                    , & ! InOut:  [real(r8) (:)]  
   		 leafc_xfer_to_leafc                 =>    cnveg_carbonflux_inst%leafc_xfer_to_leafc_patch                      , & ! InOut:  [real(r8) (:)]                                                    
   		 frootc_xfer_to_frootc               =>    cnveg_carbonflux_inst%frootc_xfer_to_frootc_patch                    , & ! InOut:  [real(r8) (:)]                                                    
   		 livestemc_xfer_to_livestemc         =>    cnveg_carbonflux_inst%livestemc_xfer_to_livestemc_patch              , & ! InOut:  [real(r8) (:)]                                                    
   		 deadstemc_xfer_to_deadstemc         =>    cnveg_carbonflux_inst%deadstemc_xfer_to_deadstemc_patch              , & ! InOut:  [real(r8) (:)]                                                    
   		 livecrootc_xfer_to_livecrootc       =>    cnveg_carbonflux_inst%livecrootc_xfer_to_livecrootc_patch            , & ! InOut:  [real(r8) (:)]                                                    
   		 deadcrootc_xfer_to_deadcrootc       =>    cnveg_carbonflux_inst%deadcrootc_xfer_to_deadcrootc_patch            , & ! InOut:  [real(r8) (:)]   
   		                                                    
   		 leafn_storage_to_xfer               =>    cnveg_nitrogenflux_inst%leafn_storage_to_xfer_patch                  , & ! InOut:  [real(r8) (:)]                                                     
   		 frootn_storage_to_xfer              =>    cnveg_nitrogenflux_inst%frootn_storage_to_xfer_patch                 , & ! InOut:  [real(r8) (:)]                                                     
   		 livestemn_storage_to_xfer           =>    cnveg_nitrogenflux_inst%livestemn_storage_to_xfer_patch              , & ! InOut:  [real(r8) (:)]                                                     
   		 deadstemn_storage_to_xfer           =>    cnveg_nitrogenflux_inst%deadstemn_storage_to_xfer_patch              , & ! InOut:  [real(r8) (:)]                                                     
   		 livecrootn_storage_to_xfer          =>    cnveg_nitrogenflux_inst%livecrootn_storage_to_xfer_patch             , & ! InOut:  [real(r8) (:)]   
   		 deadcrootn_storage_to_xfer          =>    cnveg_nitrogenflux_inst%deadcrootn_storage_to_xfer_patch             , & ! InOut:  [real(r8) (:)]                                                     
   		 leafn_xfer_to_leafn                 =>    cnveg_nitrogenflux_inst%leafn_xfer_to_leafn_patch                    , & ! InOut:  [real(r8) (:)]                                                    
   		 frootn_xfer_to_frootn               =>    cnveg_nitrogenflux_inst%frootn_xfer_to_frootn_patch                  , & ! InOut:  [real(r8) (:)]                                                    
   		 livestemn_xfer_to_livestemn         =>    cnveg_nitrogenflux_inst%livestemn_xfer_to_livestemn_patch            , & ! InOut:  [real(r8) (:)]                                                    
   		 deadstemn_xfer_to_deadstemn         =>    cnveg_nitrogenflux_inst%deadstemn_xfer_to_deadstemn_patch            , & ! InOut:  [real(r8) (:)]                                                    
   		 livecrootn_xfer_to_livecrootn       =>    cnveg_nitrogenflux_inst%livecrootn_xfer_to_livecrootn_patch          , & ! InOut:  [real(r8) (:)]                                                    
   		 deadcrootn_xfer_to_deadcrootn       =>    cnveg_nitrogenflux_inst%deadcrootn_xfer_to_deadcrootn_patch          , & ! InOut:  [real(r8) (:)]     
   		           
   		 bglfr                               => cnveg_state_inst%bglfr_patch , & ! Output: [real(r8) (:) ]  background litterfall rate (1/s)                  
   		 bgtr                                => cnveg_state_inst%bgtr_patch  , & ! Output: [real(r8) (:) ]  background transfer growth rate (1/s)             
   		 lgsf                                => cnveg_state_inst%lgsf_patch  , & ! Output: [real(r8) (:) ]  long growing season factor [0-1]                  

   		 ileafst_to_ileafxf_phc              =>    cnveg_carbonflux_inst%ileafst_to_ileafxf_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from leaf storage pool to leaf transfer pool
   		 ileafxf_to_ileaf_phc                =>    cnveg_carbonflux_inst%ileafxf_to_ileaf_ph                   , & ! Input: [integer (:)] Index of phenology related C transfer from leaf transfer pool to leaf pool
   		 ifrootst_to_ifrootxf_phc            =>    cnveg_carbonflux_inst%ifrootst_to_ifrootxf_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from fine root storage pool to fine root transfer pool
   		 ifrootxf_to_ifroot_phc              =>    cnveg_carbonflux_inst%ifrootxf_to_ifroot_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from fine root transfer pool to fine root pool
   		 ilivestemst_to_ilivestemxf_phc      =>    cnveg_carbonflux_inst%ilivestemst_to_ilivestemxf_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live stem storage pool to live stem transfer pool
   		 ilivestemxf_to_ilivestem_phc        =>    cnveg_carbonflux_inst%ilivestemxf_to_ilivestem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem transfer pool to live stem pool
   		 ideadstemst_to_ideadstemxf_phc      =>    cnveg_carbonflux_inst%ideadstemst_to_ideadstemxf_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem storage pool to dead stem transfer pool
   		 ideadstemxf_to_ideadstem_phc        =>    cnveg_carbonflux_inst%ideadstemxf_to_ideadstem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem transfer pool to dead stem pool
   		 ilivecrootst_to_ilivecrootxf_phc    =>    cnveg_carbonflux_inst%ilivecrootst_to_ilivecrootxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root storage pool to live coarse root transfer pool
   		 ilivecrootxf_to_ilivecroot_phc      =>    cnveg_carbonflux_inst%ilivecrootxf_to_ilivecroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root transfer pool to live coarse root pool
   		 ideadcrootst_to_ideadcrootxf_phc    =>    cnveg_carbonflux_inst%ideadcrootst_to_ideadcrootxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root storage pool to dead coarse root transfer pool
   		 ideadcrootxf_to_ideadcroot_phc      =>    cnveg_carbonflux_inst%ideadcrootxf_to_ideadcroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root transfer pool to dead coarse root pool
   		 ilivestem_to_ideadstem_phc          =>    cnveg_carbonflux_inst%ilivestem_to_ideadstem_ph             , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to dead stem pool
   		 ilivecroot_to_ideadcroot_phc        =>    cnveg_carbonflux_inst%ilivecroot_to_ideadcroot_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root to dead coarse root pool
   		 ileaf_to_iout_phc                   =>    cnveg_carbonflux_inst%ileaf_to_iout_ph                      , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to outside of vegetation pools
   		 ifroot_to_iout_phc                  =>    cnveg_carbonflux_inst%ifroot_to_iout_ph                     , & ! Input: [integer (:)] Index of phenology related C transfer from fine root pool to outside of vegetation pools
   		 ilivestem_to_iout_phc               =>    cnveg_carbonflux_inst%ilivestem_to_iout_ph                  , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to outside of vegetation pools
   		 igrain_to_iout_phc                  =>    cnveg_carbonflux_inst%igrain_to_iout_ph                     , & ! Input: [integer (:)] Index of phenology related C transfer from grain pool to outside of vegetation pools
   		 ileafst_to_ileafxf_phn              =>    cnveg_nitrogenflux_inst%ileafst_to_ileafxf_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from leaf storage pool to leaf transfer pool
   		 ileafxf_to_ileaf_phn                =>    cnveg_nitrogenflux_inst%ileafxf_to_ileaf_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from leaf transfer pool to leaf pool
   		 ifrootst_to_ifrootxf_phn            =>    cnveg_nitrogenflux_inst%ifrootst_to_ifrootxf_ph             , & ! Input: [integer (:)] Index of phenology related C transfer from fine root storage pool to fine root transfer pool
   		 ifrootxf_to_ifroot_phn              =>    cnveg_nitrogenflux_inst%ifrootxf_to_ifroot_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from fine root transfer pool to fine root pool
   		 ilivestemst_to_ilivestemxf_phn      =>    cnveg_nitrogenflux_inst%ilivestemst_to_ilivestemxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live stem storage pool to live stem transfer pool
   		 ilivestemxf_to_ilivestem_phn        =>    cnveg_nitrogenflux_inst%ilivestemxf_to_ilivestem_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live stem transfer pool to live stem pool
   		 ideadstemst_to_ideadstemxf_phn      =>    cnveg_nitrogenflux_inst%ideadstemst_to_ideadstemxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem storage pool to dead stem transfer pool
   		 ideadstemxf_to_ideadstem_phn        =>    cnveg_nitrogenflux_inst%ideadstemxf_to_ideadstem_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem transfer pool to dead stem pool
   		 ilivecrootst_to_ilivecrootxf_phn    =>    cnveg_nitrogenflux_inst%ilivecrootst_to_ilivecrootxf_ph     , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root storage pool to live coarse root transfer pool
   		 ilivecrootxf_to_ilivecroot_phn      =>    cnveg_nitrogenflux_inst%ilivecrootxf_to_ilivecroot_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root transfer pool to live coarse root pool
   		 ideadcrootst_to_ideadcrootxf_phn    =>    cnveg_nitrogenflux_inst%ideadcrootst_to_ideadcrootxf_ph     , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root storage pool to dead coarse root transfer pool
   		 ideadcrootxf_to_ideadcroot_phn      =>    cnveg_nitrogenflux_inst%ideadcrootxf_to_ideadcroot_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root transfer pool to dead coarse root pool
   		 ilivestem_to_ideadstem_phn          =>    cnveg_nitrogenflux_inst%ilivestem_to_ideadstem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to dead stem pool
   		 ilivecroot_to_ideadcroot_phn        =>    cnveg_nitrogenflux_inst%ilivecroot_to_ideadcroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root to dead coarse root pool
   		 ileaf_to_iout_phn                   =>    cnveg_nitrogenflux_inst%ileaf_to_iout_ph                    , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to outside of vegetation pools
   		 ifroot_to_iout_phn                  =>    cnveg_nitrogenflux_inst%ifroot_to_iout_ph                   , & ! Input: [integer (:)] Index of phenology related C transfer from fine root pool to outside of vegetation pools
   		 ilivestem_to_iout_phn               =>    cnveg_nitrogenflux_inst%ilivestem_to_iout_ph                , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to outside of vegetation pools
   		 ileaf_to_iretransn_phn              =>    cnveg_nitrogenflux_inst%ileaf_to_iretransn_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to retranslocation pool
   		 ilivestem_to_iretransn_phn          =>    cnveg_nitrogenflux_inst%ilivestem_to_iretransn_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to retranslocation pool
   		 ilivecroot_to_iretransn_phn         =>    cnveg_nitrogenflux_inst%ilivecroot_to_iretransn_ph          , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root pool to retranslocation pool
   		 igrain_to_iout_phn                  =>    cnveg_nitrogenflux_inst%igrain_to_iout_ph                     & ! Input: [integer (:)] Index of phenology related C transfer from grain pool to outside of vegetation pools
         )

      avg_dayspyr = get_average_days_per_year()

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         if (evergreen(ivt(p)) == 1._r8) then
            bglfr(p) = 1._r8/(leaf_long(ivt(p)) * avg_dayspyr * secspday)
            bgtr(p)  = 0._r8
            lgsf(p)  = 0._r8
         end if
      end do
               
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (CN_evergreen_phenology_opt == 1) then    
   do fp = 1,num_soilp    
      p = filter_soilp(fp)    
      if (evergreen(ivt(p)) == 1._r8) then    
   
         tranr=0.0002_r8   
         ! set carbon fluxes for shifting storage pools to transfer pools    
         if (use_matrixcn) then    
            leafc_storage_to_xfer(p)  = leafc_storage(p)  * matrix_update_phc(p,ileafst_to_ileafxf_phc,tranr/dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
            frootc_storage_to_xfer(p) = frootc_storage(p) * matrix_update_phc(p,ifrootst_to_ifrootxf_phc,tranr/dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
            if (woody(ivt(p)) == 1.0_r8) then   
               livestemc_storage_to_xfer(p)   = livestemc_storage(p)  * matrix_update_phc(p,ilivestemst_to_ilivestemxf_phc,tranr/dt ,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
               deadstemc_storage_to_xfer(p)   = deadstemc_storage(p)  * matrix_update_phc(p,ideadstemst_to_ideadstemxf_phc,tranr/dt ,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
               livecrootc_storage_to_xfer(p)  = livecrootc_storage(p) * matrix_update_phc(p,ilivecrootst_to_ilivecrootxf_phc,tranr/dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
               deadcrootc_storage_to_xfer(p)  = deadcrootc_storage(p) * matrix_update_phc(p,ideadcrootst_to_ideadcrootxf_phc,tranr/dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
            end if
         else
            ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
            leafc_storage_to_xfer(p)  = tranr * leafc_storage(p)/dt    
            frootc_storage_to_xfer(p) = tranr * frootc_storage(p)/dt
            if (woody(ivt(p)) == 1.0_r8) then    
               livestemc_storage_to_xfer(p)  = tranr * livestemc_storage(p)/dt    
               deadstemc_storage_to_xfer(p)  = tranr * deadstemc_storage(p)/dt    
               livecrootc_storage_to_xfer(p) = tranr * livecrootc_storage(p)/dt   
               deadcrootc_storage_to_xfer(p) = tranr * deadcrootc_storage(p)/dt   
               gresp_storage_to_xfer(p)      = tranr * gresp_storage(p)/dt        
            end if    
         end if !use_matrixcn

        ! set nitrogen fluxes for shifting storage pools to transfer pools    
        if (use_matrixcn) then    
           leafn_storage_to_xfer(p)  = leafn_storage(p)  * matrix_update_phn(p,ileafst_to_ileafxf_phn,tranr/dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
           frootn_storage_to_xfer(p) = frootn_storage(p) * matrix_update_phn(p,ifrootst_to_ifrootxf_phn,tranr/dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
           if (woody(ivt(p)) == 1.0_r8) then   
              livestemn_storage_to_xfer(p)  = livestemn_storage(p)  * matrix_update_phn(p,ilivestemst_to_ilivestemxf_phn,tranr/dt ,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
              deadstemn_storage_to_xfer(p)  = deadstemn_storage(p)  * matrix_update_phn(p,ideadstemst_to_ideadstemxf_phn,tranr/dt ,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
              livecrootn_storage_to_xfer(p) = livecrootn_storage(p) * matrix_update_phn(p,ilivecrootst_to_ilivecrootxf_phn,tranr/dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
              deadcrootn_storage_to_xfer(p) = deadcrootn_storage(p) * matrix_update_phn(p,ideadcrootst_to_ideadcrootxf_phn,tranr/dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
           end if
        else
           ! NOTE: The non matrix version of this is in CNNStateUpdate1::NStateUpdate1 EBK (11/26/2019)
           leafn_storage_to_xfer(p)  = tranr * leafn_storage(p)/dt    
           frootn_storage_to_xfer(p) = tranr * frootn_storage(p)/dt   
           if (woody(ivt(p)) == 1.0_r8) then    
               livestemn_storage_to_xfer(p)  = tranr * livestemn_storage(p)/dt    
               deadstemn_storage_to_xfer(p)  = tranr * deadstemn_storage(p)/dt    
               livecrootn_storage_to_xfer(p) = tranr * livecrootn_storage(p)/dt   
               deadcrootn_storage_to_xfer(p) = tranr * deadcrootn_storage(p)/dt   
           end if    
        end if !use_matrixcn  
                        
        t1 = 1.0_r8 / dt   

        if (use_matrixcn) then
           leafc_xfer_to_leafc(p)   = leafc_xfer(p)  * matrix_update_phc(p,ileafxf_to_ileaf_phc,t1,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
           frootc_xfer_to_frootc(p) = frootc_xfer(p) * matrix_update_phc(p,ifrootxf_to_ifroot_phc,t1,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)

           leafn_xfer_to_leafn(p)   = leafn_xfer(p)  * matrix_update_phn(p,ileafxf_to_ileaf_phn,t1,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
           frootn_xfer_to_frootn(p) = frootn_xfer(p) * matrix_update_phn(p,ifrootxf_to_ifroot_phn,t1,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
           if (woody(ivt(p)) == 1.0_r8) then
              livestemc_xfer_to_livestemc(p)   = livestemc_xfer(p)  * matrix_update_phc(p,ilivestemxf_to_ilivestem_phc,t1,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
              deadstemc_xfer_to_deadstemc(p)   = deadstemc_xfer(p)  * matrix_update_phc(p,ideadstemxf_to_ideadstem_phc,t1,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
              livecrootc_xfer_to_livecrootc(p) = livecrootc_xfer(p) * matrix_update_phc(p,ilivecrootxf_to_ilivecroot_phc,t1,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
              deadcrootc_xfer_to_deadcrootc(p) = deadcrootc_xfer(p) * matrix_update_phc(p,ideadcrootxf_to_ideadcroot_phc,t1,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)

              livestemn_xfer_to_livestemn(p)   = livestemn_xfer(p)  * matrix_update_phn(p,ilivestemxf_to_ilivestem_phn,t1,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
              deadstemn_xfer_to_deadstemn(p)   = deadstemn_xfer(p)  * matrix_update_phn(p,ideadstemxf_to_ideadstem_phn,t1,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
              livecrootn_xfer_to_livecrootn(p) = livecrootn_xfer(p) * matrix_update_phn(p,ilivecrootxf_to_ilivecroot_phn,t1,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
              deadcrootn_xfer_to_deadcrootn(p) = deadcrootn_xfer(p) * matrix_update_phn(p,ideadcrootxf_to_ideadcroot_phn,t1,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
           end if
        else
           ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
           !                                        and CNNStateUpdate1::NStateUpdate1
           leafc_xfer_to_leafc(p)   = t1 * leafc_xfer(p)    
           frootc_xfer_to_frootc(p) = t1 * frootc_xfer(p)   
            
           leafn_xfer_to_leafn(p)   = t1 * leafn_xfer(p)    
           frootn_xfer_to_frootn(p) = t1 * frootn_xfer(p)   
           if (woody(ivt(p)) == 1.0_r8) then   
               livestemc_xfer_to_livestemc(p)   = t1 * livestemc_xfer(p)   
               deadstemc_xfer_to_deadstemc(p)   = t1 * deadstemc_xfer(p)   
               livecrootc_xfer_to_livecrootc(p) = t1 * livecrootc_xfer(p)  
               deadcrootc_xfer_to_deadcrootc(p) = t1 * deadcrootc_xfer(p)  
                
               livestemn_xfer_to_livestemn(p)   = t1 * livestemn_xfer(p)   
               deadstemn_xfer_to_deadstemn(p)   = t1 * deadstemn_xfer(p)   
               livecrootn_xfer_to_livecrootn(p) = t1 * livecrootn_xfer(p)  
               deadcrootn_xfer_to_deadcrootn(p) = t1 * deadcrootn_xfer(p)  
           end if
        end if !use_matrixcn
                
      end if ! end of if (evergreen(ivt(p)) == 1._r8) then    
     
   end do ! end of pft loop 
   
   end if ! end of if (CN_evergreen_phenology_opt == 1) then    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end associate

  end subroutine CNEvergreenPhenology

  !-----------------------------------------------------------------------
  subroutine CNSeasonDecidPhenology (num_soilp, filter_soilp       , &
       temperature_inst, waterdiagnosticbulk_inst, cnveg_state_inst, dgvs_inst , &
       cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! For coupled carbon-nitrogen code (CN).
    ! This routine handles the seasonal deciduous phenology code (temperate
    ! deciduous vegetation that has only one growing season per year).
    !
    ! !USES:
    use shr_const_mod   , only: SHR_CONST_TKFRZ, SHR_CONST_PI
    use clm_varcon      , only: secspday
    use clm_varctl      , only: use_cndv
    !
    ! !ARGUMENTS:
    integer                        , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                        , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(temperature_type)         , intent(in)    :: temperature_inst
    type(waterdiagnosticbulk_type)          , intent(in)    :: waterdiagnosticbulk_inst
    type(cnveg_state_type)         , intent(inout) :: cnveg_state_inst
    type(dgvs_type)                , intent(inout) :: dgvs_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type) , intent(inout) :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: g,c,p          !indices
    integer :: fp             !filter patch index
    real(r8):: ws_flag        !winter-summer solstice flag (0 or 1)
    real(r8):: crit_onset_gdd !critical onset growing degree-day sum
    real(r8):: crit_daylat    !latitudinal light gradient in arctic-boreal 
    logical :: do_onset       ! Flag if onset should happen
    real(r8):: soilt          ! soil temperature for layer to use for seasonal phenology trigger
    !-----------------------------------------------------------------------

    associate(                                                                                                   & 
         ivt                                 =>    patch%itype                                                   , & ! Input:  [integer   (:)   ]  patch vegetation type                                
         dayl                                =>    grc%dayl                                                    , & ! Input:  [real(r8)  (:)   ]  daylength (s)
         prev_dayl                           =>    grc%prev_dayl                                               , & ! Input:  [real(r8)  (:)   ]  daylength from previous time step (s)
         
         woody                               =>    pftcon%woody                                                , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)
         season_decid                        =>    pftcon%season_decid                                         , & ! Input:  binary flag for seasonal-deciduous leaf habit (0 or 1)
         season_decid_temperate              =>    pftcon%season_decid_temperate                               , & ! Input:  binary flag for seasonal-deciduous temperate leaf habit (0 or 1)
         crit_onset_gdd_sf                   =>    pftcon%crit_onset_gdd_sf                                    , & ! Input:  scale factor for crit_onset_gdd (unitless)
         ndays_on                            =>    pftcon%ndays_on                                             , & ! Input:  number of days to complete leaf onset (days
         
         t_soisno                            =>    temperature_inst%t_soisno_col                               , & ! Input:  [real(r8)  (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
         soila10                             =>    temperature_inst%soila10_col                                , & ! Input:  [real(r8) (:)   ] 
         t_a5min                            =>    temperature_inst%t_a5min_patch                             , & ! input:  [real(r8) (:)   ]
         snow_5day                           =>    waterdiagnosticbulk_inst%snow_5day_col                      , & ! input:  [real(r8) (:)   ] 

         pftmayexist                         =>    dgvs_inst%pftmayexist_patch                                 , & ! Output: [logical   (:)   ]  exclude seasonal decid patches from tropics           

         annavg_t2m                          =>    cnveg_state_inst%annavg_t2m_patch                           , & ! Input:  [real(r8)  (:)   ]  annual average 2m air temperature (K)             
         dormant_flag                        =>    cnveg_state_inst%dormant_flag_patch                         , & ! Output: [real(r8)  (:)   ]  dormancy flag                                     
         days_active                         =>    cnveg_state_inst%days_active_patch                          , & ! Output: [real(r8)  (:)   ]  number of days since last dormancy                
         onset_flag                          =>    cnveg_state_inst%onset_flag_patch                           , & ! Output: [real(r8)  (:)   ]  onset flag                                        
         onset_counter                       =>    cnveg_state_inst%onset_counter_patch                        , & ! Output: [real(r8)  (:)   ]  onset counter (seconds)                           
         onset_gddflag                       =>    cnveg_state_inst%onset_gddflag_patch                        , & ! Output: [real(r8)  (:)   ]  onset freeze flag                                 
         onset_gdd                           =>    cnveg_state_inst%onset_gdd_patch                            , & ! Output: [real(r8)  (:)   ]  onset growing degree days                         
         offset_flag                         =>    cnveg_state_inst%offset_flag_patch                          , & ! Output: [real(r8)  (:)   ]  offset flag                                       
         offset_counter                      =>    cnveg_state_inst%offset_counter_patch                       , & ! Output: [real(r8)  (:)   ]  offset counter (seconds)                          
         bglfr                               =>    cnveg_state_inst%bglfr_patch                                , & ! Output: [real(r8)  (:)   ]  background litterfall rate (1/s)                  
         bgtr                                =>    cnveg_state_inst%bgtr_patch                                 , & ! Output: [real(r8)  (:)   ]  background transfer growth rate (1/s)             
         lgsf                                =>    cnveg_state_inst%lgsf_patch                                 , & ! Output: [real(r8)  (:)   ]  long growing season factor [0-1]                  
         
         leafc_storage                       =>    cnveg_carbonstate_inst%leafc_storage_patch                  , & ! Input:  [real(r8)  (:)   ]  (gC/m2) leaf C storage                            
         frootc_storage                      =>    cnveg_carbonstate_inst%frootc_storage_patch                 , & ! Input:  [real(r8)  (:)   ]  (gC/m2) fine root C storage                       
         livestemc_storage                   =>    cnveg_carbonstate_inst%livestemc_storage_patch              , & ! Input:  [real(r8)  (:)   ]  (gC/m2) live stem C storage                       
         deadstemc_storage                   =>    cnveg_carbonstate_inst%deadstemc_storage_patch              , & ! Input:  [real(r8)  (:)   ]  (gC/m2) dead stem C storage                       
         livecrootc_storage                  =>    cnveg_carbonstate_inst%livecrootc_storage_patch             , & ! Input:  [real(r8)  (:)   ]  (gC/m2) live coarse root C storage                
         deadcrootc_storage                  =>    cnveg_carbonstate_inst%deadcrootc_storage_patch             , & ! Input:  [real(r8)  (:)   ]  (gC/m2) dead coarse root C storage                
         gresp_storage                       =>    cnveg_carbonstate_inst%gresp_storage_patch                  , & ! Input:  [real(r8)  (:)   ]  (gC/m2) growth respiration storage                
         leafc_xfer                          =>    cnveg_carbonstate_inst%leafc_xfer_patch                     , & ! Output:  [real(r8) (:)   ]  (gC/m2) leaf C transfer                           
         frootc_xfer                         =>    cnveg_carbonstate_inst%frootc_xfer_patch                    , & ! Output:  [real(r8) (:)   ]  (gC/m2) fine root C transfer                      
         livestemc_xfer                      =>    cnveg_carbonstate_inst%livestemc_xfer_patch                 , & ! Output:  [real(r8) (:)   ]  (gC/m2) live stem C transfer                      
         deadstemc_xfer                      =>    cnveg_carbonstate_inst%deadstemc_xfer_patch                 , & ! Output:  [real(r8) (:)   ]  (gC/m2) dead stem C transfer                      
         livecrootc_xfer                     =>    cnveg_carbonstate_inst%livecrootc_xfer_patch                , & ! Output:  [real(r8) (:)   ]  (gC/m2) live coarse root C transfer               
         deadcrootc_xfer                     =>    cnveg_carbonstate_inst%deadcrootc_xfer_patch                , & ! Output:  [real(r8) (:)   ]  (gC/m2) dead coarse root C transfer               
         
         leafn_storage                       =>    cnveg_nitrogenstate_inst%leafn_storage_patch                , & ! Input:  [real(r8)  (:)   ]  (gN/m2) leaf N storage                            
         frootn_storage                      =>    cnveg_nitrogenstate_inst%frootn_storage_patch               , & ! Input:  [real(r8)  (:)   ]  (gN/m2) fine root N storage                       
         livestemn_storage                   =>    cnveg_nitrogenstate_inst%livestemn_storage_patch            , & ! Input:  [real(r8)  (:)   ]  (gN/m2) live stem N storage                       
         deadstemn_storage                   =>    cnveg_nitrogenstate_inst%deadstemn_storage_patch            , & ! Input:  [real(r8)  (:)   ]  (gN/m2) dead stem N storage                       
         livecrootn_storage                  =>    cnveg_nitrogenstate_inst%livecrootn_storage_patch           , & ! Input:  [real(r8)  (:)   ]  (gN/m2) live coarse root N storage                
         deadcrootn_storage                  =>    cnveg_nitrogenstate_inst%deadcrootn_storage_patch           , & ! Input:  [real(r8)  (:)   ]  (gN/m2) dead coarse root N storage                
         leafn_xfer                          =>    cnveg_nitrogenstate_inst%leafn_xfer_patch                   , & ! Output:  [real(r8) (:)   ]  (gN/m2) leaf N transfer                           
         frootn_xfer                         =>    cnveg_nitrogenstate_inst%frootn_xfer_patch                  , & ! Output:  [real(r8) (:)   ]  (gN/m2) fine root N transfer                      
         livestemn_xfer                      =>    cnveg_nitrogenstate_inst%livestemn_xfer_patch               , & ! Output:  [real(r8) (:)   ]  (gN/m2) live stem N transfer                      
         deadstemn_xfer                      =>    cnveg_nitrogenstate_inst%deadstemn_xfer_patch               , & ! Output:  [real(r8) (:)   ]  (gN/m2) dead stem N transfer                      
         livecrootn_xfer                     =>    cnveg_nitrogenstate_inst%livecrootn_xfer_patch              , & ! Output:  [real(r8) (:)   ]  (gN/m2) live coarse root N transfer               
         deadcrootn_xfer                     =>    cnveg_nitrogenstate_inst%deadcrootn_xfer_patch              , & ! Output:  [real(r8) (:)   ]  (gN/m2) dead coarse root N transfer               

         prev_leafc_to_litter                =>    cnveg_carbonflux_inst%prev_leafc_to_litter_patch            , & ! Output: [real(r8)  (:)   ]  previous timestep leaf C litterfall flux (gC/m2/s)
         prev_frootc_to_litter               =>    cnveg_carbonflux_inst%prev_frootc_to_litter_patch           , & ! Output: [real(r8)  (:)   ]  previous timestep froot C litterfall flux (gC/m2/s)
         leafc_xfer_to_leafc                 =>    cnveg_carbonflux_inst%leafc_xfer_to_leafc_patch             , & ! Output:  [real(r8) (:)   ]                                                    
         frootc_xfer_to_frootc               =>    cnveg_carbonflux_inst%frootc_xfer_to_frootc_patch           , & ! Output:  [real(r8) (:)   ]                                                    
         livestemc_xfer_to_livestemc         =>    cnveg_carbonflux_inst%livestemc_xfer_to_livestemc_patch     , & ! Output:  [real(r8) (:)   ]                                                    
         deadstemc_xfer_to_deadstemc         =>    cnveg_carbonflux_inst%deadstemc_xfer_to_deadstemc_patch     , & ! Output:  [real(r8) (:)   ]                                                    
         livecrootc_xfer_to_livecrootc       =>    cnveg_carbonflux_inst%livecrootc_xfer_to_livecrootc_patch   , & ! Output:  [real(r8) (:)   ]                                                    
         deadcrootc_xfer_to_deadcrootc       =>    cnveg_carbonflux_inst%deadcrootc_xfer_to_deadcrootc_patch   , & ! Output:  [real(r8) (:)   ]                                                    
         leafc_storage_to_xfer               =>    cnveg_carbonflux_inst%leafc_storage_to_xfer_patch           , & ! Output:  [real(r8) (:)   ]                                                    
         frootc_storage_to_xfer              =>    cnveg_carbonflux_inst%frootc_storage_to_xfer_patch          , & ! Output:  [real(r8) (:)   ]                                                    
         livestemc_storage_to_xfer           =>    cnveg_carbonflux_inst%livestemc_storage_to_xfer_patch       , & ! Output:  [real(r8) (:)   ]                                                    
         deadstemc_storage_to_xfer           =>    cnveg_carbonflux_inst%deadstemc_storage_to_xfer_patch       , & ! Output:  [real(r8) (:)   ]                                                    
         livecrootc_storage_to_xfer          =>    cnveg_carbonflux_inst%livecrootc_storage_to_xfer_patch      , & ! Output:  [real(r8) (:)   ]                                                    
         deadcrootc_storage_to_xfer          =>    cnveg_carbonflux_inst%deadcrootc_storage_to_xfer_patch      , & ! Output:  [real(r8) (:)   ]                                                    
         gresp_storage_to_xfer               =>    cnveg_carbonflux_inst%gresp_storage_to_xfer_patch           , & ! Output:  [real(r8) (:)   ]                                                    
         
         leafn_xfer_to_leafn                 =>    cnveg_nitrogenflux_inst%leafn_xfer_to_leafn_patch           , & ! Output:  [real(r8) (:)   ]                                                    
         frootn_xfer_to_frootn               =>    cnveg_nitrogenflux_inst%frootn_xfer_to_frootn_patch         , & ! Output:  [real(r8) (:)   ]                                                    
         livestemn_xfer_to_livestemn         =>    cnveg_nitrogenflux_inst%livestemn_xfer_to_livestemn_patch   , & ! Output:  [real(r8) (:)   ]                                                    
         deadstemn_xfer_to_deadstemn         =>    cnveg_nitrogenflux_inst%deadstemn_xfer_to_deadstemn_patch   , & ! Output:  [real(r8) (:)   ]                                                    
         livecrootn_xfer_to_livecrootn       =>    cnveg_nitrogenflux_inst%livecrootn_xfer_to_livecrootn_patch , & ! Output:  [real(r8) (:)   ]                                                    
         deadcrootn_xfer_to_deadcrootn       =>    cnveg_nitrogenflux_inst%deadcrootn_xfer_to_deadcrootn_patch , & ! Output:  [real(r8) (:)   ]                                                    
         leafn_storage_to_xfer               =>    cnveg_nitrogenflux_inst%leafn_storage_to_xfer_patch         , & ! Output:  [real(r8) (:)   ]                                                    
         frootn_storage_to_xfer              =>    cnveg_nitrogenflux_inst%frootn_storage_to_xfer_patch        , & ! Output:  [real(r8) (:)   ]                                                    
         livestemn_storage_to_xfer           =>    cnveg_nitrogenflux_inst%livestemn_storage_to_xfer_patch     , & ! Output:  [real(r8) (:)   ]                                                    
         deadstemn_storage_to_xfer           =>    cnveg_nitrogenflux_inst%deadstemn_storage_to_xfer_patch     , & ! Output:  [real(r8) (:)   ]                                                    
         livecrootn_storage_to_xfer          =>    cnveg_nitrogenflux_inst%livecrootn_storage_to_xfer_patch    , & ! Output:  [real(r8) (:)   ]                                                    
         deadcrootn_storage_to_xfer          =>    cnveg_nitrogenflux_inst%deadcrootn_storage_to_xfer_patch    , & ! Output:  [real(r8) (:)   ]   
         ileafst_to_ileafxf_phc              =>    cnveg_carbonflux_inst%ileafst_to_ileafxf_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from leaf storage pool to leaf transfer pool
         ileafxf_to_ileaf_phc                =>    cnveg_carbonflux_inst%ileafxf_to_ileaf_ph                   , & ! Input: [integer (:)] Index of phenology related C transfer from leaf transfer pool to leaf pool
         ifrootst_to_ifrootxf_phc            =>    cnveg_carbonflux_inst%ifrootst_to_ifrootxf_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from fine root storage pool to fine root transfer pool
         ifrootxf_to_ifroot_phc              =>    cnveg_carbonflux_inst%ifrootxf_to_ifroot_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from fine root transfer pool to fine root pool
         ilivestemst_to_ilivestemxf_phc      =>    cnveg_carbonflux_inst%ilivestemst_to_ilivestemxf_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live stem storage pool to live stem transfer pool
         ilivestemxf_to_ilivestem_phc        =>    cnveg_carbonflux_inst%ilivestemxf_to_ilivestem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem transfer pool to live stem pool
         ideadstemst_to_ideadstemxf_phc      =>    cnveg_carbonflux_inst%ideadstemst_to_ideadstemxf_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem storage pool to dead stem transfer pool
         ideadstemxf_to_ideadstem_phc        =>    cnveg_carbonflux_inst%ideadstemxf_to_ideadstem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem transfer pool to dead stem pool
         ilivecrootst_to_ilivecrootxf_phc    =>    cnveg_carbonflux_inst%ilivecrootst_to_ilivecrootxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root storage pool to live coarse root transfer pool
         ilivecrootxf_to_ilivecroot_phc      =>    cnveg_carbonflux_inst%ilivecrootxf_to_ilivecroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root transfer pool to live coarse root pool
         ideadcrootst_to_ideadcrootxf_phc    =>    cnveg_carbonflux_inst%ideadcrootst_to_ideadcrootxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root storage pool to dead coarse root transfer pool
         ideadcrootxf_to_ideadcroot_phc      =>    cnveg_carbonflux_inst%ideadcrootxf_to_ideadcroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root transfer pool to dead coarse root pool
         ilivestem_to_ideadstem_phc          =>    cnveg_carbonflux_inst%ilivestem_to_ideadstem_ph             , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to dead stem pool
         ilivecroot_to_ideadcroot_phc        =>    cnveg_carbonflux_inst%ilivecroot_to_ideadcroot_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root to dead coarse root pool
         ileaf_to_iout_phc                   =>    cnveg_carbonflux_inst%ileaf_to_iout_ph                      , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to outside of vegetation pools
         ifroot_to_iout_phc                  =>    cnveg_carbonflux_inst%ifroot_to_iout_ph                     , & ! Input: [integer (:)] Index of phenology related C transfer from fine root pool to outside of vegetation pools
         ilivestem_to_iout_phc               =>    cnveg_carbonflux_inst%ilivestem_to_iout_ph                  , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to outside of vegetation pools
         igrain_to_iout_phc                  =>    cnveg_carbonflux_inst%igrain_to_iout_ph                     , & ! Input: [integer (:)] Index of phenology related C transfer from grain pool to outside of vegetation pools
         ileafst_to_ileafxf_phn              =>    cnveg_nitrogenflux_inst%ileafst_to_ileafxf_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from leaf storage pool to leaf transfer pool
         ileafxf_to_ileaf_phn                =>    cnveg_nitrogenflux_inst%ileafxf_to_ileaf_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from leaf transfer pool to leaf pool
         ifrootst_to_ifrootxf_phn            =>    cnveg_nitrogenflux_inst%ifrootst_to_ifrootxf_ph             , & ! Input: [integer (:)] Index of phenology related C transfer from fine root storage pool to fine root transfer pool
         ifrootxf_to_ifroot_phn              =>    cnveg_nitrogenflux_inst%ifrootxf_to_ifroot_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from fine root transfer pool to fine root pool
         ilivestemst_to_ilivestemxf_phn      =>    cnveg_nitrogenflux_inst%ilivestemst_to_ilivestemxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live stem storage pool to live stem transfer pool
         ilivestemxf_to_ilivestem_phn        =>    cnveg_nitrogenflux_inst%ilivestemxf_to_ilivestem_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live stem transfer pool to live stem pool
         ideadstemst_to_ideadstemxf_phn      =>    cnveg_nitrogenflux_inst%ideadstemst_to_ideadstemxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem storage pool to dead stem transfer pool
         ideadstemxf_to_ideadstem_phn        =>    cnveg_nitrogenflux_inst%ideadstemxf_to_ideadstem_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem transfer pool to dead stem pool
         ilivecrootst_to_ilivecrootxf_phn    =>    cnveg_nitrogenflux_inst%ilivecrootst_to_ilivecrootxf_ph     , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root storage pool to live coarse root transfer pool
         ilivecrootxf_to_ilivecroot_phn      =>    cnveg_nitrogenflux_inst%ilivecrootxf_to_ilivecroot_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root transfer pool to live coarse root pool
         ideadcrootst_to_ideadcrootxf_phn    =>    cnveg_nitrogenflux_inst%ideadcrootst_to_ideadcrootxf_ph     , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root storage pool to dead coarse root transfer pool
         ideadcrootxf_to_ideadcroot_phn      =>    cnveg_nitrogenflux_inst%ideadcrootxf_to_ideadcroot_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root transfer pool to dead coarse root pool
         ilivestem_to_ideadstem_phn          =>    cnveg_nitrogenflux_inst%ilivestem_to_ideadstem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to dead stem pool
         ilivecroot_to_ideadcroot_phn        =>    cnveg_nitrogenflux_inst%ilivecroot_to_ideadcroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root to dead coarse root pool
         ileaf_to_iout_phn                   =>    cnveg_nitrogenflux_inst%ileaf_to_iout_ph                    , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to outside of vegetation pools
         ifroot_to_iout_phn                  =>    cnveg_nitrogenflux_inst%ifroot_to_iout_ph                   , & ! Input: [integer (:)] Index of phenology related C transfer from fine root pool to outside of vegetation pools
         ilivestem_to_iout_phn               =>    cnveg_nitrogenflux_inst%ilivestem_to_iout_ph                , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to outside of vegetation pools
         ileaf_to_iretransn_phn              =>    cnveg_nitrogenflux_inst%ileaf_to_iretransn_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to retranslocation pool
         ilivestem_to_iretransn_phn          =>    cnveg_nitrogenflux_inst%ilivestem_to_iretransn_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to retranslocation pool
         ilivecroot_to_iretransn_phn         =>    cnveg_nitrogenflux_inst%ilivecroot_to_iretransn_ph          , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root pool to retranslocation pool
         igrain_to_iout_phn                  =>    cnveg_nitrogenflux_inst%igrain_to_iout_ph                     & ! Input: [integer (:)] Index of phenology related C transfer from grain pool to outside of vegetation pools
         )

      ! start patch loop


      do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = patch%column(p)
         g = patch%gridcell(p)

         if (season_decid(ivt(p)) == 1._r8) then

            ! set background litterfall rate, background transfer rate, and
            ! long growing season factor to 0 for seasonal deciduous types
            bglfr(p) = 0._r8
            bgtr(p) = 0._r8
            lgsf(p) = 0._r8

            ! onset gdd sum from Biome-BGC, v4.1.2
            crit_onset_gdd = crit_onset_gdd_sf(ivt(p)) * exp(4.8_r8 + 0.13_r8*(annavg_t2m(p) &
                             - SHR_CONST_TKFRZ))

            ! set flag for solstice period (winter->summer = 1, summer->winter = 0)
            if (dayl(g) >= prev_dayl(g)) then
               ws_flag = 1._r8
            else
               ws_flag = 0._r8
            end if

            ! update offset_counter and test for the end of the offset period
            if (offset_flag(p) == 1.0_r8) then
               ! decrement counter for offset period
               offset_counter(p) = offset_counter(p) - dt

               ! if this is the end of the offset_period, reset phenology
               ! flags and indices
               if (offset_counter(p) < dt/2._r8) then
                  ! this code block was originally handled by call cn_offset_cleanup(p)
                  ! inlined during vectorization

                  offset_flag(p) = 0._r8
                  offset_counter(p) = 0._r8
                  dormant_flag(p) = 1._r8
                  days_active(p) = 0._r8
                  if (use_cndv) then
                     pftmayexist(p) = .true.
                  end if

                  ! reset the previous timestep litterfall flux memory
                  prev_leafc_to_litter(p) = 0._r8
                  prev_frootc_to_litter(p) = 0._r8
               end if
            end if

            ! update onset_counter and test for the end of the onset period
            if (onset_flag(p) == 1.0_r8) then
               ! decrement counter for onset period
               onset_counter(p) = onset_counter(p) - dt

               ! if this is the end of the onset period, reset phenology
               ! flags and indices
               if (onset_counter(p) < dt/2._r8) then
                  ! this code block was originally handled by call cn_onset_cleanup(p)
                  ! inlined during vectorization

                  onset_flag(p) = 0.0_r8
                  onset_counter(p) = 0.0_r8
                  ! set all transfer growth rates to 0.0
                  leafc_xfer_to_leafc(p)   = 0.0_r8
                  frootc_xfer_to_frootc(p) = 0.0_r8
                  leafn_xfer_to_leafn(p)   = 0.0_r8
                  frootn_xfer_to_frootn(p) = 0.0_r8
                  if (woody(ivt(p)) == 1.0_r8) then
                     livestemc_xfer_to_livestemc(p)   = 0.0_r8
                     deadstemc_xfer_to_deadstemc(p)   = 0.0_r8
                     livecrootc_xfer_to_livecrootc(p) = 0.0_r8
                     deadcrootc_xfer_to_deadcrootc(p) = 0.0_r8
                     livestemn_xfer_to_livestemn(p)   = 0.0_r8
                     deadstemn_xfer_to_deadstemn(p)   = 0.0_r8
                     livecrootn_xfer_to_livecrootn(p) = 0.0_r8
                     deadcrootn_xfer_to_deadcrootn(p) = 0.0_r8
                  end if
                  ! set transfer pools to 0.0
                  leafc_xfer(p) = 0.0_r8
                  leafn_xfer(p) = 0.0_r8
                  frootc_xfer(p) = 0.0_r8
                  frootn_xfer(p) = 0.0_r8
                  if (woody(ivt(p)) == 1.0_r8) then
                     livestemc_xfer(p) = 0.0_r8
                     livestemn_xfer(p) = 0.0_r8
                     deadstemc_xfer(p) = 0.0_r8
                     deadstemn_xfer(p) = 0.0_r8
                     livecrootc_xfer(p) = 0.0_r8
                     livecrootn_xfer(p) = 0.0_r8
                     deadcrootc_xfer(p) = 0.0_r8
                     deadcrootn_xfer(p) = 0.0_r8
                  end if
               end if
            end if

            ! test for switching from dormant period to growth period
            if (dormant_flag(p) == 1.0_r8) then
               soilt = t_soisno(c, phenology_soil_layer)

               do_onset = SeasonalDecidOnset( onset_gdd(p), onset_gddflag(p), soilt, soila10(c), t_a5min(p), dayl(g), &
                                              snow_5day(c), ws_flag, crit_onset_gdd, season_decid_temperate(ivt(p)) )
               ! If onset is being triggered
               if (do_onset) then
                  onset_flag(p) = 1.0_r8
                  dormant_flag(p) = 0.0_r8
                  onset_gddflag(p) = 0.0_r8
                  onset_gdd(p) = 0.0_r8
                  do_onset = .false.
                  onset_counter(p) = ndays_on(ivt(p)) * secspday


                  ! move all the storage pools into transfer pools,
                  ! where they will be transfered to displayed growth over the onset period.
                  ! this code was originally handled with call cn_storage_to_xfer(p)
                  ! inlined during vectorization

                  ! set carbon fluxes for shifting storage pools to transfer pools
                  if(use_matrixcn)then
                     leafc_storage_to_xfer(p)  = leafc_storage(p)  * matrix_update_phc(p,ileafst_to_ileafxf_phc,fstor2tran/dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                     frootc_storage_to_xfer(p) = frootc_storage(p) * matrix_update_phc(p,ifrootst_to_ifrootxf_phc,fstor2tran/dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)

                     if (woody(ivt(p)) == 1.0_r8) then
                       livestemc_storage_to_xfer(p)  = livestemc_storage(p)  * matrix_update_phc(p,ilivestemst_to_ilivestemxf_phc ,fstor2tran/dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                        deadstemc_storage_to_xfer(p)  = deadstemc_storage(p)  * matrix_update_phc(p,ideadstemst_to_ideadstemxf_phc ,fstor2tran/dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                        livecrootc_storage_to_xfer(p) = livecrootc_storage(p) * matrix_update_phc(p,ilivecrootst_to_ilivecrootxf_phc,fstor2tran/dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                        deadcrootc_storage_to_xfer(p) = deadcrootc_storage(p) * matrix_update_phc(p,ideadcrootst_to_ideadcrootxf_phc,fstor2tran/dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                        gresp_storage_to_xfer(p)      = fstor2tran * gresp_storage(p)/dt
                     end if
                     leafn_storage_to_xfer(p)  = leafn_storage(p)  * matrix_update_phn(p,ileafst_to_ileafxf_phn,fstor2tran/dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                     frootn_storage_to_xfer(p) = frootn_storage(p) * matrix_update_phn(p,ifrootst_to_ifrootxf_phn,fstor2tran/dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)

                     if (woody(ivt(p)) == 1.0_r8) then
                        livestemn_storage_to_xfer(p)  = livestemn_storage(p)  * matrix_update_phn(p,ilivestemst_to_ilivestemxf_phn ,fstor2tran/dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                        deadstemn_storage_to_xfer(p)  = deadstemn_storage(p)  * matrix_update_phn(p,ideadstemst_to_ideadstemxf_phn ,fstor2tran/dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                        livecrootn_storage_to_xfer(p) = livecrootn_storage(p) * matrix_update_phn(p,ilivecrootst_to_ilivecrootxf_phn,fstor2tran/dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                        deadcrootn_storage_to_xfer(p) = deadcrootn_storage(p) * matrix_update_phn(p,ideadcrootst_to_ideadcrootxf_phn,fstor2tran/dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                     end if
                  else
                     ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
                     !                                        and CNNStateUpdate1::NStateUpdate1
                     leafc_storage_to_xfer(p)  = fstor2tran * leafc_storage(p)/dt
                     frootc_storage_to_xfer(p) = fstor2tran * frootc_storage(p)/dt
                     if (woody(ivt(p)) == 1.0_r8) then
                        livestemc_storage_to_xfer(p)  = fstor2tran * livestemc_storage(p)/dt
                        deadstemc_storage_to_xfer(p)  = fstor2tran * deadstemc_storage(p)/dt
                        livecrootc_storage_to_xfer(p) = fstor2tran * livecrootc_storage(p)/dt
                        deadcrootc_storage_to_xfer(p) = fstor2tran * deadcrootc_storage(p)/dt
                        gresp_storage_to_xfer(p)      = fstor2tran * gresp_storage(p)/dt
                     end if

                     ! set nitrogen fluxes for shifting storage pools to transfer pools
                     leafn_storage_to_xfer(p)  = fstor2tran * leafn_storage(p)/dt
                     frootn_storage_to_xfer(p) = fstor2tran * frootn_storage(p)/dt
                     if (woody(ivt(p)) == 1.0_r8) then
                        livestemn_storage_to_xfer(p)  = fstor2tran * livestemn_storage(p)/dt
                        deadstemn_storage_to_xfer(p)  = fstor2tran * deadstemn_storage(p)/dt
                        livecrootn_storage_to_xfer(p) = fstor2tran * livecrootn_storage(p)/dt
                        deadcrootn_storage_to_xfer(p) = fstor2tran * deadcrootn_storage(p)/dt
                     end if
                  end if  ! use_matrixcn
               end if

               ! test for switching from growth period to offset period
            else if (offset_flag(p) == 0.0_r8) then
               if (use_cndv) then
                  ! If days_active > 355, then remove patch in
                  ! CNDVEstablishment at the end of the year.
                  ! days_active > 355 is a symptom of seasonal decid. patches occurring in
                  ! gridcells where dayl never drops below crit_dayl.
                  ! This results in TLAI>1e4 in a few gridcells.
                  days_active(p) = days_active(p) + fracday
                  if (days_active(p) > 355._r8) pftmayexist(p) = .false.
               end if

               crit_daylat = SeasonalCriticalDaylength( g, p )
               
               ! only begin to test for offset daylength once past the summer sol
               if (ws_flag == 0._r8 .and. dayl(g) < crit_daylat) then
                  offset_flag(p) = 1._r8
                  offset_counter(p) = ndays_off * secspday
                  prev_leafc_to_litter(p) = 0._r8
                  prev_frootc_to_litter(p) = 0._r8
               end if
            end if

         end if ! end if seasonal deciduous

      end do ! end of patch loop

    end associate
 
  end subroutine CNSeasonDecidPhenology

  !-----------------------------------------------------------------------
  function SeasonalCriticalDaylength( g, p ) result( crit_daylat )
    !
    ! !DESCRIPTION:
    ! Function to determine the critical day length needed for seasonal
    ! decidious leaf offset. When depends on latitude it's higher for
    ! high latitudes and lower for temperate regions.
    !
    ! !ARGUMENTS:
    integer, intent(IN) :: g ! Gridcell index
    integer, intent(IN) :: p ! Patch index
    real(r8) :: crit_daylat  ! Return value
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

    select case( critical_daylight_method )
    ! Critical day length depends on both vegetation type and latitude
    case(critical_daylight_depends_on_latnveg)
        ! Critical day length for offset is fixed for temperate type vegetation
        if ( pftcon%season_decid_temperate(patch%itype(p)) == 1 )then
           crit_daylat = crit_dayl
        ! For Arctic vegetation -- critical daylength is longer at high latitudes and shorter
        ! at midlatitudes
        else
           crit_daylat=params_inst%crit_dayl_at_high_lat-params_inst%crit_dayl_lat_slope* &
                       (critical_offset_high_lat-abs(grc%latdeg(g)))
           if (crit_daylat < crit_dayl) then
              crit_daylat = crit_dayl !maintain previous offset from White 2001 as minimum
           end if
        end if
    ! Critical day length depends just on vegetation type
    case(critical_daylight_depends_on_veg)
        if ( pftcon%season_decid_temperate(patch%itype(p)) == 1 )then
           crit_daylat = crit_dayl
        else
           crit_daylat=params_inst%crit_dayl_at_high_lat
        end if
    ! Critical day length depends on latitude
    case(critical_daylight_depends_on_lat)
        ! Critical daylength is higher at high latitudes and shorter
        ! for temperatre regions
        crit_daylat=params_inst%crit_dayl_at_high_lat-params_inst%crit_dayl_lat_slope* &
                    (critical_offset_high_lat-abs(grc%latdeg(g)))
        if (crit_daylat < crit_dayl) then
           crit_daylat = crit_dayl !maintain previous offset from White 2001 as minimum
        end if
    ! Critical day length is constant
    case(critical_daylight_constant)
        crit_daylat = crit_dayl
    case default
        call endrun(msg="ERROR SeasonalCriticalDaylength critical_daylight_method not implemented "//errmsg(sourcefile, __LINE__))
    end select

  end function SeasonalCriticalDaylength

  !-----------------------------------------------------------------------
  function SeasonalDecidOnset( onset_gdd, onset_gddflag, soilt, soila10, t_a5min, dayl, &
                               snow_5day, ws_flag, crit_onset_gdd, season_decid_temperate ) &
                       result( do_onset )

    !
    ! !DESCRIPTION:
    ! Function to determine if seasonal deciduous leaf onset should happen.
    !
    ! !USES:
    use shr_const_mod   , only: SHR_CONST_TKFRZ
    ! !ARGUMENTS:
    real(r8), intent(INOUT) :: onset_gdd      ! onset growing degree days 
    real(r8), intent(INOUT) :: onset_gddflag  ! Onset freeze flag
    real(r8), intent(IN)    :: soilt          ! Soil temperature at specific level for this evaluation
    real(r8), intent(IN)    :: soila10        ! 10-day running mean of the 12cm soil layer temperature (K)
    real(r8), intent(IN)    :: t_a5min        ! 5-day running mean of min 2-m temperature
    real(r8), intent(IN)    :: dayl           ! Day length
    real(r8), intent(IN)    :: snow_5day      ! 5-day average of snow
    real(r8), intent(IN)    :: ws_flag        ! winter-summer solstice flag (0 or 1)
    real(r8), intent(IN)    :: crit_onset_gdd ! critical onset growing degree-day sum
    real(r8), intent(IN)    :: season_decid_temperate  ! If this is a temperate seasonal decidious type 
    logical :: do_onset                       ! Flag if onset should happen (return value)
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: min_critical_daylength_onset = 39300._r8/2._r8 ! Minimum daylength for onset to happen
                                                                          ! NOTE above: The 39300/2(19650) value is what we've
                                                                          ! tested with, we are concerned that changing 
                                                                          ! it might change answers. This trigger was just
                                                                          ! added to make sure that onset doesn't happen in Jan/Feb.
                                                                          ! See more notes on this parameter below.
    !-----------------------------------------------------------------------

    do_onset = .false.
    ! Test to turn on growing degree-day sum, if off.
    ! switch on the growing degree day sum on the winter solstice

    if (onset_gddflag == 0._r8 .and. ws_flag == 1._r8) then
        onset_gddflag = 1._r8
        onset_gdd = 0._r8
    end if

    ! Test to turn off growing degree-day sum, if on.
    ! This test resets the growing degree day sum if it gets past
    ! the summer solstice without reaching the threshold value.
    ! In that case, it will take until the next winter solstice
    ! before the growing degree-day summation starts again.

    if (onset_gddflag == 1._r8 .and. ws_flag == 0._r8) then
        onset_gddflag = 0._r8
        onset_gdd = 0._r8
    end if

    ! if the gdd flag is set, and if the soil is above freezing
    ! then accumulate growing degree days for onset trigger

    if (onset_gddflag == 1.0_r8 .and. soilt > SHR_CONST_TKFRZ) then
        onset_gdd = onset_gdd + (soilt-SHR_CONST_TKFRZ)*fracday
    end if
    if ( onset_thresh_depends_on_veg) then
       ! separate into non-arctic seasonally deciduous pfts
       ! (temperate broadleaf deciduous
       ! tree) and arctic/boreal seasonally deciduous pfts (boreal
       ! needleleaf deciduous tree,
       ! boreal broadleaf deciduous tree, boreal broadleaf deciduous
       ! shrub, C3 arctic grass)
       if (onset_gdd > crit_onset_gdd .and.  season_decid_temperate == 1) then
           do_onset = .true.
        ! Note: The check "dayl>min_critical_daylength_onset" in the if
        ! statement was added because for some coastal
        ! points the other triggers could allow onset in January/February
        ! which isn't sustainable and is a degenerate case. To prevent this
        ! condition was added, but now the other conditions aren't triggered
        ! until much later so it's value just needs to be high enough to prevent
        ! the degenerate case of happening too early, and low enough that it
        ! doesn't restrict onset. As such the value of this parameter shouldn't
        ! matter for reasonable values between the two degenerate cases.
        else if (season_decid_temperate == 0 .and.  onset_gddflag == 1.0_r8 .and. &
                soila10 > SHR_CONST_TKFRZ .and. &
                t_a5min > SHR_CONST_TKFRZ .and. ws_flag==1.0_r8 .and. &
                dayl>min_critical_daylength_onset .and. &
                snow_5day<params_inst%snow5d_thresh_for_onset) then
           do_onset = .true.
        end if
    else
       ! set do_onset if critical growing degree-day sum is exceeded
       if (onset_gdd > crit_onset_gdd) do_onset = .true.
    end if

  end function SeasonalDecidOnset

  !-----------------------------------------------------------------------
  subroutine CNStressDecidPhenology (num_soilp, filter_soilp , &
       soilstate_inst, temperature_inst, atm2lnd_inst, wateratm2lndbulk_inst, cnveg_state_inst, &
       cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, &
       cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! This routine handles phenology for vegetation types, such as grasses and
    ! tropical drought deciduous trees, that respond to cold and drought stress
    ! signals and that can have multiple growing seasons in a given year.
    ! This routine allows for the possibility that leaves might persist year-round
    ! in the absence of a suitable stress trigger, by switching to an essentially
    ! evergreen habit, but maintaining a deciduous leaf longevity, while waiting
    ! for the next stress trigger.  This is in contrast to the seasonal deciduous
    ! algorithm (for temperate deciduous trees) that forces a single growing season
    ! per year.
    !
    ! !USES:
    use clm_time_manager , only : get_average_days_per_year
    use CNSharedParamsMod, only : use_fun
    use clm_varcon       , only : secspday
    use shr_const_mod    , only : SHR_CONST_TKFRZ, SHR_CONST_PI
    use CNSharedParamsMod, only : CNParamsShareInst
    !
    ! !ARGUMENTS:
    integer                        , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                        , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(soilstate_type)           , intent(in)    :: soilstate_inst
    type(temperature_type)         , intent(in)    :: temperature_inst
    type(atm2lnd_type)             , intent(in)    :: atm2lnd_inst
    type(wateratm2lndbulk_type)             , intent(in)    :: wateratm2lndbulk_inst
    type(cnveg_state_type)         , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type) , intent(inout) :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    real(r8),parameter :: secspqtrday = secspday / 4  ! seconds per quarter day
    integer :: g,c,p           ! indices
    integer :: fp              ! filter patch index
    real(r8):: avg_dayspyr     ! average days per year
    real(r8):: crit_onset_gdd  ! degree days for onset trigger
    real(r8):: soilt           ! temperature of top soil layer
    real(r8):: psi             ! water stress of top soil layer
    real(r8):: rain_threshold  ! rain threshold for leaf on [mm]
    logical :: additional_onset_condition ! additional condition for leaf onset
    !-----------------------------------------------------------------------

    associate(                                                                                                   & 
         ivt                                 =>    patch%itype                                                   , & ! Input:  [integer   (:)   ]  patch vegetation type                                
         dayl                                =>    grc%dayl                                                    , & ! Input:  [real(r8)  (:)   ]  daylength (s)
         
         prec10                              => wateratm2lndbulk_inst%prec10_patch                                     , & ! Input:  [real(r8) (:)     ]  10-day running mean of tot. precipitation
         leaf_long                           =>    pftcon%leaf_long                                            , & ! Input:  leaf longevity (yrs)                              
         woody                               =>    pftcon%woody                                                , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)
         stress_decid                        =>    pftcon%stress_decid                                         , & ! Input:  binary flag for stress-deciduous leaf habit (0 or 1)
         leafcn                              =>    pftcon%leafcn                                               , & ! Input:  leaf C:N (gC/gN)
         frootcn                             =>    pftcon%frootcn                                              , & ! Input:  fine root C:N (gC/gN) 

         crit_onset_gdd_sf                   =>    pftcon%crit_onset_gdd_sf                                    , & ! Input:  scale factor for crit_onset_gdd (unitless)
         ndays_on                            =>    pftcon%ndays_on                                             , & ! Input:  number of days to complete leaf onset (days)

         soilpsi                             =>    soilstate_inst%soilpsi_col                                  , & ! Input:  [real(r8)  (:,:) ]  soil water potential in each soil layer (MPa)   
         
         t_soisno                            =>    temperature_inst%t_soisno_col                               , & ! Input:  [real(r8)  (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
         dormant_flag                        =>    cnveg_state_inst%dormant_flag_patch                         , & ! Output:  [real(r8) (:)   ]  dormancy flag                                     
         days_active                         =>    cnveg_state_inst%days_active_patch                          , & ! Output:  [real(r8) (:)   ]  number of days since last dormancy                
         onset_flag                          =>    cnveg_state_inst%onset_flag_patch                           , & ! Output:  [real(r8) (:)   ]  onset flag                                        
         onset_counter                       =>    cnveg_state_inst%onset_counter_patch                        , & ! Output:  [real(r8) (:)   ]  onset counter (seconds)                           
         onset_gddflag                       =>    cnveg_state_inst%onset_gddflag_patch                        , & ! Output:  [real(r8) (:)   ]  onset freeze flag                                 
         onset_fdd                           =>    cnveg_state_inst%onset_fdd_patch                            , & ! Output:  [real(r8) (:)   ]  onset freezing degree days counter                
         onset_gdd                           =>    cnveg_state_inst%onset_gdd_patch                            , & ! Output:  [real(r8) (:)   ]  onset growing degree days                         
         onset_swi                           =>    cnveg_state_inst%onset_swi_patch                            , & ! Output:  [real(r8) (:)   ]  onset soil water index                            
         offset_flag                         =>    cnveg_state_inst%offset_flag_patch                          , & ! Output:  [real(r8) (:)   ]  offset flag                                       
         offset_counter                      =>    cnveg_state_inst%offset_counter_patch                       , & ! Output:  [real(r8) (:)   ]  offset counter (seconds)                          
         offset_fdd                          =>    cnveg_state_inst%offset_fdd_patch                           , & ! Output:  [real(r8) (:)   ]  offset freezing degree days counter               
         offset_swi                          =>    cnveg_state_inst%offset_swi_patch                           , & ! Output:  [real(r8) (:)   ]  offset soil water index                           
         lgsf                                =>    cnveg_state_inst%lgsf_patch                                 , & ! Output:  [real(r8) (:)   ]  long growing season factor [0-1]                  
         bglfr                               =>    cnveg_state_inst%bglfr_patch                                , & ! Output:  [real(r8) (:)   ]  background litterfall rate (1/s)                  
         bgtr                                =>    cnveg_state_inst%bgtr_patch                                 , & ! Output:  [real(r8) (:)   ]  background transfer growth rate (1/s)             
         annavg_t2m                          =>    cnveg_state_inst%annavg_t2m_patch                           , & ! Output:  [real(r8) (:)   ]  annual average 2m air temperature (K)             
         leafc                               =>    cnveg_carbonstate_inst%leafc_patch                          , & ! Input:  [real(r8)  (:)   ]  (gC/m2) leaf C 
     
         frootc                              =>    cnveg_carbonstate_inst%frootc_patch                         , & ! Input:  [real(r8) (:)    ]  (gC/m2) fine root C
         leafc_storage                       =>    cnveg_carbonstate_inst%leafc_storage_patch                  , & ! Input:  [real(r8)  (:)   ]  (gC/m2) leaf C storage                            
         frootc_storage                      =>    cnveg_carbonstate_inst%frootc_storage_patch                 , & ! Input:  [real(r8)  (:)   ]  (gC/m2) fine root C storage                       
         livestemc_storage                   =>    cnveg_carbonstate_inst%livestemc_storage_patch              , & ! Input:  [real(r8)  (:)   ]  (gC/m2) live stem C storage                       
         deadstemc_storage                   =>    cnveg_carbonstate_inst%deadstemc_storage_patch              , & ! Input:  [real(r8)  (:)   ]  (gC/m2) dead stem C storage                       
         livecrootc_storage                  =>    cnveg_carbonstate_inst%livecrootc_storage_patch             , & ! Input:  [real(r8)  (:)   ]  (gC/m2) live coarse root C storage                
         deadcrootc_storage                  =>    cnveg_carbonstate_inst%deadcrootc_storage_patch             , & ! Input:  [real(r8)  (:)   ]  (gC/m2) dead coarse root C storage                
         gresp_storage                       =>    cnveg_carbonstate_inst%gresp_storage_patch                  , & ! Input:  [real(r8)  (:)   ]  (gC/m2) growth respiration storage                
         leafc_xfer                          =>    cnveg_carbonstate_inst%leafc_xfer_patch                     , & ! Output:  [real(r8) (:)   ]  (gC/m2) leaf C transfer
         frootc_xfer                         =>    cnveg_carbonstate_inst%frootc_xfer_patch                    , & ! Output:  [real(r8) (:)   ]  (gC/m2) fine root C transfer                      
         livestemc_xfer                      =>    cnveg_carbonstate_inst%livestemc_xfer_patch                 , & ! Output:  [real(r8) (:)   ]  (gC/m2) live stem C transfer                      
         deadstemc_xfer                      =>    cnveg_carbonstate_inst%deadstemc_xfer_patch                 , & ! Output:  [real(r8) (:)   ]  (gC/m2) dead stem C transfer                      
         livecrootc_xfer                     =>    cnveg_carbonstate_inst%livecrootc_xfer_patch                , & ! Output:  [real(r8) (:)   ]  (gC/m2) live coarse root C transfer               
         deadcrootc_xfer                     =>    cnveg_carbonstate_inst%deadcrootc_xfer_patch                , & ! Output:  [real(r8) (:)   ]  (gC/m2) dead coarse root C transfer               
         leafn_storage                       =>    cnveg_nitrogenstate_inst%leafn_storage_patch                , & ! Input:  [real(r8)  (:)   ]  (gN/m2) leaf N storage                            
         frootn_storage                      =>    cnveg_nitrogenstate_inst%frootn_storage_patch               , & ! Input:  [real(r8)  (:)   ]  (gN/m2) fine root N storage                       
         livestemn_storage                   =>    cnveg_nitrogenstate_inst%livestemn_storage_patch            , & ! Input:  [real(r8)  (:)   ]  (gN/m2) live stem N storage                       
         deadstemn_storage                   =>    cnveg_nitrogenstate_inst%deadstemn_storage_patch            , & ! Input:  [real(r8)  (:)   ]  (gN/m2) dead stem N storage                       
         livecrootn_storage                  =>    cnveg_nitrogenstate_inst%livecrootn_storage_patch           , & ! Input:  [real(r8)  (:)   ]  (gN/m2) live coarse root N storage                
         deadcrootn_storage                  =>    cnveg_nitrogenstate_inst%deadcrootn_storage_patch           , & ! Input:  [real(r8)  (:)   ]  (gN/m2) dead coarse root N storage                
         leafn_xfer                          =>    cnveg_nitrogenstate_inst%leafn_xfer_patch                   , & ! Output:  [real(r8) (:)   ]  (gN/m2) leaf N transfer                           
         frootn_xfer                         =>    cnveg_nitrogenstate_inst%frootn_xfer_patch                  , & ! Output:  [real(r8) (:)   ]  (gN/m2) fine root N transfer                      
         livestemn_xfer                      =>    cnveg_nitrogenstate_inst%livestemn_xfer_patch               , & ! Output:  [real(r8) (:)   ]  (gN/m2) live stem N transfer                      
         deadstemn_xfer                      =>    cnveg_nitrogenstate_inst%deadstemn_xfer_patch               , & ! Output:  [real(r8) (:)   ]  (gN/m2) dead stem N transfer                      
         livecrootn_xfer                     =>    cnveg_nitrogenstate_inst%livecrootn_xfer_patch              , & ! Output:  [real(r8) (:)   ]  (gN/m2) live coarse root N transfer               
         deadcrootn_xfer                     =>    cnveg_nitrogenstate_inst%deadcrootn_xfer_patch              , & ! Output:  [real(r8) (:)   ]  (gN/m2) dead coarse root N transfer               
         
         prev_leafc_to_litter                =>    cnveg_carbonflux_inst%prev_leafc_to_litter_patch            , & ! Output:  [real(r8) (:)   ]  previous timestep leaf C litterfall flux (gC/m2/s)
         prev_frootc_to_litter               =>    cnveg_carbonflux_inst%prev_frootc_to_litter_patch           , & ! Output:  [real(r8) (:)   ]  previous timestep froot C litterfall flux (gC/m2/s)
         leafc_xfer_to_leafc                 =>    cnveg_carbonflux_inst%leafc_xfer_to_leafc_patch             , & ! Output:  [real(r8) (:)   ]                                                    
         frootc_xfer_to_frootc               =>    cnveg_carbonflux_inst%frootc_xfer_to_frootc_patch           , & ! Output:  [real(r8) (:)   ]                                                    
         livestemc_xfer_to_livestemc         =>    cnveg_carbonflux_inst%livestemc_xfer_to_livestemc_patch     , & ! Output:  [real(r8) (:)   ]                                                    
         deadstemc_xfer_to_deadstemc         =>    cnveg_carbonflux_inst%deadstemc_xfer_to_deadstemc_patch     , & ! Output:  [real(r8) (:)   ]                                                    
         livecrootc_xfer_to_livecrootc       =>    cnveg_carbonflux_inst%livecrootc_xfer_to_livecrootc_patch   , & ! Output:  [real(r8) (:)   ]                                                    
         deadcrootc_xfer_to_deadcrootc       =>    cnveg_carbonflux_inst%deadcrootc_xfer_to_deadcrootc_patch   , & ! Output:  [real(r8) (:)   ]                                                    
         leafc_storage_to_xfer               =>    cnveg_carbonflux_inst%leafc_storage_to_xfer_patch           , & ! Output:  [real(r8) (:)   ]                                                    
         frootc_storage_to_xfer              =>    cnveg_carbonflux_inst%frootc_storage_to_xfer_patch          , & ! Output:  [real(r8) (:)   ]                                                    
         livestemc_storage_to_xfer           =>    cnveg_carbonflux_inst%livestemc_storage_to_xfer_patch       , & ! Output:  [real(r8) (:)   ]                                                    
         deadstemc_storage_to_xfer           =>    cnveg_carbonflux_inst%deadstemc_storage_to_xfer_patch       , & ! Output:  [real(r8) (:)   ]                                                    
         livecrootc_storage_to_xfer          =>    cnveg_carbonflux_inst%livecrootc_storage_to_xfer_patch      , & ! Output:  [real(r8) (:)   ]                                                    
         deadcrootc_storage_to_xfer          =>    cnveg_carbonflux_inst%deadcrootc_storage_to_xfer_patch      , & ! Output:  [real(r8) (:)   ]                                                    
         gresp_storage_to_xfer               =>    cnveg_carbonflux_inst%gresp_storage_to_xfer_patch           , & ! Output:  [real(r8) (:)   ]                                                    
         
         leafn_xfer_to_leafn                 =>    cnveg_nitrogenflux_inst%leafn_xfer_to_leafn_patch           , & ! Output:  [real(r8) (:)   ]                                                    
         frootn_xfer_to_frootn               =>    cnveg_nitrogenflux_inst%frootn_xfer_to_frootn_patch         , & ! Output:  [real(r8) (:)   ]                                                    
         livestemn_xfer_to_livestemn         =>    cnveg_nitrogenflux_inst%livestemn_xfer_to_livestemn_patch   , & ! Output:  [real(r8) (:)   ]                                                    
         deadstemn_xfer_to_deadstemn         =>    cnveg_nitrogenflux_inst%deadstemn_xfer_to_deadstemn_patch   , & ! Output:  [real(r8) (:)   ]                                                    
         livecrootn_xfer_to_livecrootn       =>    cnveg_nitrogenflux_inst%livecrootn_xfer_to_livecrootn_patch , & ! Output:  [real(r8) (:)   ]                                                    
         deadcrootn_xfer_to_deadcrootn       =>    cnveg_nitrogenflux_inst%deadcrootn_xfer_to_deadcrootn_patch , & ! Output:  [real(r8) (:)   ]                                                    
         leafn_storage_to_xfer               =>    cnveg_nitrogenflux_inst%leafn_storage_to_xfer_patch         , & ! Output:  [real(r8) (:)   ]                                                    
         frootn_storage_to_xfer              =>    cnveg_nitrogenflux_inst%frootn_storage_to_xfer_patch        , & ! Output:  [real(r8) (:)   ]                                                    
         livestemn_storage_to_xfer           =>    cnveg_nitrogenflux_inst%livestemn_storage_to_xfer_patch     , & ! Output:  [real(r8) (:)   ]                                                    
         deadstemn_storage_to_xfer           =>    cnveg_nitrogenflux_inst%deadstemn_storage_to_xfer_patch     , & ! Output:  [real(r8) (:)   ]                                                    
         livecrootn_storage_to_xfer          =>    cnveg_nitrogenflux_inst%livecrootn_storage_to_xfer_patch    , & ! Output:  [real(r8) (:)   ]                                                    
         deadcrootn_storage_to_xfer          =>    cnveg_nitrogenflux_inst%deadcrootn_storage_to_xfer_patch    , & ! Output:  [real(r8) (:)		 ] 
         ileafst_to_ileafxf_phc              =>    cnveg_carbonflux_inst%ileafst_to_ileafxf_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from leaf storage pool to leaf transfer pool
         ileafxf_to_ileaf_phc                =>    cnveg_carbonflux_inst%ileafxf_to_ileaf_ph                   , & ! Input: [integer (:)] Index of phenology related C transfer from leaf transfer pool to leaf pool
         ifrootst_to_ifrootxf_phc            =>    cnveg_carbonflux_inst%ifrootst_to_ifrootxf_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from fine root storage pool to fine root transfer pool
         ifrootxf_to_ifroot_phc              =>    cnveg_carbonflux_inst%ifrootxf_to_ifroot_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from fine root transfer pool to fine root pool
         ilivestemst_to_ilivestemxf_phc      =>    cnveg_carbonflux_inst%ilivestemst_to_ilivestemxf_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live stem storage pool to live stem transfer pool
         ilivestemxf_to_ilivestem_phc        =>    cnveg_carbonflux_inst%ilivestemxf_to_ilivestem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem transfer pool to live stem pool
         ideadstemst_to_ideadstemxf_phc      =>    cnveg_carbonflux_inst%ideadstemst_to_ideadstemxf_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem storage pool to dead stem transfer pool
         ideadstemxf_to_ideadstem_phc        =>    cnveg_carbonflux_inst%ideadstemxf_to_ideadstem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem transfer pool to dead stem pool
         ilivecrootst_to_ilivecrootxf_phc    =>    cnveg_carbonflux_inst%ilivecrootst_to_ilivecrootxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root storage pool to live coarse root transfer pool
         ilivecrootxf_to_ilivecroot_phc      =>    cnveg_carbonflux_inst%ilivecrootxf_to_ilivecroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root transfer pool to live coarse root pool
         ideadcrootst_to_ideadcrootxf_phc    =>    cnveg_carbonflux_inst%ideadcrootst_to_ideadcrootxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root storage pool to dead coarse root transfer pool
         ideadcrootxf_to_ideadcroot_phc      =>    cnveg_carbonflux_inst%ideadcrootxf_to_ideadcroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root transfer pool to dead coarse root pool
         ilivestem_to_ideadstem_phc          =>    cnveg_carbonflux_inst%ilivestem_to_ideadstem_ph             , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to dead stem pool
         ilivecroot_to_ideadcroot_phc        =>    cnveg_carbonflux_inst%ilivecroot_to_ideadcroot_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root to dead coarse root pool
         ileaf_to_iout_phc                   =>    cnveg_carbonflux_inst%ileaf_to_iout_ph                      , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to outside of vegetation pools
         ifroot_to_iout_phc                  =>    cnveg_carbonflux_inst%ifroot_to_iout_ph                     , & ! Input: [integer (:)] Index of phenology related C transfer from fine root pool to outside of vegetation pools
         ilivestem_to_iout_phc               =>    cnveg_carbonflux_inst%ilivestem_to_iout_ph                  , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to outside of vegetation pools
         igrain_to_iout_phc                  =>    cnveg_carbonflux_inst%igrain_to_iout_ph                     , & ! Input: [integer (:)] Index of phenology related C transfer from grain pool to outside of vegetation pools
         ileafst_to_ileafxf_phn              =>    cnveg_nitrogenflux_inst%ileafst_to_ileafxf_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from leaf storage pool to leaf transfer pool
         ileafxf_to_ileaf_phn                =>    cnveg_nitrogenflux_inst%ileafxf_to_ileaf_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from leaf transfer pool to leaf pool
         ifrootst_to_ifrootxf_phn            =>    cnveg_nitrogenflux_inst%ifrootst_to_ifrootxf_ph             , & ! Input: [integer (:)] Index of phenology related C transfer from fine root storage pool to fine root transfer pool
         ifrootxf_to_ifroot_phn              =>    cnveg_nitrogenflux_inst%ifrootxf_to_ifroot_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from fine root transfer pool to fine root pool
         ilivestemst_to_ilivestemxf_phn      =>    cnveg_nitrogenflux_inst%ilivestemst_to_ilivestemxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live stem storage pool to live stem transfer pool
         ilivestemxf_to_ilivestem_phn        =>    cnveg_nitrogenflux_inst%ilivestemxf_to_ilivestem_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live stem transfer pool to live stem pool
         ideadstemst_to_ideadstemxf_phn      =>    cnveg_nitrogenflux_inst%ideadstemst_to_ideadstemxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem storage pool to dead stem transfer pool
         ideadstemxf_to_ideadstem_phn        =>    cnveg_nitrogenflux_inst%ideadstemxf_to_ideadstem_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem transfer pool to dead stem pool
         ilivecrootst_to_ilivecrootxf_phn    =>    cnveg_nitrogenflux_inst%ilivecrootst_to_ilivecrootxf_ph     , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root storage pool to live coarse root transfer pool
         ilivecrootxf_to_ilivecroot_phn      =>    cnveg_nitrogenflux_inst%ilivecrootxf_to_ilivecroot_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root transfer pool to live coarse root pool
         ideadcrootst_to_ideadcrootxf_phn    =>    cnveg_nitrogenflux_inst%ideadcrootst_to_ideadcrootxf_ph     , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root storage pool to dead coarse root transfer pool
         ideadcrootxf_to_ideadcroot_phn      =>    cnveg_nitrogenflux_inst%ideadcrootxf_to_ideadcroot_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root transfer pool to dead coarse root pool
         ilivestem_to_ideadstem_phn          =>    cnveg_nitrogenflux_inst%ilivestem_to_ideadstem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to dead stem pool
         ilivecroot_to_ideadcroot_phn        =>    cnveg_nitrogenflux_inst%ilivecroot_to_ideadcroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root to dead coarse root pool
         ileaf_to_iout_phn                   =>    cnveg_nitrogenflux_inst%ileaf_to_iout_ph                    , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to outside of vegetation pools
         ifroot_to_iout_phn                  =>    cnveg_nitrogenflux_inst%ifroot_to_iout_ph                   , & ! Input: [integer (:)] Index of phenology related C transfer from fine root pool to outside of vegetation pools
         ilivestem_to_iout_phn               =>    cnveg_nitrogenflux_inst%ilivestem_to_iout_ph                , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to outside of vegetation pools
         ileaf_to_iretransn_phn              =>    cnveg_nitrogenflux_inst%ileaf_to_iretransn_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to retranslocation pool
         ilivestem_to_iretransn_phn          =>    cnveg_nitrogenflux_inst%ilivestem_to_iretransn_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to retranslocation pool
         ilivecroot_to_iretransn_phn         =>    cnveg_nitrogenflux_inst%ilivecroot_to_iretransn_ph          , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root pool to retranslocation pool
         igrain_to_iout_phn                  =>    cnveg_nitrogenflux_inst%igrain_to_iout_ph                     & ! Input: [integer (:)] Index of phenology related C transfer from grain pool to outside of vegetation pools
         )

      avg_dayspyr = get_average_days_per_year()

      ! specify rain threshold for leaf onset
      rain_threshold = 20._r8

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = patch%column(p)
         g = patch%gridcell(p)

         if (stress_decid(ivt(p)) == 1._r8) then
            soilt = t_soisno(c, phenology_soil_layer)
            psi = soilpsi(c, phenology_soil_layer)

            ! onset gdd sum from Biome-BGC, v4.1.2
            crit_onset_gdd = crit_onset_gdd_sf(ivt(p)) * exp(4.8_r8 + 0.13_r8*(annavg_t2m(p) &
                             - SHR_CONST_TKFRZ))

            ! update offset_counter and test for the end of the offset period
            if (offset_flag(p) == 1._r8) then
               ! decrement counter for offset period
               offset_counter(p) = offset_counter(p) - dt

               ! if this is the end of the offset_period, reset phenology
               ! flags and indices
               if (offset_counter(p) < dt/2._r8) then
                  ! this code block was originally handled by call cn_offset_cleanup(p)
                  ! inlined during vectorization
                  offset_flag(p) = 0._r8
                  offset_counter(p) = 0._r8
                  dormant_flag(p) = 1._r8
                  days_active(p) = 0._r8

                  ! reset the previous timestep litterfall flux memory
                  prev_leafc_to_litter(p) = 0._r8
                  prev_frootc_to_litter(p) = 0._r8
               end if
            end if

            ! update onset_counter and test for the end of the onset period
            if (onset_flag(p) == 1.0_r8) then
               ! decrement counter for onset period
               onset_counter(p) = onset_counter(p) - dt

               ! if this is the end of the onset period, reset phenology
               ! flags and indices
               if (onset_counter(p) < dt/2._r8) then
                  ! this code block was originally handled by call cn_onset_cleanup(p)
                  ! inlined during vectorization
                  onset_flag(p) = 0._r8
                  onset_counter(p) = 0._r8
                  ! set all transfer growth rates to 0.0
                  leafc_xfer_to_leafc(p)   = 0._r8
                  frootc_xfer_to_frootc(p) = 0._r8
                  leafn_xfer_to_leafn(p)   = 0._r8
                  frootn_xfer_to_frootn(p) = 0._r8
                  if (woody(ivt(p)) == 1.0_r8) then
                     livestemc_xfer_to_livestemc(p)   = 0._r8
                     deadstemc_xfer_to_deadstemc(p)   = 0._r8
                     livecrootc_xfer_to_livecrootc(p) = 0._r8
                     deadcrootc_xfer_to_deadcrootc(p) = 0._r8
                     livestemn_xfer_to_livestemn(p)   = 0._r8
                     deadstemn_xfer_to_deadstemn(p)   = 0._r8
                     livecrootn_xfer_to_livecrootn(p) = 0._r8
                     deadcrootn_xfer_to_deadcrootn(p) = 0._r8
                  end if
                  ! set transfer pools to 0.0
                  leafc_xfer(p) = 0._r8
                  leafn_xfer(p) = 0._r8
                  frootc_xfer(p) = 0._r8
                  frootn_xfer(p) = 0._r8
                  if (woody(ivt(p)) == 1.0_r8) then
                     livestemc_xfer(p) = 0._r8
                     livestemn_xfer(p) = 0._r8
                     deadstemc_xfer(p) = 0._r8
                     deadstemn_xfer(p) = 0._r8
                     livecrootc_xfer(p) = 0._r8
                     livecrootn_xfer(p) = 0._r8
                     deadcrootc_xfer(p) = 0._r8
                     deadcrootn_xfer(p) = 0._r8
                  end if
               end if
            end if

            ! test for switching from dormant period to growth period
            if (dormant_flag(p) == 1._r8) then

               ! keep track of the number of freezing degree days in this
               ! dormancy period (only if the freeze flag has not previously been set
               ! for this dormancy period

               if (onset_gddflag(p) == 0._r8 .and. soilt < SHR_CONST_TKFRZ) onset_fdd(p) = onset_fdd(p) + fracday

               ! if the number of freezing degree days exceeds a critical value,
               ! then onset will require both wet soils and a critical soil
               ! temperature sum.  If this case is triggered, reset any previously
               ! accumulated value in onset_swi, so that onset now depends on
               ! the accumulated soil water index following the freeze trigger

               if (onset_fdd(p) > crit_onset_fdd) then
                  onset_gddflag(p) = 1._r8
                  onset_fdd(p) = 0._r8
                  onset_swi(p) = 0._r8
               end if

               ! if the freeze flag is set, and if the soil is above freezing
               ! then accumulate growing degree days for onset trigger

               if (onset_gddflag(p) == 1._r8 .and. soilt > SHR_CONST_TKFRZ) then
                  onset_gdd(p) = onset_gdd(p) + (soilt-SHR_CONST_TKFRZ)*fracday
               end if

              ! if soils are wet, accumulate soil water index for onset trigger
               additional_onset_condition = .true.
               if(CNParamsShareInst%constrain_stress_deciduous_onset) then
                  ! if additional constraint condition not met,  set to false
                  if ((prec10(p) * (3600.0_r8*10.0_r8*24.0_r8)) < rain_threshold) then
                     additional_onset_condition = .false.
                  endif
               endif

               if (psi >= soilpsi_on) then
                  onset_swi(p) = onset_swi(p) + fracday
               endif

               ! if critical soil water index is exceeded, set onset_flag, and
               ! then test for soil temperature criteria

               ! Adding in Kyla's rainfall trigger when fun on. RF. prec10 (mm/s) needs to be higher than 8mm over 10 days. 

               if (onset_swi(p) > crit_onset_swi.and. additional_onset_condition)  then
                  onset_flag(p) = 1._r8
              
                  ! only check soil temperature criteria if freeze flag set since
                  ! beginning of last dormancy.  If freeze flag set and growing
                  ! degree day sum (since freeze trigger) is lower than critical
                  ! value, then override the onset_flag set from soil water.

                  if (onset_gddflag(p) == 1._r8 .and. onset_gdd(p) < crit_onset_gdd) onset_flag(p) = 0._r8
               end if

               ! only allow onset if dayl > 6hrs
               if (onset_flag(p) == 1._r8 .and. dayl(g) <= secspqtrday) then
                  onset_flag(p) = 0._r8
               end if

               ! if this is the beginning of the onset period
               ! then reset the phenology flags and indices

               if (onset_flag(p) == 1._r8) then
                  dormant_flag(p) = 0._r8
                  days_active(p) = 0._r8
                  onset_gddflag(p) = 0._r8
                  onset_fdd(p) = 0._r8
                  onset_gdd(p) = 0._r8
                  onset_swi(p) = 0._r8
                  onset_counter(p) = ndays_on(ivt(p)) * secspday

                  ! call subroutine to move all the storage pools into transfer pools,
                  ! where they will be transfered to displayed growth over the onset period.
                  ! this code was originally handled with call cn_storage_to_xfer(p)
                  ! inlined during vectorization

                  ! set carbon fluxes for shifting storage pools to transfer pools
                  if (use_matrixcn) then 
                     leafc_storage_to_xfer(p)  = leafc_storage(p)  * matrix_update_phc(p,ileafst_to_ileafxf_phc,fstor2tran/dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                     frootc_storage_to_xfer(p) = frootc_storage(p) * matrix_update_phc(p,ifrootst_to_ifrootxf_phc,fstor2tran/dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                     if (woody(ivt(p)) == 1.0_r8) then
                        livestemc_storage_to_xfer(p)  = livestemc_storage(p)  * matrix_update_phc(p,ilivestemst_to_ilivestemxf_phc,fstor2tran/dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                        deadstemc_storage_to_xfer(p)  = deadstemc_storage(p)  * matrix_update_phc(p,ideadstemst_to_ideadstemxf_phc,fstor2tran/dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                        livecrootc_storage_to_xfer(p) = livecrootc_storage(p) * matrix_update_phc(p,ilivecrootst_to_ilivecrootxf_phc,fstor2tran/dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                        deadcrootc_storage_to_xfer(p) = deadcrootc_storage(p) * matrix_update_phc(p,ideadcrootst_to_ideadcrootxf_phc,fstor2tran/dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                     end if

                     leafn_storage_to_xfer(p)  = leafn_storage(p)  * matrix_update_phn(p,ileafst_to_ileafxf_phn,fstor2tran/dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                     frootn_storage_to_xfer(p) = frootn_storage(p) * matrix_update_phn(p,ifrootst_to_ifrootxf_phn,fstor2tran/dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                     if (woody(ivt(p)) == 1.0_r8) then
                        livestemn_storage_to_xfer(p)  = livestemn_storage(p)  * matrix_update_phn(p,ilivestemst_to_ilivestemxf_phn,fstor2tran/dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                        deadstemn_storage_to_xfer(p)  = deadstemn_storage(p)  * matrix_update_phn(p,ideadstemst_to_ideadstemxf_phn,fstor2tran/dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                        livecrootn_storage_to_xfer(p) = livecrootn_storage(p) * matrix_update_phn(p,ilivecrootst_to_ilivecrootxf_phn,fstor2tran/dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                        deadcrootn_storage_to_xfer(p) = deadcrootn_storage(p) * matrix_update_phn(p,ideadcrootst_to_ideadcrootxf_phn,fstor2tran/dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                     end if
                  else
                     ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
                     !                                        and CNNStateUpdate1::NStateUpdate1
                     leafc_storage_to_xfer(p)  = fstor2tran * leafc_storage(p)/dt
                     frootc_storage_to_xfer(p) = fstor2tran * frootc_storage(p)/dt
                     if (woody(ivt(p)) == 1.0_r8) then
                        livestemc_storage_to_xfer(p)  = fstor2tran * livestemc_storage(p)/dt
                        deadstemc_storage_to_xfer(p)  = fstor2tran * deadstemc_storage(p)/dt
                        livecrootc_storage_to_xfer(p) = fstor2tran * livecrootc_storage(p)/dt
                        deadcrootc_storage_to_xfer(p) = fstor2tran * deadcrootc_storage(p)/dt
                        gresp_storage_to_xfer(p)      = fstor2tran * gresp_storage(p)/dt
                     end if

                     ! set nitrogen fluxes for shifting storage pools to transfer pools
                     leafn_storage_to_xfer(p)  = fstor2tran * leafn_storage(p)/dt
                     frootn_storage_to_xfer(p) = fstor2tran * frootn_storage(p)/dt
                     if (woody(ivt(p)) == 1.0_r8) then
                        livestemn_storage_to_xfer(p)  = fstor2tran * livestemn_storage(p)/dt
                        deadstemn_storage_to_xfer(p)  = fstor2tran * deadstemn_storage(p)/dt
                        livecrootn_storage_to_xfer(p) = fstor2tran * livecrootn_storage(p)/dt
                        deadcrootn_storage_to_xfer(p) = fstor2tran * deadcrootn_storage(p)/dt
                     end if
                  end if
               end if

               ! test for switching from growth period to offset period
            else if (offset_flag(p) == 0._r8) then

               ! if soil water potential lower than critical value, accumulate
               ! as stress in offset soil water index

               if (psi <= soilpsi_off) then
                  offset_swi(p) = offset_swi(p) + fracday

                  ! if the offset soil water index exceeds critical value, and
                  ! if this is not the middle of a previously initiated onset period,
                  ! then set flag to start the offset period and reset index variables

                  if (offset_swi(p) >= crit_offset_swi .and. onset_flag(p) == 0._r8) offset_flag(p) = 1._r8

                  ! if soil water potential higher than critical value, reduce the
                  ! offset water stress index.  By this mechanism, there must be a
                  ! sustained period of water stress to initiate offset.

               else if (psi >= soilpsi_on) then
                  offset_swi(p) = offset_swi(p) - fracday
                  offset_swi(p) = max(offset_swi(p),0._r8)
               end if

               ! decrease freezing day accumulator for warm soil
               if (offset_fdd(p) > 0._r8 .and. soilt > SHR_CONST_TKFRZ) then
                  offset_fdd(p) = offset_fdd(p) - fracday
                  offset_fdd(p) = max(0._r8, offset_fdd(p))
               end if

               ! increase freezing day accumulator for cold soil
               if (soilt <= SHR_CONST_TKFRZ) then
                  offset_fdd(p) = offset_fdd(p) + fracday

                  ! if freezing degree day sum is greater than critical value, initiate offset
                  if (offset_fdd(p) > crit_offset_fdd .and. onset_flag(p) == 0._r8) offset_flag(p) = 1._r8
               end if

               ! force offset if daylength is < 6 hrs
               if (dayl(g) <= secspqtrday) then
                  offset_flag(p) = 1._r8
               end if

               ! if this is the beginning of the offset period
               ! then reset flags and indices
               if (offset_flag(p) == 1._r8) then
                  offset_fdd(p) = 0._r8
                  offset_swi(p) = 0._r8
                  offset_counter(p) = ndays_off * secspday
                  prev_leafc_to_litter(p) = 0._r8
                  prev_frootc_to_litter(p) = 0._r8
               end if
            end if

            ! keep track of number of days since last dormancy for control on
            ! fraction of new growth to send to storage for next growing season

            if (dormant_flag(p) == 0.0_r8) then
               days_active(p) = days_active(p) + fracday
            end if

            ! calculate long growing season factor (lgsf)
            ! only begin to calculate a lgsf greater than 0.0 once the number
            ! of days active exceeds days/year.
            lgsf(p) = max(min(3.0_r8*(days_active(p)-leaf_long(ivt(p))*avg_dayspyr )/avg_dayspyr, 1._r8),0._r8)
            ! RosieF. 5 Nov 2015.  Changed this such that the increase in leaf turnover is faster after
            ! trees enter the 'fake evergreen' state. Otherwise, they have a whole year of 
            ! cheating, with less litterfall than they should have, resulting in very high LAI. 
            ! Further, the 'fake evergreen' state (where lgsf>0) is entered at the end of a single leaf lifespan
            ! and not a whole year. The '3' is arbitrary, given that this entire system is quite abstract. 


            ! set background litterfall rate, when not in the phenological offset period
            if (offset_flag(p) == 1._r8) then
               bglfr(p) = 0._r8
            else
               ! calculate the background litterfall rate (bglfr)
               ! in units 1/s, based on leaf longevity (yrs) and correction for long growing season

               bglfr(p) = (1._r8/(leaf_long(ivt(p))*avg_dayspyr*secspday))*lgsf(p)
            end if

            ! set background transfer rate when active but not in the phenological onset period
            if (onset_flag(p) == 1._r8) then
               bgtr(p) = 0._r8
            else
               ! the background transfer rate is calculated as the rate that would result
               ! in complete turnover of the storage pools in one year at steady state,
               ! once lgsf has reached 1.0 (after 730 days active).

               bgtr(p) = (1._r8/(avg_dayspyr*secspday))*lgsf(p)

               ! set carbon fluxes for shifting storage pools to transfer pools

               ! reduced the amount of stored carbon flowing to display pool by only counting the delta
               ! between leafc and leafc_store in the flux. RosieF, Nov5 2015. 
               leafc_storage_to_xfer(p)  = max(0.0_r8,(leafc_storage(p)-leafc(p))) * bgtr(p)
               frootc_storage_to_xfer(p) = max(0.0_r8,(frootc_storage(p)-frootc(p))) * bgtr(p)
               if (use_matrixcn) then
                  if(leafc_storage(p) .gt. 0)then
                     leafc_storage_to_xfer(p) = leafc_storage(p) * matrix_update_phc(p,ileafst_to_ileafxf_phc,&
                                                leafc_storage_to_xfer(p) / leafc_storage(p), dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                  else
                     leafc_storage_to_xfer(p) = 0
                  end if
                  if(frootc_storage(p) .gt. 0)then
                     frootc_storage_to_xfer(p) = frootc_storage(p) * matrix_update_phc(p,ifrootst_to_ifrootxf_phc,&
                                                 frootc_storage_to_xfer(p) / frootc_storage(p), dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                  else
                     frootc_storage_to_xfer(p) = 0   
                  end if
                if (woody(ivt(p)) == 1.0_r8) then
                     livestemc_storage_to_xfer(p)  = livestemc_storage(p)  * matrix_update_phc(p,ilivestemst_to_ilivestemxf_phc ,bgtr(p),dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                     deadstemc_storage_to_xfer(p)  = deadstemc_storage(p)  * matrix_update_phc(p,ideadstemst_to_ideadstemxf_phc ,bgtr(p),dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                     livecrootc_storage_to_xfer(p) = livecrootc_storage(p) * matrix_update_phc(p,ilivecrootst_to_ilivecrootxf_phc,bgtr(p),dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                     deadcrootc_storage_to_xfer(p) = deadcrootc_storage(p) * matrix_update_phc(p,ideadcrootst_to_ideadcrootxf_phc,bgtr(p),dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                 end if
              else
                 ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
                 !                                        and CNNStateUpdate1::NStateUpdate1
                 if (woody(ivt(p)) == 1.0_r8) then
                    livestemc_storage_to_xfer(p)  = livestemc_storage(p) * bgtr(p)
                    deadstemc_storage_to_xfer(p)  = deadstemc_storage(p) * bgtr(p)
                    livecrootc_storage_to_xfer(p) = livecrootc_storage(p) * bgtr(p)
                    deadcrootc_storage_to_xfer(p) = deadcrootc_storage(p) * bgtr(p)
                    gresp_storage_to_xfer(p)      = gresp_storage(p) * bgtr(p)
                 end if
              end if !use_matrixcn

               ! set nitrogen fluxes for shifting storage pools to transfer pools
               if (use_matrixcn) then 
                  leafn_storage_to_xfer(p)  = leafn_storage(p)  * matrix_update_phn(p,ileafst_to_ileafxf_phn,bgtr(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                  frootn_storage_to_xfer(p) = frootn_storage(p) * matrix_update_phn(p,ifrootst_to_ifrootxf_phn,bgtr(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                  if (woody(ivt(p)) == 1.0_r8) then
                     livestemn_storage_to_xfer(p)  = livestemn_storage(p)  * matrix_update_phn(p,ilivestemst_to_ilivestemxf_phn,bgtr(p) ,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                     deadstemn_storage_to_xfer(p)  = deadstemn_storage(p)  * matrix_update_phn(p,ideadstemst_to_ideadstemxf_phn,bgtr(p) ,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                     livecrootn_storage_to_xfer(p) = livecrootn_storage(p) * matrix_update_phn(p,ilivecrootst_to_ilivecrootxf_phn,bgtr(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                     deadcrootn_storage_to_xfer(p) = deadcrootn_storage(p) * matrix_update_phn(p,ideadcrootst_to_ideadcrootxf_phn,bgtr(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                  end if
               else
                  ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
                  !                                        and CNNStateUpdate1::NStateUpdate1
                  leafn_storage_to_xfer(p)  = leafn_storage(p) * bgtr(p)
                  frootn_storage_to_xfer(p) = frootn_storage(p) * bgtr(p)
                  if (woody(ivt(p)) == 1.0_r8) then
                     livestemn_storage_to_xfer(p)  = livestemn_storage(p) * bgtr(p)
                     deadstemn_storage_to_xfer(p)  = deadstemn_storage(p) * bgtr(p)
                     livecrootn_storage_to_xfer(p) = livecrootn_storage(p) * bgtr(p)
                     deadcrootn_storage_to_xfer(p) = deadcrootn_storage(p) * bgtr(p)
                  end if
               end if !use_matrixcn
            end if

         end if ! end if stress deciduous

      end do ! end of patch loop

    end associate

  end subroutine CNStressDecidPhenology


  !-----------------------------------------------------------------------
  subroutine get_swindow(jday, rx_starts, rx_ends, param_start, param_end, w, start_w, end_w)
    ! !DESCRIPTION:
    ! Determine when the "next" sowing window is. This is either the sowing window we are
    ! currently in or, if not in a sowing window, the next one that will occur.

    ! !USES:
    use clm_time_manager , only : get_curr_days_per_year, is_doy_in_interval, get_doy_tomorrow
    ! !ARGUMENTS:
    integer,                          intent(in)    :: jday ! Day of year
    integer, dimension(:), intent(in)               :: rx_starts, rx_ends ! All prescribed sowing window start and end dates for this patch
    integer,                          intent(in)    :: param_start, param_end ! Sowing window start and end dates from parameter file
    integer,                          intent(out)   :: w ! Index of "next" sowing window
    integer,                          intent(out)   :: start_w, end_w ! Start and end dates of "next" sowing window
    !
    ! !LOCAL VARIABLES
    integer :: jday_tomorrow
    integer :: mxsowings_in ! Due to unit testing, we can't assume the length of the rx sowing window arrays is mxsowings as set in clm_varpar

    ! Initialize
    w = -1
    start_w = -1
    end_w   = -1

    ! Get info
    jday_tomorrow = get_doy_tomorrow(jday)
    mxsowings_in = size(rx_starts)

    ! If no sowing windows are prescribed, use the values from the parameter file.
    if (maxval(rx_starts) < 1) then
        w = 1
        start_w = param_start
        end_w   = param_end
        return

    ! Otherwise, if today is after the latest sowing window end date, use the first sowing window. This works only if sowing windows that span the new year are located at index w = 1.
    else if (jday > maxval(rx_ends)) then
        w = 1
        start_w = rx_starts(w)
        end_w   = rx_ends(w)
        return
    end if

    ! Otherwise, use the first prescribed sowing window we find whose end is >= today. This works only if sowing windows that span the new year are located at index w = 1.
    do w = 1, mxsowings_in
      ! If nothing prescribed at this w, stop looking and exit loop. Will trigger "No sowing window found" error, which we do not move here because it's possible that no start or end date is < 1.
        if (min(rx_starts(w), rx_ends(w)) < 1) then
            exit
        end if

        if (jday <= rx_ends(w)) then
            start_w = rx_starts(w)
            end_w   = rx_ends(w)
            exit
        end if
    end do

    ! Ensure that a window was found.
    ! SSR 2023-10-17: This shouldn't currently be reachable, but its being here casts the widest possible net in case code changes in future.
    if (start_w < 1 .or. end_w < 1) then
        call endrun(msg="get_swindow(): No sowing window found")
    end if

  end subroutine get_swindow


  !-----------------------------------------------------------------------
  function was_sown_in_this_window(sowing_window_startdate, sowing_window_enddate, jday, idop, sown_in_this_window)
    ! !DESCRIPTION:
    ! Determine whether the crop was sown in the current sowing window. Although sown_in_this_window is set to false in last timestep of sowing window at the end of CropPhenology(), these extra checks may be necessary if sowing windows change.
    !
    ! !USES:
    use clm_time_manager , only : is_doy_in_interval
    ! !ARGUMENTS:
    integer, intent(in)    :: sowing_window_startdate, sowing_window_enddate, jday, idop
    logical, intent(in)    :: sown_in_this_window
    ! !LOCAL VARIABLES
    logical :: is_in_sowing_window, idop_in_sowing_window
    ! !RESULT
    logical :: was_sown_in_this_window

    was_sown_in_this_window = sown_in_this_window

    ! If not in a sowing window, sown_in_this_window must be false.
    is_in_sowing_window  = is_doy_in_interval(sowing_window_startdate, sowing_window_enddate, jday)
    if (.not. is_in_sowing_window) then
        was_sown_in_this_window = .false.
        return
    end if

    ! If we're in a sowing window but the day of planting isn't in the active sowing window, we must be in a different sowing window.
    idop_in_sowing_window  = is_doy_in_interval(sowing_window_startdate, sowing_window_enddate, idop)
    if (is_in_sowing_window .and. .not. idop_in_sowing_window) then
        was_sown_in_this_window = .false.
        return
    end if

    ! Sometimes we're in an active sowing window, and the patch was sown between the start and end dates of the window, but not *the currently active* window. Note that windows with start==end are not checked here; we always trust the input value of sown_in_this_window in such cases.
    if (sowing_window_startdate < sowing_window_enddate .and. idop > jday) then
        was_sown_in_this_window = .false.
    else if (sowing_window_startdate > sowing_window_enddate) then
        if (jday <= sowing_window_enddate .and. idop <= sowing_window_enddate .and. idop > jday) then
            was_sown_in_this_window = .false.
        else if (jday >= sowing_window_startdate .and. (idop > jday .or. idop <= sowing_window_enddate)) then
            was_sown_in_this_window = .false.
        end if
    end if

  end function was_sown_in_this_window


  !-----------------------------------------------------------------------
  subroutine CropPhenology(num_pcropp, filter_pcropp                     , &
       waterdiagnosticbulk_inst, temperature_inst, crop_inst, canopystate_inst, cnveg_state_inst , &
       cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst,&
       c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst)

    ! !DESCRIPTION:
    ! Code from AgroIBIS to determine crop phenology and code from CN to
    ! handle CN fluxes during the phenological onset                       & offset periods.
    
    ! !USES:
    use clm_time_manager , only : get_prev_calday, get_curr_days_per_year, is_beg_curr_year
    use clm_time_manager , only : get_average_days_per_year
    use clm_time_manager , only : get_prev_date
    use clm_time_manager , only : is_doy_in_interval, is_end_curr_day
    use clm_time_manager , only : get_doy_tomorrow
    use pftconMod        , only : ntmp_corn, nswheat, nwwheat, ntmp_soybean
    use pftconMod        , only : nirrig_tmp_corn, nirrig_swheat, nirrig_wwheat, nirrig_tmp_soybean
    use pftconMod        , only : ntrp_corn, nsugarcane, ntrp_soybean, ncotton, nrice
    use pftconMod        , only : nirrig_trp_corn, nirrig_sugarcane, nirrig_trp_soybean
    use pftconMod        , only : nirrig_cotton, nirrig_rice
    use pftconMod        , only : nmiscanthus, nirrig_miscanthus, nswitchgrass, nirrig_switchgrass
    
    use clm_varcon       , only : spval, secspday
    use clm_varctl       , only : use_fertilizer 
    use clm_varctl       , only : use_c13, use_c14
    use clm_varcon       , only : c13ratio, c14ratio
    use clm_varctl       , only : use_cropcal_rx_swindows
    !
    ! !ARGUMENTS:
    integer                        , intent(in)    :: num_pcropp       ! number of prog crop patches in filter
    integer                        , intent(in)    :: filter_pcropp(:) ! filter for prognostic crop patches
    type(waterdiagnosticbulk_type)          , intent(in)    :: waterdiagnosticbulk_inst
    type(temperature_type)         , intent(in)    :: temperature_inst
    type(crop_type)                , intent(inout) :: crop_inst
    type(canopystate_type)         , intent(in)    :: canopystate_inst
    type(cnveg_state_type)         , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type) , intent(inout) :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: c13_cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: c14_cnveg_carbonstate_inst
    !
    ! LOCAL VARAIBLES:
    integer jday      ! julian day of the year
    integer fp,p      ! patch indices
    integer c         ! column indices
    integer g         ! gridcell indices
    integer h         ! hemisphere indices
    integer s         ! growing season indices
    integer k         ! grain pool indices
    integer w         ! sowing window index
    integer idpp      ! number of days past planting
    integer mxmat     ! maximum growing season length
    integer kyr       ! current year
    integer kmo       ! month of year  (1, ..., 12)
    integer kda       ! day of month   (1, ..., 31)
    integer mcsec     ! seconds of day (0, ..., seconds/day)
    integer sowing_window_startdate ! date (day of year) of first day of sowing window
    integer sowing_window_enddate   ! date (day of year) of last  day of sowing window
    real(r8) harvest_reason
    real(r8) dayspyr  ! days per year in this year
    real(r8) avg_dayspyr ! average number of days per year
    real(r8) crmcorn  ! comparitive relative maturity for corn
    real(r8) ndays_on ! number of days to fertilize
    logical has_rx_sowing_date ! does the crop have a single sowing date instead of a window?
    logical is_in_sowing_window ! is the crop in its sowing window?
    logical is_end_sowing_window ! is it the last day of the crop's sowing window?
    logical sowing_gdd_requirement_met ! has the gridcell historically been warm enough to support the crop?
    logical do_plant_normal ! are the normal planting rules defined and satisfied?
    logical do_plant_lastchance ! if not the above, what about relaxed rules for the last day of the planting window?
    logical do_plant_prescribed ! is today the prescribed sowing date?
    logical do_plant_prescribed_tomorrow  ! is tomorrow the prescribed sowing date?
    logical do_plant  ! are we planting in this time step for any reason?
    logical did_plant ! did we plant the crop in this time step?
    logical allow_unprescribed_planting ! should crop be allowed to be planted according to sowing window rules?
    logical do_harvest    ! Are harvest conditions satisfied?
    logical fake_harvest  ! Dealing with incorrect Dec. 31 planting
    logical did_plant_prescribed_today    ! Was the crop sown today?
    logical vernalization_forces_harvest ! Was the crop killed by freezing during vernalization?
    !------------------------------------------------------------------------

    associate(                                                                   & 
         ivt               =>    patch%itype                                     , & ! Input:  [integer  (:) ]  patch vegetation type                                
         
         leaf_long         =>    pftcon%leaf_long                                , & ! Input:  leaf longevity (yrs)                              
         leafcn            =>    pftcon%leafcn                                 , & ! Input:  leaf C:N (gC/gN)                                  
         manunitro         =>    pftcon%manunitro                              , & ! Input:  max manure to be applied in total (kgN/m2)
         minplanttemp      =>    pftcon%minplanttemp                           , & ! Input:  
         planttemp         =>    pftcon%planttemp                              , & ! Input:  
         gddmin            =>    pftcon%gddmin                                 , & ! Input:  
         lfemerg           =>    pftcon%lfemerg                                , & ! Input:  
         grnfill           =>    pftcon%grnfill                               , & ! Input:  

         t_ref2m_min       =>    temperature_inst%t_ref2m_min_patch            , & ! Input:  [real(r8) (:) ]  daily minimum of average 2 m height surface air temperature (K)
         t10               =>    temperature_inst%t_a10_patch                  , & ! Input:  [real(r8) (:) ]  10-day running mean of the 2 m temperature (K)    
         a5tmin            =>    temperature_inst%t_a5min_patch                , & ! Input:  [real(r8) (:) ]  5-day running mean of min 2-m temperature         
         a10tmin           =>    temperature_inst%t_a10min_patch               , & ! Input:  [real(r8) (:) ]  10-day running mean of min 2-m temperature        
         gdd020            =>    temperature_inst%gdd020_patch                 , & ! Input:  [real(r8) (:) ]  20 yr mean of gdd0                                
         gdd820            =>    temperature_inst%gdd820_patch                 , & ! Input:  [real(r8) (:) ]  20 yr mean of gdd8                                

         fertnitro         =>    crop_inst%fertnitro_patch                     , & ! Input:  [real(r8) (:) ]  fertilizer nitrogen
         hui               =>    crop_inst%hui_patch                           , & ! Input:  [real(r8) (:) ]  crop patch heat unit index (growing degree-days); set to 0 at sowing and accumulated until harvest
         leafout           =>    crop_inst%gddtsoi_patch                       , & ! Input:  [real(r8) (:) ]  gdd from top soil layer temperature
         harvdate          =>    crop_inst%harvdate_patch                      , & ! Output: [integer  (:) ]  harvest date                                       
         croplive          =>    crop_inst%croplive_patch                      , & ! Output: [logical  (:) ]  Flag, true if planted, not harvested
         vf                =>    crop_inst%vf_patch                            , & ! Output: [real(r8) (:) ]  vernalization factor                              
         sowing_count      =>    crop_inst%sowing_count                        , & ! Inout:  [integer  (:) ]  number of sowing events this year for this patch
         harvest_count     =>    crop_inst%harvest_count                       , & ! Inout:  [integer  (:) ]  number of harvest events this year for this patch
         peaklai           =>    cnveg_state_inst%peaklai_patch                , & ! Output: [integer  (:) ]  1: max allowed lai; 0: not at max
         tlai              =>    canopystate_inst%tlai_patch                   , & ! Input:  [real(r8) (:) ]  one-sided leaf area index, no burying by snow     
         
         idop              =>    cnveg_state_inst%idop_patch                   , & ! Output: [integer  (:) ]  date of planting (day of year)
         iyop              =>    cnveg_state_inst%iyop_patch                   , & ! Output: [integer  (:) ]  year of planting (day of year)
         gddmaturity       =>    cnveg_state_inst%gddmaturity_patch            , & ! Output: [real(r8) (:) ]  gdd needed to harvest                             
         huileaf           =>    cnveg_state_inst%huileaf_patch                , & ! Output: [real(r8) (:) ]  heat unit index needed from planting to leaf emergence
         huigrain          =>    cnveg_state_inst%huigrain_patch               , & ! Output: [real(r8) (:) ]  same to reach vegetative maturity                 
         cumvd             =>    cnveg_state_inst%cumvd_patch                  , & ! Output: [real(r8) (:) ]  cumulative vernalization d?ependence?             
         hdidx             =>    cnveg_state_inst%hdidx_patch                  , & ! Output: [real(r8) (:) ]  cold hardening index?                             
         bglfr             =>    cnveg_state_inst%bglfr_patch                  , & ! Output: [real(r8) (:) ]  background litterfall rate (1/s)                  
         bgtr              =>    cnveg_state_inst%bgtr_patch                   , & ! Output: [real(r8) (:) ]  background transfer growth rate (1/s)             
         lgsf              =>    cnveg_state_inst%lgsf_patch                   , & ! Output: [real(r8) (:) ]  long growing season factor [0-1]                  
         onset_flag        =>    cnveg_state_inst%onset_flag_patch             , & ! Output: [real(r8) (:) ]  onset flag                                        
         offset_flag       =>    cnveg_state_inst%offset_flag_patch            , & ! Output: [real(r8) (:) ]  offset flag                                       
         onset_counter     =>    cnveg_state_inst%onset_counter_patch          , & ! Output: [real(r8) (:) ]  onset counter                                     
         offset_counter    =>    cnveg_state_inst%offset_counter_patch         , & ! Output: [real(r8) (:) ]  offset counter                                    
         
         leafc_xfer        =>    cnveg_carbonstate_inst%leafc_xfer_patch       , & ! Output: [real(r8) (:) ]  (gC/m2)   leaf C transfer                           

         crop_seedc_to_leaf =>   cnveg_carbonflux_inst%crop_seedc_to_leaf_patch, & ! Output: [real(r8) (:) ]  (gC/m2/s) seed source to leaf

         fert_counter      =>    cnveg_nitrogenflux_inst%fert_counter_patch    , & ! Output: [real(r8) (:) ]  >0 fertilize; <=0 not (seconds)                   
         leafn_xfer        =>    cnveg_nitrogenstate_inst%leafn_xfer_patch     , & ! Output: [real(r8) (:) ]  (gN/m2)   leaf N transfer                           
         crop_seedn_to_leaf =>   cnveg_nitrogenflux_inst%crop_seedn_to_leaf_patch, & ! Output: [real(r8) (:) ]  (gN/m2/s) seed source to leaf
         cphase            =>    crop_inst%cphase_patch                        , & ! Output: [real(r8) (:)]   phenology phase
         fert              =>    cnveg_nitrogenflux_inst%fert_patch              & ! Output: [real(r8) (:) ]  (gN/m2/s) fertilizer applied each timestep 
         )

      ! get time info
      dayspyr = get_curr_days_per_year()
      avg_dayspyr = get_average_days_per_year()
      jday    = get_prev_calday()
      call get_prev_date(kyr, kmo, kda, mcsec)

      if (use_fertilizer) then
       ndays_on = 20._r8 ! number of days to fertilize
      else
       ndays_on = 0._r8 ! number of days to fertilize
      end if

      do fp = 1, num_pcropp
         p = filter_pcropp(fp)
         c = patch%column(p)
         g = patch%gridcell(p)
         h = inhemi(p)

         ! background litterfall and transfer rates; long growing season factor

         bglfr(p) = 0._r8 ! this value changes later in a crop's life cycle
         bgtr(p)  = 0._r8
         lgsf(p)  = 0._r8

         ! Should never be saved as zero, but including this so it's initialized just in case
         harvest_reason = 0._r8

         ! ---------------------------------
         ! from AgroIBIS subroutine planting
         ! ---------------------------------

         ! initialize other variables that are calculated for crops
         ! on an annual basis in cropresidue subroutine

         ! Second condition ensures everything is correctly set when resuming from a run with old code
         ! OR starting a run mid-year without any restart file OR handling a new crop column that just
         ! came into existence (and not at the year boundary for some reason).
         if ( is_beg_curr_year() .or. crop_inst%sdates_thisyr_patch(p,1) == spval ) then
            sowing_count(p) = 0
            harvest_count(p) = 0
            do s = 1, mxsowings
               crop_inst%sdates_thisyr_patch(p,s) = -1._r8
               crop_inst%swindow_starts_thisyr_patch(p,s) = -1._r8
               crop_inst%swindow_ends_thisyr_patch  (p,s) = -1._r8
               crop_inst%sowing_reason_thisyr_patch(p,s) = -1._r8
            end do
            do s = 1, mxharvests
               crop_inst%sdates_perharv_patch(p,s) = -1._r8
               crop_inst%syears_perharv_patch(p,s) = -1._r8
               crop_inst%hdates_thisyr_patch(p,s) = -1._r8
               cnveg_state_inst%gddmaturity_thisyr(p,s) = -1._r8
               crop_inst%gddaccum_thisyr_patch(p,s) = -1._r8
               crop_inst%hui_thisyr_patch(p,s) = -1._r8
               crop_inst%sowing_reason_perharv_patch(p,s) = -1._r8
               crop_inst%harvest_reason_thisyr_patch(p,s) = -1._r8
               do k = repr_grain_min, repr_grain_max
                  cnveg_carbonflux_inst%repr_grainc_to_food_perharv_patch(p,s,k) = 0._r8
                  cnveg_carbonflux_inst%repr_grainc_to_seed_perharv_patch(p,s,k) = 0._r8
                  cnveg_nitrogenflux_inst%repr_grainn_to_food_perharv_patch(p,s,k) = 0._r8
                  cnveg_nitrogenflux_inst%repr_grainn_to_seed_perharv_patch(p,s,k) = 0._r8
               end do
            end do
            do k = repr_grain_min, repr_grain_max
               cnveg_carbonflux_inst%repr_grainc_to_food_thisyr_patch(p,k) = 0._r8
               cnveg_carbonflux_inst%repr_grainc_to_seed_thisyr_patch(p,k) = 0._r8
               cnveg_nitrogenflux_inst%repr_grainn_to_food_thisyr_patch(p,k) = 0._r8
               cnveg_nitrogenflux_inst%repr_grainn_to_seed_thisyr_patch(p,k) = 0._r8
            end do
         end if

         ! Get dates of current or next sowing window.
         call get_swindow(jday, crop_inst%rx_swindow_starts_thisyr_patch(p,:), crop_inst%rx_swindow_ends_thisyr_patch(p,:), minplantjday(ivt(p),h), maxplantjday(ivt(p),h), w, sowing_window_startdate, sowing_window_enddate)

         ! Are we currently in a sowing window?
         ! This is outside the croplive check so that the "harvest if planting conditions were met today" conditional works.
         is_in_sowing_window  = is_doy_in_interval(sowing_window_startdate, sowing_window_enddate, jday)
         crop_inst%sown_in_this_window(p) = was_sown_in_this_window(sowing_window_startdate, sowing_window_enddate, jday, idop(p), crop_inst%sown_in_this_window(p))
         is_end_sowing_window = jday == sowing_window_enddate

         ! We only want to plant on a specific day if the prescribed sowing window starts AND ends on the same day. Also make sure we haven't planted yet today.
         has_rx_sowing_date = sowing_window_startdate == sowing_window_enddate
         do_plant_prescribed = has_rx_sowing_date .and. &
                               sowing_window_startdate == jday .and. &
                               .not. crop_inst%sown_in_this_window(p)
         do_plant_prescribed_tomorrow = &
             has_rx_sowing_date .and. &
             sowing_window_startdate == get_doy_tomorrow(jday)

         ! BACKWARDS_COMPATIBILITY(wjs/ssr, 2022-02-18)
         ! When resuming from a run with old code, may need to manually set these.
         ! Will be needed until we can rely on all restart files have been generated
         ! with CropPhenology() getting the day of the year from the START of the timestep
         ! (i.e., jday = get_prev_calday()) instead of the END of the timestep (i.e.,
         ! jday = get_calday()). See CTSM issue #1623.
         ! Once removed, can also remove the "Instead, always harvest the day before idop" bit.
         if (croplive(p) .and. idop(p) <= jday .and. sowing_count(p) == 0) then
             sowing_count(p) = 1
             crop_inst%sdates_thisyr_patch(p,1) = real(idop(p), r8)
         end if

         ! Save these diagnostic variables only on the last day of the window to ensure that windows spanning the new year aren't double-counted. Doing this on the last day ensures that outputs are ordered as inputs should be.
         if (jday == sowing_window_enddate) then
             crop_inst%swindow_starts_thisyr_patch(p,w) = sowing_window_startdate
             crop_inst%swindow_ends_thisyr_patch  (p,w) = sowing_window_enddate
         end if
         !
         ! Only allow sowing according to normal "window" rules if not using prescribed
         ! sowing dates.
         allow_unprescribed_planting = .not. has_rx_sowing_date
         if (sowing_count(p) == mxsowings) then
            do_plant_normal = .false.
            do_plant_lastchance = .false.
         else if (ivt(p) == nwwheat .or. ivt(p) == nirrig_wwheat) then
            ! winter temperate cereal : use gdd0 as a limit to plant winter cereal
            sowing_gdd_requirement_met = gdd020(p) /= spval .and. gdd020(p) >= gddmin(ivt(p))
            ! Are all the normal requirements for planting met?
            do_plant_normal = allow_unprescribed_planting           .and. &
                              a5tmin(p)   /= spval                  .and. &
                              a5tmin(p)   <= minplanttemp(ivt(p))   .and. &
                              is_in_sowing_window                   .and. &
                              sowing_gdd_requirement_met
            ! If not, but it's the last day of the planting window, what about relaxed rules?
            do_plant_lastchance = allow_unprescribed_planting           .and. &
                                  (.not. do_plant_normal)               .and. &
                                  is_end_sowing_window                  .and. &
                                  sowing_gdd_requirement_met
         else ! not winter cereal... slevis: added distinction between NH and SH
            ! slevis: The idea is that jday will equal idop sooner or later in the year
            !         while the gdd part is either true or false for the year.
            ! Are all the normal requirements for planting met?
            do_plant_normal = allow_unprescribed_planting               .and. &
                              t10(p) /= spval .and. a10tmin(p) /= spval .and. &
                              t10(p)     > planttemp(ivt(p))            .and. &
                              a10tmin(p) > minplanttemp(ivt(p))         .and. &
                              is_in_sowing_window                       .and. &
                              gdd820(p)  /= spval                       .and. &
                              gdd820(p)  >= gddmin(ivt(p))
            ! If not, but it's the last day of the planting window, what about relaxed rules?
            do_plant_lastchance = allow_unprescribed_planting    .and. &
                                  (.not. do_plant_normal)        .and. &
                                  is_end_sowing_window           .and. &
                                  gdd820(p) > 0._r8 .and. &
                                  gdd820(p) /= spval
         end if
         do_plant = do_plant_prescribed .or. do_plant_normal .or. do_plant_lastchance
         do_plant = do_plant .and. .not. crop_inst%sown_in_this_window(p)
         did_plant = .false.

         ! Once outputs can handle >1 planting per year, remove 2nd condition.
         if ( (.not. croplive(p)) .and. sowing_count(p) == 0 ) then

            ! gdd needed for * chosen crop and a likely hybrid (for that region) *
            ! to reach full physiological maturity

            ! based on accumulated seasonal average growing degree days from
            ! April 1 - Sept 30 (inclusive)
            ! for corn and soybeans in the United States -
            ! decided upon by what the typical average growing season length is
            ! and the gdd needed to reach maturity in those regions

            ! first choice is used for spring temperate cereal and/or soybeans and maize

            ! slevis: ibis reads xinpdate in io.f from control.crops.nc variable name 'plantdate'
            !         According to Chris Kucharik, the dataset of
            !         xinpdate was generated from a previous model run at 0.5 deg resolution

            if (do_plant) then

               if (ivt(p) == nwwheat .or. ivt(p) == nirrig_wwheat) then
                  cumvd(p)       = 0._r8
                  hdidx(p)       = 0._r8
                  vf(p)          = 0._r8
               end if

               call PlantCrop(p, leafcn(ivt(p)), jday, kyr, do_plant_normal, &
                              do_plant_lastchance, do_plant_prescribed,  &
                              temperature_inst, crop_inst, cnveg_state_inst, &
                              cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, &
                              cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
                              c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst)
               did_plant = .true.

            else
               gddmaturity(p) = 0._r8
            end if

            ! crop phenology (gdd thresholds) controlled by gdd needed for
            ! maturity (physiological) which is based on the average gdd
            ! accumulation and hybrids in United States from April 1 - Sept 30,
            ! unless using cultivar GDD target inputs

            ! calculate threshold from phase 1 to phase 2:
            ! threshold for attaining leaf emergence (based on fraction of
            ! gdd(i) -- climatological average)
            ! Hayhoe and Dwyer, 1990, Can. J. Soil Sci 70:493-497
            ! Carlson and Gage, 1989, Agric. For. Met., 45: 313-324
            ! J.T. Ritchie, 1991: Modeling Plant and Soil systems

            huileaf(p) = lfemerg(ivt(p)) * gddmaturity(p) ! 3-7% in cereal

            ! calculate threshhold from phase 2 to phase 3:
            ! from leaf emergence to beginning of grain-fill period
            ! this hypothetically occurs at the end of tassling, not the beginning
            ! tassel initiation typically begins at 0.5-0.55 * gddmaturity

            ! calculate linear relationship between huigrain fraction and relative
            ! maturity rating for maize

            if (ivt(p) == ntmp_corn .or. ivt(p) == nirrig_tmp_corn .or. &
                ivt(p) == ntrp_corn .or. ivt(p) == nirrig_trp_corn .or. &
                ivt(p) == nsugarcane .or. ivt(p) == nirrig_sugarcane .or. &
                ivt(p) == nmiscanthus .or. ivt(p) == nirrig_miscanthus .or. &
                ivt(p) == nswitchgrass .or. ivt(p) == nirrig_switchgrass) then
               ! the following estimation of crmcorn from gddmaturity is based on a linear
               ! regression using data from Pioneer-brand corn hybrids (Kucharik, 2003,
               ! Earth Interactions 7:1-33: fig. 2)
               crmcorn = max(73._r8, min(135._r8, (gddmaturity(p)+ 53.683_r8)/13.882_r8))

               ! the following adjustment of grnfill based on crmcorn is based on a tuning
               ! of Agro-IBIS to give reasonable results for max LAI and the seasonal
               ! progression of LAI growth (pers. comm. C. Kucharik June 10, 2010)
               huigrain(p) = -0.002_r8  * (crmcorn - 73._r8) + grnfill(ivt(p))

               huigrain(p) = min(max(huigrain(p), grnfill(ivt(p))-0.1_r8), grnfill(ivt(p)))
               huigrain(p) = huigrain(p) * gddmaturity(p)     ! Cabelguenne et
            else
               huigrain(p) = grnfill(ivt(p)) * gddmaturity(p) ! al. 1999
            end if

         end if ! crop not live nor planted

         ! ----------------------------------
         ! from AgroIBIS subroutine phenocrop
         ! ----------------------------------

         ! all of the phenology changes are based on the total number of gdd needed
         ! to change to the next phase - based on fractions of the total gdd typical
         ! for  that region based on the April 1 - Sept 30 window of development

         ! crop phenology (gdd thresholds) controlled by gdd needed for
         ! maturity (physiological) which is based on the average gdd
         ! accumulation and hybrids in United States from April 1 - Sept 30

         ! Phase 1: Planting to leaf emergence (now in CNAllocation)
         ! Phase 2: Leaf emergence to beginning of grain fill (general LAI accumulation)
         ! Phase 3: Grain fill to physiological maturity and harvest (LAI decline)
         ! Harvest: if gdd past grain fill initiation exceeds limit
         ! or number of days past planting reaches a maximum, the crop has
         ! reached physiological maturity and plant is harvested;
         ! crop could be live or dead at this stage - these limits
         ! could lead to reaching physiological maturity or determining
         ! a harvest date for a crop killed by an early frost (see next comments)
         ! --- --- ---
         ! keeping comments without the code (slevis):
         ! if minimum temperature, t_ref2m_min <= freeze kill threshold, tkill
         ! for 3 consecutive days and lai is above a minimum,
         ! plant will be damaged/killed. This function is more for spring freeze events
         ! or for early fall freeze events

         ! spring temperate cereal is affected by this, winter cereal kill function
         ! is determined in crops.f - is a more elaborate function of
         ! cold hardening of the plant

         ! currently simulates too many grid cells killed by freezing temperatures

         ! removed on March 12 2002 - C. Kucharik
         ! until it can be a bit more refined, or used at a smaller scale.
         ! we really have no way of validating this routine
         ! too difficult to implement on 0.5 degree scale grid cells
         ! --- --- ---

         onset_flag(p)  = 0._r8 ! CN terminology to trigger certain
         offset_flag(p) = 0._r8 ! carbon and nitrogen transfers

         if (croplive(p)) then
            cphase(p) = cphase_planted

            ! call vernalization if winter temperate cereal planted, living, and the
            ! vernalization factor is not 1;
            ! vf affects the calculation of gddtsoi & hui

            vernalization_forces_harvest = .false.
            if (t_ref2m_min(p) < 1.e30_r8 .and. vf(p) /= 1._r8 .and. &
               (ivt(p) == nwwheat .or. ivt(p) == nirrig_wwheat)) then
               call vernalization(p, &
                    canopystate_inst, temperature_inst, waterdiagnosticbulk_inst, cnveg_state_inst, &
                    crop_inst, vernalization_forces_harvest)
            end if

            ! days past planting may determine harvest
            idpp = DaysPastPlanting(idop(p), jday)

            ! onset_counter initialized to zero when .not. croplive
            ! offset_counter relevant only at time step of harvest

            onset_counter(p) = onset_counter(p) - dt

            ! enter phase 2 onset for one time step:
            ! transfer seed carbon to leaf emergence

            if (peaklai(p) >= 1) then
               hui(p) = max(hui(p),huigrain(p))
            endif

            do_harvest = .false.
            fake_harvest = .false.
            did_plant_prescribed_today = .false.
            if (use_cropcal_rx_swindows .and. sowing_count(p) > 0) then
                did_plant_prescribed_today = crop_inst%sdates_thisyr_patch(p,sowing_count(p)) == real(jday, r8)
            end if

            ! Optionally ignore maximum growing season length
            mxmat = pftcon%mxmat(ivt(p))
            if (.not. use_mxmat) then
                mxmat = 999
            end if

            if (jday == 1 .and. croplive(p) .and. idop(p) == 1 .and. sowing_count(p) == 0) then
                ! BACKWARDS_COMPATIBILITY(ssr, 2022-02-03): To get rid of crops incorrectly planted in last time step of Dec. 31. That was fixed in commit dadbc62 ("Call CropPhenology regardless of doalb"), but this handles restart files with the old behavior. fake_harvest ensures that outputs aren't polluted.
                do_harvest = .true.
                fake_harvest = .true.
                harvest_reason = HARVEST_REASON_SOWNBADDEC31
            else if (use_cropcal_rx_swindows .and. do_plant .and. .not. did_plant) then
                ! Today was supposed to be the planting day, but the previous crop still hasn't been harvested.
                do_harvest = .true.
                harvest_reason = HARVEST_REASON_SOWTODAY

            ! If generate_crop_gdds and this patch has prescribed sowing inputs
            else if (generate_crop_gdds .and. crop_inst%rx_swindow_starts_thisyr_patch(p,1) .gt. 0) then
               ! Harvest the day before the next prescribed sowing.
               do_harvest = do_plant_prescribed_tomorrow

               ! ... unless that will lead to growing season length 365 (or 366,
               ! if last year was a leap year). This would result in idop==jday,
               ! which would invoke the "manually setting sowing_count and
               ! sdates_thisyr" code. This would lead to crops never getting
               ! harvested. Instead, always harvest the day before idop.
               if ((.not. do_harvest) .and. &
                   (idop(p) > 1 .and. jday == idop(p) - 1) .or. &
                   (idop(p) == 1 .and. jday == dayspyr)) then
                   do_harvest = .true.
                   harvest_reason = HARVEST_REASON_IDOPTOMORROW
               else if (do_harvest) then
                   harvest_reason = HARVEST_REASON_SOWTOMORROW
               end if

            else if (did_plant_prescribed_today) then
               ! Do not harvest on the day this growing season began;
               ! would create challenges for postprocessing.
               do_harvest = .false.
            else if (vernalization_forces_harvest) then
               do_harvest = .true.
               harvest_reason = HARVEST_REASON_VERNFREEZEKILL
            else
               ! Original harvest rule
               do_harvest = hui(p) >= gddmaturity(p) .or. idpp >= mxmat

               ! Always harvest the day before the next prescribed sowing date, if still alive.
               ! WARNING: This implementation assumes that prescribed sowing dates don't change over time!
               ! In order to avoid this, you'd have to read this year's AND next year's prescribed
               ! sowing dates.
               do_harvest = do_harvest .or. do_plant_prescribed_tomorrow

               if (hui(p) >= gddmaturity(p)) then
                   harvest_reason = HARVEST_REASON_MATURE
               else if (idpp >= mxmat) then
                   harvest_reason = HARVEST_REASON_MAXSEASLENGTH
               else if (do_plant_prescribed_tomorrow) then
                   harvest_reason = HARVEST_REASON_SOWTOMORROW
               end if
            endif

            ! The following conditionals are similar to those in CropPhase. However, they
            ! differ slightly because here we are potentially setting a new crop phase,
            ! whereas CropPhase is just designed to get the current, already-determined
            ! phase. However, despite these differences: if you make changes to the
            ! following conditionals, you should also check to see if you should make
            ! similar changes in CropPhase.
            if ((.not. do_harvest) .and. leafout(p) >= huileaf(p) .and. hui(p) < huigrain(p) .and. idpp < mxmat) then
               cphase(p) = cphase_leafemerge
               if (abs(onset_counter(p)) > 1.e-6_r8) then
                  onset_flag(p)    = 1._r8
                  onset_counter(p) = dt
                    fert_counter(p)  = ndays_on * secspday
                    if (ndays_on .gt. 0) then
                       fert(p) = (manunitro(ivt(p)) * 1000._r8 + fertnitro(p))/ fert_counter(p)
                    else
                       fert(p) = 0._r8
                    end if
               else
                  ! this ensures no re-entry to onset of phase2
                  ! b/c onset_counter(p) = onset_counter(p) - dt
                  ! at every time step

                  onset_counter(p) = dt
               end if

               ! enter harvest for one time step:
               ! - transfer live biomass to litter and to crop yield
               ! - send xsmrpool to the atmosphere
               ! if onset and harvest needed to last longer than one timestep
               ! the onset_counter would change from dt and you'd need to make
               ! changes to the offset subroutine below

            else if (do_harvest) then
               ! Don't update these if you're just harvesting because of incorrect Dec.
               ! 31 planting
               if (.not. fake_harvest) then
                  if (harvdate(p) >= NOT_Harvested) harvdate(p) = jday
                  harvest_count(p) = harvest_count(p) + 1
                  crop_inst%sdates_perharv_patch(p, harvest_count(p)) = real(idop(p), r8)
                  crop_inst%syears_perharv_patch(p, harvest_count(p)) = real(iyop(p), r8)
                  crop_inst%hdates_thisyr_patch(p, harvest_count(p)) = real(jday, r8)
                  cnveg_state_inst%gddmaturity_thisyr(p,harvest_count(p)) = gddmaturity(p)
                  crop_inst%gddaccum_thisyr_patch(p, harvest_count(p)) = crop_inst%gddaccum_patch(p)
                  crop_inst%hui_thisyr_patch(p, harvest_count(p)) = hui(p)
                  crop_inst%sowing_reason_perharv_patch(p, harvest_count(p)) = real(crop_inst%sowing_reason_patch(p), r8)
                  crop_inst%sowing_reason_patch(p) = -1 ! "Reason for most recent sowing of this patch." So in the line above we save, and here we reset.
                  crop_inst%harvest_reason_thisyr_patch(p, harvest_count(p)) = harvest_reason
               endif

               croplive(p) = .false.     ! no re-entry in greater if-block
               cphase(p) = cphase_harvest
               if (tlai(p) > 0._r8) then ! plant had emerged before harvest
                  offset_flag(p) = 1._r8
                  offset_counter(p) = dt
               else                      ! plant never emerged from the ground
                  ! Revert planting transfers; this will replenish the crop seed deficit.
                  ! We subtract from any existing value in crop_seedc_to_leaf /
                  ! crop_seedn_to_leaf in the unlikely event that we enter this block of
                  ! code in the same time step where the planting transfer originally
                  ! occurred.
                  crop_seedc_to_leaf(p) = crop_seedc_to_leaf(p) - leafc_xfer(p)/dt
                  crop_seedn_to_leaf(p) = crop_seedn_to_leaf(p) - leafn_xfer(p)/dt
                  leafc_xfer(p) = 0._r8

                  leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p))
                  if (use_c13) then
                     c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = 0._r8
                  endif
                  if (use_c14) then
                     c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = 0._r8
                  endif

               end if

               ! enter phase 3 while previous criteria fail and next is true;
               ! in terms of order, phase 3 occurs before harvest, but when
               ! harvest *can* occur, we want it to have first priority.
               ! AgroIBIS uses a complex formula for lai decline.
               ! Use CN's simple formula at least as a place holder (slevis)

            else if (hui(p) >= huigrain(p)) then
               cphase(p) = cphase_grainfill
               bglfr(p) = 1._r8/(leaf_long(ivt(p))*avg_dayspyr*secspday)
            end if

            ! continue fertilizer application while in phase 2;
            ! assumes that onset of phase 2 took one time step only

              if (fert_counter(p) <= 0._r8) then
                 fert(p) = 0._r8
              else ! continue same fert application every timestep
                 fert_counter(p) = fert_counter(p) - dt
              end if

         else   ! crop not live
            ! next 2 lines conserve mass if leaf*_xfer > 0 due to interpinic.
            ! We subtract from any existing value in crop_seedc_to_leaf /
            ! crop_seedn_to_leaf in the unlikely event that we enter this block of
            ! code in the same time step where the planting transfer originally
            ! occurred.
            crop_seedc_to_leaf(p) = crop_seedc_to_leaf(p) - leafc_xfer(p)/dt
            crop_seedn_to_leaf(p) = crop_seedn_to_leaf(p) - leafn_xfer(p)/dt
            onset_counter(p) = 0._r8
            leafc_xfer(p) = 0._r8
            leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p))
            if (use_c13) then
               c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = 0._r8
            endif
            if (use_c14) then
               c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = 0._r8
            endif
         end if ! croplive

         ! At the end of the sowing window, AFTER we've done everything crop-related, set this to false
         if (is_end_sowing_window .and. is_end_curr_day()) then
            crop_inst%sown_in_this_window(p) = .false.
         end if

      end do ! prognostic crops loop

    end associate

  end subroutine CropPhenology

  !-----------------------------------------------------------------------
  subroutine CropPhase(bounds, num_pcropp, filter_pcropp, &
       crop_inst, cnveg_state_inst, crop_phase)
    !
    ! !DESCRIPTION:
    ! Get the current phase of each crop patch.
    !
    ! The returned values (in crop_phase) are from the set of cphase_* values defined in
    ! CropType. The returned values in crop_phase are only valid for patches where
    ! croplive is true; the values are undefined where croplive is false and should not be
    ! used there!
    !
    ! This has logic similar to that in CropPhenology. If you make changes here, you
    ! should also check if similar changes need to be made in CropPhenology.
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_pcropp       ! number of prog crop patches in filter
    integer                , intent(in)    :: filter_pcropp(:) ! filter for prognostic crop patches
    type(crop_type)        , intent(in)    :: crop_inst
    type(cnveg_state_type) , intent(in)    :: cnveg_state_inst
    real(r8)               , intent(inout) :: crop_phase(bounds%begp:)
    !
    ! !LOCAL VARIABLES:
    integer :: p, fp

    character(len=*), parameter :: subname = 'CropPhase'
    !-----------------------------------------------------------------------
    SHR_ASSERT_ALL_FL((ubound(crop_phase) == [bounds%endp]), sourcefile, __LINE__)

    associate( &
         croplive =>    crop_inst%croplive_patch        , & ! Input: [logical  (:) ]  Flag, true if planted, not harvested
         hui      =>    crop_inst%hui_patch             , & ! Input: [real(r8) (:) ]  gdd since planting (gddplant)
         leafout  =>    crop_inst%gddtsoi_patch         , & ! Input: [real(r8) (:) ]  gdd from top soil layer temperature
         huileaf  =>    cnveg_state_inst%huileaf_patch  , & ! Input: [real(r8) (:) ]  heat unit index needed from planting to leaf emergence
         huigrain =>    cnveg_state_inst%huigrain_patch   & ! Input: [real(r8) (:) ]  same to reach vegetative maturity
         )

    do fp = 1, num_pcropp
       p = filter_pcropp(fp)

       if (croplive(p)) then
          ! Start with cphase_planted, but this might get changed in the later
          ! conditional blocks.
          crop_phase(p) = cphase_planted
          if (leafout(p) >= huileaf(p) .and. hui(p) < huigrain(p)) then
             crop_phase(p) = cphase_leafemerge
          else if (hui(p) >= huigrain(p)) then
             ! Since we know croplive is true, any hui greater than huigrain implies that
             ! we're in the grainfill stage: if we were passt gddmaturity then croplive
             ! would be false.
             crop_phase(p) = cphase_grainfill
          end if
       end if
    end do

    end associate

  end subroutine CropPhase


  !-----------------------------------------------------------------------
  subroutine CropPhenologyInit(bounds)
    !
    ! !DESCRIPTION:
    ! Initialization of CropPhenology. Must be called after time-manager is
    ! initialized, and after pftcon file is read in.
    !
    ! !USES:
    use pftconMod       , only: npcropmin, npcropmax
    use clm_time_manager, only: get_calday
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  
    !
    ! LOCAL VARAIBLES:
    integer           :: p,g,n,i                     ! indices
    !------------------------------------------------------------------------

    allocate( inhemi(bounds%begp:bounds%endp) )

    allocate( minplantjday(0:maxveg,inSH)) ! minimum planting julian day
    allocate( maxplantjday(0:maxveg,inSH)) ! minimum planting julian day

    ! Julian day for the start of the year (mid-winter)
    jdayyrstart(inNH) =   1
    jdayyrstart(inSH) = 182

    ! Convert planting dates into julian day
    minplantjday(:,:) = huge(1)
    maxplantjday(:,:) = huge(1)
    do n = npcropmin, npcropmax
       if (pftcon%is_pft_known_to_model(n)) then
          minplantjday(n, inNH) = int( get_calday( pftcon%mnNHplantdate(n), 0 ) )
          maxplantjday(n, inNH) = int( get_calday( pftcon%mxNHplantdate(n), 0 ) )

          minplantjday(n, inSH) = int( get_calday( pftcon%mnSHplantdate(n), 0 ) )
          maxplantjday(n, inSH) = int( get_calday( pftcon%mxSHplantdate(n), 0 ) )
       end if
    end do

    ! Figure out what hemisphere each PATCH is in
    do p = bounds%begp, bounds%endp
       g = patch%gridcell(p)
       ! Northern hemisphere
       if ( grc%latdeg(g) > 0.0_r8 )then
          inhemi(p) = inNH
       else
          inhemi(p) = inSH
       end if
    end do

    !
    ! Constants for Crop vernalization
    !
    ! photoperiod factor calculation
    ! genetic constant - can be modified

    p1d = 0.004_r8  ! average for genotypes from Ritchey, 1991.
    ! Modeling plant & soil systems: Wheat phasic developmt
    p1v = 0.003_r8  ! average for genotypes from Ritchey, 1991.

    hti   = 1._r8
    tbase = 0._r8

  end subroutine CropPhenologyInit

    !-----------------------------------------------------------------------
  subroutine PlantCrop(p, leafcn_in, jday, kyr, do_plant_normal, &
       do_plant_lastchance, do_plant_prescribed,            &
       temperature_inst, crop_inst, cnveg_state_inst,               &
       cnveg_carbonstate_inst, cnveg_nitrogenstate_inst,            &
       cnveg_carbonflux_inst, cnveg_nitrogenflux_inst,              &
       c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst)
    !
    ! !DESCRIPTION:
    !
    ! subroutine calculates initializes variables to what they should be upon
    ! planting. Includes only operations that apply to all crop types; 
    ! additional operations need to happen in CropPhenology().

    ! !USES:
    use clm_varctl       , only : use_c13, use_c14
    use clm_varctl       , only : use_cropcal_rx_cultivar_gdds, adapt_cropcal_rx_cultivar_gdds
    use clm_varcon       , only : c13ratio, c14ratio
    use clm_varpar       , only : mxsowings
    use pftconMod        , only : ntmp_corn, nswheat, nwwheat, ntmp_soybean
    use pftconMod        , only : nirrig_tmp_corn, nirrig_swheat, nirrig_wwheat, nirrig_tmp_soybean
    use pftconMod        , only : ntrp_corn, nsugarcane, ntrp_soybean, ncotton, nrice
    use pftconMod        , only : nirrig_trp_corn, nirrig_sugarcane, nirrig_trp_soybean
    use pftconMod        , only : nirrig_cotton, nirrig_rice
    use pftconMod        , only : nmiscanthus, nirrig_miscanthus, nswitchgrass, nirrig_switchgrass
    !
    ! !ARGUMENTS:
    integer                , intent(in)    :: p         ! PATCH index running over
    real(r8)               , intent(in)    :: leafcn_in ! leaf C:N (gC/gN) of this patch's vegetation type (pftcon%leafcn(ivt(p)))
    integer                , intent(in)    :: jday      ! julian day of the year
    integer                , intent(in)    :: kyr       ! current year
    logical                , intent(in)    :: do_plant_normal ! Are all the normal requirements for planting met?
    logical                , intent(in)    :: do_plant_lastchance ! Are the last-chance requirements for planting met?
    logical                , intent(in)    :: do_plant_prescribed ! are we planting because it was prescribed?
    type(temperature_type)         , intent(in)    :: temperature_inst
    type(crop_type)                , intent(inout) :: crop_inst
    type(cnveg_state_type)         , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type) , intent(inout) :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: c13_cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)   , intent(inout) :: c14_cnveg_carbonstate_inst
    !
    ! LOCAL VARAIBLES:
    integer s              ! growing season index
    integer k              ! grain pool index
    real(r8) gdd20         ! GDD*20 value for this crop type
    real(r8) gdd_target    ! cultivar GDD target this growing season
    real(r8) this_sowing_reason ! number representing sowing reason(s)
    logical did_rx_gdds    ! did this patch use a prescribed harvest requirement?
    !------------------------------------------------------------------------

    associate(                                                                     & 
         ivt               =>    patch%itype                                     , & ! Input:  [integer  (:) ]  patch vegetation type
         croplive          =>    crop_inst%croplive_patch                        , & ! Output: [logical  (:) ]  Flag, true if planted, not harvested
         harvdate          =>    crop_inst%harvdate_patch                        , & ! Output: [integer  (:) ]  harvest date
         sowing_count      =>    crop_inst%sowing_count                          , & ! Inout:  [integer  (:) ]  number of sowing events this year for this patch
         sowing_reason     =>    crop_inst%sowing_reason_thisyr_patch            , & ! Output:  [real(r8)  (:) ]  reason for each sowing this year for this patch
         gddmaturity       =>    cnveg_state_inst%gddmaturity_patch            , & ! Output: [real(r8) (:) ]  gdd needed to harvest
         idop              =>    cnveg_state_inst%idop_patch                     , & ! Output: [integer  (:) ]  date of planting
         iyop              =>    cnveg_state_inst%iyop_patch                     , & ! Output: [integer  (:) ]  year of planting
         leafc_xfer        =>    cnveg_carbonstate_inst%leafc_xfer_patch         , & ! Output: [real(r8) (:) ]  (gC/m2)   leaf C transfer
         leafn_xfer        =>    cnveg_nitrogenstate_inst%leafn_xfer_patch       , & ! Output: [real(r8) (:) ]  (gN/m2)   leaf N transfer
         crop_seedc_to_leaf =>   cnveg_carbonflux_inst%crop_seedc_to_leaf_patch  , & ! Output: [real(r8) (:) ]  (gC/m2/s) seed source to leaf
         crop_seedn_to_leaf =>   cnveg_nitrogenflux_inst%crop_seedn_to_leaf_patch, & ! Output: [real(r8) (:) ]  (gN/m2/s) seed source to leaf
         hybgdd            =>    pftcon%hybgdd                                 , & ! Input:  [real(r8) (:) ]
         gddmin            =>    pftcon%gddmin                                 , & ! Input:
         gdd020            =>    temperature_inst%gdd020_patch                 , & ! Input:  [real(r8) (:) ]  20 yr mean of gdd0
         gdd820            =>    temperature_inst%gdd820_patch                 , & ! Input:  [real(r8) (:) ]  20 yr mean of gdd8
         gdd1020           =>    temperature_inst%gdd1020_patch                , & ! Input:  [real(r8) (:) ]  20 yr mean of gdd10
         aleafi            =>    cnveg_state_inst%aleafi_patch                 , & ! Output: [real(r8) (:)   ]  saved allocation coefficient from phase 2
         astemi            =>    cnveg_state_inst%astemi_patch                 , & ! Output: [real(r8) (:)   ]  saved allocation coefficient from phase 2
         aleaf             =>    cnveg_state_inst%aleaf_patch                  , & ! Output: [real(r8) (:)   ]  leaf allocation coefficient
         astem             =>    cnveg_state_inst%astem_patch                  , & ! Output: [real(r8) (:)   ]  stem allocation coefficient
         aroot             =>    cnveg_state_inst%aroot_patch                  , & ! Output: [real(r8) (:)   ]  root allocation coefficient
         arepr             =>    cnveg_state_inst%arepr_patch                    & ! Output: [real(r8) (:,:) ]  reproductive allocation coefficient(s)
         )

      ! impose limit on growing season length needed
      ! for crop maturity - for cold weather constraints
      croplive(p)  = .true.
      crop_inst%sown_in_this_window(p) = .true.
      idop(p)      = jday
      iyop(p)      = kyr
      harvdate(p)  = NOT_Harvested
      sowing_count(p) = sowing_count(p) + 1

      crop_inst%sdates_thisyr_patch(p,sowing_count(p)) = real(jday, r8)

      this_sowing_reason = 0._r8
      if (do_plant_prescribed) then
          this_sowing_reason = 10._r8
      end if
      if (do_plant_normal) then
          this_sowing_reason = this_sowing_reason + 1._r8
      else if (do_plant_lastchance) then
          this_sowing_reason = this_sowing_reason + 2._r8
      end if
      sowing_reason(p,sowing_count(p)) = this_sowing_reason
      crop_inst%sowing_reason_patch(p) = this_sowing_reason

      leafc_xfer(p)  = initial_seed_at_planting
      leafn_xfer(p) = leafc_xfer(p) / leafcn_in ! with onset
      crop_seedc_to_leaf(p) = leafc_xfer(p)/dt
      crop_seedn_to_leaf(p) = leafn_xfer(p)/dt

      ! because leafc_xfer is set above rather than incremneted through the normal process, must also set its isotope
      ! pools here.  use totvegc_patch as the closest analogue if nonzero, and use initial value otherwise
      if (use_c13) then
         if ( cnveg_carbonstate_inst%totvegc_patch(p) .gt. 0._r8) then
            c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * &
            c13_cnveg_carbonstate_inst%totvegc_patch(p) / cnveg_carbonstate_inst%totvegc_patch(p)
         else
            c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * c13ratio
         endif
      endif
      if (use_c14) then
         if ( cnveg_carbonstate_inst%totvegc_patch(p) .gt. 0._r8) then
            c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * &
               c14_cnveg_carbonstate_inst%totvegc_patch(p) / cnveg_carbonstate_inst%totvegc_patch(p)
         else
            c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * c14ratio
         endif
      endif

      ! which GDD*20 variable does this crop use?
      if (ivt(p) == ntmp_soybean .or. ivt(p) == nirrig_tmp_soybean .or. &
            ivt(p) == ntrp_soybean .or. ivt(p) == nirrig_trp_soybean) then
         gdd20 = gdd1020(p)
      else if (ivt(p) == ntmp_corn .or. ivt(p) == nirrig_tmp_corn .or. &
            ivt(p) == ntrp_corn .or. ivt(p) == nirrig_trp_corn .or. &
            ivt(p) == nsugarcane .or. ivt(p) == nirrig_sugarcane .or. &
            ivt(p) == nmiscanthus .or. ivt(p) == nirrig_miscanthus .or. &
            ivt(p) == nswitchgrass .or. ivt(p) == nirrig_switchgrass) then
         gdd20 = gdd820(p)
      else if (ivt(p) == nswheat .or. ivt(p) == nirrig_swheat .or. &
            ivt(p) == nwwheat .or. ivt(p) == nirrig_wwheat .or. &
            ivt(p) == ncotton .or. ivt(p) == nirrig_cotton .or. &
            ivt(p) == nrice   .or. ivt(p) == nirrig_rice) then
         gdd20 = gdd020(p)
      else
         write(iulog, *) 'ERROR: PlantCrop(): unrecognized ivt for gdd20: ',ivt(p)
         call endrun(msg="Stopping")
      end if

      ! set GDD target
      did_rx_gdds = .false.
      if (use_cropcal_rx_cultivar_gdds .and. crop_inst%rx_cultivar_gdds_thisyr_patch(p,sowing_count(p)) .ge. 0._r8) then
         gddmaturity(p) = crop_inst%rx_cultivar_gdds_thisyr_patch(p,sowing_count(p))
         did_rx_gdds = .true.
         if (adapt_cropcal_rx_cultivar_gdds .and. crop_inst%gdd20_baseline_patch(p) > min_gdd20_baseline) then
            gddmaturity(p) = gddmaturity(p) * gdd20 / crop_inst%gdd20_baseline_patch(p)
            !TODO Sam Rabin: Set maximum and minimum gddmaturity
         end if
      else if (ivt(p) == nwwheat .or. ivt(p) == nirrig_wwheat) then
         gddmaturity(p) = hybgdd(ivt(p))
      else
         if (ivt(p) == ntmp_soybean .or. ivt(p) == nirrig_tmp_soybean .or. &
               ivt(p) == ntrp_soybean .or. ivt(p) == nirrig_trp_soybean .or. &
               ivt(p) == nswheat .or. ivt(p) == nirrig_swheat .or. &
               ivt(p) == ncotton .or. ivt(p) == nirrig_cotton .or. &
               ivt(p) == nrice   .or. ivt(p) == nirrig_rice) then
            gddmaturity(p) = min(gdd20, hybgdd(ivt(p)))
         else if (ivt(p) == ntmp_corn .or. ivt(p) == nirrig_tmp_corn .or. &
               ivt(p) == ntrp_corn .or. ivt(p) == nirrig_trp_corn .or. &
               ivt(p) == nsugarcane .or. ivt(p) == nirrig_sugarcane .or. &
               ivt(p) == nmiscanthus .or. ivt(p) == nirrig_miscanthus .or. &
               ivt(p) == nswitchgrass .or. ivt(p) == nirrig_switchgrass) then
            gddmaturity(p) = max(950._r8, min(gdd20*0.85_r8, hybgdd(ivt(p))))
            if (do_plant_normal) then
               gddmaturity(p) = max(950._r8, min(gddmaturity(p)+150._r8, 1850._r8))
            end if
         else
            write(iulog, *) 'ERROR: PlantCrop(): unrecognized ivt for GDD target: ',ivt(p)
            call endrun(msg="Stopping")
         end if

      endif

      if (gddmaturity(p) < min_gddmaturity) then
         if (use_cropcal_rx_cultivar_gdds .or. generate_crop_gdds) then
             if (did_rx_gdds) then
                 write(iulog,*) 'Some patch with ivt ',ivt(p),' has rx gddmaturity ',gddmaturity(p),'; using min_gddmaturity instead (',min_gddmaturity,')'
             end if
             gddmaturity(p) = min_gddmaturity
         else
            write(iulog, *) 'ERROR: PlantCrop(): gddmaturity < minimum for ivt ',ivt(p)
            call endrun(msg="Stopping")
         end if
      end if

      ! Initialize allocation coefficients.
      ! Because crops have no live carbon pools when planted but not emerged, this shouldn't
      ! matter unless they skip the vegetative phase (which only happens in very weird run
      ! setups).
      aleaf(p) = 1._r8
      aleafi(p) = 1._r8
      astem(p) = 0._r8
      astemi(p) = 0._r8
      aroot(p) = 0._r8
      do k = 1, nrepr
         arepr(p,k) = 0._r8
      end do

    end associate

  end subroutine PlantCrop

  !-----------------------------------------------------------------------
  function DaysPastPlanting(idop, jday_in)
    ! !USES:
    use clm_time_manager, only : get_prev_calday, get_curr_days_per_year
    !
    ! !ARGUMENTS:
    integer,           intent(in) :: idop ! patch day of planting
    integer, optional, intent(in) :: jday_in ! julian day of the year
    !
    ! !LOCAL VARIABLES
    integer :: DaysPastPlanting
    integer :: jday

    ! Must use separate jday_in and jday because we can't redefine an intent(in)
    ! variable, even if it wasn't provided in the function call.
    if (present(jday_in)) then
       jday = jday_in
    else
       ! Use prev instead of curr to avoid jday=1 in last timestep of year
       jday = get_prev_calday()
    end if

    if (jday >= idop) then
       DaysPastPlanting = jday - idop
    else
      ! As long as crops have at most a 365-day growing season, using get_curr_days_per_year()
      ! should give the same result of this function as using get_prev_days_per_year().
       DaysPastPlanting = jday - idop + get_curr_days_per_year()
    end if

  end function DaysPastPlanting

  !-----------------------------------------------------------------------
  subroutine vernalization(p, &
       canopystate_inst, temperature_inst, waterdiagnosticbulk_inst, cnveg_state_inst, crop_inst, &
       force_harvest)
    !
    ! !DESCRIPTION:
    !
    ! * * * only call for winter temperate cereal * * *
    !
    ! subroutine calculates vernalization and photoperiod effects on
    ! gdd accumulation in winter temperate cereal varieties. Thermal time accumulation
    ! is reduced in 1st period until plant is fully vernalized. During this
    ! time of emergence to spikelet formation, photoperiod can also have a
    ! drastic effect on plant development.
    !
    ! !ARGUMENTS:
    integer                , intent(in)    :: p    ! PATCH index running over
    type(canopystate_type) , intent(in)    :: canopystate_inst
    type(temperature_type) , intent(in)    :: temperature_inst
    type(waterdiagnosticbulk_type)  , intent(in)    :: waterdiagnosticbulk_inst
    type(cnveg_state_type) , intent(inout) :: cnveg_state_inst
    type(crop_type)        , intent(inout) :: crop_inst
    logical                , intent(inout) :: force_harvest ! "harvest" forced if freeze-killed
    !
    ! LOCAL VARAIBLES:
    real(r8) tcrown                     ! ?
    real(r8) vd, vd1, vd2               ! vernalization dependence
    real(r8) tkil                       ! Freeze kill threshold
    integer  c,g                        ! indices
    !------------------------------------------------------------------------

    associate(                                               & 
         tlai        => canopystate_inst%tlai_patch        , & ! Input:  [real(r8) (:) ]  one-sided leaf area index, no burying by snow     

         t_ref2m     => temperature_inst%t_ref2m_patch     , & ! Input:  [real(r8) (:) ]  2 m height surface air temperature (K)            
         t_ref2m_min => temperature_inst%t_ref2m_min_patch , & ! Input:  [real(r8) (:) ] daily minimum of average 2 m height surface air temperature (K)
         t_ref2m_max => temperature_inst%t_ref2m_max_patch , & ! Input:  [real(r8) (:) ] daily maximum of average 2 m height surface air temperature (K)

         snow_depth  => waterdiagnosticbulk_inst%snow_depth_col     , & ! Input:  [real(r8) (:) ]  snow height (m)                                   

         hdidx       => cnveg_state_inst%hdidx_patch       , & ! Output: [real(r8) (:) ]  cold hardening index?                             
         cumvd       => cnveg_state_inst%cumvd_patch       , & ! Output: [real(r8) (:) ]  cumulative vernalization d?ependence?             

         vf          => crop_inst%vf_patch                   & ! Output: [real(r8) (:) ]  vernalization factor for cereal
         )

      c = patch%column(p)

      ! for all equations - temperatures must be in degrees (C)
      ! calculate temperature of crown of crop (e.g., 3 cm soil temperature)
      ! snow depth in centimeters
      
      if (t_ref2m(p) < tfrz) then !slevis: t_ref2m inst of td=daily avg (K)
         tcrown = 2._r8 + (t_ref2m(p) - tfrz) * (0.4_r8 + 0.0018_r8 * &
              (min(snow_depth(c)*100._r8, 15._r8) - 15._r8)**2)
      else !slevis: snow_depth inst of adsnod=daily average (m)
         tcrown = t_ref2m(p) - tfrz
      end if

      ! vernalization factor calculation
      ! if vf(p) = 1.  then plant is fully vernalized - and thermal time
      ! accumulation in phase 1 will be unaffected
      ! refers to gddtsoi & hui, defined in the accumulation routines (slevis)
      ! reset vf, cumvd, and hdidx to 0 at planting of crop (slevis)

      if (t_ref2m_max(p) > tfrz) then
         if (t_ref2m_min(p) <= tfrz+15._r8) then
            vd1      = 1.4_r8 - 0.0778_r8 * tcrown
            vd2      = 0.5_r8 + 13.44_r8 / ((t_ref2m_max(p)-t_ref2m_min(p)+3._r8)**2) * tcrown
            vd       = max(0._r8, min(1._r8, vd1, vd2))
            cumvd(p) = cumvd(p) + vd
         end if

         if (cumvd(p) < 10._r8 .and. t_ref2m_max(p) > tfrz+30._r8) then
            cumvd(p) = cumvd(p) - 0.5_r8 * (t_ref2m_max(p) - tfrz - 30._r8)
         end if
         cumvd(p) = max(0._r8, cumvd(p))       ! must be > 0

         vf(p) = 1._r8 - p1v * (50._r8 - cumvd(p))
         vf(p) = max(0._r8, min(vf(p), 1._r8)) ! must be between 0 - 1
      end if

      ! calculate cold hardening of plant
      ! determines for winter cereal varieties whether the plant has completed
      ! a period of cold hardening to protect it from freezing temperatures. If
      ! not, then exposure could result in death or killing of plants.

      ! there are two distinct phases of hardening

      if (t_ref2m_min(p) <= tfrz-3._r8 .or. hdidx(p) /= 0._r8) then
         if (hdidx(p) >= hti) then   ! done with phase 1
            hdidx(p) = hdidx(p) + 0.083_r8
            hdidx(p) = min(hdidx(p), hti*2._r8)
         end if

         if (t_ref2m_max(p) >= tbase + tfrz + 10._r8) then
            hdidx(p) = hdidx(p) - 0.02_r8 * (t_ref2m_max(p)-tbase-tfrz-10._r8)
            if (hdidx(p) > hti) hdidx(p) = hdidx(p) - 0.02_r8 * (t_ref2m_max(p)-tbase-tfrz-10._r8)
            hdidx(p) = max(0._r8, hdidx(p))
         end if

      else if (tcrown >= tbase-1._r8) then
         if (tcrown <= tbase+8._r8) then
            hdidx(p) = hdidx(p) + 0.1_r8 - (tcrown-tbase+3.5_r8)**2 / 506._r8
            if (hdidx(p) >= hti .and. tcrown <= tbase + 0._r8) then
               hdidx(p) = hdidx(p) + 0.083_r8
               hdidx(p) = min(hdidx(p), hti*2._r8)
            end if
         end if

         if (t_ref2m_max(p) >= tbase + tfrz + 10._r8) then
            hdidx(p) = hdidx(p) - 0.02_r8 * (t_ref2m_max(p)-tbase-tfrz-10._r8)
            if (hdidx(p) > hti) hdidx(p) = hdidx(p) - 0.02_r8 * (t_ref2m_max(p)-tbase-tfrz-10._r8)
            hdidx(p) = max(0._r8, hdidx(p))
         end if
      end if

      ! calculate what the cereal killing temperature
      ! there is a linear inverse relationship between
      ! hardening of the plant and the killing temperature or
      ! threshold that the plant can withstand
      ! when plant is fully-hardened (hdidx = 2), the killing threshold is -18 C

      ! will have to develop some type of relationship that reduces LAI and
      ! biomass pools in response to cold damaged crop

      if (t_ref2m_min(p) <= tfrz - 6._r8) then
         tkil = (tbase - 6._r8) - 6._r8 * hdidx(p)
         if (tkil >= tcrown) then
            if ((0.95_r8 - 0.02_r8 * (tcrown - tkil)**2) >= 0.02_r8) then
               write (iulog,*)  'crop damaged by cold temperatures at p,c =', p,c
            else if (tlai(p) > 0._r8) then
               ! slevis: kill if past phase1 by forcing through harvest
               ! srabin: do this with force_harvest instead of setting
               !         gddmaturity = huigrain = 0, since gddmaturity==0 can
               !         lead to 0/0 when the crop isn't actually harvested based
               !         on "maturity." This can occur when generate_crop_gdds
               !         is true.
               force_harvest = .true.
               write (iulog,*)  '95% of crop killed by cold temperatures at p,c =', p,c
            end if
         end if
      end if

    end associate 

  end subroutine vernalization

  !-----------------------------------------------------------------------
  subroutine CNOnsetGrowth (num_soilp, filter_soilp, &
       cnveg_state_inst, &
       cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! Determines the flux of stored C and N from transfer pools to display
    ! pools during the phenological onset period.
    !
    ! !ARGUMENTS:
    integer                        , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                        , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_state_type)             , intent(in)    :: cnveg_state_inst
    type(cnveg_carbonstate_type)   , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type) , intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: p            ! indices
    integer :: fp           ! filter patch index
    real(r8):: t1           ! temporary variable
    !-----------------------------------------------------------------------

    associate(                                                                                             & 
         ivt                                 =>    patch%itype                                                   , & ! Input:  [integer   (:) ]  patch vegetation type                                

         woody                               =>    pftcon%woody                                                , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)
         
         onset_flag                          =>    cnveg_state_inst%onset_flag_patch                           , & ! Input:  [real(r8)  (:) ]  onset flag                                        
         onset_counter                       =>    cnveg_state_inst%onset_counter_patch                        , & ! Input:  [real(r8)  (:) ]  onset days counter                                
         bgtr                                =>    cnveg_state_inst%bgtr_patch                                 , & ! Input:  [real(r8)  (:) ]  background transfer growth rate (1/s)             
         
         leafc_xfer                          =>    cnveg_carbonstate_inst%leafc_xfer_patch                     , & ! Input:  [real(r8)  (:) ]  (gC/m2) leaf C transfer                           
         frootc_xfer                         =>    cnveg_carbonstate_inst%frootc_xfer_patch                    , & ! Input:  [real(r8)  (:) ]  (gC/m2) fine root C transfer                      
         livestemc_xfer                      =>    cnveg_carbonstate_inst%livestemc_xfer_patch                 , & ! Input:  [real(r8)  (:) ]  (gC/m2) live stem C transfer                      
         deadstemc_xfer                      =>    cnveg_carbonstate_inst%deadstemc_xfer_patch                 , & ! Input:  [real(r8)  (:) ]  (gC/m2) dead stem C transfer                      
         livecrootc_xfer                     =>    cnveg_carbonstate_inst%livecrootc_xfer_patch                , & ! Input:  [real(r8)  (:) ]  (gC/m2) live coarse root C transfer               
         deadcrootc_xfer                     =>    cnveg_carbonstate_inst%deadcrootc_xfer_patch                , & ! Input:  [real(r8)  (:) ]  (gC/m2) dead coarse root C transfer               
         
         leafn_xfer                          =>    cnveg_nitrogenstate_inst%leafn_xfer_patch                   , & ! Input:  [real(r8)  (:) ]  (gN/m2) leaf N transfer                           
         frootn_xfer                         =>    cnveg_nitrogenstate_inst%frootn_xfer_patch                  , & ! Input:  [real(r8)  (:) ]  (gN/m2) fine root N transfer                      
         livestemn_xfer                      =>    cnveg_nitrogenstate_inst%livestemn_xfer_patch               , & ! Input:  [real(r8)  (:) ]  (gN/m2) live stem N transfer                      
         deadstemn_xfer                      =>    cnveg_nitrogenstate_inst%deadstemn_xfer_patch               , & ! Input:  [real(r8)  (:) ]  (gN/m2) dead stem N transfer                      
         livecrootn_xfer                     =>    cnveg_nitrogenstate_inst%livecrootn_xfer_patch              , & ! Input:  [real(r8)  (:) ]  (gN/m2) live coarse root N transfer               
         deadcrootn_xfer                     =>    cnveg_nitrogenstate_inst%deadcrootn_xfer_patch              , & ! Input:  [real(r8)  (:) ]  (gN/m2) dead coarse root N transfer               
         
         leafc_xfer_to_leafc                 =>    cnveg_carbonflux_inst%leafc_xfer_to_leafc_patch             , & ! Output:  [real(r8) (:) ]                                                    
         frootc_xfer_to_frootc               =>    cnveg_carbonflux_inst%frootc_xfer_to_frootc_patch           , & ! Output:  [real(r8) (:) ]                                                    
         livestemc_xfer_to_livestemc         =>    cnveg_carbonflux_inst%livestemc_xfer_to_livestemc_patch     , & ! Output:  [real(r8) (:) ]                                                    
         deadstemc_xfer_to_deadstemc         =>    cnveg_carbonflux_inst%deadstemc_xfer_to_deadstemc_patch     , & ! Output:  [real(r8) (:) ]                                                    
         livecrootc_xfer_to_livecrootc       =>    cnveg_carbonflux_inst%livecrootc_xfer_to_livecrootc_patch   , & ! Output:  [real(r8) (:) ]                                                    
         deadcrootc_xfer_to_deadcrootc       =>    cnveg_carbonflux_inst%deadcrootc_xfer_to_deadcrootc_patch   , & ! Output:  [real(r8) (:) ]                                                    
         
         leafn_xfer_to_leafn                 =>    cnveg_nitrogenflux_inst%leafn_xfer_to_leafn_patch           , & ! Output:  [real(r8) (:) ]                                                    
         frootn_xfer_to_frootn               =>    cnveg_nitrogenflux_inst%frootn_xfer_to_frootn_patch         , & ! Output:  [real(r8) (:) ]                                                    
         livestemn_xfer_to_livestemn         =>    cnveg_nitrogenflux_inst%livestemn_xfer_to_livestemn_patch   , & ! Output:  [real(r8) (:) ]                                                    
         deadstemn_xfer_to_deadstemn         =>    cnveg_nitrogenflux_inst%deadstemn_xfer_to_deadstemn_patch   , & ! Output:  [real(r8) (:) ]                                                    
         livecrootn_xfer_to_livecrootn       =>    cnveg_nitrogenflux_inst%livecrootn_xfer_to_livecrootn_patch , & ! Output:  [real(r8) (:) ]                                                    
         deadcrootn_xfer_to_deadcrootn       =>    cnveg_nitrogenflux_inst%deadcrootn_xfer_to_deadcrootn_patch , & ! Output:  [real(r8) (:) ]
         ileafst_to_ileafxf_phc              =>    cnveg_carbonflux_inst%ileafst_to_ileafxf_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from leaf storage pool to leaf transfer pool
         ileafxf_to_ileaf_phc                =>    cnveg_carbonflux_inst%ileafxf_to_ileaf_ph                   , & ! Input: [integer (:)] Index of phenology related C transfer from leaf transfer pool to leaf pool
         ifrootst_to_ifrootxf_phc            =>    cnveg_carbonflux_inst%ifrootst_to_ifrootxf_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from fine root storage pool to fine root transfer pool
         ifrootxf_to_ifroot_phc              =>    cnveg_carbonflux_inst%ifrootxf_to_ifroot_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from fine root transfer pool to fine root pool
         ilivestemst_to_ilivestemxf_phc      =>    cnveg_carbonflux_inst%ilivestemst_to_ilivestemxf_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live stem storage pool to live stem transfer pool
         ilivestemxf_to_ilivestem_phc        =>    cnveg_carbonflux_inst%ilivestemxf_to_ilivestem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem transfer pool to live stem pool
         ideadstemst_to_ideadstemxf_phc      =>    cnveg_carbonflux_inst%ideadstemst_to_ideadstemxf_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem storage pool to dead stem transfer pool
         ideadstemxf_to_ideadstem_phc        =>    cnveg_carbonflux_inst%ideadstemxf_to_ideadstem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem transfer pool to dead stem pool
         ilivecrootst_to_ilivecrootxf_phc    =>    cnveg_carbonflux_inst%ilivecrootst_to_ilivecrootxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root storage pool to live coarse root transfer pool
         ilivecrootxf_to_ilivecroot_phc      =>    cnveg_carbonflux_inst%ilivecrootxf_to_ilivecroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root transfer pool to live coarse root pool
         ideadcrootst_to_ideadcrootxf_phc    =>    cnveg_carbonflux_inst%ideadcrootst_to_ideadcrootxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root storage pool to dead coarse root transfer pool
         ideadcrootxf_to_ideadcroot_phc      =>    cnveg_carbonflux_inst%ideadcrootxf_to_ideadcroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root transfer pool to dead coarse root pool
         ilivestem_to_ideadstem_phc          =>    cnveg_carbonflux_inst%ilivestem_to_ideadstem_ph             , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to dead stem pool
         ilivecroot_to_ideadcroot_phc        =>    cnveg_carbonflux_inst%ilivecroot_to_ideadcroot_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root to dead coarse root pool
         ileaf_to_iout_phc                   =>    cnveg_carbonflux_inst%ileaf_to_iout_ph                      , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to outside of vegetation pools
         ifroot_to_iout_phc                  =>    cnveg_carbonflux_inst%ifroot_to_iout_ph                     , & ! Input: [integer (:)] Index of phenology related C transfer from fine root pool to outside of vegetation pools
         ilivestem_to_iout_phc               =>    cnveg_carbonflux_inst%ilivestem_to_iout_ph                  , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to outside of vegetation pools
         igrain_to_iout_phc                  =>    cnveg_carbonflux_inst%igrain_to_iout_ph                     , & ! Input: [integer (:)] Index of phenology related C transfer from grain pool to outside of vegetation pools
         ileafst_to_ileafxf_phn              =>    cnveg_nitrogenflux_inst%ileafst_to_ileafxf_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from leaf storage pool to leaf transfer pool
         ileafxf_to_ileaf_phn                =>    cnveg_nitrogenflux_inst%ileafxf_to_ileaf_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from leaf transfer pool to leaf pool
         ifrootst_to_ifrootxf_phn            =>    cnveg_nitrogenflux_inst%ifrootst_to_ifrootxf_ph             , & ! Input: [integer (:)] Index of phenology related C transfer from fine root storage pool to fine root transfer pool
         ifrootxf_to_ifroot_phn              =>    cnveg_nitrogenflux_inst%ifrootxf_to_ifroot_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from fine root transfer pool to fine root pool
         ilivestemst_to_ilivestemxf_phn      =>    cnveg_nitrogenflux_inst%ilivestemst_to_ilivestemxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live stem storage pool to live stem transfer pool
         ilivestemxf_to_ilivestem_phn        =>    cnveg_nitrogenflux_inst%ilivestemxf_to_ilivestem_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live stem transfer pool to live stem pool
         ideadstemst_to_ideadstemxf_phn      =>    cnveg_nitrogenflux_inst%ideadstemst_to_ideadstemxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem storage pool to dead stem transfer pool
         ideadstemxf_to_ideadstem_phn        =>    cnveg_nitrogenflux_inst%ideadstemxf_to_ideadstem_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem transfer pool to dead stem pool
         ilivecrootst_to_ilivecrootxf_phn    =>    cnveg_nitrogenflux_inst%ilivecrootst_to_ilivecrootxf_ph     , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root storage pool to live coarse root transfer pool
         ilivecrootxf_to_ilivecroot_phn      =>    cnveg_nitrogenflux_inst%ilivecrootxf_to_ilivecroot_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root transfer pool to live coarse root pool
         ideadcrootst_to_ideadcrootxf_phn    =>    cnveg_nitrogenflux_inst%ideadcrootst_to_ideadcrootxf_ph     , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root storage pool to dead coarse root transfer pool
         ideadcrootxf_to_ideadcroot_phn      =>    cnveg_nitrogenflux_inst%ideadcrootxf_to_ideadcroot_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root transfer pool to dead coarse root pool
         ilivestem_to_ideadstem_phn          =>    cnveg_nitrogenflux_inst%ilivestem_to_ideadstem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to dead stem pool
         ilivecroot_to_ideadcroot_phn        =>    cnveg_nitrogenflux_inst%ilivecroot_to_ideadcroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root to dead coarse root pool
         ileaf_to_iout_phn                   =>    cnveg_nitrogenflux_inst%ileaf_to_iout_ph                    , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to outside of vegetation pools
         ifroot_to_iout_phn                  =>    cnveg_nitrogenflux_inst%ifroot_to_iout_ph                   , & ! Input: [integer (:)] Index of phenology related C transfer from fine root pool to outside of vegetation pools
         ilivestem_to_iout_phn               =>    cnveg_nitrogenflux_inst%ilivestem_to_iout_ph                , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to outside of vegetation pools
         ileaf_to_iretransn_phn              =>    cnveg_nitrogenflux_inst%ileaf_to_iretransn_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to retranslocation pool
         ilivestem_to_iretransn_phn          =>    cnveg_nitrogenflux_inst%ilivestem_to_iretransn_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to retranslocation pool
         ilivecroot_to_iretransn_phn         =>    cnveg_nitrogenflux_inst%ilivecroot_to_iretransn_ph          , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root pool to retranslocation pool
         igrain_to_iout_phn                  =>    cnveg_nitrogenflux_inst%igrain_to_iout_ph                     & ! Input: [integer (:)] Index of phenology related C transfer from grain pool to outside of vegetation pools
         )

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! only calculate these fluxes during onset period
         if (onset_flag(p) == 1._r8) then

            ! The transfer rate is a linearly decreasing function of time,
            ! going to zero on the last timestep of the onset period

            if (abs(onset_counter(p) - dt) <= dt/2._r8) then
               t1 = 1.0_r8 / dt
            else
               t1 = 2.0_r8 / (onset_counter(p))
            end if
            if (use_matrixcn)then
               leafc_xfer_to_leafc(p)   = leafc_xfer(p) * matrix_update_phc(p,ileafxf_to_ileaf_phc,t1,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
               frootc_xfer_to_frootc(p) = frootc_xfer(p) * matrix_update_phc(p,ifrootxf_to_ifroot_phc,t1,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
               leafn_xfer_to_leafn(p)   = leafn_xfer(p) * matrix_update_phn(p,ileafxf_to_ileaf_phn,t1,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
               frootn_xfer_to_frootn(p) = frootn_xfer(p) * matrix_update_phn(p,ifrootxf_to_ifroot_phn,t1,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
               if (woody(ivt(p)) == 1.0_r8) then

                  livestemc_xfer_to_livestemc(p)   = livestemc_xfer(p)  * matrix_update_phc(p,ilivestemxf_to_ilivestem_phc,t1,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                  deadstemc_xfer_to_deadstemc(p)   = deadstemc_xfer(p)  * matrix_update_phc(p,ideadstemxf_to_ideadstem_phc,t1,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                  livecrootc_xfer_to_livecrootc(p) = livecrootc_xfer(p) * matrix_update_phc(p,ilivecrootxf_to_ilivecroot_phc,t1,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                  deadcrootc_xfer_to_deadcrootc(p) = deadcrootc_xfer(p) * matrix_update_phc(p,ideadcrootxf_to_ideadcroot_phc,t1,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)

                  livestemn_xfer_to_livestemn(p)   = livestemn_xfer(p)  * matrix_update_phn(p,ilivestemxf_to_ilivestem_phn,t1,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                  deadstemn_xfer_to_deadstemn(p)   = deadstemn_xfer(p)  * matrix_update_phn(p,ideadstemxf_to_ideadstem_phn,t1,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                  livecrootn_xfer_to_livecrootn(p) = livecrootn_xfer(p) * matrix_update_phn(p,ilivecrootxf_to_ilivecroot_phn,t1,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                  deadcrootn_xfer_to_deadcrootn(p) = deadcrootn_xfer(p) * matrix_update_phn(p,ideadcrootxf_to_ideadcroot_phn,t1,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
               end if
            else
               ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
               !                                        and CNNStateUpdate1::NStateUpdate1
               leafc_xfer_to_leafc(p)   = t1 * leafc_xfer(p)
               frootc_xfer_to_frootc(p) = t1 * frootc_xfer(p)
               leafn_xfer_to_leafn(p)   = t1 * leafn_xfer(p)
               frootn_xfer_to_frootn(p) = t1 * frootn_xfer(p)
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemc_xfer_to_livestemc(p)   = t1 * livestemc_xfer(p)
                  deadstemc_xfer_to_deadstemc(p)   = t1 * deadstemc_xfer(p)
                  livecrootc_xfer_to_livecrootc(p) = t1 * livecrootc_xfer(p)
                  deadcrootc_xfer_to_deadcrootc(p) = t1 * deadcrootc_xfer(p)
                  livestemn_xfer_to_livestemn(p)   = t1 * livestemn_xfer(p)
                  deadstemn_xfer_to_deadstemn(p)   = t1 * deadstemn_xfer(p)
                  livecrootn_xfer_to_livecrootn(p) = t1 * livecrootn_xfer(p)
                  deadcrootn_xfer_to_deadcrootn(p) = t1 * deadcrootn_xfer(p)
               end if
            end if !use_matrixcn

         end if ! end if onset period

         ! calculate the background rate of transfer growth (used for stress
         ! deciduous algorithm). In this case, all of the mass in the transfer
         ! pools should be moved to displayed growth in each timestep.

         if (bgtr(p) > 0._r8) then
            if(use_matrixcn)then
               leafc_xfer_to_leafc(p)   = leafc_xfer(p)  * matrix_update_phc(p,ileafxf_to_ileaf_phc,1._r8 / dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
               frootc_xfer_to_frootc(p) = frootc_xfer(p) * matrix_update_phc(p,ifrootxf_to_ifroot_phc,1._r8 / dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
               leafn_xfer_to_leafn(p)   = leafn_xfer(p)  * matrix_update_phn(p,ileafxf_to_ileaf_phn,1._r8 / dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
               frootn_xfer_to_frootn(p) = frootn_xfer(p) * matrix_update_phn(p,ifrootxf_to_ifroot_phn,1._r8 / dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
               if (woody(ivt(p)) == 1.0_r8) then

                  livestemc_xfer_to_livestemc(p)   = livestemc_xfer(p)  * matrix_update_phc(p,ilivestemxf_to_ilivestem_phc,1._r8 / dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                  deadstemc_xfer_to_deadstemc(p)   = deadstemc_xfer(p)  * matrix_update_phc(p,ideadstemxf_to_ideadstem_phc,1._r8 / dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                  livecrootc_xfer_to_livecrootc(p) = livecrootc_xfer(p)  * matrix_update_phc(p,ilivecrootxf_to_ilivecroot_phc,1._r8 / dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                  deadcrootc_xfer_to_deadcrootc(p) = deadcrootc_xfer(p)  * matrix_update_phc(p,ideadcrootxf_to_ideadcroot_phc,1._r8 / dt,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)

                  livestemn_xfer_to_livestemn(p)   = livestemn_xfer(p)  * matrix_update_phn(p,ilivestemxf_to_ilivestem_phn,1._r8 / dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                  deadstemn_xfer_to_deadstemn(p)   = deadstemn_xfer(p)  * matrix_update_phn(p,ideadstemxf_to_ideadstem_phn,1._r8 / dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                  livecrootn_xfer_to_livecrootn(p) = livecrootn_xfer(p)  * matrix_update_phn(p,ilivecrootxf_to_ilivecroot_phn,1._r8 / dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                  deadcrootn_xfer_to_deadcrootn(p) = deadcrootn_xfer(p)  * matrix_update_phn(p,ideadcrootxf_to_ideadcroot_phn,1._r8 / dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
               end if
            else
               ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
               !                                        and CNNStateUpdate1::NStateUpdate1
               leafc_xfer_to_leafc(p)   = leafc_xfer(p) / dt
               frootc_xfer_to_frootc(p) = frootc_xfer(p) / dt
               leafn_xfer_to_leafn(p)   = leafn_xfer(p) / dt
               frootn_xfer_to_frootn(p) = frootn_xfer(p) / dt
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemc_xfer_to_livestemc(p)   = livestemc_xfer(p) / dt
                  deadstemc_xfer_to_deadstemc(p)   = deadstemc_xfer(p) / dt
                  livecrootc_xfer_to_livecrootc(p) = livecrootc_xfer(p) / dt
                  deadcrootc_xfer_to_deadcrootc(p) = deadcrootc_xfer(p) / dt
                  livestemn_xfer_to_livestemn(p)   = livestemn_xfer(p) / dt
                  deadstemn_xfer_to_deadstemn(p)   = deadstemn_xfer(p) / dt
                  livecrootn_xfer_to_livecrootn(p) = livecrootn_xfer(p) / dt
                  deadcrootn_xfer_to_deadcrootn(p) = deadcrootn_xfer(p) / dt
               end if
            end if !use_matrixcn
         end if ! end if bgtr

      end do ! end patch loop

    end associate

  end subroutine CNOnsetGrowth

  !-----------------------------------------------------------------------
  subroutine CNOffsetLitterfall (num_soilp, filter_soilp, &
       cnveg_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
       crop_inst)
    !
    ! !DESCRIPTION:
    ! Determines the flux of C and N from displayed pools to litter
    ! pools during the phenological offset period.
    !
    ! !USES:
    use pftconMod        , only : npcropmin
    use pftconMod        , only : nmiscanthus, nirrig_miscanthus, nswitchgrass, nirrig_switchgrass
    
    use CNSharedParamsMod, only : use_fun
    use clm_varctl       , only : CNratio_floating, crop_residue_removal_frac
    !
    ! !ARGUMENTS:
    integer                       , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                       , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_state_type)        , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)  , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type), intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)   , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type) , intent(inout) :: cnveg_nitrogenflux_inst
    type(crop_type)               , intent(in)    :: crop_inst
    !
    ! !LOCAL VARIABLES:
    integer :: p, c, k, h   ! indices
    integer :: fp           ! filter patch index
    real(r8):: t1           ! temporary variable
    real(r8):: denom        ! temporary variable for divisor
    real(r8) :: ntovr_leaf  
    real(r8) :: fr_leafn_to_litter ! fraction of the nitrogen turnover that goes to litter; remaining fraction is retranslocated
    real(r8) :: grainc_to_out, grainn_to_out ! Temporary for grain Carbon and grain Nitrogen output
    real(r8) :: cropseedc_deficit_remaining  ! remaining amount of crop seed C deficit that still needs to be restored (gC/m2) (positive, in contrast to the negative cropseedc_deficit)
    real(r8) :: cropseedn_deficit_remaining  ! remaining amount of crop seed N deficit that still needs to be restored (gN/m2) (positive, in contrast to the negative cropseedn_deficit)
    real(r8) :: cropseedc_deficit_to_restore ! amount of crop seed C deficit that will be restored from this grain pool (gC/m2)
    real(r8) :: cropseedn_deficit_to_restore ! amount of crop seed N deficit that will be restored from this grain pool (gN/m2)
    real(r8) :: repr_grainc_to_food_thispool ! amount added to / subtracted from repr_grainc_to_food for the pool in question (gC/m2/s)
    real(r8) :: repr_grainn_to_food_thispool ! amount added to / subtracted from repr_grainn_to_food for the pool in question (gN/m2/s)
    real(r8) :: leafc_remaining, livestemc_remaining
    real(r8) :: leafn_remaining, livestemn_remaining
    real(r8) :: removedresidue_fraction
    !-----------------------------------------------------------------------

    associate(                                                                           & 
         ivt                   =>    patch%itype                                       , & ! Input:  [integer  (:) ]  patch vegetation type                                

         leafcn                =>    pftcon%leafcn                                     , & ! Input:  leaf C:N (gC/gN) 
         
         biofuel_harvfrac      =>    pftcon%biofuel_harvfrac                           , & ! Input:  cut a fraction of leaf & stem for biofuel (-)
         repr_structure_harvfrac =>  pftcon%repr_structure_harvfrac                    , & ! Input:  fraction of each reproductive structure component that is harvested and sent to the crop products pool

         lflitcn               =>    pftcon%lflitcn                                    , & ! Input:  leaf litter C:N (gC/gN)                           
         frootcn               =>    pftcon%frootcn                                    , & ! Input:  fine root C:N (gC/gN)                             
         graincn               =>    pftcon%graincn                                    , & ! Input:  grain C:N (gC/gN)                                 

         offset_flag           =>    cnveg_state_inst%offset_flag_patch                , & ! Input:  [real(r8) (:) ]  offset flag                                       
         offset_counter        =>    cnveg_state_inst%offset_counter_patch             , & ! Input:  [real(r8) (:) ]  offset days counter                               

         leafc                 =>    cnveg_carbonstate_inst%leafc_patch                , & ! Input:  [real(r8) (:) ]  (gC/m2) leaf C                                    
         frootc                =>    cnveg_carbonstate_inst%frootc_patch               , & ! Input:  [real(r8) (:) ]  (gC/m2) fine root C                               
         reproductivec                =>    cnveg_carbonstate_inst%reproductivec_patch               , & ! Input:  [real(r8) (:,:) ]  (gC/m2) grain C
         cropseedc_deficit     =>    cnveg_carbonstate_inst%cropseedc_deficit_patch    , & ! Input:  [real(r8) (:) ]  (gC/m2) crop seed C deficit
         livestemc             =>    cnveg_carbonstate_inst%livestemc_patch            , & ! Input:  [real(r8) (:) ]  (gC/m2) livestem C                                
         cropseedn_deficit     =>    cnveg_nitrogenstate_inst%cropseedn_deficit_patch  , & ! Input:  [real(r8) (:) ]  (gC/m2) crop seed N deficit
         livestemn             =>    cnveg_nitrogenstate_inst%livestemn_patch          , & ! Input:  [real(r8) (:) ]  (gN/m2) livestem N

         cpool_to_reproductivec       =>    cnveg_carbonflux_inst%cpool_to_reproductivec_patch       , & ! Input:  [real(r8) (:,:) ]  allocation to grain C (gC/m2/s)
         npool_to_reproductiven       =>    cnveg_nitrogenflux_inst%npool_to_reproductiven_patch     , & ! Input: [real(r8) (:,:)  ]  allocation to grain N (gN/m2/s)
         reproductiven                =>    cnveg_nitrogenstate_inst%reproductiven_patch             , & ! Input: [real(r8) (:,:)  ]  (kgN/m2) grain N
         cpool_to_livestemc    =>    cnveg_carbonflux_inst%cpool_to_livestemc_patch    , & ! Input:  [real(r8) (:) ]  allocation to live stem C (gC/m2/s)               
         cpool_to_leafc        =>    cnveg_carbonflux_inst%cpool_to_leafc_patch        , & ! Input:  [real(r8) (:) ]  allocation to leaf C (gC/m2/s)                    
         cpool_to_frootc       =>    cnveg_carbonflux_inst%cpool_to_frootc_patch       , & ! Input:  [real(r8) (:) ]  allocation to fine root C (gC/m2/s)               
         prev_leafc_to_litter  =>    cnveg_carbonflux_inst%prev_leafc_to_litter_patch  , & ! Output: [real(r8) (:) ]  previous timestep leaf C litterfall flux (gC/m2/s)
         prev_frootc_to_litter =>    cnveg_carbonflux_inst%prev_frootc_to_litter_patch , & ! Output: [real(r8) (:) ]  previous timestep froot C litterfall flux (gC/m2/s)
         leafc_to_litter       =>    cnveg_carbonflux_inst%leafc_to_litter_patch       , & ! Output: [real(r8) (:) ]  leaf C litterfall (gC/m2/s)                       
         frootc_to_litter      =>    cnveg_carbonflux_inst%frootc_to_litter_patch      , & ! Output: [real(r8) (:) ]  fine root C litterfall (gC/m2/s)                  
         livestemc_to_litter   =>    cnveg_carbonflux_inst%livestemc_to_litter_patch   , & ! Output: [real(r8) (:) ]  live stem C litterfall (gC/m2/s)                  
         repr_grainc_to_food   =>    cnveg_carbonflux_inst%repr_grainc_to_food_patch   , & ! Output: [real(r8) (:,:) ]  grain C to food (gC/m2/s)
         repr_grainc_to_food_perharv => cnveg_carbonflux_inst%repr_grainc_to_food_perharv_patch, & ! Output: [real(r8) (:,:,:) ]  grain C to food per harvest (gC/m2)
         repr_grainc_to_food_thisyr => cnveg_carbonflux_inst%repr_grainc_to_food_thisyr_patch, & ! Output: [real(r8) (:,:) ]  grain C to food harvested this calendar year (gC/m2)
         repr_grainc_to_seed   =>    cnveg_carbonflux_inst%repr_grainc_to_seed_patch   , & ! Output: [real(r8) (:,:) ]  grain C to seed (gC/m2/s)
         repr_grainc_to_seed_perharv => cnveg_carbonflux_inst%repr_grainc_to_seed_perharv_patch, & ! Output: [real(r8) (:,:,:) ]  grain C to seed per harvest (gC/m2)
         repr_grainc_to_seed_thisyr => cnveg_carbonflux_inst%repr_grainc_to_seed_thisyr_patch, & ! Output: [real(r8) (:,:) ]  grain C to seed harvested this calendar year (gC/m2)
         repr_structurec_to_cropprod => cnveg_carbonflux_inst%repr_structurec_to_cropprod_patch, & ! Output: [real(r8) (:,:) ] reproductive structure C to crop product pool (gC/m2/s)
         repr_structurec_to_litter   => cnveg_carbonflux_inst%repr_structurec_to_litter_patch,   & ! Output: [real(r8) (:,:) ] reproductive structure C to litter (gC/m2/s)
         leafc_to_biofuelc     =>    cnveg_carbonflux_inst%leafc_to_biofuelc_patch     , & ! Output: [real(r8) (:) ]  leaf C to biofuel C (gC/m2/s)
         livestemc_to_biofuelc =>    cnveg_carbonflux_inst%livestemc_to_biofuelc_patch , & ! Output: [real(r8) (:) ]  livestem C to biofuel C (gC/m2/s)
         leafc_to_removedresiduec     => cnveg_carbonflux_inst%leafc_to_removedresiduec_patch     , & ! Output: [real(r8) (:) ]  leaf C to removed residue C (gC/m2/s)
         livestemc_to_removedresiduec => cnveg_carbonflux_inst%livestemc_to_removedresiduec_patch , & ! Output: [real(r8) (:) ]  livestem C to removed residue C (gC/m2/s)
         leafn                 =>    cnveg_nitrogenstate_inst%leafn_patch              , & ! Input:  [real(r8) (:) ]  (gN/m2) leaf N      
         frootn                =>    cnveg_nitrogenstate_inst%frootn_patch             , & ! Input:  [real(r8) (:) ]  (gN/m2) fine root N                        

         livestemn_to_litter   =>    cnveg_nitrogenflux_inst%livestemn_to_litter_patch , & ! Output: [real(r8) (:) ]  livestem N to litter (gN/m2/s)                    
         repr_grainn_to_food        =>    cnveg_nitrogenflux_inst%repr_grainn_to_food_patch      , & ! Output: [real(r8) (:,:) ]  grain N to food (gN/m2/s)
         repr_grainn_to_food_perharv => cnveg_nitrogenflux_inst%repr_grainn_to_food_perharv_patch, & ! Output: [real(r8) (:,:,:) ]  grain N to food per harvest (gN/m2)
         repr_grainn_to_food_thisyr => cnveg_nitrogenflux_inst%repr_grainn_to_food_thisyr_patch, & ! Output: [real(r8) (:,:) ]  grain N to food harvested this calendar year (gN/m2)
         repr_grainn_to_seed        =>    cnveg_nitrogenflux_inst%repr_grainn_to_seed_patch      , & ! Output: [real(r8) (:,:) ]  grain N to seed (gN/m2/s)
         repr_grainn_to_seed_perharv => cnveg_nitrogenflux_inst%repr_grainn_to_seed_perharv_patch, & ! Output: [real(r8) (:,:,:) ]  grain N to seed per harvest (gN/m2)
         repr_grainn_to_seed_thisyr => cnveg_nitrogenflux_inst%repr_grainn_to_seed_thisyr_patch, & ! Output: [real(r8) (:,:) ]  grain N to seed harvested this calendar year (gN/m2)
         repr_structuren_to_cropprod => cnveg_nitrogenflux_inst%repr_structuren_to_cropprod_patch, & ! Output: [real(r8) (:,:) ] reproductive structure N to crop product pool (gN/m2/s)
         repr_structuren_to_litter   => cnveg_nitrogenflux_inst%repr_structuren_to_litter_patch,   & ! Output: [real(r8) (:,:) ] reproductive structure N to litter (gN/m2/s)
         leafn_to_biofueln     =>    cnveg_nitrogenflux_inst%leafn_to_biofueln_patch   , & ! Output: [real(r8) (:) ]  leaf N to biofuel N (gN/m2/s)
         livestemn_to_biofueln =>    cnveg_nitrogenflux_inst%livestemn_to_biofueln_patch, & ! Output: [real(r8) (:) ]  livestem N to biofuel N (gN/m2/s)     
         leafn_to_removedresiduen     => cnveg_nitrogenflux_inst%leafn_to_removedresiduen_patch    , & ! Output: [real(r8) (:) ]  leaf N to removed residue N (gN/m2/s)
         livestemn_to_removedresiduen => cnveg_nitrogenflux_inst%livestemn_to_removedresiduen_patch, & ! Output: [real(r8) (:) ]  livestem N to removed residue N (gN/m2/s)     
         leafn_to_litter       =>    cnveg_nitrogenflux_inst%leafn_to_litter_patch     , & ! Output: [real(r8) (:) ]  leaf N litterfall (gN/m2/s)                       
         leafn_to_retransn     =>    cnveg_nitrogenflux_inst%leafn_to_retransn_patch   , & ! Input: [real(r8) (:) ]  leaf N to retranslocated N pool (gN/m2/s)         
         free_retransn_to_npool=>    cnveg_nitrogenflux_inst%free_retransn_to_npool_patch  , & ! Input: [real(r8) (:) ] free leaf N to retranslocated N pool (gN/m2/s)          
         paid_retransn_to_npool=>    cnveg_nitrogenflux_inst%retransn_to_npool_patch, & ! Input: [real(r8) (:) ] free leaf N to retranslocated N pool (gN/m2/s)          
         frootn_to_litter      =>    cnveg_nitrogenflux_inst%frootn_to_litter_patch    , & ! Output: [real(r8) (:) ]  fine root N litterfall (gN/m2/s)                  
         leafc_to_litter_fun   =>    cnveg_carbonflux_inst%leafc_to_litter_fun_patch   , & ! Output:  [real(r8) (:) ]  leaf C litterfall used by FUN (gC/m2/s)
         leafcn_offset         =>    cnveg_state_inst%leafcn_offset_patch              , & ! Output:  [real(r8) (:) ]  Leaf C:N used by FUN

         ileafst_to_ileafxf_phc              =>    cnveg_carbonflux_inst%ileafst_to_ileafxf_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from leaf storage pool to leaf transfer pool
         ileafxf_to_ileaf_phc                =>    cnveg_carbonflux_inst%ileafxf_to_ileaf_ph                   , & ! Input: [integer (:)] Index of phenology related C transfer from leaf transfer pool to leaf pool
         ifrootst_to_ifrootxf_phc            =>    cnveg_carbonflux_inst%ifrootst_to_ifrootxf_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from fine root storage pool to fine root transfer pool
         ifrootxf_to_ifroot_phc              =>    cnveg_carbonflux_inst%ifrootxf_to_ifroot_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from fine root transfer pool to fine root pool
         ilivestemst_to_ilivestemxf_phc      =>    cnveg_carbonflux_inst%ilivestemst_to_ilivestemxf_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live stem storage pool to live stem transfer pool
         ilivestemxf_to_ilivestem_phc        =>    cnveg_carbonflux_inst%ilivestemxf_to_ilivestem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem transfer pool to live stem pool
         ideadstemst_to_ideadstemxf_phc      =>    cnveg_carbonflux_inst%ideadstemst_to_ideadstemxf_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem storage pool to dead stem transfer pool
         ideadstemxf_to_ideadstem_phc        =>    cnveg_carbonflux_inst%ideadstemxf_to_ideadstem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem transfer pool to dead stem pool
         ilivecrootst_to_ilivecrootxf_phc    =>    cnveg_carbonflux_inst%ilivecrootst_to_ilivecrootxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root storage pool to live coarse root transfer pool
         ilivecrootxf_to_ilivecroot_phc      =>    cnveg_carbonflux_inst%ilivecrootxf_to_ilivecroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root transfer pool to live coarse root pool
         ideadcrootst_to_ideadcrootxf_phc    =>    cnveg_carbonflux_inst%ideadcrootst_to_ideadcrootxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root storage pool to dead coarse root transfer pool
         ideadcrootxf_to_ideadcroot_phc      =>    cnveg_carbonflux_inst%ideadcrootxf_to_ideadcroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root transfer pool to dead coarse root pool
         ilivestem_to_ideadstem_phc          =>    cnveg_carbonflux_inst%ilivestem_to_ideadstem_ph             , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to dead stem pool
         ilivecroot_to_ideadcroot_phc        =>    cnveg_carbonflux_inst%ilivecroot_to_ideadcroot_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root to dead coarse root pool
         ileaf_to_iout_phc                   =>    cnveg_carbonflux_inst%ileaf_to_iout_ph                      , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to outside of vegetation pools
         ileaf_to_iout_gmc                   =>    cnveg_carbonflux_inst%ileaf_to_iout_gm                      , & ! Input: [integer (:)] Index of gap mortality related C transfer from leaf pool to outside of vegetation pools
         ileaf_to_iout_gmn                   =>    cnveg_nitrogenflux_inst%ileaf_to_iout_gm                    , & ! Input: [integer (:)] Index of gap mortality related N transfer from leaf pool to outside of vegetation pools
         ifroot_to_iout_phc                  =>    cnveg_carbonflux_inst%ifroot_to_iout_ph                     , & ! Input: [integer (:)] Index of phenology related C transfer from fine root pool to outside of vegetation pools
         ilivestem_to_iout_phc               =>    cnveg_carbonflux_inst%ilivestem_to_iout_ph                  , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to outside of vegetation pools
         igrain_to_iout_phc                  =>    cnveg_carbonflux_inst%igrain_to_iout_ph                     , & ! Input: [integer (:)] Index of phenology related C transfer from grain pool to outside of vegetation pools
         ileafst_to_ileafxf_phn              =>    cnveg_nitrogenflux_inst%ileafst_to_ileafxf_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from leaf storage pool to leaf transfer pool
         ileafxf_to_ileaf_phn                =>    cnveg_nitrogenflux_inst%ileafxf_to_ileaf_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from leaf transfer pool to leaf pool
         ifrootst_to_ifrootxf_phn            =>    cnveg_nitrogenflux_inst%ifrootst_to_ifrootxf_ph             , & ! Input: [integer (:)] Index of phenology related C transfer from fine root storage pool to fine root transfer pool
         ifrootxf_to_ifroot_phn              =>    cnveg_nitrogenflux_inst%ifrootxf_to_ifroot_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from fine root transfer pool to fine root pool
         ilivestemst_to_ilivestemxf_phn      =>    cnveg_nitrogenflux_inst%ilivestemst_to_ilivestemxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live stem storage pool to live stem transfer pool
         ilivestemxf_to_ilivestem_phn        =>    cnveg_nitrogenflux_inst%ilivestemxf_to_ilivestem_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live stem transfer pool to live stem pool
         ideadstemst_to_ideadstemxf_phn      =>    cnveg_nitrogenflux_inst%ideadstemst_to_ideadstemxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem storage pool to dead stem transfer pool
         ideadstemxf_to_ideadstem_phn        =>    cnveg_nitrogenflux_inst%ideadstemxf_to_ideadstem_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem transfer pool to dead stem pool
         ilivecrootst_to_ilivecrootxf_phn    =>    cnveg_nitrogenflux_inst%ilivecrootst_to_ilivecrootxf_ph     , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root storage pool to live coarse root transfer pool
         ilivecrootxf_to_ilivecroot_phn      =>    cnveg_nitrogenflux_inst%ilivecrootxf_to_ilivecroot_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root transfer pool to live coarse root pool
         ideadcrootst_to_ideadcrootxf_phn    =>    cnveg_nitrogenflux_inst%ideadcrootst_to_ideadcrootxf_ph     , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root storage pool to dead coarse root transfer pool
         ideadcrootxf_to_ideadcroot_phn      =>    cnveg_nitrogenflux_inst%ideadcrootxf_to_ideadcroot_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root transfer pool to dead coarse root pool
         ilivestem_to_ideadstem_phn          =>    cnveg_nitrogenflux_inst%ilivestem_to_ideadstem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to dead stem pool
         ilivecroot_to_ideadcroot_phn        =>    cnveg_nitrogenflux_inst%ilivecroot_to_ideadcroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root to dead coarse root pool
         ileaf_to_iout_phn                   =>    cnveg_nitrogenflux_inst%ileaf_to_iout_ph                    , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to outside of vegetation pools
         ifroot_to_iout_phn                  =>    cnveg_nitrogenflux_inst%ifroot_to_iout_ph                   , & ! Input: [integer (:)] Index of phenology related C transfer from fine root pool to outside of vegetation pools
         ilivestem_to_iout_phn               =>    cnveg_nitrogenflux_inst%ilivestem_to_iout_ph                , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to outside of vegetation pools
         ilivestem_to_iout_gmc               =>    cnveg_carbonflux_inst%ilivestem_to_iout_gm                  , & ! Input: [integer (:)] Index of gap mortality related C transfer from live stem pool to outside of vegetation pools
         ilivestem_to_iout_gmn               =>    cnveg_nitrogenflux_inst%ilivestem_to_iout_gm                , & ! Input: [integer (:)] Index of gap mortality related N transfer from live stem pool to outside of vegetation pools
         ileaf_to_iretransn_phn              =>    cnveg_nitrogenflux_inst%ileaf_to_iretransn_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to retranslocation pool
         ilivestem_to_iretransn_phn          =>    cnveg_nitrogenflux_inst%ilivestem_to_iretransn_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to retranslocation pool
         ilivecroot_to_iretransn_phn         =>    cnveg_nitrogenflux_inst%ilivecroot_to_iretransn_ph          , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root pool to retranslocation pool
         igrain_to_iout_phn                  =>    cnveg_nitrogenflux_inst%igrain_to_iout_ph                     & ! Input: [integer (:)] Index of phenology related C transfer from grain pool to outside of vegetation pools
         )

      ! The litterfall transfer rate starts at 0.0 and increases linearly
      ! over time, with displayed growth going to 0.0 on the last day of litterfall

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! only calculate fluxes during offset period
         if (offset_flag(p) == 1._r8) then

            if (abs(offset_counter(p) - dt) <= dt/2._r8) then
               t1 = 1.0_r8 / dt
               frootc_to_litter(p) = t1 * frootc(p) + cpool_to_frootc(p)

               ! frootc_to_litter for matrix
               if (use_matrixcn) then
                  if(frootc(p) .gt. 0)then
                     frootc_to_litter(p) = frootc(p) * matrix_update_phc(p,ifroot_to_iout_phc,frootc_to_litter(p) / frootc(p),dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                  else
                     frootc_to_litter(p) = 0
                  end if
               else
                  ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
                  !                                        and CNNStateUpdate1::NStateUpdate1
               end if ! use_matrixcn
               ! this assumes that offset_counter == dt for crops
               ! if this were ever changed, we'd need to add code to the "else"
               if (ivt(p) >= npcropmin) then

                  ! How many harvests have occurred?
                  h = crop_inst%harvest_count(p)

                  ! Replenish the seed deficits from grain, if there is enough available
                  ! grain. (If there is not enough available grain, the seed deficits will
                  ! accumulate until there is eventually enough grain to replenish them.)
                  ! Note that, if there are multiple grain pools, we arbitrarily pull
                  ! first from grain pool 1, then from grain pool 2, etc., until we have
                  ! fully replenished the seed deficit.
                  if (for_testing_no_crop_seed_replenishment) then
                     cropseedc_deficit_remaining = 0._r8
                     cropseedn_deficit_remaining = 0._r8
                  else
                     cropseedc_deficit_remaining = -cropseedc_deficit(p)
                     cropseedn_deficit_remaining = -cropseedn_deficit(p)
                  end if
                  do k = repr_grain_min, repr_grain_max
                     cropseedc_deficit_to_restore = min(cropseedc_deficit_remaining, reproductivec(p,k))
                     cropseedc_deficit_remaining = cropseedc_deficit_remaining - cropseedc_deficit_to_restore
                     repr_grainc_to_seed(p,k) = t1 * cropseedc_deficit_to_restore
                     if (cropseedc_deficit_to_restore > 0._r8) then
                         repr_grainc_to_seed_perharv(p,h,k) = cropseedc_deficit_to_restore
                         repr_grainc_to_seed_thisyr(p,k) = repr_grainc_to_seed_thisyr(p,k) &
                             + repr_grainc_to_seed_perharv(p,h,k)
                     end if
                     cropseedn_deficit_to_restore = min(cropseedn_deficit_remaining, reproductiven(p,k))
                     cropseedn_deficit_remaining = cropseedn_deficit_remaining - cropseedn_deficit_to_restore
                     repr_grainn_to_seed(p,k) = t1 * cropseedn_deficit_to_restore
                     if (cropseedn_deficit_to_restore > 0._r8) then
                         repr_grainn_to_seed_perharv(p,h,k) = cropseedn_deficit_to_restore
                         repr_grainn_to_seed_thisyr(p,k) = repr_grainn_to_seed_thisyr(p,k) &
                             + repr_grainn_to_seed_perharv(p,h,k)
                     end if

                     ! Send the remaining grain to the food product pool
                     repr_grainc_to_food_thispool = cpool_to_reproductivec(p,k) - repr_grainc_to_seed(p,k)
                     repr_grainc_to_food(p,k) = t1 * reproductivec(p,k) &
                          + repr_grainc_to_food_thispool
                     if (reproductivec(p,k) + repr_grainc_to_food_thispool * dt > 0._r8) then
                         repr_grainc_to_food_perharv(p,h,k) = reproductivec(p,k) &
                             + repr_grainc_to_food_thispool * dt
                         repr_grainc_to_food_thisyr(p,k) = repr_grainc_to_food_thisyr(p,k) &
                             + repr_grainc_to_food_perharv(p,h,k)
                     end if
                     repr_grainn_to_food_thispool = npool_to_reproductiven(p,k) - repr_grainn_to_seed(p,k)
                     repr_grainn_to_food(p,k) = t1 * reproductiven(p,k) &
                          + npool_to_reproductiven(p,k) - repr_grainn_to_seed(p,k)
                     if (reproductiven(p,k) + repr_grainn_to_food_thispool * dt > 0._r8) then
                         repr_grainn_to_food_perharv(p,h,k) = reproductiven(p,k) &
                             + repr_grainn_to_food_thispool * dt
                         repr_grainn_to_food_thisyr(p,k) = repr_grainn_to_food_thisyr(p,k) &
                             + repr_grainn_to_food_perharv(p,h,k)
                     end if
                  end do

                  do k = repr_structure_min, repr_structure_max
                     repr_structurec_to_cropprod(p,k) = (t1 * reproductivec(p,k) + cpool_to_reproductivec(p,k)) &
                          * repr_structure_harvfrac(ivt(p), k)
                     repr_structurec_to_litter(p,k)   = (t1 * reproductivec(p,k) + cpool_to_reproductivec(p,k)) &
                          * (1._r8 - repr_structure_harvfrac(ivt(p), k))
                     repr_structuren_to_cropprod(p,k) = (t1 * reproductiven(p,k) + npool_to_reproductiven(p,k)) &
                          * repr_structure_harvfrac(ivt(p), k)
                     repr_structuren_to_litter(p,k)   = (t1 * reproductiven(p,k) + npool_to_reproductiven(p,k)) &
                          * (1._r8 - repr_structure_harvfrac(ivt(p), k))
                  end do
                  
                  ! Cut a certain fraction (i.e., biofuel_harvfrac(ivt(p))) (e.g., biofuel_harvfrac(ivt(p)=70% for bioenergy crops) of leaf C
                  ! and move this fration of leaf C to biofuel C, rather than move it to litter
                  leafc_to_biofuelc(p) = t1 * leafc(p) * biofuel_harvfrac(ivt(p))
                  leafc_remaining = leafc(p)*(1._r8-biofuel_harvfrac(ivt(p)))
                  leafn_to_biofueln(p) = t1 * leafn(p) * biofuel_harvfrac(ivt(p))
                  leafn_remaining = leafn(p)*(1._r8-biofuel_harvfrac(ivt(p)))

                  ! Cut a certain fraction (i.e., biofuel_harvfrac(ivt(p))) (e.g., biofuel_harvfrac(ivt(p)=70% for bioenergy crops) of livestem C
                  ! and move this fration of leaf C to biofuel C, rather than move it to litter
                  livestemc_to_biofuelc(p) = t1 * livestemc(p) * biofuel_harvfrac(ivt(p))
                  livestemn_to_biofueln(p) = t1 * livestemn(p) * biofuel_harvfrac(ivt(p))
                  livestemc_remaining = livestemc(p)*(1._r8-biofuel_harvfrac(ivt(p)))
                  livestemn_remaining = livestemn(p)*(1._r8-biofuel_harvfrac(ivt(p)))

                  ! Remove residues
                  leafc_to_removedresiduec(p) = t1 * leafc_remaining * crop_residue_removal_frac
                  leafn_to_removedresiduen(p) = t1 * leafn_remaining * crop_residue_removal_frac
                  livestemc_to_removedresiduec(p) = t1 * livestemc_remaining * crop_residue_removal_frac
                  livestemn_to_removedresiduen(p) = t1 * livestemn_remaining * crop_residue_removal_frac
                  leafc_remaining     = leafc_remaining     * (1._r8 - crop_residue_removal_frac)
                  leafn_remaining     = leafn_remaining     * (1._r8 - crop_residue_removal_frac)
                  livestemc_remaining = livestemc_remaining * (1._r8 - crop_residue_removal_frac)
                  livestemn_remaining = livestemn_remaining * (1._r8 - crop_residue_removal_frac)
                  
                  leafc_to_litter(p)  = t1 * leafc_remaining  + cpool_to_leafc(p)
                  livestemc_to_litter(p)   = t1 * livestemc_remaining  + cpool_to_livestemc(p)
                  livestemn_to_litter(p)   = t1 * livestemn_remaining
                  ! Sam Rabin 2023-09-11:
                  ! leafn_to_litter is calculated below based on leafc_to_litter (updated above)
                  ! as well as leaf C:N ratio (unaffected by biofuel harvest). It thus does not
                  ! need to be updated here.

                  ! Matrix for grain
                  ! Matrix for livestem/leaf to litter, biofuel, and removed residue
                  if(use_matrixcn)then
                     if(reproductivec(p,1) > 0._r8)then
                        grainc_to_out = reproductivec(p,1) * matrix_update_phc(p,igrain_to_iout_phc,(repr_grainc_to_seed(p,1) + repr_grainc_to_food(p,1)) / reproductivec(p,1),dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                     else
                        repr_grainc_to_seed(p,1) = 0._r8
                        repr_grainc_to_food(p,1) = 0._r8
                     end if
                     if(reproductiven(p,1) > 0._r8)then
                        grainn_to_out = reproductiven(p,1) * matrix_update_phn(p,igrain_to_iout_phn,(repr_grainn_to_seed(p,1) + repr_grainn_to_food(p,1)) / reproductiven(p,1),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                     else
                        repr_grainn_to_seed(p,1) = 0._r8
                        repr_grainn_to_food(p,1) = 0._r8
                     end if
                     if(livestemc(p) > 0._r8)then
                        livestemc_to_litter(p) = livestemc(p) * matrix_update_phc(p,ilivestem_to_iout_phc,livestemc_to_litter(p) / livestemc(p),dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                        livestemc_to_removedresiduec(p) = livestemc(p) * matrix_update_gmc(p,ilivestem_to_iout_gmc,livestemc_to_removedresiduec(p) / livestemc(p),dt,cnveg_carbonflux_inst,matrixcheck_ph,.True.)
                        livestemc_to_biofuelc(p) = livestemc(p) * matrix_update_gmc(p,ilivestem_to_iout_gmc,livestemc_to_biofuelc(p) / livestemc(p),dt,cnveg_carbonflux_inst,matrixcheck_ph,.True.)
                     else
                        livestemc_to_litter(p) = 0._r8
                        livestemc_to_removedresiduec(p) = 0._r8
                        livestemc_to_biofuelc(p) = 0._r8
                     end if
                     if(livestemn(p) > 0._r8)then
                        livestemn_to_biofueln(p) = livestemn(p) * matrix_update_gmn(p,ilivestem_to_iout_gmn,livestemn_to_biofueln(p) / livestemn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,.True.)
                        livestemn_to_removedresiduen(p) = livestemn(p) * matrix_update_gmn(p,ilivestem_to_iout_gmn,livestemn_to_removedresiduen(p) / livestemn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,.True.)
                        livestemn_to_litter(p) = livestemn(p) * matrix_update_phn(p,ilivestem_to_iout_phn, (1._r8- biofuel_harvfrac(ivt(p)))/dt, dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                     else
                        livestemn_to_biofueln(p) = 0._r8
                        livestemn_to_removedresiduen(p) = 0._r8
                        livestemn_to_litter(p) = 0._r8
                     end if
                     if(leafn(p) > 0._r8)then
                        leafn_to_biofueln(p) =  leafn(p) * matrix_update_gmn(p,ileaf_to_iout_gmn,leafn_to_biofueln(p) / leafn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,.True.)
                        leafn_to_removedresiduen(p) = leafn(p) * matrix_update_gmn(p,ileaf_to_iout_gmn,leafn_to_removedresiduen(p) / leafn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,.True.)
                     else
                        leafn_to_biofueln(p) = 0._r8
                        leafn_to_removedresiduen(p) = 0._r8
                     end if
                     if (leafc(p) > 0._r8)then
                        leafc_to_biofuelc(p) =  leafc(p) * matrix_update_gmc(p,ileaf_to_iout_gmc,leafc_to_biofuelc(p) / leafc(p),dt,cnveg_carbonflux_inst,matrixcheck_ph,.True.)
                        leafc_to_removedresiduec(p) = leafc(p) * matrix_update_gmc(p,ileaf_to_iout_gmc,leafc_to_removedresiduec(p) / leafc(p),dt,cnveg_carbonflux_inst,matrixcheck_ph,.True.)
                        leafc_to_litter(p) = leafc(p) * matrix_update_phc(p,ileaf_to_iout_phc,leafc_to_litter(p) / leafc(p),dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                     else
                        leafc_to_biofuelc(p) = 0._r8
                        leafc_to_removedresiduec(p) = 0._r8
                        leafc_to_litter(p) = 0._r8
                     end if
                  else
                     ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
                     !                                        and CNNStateUpdate1::NStateUpdate1
                  end if ! use_matrixcn
               end if
            else
               t1 = dt * 2.0_r8 / (offset_counter(p) * offset_counter(p))
               leafc_to_litter(p)  = prev_leafc_to_litter(p)  + t1*(leafc(p)  - prev_leafc_to_litter(p)*offset_counter(p))
               frootc_to_litter(p) = prev_frootc_to_litter(p) + t1*(frootc(p) - prev_frootc_to_litter(p)*offset_counter(p))

               if (use_matrixcn) then
                  if(leafc(p) .gt. 0)then
                     leafc_to_litter(p) = leafc(p) * matrix_update_phc(p,ileaf_to_iout_phc,leafc_to_litter(p) / leafc(p),dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                  else
                     leafc_to_litter(p) = 0
                  end if
                  if(frootc(p) .gt. 0)then
                     frootc_to_litter(p) = frootc(p) * matrix_update_phc(p,ifroot_to_iout_phc,frootc_to_litter(p) / frootc(p),dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                  else
                     frootc_to_litter(p) = 0  ! TODO slevis here and elsewhere
                  end if
               else
                  ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
                  !                                        and CNNStateUpdate1::NStateUpdate1
               end if !use_matrixcn
            end if
            
            if ( use_fun ) then
               if(leafc_to_litter(p)*dt.gt.leafc(p))then
                   leafc_to_litter(p) = leafc(p)/dt + cpool_to_leafc(p)
                  if (use_matrixcn) then
                     if(leafc(p) .gt. 0)then
                        leafc_to_litter(p) = leafc(p) * matrix_update_phc(p,ileaf_to_iout_phc,leafc_to_litter(p) / leafc(p),dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                     else
                        leafc_to_litter(p) = 0
                     end if
                  else
                     ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
                  end if !use_matrixcn 
               endif
               if(frootc_to_litter(p)*dt.gt.frootc(p))then
                   frootc_to_litter(p) = frootc(p)/dt + cpool_to_frootc(p)
                  if (use_matrixcn) then
                     if(frootc(p) .gt. 0)then
                        frootc_to_litter(p) = frootc(p) * matrix_update_phc(p,ifroot_to_iout_phc,frootc_to_litter(p) / frootc(p),dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
                     else
                        frootc_to_litter(p) = 0
                     end if
                  else
                     ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
                  end if !use_matrixcn
               endif
            end if
            
            
            if ( use_fun ) then
               leafc_to_litter_fun(p)      =  leafc_to_litter(p)
               leafn_to_retransn(p)        =  paid_retransn_to_npool(p) + free_retransn_to_npool(p)
               if (leafn(p).gt.0._r8) then
                  if (leafn(p)-leafn_to_retransn(p)*dt.gt.0._r8) then
                      leafcn_offset(p)     =  leafc(p)/(leafn(p)-leafn_to_retransn(p)*dt)
                  else
                      leafcn_offset(p)     =  leafc(p)/leafn(p)
                  end if
               else
                  leafcn_offset(p)         =  leafcn(ivt(p))
               end if
               leafn_to_litter(p)          =  leafc_to_litter(p)/leafcn_offset(p) - leafn_to_retransn(p)
               leafn_to_litter(p)          =  max(leafn_to_litter(p),0._r8)
               if (use_matrixcn) then   
                  if(leafn(p) .gt. 0)then
                      leafn_to_litter(p)   = leafn(p) * matrix_update_phn(p,ileaf_to_iout_phn,leafn_to_litter(p) / leafn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                      leafn_to_retransn(p) = leafn(p) * matrix_update_phn(p,ileaf_to_iretransn_phn,leafn_to_retransn(p) / leafn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                  else
                      leafn_to_litter(p)   = 0
                      leafn_to_retransn(p) = 0
                  end if
               else
                  ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
                  !                                        and CNNStateUpdate1::NStateUpdate1
               end if !use_matrixcn
               
               denom = ( leafn_to_retransn(p) + leafn_to_litter(p) )
               if ( denom /= 0.0_r8 ) then
                  fr_leafn_to_litter =  leafn_to_litter(p) / ( leafn_to_retransn(p) + leafn_to_litter(p) )
               else if ( leafn_to_litter(p) == 0.0_r8 ) then
                  fr_leafn_to_litter =  0.0_r8
               else
                  fr_leafn_to_litter =  1.0_r8
               end if

            else
               if (CNratio_floating .eqv. .true.) then    
                  fr_leafn_to_litter = 0.5_r8    ! assuming 50% of nitrogen turnover goes to litter
               end if
               ! calculate the leaf N litterfall and retranslocation
               leafn_to_litter(p)   = leafc_to_litter(p)  / lflitcn(ivt(p))
               leafn_to_retransn(p) = (leafc_to_litter(p) / leafcn(ivt(p))) - leafn_to_litter(p)

               if (use_matrixcn) then   
                  if(leafn(p) .gt. 0)then
                      leafn_to_litter(p)   = leafn(p) * matrix_update_phn(p,ileaf_to_iout_phn,leafn_to_litter(p) / leafn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                      leafn_to_retransn(p) = leafn(p) * matrix_update_phn(p,ileaf_to_iretransn_phn,leafn_to_retransn(p) / leafn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                  else
                      leafn_to_litter(p)   = 0
                      leafn_to_retransn(p) = 0
                  end if
               else
                  ! NOTE: The non matrix version of this is in CNNStateUpdate1::NStateUpdate1 EBK (11/26/2019)
               end if !use_matrixcn
            end if    

            ! calculate fine root N litterfall (no retranslocation of fine root N)
            frootn_to_litter(p) = frootc_to_litter(p) / frootcn(ivt(p))
            if (use_matrixcn) then   
               if(frootn(p) .gt. 0)then
                  frootn_to_litter(p) = frootn(p) * matrix_update_phn(p,ifroot_to_iout_phn,frootn_to_litter(p) / frootn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
               else
                  frootn_to_litter(p) = 0
               end if
            else
               ! NOTE: The non matrix version of this is in CNNStateUpdate1::NStateUpdate1 EBK (11/26/2019)
            end if !use_matrixcn
            
            if (CNratio_floating .eqv. .true.) then    
               if (leafc(p) == 0.0_r8) then    
                  ntovr_leaf = 0.0_r8    
               else    
                  ntovr_leaf = leafc_to_litter(p) * (leafn(p) / leafc(p))   
               end if   
           
               leafn_to_litter(p)   = fr_leafn_to_litter * ntovr_leaf
               leafn_to_retransn(p) = ntovr_leaf - leafn_to_litter(p)
               if (use_matrixcn) then   
                  if(leafn(p) .gt. 0)then
                      leafn_to_litter(p) = leafn(p) * matrix_update_phn(p,ileaf_to_iout_phn,leafn_to_litter(p) / leafn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                      leafn_to_retransn(p) = leafn(p) * matrix_update_phn(p,ileaf_to_iretransn_phn,leafn_to_retransn(p) / leafn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                  else
                      leafn_to_litter(p)   = 0
                      leafn_to_retransn(p) = 0
                  end if
               end if !use_matrixcn
               if (frootc(p) == 0.0_r8) then    
                   frootn_to_litter(p) = 0.0_r8    
                else    
                   frootn_to_litter(p) = frootc_to_litter(p) * (frootn(p) / frootc(p))   
                end if   
               if (use_matrixcn) then   
                  if(frootn(p) .gt. 0)then
                     frootn_to_litter(p) = frootn(p) * matrix_update_phn(p,ifroot_to_iout_phn,frootn_to_litter(p) / frootn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                  else
                     frootn_to_litter(p) = 0
                  end if
               else
                  ! NOTE: The non matrix version of this is in CNNStateUpdate1::NStateUpdate1 EBK (11/26/2019)
               end if !use_matrixcn
            end if  
            
            if ( use_fun ) then
               if(frootn_to_litter(p)*dt.gt.frootn(p))then
                  if (.not. use_matrixcn) then   
                     frootn_to_litter(p) = frootn(p)/dt
                  else
                     frootn_to_litter(p) = frootn(p) * matrix_update_phn(p,ifroot_to_iout_phn,1._r8/dt,dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                  end if
               endif    
            end if

            ! save the current litterfall fluxes
            prev_leafc_to_litter(p)  = leafc_to_litter(p)
            prev_frootc_to_litter(p) = frootc_to_litter(p)
                
         end if ! end if offset period

      end do ! end patch loop
      !matrix for leafn_to_retran will be added in allocation subroutine

    end associate 

  end subroutine CNOffsetLitterfall

  !-----------------------------------------------------------------------
  subroutine CNBackgroundLitterfall (num_soilp, filter_soilp, &
       cnveg_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! Determines the flux of C and N from displayed pools to litter
    ! pools as the result of background litter fall.
    !
    ! !USES:
    use CNSharedParamsMod   , only : use_fun
    use clm_varctl          , only : CNratio_floating    
    ! !ARGUMENTS:
    implicit none
    integer                       , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                       , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_state_type)        , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)  , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type), intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)   , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type) , intent(inout) :: cnveg_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: p            ! indices
    integer :: fp           ! filter patch index
    real(r8) :: fr_leafn_to_litter ! fraction of the nitrogen turnover that goes to litter; remaining fraction is retranslocated
    real(r8) :: ntovr_leaf  
    real(r8) :: denom       
    !-----------------------------------------------------------------------

    associate(                                                                     & 
         ivt               =>    patch%itype                                     , & ! Input:  [integer  (:) ]  patch vegetation type                                

         leafcn            =>    pftcon%leafcn                                   , & ! Input:  leaf C:N (gC/gN)                                  
         lflitcn           =>    pftcon%lflitcn                                  , & ! Input:  leaf litter C:N (gC/gN)                           
         frootcn           =>    pftcon%frootcn                                  , & ! Input:  fine root C:N (gC/gN)                             

         bglfr             =>    cnveg_state_inst%bglfr_patch                    , & ! Input:  [real(r8) (:) ]  background litterfall rate (1/s)                  

         leafc             =>    cnveg_carbonstate_inst%leafc_patch              , & ! Input:  [real(r8) (:) ]  (gC/m2) leaf C                                    
         frootc            =>    cnveg_carbonstate_inst%frootc_patch             , & ! Input:  [real(r8) (:) ]  (gC/m2) fine root C                               
         
         leafc_to_litter   =>    cnveg_carbonflux_inst%leafc_to_litter_patch     , & ! Output: [real(r8) (:) ]                                                    
         frootc_to_litter  =>    cnveg_carbonflux_inst%frootc_to_litter_patch    , & ! Output: [real(r8) (:) ]                                                    

         leafn             =>    cnveg_nitrogenstate_inst%leafn_patch            , & ! Input:  [real(r8) (:) ]  (gN/m2) leaf N  
         frootn            =>    cnveg_nitrogenstate_inst%frootn_patch           , & ! Input:  [real(r8) (:) ]  (gN/m2) fine root N 
         leafn_to_litter   =>    cnveg_nitrogenflux_inst%leafn_to_litter_patch   , & ! Output: [real(r8) (:) ]                                                    
         leafn_to_retransn =>    cnveg_nitrogenflux_inst%leafn_to_retransn_patch , & ! Output: [real(r8) (:) ]                                                    
         frootn_to_litter  =>    cnveg_nitrogenflux_inst%frootn_to_litter_patch  , & ! Output: [real(r8) (:) ]                                                    
         leafc_to_litter_fun   => cnveg_carbonflux_inst%leafc_to_litter_fun_patch, & ! Output:  [real(r8) (:) ] leaf C litterfall used by FUN (gC/m2/s)
         leafcn_offset         => cnveg_state_inst%leafcn_offset_patch           , & ! Output:  [real(r8) (:) ] Leaf C:N used by FUN
         free_retransn_to_npool=>    cnveg_nitrogenflux_inst%free_retransn_to_npool_patch  , & ! Input: [real(r8) (:) ] free leaf N to retranslocated N pool (gN/m2/s)          
         paid_retransn_to_npool=>    cnveg_nitrogenflux_inst%retransn_to_npool_patch   , & ! Input: [real(r8) (:) ] free leaf N to retranslocated N pool (gN/m2/s)

         ileafst_to_ileafxf_phc              =>    cnveg_carbonflux_inst%ileafst_to_ileafxf_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from leaf storage pool to leaf transfer pool
         ileafxf_to_ileaf_phc                =>    cnveg_carbonflux_inst%ileafxf_to_ileaf_ph                   , & ! Input: [integer (:)] Index of phenology related C transfer from leaf transfer pool to leaf pool
         ifrootst_to_ifrootxf_phc            =>    cnveg_carbonflux_inst%ifrootst_to_ifrootxf_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from fine root storage pool to fine root transfer pool
         ifrootxf_to_ifroot_phc              =>    cnveg_carbonflux_inst%ifrootxf_to_ifroot_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from fine root transfer pool to fine root pool
         ilivestemst_to_ilivestemxf_phc      =>    cnveg_carbonflux_inst%ilivestemst_to_ilivestemxf_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live stem storage pool to live stem transfer pool
         ilivestemxf_to_ilivestem_phc        =>    cnveg_carbonflux_inst%ilivestemxf_to_ilivestem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem transfer pool to live stem pool
         ideadstemst_to_ideadstemxf_phc      =>    cnveg_carbonflux_inst%ideadstemst_to_ideadstemxf_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem storage pool to dead stem transfer pool
         ideadstemxf_to_ideadstem_phc        =>    cnveg_carbonflux_inst%ideadstemxf_to_ideadstem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem transfer pool to dead stem pool
         ilivecrootst_to_ilivecrootxf_phc    =>    cnveg_carbonflux_inst%ilivecrootst_to_ilivecrootxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root storage pool to live coarse root transfer pool
         ilivecrootxf_to_ilivecroot_phc      =>    cnveg_carbonflux_inst%ilivecrootxf_to_ilivecroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root transfer pool to live coarse root pool
         ideadcrootst_to_ideadcrootxf_phc    =>    cnveg_carbonflux_inst%ideadcrootst_to_ideadcrootxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root storage pool to dead coarse root transfer pool
         ideadcrootxf_to_ideadcroot_phc      =>    cnveg_carbonflux_inst%ideadcrootxf_to_ideadcroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root transfer pool to dead coarse root pool
         ilivestem_to_ideadstem_phc          =>    cnveg_carbonflux_inst%ilivestem_to_ideadstem_ph             , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to dead stem pool
         ilivecroot_to_ideadcroot_phc        =>    cnveg_carbonflux_inst%ilivecroot_to_ideadcroot_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root to dead coarse root pool
         ileaf_to_iout_phc                   =>    cnveg_carbonflux_inst%ileaf_to_iout_ph                      , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to outside of vegetation pools
         ifroot_to_iout_phc                  =>    cnveg_carbonflux_inst%ifroot_to_iout_ph                     , & ! Input: [integer (:)] Index of phenology related C transfer from fine root pool to outside of vegetation pools
         ilivestem_to_iout_phc               =>    cnveg_carbonflux_inst%ilivestem_to_iout_ph                  , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to outside of vegetation pools
         igrain_to_iout_phc                  =>    cnveg_carbonflux_inst%igrain_to_iout_ph                     , & ! Input: [integer (:)] Index of phenology related C transfer from grain pool to outside of vegetation pools
         ileafst_to_ileafxf_phn              =>    cnveg_nitrogenflux_inst%ileafst_to_ileafxf_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from leaf storage pool to leaf transfer pool
         ileafxf_to_ileaf_phn                =>    cnveg_nitrogenflux_inst%ileafxf_to_ileaf_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from leaf transfer pool to leaf pool
         ifrootst_to_ifrootxf_phn            =>    cnveg_nitrogenflux_inst%ifrootst_to_ifrootxf_ph             , & ! Input: [integer (:)] Index of phenology related C transfer from fine root storage pool to fine root transfer pool
         ifrootxf_to_ifroot_phn              =>    cnveg_nitrogenflux_inst%ifrootxf_to_ifroot_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from fine root transfer pool to fine root pool
         ilivestemst_to_ilivestemxf_phn      =>    cnveg_nitrogenflux_inst%ilivestemst_to_ilivestemxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live stem storage pool to live stem transfer pool
         ilivestemxf_to_ilivestem_phn        =>    cnveg_nitrogenflux_inst%ilivestemxf_to_ilivestem_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live stem transfer pool to live stem pool
         ideadstemst_to_ideadstemxf_phn      =>    cnveg_nitrogenflux_inst%ideadstemst_to_ideadstemxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem storage pool to dead stem transfer pool
         ideadstemxf_to_ideadstem_phn        =>    cnveg_nitrogenflux_inst%ideadstemxf_to_ideadstem_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem transfer pool to dead stem pool
         ilivecrootst_to_ilivecrootxf_phn    =>    cnveg_nitrogenflux_inst%ilivecrootst_to_ilivecrootxf_ph     , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root storage pool to live coarse root transfer pool
         ilivecrootxf_to_ilivecroot_phn      =>    cnveg_nitrogenflux_inst%ilivecrootxf_to_ilivecroot_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root transfer pool to live coarse root pool
         ideadcrootst_to_ideadcrootxf_phn    =>    cnveg_nitrogenflux_inst%ideadcrootst_to_ideadcrootxf_ph     , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root storage pool to dead coarse root transfer pool
         ideadcrootxf_to_ideadcroot_phn      =>    cnveg_nitrogenflux_inst%ideadcrootxf_to_ideadcroot_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root transfer pool to dead coarse root pool
         ilivestem_to_ideadstem_phn          =>    cnveg_nitrogenflux_inst%ilivestem_to_ideadstem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to dead stem pool
         ilivecroot_to_ideadcroot_phn        =>    cnveg_nitrogenflux_inst%ilivecroot_to_ideadcroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root to dead coarse root pool
         ileaf_to_iout_phn                   =>    cnveg_nitrogenflux_inst%ileaf_to_iout_ph                    , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to outside of vegetation pools
         ifroot_to_iout_phn                  =>    cnveg_nitrogenflux_inst%ifroot_to_iout_ph                   , & ! Input: [integer (:)] Index of phenology related C transfer from fine root pool to outside of vegetation pools
         ilivestem_to_iout_phn               =>    cnveg_nitrogenflux_inst%ilivestem_to_iout_ph                , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to outside of vegetation pools
         ileaf_to_iretransn_phn              =>    cnveg_nitrogenflux_inst%ileaf_to_iretransn_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to retranslocation pool
         ilivestem_to_iretransn_phn          =>    cnveg_nitrogenflux_inst%ilivestem_to_iretransn_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to retranslocation pool
         ilivecroot_to_iretransn_phn         =>    cnveg_nitrogenflux_inst%ilivecroot_to_iretransn_ph          , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root pool to retranslocation pool
         igrain_to_iout_phn                  =>    cnveg_nitrogenflux_inst%igrain_to_iout_ph                     & ! Input: [integer (:)] Index of phenology related C transfer from grain pool to outside of vegetation pools
         )

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! only calculate these fluxes if the background litterfall rate is non-zero
         if (bglfr(p) > 0._r8) then
            ! units for bglfr are already 1/s
            leafc_to_litter(p)  = bglfr(p) * leafc(p)
            frootc_to_litter(p) = bglfr(p) * frootc(p)
            if (use_matrixcn) then
               leafc_to_litter(p)  = leafc(p)  * matrix_update_phc(p,ileaf_to_iout_phc,bglfr(p),dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
               frootc_to_litter(p) = frootc(p) * matrix_update_phc(p,ifroot_to_iout_phc,bglfr(p),dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
            end if
            if ( use_fun ) then
               leafc_to_litter_fun(p)     = leafc_to_litter(p)
               leafn_to_retransn(p)       = paid_retransn_to_npool(p) + free_retransn_to_npool(p)
               if (leafn(p).gt.0._r8) then
                  if (leafn(p)-leafn_to_retransn(p)*dt.gt.0._r8) then
                     leafcn_offset(p)     = leafc(p)/(leafn(p)-leafn_to_retransn(p)*dt)
                  else
                     leafcn_offset(p)     = leafc(p)/leafn(p)
                  end if
               else
                  leafcn_offset(p)        = leafcn(ivt(p))
               end if
               leafn_to_litter(p)         = leafc_to_litter(p)/leafcn_offset(p) - leafn_to_retransn(p)
               leafn_to_litter(p)         = max(leafn_to_litter(p),0._r8)
               if(use_matrixcn)then
                  if(leafn(p) .ne. 0._r8)then
                     leafn_to_litter(p)   = leafn(p) * matrix_update_phn(p,ileaf_to_iout_phn,leafn_to_litter(p)        / leafn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                     leafn_to_retransn(p) = leafn(p) * matrix_update_phn(p,ileaf_to_iretransn_phn,leafn_to_retransn(p) / leafn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                  end if
               else
                  ! NOTE: The non matrix version of this is in CNNStateUpdate1::NStateUpdate1 EBK (11/26/2019)
               end if !use_matrixcn

               denom = ( leafn_to_retransn(p) + leafn_to_litter(p) )
               if ( denom /= 0.0_r8 ) then
                  fr_leafn_to_litter =  leafn_to_litter(p) / ( leafn_to_retransn(p) + leafn_to_litter(p) )
               else if ( leafn_to_litter(p) == 0.0_r8 ) then
                  fr_leafn_to_litter =  0.0_r8
               else
                  fr_leafn_to_litter =  1.0_r8
               end if


            else
               if (CNratio_floating .eqv. .true.) then    
                  fr_leafn_to_litter = 0.5_r8    ! assuming 50% of nitrogen turnover goes to litter
               end if
               ! calculate the leaf N litterfall and retranslocation
               leafn_to_litter(p)   = leafc_to_litter(p)  / lflitcn(ivt(p))
               leafn_to_retransn(p) = (leafc_to_litter(p) / leafcn(ivt(p))) - leafn_to_litter(p)

               if (use_matrixcn) then   
                  if(leafn(p) .ne. 0)then
                     leafn_to_litter(p)   = leafn(p) * matrix_update_phn(p,ileaf_to_iout_phn,leafn_to_litter(p)        / leafn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                     leafn_to_retransn(p) = leafn(p) * matrix_update_phn(p,ileaf_to_iretransn_phn,leafn_to_retransn(p) / leafn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                  end if
               else
                  ! NOTE: The non matrix version of this is in CNNStateUpdate1::NStateUpdate1 EBK (11/26/2019)
               end if !use_matrixcn
            end if    

            ! calculate fine root N litterfall (no retranslocation of fine root N)
            frootn_to_litter(p) = frootc_to_litter(p) / frootcn(ivt(p))
            
            if (CNratio_floating .eqv. .true.) then    
               if (leafc(p) == 0.0_r8) then    
                  ntovr_leaf = 0.0_r8    
               else    
                  ntovr_leaf = leafc_to_litter(p) * (leafn(p) / leafc(p))   
               end if   
           
               leafn_to_litter(p)   = fr_leafn_to_litter * ntovr_leaf
               leafn_to_retransn(p) = ntovr_leaf - leafn_to_litter(p)
               if (use_matrixcn) then   
                  if(leafn(p) .gt. 0)then
                     leafn_to_litter(p)   = leafn(p) * matrix_update_phn(p,ileaf_to_iout_phn,leafn_to_litter(p)        / leafn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                     leafn_to_retransn(p) = leafn(p) * matrix_update_phn(p,ileaf_to_iretransn_phn,leafn_to_retransn(p) / leafn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                  else
                     leafn_to_litter(p)   = 0
                     leafn_to_retransn(p) = 0
                  end if
               else
                  ! NOTE: The non matrix version of this is in CNNStateUpdate1::NStateUpdate1 EBK (11/26/2019)
               end if !use_matrixcn
               if (frootc(p) == 0.0_r8) then    
                   frootn_to_litter(p) = 0.0_r8    
                else    
                   frootn_to_litter(p) = frootc_to_litter(p) * (frootn(p) / frootc(p))   
                end if   
            end if    

            if ( use_fun ) then
               if(frootn_to_litter(p)*dt.gt.frootn(p))then
                    frootn_to_litter(p) = frootn(p)/dt
               endif
            end if

            if (use_matrixcn) then   
               if(frootn(p) .ne. 0)then
                  frootn_to_litter(p) = frootn(p) * matrix_update_phn(p,ifroot_to_iout_phn,frootn_to_litter(p) / frootn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
               end if
            else
               ! NOTE: The non matrix version of this is in CNNStateUpdate1::NStateUpdate1 EBK (11/26/2019)
            end if !use_matrixcn
         end if

      end do
      !matrix for retransn_to_leafn will be added in allocation subroutine
    end associate 

  end subroutine CNBackgroundLitterfall

  !-----------------------------------------------------------------------
  subroutine CNLivewoodTurnover (num_soilp, filter_soilp, &
       cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! Determines the flux of C and N from live wood to
    ! dead wood pools, for stem and coarse root.
    !
    use CNSharedParamsMod, only: use_fun
    use clm_varctl          , only : CNratio_floating    
    ! !ARGUMENTS:
    integer                        , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                        , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_carbonstate_type)   , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type) , intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: p            ! indices
    integer :: fp           ! filter patch index
    real(r8):: ctovr        ! temporary variable for carbon turnover
    real(r8):: ntovr        ! temporary variable for nitrogen turnover
    !-----------------------------------------------------------------------

    associate(                                                                                   & 
         ivt                      =>    patch%itype                                              , & ! Input:  [integer  (:) ]  patch vegetation type                                

         woody                    =>    pftcon%woody                                           , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)
         livewdcn                 =>    pftcon%livewdcn                                        , & ! Input:  live wood (phloem and ray parenchyma) C:N (gC/gN) 
         deadwdcn                 =>    pftcon%deadwdcn                                        , & ! Input:  dead wood (xylem and heartwood) C:N (gC/gN)       

         livestemc                =>    cnveg_carbonstate_inst%livestemc_patch                 , & ! Input:  [real(r8) (:) ]  (gC/m2) live stem C                               
         livecrootc               =>    cnveg_carbonstate_inst%livecrootc_patch                , & ! Input:  [real(r8) (:) ]  (gC/m2) live coarse root C                        

         livestemn                =>    cnveg_nitrogenstate_inst%livestemn_patch               , & ! Input:  [real(r8) (:) ]  (gN/m2) live stem N                               
         livecrootn               =>    cnveg_nitrogenstate_inst%livecrootn_patch              , & ! Input:  [real(r8) (:) ]  (gN/m2) live coarse root N                        
         retransn                 =>    cnveg_nitrogenstate_inst%retransn_patch                , & ! Output: [real(r8)  (:)]
         
         livestemc_to_deadstemc   =>    cnveg_carbonflux_inst%livestemc_to_deadstemc_patch     , & ! Output: [real(r8) (:) ]                                                    
         livecrootc_to_deadcrootc =>    cnveg_carbonflux_inst%livecrootc_to_deadcrootc_patch   , & ! Output: [real(r8) (:) ]                                                    

         livestemn_to_deadstemn   =>    cnveg_nitrogenflux_inst%livestemn_to_deadstemn_patch   , & ! Output: [real(r8) (:) ]                                                    
         livestemn_to_retransn    =>    cnveg_nitrogenflux_inst%livestemn_to_retransn_patch    , & ! Output: [real(r8) (:) ]                                                    
         livecrootn_to_deadcrootn =>    cnveg_nitrogenflux_inst%livecrootn_to_deadcrootn_patch , & ! Output: [real(r8) (:) ]                                                    
         livecrootn_to_retransn   =>    cnveg_nitrogenflux_inst%livecrootn_to_retransn_patch   , & ! Output: [real(r8) (:) ] 
         free_retransn_to_npool   =>    cnveg_nitrogenflux_inst%free_retransn_to_npool_patch   , & ! Input:  [real(r8) (:) ] free leaf N to retranslocated N pool (gN/m2/s)          
         ileafst_to_ileafxf_phc              =>    cnveg_carbonflux_inst%ileafst_to_ileafxf_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from leaf storage pool to leaf transfer pool
         ileafxf_to_ileaf_phc                =>    cnveg_carbonflux_inst%ileafxf_to_ileaf_ph                   , & ! Input: [integer (:)] Index of phenology related C transfer from leaf transfer pool to leaf pool
         ifrootst_to_ifrootxf_phc            =>    cnveg_carbonflux_inst%ifrootst_to_ifrootxf_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from fine root storage pool to fine root transfer pool
         ifrootxf_to_ifroot_phc              =>    cnveg_carbonflux_inst%ifrootxf_to_ifroot_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from fine root transfer pool to fine root pool
         ilivestemst_to_ilivestemxf_phc      =>    cnveg_carbonflux_inst%ilivestemst_to_ilivestemxf_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live stem storage pool to live stem transfer pool
         ilivestemxf_to_ilivestem_phc        =>    cnveg_carbonflux_inst%ilivestemxf_to_ilivestem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem transfer pool to live stem pool
         ideadstemst_to_ideadstemxf_phc      =>    cnveg_carbonflux_inst%ideadstemst_to_ideadstemxf_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem storage pool to dead stem transfer pool
         ideadstemxf_to_ideadstem_phc        =>    cnveg_carbonflux_inst%ideadstemxf_to_ideadstem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem transfer pool to dead stem pool
         ilivecrootst_to_ilivecrootxf_phc    =>    cnveg_carbonflux_inst%ilivecrootst_to_ilivecrootxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root storage pool to live coarse root transfer pool
         ilivecrootxf_to_ilivecroot_phc      =>    cnveg_carbonflux_inst%ilivecrootxf_to_ilivecroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root transfer pool to live coarse root pool
         ideadcrootst_to_ideadcrootxf_phc    =>    cnveg_carbonflux_inst%ideadcrootst_to_ideadcrootxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root storage pool to dead coarse root transfer pool
         ideadcrootxf_to_ideadcroot_phc      =>    cnveg_carbonflux_inst%ideadcrootxf_to_ideadcroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root transfer pool to dead coarse root pool
         ilivestem_to_ideadstem_phc          =>    cnveg_carbonflux_inst%ilivestem_to_ideadstem_ph             , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to dead stem pool
         ilivecroot_to_ideadcroot_phc        =>    cnveg_carbonflux_inst%ilivecroot_to_ideadcroot_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root to dead coarse root pool
         ileaf_to_iout_phc                   =>    cnveg_carbonflux_inst%ileaf_to_iout_ph                      , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to outside of vegetation pools
         ifroot_to_iout_phc                  =>    cnveg_carbonflux_inst%ifroot_to_iout_ph                     , & ! Input: [integer (:)] Index of phenology related C transfer from fine root pool to outside of vegetation pools
         ilivestem_to_iout_phc               =>    cnveg_carbonflux_inst%ilivestem_to_iout_ph                  , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to outside of vegetation pools
         igrain_to_iout_phc                  =>    cnveg_carbonflux_inst%igrain_to_iout_ph                     , & ! Input: [integer (:)] Index of phenology related C transfer from grain pool to outside of vegetation pools
         ileafst_to_ileafxf_phn              =>    cnveg_nitrogenflux_inst%ileafst_to_ileafxf_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from leaf storage pool to leaf transfer pool
         ileafxf_to_ileaf_phn                =>    cnveg_nitrogenflux_inst%ileafxf_to_ileaf_ph                 , & ! Input: [integer (:)] Index of phenology related C transfer from leaf transfer pool to leaf pool
         ifrootst_to_ifrootxf_phn            =>    cnveg_nitrogenflux_inst%ifrootst_to_ifrootxf_ph             , & ! Input: [integer (:)] Index of phenology related C transfer from fine root storage pool to fine root transfer pool
         ifrootxf_to_ifroot_phn              =>    cnveg_nitrogenflux_inst%ifrootxf_to_ifroot_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from fine root transfer pool to fine root pool
         ilivestemst_to_ilivestemxf_phn      =>    cnveg_nitrogenflux_inst%ilivestemst_to_ilivestemxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live stem storage pool to live stem transfer pool
         ilivestemxf_to_ilivestem_phn        =>    cnveg_nitrogenflux_inst%ilivestemxf_to_ilivestem_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live stem transfer pool to live stem pool
         ideadstemst_to_ideadstemxf_phn      =>    cnveg_nitrogenflux_inst%ideadstemst_to_ideadstemxf_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem storage pool to dead stem transfer pool
         ideadstemxf_to_ideadstem_phn        =>    cnveg_nitrogenflux_inst%ideadstemxf_to_ideadstem_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from dead stem transfer pool to dead stem pool
         ilivecrootst_to_ilivecrootxf_phn    =>    cnveg_nitrogenflux_inst%ilivecrootst_to_ilivecrootxf_ph     , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root storage pool to live coarse root transfer pool
         ilivecrootxf_to_ilivecroot_phn      =>    cnveg_nitrogenflux_inst%ilivecrootxf_to_ilivecroot_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root transfer pool to live coarse root pool
         ideadcrootst_to_ideadcrootxf_phn    =>    cnveg_nitrogenflux_inst%ideadcrootst_to_ideadcrootxf_ph     , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root storage pool to dead coarse root transfer pool
         ideadcrootxf_to_ideadcroot_phn      =>    cnveg_nitrogenflux_inst%ideadcrootxf_to_ideadcroot_ph       , & ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root transfer pool to dead coarse root pool
         ilivestem_to_ideadstem_phn          =>    cnveg_nitrogenflux_inst%ilivestem_to_ideadstem_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to dead stem pool
         ilivecroot_to_ideadcroot_phn        =>    cnveg_nitrogenflux_inst%ilivecroot_to_ideadcroot_ph         , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root to dead coarse root pool
         ileaf_to_iout_phn                   =>    cnveg_nitrogenflux_inst%ileaf_to_iout_ph                    , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to outside of vegetation pools
         ifroot_to_iout_phn                  =>    cnveg_nitrogenflux_inst%ifroot_to_iout_ph                   , & ! Input: [integer (:)] Index of phenology related C transfer from fine root pool to outside of vegetation pools
         ilivestem_to_iout_phn               =>    cnveg_nitrogenflux_inst%ilivestem_to_iout_ph                , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to outside of vegetation pools
         ileaf_to_iretransn_phn              =>    cnveg_nitrogenflux_inst%ileaf_to_iretransn_ph               , & ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to retranslocation pool
         ilivestem_to_iretransn_phn          =>    cnveg_nitrogenflux_inst%ilivestem_to_iretransn_ph           , & ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to retranslocation pool
         ilivecroot_to_iretransn_phn         =>    cnveg_nitrogenflux_inst%ilivecroot_to_iretransn_ph          , & ! Input: [integer (:)] Index of phenology related C transfer from live coarse root pool to retranslocation pool
         iretransn_to_iout                   =>    cnveg_nitrogenflux_inst%iretransn_to_iout_ph                , & ! Input: [integer    ]
         igrain_to_iout_phn                  =>    cnveg_nitrogenflux_inst%igrain_to_iout_ph                     & ! Input: [integer (:)] Index of phenology related C transfer from grain pool to outside of vegetation pools
         )



      ! patch loop
ptch: do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! only calculate these fluxes for woody types
         if (woody(ivt(p)) > 0._r8) then

            ! live stem to dead stem turnover

            ctovr = livestemc(p) * lwtop
            ntovr = ctovr / livewdcn(ivt(p))
            livestemc_to_deadstemc(p) = ctovr
            livestemn_to_deadstemn(p) = ctovr / deadwdcn(ivt(p))
            if( use_matrixcn)then
               livestemc_to_deadstemc(p) = livestemc(p) * matrix_update_phc(p,ilivestem_to_ideadstem_phc,lwtop,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
               if (livestemn(p) .gt. 0.0_r8) then
                  livestemn_to_deadstemn(p) = livestemn(p) * matrix_update_phn(p,ilivestem_to_ideadstem_phn,livestemn_to_deadstemn(p)/livestemn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
               else
                  livestemn_to_deadstemn(p) = 0
               end if
            else
               ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
               !                                        and CNNStateUpdate1::NStateUpdate1
            end if
            if (CNratio_floating .eqv. .true.) then    
               if (livestemc(p) == 0.0_r8) then    
                   ntovr = 0.0_r8    
                   livestemn_to_deadstemn(p) = 0.0_r8 
                else    
                   ntovr = ctovr * (livestemn(p) / livestemc(p))   
                   livestemn_to_deadstemn(p) = ctovr / deadwdcn(ivt(p)) 
                end if   

               if (use_matrixcn)then 
                  if (livestemn(p) .gt. 0.0_r8) then
                     livestemn_to_deadstemn(p) = livestemn(p) * matrix_update_phn(p,ilivestem_to_ideadstem_phn,&
                                                 livestemn_to_deadstemn(p) / livestemn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
                  else
                     livestemn_to_deadstemn(p) = 0
                  end if
               else
                  ! NOTE: The non matrix version of this is in CNNStateUpdate1::NStateUpdate1 EBK (11/26/2019)
               end if
            end if
            
            livestemn_to_retransn(p)  = ntovr - livestemn_to_deadstemn(p)
            !matrix for livestemn_to_retransn will be added in allocation subroutine

            ! live coarse root to dead coarse root turnover

            ctovr = livecrootc(p) * lwtop
            ntovr = ctovr / livewdcn(ivt(p))
            if(.not. use_matrixcn)then
               ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
               !                                        and CNNStateUpdate1::NStateUpdate1
               livecrootc_to_deadcrootc(p) = ctovr
               livecrootn_to_deadcrootn(p) = ctovr / deadwdcn(ivt(p))
            else
               livecrootc_to_deadcrootc(p) = livecrootc(p) * matrix_update_phc(p,ilivecroot_to_ideadcroot_phc,lwtop,dt,cnveg_carbonflux_inst,matrixcheck_ph,acc_ph)
               livecrootn_to_deadcrootn(p) = livecrootn(p) * matrix_update_phn(p,ilivecroot_to_ideadcroot_phn,lwtop/deadwdcn(ivt(p)),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph)
            end if !use_matrixcn
            
            if (CNratio_floating .eqv. .true.) then    
               if (livecrootc(p) == 0.0_r8) then    
                  ntovr = 0.0_r8    
                  livecrootn_to_deadcrootn(p) = 0.0_r8 
               else    
                  ntovr = ctovr * (livecrootn(p) / livecrootc(p))   
                  livecrootn_to_deadcrootn(p) = ctovr / deadwdcn(ivt(p)) 
               end if   

               if (use_matrixcn)then 
                  if (livecrootn(p) .ne.0.0_r8 )then
                     livecrootn_to_deadcrootn(p) = matrix_update_phn(p,ilivecroot_to_ideadcroot_phn,&
                                                   livecrootn_to_deadcrootn(p) / livecrootn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph) * livecrootn(p)
                  end if
               else
                  ! NOTE: The non matrix version of this is in CNNStateUpdate1::NStateUpdate1 EBK (11/26/2019)
               end if !use_matrixcn
            end if

            livecrootn_to_retransn(p)  = ntovr - livecrootn_to_deadcrootn(p)
               if(use_matrixcn)then
                  if(livecrootn(p) .gt. 0.0_r8) then
                     livecrootn_to_retransn(p) = matrix_update_phn(p,ilivecroot_to_iretransn_phn,&
                                                 livecrootn_to_retransn(p) / livecrootn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph) * livecrootn(p)
                  else
                     livecrootn_to_retransn(p) = 0
                  end if
                  if(livestemn(p) .gt. 0.0_r8) then
                     livestemn_to_retransn(p) = matrix_update_phn(p,ilivestem_to_iretransn_phn,&
                                                livestemn_to_retransn(p) / livestemn(p),dt,cnveg_nitrogenflux_inst,matrixcheck_ph,acc_ph) * livestemn(p)
                  else
                     livestemn_to_retransn(p)  = 0
                  end if
                  ! WW change logic so livestem_retrans goes to npool (via
                  ! free_retrans flux)
                  ! this should likely be done more cleanly if it works, i.e. not
                  ! update fluxes w/ states
                  ! additional considerations for crop?
                  ! The non-matrix version of this is in NStateUpdate1
                  if (use_fun) then
                     if (retransn(p) .gt. 0._r8) then
                        ! The acc matrix check MUST be turned on, or this will
                        ! fail with Nitrogen balance error EBK 03/11/2021
                        free_retransn_to_npool(p) = free_retransn_to_npool(p) + retransn(p) * matrix_update_phn(p,iretransn_to_iout, &
                                                    (livestemn_to_retransn(p) + livecrootn_to_retransn(p)) / retransn(p),dt,         &
                                                    cnveg_nitrogenflux_inst, matrixcheck_ph, acc=.true.)
                     else
                        free_retransn_to_npool(p) = 0._r8
                     end if
                  end if
               end if !use_matrixcn

         end if

      end do ptch

    end associate 

  end subroutine CNLivewoodTurnover

  !-----------------------------------------------------------------------
  subroutine CNCropHarvestToProductPools(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, &
       cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! If using prognostic crop, then move any necessary harvested amounts into fluxes
    ! destined for the product pools.
    !
    ! !USES:
    use clm_varctl    , only : use_crop
    use clm_varctl    , only : use_grainproduct
    use subgridAveMod , only : p2c
    !
    ! !ARGUMENTS:
    type(bounds_type)             , intent(in)    :: bounds
    integer                       , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                       , intent(in)    :: filter_soilp(:) ! filter for soil patches
    integer                       , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                       , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(cnveg_carbonflux_type)   , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type) , intent(inout) :: cnveg_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: fp, p, k

    character(len=*), parameter :: subname = 'CNCropHarvestToProductPools'
    !-----------------------------------------------------------------------

    if (use_crop) then
       do fp = 1, num_soilp
          p = filter_soilp(fp)
          cnveg_carbonflux_inst%crop_harvestc_to_cropprodc_patch(p) = &
               cnveg_carbonflux_inst%leafc_to_biofuelc_patch(p) + &
               cnveg_carbonflux_inst%livestemc_to_biofuelc_patch(p) + &
               cnveg_carbonflux_inst%leafc_to_removedresiduec_patch(p) + &
               cnveg_carbonflux_inst%livestemc_to_removedresiduec_patch(p)
          cnveg_nitrogenflux_inst%crop_harvestn_to_cropprodn_patch(p) = &
               cnveg_nitrogenflux_inst%leafn_to_biofueln_patch(p) + &
               cnveg_nitrogenflux_inst%livestemn_to_biofueln_patch(p) + &
               cnveg_nitrogenflux_inst%leafn_to_removedresiduen_patch(p) + &
               cnveg_nitrogenflux_inst%livestemn_to_removedresiduen_patch(p)
       end do

       if (use_grainproduct) then
          do k = repr_grain_min, repr_grain_max
             do fp = 1, num_soilp
                p = filter_soilp(fp)
                cnveg_carbonflux_inst%crop_harvestc_to_cropprodc_patch(p) = &
                     cnveg_carbonflux_inst%crop_harvestc_to_cropprodc_patch(p) + &
                     cnveg_carbonflux_inst%repr_grainc_to_food_patch(p,k)
                cnveg_nitrogenflux_inst%crop_harvestn_to_cropprodn_patch(p) = &
                     cnveg_nitrogenflux_inst%crop_harvestn_to_cropprodn_patch(p) + &
                     cnveg_nitrogenflux_inst%repr_grainn_to_food_patch(p,k)
             end do
          end do
       end if

       do k = repr_structure_min, repr_structure_max
          do fp = 1, num_soilp
             p = filter_soilp(fp)
             cnveg_carbonflux_inst%crop_harvestc_to_cropprodc_patch(p) = &
                  cnveg_carbonflux_inst%crop_harvestc_to_cropprodc_patch(p) + &
                  cnveg_carbonflux_inst%repr_structurec_to_cropprod_patch(p,k)
             cnveg_nitrogenflux_inst%crop_harvestn_to_cropprodn_patch(p) = &
                  cnveg_nitrogenflux_inst%crop_harvestn_to_cropprodn_patch(p) + &
                  cnveg_nitrogenflux_inst%repr_structuren_to_cropprod_patch(p,k)
          end do
       end do

       call p2c (bounds, num_soilc, filter_soilc, &
            cnveg_carbonflux_inst%crop_harvestc_to_cropprodc_patch(bounds%begp:bounds%endp), &
            cnveg_carbonflux_inst%crop_harvestc_to_cropprodc_col(bounds%begc:bounds%endc))

       call p2c (bounds, num_soilc, filter_soilc, &
            cnveg_nitrogenflux_inst%crop_harvestn_to_cropprodn_patch(bounds%begp:bounds%endp), &
            cnveg_nitrogenflux_inst%crop_harvestn_to_cropprodn_col(bounds%begc:bounds%endc))

    end if

  end subroutine CNCropHarvestToProductPools

  !-----------------------------------------------------------------------
  subroutine CNLitterToColumn (bounds, num_bgc_vegp, filter_bgc_vegp,         &
       cnveg_state_inst,cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
       leaf_prof_patch, froot_prof_patch)
    !
    ! !DESCRIPTION:
    ! called at the end of cn_phenology to gather all patch-level litterfall fluxes
    ! to the column level and assign them to the three litter pools
    !
    ! !USES:
    use clm_varpar , only : nlevdecomp
    use pftconMod  , only : npcropmin
    use clm_varctl , only : use_grainproduct
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_bgc_vegp       ! number of bgc veg patches
    integer                         , intent(in)    :: filter_bgc_vegp(:) ! filter for bgc veg patches
    type(cnveg_state_type)          , intent(in)    :: cnveg_state_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    real(r8)                        , intent(in)    :: leaf_prof_patch(bounds%begp:,1:)
    real(r8)                        , intent(in)    :: froot_prof_patch(bounds%begp:,1:)
    !
    ! !LOCAL VARIABLES:
    integer :: fp,c,p,k,j,i  ! indices
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(leaf_prof_patch)   == (/bounds%endp,nlevdecomp_full/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(froot_prof_patch)  == (/bounds%endp,nlevdecomp_full/)), sourcefile, __LINE__)

    associate(                                                                                & 
         leaf_prof                 => leaf_prof_patch                                       , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
         froot_prof                => froot_prof_patch                                      , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     

         ivt                       => patch%itype                                             , & ! Input:  [integer  (:)   ]  patch vegetation type                                
         wtcol                     => patch%wtcol                                             , & ! Input:  [real(r8) (:)   ]  weight (relative to column) for this patch (0-1)    

         lf_f                      => pftcon%lf_f                                           , & ! Input:  leaf litter fractions
         fr_f                      => pftcon%fr_f                                           , & ! Input:  fine root litter fractions

         leafc_to_litter           => cnveg_carbonflux_inst%leafc_to_litter_patch           , & ! Input:  [real(r8) (:)   ]  leaf C litterfall (gC/m2/s)                       
         frootc_to_litter          => cnveg_carbonflux_inst%frootc_to_litter_patch          , & ! Input:  [real(r8) (:)   ]  fine root N litterfall (gN/m2/s)                  
         livestemc_to_litter       => cnveg_carbonflux_inst%livestemc_to_litter_patch       , & ! Input:  [real(r8) (:)   ]  live stem C litterfall (gC/m2/s)                  
         repr_grainc_to_food       => cnveg_carbonflux_inst%repr_grainc_to_food_patch       , & ! Input:  [real(r8) (:,:) ]  grain C to food (gC/m2/s)
         repr_structurec_to_litter => cnveg_carbonflux_inst%repr_structurec_to_litter_patch,  & ! Input:  [real(r8) (:,:) ] reproductive structure C to litter (gC/m2/s)
         phenology_c_to_litr_c     => cnveg_carbonflux_inst%phenology_c_to_litr_c_col       , & ! Output: [real(r8) (:,:,:) ]  C fluxes associated with phenology (litterfall and crop) to litter pools (gC/m3/s)

         livestemn_to_litter       => cnveg_nitrogenflux_inst%livestemn_to_litter_patch     , & ! Input:  [real(r8) (:)   ]  livestem N to litter (gN/m2/s)                    
         repr_grainn_to_food       => cnveg_nitrogenflux_inst%repr_grainn_to_food_patch     , & ! Input:  [real(r8) (:,:) ]  grain N to food (gN/m2/s)
         repr_structuren_to_litter => cnveg_nitrogenflux_inst%repr_structuren_to_litter_patch,& ! Input:  [real(r8) (:,:) ] reproductive structure N to litter (gN/m2/s)
         leafn_to_litter           => cnveg_nitrogenflux_inst%leafn_to_litter_patch         , & ! Input:  [real(r8) (:)   ]  leaf N litterfall (gN/m2/s)                       
         frootn_to_litter          => cnveg_nitrogenflux_inst%frootn_to_litter_patch        , & ! Input:  [real(r8) (:)   ]  fine root N litterfall (gN/m2/s)                  
         phenology_n_to_litr_n     => cnveg_nitrogenflux_inst%phenology_n_to_litr_n_col       & ! Output: [real(r8) (:,:,:) ]  N fluxes associated with phenology (litterfall and crop) to litter pools (gN/m3/s)
         )
      
      do_nlev: do j = 1, nlevdecomp

         do_vegp: do fp = 1,num_bgc_vegp
            p = filter_bgc_vegp(fp)
            c = patch%column(p)

            do_ilit: do i = i_litr_min, i_litr_max
               ! leaf litter carbon fluxes
               phenology_c_to_litr_c(c,j,i) = &
                    phenology_c_to_litr_c(c,j,i) + &
                    leafc_to_litter(p) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j)
               
               ! leaf litter nitrogen fluxes
               phenology_n_to_litr_n(c,j,i) = &
                    phenology_n_to_litr_n(c,j,i) + &
                    leafn_to_litter(p) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j)
               
               ! fine root litter carbon fluxes
               phenology_c_to_litr_c(c,j,i) = &
                    phenology_c_to_litr_c(c,j,i) + &
                    frootc_to_litter(p) * fr_f(ivt(p),i) * wtcol(p) * froot_prof(p,j)
               
               ! fine root litter nitrogen fluxes
               phenology_n_to_litr_n(c,j,i) = &
                    phenology_n_to_litr_n(c,j,i) + &
                    frootn_to_litter(p) * fr_f(ivt(p),i) * wtcol(p) * froot_prof(p,j)
            end do do_ilit
            
            ! agroibis puts crop stem litter together with leaf litter
            ! so I've used the leaf lf_f* parameters instead of making
            ! new ones for now (slevis)
            ! also for simplicity I've put "food" into the litter pools
            
            if (ivt(p) >= npcropmin) then ! add livestemc to litter
               do i = i_litr_min, i_litr_max
                  ! stem litter carbon fluxes
                  phenology_c_to_litr_c(c,j,i) = &
                       phenology_c_to_litr_c(c,j,i) + &
                       livestemc_to_litter(p) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j)
                  
                  ! stem litter nitrogen fluxes
                  phenology_n_to_litr_n(c,j,i) = &
                       phenology_n_to_litr_n(c,j,i) + &
                       livestemn_to_litter(p) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j)
               end do
               
               if (.not. use_grainproduct) then
                  do i = i_litr_min, i_litr_max
                     do k = repr_grain_min, repr_grain_max
                        ! grain litter carbon fluxes
                        phenology_c_to_litr_c(c,j,i) = &
                             phenology_c_to_litr_c(c,j,i) + &
                             repr_grainc_to_food(p,k) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j)

                        ! grain litter nitrogen fluxes
                        phenology_n_to_litr_n(c,j,i) = &
                             phenology_n_to_litr_n(c,j,i) + &
                             repr_grainn_to_food(p,k) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j)
                     end do
                  end do
               end if

               do i = i_litr_min, i_litr_max
                  do k = repr_structure_min, repr_structure_max
                     ! reproductive structure litter carbon fluxes
                     phenology_c_to_litr_c(c,j,i) = &
                          phenology_c_to_litr_c(c,j,i) + &
                          repr_structurec_to_litter(p,k) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j)

                     ! reproductive structure litter nitrogen fluxes
                     phenology_n_to_litr_n(c,j,i) = &
                          phenology_n_to_litr_n(c,j,i) + &
                          repr_structuren_to_litter(p,k) * lf_f(ivt(p),i) * wtcol(p) * leaf_prof(p,j)
                  end do
               end do
            end if
         end do do_vegp
      end do do_nlev

    end associate

  end subroutine CNLitterToColumn

end module CNPhenologyMod
