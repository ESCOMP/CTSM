module DistParamType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Spatially distributed parameter data type allocation and initialization
  ! -------------------------------------------------------- 
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_sys_mod    , only : shr_sys_abort
  use clm_varcon     , only : spval, ispval, grlnd
  use clm_varctl     , only : iulog
  use clm_varctl     , only : paramfile
  use ColumnType     , only : col
  use spmdMod        , only : masterproc, mpicom
  use fileutils      , only : getfil, opnfil, getavu, relavu
  use shr_nl_mod     , only : shr_nl_find_group_name
  use shr_mpi_mod    , only : shr_mpi_bcast
  use decompMod      , only : bounds_type
  use ncdio_pio      , only : ncd_pio_closefile, ncd_pio_openfile
  use ncdio_pio      , only : file_desc_t, ncd_inqdid, ncd_inqdlen

  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  character(len=*), parameter, private :: sourcefile = __FILE__
  !

  type, public :: distparam_class
     logical               :: is_distributed = .false. ! is parameter spatially distributed?
     character(len=256)    :: name                     ! name on parameter files 
     real(r8), allocatable :: val(:)                   ! parameter value array (either distributed or length 1)
   contains
     procedure, public :: param_val                    ! return parameter value
  end type distparam_class
  
  type, public :: distributed_parameter_type
     ! all parameters that could *potentially* be spatially distributed

     ! SoilHydrologyMod
     class(distparam_class), pointer :: aq_sp_yield_min => NULL()                 ! Minimum aquifer specific yield (unitless)
     class(distparam_class), pointer :: n_baseflow => NULL()                      ! Drainage power law exponent (unitless)
     class(distparam_class), pointer :: perched_baseflow_scalar => NULL()         ! Scalar multiplier for perched base flow rate (kg/m2/s)
     class(distparam_class), pointer :: e_ice => NULL()                           ! Soil ice impedance factor (unitless)
     class(distparam_class), pointer :: baseflow_scalar => NULL()                 ! Scalar multiplier for base flow rate ()

     ! SaturatedExcessRunoff
     class(distparam_class), pointer :: fff => NULL()                             ! Decay factor for fractional saturated area (1/m)
     
     ! initVerticalMod
     class(distparam_class), pointer :: slopebeta => NULL()                       ! exponent for microtopography pdf sigma (unitless)
     class(distparam_class), pointer :: slopemax => NULL()                        ! max topographic slope for microtopography pdf sigma (unitless)

     ! SnowCoverFractionSwensonLawrence2012Mod
     class(distparam_class), pointer :: n_melt_coef => NULL()                     ! SCA shape parameter
     class(distparam_class), pointer :: accum_factor => NULL()                    ! Accumulation constant for fractional snow covered area (unitless)

     ! SnowHydrologyMod
     class(distparam_class), pointer :: upplim_destruct_metamorph => NULL()       ! Upper limit on destructive metamorphism compaction (kg/m3)

     ! SoilStateInitTimeConstMod
     class(distparam_class), pointer :: bsw_sf => NULL()                          ! Scale factor for bsw (unitless)
     class(distparam_class), pointer :: hksat_sf => NULL()                        ! Scale factor for hksat (unitless)
     class(distparam_class), pointer :: sucsat_sf => NULL()                       ! Scale factor for sucsat (unitless)
     class(distparam_class), pointer :: watsat_sf => NULL()                       ! Scale factor for watsat (unitless)

     ! SurfaceResistanceMod
     class(distparam_class), pointer :: d_max => NULL()                           ! Dry surface layer parameter (mm)
     class(distparam_class), pointer :: frac_sat_soil_dsl_init => NULL()          ! Fraction of saturated soil for moisture value at which DSL initiates (unitless)

     ! SurfaceWaterMod
     class(distparam_class), pointer :: pc => NULL()                              ! Threshold probability for surface water (unitless)
     class(distparam_class), pointer :: mu => NULL()                              ! Connectivity exponent for surface water (unitless)

     ! WaterDiagnosticBulkType
     class(distparam_class), pointer :: zlnd => NULL()                            ! Momentum roughness length for soil, glacier, wetland (m)
     class(distparam_class), pointer :: snw_rds_min => NULL()                     ! minimum allowed snow effective radius (also cold "fresh snow" value) [microns]

     ! atm2lndType
     class(distparam_class), pointer :: precip_repartition_nonglc_all_rain_t_celsius => NULL()      ! Rain temperature threshold for non-glacier landunits (C)
     class(distparam_class), pointer :: precip_repartition_nonglc_all_snow_t_celsius => NULL()      ! Snow temperature threshold for non-glacier landunits (C)
     class(distparam_class), pointer :: precip_repartition_glc_all_rain_t_celsius => NULL()      ! Rain temperature threshold for glacier landunits (C)
     class(distparam_class), pointer :: precip_repartition_glc_all_snow_t_celsius => NULL()      ! Snow temperature threshold for glacier landunits (C)

   contains

     procedure, public :: Init
     procedure, public :: Clean

  end type distributed_parameter_type

  type(distributed_parameter_type), public, target        :: distributed_parameters 

  public :: InitGlobalParameters

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  function param_val(this,g)
    !
    ! !DESCRIPTION:
    ! Return distributed parameter value if allocated, otherwise return scalar value
    !
    ! !ARGUMENTS:
    implicit none
    class(distparam_class) :: this
    real(r8)                          :: param_val
    integer                           :: g
    !
    if ( this%is_distributed )then
       param_val = this%val(g)
    else
       param_val = this%val(1)
    end if
  end function param_val

  !------------------------------------------------------------------------
  subroutine ReadScalarParameter(this,ncid)
    !
    ! !USES:
    use paramUtilMod    , only : readNcdioScalar
    ! !ARGUMENTS:
    implicit none
    class(distparam_class), intent(inout) :: this
    type(file_desc_t)     , intent(inout) :: ncid  ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    real(r8)           :: fscalar_in               ! read in - scalar - float
    character(len=*), parameter :: subname = 'ReadScalarParameter'
    !--------------------------------------------------------------------

    if ( .not. this%is_distributed ) then
       allocate(this%val(1))
       call readNcdioScalar(ncid, this%name, subname, fscalar_in)
       this%val(1) = fscalar_in
    endif

  end subroutine ReadScalarParameter
  
  !------------------------------------------------------------------------
  subroutine InitGlobalParameters(bounds)
    !
    ! !USES:
    use ncdio_pio       , only : check_var, ncd_io
    use paramUtilMod    , only : readNcdioScalar
    use abortutils      , only : endrun
    use shr_log_mod     , only : errMsg => shr_log_errMsg
    use domainMod       , only : ldomain
    use clm_varctl      , only : NLFilename => NLFilename_in

    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in)      :: bounds
    !
    ! !LOCAL VARIABLES:
    logical            :: readvar               ! whether the variable was found
    integer            :: c, g
    real(r8)           :: fscalar_in            ! read in - scalar - float
    real(r8), pointer  :: fparam_in(:)          ! read in - 1D - float
    type(file_desc_t)  :: ncid                  ! pio netCDF file id
    character(len=256) :: locfn                 ! local filename
    character(len=*), parameter :: subname = 'InitGlobalParameters'
    !--------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) trim(subname)//' :: reading CLM parameters'
    end if

    !distributed_parameter_stream%Init must be called before this routine is called.

    ! global (scalar) parameters
    call getfil (paramfile, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    ! Any parameters not read from stream will be read here from parameter file or namelist
    !-----------------------------------------------------------
    ! SoilHydrology !
    !-----------------------------------------------------------

    ! Minimum aquifer specific yield (unitless)
    call ReadScalarParameter(distributed_parameters%aq_sp_yield_min,ncid)

    ! Drainage power law exponent (unitless)
    call ReadScalarParameter(distributed_parameters%n_baseflow,ncid)
    
    ! Scalar multiplier for perched base flow rate ()
    call ReadScalarParameter(distributed_parameters%baseflow_scalar,ncid)
    
    ! Scalar multiplier for perched base flow rate (kg/m2/s)
    call ReadScalarParameter(distributed_parameters%perched_baseflow_scalar,ncid)
    
    ! Soil ice impedance factor (unitless)
    call ReadScalarParameter(distributed_parameters%e_ice,ncid)

    !-----------------------------------------------------------
    ! SaturatedExcessRunoff !
    !-----------------------------------------------------------
    
    ! Decay factor for fractional saturated area (1/m)
    call ReadScalarParameter(distributed_parameters%fff,ncid)

    !-----------------------------------------------------------
    ! initVertical !
    !-----------------------------------------------------------

    ! exponent for microtopography pdf sigma (unitless)
    call ReadScalarParameter(distributed_parameters%slopebeta,ncid)

    ! max topographic slope for microtopography pdf sigma (unitless)
    call ReadScalarParameter(distributed_parameters%slopemax,ncid)

    !-----------------------------------------------------------
    ! SnowCoverFractionSwensonLawrence2012 !
    !-----------------------------------------------------------

    ! SCA shape parameter
    call ReadScalarParameter(distributed_parameters%n_melt_coef,ncid)
       
    ! Accumulation constant for fractional snow covered area (unitless) 
    call ReadScalarParameter(distributed_parameters%accum_factor,ncid)

    !-----------------------------------------------------------
    ! SnowHydrologyMod !
    !-----------------------------------------------------------

    ! Upper limit on destructive metamorphism compaction (kg/m3)
    call ReadScalarParameter(distributed_parameters%upplim_destruct_metamorph,ncid)

    !-----------------------------------------------------------
    ! SoilStateInitTimeConstMod !
    !-----------------------------------------------------------

    call ReadScalarParameter(distributed_parameters%bsw_sf,ncid)

    call ReadScalarParameter(distributed_parameters%hksat_sf,ncid)

    call ReadScalarParameter(distributed_parameters%sucsat_sf,ncid)

    call ReadScalarParameter(distributed_parameters%watsat_sf,ncid)

    !-----------------------------------------------------------
    ! SurfaceResistanceMod !
    !-----------------------------------------------------------

    ! Dry surface layer parameter (mm)
    call ReadScalarParameter(distributed_parameters%d_max,ncid)

    ! Fraction of saturated soil for moisture value at which DSL initiates (unitless)
    call ReadScalarParameter(distributed_parameters%frac_sat_soil_dsl_init,ncid)

    !-----------------------------------------------------------
    ! SurfaceWaterMod !
    !-----------------------------------------------------------

    ! Threshold probability for surface water (unitless)
    call ReadScalarParameter(distributed_parameters%pc,ncid)

    ! Connectivity exponent for surface water (unitless)
    call ReadScalarParameter(distributed_parameters%mu,ncid)

    !-----------------------------------------------------------
    ! WaterDiagnosticBulkType !
    !-----------------------------------------------------------

    ! Momentum roughness length for soil, glacier, wetland (m) 
    call ReadScalarParameter(distributed_parameters%zlnd,ncid)

    ! Minimum allowed snow effective radius (also cold "fresh snow" value) [microns] 
    call ReadScalarParameter(distributed_parameters%snw_rds_min,ncid)
    
    !-----------------------------------------------------------
    ! atm2lndType !
    !-----------------------------------------------------------

    ! non-glacier all rain temperature (degrees C)
    call ReadScalarParameter(distributed_parameters%precip_repartition_nonglc_all_rain_t_celsius,ncid)
    ! non-glacier all snow temperature (degrees C)
    call ReadScalarParameter(distributed_parameters%precip_repartition_nonglc_all_snow_t_celsius,ncid)
    ! glacier all rain temperature (degrees C)
    call ReadScalarParameter(distributed_parameters%precip_repartition_glc_all_rain_t_celsius,ncid)
    ! glacier all snow temperature (degrees C)
    call ReadScalarParameter(distributed_parameters%precip_repartition_glc_all_snow_t_celsius,ncid)

    ! close parameter file
    call ncd_pio_closefile(ncid)

  end subroutine InitGlobalParameters

  !------------------------------------------------------------------------
  subroutine Init(this)
    !
    ! !ARGUMENTS:
    class(distributed_parameter_type) :: this
    !------------------------------------------------------------------------

    allocate(this%aq_sp_yield_min)
    allocate(this%n_baseflow)
    allocate(this%perched_baseflow_scalar)
    allocate(this%e_ice)
    allocate(this%baseflow_scalar)
    allocate(this%fff)
    allocate(this%slopebeta)
    allocate(this%slopemax)
    allocate(this%n_melt_coef)
    allocate(this%accum_factor)
    allocate(this%upplim_destruct_metamorph)
    allocate(this%pc)
    allocate(this%mu)
    allocate(this%bsw_sf)
    allocate(this%hksat_sf)
    allocate(this%sucsat_sf)
    allocate(this%watsat_sf)
    allocate(this%d_max)
    allocate(this%frac_sat_soil_dsl_init)
    allocate(this%zlnd)
    allocate(this%snw_rds_min)
    allocate(this%precip_repartition_nonglc_all_rain_t_celsius)
    allocate(this%precip_repartition_nonglc_all_snow_t_celsius)
    allocate(this%precip_repartition_glc_all_rain_t_celsius)
    allocate(this%precip_repartition_glc_all_snow_t_celsius)

    this%aq_sp_yield_min%name = 'aq_sp_yield_min'
    this%aq_sp_yield_min%is_distributed = .false.

    this%n_baseflow%name = 'n_baseflow'
    this%n_baseflow%is_distributed = .false.

    this%perched_baseflow_scalar%name = 'perched_baseflow_scalar'
    this%perched_baseflow_scalar%is_distributed = .false.

    this%e_ice%name = 'e_ice'
    this%e_ice%is_distributed = .false.

    this%baseflow_scalar%name = 'baseflow_scalar'
    this%baseflow_scalar%is_distributed = .false.

    this%fff%name = 'fff'
    this%fff%is_distributed = .false.

    this%slopebeta%name = 'slopebeta'
    this%slopebeta%is_distributed = .false.

    this%slopemax%name = 'slopemax'
    this%slopemax%is_distributed = .false.

    this%n_melt_coef%name = 'n_melt_coef'
    this%n_melt_coef%is_distributed = .false.

    this%accum_factor%name = 'accum_factor'
    this%accum_factor%is_distributed = .false.

    this%upplim_destruct_metamorph%name = 'upplim_destruct_metamorph'
    this%upplim_destruct_metamorph%is_distributed = .false.

    this%pc%name = 'pc'
    this%pc%is_distributed = .false.

    this%mu%name = 'mu'
    this%mu%is_distributed = .false.

    this%bsw_sf%name = 'bsw_sf'
    this%bsw_sf%is_distributed = .false.

    this%hksat_sf%name = 'hksat_sf'
    this%hksat_sf%is_distributed = .false.

    this%sucsat_sf%name = 'sucsat_sf'
    this%sucsat_sf%is_distributed = .false.

    this%watsat_sf%name = 'watsat_sf'
    this%watsat_sf%is_distributed = .false.

    this%d_max%name = 'd_max'
    this%d_max%is_distributed = .false.

    this%frac_sat_soil_dsl_init%name = 'frac_sat_soil_dsl_init'
    this%frac_sat_soil_dsl_init%is_distributed = .false.

    this%zlnd%name = 'zlnd'
    this%zlnd%is_distributed = .false.

    this%snw_rds_min%name = 'snw_rds_min'
    this%snw_rds_min%is_distributed = .false.

    this%precip_repartition_nonglc_all_rain_t_celsius%name = &
         'precip_repartition_nonglc_all_rain_t'
    this%precip_repartition_nonglc_all_rain_t_celsius%is_distributed = .false.

    this%precip_repartition_nonglc_all_snow_t_celsius%name = &
         'precip_repartition_nonglc_all_snow_t'
    this%precip_repartition_nonglc_all_snow_t_celsius%is_distributed = .false.

    this%precip_repartition_glc_all_rain_t_celsius%name = &
         'precip_repartition_glc_all_rain_t'
    this%precip_repartition_glc_all_rain_t_celsius%is_distributed = .false.

    this%precip_repartition_glc_all_snow_t_celsius%name = &
         'precip_repartition_glc_all_snow_t'
    this%precip_repartition_glc_all_snow_t_celsius%is_distributed = .false.

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine Clean(this)
    !
    ! !ARGUMENTS:
    class(distributed_parameter_type) :: this
    !------------------------------------------------------------------------

    deallocate(this%aq_sp_yield_min%val)
    deallocate(this%n_baseflow%val)
    deallocate(this%perched_baseflow_scalar%val)
    deallocate(this%e_ice%val)
    deallocate(this%baseflow_scalar%val)
    deallocate(this%fff%val)
    deallocate(this%slopebeta%val)
    deallocate(this%slopemax%val)
    deallocate(this%n_melt_coef%val)
    deallocate(this%accum_factor%val)
    deallocate(this%upplim_destruct_metamorph%val)
    deallocate(this%pc%val)
    deallocate(this%mu%val)
    deallocate(this%bsw_sf%val)
    deallocate(this%hksat_sf%val)
    deallocate(this%sucsat_sf%val)
    deallocate(this%watsat_sf%val)
    deallocate(this%d_max%val)
    deallocate(this%frac_sat_soil_dsl_init%val)
    deallocate(this%zlnd%val)
    deallocate(this%snw_rds_min%val)
    deallocate(this%precip_repartition_nonglc_all_rain_t_celsius%val)
    deallocate(this%precip_repartition_nonglc_all_snow_t_celsius%val)
    deallocate(this%precip_repartition_glc_all_rain_t_celsius%val)
    deallocate(this%precip_repartition_glc_all_snow_t_celsius%val)

  end subroutine Clean

end module Distparamtype
